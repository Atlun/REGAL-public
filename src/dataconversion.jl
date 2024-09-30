using MAT, HDF5, H5Zzstd, Random, StatsBase, Shapefile, Rasters, ArchGDAL, GDAL_jll, Rasters.LookupArrays, Polynomials
import GeoDataFrames as GDF
import GeoFormatTypes as GFT
GI = Shapefile.GeoInterface
AG = ArchGDAL

export readgeodata, savegeodata_DESO, write_cellcompanydata, loadHDF5, GDF, GFT, GI, AG, CSV, Shapefile, Rasters, mean, openmap

"Read 1 km population data into a GeoDataFrame."
function getpop()
    # Geopackage of population in 1 km squares (114800 polygons of vector data):
    # https://www.scb.se/vara-tjanster/oppna-data/oppna-geodata/statistik-pa-rutor/
    # gdf = GDF.read("$DATAFOLDER/TotRut_SweRef.gpkg")              # 534712 polygons and SWEREF99 cell ID ("Ruta")
    # gdf = GDF.read("$DATAFOLDER/totalbefolkning_1km_221231.gpkg")   # 114800 polygons, SWEREF99 cell ID, cellsize (=1000), population
    gdf = GDF.read("$DATAFOLDER/Totalbefolkning_1km_211231.gpkg")   # 114800 polygons, SWEREF99 cell ID, cellsize (=1000), population
    select!(gdf, [:geom, :Ruta, :Pop])
    rename!(gdf, [:geometry, :cellid, :pop])

    # Ruta/cellid refers to the lower left corner of each cell, so use the center to determine parent DeSO & municipality.
    #
    # julia> Gridmodel.Rasters.GeoInterface.coordinates(pp[1000,1])
    # 1-element Vector{Vector{Vector{Float64}}}:
    #  [[749000.0, 7.204e6], [749000.0, 7.205e6], [750000.0, 7.205e6], [750000.0, 7.204e6], [749000.0, 7.204e6]]
    #
    # julia> pp[1000,:].Ruta
    # "7490007204000"
    #
    # julia> Gridmodel.Rasters.GeoInterface.extent(pp[1000,1])
    # Extent(X = (749000.0, 750000.0), Y = (7.204e6, 7.205e6))
end

"""Read and merge GeoDataFrames of 1 km population and DeSO boundaries.

    When a 1 km cell is not inside a DeSO, assign it to DeSO in neighboring cells instead.
    
    Return (dfpop, dfdeso), where dfpop is a GeoDataFrame of 114799 (1 km) cells of population, coordinates and DeSO ID  
    and dfdeso is the original GeoDataFrame of 5984 DeSO polygons, names, ID and parent municipalities."""
function merge_pop_coords_deso()
    swe, dfdeso = getdesoraster()
    dfpop = getpop()
    dfpop.deso .= ""
    for r in eachrow(dfpop)
        x, y = decode_coords(r.cellid)
        ndx_deso = swe[X(At(x)),Y(At(y))]
        if ndx_deso > 0
            r.deso = dfdeso.deso[ndx_deso]
        else
            # 76 cells are not inside a DeSO, but 75 are immediately adjacent to one
            neighbors = [swe[X(At(x-1000)),Y(At(y))], swe[X(At(x+1000)),Y(At(y))],
                            swe[X(At(x)),Y(At(y-1000))], swe[X(At(x)),Y(At(y+1000))]]
            realneighbors = neighbors[neighbors .> 0]
            r.deso = isempty(realneighbors) ? "" : dfdeso.deso[mostfrequentelement(realneighbors)]
        end
    end
    ii = findfirst(isempty.(dfpop.deso))    # one isolated cell left, "Svenska Högarna" outside Stockholm archepelago
    delete!(dfpop, ii)     # No DeSO, so just delete it (only 7 pop)
    # east, north = decode_coords(r.cellid)
    # ArchGDAL.reproject([north, east], G.EPSG(3006), G.EPSG(4326))
    # https://www.google.com/maps/@59.44709028024455,19.50796402393453,12z
    dfpop.cellid .= parse.(Int, dfpop.cellid)   # Convert cellid to integers
    dfpop, dfdeso
end

decode_coords(v::Vector) = permutedims(hcat(collect.(decode_coords.(v))...))
decode_coords(ruta::Int) = decode_coords(string(ruta))
decode_coords(ruta::String) = (parse(Int, ruta[1:6]), parse(Int, ruta[7:end])) .+ 500       # add 500 for center instead of lower left point

"""Rasterize DeSO boundaries using center points of our 1 km grid squares.
    IMPORTANT: Many DeSO areas enclose and (incorrectly) overlap smaller DeSOs. We account for
    this by calculating their areas and rasterizing them in order of decreasing area (so larger
    regions are burned into the resulting raster first and are then overwritten by smaller regions). 
    Return (swe, dfdeso), where swe is a Raster and dfdeso a GeoDataFrame of DeSO areas."""
function getdesoraster()
    # Geopackage of DeSO boundaries (vector data):
    # https://www.scb.se/vara-tjanster/oppna-data/oppna-geodata/deso--demografiska-statistikomraden/
    dfdeso = GDF.read("$DATAFOLDER/DeSO_2018_v2.gpkg")
    dfdeso.area_km2 = Rasters.GeoInterface.area.(dfdeso.geom) / 1e6
    swe = rasterizeSWEREF99(dfdeso, :geom)
    return swe, dfdeso
    # colorindex = repeat(1:8, 748)
    # colors = copy(swe)
    # colors[swe .> 0] .= colorindex[colors[swe .>0]]
    # return swe, colors, dfdeso
end

function getcompanydata()
    dfcompany = GDF.read("$DATAFOLDER/Elnätsområden Therese/omraden.shp")
    dfcompany.bolag[276] = "Hedemorahyttorna"
    dfcompany.snitt = parse.(Int, dfcompany.snitt)
    disallowmissing!(dfcompany)
    return dfcompany
end

"""Rasterize Swedish power company areas using center points of our 1 km grid squares.
    Return (company, dfcompany), where company is a Raster and dfcompany a GeoDataFrame of company areas."""
function getcompanyraster()
    # Shapefile of Swedish power companies (purchased by Therese)
    # Is this the source?   https://www.natomraden.se/
    dfcompany = getcompanydata()
    company = rasterizeSWEREF99(dfcompany, :geometry)
    return company, dfcompany
end

"""Rasterize a SWEREF99 GeoDataFrame using center points of our 1 km grid squares.
    Return (swe, dfdeso), where swe is a Raster and dfdeso a GeoDataFrame of DeSO areas."""
function rasterizeSWEREF99(gdf, column)
    # Build a raster by sampling the vector data every 1000 m in the *center* of our 1 km grid squares.
    # The SWEREF99 cell ID ("Ruta") refers to the lower left corner of each cell, not the center.
    # (outer coordinate limits: x: 181896.33 - 1086312.94, y: 6090353.78 - 7689478.31)
    rasterize_geovector(gdf[!, column], (181500,1086500), (6089500,7689500), 1000, EPSG(3006))
end

"""rasterize_geovector(gv, xlim, ylim, interval, crs; T=UInt32, missingval=UInt32(0))

    Rasterize a geovector `gv` (i.e. the geometry column containing polygons of a GeoDataFrame)
    between coordinate limits `xlim` and `ylim` (i.e. Tuples in format (xmin, xmax)),
    spaced by `interval` and using coordinate system `crs`. Return a Raster."""
function rasterize_geovector(gv, xlim, ylim, interval, crs; T=UInt32, missingval=UInt32(0))
    # x, y = values(ext)
    dimz = X(xlim[1]:interval:xlim[2]; mode=Projected(; sampling=Points(), crs)),
        Y(ylim[1]:interval:ylim[2]; mode=Projected(; sampling=Points(), crs))
    raster = Raster(zeros(T, dimz); missingval)
    for (i, g) in enumerate(gv)
        rasterize!(raster, g, fill=i)
    end
    return raster
end

"Area of the intersection of two polygons"
intersectedarea(p1,p2) = GI.intersects(p1,p2) ? GI.area(GI.intersection(p1,p2)) : 0.0

"Calculate the intersected area matrix of all 1 km cells and DeSO areas (9 minute run time)."
cell_deso_intersection(dfpop, dfdeso, T=Float64) = T[intersectedarea(pop, deso) for pop in dfpop.geometry, deso in dfdeso.geom]

"Save indices and values of nonzero elements as a CSV file (only 0.2% of elements are > 0)."
function save_intersection_matrix()
    dfpop, dfdeso = merge_pop_coords_deso()
    sort!(dfdeso, :deso)
    areamat = cell_deso_intersection(dfpop, dfdeso)
    areamat ./= sum(areamat, dims=2)    # normalize so cell row sum = 1 (ensures even border cells belong 100% to a DeSO)
    ii = findall(areamat .> 0)
    df = DataFrame(:row => getindex.(ii,1), :col => getindex.(ii,2), :val => areamat[ii])
    CSV.write("$DATAFOLDER/intersectionmatrix.csv", df)
end

"Read back the sparse matrix saved by `save_intersection_matrix`."
function read_intersection_matrix(T=Float64)
    df = CSV.read("$DATAFOLDER/intersectionmatrix.csv", DataFrame)
    areamat = zeros(T, 114799, 5984)
    ii = CartesianIndex.(df.row, df.col)
    areamat[ii] .= df.val
    return areamat
end

# getcoords(point) = (GI.getcoord(point,1), GI.getcoord(point,2))
# centercoords(geometry) = getcoords(GI.getpoint(GI.getgeom(geometry, 1), 1)) .+ 500
# centerpoint(geometry) = AG.createpoint(centercoords(geometry))
# centerintersects(p1, p2) = GI.intersects(centerpoint(p1), p2)

function write_cellcompanydata(outputdir="$DATAFOLDER/output")
    dfpop, dfdeso = merge_pop_coords_deso()
    println("Identifying power company for each cell...")
    dfcellcompany = getcellcompanydata(dfpop)
    
    dfcells = hcat(dfpop[!, [1,3,4]], dfcellcompany[!, 2:end])
    rename!(dfcells, Dict("snitt" => "elomrade"))
    mkpath(outputdir)
    shapefilename = "$outputdir/bolagsinfo_celler.shp"
    println("Writing Shapefile to $shapefilename...")
    GDF.write(shapefilename, dfcells; layer_name="celler", crs=EPSG(3006))

    dfcells = hcat(dfpop[!, [3,4,6,7]], dfcellcompany[!, 2:end])
    rename!(dfcells, Dict("snitt" => "elomrade"))
    csvfilename = "$outputdir/bolagsinfo_celler.csv"
    println("Writing CSV file to $csvfilename...")
    CSV.write(csvfilename, dfcells)
    nothing
end

"""aggregate_extents(gv)

    Calculate extent of the geovector `gv` by aggregating extents of its individual polygons."""
function aggregate_extents(gv)
    GI = Shapefile.GeoInterface
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    for g in gv
        x, y = values(GI.extent(g)) 
        xmin = min(xmin, x[1])
        xmax = max(xmax, x[2])
        ymin = min(ymin, y[1])
        ymax = max(ymax, y[2])
    end
    return (X=(xmin,xmax), Y=(ymin,ymax))
end

function gridmap(df, rows, z)
    dd, zz = df[rows, :], z[rows]
    coords = decode_coords(dd.cellid)
    xlim, ylim = extrema(coords[:,1]), extrema(coords[:,2])
    # xx, yy = range(xlim..., step=1000), range(ylim..., step=1000)
    xlen, ylen = coordindex(xlim[2], xlim), coordindex(ylim[2], ylim)
    rr = zeros(ylen, xlen)
    for (i, row) in enumerate(eachrow(dd))
        x, y = coordindex(coords[i,1], xlim), coordindex(coords[i,2], ylim)
        rr[y, x] += zz[i]
    end
    rr[.!isfinite.(rr)] .= 0
    return rr
end

coordindex(x, xlim) = (x - xlim[1]) ÷ 1000 + 1

"""Read statistics given per DeSO area och check compatibility of mutual DeSO row order.

    Return `(cars, persons_per_house, persons_per_contract, households, dwellings)`, where each of these is a DataFrame of those DeSO statistics."""
function readdesostats(year=2021)
    yearstr = string(year)
    # https://www.scb.se/hitta-statistik/regional-statistik-och-kartor/regionala-indelningar/deso---demografiska-statistikomraden/deso-tabellerna-i-ssd--information-och-instruktioner/
    cars0 = CSV.read("$DATAFOLDER/Personbilar 2015-2021.csv", DataFrame)
    persons_per_house0 = CSV.read("$DATAFOLDER/Personer per hustyp 2018-2021.csv", DataFrame)
    persons_per_contract0 = CSV.read("$DATAFOLDER/Personer per upplåtelseform 2018-2021.csv", DataFrame)
    households0 = CSV.read("$DATAFOLDER/Antal hushåll 2018-2021.csv", DataFrame)
    dwellings0 = CSV.read("$DATAFOLDER/Antal lägenheter 2018-2021.csv", DataFrame)

    # number of cars per DeSO in traffic, not in traffic, and total
    cars = unstack(cars0[!, ["region", "status", yearstr]], :status, yearstr)
    rename!(cars, [:deso, :traffic, :dormant, :total])
    disallowmissing!(cars)

    # number of *PERSONS* per type of dwelling (single or multi-family), in each DeSO
    persons_per_house = unstack(persons_per_house0[!, ["region", "hustyp", yearstr]], :hustyp, yearstr)
    rename!(persons_per_house, [:deso, :single, :multi, :other, :special, :missing, :total])
    disallowmissing!(persons_per_house)

    # number of *PERSONS* per type of contract (rental, condo, owner), in each DeSO
    persons_per_contract = unstack(persons_per_contract0[!, ["region", "upplåtelseform", yearstr]], :upplåtelseform, yearstr)
    rename!(persons_per_contract, [:deso, :owner, :tenantowned, :rental, :missing, :total])
    select!(persons_per_contract, [:deso, :rental, :tenantowned, :owner, :missing, :total])
    disallowmissing!(persons_per_contract)

    # number of households per DeSO by type (single/couple, kids/no kids)
    households = unstack(households0[!, ["region", "hushållstyp", yearstr]], :hushållstyp, yearstr)
    rename!(households, [:deso, :couple_kids, :couple_nokids, :single_kids, :single_nokids, :other, :total])
    disallowmissing!(households)

    # number of dwellings per DeSO by ownership type 
    dwellings = unstack(dwellings0[!, ["region", "upplåtelseform", yearstr]], :upplåtelseform, yearstr)
    rename!(dwellings, [:deso, :rental, :tenantowned, :owner, :missing])
    dwellings.total .= sum(Matrix(dwellings[!, 2:5]), dims=2)
    disallowmissing!(dwellings)

    @assert all(cars.deso .== persons_per_house.deso)
    @assert all(cars.deso .== persons_per_contract.deso)
    @assert all(cars.deso .== households.deso)
    @assert all(cars.deso .== dwellings.deso)
    return cars, persons_per_house, persons_per_contract, households, dwellings
end

"Identify power company for each cell and return a large DataFrame with company data."
function getcellcompanydata(dfpop)
    dfcompany = getcompanydata()
    ncells = size(dfpop,1)
    companyindex = zeros(Int, ncells)
    progressbar = Progress(ncells, 1)
    Threads.@threads for i = 1:ncells
        companyindex[i] = argmax(intersectedarea(dfpop.geometry[i], grid) for grid in dfcompany.geometry)
        next!(progressbar)
    end
    return dfcompany[companyindex, :]
end

"Read and merge all DeSO level statistics together into a single DataFrame."
function buildDESOdata()
    println("Merging cells with DeSO areas...")
    dfpop, dfdeso = merge_pop_coords_deso()
    println("Inserting gridarea data...")
    dfcellcompany = getcellcompanydata(dfpop)
    dfpop.gridarea .= dfcellcompany.snitt
    dfpop.cpy .= dfcellcompany.omr
    dfpop.company .= dfcellcompany.bolag
    select!(dfpop, [:cellid, :deso, :gridarea, :cpy, :company, :pop])
    select!(dfdeso, [:deso, :kommunnamn, :lannamn, :area_km2])
    sort!(dfdeso, :deso)

    println("Reading DeSO statistics...")
    cars, pph, ppc, households, dwellings = readdesostats(2021)
    @assert all(dfdeso.deso .== cars.deso)

    # sanity checks:
    # sum(dwellings.missing) / sum(dwellings.total) # 0.00012, so very little dwelling data missing
    # sum(pph.missing) / sum(pph.total)             # 0.02469, but enough missing here to throw off house/apt stats
    # sum(pph.other + pph.special) / sum(pph.total) # 0.031, so more people in weird housing types
    # sum(ppc.missing) / sum(ppc.total)             # 0.02477, very similar but not identical
    # all(pph.total .== ppc.total)                  # true, so totals per DeSO add up
    # qq = households.total ./ dwellings.total      # these should be quite similar
    # extrema(qq)                                   # (0.40, 2.43), so some DeSO differ a lot
    # sum(qq .< 0.8 .|| qq .> 1.25)                 # 101, so only 101/5984 differ more than 20%

    dfdeso.pop .= pph.total
    dfdeso.cars .= cars.traffic
    dfdeso.houses .= dwellings.owner                            # "owner" not necessarily a house
    dfdeso.apts .= dwellings.rental .+ dwellings.tenantowned    # "rental" and "tenantowned" not necessarily apartments

    dfdeso.popdens .= dfdeso.pop ./ dfdeso.area_km2        # some DeSO polygons are wrong! (when large DeSOs enclose smaller ones area is double-counted)
    dfdeso.carspercapita .= dfdeso.cars ./ dfdeso.pop
    dfdeso.carsperHH .= dfdeso.cars ./ dwellings.total
    dfdeso.popperHH .= dfdeso.pop ./ dwellings.total

    missingfix = ppc.total ./ (ppc.total .- ppc.missing)        # adjust for missing pop in ppc table
    dfdeso.popperhouse .= ppc.owner ./ dwellings.owner .* missingfix
    dfdeso.popperapt .= (ppc.rental + ppc.tenantowned) ./ (dwellings.rental + dwellings.tenantowned) .* missingfix
    dfdeso.popperhouse2 .= pph.single ./ dfdeso.houses
    dfdeso.popperapt2 .= pph.multi ./ dfdeso.apts

    # fixes to avoid NaN issues for DeSOs with no apts or houses
    dfdeso.popperapt[dfdeso.apts.==0] .= 0
    dfdeso.popperhouse[dfdeso.houses.==0] .= 0
    dfdeso.popperapt2[dfdeso.apts.==0] .= 0
    dfdeso.popperhouse2[dfdeso.houses.==0] .= 0
    # mean_aptpop = mean(dfdeso.popperapt[dfdeso.popperapt .> 0])
    # mean_housepop = mean(dfdeso.popperhouse[dfdeso.popperhouse .> 0])
    # dfdeso.popperapt[dfdeso.popperapt.==0 .&& dfdeso.apts .> 0] .= mean_aptpop
    # dfdeso.popperhouse[dfdeso.popperhouse.==0 .&& dfdeso.houses .> 0] .= mean_housepop

    dfdeso.fracAP .= dfdeso.apts ./ (dfdeso.houses + dfdeso.apts)
    dfdeso.fracAP_pw .= (ppc.rental + ppc.tenantowned) ./ (ppc.owner + ppc.rental + ppc.tenantowned)
    dfdeso.fracAP_pw2 .= pph.multi ./ (pph.single + pph.multi)

    dfdeso.pop_in_houses = dfdeso.pop .* (1 .- dfdeso.fracAP_pw)
    dfdeso.pop_in_apts = dfdeso.pop .* dfdeso.fracAP_pw

    # more checks:
    # sum(ppc.owner.*missingfix) / sum(dwellings.owner)     # roughly 2.69 persons/house nationally
    # sum(pph.single) / sum(dfdeso.houses)                  # ... but 2.8 persons/house using other stats
    # (sum(ppc.rental.*missingfix) + sum(ppc.tenantowned.*missingfix)) / (sum(dwellings.rental) + sum(dwellings.tenantowned))
    #                                                       # roughly 1.67 persons/apt nationally
    # sum(pph.multi) / sum(dfdeso.apts)                     # ... but 1.4 persons/house using other stats
    #
    # n_apts1 = dfdeso.fracAP.*(dfdeso.apts + dfdeso.houses)
    # n_apts2 = dfdeso.fracAP_pw.*dfdeso.pop./dfdeso.popperapt      # identical except for NaNs

    println("There are also DeSO stats for population sliced by age, gender, birthplace, etc.")
    println("There are also DeSO stats for income and type of work (sliced by sector).")
    println("\nShould we use fracAP_pw or fracAP_pw2?")

    # colnames = [:deso, :kommunnamn, :lannamn, :cars, :pop_in_houses, :pop_in_apts]
    # newnames = [:deso, :munic, :region, :desocars, :desopop_in_houses, :desopop_in_apts]
    # df = innerjoin(dfpop, dfdeso[!, colnames], on=:deso)
    # rename!(df, Dict(c => n for (c,n) in zip(colnames, newnames)))

    df = innerjoin(dfpop, dfdeso[!, [:deso, :kommunnamn, :lannamn, :popperapt, :popperhouse]], on=:deso)
    rename!(df, Dict(:kommunnamn => :munic, :lannamn => :region))

    # Adjust DeSO level stats for cells with multiple DeSOs (maybe fails due to type conversions)
    # areamat = read_intersection_matrix()
    # df.pop .= areamat * dfdeso.pop / 1e6
    # df.cars .= areamat * dfdeso.cars / 1e6
    # df.carspercapita .= areamat * dfdeso.carspercapita / 1e6
    # df.fracAP_pw .= areamat * dfdeso.fracAP_pw / 1e6

    # THESE COMBINED CELL STATS ARE APPROXIMATE ONLY!! (Because they assume each cell is allocated to a single DeSO)
    # calculate and insert column :area_ncells in dfdeso
    gdf = groupby(df, :deso, sort=true)
    gdf_rows = DataFrames.combine(gdf, nrow => :area_ncells, :pop => sum)
    leftjoin!(dfdeso, gdf_rows, on=:deso)
    replace!(dfdeso.area_ncells, missing => 0)
    replace!(dfdeso.pop_sum, missing => 0)
    select!(dfdeso, 1:4, :area_ncells, 5:size(dfdeso,2))
    disallowmissing!(dfdeso)

    # Adjust DeSO stats for population in cells allocated to each DeSO
    popfix = dfdeso.pop_sum ./ dfdeso.pop
    dfdeso.cars .*= popfix
    dfdeso.pop_in_houses .*= popfix
    dfdeso.pop_in_apts .*= popfix

    distribute_fracAP_and_cars_in_cells!(df, dfdeso)
    df.apts .= df.pop_in_apts ./ df.popperapt
    df.houses .= (df.pop .- df.pop_in_apts) ./ df.popperhouse
    df.apts[df.popperapt.==0] .= 0              # avoid NaN issue for DeSOs with no apts or houses
    df.houses[df.popperhouse.==0] .= 0          # avoid NaN issue for DeSOs with no apts or houses
    households = df.apts .+ df.houses
    
    rows = households .== 0 .&& df.pop .> 0                 # cells with missing data (but all in DeSOs with fracAP==1)
    df.apts[rows] .= df.pop[rows] ./ df.popperapt[rows]     # so add back apts data
    households[rows] .= df.apts[rows]
    df.pop_in_apts[rows] .= df.pop[rows]

    df.carsperHH .= df.cars ./ households
    df.fracAP .= df.apts ./ households
    df.fracAP_pw .= df.pop_in_apts ./ df.pop

    select!(df, Not(:pop_in_apts))
    rename!(df, Dict(:popperapt => :popperapt_deso, :popperhouse => :popperhouse_deso))

    return df, dfdeso
end

"""distribute_fracAP_and_cars_in_cells!(dfcell, dfdeso)

Distributes DeSO-level data on fracAP and cars to cells based on population density
using empirical relationships (see gridcharts.jl), specifically a logistic curve (fracAP)
and 2nd-degree polynomial (cars). Writes new columns to dfcell and dfdeso."""
function distribute_fracAP_and_cars_in_cells!(dfcell, dfdeso)
    log_pd = log10.(dfcell.pop)
    fracAP = fracAPformula_pw.(log_pd)     # fraction *persons* in apts or houses, logistic curve fitted to DeSO data
    carspercap = carformula.(log_pd)

    # now interpret these fracAP as weights and match DeSO-level AP to total cell-level AP
    ncells, ndesos = size(dfcell,1), size(dfdeso,1)
    dfcell.cars = zeros(ncells)
    dfcell.pop_in_apts = zeros(ncells)
    dfdeso.totcars = zeros(ndesos)
    dfdeso.totpop_in_apts = zeros(ndesos)
    for desorow in eachrow(dfdeso)
        rows = dfcell.deso .== desorow.deso
        !any(rows) && continue
        cellpop = dfcell.pop[rows]

        aptpop = cellpop .* fracAP[rows]
        desorow.totpop_in_apts = sum(aptpop)
        fudge_apts = (desorow.totpop_in_apts > 0) ? desorow.pop_in_apts / desorow.totpop_in_apts : 1.0
        new_aptpop = aptpop * fudge_apts
        new_aptpop[new_aptpop .> cellpop] .= cellpop[new_aptpop .> cellpop]
        dfcell.pop_in_apts[rows] .= new_aptpop

        cars = cellpop .* carspercap[rows]
        desorow.totcars = sum(cars)
        fudge_cars = (desorow.totcars > 0) ? desorow.cars / desorow.totcars : 1.0
        dfcell.cars[rows] .= cars * fudge_cars
    end
end

# empirical curves fitted to DeSO data, see gridcharts.jl
logistic(x, x0, k) = 1 ./ (1 + exp(-k * (x - x0)))
fracAPformula(log_pd) = log_pd < 1.7 ? 0 : logistic(log_pd, 2.5, 2.75)
fracAPformula_pw(log_pd) = log_pd < 1.7 ? 0 : logistic(log_pd, 2.75, 2.75)

carpoly(x) = Polynomial(0.6696887920332156 - 0.02658599714105161*x - 0.01884416753214118*x^2)
carformula(log_pd) = log_pd < 0 ? carpoly(0.0) : carpoly(log_pd)

function savegeodata_DESO()
    dfgeo, dfdeso = buildDESOdata()
    println("Inserting borehole data...")
    insert_boreholes!(dfgeo)
    println("Inserting district heating data...")
    insert_district_heating_areas!(dfgeo)
    CSV.write("$DATAFOLDER/geodata_DESO.csv", dfgeo)
end

function getverk()
    # Verksamhetsområden
    # Statistiken visar antal anställda, antal arbetsställen, branscher och arealer för verksamhetsområden. Ett verksamhetsområde kan bestå antingen av en koncentration av flera arbetsställen, som tillsammans bildar ett geografiskt avgränsat område, eller ett område som karaktäriseras av storskalig och industriellt präglad markanvändning och har minst ett arbetsställe. Ett verksamhetsområde är minst 3 hektar stort. Exempel på verksamhetsområden är flygplatser, gruvor och täkter, områden med industri, företagsparker och logistik eller områden med högskolor, universitet och vårdinrättningar.
    # https://www.scb.se/vara-tjanster/oppna-data/oppna-geodata/verksamhetsomraden/
    # ogrinfo WFS:"http://geodata.scb.se/geoserver/stat/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=stat%3AVerksamhetsomraden.2015&outputFormat=SHAPE-ZIP&format_options=charset:UTF-8"
    # ogr2ogr -f "ESRI Shapefile" -lco ENCODING=UTF-8 Verksamhetsomraden.shp WFS:"http://geodata.scb.se/geoserver/stat/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=stat%3AVerksamhetsomraden.2015&outputFormat=SHAPE-ZIP&format_options=charset:UTF-8"
    df = Shapefile.Table("$DATAFOLDER/Verksamhetsomraden.shp") |> DataFrame
end

# Geodataportalen
# https://www.geodata.se/geodataportalen/srv/swe/catalog.search
#     SE.PF Energiproduktionsanläggningar (NACE D) 
#     SE.PF Tillverkningsindustrier (NACE C) 

function testdeso(df)
    desos = unique(df.deso)
    totbad = 0
    for d in desos
        ii = (df.deso .== d)
        munis = df.muni_id[ii]
        main = mostfrequentelement(munis)
        totbad += sum(munis .!= main)
    end
    totbad
end

function unionall(geovector)
    len = length(geovector)
    len <= 1 && return geovector
    geom = geovector[1]
    for g in geovector[2:end]
        geom = GI.union(geom, g)
    end
    return geom
end

uniqueindexes(x) = unique(z -> x[z], 1:length(x))

function mostfrequentelement(x)
    d = Dict{eltype(x), Int}()
    for e in x
        d[e] = get(d, e, 0) + 1
    end
    return argmax(d)
end

"Get indexes of the n smallest elements of A"
function argsmallest(A::AbstractArray{T,N}, n::Integer) where {T,N}
    # should someone ask more elements than array size, just sort array
    if n>= length(vec(A))
      ind=collect(1:length(vec(A)))
      ind=sortperm(A[ind])
      return CartesianIndices(A)[ind]
    end
    # otherwise 
    ind=collect(1:n)
    mymax=maximum(A[ind])
    for j=n+1:length(vec(A))
    if A[j]<mymax
     getout=findmax(A[ind])[2]
     ind[getout]=j
     mymax=maximum(A[ind])
    end
    end
    ind=ind[sortperm(A[ind])]
    
    return CartesianIndices(A)[ind]
end

function openmap(ruta::Int; bing=false)
    east, north = decode_coords(ruta)
    lat, lon = AG.reproject((north, east), EPSG(3006), EPSG(4326))
    openmap(lon, lat; bing)
end

function openmap(lon::Real, lat::Real; bing=false)
    if bing
        url = "https://bing.com/maps/default.aspx?cp=$lat~$lon\"&\"lvl=18\"&\"style=a\"&\"sp=point.$(lat)_$(lon)_"
    else
        url = "http://maps.google.com/maps?t=k\"&\"q=loc:$lat+$lon"
    end
    # Extra quotes to avoid errors with special chars ? and &:
    # https://superuser.com/questions/36728/can-i-launch-urls-from-command-line-in-windows
    # https://discourse.julialang.org/t/quoting-special-characters-of-a-url-in-cmd-objects-on-windows/44324
    c = Cmd(`cmd /c start \"\" $url`, windows_verbatim=true)
    run(c)
end

"""readgeodata(selectcells, gridarea):

Returns a DataFrame (114799 x 15), Swedish 1 km cells in rows, important columns: pop, carsperHH, fracAP."""
function readgeodata(selectcells::String, gridarea="", use_commercial=false)
    if isempty(gridarea)
        text = ""
        gridnumber = 0
    else
        text = " for grid area $gridarea"
        gridnumber = parse(Int, gridarea[3])
    end

    println("    loading cell geodata$text...")
    df = CSV.read("$DATAFOLDER/geodata_DESO.csv", DataFrame)
    
    if selectcells == "all"
        lowerpop, upperpop = 1, Inf
    elseif selectcells == "rural"
        lowerpop, upperpop = 1, 200
    elseif selectcells == "urban"
        lowerpop, upperpop = 200, 1000
    elseif selectcells == "city"
        lowerpop, upperpop = 1000, Inf
    elseif selectcells == "fast"
        lowerpop, upperpop = 6000, Inf
    else
        error("Unrecognized selectcells, must be one of [all,rural,urban,city].")
    end
    rows = df.pop .>= lowerpop .&& df.pop .< upperpop
    if gridnumber > 0
        rows = rows .&& df.gridarea .== gridnumber
    end
    dfgeo = df[rows, :]

    return dfgeo
end

"insertcoords!(df): Calculate X & Y coordinates using the polygons in column 1 and insert into the dataframe."
function insertcoords!(df)
    len = nrow(df)
    insertcols!(df, 2, :X => zeros(Int, len), :Y => zeros(Int, len)) 
    for row in eachrow(df)
        coordvect = Shapefile.GeoInterface.coordinates(row.geometry)[1][1]
        coords = reduce(hcat, coordvect)'       # combine into matrix
        row.X, row.Y = round.(Int, mean(coords, dims=1))
    end
    return df
end

"loadHDF5(filename, varnames...): Reads variables [varname1, varname2, ...] of the HDF5 file $DATAFOLDER/[filename].h5."
function loadHDF5(filename, varnames...)
    prefixedvarnames = Tuple("/$x" for x in varnames)
    h5open("$DATAFOLDER/$filename.h5", "r") do file
        read(file, prefixedvarnames...)
    end
end

function saveHDF5(x, filename="coincidence", varname="coincidence")
    h5open("$filename.h5", "w") do file
        file["/$varname", chunk=HDF5.heuristic_chunk(x), filters=ZstdFilter()] = x
    end
end

# Convert the matlab file to Float32 and save as HDF5 using ZSTD compression.
# Saving as Float32 saves half the disk space (and load/save time), permuting dimensions
# so time steps are in the first dimension improves compression which halves disk space again.
function convert_coincidences()
    folder = "D:/Marre/Filer_coincidence_till_niclas"
    dir = readdir(folder)
    targetfolder = "$DATAFOLDER"
    for (i, matfile) in enumerate(dir)
        println("Reading $matfile (file $i of $(length(dir)))...")
        @time x = loadMAT("$folder/$matfile")
        filename, ext = splitext(matfile)
        outfile = "$targetfolder/$(lowercasefirst(filename))"
        println("    ...converting to Float32 and saving to $outfile.h5...")
        @time saveHDF5(Float32.(permutedims(x, (3,2,1))), outfile)
    end
    nothing
end

function convert_demand()
    vars = matread("C:/Griddata/ResidentialLP.mat")
    h5open("C:/Griddata/residentialdemand.h5", "w") do file
        x = permutedims(vars["HHLoadProfile"])
        file["/HHLoadProfile", chunk=HDF5.heuristic_chunk(x), filters=ZstdFilter()] = x
        x = permutedims(vars["APTLoadProfile"])
        file["/APTLoadProfile", chunk=HDF5.heuristic_chunk(x), filters=ZstdFilter()] = x
    end
end

function convert_coincidence_MAT_Float32()
    @time x = loadMAT("coincidence50100kWh")
    @time saveMAT_compressed(Float32.(x), "coincidence50_float32")
    x = nothing
    GC.gc()
    # @time x = loadMAT("CoincidenceLarge100kWh", "coincidence_ED017_CR69")
    # @time saveHDF5(Float32.(permutedims(x, (3,2,1))), "coincidence_full");
    @time x = loadMAT("CoincidenceLarge100kWh", "coincidence_ED017_CR69")
    @time saveMAT_compressed(Float32.(x), "coincidence100_float32")
    nothing
end

function loadMAT(filename, varname="coincidence")
    # coincidence50100kWh, CoincidenceLarge100kWh, coincidence_compressed, coincidence_Float32
    matopen(filename) do file
        vars = keys(file)
        if !(varname in vars)
            varname = first(vars)
        end
        read(file, varname)
    end
end

function saveMAT_compressed(x, filename="coincidence_compressed", varname="coincidence")
    matopen("$DATAFOLDER/$filename.mat", "w"; compress=true) do file
        write(file, varname, x)
    end
end

function saveMAT_float32(x)
    matopen("$DATAFOLDER/coincidence_float32.mat", "w"; compress=true) do file
        write(file, "coincidence", Float32.(x))
    end
end

function saveHDF5_c4(x, filename="coincidence_c4", varname="coincidence")
    h5open("$DATAFOLDER/$filename.h5", "w") do file
        file["/$varname", compress=4] = x
        # file["/coincidence", shuffle=(), deflate=3] = x
        # g = create_group(file, "mygroup")
        # dset = create_dataset(g, "myvector", Float64, (10,))
        # write(dset,rand(10))
    end
end

# b = matread("$DATAFOLDER/BRDChargingAnnualHomeED17CR69100kWh.mat")["BRDChargingAnnualHomeED17CR69"]
# @time saveHDF5(Float32.(b), "BRD_chargedata", "chargedata");

function EV_coincidence_maker()
    Random.seed!(23)
    charge = loadHDF5("BRD_chargedata", "chargedata")   # (52704, 429)
    coincidence = zeros(Float32, 31, 100, 52704)
    average_charge = zeros(Float32, 52704, 1000)
    for k = 9:31
        println(k)
        coincidence[k, :, :] = get_charge_samples!(k, average_charge, charge)'
    end
    return coincidence
end

function get_charge_samples!(k, average_charge, charge)
    for j = 1:1000
        chargingcars = sample(1:429, k, replace=false)      # random sample, select k of the 429 cars for charging
        average_charge[:, j] = sum(charge[:, chargingcars], dims=2) / k
    end
    ndx = sortperm(vec(maximum(average_charge, dims=1)))    # sort the 1000 samples by max charge over the year (save index only)
    return average_charge[:, ndx[1:10:1000]]                # choose every 10th sample, return avg charge for all 10 minute periods
end

function save_profile_aggregations(scenarios=["BaseDirect", "BaseSmart", "BaseV2G", "BaseMix"]; n=100)
    # Runtime 22 + 17 + 51 seconds on laptop
    (; profilenumbers) = getparams(2, 3, false)

    for scenarioEV in scenarios
        regions = (scenarioEV == "BaseDirect") ? [""] : ["_SE1", "_SE2", "_SE3", "_SE4"]
        for region in regions
            if scenarioEV == "BaseDirect"
                cc = loadHDF5("BRD_chargedata", "chargedata")   # (52704, 429)
            elseif scenarioEV == "BaseMix"
                cc = loadHDF5("BRD_chargedata", "chargedata")[1:52560, :]   # (52704, 429)
                cc3 = CSV.read("$DATAFOLDER/EVcharging_percar_kWh$region.csv", DataFrame)[:, 2:end] |> Matrix
                cc3 = repeat(cc3 / 6.9, inner=(6,1))      # normalize EV load to 1.0 (max charge power = 6.9 kW) and expand hourly data to 10 minute periods
                cc[:, 3:3:426] .= cc3[:, 3:3:426]
            else
                filename = (scenarioEV == "BaseSmart") ? "Optimal_EVcharging_percar_kWh" : "V2G_EVcharging_percar_kWh"
                cc = CSV.read("$DATAFOLDER/$filename$region.csv", DataFrame) |> Matrix
                cc = repeat(cc / 6.9, inner=(6,1))      # normalize EV load to 1.0 (max charge power = 6.9 kW) and expand hourly data to 10 minute periods 
            end

            hh, ap = loadHDF5("residentialdemand_randomized", "HHLoadProfile", "APTLoadProfile")
            DemandHH = add_shifted_profiles(hh; periods=12, weeks=1)    # makes 20*(2*12+1)*3 = 1500 profiles > max 1120 CustomersPerTransformer
            DemandAP = add_shifted_profiles(ap; periods=12, weeks=1)    # makes 15*(2*12+1)*3 = 1125 profiles > max 1120 CustomersPerTransformer
            DemandEV = add_shifted_profiles(cc; periods=0, weeks=1)     # makes 426*3 = 1278 profiles > max 962 carsperHH*CustomersPerTransformer

            h5open("$DATAFOLDER/profile_coincidences_$scenarioEV$region.h5", "w") do file
                println("$scenarioEV $region profilesHH...")
                @time x = randomize_profile_aggregations(DemandHH, profilenumbers, n; seed=23)
                file["/profilesHH", chunk=HDF5.heuristic_chunk(x), filters=ZstdFilter()] = x    # total load
                x = 0; GC.gc()
                println("$scenarioEV $region profilesAP...")
                @time x = randomize_profile_aggregations(DemandAP, profilenumbers, n; seed=24)
                file["/profilesAP", chunk=HDF5.heuristic_chunk(x), filters=ZstdFilter()] = x    # total load
                x = 0; GC.gc()
                println("$scenarioEV $region profilesEV...")
                @time x = randomize_profile_aggregations(DemandEV, profilenumbers, n; seed=25)
                file["/profilesEV", chunk=HDF5.heuristic_chunk(x), filters=ZstdFilter()] = x ./ profilenumbers'   # coincidence (so divide by profilenumbers after summing to get mean)
            end
        end
    end
    nothing
end

random_indices(nindices, nsamples, niterations) = stack(sample(1:nindices, nsamples, replace=false) for i = 1:niterations)

function randomize_profile_aggregations(profiles, num_aggregates, n; seed=23)
    Random.seed!(seed)
    ntime, nprof = size(profiles)
    num_aggregates = num_aggregates[num_aggregates .<= nprof]
    len_agg = length(num_aggregates)
    sample_indexes = [random_indices(nprof, agg, n) for agg in num_aggregates]
    aggprofiles = zeros(Float32, ntime, len_agg, n)
    Threads.@threads for i = 1:len_agg
        for (j, samples) in enumerate(eachcol(sample_indexes[i]))
            @views aggprofiles[:, i, j] .= sum(profiles[:, samples], dims=2)
        end
    end
    return aggprofiles
end

"Expand the number of available 10-minute profiles by adding timeshifted profiles (e.g. +/- 2 hours, +/- 1 week)."
function add_shifted_profiles(profilemat; periods, weeks)
    nperiods, nprofiles = size(profilemat)
    combos = [(profile, period, week) for profile = 1:nprofiles for period = -periods:periods for week = -weeks:weeks]
    new_profilemat = zeros(nperiods, length(combos))
    for (i, combo) in enumerate(combos)
        profile, period, week = combo
        shift = period + week * 24 * 7 * 6
        indexes = circshift(1:nperiods, shift)
        new_profilemat[:, i] .= profilemat[indexes, profile]
    end
    return new_profilemat
end

function aggregate_borehole_files()
    # https://www.sgu.se/produkter-och-tjanster/geologiska-data/oppna-data/grundvatten-oppna-data/brunnar/
    # https://resource.sgu.se/dokument/produkter/oppnadata/brunnar-oppnadata-beskrivning.pdf
    # https://apps.sgu.se/kartvisare/kartvisare-brunnar.html
    folder = "C:/Griddata/Brunnar"
    dir = readdir(folder)
    df = CSV.read("$folder/$(dir[1])", DataFrame)
    for file in dir[2:end]
        println("Reading $file...")
        dfnew = CSV.read("$folder/$file", DataFrame)
        df = [df; dfnew]
    end
    badrows = ismissing.(df.n .+ df.e) .|| ismissing.(df.anvandning)
    df = df[.!badrows, [:id, :n, :e, :kommunkod, :lagesnoggrannhet, :fastighet, :ort, :totaldjup, :anvandning]]
    df.n .= Int.(df.n)
    df.e .= Int.(df.e)
    df.cellid .= [string(Int(floor(e, digits=-3))) * string(Int(floor(n, digits=-3))) for (n,e) in zip(df.n, df.e)]
    df.energibrunn .= (df.anvandning .== "ENERGIBRUNN (VÄRME OCH/ELLER KYLA)")
    CSV.write("C:/Griddata/boreholes.csv", df)
end

function insert_boreholes!(dfgeo)
    df = CSV.read("C:/Griddata/boreholes.csv", DataFrame)
    df = df[df.energibrunn, :]
    counts = countmap(df.cellid)
    dfgeo.energibrunnar .= get.(Ref(counts), dfgeo.cellid, 0)
    nothing
end

function get_district_heating_areas()
    # https://www.euroheat.org/resource/paneuropean-thermal-atlas-5.html
    # https://s-eenergies-open-data-euf.hub.arcgis.com/search?categories=d5.1
    # https://s-eenergies-open-data-euf.hub.arcgis.com/maps/b62b8ad79f0e4ae38f032ad6aadb91a0
    # https://www.seenergies.eu/wp-content/uploads/sites/25/2020/04/sEEnergies-WP5_D5.1-Excess_heat_potentials_of_industrial_sites_in_Europe.pdf
    shapefile = "C:/Griddata/D5_1_District_Heating_Areas/D5_1_District_Heating_Areas.shp"
    # df = Shapefile.Table(shapefile) |> DataFrame
    gdf = GDF.read(shapefile)
    gdf_swe = gdf[gdf.CNTR_CODE .== "SE", :]
    gdf_swe.geometry .= GDF.reproject(gdf_swe.geometry, EPSG(4326), EPSG(3006), order=:trad)
    dh_swe = rasterizeSWEREF99(gdf_swe, :geometry)
    return gdf_swe, dh_swe
end

function insert_district_heating_areas!(dfgeo)
    gdf_swe, dh_swe = get_district_heating_areas()
    cellids = dfgeo.cellid
    ndx_DH = similar(cellids)
    for (i, id) in enumerate(cellids)
        x, y = decode_coords(id)
        ndx_DH[i] = dh_swe[X(At(x)),Y(At(y))]
    end
    dfgeo.index_DHarea .= ndx_DH
    nothing
end
