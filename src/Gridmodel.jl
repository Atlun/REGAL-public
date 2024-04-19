module Gridmodel

using ProgressMeter, TimerOutputs, CSV, MAT

export runmodel, readrun, readbatch, daymatrix, DataFrame

include("dataconversion.jl")
include("inputdata.jl")
include("networkmodel.jl")

# In a run like df, dd, hh = runmodel():
#
# df[:, :viol_XX] is expected violations (in a 10-minute period) per cell and year
# dd (or hh) is total expected violations (in a 10-minute period) in all cells per day (or per hour)
#
# With vv = df[:, [:viol_voltagelow, :viol_voltagehigh, :viol_transformer, :viol_cable, :viol_any]] |> Matrix
# sum(vv, dims=1) .== sum(dd, dims=1) .== sum(hh, dims=1)
#
# To get probabilities (not likelihoods!!!) of a violation in a 10-minute period:
# sum(vv, dims=1)/8760/6 .== mean(dd, dims=1)/24/6 .== mean(hh, dims=1)/365/6

function runmodel(selectcells="city", coincidence="", shareEV=1;
                runs=1, iterations_per_cell=100, save=true, seed=0, use_commercial=false, batchname="")
    if isempty(coincidence)
        coincidencesuffix, usegridareas = "", false
    elseif coincidence in ["BaseDirect"]
        coincidencesuffix, usegridareas = "_$coincidence", false
    elseif coincidence in ["BaseV2G", "BaseSmart", "BaseMix"]
        coincidencesuffix, usegridareas = "_$coincidence", true
    end

    batchname = isempty(batchname) ? "" : "$batchname "
    
    println("Loading datasets into memory...")
    @time datasets = readdatasets(selectcells, "_BaseDirect", "", use_commercial)   # Use BaseDirect to get all dfgeo cells (SE1-4)
    ncells = length(datasets.dfgeo.pop)

    seed = (seed > 0) ? seed : rand(1:999999)
    println("\nRandom seed to reproduce this run: $seed\n")

    if runs > 1
        avgbatchviolations_days, avgbatchviolations_hours = zeros(Float64, 365, 5), zeros(Float64, 24, 5)
        dfs = summarize()
        dft = DataFrame(zeros(Float32, ncells, 9), names(dfs))
        dft.min_minVoltage .= Inf
        columns = names(dfs)[3:end]
        Random.seed!(seed)
        seeds = rand(1:999999, runs)

        local df
        for i = 1:runs
            gridtext = usegridareas ? " with separate grid areas" : ""
            println("Run $i/$runs$gridtext:")
            df, dayviols, hourviols = startrun(selectcells, shareEV, coincidencesuffix, usegridareas;
                                        iterations_per_cell, datasets, save, seed=seeds[i], use_commercial, batchname)
            avgbatchviolations_days .+= dayviols
            avgbatchviolations_hours .+= hourviols
            dfs = [dfs; summarize(df, seeds[i])]
            dft[:, 3:end] .+= df[:, columns]
            dft.min_minVoltage .= min.(dft.min_minVoltage, df.minVoltage)
            println()
        end
        avgbatchviolations_days ./= runs
        avgbatchviolations_hours ./= runs
        dft[:, 3:end] ./= runs

        if save
            CSV.write("Gridmodelsummary $batchname$selectcells$coincidencesuffix $seed table.csv", dfs)
            CSV.write("Gridmodelmeancells $batchname$selectcells$coincidencesuffix $seed table.csv", dft)
        end

        return dft, dfs, datasets.dfgeo, avgbatchviolations_days, avgbatchviolations_hours
    else
        df, dayviols, hourviols = startrun(selectcells, shareEV, coincidencesuffix, usegridareas;
                                            iterations_per_cell, datasets, save, seed, use_commercial, batchname)
        dfgeo = readgeodata(selectcells, "", use_commercial)
        return df, summarize(df, seed), dfgeo, dayviols, hourviols
    end
end

function summarize(df=nothing, seed=0)
    columns = [:minVoltage, :max_deltaPower, :viol_voltagelow, :viol_voltagehigh, :viol_transformer, :viol_cable, :viol_any]
    if df === nothing
        df1 = DataFrame(min_minVoltage=Float32[], max_maxPower=Float32[])
        df2 = DataFrame(zeros(Float32, 0, 7), columns)
    else
        populated = (df.Pop_density .> 0)
        df1 = DataFrame(min_minVoltage=minimum(df.minVoltage[populated]), max_maxPower=maximum(df.max_deltaPower[populated]))
        df2 = DataFrame(mean(Matrix(df[populated, columns]), dims=1), columns)
    end
    return [df1 df2]
end

function startrun(selectcells="", shareEV=1, coincidencesuffix="", usegridareas=false;
                    iterations_per_cell=100, datasets=nothing, save=true, seed=0, use_commercial=false, batchname="")

    resultnames = [:Pop_density, :Pop_type, :nBranches, :CustomersCalc, :CustomersInitial, :CustomersPerKm,
                    :CustomersPerTransformer, :NumberOfTransformers, :TrCap, :AVGEnergy, :PowerDemand, :LLVMax,
                    :maxZearth, :meanCable, :TransformerDemand, :meanCableDemand, :minVoltage, :max_deltaPower,
                    :viol_voltagelow, :viol_voltagehigh, :viol_transformer, :viol_cable, :viol_any]

    nresults = length(resultnames)

    if usegridareas
        totalresults = zeros(Float32, nresults, 0)
        cellids = zeros(Int, 0)
        viol_days, viol_hours = zeros(Int, 365, 5), zeros(Int, 24, 5)
        for (i, gridarea) in enumerate(["SE1", "SE2", "SE3", "SE4"])
            println("Loading datasets into memory...")
            datasets = readdatasets(selectcells, coincidencesuffix, gridarea)
            results, vdays, vhours = singlerun(datasets, nresults, shareEV, iterations_per_cell, seed+i-1, use_commercial)
            viol_days .+= vdays
            viol_hours .+= vhours
            totalresults = [totalresults results]
            cellids = [cellids; datasets.dfgeo.cellid]
        end
        cellid_correctorder = readgeodata(selectcells).cellid
        ii = indexin(cellid_correctorder, cellids)
        totalresults = totalresults[:, ii]
    else
        totalresults, viol_days, viol_hours = singlerun(datasets, nresults, shareEV, iterations_per_cell, seed, use_commercial)
    end

    df = DataFrame(totalresults', resultnames)
    avgviolations_days, avgviolations_hours = viol_days/iterations_per_cell, viol_hours/iterations_per_cell
    if save
        savefiles("$batchname$selectcells$coincidencesuffix $seed", df, avgviolations_days, avgviolations_hours)
    end

    return df, avgviolations_days, avgviolations_hours
end

function singlerun(datasets, nresults, shareEV, iterations_per_cell, seed, use_commercial)
    params = getparams(2, 3, use_commercial)
    (; dfgeo) = datasets      # unpack these variables from datasets

    totalcells = length(dfgeo.pop)
    results = zeros(Float32, nresults, totalcells)

    ToTCost = [zeros(18000, length(params.TransformerCap)) for i = 1:Threads.nthreads()]
    violations_days = [zeros(Int, 365, 5) for i = 1:Threads.nthreads()]       # in a single cell only
    violations_hours = [zeros(Int, 24, 5) for i = 1:Threads.nthreads()]         # in a single cell only
    totalviolations_days = [zeros(Int, 365, 5) for i = 1:Threads.nthreads()]  # total over all cells
    totalviolations_hours = [zeros(Int, 24, 5) for i = 1:Threads.nthreads()]    # total over all cells

    Random.seed!(seed)
    random_index_100 = stack(random_indices(100, 6, iterations_per_cell) for i=1:3; dims=3)   # (nprofiles in profile_coincidences, nsegments+1, iterations_per_cell), +1 segment for load from remaining customers at transformer

    ncells = size(dfgeo, 1)
    println("Iterating over $ncells cells, with $iterations_per_cell random combinations of aggregated profiles per cell.")
    progressbar = Progress(ncells, 1)
    @time Threads.@threads for cell in 1:ncells
        thr = Threads.threadid()
        preallocated = (; ToTCost=ToTCost[thr], violations_days=violations_days[thr], violations_hours=violations_hours[thr], random_index_100)
        results[:, cell] = networkmodel!(cell, shareEV, iterations_per_cell, preallocated, datasets, params)
        totalviolations_days[thr] .+= violations_days[thr]
        totalviolations_hours[thr] .+= violations_hours[thr]
        next!(progressbar)
    end

    return results, sum(totalviolations_days), sum(totalviolations_hours)
end

function savefiles(selectcells, df, dayviolations, hourviolations)
    profilename = "$selectcells"
    CSV.write("Gridmodelrun $profilename table.csv", df)
    matopen("Gridmodelrun $profilename timeviolations.mat", "w"; compress=true) do file
        write(file, "dayviolations", dayviolations)
        write(file, "hourviolations", hourviolations)
    end
end

function readrun(selectcells, hh, apt)
    profilename = "$selectcells $(hh)_$apt"
    df = CSV.read("Gridmodelrun $profilename table.csv", DataFrame)
    matopen("Gridmodelrun $profilename timeviolations.mat") do file
        dayviolations = read(file, "dayviolations")
        hourviolations = read(file, "hourviolations")
        return df, dayviolations, hourviolations
    end
end

function readbatch(batchname, selectcells, coincidence)
    batchname = isempty(batchname) ? "" : "$batchname "
    runs = 0
    avgbatchviolations_days, avgbatchviolations_hours = zeros(Float64, 366, 5), zeros(Float64, 24, 5)
    for fname in readdir(".")
        if startswith(fname, "Gridmodelrun $batchname$selectcells") && endswith(fname, "timeviolations.mat")
            runs += 1
            matopen(fname) do file
                dayviolations = read(file, "dayviolations")
                hourviolations = read(file, "hourviolations")
                if runs == 1 && size(dayviolations, 1) == 365
                    avgbatchviolations_days = avgbatchviolations_days[1:365, :]
                end
                avgbatchviolations_days .+= dayviolations
                avgbatchviolations_hours .+= hourviolations
            end
        end
    end
    runs == 0 && error("No files with batchname=$batchname and selectcells=$selectcells found.")
    println("Aggregated results from $runs runs.")
    avgbatchviolations_days ./= runs
    avgbatchviolations_hours ./= runs
    coincidencesuffix = isempty(coincidence) ? "" : "_$coincidence"
    dfs = CSV.read("Gridmodelsummary $batchname$selectcells$coincidencesuffix table.csv", DataFrame)
    dft = CSV.read("Gridmodelmeancells $batchname$selectcells$coincidencesuffix table.csv", DataFrame)
    return dft, dfs, avgbatchviolations_days, avgbatchviolations_hours
end

function daymatrix(dailyviolations::AbstractVector)
    violmat = zeros(31,12)
    dayofyear = 0
    for m in 1:12
        ndays = daysinmonth(Date(2021,m))
        violmat[1:ndays, m] = dailyviolations[dayofyear+1:dayofyear+ndays]
        dayofyear += ndays
    end
    return violmat
end

end
