# Matlab always rounds ties upward (e.g. 2.5 -> 3), but Julia rounds to the nearest even number according to the ISO standard.
# So let's do it the Matlab way for now, while comparing exact results to the Matlab model. 
roundMAT(x::AbstractFloat) = round(x, RoundNearestTiesUp)
roundMAT(::Type{T}, x::AbstractFloat) where {T<:Integer} = round(T, x, RoundNearestTiesUp)

velander(e, k1, k2) = k1*e + k2*sqrt(e)

function networkmodel!(cell, shareEV, iterations_per_cell, preallocated, datasets, params)
    zero_preallocated!(preallocated)
    
    # SECTION 1 - Demand estimation
    celldata = get_celldata(cell, shareEV, datasets)
    demandvars = calc_load_profiles(celldata, params)

    # SECTION 2 - Feeder length calculation
    transvars = optimize_transformers(preallocated, demandvars, celldata, params)
    feedervars = calc_longest_feeder(transvars, celldata)

    # SECTION 3 - Feeder sizing, capacity, voltage and tripping criteria (voltage and thermal limits)
    branchvars = calc_branch_vectors(feedervars)
    trippingvars = check_fuse_size_and_tripping_criteria(branchvars, demandvars, transvars, celldata, params)
    cablevars = cable_capacity_earth_impedance(branchvars, trippingvars, transvars, params)
    update_cables_using_voltage_profile!(cablevars, trippingvars, transvars, celldata, params)
    update_cables_using_loop_impedance!(cablevars, params)

    # SECTION 4 - Adding new loads (EV, solar PV etc.)
    violationvars = check_violations!(preallocated, celldata, branchvars, cablevars, transvars, iterations_per_cell, datasets, params)
    (; TransformerDemand, CableDemand, min_voltage, max_deltaPower, violations) = violationvars

    (; NumberOfTransformers, CustomersPerTransformer, TrCap) = transvars
    CustomersPerKm = CustomersPerTransformer / feedervars.LengthLVPerTransformer
    CustomersCalc = CustomersPerTransformer * NumberOfTransformers

    return [celldata.pop, celldata.type, branchvars.p, CustomersCalc, celldata.NumberOfCustomers, CustomersPerKm, CustomersPerTransformer,
        NumberOfTransformers, TrCap, demandvars.AVGEnergy, demandvars.PowerDemand, feedervars.LLVMax, maximum(cablevars.z_earth),
        mean(cablevars.Cable), TransformerDemand, CableDemand, min_voltage, max_deltaPower, violations...]
end

function zero_preallocated!(preallocated)
    (; ToTCost, violations_days, violations_hours) = preallocated
    fill!(ToTCost, 0.0)
    fill!(violations_days, 0)
    fill!(violations_hours, 0)
end

function get_celldata(cell, shareEV, datasets)
    geodata = datasets.dfgeo[cell, :]
    (; pop, carsperHH, fracAP, apts, houses) = geodata
    carsperHH = carsperHH * shareEV

    # type and LLAF_LV are feeder parameters, cable costs are taken from normvardeslistan from EI.se
    type, LLAF_LV, CostPerKmLineLV, CostPerKmLineMV = 
        if pop <= 200
            1, 1.0, 177000.0, 512337.0      # rural (0-200)
        elseif pop <= 1000
            2, 1.1, 540000.0, 887790.0      # town (200-1000)
        else
            3, 1.2, 827000.0, 1140746.0     # city (>1000)
        end

    fracHH = 1 - fracAP
    fuse = 10 * fracAP + 20 * fracHH
    NumberOfCustomers = roundMAT(Int, apts + houses)

    return (; cell, pop, carsperHH, type, LLAF_LV, CostPerKmLineLV, CostPerKmLineMV, fracAP, fracHH,
                NumberOfCustomers, fuse)
end

function calc_load_profiles(celldata, params)
    (; fracAP, fracHH, NumberOfCustomers, pop) = celldata
    (; AVGEnergyAptFuture, AVGEnergyHouseFuture, AVGEnergyHouseOHFuture, margin, k1_House, k2_House, k1_Apt, k2_Apt,
        TransformerCap, poplim1, poplim2, PowerDemandHH, PowerDemandAPT, alphaL, alphaH, use_commercial) = params

    PowerDemand = fracHH .* PowerDemandHH .+ fracAP .* PowerDemandAPT
    AVGEnergy = fracAP * AVGEnergyAptFuture + fracHH * AVGEnergyHouseFuture

    if pop >= 200
        AVGEnergyHouseFuture = AVGEnergyHouseOHFuture #lower household load for highly populated areas due to availability of alternative heating 
    elseif pop >= 100
        AVGEnergyHouseFuture = AVGEnergyHouseFuture-((AVGEnergyHouseFuture-AVGEnergyHouseOHFuture)/(200-100)*(pop-100))
    end


    if pop <= poplim1
        margin = alphaH  #higher alpha value for low population density
    elseif pop <= poplim2
        margin = alphaH+((alphaL-alphaH)/(poplim2-poplim1)*(pop-poplim1)) #transition area between lower and higher population density
    else
        margin = alphaH    #lower alpha value for high population density
    end

    CustomersPerTransformer = NumberOfCustomers
    PeakLoad = margin * (velander(CustomersPerTransformer * fracHH * AVGEnergyHouseFuture, k1_House, k2_House) +
                        velander(CustomersPerTransformer * fracAP * AVGEnergyAptFuture, k1_Apt, k2_Apt))

    n0 = max(1, ceil(Int, PeakLoad / TransformerCap[1]))
    if n0 >= NumberOfCustomers + 1
        n0 -= 1
    end

    return (; PowerDemand, AVGEnergy, n0)
end

function optimize_transformers(preallocated, demandvars, celldata, params)
    (; ToTCost) = preallocated
    (; n0) = demandvars
    (; NumberOfCustomers) = celldata
    (; TransformerCap, Transformer_Impedance, Z_TransformerR, Z_TransformerX) = params
    # Iterates over all possible numbers of transformers that can supply each cell and calculate costs
    # of power lines and transformers. Select the cheapest option.
    # ToTCost = fill(Inf, NumberOfCustomers-n0+1, length(TransformerCap))    # NM: do we need the columns? Just return TrCap from calc_transformer_costs and store it. And eliminate the fill().
    costview = @view ToTCost[1:NumberOfCustomers-n0+1, :]
    costview .= Inf

    for ntf = n0:NumberOfCustomers                      # NM: can't this be reduced quite a lot?
        totalcost, indexTR = calc_transformer_costs(ntf, celldata, params)
        ToTCost[ntf-n0+1, indexTR] = totalcost
    end

    # Select option with lowest cost.
    xt, TrSize = Tuple(argmin(costview))     # Tuple() needed because argmin returns a CartesianIndex which doesn't support multiple outputs: x,y = CartesianIndex(1,2)
    # xt, TrSize = Tuple(argmin(ToTCost))     # Tuple() needed because argmin returns a CartesianIndex which doesn't support multiple outputs: x,y = CartesianIndex(1,2)
    NumberOfTransformers = n0 + xt - 1                  # NM: everything below repeats code in calc_transformer_costs(), consider reuse

    TrCap = TransformerCap[TrSize]                      # NM: Not same as TR_Cap,i in the old paper? See also PeakLoadTR below.
    Z_transformer = Transformer_Impedance[TrSize]
    Transformer_R = Z_TransformerR[TrSize]
    Transformer_X = Z_TransformerX[TrSize]

    A_TS, CustomersPerTransformer, coincidenceTR = calc_transformer_stats(NumberOfTransformers, demandvars, celldata, params)
    while !transformer_voltage_ok(coincidenceTR, CustomersPerTransformer, Transformer_R, Transformer_X, celldata, params)
        NumberOfTransformers += 1
        A_TS, CustomersPerTransformer, coincidenceTR = calc_transformer_stats(NumberOfTransformers, demandvars, celldata, params)
    end

    return (; NumberOfTransformers, CustomersPerTransformer, TrCap, TrSize, A_TS, coincidenceTR, Z_transformer, Transformer_R, Transformer_X)
end

function calc_transformer_costs(ntf, celldata, params)
    (; fracAP, fracHH, CostPerKmLineLV, CostPerKmLineMV, NumberOfCustomers, pop) = celldata
    (; AVGEnergyAptFuture, AVGEnergyHouseFuture, margin, k1_House, k2_House, k1_Apt, k2_Apt,
        TransformerCap, TransformerCost, AVGEnergyHouseOHFuture, poplim1, poplim2, alphaL, alphaH) = params
    # ntf: number of transformers to consider
    CustomersPerTransformer = roundMAT(Int, NumberOfCustomers / ntf)
    if pop >= 200
        AVGEnergyHouseFuture = AVGEnergyHouseOHFuture #lower household load for highly populated areas due to availability of alternative heating 
    elseif pop >= 100
        AVGEnergyHouseFuture = AVGEnergyHouseFuture-((AVGEnergyHouseFuture-AVGEnergyHouseOHFuture)/(200-100)*(pop-100))
    end

    if pop <= poplim1
        margin = alphaH  #higher alpha value for low population density
    elseif pop <= poplim2
        margin = alphaH+((alphaL-alphaH)/(poplim2-poplim1)*(pop-poplim1)) #transition area between lower and higher population density
    else
        margin = alphaL    #lower alpha value for high population density
    end
    PeakLoadTR = margin * (velander(CustomersPerTransformer * fracHH * AVGEnergyHouseFuture, k1_House, k2_House) +
                            velander(CustomersPerTransformer * fracAP * AVGEnergyAptFuture, k1_Apt, k2_Apt))

    indexTR = searchsortedlast(TransformerCap, PeakLoadTR, rev=true, lt=(<=))   # NM: make sure this can't error or create type instability
    indexTR == 0 && return Inf, 1
    TransformerCostArea = TransformerCost[indexTR] * ntf
    
    A_TS = 1 / ntf      # Area per transformer in each 1 km^2 cell
    NumberOfConnectionsTransformer = roundMAT(Int, NumberOfCustomers / ntf)     # NM: So exactly the same as CustomersPerTransformer above? 

    sqrt_NumberOfConnectionsTransformer = sqrt(NumberOfConnectionsTransformer)  # NM: why sometimes float and sometimes rounded???
    sqrt_ntf = sqrt(ntf)

    # Calculate distance between low (d) and medium (dMV) voltage supply points
    d = sqrt(A_TS) / (sqrt_NumberOfConnectionsTransformer + 1) + 0.02          # NM:* But +1 is inside sqrt() in the old paper. Also, what is 0.02? Adjustment factor γ?
    
    dMV = 1 / (sqrt_ntf + 1)   # Length of medium voltage cables
    if isodd(roundMAT(sqrt_NumberOfConnectionsTransformer))    # NM:* Very similar to nMV equations below, follows from the Manhattan geometry. When odd, the transformer blocks a central customer. When even, need to add extra lengths for the first smaller "ring" of customers oround the transformer.
        n = NumberOfConnectionsTransformer - 1
    else
        n = NumberOfConnectionsTransformer + sqrt_NumberOfConnectionsTransformer - 2
    end
    
    if roundMAT(sqrt_ntf) == 1     # NM: equivalently ntf <= 2, assuming n0 != 0 always [check this]
        nMV = 0.5
    elseif isodd(roundMAT(sqrt_ntf))
        nMV = ntf - 1.0
    else
        nMV = ntf + sqrt_ntf - 2.0
    end

    LengthLVPerTransformer = n * d      # NM:* this means that n = NC_Tr,i  but then NC_Tr,i can't also be NumberOfConnectionsTransformer in d=[..] above. Also definitely different n from calc_longest_feeder(). 
    LengthMV = nMV * dMV                     # NM:* But nMV differs from ntf?
    MV_cost = LengthMV * CostPerKmLineMV                        # Medium-voltage costs
    LV_cost = ntf * LengthLVPerTransformer * CostPerKmLineLV    # Low-voltage costs
    
    LineCost = LV_cost + MV_cost

    return TransformerCostArea + LineCost, indexTR
end

function calc_transformer_stats(NumberOfTransformers, demandvars, celldata, params)
    (; AVGEnergy) = demandvars
    (; NumberOfCustomers) = celldata
    (; k1_House, k2_House) = params

    A_TS = 1 / NumberOfTransformers   # Area per transformer in each 1 km^2 cell

    CustomersPerTransformer = roundMAT(Int, NumberOfCustomers / NumberOfTransformers)

    # Coincidence for each transformer.
    coincidenceTR = velander(CustomersPerTransformer * AVGEnergy, k1_House, k2_House) /
                        (CustomersPerTransformer * velander(AVGEnergy, k1_House, k2_House))

    return A_TS, CustomersPerTransformer, coincidenceTR
end

function transformer_voltage_ok(coincidenceTR, CustomersPerTransformer, Transformer_R, Transformer_X, celldata, params)
    (; Vn, Design_voltage) = params
    tfvars = (; coincidenceTR, CustomersPerTransformer, Transformer_R, Transformer_X)
    vdropTR = calc_initial_voltage_drop(tfvars, celldata, params)
    return Vn - vdropTR >= Vn * Design_voltage.lower
end

function calc_longest_feeder(transvars, celldata)
    (; CustomersPerTransformer, A_TS) = transvars
    LLAF_LV = celldata.LLAF_LV
    # Special cases for 1 or 2 customers per transformer
    # d: distance between low (d) and medium (dMV) voltage supply points
    # NM: n, nc: What are they exactly? Roughly number of customers, but how precisely? See below.
    if CustomersPerTransformer == 1
        LLVMax = 0.1                                    
        d, n, nc = LLVMax, 1, 1
    elseif CustomersPerTransformer == 2
        LLVMax = LLAF_LV * sqrt(A_TS) / 4 + 0.02            # NM: 0.02 again, "last meter", från kabelbox till hus
        d, n, nc = LLVMax, 2, 2
    else
        n = roundMAT(Int, sqrt(CustomersPerTransformer))
        nc = roundMAT(Int, CustomersPerTransformer / 4)
        if n == 2                                       # NM: i.e. CustomersPerTransformer in (3, 4, 5, 6), then nc = (1, 1, 1, 2)
            d = LLAF_LV * sqrt(A_TS) / n                # NM: same as d in else clause, except for rounding of sqrt(CustomersPerTransformer)
            LLVMax = d                                  # NM: same as LLVMax in else clause, except for 0.02
            nc = 2                                      # NM: not same as nc = (1, 1, 1, 2) for values of CustomersPerTransformer above
        else
            LLVMax = LLAF_LV * sqrt(A_TS) * (n - 1)/n + 0.02          # NM: 0.02 again
            d = LLAF_LV * sqrt(A_TS) / sqrt(CustomersPerTransformer)    # moved from above, this is where it is relevant    # NM: should sqrt(CustomersPerTransformer) be n instead?
        end
    end
    LengthLVPerTransformer = CustomersPerTransformer * d

    return (; LLVMax, d, n, nc, LengthLVPerTransformer)
end

"calc_branch_vectors(n, nc): create branch vectors m and mc according to this pattern: n = 11 => m = [2, 2, 2, 2, 3], n = 8 => m = [2, 2, 4], n = 5 => m = [2, 3]." 
function calc_branch_vectors(feedervars)
    (; n, nc, d) = feedervars
    # Set number of branches based on number of customers per low-voltage grid. (NM: distribute customers per branch?)
    p = n in (3,4) ? 2 : clamp((n-1) ÷ 2, 1, 5)    # number of branches (?)
    m1, mc1 = roundMAT(Int, (n - 1)/p), roundMAT(Int, (nc - 1)/p)
    m, mc = fill(m1, p), fill(mc1, p)
    m[p], mc[p] = n - (p - 1) * m1, nc - (p - 1) * mc1

    if n == 1
        d_long = m * d * 1000
    else
        d_long = (m .- 1/p) * d * 1000      # NM: subtract 1/p is equivalent to subtracting 1 after multiplying by p. Needed because we want to take total distance between n customers in a row, which is (n-1)*d.  
    end
    # Flip branch vector (so first element represents the branch closest to the transformer).   # NM: put them in order from the start instead
    return (; p, mc=reverse(mc), d_long=reverse(d_long))
end

function check_fuse_size_and_tripping_criteria(branchvars, demandvars, transvars, celldata, params)
    (; p, mc, d_long) = branchvars
    (; targetFuseRatings, targetLineImpedance, CableCapacity) = params

    mp, I_line_fuse = zeros(Int, p), similar(mc)     # Int
    RX_multiplier, coincidenceLine, P_demand = similar(d_long), similar(d_long), similar(d_long)   # Float64   # NM: maybe Float32 later
    L_max = zeros(7, p)
    maxfuse = maximum(targetFuseRatings)
    Z_loop = [transvars.Z_transformer*1000]
    for rt = 1:p
        rr = sum(mc[rt:end])
        mp[rt] = rr             # store cumulative customers in mp (and don't update it with rr anymore)
        _RX_multiplier = 1.0

        _cablevars = check_cable_current_vs_fuse(rr, _RX_multiplier, maxfuse, demandvars, celldata, params)
        _trippingvars = check_tripping_criteria!(_cablevars, Z_loop, d_long[rt], demandvars, celldata, params)

        (; _L_max, _I_line_fuse) = _trippingvars
        _L_max[_L_max .< d_long[rt]] .= Inf                # NM: Clean this up.
        _L_max[CableCapacity .< _I_line_fuse] .= Inf
        ixZloop = argmin(_L_max)
        Z_loopNew = targetLineImpedance[ixZloop] * d_long[rt] # Z_loop(ixZloop,rt);
        push!(Z_loop, Z_loopNew)
        # NM: Store output in vectors here.
        _L_max[.!isfinite.(_L_max)] .= 0
        RX_multiplier[rt], coincidenceLine[rt], P_demand[rt], I_line_fuse[rt], L_max[:, rt] = (_trippingvars...,)
    end
    return (; mp, RX_multiplier, coincidenceLine, P_demand, I_line_fuse, L_max)
end

function check_cable_current_vs_fuse(rr, _RX_multiplier, maxfuse, demandvars, celldata, params)
    _coincidenceLine, _P_demand, _I_line = calc_line_current(rr, demandvars, celldata, params)
    
    # Check if cable current is compatible with fuse size. Otherwise reduce demand through RX_multiplier.
    while _I_line > maxfuse
        rr = roundMAT(Int, rr/2)
        # RX_multiplier is a factor that consider multiple cables laid in parallell if needed.
        # If it is 0.5, then two cables are used, 0.25 four cables etc. 
        _RX_multiplier = _RX_multiplier * 0.5 
        _coincidenceLine, _P_demand, _I_line = calc_line_current(rr, demandvars, celldata, params)
    end
    return (; rr, _RX_multiplier, _coincidenceLine, _P_demand, _I_line)
end

function calc_line_current(rr, demandvars, celldata, params)
    (; AVGEnergy, PowerDemand) = demandvars
    (; k1_House, k2_House) = params
    coincidenceLine = velander(rr * AVGEnergy, k1_House, k2_House) / (rr * velander(AVGEnergy, k1_House, k2_House))
    P_demand = PowerDemand * rr * coincidenceLine       # NM: We don't actually need this most of the time.
    I_line = coincidenceLine * celldata.fuse * rr
    return coincidenceLine, P_demand, I_line
end

function check_tripping_criteria!(_cablevars, Z_loop, _d_long, demandvars, celldata, params)
    (; rr, _RX_multiplier, _coincidenceLine, _P_demand, _I_line) = _cablevars
    sum_Z_loop = sum(Z_loop)
    _I_line_fuse, _L_max = calc_fuse_and_impedance(_I_line, sum_Z_loop, params)
    
    # Check if the tripping critera is fullfilled given length and capacity of cables and fuse ratings.
    # If needed use multiple parallell cables (e.g. change RX_multiplier).
    while all(_L_max .<= _d_long)
        rr = roundMAT(Int, rr/2) 
        _RX_multiplier = _RX_multiplier * 0.5
        _coincidenceLine, _P_demand, _I_line = calc_line_current(rr, demandvars, celldata, params)
        _I_line_fuse, _L_max = calc_fuse_and_impedance(_I_line, sum_Z_loop, params)
    end
    return (; _RX_multiplier, _coincidenceLine, _P_demand, _I_line_fuse, _L_max)
end

function calc_fuse_and_impedance(I_line, sum_Z_loop, params)
    (; targetFuseRatings, fuses_tripping_size, c, Ufn, targetLineImpedance) = params
    index = searchsortedfirst(targetFuseRatings, I_line)    # NM: make sure this can't error or create type instability
    I_line_fuse, Iu = targetFuseRatings[index], fuses_tripping_size[index]
    ZmaxThick = 1000 * c * Ufn / Iu
    L_max = (ZmaxThick - sum_Z_loop) ./ targetLineImpedance
    return I_line_fuse, L_max
end

function cable_capacity_earth_impedance(branchvars, trippingvars, transvars, params)
    (; p, d_long) = branchvars
    (; L_max, I_line_fuse, RX_multiplier) = trippingvars
    (; CableCapacity, Z_lineR_list, Z_lineX_list, targetLineImpedance) = params
    z_earth = zeros(p+1)
    ixCable, Cable = zeros(Int, p), zeros(Int, p)
    R, X = similar(d_long), similar(d_long)
    z_earth[1] = transvars.Z_transformer * 1000
    # Check if cable capacity is sufficient and calcualte earth impedance.
    # Earth impedance is used for validation purposes.
    for br = 1:p
        ix = findfirst(L_max[:, br] .> 0)       # NM: keep Inf and switch to isfinite later, make sure this can't error or create type instability
        cablecapac = CableCapacity[ix]
        if cablecapac < I_line_fuse[br]    # && !(I_line_fuse[br] == 250 || I_line_fuse[br] == 315)     # can't happen, targetFuseRatings is max 230 now
            ix = searchsortedfirst(CableCapacity, I_line_fuse[br])  # NM: make sure this can't error or create type instability
            cablecapac = CableCapacity[ix]
        end
        R[br] = Z_lineR_list[ix] * d_long[br] * RX_multiplier[br]
        X[br] = Z_lineX_list[ix] * d_long[br] * RX_multiplier[br]
        z_earth[br+1] = z_earth[br] + targetLineImpedance[ix] * d_long[br] * RX_multiplier[br]
        ixCable[br], Cable[br] = ix, cablecapac
    end

    # Update cable resistance and reactance if cable capacity has been changed.
    # Elias: we can comment this whole loop out.
    # ix, id = findmax(ixCable)
    # for br = 1:id
    #     Cable[br] = CableCapacity[ix]
    #     R[br] = Z_lineR_list[ix] * d_long[br]       # NM: No RX_multiplier here? Elias: It's complicated.
    #     X[br] = Z_lineX_list[ix] * d_long[br]
    #     z_earth[br+1] = z_earth[br] + targetLineImpedance[ixCable[br]] * d_long[br] * RX_multiplier[br]     # NM: check this, seems different from just above
    # end
    R, X = R/1000, X/1000

    return (; p, R, X, RX_multiplier, d_long, z_earth, ixCable, Cable)
end

function update_cables_using_voltage_profile!(cablevars, trippingvars, transvars, celldata, params)
    (; R, X, RX_multiplier, d_long, z_earth, ixCable, Cable) = cablevars
    (; Vn, Design_voltage, CableCapacity, Z_lineR_list, Z_lineX_list, targetLineImpedance) = params
    len = length(CableCapacity)
    V = calc_voltage_drop(cablevars, trippingvars, transvars, celldata, params)
    # Update the cable capacity if voltage profile is outside limits.
    while minimum(V) < Vn * Design_voltage.lower        # NM: isn't min always at V[end]?
        val, pos = findmin(ixCable)
        if val == len
            pos = argmax(RX_multiplier)
            RX_multiplier[pos] = RX_multiplier[pos] * 0.5
            R[pos] = Z_lineR_list[val] * d_long[pos] * RX_multiplier[pos] / 1000    # NM: Repeated calculations below and in previous function, break out to separate function
            X[pos] = Z_lineX_list[val] * d_long[pos] * RX_multiplier[pos] / 1000
            ixCable[pos] = val
            Cable[pos] = CableCapacity[val]
            z_earth[pos+1] = z_earth[pos] + targetLineImpedance[val] * d_long[pos] * RX_multiplier[pos]
        else
            R[pos] = Z_lineR_list[val+1] * d_long[pos] * RX_multiplier[pos] / 1000
            X[pos] = Z_lineX_list[val+1] * d_long[pos] * RX_multiplier[pos] / 1000
            ixCable[pos] = val+1
            Cable[pos] = CableCapacity[val+1]
            z_earth[pos+1] = z_earth[pos] + targetLineImpedance[val+1] * d_long[pos] * RX_multiplier[pos]
        end
        V = calc_voltage_drop(cablevars, trippingvars, transvars, celldata, params)
    end
end

function calc_voltage_drop(cablevars, trippingvars, transvars, celldata, params)
    (; p, R, X) = cablevars
    (; P_demand) = trippingvars
    (; PowerFactorX, Vn) = params
    V = zeros(p+1)
    vdropTR = calc_initial_voltage_drop(transvars, celldata, params)
    V[1] = Vn - vdropTR
    # Calculate voltage profile along the cables.
    for br = 1:p
        V[br+1] = V[br] - (R[br] + X[br] * PowerFactorX) * P_demand[br] / Vn
    end
    return V
end

"Calculate the initial voltage drop at the transformer."
function calc_initial_voltage_drop(transvars, celldata, params)
    (; coincidenceTR, CustomersPerTransformer, Transformer_R, Transformer_X) = transvars
    (; fuse) = celldata
    (; PowerFactorX, Vn, Ufn) = params
    return coincidenceTR * CustomersPerTransformer * (Transformer_R + Transformer_X * PowerFactorX) * fuse * Ufn * 3 / Vn
end

function update_cables_using_loop_impedance!(cablevars, params)
    (; Design_Z, CableCapacity, Z_lineR_list, Z_lineX_list, targetLineImpedance) = params
    (; R, X, RX_multiplier, d_long, z_earth, ixCable, Cable) = cablevars
    len = length(CableCapacity)
    # Checking that loop impedance is lower than maximum value. Max value is
    # taken from a CIRED paper written by Per Norborg at Vattenfall.
    z_earth .= z_earth / 1000     # Setting correct unit for impedance (per m)
    while maximum(z_earth) > Design_Z
        val, pos = findmin(ixCable)
        if val == len
            pos = argmax(RX_multiplier)
            RX_multiplier[pos] = RX_multiplier[pos] * 0.5
            
            for gg = pos:len
                z_earth[gg+1] = z_earth[gg] + targetLineImpedance[val] * d_long[gg] * RX_multiplier[gg] / 1000
            end
            R[pos] = Z_lineR_list[val] * d_long[pos] * RX_multiplier[pos] / 1000
            X[pos] = Z_lineX_list[val] * d_long[pos] * RX_multiplier[pos] / 1000
            ixCable[pos] = val
            Cable[pos] = CableCapacity[val] * (1 / RX_multiplier[pos])  # Why 1/RX here? Because we only double cables when we exceed the largest one.
        else
            for gg = pos:length(ixCable)
                z_earth[gg+1] = z_earth[gg] + targetLineImpedance[val+1] * d_long[gg] * RX_multiplier[gg] / 1000 
            end
            R[pos] = Z_lineR_list[val+1] * d_long[pos] * RX_multiplier[pos] / 1000
            X[pos] = Z_lineX_list[val+1] * d_long[pos] * RX_multiplier[pos] / 1000
            ixCable[pos] = val+1
            Cable[pos] = CableCapacity[val+1]                           # correct, no 1/RX here
        end
    end
end

function check_violations!(preallocated, celldata, branchvars, cablevars, transvars, iterations_per_cell, datasets, params)
    (; random_index_100, violations_days, violations_hours) = preallocated
    (; fracHH, fracAP, carsperHH) = celldata
    (; mc) = branchvars
    (; R, X, RX_multiplier, Cable) = cablevars
    (; CustomersPerTransformer, Transformer_R, Transformer_X, TrCap) = transvars
    (; profilesHH, profilesAP, profilesEV, hourindex, dayindex) = datasets
    (; ChargePower, Vn, PowerFactorX, voltageLimit, thermallimit, profilenumbers) = params

    len_t, n_segments = size(profilesHH, 1), length(mc)                         # throw away EV coincidences for day 366

    remaining_customers = CustomersPerTransformer - sum(mc)                     # customers in other branches that connect to the same transformer
    push!(mc, remaining_customers)                                              # add remaining_customers to index 6 of mc (to reuse some of the code instead of duplicating it for extra customers)

    EVs_per_segment = max.(1, ceil.(Int, carsperHH * mc))                       # (5,) NMnew: switch mp -> mc (cumulative customers per segment -> customers)
    evchargepower = EVs_per_segment * ChargePower                               # (5,)

    len_p = length(profilenumbers)
    profile_index = min.(len_p, searchsortedfirst.(Ref(profilenumbers), mc))    # for both HH and AP (weighted with fracAP later)
    profile_indexEV = min.(len_p, searchsortedfirst.(Ref(profilenumbers), EVs_per_segment))
    demandfudge = mc ./ profilenumbers[profile_index]

    # Precalculate coefficients for the inner loop.                                  
    vd0_TR, vd0 = Transformer_X*PowerFactorX/Vn, Transformer_R/Vn               #       NMnew: removed CustomersPerTransformer from vd0 since we now use total loads TL: added PowerFactorX
    vd_line, vd = X * PowerFactorX/Vn, R/Vn                                     # (5,)  NMnew: removed mp here for the same reason
    dcc = RX_multiplier / Vn / sqrt(3)                                          # (5,)  NMnew: mpRX -> RX_multiplier, same reason
    thermalcapTR = TrCap * 1000 * thermallimit
    lowerVoltage, upperVoltage = Vn * voltageLimit.lower, Vn * voltageLimit.upper

    max_deltaPower, min_voltage, max_deltaCurrentCable = -Inf, Inf, fill(-Inf, n_segments)

    @fastmath @inbounds for i = 1:iterations_per_cell
        for t = 1:len_t
            voltagedrop = lineload = EVload = load = 0.0
            violation_cable = false
            # Most program time spent in this "hot inner loop". Do everything to make it fast!
            for s = n_segments:-1:1                                             # reverse segment order to accumulate loads & voltage drops towards transformer
                iHH, iAP, iEV = @view random_index_100[s, i, :]                 # prerandomized indexes 1..100
                loadHH = profilesHH[t, profile_index[s], iHH] * demandfudge[s]
                loadAP = profilesAP[t, profile_index[s], iAP] * demandfudge[s]
                loadEV = profilesEV[t, profile_indexEV[s], iEV] * evchargepower[s]

                lineload += fracHH*loadHH + fracAP*loadAP                       # Accumulate lineload & EVload separately along the segments
                EVload += loadEV
                load = lineload + EVload

                # Calculate voltage drop over the cable (EV charging is a pure resistive load).
                voltagedrop += vd_line[s] * lineload + vd[s] * load

                # Calculate power demand in cable.
                deltaCurrentCable = dcc[s] * load                               # NM: abs() in orig code, but not needed since nothing is negative or complex.
                max_deltaCurrentCable[s] = max(max_deltaCurrentCable[s], deltaCurrentCable)
                violation_cable = violation_cable || deltaCurrentCable > Cable[s]
            end

            # finally get random profiles for load from the remaining customers of the transformer, stored in index 6 (end)
            iHH, iAP, iEV = @view random_index_100[end, i, :]                 # prerandomized indexes 1..100
            loadHH = profilesHH[t, profile_index[end], iHH] * demandfudge[end]
            loadAP = profilesAP[t, profile_index[end], iAP] * demandfudge[end]
            loadEV = profilesEV[t, profile_indexEV[end], iEV] * evchargepower[end]

            lineload += fracHH*loadHH + fracAP*loadAP
            EVload += loadEV
            load = lineload + EVload

            V0 = Vn - vd0_TR * lineload - vd0 * load
            voltage = V0 - voltagedrop

            # Calculate power demand on transformer.
            max_deltaPower = max(max_deltaPower, load)
            violation_transformer = load > thermalcapTR

            min_voltage = min(min_voltage, voltage)

            violation_voltage_low = voltage < lowerVoltage
            violation_voltage_high = voltage > upperVoltage
            violation_any = violation_voltage_low || violation_voltage_high || violation_transformer || violation_cable

            allviolations = (violation_voltage_low, violation_voltage_high, violation_transformer, violation_cable, violation_any)
            viol_days = @view violations_days[dayindex[t], :]
            viol_hours = @view violations_hours[hourindex[t], :]
            viol_days .+= allviolations
            viol_hours .+= allviolations
        end
    end
    violations = dropdims(sum(violations_hours, dims=1), dims=1) / iterations_per_cell

    # Calculate how much extra power the violations require (e.g. what is needed to avoid them), and how long duration
    # this occurs over (in 10 min blocks). This can be used to identify reinforcements. VOLTAGE VIOLATIONS ARE NOT COVERED.
    TransformerDemand = TrCap * max(0, max_deltaPower / thermalcapTR - 1)               # kVA (scalar)
    CableDemand = Vn * sqrt(3) * Cable .* max.(0, max_deltaCurrentCable ./ Cable .- 1)  # W (5,)

    return (; violations, TransformerDemand, CableDemand=mean(CableDemand), min_voltage, max_deltaPower=max_deltaPower/thermalcapTR)
end
