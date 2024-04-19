using Shapefile, DataFrames, HDF5, H5Zzstd, Statistics, Dates

export readdatasets, getparams

"""readdatasets(selectcells="all"): Read large input datasets from C:/Griddata folder to a NamedTuple."""
function readdatasets(selectcells="all", coincidencesuffix="", gridarea="", use_commercial=false)
    dfgeo = readgeodata(selectcells, gridarea, use_commercial)
    
    gridmap = Dict("" => "", "SE1" => "_SE4", "SE2" => "_SE3", "SE3" => "_SE2", "SE4" => "_SE1")
    gridsuffix = gridmap[gridarea]
    profilesHH, profilesAP, profilesEV = loadHDF5("profile_coincidences$coincidencesuffix$gridsuffix", "profilesHH", "profilesAP", "profilesEV");

    # Store index lookups for hour of day and day of year (assume not a leapyear and first period at midnight)
    # year = len_coin == 52560 ? 2021 : 2020
    year = 2021     # Hardcode year to non-leapyear to avoid leapyear issues (we'll toss away data for day 366 later)
    timeperiods = DateTime(year):Minute(10):DateTime(year+1)-Minute(1)
    hourindex = hour.(timeperiods) .+ 1
    dayindex = dayofyear.(timeperiods)

    return (; dfgeo, profilesHH, profilesAP, profilesEV, hourindex, dayindex, selectcells)
end

"getparams(voltagelevel, transformersize): Read all input parameters into a convenient NamedTuple for passing numerical data around."
function getparams(voltagelevel, transformersize, use_commercial)
    # Control parameters
    Zgrid = [0.3, 0.4, 0.55][voltagelevel]
    RXratio = [0.5, 2, 5][voltagelevel]
    voltageLimit = (lower = [0.9, 0.95, 0.95][voltagelevel], upper = [1.1, 1.05, 1.05][voltagelevel])
    thermallimit = [0.8, 1, 1.2][voltagelevel]
    factor = [1, 1, 1.4][voltagelevel]    # EnergyScaleFactor

    Xgrid = Zgrid / sqrt(1 + RXratio^2)
    Rgrid = Xgrid * RXratio

    # Transformer sizing
    alpha = [1.2, 1.5, 1.8, 2.1, 2.4][transformersize]
    margin = alpha
    poplim1 = 50 #Upper population limit for using higher alpha value
    poplim2 = 500 #Lower population limit for starting to use lower alpha value 
    alphaL = 1.15
    alphaH = 1.65

    # Parameters originally defined in network_model_EV.m 
    Design_Z = 0.65         # Maximum earth loop impedance in Swedish distribution grids
    Design_voltage = (lower=0.92, upper=1.06)    # Design voltage

    AptBuilding_share = 1.25                # https://www.sciencedirect.com/science/article/pii/S0378778811005019#fig0010
    AVGEnergyHouseEH = factor * 18.558 * 1000   # Average annual electricity consumption for house in Sweden with direct electric heating, Zimmermann
    AVGEnergyHouseOH = factor * 8.416 * 1000    #Average annual electricity consumption for house in Sweden without direct electric heating, Zimmermann
    AVGEnergyApt = factor * 2.404 * 1000     # Zimmermann, 6 high

    AVGEnergyHouseFuture = AVGEnergyHouseEH
    AVGEnergyHouseOHFuture = AVGEnergyHouseOH
    AVGEnergyAptFuture = AVGEnergyApt * AptBuilding_share

    # Power demands for households (used for cable sizing)
    PowerDemandHH = 13.8 * 1000 #kW, Fuse size 20 A, most common fuse size 
    PowerDemandAPT = 6.9 * 1000 #kW, Fuse size 10 A

    # Velander coefficients
    k2_House = (0.025 + 0.05) / 2
    k1_House = 0.0003
    k2_Apt = 0.014
    k1_Apt = 0.000264

    TransformerCap = [1600, 1500, 1250, 1125, 1000, 900, 800,700, 600, 500, 400, 315, 200, 150, 100, 70, 50, 30]
    TransformerCost = [234551, 223450, 195272, 179381, 163835, 151090, 134751, 124778, 111211, 101565, 83255, 70501, 53509, 46769, 38446, 34731, 32140, 28647]
        # Transformer cost for additional transformer capacities calculated according to: C=-0,0137x^2 + 153,48x + 24055, where x is the transf capacity in kVA

    # Set power factor to 0.9 and calculate the corresponding reactive power share .
    PowerFactorX = sind(acosd(0.9))/0.9

    c = 0.85    # dimensionless country-specific tripping factor (0.85 for Sweden, typically 0.95 elsewhere?)
    Ufn = 230   # V
    Vn = 400    # V, nominal voltage of LV grid

    # Set data for fuse sizes, cable impedance and transformer impedance.
    Transformer_Impedance = [4.85, 5.15, 6.5,6.7, 7.5, 8.29, 10, 10.5, 12.1, 13, 17.6, 20, 32, 43.9, 65, 89, 130, 195] / 1000               # Ω (after dividing by 1000), 
    CableCapacity = [52, 67, 80, 94, 138, 178, 230]                             # A
    # Table 3 in EN60269, for fuses less than 16A from BS3036 (UK standard)   (1490 through linear interpolation)
    targetFuseRatings = [6, 10, 16, 20, 25, 32, 40, 50, 63, 80, 100, 125, 160, 200, 230] # 250, 315]                    # A
    fuses_tripping_size = [18, 32, 65, 85, 110, 150, 190, 250, 320, 425, 580, 715, 950, 1250, 1490] # 1650, 2200]       # A
    Z_lineR_list = [1.83, 1.15, 0.76, 0.641, 0.32, 0.206, 0.125]                # Ω/km
    Z_lineX_list = [0.0861, 0.0817, 0.0783, 0.0779, 0.0762, 0.0745, 0.0752]     # Ω/km
    targetLineImpedance = [4.18, 2.63, 1.91, 1.47, 0.746, 0.495, 0.324]         # Ω/km

    Z_TransformerR = [0.000635481, 0.000687407, 0.000848114120254821, 0.000987029, 0.001125944, 0.001279981, 0.00148841681507050, 0.00173793, 0.002096557, 0.00251028652976876,0.003434074, 0.00492770793724613, 0.00776114000116266, 0.011329575, 0.0202385770250776,0.028643915, 0.0404771540501553, 0.080327]
    Z_TransformerX = [0.005616973, 0.005935247, 0.00763302708229339, 0.008012082, 0.008391136, 0.009181162, 0.0119073345205640, 0.012122, 0.012337, 0.0125514326488438, 0.018351076, 0.0197108317489845, 0.0310445600046506, 0.042407125, 0.0607157310752329, 0.081303087, 0.121431462150466, 0.167632882]

    ChargePower = 6900     # W (max EV charge power)

    profilenumbers = [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200]
   


    # Return all parameters as a NamedTuple, so they can be accessed by e.g. params.voltageLimit, or by unpacking to variables.
    # (Why export both alpha & margin?)
    return (; voltageLimit, thermallimit, factor, Xgrid, Rgrid, alpha, Design_Z, Design_voltage, margin, use_commercial, profilenumbers,
        AVGEnergyHouseFuture, AVGEnergyHouseOHFuture, AVGEnergyAptFuture, k2_House, k1_House, k2_Apt, k1_Apt, TransformerCap, TransformerCost,
        PowerFactorX, c, Ufn, Vn, Transformer_Impedance, CableCapacity, fuses_tripping_size, poplim1, poplim2, PowerDemandHH, PowerDemandAPT,
        targetFuseRatings, Z_lineR_list, Z_lineX_list, targetLineImpedance, Z_TransformerR, Z_TransformerX, ChargePower, alphaL, alphaH)
end
