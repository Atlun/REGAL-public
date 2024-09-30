# REGAL

A model of the Swedish (residential) low-voltage power grid estimated using public datasets. The model first dimensions
the grid based on historical demand and user-selected 10-minute full-year electricity demand profiles for houses and
apartments, then adds additional demand profiles for electric vehicle charging based on user-selected coincidence
matrices, and finally calculates resulting thermal, current and voltage violations for cables and transformers.

## Initial setup

To minimize run times, the model requires precalculated coincidence matrices based on random combinations of EV charging and household demand profiles. Choose one of the following two options, depending on whether you intend to run the full set of scenarios (BaseDirect, BaseSmart, BaseV2G) or just the BaseDirect scenario.

```julia
setup("full")    # To create data for the full set of scenarios (14 GB of free disk space required)
setup("minimal") # Only creates data for the BaseDirect scenario (1.6 GB of free disk space required)
```

Also, to speed up repeated runs, the setup function will also combine open GIS data from multiple sources and save a single dataframe of geodata used in the model.

## Running: arguments and defaults

Here is the full set of arguments, along with their default values (if omitted).

```julia
df, dfs, dfgeo, dd, hh = runmodel(selectcells="city", coincidence="", hh=0, apt=0, shareEV=1; runs=1, save=true, batchname="")
```

The outputs: `df` is a DataFrame with one row per 1x1 km cell and one column for each variable. Run `names(df)`
to see a list of column names. The next DataFrame `dfs` will contain a summary of all runs, with one run per row.
The third DataFrame `dfgeo` contains cell-level geodata and transformed DeSO-level statistics.
The outputs `dd` and `hh` are daily and hourly matrices of violations, with
days of the year or hours of the day in rows and different types of violations in columns, in this order:
`[voltage_low, voltage_high, transformer, cable, any]`. More specifically, the matrices contain total violations per type
(in a 10-minute period) over a full year over all 1 km cells, averaged over lambda (dimension 2 of the EV profile coincidence matrix).
In other words, the sum of violations in all time steps divided by lambda (50 or 100).

| Argument | Valid values | Comment |
| --- | --- | --- |
| selectcells | "all", "rural", "urban", "city" | population density ranges (respectively) 0-Inf, 0-200, 200-1000, 1000-Inf |
| coincidence | "", "BaseDirect", "BaseSmart", "BaseV2G", "BaseMix", "ReducedDirect", "ReducedSmart", "ReducedV2G" | Empty quotes use the old 50-lambda coincidence |
| hh | 0 - 20 | set 0 to select a random profile of the 20 presets |
| apt | 0 - 15 | set 0 to select a random profile of the 15 presets |
| shareEV | 0.0 - 1.0 | EV market share of all cars |
| runs= | 1 - ?? | number of runs in batch, also forces `hh = apt = 0` (regardless of their given values) |
| save= | true, false | whether to save each run to disk in Matlab format |
| batchname= | any quoted string | placed in the filename to identify separate batch runs so they can be reloaded |

## Examples

Keyword arguments (`save`, `runs` and `batchname`) must be named explicitly but can appear in any order, or not at all.
Other (ordinary) arguments also have default values so they are also optional, but can only omitted in sequence from right to left.

```julia
df, dfs, dfgeo, dd, hh = runmodel("city", "ReducedSmart", 19, 8, 1.0, save=false); # single run (runs=1 and batchname="")
df, dfs, dfgeo, dd, hh = runmodel("urban", "ReducedV2G", 0, 0, 1.0, save=true); # single run with random demand profiles
df, dfs, dfgeo, dd, hh = runmodel("urban", "ReducedV2G"); # exactly the same (other arguments at defaults)
df, dfs, dfgeo, dd, hh = runmodel("all", "BaseDirect", 0, 0, 1.0, save=true, runs=10, batchname="test1"); # batch of 10 runs with random profiles
```

## Reloading saved runs

This syntax needs a bit of polish, but for now this should work...

```julia
df, dd, hh = readrun("city_ReducedSmart", 19, 8) # load the first single run in Examples above 
df, dd, hh = readrun("test1 all_BaseDirect", 1, 14) # loads one of the batch runs in the 10-run batch in Examples
```

Runs are saved to and loaded from the Julia work folder (run `pwd()` to print the current working directory).

To get the input DataFrame `dfgeo`, run `readgeodata(selectcells)` with `selectcells` one of "city", "urban", "rural" or "all" as above.

```julia
dfgeo = readgeodata("city")
```

## Reloading a batch run

This loads the summary DataFrames `dft` and `dfs` and the daily and hourly matrices from a batch run.
In a batch run, `dd` and `hh` get averaged over all runs.

```julia
dft, dfs, dd, hh = readbatch("test1", "all", "BaseDirect") # reload and recalculate dd and hh from the test1 batch run
```

## Plotting

Some examples below. Copy to a file and save as "gridcharts.jl", then run `include("gridcharts.jl")`.

```julia
using Gridmodel, Plots, Dates

plotly()

function plotdays(dayviolations; attr...)
    violmat = daymatrix(dayviolations[:,5])
    heatmap(violmat; yflip=true, yticks=1:31, xticks=1:12, color=:viridis, attr...)
end

function plothours(hourviolations; attr...)
    bar(hourviolations[:,5]; attr...)
end

function plotcellviolations(df; attr...)
    viols = df[:, [:viol_voltagelow, :viol_voltagehigh, :viol_transformer, :viol_cable, :viol_any]] |> Matrix
    plot(sort(viols, dims=1); legend=:outertopright, size=(1200,600), tickfont=16, legendfont=16, linewidth=3,
            labels=["voltagelow" "voltagehigh" "transformer" "cable" "all"], attr...)
end

function histvoltage(df; attr...)
    histogram(df[:, :minVoltage]; attr...)
end
```
