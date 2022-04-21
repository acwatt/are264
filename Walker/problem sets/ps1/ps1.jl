#==============================================================================
file: ps1.jl
description: work through questions in PS1 of ARE 264 B
focus: PS1, ARE 264 (part B) course, spring 2022
instructor: Reed Walker, UC Berkeley
author: Aaron C Watt (UCB Grad Student, Ag & Resource Econ)
notes:
    - Setting up structure of problem set
==============================================================================#

#==============================================================================
                                    PACKAGES
==============================================================================#
using Pkg
# Pkg.add(["Plots", "ZipFile", "StatFiles", "GLM", "DataFrames", "Chain", "CSV", "StatsModels", "Econometrics", "Latexify", "FixedEffectModels", "CUDA", "RDatasets", "StatsPlots", "LinearAlgebra", "TimeSeries"])
using ZipFile  # unzip compressed .zip folders
using StatFiles  # read Stata files
using DataFrames, Chain, DataFramesMeta
using CSV
# using Econometrics
using CovarianceMatrices
using CategoricalArrays
using Statistics
using StatsModels
using StatsBase:sample
using CUDA  # must be run before FixedEffectModels
using FixedEffectModels
using RegressionTables
using Dates
using Plots
using Plots:plot,plot!  # just to get rid of wavey syntax lines (VCcode bug)
using GLM
using Latexify
using RDatasets
using StatsPlots
using LinearAlgebra
using Distributions
using TimeSeries
Latexify.set_default(fmt = "%.3f")



#==============================================================================
                                    GLOBALS
==============================================================================#
# Set local directory paths and filenames
root = dirname(@__FILE__)
zip_fn = "Walker-ProblemSet1-Data.zip"
county_fn = "fips1001.dta"
employment_fn = "reis_combine.dta"
example_fn = "CountyAnnualTemperature1950to2012.dta"
pollution_fn = "poll7080.dta"
example_url = "https://www.dropbox.com/s/fnl1u0ix4e493vv/CountyAnnualTemperature1950to2012.dta?dl=0"



#==============================================================================
KEYBOARD SHORTCUTS
Fold current function: ctrl + shift + [
Fold all: ctrl+k ctrl+0
==============================================================================#


#==============================================================================
ADMIN FUNCTIONS
==============================================================================#
"""Extracts file_name from zip folder zip_name in directory root_dir."""
function extract_file_from_zip(root_dir, zip_name, file_name)
    zip_path = joinpath(root_dir, zip_name)
    save_path = joinpath(root_dir, file_name)
    # If the file is already extracted, done.
    if isfile(save_path)
        return
    end

    println("Extracting $file_name from $root_dir")
    # Open zip file
    zarchive = ZipFile.Reader(zip_path)
    # Find file_name in zip archive
    files = [f.name for f in zarchive.files]
    idx = findfirst(x -> x == file_name, files)
    # Save file_name to disk
    write(save_path, read(zarchive.files[idx]))
end


"""Return dataframe, after unzipping and saving the file."""
function df_from_zip(root_dir, zip_name, filename)
    extract_file_from_zip(root_dir, zip_name, filename)
    println("Loading $filename from file.")
    return DataFrame(load(joinpath(root_dir, filename)))
end


"""Read file at save_path to dataframe, if not present, download from url."""
function df_from_url(url, save_path)
    if isfile(save_path)
        println("Loading $save_path from file.")
        return DataFrame(load(joinpath(root, example_fn)))
    else
        println("Downloading $save_path from url.")
        download(url, joinpath(root, example_fn))
        println("Loading $save_path from file.")
        return DataFrame(load(joinpath(root, example_fn)))
    end
end


"""Return a dataframe, with types Float64 and no missing for given columns."""
function convert_float32(df, cols)
    df1 = dropmissing(df, cols=cols; disallowmissing=true)
    df1 = transform!(df1, cols .=> ByRow(Float64), renamecols=false)
    return df1
end



#==============================================================================
TEMPERATURE FUNCTIONS
==============================================================================#
"""Degree day calculation from Snyder 1985, using integral of sine."""
function degree_day_calc(MAX, MIN, THR)
    M = (MAX + MIN) / 2  # average temp
    W = (MAX - MIN) / 2  # half temp range
    θ = asin((THR - M) / W)  # radian between -π/2 and π/2 where sine crosses THR
    DD = ((M - THR)*(π / 2 - θ) + W * cos(θ)) / π
    return DD
end


"""Return degree days above threshold temperature for single day based on max and min temperatures."""
function degree_day_single_day(Tmax, Tmin, Threshold)
    if Threshold ≥ Tmax
        return 0
    elseif Threshold ≤ Tmin
        return (Tmax + Tmin) / 2 - Threshold
    else
        return degree_day_calc(Tmax, Tmin, Threshold)
    end

end


"""Return dataframe with degree day columns, degree days over each threshold in t."""
function add_degree_days(df::DataFrame, t::Vector)
    println("Adding degree days above thresholds in t.")
    for k ∈ 1:length(t)
        df[!, "dday$(t[k])C"] = degree_day_single_day.(df[!, :tMax], df[!, :tMin], t[k])
    end
    return df
end


"""Return temperature threshold after converting to -1,1 interval (for sine function)"""
function convert_temp_threshold(Tmax, Tmin, Thresh)
    proportion = (Thresh - Tmin) / (Tmax - Tmin)
    return (2 * proportion - 1)
end


"""Return portion (0 to 1) of day's temperature above temperature Threshold,
    given max and min temperature, 
    based on sinusoidal curve to model diurnal temperature cycle.
"""
function fraction_above_threshold(Tmax, Tmin, Threshold)
    α = convert_temp_threshold(Tmax, Tmin, Threshold)
    if α ≥ 1
        return 0
    elseif α ≤ -1
        return 1
    else
        return (π - 2 * asin(α)) / (2π)
    end
end


"""Return portion (0 to 1) of day's temperature between temperature Thresholds,
    given max and min temperature, 
    based on sinusoidal curve to model diurnal temperature cycle.
"""
function fraction_between_thresholds(Tmax, Tmin, Threshold_lower, Threshold_upper)
    # Fraction between = (fraction above lower) - (fraction above upper)
    lower = fraction_above_threshold(Tmax, Tmin, Threshold_lower)
    upper = fraction_above_threshold(Tmax, Tmin, Threshold_upper)
    return lower - upper
end


"""Return dataframe with binned temperature variables columns (indicator for tAvg in bin)."""
function add_bin_indicator_days(df::DataFrame, edges::Vector)
    println("Adding inidicator for average daily temperature in each bin")
    for k ∈ 1:length(edges)
        if k == 1
            # tAvg below lowest threshold
            df[!, "tempB$(edges[k])"] = (df[!, :tAvg] .< edges[k])
        else
            #  tAvg between this threshold and previous threshold
            df[!, "temp$(edges[k-1])to$(edges[k])"] = (edges[k-1] .≤ df[!, :tAvg] .< edges[k])

            if k == length(edges)
                # tAvg above threshold
                df[!, "tempA$(edges[k])"] = (df[!, :tAvg] .≥ edges[k])
            end
        end
    end
    return df
end


"""Return dataframe with binned temperature variables columns (portion of day in bin)."""
function add_bin_portion_days(df::DataFrame, bins::Vector)
    println("Adding portion of days between bin edges in bins.")
    for k ∈ 1:length(bins)
        if k == 1
            # Portion of day below lowest threshold
            df[!, "tempPropB$(bins[k])"] = 1 .- fraction_above_threshold.(df[!, :tMax], df[!, :tMin], bins[k])
        else
            # Portion of day between this threshold and previous threshold
            df[!, "tempProp$(bins[k-1])to$(bins[k])"] = fraction_between_thresholds.(df[!, :tMax], df[!, :tMin], bins[k-1], bins[k])

            if k == length(bins)
                # Portion of day above threshold
                df[!, "tempPropA$(bins[k])"] = fraction_above_threshold.(df[!, :tMax], df[!, :tMin], bins[k])
            end
        end
    end
    return df
end


"""(⋅)₊ function from Harrell 2001"""
m(a) = maximum([0, a])


"""Return jth restricted cubic spline variable using the list of knots t.
    From Frank Harrell's Stats Textbook: Regression Modeling Strategies, 2001, sect. 2.4.4
"""
function cubic_spline_var_harrell(x, j, t)
    k = length(t)
    if (j < 1) | (j > k - 1)
        println("Out of Bounds (j=$j): j is < 1 or > k-1. Not a valid variable.")
    elseif j == 1
        return x
    else
        xⱼ = m.(x .- t[j-1]) .^ 3 .-
             m.(x .- t[k-1]) .^ 3 * (t[k] - t[j-1]) / (t[k] - t[k-1]) +
             m.(x .- t[k]) .^ 3 * (t[k-1] - t[j-1]) / (t[k] - t[k-1])
        # Divide by the scale factor from the stata formula (range of the knots)
        return xⱼ / (t[k] - t[1])^2
    end
end


"""Return ith restricted cubic spline variable using the list of knots k.
    From stata's mkspline documentation
"""
function cubic_spline_var_stata(V, i, k)
    n = length(k)
    if (i < 1) | (i > n - 1)
        println("Out of Bounds (i=$i, n=$n): i < 1 or i > n-1. Not a valid variable.")
    elseif i == 1
        return V
    else
        Vᵢ = (m.(V .- k[i-1]) .^ 3 .- (k[n] - k[n-1])^(-1) *
                                      (m.(V .- k[n-1]) .^ 3 * (k[n] - k[i-1]) .- (V .- k[n]) .^ 3 * (k[n-1] - k[i-1]))) ./
             (k[n] - k[1])^2
        return Vᵢ
    end
end


"""Add restricted cubic spline basis variables to dataframe df.
    Using continuous variable varname from dataframe and list of knots.

    After trying both the stata and Harrell formulas, neither matched the example file given. After comparing the formulas, I recognized the only difference was a scaling factor of 
        `julia 1/(knots[last] - knots[first])^2`
    so I added that scaling factor to the Harrell formula and it produces results that match the example file. So I probably just messed up the stata formula somehow. I assume the example file was created in stata with mksplines.
"""
function add_cubic_spline_vars(df, varname, knots; newbasename="splineC")
    for i in 1:(length(knots)-1)
        newname = "$newbasename$i"
        println("Making $newname")
        df[!, newname] = cubic_spline_var_harrell(df[!, varname], i, knots)
    end
    return df
end


"""Return vector of element-wise minimum between a vector and constant."""
min_vecs(vec, c) = minimum.(zip(vec, repeat([c], length(vec))))
"""Return vector of element-wise maximum between a vector and constant."""
max_vecs(vec, c) = maximum.(zip(vec, repeat([c], length(vec))))


"""Return ith linear piecewise variable using the list of knots k.
    From stata's mkspline documentation
"""
function piecewise_linear_var_stata(V, i, k)
    n = length(k) + 1
    if i == 1
        return min_vecs(V, k[i])
    elseif i > 1 && i < n
        return max_vecs(min_vecs(V, k[i]), k[i-1]) .- k[i-1]
    elseif i == n
        return max_vecs(V, k[n-1]) .- k[n-1]
    else
        println("Out of Bounds (i=$i, n=$n): i < 1 or i > n-1. Not a valid variable.")
    end
end


"""Return dataframe with piecewise linear basis variable columns using given knots."""
function add_piecewise_linear_vars(df, varname, knots; newbasename="piece")
    println("Adding piecewise linear basis variable columns with breakpoints.")
    n = length(knots)+1
    for i in 1:n
        e = i == 1 ? "B$(knots[i])" :
            i == n ? "A$(knots[n-1])" :
            #= o.w =# "$(knots[i-1])to$(knots[i])"
        newname = "$newbasename$e"
        println("Making $newname")
        df[!, newname] = piecewise_linear_var_stata(df[!, varname], i, knots)
    end
    return df
end


"""Return grid-day data aggregated to county-year level"""
function aggregate_df(df, tags)
    # List of columns to average over
    avg_list = [:tMin, :tMax, :tAvg, :prec]
    # List of columns to sum over the year (all cols with a substring from tags)
    sum_list = [n for n in names(df) if any(occursin.(tags, n))]
    # sum_list = [:splineC1, :splineC2, :splineC3, :splineC4]

    # Sum over all days in each year, for each grid point
    df[!, "date"] = Date(1960, 1, 1) + Day.(df.dateNum)
    df[!, "year"] = Year.(df.date)
    gd_gy = groupby(df, [:gridNumber, :year])
    df_gy = combine(gd_gy,
                    avg_list .=> mean .=> avg_list,
                    sum_list .=> sum .=> sum_list)
    # Take the average over all grid points in the county
    gd_y = groupby(df_gy, :year)
    df_y = combine(gd_y, valuecols(gd_y) .=> mean .=> valuecols(gd_y))

    return df_y
end


"""Save temperature variables
    1.1.1 Construct 4 temperature response variables
    1.1.2 Aggregate to year-county
    - sum each new temperaature variable in each grid point over each year
    - average tMin, tMax, tAvg, perc in each gridpoint over each year
    - take a simple average over all grid points in the county to get county-year level
"""
function create_temperature_vars()
    # Extract county fips file from zip if not present, and read
    extract_file_from_zip(root, zip_fn, county_fn)

    # Load data and create average daily temp
    df_temp1 = DataFrame(load(joinpath(root, county_fn)))
    df_temp1[!, :tAvg] = (df_temp1[!, :tMax] .+ df_temp1[!, :tMin]) ./ 2

    # Add degree days for each day-gridpoint
    thresholds = [30, 32, 34]
    df_temp1 = add_degree_days(df_temp1, thresholds)

    # Add binned temperature variables
    bin_edges = [0, 4, 8, 12, 16, 20, 24, 28, 32]
    df_temp1 = add_bin_indicator_days(df_temp1, bin_edges)
    df_temp1 = add_bin_portion_days(df_temp1, bin_edges)

    # Add cubic spline basis variables for each day-gridpoint
    knots = [0 8 16 24 32]
    df_temp1 = add_cubic_spline_vars(df_temp1, :tAvg, knots)

    # Add piecewise linear basis variables for each day-gridpoint
    knots = [28, 32]
    df_temp1 = add_piecewise_linear_vars(df_temp1, :tAvg, knots)

    # sum over year, then average over all grid points in county
    column_tags = ["dday", "temp", "splineC", "piece"]
    df_temp = aggregate_df(df_temp1, column_tags)

    # Save the new file
    CSV.write(joinpath(root, "CountyAnnualTemperature_aaron.csv"), df_temp)
end








#==============================================================================
ANALYSIS FUNCTIONS
==============================================================================#
function apply_vcov(reg::StatsModels.TableRegressionModel, SE::AbstractVector)
    df = DataFrame(coeftable(reg))[!, ["Name", "Coef.", "Std. Error"]]
    @assert length(SE) == nrow(df)
    df[!, "Std. Error"] = SE
    return df
end

"""Return OLS or WLS regression using dataframe df and formula.

    Weights specified by wts can be string (column name in df) or vector of weights
    of same length as df. If wts is unspecified, OLS is assumed.

    # Arguments:
    - `df::DataFrame`: dataframe of data
    - `formula::FormulaTerm`: @formula with names of columns in df
    - `weights::Union{String, Vector, Nothing}=nothing`: string columnname of df, or vector, used as weights in WLS
    - `vcov::Union{String, Nothing}=nothing`: a string indicating which Heteroskedasticity-robust covariance matrix to report and use in standard errors.
    - `cluster::Union{String, Nothing}=nothing`: string columnname of df to cluster on
    Long & Ervin (2000) conduct a simulation study of HC estimators (HC0 to HC3) 
    in the linear regression model, recommending to use HC3 which is thus the default

    Currently only works if there are no missing values in df.
"""
function regression(df::DataFrame,
                    formula::FormulaTerm; 
                    weights::Union{String, Vector, Nothing}=nothing,
                    varcov::String="HC3",
                    cluster::Union{String, Nothing}=nothing,
                    filestub::Union{String, Nothing}=nothing,
                    append_time::Bool=true)
    println("Beginning regression of $formula")
    # If weights is a string, convert to column vector
    if typeof(weights) <: String
        try
            weights = df[!, weights]
        catch e
            if isa(e, ArgumentError)
                println("ERROR: weights argument ($weights) is not a column in dataframe.\n")
            end
            throw(e)
        end
    end

    # If weights is given, run WLS
    if typeof(weights) <: Vector  # run WLS
        reg = glm(formula, df, Normal(), IdentityLink(), wts=weights)
    # No weights, run OLS
    else
        reg = lm(formula, df)
    end

    # If not clustering, use Heteroskedastic-consistent
    if cluster === nothing
        hc = eval(Symbol(varcov))
        V = vcov(reg, hc())
        SE = stderror(reg, hc)

    else  # clustering, use cluster-robust heteroskedasticity-consistent
        V = vcov(CRHC1(cluster, df), reg)
        SE = stderror(CRHC1(cluster, df),reg)
    end

    # Combine new SEs into dataframe with coefficients
    new_reg = apply_vcov(reg, SE)


    latexify(new_reg)
    filestub = filestub === nothing ? "reg" : filestub
    filename = append_time ? "$(filestub)_$(Dates.now()).tex" : "$filestub.tex"
    write(filename, string(latexify(new_reg, env=:tabular)))

    return Dict(:reg => new_reg,
                :β => coef(reg),
                :V => V,
                :SE => SE)
end


function reg2(df::DataFrame,
    formula::FormulaTerm; 
    weights::Union{String, Vector, Nothing}=nothing,
    varcov::String="HC3",
    cluster::Union{String, Nothing}=nothing,
    filestub::Union{String, Nothing}=nothing,
    append_time::Bool=true)
    #! can I implement Econometrics.jl fit? 
    # issues:
    #   - HC3 doesn't seem to work (OutOfMemory() error)
    #   - latexify(model) does not work on the econometrics.jl model. Can I convert to statsmodel?
end





#==============================================================================
SETUP
    key:
    fn = filename
    dir = directory
    df = dataframe
==============================================================================#
# Save single county (01001) to compare to built results
f = joinpath(root, "CountyAnnualTemperature1950to2012.csv")
if !isfile(f)
    df_ex = subset(df_ex, :fips => ByRow(==(01001)), skipmissing=true)
    CSV.write(joinpath(root, "CountyAnnualTemperature1950to2012.csv"), df_ex)
end


# Open linear piecewise stata results
# df3 = DataFrame(load(joinpath(root, "bestLinearModel", "corn_year1950_2020_month3_8.dta")))














#==============================================================================
1.1 TEMPERATURE AGGREGATION
==============================================================================#
#=================================================
1.1.1 Construct 4 temperature response variables
1.1.2 Aggregate to year-county
=================================================#
# Create a new file with temperature variables for example county (01001)
# This will be compared to the file given by Reed (example_fn)
if !isfile(joinpath(root, "CountyAnnualTemperature_aaron.csv"))
    create_temperature_vars()
end
# This file will NOT be used for the rest of the problem set
# This was just an exercise to create the variables.





#==============================================================================
1.2 US Climate Impacts
==============================================================================#
#=================================================
1.2.0 Merge Income and Employment onto climate data
=================================================#

# Extract and read employment file from zip
df_employ = df_from_zip(root, zip_fn, employment_fn)
df_employ[!, :fips] = parse.(Int32, df_employ[!, :fips])

# Download and load example county file to compare variables
df_temp = df_from_url(example_url, joinpath(root, example_fn))

# drop first row because it's missing
df_temp = last(df_temp, nrow(df_temp)-1)

# Join onto climate data
df = innerjoin(df_temp, df_employ, on=[:fips, :year])

# Make year and fips categegorical for regressions
df[!,:year] = categorical(df[!,:year])
df[!,:fips] = categorical(df[!,:fips])

# Create log(farm employees in county) and replace log(0) with missing
df[!, :emp_farm_ln] = replace(log.(df[!, :emp_farm]), -Inf => missing)

# Create log tranformed inc_farm_prop_income/population 
logfn(x) = x===missing ? missing : (x ≤ 0 ? missing : log(x))
df[!, :inc_farm_prop_inc_lpc] = logfn.(df[!,:inc_farm_prop_income] ./ df[!,:pop_population])
# 40,385 missing

# Convert the year into a Date
df[!, :Date] = Date.(convert(Vector{Int}, df[!,:year]))


#=================================================
1.2.1 emp_farm vs binned temperature
    Explore the relationship between log tranformed emp_farm 
    and the vector of binned temperature controls.
    Include additional controls for county FE, year FE.
    Provide an interpretation of the coefficient of the 32+ bin.
=================================================#
formula121 = @formula(emp_farm_ln ~ tempB0 + temp0to4 + temp4to8 + temp8to12 + temp12to16 + temp16to20 + temp24to28 + temp28to32 + tempA32 +
    fe(fips) + fe(year)
)
# Reference category = temp20to24 (days in county-year where tAvg ∈ (20, 24]C)

reg121 = @time reg(df, formula121, Vcov.cluster(:fips))
regtable(reg121, 
         renderSettings = latexOutput("reg121_$(Dates.now()).tex")
)
se121 = .√[reg121.vcov[i,i] for i in 1:length(reg121.coef)]
ub121 = reg121.coef + 1.96*se121
lb121 = reg121.coef - 1.96*se121
df121 = DataFrame(tAvg=[-2,2,6,10,14,18,26,30,34],
                  yhat=reg121.coef,
                  yerror=1.96*se121,
                  y_ub=ub121,
                  y_lb=lb121)
df121 = reduce(vcat, [[df121]; [DataFrame(tAvg=22,yhat=0,yerror=0,y_lb=0,y_ub=0)]])
df121 = sort!(df121, :tAvg)

plot(df121[!,:tAvg], repeat([0], nrow(df121)), c="gray", s=:dash, label="")
@df df121 plot!(:tAvg, :y_lb, w=0, msw = 0, ms = 0, c=1,
    fillrange=:y_ub, fillalpha=0.35,
    label="95% CI")
xlabel!("Average Temperature")
ylabel!("log(Farm Employment)")
title!("Binned Temperature Response: Farm Employment      ")
p121 = @df df121 plot!(:tAvg, :yhat, 
    label="Predicted Response", 
    legend=:bottomright, c=2, shape=:circle, markerstrokewidth=0)
savefig(p121, "plot121-bins.svg")

#! % change in employees for the average county for each day above 32.














#=================================================
1.2.2 per capita farm prop income vs cubic spline temperature
    Explore the relationship between log tranformed 
    inc_farm_prop_income/population (i.e. log(per capita 
    farm prop income) using the restricted cubic spline. 
    Include additional controls for county FE and year FE.
    Plot the predicted marginal effects with associated 
    confidence intervals, and compare to the binned 
    temperature response function above.

    We want to construct the estimated effect of average 
    temperature of a day on inc_farm_prop_inc_lpc
    Calculate point estimates and Conf Intervals 
    for linear combination of parameters
=================================================#
# Regress inc_farm_prop_inc_lpc on spline variables and fixed effects
formula122 = @formula(inc_farm_prop_inc_lpc ~ splineC1 + splineC2 + splineC3 + splineC4 +
                                    fe(year) + fe(fips))
reg122 = @time reg(df, formula122, Vcov.cluster(:fips))
regtable(reg122,renderSettings = latexOutput("reg122_$(Dates.now()).tex"))


"""Return DF of avg temp, cubic spline, and predicted outcome variables needed for regressions and plotting.
    
    @param: reg_output: output from a regression on 4 cubic spline variables
    
    - create an x-axis for plotting: define an array of avg temps from 0 to 40°C, 0.25 steps
    - define the spline knots/breakpoints
    - Create spline basis variables. After regressing, the absolute level of the outcome (y-axis)
        is not identified because the county & year fixed effects shift the intercept around.
        So these basis will be compared to the spline basis variables created for X°C.
    - Create spline basis variables for X°C
    - Create relative basis variables (base - X°C variables) which is equivalent to shifting
        the base prediction by the X°C prediction, but allows for adjusted standard errors
        in the delta method calculation.
"""
function plotdf_init(reg_output; X=20)
    df = DataFrame(tAvg=0:0.25:40)      # Initialize x-axis
    knots = [0 8 16 24 32]              # breakpoints in spline
    # Create base basis spline variables
    df = add_cubic_spline_vars(df, :tAvg, knots, newbasename="splinebaseC")
    # Create X°C basis spline variables
    df[!, "XC"] = repeat([X], nrow(df))
    df = add_cubic_spline_vars(df, "XC", knots, newbasename="XC")
    # Create relative basis spline variables (relative to X°C)
    for k in 1:4
        df[!, "splinerelativeC$k"] = df[!, "splinebaseC$k"] - df[!, "XC$k"]
    end
    # Construct the predicted spline outcome
    df[!, :yhat] = Matrix(df[!,r"splinerelativeC"]) * reg_output.coef
    # Return the x-axis, predicted outcome, splinerelCnd relative spline variables (for delta method SE estimation)
    return df[!, [["tAvg", "yhat"]; ["splinerelativeC$k" for k ∈ 1:4]]]
end



#=
Standard Errors can be calculated using:
- the delta method
- analytically using additivity of variance and covariance (b/c linear)
- bootstrap resampling
=#




#= Let's try the analytical method!
    Define x = a single day's average temperature
    Define f(x) = the effect of moving a single day's average temperature from 20°C to x
    We want to estimate the 95% confidence interval of f(x) for each x ∈ [0°C, 40°C]
    Each of our spline basis variables is a function of the average temperature x
    Let Zₖ = gₖ(X) to be the spline basis variables for k ∈ {1,2,3,4}
    We run a regression on the spline basis variables to fit the spline to farm employment:
        y = ΣₖZₖβₖ + year + county (year and county fixed effects)
    After running our regression on the spline basis variables and fixed effects,
        our predicted response of percapita farm employment from moving a day 
        with average temp 0°C to average temp x is:
        f(x) = ΣₖZₖβₖ
    Then we would want to calculate the variance of the prediction f(x) at each value of x
    So Var(f(x)) = Var(ΣₖZₖβₖ) = ΣₖZₖ² Var(βₖ) + 2 Σ{k}Σ{j>k} ZₖZⱼ Cov(βₖ,βⱼ)
    ... turns out this is identical to the delta method formula, even after making f(x)
    relative to a specific temperature (20°C)
=#


# Let's try the Delta Method!
function get_diagonal(M::Matrix)
    n,_ = size(M)
    return [M[i,i] for i∈1:n]
end


"""Return DF with upper and lower bounds created from delta method.

    Delta method: SE =  √(∂f(β)/∂β ⋅ VCOV ⋅ ∂f(β)/∂β)
    Calculate the β-derivative of the function of parameters
    f(β) = Σxₖ⋅βₖ ⟹ ∂f∂βₖ = xₖ
    ∂fβ_∂β = Matrix(df_plot[!,r"splineC"])
    SEs are the diagonal of .√(∂fβ_∂β * VCOV * ∂fβ_∂β'))
    where VCOV is the variance-covariance matrix of parameter estimates
    from the regression results
"""
function deltamethod_CIbounds!(df, reg_output; column_stem="splinerelativeC", α=0.05)
    VCOV = reg_output.vcov
    ∂fβ_∂β = Matrix(df[!, Regex(column_stem)])
    n = reg_output.nobs
    df[!, :SE] = .√(get_diagonal(∂fβ_∂β * VCOV * ∂fβ_∂β'))
    # Use a t-distribution just in case the sample is small
    tcrit = quantile(TDist(reg_output.dof_residual), 1-α/2)
    # Use √(1 + 1/n) correction because both the mean (predicted value) and the variance are unknown estimated parameters
    width = tcrit * df[!, :SE] * √(1 + 1/n)
    df[!, :y_lb] = df[!, :yhat] - width
    df[!, :y_ub] = df[!, :yhat] + width
    #! If the sample size is small, consider using the second-order delta method
    return df
end


# Let's try bootstrapping!
"""Given a dataframe sample and dataframe with tAvg values to predict on,
    estimate the prediction vector"""
function prediction_sample(df, df_plot)
    reg_ = reg(df, formula122, Vcov.cluster(:fips))
    df2 = DataFrame(tAvg=0:0.25:40)
    df2[!, :yhat] = Matrix(df_plot[!,r"splinerelativeC"]) * reg_.coef
    return df2[!, [:tAvg, :yhat]]
end


"""Return DF of confidence intervals for each temperature for bootstraped simulation results"""
function confidence_interval(df, α)
    df_bs = @chain df begin
        groupby(:tAvg)
        combine(:yhat => mean => :y_mean,
                :yhat => (x -> quantile!(x, α/2)) => :y_lb,
                :yhat => (x -> quantile!(x, 1-α/2)) => :y_ub,
                :yhat => minimum => :y_min, :yhat => maximum => :y_max)
    end
    return df_bs
end


"""Return dataframe of bootstrapped predictions with α confidence intervals."""
function bootstrap_predicted_outcome(df::DataFrame,
                                     prediction_function::Function;
                                     screenwidth=50,
                                     nsimulations=10000,
                                     α=0.05
                                     )
    println("\nStarting simple bootstrap")
    cols = [:year, :fips, :tAvg, :inc_farm_prop_inc_lpc, :splineC1, :splineC2, :splineC3, :splineC4]
    preds = []
    t = @timed for k ∈ 1:nsimulations
        k%ceil(nsimulations/screenwidth) == 0 ? print("⋅") : nothing
        sample_rows = sample(1:nrow(df), nrow(df), replace=true)
        df_sample = df[sample_rows, cols]
        append!(preds, [prediction_function(df_sample, df_plot)])
    end
    # Vertically concatenate the bootstrap predictions
    preds_comb = reduce(vcat, preds)
    # Get mean and 95% CI for each value of tAvg
    df_bs = confidence_interval(preds_comb, α)
    print("$(t.time) seconds\n")
    return df_bs
end


"""Return dataframe of cluster-bootstrapped predictions with α confidence intervals.

    Clustered Bootstrap (over each fips county code)
    Need to sample 3055 counties of 3055 counties (with replacement)
    Requires repeating some counties in the data and keeping all years
    of their data
"""
function bootstrap_predicted_outcome(df::DataFrame,
                                     prediction_function::Function,
                                     cluster::Union{Symbol, String};
                                     screenwidth=50,
                                     nsimulations=10000,
                                     α=0.05
                                     )
    println("\nStarting clustered bootstrap")
    cols = [:year, :fips, :tAvg, :inc_farm_prop_inc_lpc, :splineC1, :splineC2, :splineC3, :splineC4]
    preds = []
    df_grouped = groupby(df[!, cols], cluster)
    ngroups = length(df_grouped)
    t = @timed for k ∈ 1:nsimulations
        k%ceil(nsimulations/screenwidth) == 0 ? print("⋅") : nothing
        sample_idx = sample(1:ngroups, ngroups; replace = true, ordered = false)
        df_sample = reduce(vcat, [df_grouped[i] for i in sample_idx])
        append!(preds, [prediction_function(df_sample, df_plot)])
    end
    # Vertically concatenate the bootstrap predictions
    preds_comb = reduce(vcat, preds)
    # Get mean and 95% CI for each value of tAvg
    df_bs = confidence_interval(preds_comb, α)
    print("$(t.time) seconds\n")
    return df_bs
end


# Initialize dataframe used for plotting (x-axis with temperatures, spline variables)
df_plot = plotdf_init(reg122)
# Delta method on recentered data (outcome centered on predicted outcome at 20°C)
df_plot = deltamethod_CIbounds!(df_plot, reg122)


# Get mean and 95% CI for each value of tAvg using bootstraps
if 1==0
    nsimulations = 100
    # Simple bootstrap (over all individuals)
    df_bs = bootstrap_predicted_outcome(df, prediction_sample, nsimulations=nsimulations)
    # Clustered bootstrap
    df_bs_cluster = bootstrap_predicted_outcome(df, prediction_sample, :fips; nsimulations=nsimulations)   


    # Plot Delta method and bootstraps confidence intervals and predicted outcome
    s = 0; dpi=400;
    @df df_plot plot(:tAvg, :y_lb, w=0, msw = 0, ms = s, c=1,
        fillrange=:y_ub, fillalpha=0.35,
        label="95% CI Delta Method")
    @df df_bs plot!(:tAvg, :y_lb, w=0, msw = 0, ms = s, c=2,
        fillrange=:y_ub, fillalpha=0.35,
        label="95% CI BootStrap ($nsimulations)")
    @df df_bs_cluster plot!(:tAvg, :y_lb, w=0, msw = 0, ms = s, c=3,
        fillrange=:y_ub, fillalpha=0.35,
        label="95% CI BootStrap-clustered ($nsimulations, FIPS)")
    p8 = @df df_plot plot!(:tAvg, :yhat, label="Predicted Response", legend=:bottomleft)
    savefig(p8, "plot122-delta-boot-bootcluster.svg")
end


# Plot Delta method confidence intervals and predicted outcome
plot(df_plot[!,:tAvg], repeat([0], nrow(df_plot)), c="gray", s=:dash, label="")
@df df_plot plot!(:tAvg, :y_lb, w=0, msw = 0, ms = 0, c=1,
    fillrange=:y_ub, fillalpha=0.35,
    label="95% CI")
xlabel!("Average Temperature")
ylabel!("log(Farm Income per capita)")
title!("Spline Temperature Response: Farm Employment      ")
p9 = @df df_plot plot!(:tAvg, :yhat, 
    label="Predicted Response", legend=:bottomleft, c=2)
savefig(p9, "plot122-spline.svg")



















#=================================================
1.2.3 Treatment heterogeneity
    Use the binned temperature estimator to design a 
    test for whether we observe treatment effect
    heterogeneity. One possibility for this test is 
    to interact the temperature bins with the average
    # of days in a county for which the temperature 
    falls into each respective bin (i.e. are counties
    that experience more 32+ days more/less responsive 
    to temperature).
=================================================#
#! the temperature bin variables are degree-days in the bin
#! which are averaged over all grid points in the county
#! so aren't they already "the average # of days in a county
#! for which the temperature falls inot each bin"?
#! Simon ==> use 10-year rolling averages


function add_moving_avgs!(df; windowsize=10)
    # group each county
    cols = [["fips", "Date"]; names(df, r"temp[^P]")]
    gdf = groupby(df[!, cols], :fips)
    dfs = []
    # for each county, calculate a 10-year backward-looking moving average
    for group in gdf
        # Create TimeArray
        temp = TimeArray(group, timestamp = :Date)
        # create moving average of last 10 years
        temp = moving(mean, temp, windowsize)
        # convert back to dataframe
        temp = DataFrame(temp)
        append!(dfs, [temp])
    end
    # concatenate county dataframes and leftjoin onto df using fips and Date
    df_avg = reduce(vcat, dfs)
    renamingdict = [names(df, r"temp[^P]") .=> names(df, r"temp[^P]") .* "_avg"; ["timestamp" => "Date"]]
    df_avg = DataFrames.rename(df_avg, renamingdict)
    df = leftjoin(df, df_avg, on=[:fips, :Date])
    return df
end

# For each bin, add the 10-year moving average
df = add_moving_avgs!(df)


# Calculate the deviance from the running mean
for col in names(df, r"temp[^P].*(?<!_avg)$")
    df[!, "$(col)_dev"] = df[!, col] - df[!, "$(col)_avg"]
end


# Regress the outcomes on the deviance from the running average bins
function reg123(df, y; cluster=:fips, regex=r"temp(?!20to24).*_dev")
    formula_ = term(y) ~ sum([term.(names(df, regex)); [term(fe(:fips)), term(fe(:year))]])
    reg_ = @time reg(df, formula_, Vcov.cluster(cluster))
    se = .√[reg_.vcov[i,i] for i in 1:length(reg_.coef)]
    ub = reg_.coef + 1.96*se
    lb = reg_.coef - 1.96*se
    df2 = DataFrame(tAvg=[-2,2,6,10,14,18,26,30,34],
                    yhat=reg_.coef,
                    yerror=1.96*se,
                    y_ub=ub,
                    y_lb=lb)
    df2 = reduce(vcat, [[df2]; [DataFrame(tAvg=22,yhat=0,yerror=0,y_lb=0,y_ub=0)]])
    df2 = sort!(df2, :tAvg)
    return df2
end


# Plot deviance bins: above avg # of days this year that were in this bin (days above the 10-year average # of days)
labels = Dict(:emp_farm_ln => "Farm Employment", :inc_farm_prop_inc_lpc => "Farm Income per capita")
Dict(:emp_farm_ln => "log(Farm Employment)", :inc_farm_prop_inc_lpc => "log(Farm Income per capita)")
for y ∈ [:emp_farm_ln, :inc_farm_prop_inc_lpc]
    temp = reg123(df, y)
    plot(temp[!,:tAvg], repeat([0], nrow(temp)), c="gray", s=:dash, label="")
    @df temp plot!(:tAvg, :y_lb, w=0, msw = 0, ms = 0, c=1,
    fillrange=:y_ub, fillalpha=0.35,
    label="95% CI")
    xlabel!("Average Temperature")
    ylabel!("log($(labels[y]))")
    title!("Binned DDay Deviance from Avg: $(labels[y])          ")
    p121 = @df temp plot!(:tAvg, :yhat, 
        label="Predicted Response", 
        legend=:bottomright, c=2, shape=:circle, markerstrokewidth=0)
    savefig(p121, "plot123 $(labels[y]) deviance.svg")
end


# Plot basic bins again: # of days in a county-year that day's avg temp was in this bin
for y ∈ [:emp_farm_ln, :inc_farm_prop_inc_lpc]
    temp = reg123(df, y; regex=r"^temp(?!20to24)[^P].*[^vg]$")
    plot(temp[!,:tAvg], repeat([0], nrow(temp)), c="gray", s=:dash, label="")
    @df temp plot!(:tAvg, :y_lb, w=0, msw = 0, ms = 0, c=1,
    fillrange=:y_ub, fillalpha=0.35,
    label="95% CI")
    xlabel!("Average Temperature")
    ylabel!("log($(labels[y]))")
    title!("Binned Temperature Response: $(labels[y])         ")
    p121 = @df temp plot!(:tAvg, :yhat, 
        label="Predicted Response", 
        legend=:bottomright, c=2, shape=:circle, markerstrokewidth=0)
    savefig(p121, "plot123 $(labels[y]) bins.svg")
end

























#==============================================================================
2.2 Hedonic Air Quality Analysis
    Does Air Quality Get Capitalized into Housing Prices? The outcome of interest
    is  the change in county housing prices during the 1970s. We want to estimate
    the  “causal” effect of air pollution changes on housing price changes.
    According  to hedonic price theory, the housing market may be used to estimate
    the  implicit prices of clean air and the economic value of pollution
    reductions to individuals. A statistically significant negative relationship
    between changes in property values and pollution levels across counties is
    interpreted as  evidence that clean air has economic benefits. (For a summary
    of the theory,  you could read: Rosen, Sherwin, “The Theory of Equalizing
    Differences,” Chapter 12 in Handbook of Labor Economics, Volume 1, 1986, pp.
    641-92.) A basic model for the change in housing prices at the county level
    could be:  Change in housing price = g(economic shocks, changes in county
    characteristics,  change in air pollution)
==============================================================================#
#=================================================
2.2.0 Load Data
=================================================#
println("Reading $pollution_fn from zip folder.")
extract_file_from_zip(root, zip_fn, pollution_fn)
pollution = DataFrame(load(joinpath(root, pollution_fn)))

# variables
labels = Dict(
    :dlhouse => "\\\Delta log-housing values",
    :dgtsp => "\\\Delta"
)
# dlhouse = change in log-housing values from 1970 to 1980
# dgtsp = change in the annual mean of TSPs from 69-72 to 77-80
# tsp7576 = indicator: county was regulated EPA in 1975 or 1976
# mtspgm74 = annual geometric mean of TSPs in 1974.
# pop7080 = 1970 + 1980 county populations; weight for regressions

# dincome = change in income per-capita.
# dunemp = change in unemployment rate.
# dmnfcg = change in % manufacturing employment.

# ddens = 1970-80 change in population density.
# dwhite = change in fraction of population that is white.
# dage65 = change in fraction over 65 years old.
# dhs = change in fraction with at least a high school degree.
# dcoll = change in fraction with at least a college degree.
# durban = change in fraction living in urban area.
# downer = change in fraction of houses that are owner-occupied.
# dplumb = change in fraction of houses with plumbing.
# drevenue = change in government revenue per-capita.
# dtaxprop = change in property taxes per-capita,
# depend = change in general expenditures per-capita.

# Other varaible definitions in Walker-ProblemSet1.pdf



#=================================================
2.2.1  Housing prices vs air pollution

Estimate the relationship between changes in air
pollution and housing prices, both not adjusting
and adjusting for other housing price
determinants (use pop7080 as weights). What do
your estimates imply and do they make sense?
Describe the potential omitted variables biases.
What is the likely relationship between economic
shocks and pollution and housing price changes?
Using the observable measures of economic shocks
(dincome, dunemp, dmnfcg), provide evidence on
this.
=================================================#
# Remove missing and convert to Float64
# df1 = convert_float32(pollution, [:dlhouse, :dgtsp, :pop7080])
# df212 = convert_float32(pollution, [[:dlhouse, :pop7080]; cols212])
# reg211 = regression(df1, @formula(dlhouse ~ dgtsp); weights="pop7080", filestub="reg211")

# Reg: dlhouse on dgtsp, weight with pop7080
reg211 = @time reg(pollution, @formula(dlhouse ~ dgtsp), Vcov.robust(); weights=:pop7080)
regtable(reg211, renderSettings = latexOutput("reg211base_$(Dates.now()).tex"))

# Reg: dlhouse on dgtsp, weight with pop7080,
# controlling for ddens dwhite dfeml dage65
# dhs dcoll durban dpoverty dvacant downer
# dplumb dtaxprop (depend?)
cols221 = [:dgtsp, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
formula221 = (term(:dlhouse) ~ sum(term.(Symbol.(cols221))))
reg221 = @time reg(pollution, formula221, Vcov.robust(); weights=:pop7080)
regtable(reg221, renderSettings = latexOutput("reg221controls_$(Dates.now()).tex"))

#! What do your estimates imply and do they make sense?
#! Describe the potential omitted variables biases
#=
both these seem to indicate that increased pollution increases
housing values! This doesn't makes sense since air pollution
should be a disamenity. This result could be caused by OVB,
if the omitted variable has correlation with housing price in
the same direction as the correlation with pollution (ie,
if x2 is the OV, cov(x2, pollution) and cov(x2, houseprice) are both >0 or both <0).
For example, if there is a recession, and the economy slows, 
this generally decreases the amount of pollution being produced
by firms supplying goods to the economy. So an indicator of positive
economic activity would be positively correlated with pollution. At the same
time, a recession would likely cause the average price of housing to decrease.
So an indicator of positive
economic activity would be positively correlated with housing price.
Then bias from an omitted economic variable that is positively correlated with
both pollution and housing prices would cause an upward bias in our estimate
of the partial derivative of house price changes w.r.t. pollution changes.
If the bias is large enough, it could flip the sign of the coefficient,
and since our coeffience on pollution is nearly statistically indistinguishable from 0
after controlling for other observed housing determinants, upward OVB could
be likely here.

There was a large recession in the 70s, so variables describing the impact of 
the recession on the supply and demand for housing seem like they
would be important to control for. For example, many manufacturing jobs were lost
during the recession (more relative to other sectors)
which affects the housing market. We can test three of our economic
indicators to see how correlated with pollution they are in order to see if
we should be concerned about OVB.
=#




#! What is the likely relationship between economic
#! shocks and pollution and housing price changes?
#! Using the observable measures of economic shocks
#! (dincome, dunemp, dmnfcg), provide evidence on
#! this.
#=
Economic shocks can change pollution because pollution
is generated by economic activity (i.e., if there is 
a recession, there will be less pollution) and economic shocks
can change what people are able/willing to pay for a house 
(ie, recession -> house prices decrease), which creates
endogeniety (house price is correlated with the error term)=#
# Regress economic shocks on pollution
reg221b = @time reg(pollution, @formula(dgtsp ~ dincome + dunemp + dmnfcg), Vcov.robust(); weights=:pop7080)
regtable(reg221b, renderSettings = latexOutput("reg221pollution-econshocks_$(Dates.now()).tex"))

#! Change in manufacturing employment has a strong effect on change in pollution
#! A 10% increase in manufacturing employment between 1970 and 1980 is correlated
#! with an increase in the 1970-1980 change in mean TSPs by about 10.6 units.
#! the change 





#=
We could further
control for unemployment more generally and average county
income changes to ensure we are comparing willingness to pay across
similar counties.  After controlling for these three economic 
indicators, we see the hedonic estimate of willingness to pay for
changes in pollution become statistically indistinguishable from zero.
=#
cols221c = [:dgtsp, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop, :dincome, :dunemp, :dmnfcg]
formula221c = (term(:dlhouse) ~ sum(term.(Symbol.(cols221c))))
reg221c = @time reg(pollution, formula221c, Vcov.robust(); weights=:pop7080)
regtable(reg221c, renderSettings = latexOutput("reg221controls-econshocks_$(Dates.now()).tex"))







#=================================================
2.2.2  Attainment status as IV for pollution changes

Suppose that federal EPA pollution regulation is
a potential instrumental variable for pollution
changes during the 1970s. What are the
assumptions required for mid-decade regulatory
status (tsp7576) to be a valid instrument for
pollution changes when the outcome of interest is
housing price changes? Provide evidence on the
relationship between the regulatory status
indicator and the observable economic shock
measures. Interpret your findings.
=================================================#

reg222 = @time reg(pollution, @formula(tsp7576 ~ dincome + dunemp + dmnfcg), Vcov.robust(); weights=:pop7080)
regtable(reg222, renderSettings = latexOutput("reg222status-econshocks_$(Dates.now()).tex"))

#=
1. Relevance - first stage, next question
2. Exogeneity - not correlated with the economic shocks
=#












#=================================================
2.2.3 First stage, reduced form

Document the first-stage relationship between
regulation and air pollution changes and the
reduced form relationship between regulation and
housing price changes, both not adjusting and
adjusting for other covariates. Interpret your
findings. How does two-stage least squares use
these two equations? Now estimate the effect of
air quality changes on housing price changes
using two-stage least squares and the tsp7576
indicator as an instrument (not conditioning and
conditioning on other observables). Interpret the
results.
=================================================#

#! first-stage relationship between regulation and air pollution changes
#! Relevance

# Reg: dgtsp on tsp7576, weight with pop7080
reg223a = @time reg(pollution, @formula(dgtsp ~ tsp7576), Vcov.robust(); weights=:pop7080)
regtable(reg223a, renderSettings = latexOutput("reg223afirst-base_$(Dates.now()).tex"))

# Reg: dgtsp on tsp7576, weight with pop7080, with controls
cols223b = [:tsp7576, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
formula223b = (term(:dgtsp) ~ sum(term.(Symbol.(cols223b))))
reg223b = @time reg(pollution, formula223b, Vcov.robust(); weights=:pop7080)
regtable(reg223b, renderSettings = latexOutput("reg223bfirst-controls_$(Dates.now()).tex"))

# Don't need to include econ shocks because I showed in previous question it was exogenous to the shocks


#! reduced form relationship between regulation and housing price changes
#! Reduced form estimate

# Reg: dlhouse on tsp7576, weight with pop7080
reg223c = @time reg(pollution, @formula(dlhouse ~ tsp7576), Vcov.robust(); weights=:pop7080)
regtable(reg223c, renderSettings = latexOutput("reg223creduced-base_$(Dates.now()).tex"))

# Reg: dlhouse on tsp7576, weight with pop7080, with controls
cols223d = [:tsp7576, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
formula223d = (term(:dlhouse) ~ sum(term.(Symbol.(cols223d))))
reg223d = @time reg(pollution, formula223d, Vcov.robust(); weights=:pop7080)
regtable(reg223d, renderSettings = latexOutput("reg223dreduced-controls_$(Dates.now()).tex"))

#! How does two-stage least squares use these two equations?
#=
These are the first two coef. ratio: e = c / a,  f = d / b
=#




#! air quality changes on housing price changes using two-stage least squares and the tsp7576 indicator as an instrument
#! IV estimate

# 2SLS Reg: dlhouse on dgtsp, weight with pop7080
reg223e = @time reg(pollution, @formula(dlhouse ~ (dgtsp ~ tsp7576)), Vcov.robust(); weights=:pop7080)
regtable(reg223e, renderSettings = latexOutput("reg223eIV-base_$(Dates.now()).tex"))

# 2SLS Reg: dlhouse on dgtsp, weight with pop7080, with controls
cols223f = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
formula223f = (term(:dlhouse) ~ sum(term.(Symbol.(cols223f))) + (term(:dgtsp) ~ term(:tsp7576)))
reg223f = @time reg(pollution, formula223f, Vcov.robust(); weights=:pop7080)
regtable(reg223f, renderSettings = latexOutput("reg223fIV-controls_$(Dates.now()).tex"))


#! Interpret the results.
#=
A 1-unit increase in mean TSPs from 1970 to 1980 results in an average
housing price change from 1970 to 1980 being 0.7 percentage points larger (more positive).
=#










#=================================================
2.2.4 Mean of TSPs as IV for pollution changes

In principle, the regulation indicator variable
should be a discrete function of pollution levels
in 1974. Specifically, the EPA is supposed to
regulate those counties in 1975 who had either an
annual geometric mean of TSPs above 75 units
(mg/m3) or a 2nd highest daily concentration
above 260 units in 1974. Based on this notion,
redo part (3) using mtspgm74 as an instrument for
pollution changes. Interpret your findings.
=================================================#

# Create an indicator for mtspgm74<75
pollution[!, :above75] = pollution[!, :mtspgm74] .> 75


#! Exogeniety
reg224_ = @time reg(pollution, @formula(above75 ~ dincome + dunemp + dmnfcg), Vcov.robust(); weights=:pop7080)
regtable(reg224_, renderSettings = latexOutput("reg224_status-econshocks_$(Dates.now()).tex"))



#! first-stage relationship between regulation and air pollution changes
#! Relevance
# Reg: dgtsp on above75, weight with pop7080
reg224a = @time reg(pollution, @formula(dgtsp ~ above75), Vcov.robust(); weights=:pop7080)
regtable(reg224a, renderSettings = latexOutput("reg224afirst-base_$(Dates.now()).tex"))

# Reg: dgtsp on above75, weight with pop7080, with controls
cols224b = [:above75, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
formula224b = (term(:dgtsp) ~ sum(term.(Symbol.(cols224b))))
reg224b = @time reg(pollution, formula224b, Vcov.robust(); weights=:pop7080)
regtable(reg224b, renderSettings = latexOutput("reg224bfirst-controls_$(Dates.now()).tex"))


#! reduced form relationship between regulation and housing price changes
#! Reduced form estimate
# Reg: dlhouse on above75, weight with pop7080
reg224c = @time reg(pollution, @formula(dlhouse ~ above75), Vcov.robust(); weights=:pop7080)
regtable(reg224c, renderSettings = latexOutput("reg224creduced-base_$(Dates.now()).tex"))

# Reg: dlhouse on above75, weight with pop7080, with controls
cols224d = [:above75, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
formula224d = (term(:dlhouse) ~ sum(term.(Symbol.(cols224d))))
reg224d = @time reg(pollution, formula224d, Vcov.robust(); weights=:pop7080)
regtable(reg224d, renderSettings = latexOutput("reg224dreduced-controls_$(Dates.now()).tex"))



#! air quality changes on housing price changes using two-stage least squares and the above75 indicator as an instrument
#! IV estimate
# 2SLS Reg: dlhouse on dgtsp (above75 as IV), weight with pop7080
reg224e = @time reg(pollution, @formula(dlhouse ~ (dgtsp ~ above75)), Vcov.robust(); weights=:pop7080)
regtable(reg224e, renderSettings = latexOutput("reg224eIV-base_$(Dates.now()).tex"))

# 2SLS Reg: dlhouse on dgtsp (above75 as IV), weight with pop7080, with controls
cols224f = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
formula224f = (term(:dlhouse) ~ sum(term.(Symbol.(cols224f))) + (term(:dgtsp) ~ term(:above75)))
reg224f = @time reg(pollution, formula224f, Vcov.robust(); weights=:pop7080)
regtable(reg224f, renderSettings = latexOutput("reg224fIV-controls_$(Dates.now()).tex"))















#=================================================
2.2.5 Regression Discontinuity, Diff-in-Diffs?

Describe how one could use this discontinuity in
treatment assignment to derive an alternative
estimate of the capitalization of pollution
changes. Under what conditions will this estimate
be valid? Based on this concept, estimate the
nonparametric bivariate relationships between
pollution changes and 1974 TSPs levels and
housing price changes and 1974 TSPs levels. Do
this using the lowess STATA command or a similar
command in R (experiment with relatively small
bandwidths between 2-4). Plot the two conditional
mean functions on either side of the
discontinutiy (changes in air quality for
discrete 1974 values on either side of
threshold). Interpret your findings and relate
them to the results in (4).
=================================================#















#=================================================
2.2.6 

Now use linear regression to estimate the
predicted change in housing prices based on the
control variables (excluding the change in TSPs).
This provides a single-index measure of the
housing price changes predicted to occur due to
other variables changing. Estimate the
nonparametric bivariate relation between this
index and 1974 TSPs levels and plot it against
the smoothed housing price changes from (5).
Explain how this provides an indirect test of the
smoothness condition required by the regression
discontinuity design and interpret your findings.
=================================================#















#=================================================
2.2.7 

A number of counties with annual geometric mean
TSPs below 75 units in 1974 were regulated due to
having as few as 2 bad days. This implies that we
can compare regulated and unregulated counties
with identical average TSPs levels in the
regulation selection year (be- low 75 units). For
those counties with 1974 mean TSPs between 50 and
75 units, estimate the bivariate relation between
TSPs changes and 1974 TSPs levels separately for
regulated and unregulated counties and plot the
results. Now do the same for the relation between
log-housing price changes and 1974 TSPs levels.
Interpret your findings. Since there are fewer
observations, you may need to use bigger
bandwidths than those in part (5) (e.g.,
bandwidths between 6-9).
=================================================#















#=================================================
2.2.8 

Under what assumptions will two-stage least
squares identify the average treatment effect
(ATE)? If ATE is not identified, describe what
may be identifiable with two-stage least squares
estimation. Under what conditions is this effect
identified? Give some intuition on what this
effect may represent when one uses EPA regulation
as an instrument variable.
=================================================#















#=================================================
2.2.9 

Provide a concise synthesis/summary of your
results. Discuss the credibility of the various
research designs underlying the results.
=================================================#




















#==============================================================================
Extra Notes / leftovers
==============================================================================#


#=
Merge gridNumber with latlon then with Decennial Census Block population counts.
Or are they weighted using crop weights like in degreeDays.do? 

# lat-lon from degreeDays.do
use ../metaData/cropArea, clear;
gen longitude = -125 + mod(gridNumber-1,1405)/24;
label var longitude "longitude of grid centroid (decimal degrees)";
gen latitude  = 49.9375+1/48 - ceil(gridNumber/1405)/24;
label var latitude  "latitude of grid centroid (decimal degrees)";
merge 1:1 gridNumber using ../metaData/linkGridnumberFIPS;


# Decennial census block population counts
=#





    # combine(groupby(combine(gd, :tMax => mean => :tMax), :year), :tMax => mean)
    # Weighted average using Decennial Census Block population counts as weights
    # gridnum_fips_fn = "linkGridnumberFIPS.dta"
    # lookup = DataFrame(load(joinpath(root, gridnum_fips_fn)))
