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
# Pkg.add(["ZipFile", "StatFiles"])
using ZipFile  # unzip compressed .zip folders
using StatFiles  # read Stata files
using DataFrames
using CSV
using Statistics
using Dates
using Plots
using GLM
using CovarianceMatrices
using Latexify
Latexify.set_default(fmt = "%.3f")




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
    # Open zip file
    zarchive = ZipFile.Reader(zip_path)
    # Find file_name in zip archive
    files = [f.name for f in zarchive.files]
    idx = findfirst(x -> x == file_name, files)
    # Save file_name to disk
    write(save_path, read(zarchive.files[idx]))
end


"""Read file at save_path to dataframe, if not present, download from url."""
function df_from_url(url, save_path)
    if isfile(save_path)
        return DataFrame(load(joinpath(root, example_fn)))
    else
        download(example_url, joinpath(root, example_fn))
        return DataFrame(load(joinpath(root, example_fn)))
    end
end


"""Return a dataframe, with types Float64 and no missing for given columns."""
function convert_float32(df, cols)
    df1 = dropmissing(df[!, cols], disallowmissing=true)
    df1 = transform!(df1, names(df1) .=> ByRow(Float64), renamecols=false)
    return df1
end



#==============================================================================
                            TEMPERATURE FUNCTIONS
==============================================================================#
"""Return temperature threshold after converting to -1,1 interval (for sine function)"""
function convert_temp_threshold(Tmax, Tmin, Thresh)
    proportion = (Thresh - Tmin) / (Tmax - Tmin)
    return (2 * proportion - 1)
end


"""
Return portion (0 to 1) of day's temperature above temperature Threshold,
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


"""Degree day calculation from Snyder 1985, using integral of sine."""
function degree_day_calc(MAX, MIN, THR)
    M = (MAX + MIN) / 2  # average temp
    W = (MAX - MIN) / 2  # half temp range
    θ = asin((THR - M) / W)  # radian between -π/2 and π/2 where sine crosses THR
    DD = ((M - THR)(π / 2 - θ) + W * cos(θ)) / π
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


"""Return jth restricted cubic spline variable using the list of knots t.
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

After trying both the stata and Harrell formulas, neither matched the example file given.
After comparing the formulas, I recognized the only difference was a scaling factor of 
    1/(knots[last] - knots[first])^2
so I added that scaling factor to the Harrell formula and it produces results that match
the example file. So I probably just messed up the stata formula somehow. I assume the
example file was created in stata with mksplines.
"""
function cubic_spline_vars(df, varname, knots, newbasename="splineC")
    for i in 1:(length(knots)-1)
        newname = "$newbasename$i"
        println("Making $newname")
        df[!, newname] = cubic_spline_var_harrell(df[!, varname], i, knots)
    end
    return df
end


#! need to apply degree_day_single_day() to all day-gridpoints in df1.

#=
Using data on exposure to each 1-degree Celsius temperature interval, we approximate
the above integral... pg 7, appendix

First, we approximate g(h) using dummy variables for each three-degree 
temperature interval. This step function effectively regresses yield on 
season-total time within each temperature interval

The second speciﬁcation assumes g(h) is an m-th order Chebychev polynomial
=#


"""Return grid-day data aggregated to county-year level"""
function aggregate_df(df)
    # List of variables to average or sum over
    avg_list = [:tMin, :tMax, :tAvg, :prec]
    sum_list = [:splineC1, :splineC2, :splineC3, :splineC4]
    # Sum over all days in each year, for each grid point
    df[!, "date"] = Date(1960, 1, 1) + Day.(df.dateNum)
    df[!, "year"] = Year.(df.date)
    gd_gy = groupby(df, [:gridNumber, :year])
    df_gy = combine(gd_gy,
        avg_list .=> mean .=> avg_list,
        sum_list .=> sum .=> sum_list)
    gd_y = groupby(df_gy, :year)
    df_y = combine(gd_y, valuecols(gd_y) .=> mean .=> valuecols(gd_y))

    # combine(groupby(combine(gd, :tMax => mean => :tMax), :year), :tMax => mean)
    # Weighted average using Decennial Census Block population counts as weights
    gridnum_fips_fn = "linkGridnumberFIPS.dta"
    lookup = DataFrame(load(joinpath(root, gridnum_fips_fn)))

    return df_y
end





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
        V = vcov(hc(), reg)
        SE = stderror(hc(), reg)

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






#==============================================================================
                                    SETUP
key:
fn = filename
dir = directory
df = dataframe
==============================================================================#
# Set local directory paths and filenames
root = dirname(@__FILE__)
zip_fn = "Walker-ProblemSet1-Data.zip"
county_fn = "fips1001.dta"
employment_fn = "reis_combine.dta"
example_fn = "CountyAnnualTemperature1950to2012.dta"
pollution_fn = "poll7080.dta"
example_url = "https://www.dropbox.com/s/fnl1u0ix4e493vv/CountyAnnualTemperature1950to2012.dta?dl=0"


# Extract and read county fips file from zip
println("Reading $county_fn from zip folder.")
extract_file_from_zip(root, zip_fn, county_fn)
df_temp1 = DataFrame(load(joinpath(root, county_fn)))
# Extract and read emplyment file from zip
println("Reading $employment_fn from zip folder.")
extract_file_from_zip(root, zip_fn, employment_fn)
df_employ = DataFrame(load(joinpath(root, employment_fn)))
# Download and load example county file to compare variables
println("Downloading $example_fn from dropbox folder.")
df_ex = df_from_url(example_url, joinpath(root, example_fn))
df_ex = subset(df2, :fips => ByRow(==(01001)), skipmissing=true)
CSV.write(joinpath(root, "CountyAnnualTemperature1950to2012.csv"), df2_ex)


# Add cubic spline basis variables
knots = [0 8 16 24 32]
df_temp1 = DataFrame(load(joinpath(root, county_fn)))
df_temp1[!, :tAvg] = (df_temp1[!, :tMax] .+ df_temp1[!, :tMin]) ./ 2
df_temp1 = cubic_spline_vars(df_temp1, :tAvg, knots)

# sum over year, then average over all grid points
df_temp = aggregate_df(df_temp1)
df_temp[20, [:year, :splineC1, :splineC2, :splineC3, :splineC4]]

# Open linear piecewise stata results
df3 = DataFrame(load(joinpath(root, "bestLinearModel", "corn_year1950_2020_month3_8.dta")))






#==============================================================================
                        1.1 TEMPERATURE AGGREGATION
==============================================================================#
#=================================================
1.1.1 Construct 4 temperature response variables
=================================================#





#=================================================
1.1.2 Aggregate to year-county
- sum each new temperaature variable in each grid point over each year
- average tMin, tMax, tAvg, perc in each gridpoint over each year
- take a simple average over all grid points in the county to get county-year level
=================================================#















#==============================================================================
                        1.2 US Climate Impacts
==============================================================================#
#=================================================
1.2.0 Merge Income and Employment onto climate data
=================================================#







#=================================================
1.2.1 emp_farm vs binned temperature
Explore the relationship between log tranformed emp_farm 
and the vector of binned temperature controls.
Include additional controls for county FE, year FE.
Provide an interpretation of the coefficient of the 32+ bin.
=================================================#







#=================================================
1.2.2 per capita farm prop income vs cubic spline temperature
Explore the relationship between log tranformed 
inc_farm_prop_income/population (i.e. log(per capita 
farm prop income) using the restricted cubic spline. 
Include additional controls for county FE and year FE.
Plot the predicted marginal effects with associated 
confidence intervals, and compare to the binned 
temperature response function above.
=================================================#







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
df1 = convert_float32(pollution, [:dlhouse, :dgtsp, :pop7080])
# Reg: dlhouse on dgtsp, weight with pop7080
reg211 = regression(df1, @formula(dlhouse ~ dgtsp); weights="pop7080", filestub="reg211")

# Reg: dlhouse on dgtsp, weight with pop7080,
# controlling for ddens dwhite dfeml dage65
# dhs dcoll durban dpoverty dvacant downer
# dplumb dtaxprop (depend?)
cols212 = [:dgtsp, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
df212 = convert_float32(pollution, [[:dlhouse, :pop7080]; cols212])
formula212 = (term(:dlhouse) ~ sum(term.(Symbol.(cols212))))
reg212 = regression(df212, formula212; weights="pop7080", filestub="reg212")
reg212[:reg]

#! Save output to CSV, then use web to create latex tables for now











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




















