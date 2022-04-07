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
    idx = findfirst(x->x==file_name, files)
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




#==============================================================================
                            TEMPERATURE FUNCTIONS
==============================================================================#
"""Return temperature threshold after converting to -1,1 interval (for sine function)"""
function convert_temp_threshold(Tmax, Tmin, Thresh)
    proportion = (Thresh - Tmin) / (Tmax - Tmin)
    return (2*proportion - 1)
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
        return (π - 2*asin(α)) / (2π)
    end
end


"""Degree day calculation from Snyder 1985, using integral of sine."""
function degree_day_calc(MAX, MIN, THR)
    M = (MAX + MIN) / 2  # average temp
    W = (MAX - MIN) / 2  # half temp range
    θ = asin( (THR - M) / W )  # radian between -π/2 and π/2 where sine crosses THR
    DD = ( (M - THR)(π/2 - θ) + W*cos(θ) ) / π
    return DD
end


"""Return degree days above threshold temperature for single day based on max and min temperatures."""
function degree_day_single_day(Tmax, Tmin, Threshold)
    if Threshold ≥ Tmax
        return 0
    elseif Threshold ≤ Tmin
        return (Tmax + Tmin)/2 - Threshold
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
    if (j < 1) | (j > k-1)
        println("Out of Bounds (j=$j): j is < 1 or > k-1. Not a valid variable.")
    elseif j == 1
        return x
    else
        xⱼ = m.(x .- t[j-1]).^3 .-
             m.(x .- t[k-1]).^3 * (t[k]-t[j-1])/(t[k]-t[k-1]) +
             m.(x .- t[k]).^3 * (t[k-1]-t[j-1])/(t[k]-t[k-1])
        # Divide by the scale factor from the stata formula (range of the knots)
        return xⱼ / (t[k] - t[1])^2
    end
end


"""Return jth restricted cubic spline variable using the list of knots t.
From stata's mkspline documentation
"""
function cubic_spline_var_stata(V, i, k)
    n = length(k)
    if (i < 1) | (i > n-1)
        println("Out of Bounds (i=$i, n=$n): i < 1 or i > n-1. Not a valid variable.")
    elseif i == 1
        return V
    else
        Vᵢ = (m.(V .- k[i-1]).^3 .- (k[n]-k[n-1])^(-1) * 
                ( m.(V .- k[n-1]).^3 * (k[n]-k[i-1]) .- (V .- k[n]).^3 * (k[n-1]-k[i-1]) ) ) ./ 
             (k[n]-k[1])^2
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
function cubic_spline_vars(df, varname, knots, newbasename = "splineC")
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
    df[!,"date"] = Date(1960,1,1) + Day.(df.dateNum)
    df[!,"year"] = Year.(df.date)
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
                               ANALYSIS FUNCTIONS
==============================================================================#





#==============================================================================
                                    ANALYSIS
==============================================================================#


