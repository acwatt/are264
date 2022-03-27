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
Pkg.add(["ZipFile", "StatFiles"])
using ZipFile  # unzip compressed .zip folders
using StatFiles  # read Stata files
using DataFrames


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
example_fn = "CountyAnnualTemperature1950to2012.dta"

# Extract and read county fips file from zip
extract_file_from_zip(root, zip_fn, county_fn)
df1 = DataFrame(load(joinpath(root, county_fn)))
# Extract and read example county file to compare variables
extract_file_from_zip(root, zip_fn, example_fn)
df2 = DataFrame(load(joinpath(root, example_fn)))





#==============================================================================
                                    ADMIN FUNCTIONS
==============================================================================#
"""Extracts file_name from zip folder zip_name in directory root_dir."""
function extract_file_from_zip(root_dir, zip_name, file_name)
    zip_path = joinpath(root_dir, zip_name)
    save_path = joinpath(root_dir, file_name)
    # Open zip file
    zarchive = ZipFile.Reader(zip_path)
    # Find file_name in zip archive
    files = [f.name for f in zarchive.files]
    idx = findfirst(x->x==file_name, files)
    # Save file_name to disk
    write(save_path, read(zarchive.files[idx]))
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








#==============================================================================
                               ANALYSIS FUNCTIONS
==============================================================================#





#==============================================================================
                                    ANALYSIS
==============================================================================#


