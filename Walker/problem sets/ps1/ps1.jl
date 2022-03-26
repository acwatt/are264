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
root = "/media/a/E/Programming/github/are264/Walker/problem sets/ps1/"
zip_fn = "Walker-ProblemSet1-Data.zip"
county_fn = "fips1001.dta"

# Extract county fips file from zip
extract_file_from_zip(root, zip_fn, county_fn)
# Read Autauga temperature data
df1 = DataFrame(load(string(root, county_fn)))





#==============================================================================
                                    ADMIN FUNCTIONS
==============================================================================#
"""Extracts file_name from zip folder zip_name in directory root_dir."""
function extract_file_from_zip(root_dir, zip_name, file_name)
    zip_path = string(root_dir, zip_name)
    save_path = string(root_dir, file_name)
    # Open zip file
    zarchive = ZipFile.Reader(zfile)
    # Find file_name in zip archive
    files = [f.name for f in zarchive.files]
    idx = findfirst(x->x==file_name, files)
    # Save file_name to disk
    write(save_path, read(zarchive.files[idx]))
end




#==============================================================================
                                    ANALYSIS FUNCTIONS
==============================================================================#









#==============================================================================
                                    ANALYSIS
==============================================================================#


