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
        return DataFrame(load(save_path))
    else
        println("Downloading $save_path from url.")
        download(url, joinpath(root, temp_fn))
        println("Loading $save_path from file.")
        return DataFrame(load(save_path))
    end
end


"""Return a dataframe, with types Float64 and no missing for given columns."""
function convert_float32(df, cols)
    df1 = dropmissing(df, cols=cols; disallowmissing=true)
    df1 = transform!(df1, cols .=> ByRow(Float64), renamecols=false)
    return df1
end


"""Return a merged dataframe for problem 2 with necessary variables created."""
function merge_prep_est_data()
    # Extract and read employment file from zip
    df_employ = df_from_zip(root, zip_fn, "reis_combine.dta")
    df_employ[!, :fips] = parse.(Int32, df_employ[!, :fips])
    # Download and load example county file to compare variables
    df_temp = df_from_url(example_url, joinpath(root, temp_fn))
    # drop first row because it's missing from dta transformation
    df_temp = last(df_temp, nrow(df_temp)-1)
    # Join onto climate data
    df = innerjoin(df_temp, df_employ, on=[:fips, :year])
    # Make year and fips categegorical for regressions
    df[!,:year] = categorical(df[!,:year])
    df[!,:fips] = categorical(df[!,:fips])
    # Convert the year into a Date for rolling average calculation
    df[!, :Date] = Date.(convert(Vector{Int}, df[!,:year]))
    # Create log(farm employees in county) and replace log(0) with missing
    df[!, :emp_farm_ln] = replace(log.(df[!, :emp_farm]), -Inf => missing)
    # Create log tranformed inc_farm_prop_income/population 
    logfn(x) = x===missing ? missing : (x ≤ 0 ? missing : log(x))
    df[!, :inc_farm_prop_inc_lpc] = logfn.(df[!,:inc_farm_prop_income] ./ df[!,:pop_population])
    # 40,385 missing
    return df
end
