using ArgParse, CSV, DataFrames, Tables, Logging, Dates, LoggingExtras, Tidier, Printf, Distributions

struct Duration
    value::String
end

Base.show(io::IO, d::Duration) = print(io, d.value)


function get_time()::Duration
    global start_time
    total_seconds = round(Int, (now() - start_time).value / 1000)

    minutes = total_seconds ÷ 60
    seconds = total_seconds % 60

    return Duration("$(minutes):$(lpad(seconds, 2, "0"))")
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--in"
            help = "Input sumstats filename"
            arg_type = String
        "--out"
            help = "Output dir"
            arg_type = String 
        "--info-min"
            help = "Lower Info threshold"
            default = 0.9
            arg_type = Real
        "--merge-hapmap"
            help = "Toggle whether to merge with hapmap, and only keep overlapping SNPs"
            action = :store_true
        "--keep-biallelic"
            help = "Toggle whether to remove non-biallelic SNPS"
            action = :store_true
    end

    return parse_args(s)
end


function rename_columns(df)
    replacements_df = CSV.File("replacements.csv", delim=';') |> DataFrame

    # build a dictionary mapping old names to new ones
    replacements = Dict(lowercase(string(row.old)) => string(row.new) for row in eachrow(replacements_df))

    old_names = names(df)
    new_names = [get(replacements, lowercase(string(name)), string(name)) for name in old_names]
    pairs = Pair.(old_names, Symbol.(new_names))
    rename!(df, pairs)

    return df
end


function read_sumstats(file_name::String) 
    df = CSV.File(file_name, limit = 1000) |> 
        DataFrame

    if size(df, 2) == 1
        df = CSV.File(file_name, delim=' ', ignorerepeated=true, missingstring="NA") |> DataFrame
    end

    df = df |> rename_columns

    @info "Sumstats read." SNPs = nrow(df) time = get_time()

    df
end


function filter_info(df, info_min::Real)
    if !("Info" in names(df))
        @info "No Info column found in file"
        return(df)
    end

    original_size = size(df, 1)
    filter!(row -> row[:Info] >= info_min, df)

    num_invalid = count(row -> row[:Info] < 0 || row[:Info] > 1, eachrow(df))
    prop_invalid = percent(num_invalid / nrow(df))
    if num_invalid > 0
        @warn "$prop_invalid of rows with Info values outside the range [0, 1]"
    end

    num_removed = original_size - size(df, 1)
    @info "Filter Info < $info_min." "SNPs removed" = num_removed SNPs = nrow(df) time = get_time()

    df
end

function check_cols(df)
    required_columns = ["n", "A1", "A2", "MarkerName"]
    missing_columns = setdiff(required_columns, names(df))


    if !isempty(missing_columns)
        missing_columns_str = join(missing_columns, ", ")
            @warn "The following columns are missing from the input file: $missing_columns_str"
    end

    required_cols = [["Zb"], ["Zp"], ["Z"], ["Direction", "Stat"], ["Direction", "P"]]

    for (i, cols) in enumerate(required_cols)
        if all(x -> x in names(df), cols)
            break
        elseif i == length(required_cols)
            @warn "The file does not contain Z or Direction,Stat, or Direction,P"
        end
    end

    df
end


function merge_hapmap(df, merge_hapmap::Bool)
    if !merge_hapmap
        if !("SNP" in names(df))
            if "Chr" in names(df) && "Pos" in names(df)
                df = transform(df, [:Chr, :Pos] => ByRow((Chr, Pos) -> "$Chr:$Pos") => :SNP)
            else
                @warn "Data not merged with hapmap. SNP positions missing"
            end
        end
            
        return(df)
    end

    merge_df = CSV.File("files/hapmap3.snps") |> DataFrame
    original_size = size(df, 1)
    merged = innerjoin(df, select(merge_df, [:MarkerName, :SNP]), on = :MarkerName)
    new_size = size(merged, 1)
    num_removed = original_size - new_size

    @info "SNP positions merged from hapmap. SNPs removed if not present in hapmap." "SNPs removed" = num_removed SNPs = new_size time = get_time()
    if num_removed / original_size > 0.3333
        @warn "more than 33% of SNPs removed"
    end

    merged
end

function split_alleles(df)
    if "Alleles" in names(df)
        start_time = now()
        alleles_split = split.(df.Alleles, '/')
        df.A1 = first.(alleles_split)
        df.A2 = last.(alleles_split)
        @info "A1 and A2 created from Alleles column." time = get_time()
    end
    df
end


function remove_alleles(df, keep_biallelic::Bool)
    original_size = size(df, 1)

    filter!(row -> !((row[:A1] == "A" && row[:A2] == "T") || 
                     (row[:A1] == "T" && row[:A2] == "A") || 
                     (row[:A1] == "G" && row[:A2] == "C") || 
                     (row[:A1] == "C" && row[:A2] == "G")), df)

    num_removed = original_size - size(df, 1)
    @info "Removing ambiguous alleles (A/T and G/C)." "SNPs removed" = num_removed SNPs = (size(df, 1)) time = get_time()

    if !(keep_biallelic)
        original_size = nrow(df)
        df = filter(row -> row.A1 in ["A", "T", "G", "C"] && row.A2 in ["A", "T", "G", "C"], df)
        new_size = nrow(df)
        num_removed = original_size - new_size
        @info "Removing non-biallelic SNPs." "SNPs removed" = num_removed SNPs = new_size time = get_time()
    end

    df  
end


function remove_duplicates(df)
    original_size = nrow(df)
    unique!(df, :MarkerName)
    num_removed = original_size - nrow(df)

    @info "Remove duplicate markers." "SNPs removed" = num_removed SNPs = nrow(df) time=time()

    df
end


function calculate_statistics(df)
    if "NeglogPval" in names(df) && !("Pval" in names(df))
        df = transform(df, :NeglogPval => (x -> 10 .^ -x) => :Pval)
        @info "Pval calculated from NeglogPval." time = get_time()
    end
    
    if "Pval" in names(df)
        df = transform(df, :Pval => (Pval -> sqrt.(quantile.(Chisq(1), 1 .- Pval))) => :Zp)
        @info "Zp calculated from Pval." time = get_time()
    end

    if "ndiv2" in names(df) && !("n" in names(df))
        df.n = df.ndiv2 .* 2
        @info "N calculated from ndiv2" time = get_time()
    end

    if "n" ∉ names(df) && "neff" in names(df)
        rename!(df, :neff => :n)
        @info "neff used for n"
    end

    if "n" ∉ names(df) && "nmax" in names(df)
        rename!(df, :nmax => :n)
        @info "nmax used for n"
    end

    df
end


percent(n) = @sprintf("%.1f%%", n * 100)


function check_statistics(df)
    if "Pval" in names(df)
        original_size = nrow(df)
        filter!(row -> 0 < row[:Pval] < 1, df)
        num_removed = original_size - nrow(df)
        if num_removed > 0
            @warn "Impossible p values detected." "SNPs removed" = num_removed SNPs = nrow(df) time = time()
        end
    end

    if "Zp" in names(df) && "Z" in names(df)
        diff = (abs.(df.Zp) .- abs.(df.Z))
        large_diff_rows = diff .> 0.05
        large_diff_df = df[large_diff_rows, :]
        println(first(large_diff_df, 10))
        num_large_diff = sum(large_diff_rows)
        prop_diff = percent(num_large_diff / nrow(df))
        if num_large_diff > 0
            @warn "$prop_diff of rows have more than 5% difference between |Zp| and |Z|."
        end        
    end

    df
end


function identify_effect(df)
    m = median(df.Effect)
    if abs(m) < 0.3
    rename!(df, :Effect => :Beta)
        @info "Effect is beta"
        if "Stderr" in names(df)
            df.Z = df.Beta ./ df.Stderr
            @info "Z calculated from Beta and Stderr" time = get_time()
        else
            @warn "As Stderr is missing, Z can not be calculated"
        end
    elseif abs(m) > 0.7 && abs(m) < 1.3
        rename!(df, :Effect => :OR)
        df.Beta = log.(df.OR)
        @info "Effect is OR. Beta calculated." time = get_time()
        
        if "Stderr" in names(df)
            df.Z = df.Beta ./ df.Stderr
            @info "Z calculated from Beta and Stderr" time = get_time()
        else
            @warn "As Stderr is missing, Z can not be calculated"
        end
    else
        @assert abs(m) < 0.1 "Effect column has unusual values. Check data."
    end

    df
end
 

function allele_columns(df)
    @assert ("A1" in names(df) && "A2" in names(df)) "Allele columns A1 and A2 must be present."
    df
end


function main()
    args = parse_commandline()
    inp = basename(args["in"])
    output_file_path = joinpath(args["out"], "$(inp).processed")
    use_hapmap = get(args, :merge_hapmap, false)
    keep_biallelic = get(args, :keep_biallelic, false)
    
    log_file = joinpath(args["out"], "$(inp).log")
    file_logger = Logging.SimpleLogger(open(log_file, "a"))
    console_logger = Logging.ConsoleLogger(stderr)
    Logging.global_logger(LoggingExtras.TeeLogger(console_logger, file_logger))
    @info "Program started"
    global start_time = Dates.now()

    df = read_sumstats(args["in"])
    println(first(df, 3))

    df = @chain df begin
        filter_info(args["info-min"])
        merge_hapmap(use_hapmap)
        split_alleles()
        allele_columns()
        remove_alleles(keep_biallelic)
        calculate_statistics()
        identify_effect()
        check_statistics()
        check_cols()

    end

    println(first(df, 3))

# Specify the path to your directory
path = "../rawdata/sumstats1"

# Get a list of all files in the directory
all_files = readdir(path)

# Filter out files that end with .zip
non_zip_files = filter(f -> !endswith(f, ".zip"), all_files)

# Convert to a DataFrame and save as a CSV
df = DataFrame(study = non_zip_files)
CSV.write("studies.csv", df)

    @info "Program finished." time = get_time()
end

main()

