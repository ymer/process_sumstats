using ArgParse, CSV, DataFrames, Logging, LoggingExtras, Distributions

include("Utils.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--in"
            help = "Input sumstats filename"
            arg_type = String
            required = true
        "--out"
            help = "Output dir"
            arg_type = String 
            required = true
        "--info-min"
            help = "Minimum Info threshold"
            default = 0.9
            arg_type = Real
        "--keep-biallelic"
            help = "Toggle whether to remove non-biallelic SNPS"
            action = :store_true
        "--reference"
            help = "Path to a reference data set"
            default = ""
            arg_type = String
        "--check-with-reference"
            help = "Checks if the build is hg19 using a reference data set"
            action = :store_true
        "--filter-with-reference"
            help = "Keep only variants that are present in the reference data set."
            action = :store_true
    end

    return parse_args(s)
end


function rename_columns(df)
    replacements_df = CSV.File("replacements.csv", delim=';') |> DataFrame

    replacements = Dict(lowercase(string(row.old)) => string(row.new) for row in eachrow(replacements_df))

    old_names = names(df)
    new_names = [get(replacements, lowercase(string(name)), string(name)) for name in old_names]
    pairs = Pair.(old_names, Symbol.(new_names))
    rename!(df, pairs)

    return df
end


function read_sumstats(file_name::String) 
    df = CSV.File(file_name, limit = 1000) |> DataFrame
    
    if size(df, 2) == 1
        df = CSV.File(file_name, delim=' ', ignorerepeated=true, missingstring="NA") |> DataFrame
        @warn "Space used as separator instead of tab"
    end
 
    df = df |> rename_columns

    @info "Sumstats read." SNPs = nrow(df) time = get_time()

    df
end


function filter_info(df)
    info_min = args["info-min"]
    if :Info ∉ df
        @info "No Info column found in file"
        return(df)
    end

    original_size = size(df, 1)
    filter!(row -> row[:Info] >= info_min, df)

    num_invalid = count(row -> row[:Info] < 0 || row[:Info] > 1, eachrow(df))
    prop_invalid = Percent(num_invalid / nrow(df))
    if num_invalid > 0
        @warn "$prop_invalid of rows with Info values outside the range [0, 1]"
    end

    num_removed = original_size - size(df, 1)
    @info "Remove SNPs with Info < $info_min." "SNPs removed" = num_removed SNPs = nrow(df) time = get_time()

    df
end


function use_reference(df)

    if args["check-with-reference"] || args["filter-with-reference"]
        @assert args["reference"] != "" "--reference must be supplied if you want to filter or check with reference."
    else
        @info "SNP locations not checked with hapmap. Could be from unknown build."
        return(df)
    end

    ref = CSV.File(args["reference"]) |> DataFrame
    original_size = size(df, 1)
    merged = innerjoin(df, select(ref, [:Rsid, :SNP_ref]), on = :Rsid)
    new_size = size(merged, 1)

    @info "Merged with reference." time = get_time() "SNPs found in reference" = Percent(new_size / original_size)
    
    #@assert SNP ∈ df && SNP_ref ∈ df "SNP and SNP_ref must be present in DataFrame"
    df
end


function split_alleles(df)
    if :Alleles ∈ df
        start_time = now()
        alleles_split = split.(df.Alleles, '/')
        df.A1 = first.(alleles_split)
        df.A2 = last.(alleles_split)
        @info "A1 and A2 created from Alleles column." time = get_time()
    end
    df
end


function filter_alleles(df)
    original_size = size(df, 1)

    filter!(row -> !((row[:A1] == "A" && row[:A2] == "T") || 
                     (row[:A1] == "T" && row[:A2] == "A") || 
                     (row[:A1] == "G" && row[:A2] == "C") || 
                     (row[:A1] == "C" && row[:A2] == "G")), df)

    num_removed = original_size - size(df, 1)
    @info "Removing ambiguous alleles (A/T and G/C)." "SNPs removed" = num_removed SNPs = (size(df, 1)) time = get_time()

    if !(args["keep-biallelic"])
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
    unique!(df, :Rsid)
    num_removed = original_size - nrow(df)

    @info "Remove duplicate markers." "SNPs removed" = num_removed SNPs = nrow(df) time=time()

    df
end


function calculate_statistics(df)
    if :NeglogPval in propertynames(df) && !(:Pval in propertynames(df))
        df = transform(df, :NeglogPval => (x -> 10 .^ -x) => :Pval)
        @info "Pval calculated from NeglogPval." time = get_time()
    end
    
    if :Pval ∈ df
        df = transform(df, :Pval => (Pval -> sqrt.(quantile.(Chisq(1), 1 .- Pval))) => :Zp)
        @info "Zp calculated from Pval." time = get_time()
    end

    if :ndiv2 ∈ df && :n ∉ df
        df.n = df.ndiv2 .* 2
        @info "N calculated from ndiv2" time = get_time()
    end

    if :n ∉ df && :neff ∈ df
        rename!(df, :neff => :n)
        @info "neff used for n"
    end

    if :n ∉ df && :nmax ∈ df
        rename!(df, :nmax => :n)
        @info "nmax used for n"
    end

    df
end


function filter_statistics(df)
    if :Pval ∈ df
        original_size = nrow(df)
        invalid_rows = filter(row -> row[:Pval] < 0 || row[:Pval] > 1, df)
        filter!(row -> 0 <= row[:Pval] <= 1, df)
        num_removed = original_size - nrow(df)
        if num_removed > 0
            @warn "Impossible p values detected." "SNPs removed" = num_removed SNPs = nrow(df) time = get_time()
            cols_to_select = intersect([:SNP, :Rsid, :Pval], propertynames(large_diff_df))
            selected_df = select(large_diff_df, cols_to_select)
            print_header(selected_df, "Invalid p-values")
        end
    end

    if :Zp ∈ df && :Z ∈ df
        diff = (abs.(df.Zp) .- abs.(df.Z))
        large_diff_rows = diff .> 0.1
        large_diff_df = df[large_diff_rows, :]
        num_large_diff = sum(large_diff_rows)
        if num_large_diff > 0
            @warn "$num_large_diff rows have more than 0.1 difference between |Zp| and |Z|."
            cols_to_select = intersect([:SNP, :Rsid, :Pval, :OR, :Beta, :OR, :Z, :Zp], propertynames(large_diff_df))
            selected_df = select(large_diff_df, cols_to_select)
            print_header(selected_df, "Large diff between Pz and Z")
        end        
    end

    df
end


function identify_effect(df)
    m = median(df.Effect)
    if abs(m) < 0.3
    rename!(df, :Effect => :Beta)
        @info "Effect is beta"
        if :Stderr ∈ df
            df.Z = df.Beta ./ df.Stderr
            @info "Z calculated from Beta and Stderr" time = get_time()
        else
            @warn "As Stderr is missing, Z can not be calculated"
        end
    elseif abs(m) > 0.7 && abs(m) < 1.3
        rename!(df, :Effect => :OR)
        df.Beta = log.(df.OR)
        @info "Effect is OR. Beta calculated." time = get_time()
        
        if :Stderr ∈ df
            df.Z = df.Beta ./ df.Stderr
            @info "Z calculated from Beta and Stderr" time = get_time()
        else
            @warn "As Stderr is missing, Z can not be calculated"
        end
    else
        @assert abs(m) < 0.3 "Effect column has unusual values. Check data."
    end

    df
end
 

function check_allele_columns(df)
    @assert :A1 ∈ df && :A2 ∈ df "Allele columns A1 and A2 must be present."
end


function create_SNP(df)
    if :SNP ∉ names(df)
        if "Chr" in names(df) && "Pos" in names(df)
            df = transform(df, [:Chr, :Pos] => ByRow((Chr, Pos) -> "$Chr:$Pos") => :SNP)
            @info "SNP column created from Chr and Pos"
        else
            @warn "SNP column not found. Neither is Chr + Pos"
        end
    end
    df
end


function select_cols(df)
    missing_columns = setdiff([:SNP, :n, :A1, :A2, :Rsid], propertynames(df))

    if !isempty(missing_columns)
        missing_columns_str = join(missing_columns, ", ")
            @warn "The following columns are missing from the input file: $missing_columns_str"
    end

    required_cols = [[:Z], [:Direction, :Stat], [:Direction, :P]]

    for (i, cols) in enumerate(required_cols)
        if all(x -> x in propertynames(df), cols)
            break
        elseif i == length(required_cols)
            @warn "The file does not contain Z or Direction,Stat, or Direction,P"
        end
    end

    select(df, intersect([:SNP, :SNP_ref, :Rsid, :A1, :A2, :Stat, :P, :OR, :Beta, :Z, :n], propertynames(df)))
end

function create_logger(log_file)
    file_logger = Logging.SimpleLogger(open(log_file, "a"))
    console_logger = Logging.ConsoleLogger(stderr)
    Logging.global_logger(LoggingExtras.TeeLogger(console_logger, file_logger))
    global start_time = Dates.now()
end

static2(f, arg) = x -> (f(x, arg); x)
static(f) = x -> (f(x); x)

function main()
    global args = parse_commandline()
    inp = basename(args["in"])
    outp = joinpath(args["out"], "$(inp)")
    create_logger(outp * ".log")

    @info "Program started"

    df = read_sumstats(args["in"]) |>
        static2(print_header, "Original data") |>
        filter_info |>
        create_SNP |>
        use_reference |>
        split_alleles |>
        static(check_allele_columns) |>
        filter_alleles |>
        calculate_statistics |>
        identify_effect |>
        filter_statistics |>
        select_cols |>
        static2(print_header, "After processing")

    @info "Program finished." time = get_time()
end

main()

