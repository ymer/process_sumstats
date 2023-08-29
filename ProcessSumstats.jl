using ArgParse, CSV, DataFrames, Logging, LoggingExtras, Distributions, CodecZlib, Formatting

include("Utils.jl"); include("NameChanges.jl")


macro log(msg, df)
    return quote
        formatted_snps = format(nrow($(esc(df))), commas = true)
        @info $(esc(msg)) * "\nSNPs: " * formatted_snps * "\ntime: " * string(get_time())
    end
end


function parse_commandline()
    s = ArgParseSettings(
            description = "This program is used to process gwas sumstat files in various formats to a uniform format.",
            epilog = 
        "The ID of a SNP can be in two formats, called 'SNP' and 'Rsid'. SNP is in the form of chr:position, while Rsid is a rs-number. Examples:\n
        'SNP': 10:157346\n
        'Rsid': rs162373\n
        \n
        A .bim file is a file in the plink .bim format. It has the following unlabelled columns: Chr, SNP/Rsid, GPos, Pos, A1, A2\n
        \n
        A .snpid file is a custom file format. It has the following labelled columns: SNP, Rsid, A1, A2."
       )

    @add_arg_table! s begin
        "--in"
            help = "Input sumstats filename"
            arg_type = String
            required = true
        "--outdir"
            help = "Output directory"
            arg_type = String 
            required = true
        "--info-min"
            help = "Minimum Info thresholdi"
            default = 0.9
            arg_type = Real
        "--keep-biallelic"
            help = "Non-biallelic SNPS are removed by default. Use this option to keep them in."
            action = :store_true
        "--filter"
            help = "Path(s) to one or more .bim / .snpid files (separated by comma). Keep only variants that are present here."
            #arg_type = String
            default = ""
        "--write-snpid"
            help = "Write a .snpid file"
            action = :store_true
        "--check-build"
            help = "Path to a .snpid file. Used to check if the sumstats file SNP positions are in the same build."
            default = ""
            arg_type = String
        "--pos-from-snpid"
            help = "Path to a .snpid file. Used to find the SNP position if the sumstat file only has Rsid."
            default = ""
        "--head"
            help = "Read only the first n columns"
            default = -1
            arg_type = Int64
        "--convert-to-tabs"
            help = "Convert spaces in file to tabs"
            action = :store_true
    end

    return parse_args(s)
end


function read_sumstats(file_name::String)
    function rename_columns(df)
        old_names = names(df)
        new_names = [get(replacements, lowercase(string(name)), string(name)) for name in old_names]
        pairs = Pair.(old_names, Symbol.(new_names))
        rename!(df, pairs)

        if :Rsid ∈ df && :MarkerName ∈ df && :SNP ∉ df
            rename!(df, :MarkerName => :SNP)
            @info "MarkerName renamed to SNP"
        end

        if :Rsid ∉ df && :MarkerName ∈ df
            first_value = df[1, :MarkerName]
            if first_value isa AbstractString && startswith(first_value, "rs")
                rename!(df, :MarkerName => :Rsid)
                @info "MarkerName renamed to Rsid"
            end
        end

        df
    end


    function convert_spaces_to_tabs(file_name, new_filename)
        _, ext = splitext(file_name)

        if ext == ".gz"
            stream = GzipDecompressorStream(open(file_name, "r")) 
            contents = read(stream, String)
            close(stream)
        else
            contents = read(file_name, String)
        end

        contents = replace(contents, " " => "\t")
        write(new_filename, contents)
        @info "Spaces converted to tabs"
    end
     
    if args["convert-to-tabs"]
        new_filename = file_name * "_tabs"
        convert_spaces_to_tabs(file_name, new_filename)
        file_name = new_filename
    end

    n = args["head"]
    if n != -1
        @info "Only first $n rows read from sumstats file"
        df = CSV.File(file_name, limit = n, missingstring = ["NA", "N/A"]) |> DataFrame
    else
        df = CSV.File(file_name, missingstring = ["NA", "N/A"]) |> DataFrame
    end
    
    if size(df, 2) < 3
        new_filename = file_name * "_tabs"
        convert_spaces_to_tabs(file_name, new_filename)
        df = read_sumstats(new_filename)
    end

    df = df |> rename_columns

    @log "Sumstats read" df

    df
end


function check_data_types(df)
    if :Effect ∈ names(df)
        non_numeric = [x for x in df.Effect if isa(x, String) && !isnumber(x)]
        if !isempty(non_numeric)
            println.(non_numeric[1:5])
            @error "Effect column has the above values which are not numeric"
        end
    end
end


function filter_info(df)
    info_min = args["info-min"]
    if :Info ∉ df
        @info "No Info column found in file"
        return(df)
    end

    invalid_df = filter(row -> row[:Info] < 0 || row[:Info] > 1.02, df)
    if nrow(invalid_df) > 0
        @warn "$(nrow(invalid_df)) of rows with Info values outside the range [0, 1.02]"
        print_sample(invalid_df, [:SNP, :Rsid, :Info], "Invalid Info values")
    end

    filter!(row -> row[:Info] >= info_min, df)
    @log "Remove SNPs with Info < $info_min." df

    df
end


function create_SNP(df)
    args["pos-from-snpid"] != "" && return df
    if :SNP ∉ df
        if :Chr in df && :Pos in df
            df = transform(df, [:Chr, :Pos] => ByRow((Chr, Pos) -> "$Chr:$Pos") => :SNP)
            @info "SNP column created from Chr and Pos"
        else
            @warn "SNP column not found. Neither is Chr + Pos"
        end
    end
    df
end


function pos_from_snpid(df)
    fn = args["pos-from-snpid"]
    fn == "" && return df
    
    snpid = CSV.File(fn) |> DataFrame
    merged = innerjoin(df, snpid, on = [:Rsid, :A1, :A2])

    @log "SNP position acquired from reference. " * "Reference file = $fn" merged
    
    merged
end


function check_build(df)

    if args["check-build"] == ""
        @info "SNP positions not checked. Could be from unknown build."
        return(df)
    end

    @assert :Rsid ∈ df "Rsid must be in dataframe to check with reference."
    @assert :SNP ∈ df "SNP must be in dataframe to check with reference."

    snpid = CSV.File(args["check-build"]) |> DataFrame
    rename!(snpid, :SNP => :SNP_ref)

    original_size = size(df, 1)
    merged = innerjoin(df, select(snpid, [:SNP_ref, :Rsid]), on = :Rsid, matchmissing = :notequal)
    new_size = size(merged, 1)

    @info "Build check" "SNPs found in reference" = Percent(new_size / original_size) time = get_time()
    
    not_equal_df = merged[merged.SNP .!= merged.SNP_ref, :]
    if nrow(not_equal_df) > 0
        @warn "Not all SNP positions matched with reference" "SNPs not matched" = nrow(not_equal_df)
        cols_to_select = intersect([:SNP, :SNP_ref, :Rsid], propertynames(not_equal_df))
        selected_df = select(not_equal_df, cols_to_select)
        print_header(selected_df, "SNPs not matching reference")
    else
        @info "All SNP positions match the reference."
    end
end


function filter_ids(df)
    args["filter"] == "" && return df

    for fn in split(args["filter"], ',')
        ext = splitext(fn)[2]
        if ext == ".bim"
            ref = CSV.File(fn, header=["Chr", "Rsid", "GPos", "Pos", "A1", "A2"]) |> DataFrame
            ref[!, "SNP"] = string.(ref[!, "Chr"]) .* ":" .* string.(ref[!, "Pos"])
        elseif ext == ".snpid"
            ref = CSV.File(fn) |> DataFrame
        else
            @error "Unsupported file extension" ext
            return df
        end

        df = innerjoin(df, select(ref, [:SNP, :A1, :A2]), on = [:SNP, :A1, :A2])

        @log "Filtered with $fn" df
    end
    
    df
end



function split_alleles(df)
    if :Alleles ∈ df
        alleles_split = split.(df.Alleles, '/')
        df.A1 = first.(alleles_split)
        df.A2 = last.(alleles_split)
        @info "A1 and A2 created from Alleles column." time = get_time()
    end
    df
end


function filter_alleles(df)
    @assert :A1 ∈ df && :A2 ∈ df "Allele columns A1 and A2 must be present."
    
    df.A1 = uppercase.(df.A1)
    df.A2 = uppercase.(df.A2)
    
    pairs = [("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")]
    filter!(row -> !any([row[:A1] == p[1] && row[:A2] == p[2] for p in pairs]), df)

    @log "Removing ambiguous alleles (A/T and G/C)" df

    if !(args["keep-biallelic"])
        df = filter(row -> row.A1 in ["A", "T", "G", "C"] && row.A2 in ["A", "T", "G", "C"], df)
        @log "Removing non-biallelic SNPs" df
    end

    df  
end


function remove_duplicates(df)
    unique!(df, :Rsid)
    unique!(df, :SNP)
    @log "Removing duplicate markers" df
    df
end


function identify_effect(df)
    if :Effect ∉ df
        @assert :Z ∈ df "Neither effect nor Z found"
        @info "No Effect column, but Z is present"
        return(df)
    end

    m = median(df.Effect)
    if abs(m) < 0.3
        rename!(df, :Effect => :Beta)
        @info "Effect is Beta"
    elseif m > 0.7 && m < 1.3
        rename!(df, :Effect => :OR)
        df.Beta = log.(df.OR)
        @info "Effect is OR. Beta calculated." time = get_time()
    else
        error("Effect column has unusual values. Check data.")
    end
    
    if :Stderr ∈ df
        df.Z = df.Beta ./ df.Stderr
        @info "Z calculated from Beta and Stderr" time = get_time()
    else
        df.Direction = ifelse.(df.Beta .> 0, 1, -1)
        @warn "As Stderr is missing, Z can not be calculated. Direction calculated from Beta."
    end

    df
end


function drop_missing(df)
    df = dropmissing(df, [:SNP, :Rsid, :A1, :A2, :Effect])
    @log "Drop rows with missing values" df
    df
end


function calculate_statistics(df)
    if :NeglogP in propertynames(df) && !(:P in propertynames(df))
        df = transform(df, :NeglogP => (x -> 10 .^ -x) => :P)
        @info "P calculated from NeglogP." time = get_time()
    end
    
    if :P ∈ df
        df = transform(df, :P => (P -> sqrt.(quantile.(Chisq(1), 1 .- P))) => :Zp)
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
    if :P ∈ df
        original_size = nrow(df)
        invalid_rows = filter(row -> row[:P] < 0 || row[:P] > 1, df)
        filter!(row -> 0 <= row[:P] <= 1, df)
        num_removed = original_size - nrow(df)
        if num_removed > 0
            @warn "Impossible p values detected." "SNPs removed" = num_removed SNPs = nrow(df) time = get_time()
            cols_to_select = intersect([:SNP, :Rsid, :P], propertynames(large_diff_df))
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
            cols_to_select = intersect([:SNP, :Rsid, :P, :OR, :Beta, :OR, :Z, :Zp], propertynames(large_diff_df))
            selected_df = select(large_diff_df, cols_to_select)
            print_header(selected_df, "Large diff between Pz and Z")
        end        
    end

    df
end


function select_cols(df)
    missing_columns = setdiff([:SNP, :n, :A1, :A2], propertynames(df))

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

    select(df, intersect([:SNP, :SNP_ref, :Rsid, :A1, :A2, :Stat, :P, :OR, :Beta, :Z, :Direction, :n], propertynames(df)))
end


function write_output(df)
    fn = basename(args["in"])
    if args["write-snpid"] != ""
        CSV.write(joinpath(args["outdir"], fn * ".snpid"), df[:, [:SNP, :Rsid, :A1, :A2]], delim = "\t")
    end
    CSV.write(joinpath(args["outdir"], fn), df, delim = "\t")
end


function create_logger(log_file)
    file_logger = Logging.SimpleLogger(open(log_file, "a"))
    console_logger = Logging.ConsoleLogger(stderr)
    Logging.global_logger(LoggingExtras.TeeLogger(console_logger, file_logger))
    global start_time = Dates.now()
end



function print_sample(df, cols, title)
    print_header(select(df, intersect(cols, propertynames(df))), title)
end
 

function main()
    global args = parse_commandline()
    create_logger(joinpath(args["outdir"], basename(args["in"]) * ".log"))

    @info "Program started"

    df = read_sumstats(args["in"]) |>
        static(print_header, "Original data") |>
        static(check_data_types) |>
        filter_info |>
        create_SNP |>
        pos_from_snpid |>
        static(check_build) |>
        split_alleles |>
        filter_alleles |>
        filter_ids |>
        remove_duplicates |>
        drop_missing |>
        identify_effect |>
        calculate_statistics |>
        filter_statistics |>
        select_cols |>
        static(write_output) |>
        static(print_header, "After processing")

    @info "Program finished" time = get_time()
end

main()

