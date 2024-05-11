using ArgParse, CSV, DataFrames, Logging, LoggingExtras, Distributions, CodecZlib, Formatting

include("Utils.jl"); include("NameChanges.jl")


function parse_commandline()
    s = ArgParseSettings(
            description = "This program is used to process gwas sumstat files in various formats to a uniform format.",
            epilog = 
        "The ID of a SNP can be in two formats, called 'SNP' and 'RSID'. SNP is in the form of chr:position, while RSID is a rs-number. Examples:\n
        'SNP': 10:157346\n
        'RSID': rs162373\n
        \n
        A .bim file is a file in the plink .bim format. It has the following unlabelled columns: CHR, SNP/RSID, GPOS, POS, A1, A2\n
        \n
        A .snpid file is a custom file format. It has the following labelled columns: SNP, RSID, A1, A2."
       )

    @add_arg_table! s begin
        "in"
            help = "Input sumstats filename"
            arg_type = String
            required = true
        "out"
            help = "Output filename"
            arg_type = String 
            required = true
        "--info-min"
            help = "Minimum Info threshold"
            default = 0.9
            arg_type = Real
        "--keep-biallelic"
            help = "Non-biallelic SNPS are removed by default. Use this option to keep them in."
            action = :store_true
        "--filter"
            help = "Path(s) to one or more .bim / .snpid files (separated by comma). Keep only variants that are present here."
            arg_type = String
            default = ""
        "--write-snpid"
            help = "Write a .snpid file"
            action = :store_true
        "--check-with-ref"
            help = "Check if the positions are the same as in the reference. (Requires --ref.)"
            action = :store_true
        "--pos-from-ref"
            help = "Match on RSID, and use the position from the reference instead. (Requires --ref.)"
            action = :store_true
        "--ref"
            help = "Path to a reference .bim file."
            default = ""
            arg_type = String
        "--head"
            help = "Read only the first n columns"
            default = -1
            arg_type = Int64
        "--convert-to-tabs"
            help = "Convert spaces in file to tabs"
            action = :store_true
        "--convert-to-dots"
            help = "Convert commas in numerical columns to dots"
            action = :store_true
        "--assign-n"
            help = "Assign an n value"
            default = -1
            arg_type = Int64
        "--skip"
            help = "Skip first n rows when reading the file"
            default = -1
            arg_type = Int64
        "--no-error"
            help = "Do not prompt an error is important columns are missing"
            action = :store_true
        "--effect_col"
            help = "The name of the effect column (eg OR / beta) (if unusual)"
            default = nothing
            arg_type = String
        "--a1_col"
            help = "The name of the A1 column (if unusual)"
            default = nothing
            arg_type = String
        "--a2_col"
            help = "The name of the A2 column (if unusual)"
            default = nothing
            arg_type = String
         "--se_col"
            help = "The name of the SE column (if unusual)"
            default = nothing
            arg_type = String
         "--p_col"
            help = "The name of the p-value column (if unusual)"
            default = nothing
            arg_type = String
         "--rsid_col"
            help = "The name of the RSID column (if unusual)"
            default = nothing
            arg_type = String
         "--snp_col"
            help = "The name of the SNP column (if unusual)"
            default = nothing         
         "--n_col"
            help = "The name of the N column (if unusual)"
            default = nothing
         "--rsid-as-snp"
            help = "In the output file, have RSID in the column called SNP"
            action = :store_true
        end

    return parse_args(s)
end


function read_sumstats(file_name::String)
    function rename_columns(df)

        replacements2 = Dict{String, String}()
        for (arg_key, replacement) in [("effect_col", "Effect"), ("p_col", "P"), ("rsid_col", "RSID"), ("a1_col", "A1"), ("a2_col", "A2"), ("se_col", "SE"), ("snp_col", "SNP"), ("n_col", "N")]
            args[arg_key] !== nothing && (replacements2[args[arg_key]] = replacement)
        end
        
        old_names = names(df)
        new_names = [get(replacements2, string(name), get(replacements, lowercase(string(name)), string(name))) for name in old_names]
        unique_new_names = Symbol[]
        seen = Dict{Symbol, Int}()

        for name in new_names
            count = get!(seen, Symbol(name), 0)
            seen[Symbol(name)] += 1
            unique_name = count == 0 ? Symbol(name) : Symbol("$(name)_$(count)")
            push!(unique_new_names, unique_name)
        end

        pairs = Pair.(old_names, unique_new_names)
        rename!(df, pairs)

        if :RSID ∈ df && :MarkerName ∈ df && :SNP ∉ df
            rename!(df, :MarkerName => :SNP)
            @info "MarkerName renamed to SNP"
        end

        if :RSID ∉ df && :MarkerName ∈ df
            first_value = df[1, :MarkerName]
            if first_value isa AbstractString && startswith(first_value, "rs")
                rename!(df, :MarkerName => :RSID)
                @info "MarkerName renamed to RSID"
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
        contents = replace(contents, r"[ \t]+" => "\t")
        #contents = replace(contents, " " => "\t")
        write(new_filename, contents)
        @info "Spaces converted to tabs"
    end

    t1 = now()

    if args["convert-to-tabs"]
        new_filename = file_name * "_tabs"
        convert_spaces_to_tabs(file_name, new_filename)
        file_name = new_filename
    end

    csv_args = Dict{String, Any}("missingstring" => ["NA", "N/A", ".", ""], "comment" => "#")


    head = args["head"]
    if head != -1
        csv_args["limit"] = head
        @info "Only first $head rows read from sumstats file"
    end

    skip = args["skip"]
    if skip != -1
        csv_args["skipto"] = skip + 2
        csv_args["header"] = skip + 1
    end

    if args["convert-to-dots"]
        csv_args["decimal"] = ','
    end

    symbol_csv_args = Dict(Symbol(k)=>v for (k,v) in csv_args)
    df = CSV.File(file_name; symbol_csv_args...) |> DataFrame

    if args["convert-to-dots"]
        for col in names(df)
            if any(x -> occursin(",", string(x)), df[!, col])
                df[!, col] = map(x -> isempty(string(x)) ? missing : parse(Float64, replace(string(x), "," => ".")), df[!, col])
            end
        end
    end
     
    for col in names(df)
        if all(ismissing, df[!, col])
            select!(df, Not(col))
        end
    end
 
    if size(df, 2) < 3
        new_filename = file_name * "_tabs"
        convert_spaces_to_tabs(file_name, new_filename)
        df = read_sumstats(new_filename)
    end
   
    df = df |> rename_columns

    function is_number(s)
        try
            parse(Float64, s)
            true
        catch
            false
        end
    end


    numcols = [:Info, :P, :Effect, :Stderr]
    for col in intersect(numcols, names(df))
        if df[!, col] isa AbstractVector{<:AbstractString} && any(occursin(",", x) for x in df[!, col])
            df[!, col] = parse.(Float64, replace.(df[!, col], "," => "."))
        end
    end

    @log "Sumstats read" file_name df t1
    global initial_snps = nrow(df)

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
        return df
    end

    invalid_df = filter(row -> row[:Info] < 0 || row[:Info] > 1.1, df)
    if nrow(invalid_df) > 0
        @warn "$(nrow(invalid_df)) of rows with Info values outside the range [0, 1.1]"
        print_sample(invalid_df, [:SNP, :RSID, :Info], "Invalid Info values")
    end

    filter!(row -> row[:Info] >= info_min, df)
    @log "Remove SNPs with Info < $info_min." df
    
    check_snps_left(df)
    df
end


function create_SNP(df)
    args["pos-from-ref"] && return df
    if :SNP ∉ df
        if :CHR in df && :POS in df
            df = transform(df, [:CHR, :POS] => ByRow((Chr, POS) -> "$Chr:$POS") => :SNP)
            @info "SNP column created from CHR and POS"
        else
            @warn "SNP column not found. Neither is CHR + POS"
        end
    end
    df
end


function pos_from_snpid(df)
    fn = args["pos-from-snpid"]
    fn == "" && return df
    
    snpid = CSV.File(fn) |> DataFrame
    merged = innerjoin(df, snpid, on = [:RSID, :A1, :A2])

    @log "SNP position acquired from reference. " * "Reference file = $fn" merged
    
    check_snps_left(merged)
    merged
end


function read_bim(fn) 
    bim = CSV.File(fn, header=["CHR", "RSID", "GPOS", "POS", "A1", "A2"]) |> DataFrame
    if all(occursin(":", x) for x in bim[1:5, "RSID"])
        rename!(bim, "RSID" => "SNP")
    else      
        bim[!, "SNP"] = string.(bim[!, "CHR"]) .* ":" .* string.(bim[!, "POS"])
    end

    bim
end

function use_ref(df)

    if !args["check-with-ref"] & !args["pos-from-ref"]
        @info "SNP positions not checked. Could be from unknown build."
        return df
    end

    t1 = now()
    fn = args["ref"]
    snpid = read_bim(fn)
    print_header(snpid)
    print_header(df)
    @log "Reference read"

    if args["pos-from-ref"]
        df = "SNP" in names(df) ? select(df, Not(:SNP)) : df
        merged = innerjoin(df, select(snpid, [:SNP, :RSID]), on = :RSID, matchmissing = :notequal)
        @log "SNP position changed to those from reference" merged
        check_snps_left(merged)
        return merged
    end

    @assert :RSID ∈ df "RSID must be in dataframe to check with reference."
    @assert :SNP ∈ df "SNP must be in dataframe to check with reference."

    rename!(snpid, :SNP => :SNP_ref)

    original_size = size(df, 1)
    merged = innerjoin(df, select(snpid, [:SNP_ref, :RSID]), on = :RSID, matchmissing = :notequal)
    new_size = size(merged, 1)
   
    @info "Build check" "SNPs found in reference" = Percent(new_size / original_size)
    
    not_equal_df = merged[merged.SNP .!= merged.SNP_ref, :]
    if nrow(not_equal_df) > 0
        @warn "Not all SNP positions matched with reference" "SNPs not matched" = nrow(not_equal_df)
        cols_to_select = intersect([:SNP, :SNP_ref, :RSID], propertynames(not_equal_df))
        selected_df = select(not_equal_df, cols_to_select)
        print_header(selected_df, "SNPs not matching reference")
    else
        @info "All SNP positions match the reference."
    end

    df
end


function filter_ids(df)
    args["filter"] == "" && return df

    for fn in split(args["filter"], ',')
        t = now()
        ext = splitext(fn)[2]
        if ext == ".bim"
            filter = read_bim(fn)
        elseif ext == ".snpid"
            filter = CSV.File(fn) |> DataFrame
        else
            @error "Unsupported file extension" ext
            return df
        end

        @log "Filter list read" fn filter t

        df = innerjoin(df, select(filter, [:SNP, :A1, :A2]), on = [:SNP, :A1, :A2])
        @log "Perform filtering" df
    end
    
    check_snps_left(df)
    df
end



function split_alleles(df)
    if :Alleles ∈ df
        alleles_split = split.(df.Alleles, '/')
        df.A1 = first.(alleles_split)
        df.A2 = last.(alleles_split)
        @info "A1 and A2 created from Alleles column."
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

    initial_row_count = size(df, 1)
    filter!(row -> row[:A1] != row[:A2], df)
    size(df, 1) < initial_row_count && @log "Removing equal allele pairs" df

    check_snps_left(df)
    df  
end


function remove_duplicates(df)
    :RSID in df && unique!(df, :RSID)
    :SNP in df && unique!(df, :SNP)
    @log "Removing duplicate markers" df
    check_snps_left(df)
    df
end


function identify_effect(df)
    if :Effect ∉ df
        if :Z in df
            if :P in df 
                if :SE ∉ df
                    @info "Calculate SE from P and Z"
                    df.SE = abs.(df.Z) ./ sqrt.(-2 * log.(df.P))
                    @info "Calculate BETA from SE and Z"
                    df.BETA = df.Z .* df.SE
                end
            end
        else
            @warn "Neither effect nor Z found"
        end
        return df
    end

    m = median(df.Effect)
    
    if abs(m) < 0.3
        rename!(df, :Effect => :BETA)
        @info "Effect is BETA"
    
    elseif m > 0.7 && m < 1.3
        rename!(df, :Effect => :OR)
        if "Direction" in names(df) 
            if Set(df.Direction) == Set([-1, 1])
                df.BETA = log.(df.OR) .* df.Direction * -1
            else
                @warn "Uninterpretable direction"
                df.BETA = log.(df.OR)
            end
        else
            df.BETA = log.(df.OR)
        end
        @info "Effect is OR. BETA calculated."
    else
        error("Effect column has unusual values. Check data.")
    end
    
    if :SE ∈ df
        df.Z = df.BETA ./ df.SE
        @info "Z calculated from BETA and Standard error"
    else
        @warn "As Standard error is missing, Z can not be calculated."
        if !(:Direction in df)
            df.Direction = ifelse.(df.BETA .> 0, 1, -1)
            @info "Direction calculated from BETA"
        end
    end

    df
end


function drop_missing(df)
    cols_to_check = Symbol.(intersect(string.([:SNP, :RSID, :A1, :A2, :Effect, :Z, :BETA]), string.(names(df))))
    df = dropmissing(df, cols_to_check)
    @log "Drop rows with missing values" df
    check_snps_left(df)
    df
end


function calculate_statistics(df)
    if :NeglogP ∈ df && :P ∉ df
        df = transform(df, :NeglogP => (x -> 10 .^ -x) => :P)
        @info "P calculated from NeglogP."
    end

    #if :P ∈ df && :Z ∉ df
    #    df.Z = [quantile(Normal(), 1 - p / 2) for p in df.P]
    #    @info "Z calculated from P"
    #end

    if :ndiv2 ∈ df && :N ∉ df
        df.N = df.ndiv2 .* 2
        @info "N calculated from ndiv2"
    end

    if :N ∉ df && :neff ∈ df
        rename!(df, :neff => :N)
        @info "neff used for n"
    end

    if :N ∉ df && :nmax ∈ df
        rename!(df, :nmax => :N)
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
            @warn "Impossible p values detected." "SNPs removed" = num_removed SNPs = nrow(df)
            cols_to_select = intersect([:SNP, :RSID, :P], propertynames(large_diff_df))
            selected_df = select(large_diff_df, cols_to_select)
            print_header(selected_df, "Invalid p-values")
        end
    end

    df
end


function check_snps_left(df)
    if nrow(df) < 1000
        error("Fewer than 1000 SNPs remaining")
    end
end


function assign_n(df)
    n = args["assign-n"]
    if n != -1
        if :N in df
            @warn "assign-n selected, but N is already present. Ignored"
            return df
        end
        df.N = fill(n, size(df, 1))
    end

    df
end


function select_cols(df)
    missing_columns = setdiff([:SNP, :N, :A1, :A2], propertynames(df))

    if !isempty(missing_columns)
        missing_columns_str = join(missing_columns, ", ")
            @warn "The following columns are missing from the input file: $missing_columns_str"
    end

    required_cols = [[:Z], [:Direction, :BETA], [:Direction, :OR], [:Direction, :P]]

    for (i, cols) in enumerate(required_cols)
        if all(x -> x in propertynames(df), cols)
            break
        elseif i == length(required_cols)
            msg ="The file does not contain Z or Direction,Stat, or Direction,P"
            if args["no-error"]
                @warn msg
            else
                @error msg
                error(msg)
            end 
        end
    end

    df = select(df, intersect([:SNP, :SNP_ref, :RSID, :CHR, :POS, :A1, :A2, :Stat, :P, :OR, :BETA, :SE, :Z, :Direction, :N], propertynames(df)))
    if args["rsid-as-snp"]
        select!(df, Not(:SNP))
        rename!(df, :RSID => :SNP)
    end
    df
end


function write_output(df)
    args["write-snpid"] && CSV.write(args["out"] * ".snpid", df[:, [:SNP, :RSID, :A1, :A2]], delim = "\t")
    t1 = now()
    CSV.write(args["out"]  * ".sumstats", df, delim = "\t")

    @log "Write formatted sumstat file"
end


function create_logger(log_file)
    file_logger = Logging.SimpleLogger(open(log_file, "w"))
    console_logger = Logging.ConsoleLogger(stderr)
    Logging.global_logger(LoggingExtras.TeeLogger(console_logger, file_logger))
    global start_time = Dates.now()
end



function print_sample(df, cols, title)
    print_header(select(df, intersect(cols, propertynames(df))), title)
end
 

function main()
    global args = parse_commandline()
    create_logger(args["out"] * ".log")

    @info "\n=========================== Program started ===========================\n\n"

    df = read_sumstats(args["in"]) |>
        static(print_header, "Original data") |>
        static(check_data_types) |>
        drop_missing |>
        filter_info |>
        create_SNP |>
        use_ref |>
        split_alleles |>
        filter_alleles |>
        filter_ids |>
        remove_duplicates |>
        identify_effect |>
        calculate_statistics |>
        filter_statistics |>
        assign_n |>
        select_cols |>
        static(write_output) |>
        static(print_header, "After processing")

    @info "Program finished" time=get_time()
    
    data = DataFrame(Type = ["initial SNPs", "final SNPs"], Count = [initial_snps, nrow(df)])
    pretty_table(data, header = ["Type", "Count"], alignment = :r, show_header = false)    
end

main()

