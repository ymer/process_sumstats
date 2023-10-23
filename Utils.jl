using PrettyTables, Printf, Dates

import Base: ∉
∉(c::Symbol, df::DataFrame) = !(c in propertynames(df))
∉(c::String, df::DataFrame) = !(c in names(df))

import Base: ∈
∈(c::Symbol, df::DataFrame) = c in propertynames(df)
∈(c::String, df::DataFrame) = c in names(df)


struct Percent
    value::Float64
end

Base.show(io::IO, p::Percent) = print(io, @sprintf("%.1f%%", p.value * 100))


function static(f, arg=nothing)
    return x -> begin
        if arg === nothing
            f(x)
        else
            f(x, arg)
        end
        x
    end
end


function print_header(df, title = "")
    pretty_table(first(df, 4), header_crayon = crayon"yellow bold", title = title, show_omitted_cell_summary = true)
    col_widths = [length(string(col)) for col in names(df)]
    total_width = sum(col_widths) + length(names(df)) * 3
    term_width = displaysize(stdout)[2]

    if total_width > term_width
        println("The table does not fit the screen. Column names are:")
        println(names(df))
    end
end


function print_header2(df, title = "")
    pretty_table(
        first(df, 4),
        header_crayon = crayon"yellow bold",
        title = title,
        #display_size = (13, 300),
        show_omitted_cell_summary = true
        )
end


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


function get_time(t1)::Duration

    total_seconds = round(Int, (now() - t1).value / 1000)

    minutes = total_seconds ÷ 60
    seconds = total_seconds % 60

    return Duration("$(minutes):$(lpad(seconds, 2, "0"))")
end


macro log(msg, df)
    return quote
        formatted_snps = format(nrow($(esc(df))), commas=true)
        @info $(esc(msg)) * "\n SNPs: " * formatted_snps
    end
end

macro log(msg, df, t)
    return quote
        formatted_time = get_time($(esc(t)))
        formatted_snps = format(nrow($(esc(df))), commas=true)
        @info $(esc(msg)) * "\n SNPs: " * formatted_snps * "\n time: " * formatted_time
    end
end

macro log(msg, fn, df, t)
    return quote
        local formatted_time = get_time($(esc(t)))
        local formatted_snps = format(nrow($(esc(df))), commas=true)
        local log_msg = $(esc(msg)) * "\n file: " * $(esc(fn)) * "\n SNPs in file: " * formatted_snps * "\n time: " * string(formatted_time)
        @info log_msg
    end
end


