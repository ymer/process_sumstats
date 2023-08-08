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
    pretty_table(
        first(df, 3),
        header_crayon = crayon"yellow bold",
        title = title,
        #display_size = (13, 300),
        show_omitted_cell_summary = false
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


