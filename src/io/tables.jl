# Delimited tables (CSV / TSV).
#
# Written tables have one row per resource: a `class` column, an `id` column and one
# column per field (multi-valued cells joined with `|`), so they read back losslessly.
# Any other table can be imported by naming the class with `type` and renaming
# columns with `mapping`.

function read_delim(io::IO, delim::Char)
    rows, row, cell = Vector{String}[], String[], IOBuffer()
    s = read(io, String)
    startswith(s, '﻿') && (s = s[nextind(s, 1):end])     # BOM
    i, inquote = firstindex(s), false
    while i <= lastindex(s)
        c = s[i]
        if inquote
            if c == '"'
                nxt = nextind(s, i)
                nxt <= lastindex(s) && s[nxt] == '"' ? (print(cell, '"'); i = nxt) : (inquote = false)
            else
                print(cell, c)
            end
        elseif c == '"'
            inquote = true
        elseif c == delim
            push!(row, String(take!(cell)))
        elseif c == '\n' || c == '\r'
            c == '\r' && nextind(s, i) <= lastindex(s) && s[nextind(s, i)] == '\n' && (i = nextind(s, i))
            push!(row, String(take!(cell))); push!(rows, row); row = String[]
        else
            print(cell, c)
        end
        i = nextind(s, i)
    end
    (position(cell) > 0 || !isempty(row)) && (push!(row, String(take!(cell))); push!(rows, row))
    return filter!(r -> !(length(r) == 1 && isempty(strip(r[1]))), rows)
end

function cell_string(v, list_sep)
    v === nothing && return ""
    v isa AbstractVector && return join((string(plain_value(x)) for x in v), list_sep)
    return string(plain_value(v))
end

quote_cell(s, delim) = any(c -> c == delim || c == '"' || c == '\n' || c == '\r', s) ?
    "\"" * replace(s, "\"" => "\"\"") * "\"" : s

"""
    write_table(io_or_path, kg; delim = '\\t', types = Resource, list_sep = "|")

One row per resource of type `types` (e.g. `types = EdgeEntity` for an edge list).
Columns: `class`, `id`, then every field used by those types.
"""
function write_table(io::IO, kg::KnowledgeGraph; delim::Char = '\t', types::Type = Resource, list_sep::AbstractString = "|")
    rs = [r for r in kg if r isa types]
    cols = Symbol[]
    for T in unique(typeof.(rs)), f in fieldnames(T)
        f === :id || f in cols || push!(cols, f)
    end
    println(io, join(quote_cell.(["class"; "id"; string.(cols)], delim), delim))
    for r in rs
        cells = [class_name(r); r.id; [hasfield(typeof(r), f) ? cell_string(getfield(r, f), list_sep) : "" for f in cols]]
        println(io, join(quote_cell.(cells, delim), delim))
    end
end

"""
    read_table(io_or_path; delim = '\\t', type = nothing, mapping = Dict(), list_sep = "|", into = nothing)

Import a delimited table. Each row becomes (or updates) one resource.

- the class comes from a `class` column, or from `type` (e.g. `type = Protein`)
- `mapping` renames columns to fields, e.g. `Dict("Gene symbol" => :name, "source" => :left)`;
  other columns are matched to fields by name, and unknown columns are ignored
- rows without an `id` column get generated ids (`protein-1`, ...)
- rows whose id already exists with the same class update that resource
- `into` adds to an existing graph instead of a new one
"""
function read_table(io::IO; delim::Char = '\t', type::Union{Type,Nothing} = nothing,
                    mapping = Dict{String,Symbol}(), list_sep::AbstractString = "|",
                    into::Union{KnowledgeGraph,Nothing} = nothing)
    rows = read_delim(io, delim)
    kg = into === nothing ? KnowledgeGraph() : into
    isempty(rows) && return kg
    header = strip.(rows[1])
    fields = [Symbol(get(mapping, h, h)) for h in header]
    classcol = findfirst(==(:class), fields)
    idcol = findfirst(==(:id), fields)
    classcol === nothing && type === nothing &&
        throw(ArgumentError("table has no `class` column; pass `type = ...`"))
    for (n, row) in enumerate(rows[2:end])
        cell(k) = k <= length(row) ? strip(row[k]) : ""
        T = classcol !== nothing && !isempty(cell(classcol)) ? class_type(cell(classcol)) : type
        T === nothing && throw(ArgumentError("row $(n + 1): no class"))
        id = idcol === nothing || isempty(cell(idcol)) ? nothing : String(cell(idcol))
        if id === nothing
            k = n
            while haskey(kg, "$(lowercase(class_name(T)))-$k"); k += 1; end
            id = "$(lowercase(class_name(T)))-$k"
        end
        existing = get(kg, id, nothing)
        r = existing isa T ? existing : new_resource(T, id)
        for (k, f) in enumerate(fields)
            (f === :id || f === :class || !hasfield(T, f) || isempty(cell(k))) && continue
            v = fieldtype(T, f) <: AbstractVector ? [strip(x) for x in split(cell(k), list_sep) if !isempty(strip(x))] : cell(k)
            try
                assign!(r, f, v)
            catch e
                throw(ArgumentError("row $(n + 1), column \"$(header[k])\": $(sprint(showerror, e))"))
            end
        end
        existing === r || add!(kg, r)
    end
    return kg
end
