# Minimal JSON reader/writer (no dependencies), plus the native JSON and JSON-LD formats.

# --- JSON reader -------------------------------------------------------------

mutable struct JSONReader
    s::String
    i::Int
end

json_error(r::JSONReader, msg) = throw(ArgumentError("JSON parse error at byte $(r.i): $msg"))

function json_ws!(r::JSONReader)
    while r.i <= ncodeunits(r.s) && r.s[r.i] in (' ', '\t', '\n', '\r')
        r.i += 1
    end
end

function json_value!(r::JSONReader)
    json_ws!(r)
    r.i > ncodeunits(r.s) && json_error(r, "unexpected end of input")
    c = r.s[r.i]
    if c == '{'
        r.i += 1
        obj = Dict{String,Any}()
        json_ws!(r)
        r.s[r.i] == '}' && (r.i += 1; return obj)
        while true
            json_ws!(r)
            r.s[r.i] == '"' || json_error(r, "expected object key")
            k = json_string!(r)
            json_ws!(r)
            r.s[r.i] == ':' || json_error(r, "expected ':'")
            r.i += 1
            obj[k] = json_value!(r)
            json_ws!(r)
            r.s[r.i] == ',' ? (r.i += 1) : r.s[r.i] == '}' ? (r.i += 1; return obj) : json_error(r, "expected ',' or '}'")
        end
    elseif c == '['
        r.i += 1
        arr = Any[]
        json_ws!(r)
        r.s[r.i] == ']' && (r.i += 1; return arr)
        while true
            push!(arr, json_value!(r))
            json_ws!(r)
            r.s[r.i] == ',' ? (r.i += 1) : r.s[r.i] == ']' ? (r.i += 1; return arr) : json_error(r, "expected ',' or ']'")
        end
    elseif c == '"'
        return json_string!(r)
    elseif (m = match(r"\G-?\d+(\.\d+)?([eE][+-]?\d+)?", r.s, r.i)) !== nothing
        r.i += ncodeunits(m.match)
        return m[1] === nothing && m[2] === nothing ? parse(Int, m.match) : parse(Float64, m.match)
    end
    for (word, v) in (("true", true), ("false", false), ("null", nothing))
        startswith(SubString(r.s, r.i), word) && (r.i += length(word); return v)
    end
    json_error(r, "unexpected character '$c'")
end

function json_string!(r::JSONReader)
    j = r.i + 1
    while r.s[j] != '"'
        j = nextind(r.s, r.s[j] == '\\' ? nextind(r.s, j) : j)
        j > ncodeunits(r.s) && json_error(r, "unterminated string")
    end
    raw = r.s[r.i+1:prevind(r.s, j)]
    r.i = j + 1
    # JSON escapes are a subset of the N-Triples ones, except '\/' and UTF-16 surrogate pairs.
    raw = replace(raw, "\\/" => "/")
    raw = replace(raw, r"\\u(D[89AB][0-9A-Fa-f]{2})\\u(D[C-F][0-9A-Fa-f]{2})"i => function (m)
        hi, lo = parse(UInt32, m[3:6]; base = 16), parse(UInt32, m[9:12]; base = 16)
        string(Char(0x10000 + ((hi - 0xD800) << 10) + (lo - 0xDC00)))
    end)
    return unescape_literal(raw)
end

"""    parse_json(io_or_string) -> Any"""
parse_json(io::IO) = parse_json(read(io, String))
function parse_json(s::AbstractString)
    r = JSONReader(String(s), 1)
    v = json_value!(r)
    json_ws!(r)
    r.i <= ncodeunits(r.s) && json_error(r, "trailing characters")
    return v
end

# --- JSON writer -------------------------------------------------------------
# Objects are written from Vector{Pair} (order kept) or AbstractDict (keys sorted).

function json_print(io::IO, x, indent::Int = 0)
    pad(n) = "\n" * " "^n
    if x === nothing
        print(io, "null")
    elseif x isa Bool || x isa Integer
        print(io, x)
    elseif x isa AbstractFloat
        print(io, isfinite(x) ? x : "null")
    elseif x isa AbstractString || x isa Symbol
        print(io, '"', escape_json(string(x)), '"')
    elseif x isa AbstractDict || (x isa AbstractVector && eltype(x) <: Pair)
        pairs_ = x isa AbstractDict ? sort!(collect(x); by = first) : x
        isempty(pairs_) && return print(io, "{}")
        print(io, "{")
        for (k, (key, v)) in enumerate(pairs_)
            print(io, k == 1 ? "" : ",", pad(indent + 2), '"', escape_json(string(key)), "\": ")
            json_print(io, v, indent + 2)
        end
        print(io, pad(indent), "}")
    elseif x isa AbstractVector || x isa Tuple
        isempty(x) && return print(io, "[]")
        if all(v -> v isa Union{AbstractString,Number,Symbol,Nothing}, x)
            print(io, "[")
            join(io, (sprint(json_print, v) for v in x), ", ")
            print(io, "]")
        else
            print(io, "[")
            for (k, v) in enumerate(x)
                print(io, k == 1 ? "" : ",", pad(indent + 2))
                json_print(io, v, indent + 2)
            end
            print(io, pad(indent), "]")
        end
    else
        json_print(io, string(x), indent)
    end
end

function escape_json(s::AbstractString)
    io = IOBuffer()
    for c in s
        c == '"' ? print(io, "\\\"") : c == '\\' ? print(io, "\\\\") :
        c == '\n' ? print(io, "\\n") : c == '\r' ? print(io, "\\r") : c == '\t' ? print(io, "\\t") :
        c < ' ' ? print(io, "\\u", string(UInt32(c); base = 16, pad = 4)) : print(io, c)
    end
    return String(take!(io))
end

# --- converting plain values into fields (shared by JSON and tables) ----------

convert_value(::Type{String}, x) = string(x)
convert_value(::Type{Symbol}, x) = Symbol(x)
convert_value(::Type{Float64}, x) = x isa Number ? Float64(x) : parse(Float64, strip(string(x)))
convert_value(::Type{Int}, x) = x isa Integer ? Int(x) : x isa Number ? Int(x) : parse(Int, strip(string(x)))
convert_value(::Type{DateTime}, x) = x isa DateTime ? x : DateTime(first(strip(string(x)), 19))

"""
    assign!(r, field, value)

Set `field` of resource `r` from a plain value (string, number, vector), converting to
the field's type. Vector fields accept a vector or a single value; `nothing` clears.
"""
function assign!(r::Resource, f::Symbol, x)
    FT = fieldtype(typeof(r), f)
    if FT <: AbstractVector
        vals = x === nothing ? [] : x isa AbstractVector ? x : [x]
        setfield!(r, f, eltype(FT)[convert_value(eltype(FT), v) for v in vals])
    else
        T = Base.nonnothingtype(FT)
        setfield!(r, f, x === nothing ? (Nothing <: FT ? nothing : throw(ArgumentError("$f is required"))) :
                        convert_value(T, x))
    end
    return r
end

new_resource(::Type{T}, id; type = :Unknown) where {T<:Resource} =
    T <: EdgeEntity ? T(; id = String(id), type = Symbol(type)) : T(; id = String(id))

function class_type(name::AbstractString)
    T = get(CLASS_TYPES, name, nothing)
    T !== nothing && haskey(SCHEMA, T) || throw(ArgumentError("unknown BioNExt class \"$name\""))
    return T
end

plain_value(v::DateTime) = Dates.format(v, dateformat"yyyy-mm-ddTHH:MM:SS")
plain_value(v) = v

# --- native JSON ---------------------------------------------------------------

"""
    write_json(io_or_path, kg)

Plain JSON document: `{"format": "bionext-json", "base", "prefixes", "resources": [...]}`,
one object per resource with `"class"`, `"id"` and every non-empty field by its Julia name.
"""
function write_json(io::IO, kg::KnowledgeGraph)
    rs = map(collect(kg)) do r
        obj = Pair{String,Any}["class" => class_name(r), "id" => r.id]
        for f in fieldnames(typeof(r))
            f === :id && continue
            v = getfield(r, f)
            (v === nothing || (v isa AbstractVector && isempty(v))) && continue
            push!(obj, string(f) => v isa AbstractVector ? map(plain_value, v) : plain_value(v))
        end
        obj
    end
    json_print(io, Pair{String,Any}["format" => "bionext-json", "version" => 1, "base" => kg.base,
                                    "prefixes" => kg.prefixes, "resources" => rs])
    println(io)
end

"""    read_json(io_or_path) -> KnowledgeGraph"""
function read_json(src)
    doc = src isa IO ? parse_json(src) : parse_json(read(src, String))
    doc isa AbstractDict && haskey(doc, "resources") || throw(ArgumentError("not a BioNExt JSON document (no \"resources\")"))
    kg = KnowledgeGraph(; base = get(doc, "base", BNR), prefixes = Dict{String,String}(get(doc, "prefixes", Dict())))
    for obj in doc["resources"]
        T = class_type(obj["class"])
        r = new_resource(T, obj["id"]; type = get(obj, "type", :Unknown))
        for (k, v) in obj
            k in ("class", "id") && continue
            f = Symbol(k)
            hasfield(T, f) ? assign!(r, f, v) : @warn "ignoring unknown field \"$k\" on $(obj["id"])"
        end
        add!(kg, r)
    end
    return kg
end

# --- JSON-LD -------------------------------------------------------------------

function jsonld_object(x::Literal, ctx_iri)
    x.lang !== nothing && return Pair{String,Any}["@value" => x.value, "@language" => x.lang]
    x.datatype === nothing || x.datatype == XSD * "string" ? x.value :
        Pair{String,Any}["@value" => x.value, "@type" => ctx_iri(x.datatype)]
end
jsonld_object(x::IRI, ctx_iri) = Pair{String,Any}["@id" => ctx_iri(x.value)]

"""
    write_jsonld(io_or_path, kg; inverses = true)

Compact JSON-LD with an `@context` of the graph's prefixes and an `@graph` array.
"""
function write_jsonld(io::IO, kg::KnowledgeGraph; inverses::Bool = true)
    pfx = prefixes(kg)
    ctx_iri(iri) = (t = turtle_iri(iri, pfx); startswith(t, '<') ? iri : t)
    order, nodes_ = String[], Dict{String,Vector{Pair{String,Any}}}()
    for t in unique(triples(kg; inverses))
        s = ctx_iri(t.subject.value)
        haskey(nodes_, s) || (push!(order, s); nodes_[s] = Pair{String,Any}["@id" => s])
        node = nodes_[s]
        key = t.predicate.value == RDF_TYPE ? "@type" : ctx_iri(t.predicate.value)
        val = key == "@type" ? ctx_iri(t.object.value) : jsonld_object(t.object, ctx_iri)
        i = findfirst(p -> p.first == key, node)
        i === nothing ? push!(node, key => Any[val]) : push!(node[i].second, val)
    end
    graph = [Pair{String,Any}[k => (v isa Vector && length(v) == 1 ? v[1] : v) for (k, v) in nodes_[s]] for s in order]
    json_print(io, Pair{String,Any}["@context" => pfx, "@graph" => graph])
    println(io)
end

"""
    parse_jsonld(io_or_path) -> (triples, prefixes)

Read JSON-LD in compacted or expanded form: `@context` with prefix/term
definitions (string or `{"@id", "@type": "@id"}`), `@graph`, `@id`, `@type`,
value objects, node references and nested node objects. Remote contexts are not fetched.
"""
function parse_jsonld(src)
    doc = src isa IO ? parse_json(src) : parse_json(read(src, String))
    ctx = Dict{String,Any}()
    collect_ctx!(c) = c isa AbstractDict ? merge!(ctx, c) : c isa AbstractVector ? foreach(collect_ctx!, c) : nothing
    doc isa AbstractDict && collect_ctx!(get(doc, "@context", nothing))
    vocab = get(ctx, "@vocab", nothing)
    ts, bn = Triple[], Ref(0)

    function iri(s::AbstractString; term = false)
        if term && haskey(ctx, s)
            d = ctx[s]
            return iri(d isa AbstractDict ? d["@id"] : d)
        end
        m = match(r"^([^:/?#]*):(.*)$", s)
        if m !== nothing && !startswith(m[2], "//") && haskey(ctx, m[1])
            d = ctx[m[1]]
            return (d isa AbstractDict ? d["@id"] : d) * m[2]
        end
        m !== nothing && return String(s)
        term && vocab !== nothing && return vocab * s
        return String(s)
    end
    coerced(key) = (d = get(ctx, key, nothing); d isa AbstractDict ? get(d, "@type", nothing) : nothing)

    function value(v, key)
        if v isa AbstractDict
            haskey(v, "@value") || return node!(v)
            lang = get(v, "@language", nothing)
            dt = haskey(v, "@type") ? iri(v["@type"]; term = true) : nothing
            return Literal(string(v["@value"]), dt, lang)
        end
        ct = coerced(key)
        ct == "@id" && v isa AbstractString && return IRI(iri(v))
        v isa Bool && return Literal(string(v), XSD_BOOLEAN, nothing)
        v isa Integer && return Literal(string(v), XSD_INTEGER, nothing)
        v isa AbstractFloat && return Literal(string(v), XSD_DOUBLE, nothing)
        return Literal(string(v), ct === nothing ? nothing : iri(ct; term = true), nothing)
    end

    function node!(obj)
        s = IRI(haskey(obj, "@id") ? iri(obj["@id"]) : "_:j" * string(bn[] += 1))
        for (k, v) in obj
            k in ("@id", "@context") && continue
            vals = v isa AbstractVector ? v : [v]
            if k == "@type"
                foreach(t -> push!(ts, Triple(s, IRI(RDF_TYPE), IRI(iri(t; term = true)))), vals)
            elseif k == "@graph"
                foreach(node!, vals)
            elseif !startswith(k, "@")
                p = IRI(iri(k; term = true))
                foreach(x -> push!(ts, Triple(s, p, value(x, k))), vals)
            end
        end
        return s
    end

    tops = doc isa AbstractVector ? doc : haskey(doc, "@graph") ? doc["@graph"] : [doc]
    foreach(node!, tops isa AbstractVector ? tops : [tops])
    pfx = Dict{String,String}(k => v for (k, v) in ctx if v isa AbstractString && !startswith(k, "@") && occursin(r"[/#:]$", v))
    return ts, pfx
end

"""    read_jsonld(io_or_path; base = nothing, prefixes = Dict()) -> KnowledgeGraph"""
function read_jsonld(src; base = nothing, prefixes = Dict{String,String}())
    ts, pfx = parse_jsonld(src)
    return from_triples(ts; graph_options(pfx; base, prefixes)...)
end
