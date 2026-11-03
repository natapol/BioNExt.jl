# RDF serialization (Turtle, N-Triples) and N-Triples import.

const RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
const RDFS = "http://www.w3.org/2000/01/rdf-schema#"
const XSD = "http://www.w3.org/2001/XMLSchema#"
const IDENTIFIERS = "https://identifiers.org/"
const RDF_TYPE = RDF * "type"

"""An RDF IRI (or a blank node written `_:label`)."""
struct IRI
    value::String
end

"""An RDF literal with an optional datatype IRI or language tag."""
struct Literal
    value::String
    datatype::Union{String,Nothing}
    lang::Union{String,Nothing}
end
Literal(value::AbstractString) = Literal(String(value), nothing, nothing)

"""An RDF triple."""
struct Triple
    subject::IRI
    predicate::IRI
    object::Union{IRI,Literal}
end

prefixes(kg::KnowledgeGraph) =
    merge(Dict("bn" => BN, "bnr" => kg.base, "rdf" => RDF, "rdfs" => RDFS, "xsd" => XSD), kg.prefixes)

const CURIE = r"^([A-Za-z][\w.\-]*):(.+)$"
isabsolute(s) = startswith(s, "http://") || startswith(s, "https://") || startswith(s, "urn:")

"""
    expand(kg, id) -> String

IRI for an identifier: full IRIs pass through, CURIEs with a known prefix are
expanded (unknown prefixes go to identifiers.org), bare names are minted under `kg.base`.
"""
function expand(kg::KnowledgeGraph, id::AbstractString)
    isabsolute(id) && return String(id)
    m = match(CURIE, id)
    m === nothing && return kg.base * id
    ns = get(prefixes(kg), m[1], nothing)
    return ns === nothing ? IDENTIFIERS * id : ns * m[2]
end

"""    compact(kg, iri) -> String

Inverse of `expand`."""
function compact(kg::KnowledgeGraph, iri::AbstractString)
    if startswith(iri, kg.base)
        local_ = iri[length(kg.base)+1:end]
        return occursin(':', local_) ? String(iri) : local_
    end
    startswith(iri, IDENTIFIERS) && return iri[length(IDENTIFIERS)+1:end]
    for (p, ns) in prefixes(kg)
        p == "bnr" && continue
        startswith(iri, ns) && return p * ":" * iri[length(ns)+1:end]
    end
    return String(iri)
end

to_literal(v::AbstractString, _) = Literal(v)
to_literal(v::AbstractFloat, range) = Literal(string(v), XSD * range, nothing)
to_literal(v::Integer, range) = Literal(string(v), XSD * range, nothing)
to_literal(v::DateTime, range) = Literal(Dates.format(v, dateformat"yyyy-mm-ddTHH:MM:SS"), XSD * range, nothing)

const INVERSE_OF = Dict(v => k for (k, v) in INVERSES)

"""
    triples(kg; inverses = true) -> Vector{Triple}

The graph as RDF triples. With `inverses = true`, inverse properties that are
stored only one way (e.g. `bn:hasVariant`, `bn:participatesInPathway`,
`bn:participatesInEdge`) are materialized as well.
"""
function triples(kg::KnowledgeGraph; inverses::Bool = true)
    out = Triple[]
    bn(p) = IRI(BN * p)
    for r in kg
        s = IRI(expand(kg, r.id))
        push!(out, Triple(s, IRI(RDF_TYPE), bn(class_name(r))))
        for p in schema(r), v in field_values(r, p.field)
            pred = p.pred
            if r isa OrthologGroup && p.field === :members && get(kg, v, nothing) isa Protein
                pred = "hasProteinMember"
            end
            o = p.kind === :literal ? to_literal(v, p.range) :
                p.kind === :term ? bn(string(v)) : IRI(expand(kg, v))
            push!(out, Triple(s, bn(pred), o))
            inverses && p.kind === :ref || continue
            haskey(INVERSE_OF, pred) && push!(out, Triple(o, bn(INVERSE_OF[pred]), s))
            if r isa EdgeEntity && p.field in (:left, :right) && get(kg, v, nothing) isa NodeEntity
                push!(out, Triple(o, bn("participatesInEdge"), s))
            end
        end
    end
    return out
end

# --- writing ---------------------------------------------------------------

function escape_string_literal(s::AbstractString)
    io = IOBuffer()
    for c in s
        c == '\\' ? print(io, "\\\\") : c == '"' ? print(io, "\\\"") :
        c == '\n' ? print(io, "\\n") : c == '\r' ? print(io, "\\r") :
        c == '\t' ? print(io, "\\t") : print(io, c)
    end
    return String(take!(io))
end

const PN_LOCAL = r"^[A-Za-z0-9_]([A-Za-z0-9_.\-]*[A-Za-z0-9_\-])?$"

function turtle_iri(iri::String, pfx)
    best = ""
    for (p, ns) in pfx
        if startswith(iri, ns) && occursin(PN_LOCAL, iri[length(ns)+1:end]) &&
           (isempty(best) || length(ns) > length(pfx[best]))
            best = p
        end
    end
    return isempty(best) ? "<$iri>" : best * ":" * iri[length(pfx[best])+1:end]
end

function format_term(x::Literal, fmt_iri)
    s = "\"" * escape_string_literal(x.value) * "\""
    x.lang !== nothing && return s * "@" * x.lang
    x.datatype === nothing || x.datatype == XSD * "string" ? s : s * "^^" * fmt_iri(x.datatype)
end
format_term(x::IRI, fmt_iri) = fmt_iri(x.value)

"""
    write_turtle(io_or_path, kg; inverses = true)

Serialize the graph as Turtle, grouped by subject.
"""
function write_turtle(io::IO, kg::KnowledgeGraph; inverses::Bool = true)
    pfx = prefixes(kg)
    fmt(iri) = turtle_iri(iri, pfx)
    for p in sort!(collect(keys(pfx)))
        println(io, "@prefix ", p, ": <", pfx[p], "> .")
    end
    order, groups = IRI[], Dict{IRI,Vector{Triple}}()
    for t in triples(kg; inverses)
        haskey(groups, t.subject) || (push!(order, t.subject); groups[t.subject] = Triple[])
        t in groups[t.subject] || push!(groups[t.subject], t)
    end
    for s in order
        preds, objs = IRI[], Dict{IRI,Vector{String}}()
        for t in groups[s]
            haskey(objs, t.predicate) || (push!(preds, t.predicate); objs[t.predicate] = String[])
            push!(objs[t.predicate], format_term(t.object, fmt))
        end
        lines = [(p.value == RDF_TYPE ? "a" : fmt(p.value)) * " " * join(objs[p], ", ") for p in preds]
        print(io, "\n", fmt(s.value), " ", join(lines, " ;\n    "), " .\n")
    end
end

"""
    write_ntriples(io_or_path, kg; inverses = true)

Serialize the graph as N-Triples (one triple per line, full IRIs).
"""
function write_ntriples(io::IO, kg::KnowledgeGraph; inverses::Bool = true)
    fmt(iri) = "<$iri>"
    for t in unique(triples(kg; inverses))
        println(io, fmt(t.subject.value), " ", fmt(t.predicate.value), " ", format_term(t.object, fmt), " .")
    end
end

for f in (:write_turtle, :write_ntriples)
    @eval $f(path::AbstractString, kg::KnowledgeGraph; kw...) = open(io -> $f(io, kg; kw...), path, "w")
end

# --- reading ---------------------------------------------------------------

function unescape_literal(s::AbstractString)
    io, i = IOBuffer(), firstindex(s)
    while i <= lastindex(s)
        c = s[i]
        if c == '\\' && i < lastindex(s)
            i = nextind(s, i); e = s[i]
            if e == 'u' || e == 'U'
                n = e == 'u' ? 4 : 8
                print(io, Char(parse(UInt32, s[nextind(s, i):nextind(s, i, n)]; base = 16)))
                i = nextind(s, i, n)
            else
                print(io, get(Dict('t' => '\t', 'n' => '\n', 'r' => '\r', 'b' => '\b', 'f' => '\f'), e, e))
            end
        else
            print(io, c)
        end
        i = nextind(s, i)
    end
    return String(take!(io))
end

function parse_nt_term(s::AbstractString, i::Int)
    while isspace(s[i]); i = nextind(s, i); end
    c = s[i]
    if c == '<'
        j = findnext('>', s, i)
        return IRI(s[nextind(s, i):prevind(s, j)]), nextind(s, j)
    elseif c == '_'
        j = something(findnext(isspace, s, i), lastindex(s) + 1)
        return IRI(s[i:prevind(s, j)]), j
    elseif c == '"'
        j = nextind(s, i)
        while s[j] != '"'
            j = nextind(s, s[j] == '\\' ? nextind(s, j) : j)
        end
        value = unescape_literal(s[nextind(s, i):prevind(s, j)])
        k = nextind(s, j)
        if k <= lastindex(s) && s[k] == '@'
            m = match(r"^@([A-Za-z0-9\-]+)", SubString(s, k))
            return Literal(value, nothing, m[1]), k + ncodeunits(m.match)
        elseif startswith(SubString(s, k), "^^<")
            e = findnext('>', s, k)
            return Literal(value, s[k+3:prevind(s, e)], nothing), nextind(s, e)
        end
        return Literal(value), k
    end
    error("unexpected character '$c' at position $i")
end

"""    parse_ntriples(io_or_path) -> Vector{Triple}"""
function parse_ntriples(io::IO)
    out = Triple[]
    for (n, line) in enumerate(eachline(io))
        line = strip(line)
        (isempty(line) || startswith(line, '#')) && continue
        try
            s, i = parse_nt_term(line, 1)
            p, i = parse_nt_term(line, i)
            o, _ = parse_nt_term(line, i)
            push!(out, Triple(s, p, o))
        catch e
            throw(ArgumentError("N-Triples parse error on line $n: $(sprint(showerror, e))"))
        end
    end
    return out
end
parse_ntriples(path::AbstractString) = open(parse_ntriples, path)

parse_value(::Type{String}, s) = s
parse_value(::Type{Float64}, s) = parse(Float64, s)
parse_value(::Type{Int}, s) = parse(Int, s)
parse_value(::Type{DateTime}, s) = DateTime(first(s, 19))

function set_value!(r::Resource, f::Symbol, v)
    T = fieldtype(typeof(r), f)
    if T <: AbstractVector
        vec = getfield(r, f)
        v in vec || push!(vec, v)
    else
        setfield!(r, f, v)
    end
end

# Earlier namespaces of the ontology; data using them is read as the current ones.
const LEGACY_NAMESPACES = [
    "https://jamiesilakom618.github.io/BioNExt/ontology/" => BN,
    "https://jamiesilakom618.github.io/BioNExt/resource/" => BNR,
]

function upgrade_iri(s::AbstractString)
    for (old, new) in LEGACY_NAMESPACES
        startswith(s, old) && return new * s[length(old)+1:end]
    end
    return String(s)
end
upgrade(x::IRI) = IRI(upgrade_iri(x.value))
upgrade(x::Literal) = x.datatype === nothing ? x : Literal(x.value, upgrade_iri(x.datatype), x.lang)
upgrade(t::Triple) = Triple(upgrade(t.subject), upgrade(t.predicate), upgrade(t.object))

"""
    from_triples(triples; base = BNR, prefixes = Dict()) -> KnowledgeGraph

Build a graph from RDF triples. Subjects typed with a concrete BioNExt class
become resources; `bn:` predicates are mapped back through the schema (inverse
and sub-properties are folded into their canonical field). Other triples are ignored.
IRIs in a legacy BioNExt namespace are upgraded to the current one.
"""
function from_triples(ts::AbstractVector{Triple}; base::AbstractString = BNR, prefixes = Dict{String,String}())
    ts = map(upgrade, ts)
    kg = KnowledgeGraph(; base = upgrade_iri(base), prefixes = Dict(k => upgrade_iri(v) for (k, v) in prefixes))
    for t in ts
        t.predicate.value == RDF_TYPE && t.object isa IRI && startswith(t.object.value, BN) || continue
        T = get(CLASS_TYPES, t.object.value[length(BN)+1:end], nothing)
        (T === nothing || !haskey(SCHEMA, T)) && continue
        id = compact(kg, t.subject.value)
        haskey(kg, id) && continue
        add!(kg, T <: EdgeEntity ? T(; id, type = :Unknown) : T(; id))
    end
    for t in ts
        pred = t.predicate.value
        startswith(pred, BN) || continue
        pred = pred[length(BN)+1:end]
        subj, obj = t.subject, t.object
        if haskey(INVERSES, pred)
            obj isa IRI || continue
            pred, subj, obj = INVERSES[pred], obj, subj
        end
        pred = get(SUBPROPERTIES, pred, pred)
        r = get(kg, compact(kg, subj.value), nothing)
        r === nothing && continue
        i = findfirst(p -> p.pred == pred, schema(r))
        i === nothing && continue
        p = schema(r)[i]
        if p.kind === :literal
            obj isa Literal || continue
            FT = fieldtype(typeof(r), p.field)
            v = parse_value(FT <: AbstractVector ? eltype(FT) : Base.nonnothingtype(FT), obj.value)
        elseif p.kind === :term
            obj isa IRI && startswith(obj.value, BN) || continue
            v = Symbol(obj.value[length(BN)+1:end])
        else
            obj isa IRI || continue
            v = compact(kg, obj.value)
        end
        set_value!(r, p.field, v)
    end
    return kg
end

"""    read_ntriples(io_or_path; base = BNR, prefixes = Dict()) -> KnowledgeGraph"""
read_ntriples(src; kw...) = from_triples(parse_ntriples(src); kw...)
