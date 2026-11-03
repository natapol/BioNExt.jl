# BN and BNR (the ontology and default resource namespaces) come from vocabulary_data.jl.

"""
    KnowledgeGraph(; base = BNR, prefixes = Dict())

A collection of BioNExt resources keyed by identifier, kept in insertion order.

Identifiers are plain strings: a bare local name (`"TP53"`, minted under `base`),
a CURIE (`"UniProt:P04637"`, expanded with `prefixes` or identifiers.org), or a full IRI.

Queries scan the graph unless it is indexed; see `index!`.
"""
mutable struct KnowledgeGraph
    resources::Dict{String,Resource}
    order::Vector{String}
    base::String
    prefixes::Dict{String,String}
    seq::Dict{String,Int}          # insertion sequence number of each id, for ordering
    counter::Int
    # Reverse reference index (target id => [referring id => field]); `nothing` when off.
    index::Union{Nothing,Dict{String,Vector{Pair{String,Symbol}}}}
end

function KnowledgeGraph(; base::AbstractString = BNR, prefixes = Dict{String,String}())
    KnowledgeGraph(Dict{String,Resource}(), String[], String(base), Dict{String,String}(prefixes),
                   Dict{String,Int}(), 0, nothing)
end

Base.length(kg::KnowledgeGraph) = length(kg.order)
Base.haskey(kg::KnowledgeGraph, id::AbstractString) = haskey(kg.resources, id)
Base.getindex(kg::KnowledgeGraph, id::AbstractString) = kg.resources[id]
Base.get(kg::KnowledgeGraph, id::AbstractString, default) = get(kg.resources, id, default)
Base.iterate(kg::KnowledgeGraph, i = 1) =
    i > length(kg.order) ? nothing : (kg.resources[kg.order[i]], i + 1)
Base.eltype(::Type{KnowledgeGraph}) = Resource
Base.keys(kg::KnowledgeGraph) = kg.order

"""
    add!(kg, resources...; replace = false) -> kg

Add resources. Adding an id that already exists throws unless `replace = true`.
"""
function add!(kg::KnowledgeGraph, rs::Resource...; replace::Bool = false)
    for r in rs
        isempty(r.id) && throw(ArgumentError("resource of type $(typeof(r)) has an empty id"))
        if haskey(kg.resources, r.id)
            replace || throw(ArgumentError("duplicate id \"$(r.id)\" (use replace = true to overwrite)"))
            index_remove!(kg, kg.resources[r.id])
        else
            push!(kg.order, r.id)
            kg.seq[r.id] = (kg.counter += 1)
        end
        kg.resources[r.id] = r
        index_add!(kg, r)
    end
    return kg
end
Base.push!(kg::KnowledgeGraph, rs::Resource...) = add!(kg, rs...)

"""    delete!(kg, id) -> kg

Remove a resource. References to it from other resources are left in place (see `validate`)."""
function Base.delete!(kg::KnowledgeGraph, id::AbstractString)
    haskey(kg.resources, id) && delete_many!(kg, Set([String(id)]))
    return kg
end

# Remove several resources with a single pass over the insertion order.
function delete_many!(kg::KnowledgeGraph, ids::AbstractSet{String})
    for id in ids
        r = get(kg.resources, id, nothing)
        r === nothing && continue
        index_remove!(kg, r)
        delete!(kg.resources, id)
        delete!(kg.seq, id)
    end
    filter!(id -> !(id in ids), kg.order)
    return kg
end

# --- reverse reference index ---------------------------------------------------

"""
    index!(kg) -> kg

Build a reverse index of all references and keep it up to date. Afterwards, `edges_of`,
`neighbors`, `pathways_of`, `complexes_of`, `variants`, `ortholog_groups_of` and
`references` look up the index instead of scanning the whole graph.

`add!`, `delete!`, `merge!`, `load!`, `resolve!` and `integrate!` keep the index current.
Editing reference fields directly (e.g. `push!(kg["e1"].right, "x")`) does **not**:
call `reindex!(kg)` after such edits. Results are verified against the current fields,
so a stale index can miss new references but never reports removed ones.
"""
function index!(kg::KnowledgeGraph)
    kg.index = Dict{String,Vector{Pair{String,Symbol}}}()
    foreach(r -> index_add!(kg, r), kg)
    return kg
end

"""    reindex!(kg) -> kg

Rebuild the index after editing reference fields directly. Does nothing if `kg` is not indexed."""
reindex!(kg::KnowledgeGraph) = kg.index === nothing ? kg : index!(kg)

"""    unindex!(kg) -> kg

Drop the index; queries scan again."""
unindex!(kg::KnowledgeGraph) = (kg.index = nothing; kg)

"""    isindexed(kg) -> Bool"""
isindexed(kg::KnowledgeGraph) = kg.index !== nothing

function index_add!(kg::KnowledgeGraph, r::Resource)
    kg.index === nothing && return
    for p in schema(r)
        p.kind === :ref || continue
        for v in field_values(r, p.field)
            push!(get!(kg.index, v, Pair{String,Symbol}[]), r.id => p.field)
        end
    end
end

function index_remove!(kg::KnowledgeGraph, r::Resource)
    kg.index === nothing && return
    for p in schema(r)
        p.kind === :ref || continue
        for v in field_values(r, p.field)
            lst = get(kg.index, v, nothing)
            lst === nothing || filter!(e -> e.first != r.id, lst)
        end
    end
end

"""    referrers(kg, x) -> Vector{Pair{String,Symbol}}

Every `(resource id => field)` whose reference field contains `x`, in insertion order.
Uses the index when there is one."""
function referrers(kg::KnowledgeGraph, x)
    id = _id(x)
    if kg.index === nothing
        out = Pair{String,Symbol}[]
        for r in kg, p in schema(r)
            p.kind === :ref && id in field_values(r, p.field) && push!(out, r.id => p.field)
        end
        return out
    end
    lst = get(kg.index, id, nothing)
    lst === nothing && return Pair{String,Symbol}[]
    out = unique!([e for e in lst if haskey(kg.resources, e.first) && id in field_values(kg[e.first], e.second)])
    return sort!(out; by = e -> kg.seq[e.first], alg = Base.Sort.DEFAULT_STABLE)
end

"""    merge!(kg, other; replace = false) -> kg"""
function Base.merge!(kg::KnowledgeGraph, other::KnowledgeGraph; replace::Bool = false)
    kg === other && return kg
    merge!(kg.prefixes, other.prefixes)
    for r in other
        (replace || !haskey(kg, r.id)) && add!(kg, r; replace = true)
    end
    return kg
end

"""    resources(kg, T = Resource) -> Vector{T}

All resources of type `T` (abstract types allowed), e.g. `resources(kg, NodeEntity)`."""
resources(kg::KnowledgeGraph, ::Type{T} = Resource) where {T<:Resource} =
    T[r for r in kg if r isa T]

"""    nodes(kg) -> Vector{NodeEntity}

All node entities, in insertion order."""
nodes(kg::KnowledgeGraph) = resources(kg, NodeEntity)
"""    edges(kg) -> Vector{EdgeEntity}

All edge entities (interactions, influences, reactions), in insertion order."""
edges(kg::KnowledgeGraph) = resources(kg, EdgeEntity)
"""    pathways(kg) -> Vector{PathwayModule}"""
pathways(kg::KnowledgeGraph) = resources(kg, PathwayModule)

_id(x::AbstractString) = String(x)
_id(r::Resource) = r.id

"""
    edges_of(kg, x; side = :any, type = nothing) -> Vector{EdgeEntity}

Edges in which `x` (id or resource) participates — the `bn:participatesInEdge` view.
`side` is `:left`, `:right` or `:any`; `type` optionally restricts to edges whose
relationship type `is_a` the given term.
"""
function edges_of(kg::KnowledgeGraph, x; side::Symbol = :any, type::Union{Symbol,Nothing} = nothing)
    if !isindexed(kg)       # direct scan of the edges: much cheaper than a generic `referrers` scan
        id = _id(x)
        return filter(edges(kg)) do e
            (type === nothing || is_a(e.type, type)) &&
                ((side !== :right && id in e.left) || (side !== :left && id in e.right))
        end
    end
    fields = side === :left ? (:left,) : side === :right ? (:right,) : (:left, :right)
    out, seen = EdgeEntity[], Set{String}()
    for (s, f) in referrers(kg, x)
        f in fields && !(s in seen) || continue
        e = kg[s]
        e isa EdgeEntity && (type === nothing || is_a(e.type, type)) || continue
        push!(seen, s); push!(out, e)
    end
    return out
end

function _referring(kg, x, field::Symbol, ::Type{T}) where {T}
    isindexed(kg) || return (id = _id(x); filter(r -> hasfield(typeof(r), field) && id in getfield(r, field), resources(kg, T)))
    return unique!(T[kg[s] for (s, f) in referrers(kg, x) if f === field && kg[s] isa T])
end

"""
    neighbors(kg, x; direction = :both, type = nothing) -> Vector{String}

Ids connected to `x` through edges. `:out` follows left→right, `:in` right→left.
"""
function neighbors(kg::KnowledgeGraph, x; direction::Symbol = :both, type::Union{Symbol,Nothing} = nothing)
    id = _id(x)
    out = String[]
    for e in edges_of(kg, id; type)
        direction !== :in && id in e.left && append!(out, e.right)
        direction !== :out && id in e.right && append!(out, e.left)
    end
    return unique!(filter!(!=(id), out))
end

"""    pathways_of(kg, x) -> Vector{PathwayModule}

Pathways listing `x` as a component (`bn:participatesInPathway`)."""
pathways_of(kg::KnowledgeGraph, x) = _referring(kg, x, :components, PathwayModule)

"""    components(kg, pathway; recursive = false) -> Vector{String}

Component ids of a pathway; with `recursive = true`, sub-pathways are expanded too."""
function components(kg::KnowledgeGraph, x; recursive::Bool = false)
    p = kg[_id(x)]::PathwayModule
    recursive || return copy(p.components)
    out, seen, stack = String[], Set{String}(), copy(p.components)
    while !isempty(stack)
        c = popfirst!(stack)
        c in seen && continue
        push!(seen, c); push!(out, c)
        sub = get(kg, c, nothing)
        sub isa PathwayModule && append!(stack, sub.components)
    end
    return out
end

"""    complexes_of(kg, x) -> Vector{MolecularComplex}

Complexes having `x` as a subunit (`bn:memberOfComplex`)."""
complexes_of(kg::KnowledgeGraph, x) = _referring(kg, x, :subunits, MolecularComplex)

"""    variants(kg, x) -> Vector{NodeEntity}

Resources declared `isVariantOf` `x` (`bn:hasVariant`)."""
variants(kg::KnowledgeGraph, x) = _referring(kg, x, :variant_of, NodeEntity)

"""    ortholog_groups_of(kg, x) -> Vector{OrthologGroup}"""
ortholog_groups_of(kg::KnowledgeGraph, x) = _referring(kg, x, :members, OrthologGroup)

"""    references(kg, x) -> Vector{Pair{String,Symbol}}

Every `(resource id => field)` that refers to `x`; useful before deleting it."""
references(kg::KnowledgeGraph, x) = referrers(kg, x)

function Base.show(io::IO, ::MIME"text/plain", kg::KnowledgeGraph)
    counts = Dict{String,Int}()
    for r in kg
        counts[class_name(r)] = get(counts, class_name(r), 0) + 1
    end
    print(io, "KnowledgeGraph with ", length(kg), " resources", isindexed(kg) ? " (indexed)" : "")
    for (k, v) in sort!(collect(counts))
        print(io, "\n  ", rpad(k, 15), v)
    end
end
Base.show(io::IO, kg::KnowledgeGraph) = print(io, "KnowledgeGraph(", length(kg), " resources)")
