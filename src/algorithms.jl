# Graph algorithms that understand BioNExt edges: relationship types (with the term
# hierarchy), direction, and edges whose participants are other edges or pathways.
# Names avoid Graphs.jl's exports (`degree`, `connected_components`, `neighborhood`)
# so both packages can be loaded together; `to_graph` hands the graph to Graphs.jl.

"""
    degree_of(kg, x; direction = :both, type = nothing) -> Int

Number of edges `x` takes part in. `direction = :out` counts edges where `x` is on the
left, `:in` where it is on the right; `type` restricts to a relationship type and its subtypes.
"""
function degree_of(kg::KnowledgeGraph, x; direction::Symbol = :both, type::Union{Symbol,Nothing} = nothing)
    side = direction === :out ? :left : direction === :in ? :right : :any
    return length(edges_of(kg, x; side, type))
end

# One traversal step, matching the arcs of `network(kg; reify, type)`.
function steps(kg, x, direction, type, reify)
    reify || return neighbors(kg, x; direction, type)
    out = String[]
    # node -> edge: edges having x on the left (out) / right (in)
    direction !== :in && append!(out, (e.id for e in edges_of(kg, x; side = :left, type)))
    direction !== :out && append!(out, (e.id for e in edges_of(kg, x; side = :right, type)))
    # edge -> node: an edge leads on to its right participants (out) / back to its left ones (in)
    e = get(kg, x, nothing)
    if e isa EdgeEntity && (type === nothing || is_a(e.type, type))
        direction !== :in && append!(out, e.right)
        direction !== :out && append!(out, e.left)
    end
    return unique!(filter!(!=(x), out))
end

"""
    shortest_path(kg, from, to; direction = :out, type = nothing, reify = false,
                  max_hops = typemax(Int)) -> Union{Vector{String}, Nothing}

Fewest-hop path of ids from `from` to `to`, or `nothing` if there is none. It follows
the same arcs as `network(kg; reify, type)` (and so `to_graph`): left → right for
`direction = :out`, right → left for `:in`, either way for `:both`.

With `reify = false`, each step goes from a participant to a participant of the same
edge; an edge that is itself a participant (the target of an `Influence`) is an end point.
With `reify = true`, edges are steps of their own (`A → e1 → B`), so a path can run
through an influence on an edge to that edge's participants.
"""
function shortest_path(kg::KnowledgeGraph, from, to; direction::Symbol = :out,
                       type::Union{Symbol,Nothing} = nothing, reify::Bool = false,
                       max_hops::Integer = typemax(Int))
    s, t = _id(from), _id(to)
    s == t && return [s]
    prev = Dict{String,String}(s => s)
    frontier = [s]
    for _ in 1:max_hops
        next = String[]
        for x in frontier, y in steps(kg, x, direction, type, reify)
            haskey(prev, y) && continue
            prev[y] = x
            if y == t
                path = [t]
                while path[1] != s
                    pushfirst!(path, prev[path[1]])
                end
                return path
            end
            push!(next, y)
        end
        isempty(next) && return nothing
        frontier = next
    end
    return nothing
end

"""
    network_components(kg; type = nothing, reify = false) -> Vector{Vector{String}}

Weakly connected components (ignoring direction) of `network(kg; reify, type)`, largest
first, each in insertion order. Every node entity is in exactly one component (isolated
nodes as singletons), as is every id that occurs as an edge participant. With
`reify = true`, every edge is a member too and is joined to its participants, so an
influence on an edge connects its source to that edge's participants.
"""
function network_components(kg::KnowledgeGraph; type::Union{Symbol,Nothing} = nothing, reify::Bool = false)
    parent = Dict{String,String}()
    find(x) = (p = get!(parent, x, x); p == x ? x : (parent[x] = find(p)))
    unite(a, b) = (ra = find(a); rb = find(b); ra != rb && (parent[rb] = ra))
    foreach(n -> find(n.id), nodes(kg))
    for e in edges(kg)
        type === nothing || is_a(e.type, type) || continue
        ids = reify ? [e.id; e.left; e.right] : [e.left; e.right]
        isempty(ids) && continue
        foreach(i -> unite(ids[1], i), ids)
    end
    groups = Dict{String,Vector{String}}()
    order = Dict(id => i for (i, id) in enumerate(keys(kg)))
    for x in sort!(collect(keys(parent)); by = x -> get(order, x, typemax(Int)))
        push!(get!(groups, find(x), String[]), x)
    end
    return sort!(collect(values(groups)); by = g -> (-length(g), get(order, g[1], typemax(Int))))
end

"""
    subgraph(kg, seeds; hops = 0, type = nothing, include_context = true) -> KnowledgeGraph

A new graph containing copies of `seeds` and everything within `hops` edge steps of them.
Each step adds the edges of the current ids and all their participants. With
`include_context = true`, resources the kept ones refer to without being part of the
network (taxa, provenance, cofactors, complex subunits, ...) are copied as well, so the
result validates without dangling-reference warnings for them. Pathways are kept only
if they are seeds; their `components` are trimmed to what was kept.
"""
function subgraph(kg::KnowledgeGraph, seeds; hops::Integer = 0, type::Union{Symbol,Nothing} = nothing,
                  include_context::Bool = true)
    keep = Set{String}(_id(s) for s in (seeds isa Union{AbstractString,Resource} ? [seeds] : seeds))
    frontier = collect(keep)
    for _ in 1:hops
        next = String[]
        for x in frontier, e in edges_of(kg, x; type), y in [e.id; e.left; e.right]
            y in keep || (push!(keep, y); push!(next, y))
        end
        isempty(next) && break
        frontier = next
    end
    if include_context
        stack = collect(keep)
        while !isempty(stack)
            r = get(kg, pop!(stack), nothing)
            r === nothing && continue
            for p in schema(r), v in field_values(r, p.field)
                p.kind === :ref && p.field !== :components || continue
                target = get(kg, v, nothing)
                target === nothing || target isa EdgeEntity || target isa PathwayModule || v in keep ||
                    (push!(keep, v); push!(stack, v))
            end
        end
    end
    sub = KnowledgeGraph(; base = kg.base, prefixes = copy(kg.prefixes))
    for id in keys(kg)
        id in keep || continue
        r = deepcopy(kg[id])
        r isa PathwayModule && filter!(in(keep), r.components)
        add!(sub, r)
    end
    isindexed(kg) && index!(sub)
    return sub
end

"""
    to_graph(kg; reify = false, directed = true, type = nothing)
        -> (graph, ids::Vector{String}, vertex::Dict{String,Int})

Convert to a Graphs.jl graph so its algorithms (centrality, communities, flows, ...)
can run on the network. Requires `using Graphs` (or `import Graphs`), which loads the
package extension. Vertices are the nodes of `network(kg; reify)`: vertex `i` is
`ids[i]`, and `vertex[id]` is the reverse lookup. `type` keeps only edges of that
relationship type (and its subtypes). The graph is a `SimpleDiGraph`, or a
`SimpleGraph` with `directed = false`; parallel arcs collapse into one.
"""
to_graph(args...; kwargs...) =
    error("to_graph needs Graphs.jl: run `using Graphs` (or `import Graphs`) first")
