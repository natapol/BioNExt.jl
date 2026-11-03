# 4. Querying

All query functions take either an id string or a resource, and return results in insertion order.
By default they scan the graph on every call. For anything beyond a few thousand resources,
call `index!(kg)` once and they become lookups (see [Indexing](#indexing)).

The examples on this page use the following graph:

```julia
using BioNExt

kg = KnowledgeGraph()
add!(kg,
    Gene(id = "HGNC:11998", name = "TP53"),
    Gene(id = "HGNC:6973", name = "MDM2"),
    Gene(id = "HGNC:1784", name = "CDKN1A"),
    Protein(id = "UniProt:P04637", name = "p53", cellular_location = [:Nucleus, :Cytoplasm]),
    Protein(id = "UniProt:Q00987", name = "MDM2", cellular_location = [:Nucleus]),
    Protein(id = "UniProt:Q13315", name = "ATM", cellular_location = [:Nucleus]),
    Protein(id = "UniProt:P04637_p.R175H", variant_of = ["UniProt:P04637"], hgvs = "p.Arg175His"),
    MolecularComplex(id = "CPX-p53-MDM2", subunits = ["UniProt:P04637", "UniProt:Q00987"]),
    Interaction(id = "e1", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
    Influence(id = "e2", type = :Ubiquitination, left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
    Influence(id = "e3", type = :RegulateExpression, left = ["UniProt:P04637"], right = ["HGNC:6973"]),
    Influence(id = "e4", type = :RegulateExpression, left = ["UniProt:P04637"], right = ["HGNC:1784"]),
    Influence(id = "e5", type = :Phosphorylation, left = ["UniProt:Q13315"], right = ["UniProt:P04637"]),
    Interaction(id = "e6", type = :PositiveCorrelation, left = ["HGNC:6973"], right = ["HGNC:1784"]),
    Influence(id = "e7", type = :Inhibition, left = ["UniProt:Q13315"], right = ["e1"]),
    PathwayModule(id = "PW-feedback", components = ["e1", "e2", "e3"]),
    PathwayModule(id = "PW-p53", components = ["PW-feedback", "e4", "e5"]),
)
```

## Listing resources by class

`resources(kg, T)` returns every resource of type `T`. `T` may be abstract. `nodes`, `edges` and
`pathways` are shortcuts:

```julia
length(resources(kg, Protein))          # -> 4
length(resources(kg, NodeEntity))       # -> 8
length(edges(kg))                       # -> 7
[p.id for p in pathways(kg)]            # -> ["PW-feedback", "PW-p53"]
```

Filter them with ordinary Julia:

```julia
nuclear = [p.name for p in resources(kg, Protein) if :Nucleus in p.cellular_location]
# -> ["p53", "MDM2", "ATM"]

canonical = filter(p -> isempty(p.variant_of), resources(kg, Protein))
by_name = Dict(r.name => r.id for r in resources(kg, Entity) if r.name !== nothing)
by_name["ATM"]                          # -> "UniProt:Q13315"
```

## Edges of a resource

`edges_of(kg, x)` returns the edges in which `x` takes part, which is the ontology's
`bn:participatesInEdge`:

```julia
[e.id for e in edges_of(kg, "UniProt:P04637")]                 # -> ["e1", "e2", "e3", "e4", "e5"]
[e.id for e in edges_of(kg, "UniProt:P04637"; side = :left)]   # -> ["e3", "e4"]
[e.id for e in edges_of(kg, "UniProt:P04637"; side = :right)]  # -> ["e1", "e2", "e5"]
```

`type` keeps only edges whose relationship type **is, or descends from**, the given term. That
lets you ask for a whole family at once:

```julia
[e.id for e in edges_of(kg, "HGNC:6973"; type = :CoExpression)]   # -> ["e6"]  (PositiveCorrelation is a CoExpression)
```

Edges can also be targets. Here `e7` inhibits the binding event `e1`:

```julia
[e.id for e in edges_of(kg, "e1"; side = :right)]                 # -> ["e7"]
```

## Neighbours

`neighbors(kg, x)` returns the ids connected to `x` through any edge. `direction` follows the
edges left → right (`:out`), right → left (`:in`), or both ways (`:both`, the default):

```julia
neighbors(kg, "UniProt:P04637"; direction = :out)            # -> ["HGNC:6973", "HGNC:1784"]
neighbors(kg, "UniProt:P04637"; direction = :in)             # -> ["UniProt:Q00987", "UniProt:Q13315"]
neighbors(kg, "UniProt:Q13315")                              # -> ["UniProt:P04637", "e1"]
neighbors(kg, "UniProt:P04637"; type = :RegulateExpression)  # -> ["HGNC:6973", "HGNC:1784"]
```

## Pathways

```julia
[p.id for p in pathways_of(kg, "e1")]           # -> ["PW-feedback"]
components(kg, "PW-p53")                        # -> ["PW-feedback", "e4", "e5"]
components(kg, "PW-p53"; recursive = true)      # -> ["PW-feedback", "e4", "e5", "e1", "e2", "e3"]
```

`pathways_of` returns only the pathways that list `x` directly. To find every enclosing pathway,
walk upwards:

```julia
function enclosing_pathways(kg, x)
    out, frontier = String[], [x]
    while !isempty(frontier)
        for p in pathways_of(kg, popfirst!(frontier))
            p.id in out || (push!(out, p.id); push!(frontier, p.id))
        end
    end
    return out
end
enclosing_pathways(kg, "e1")                    # -> ["PW-feedback", "PW-p53"]
```

## Complexes, variants and ortholog groups

```julia
[c.id for c in complexes_of(kg, "UniProt:Q00987")]   # -> ["CPX-p53-MDM2"]
[v.id for v in variants(kg, "UniProt:P04637")]       # -> ["UniProt:P04637_p.R175H"]
ortholog_groups_of(kg, "UniProt:P04637")             # -> OrthologGroup[]
```

## Who refers to a resource?

`references(kg, x)` lists every `(id => field)` that holds `x` in a reference field. Use it before
you delete or rename `x`:

```julia
references(kg, "UniProt:Q00987")
# -> ["CPX-p53-MDM2" => :subunits, "e1" => :left, "e2" => :left]
```

## Indexing

Without an index, each query scans the graph. `index!(kg)` builds a reverse index of every
reference, and from then on `edges_of`, `neighbors`, `degree_of`, `pathways_of`,
`complexes_of`, `variants`, `ortholog_groups_of` and `references` are direct lookups. On a
graph of 100,000 proteins and 300,000 interactions, a `neighbors` call drops from about
80 ms to about 0.04 ms. Building the index takes about 0.4 s.

```julia
index!(kg)
isindexed(kg)                                    # -> true
[e.id for e in edges_of(kg, "UniProt:P04637")]   # -> ["e1", "e2", "e3", "e4", "e5"]
```

The index updates itself when the graph changes through `add!`, `delete!`, `merge!`,
`load!`, `resolve!` or `integrate!`. It **can't see edits made directly to fields**. After
changing a reference field in place, call `reindex!`:

```julia
add!(kg, Protein(id = "UniProt:P38936", name = "p21"))
push!(kg["e4"].right, "UniProt:P38936")          # direct edit: the index doesn't see it
[e.id for e in edges_of(kg, "UniProt:P38936")]   # -> String[]
reindex!(kg)
[e.id for e in edges_of(kg, "UniProt:P38936")]   # -> ["e4"]
```

A stale index can miss new references, but it never returns references that have been
removed, because every result is checked against the current fields. `unindex!(kg)` turns
indexing off again. A good habit: build or load the graph, call `index!` once, and keep
direct edits to non-reference fields such as `name` or `confidence`, which don't affect the index.

## Graph algorithms

Four algorithms understand BioNExt edges: relationship types (including subtypes), direction,
and edges whose participants are other edges. Each accepts `type` to restrict the edges used.

```julia
degree_of(kg, "UniProt:P04637")                          # -> 5
degree_of(kg, "UniProt:P04637"; direction = :out)        # -> 2
degree_of(kg, "UniProt:P04637"; type = :RegulateExpression)  # -> 2

shortest_path(kg, "UniProt:Q13315", "HGNC:1784")
# -> ["UniProt:Q13315", "UniProt:P04637", "HGNC:1784"]
shortest_path(kg, "HGNC:1784", "UniProt:Q13315")         # -> nothing  (edges are directed)
shortest_path(kg, "HGNC:1784", "UniProt:Q13315"; direction = :both) |> length   # -> 3

network_components(kg)
# -> [["HGNC:6973", "HGNC:1784", "UniProt:P04637", "UniProt:Q00987", "UniProt:Q13315", "e1", "UniProt:P38936"],
#     ["HGNC:11998"], ["UniProt:P04637_p.R175H"], ["CPX-p53-MDM2"]]
```

`network_components` returns weakly connected components, largest first, each in insertion
order. `e1` appears because `e7` targets it. Resources that aren't connected by any edge
(the TP53 gene, the variant, the complex) are singletons.

### Edges that point at edges: `reify`

`e7` says that ATM inhibits the *binding* `e1` between MDM2 and p53. By default, an edge that is
a participant is an end point, as in the network that `to_graph(kg)` hands to Graphs.jl. So
there is no path from ATM through `e1` to p53. With `reify = true`, every edge becomes a step of
its own (`left → edge → right`), so the path exists:

```julia
shortest_path(kg, "UniProt:Q13315", "UniProt:P04637"; type = :Inhibition)                 # -> nothing
shortest_path(kg, "UniProt:Q13315", "UniProt:P04637"; type = :Inhibition, reify = true) === nothing  # -> true
shortest_path(kg, "UniProt:Q13315", "UniProt:P04637"; reify = true)
# -> ["UniProt:Q13315", "e5", "UniProt:P04637"]
p = shortest_path(kg, "UniProt:Q13315", "UniProt:Q00987"; reify = true, direction = :both)
p                                                # -> ["UniProt:Q13315", "e7", "e1", "UniProt:Q00987"]
```

The second line finds nothing because `type = :Inhibition` also filters out `e1` (a `:Binding`)
as a step. `network_components(kg; reify = true)` likewise joins every edge to its participants.

### Sub-graphs

`subgraph(kg, seeds; hops)` copies the seeds and everything within `hops` edge steps into a new
graph. Each step adds the edges of the current ids and all their participants. The taxa,
provenance records and other resources that the kept ones refer to come along as well (turn
this off with `include_context = false`), so the result validates cleanly. Here ATM's edges are
`e5` and `e7`; `e7` targets `e1`, whose other participant, MDM2, comes along as context:

```julia
sub = subgraph(kg, "UniProt:Q13315"; hops = 1)
sort(collect(keys(sub)))
# -> ["UniProt:P04637", "UniProt:Q00987", "UniProt:Q13315", "e1", "e5", "e7"]
sort(collect(keys(subgraph(kg, "PW-feedback"))))   # a pathway seed, trimmed to what is kept
# -> ["PW-feedback"]
```

`hops = 0` copies exactly the seeds, plus their context. The copies are independent of `kg`.

### Graphs.jl

For anything else (centrality, communities, flows, layouts, isomorphism), hand the network to
[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl). `to_graph` is available once Graphs.jl
is loaded. It's an optional extension, and BioNExt.jl doesn't depend on Graphs.jl.

<!-- norun -->
```julia
import Graphs                       # `import`, not `using`: see the note below

g = to_graph(kg)                    # (graph, ids, vertex)
g.graph                             # a Graphs.SimpleDiGraph
bc = Graphs.betweenness_centrality(g.graph)
top = g.ids[argmax(bc)]             # the most central id
Graphs.has_path(g.graph, g.vertex["UniProt:Q13315"], g.vertex["HGNC:1784"])

to_graph(kg; reify = true, directed = false, type = :RegulateExpression)
```

Vertex `i` is `g.ids[i]`, and `g.vertex[id]` gives the vertex number. The vertices and arcs are
exactly those of `network(kg; reify, type)`. Parallel arcs collapse into one.

**Name clash.** Both packages export `neighbors` and `edges`. If you load them with
`using BioNExt, Graphs`, Julia makes you qualify those two names
(`BioNExt.neighbors`, `Graphs.neighbors`). Write `import Graphs` and prefix Graphs.jl's
functions to avoid this. BioNExt.jl's own algorithm names (`degree_of`, `shortest_path`,
`network_components`, `subgraph`) were chosen not to clash.

## Statistics

```julia
hubs = sort([(n.id, degree_of(kg, n)) for n in nodes(kg)]; by = last, rev = true)
first(hubs)                                     # -> ("UniProt:P04637", 5)

type_counts = Dict{Symbol,Int}()
for e in edges(kg)
    type_counts[e.type] = get(type_counts, e.type, 0) + 1
end
type_counts[:RegulateExpression]                # -> 2
```
