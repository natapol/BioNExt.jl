# Graph statistics: composition, degrees, connectivity, structure, sources and annotation coverage.

"""
    GraphStats

Statistics of a knowledge graph, computed by `graph_stats`. Fields:

- `resources`, `classes` (class ⇒ count), `relationship_types` (type ⇒ count)
- `nodes`, `edges`, `isolated` (nodes in no edge)
- `degree`, `in_degree`, `out_degree`: `(min, max, mean, median)` over nodes; a node's
  degree is the number of edges it takes part in (`degree_of`); in = as a target (right),
  out = as a source (left)
- `degree_distribution` (degree ⇒ number of nodes), `hubs` (id ⇒ degree, highest first)
- `components`, `largest_component`, `largest_fraction`: weakly connected components of
  the flattened network (`network_components`)
- `arcs` (distinct directed pairs in `network(kg)`), `density`, `reciprocity` (fraction of
  arcs whose reverse also exists), `self_loops`, `clustering` (mean local clustering
  coefficient over nodes with ≥ 2 neighbours, undirected; `nothing` if not computed)
- `pathways` and `complexes`: `(count, mean, max)` of their direct components / subunits;
  `pathway_depth`: deepest nesting of sub-pathways
- `organisms` (taxon ⇒ entities), `locations` (term ⇒ nodes), `id_prefixes` (prefix ⇒ resources)
- `sources` (provenance ⇒ entities), `multi_source` (entities with ≥ 2 provenance records),
  `source_overlap` ((a, b) ⇒ entities shared by both)
- `publications` (distinct references), `completeness` (annotation ⇒ fraction of applicable resources)
"""
struct GraphStats
    resources::Int
    classes::Vector{Pair{String,Int}}
    relationship_types::Vector{Pair{Symbol,Int}}
    nodes::Int
    edges::Int
    isolated::Int
    degree::NamedTuple
    in_degree::NamedTuple
    out_degree::NamedTuple
    degree_distribution::Vector{Pair{Int,Int}}
    hubs::Vector{Pair{String,Int}}
    components::Int
    largest_component::Int
    largest_fraction::Float64
    arcs::Int
    density::Float64
    reciprocity::Float64
    self_loops::Int
    clustering::Union{Float64,Nothing}
    pathways::NamedTuple
    pathway_depth::Int
    complexes::NamedTuple
    organisms::Vector{Pair{String,Int}}
    locations::Vector{Pair{Symbol,Int}}
    id_prefixes::Vector{Pair{String,Int}}
    sources::Vector{Pair{String,Int}}
    multi_source::Int
    source_overlap::Vector{Pair{Tuple{String,String},Int}}
    publications::Int
    completeness::Vector{Pair{String,Float64}}
end

summarize(v) = isempty(v) ? (min = 0, max = 0, mean = 0.0, median = 0.0) :
    (s = sort(v); (min = s[1], max = s[end], mean = sum(s) / length(s),
                   median = isodd(length(s)) ? float(s[(end + 1) ÷ 2]) : (s[end ÷ 2] + s[end ÷ 2 + 1]) / 2))
sizes(v) = isempty(v) ? (count = 0, mean = 0.0, max = 0) : (count = length(v), mean = sum(v) / length(v), max = maximum(v))
sorted_counts(d) = sort!(collect(d); by = p -> (-p.second, string(p.first)))

function id_prefix(id)
    isabsolute(id) && return "(IRI)"
    m = match(CURIE, id)
    m === nothing && return (occursin('/', id) ? "(minted) " * first(split(id, '/')) * "/" : "(bare)")
    return m[1]
end

# Mean local clustering coefficient of the undirected projection (nodes with ≥ 2 neighbours).
function mean_clustering(arcs)
    nb = Dict{String,Set{String}}()
    for (a, b) in arcs
        a == b && continue
        push!(get!(nb, a, Set{String}()), b)
        push!(get!(nb, b, Set{String}()), a)
    end
    tri = Dict{String,Int}()
    for (u, nu) in nb, v in nu
        u < v || continue
        small, large = length(nu) <= length(nb[v]) ? (nu, nb[v]) : (nb[v], nu)
        c = count(w -> w in large, small)
        tri[u] = get(tri, u, 0) + c
        tri[v] = get(tri, v, 0) + c
    end
    cs = [(get(tri, u, 0) / 2) / (length(n) * (length(n) - 1) / 2) for (u, n) in nb if length(n) >= 2]
    return isempty(cs) ? 0.0 : sum(cs) / length(cs)
end

"""
    graph_stats(kg; top = 10, clustering = :auto) -> GraphStats

Collect statistics of `kg` (see `GraphStats` for every field). `top` is the number of hubs
listed. The clustering coefficient is computed when `clustering = true`, or with `:auto`
when the graph has at most 200,000 nodes. Display the result for a report, or turn it into
rows with `stats_table` and save it with `write_stats`.
"""
function graph_stats(kg::KnowledgeGraph; top::Integer = 10, clustering = :auto)
    classes = Dict{String,Int}()
    foreach(r -> classes[class_name(r)] = get(classes, class_name(r), 0) + 1, kg)
    es, ns = edges(kg), nodes(kg)
    types = Dict{Symbol,Int}()
    deg, indeg, outdeg = Dict{String,Int}(), Dict{String,Int}(), Dict{String,Int}()
    for e in es
        types[e.type] = get(types, e.type, 0) + 1
        for x in unique([e.left; e.right])
            deg[x] = get(deg, x, 0) + 1
        end
        foreach(x -> outdeg[x] = get(outdeg, x, 0) + 1, unique(e.left))
        foreach(x -> indeg[x] = get(indeg, x, 0) + 1, unique(e.right))
    end
    nd = [get(deg, n.id, 0) for n in ns]
    dist = Dict{Int,Int}()
    foreach(d -> dist[d] = get(dist, d, 0) + 1, nd)
    hubs = first(sort!([n.id => get(deg, n.id, 0) for n in ns]; by = p -> (-p.second, p.first)), top)

    net = network(kg)
    pairs = Set{Tuple{String,String}}((a.source, a.target) for a in net.arcs)
    loops = count(((a, b),) -> a == b, pairs)
    proper = [(a, b) for (a, b) in pairs if a != b]
    n = length(net.nodes)
    recip = isempty(proper) ? 0.0 : count(((a, b),) -> (b, a) in pairs, proper) / length(proper)
    comps = network_components(kg)
    do_clust = clustering === true || (clustering === :auto && n <= 200_000)

    pws = pathways(kg)
    depth_memo = Dict{String,Int}()
    function depth(id, seen = Set{String}())
        haskey(depth_memo, id) && return depth_memo[id]
        id in seen && return 0
        push!(seen, id)
        subs = [c for c in kg[id].components if get(kg, c, nothing) isa PathwayModule]
        return depth_memo[id] = 1 + (isempty(subs) ? 0 : maximum(depth(c, seen) for c in subs))
    end
    cxs = resources(kg, MolecularComplex)

    orgs, locs, prefixes, srcs = Dict{String,Int}(), Dict{Symbol,Int}(), Dict{String,Int}(), Dict{String,Int}()
    overlap = Dict{Tuple{String,String},Int}()
    pubs = Set{String}()
    multi = 0
    for r in kg
        prefixes[id_prefix(r.id)] = get(prefixes, id_prefix(r.id), 0) + 1
        r isa Entity || continue
        if r.organism !== nothing
            label = (t = get(kg, r.organism, nothing)) isa Taxon && t.name !== nothing ? "$(r.organism) ($(t.name))" : r.organism
            orgs[label] = get(orgs, label, 0) + 1
        end
        r isa NodeEntity && foreach(l -> locs[l] = get(locs, l, 0) + 1, r.cellular_location)
        ps = unique(r.provenance)
        foreach(p -> srcs[p] = get(srcs, p, 0) + 1, ps)
        length(ps) >= 2 && (multi += 1)
        for i in eachindex(ps), j in i+1:length(ps)
            key = ps[i] < ps[j] ? (ps[i], ps[j]) : (ps[j], ps[i])
            overlap[key] = get(overlap, key, 0) + 1
        end
        union!(pubs, r.publications)
    end

    frac(xs, f) = isempty(xs) ? NaN : count(f, xs) / length(xs)
    ents = resources(kg, Entity)
    named = filter(r -> !(r isa EdgeEntity), ents)
    molecules = filter(r -> r isa Union{Gene,Protein,Dna,Rna}, ents)
    prots = resources(kg, Protein)
    completeness = [
        "name (nodes, pathways, taxa)" => frac(named, r -> r.name !== nothing),
        "organism (genes, proteins, DNA, RNA)" => frac(molecules, r -> r.organism !== nothing),
        "sequence (proteins)" => frac(prots, r -> r.sequence !== nothing),
        "cellular location (nodes)" => frac(ns, r -> !isempty(r.cellular_location)),
        "provenance (entities)" => frac(ents, r -> !isempty(r.provenance)),
        "publications (edges)" => frac(es, r -> !isempty(r.publications)),
        "confidence (edges)" => frac(es, r -> r.confidence !== nothing),
        "external ids (nodes)" => frac(ns, r -> !isempty(r.external_ids)),
    ]

    return GraphStats(
        length(kg), sorted_counts(classes), sorted_counts(types), length(ns), length(es), count(==(0), nd),
        summarize(nd), summarize([get(indeg, x.id, 0) for x in ns]), summarize([get(outdeg, x.id, 0) for x in ns]),
        sort!(collect(dist)), hubs,
        length(comps), isempty(comps) ? 0 : length(comps[1]), isempty(comps) || n == 0 ? 0.0 : length(comps[1]) / n,
        length(proper), n > 1 ? length(proper) / (n * (n - 1)) : 0.0, recip, loops,
        do_clust ? mean_clustering(pairs) : nothing,
        sizes([length(p.components) for p in pws]), isempty(pws) ? 0 : maximum(depth(p.id) for p in pws),
        sizes([length(c.subunits) + length(c.subcomplexes) for c in cxs]),
        sorted_counts(orgs), sorted_counts(locs), sorted_counts(prefixes), sorted_counts(srcs), multi,
        sorted_counts(overlap), length(pubs), completeness)
end

fmt(x::AbstractFloat) = isnan(x) ? "–" : x == round(x) && abs(x) < 1e6 ? string(Int(x)) : string(round(x; sigdigits = 4))
fmt(x) = string(x)
pct(x) = isnan(x) ? "–" : string(round(100x; digits = 1), "%")

"""
    stats_table(s::GraphStats) -> Vector{NamedTuple{(:section, :metric, :value)}}

The statistics as rows, for printing, saving (`write_stats`) or a DataFrame.
"""
function stats_table(s::GraphStats)
    rows = NamedTuple{(:section, :metric, :value),Tuple{String,String,Any}}[]
    add(sec, m, v) = push!(rows, (section = sec, metric = string(m), value = v))
    add("overview", "resources", s.resources)
    add("overview", "nodes", s.nodes)
    add("overview", "edges", s.edges)
    foreach(p -> add("classes", p.first, p.second), s.classes)
    foreach(p -> add("relationship types", p.first, p.second), s.relationship_types)
    add("degree", "isolated nodes", s.isolated)
    for (name, d) in (("degree", s.degree), ("in-degree", s.in_degree), ("out-degree", s.out_degree)), k in keys(d)
        add("degree", "$name $k", d[k])
    end
    foreach(p -> add("hubs", p.first, p.second), s.hubs)
    foreach(p -> add("degree distribution", p.first, p.second), s.degree_distribution)
    add("network", "components", s.components)
    add("network", "largest component", s.largest_component)
    add("network", "largest component fraction", s.largest_fraction)
    add("network", "arcs", s.arcs)
    add("network", "density", s.density)
    add("network", "reciprocity", s.reciprocity)
    add("network", "self-loops", s.self_loops)
    add("network", "mean clustering coefficient", something(s.clustering, NaN))
    for (name, t) in (("pathways", s.pathways), ("complexes", s.complexes)), k in keys(t)
        add("structure", "$name $k", t[k])
    end
    add("structure", "pathway nesting depth", s.pathway_depth)
    foreach(p -> add("organisms", p.first, p.second), s.organisms)
    foreach(p -> add("cellular locations", p.first, p.second), s.locations)
    foreach(p -> add("identifier prefixes", p.first, p.second), s.id_prefixes)
    foreach(p -> add("sources", p.first, p.second), s.sources)
    add("sources", "entities from ≥ 2 sources", s.multi_source)
    foreach(p -> add("source overlap", "$(p.first[1]) & $(p.first[2])", p.second), s.source_overlap)
    add("annotation", "distinct publications", s.publications)
    foreach(p -> add("annotation", p.first, p.second), s.completeness)
    return rows
end

function Base.show(io::IO, ::MIME"text/plain", s::GraphStats)
    line(k, v) = println(io, "  ", rpad(k, 37), " ", v)
    lst(ps; n = 8) = join(("$(p.first) $(p.second)" for p in first(ps, n)), ", ") * (length(ps) > n ? ", …" : "")
    println(io, "Knowledge graph statistics")
    println(io, "Composition")
    line("resources / nodes / edges", "$(s.resources) / $(s.nodes) / $(s.edges)")
    line("classes", lst(s.classes; n = 20))
    line("relationship types", lst(s.relationship_types))
    println(io, "Degree (edges per node)")
    d = s.degree
    line("min / median / mean / max", "$(d.min) / $(fmt(d.median)) / $(fmt(d.mean)) / $(d.max)")
    line("isolated nodes", s.isolated)
    line("hubs", lst(s.hubs; n = 5))
    println(io, "Network")
    line("components (largest)", "$(s.components) ($(s.largest_component) = $(pct(s.largest_fraction)) of the network)")
    line("arcs / density / reciprocity", "$(s.arcs) / $(fmt(s.density)) / $(pct(s.reciprocity))")
    line("self-loops", s.self_loops)
    line("mean clustering coefficient", s.clustering === nothing ? "not computed" : fmt(s.clustering))
    println(io, "Structure")
    line("pathways (mean / max size, depth)", "$(s.pathways.count) ($(fmt(s.pathways.mean)) / $(s.pathways.max), $(s.pathway_depth))")
    line("complexes (mean / max subunits)", "$(s.complexes.count) ($(fmt(s.complexes.mean)) / $(s.complexes.max))")
    println(io, "Context")
    line("organisms", isempty(s.organisms) ? "–" : lst(s.organisms; n = 4))
    line("cellular locations", isempty(s.locations) ? "–" : lst(s.locations; n = 5))
    line("identifier prefixes", lst(s.id_prefixes; n = 6))
    println(io, "Sources")
    line("entities per source", isempty(s.sources) ? "–" : lst(s.sources))
    line("entities from ≥ 2 sources", s.multi_source)
    isempty(s.source_overlap) || line("largest overlaps", lst([("$(p.first[1]) & $(p.first[2])" => p.second) for p in s.source_overlap]; n = 3))
    println(io, "Annotation coverage")
    line("distinct publications", s.publications)
    for (k, v) in s.completeness
        line(k, pct(v))
    end
end
Base.show(io::IO, s::GraphStats) = print(io, "GraphStats(", s.nodes, " nodes, ", s.edges, " edges)")

"""
    write_stats(io_or_path, s::GraphStats; format = :markdown)

Save statistics as a Markdown report (`:markdown`, `.md`), a table (`:tsv`, `.tsv`;
`:csv`, `.csv`) or JSON (`:json`, `.json`). With a path, the format follows the extension.
"""
function write_stats(io::IO, s::GraphStats; format::Symbol = :markdown)
    rows = stats_table(s)
    val(v) = v isa AbstractFloat ? fmt(v) : string(v)
    if format === :markdown
        println(io, "# Knowledge graph statistics")
        for sec in unique(r.section for r in rows)
            println(io, "\n## ", uppercasefirst(sec), "\n\n| Metric | Value |\n|---|---|")
            for r in rows
                r.section == sec && println(io, "| ", replace(r.metric, "|" => "\\|"), " | ",
                                           occursin("fraction", r.metric) || r.section == "annotation" && r.value isa AbstractFloat ?
                                           pct(r.value) : val(r.value), " |")
            end
        end
    elseif format in (:tsv, :csv)
        d = format === :tsv ? '\t' : ','
        println(io, join(["section", "metric", "value"], d))
        foreach(r -> println(io, join(quote_cell.([r.section, r.metric, val(r.value)], d), d)), rows)
    elseif format === :json
        sections = Pair{String,Any}[]
        for sec in unique(r.section for r in rows)
            push!(sections, sec => Pair{String,Any}[r.metric => (r.value isa AbstractFloat && !isfinite(r.value) ? nothing : r.value)
                                                    for r in rows if r.section == sec])
        end
        json_print(io, sections)
        println(io)
    else
        throw(ArgumentError("unknown statistics format :$format (use :markdown, :tsv, :csv or :json)"))
    end
end

function write_stats(path::AbstractString, s::GraphStats;
                     format::Symbol = get(Dict(".md" => :markdown, ".tsv" => :tsv, ".csv" => :csv, ".json" => :json),
                                          lowercase(splitext(path)[2]), :markdown))
    open(io -> write_stats(io, s; format), path, "w")
    return path
end
