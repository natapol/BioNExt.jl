# Network views of the graph for network tools: SIF (Cytoscape), GraphML, Graphviz DOT.
#
# BioNExt edges are n-ary (several left and right participants) and may target other
# edges or pathways. Two flattenings are offered:
#   reify = false  one arc per (left, right) pair, labelled with the relationship type
#   reify = true   each edge entity becomes a node, with left -> edge -> right arcs (lossless)

"""One arc of a flattened network; `role` is `:pair`, `:left` or `:right`."""
struct Arc
    source::String
    target::String
    type::Symbol
    edge::String
    role::Symbol
end

"""
    network(kg; reify = false, type = nothing) -> (nodes = Vector{String}, arcs = Vector{Arc})

Flatten the graph to a simple directed network (see `reify` above). Nodes are all
node entities plus anything that appears as an arc endpoint. `type` keeps only edges
whose relationship type `is_a` that term.
"""
function network(kg::KnowledgeGraph; reify::Bool = false, type::Union{Symbol,Nothing} = nothing)
    arcs = Arc[]
    for e in edges(kg)
        type === nothing || is_a(e.type, type) || continue
        if reify
            foreach(l -> push!(arcs, Arc(l, e.id, e.type, e.id, :left)), e.left)
            foreach(r -> push!(arcs, Arc(e.id, r, e.type, e.id, :right)), e.right)
        else
            for l in e.left, r in e.right
                push!(arcs, Arc(l, r, e.type, e.id, :pair))
            end
        end
    end
    ids = [n.id for n in nodes(kg)]
    seen = Set(ids)
    for a in arcs, x in (a.source, a.target)
        x in seen || (push!(seen, x); push!(ids, x))
    end
    return (nodes = ids, arcs = arcs)
end

arc_label(a::Arc) = a.role === :pair ? string(a.type) : string(a.role)

# --- SIF -----------------------------------------------------------------------

"""
    write_sif(io_or_path, kg; reify = false)

Simple Interaction Format: `source<TAB>relationship<TAB>target`, isolated nodes on their own line.
"""
function write_sif(io::IO, kg::KnowledgeGraph; reify::Bool = false)
    net = network(kg; reify)
    connected = Set{String}()
    for a in net.arcs
        println(io, a.source, '\t', arc_label(a), '\t', a.target)
        push!(connected, a.source, a.target)
    end
    foreach(n -> n in connected || println(io, n), net.nodes)
end

function match_relationship(rel::AbstractString, relationships)
    haskey(relationships, rel) && return Symbol(relationships[rel])
    s = Symbol(rel)
    is_term(s) && return s
    key = lowercase(replace(rel, r"[\s_\-]" => ""))
    for t in terms(:RelationshipType)
        (lowercase(string(t)) == key || lowercase(replace(term(t).label, r"[\s_\-]" => "")) == key) && return t
    end
    return s   # left as-is; `validate` reports unknown terms
end

"""
    read_sif(io_or_path; node_type = Protein, edge_type = Interaction,
             relationships = Dict(), into = nothing) -> KnowledgeGraph

Import a SIF file (`source rel target1 target2 ...`, tab- or space-separated).
Unknown nodes are created as `node_type`, each arc becomes an `edge_type` edge.
Relationship names are matched to vocabulary terms by id or label, case-insensitively
(`"activation"` → `:Activation`); `relationships` overrides, e.g. `Dict("pp" => :Binding)`.
"""
function read_sif(io::IO; node_type::Type{<:NodeEntity} = Protein, edge_type::Type{<:EdgeEntity} = Interaction,
                  relationships = Dict{String,Symbol}(), into::Union{KnowledgeGraph,Nothing} = nothing)
    kg = into === nothing ? KnowledgeGraph() : into
    ensure!(id) = haskey(kg, id) || add!(kg, node_type(; id))
    n = 0
    for line in eachline(io)
        isempty(strip(line)) && continue
        tok = occursin('\t', line) ? strip.(split(line, '\t')) : split(strip(line))
        tok = String.(filter(!isempty, tok))
        ensure!(tok[1])
        length(tok) == 2 && throw(ArgumentError("SIF line without target: \"$line\""))
        length(tok) >= 3 || continue
        type = match_relationship(tok[2], relationships)
        for target in tok[3:end]
            ensure!(target)
            while haskey(kg, "edge-$(n += 1)"); end
            add!(kg, edge_type(; id = "edge-$n", type, left = [tok[1]], right = [target]))
        end
    end
    return kg
end

# --- GraphML -------------------------------------------------------------------

"""
    write_graphml(io_or_path, kg; reify = false)

GraphML for Cytoscape, Gephi, yEd or Graphs.jl. Node attributes: `class`, `name`;
edge attributes: `type`, `edge_id`, `class`, `role`, `weight`, `confidence`.
"""
function write_graphml(io::IO, kg::KnowledgeGraph; reify::Bool = false)
    net = network(kg; reify)
    x = xml_escape
    println(io, """<?xml version="1.0" encoding="UTF-8"?>
    <graphml xmlns="http://graphml.graphdrawing.org/xmlns">
      <key id="class" for="node" attr.name="class" attr.type="string"/>
      <key id="name" for="node" attr.name="name" attr.type="string"/>
      <key id="type" for="edge" attr.name="type" attr.type="string"/>
      <key id="edge_id" for="edge" attr.name="edge_id" attr.type="string"/>
      <key id="edge_class" for="edge" attr.name="class" attr.type="string"/>
      <key id="role" for="edge" attr.name="role" attr.type="string"/>
      <key id="weight" for="edge" attr.name="weight" attr.type="double"/>
      <key id="confidence" for="edge" attr.name="confidence" attr.type="double"/>
      <graph id="BioNExt" edgedefault="directed">""")
    for id in net.nodes
        r = get(kg, id, nothing)
        print(io, "    <node id=\"", x(id), "\">")
        r === nothing || print(io, "<data key=\"class\">", class_name(r), "</data>")
        r isa Entity && r.name !== nothing && print(io, "<data key=\"name\">", x(r.name), "</data>")
        println(io, "</node>")
    end
    for a in net.arcs
        e = kg[a.edge]
        print(io, "    <edge source=\"", x(a.source), "\" target=\"", x(a.target), "\">",
              "<data key=\"type\">", a.type, "</data><data key=\"edge_id\">", x(a.edge), "</data>",
              "<data key=\"edge_class\">", class_name(e), "</data><data key=\"role\">", a.role, "</data>")
        e.weight === nothing || print(io, "<data key=\"weight\">", e.weight, "</data>")
        e.confidence === nothing || print(io, "<data key=\"confidence\">", e.confidence, "</data>")
        println(io, "</edge>")
    end
    println(io, "  </graph>\n</graphml>")
end

# --- DOT -----------------------------------------------------------------------

const DOT_SHAPES = Dict("Gene" => "box", "Dna" => "box", "Rna" => "note", "Protein" => "ellipse",
                        "SmallMolecule" => "diamond", "Complex" => "octagon", "OrthologGroup" => "folder",
                        "PathwayModule" => "tab", "Interaction" => "point", "Influence" => "point",
                        "Reaction" => "point")

dot_quote(s) = "\"" * replace(string(s), "\\" => "\\\\", "\"" => "\\\"") * "\""

"""
    write_dot(io_or_path, kg; reify = false)

Graphviz DOT; render with e.g. `dot -Tsvg graph.dot -o graph.svg`.
"""
function write_dot(io::IO, kg::KnowledgeGraph; reify::Bool = false)
    net = network(kg; reify)
    println(io, "digraph BioNExt {\n  rankdir=LR;\n  node [fontname=\"Helvetica\", fontsize=10];\n  edge [fontname=\"Helvetica\", fontsize=9];")
    for id in net.nodes
        r = get(kg, id, nothing)
        cls = r === nothing ? "" : class_name(r)
        label = r isa Entity && r.name !== nothing ? r.name : id
        shape = get(DOT_SHAPES, cls, "ellipse")
        attrs = shape == "point" ? "shape=point, xlabel=$(dot_quote(r.type))" :
                "label=$(dot_quote(label)), shape=$shape, tooltip=$(dot_quote(id))"
        println(io, "  ", dot_quote(id), " [", attrs, "];")
    end
    for a in net.arcs
        style = is_a(a.type, :Inhibition) || is_a(a.type, :Repression) || is_a(a.type, :Suppression) ?
                ", arrowhead=tee" : ""
        label = a.role === :pair ? "label=$(dot_quote(a.type))" : "label=\"\""
        println(io, "  ", dot_quote(a.source), " -> ", dot_quote(a.target), " [", label, style, "];")
    end
    println(io, "}")
end
