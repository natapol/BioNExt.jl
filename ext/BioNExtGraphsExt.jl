module BioNExtGraphsExt

using BioNExt: BioNExt, KnowledgeGraph, network
using Graphs: SimpleDiGraph, SimpleGraph, add_edge!

function BioNExt.to_graph(kg::KnowledgeGraph; reify::Bool = false, directed::Bool = true,
                          type::Union{Symbol,Nothing} = nothing)
    net = network(kg; reify, type)
    vertex = Dict(id => i for (i, id) in enumerate(net.nodes))
    g = directed ? SimpleDiGraph(length(net.nodes)) : SimpleGraph(length(net.nodes))
    for a in net.arcs
        add_edge!(g, vertex[a.source], vertex[a.target])
    end
    return (graph = g, ids = net.nodes, vertex = vertex)
end

end
