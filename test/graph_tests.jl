using Random

@testset "indexing" begin
    @testset "indexed and scanned queries agree under random edits" begin
        rng = Random.MersenneTwister(7)
        kg = KnowledgeGraph()
        ids = ["P$i" for i in 1:40]
        foreach(id -> add!(kg, Protein(id = id)), ids)
        add!(kg, MolecularComplex(id = "C1", subunits = ids[1:3]),
                 OrthologGroup(id = "OG", members = ids[4:6]),
                 Protein(id = "V1", hgvs = "p.X1Y", variant_of = ["P1"]))
        types = [:Binding, :Activation, :Inhibition, :PositiveCorrelation]
        queries(kg, x) = (
            [e.id for e in edges_of(kg, x)], [e.id for e in edges_of(kg, x; side = :left)],
            [e.id for e in edges_of(kg, x; type = :CoExpression)], neighbors(kg, x),
            neighbors(kg, x; direction = :out), [p.id for p in pathways_of(kg, x)],
            [c.id for c in complexes_of(kg, x)], [v.id for v in variants(kg, x)],
            [g.id for g in ortholog_groups_of(kg, x)], references(kg, x), degree_of(kg, x))
        index!(kg)
        n = 0
        for step in 1:300
            r = rand(rng)
            if r < 0.6
                n += 1
                T = rand(rng, (Interaction, Influence))
                add!(kg, T(id = "e$n", type = rand(rng, types),
                           left = [rand(rng, ids)], right = [rand(rng, [ids; ["e$k" for k in 1:n-1]])]))
            elseif r < 0.7 && n > 0
                delete!(kg, "e$(rand(rng, 1:n))")
            elseif r < 0.8
                add!(kg, PathwayModule(id = "PW$(rand(rng, 1:5))", components = rand(rng, collect(keys(kg)), 3));
                     replace = true)
            elseif r < 0.9 && n > 0
                id = "e$(rand(rng, 1:n))"
                haskey(kg, id) && (push!(kg[id].left, rand(rng, ids)); reindex!(kg))
            end
            if step % 25 == 0
                for x in rand(rng, collect(keys(kg)), 10)
                    indexed = queries(kg, x)
                    unindex!(kg); scanned = queries(kg, x); index!(kg)
                    @test indexed == scanned
                end
            end
        end
    end

    @testset "lifecycle" begin
        kg = KnowledgeGraph()
        add!(kg, Protein(id = "A"), Protein(id = "B"), Interaction(id = "e", type = :Binding, left = ["A"], right = ["B"]))
        @test !isindexed(kg)
        index!(kg)
        @test isindexed(kg) && contains(sprint(show, MIME"text/plain"(), kg), "(indexed)")
        # direct edit without reindex: removed references disappear at once, new ones after reindex!
        kg["e"].right = ["C"]
        add!(kg, Protein(id = "C"))
        @test isempty(edges_of(kg, "B"))
        @test isempty(edges_of(kg, "C"))
        reindex!(kg)
        @test [e.id for e in edges_of(kg, "C")] == ["e"]
        # merge! keeps the index current
        other = KnowledgeGraph()
        add!(other, Interaction(id = "f", type = :Binding, left = ["C"], right = ["A"]))
        merge!(kg, other)
        @test [e.id for e in edges_of(kg, "C")] == ["e", "f"]
        unindex!(kg)
        @test !isindexed(kg) && [e.id for e in edges_of(kg, "C")] == ["e", "f"]
    end
end

@testset "algorithms" begin
    kg = KnowledgeGraph()
    add!(kg, (Protein(id = x) for x in ["A", "B", "C", "D", "E", "Z"])...,
        Taxon(id = "NCBITaxon:9606"), Provenance(id = "prov"),
        Interaction(id = "e1", type = :Binding, left = ["A"], right = ["B"], provenance = ["prov"]),
        Influence(id = "e2", type = :Activation, left = ["B"], right = ["C"]),
        Influence(id = "e3", type = :Inhibition, left = ["C"], right = ["D"]),
        Influence(id = "e4", type = :Inhibition, left = ["E"], right = ["e1"]),
        PathwayModule(id = "PW", components = ["e1", "e2", "e3", "e4"]))
    kg["A"].organism = "NCBITaxon:9606"

    for indexed in (false, true)
        indexed ? index!(kg) : unindex!(kg)
        @test degree_of(kg, "B") == 2
        @test degree_of(kg, "B"; direction = :out) == 1
        @test degree_of(kg, "C"; type = :Inhibition) == 1
        @test shortest_path(kg, "A", "D") == ["A", "B", "C", "D"]
        @test shortest_path(kg, "D", "A") === nothing
        @test shortest_path(kg, "D", "A"; direction = :both) == ["D", "C", "B", "A"]
        @test shortest_path(kg, "A", "D"; max_hops = 2) === nothing
        @test shortest_path(kg, "A", "C"; type = :Binding) === nothing
        # E inhibits the edge e1 (A binds B): an end point unless edges are reified
        @test shortest_path(kg, "E", "B") === nothing
        @test shortest_path(kg, "E", "e1") == ["E", "e1"]
        @test shortest_path(kg, "E", "B"; reify = true) == ["E", "e4", "e1", "B"]
        @test shortest_path(kg, "B", "E"; reify = true, direction = :in) == ["B", "e1", "e4", "E"]
        @test shortest_path(kg, "A", "D"; reify = true) == ["A", "e1", "B", "e2", "C", "e3", "D"]
        @test shortest_path(kg, "A", "A") == ["A"]
        comps = network_components(kg)
        @test comps == [["A", "B", "C", "D"], ["E", "e1"], ["Z"]]
        @test network_components(kg; reify = true) ==
              [["A", "B", "C", "D", "E", "e1", "e2", "e3", "e4"], ["Z"]]
        @test network_components(kg; type = :Inhibition) == [["C", "D"], ["E", "e1"], ["A"], ["B"], ["Z"]]
    end

    sub = subgraph(kg, "B"; hops = 1)
    @test Set(keys(sub)) == Set(["A", "B", "C", "e1", "e2", "prov", "NCBITaxon:9606"])
    @test isvalid(sub; strict = true)
    @test sub["B"] !== kg["B"]                                     # copies
    sub0 = subgraph(kg, ["PW", "e3"])
    @test sub0["PW"].components == ["e3"]
    @test Set(keys(sub0)) == Set(["PW", "e3", "C", "D"])
    @test Set(keys(subgraph(kg, "A"; hops = 1, include_context = false))) == Set(["A", "B", "e1"])
    @test length(subgraph(kg, "A"; hops = 10)) == 11      # A..E, e1..e4 (e4 targets e1), taxon, provenance
    @test_throws ErrorException (Base.get_extension(BioNExt, :BioNExtGraphsExt) === nothing ? to_graph(kg) : error())
end

@testset "Graphs.jl extension" begin
    if Base.find_package("Graphs") === nothing
        @info "Graphs.jl not installed; skipping extension tests"
    else
        @eval import Graphs
        kg = KnowledgeGraph()
        add!(kg, (Protein(id = x) for x in ["A", "B", "C"])...,
             Interaction(id = "e1", type = :Binding, left = ["A"], right = ["B"]),
             Influence(id = "e2", type = :Inhibition, left = ["C"], right = ["e1"]))
        g = Base.invokelatest(to_graph, kg)
        @test Graphs.nv(g.graph) == 4 && Graphs.ne(g.graph) == 2
        @test g.ids[1:3] == ["A", "B", "C"] && g.vertex["e1"] == 4
        @test Graphs.has_edge(g.graph, g.vertex["A"], g.vertex["B"])
        r = Base.invokelatest(to_graph, kg; reify = true, directed = false)
        @test r.graph isa Graphs.SimpleGraph && Graphs.ne(r.graph) == 4
        t = Base.invokelatest(to_graph, kg; type = :Binding)
        @test Graphs.ne(t.graph) == 1
        @test Graphs.is_connected(r.graph)
    end
end
