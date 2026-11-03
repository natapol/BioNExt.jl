same_graph(a, b) = keys(a) == keys(b) &&
    all(getfield(b[r.id], f) == getfield(r, f) for r in a for f in fieldnames(typeof(r)))

@testset "I/O formats" begin
    kg = sample_graph()
    dir = mktempdir()

    @testset "round trip via $ext" for ext in (".ttl", ".nt", ".jsonld", ".json", ".tsv", ".csv")
        path = joinpath(dir, "graph" * ext)
        save(path, kg)
        @test same_graph(kg, load(path))
    end

    @testset "write-only formats" begin
        for ext in (".rdf", ".graphml", ".dot", ".sif")
            save(joinpath(dir, "g" * ext), kg)
        end
        xml = read(joinpath(dir, "g.rdf"), String)
        @test contains(xml, "<rdf:Description rdf:about=\"https://identifiers.org/UniProt:P04637\">")
        @test contains(xml, "<bn:hasHGVS>p.Arg273His</bn:hasHGVS>")
        @test contains(xml, "<bn:hasConfidenceScore rdf:datatype=\"http://www.w3.org/2001/XMLSchema#float\">0.8</bn:hasConfidenceScore>")
        @test contains(xml, "TP53 &quot;tumor&quot; protein")
        gml = read(joinpath(dir, "g.graphml"), String)
        @test contains(gml, "<edge source=\"UniProt:Q00987\" target=\"UniProt:P04637\"><data key=\"type\">Binding</data>")
        dot = read(joinpath(dir, "g.dot"), String)
        @test contains(dot, "\"cx1\" -> \"r1\" [label=\"Inhibition\", arrowhead=tee];")
        @test_throws ArgumentError load(joinpath(dir, "g.graphml"))
        @test_throws ArgumentError save(joinpath(dir, "g.xyz"), kg)
    end

    @testset "network flattening" begin
        flat = network(kg)
        @test Arc("UniProt:Q00987", "UniProt:P04637", :Binding, "e1", :pair) in flat.arcs
        @test "r1" in flat.nodes                      # e3 targets an edge, kept as a node
        reified = network(kg; reify = true)
        @test Arc("cx1", "e3", :Inhibition, "e3", :left) in reified.arcs
        @test Arc("e3", "r1", :Inhibition, "e3", :right) in reified.arcs
    end

    @testset "SIF import" begin
        sif = """
        TP53\tactivation\tMDM2\tCDKN1A
        MDM2\tpp\tTP53
        ATM
        """
        g = read_sif(IOBuffer(sif); relationships = Dict("pp" => :Binding))
        @test Set(keys(g)) ⊇ Set(["TP53", "MDM2", "CDKN1A", "ATM"])
        @test g["ATM"] isa Protein
        @test [e.type for e in edges(g)] == [:Activation, :Activation, :Binding]
        @test isvalid(g)
        g2 = read_sif(IOBuffer("A\tsomething_odd\tB\n"); node_type = Gene, edge_type = Influence)
        @test g2["A"] isa Gene && edges(g2)[1] isa Influence
        @test any(i -> contains(i.message, "unknown vocabulary term :something_odd"), validate(g2))
    end

    @testset "table import with mapping" begin
        genes = "symbol,hgnc,aliases,taxon\nTP53,HGNC:11998,p53|LFS1,NCBITaxon:9606\nMDM2,HGNC:6973,,NCBITaxon:9606\n"
        g = read_table(IOBuffer(genes); delim = ',', type = Gene,
                       mapping = Dict("symbol" => :name, "hgnc" => :id, "aliases" => :synonyms, "taxon" => :organism))
        @test g["HGNC:11998"].synonyms == ["p53", "LFS1"]
        @test g["HGNC:6973"].organism == "NCBITaxon:9606"

        # an edge list without ids, read into the same graph; a quoted cell with the delimiter
        ppi = "source\ttarget\tinteraction\tscore\nHGNC:11998\tHGNC:6973\tBinding\t0.7\n\"HGNC:6973\"\tHGNC:11998\tInhibition\t\n"
        read_table(IOBuffer(ppi); type = Influence, into = g,
                   mapping = Dict("source" => :left, "target" => :right, "interaction" => :type, "score" => :confidence))
        @test [e.id for e in edges(g)] == ["influence-1", "influence-2"]
        @test g["influence-1"].confidence == 0.7 && g["influence-2"].type == :Inhibition

        # updating an existing row by id
        read_table(IOBuffer("class\tid\tname\nGene\tHGNC:6973\tMDM2 proto-oncogene\n"); into = g)
        @test g["HGNC:6973"].name == "MDM2 proto-oncogene" && g["HGNC:6973"].synonyms == String[]

        @test_throws ArgumentError read_table(IOBuffer("name\nx\n"))
        @test_throws ArgumentError read_table(IOBuffer("class\tid\tnet_charge\nSmallMolecule\tm1\tnot-a-number\n"))

        # partial export: just the edges
        tsv = sprint(io -> write_table(io, kg; types = EdgeEntity))
        @test startswith(tsv, "class\tid\t")
        @test count(==('\n'), tsv) == 1 + length(edges(kg))
    end

    @testset "load! merges into an existing graph" begin
        g = KnowledgeGraph()
        add!(g, Protein(id = "UniProt:P04637", name = "kept"))
        load!(g, joinpath(dir, "graph.ttl"))
        @test g["UniProt:P04637"].name == "kept" && length(g) == length(kg)
        load!(g, joinpath(dir, "graph.ttl"); replace = true)
        @test g["UniProt:P04637"].name == "p53"
    end

    @testset "Turtle syntax coverage" begin
        ttl = """
        PREFIX bn: <https://w3id.org/bionext/ontology/>
        @prefix up: <http://purl.uniprot.org/uniprot/> .
        @base <http://example.org/data/> .
        # comment line
        up:P04637 a bn:Protein ;
            bn:name '''p53
        tumor suppressor''' ;
            bn:synonyms "TP53"@en, "p53\\ttab" ;
            bn:cellularLocation bn:Nucleus ; .
        <rel1> a bn:Interaction ; bn:hasType bn:Binding ;
            bn:hasLeft up:Q00987 ; bn:hasRight up:P04637 ; bn:hasWeight 0.5 ; bn:hasConfidenceScore 1e-1 .
        up:Q00987 a bn:Protein ; bn:name "MDM2" .
        [] bn:note ( 1 2.5 true [ bn:x "y" ] ) .
        _:b1 bn:name "blank" .
        """
        ts, pfx = parse_turtle(IOBuffer(ttl))
        @test pfx["up"] == "http://purl.uniprot.org/uniprot/"
        g = read_turtle(IOBuffer(ttl))
        @test g["up:P04637"].name == "p53\ntumor suppressor"
        @test g["up:P04637"].synonyms == ["TP53", "p53\ttab"]
        @test g["http://example.org/data/rel1"].weight == 0.5
        @test g["http://example.org/data/rel1"].confidence == 0.1
        @test count(t -> endswith(t.predicate.value, "#first"), ts) == 4
        @test any(t -> t.object isa Literal && t.object.datatype == BioNExt.XSD_BOOLEAN, ts)
        @test_throws ArgumentError parse_turtle(IOBuffer("nope:x a bn:Gene ."))
        @test_throws ArgumentError parse_turtle(IOBuffer("<a> <b> \"unterminated ."))
    end

    @testset "Turtle parser on the BioNExt ontology" begin
        path = joinpath(@__DIR__, "..", "ontology", "bne_ontology.ttl")
        if isfile(path)
            ts, pfx = parse_turtle(path)
            # the generated vocabulary must be in sync with the ontology file
            @test pfx["bn"] == BioNExt.BN && pfx["bnr"] == BioNExt.BNR
            onto = BioNExt.IRI(BioNExt.BN)
            @test any(t -> t.subject == onto && endswith(t.predicate.value, "versionInfo") &&
                           t.object.value == BioNExt.ONTOLOGY_VERSION, ts)
            owl(x) = "http://www.w3.org/2002/07/owl#" * x
            ntype(o) = count(t -> t.predicate.value == BioNExt.RDF_TYPE && t.object == IRI(owl(o)), ts)
            @test ntype("Class") == 141 && ntype("ObjectProperty") == 36 && ntype("DatatypeProperty") == 20
        end
    end

    @testset "legacy namespace is upgraded on read" begin
        old = """
        @prefix bn: <https://jamiesilakom618.github.io/BioNExt/ontology/> .
        @prefix bnr: <https://jamiesilakom618.github.io/BioNExt/resource/> .
        bnr:e1 a bn:Interaction ; bn:hasType bn:Binding ; bn:hasLeft bnr:a ; bn:hasRight bnr:b .
        bnr:a a bn:Protein ; bn:name "A" .
        bnr:b a bn:Protein .
        """
        g = read_turtle(IOBuffer(old))
        @test g.base == BioNExt.BNR
        @test g["e1"].type == :Binding && g["e1"].left == ["a"] && g["a"].name == "A"
        nt = "<https://jamiesilakom618.github.io/BioNExt/resource/x> " *
             "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://jamiesilakom618.github.io/BioNExt/ontology/Gene> .\n"
        @test read_ntriples(IOBuffer(nt))["x"] isa Gene
        @test !contains(sprint(write_turtle, g), "jamiesilakom618")
    end

    @testset "JSON-LD from other tools" begin
        doc = """
        {"@context": {"bn": "https://w3id.org/bionext/ontology/",
                      "name": "bn:name",
                      "left": {"@id": "bn:hasLeft", "@type": "@id"},
                      "right": {"@id": "bn:hasRight", "@type": "@id"},
                      "kind": {"@id": "bn:hasType", "@type": "@id"}},
         "@graph": [
           {"@id": "http://ex.org/a", "@type": "bn:Gene", "name": "A \\u00e9 \\ud83e\\uddec"},
           {"@id": "http://ex.org/b", "@type": ["bn:Gene"], "name": {"@value": "B", "@language": "en"}},
           {"@id": "http://ex.org/e", "@type": "bn:Influence", "kind": "bn:Activation",
            "left": "http://ex.org/a", "right": ["http://ex.org/b"], "bn:hasWeight": 2}
         ]}
        """
        g = read_jsonld(IOBuffer(doc))
        @test g["http://ex.org/a"].name == "A é 🧬"
        @test g["http://ex.org/e"].type == :Activation
        @test g["http://ex.org/e"].right == ["http://ex.org/b"]
        @test g["http://ex.org/e"].weight == 2.0
        @test isvalid(g)
    end

    @testset "JSON" begin
        @test parse_json("""{"a": [1, 2.5, -3e2, true, null, "x\\"y"]}""") ==
              Dict("a" => Any[1, 2.5, -300.0, true, nothing, "x\"y"])
        @test_throws ArgumentError parse_json("{\"a\": }")
        @test_throws ArgumentError read_json(IOBuffer("""{"resources": [{"class": "Nope", "id": "x"}]}"""))
    end
end
