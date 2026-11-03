using BioNExt, Dates, Test

function sample_graph()
    kg = KnowledgeGraph()
    add!(kg,
        Taxon(id = "NCBITaxon:9606", name = "Homo sapiens"),
        Taxon(id = "NCBITaxon:4932", name = "Saccharomyces cerevisiae"),
        Provenance(id = "prov1", collected_date = DateTime(2026, 10, 1, 12, 30),
                   sources = ["https://www.uniprot.org"]),
        Dna(id = "chr17-TP53", chromosome = "17", start_pos = 7661779, end_pos = 7687538,
            regions = [:Promoter, :Exon]),
        Gene(id = "HGNC:11998", name = "TP53", templates = ["chr17-TP53"], organism = "NCBITaxon:9606"),
        Gene(id = "HGNC:6973", name = "MDM2"),
        Rna(id = "ENST00000269305", encoded_by = ["chr17-TP53"], regions = [:MessengerRNA]),
        Protein(id = "UniProt:P04637", name = "p53", sequence = "MEEPQSDPSV",
                cellular_location = [:Nucleus], provenance = ["prov1"],
                modification_sites = ["Ser15"], synonyms = ["TP53 \"tumor\" protein"]),
        Protein(id = "UniProt:P04637_p.R273H", hgvs = "p.Arg273His", variant_of = ["UniProt:P04637"]),
        Protein(id = "UniProt:Q00987", name = "MDM2"),
        Protein(id = "UniProt:P38123", name = "yeast homolog", organism = "NCBITaxon:4932"),
        SmallMolecule(id = "CHEBI:15422", name = "ATP", net_charge = -4, smiles = "Nc1ncnc2..."),
        SmallMolecule(id = "CHEBI:16761", name = "ADP"),
        MolecularComplex(id = "cx1", subunits = ["UniProt:P04637", "UniProt:Q00987"]),
        OrthologGroup(id = "OG1", members = ["UniProt:P04637", "UniProt:P38123", "HGNC:11998"],
                      species = ["NCBITaxon:9606", "NCBITaxon:4932"], clustering_method = "OrthoDB",
                      confidence = 0.8),
        Interaction(id = "e1", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
        Influence(id = "e2", type = :RegulateExpression, left = ["UniProt:P04637"], right = ["HGNC:6973"]),
        Reaction(id = "r1", type = :Irreversible, left = ["CHEBI:15422"], right = ["CHEBI:16761"],
                 ec_numbers = ["EC:2.7.11.1"]),
        Influence(id = "e3", type = :Inhibition, left = ["cx1"], right = ["r1"]),
        PathwayModule(id = "sub", components = ["r1"]),
        PathwayModule(id = "pw", name = "p53 signalling", components = ["e1", "e2", "e3", "sub"]),
    )
    return kg
end

@testset "BioNExt" begin
    @testset "vocabulary" begin
        @test is_a(:NegativeCorrelation, :RelationshipType)
        @test is_a(:ThreePrimeUTR, :DnaRegion)
        @test !is_a(:Nucleus, :RelationshipType)
        @test :Binding in terms(:RelationshipType)
        @test term(:Activation).expected_right == [:EdgeEntity, :NodeEntity, :PathwayModule]
        @test_throws KeyError term(:NotATerm)
    end

    @testset "graph basics and queries" begin
        kg = sample_graph()
        @test length(kg) == 21
        @test kg["UniProt:P04637"] isa Protein
        @test_throws ArgumentError add!(kg, Gene(id = "HGNC:6973"))
        @test length(resources(kg, EdgeEntity)) == 4
        @test Set(e.id for e in edges_of(kg, "UniProt:P04637")) == Set(["e1", "e2"])
        @test [e.id for e in edges_of(kg, "UniProt:P04637"; side = :left)] == ["e2"]
        @test neighbors(kg, "UniProt:Q00987"; direction = :out) == ["UniProt:P04637"]
        @test isempty(edges_of(kg, "UniProt:P04637"; type = :CoExpression))
        @test [c.id for c in complexes_of(kg, "UniProt:Q00987")] == ["cx1"]
        @test [v.id for v in variants(kg, "UniProt:P04637")] == ["UniProt:P04637_p.R273H"]
        @test [p.id for p in pathways_of(kg, "r1")] == ["sub"]
        @test Set(components(kg, "pw"; recursive = true)) == Set(["e1", "e2", "e3", "sub", "r1"])
        @test ("cx1" => :subunits) in references(kg, "UniProt:P04637")
        @test [g.id for g in ortholog_groups_of(kg, "UniProt:P38123")] == ["OG1"]
        delete!(kg, "sub")
        @test !haskey(kg, "sub") && length(kg) == 20
    end

    @testset "validation" begin
        kg = sample_graph()
        @test isempty(validate(kg))
        @test isvalid(kg)

        bad = sample_graph()
        add!(bad,
            Interaction(id = "x1", type = :Binding, left = ["UniProt:P04637"], right = ["pw"]),  # Interaction right must be NodeEntity
            Interaction(id = "x2", type = :Nucleus, left = ["UniProt:P04637"], right = ["UniProt:Q00987"]),
            Influence(id = "x3", type = :RegulateExpression, left = ["UniProt:P04637"], right = ["UniProt:Q00987"]),
            Influence(id = "x4", type = :Activation, left = String[], right = ["missing-id"]),
            Protein(id = "v-nohgvs", variant_of = ["HGNC:11998"]),   # variant of wrong class, no HGVS
            Rna(id = "rna-bad", regions = [:Promoter], start_pos = 0),
        )
        msgs = Dict{String,Vector{String}}()
        for i in validate(bad)
            push!(get!(msgs, i.id, String[]), i.message)
        end
        @test any(contains("expected NodeEntity"), msgs["x1"])
        @test any(contains("is not a RelationshipType"), msgs["x2"])
        @test any(contains(":RegulateExpression expects Gene"), msgs["x3"])
        @test any(contains("no left participant"), msgs["x4"])
        @test any(i -> i.id == "x4" && i.severity === :warning, validate(bad))
        @test any(contains("HGVS"), msgs["v-nohgvs"])
        @test any(contains("expected Protein"), msgs["v-nohgvs"])
        @test any(contains("not a RnaRegion"), msgs["rna-bad"])
        @test any(contains("positive integer"), msgs["rna-bad"])
        @test !isvalid(bad)
    end

    @testset "RDF round trip" begin
        kg = sample_graph()
        ttl = sprint(write_turtle, kg)
        @test contains(ttl, "bnr:e1 a bn:Interaction")
        @test contains(ttl, "bn:hasProteinMember")
        @test contains(ttl, "bn:hasVariant")
        @test contains(ttl, "\\\"tumor\\\"")

        for inverses in (true, false)
            nt = sprint(io -> write_ntriples(io, kg; inverses))
            kg2 = read_ntriples(IOBuffer(nt))
            @test keys(kg2) == keys(kg)
            for r in kg, f in fieldnames(typeof(r))
                @test getfield(kg2[r.id], f) == getfield(r, f)
            end
        end

        # inverse-only input is folded into the canonical field
        nt = """
        <$(BioNExt.BNR)g1> <$(BioNExt.RDF_TYPE)> <$(BioNExt.BN)Gene> .
        <$(BioNExt.BNR)d1> <$(BioNExt.RDF_TYPE)> <$(BioNExt.BN)Dna> .
        <$(BioNExt.BNR)d1> <$(BioNExt.BN)annotates> <$(BioNExt.BNR)g1> .
        <$(BioNExt.BNR)d1> <$(BioNExt.BN)name> "DNA \\u00e9\\"x\\""@en .
        """
        kg3 = read_ntriples(IOBuffer(nt))
        @test kg3["g1"].templates == ["d1"]
        @test kg3["d1"].name == "DNA é\"x\""
    end
end

include("io_tests.jl")
include("graph_tests.jl")
include("resolve_tests.jl")
include("quality_tests.jl")
include("docs_tests.jl")
