@testset "audit" begin
    # the shared sample graph has a real inconsistency: Ser15 on a 10-residue sequence
    sample = audit(sample_graph())
    @test count_findings(sample; severity = :error) == 1
    @test only(findings(sample; severity = :error)).message == "modification site Ser15 is beyond the sequence (10 residues)"
    @test sample.checks == first.(AUDIT_CHECKS)
    good = sample_graph()
    good["UniProt:P04637"].sequence = nothing
    r = audit(good)
    @test passed(r) && count_findings(r; severity = :error) == 0
    @test contains(sprint(show, MIME"text/plain"(), r), "Audit of 21 resources: passed")

    kg = KnowledgeGraph()
    add!(kg,
        Taxon(id = "NCBITaxon:9606"), Taxon(id = "NCBITaxon:10090"), Taxon(id = "NCBITaxon:7227"),
        Provenance(id = "prov-unused"),
        Protein(id = "UniProt:P04637", name = "p53", organism = "NCBITaxon:9606", sequence = "MEEPQSDPSVEPPLSQ",
                modification_sites = ["Ser7", "Ser9", "Thr99", "phospho"], provenance = ["prov-unused"]),
        Protein(id = "UniProt:NOT_AN_ACC", name = "bad id", organism = "NCBITaxon:10090", sequence = "MEEP1Q"),
        Protein(id = "Ensembl:ENSP00000269305", name = "p53 again", external_ids = ["UniProt:P04637"]),
        Protein(id = "lonely"),
        Dna(id = "d1", sequence = "ACGTU", length = 10, start_pos = 50, end_pos = 40),
        Dna(id = "d2", start_pos = 1, end_pos = 10, length = 12),
        SmallMolecule(id = "CHEBI:15422", name = "ATP", inchikey = "not-a-key", cas_number = "56-65-5", formula = "C10H16N5O13P3"),
        SmallMolecule(id = "CHEBI:x", name = "bad", cas_number = "56-65-4", formula = "hello world"),
        Protein(id = "UniProt:P04637_v", name = "v", hgvs = "R175H", variant_of = ["UniProt:P04637"]),
        MolecularComplex(id = "C-empty"), MolecularComplex(id = "C-loop", subcomplexes = ["C-loop2"]),
        MolecularComplex(id = "C-loop2", subcomplexes = ["C-loop"]),
        Interaction(id = "e1", type = :Binding, left = ["UniProt:P04637"], right = ["UniProt:NOT_AN_ACC"],
                    confidence = 1.5, publications = ["PMID:abc"]),
        Interaction(id = "e2", type = :Binding, left = ["UniProt:P04637"], right = ["UniProt:NOT_AN_ACC"], weight = NaN),
        Interaction(id = "e3", type = :Binding, left = ["UniProt:P04637"], right = ["UniProt:P04637"]),
        Influence(id = "up", type = :Activation, left = ["UniProt:P04637"], right = ["d1"]),
        Influence(id = "down", type = :Inhibition, left = ["UniProt:P04637"], right = ["d1"]),
        Influence(id = "badtype", type = :RegulateExpression, left = ["UniProt:P04637"], right = ["lonely2"]),
        Protein(id = "lonely2", name = "x"),
        PathwayModule(id = "PW-empty"), PathwayModule(id = "PW-a", components = ["PW-b"]),
        PathwayModule(id = "PW-b", components = ["PW-a"]),
    )
    r = audit(kg)
    has(check, sev, id, text) = any(f -> f.check === check && f.severity === sev && f.id == id && occursin(text, f.message), r.findings)
    @test !passed(r) && !passed(r; fail_on = :warning)
    @test has(:ontology, :error, "badtype", ":RegulateExpression expects Gene")
    @test has(:identifiers, :warning, "UniProt:NOT_AN_ACC", "not a well-formed UniProt")
    @test !has(:identifiers, :warning, "UniProt:P04637_v", "")                  # variant suffix allowed
    @test has(:identifiers, :warning, "CHEBI:x", "id")
    @test has(:identifiers, :warning, "e1", "PMID:abc")
    @test has(:sequences, :error, "UniProt:NOT_AN_ACC", "'1'")
    @test has(:sequences, :error, "d1", "'U'") && has(:sequences, :warning, "d1", "length is 10 but the sequence has 5")
    @test has(:coordinates, :error, "d1", "end position 40 is before start position 50")
    @test has(:coordinates, :warning, "d2", "length 12 does not match the coordinates 1–10 (10)")
    @test has(:modification_sites, :warning, "UniProt:P04637", "Ser7: the sequence has D")
    @test has(:modification_sites, :error, "UniProt:P04637", "Thr99 is beyond the sequence")
    @test has(:modification_sites, :info, "UniProt:P04637", "\"phospho\"")
    @test !has(:modification_sites, :warning, "UniProt:P04637", "Ser9")          # sequence[9] is S
    @test has(:values, :warning, "e1", "confidence 1.5") && has(:values, :error, "e2", "weight is NaN")
    @test has(:values, :warning, "CHEBI:15422", "InChIKey") && !has(:values, :warning, "CHEBI:15422", "CAS")
    @test has(:values, :warning, "CHEBI:x", "CAS number") && has(:values, :warning, "CHEBI:x", "formula")
    @test has(:values, :warning, "UniProt:P04637_v", "HGVS")
    @test has(:duplicates, :warning, "UniProt:P04637", "Ensembl:ENSP00000269305")
    @test has(:duplicates, :warning, "e1", "e2")
    @test has(:edges, :info, "e3", "self-loop")
    @test has(:edges, :info, "e1", "different organisms")
    @test has(:edges, :info, "up", "inhibiting edge(s) down")
    @test has(:structure, :warning, "C-empty", "no subunits") && has(:structure, :warning, "PW-empty", "no components")
    @test has(:structure, :warning, "C-loop", "contains itself") && has(:structure, :warning, "PW-a", "contains itself")
    @test has(:connectivity, :info, "lonely", "not connected") && has(:connectivity, :info, "NCBITaxon:7227", "not used")
    @test has(:annotation, :info, "lonely", "no name") && has(:annotation, :info, "e2", "no publication")

    @test all(f -> f.check === :values, findings(r; check = :values))
    @test count_findings(r; check = :sequences, severity = :error) == 2
    small = audit(kg; checks = [:sequences, :values], max_per_check = 1)
    @test small.checks == [:sequences, :values] && length(small.findings) == 2
    @test count_findings(small; check = :values) > 1                         # counted beyond the cap
    @test audit(kg; skip = [:ontology, :annotation]).checks == setdiff(first.(AUDIT_CHECKS), [:ontology, :annotation])
    @test_throws ArgumentError audit(kg; checks = [:nonsense])
    @test_throws QualityError assert_quality(kg)
    @test assert_quality(good) isa AuditReport
    @test_throws QualityError assert_quality(good; fail_on = :info)
    @test contains(sprint(showerror, QualityError(r, :error)), "failed its audit")
    @test audit(kg; checks = [:identifiers], patterns = Dict("UniProt" => r".*")).findings |>
          fs -> !any(f -> f.id == "UniProt:NOT_AN_ACC", fs)
end

@testset "graph_stats" begin
    kg = sample_graph()
    s = graph_stats(kg; top = 3)
    @test s.resources == 21 && s.nodes == 12 && s.edges == 4
    @test first(s.classes) == ("Protein" => 4)
    @test Dict(s.relationship_types)[:Binding] == 1
    @test s.hubs[1] == ("UniProt:P04637" => 2) && length(s.hubs) == 3
    @test s.degree.max == 2 && s.isolated == count(n -> isempty(edges_of(kg, n)), nodes(kg))
    @test sum(last, s.degree_distribution) == s.nodes
    @test s.pathways.count == 2 && s.pathway_depth == 2
    @test s.complexes == (count = 1, mean = 2.0, max = 2)
    @test s.components == length(network_components(kg))
    @test 0 <= s.density <= 1 && 0 <= s.reciprocity <= 1 && s.clustering !== nothing
    @test graph_stats(kg; clustering = false).clustering === nothing
    @test ("NCBITaxon:9606 (Homo sapiens)" => 1) in s.organisms
    @test Dict(s.id_prefixes)["UniProt"] == 4
    @test Dict(s.completeness)["provenance (entities)"] ≈ 1 / count(r -> r isa Entity, kg)

    # sources and overlap
    g = KnowledgeGraph()
    add!(g, Provenance(id = "A"), Provenance(id = "B"),
         Protein(id = "p", name = "p", provenance = ["A", "B"]), Protein(id = "q", provenance = ["A"]),
         Protein(id = "r", provenance = ["B"]),
         Interaction(id = "x", type = :Binding, left = ["p"], right = ["q"], publications = ["PMID:1"]),
         Interaction(id = "y", type = :Binding, left = ["q"], right = ["r"]),
         Interaction(id = "z", type = :Binding, left = ["r"], right = ["p"], confidence = 0.5))
    t = graph_stats(g)
    @test t.sources == ["A" => 2, "B" => 2] && t.multi_source == 1 && t.source_overlap == [("A", "B") => 1]
    @test t.clustering ≈ 1.0 && t.components == 1 && t.reciprocity == 0.0
    @test Dict(t.completeness)["publications (edges)"] ≈ 1 / 3
    @test t.publications == 1

    rows = stats_table(t)
    @test (section = "overview", metric = "edges", value = 3) in rows
    for (ext, needle) in ((".md", "| edges | 3 |"), (".tsv", "overview\tedges\t3"), (".csv", "overview,edges,3"), (".json", "\"edges\": 3"))
        path = joinpath(mktempdir(), "stats" * ext)
        write_stats(path, t)
        @test contains(read(path, String), needle)
    end
    @test BioNExt.parse_json(sprint(io -> write_stats(io, t; format = :json)))["overview"]["nodes"] == 3
    @test_throws ArgumentError write_stats(IOBuffer(), t; format = :xml)
    @test contains(sprint(show, MIME"text/plain"(), t), "entities from ≥ 2 sources")
    @test graph_stats(KnowledgeGraph()).nodes == 0
end
