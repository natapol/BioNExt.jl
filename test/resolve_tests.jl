@testset "entity resolution" begin
    function sources()
        kg = KnowledgeGraph()
        add!(kg,
            Provenance(id = "prov-uniprot"), Provenance(id = "prov-string"),
            Protein(id = "UniProt:P04637", name = "p53", cellular_location = [:Nucleus], provenance = ["prov-uniprot"]),
            Protein(id = "Ensembl:ENSP00000269305", name = "TP53", external_ids = ["uniprot:P04637"],
                    cellular_location = [:Cytoplasm], sequence = "MEEP"),
            Protein(id = "STRING:9606.ENSP00000269305", external_ids = ["Ensembl:ENSP00000269305"],
                    provenance = ["prov-string"]),
            Protein(id = "UniProt:Q00987", name = "MDM2"),
            Gene(id = "HGNC:11998", external_ids = ["UniProt:P04637"]),          # other class: never merged
            Protein(id = "UniProt:P04637_p.R175H", hgvs = "p.Arg175His", variant_of = ["Ensembl:ENSP00000269305"]),
            Interaction(id = "a1", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"], confidence = 0.9),
            Interaction(id = "s1", type = :Binding, left = ["UniProt:Q00987"], right = ["STRING:9606.ENSP00000269305"],
                        confidence = 0.7, provenance = ["prov-string"]),
            Influence(id = "i1", type = :Inhibition, left = ["UniProt:Q00987"], right = ["a1"]),
            Influence(id = "i2", type = :Inhibition, left = ["UniProt:Q00987"], right = ["s1"]),
            Interaction(id = "a2", type = :Association, left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
            PathwayModule(id = "PW", components = ["a1", "s1", "i1", "i2", "a2"]),
        )
        return kg
    end

    @testset "duplicates and identifiers" begin
        kg = sources()
        @test duplicates(kg) == [["UniProt:P04637", "Ensembl:ENSP00000269305", "STRING:9606.ENSP00000269305"]]
        @test normalize_identifier(kg, "UniProt:P04637") == normalize_identifier(kg, "https://identifiers.org/uniprot:P04637")
        @test normalize_identifier(kg, "UniProtKB:P04637"; prefix_aliases = Dict("UniProtKB" => "uniprot")) ==
              normalize_identifier(kg, "UniProt:P04637")
        @test find_entity(kg, "ensembl:ENSP00000269305").id == "Ensembl:ENSP00000269305"
        @test find_entity(kg, "nothing:here") === nothing
        @test identifier_map(kg)[normalize_identifier(kg, "Ensembl:ENSP00000269305")] == "Ensembl:ENSP00000269305"

        xr = KnowledgeGraph()
        add!(xr, Protein(id = "A", db_xrefs = ["PDB:1TUP"]), Protein(id = "B", db_xrefs = ["PDB:1TUP"]))
        @test isempty(duplicates(xr))                      # xrefs do not assert equivalence...
        @test duplicates(xr; use_xrefs = true) == [["A", "B"]]  # ...unless asked
    end

    @testset "resolve!" begin
        for indexed in (false, true)
            kg = sources()
            indexed && index!(kg)
            res = resolve!(kg; prefer = ["UniProt"], on_conflict = (f, a, b) -> f === :confidence ? max(a, b) : a)
            @test res.groups == [["UniProt:P04637", "Ensembl:ENSP00000269305", "STRING:9606.ENSP00000269305"],
                                 ["a1", "s1"], ["i1", "i2"]]
            @test res.mapping["STRING:9606.ENSP00000269305"] == "UniProt:P04637"
            @test res.mapping["s1"] == "a1"
            @test contains(sprint(show, res), "4 resources merged into 3")

            p = kg["UniProt:P04637"]
            @test p.name == "p53"                                      # canonical kept its value
            @test p.cellular_location == [:Nucleus, :Cytoplasm]         # vectors united
            @test p.sequence == "MEEP"                                  # empty scalar filled in
            @test p.provenance == ["prov-uniprot", "prov-string"]
            @test p.external_ids == ["Ensembl:ENSP00000269305", "STRING:9606.ENSP00000269305"]
            @test any(c -> c.field === :name && c.kept == "p53" && c.discarded == "TP53", res.conflicts)
            @test kg["a1"].confidence == 0.9                            # custom on_conflict
            @test kg["a1"].provenance == ["prov-string"]

            @test kg["UniProt:P04637_p.R175H"].variant_of == ["UniProt:P04637"]   # references rewritten
            @test kg["i1"].right == ["a1"]
            @test kg["PW"].components == ["a1", "i1", "a2"]
            @test haskey(kg, "HGNC:11998") && haskey(kg, "a2")          # other class / other type untouched
            @test !haskey(kg, "s1") && !haskey(kg, "Ensembl:ENSP00000269305")
            @test find_entity(kg, "STRING:9606.ENSP00000269305").id == "UniProt:P04637"
            @test [e.id for e in edges_of(kg, "UniProt:P04637")] == ["a1", "a2"]
            @test [v.id for v in variants(kg, "UniProt:P04637")] == ["UniProt:P04637_p.R175H"]
            @test isvalid(kg; strict = true)
            @test isempty(duplicates(kg))
        end

        kg = sources()
        res = resolve!(kg; merge_edges = false)
        @test res.groups == [["UniProt:P04637", "Ensembl:ENSP00000269305", "STRING:9606.ENSP00000269305"]]
        @test haskey(kg, "s1") && kg["s1"].right == ["UniProt:P04637"]

        kg = sources()
        resolve!(kg; prefer = ["Ensembl"])
        @test haskey(kg, "Ensembl:ENSP00000269305") && !haskey(kg, "UniProt:P04637")
        @test kg["Ensembl:ENSP00000269305"].name == "TP53"
        kg = sources()
        resolve!(kg; prefer = g -> last(g))
        @test haskey(kg, "STRING:9606.ENSP00000269305")

        # nothing to do
        kg = KnowledgeGraph(); add!(kg, Gene(id = "g"))
        @test isempty(resolve!(kg).groups)
    end

    @testset "integrate!" begin
        base = KnowledgeGraph()
        add!(base, Protein(id = "UniProt:P04637", name = "p53"), Protein(id = "UniProt:Q00987", name = "MDM2"),
                   Interaction(id = "lit1", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"],
                               publications = ["PMID:1535557"]))
        index!(base)
        string_db = KnowledgeGraph(prefixes = Dict("string" => "https://string-db.org/network/"))
        add!(string_db,
            Protein(id = "string:9606.ENSP00000269305", external_ids = ["UniProt:P04637"]),
            Protein(id = "string:9606.ENSP00000258149", external_ids = ["UniProt:Q00987"]),
            Protein(id = "UniProt:Q00987", synonyms = ["HDM2"]),                          # same id: merged field-wise
            Interaction(id = "str1", type = :Binding, left = ["string:9606.ENSP00000258149"],
                        right = ["string:9606.ENSP00000269305"], confidence = 0.999))
        res = integrate!(base, string_db; prefer = ["UniProt"])
        @test sort(collect(keys(base))) == ["UniProt:P04637", "UniProt:Q00987", "lit1"]
        @test base["lit1"].confidence == 0.999 && base["lit1"].publications == ["PMID:1535557"]
        @test base["UniProt:Q00987"].synonyms == ["HDM2"]
        @test base.prefixes["string"] == "https://string-db.org/network/"
        @test res.mapping["str1"] == "lit1"
        @test [e.id for e in edges_of(base, "UniProt:P04637")] == ["lit1"]   # index still current
        @test length(string_db) == 4 && string_db["UniProt:Q00987"].synonyms == ["HDM2"]  # source untouched

        clash = KnowledgeGraph(); add!(clash, Gene(id = "UniProt:P04637"))
        @test_throws ArgumentError integrate!(base, clash)
    end
end
