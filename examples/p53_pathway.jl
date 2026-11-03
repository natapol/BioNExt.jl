# A small human p53–MDM2 module built with BioNExt.jl.
#   julia --project=. examples/p53_pathway.jl

using BioNExt, Dates

kg = KnowledgeGraph()

add!(kg,
    Taxon(id = "NCBITaxon:9606", name = "Homo sapiens"),
    Provenance(id = "prov-uniprot-2026", sources = ["https://www.uniprot.org"],
               collected_date = DateTime(2026, 10, 1)),

    Gene(id = "HGNC:11998", name = "TP53", organism = "NCBITaxon:9606",
         functions = ["GO:0003700"], synonyms = ["p53", "LFS1"]),
    Gene(id = "HGNC:6973", name = "MDM2", organism = "NCBITaxon:9606"),
    Protein(id = "UniProt:P04637", name = "Cellular tumor antigen p53",
            organism = "NCBITaxon:9606", cellular_location = [:Nucleus, :Cytoplasm],
            provenance = ["prov-uniprot-2026"], publications = ["PMID:2046748"]),
    Protein(id = "UniProt:P04637_p.R273H", name = "p53 R273H", hgvs = "p.Arg273His",
            variant_of = ["UniProt:P04637"], variant_type = ["SO:1000002"]),
    Protein(id = "UniProt:Q00987", name = "E3 ubiquitin-protein ligase Mdm2",
            organism = "NCBITaxon:9606", cellular_location = [:Nucleus]),
    SmallMolecule(id = "CHEBI:29105", name = "zinc(2+)", formula = "Zn", net_charge = 2),

    Interaction(id = "e-p53-mdm2-binding", type = :Binding,
                left = ["UniProt:Q00987"], right = ["UniProt:P04637"], confidence = 0.98),
    Influence(id = "e-mdm2-ubq-p53", type = :Ubiquitination,
              left = ["UniProt:Q00987"], right = ["UniProt:P04637"], weight = 0.9),
    Influence(id = "e-p53-expr-mdm2", type = :RegulateExpression,
              left = ["UniProt:P04637"], right = ["HGNC:6973"]),

    PathwayModule(id = "p53-mdm2-loop", name = "p53–MDM2 negative feedback loop",
                  components = ["e-p53-mdm2-binding", "e-mdm2-ubq-p53", "e-p53-expr-mdm2"]),
)
kg["UniProt:P04637"].cofactors = ["CHEBI:29105"]

display(kg); println()

println("\nIssues:")
foreach(println, validate(kg))

println("\nWhat does MDM2 act on? ", neighbors(kg, "UniProt:Q00987"; direction = :out))
println("p53 edges of type Binding: ", [e.id for e in edges_of(kg, "UniProt:P04637"; type = :Binding)])
println("Variants of p53: ", [v.id for v in variants(kg, "UniProt:P04637")])

out = joinpath(@__DIR__, "output")
mkpath(out)
for ext in (".ttl", ".nt", ".jsonld", ".rdf", ".json", ".tsv", ".sif", ".graphml", ".dot")
    save(joinpath(out, "p53_pathway" * ext), kg)
end
save(joinpath(out, "p53_edges.csv"), kg; types = EdgeEntity)   # just the edge list
println("\nWrote examples/output/p53_pathway.{ttl,nt,jsonld,rdf,json,tsv,sif,graphml,dot}")

# Any readable format loads back to the same graph.
@assert keys(load(joinpath(out, "p53_pathway.jsonld"))) == keys(kg)
