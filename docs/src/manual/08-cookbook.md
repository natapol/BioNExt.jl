# 8. Cookbook

Complete recipes for common jobs. Each one writes its own small input files, so you can copy
it as-is and then point it at your real data.

```julia
using BioNExt, Dates
```

## Importing a scored protein–protein interaction table

A STRING-style table with scores from 0 to 1000. This recipe keeps confident pairs only,
scales the scores to 0–1, and records where the data came from.

```julia
write("string_links.tsv", """
protein1\tprotein2\tcombined_score
9606.ENSP00000269305\t9606.ENSP00000258149\t999
9606.ENSP00000269305\t9606.ENSP00000278616\t960
9606.ENSP00000258149\t9606.ENSP00000340989\t310
""")

kg = KnowledgeGraph(prefixes = Dict("string" => "https://string-db.org/network/"))
add!(kg,
    Taxon(id = "NCBITaxon:9606", name = "Homo sapiens"),
    Provenance(id = "prov-string-v12", sources = ["https://string-db.org"],
               collected_date = DateTime(2026, 10, 30)))

for (n, line) in enumerate(Iterators.drop(eachline("string_links.tsv"), 1))
    a, b, score = split(line, '\t')
    s = parse(Int, score) / 1000
    s < 0.7 && continue                                   # high-confidence links only
    ids = ["string:" * a, "string:" * b]
    for id in ids
        haskey(kg, id) || add!(kg, Protein(id = id, organism = "NCBITaxon:9606"))
    end
    add!(kg, Interaction(id = "string-$n", type = :Association, left = ids[1:1], right = ids[2:2],
                         confidence = s, provenance = ["prov-string-v12"]))
end
length(edges(kg))          # -> 2
isvalid(kg)                # -> true
```

`:Association` is the honest type for co-evidence scores. Use `:Binding` only for physical
interactions.

## A reaction, its catalyst and an inhibitor

Reactions transform `left` (substrates) into `right` (products). A catalyst is an `Influence` of
type `:Catalysis` whose **target is the reaction**. A drug that inhibits the enzyme is an
`Influence` on the enzyme:

```julia
kg = KnowledgeGraph()
add!(kg,
    Protein(id = "UniProt:P00533", name = "EGFR"),
    SmallMolecule(id = "CHEBI:15422", name = "ATP"),
    SmallMolecule(id = "CHEBI:456216", name = "ADP"),
    Protein(id = "UniProt:P00533_pY1068", name = "EGFR pY1068", modification_sites = ["Tyr1068"]),
    SmallMolecule(id = "CHEBI:49668", name = "gefitinib"),
    Reaction(id = "R-autophos", type = :Irreversible, ec_numbers = ["EC:2.7.10.1"],
             left = ["UniProt:P00533", "CHEBI:15422"], right = ["UniProt:P00533_pY1068", "CHEBI:456216"]),
    Influence(id = "C-egfr", type = :Catalysis, left = ["UniProt:P00533"], right = ["R-autophos"]),
    Influence(id = "D-gef", type = :Inhibition, left = ["CHEBI:49668"], right = ["UniProt:P00533"]),
    PathwayModule(id = "PW-egfr", name = "EGFR activation",
                  components = ["R-autophos", "C-egfr", "D-gef"]),
)
isvalid(kg)                                            # -> true
[e.id for e in edges_of(kg, "R-autophos"; side = :right)]   # -> ["C-egfr"]
```

Which reactions does a drug affect, directly or through the enzyme that catalyses them?

```julia
function affected_reactions(kg, drug)
    hit = Set{String}()
    for inh in edges_of(kg, drug; side = :left, type = :Inhibition), target in inh.right
        kg[target] isa Reaction && push!(hit, target)
        for cat in edges_of(kg, target; side = :left, type = :Catalysis)
            union!(hit, cat.right)
        end
    end
    return sort!(collect(hit))
end
affected_reactions(kg, "CHEBI:49668")                  # -> ["R-autophos"]
```

## Yeast–human orthology and interologs

Group orthologous proteins, then transfer interactions known in yeast to the human orthologs.
These transferred interactions are called *interologs*.

```julia
kg = KnowledgeGraph()
add!(kg,
    Taxon(id = "NCBITaxon:559292", name = "Saccharomyces cerevisiae S288C"),
    Taxon(id = "NCBITaxon:9606", name = "Homo sapiens"),
    Protein(id = "UniProt:P32598", name = "GLC7", organism = "NCBITaxon:559292"),
    Protein(id = "UniProt:P36047", name = "SDS22", organism = "NCBITaxon:559292"),
    Protein(id = "UniProt:P62136", name = "PPP1CA", organism = "NCBITaxon:9606"),
    Protein(id = "UniProt:Q15435", name = "PPP1R7", organism = "NCBITaxon:9606"),
    OrthologGroup(id = "OG-PP1", members = ["UniProt:P32598", "UniProt:P62136"],
                  species = ["NCBITaxon:559292", "NCBITaxon:9606"], clustering_method = "OrthoDB"),
    OrthologGroup(id = "OG-SDS22", members = ["UniProt:P36047", "UniProt:Q15435"],
                  species = ["NCBITaxon:559292", "NCBITaxon:9606"], clustering_method = "OrthoDB"),
    Interaction(id = "y1", type = :Binding, left = ["UniProt:P32598"], right = ["UniProt:P36047"]),
)

orthologs(kg, p, taxon) = [m for g in ortholog_groups_of(kg, p) for m in g.members
                           if m != p && kg[m].organism == taxon]

n = 0
for e in collect(edges(kg)), a in e.left, b in e.right
    kg[a].organism == "NCBITaxon:559292" || continue
    for a2 in orthologs(kg, a, "NCBITaxon:9606"), b2 in orthologs(kg, b, "NCBITaxon:9606")
        global n += 1
        add!(kg, Interaction(id = "interolog-$n", type = e.type, left = [a2], right = [b2],
                             confidence = 0.5, db_xrefs = [expand(kg, e.id)]))
    end
end
[(e.left[1], e.right[1]) for e in edges(kg) if startswith(e.id, "interolog")]
# -> [("UniProt:P62136", "UniProt:Q15435")]
```

## Merging two sources field by field

`integrate!` merges the fields of resources that share an id (and of those linked through
`external_ids`): it unites vector fields and fills in empty ones. See
[Integrating data sources](10-integration.md) for the full story.

```julia
uniprot = KnowledgeGraph(); add!(uniprot, Protein(id = "UniProt:P04637", name = "p53",
    cellular_location = [:Nucleus], provenance = ["prov-uniprot"]))
curated = KnowledgeGraph(); add!(curated, Protein(id = "UniProt:P04637",
    cellular_location = [:Cytoplasm, :Nucleus], sequence = "MEEPQSDPSV", provenance = ["prov-lab"]))
merged = deepcopy(uniprot)
integrate!(merged, curated)
merged["UniProt:P04637"].cellular_location     # -> [:Nucleus, :Cytoplasm]
merged["UniProt:P04637"].provenance            # -> ["prov-uniprot", "prov-lab"]
merged["UniProt:P04637"].sequence              # -> "MEEPQSDPSV"
```

Keeping one `Provenance` id per source, as here, means you can always trace which source said what.

## Exporting to Cytoscape

```julia
write("demo.sif", "TP53\tactivation\tMDM2\nMDM2\tinhibition\tTP53\nATM\tphosphorylation\tTP53\n")
kg = load("demo.sif")
foreach(n -> n.name = n.id, nodes(kg))

save("network.graphml", kg)                    # File → Import → Network from File
save("node_attributes.csv", kg; types = NodeEntity)   # File → Import → Table from File
```

In Cytoscape, style edges by the `type` attribute, for example red for `Inhibition` and green for
`Activation`. The node table joins on `id`. If some edges target other edges, export with
`reify = true` so that nothing is lost.

## Cleaning up dangling references

After deleting resources or importing partial data, remove references that point nowhere:

```julia
kg = KnowledgeGraph()
add!(kg, Protein(id = "A"), Interaction(id = "e", type = :Binding, left = ["A"], right = ["B", "C"]),
         PathwayModule(id = "PW", components = ["e", "gone"]))

function prune_dangling!(kg)
    for r in kg, p in schema(r)
        p.kind === :ref || continue
        v = getfield(r, p.field)
        if v isa AbstractVector
            filter!(id -> haskey(kg, id), v)
        elseif v !== nothing && !haskey(kg, v)
            setfield!(r, p.field, nothing)
        end
    end
    return kg
end
prune_dangling!(kg)
kg["e"].right                       # -> String[]
kg["PW"].components                 # -> ["e"]
```

`e` now has no right participant. Run `validate` afterwards and decide whether to delete such edges.

## Converting files in bulk

```julia
mkpath("in"); mkpath("out")
save("in/a.ttl", kg); save("in/b.json", kg)
for f in readdir("in")
    g = load(joinpath("in", f))
    save(joinpath("out", splitext(f)[1] * ".jsonld"), g)
end
readdir("out")                      # -> ["a.jsonld", "b.jsonld"]
```

## Summarising a graph

```julia
function summary_table(kg)
    counts = Dict{String,Int}()
    foreach(r -> counts[class_name(r)] = get(counts, class_name(r), 0) + 1, kg)
    types = Dict{Symbol,Int}()
    foreach(e -> types[e.type] = get(types, e.type, 0) + 1, edges(kg))
    issues = validate(kg)
    return (resources = sort!(collect(counts)), relationship_types = sort!(collect(types)),
            errors = count(i -> i.severity === :error, issues),
            warnings = count(i -> i.severity === :warning, issues))
end
summary_table(kg).warnings          # -> 0
```
