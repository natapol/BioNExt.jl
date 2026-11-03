# 3. Building graphs

```julia
using BioNExt, Dates
```

## Creating a graph

```julia
kg = KnowledgeGraph()
```

`KnowledgeGraph` takes two optional keyword arguments:

| Keyword | Default | Purpose |
|---|---|---|
| `base` | `"https://w3id.org/bionext/resource/"` | the namespace under which bare ids are minted |
| `prefixes` | `Dict()` | extra CURIE prefixes, e.g. `Dict("up" => "http://purl.uniprot.org/uniprot/")` |

## Adding resources

`add!` takes any number of resources and returns the graph. Resources are kept in the order
they were added, and every listing and every file export follows that order.

```julia
add!(kg,
    Gene(id = "HGNC:11998", name = "TP53"),
    Protein(id = "UniProt:P04637", name = "p53"),
)
push!(kg, Protein(id = "UniProt:Q00987", name = "MDM2"))   # same as add!
length(kg)                                                 # -> 3
```

Ids must be unique and non-empty. Adding an id that already exists throws an `ArgumentError`
unless you ask for replacement:

```julia
try
    add!(kg, Gene(id = "HGNC:11998"))
catch err
    println(err.msg)      # duplicate id "HGNC:11998" (use replace = true to overwrite)
end
add!(kg, Gene(id = "HGNC:11998", name = "TP53", synonyms = ["p53"]); replace = true)
```

You can add references before their targets exist. `validate` reports anything still
missing at the end as a warning.

## Reading and editing

A graph behaves like an ordered dictionary of resources:

```julia
kg["UniProt:P04637"]                    # the Protein
haskey(kg, "UniProt:P04637")            # -> true
get(kg, "nope", nothing)                # -> nothing
keys(kg)                                # ids in insertion order
for r in kg                             # iterate resources in insertion order
    println(class_name(r), "\t", r.id)
end
```

Resources are mutable, so edit them in place:

```julia
p = kg["UniProt:P04637"]
p.cellular_location = [:Nucleus]
push!(p.publications, "PMID:2046748")
p.sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
```

Don't change `id` on a resource that is already in a graph, because the graph indexes resources by id.
To rename, remove the resource and add it again, then update every reference (see `references` below).

To set a field from untyped data such as strings from a file, use `assign!`. It converts the
value to the field's type, and is what the table and JSON readers use:

```julia
sm = SmallMolecule(id = "CHEBI:15422")
assign!(sm, :name, "ATP")
assign!(sm, :net_charge, "-4")               # "-4" -> -4
assign!(sm, :synonyms, ["adenosine 5'-triphosphate"])
sm.net_charge                                # -> -4
```

## Building every kind of resource

### Genes, DNA and RNA

```julia
add!(kg,
    Dna(id = "chr17:7661779-7687538", chromosome = "17", start_pos = 7661779,
        end_pos = 7687538, regions = [:Promoter, :Exon, :Intron]),
    Rna(id = "ENST00000269305", name = "TP53-201", encoded_by = ["chr17:7661779-7687538"],
        regions = [:MessengerRNA], length = 2512),
)
kg["HGNC:11998"].templates = ["chr17:7661779-7687538"]   # gene --hasTemplate--> DNA
kg["HGNC:11998"].functions = ["GO:0003700"]               # DNA-binding TF activity
```

`regions` holds region-type terms (`:Promoter`, `:CDS`, `:ThreePrimeUTR`, … for DNA and genes;
`:MessengerRNA`, `:MicroRNA`, … for RNA). Coordinates must be positive integers.

### Proteins, isoforms and modifications

```julia
add!(kg,
    Protein(id = "UniProt:P04637-2", name = "p53 isoform 2", isoform_of = ["UniProt:P04637"]),
)
kg["UniProt:P04637"].modification_sites = ["Ser15", "Ser20"]
kg["UniProt:P04637"].modification_type = ["MOD:00046"]      # O-phospho-L-serine
```

### Small molecules and cofactors

```julia
add!(kg,
    SmallMolecule(id = "CHEBI:29105", name = "zinc(2+)", formula = "Zn", net_charge = 2),
    SmallMolecule(id = "CHEBI:15422", name = "ATP", formula = "C10H16N5O13P3",
                  inchikey = "ZKHQWZAMYRWXGA-KQYNXXCUSA-N", smiles = "Nc1ncnc2...", net_charge = -4),
)
kg["UniProt:P04637"].cofactors = ["CHEBI:29105"]
```

### Complexes

A complex lists its monomer `subunits` (any node entities) and its `subcomplexes` (other complexes):

```julia
add!(kg, MolecularComplex(id = "CPX-1", name = "p53–MDM2 complex",
                          subunits = ["UniProt:P04637", "UniProt:Q00987"]))
```

### Ortholog groups and taxa

```julia
add!(kg,
    Taxon(id = "NCBITaxon:9606", name = "Homo sapiens"),
    Taxon(id = "NCBITaxon:10090", name = "Mus musculus"),
    Protein(id = "UniProt:P02340", name = "p53 (mouse)", organism = "NCBITaxon:10090"),
    OrthologGroup(id = "OrthoDB:12345", clustering_method = "OrthoDB v11", confidence = 0.95,
                  members = ["UniProt:P04637", "UniProt:P02340"],
                  species = ["NCBITaxon:9606", "NCBITaxon:10090"]),
)
```

`members` may contain genes and proteins. Protein members are written to RDF as
`bn:hasProteinMember`, and gene members as `bn:hasMember`.

### Edges

All three edge types share the same fields. `left` and `right` are vectors because the ontology
allows several participants on each side, as in a reaction with several substrates:

```julia
add!(kg,
    SmallMolecule(id = "CHEBI:456216", name = "ADP"),
    Reaction(id = "R-phospho", type = :Phosphorylation, ec_numbers = ["EC:2.7.11.1"],
             left = ["UniProt:P04637", "CHEBI:15422"], right = ["UniProt:P04637", "CHEBI:456216"]),
    Interaction(id = "I-bind", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"],
                confidence = 0.98, publications = ["PMID:1535557"]),
    Influence(id = "F-inh", type = :Inhibition, left = ["CPX-1"], right = ["R-phospho"], weight = 0.7),
)
```

`F-inh` is an influence whose target is another edge. Only `Influence` may do this; an
`Interaction` or `Reaction` with an edge on the right fails validation.

### Pathways

A pathway lists its `components`: nodes, edges and sub-pathways.

```julia
add!(kg,
    PathwayModule(id = "PW-core", name = "p53 core regulation", components = ["I-bind", "F-inh"]),
    PathwayModule(id = "PW-p53", name = "p53 signalling", components = ["PW-core", "R-phospho"]),
)
components(kg, "PW-p53"; recursive = true)   # -> ["PW-core", "R-phospho", "I-bind", "F-inh"]
```

### Provenance

```julia
add!(kg, Provenance(id = "prov-curation-2026", sources = ["PMID:1535557"],
                    created_by = ["https://orcid.org/0009-0003-1344-2449"],
                    collected_date = DateTime(2026, 10, 30), modified_date = now()))
push!(kg["I-bind"].provenance, "prov-curation-2026")
```

## Removing resources

`delete!` removes one resource. References to it from other resources are **not** removed, so
call `references` first to see what points at it:

```julia
references(kg, "CHEBI:29105")      # -> ["UniProt:P04637" => :cofactors]
delete!(kg, "CHEBI:29105")
filter!(!=("CHEBI:29105"), kg["UniProt:P04637"].cofactors)
```

## Merging graphs

`merge!(kg, other)` copies every resource of `other` into `kg`, along with its prefixes. When an
id exists in both, the existing resource is kept, unless you pass `replace = true`:

```julia
other = KnowledgeGraph()
add!(other, Protein(id = "UniProt:P04637", name = "from another source"),
            Protein(id = "UniProt:Q9Y6K9", name = "NEMO"))
merge!(kg, other)
kg["UniProt:P04637"].name            # -> "p53"   (kept)
merge!(kg, other; replace = true)
kg["UniProt:P04637"].name            # -> "from another source"
```

Merging copies the resource objects, so after a merge both graphs share them. If you need
independent copies, `deepcopy` the source graph first.

`merge!` matches resources by id only. To combine sources field by field, and to merge
entities that appear under different ids (UniProt vs Ensembl vs STRING), use `integrate!`
(see [Integrating data sources](10-integration.md)).

## Namespaces and prefixes

Bare ids are minted under `kg.base`. To publish your graph under your own namespace, set the base:

```julia
lab = KnowledgeGraph(base = "https://example.org/mylab/")
add!(lab, Gene(id = "g1"))
expand(lab, "g1")                    # -> "https://example.org/mylab/g1"
```

To make CURIEs resolve to specific namespaces instead of identifiers.org, register prefixes. They
are also used to shorten IRIs in Turtle and JSON-LD output, and to compact IRIs back into ids on input:

```julia
up = KnowledgeGraph(prefixes = Dict("up" => "http://purl.uniprot.org/uniprot/"))
add!(up, Protein(id = "up:P04637"))
expand(up, "up:P04637")              # -> "http://purl.uniprot.org/uniprot/P04637"
up.prefixes["taxon"] = "http://purl.uniprot.org/taxonomy/"
```

The prefixes `bn`, `bnr`, `rdf`, `rdfs` and `xsd` are always defined. Ids that are full IRIs are
left untouched.
