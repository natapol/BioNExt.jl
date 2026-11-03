# BioNExt.jl

[![CI](https://github.com/natapol/BioNExt.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/natapol/BioNExt.jl/actions/workflows/CI.yml)

Build, query, validate and serialize biological knowledge graphs structured by the
[BioNExt ontology](https://w3id.org/bionext/ontology/) (v1.0.0, namespace `https://w3id.org/bionext/`).
Only depends on the `Dates` standard library.

## Installation

BioNExt.jl needs Julia 1.10 or newer. It isn't registered yet, so add it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/natapol/BioNExt.jl")
```

To import data from UniProt, KEGG, IntAct, BioGRID, SIGNOR, Reactome and BioCyc, see the companion
package [BioNExtSources.jl](https://github.com/natapol/BioNExtSources.jl).

## Documentation

The full manual is in [`docs/src/index.md`](docs/src/index.md). It has:

- **Guide chapters:** getting started, concepts, building, querying, validation, input/output, RDF mapping, a cookbook, and maintenance.
- **Generated reference pages:** the [schema](docs/src/reference/schema.md), the [vocabulary](docs/src/reference/vocabulary.md) and the [API](docs/src/reference/api.md).

The test suite runs every code example in the manual.

## Model

Each ontology class has a Julia type; each ontology property is a field.

| Julia type | Ontology class | Own fields (besides `id` and the inherited ones) |
|---|---|---|
| `Taxon` | `bn:Taxon` | — |
| `Provenance` | `bn:Provenance` | `collected_date`, `modified_date`, `created_by`, `derived_from`, `sources` |
| `Gene` | `bn:Gene` | `functions`, `templates`→Dna, `regions`, `variant_of`, `variant_type` |
| `Dna` | `bn:Dna` | `chromosome`, `sequence`, `start_pos`, `end_pos`, `length`, `regions`, `variant_of`, `variant_type` |
| `Rna` | `bn:Rna` | `sequence`, coordinates, `regions`, `encoded_by`→Dna, `modification_sites`, `modification_type`, `variant_of`, `variant_type` |
| `Protein` | `bn:Protein` | `sequence`, `cofactors`→SmallMolecule, `isoform_of`, `modification_sites`, `modification_type`, `variant_of`, `variant_type` |
| `MolecularComplex` | `bn:Complex` | `subunits`, `subcomplexes` |
| `SmallMolecule` | `bn:SmallMolecule` | `cas_number`, `formula`, `inchi`, `inchikey`, `net_charge`, `smiles` |
| `OrthologGroup` | `bn:OrthologGroup` | `clustering_method`, `confidence`, `members`→Gene/Protein, `species`→Taxon |
| `Interaction`, `Influence`, `Reaction` | `bn:EdgeEntity` subclasses | `Reaction`: `ec_numbers` |
| `PathwayModule` | `bn:PathwayModule` | `components` |

Fields inherited from the ontology superclasses:

- **All entities:** `name`, `synonyms`, `organism`, `db_xrefs`, `external_ids`, `provenance`, `publications`.
- **Node entities:** `cellular_location` and `hgvs`.
- **Edge entities:** `type` (required, a `RelationshipType` term), `left`, `right`, `weight` and `confidence`.

Run `schema(T)` to list the exact field → `bn:` predicate mapping. That one table drives validation, export and import.

- **References** are identifier strings. An id can be a bare name (`"TP53"`, minted under `bnr:`), a CURIE (`"UniProt:P04637"`, expanded through `kg.prefixes` or identifiers.org) or a full IRI.
- **Controlled-vocabulary terms** are `Symbol`s: `:Binding`, `:Nucleus`, `:Promoter`, and so on. They are generated from the ontology, and `terms(:RelationshipType)` lists them.
- **Inverse properties** are stored in one direction only, and the other direction is derived. For example, `bn:hasVariant` comes from `variant_of`, `bn:participatesInPathway` from `components`, and `bn:memberOfComplex` from `subunits`.

## Usage

```julia
using BioNExt

kg = KnowledgeGraph()
add!(kg,
    Gene(id = "HGNC:6973", name = "MDM2"),
    Protein(id = "UniProt:P04637", name = "p53", cellular_location = [:Nucleus]),
    Protein(id = "UniProt:Q00987", name = "MDM2"),
    Interaction(id = "e1", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
    Influence(id = "e2", type = :RegulateExpression, left = ["UniProt:P04637"], right = ["HGNC:6973"]),
    PathwayModule(id = "p53-loop", components = ["e1", "e2"]),
)

validate(kg)                       # Vector{Issue}; isvalid(kg) for a Bool
edges_of(kg, "UniProt:P04637")     # bn:participatesInEdge view
neighbors(kg, "UniProt:Q00987"; direction = :out, type = :Binding)
pathways_of(kg, "e1"); components(kg, "p53-loop"; recursive = true)

save("graph.ttl", kg)              # format chosen by extension
kg2 = load("graph.jsonld")
```

### Speed and graph algorithms

```julia
index!(kg)                         # queries become lookups instead of scans
shortest_path(kg, "UniProt:Q00987", "HGNC:6973")      # typed, directed BFS
degree_of(kg, "UniProt:P04637"; type = :Binding)
network_components(kg); subgraph(kg, "UniProt:P04637"; hops = 2)

import Graphs                      # optional: loads the Graphs.jl extension
g = to_graph(kg)                   # (graph, ids, vertex) for any Graphs.jl algorithm
```

On 100,000 proteins with 300,000 interactions, a `neighbors` query takes about 80 ms by scanning
and about 0.04 ms with the index. The index updates itself through `add!`, `delete!`, `merge!` and
`integrate!`. After editing reference fields directly, call `reindex!(kg)`.

### Data quality and statistics

```julia
report = audit(kg)                 # beyond validate: identifiers, sequences, modification sites,
passed(report)                     # scores, duplicates, cycles, orphans, missing annotation
assert_quality(kg)                 # throws a QualityError for pipelines and tests
s = graph_stats(kg)                # composition, degrees, hubs, components, clustering,
write_stats("stats.md", s)         # sources and their overlap, annotation coverage
```

### Merging data from several databases

```julia
integrate!(kg, string_graph; prefer = ["UniProt", "HGNC"])   # copy in, then merge duplicates
duplicates(kg); resolve!(kg)                                 # or clean up one graph
find_entity(kg, "Ensembl:ENSP00000269305")                   # look up by any known identifier
```

Entities of the same class that share an id or an `external_ids` entry (`bn:hasExternalId`)
are merged:

- Vector fields are united.
- Conflicting values are settled by a policy you choose, and every conflict is recorded.
- Every reference in the graph is rewritten to the surviving id.
- Edges that become identical are merged too.

`examples/p53_pathway.jl` contains a fuller example.

## Input and output formats

`load(path)`, `load!(kg, path)` and `save(path, kg)` pick the format from the file extension. Pass `format = :tsv` (or another format name) to override it, and run `formats()` to list what's supported. Each format also has its own `read_*`/`write_*` function, which accepts either an `IO` or a path.

| Format | Extensions | Read | Write | Notes |
|---|---|:-:|:-:|---|
| Turtle | `.ttl` | ✓ | ✓ | Full Turtle 1.1 parser: prefixes, blank nodes, collections, long strings, numbers |
| N-Triples | `.nt` | ✓ | ✓ | |
| JSON-LD | `.jsonld` | ✓ | ✓ | Reads compacted or expanded JSON-LD with an inline `@context`; remote contexts are not fetched |
| RDF/XML | `.rdf` `.owl` `.xml` | | ✓ | For Protégé, Jena, rdflib |
| JSON | `.json` | ✓ | ✓ | Plain, non-RDF; one object per resource, using the Julia field names |
| CSV / TSV | `.csv` `.tsv` `.txt` | ✓ | ✓ | One row per resource, with `class` and `id` columns; multiple values are joined by `\|` |
| SIF | `.sif` | ✓ | ✓ | Cytoscape Simple Interaction Format |
| GraphML | `.graphml` | | ✓ | For Cytoscape, Gephi, yEd and Graphs.jl |
| DOT | `.dot` `.gv` | | ✓ | Graphviz; inhibitory edges are drawn with a tee arrowhead |

- **RDF output** (Turtle, N-Triples, JSON-LD, RDF/XML) writes inverse properties in both directions by default. Pass `inverses = false` to write each one only once.
- **RDF input** keeps the prefixes declared in the file, so ids come back as `up:P04637` and not as full IRIs.
- **Network formats** turn BioNExt's n-ary edges into simple arcs:
  - By default, each edge becomes one arc per left/right pair.
  - With `reify = true`, each edge becomes a node of its own, linked as left → edge → right. Nothing is lost this way, including edges that target other edges or pathways.
  - `network(kg; reify)` returns the same arcs as a Julia value.

### Importing your own tables

```julia
# a gene list: name the class and map columns to fields
kg = load("genes.csv"; type = Gene,
          mapping = Dict("symbol" => :name, "hgnc_id" => :id, "aliases" => :synonyms))

# an edge list read into the same graph; rows without an id get one (influence-1, ...)
load!(kg, "regulation.tsv"; type = Influence,
      mapping = Dict("source" => :left, "target" => :right, "effect" => :type, "score" => :confidence))

# a SIF file; relationship names are matched to vocabulary terms ("activation" → :Activation)
load!(kg, "ppi.sif"; node_type = Protein, relationships = Dict("pp" => :Binding))

save("edges.tsv", kg; types = EdgeEntity)    # export only the edges
validate(kg)                                 # check the imported data against the ontology
```

### What `validate` checks

- **Property ranges and `owl:allValuesFrom` restrictions.** For example, the right side of an `Interaction` or `Reaction` must be NodeEntities, and a variant must belong to the same class as the entity it varies.
- **Vocabulary membership.** Relationship types, cellular locations, and DNA or RNA regions must be valid terms.
- **The `bn:expectedLeft` and `bn:expectedRight` annotations** of each relationship type. For example, `:RegulateExpression` expects a `Gene` on the right.
- **Required parts.** Every edge needs a left and a right participant, and every variant needs an HGVS string.
- **Coordinates** must be positive integers.
- **Dangling references** are reported as warnings.

## Development

```sh
julia --project=. test/runtests.jl
julia scripts/generate_vocabulary.jl ontology/bne_ontology.ttl          # after the ontology changes
julia --project=. scripts/check_ontology.jl ontology/bne_ontology.ttl   # report schema/ontology differences
julia --project=. docs/generate_reference.jl                      # refresh the reference pages
```

See [Maintaining and extending](docs/src/manual/09-maintenance.md) for the full procedure.

## License

The package code is under the MIT license (see [LICENSE](LICENSE)). The ontology files in
[`ontology/`](ontology) are the BioNExt Ontology 1.0.0, under CC BY 4.0; see
[ontology/README.md](ontology/README.md) for attribution.
