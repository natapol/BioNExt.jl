# 1. Getting started

## Installation

BioNExt.jl needs Julia 1.10 or newer. It has no dependencies outside the standard library.

The package is in Julia's General registry:

```
pkg> add BioNExt
```

To work on the package itself, clone it and start Julia in its folder:

```sh
git clone https://github.com/natapol/BioNExt.jl
cd BioNExt.jl
julia --project=.
```

Then load it:

```julia
using BioNExt
```

## Your first graph

A `KnowledgeGraph` is a container of resources, each identified by an `id` string. This example
models part of the p53–MDM2 feedback loop: MDM2 binds p53 and marks it for degradation, and
p53 switches on expression of the *MDM2* gene.

```julia
kg = KnowledgeGraph()

add!(kg,
    Taxon(id = "NCBITaxon:9606", name = "Homo sapiens"),
    Gene(id = "HGNC:6973", name = "MDM2", organism = "NCBITaxon:9606"),
    Protein(id = "UniProt:P04637", name = "p53", organism = "NCBITaxon:9606",
            cellular_location = [:Nucleus, :Cytoplasm]),
    Protein(id = "UniProt:Q00987", name = "MDM2", organism = "NCBITaxon:9606"),
)
```

Each constructor takes keyword arguments named after its fields. Only `id` is required, plus
`type` for edges. Every other field is optional and starts empty. The full list of fields for each type
is in the [schema reference](../reference/schema.md).

Now the relationships. An edge has a relationship `type` from the controlled vocabulary,
participants on the `left` (the source) and on the `right` (the target), and optional
`weight` and `confidence`:

```julia
add!(kg,
    Interaction(id = "e-bind", type = :Binding,
                left = ["UniProt:Q00987"], right = ["UniProt:P04637"], confidence = 0.98),
    Influence(id = "e-ubq", type = :Ubiquitination,
              left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
    Influence(id = "e-expr", type = :RegulateExpression,
              left = ["UniProt:P04637"], right = ["HGNC:6973"]),
    PathwayModule(id = "p53-mdm2", name = "p53–MDM2 feedback loop",
                  components = ["e-bind", "e-ubq", "e-expr"]),
)
```

Show a summary:

```julia
show(stdout, MIME"text/plain"(), kg)
# KnowledgeGraph with 8 resources
#   Gene           1
#   Influence      2
#   Interaction    1
#   PathwayModule  1
#   Protein        2
#   Taxon          1
```

## Check it against the ontology

```julia
issues = validate(kg)
@assert isempty(issues)
```

`validate` returns a list of `Issue`s. The list is empty here, but the validator knows, for instance,
that `:RegulateExpression` must target a `Gene`. Point it at the protein and you get an error:

```julia
kg["e-expr"].right = ["UniProt:Q00987"]
validate(kg)
# -> [ERROR [e-expr] right participant "UniProt:Q00987" is a Protein; :RegulateExpression expects Gene]
kg["e-expr"].right = ["HGNC:6973"]        # put it back
```

## Ask questions

```julia
neighbors(kg, "UniProt:Q00987"; direction = :out)      # -> ["UniProt:P04637"]
[e.id for e in edges_of(kg, "UniProt:P04637")]         # -> ["e-bind", "e-ubq", "e-expr"]
[p.id for p in pathways_of(kg, "e-bind")]              # -> ["p53-mdm2"]
```

## Save and load

```julia
save("p53.ttl", kg)          # RDF Turtle
save("p53.jsonld", kg)       # JSON-LD
save("p53.tsv", kg)          # a spreadsheet-friendly table
save("p53.graphml", kg)      # for Cytoscape / Gephi

kg2 = load("p53.ttl")
@assert keys(kg2) == keys(kg)
```

## Where next

- [Concepts](02-concepts.md) explains how the model maps onto the ontology. Read it before
  modelling your own data.
- [Input and output](06-io.md) covers importing your existing tables and files.
