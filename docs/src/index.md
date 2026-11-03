# BioNExt.jl manual

BioNExt.jl is a Julia library for building, querying, validating and exchanging biological
knowledge graphs that follow the [BioNExt ontology](https://w3id.org/bionext/ontology/).
The ontology describes genes, DNA, RNA, proteins, small molecules, complexes, ortholog
groups, the interactions, influences and reactions between them, the pathways that group them,
and the provenance of the data. BioNExt.jl gives each of these a Julia type and a set of tools to work with them.

## Contents

**Guide**

1. [Getting started](manual/01-getting-started.md): install the package and build a first graph in five minutes
2. [Concepts](manual/02-concepts.md): the data model, identifiers, references, vocabulary terms, inverse properties
3. [Building graphs](manual/03-building.md): create, edit, delete and merge resources; namespaces and prefixes
4. [Querying](manual/04-querying.md): find resources, traverse edges, indexing, graph algorithms, Graphs.jl
5. [Validation](manual/05-validation.md): what is checked against the ontology and how to fix each kind of issue
6. [Input and output](manual/06-io.md): every file format, with its options and limitations
7. [RDF mapping](manual/07-rdf.md): how the graph becomes triples, and how to use it with triple stores and SPARQL
8. [Cookbook](manual/08-cookbook.md): complete recipes for common tasks
9. [Maintaining and extending](manual/09-maintenance.md): ontology updates, the vocabulary files, adding properties or classes, the project layout, tests
10. [Integrating data sources](manual/10-integration.md): find and merge duplicate entities across databases, settle conflicts, look up identifiers
11. [Data quality and statistics](manual/11-quality.md): audit the data beyond the ontology rules, gate pipelines on quality, collect and export graph statistics

**Reference** (generated from the code)

- [Schema](reference/schema.md): every type, field, ontology property and allowed range
- [Vocabulary](reference/vocabulary.md): all relationship types, cellular locations and region types
- [API](reference/api.md): every exported function and type

## At a glance

```julia
using BioNExt

kg = KnowledgeGraph()
add!(kg,
    Protein(id = "UniProt:P04637", name = "p53", cellular_location = [:Nucleus]),
    Protein(id = "UniProt:Q00987", name = "MDM2"),
    Interaction(id = "e1", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
)
validate(kg)                         # -> Issue[] (no problems)
neighbors(kg, "UniProt:Q00987")      # -> ["UniProt:P04637"]
save("p53.ttl", kg)                  # Turtle; also .jsonld .json .csv .tsv .sif .graphml .dot ...
```

## Conventions in this manual

- Code blocks are runnable. They are executed by the test suite (`test/docs_tests.jl`), in
  order, page by page, in an empty temporary directory. Comments of the form `# -> value`
  show what an expression returns.
- `kg` is always a `KnowledgeGraph`.
- "Resource" means anything stored in a graph: every entity, plus provenance records.
- Ontology names are written with the `bn:` prefix (`bn:hasLeft`), and Julia names in code style (`left`).
