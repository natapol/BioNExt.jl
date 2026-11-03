# 9. Maintaining and extending

This chapter is for the people who maintain BioNExt.jl: keeping it in step with the ontology,
adding properties, classes and file formats, and running the tests.

## Project layout

```
BioNExt/
├── Project.toml               package metadata (depends only on Dates)
├── src/
│   ├── BioNExt.jl             module: includes and exports
│   ├── vocabulary.jl          Term, term, terms, is_a
│   ├── vocabulary_data.jl     GENERATED: namespaces, version, all vocabulary terms
│   ├── types.jl               the @bionext macro and every resource type
│   ├── schema.jl              SCHEMA: field ⇄ ontology property table, inverse pairs
│   ├── graph.jl               KnowledgeGraph, add!, queries
│   ├── validate.jl            validate, Issue
│   ├── resolve.jl             entity resolution: duplicates, resolve!, integrate!, find_entity
│   ├── algorithms.jl          degree_of, shortest_path, network_components, subgraph, to_graph stub
│   ├── audit.jl               data-quality audit: checks, ID_PATTERNS, assert_quality
│   ├── stats.jl               graph_stats, stats_table, write_stats
│   ├── rdf.jl                 IRI/Literal/Triple, triples, from_triples, Turtle and N-Triples writers, N-Triples reader
│   └── io/
│       ├── turtle.jl          Turtle reader
│       ├── json.jl            JSON reader/writer, native JSON, JSON-LD, assign!
│       ├── rdfxml.jl          RDF/XML writer
│       ├── tables.jl          CSV/TSV
│       ├── network.jl         network(), SIF, GraphML, DOT
│       └── io.jl              load/save/load!, format registry
├── ext/
│   └── BioNExtGraphsExt.jl    Graphs.jl extension (to_graph), loaded when Graphs.jl is
├── ontology/                  the BioNExt ontology this version is generated from (.ttl and .owl)
├── scripts/
│   ├── generate_vocabulary.jl regenerate src/vocabulary_data.jl from the ontology
│   └── check_ontology.jl      report differences between the schema and the ontology
├── docs/
│   ├── generate_reference.jl  regenerate docs/src/reference/*.md
│   └── src/                   this manual
├── test/
│   ├── runtests.jl            core tests (includes the two below)
│   ├── io_tests.jl            every file format
│   ├── graph_tests.jl         indexing (randomized indexed-vs-scan check), algorithms, Graphs.jl extension
│   ├── resolve_tests.jl       entity resolution
│   ├── quality_tests.jl       audit (one broken example per check) and statistics
│   └── docs_tests.jl          runs every code example in the manual
└── examples/p53_pathway.jl
```

The design rests on one idea: **`SCHEMA` in `src/schema.jl` is the single source of truth**
for how fields map to the ontology. Validation, RDF output and input, JSON, tables and the
reference docs all read from it. Add a property there and every format supports it with no
further changes.

## Running the tests

The tests cover the data model and queries, indexing, algorithms, entity resolution, validation,
every file format, a parse of `ontology/bne_ontology.ttl`, a check that the generated vocabulary matches
it, and every code example in this manual. `Pkg.test()` runs everything, including the Graphs.jl
extension tests (Graphs.jl is a test-only dependency). Running `test/runtests.jl` directly skips
those tests, with a message, unless Graphs.jl can be loaded.

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The companion package BioNExtSources (a separate repository) uses a few internal functions
(`BioNExt.merge_into!`, `BioNExt.assign!`, `BioNExt.CLASS_TYPES`). Run its tests too after changing them.

### Keeping the index correct

Any new function that changes reference fields must keep the index current. Use `add!` and
`delete_many!`, which maintain it. If you edit fields in place, call `index_remove!(kg, r)`
before and `index_add!(kg, r)` after, or `reindex!(kg)` at the end, as `resolve!` does. The
randomized test in `graph_tests.jl` compares indexed and scanned results after random edits.
Extend it when you add a query.

## When the ontology changes

Run these steps from the package folder after replacing `ontology/bne_ontology.ttl` (and `ontology/bne_ontology.owl`) with the new release:

```sh
# 1. regenerate namespaces, version and vocabulary terms
julia scripts/generate_vocabulary.jl ontology/bne_ontology.ttl

# 2. compare the ontology with the schema
julia --project=. scripts/check_ontology.jl ontology/bne_ontology.ttl

# 3. fix whatever step 2 reports (see the sections below), then repeat step 2

# 4. regenerate the reference pages of this manual, and run the tests
julia --project=. docs/generate_reference.jl
julia --project=. test/runtests.jl
```

`check_ontology.jl` reports:

- a changed `bn:`/`bnr:` namespace or `owl:versionInfo` (fixed by step 1)
- an ontology class with no Julia type ([add a class](#adding-a-class))
- an ontology property whose domain class has no field for it ([add a property](#adding-a-property))
- a field whose ontology property no longer exists ([remove or rename it](#renaming-or-removing-a-property))
- a datatype mismatch, such as `xsd:integer` in the package but `xsd:float` in the ontology
- a different number of vocabulary terms (fixed by step 1)

The script exits with status 1 when it finds problems, so you can run it in CI.

Changes in **vocabulary terms** (new relationship types, cellular locations, region types, new
`expectedLeft`/`expectedRight` annotations) need step 1 only. The validator reads them from the
generated data.

### The two vocabulary files

The vocabulary lives in two files, which are maintained differently.

**`src/vocabulary_data.jl` is generated; never edit it by hand.** Step 1 rewrites the whole file,
so hand edits would be lost (it starts with `# AUTO-GENERATED … do not edit`). From the ontology
the script reads:

- the `bn:` and `bnr:` namespaces from the `@prefix` lines, which become `BN` and `BNR`;
- `owl:versionInfo`, which becomes `ONTOLOGY_VERSION`;
- every controlled-vocabulary term, that is, each `bn:` block typed `owl:NamedIndividual`. For each
  term it takes the label, comment, parent class (`rdfs:subClassOf`), `rdfs:seeAlso`, and
  `bn:expectedLeft`/`bn:expectedRight`. These become `VOCABULARY_DATA`.

**`src/vocabulary.jl` is hand-written, and ordinary ontology changes don't touch it.** It defines
the `Term` struct and the functions built on it: `term`, `terms`, `is_term` and `is_a`. It reads
`VOCABULARY_DATA` from the generated file.

You only change it when the ontology gains a new *kind* of annotation on terms that you want in
Julia, say a `bn:synonym` on each term. That takes three steps:

1. **`src/vocabulary.jl`:** add a field to `Term`, for example `synonyms::Vector{String}`.
2. **`scripts/generate_vocabulary.jl`:** in the loop over blocks, extract the annotation with the
   existing helpers (`curie_list` for IRIs, `first_string` for literals), and write it into the
   `Term(...)` call the script prints. The arguments must be in the same order as the struct's fields.
3. **Regenerate** with step 1, then run the tests. They include a check that the generated data match
   the ontology.

Two caveats:

- **What the generator can read:** it is not a general Turtle parser. It splits the file into
  blank-line-separated blocks in the layout that rdflib and Protégé write. If you save the TTL from a
  tool that lays it out differently, terms could be missed. The vocabulary test would then fail.
- **What the generator doesn't cover:** new classes and properties. Those go in `types.jl` and
  `schema.jl` by hand. `check_ontology.jl` lists them, and the sections below walk through each case.

### When the namespace changes

Step 1 updates `BN` and `BNR`, and every output then uses the new namespace. To keep reading
files written with the old namespace, add the old one to `LEGACY_NAMESPACES` in `src/rdf.jl`:

<!-- norun -->

```julia
const LEGACY_NAMESPACES = [
    "https://jamiesilakom618.github.io/BioNExt/ontology/" => BN,
    "https://jamiesilakom618.github.io/BioNExt/resource/" => BNR,
    # "https://old.example.org/bionext/ontology/" => BN,     # add older namespaces here
]
```

## Adding a property

Suppose the ontology gains a datatype property `bn:hasHalfLife` with domain `bn:Protein` and
range `xsd:float`.

**1. Add a field** to the struct in `src/types.jl`, with a comment naming the property:

<!-- norun -->

```julia
@bionext struct Protein <: NodeEntity
    # ... existing fields ...
    half_life::OptFloat = nothing             # bn:hasHalfLife
end
```

Use `OptString`/`OptFloat`/`OptInt` or `Union{T,Nothing} = nothing` for single values, and
`Vector{T} = T[]` for multi-valued properties. If the domain is a superclass (`NodeEntity`,
`EdgeEntity`, `Entity`), add the field to `NODE_FIELDS`, `EDGE_FIELDS` or `ENTITY_FIELDS` at the
top of `types.jl` instead. Every subtype then inherits it.

**2. Map it** in `SCHEMA` (`src/schema.jl`):

<!-- norun -->

```julia
Protein => [
    # ... existing props ...
    Prop(:half_life, "hasHalfLife", :literal, "float");
],
```

Choose the kind as described in [Concepts](02-concepts.md#properties-become-fields):

- `:ref` with a tuple of allowed Julia types, for properties that point at resources in the graph
- `:term` with the root term, for controlled vocabulary
- `:iri` with `nothing`, for external identifiers
- `:literal` with the XSD local name, for data values

For superclass properties, add the `Prop` to `ENTITY_PROPS`, `NODE_PROPS` or `edge_props`.

**3. If it is an inverse pair**, store only one direction and register the other in `INVERSES`
(`"inverseName" => "storedName"`). Readers then accept both directions and writers emit both.
Add a query function in `graph.jl` if users need the derived direction.

**4. Add a test**, regenerate the reference docs, and run the tests. Every format now supports
the field.

## Renaming or removing a property

When the ontology renames a property, change the predicate string in `SCHEMA` and keep the Julia
field name, so user code doesn't break. To stay able to read old files, map the old predicate in
`from_triples`. The simplest way is an entry in `SUBPROPERTIES` (`"oldName" => "newName"`), which
the reader already uses to redirect predicates.

When the ontology removes a property, remove the `Prop` from `SCHEMA`, and, after a deprecation
period, the field.

## Adding a class

Suppose the ontology gains `bn:Lipid`, a subclass of `bn:NodeEntity` with a property `bn:hasLipidMaps`.

1. Define the type in `src/types.jl`. `@bionext` adds `id` and the inherited fields:

   <!-- norun -->

   ```julia
   """A lipid (`bn:Lipid`)."""
   @bionext struct Lipid <: NodeEntity
       lipidmaps_id::OptString = nothing     # bn:hasLipidMaps
   end
   ```

2. Register its ontology name in `CLASS_NAMES` (`Lipid => "Lipid"`).
3. Add `Lipid => [NODE_PROPS; Prop(:lipidmaps_id, "hasLipidMaps", :literal, "string")]` to `SCHEMA`.
4. Export `Lipid` from `src/BioNExt.jl`, and add it to `ORDER` and `GROUPS` in
   `docs/generate_reference.jl`. The generator stops with an error if an exported name is
   missing from `GROUPS`.
5. If the class appears in the range of an existing property (for example, `members` may now
   include lipids), extend that `Prop`'s range tuple.
6. If the class carries special rules (as `EdgeEntity` and variants do), add a `_validate(kg, r::Lipid, err)` method in `validate.jl`.
7. Optionally, give it a shape in `DOT_SHAPES` (`src/io/network.jl`).

## Adding a validation rule

Rules that apply to a type go into a `_validate(kg, r::T, err)` method in `src/validate.jl`.
Call `err(r, "message")` for each violation. Messages should start with the field name,
explain what was expected, and include the offending value. Add a test in the "validation"
test set and an example in [Validation](05-validation.md).

## Adding an audit check

A check is a function `(kg, emit, opts) -> nothing` that calls `emit(severity, id, message)`
for each problem it finds. Add it to `AUDIT_CHECKS` in `src/audit.jl` with a one-line
description, add a row to the tables in `audit`'s docstring and in [chapter 11](11-quality.md),
and add tests with broken data to `test/quality_tests.jl`. Use `:error` only for data that are
certainly wrong. To support another database's identifier format, add its pattern to `ID_PATTERNS`.

## Adding a file format

1. Write `write_<fmt>(io::IO, kg; kwargs...)`, `read_<fmt>(io::IO; kwargs...)`, or both, in a
   new file under `src/io/`. Build writers on `triples(kg)` (RDF-like formats), `network(kg)`
   (network formats) or the struct fields. Build readers on `from_triples` or `new_resource` +
   `assign!`, so that type conversion matches the other formats.
2. Include the file in `src/BioNExt.jl` **before** `io/io.jl`, and export its functions.
3. Register it in `src/io/io.jl`: add the extensions to `EXTENSIONS`, the functions to `READERS` and
   `WRITERS`, and the writer to the loop that adds path-accepting methods.
4. Add a round-trip test to `test/io_tests.jl`, a section to [Input and output](06-io.md), and
   the functions to `GROUPS` in `docs/generate_reference.jl`.

## Updating this manual

- **Hand-written pages:** `docs/src/index.md` and `docs/src/manual/*.md`. Code in fenced
  ```` ```julia ```` blocks is executed by `test/docs_tests.jl`, page by page, in an empty
  temporary directory. Write `expr  # -> value` on a single line to have the result checked too.
  To show Julia code that shouldn't run (for example an edit to the source), put the line
  `<!-- norun -->` directly above the block. It is invisible when the page is rendered.
- **Generated pages:** `docs/src/reference/*.md`. Never edit them by hand; run
  `docs/generate_reference.jl`. They are built from docstrings, `SCHEMA` and `VOCABULARY`, so
  improving a docstring improves the API page.
- The pages are plain Markdown, readable on GitHub/GitLab or in any editor. They also follow
  Documenter.jl's layout (`docs/src`), so a `docs/make.jl` with `makedocs(sitename = "BioNExt.jl", pages = [...])`
  can turn them into a website.

## Known issues in the ontology (version 1.0.0)

These are recorded here so they can be raised with the ontology authors:

- Label typos: `bn:RegulateExpression` "Regulates xpression", `bn:ControlledVocabulary`
  "Contolled vocabulary", `bn:hasChemFormula` "has chemical formular".
- `bn:hasOrthologGroup` (Protein → OrthologGroup) is not declared `owl:inverseOf bn:hasProteinMember`.
  The package treats the two as inverses.
- `bn:variantType` and `bn:modificationType` have range `bn:ControlledVocabulary`, but the
  vocabulary has no mutation or modification terms. The package accepts external ids (SO, PSI-MOD) for them.
- `bn:hasHGVS` has domain `bn:NodeEntity`, but HGVS is only meaningful for variants of
  DNA, RNA, genes and proteins. The package lets any node carry it and requires it only on variants.
