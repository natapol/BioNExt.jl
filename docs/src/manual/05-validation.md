# 5. Validation

Constructors only check Julia types: a `Gene` can't have a `Float64` name. They don't check
whether the data makes sense under the ontology. That is the job of `validate(kg)`, which checks the
whole graph and returns a vector of `Issue`s. Validate after building or importing a graph,
before you publish it.

```julia
using BioNExt

kg = KnowledgeGraph()
add!(kg,
    Gene(id = "HGNC:6973", name = "MDM2"),
    Protein(id = "UniProt:P04637", name = "p53"),
    Protein(id = "UniProt:Q00987", name = "MDM2"),
    Interaction(id = "e1", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
)
validate(kg)               # -> Issue[]
isvalid(kg)                # -> true
```

## Issues

Each `Issue` has three fields:

| Field | Meaning |
|---|---|
| `severity` | `:error`, meaning the data contradicts the ontology, or `:warning`, meaning a reference to something not in the graph |
| `id` | the resource the issue is about |
| `message` | a human-readable explanation that starts with the field name |

`isvalid(kg)` is `true` when there are no errors. Pass `strict = true` to also fail on warnings.

```julia
add!(kg, Influence(id = "e2", type = :Activation, left = ["UniProt:Q13315"], right = ["UniProt:P04637"]))
issues = validate(kg)
issues[1].severity          # -> :warning
issues[1].message           # -> "left: \"UniProt:Q13315\" is not in the graph"
isvalid(kg)                 # -> true
isvalid(kg; strict = true)  # -> false
```

Print the issues for a readable report:

```julia
foreach(println, issues)
# WARNING [e2] left: "UniProt:Q13315" is not in the graph
```

## What is checked

### 1. References point at the right class

Each reference field has a range taken from the ontology's `rdfs:range` and `owl:allValuesFrom`
restrictions. The [schema reference](../reference/schema.md) lists the range of every field.

```julia
add!(kg, Interaction(id = "bad1", type = :Binding, left = ["UniProt:P04637"], right = ["e1"]))
foreach(println, filter(i -> i.id == "bad1", validate(kg)))
# ERROR [bad1] right: "e1" is an Interaction, expected NodeEntity
# ERROR [bad1] right participant "e1" is an Interaction; :Binding expects NodeEntity
```

One mistake can produce more than one issue. Here the edge breaks both the range of
`Interaction.right` and the expected right participant of `:Binding`.

Ranges you are likely to run into:

| Field | Allowed |
|---|---|
| `Interaction.right`, `Reaction.right` | node entities only |
| `Influence.right`, any edge's `left` | nodes, edges or pathways |
| `variant_of` | the **same class** as the variant (a protein variant of a protein) |
| `subunits` | node entities; `subcomplexes`: complexes |
| `templates` (Gene) / `encoded_by` (Rna) | `Dna` |
| `cofactors` | `SmallMolecule` |
| `members` (OrthologGroup) | `Gene` or `Protein` |
| `organism`, `species` | `Taxon` |
| `provenance` | `Provenance` |

### 2. Vocabulary terms exist and have the right kind

```julia
add!(kg, Protein(id = "bad2", cellular_location = [:Nucleus, :Binding, :Nucleous]))
foreach(println, filter(i -> i.id == "bad2", validate(kg)))
# ERROR [bad2] cellular_location: :Binding is not a CellularLocation
# ERROR [bad2] cellular_location: unknown vocabulary term :Nucleous
```

`type` must be under `:RelationshipType`, `cellular_location` under `:CellularLocation`, the
`regions` of `Dna` and `Gene` under `:DnaRegion`, and the `regions` of `Rna` under `:RnaRegion`.

### 3. Participants match the relationship type

Many relationship types declare which classes they expect on each side (`bn:expectedLeft`
and `bn:expectedRight`). The validator applies them:

```julia
add!(kg, Influence(id = "bad3", type = :RegulateExpression,
                   left = ["UniProt:P04637"], right = ["UniProt:Q00987"]))
last(validate(kg)).message
# -> "right participant \"UniProt:Q00987\" is a Protein; :RegulateExpression expects Gene"
```

The fix here is to point the edge at the gene (`HGNC:6973`) instead of the protein. Run
`term(:X).expected_left` and `term(:X).expected_right` to look these constraints up, or see the
[vocabulary reference](../reference/vocabulary.md#relationship-types-relationshiptype). Types
with no expectation, such as `:Phosphorylation`, accept any participant.

### 4. Edges are complete

Every edge needs at least one `left` and one `right` participant. The ontology states this with
`owl:someValuesFrom`.

```julia
add!(kg, Influence(id = "bad4", type = :Activation, right = ["UniProt:P04637"]))
last(validate(kg)).message     # -> "edge has no left participant (hasLeft)"
```

### 5. Variants are well-formed

A resource with `variant_of` must have an `hgvs` description, and it can't be a variant of itself:

```julia
add!(kg, Protein(id = "bad5", variant_of = ["UniProt:P04637"]))
last(validate(kg)).message     # -> "variant (isVariantOf) must have an HGVS description (hasHGVS)"
```

### 6. Coordinates are positive

`start_pos`, `end_pos` and `length` are `xsd:positiveInteger` in the ontology:

```julia
add!(kg, Dna(id = "bad6", start_pos = 0))
last(validate(kg)).message     # -> "start_pos: must be a positive integer, got 0"
```

### 7. References resolve (warning)

A reference to an id that isn't in the graph is a warning, not an error. This is common
while importing data in several steps, or when a graph deliberately points at resources
held elsewhere.

## Working with many issues

`validate` returns plain data, so slice it however you like:

```julia
issues = validate(kg)
errors = filter(i -> i.severity === :error, issues)
length(errors)                                  # -> 8

by_resource = Dict{String,Vector{String}}()
for i in issues
    push!(get!(by_resource, i.id, String[]), i.message)
end
sort(collect(keys(by_resource)))
# -> ["bad1", "bad2", "bad3", "bad4", "bad5", "bad6", "e2"]

# drop everything that has an error
foreach(id -> delete!(kg, id), unique(i.id for i in errors))
isvalid(kg)                                     # -> true
```

To fail fast in a pipeline:

```julia
function check(kg)
    errs = filter(i -> i.severity === :error, validate(kg))
    isempty(errs) || error("graph has $(length(errs)) errors:\n" * join(string.(errs), "\n"))
    return kg
end
check(kg);
```

## What is not checked

- **Disjointness** of classes (a resource can't be both a `Gene` and a `Protein`) is guaranteed
  by construction, because every resource has exactly one Julia type.
- **External identifiers** (`db_xrefs`, `publications`, `functions`, …) are not looked up.
- **Cardinality** beyond the cases above is not enforced. For example, the ontology doesn't
  limit how many organisms an entity has, and neither does the library.
- **Cycles** in `subcomplexes` or pathway `components` are allowed. `components(...; recursive = true)`
  handles them safely.
