# 2. Concepts

This chapter explains how BioNExt.jl represents the ontology. Once you know these ideas,
everything else follows from them.

```julia
using BioNExt
```

## Classes become types

Each class of the ontology has a Julia type. Abstract ontology classes become abstract types,
so you can dispatch and filter on them:

```
Resource
├── Provenance                       bn:Provenance
└── Entity                           bn:Entity
    ├── Taxon                        bn:Taxon
    ├── PathwayModule                bn:PathwayModule
    ├── NodeEntity                   bn:NodeEntity
    │   ├── Gene                     bn:Gene
    │   ├── Dna                      bn:Dna
    │   ├── Rna                      bn:Rna
    │   ├── Protein                  bn:Protein
    │   ├── SmallMolecule            bn:SmallMolecule
    │   ├── OrthologGroup            bn:OrthologGroup
    │   └── MolecularComplex         bn:Complex
    └── EdgeEntity                   bn:EdgeEntity
        ├── Interaction              bn:Interaction
        ├── Influence                bn:Influence
        └── Reaction                 bn:Reaction
```

`bn:Complex` is called `MolecularComplex` because `Complex` already means complex numbers in
Julia. `class_name` gives the ontology name of any type or resource:

```julia
class_name(MolecularComplex)        # -> "Complex"
class_name(Gene(id = "g1"))         # -> "Gene"
Gene <: NodeEntity <: Entity        # -> true
```

The three edge classes are disjoint in the ontology, and they differ in what they may point at:

| Type | Meaning | `right` may contain |
|---|---|---|
| `Interaction` | an association or contact between entities (binding, co-expression, …) | node entities only |
| `Reaction` | a transformation of substrates (`left`) into products (`right`) | node entities only |
| `Influence` | a directed effect of a source on a target (activation, inhibition, …) | nodes, **other edges**, or **pathways** |

That makes `Influence` the right class for statements such as "this complex inhibits that
reaction" or "this drug suppresses that pathway".

## Properties become fields

Each ontology property whose domain is a class becomes a field of that type. Fields that the
ontology defines on a superclass are inherited:

- **Every entity** has `name`, `synonyms`, `organism`, `db_xrefs`, `external_ids`, `provenance`
  and `publications`.
- **Every node entity** also has `cellular_location` and `hgvs`.
- **Every edge** also has `type`, `left`, `right`, `weight` and `confidence`.

```julia
p = Protein(id = "UniProt:P04637", name = "p53", sequence = "MEEPQSDPSV")
p.name                       # -> "p53"
p.cellular_location          # -> Symbol[]   (empty, can hold several)
p.hgvs                       # -> nothing    (optional single value)
```

Fields that hold a single value are typed `Union{T, Nothing}` and default to `nothing`. Fields
that can hold several values are `Vector`s and default to empty. All resources are mutable:

```julia
p.name = "Cellular tumor antigen p53"
push!(p.synonyms, "TP53")
```

`schema(T)` lists the mapping from fields to ontology properties that the library uses for
validation and file conversion:

```julia
first(schema(Protein), 3)
# -> [Prop(:name, "name", :literal, "string"),
#     Prop(:synonyms, "synonyms", :literal, "string"),
#     Prop(:organism, "organism", :ref, (Taxon,))]
```

Every field has one of four kinds:

| Kind | Stored as | Example fields | Example value |
|---|---|---|---|
| `:ref` | the id of another resource in the graph | `left`, `organism`, `subunits`, `variant_of` | `"UniProt:P04637"` |
| `:term` | a controlled-vocabulary `Symbol` | `type`, `cellular_location`, `regions` | `:Binding` |
| `:iri` | an identifier outside the graph | `db_xrefs`, `publications`, `functions`, `ec_numbers` | `"PMID:2046748"`, `"GO:0003700"` |
| `:literal` | a data value | `name`, `sequence`, `weight`, `start_pos` | `"p53"`, `0.9`, `7661779` |

## Identifiers

Every resource has an `id` string. References between resources use these ids, not Julia
object references. That keeps graphs easy to serialize and lets you refer to something
before you have added it. An id can take one of three forms:

| Form | Example | Its IRI |
|---|---|---|
| bare name | `"e-bind"` | minted under the graph's base namespace: `https://w3id.org/bionext/resource/e-bind` |
| CURIE with a known prefix | `"up:P04637"` (graph has prefix `up`) | the prefix's namespace + local part |
| CURIE with an unknown prefix | `"UniProt:P04637"` | `https://identifiers.org/UniProt:P04637` |
| full IRI | `"http://purl.uniprot.org/uniprot/P04637"` | itself |

CURIEs such as `HGNC:6973`, `UniProt:P04637`, `CHEBI:15422`, `NCBITaxon:9606` and `GO:0005634`
resolve through [identifiers.org](https://identifiers.org), so they are globally meaningful
without any setup. Fields of kind `:iri` (cross-references, publications, GO functions, …)
use the same rules. See [Building graphs](03-building.md#namespaces-and-prefixes) to configure
prefixes, and `expand`/`compact` to convert by hand:

```julia
kg = KnowledgeGraph()
expand(kg, "e-bind")             # -> "https://w3id.org/bionext/resource/e-bind"
expand(kg, "UniProt:P04637")     # -> "https://identifiers.org/UniProt:P04637"
compact(kg, "https://identifiers.org/UniProt:P04637")   # -> "UniProt:P04637"
```

## Controlled vocabulary

Relationship types, cellular locations and DNA/RNA region types come from the ontology's
controlled vocabulary. They are written as `Symbol`s and organised as a hierarchy:

```julia
term(:Binding)                              # -> Term(:Binding, "Binding")
term(:Binding).parent                       # -> :RelationshipType
is_a(:NegativeCorrelation, :CoExpression)   # -> true
is_a(:NegativeCorrelation, :RelationshipType) # -> true
terms(:RnaRegion)[1:3]                      # -> [:AntisenseRNA, :CircularRNA, :DoubleStrandedRNA]
```

A `Term` also carries its `label`, its `comment` (definition), its `see_also` links to external
ontologies (MI, GO, SO, …), and for relationship types the **expected participants**: which
classes may appear on the left and right.

```julia
term(:RegulateExpression).expected_left    # -> [:Gene, :Protein, :Rna]
term(:RegulateExpression).expected_right   # -> [:Gene]
```

The [vocabulary reference](../reference/vocabulary.md) lists every term. Some fields whose range
in the ontology is `ControlledVocabulary` take external identifiers in practice, because the
ontology has no terms for them yet: `variant_type` (e.g. Sequence Ontology `"SO:1000002"`) and
`modification_type` (e.g. PSI-MOD `"MOD:00696"`).

## Inverse properties are stored once

The ontology has several pairs of inverse properties: `hasVariant`/`isVariantOf`,
`pathwayComponent`/`participatesInPathway`, `hasSubunit`/`memberOfComplex`, and more. Storing
both directions would let them disagree, so each pair is stored on **one side only**, and the
other side is computed:

| Stored field | Ontology property | The inverse, and how to get it |
|---|---|---|
| `variant.variant_of` | `bn:isVariantOf` | `bn:hasVariant`: `variants(kg, x)` |
| `pathway.components` | `bn:pathwayComponent` | `bn:participatesInPathway`: `pathways_of(kg, x)` |
| `complex.subunits` | `bn:hasSubunit` | `bn:memberOfComplex`: `complexes_of(kg, x)` |
| `group.members` | `bn:hasMember` / `bn:hasProteinMember` | `bn:hasOrthologGroup`: `ortholog_groups_of(kg, x)` |
| `gene.templates` | `bn:hasTemplate` | `bn:annotates` |
| `rna.encoded_by` | `bn:encodedBy` | `bn:encodes` |
| `edge.left`, `edge.right` | `bn:hasLeft`, `bn:hasRight` | `bn:participatesInEdge`: `edges_of(kg, x)` |

So, to say that a protein takes part in a pathway, add it to the pathway's `components`. To say
that a variant belongs to its canonical protein, set the variant's `variant_of`. When you
export to RDF, both directions are written by default, and when you import, either direction is
accepted (see [RDF mapping](07-rdf.md)).

## Variants

A variant is an ordinary node of the same class as its reference, with `variant_of` pointing at
the reference and an HGVS description in `hgvs`. The ontology requires exactly one HGVS string
per variant:

```julia
add!(kg,
    Protein(id = "UniProt:P04637", name = "p53"),
    Protein(id = "UniProt:P04637_p.R273H", name = "p53 R273H",
            variant_of = ["UniProt:P04637"], hgvs = "p.Arg273His", variant_type = ["SO:1000002"]),
)
[v.id for v in variants(kg, "UniProt:P04637")]   # -> ["UniProt:P04637_p.R273H"]
```

## Provenance

`Provenance` records describe where data came from. They aren't entities; they are attached
to entities by id through the `provenance` field:

```julia
using Dates
add!(kg, Provenance(id = "prov-uniprot", sources = ["https://www.uniprot.org"],
                    created_by = ["https://orcid.org/0000-0000-0000-0000"],
                    collected_date = DateTime(2026, 10, 1)))
push!(kg["UniProt:P04637"].provenance, "prov-uniprot")
```

## What the library does *not* do

- It isn't an OWL reasoner. It checks the restrictions that matter for data quality (see
  [Validation](05-validation.md)), but it doesn't infer new class memberships.
- It doesn't fetch remote data. Identifiers are kept as given, and nothing is looked up online.
- It doesn't store arbitrary extra properties. Triples whose predicate isn't part of the schema are
  ignored on import. To keep extra data, extend the schema (see [Maintaining](09-maintenance.md)).
