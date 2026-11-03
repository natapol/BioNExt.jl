# 10. Integrating data sources

Data about one molecule arrives under many names. UniProt calls p53 `P04637`, Ensembl calls it
`ENSP00000269305`, and STRING calls it `9606.ENSP00000269305`. Load all three and the graph
holds three proteins, each with a third of the interactions. This chapter shows how to find
such duplicates and merge them into one entity, without losing information or the record of
where it came from.

```julia
using BioNExt
```

## How sameness is decided

The ontology distinguishes two kinds of cross-reference. The library follows that distinction:

| Field | Property | Meaning | Used for matching |
|---|---|---|---|
| `external_ids` | `bn:hasExternalId` | "this identifier refers to the **same** entity" | yes |
| `db_xrefs` | `bn:hasDbXref` | "this record is related to that one" (no equivalence asserted) | only with `use_xrefs = true` |

Two resources are duplicates when they:

1. belong to the **same class**, because a gene and its protein product are never merged, even if
   one lists the other's id; and
2. share **an identifier**: their `id` or any of their `external_ids`. This works
   transitively, so if A lists B and B lists C, all three are one entity.

Identifiers are compared in a normalized form. `"UniProt:P04637"`, `"uniprot:P04637"` and
`"https://identifiers.org/UniProt:P04637"` are equal, because the identifiers.org prefix is
compared ignoring case. Different names for the same database can be declared with
`prefix_aliases`:

```julia
kg = KnowledgeGraph()
normalize_identifier(kg, "uniprot:P04637") == normalize_identifier(kg, "UniProt:P04637")  # -> true
normalize_identifier(kg, "UniProtKB:P04637"; prefix_aliases = Dict("UniProtKB" => "UniProt")) ==
    normalize_identifier(kg, "UniProt:P04637")                                            # -> true
```

Accession *versions*, as in `ENSP00000269305.4` versus `ENSP00000269305`, are different strings.
Strip versions when you load the data if they should match.

## An example: three sources

```julia
uniprot = KnowledgeGraph()
add!(uniprot,
    Provenance(id = "prov-uniprot", sources = ["https://www.uniprot.org"]),
    Protein(id = "UniProt:P04637", name = "Cellular tumor antigen p53", cellular_location = [:Nucleus],
            external_ids = ["Ensembl:ENSP00000269305"], provenance = ["prov-uniprot"]),
    Protein(id = "UniProt:Q00987", name = "E3 ubiquitin-protein ligase Mdm2",
            external_ids = ["Ensembl:ENSP00000258149"], provenance = ["prov-uniprot"]),
)

string_db = KnowledgeGraph()
add!(string_db,
    Provenance(id = "prov-string", sources = ["https://string-db.org"]),
    Protein(id = "STRING:9606.ENSP00000269305", name = "TP53",
            external_ids = ["Ensembl:ENSP00000269305"], provenance = ["prov-string"]),
    Protein(id = "STRING:9606.ENSP00000258149", name = "MDM2",
            external_ids = ["Ensembl:ENSP00000258149"], provenance = ["prov-string"]),
    Interaction(id = "string-1", type = :Binding, confidence = 0.999, provenance = ["prov-string"],
                left = ["STRING:9606.ENSP00000258149"], right = ["STRING:9606.ENSP00000269305"]),
)

literature = KnowledgeGraph()
add!(literature,
    Protein(id = "UniProt:P04637", sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"),
    Protein(id = "UniProt:Q00987"),
    Interaction(id = "pmid-1535557", type = :Binding, publications = ["PMID:1535557"],
                left = ["UniProt:Q00987"], right = ["UniProt:P04637"]),
)
```

## `integrate!`: bringing a source in

`integrate!(kg, other)` copies `other` into `kg` and resolves duplicates in one step:

```julia
kg = KnowledgeGraph()
integrate!(kg, uniprot)
res = integrate!(kg, string_db; prefer = ["UniProt"])
show(res)                              # Resolution: 2 resources merged into 2, 2 conflicts
res.groups
# -> [["UniProt:P04637", "STRING:9606.ENSP00000269305"],
#     ["UniProt:Q00987", "STRING:9606.ENSP00000258149"]]
res2 = integrate!(kg, literature; prefer = ["UniProt"])
res2.groups
# -> [["string-1", "pmid-1535557"]]
```

Here is what happened:

- The STRING proteins share an Ensembl id with the UniProt ones, so they were merged into the
  UniProt records. `prefer = ["UniProt"]` chose which id survives.
- The two sources disagree on the names (`"Cellular tumor antigen p53"` and `"TP53"`).
  The canonical record kept its own name, and each disagreement was recorded:

  ```julia
  [(c.id, c.kept, c.discarded) for c in res.conflicts][1]
  # -> ("UniProt:P04637", "Cellular tumor antigen p53", "TP53")
  ```

  The discarded name isn't lost if you want it: add it to `synonyms`.
- The literature graph uses the same UniProt ids, so its proteins were merged field by field
  straight away. The p53 sequence was added.
- Both interaction records now describe "MDM2 binds p53". After the proteins were merged they
  became identical, so they were merged as well. The surviving edge carries the STRING
  confidence *and* the literature citation:

```julia
length(kg)                             # -> 5
e = only(edges(kg))
e.id, e.confidence, e.publications     # -> ("string-1", 0.999, ["PMID:1535557"])
e.left, e.right                        # -> (["UniProt:Q00987"], ["UniProt:P04637"])
p = kg["UniProt:P04637"]
p.provenance                           # -> ["prov-uniprot", "prov-string"]
p.external_ids                         # -> ["Ensembl:ENSP00000269305", "STRING:9606.ENSP00000269305"]
isvalid(kg; strict = true)             # -> true
```

The merged ids are kept in `external_ids`, so nothing is lost and old identifiers still find the entity:

```julia
find_entity(kg, "STRING:9606.ENSP00000269305").id     # -> "UniProt:P04637"
find_entity(kg, "ensembl:ENSP00000258149").id         # -> "UniProt:Q00987"
```

`integrate!` doesn't modify `other`, and it keeps an index of `kg` up to date. Resources with
the same id but different classes can't be merged, so `integrate!` throws an `ArgumentError`
that names them.

## `resolve!`: cleaning a single graph

If the duplicates are already in one graph, for example after loading several files with
`load!`, call `resolve!` directly. `duplicates` shows what it would merge without changing anything:

```julia
kg = KnowledgeGraph()
add!(kg,
    Protein(id = "UniProt:P04637", name = "p53"),
    Protein(id = "Ensembl:ENSP00000269305", name = "TP53", external_ids = ["UniProt:P04637"]),
    Gene(id = "HGNC:11998", external_ids = ["UniProt:P04637"]),       # a gene: not merged with proteins
    Protein(id = "V1", hgvs = "p.Arg175His", variant_of = ["Ensembl:ENSP00000269305"]),
)
duplicates(kg)                         # -> [["UniProt:P04637", "Ensembl:ENSP00000269305"]]
res = resolve!(kg)
kg["V1"].variant_of                    # -> ["UniProt:P04637"]   (references are rewritten)
res.conflicts[1]
# Conflict(UniProt:P04637.name: kept "p53", discarded "TP53" from Ensembl:ENSP00000269305)
```

## What a merge does

For each group, one member becomes **canonical** and the others are folded into it.

| Step | Rule |
|---|---|
| Canonical id | the member whose prefix comes first in `prefer` (e.g. `["UniProt", "HGNC", "Ensembl"]`, ignoring case), otherwise the earliest added. `prefer` can also be a function `group -> id`. |
| Vector fields | united, in order, without repeats (synonyms, locations, provenance, publications, …) |
| Empty single-valued fields | filled in from the others (sequence, formula, …) |
| Conflicting single values | settled by `on_conflict` and recorded as a `Conflict` |
| Merged ids | added to the canonical's `external_ids`; repeats of the canonical id itself are dropped |
| References | every reference to a merged id, anywhere in the graph, now points at the canonical id; repeats are dropped |
| Edges | edges that become identical (same class, type and participant sets) are merged too, repeatedly, because merging two edges can make the influences that target them identical |
| Index | kept current |

The result is a `Resolution` with three fields: `mapping` (merged id ⇒ canonical id), `groups`
(each group, canonical id first) and `conflicts`. Keep `mapping` if other data still uses the
old ids:

```julia
res.mapping                            # -> Dict("Ensembl:ENSP00000269305" => "UniProt:P04637")
```

## Settling conflicts

By default the canonical resource keeps its value (`on_conflict = keep_first`). Pass a function
`(field, kept, other) -> value` to decide per field, for example to keep the highest confidence
and the longest sequence:

```julia
policy(field, kept, other) =
    field in (:confidence, :weight) ? max(kept, other) :
    field === :sequence ? (length(other) > length(kept) ? other : kept) : kept

kg = KnowledgeGraph()
add!(kg, Protein(id = "A"), Protein(id = "B"),
    Interaction(id = "x", type = :Binding, left = ["A"], right = ["B"], confidence = 0.4),
    Interaction(id = "y", type = :Binding, left = ["A"], right = ["B"], confidence = 0.9))
res = resolve!(kg; on_conflict = policy)
kg["x"].confidence                     # -> 0.9
res.conflicts[1].discarded             # -> 0.4
```

Every conflict is recorded even when the policy picks the other value. Review
`res.conflicts` after each integration, since they are often the first sign of a mapping error.

## Options

`resolve!` and `integrate!` take the same keywords, and `duplicates` takes the matching ones:

| Keyword | Default | Effect |
|---|---|---|
| `prefer` | `String[]` | prefix priority for the canonical id, or a function `group -> id` |
| `merge_edges` | `true` | also merge edges that become identical |
| `use_xrefs` | `false` | also treat `db_xrefs` as equivalence (use with care) |
| `prefix_aliases` | `Dict()` | prefix ⇒ prefix, for databases known under several names |
| `on_conflict` | `keep_first` | `(field, kept, other) -> value` |

## Looking up identifiers

`find_entity(kg, identifier)` returns the resource whose id or external ids match. It's
handy for annotating rows of an external table. For many lookups, build the map once:

```julia
kg = KnowledgeGraph()
add!(kg, Protein(id = "UniProt:P04637", external_ids = ["Ensembl:ENSP00000269305", "RefSeq:NP_000537"]))
idmap = identifier_map(kg)
ids = ["refseq:NP_000537", "Ensembl:ENSP00000269305", "UniProt:P99999"]
[get(idmap, normalize_identifier(kg, x), missing) for x in ids]
# -> ["UniProt:P04637", "UniProt:P04637", missing]
```

## A workflow for real data

1. Load each source into its own graph, with its own `Provenance` record attached to everything.
2. When loading, fill `external_ids` with the cross-references the source **asserts as the same
   molecule**, such as UniProt's Ensembl mappings or Ensembl's UniProt mappings. Put weaker links
   in `db_xrefs`.
3. `integrate!` the sources one by one into a main graph, with a fixed `prefer` order and a
   conflict policy.
4. Look at `res.conflicts`, and check `validate(kg)`.
5. `index!(kg)` before heavy querying (see [Querying](04-querying.md#indexing)).
