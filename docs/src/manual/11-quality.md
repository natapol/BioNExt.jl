# 11. Data quality and statistics

`validate` checks that a graph obeys the ontology's rules. That is necessary, but it can't
tell whether the data themselves are right. A protein sequence with a digit in it, a
phosphosite on a residue the sequence doesn't have, or the same protein stored under two
identifiers all pass `validate`. This chapter covers two tools for looking at the data:

- `audit(kg)` checks data quality and returns a report of findings.
- `graph_stats(kg)` collects statistics describing the graph.

```julia
using BioNExt

kg = KnowledgeGraph()
add!(kg,
    Taxon(id = "NCBITaxon:9606", name = "Homo sapiens"),
    Provenance(id = "prov-uniprot"), Provenance(id = "prov-lit"),
    Protein(id = "UniProt:P04637", name = "p53", organism = "NCBITaxon:9606",
            sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPEN", modification_sites = ["Ser15", "Ser20", "Thr18"],
            provenance = ["prov-uniprot", "prov-lit"]),
    Protein(id = "UniProt:Q00987", name = "MDM2", organism = "NCBITaxon:9606", provenance = ["prov-uniprot"]),
    Protein(id = "UniProt:Q13315", name = "ATM", organism = "NCBITaxon:9606", provenance = ["prov-lit"]),
    Interaction(id = "e1", type = :Binding, left = ["UniProt:Q00987"], right = ["UniProt:P04637"],
                confidence = 0.98, publications = ["PMID:1535557"], provenance = ["prov-lit"]),
    Influence(id = "e2", type = :Phosphorylation, left = ["UniProt:Q13315"], right = ["UniProt:P04637"],
              publications = ["PMID:9733515"], provenance = ["prov-lit"]),
)
```

## Auditing

```julia
report = audit(kg)
passed(report)                         # -> true
count_findings(report; severity = :error)   # -> 0
```

Show the report to read it. Each check is listed with its counts and an example finding:

```julia
show(stdout, MIME"text/plain"(), report)
# Audit of 7 resources: passed — 0 errors, 0 warnings, 1 notes
#   ontology            ok
#   identifiers         ok
#   ...
#   annotation          1 info
#       e.g. UniProt:Q00987: ...
```

Now introduce some typical problems: a modification site that isn't in the sequence, a
malformed accession, a confidence of 98 that should be 0.98, and the same protein again
under its Ensembl id.

```julia
push!(kg["UniProt:P04637"].modification_sites, "Ser99", "Ser16")
add!(kg,
    Protein(id = "UniProt:Q0098", name = "MDM2 typo"),
    Protein(id = "Ensembl:ENSP00000269305", name = "TP53", external_ids = ["UniProt:P04637"]),
    Interaction(id = "e3", type = :Binding, left = ["UniProt:Q13315"], right = ["UniProt:P04637"], confidence = 98.0),
)
report = audit(kg)
passed(report)                         # -> false
[(f.check, f.severity, f.id) for f in findings(report; severity = :error)]
# -> [(:modification_sites, :error, "UniProt:P04637")]
[f.message for f in findings(report; check = :modification_sites)]
# -> ["modification site Ser99 is beyond the sequence (29 residues)",
#     "modification site Ser16: the sequence has Q at position 16"]
```

The other problems are warnings, listed in the order the checks run:

```julia
[(f.check, f.id) for f in findings(report; severity = :warning)]
# -> [(:identifiers, "UniProt:Q0098"), (:modification_sites, "UniProt:P04637"),
#     (:values, "e3"), (:duplicates, "UniProt:P04637")]
```

Each finding is a `Finding` with four fields: `check`, `severity`,
`id` and `message`. Severities mean:

| Severity | Meaning | Examples |
|---|---|---|
| `:error` | the data are wrong | invalid sequence characters, end before start, a site beyond the sequence, a NaN score, anything `validate` calls an error |
| `:warning` | probably wrong, or needs action | malformed identifiers, a site on the wrong residue, confidence outside [0, 1], duplicates, empty complexes, containment cycles |
| `:info` | worth knowing, often fine | self-loops (homodimers), cross-species edges, opposite effects between the same participants, isolated nodes, missing names or publications |

### The checks

| Check | What it looks for |
|---|---|
| `:ontology` | everything `validate` reports |
| `:identifiers` | ids and cross-references whose prefix is in `ID_PATTERNS` (UniProt, ChEBI, NCBI Gene, GO, PubMed, EC, Reactome, IntAct, KEGG, …) but whose local part has the wrong format |
| `:sequences` | characters that aren't valid for the molecule type (protein, DNA or RNA alphabets, IUPAC codes allowed); `length` that differs from the sequence |
| `:coordinates` | `end_pos` before `start_pos`; `length` inconsistent with the coordinates |
| `:modification_sites` | sites written as `Ser15` that lie beyond the sequence or on a different residue; sites not in that form (info) |
| `:values` | NaN or infinite scores; confidence outside [0, 1]; malformed InChIKey, CAS number (checksum included), formula or HGVS |
| `:duplicates` | entities that `resolve!` would merge; edges with the same class, type and participants |
| `:edges` | self-loops; edges between different organisms; activating and inhibiting edges between the same participants |
| `:structure` | complexes, pathways and ortholog groups without members; a complex or pathway that contains itself |
| `:connectivity` | nodes in no edge, complex, group or pathway; taxa and provenance records nobody uses |
| `:annotation` | entities without a name, molecules without an organism, entities without provenance, edges without publications |

Run only some checks, or leave some out. For large graphs, `max_per_check` limits how many
findings are kept; the counts are always complete:

```julia
audit(kg; checks = [:sequences, :modification_sites]).checks   # -> [:sequences, :modification_sites]
r = audit(kg; skip = [:annotation, :connectivity], max_per_check = 1)
length(findings(r; check = :modification_sites))           # -> 1
count_findings(r; check = :modification_sites)             # -> 2
```

### Identifier formats

`ID_PATTERNS` maps a prefix to a pattern for the part after the colon. Add your own databases,
or relax a pattern, by passing a modified copy:

```julia
patterns = copy(ID_PATTERNS)
patterns["SGD"] = r"^S\d{9}$"
add!(kg, Gene(id = "SGD:S000001855", name = "CDC28"), Gene(id = "SGD:CDC28", name = "CDC28 typo"))
[f.id for f in findings(audit(kg; checks = [:identifiers], patterns))]
# -> ["UniProt:Q0098", "SGD:CDC28"]
```

### As a gate in a pipeline

`assert_quality` runs the audit and throws a `QualityError`, which shows the report, when the graph
doesn't pass. Choose the level that fails: `:error` (the default), `:warning`, or `:info`.

```julia
try
    assert_quality(kg)
catch err
    err isa QualityError                                   # -> true
end
fixed = deepcopy(kg)
filter!(!in(["Ser99", "Ser16"]), fixed["UniProt:P04637"].modification_sites)   # remove the bad sites
assert_quality(fixed) |> passed                            # -> true
```

In a test suite, write `@test passed(audit(kg))`, or `@test passed(audit(kg); fail_on = :warning)`.

## Statistics

The statistics below describe the graph as the audit section left it, with its extra records included.

```julia
s = graph_stats(kg)
s.nodes, s.edges                          # -> (7, 3)
s.hubs[1]                                 # -> "UniProt:P04637" => 3
s.multi_source                            # -> 1
```

Displaying the result gives a report like this:

```julia
show(stdout, MIME"text/plain"(), s)
```

```
Knowledge graph statistics
Composition
  resources / nodes / edges             13 / 7 / 3
  classes                               Protein 5, Gene 2, Interaction 2, Provenance 2, Influence 1, Taxon 1
  relationship types                    Binding 2, Phosphorylation 1
Degree (edges per node)
  min / median / mean / max             0 / 0 / 0.8571 / 3
  isolated nodes                        4
  hubs                                  UniProt:P04637 3, UniProt:Q13315 2, UniProt:Q00987 1, …
Network
  components (largest)                  5 (3 = 42.9% of the network)
  arcs / density / reciprocity          2 / 0.04762 / 0.0%
  ...
Sources
  entities per source                   prov-lit 4, prov-uniprot 2
  entities from ≥ 2 sources             1
  largest overlaps                      prov-lit & prov-uniprot 1
Annotation coverage
  distinct publications                 2
  name (nodes, pathways, taxa)          100.0%
  sequence (proteins)                   20.0%
  publications (edges)                  66.7%
  ...
```

There are 2 arcs for 3 edges because `e2` and `e3` connect the same pair of proteins.

`GraphStats` holds the following (`?GraphStats` lists every field):

| Group | Fields |
|---|---|
| composition | `resources`, `classes`, `relationship_types`, `nodes`, `edges` |
| degree | `degree`, `in_degree`, `out_degree` (min, max, mean, median), `degree_distribution`, `hubs`, `isolated` |
| network | `components`, `largest_component`, `largest_fraction`, `arcs`, `density`, `reciprocity`, `self_loops`, `clustering` |
| structure | `pathways`, `complexes` (count, mean and max size), `pathway_depth` |
| context | `organisms`, `locations`, `id_prefixes` |
| sources | `sources` (entities per provenance record), `multi_source`, `source_overlap` (entities shared by each pair of sources) |
| annotation | `publications`, `completeness` (the fraction of resources that have each annotation) |

Some notes on how these are computed:

- **Degree** is the number of edges a node takes part in, as `degree_of` counts it.
- **Network measures** use the flattened network from `network(kg)`: one arc per left/right pair.
- **The clustering coefficient** is the mean local coefficient over nodes with at least two
  neighbours, on the undirected network. With the default `clustering = :auto`, it is skipped
  for graphs of more than 200,000 nodes.

```julia
s.degree.max, s.largest_component         # -> (3, 3)
Dict(s.completeness)["publications (edges)"]   # -> 2 / 3
s.source_overlap                          # -> [("prov-lit", "prov-uniprot") => 1]
```

### Saving

`stats_table` turns the statistics into rows (section, metric, value). `write_stats` saves them
as Markdown, TSV, CSV or JSON:

```julia
first(stats_table(s), 3)
# -> [(section = "overview", metric = "resources", value = 13),
#     (section = "overview", metric = "nodes", value = 7),
#     (section = "overview", metric = "edges", value = 3)]
write_stats("stats.md", s)                # a report with one table per section
write_stats("stats.tsv", s)               # for spreadsheets, R, pandas
write_stats("stats.json", s)              # nested by section
BioNExt.parse_json(read("stats.json", String))["overview"]["edges"]   # -> 3
```

## Typical workflow

After building or importing a graph:

```julia
report = audit(kg; skip = [:annotation])
for f in findings(report; severity = :error)
    println(f)                            # fix these first
end
res = resolve!(kg)                        # merge the duplicates the audit found
count_findings(audit(kg); check = :duplicates)   # -> 0
write_stats("graph_stats.md", graph_stats(kg))
```

Audit and statistics are fast. On a BioGRID release with 156,669 interactions and 16,695
proteins, the full audit takes about 4 s and `graph_stats` (clustering included) about 1.5 s.
