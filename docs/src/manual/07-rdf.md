# 7. RDF mapping

This chapter describes exactly how a graph becomes RDF, and how RDF becomes a graph. Read it
if you plan to query exported data with SPARQL, load it into a triple store, or import RDF
produced by other tools.

```julia
using BioNExt, Dates

kg = KnowledgeGraph()
add!(kg,
    Protein(id = "UniProt:P04637", name = "p53", cellular_location = [:Nucleus]),
    Protein(id = "UniProt:P04637_p.R175H", hgvs = "p.Arg175His", variant_of = ["UniProt:P04637"]),
    Interaction(id = "e1", type = :Binding, left = ["UniProt:P04637"], right = ["UniProt:P04637"],
                weight = 0.5),
    PathwayModule(id = "PW", components = ["e1"]),
    Provenance(id = "prov1", collected_date = DateTime(2026, 10, 1)),
)
```

## Triples

`triples(kg)` returns the graph as `Triple(subject, predicate, object)` values. Subjects and
predicates are `IRI`s. Objects are `IRI`s or `Literal`s (`value`, `datatype`, `lang`). Every
writer is built on this function.

```julia
ts = triples(kg)
ts[1]
# -> Triple(IRI("https://identifiers.org/UniProt:P04637"),
#           IRI("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"),
#           IRI("https://w3id.org/bionext/ontology/Protein"))
```

## Mapping rules

| Julia | RDF |
|---|---|
| a resource with id `x` | subject `expand(kg, x)` (see [identifiers](02-concepts.md#identifiers)) |
| its type | `rdf:type bn:<Class>` (`MolecularComplex` → `bn:Complex`) |
| a field | the predicate `bn:<property>` from the [schema](../reference/schema.md) |
| a vector field | one triple per element |
| `nothing` / empty vector | no triple |
| `:ref` or `:iri` value | an IRI, `expand(kg, value)` |
| `:term` value `:Nucleus` | the IRI `bn:Nucleus` |
| `String` | a plain literal (`xsd:string`) |
| `Float64` | a literal typed `xsd:float` |
| `Int` | a literal typed `xsd:integer`, or `xsd:positiveInteger` for coordinates |
| `DateTime` | a literal typed `xsd:dateTime`, `yyyy-mm-ddTHH:MM:SS` |
| `OrthologGroup.members` | `bn:hasProteinMember` for proteins, `bn:hasMember` otherwise |

```julia
w = only(t for t in ts if t.predicate == IRI(BioNExt.BN * "hasWeight"))
w.object        # -> Literal("0.5", "http://www.w3.org/2001/XMLSchema#float", nothing)
```

### Inverse triples

With `inverses = true` (the default for `triples` and every RDF writer), the inverse of each
stored one-way property is added as well:

| Stored | Also written |
|---|---|
| `V bn:isVariantOf P` | `P bn:hasVariant V` |
| `PW bn:pathwayComponent X` | `X bn:participatesInPathway PW` |
| `C bn:hasSubunit X` | `X bn:memberOfComplex C` |
| `G bn:hasProteinMember P` | `P bn:hasOrthologGroup G` |
| `Gene bn:hasTemplate D` | `D bn:annotates Gene` |
| `R bn:encodedBy D` | `D bn:encodes R` |
| `E bn:hasLeft N` / `E bn:hasRight N` (node `N`) | `N bn:participatesInEdge E` |

```julia
has(s, p, o) = Triple(IRI(expand(kg, s)), IRI(BioNExt.BN * p), IRI(expand(kg, o))) in ts
has("UniProt:P04637", "hasVariant", "UniProt:P04637_p.R175H")      # -> true
has("e1", "participatesInPathway", "PW")                           # -> true
```

This lets SPARQL queries follow both directions without a reasoner. The ontology declares
most of these pairs with `owl:inverseOf`, so an OWL reasoner would infer the same triples.

## From RDF back to a graph

`from_triples(ts; base, prefixes)` is the inverse of `triples`, and every RDF reader uses it.

1. **Resources.** Every subject with `rdf:type bn:<C>`, where `C` is a concrete class, becomes a
   resource. Its id is `compact(kg, subject)`.
2. **Fields.** Every `bn:` predicate on a known subject is looked up in that type's schema, and
   its value is converted to the field type.
3. **Inverse predicates** (`bn:hasVariant`, `bn:participatesInPathway`, …) are turned around and
   stored on the canonical side. `bn:hasProteinMember` goes into `members`.
4. **Ignored triples:** subjects without a BioNExt type, predicates outside `bn:`,
   `bn:participatesInEdge` (it is implied by edges), and predicates not in the schema for that
   class.
5. **Edges** read without `bn:hasType` keep the placeholder type `:Unknown`, which `validate` reports.
6. **Old namespaces:** IRIs in the BioNExt namespace before `w3id.org` are upgraded.

```julia
back = from_triples(ts)
back["UniProt:P04637_p.R175H"].variant_of     # -> ["UniProt:P04637"]
```

So you can build a graph from RDF that another tool produced, as long as it uses the BioNExt
vocabulary:

```julia
nt = """
<http://example.org/g1> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/bionext/ontology/Gene> .
<http://example.org/g1> <https://w3id.org/bionext/ontology/name> "BRCA1" .
<http://example.org/g1> <http://www.w3.org/2000/01/rdf-schema#label> "ignored" .
"""
g = read_ntriples(IOBuffer(nt))
g["http://example.org/g1"].name                # -> "BRCA1"
g2 = read_ntriples(IOBuffer(nt); prefixes = Dict("ex" => "http://example.org/"))
collect(keys(g2))                              # -> ["ex:g1"]
```

## Querying exported data with SPARQL

Export to N-Triples or Turtle, then load the file into any SPARQL store: Apache Jena Fuseki,
GraphDB, Virtuoso, Oxigraph, or rdflib in Python. Load `bne_ontology.ttl` alongside if you
want labels and class hierarchy in your queries.

```julia
save("graph.nt", kg)
```

Every protein–protein binding, with names:

```sparql
PREFIX bn: <https://w3id.org/bionext/ontology/>
SELECT ?a ?aName ?b ?bName WHERE {
  ?e a bn:Interaction ; bn:hasType bn:Binding ; bn:hasLeft ?a ; bn:hasRight ?b .
  ?a a bn:Protein . ?b a bn:Protein .
  OPTIONAL { ?a bn:name ?aName } OPTIONAL { ?b bn:name ?bName }
}
```

Everything in a pathway, including sub-pathways (this relies on the inverse triples written by default):

```sparql
PREFIX bn: <https://w3id.org/bionext/ontology/>
PREFIX bnr: <https://w3id.org/bionext/resource/>
SELECT DISTINCT ?x WHERE { ?x bn:participatesInPathway+ bnr:PW }
```

Every co-expression edge, including positive and negative correlations, using the ontology's class hierarchy (load
the ontology too):

```sparql
PREFIX bn: <https://w3id.org/bionext/ontology/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
SELECT ?e ?type WHERE { ?e bn:hasType ?type . ?type rdfs:subClassOf* bn:CoExpression }
```

Variants and their HGVS descriptions:

```sparql
PREFIX bn: <https://w3id.org/bionext/ontology/>
SELECT ?ref ?variant ?hgvs WHERE { ?ref bn:hasVariant ?variant . ?variant bn:hasHGVS ?hgvs }
```

## Notes for publishing

- Set a `base` of your own (`KnowledgeGraph(base = "https://example.org/mydata/")`) so that bare
  ids get stable, resolvable IRIs.
- The output contains only data triples. If your publication guidelines ask for a dataset
  description (`void:Dataset`, `dct:license`, `owl:imports` of the ontology), add those triples
  to the exported file with a separate step.
- External ids written as CURIEs become identifiers.org IRIs, which resolve on the web and are
  widely understood. Register prefixes if you prefer a database's native IRIs, such as
  `http://purl.uniprot.org/uniprot/`.
