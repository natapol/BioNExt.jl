# Compare the package's schema with an ontology file and report what is out of sync.
#
#   julia --project=. scripts/check_ontology.jl [path/to/bne_ontology.ttl]
#
# Reports ontology classes without a Julia type, properties whose domain class lacks a
# field for them, schema predicates the ontology no longer defines, datatype mismatches,
# and namespace/version drift. Exits with status 1 if anything needs attention.

using BioNExt
using BioNExt: SCHEMA, INVERSES, SUBPROPERTIES, CLASS_TYPES, BN, BNR, ONTOLOGY_VERSION, RDF, RDFS, XSD

const OWL = "http://www.w3.org/2002/07/owl#"
path = isempty(ARGS) ? joinpath(@__DIR__, "..", "ontology", "bne_ontology.ttl") : ARGS[1]
ts, pfx = parse_turtle(path)

problems = String[]
report(msg) = (push!(problems, msg); println("  ✗ ", msg))

bysubj = Dict{String,Vector{Triple}}()
foreach(t -> push!(get!(bysubj, t.subject.value, Triple[]), t), ts)
objects(s, p) = [t.object for t in get(bysubj, s, Triple[]) if t.predicate.value == p]
local_name(iri) = startswith(iri, BN) ? iri[length(BN)+1:end] : iri

# A class expression: a named class, or a blank node with owl:unionOf ( ... ).
function class_list(o)
    o isa Literal && return String[]
    startswith(o.value, "_:") || return [local_name(o.value)]
    out = String[]
    for u in objects(o.value, OWL * "unionOf")
        n = u.value
        while n != RDF * "nil"
            push!(out, local_name(only(objects(n, RDF * "first")).value))
            n = only(objects(n, RDF * "rest")).value
        end
    end
    return out
end

println("Namespaces and version")
get(pfx, "bn", "") == BN || report("bn: is $(get(pfx, "bn", "missing")) in the ontology but $BN in the package")
get(pfx, "bnr", "") == BNR || report("bnr: is $(get(pfx, "bnr", "missing")) in the ontology but $BNR in the package")
version = [t.object.value for t in ts if endswith(t.predicate.value, "versionInfo") && t.subject.value == BN]
isempty(version) || version[1] == ONTOLOGY_VERSION ||
    report("ontology version is $(version[1]) but the package was generated from $ONTOLOGY_VERSION")

println("Classes")
istype(s, c) = IRI(c) in objects(s, RDF * "type")
vocab_term(s) = istype(s, OWL * "NamedIndividual") || any(o -> o isa IRI && startswith(o.value, BN), objects(s, RDF * "type"))
classes = [local_name(s) for s in keys(bysubj) if startswith(s, BN) && istype(s, OWL * "Class") && !vocab_term(s)]
for c in sort(classes)
    haskey(CLASS_TYPES, c) || report("class bn:$c has no Julia type")
end

println("Properties")
derived = Set([collect(keys(INVERSES)); "participatesInEdge"; "participatesIn"])
schema_preds = Dict{Type,Set{String}}(T => Set(p.pred for p in props) for (T, props) in SCHEMA)
foreach(s -> union!(s, keys(SUBPROPERTIES)), values(schema_preds))
ontology_props = Set{String}()
for s in keys(bysubj)
    (istype(s, OWL * "ObjectProperty") || istype(s, OWL * "DatatypeProperty")) && startswith(s, BN) || continue
    p = local_name(s)
    push!(ontology_props, p)
    p in derived && continue
    for d in vcat(class_list.(objects(s, RDFS * "domain"))...)
        D = get(CLASS_TYPES, d, nothing)
        D === nothing && (report("bn:$p has domain bn:$d, which has no Julia type"); continue)
        for (T, preds) in schema_preds
            T <: D && !(p in preds) && report("bn:$p (domain bn:$d) has no field on $T")
        end
    end
    ranges = [o.value for o in objects(s, RDFS * "range") if o isa IRI && startswith(o.value, XSD)]
    for (T, props) in SCHEMA, q in props
        q.pred == p && q.kind === :literal && !isempty(ranges) && XSD * q.range != ranges[1] &&
            report("$T.$(q.field): xsd:$(q.range) in the package, $(local_name(replace(ranges[1], XSD => "xsd:"))) in the ontology")
    end
end
for (T, props) in SCHEMA, q in props
    q.pred in ontology_props || report("$T.$(q.field) maps to bn:$(q.pred), which the ontology does not define")
end

println("Vocabulary")
n_terms = count(s -> startswith(s, BN) && istype(s, OWL * "NamedIndividual"), keys(bysubj))
n_terms == length(VOCABULARY) ||
    report("ontology has $n_terms vocabulary terms, package has $(length(VOCABULARY)); run scripts/generate_vocabulary.jl")

if isempty(problems)
    println("\nThe package is in sync with $(basename(path)) (version $ONTOLOGY_VERSION).")
else
    println("\n$(length(problems)) problem(s) found.")
    exit(1)
end
