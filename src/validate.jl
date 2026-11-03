"""
    Issue(severity, id, message)

A validation finding. `severity` is `:error` (violates the ontology) or
`:warning` (e.g. a reference to an id not present in the graph).
"""
struct Issue
    severity::Symbol
    id::String
    message::String
end

Base.show(io::IO, i::Issue) = print(io, uppercase(string(i.severity)), " [", i.id, "] ", i.message)

_names(types) = join(map(class_name, types), "|")
_a(name) = (occursin(r"^[AEIOU]", name) ? "an " : "a ") * name

"""
    validate(kg) -> Vector{Issue}

Check the graph against the BioNExt ontology:
- reference ranges (`rdfs:range` and `owl:allValuesFrom` restrictions, e.g. a
  `Reaction`'s right side must be NodeEntities; variants must be of the same class)
- controlled-vocabulary membership (relationship types, cellular locations, regions)
- `bn:expectedLeft` / `bn:expectedRight` of each edge's relationship type
- edges need a left and a right participant; variants need exactly one HGVS string
- positive integer coordinates
- dangling references (warning)
"""
function validate(kg::KnowledgeGraph)
    issues = Issue[]
    err(r, msg) = push!(issues, Issue(:error, r.id, msg))
    for r in kg
        for p in schema(r), v in field_values(r, p.field)
            if p.kind === :ref
                if !haskey(kg, v)
                    push!(issues, Issue(:warning, r.id, "$(p.field): \"$v\" is not in the graph"))
                elseif !any(T -> kg[v] isa T, p.range)
                    err(r, "$(p.field): \"$v\" is $(_a(class_name(kg[v]))), expected $(_names(p.range))")
                end
            elseif p.kind === :term
                is_term(v) || (err(r, "$(p.field): unknown vocabulary term :$v"); continue)
                is_a(v, p.range) || err(r, "$(p.field): :$v is not a $(p.range)")
            elseif p.kind === :literal && p.range == "positiveInteger"
                v > 0 || err(r, "$(p.field): must be a positive integer, got $v")
            end
        end
        _validate(kg, r, err)
    end
    return issues
end

_validate(kg, r, err) = nothing

function _validate(kg, e::EdgeEntity, err)
    isempty(e.left) && err(e, "edge has no left participant (hasLeft)")
    isempty(e.right) && err(e, "edge has no right participant (hasRight)")
    is_term(e.type) && is_a(e.type, :RelationshipType) || return
    t = term(e.type)
    for (side, ids, expected) in ((:left, e.left, t.expected_left), (:right, e.right, t.expected_right))
        isempty(expected) && continue
        allowed = [CLASS_TYPES[string(c)] for c in expected]
        for id in ids
            haskey(kg, id) || continue
            any(T -> kg[id] isa T, allowed) ||
                err(e, "$side participant \"$id\" is $(_a(class_name(kg[id]))); :$(e.type) expects $(_names(allowed))")
        end
    end
end

function _validate(kg, n::NodeEntity, err)
    hasfield(typeof(n), :variant_of) || return
    !isempty(n.variant_of) && n.hgvs === nothing &&
        err(n, "variant (isVariantOf) must have an HGVS description (hasHGVS)")
    n.id in n.variant_of && err(n, "is declared a variant of itself")
end

"""    isvalid(kg; strict = false) -> Bool

True if `validate` reports no errors (and no warnings when `strict`)."""
Base.isvalid(kg::KnowledgeGraph; strict::Bool = false) =
    all(i -> i.severity === :warning && !strict, validate(kg))
