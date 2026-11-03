"""
    Term

A controlled-vocabulary term of the BioNExt ontology (a subclass of
`bn:ControlledVocabulary`, e.g. a relationship type, cellular component or region type).
`expected_left`/`expected_right` hold the `bn:expectedLeft`/`bn:expectedRight`
annotations of relationship types (empty = unconstrained).
"""
struct Term
    id::Symbol
    label::String
    comment::String
    parent::Union{Symbol,Nothing}
    see_also::Vector{String}
    expected_left::Vector{Symbol}
    expected_right::Vector{Symbol}
end

include("vocabulary_data.jl")

"""All controlled-vocabulary terms of the ontology, keyed by id (e.g. `VOCABULARY[:Binding]`)."""
const VOCABULARY = Dict{Symbol,Term}(t.id => t for t in VOCABULARY_DATA)

"""    term(id::Symbol) -> Term

Look up a controlled-vocabulary term, e.g. `term(:Phosphorylation)`."""
term(id::Symbol) = haskey(VOCABULARY, id) ? VOCABULARY[id] :
    throw(KeyError("unknown BioNExt vocabulary term :$id"))

"""    is_term(id::Symbol) -> Bool"""
is_term(id::Symbol) = haskey(VOCABULARY, id)

"""    is_a(id::Symbol, ancestor::Symbol) -> Bool

True if term `id` equals `ancestor` or is a (transitive) subclass of it."""
function is_a(id::Symbol, ancestor::Symbol)
    cur = id
    while cur !== nothing
        cur === ancestor && return true
        haskey(VOCABULARY, cur) || return false
        cur = VOCABULARY[cur].parent
    end
    return false
end

"""    terms(root::Symbol) -> Vector{Symbol}

All terms below `root` (excluding `root` itself), sorted, e.g. `terms(:RelationshipType)`."""
terms(root::Symbol) = sort!([k for k in keys(VOCABULARY) if k !== root && is_a(k, root)])

Base.show(io::IO, t::Term) = print(io, "Term(:", t.id, ", \"", t.label, "\")")
