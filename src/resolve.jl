# Entity resolution: find resources that describe the same thing and merge them.
#
# Two resources are the same entity when they have the same class and share an
# identifier: their ids and `external_ids` (bn:hasExternalId, which the ontology defines
# as an accepted equivalent identifier). `db_xrefs` (bn:hasDbXref) do not assert
# equivalence and are only used with `use_xrefs = true`.

"""
    Conflict(id, field, kept, discarded, source)

Two merged resources had different values for a single-valued field. `kept` is the value
the canonical resource `id` ended up with; `discarded` came from (or stayed with) `source`.
"""
struct Conflict
    id::String
    field::Symbol
    kept::Any
    discarded::Any
    source::String
end

Base.show(io::IO, c::Conflict) =
    print(io, "Conflict(", c.id, ".", c.field, ": kept ", repr(c.kept), ", discarded ", repr(c.discarded),
          " from ", c.source, ")")

"""
    Resolution

Result of `resolve!` / `integrate!`:
- `mapping`: merged id ⇒ canonical id (canonical ids are not keys)
- `groups`: each merged group, canonical id first
- `conflicts`: single-valued fields whose values disagreed
"""
struct Resolution
    mapping::Dict{String,String}
    groups::Vector{Vector{String}}
    conflicts::Vector{Conflict}
end

function Base.show(io::IO, r::Resolution)
    print(io, "Resolution: ", length(r.mapping), " resources merged into ", length(r.groups),
          ", ", length(r.conflicts), " conflict", length(r.conflicts) == 1 ? "" : "s")
end

"""
    normalize_identifier(kg, id; prefix_aliases = Dict()) -> String

Canonical comparison key for an identifier: its full IRI, with the identifiers.org
prefix lower-cased and renamed through `prefix_aliases` (e.g. `Dict("uniprotkb" => "uniprot")`),
so `"UniProt:P04637"`, `"uniprot:P04637"` and `"https://identifiers.org/UniProt:P04637"` match.
"""
function normalize_identifier(kg::KnowledgeGraph, id::AbstractString; prefix_aliases = Dict{String,String}())
    iri = expand(kg, id)
    startswith(iri, IDENTIFIERS) || return iri
    m = match(r"^([^:]+):(.*)$", iri[length(IDENTIFIERS)+1:end])
    m === nothing && return iri
    prefix = lowercase(m[1])
    prefix = lowercase(get(Dict(lowercase(k) => v for (k, v) in prefix_aliases), prefix, prefix))
    return IDENTIFIERS * prefix * ":" * m[2]
end

identifiers_of(r::Resource; use_xrefs::Bool = false) =
    r isa Entity ? [r.id; r.external_ids; use_xrefs ? r.db_xrefs : String[]] : [r.id]

"""
    identifier_map(kg; use_xrefs = false, prefix_aliases = Dict()) -> Dict{String,String}

Normalized identifier ⇒ resource id, covering every id and external id. Build it once
for many lookups; `find_entity` does a single lookup.
"""
function identifier_map(kg::KnowledgeGraph; use_xrefs::Bool = false, prefix_aliases = Dict{String,String}())
    out = Dict{String,String}()
    for r in kg, x in identifiers_of(r; use_xrefs)
        get!(out, normalize_identifier(kg, x; prefix_aliases), r.id)
    end
    return out
end

"""
    find_entity(kg, identifier; use_xrefs = false, prefix_aliases = Dict()) -> Union{Resource,Nothing}

The resource whose id or `external_ids` match `identifier` (after normalization), e.g.
`find_entity(kg, "Ensembl:ENSP00000269305")` finds the protein stored as `UniProt:P04637`
if that Ensembl id is among its external ids. Returns the first match.
"""
function find_entity(kg::KnowledgeGraph, x::AbstractString; use_xrefs::Bool = false, prefix_aliases = Dict{String,String}())
    haskey(kg, x) && return kg[x]
    key = normalize_identifier(kg, x; prefix_aliases)
    for r in kg, y in identifiers_of(r; use_xrefs)
        normalize_identifier(kg, y; prefix_aliases) == key && return r
    end
    return nothing
end

"""
    duplicates(kg; use_xrefs = false, prefix_aliases = Dict()) -> Vector{Vector{String}}

Groups of resource ids that describe the same entity: same class and at least one shared
identifier (id or external id, transitively). Only groups of two or more are returned,
each in insertion order.
"""
function duplicates(kg::KnowledgeGraph; use_xrefs::Bool = false, prefix_aliases = Dict{String,String}())
    parent = Dict{String,String}()
    find(x) = (p = get!(parent, x, x); p == x ? x : (parent[x] = find(p)))
    owner = Dict{Tuple{Type,String},String}()      # (class, identifier key) => first resource id
    for r in kg
        r isa Entity || continue
        find(r.id)
        for x in identifiers_of(r; use_xrefs)
            key = (typeof(r), normalize_identifier(kg, x; prefix_aliases))
            first_ = get!(owner, key, r.id)
            a, b = find(first_), find(r.id)
            a != b && (kg.seq[a] < kg.seq[b] ? (parent[b] = a) : (parent[a] = b))
        end
    end
    groups = Dict{String,Vector{String}}()
    foreach(id -> haskey(parent, id) && push!(get!(groups, find(id), String[]), id), keys(kg))
    return sort!([g for g in values(groups) if length(g) > 1]; by = g -> kg.seq[g[1]])
end

keep_first(field, kept, other) = kept

# Pick the canonical id of a group: by prefix preference, then insertion order.
function choose_canonical(kg, group, prefer)
    prefer isa Function && return prefer(group)
    ranks = Dict(lowercase(p) => i for (i, p) in enumerate(prefer))
    rank(id) = (m = match(CURIE, id); m === nothing ? length(ranks) + 1 : get(ranks, lowercase(m[1]), length(ranks) + 1))
    return argmin(id -> (rank(id), kg.seq[id]), group)
end

# Fold `other` into `canon`: union vector fields, fill empty scalars, record conflicts.
function merge_into!(canon::Resource, other::Resource, conflicts, on_conflict)
    typeof(canon) == typeof(other) ||
        throw(ArgumentError("cannot merge $(other.id) ($(class_name(other))) into $(canon.id) ($(class_name(canon)))"))
    for f in fieldnames(typeof(canon))
        f === :id && continue
        a, b = getfield(canon, f), getfield(other, f)
        if a isa AbstractVector
            foreach(x -> x in a || push!(a, x), b)
        elseif a === nothing || (f === :type && a === :Unknown)
            setfield!(canon, f, b)
        elseif b !== nothing && a != b
            v = on_conflict(f, a, b)
            v == a || setfield!(canon, f, v)
            push!(conflicts, Conflict(canon.id, f, v, v == a ? b : a, other.id))
        end
    end
    canon isa Entity && other.id != canon.id && !(other.id in canon.external_ids) && push!(canon.external_ids, other.id)
    return canon
end

# Point every reference at canonical ids; drop duplicates and self-references.
function rewrite_references!(kg, mapping)
    isempty(mapping) && return
    for r in kg, p in schema(r)
        p.kind === :ref || continue
        v = getfield(r, p.field)
        if v isa AbstractVector
            new = unique!([get(mapping, x, x) for x in v])
            filter!(!=(r.id), new)
            setfield!(r, p.field, new)
        elseif v !== nothing
            setfield!(r, p.field, get(mapping, v, v) == r.id ? nothing : get(mapping, v, v))
        end
    end
end

# Merge each group into its canonical member; returns the set of ids that were merged away.
function merge_groups!(kg, groups, prefer, conflicts, on_conflict, mapping, done)
    merged = Set{String}()
    for g in groups
        canon = choose_canonical(kg, g, prefer)
        rest = filter(!=(canon), g)
        for id in rest
            merge_into!(kg[canon], kg[id], conflicts, on_conflict)
            mapping[id] = canon
            push!(merged, id)
        end
        push!(done, [canon; rest])
    end
    # a mapping that pointed at an id merged in this round follows it to its new canonical id
    for (k, v) in mapping
        while haskey(mapping, v)
            v = mapping[v]
        end
        mapping[k] = v
    end
    return merged
end

# Edges that are identical after node merging: same class, type and participant sets.
function duplicate_edges(kg)
    groups = Dict{Tuple{Type,Symbol,Set{String},Set{String}},Vector{String}}()
    for e in edges(kg)
        push!(get!(groups, (typeof(e), e.type, Set(e.left), Set(e.right)), String[]), e.id)
    end
    return sort!([g for g in values(groups) if length(g) > 1]; by = g -> kg.seq[g[1]])
end

"""
    resolve!(kg; prefer = String[], merge_edges = true, use_xrefs = false,
             prefix_aliases = Dict(), on_conflict = keep_first) -> Resolution

Merge every group found by `duplicates` into one canonical resource, then rewrite all
references in the graph to the canonical ids.

- **Canonical id**: the member whose id prefix comes first in `prefer`
  (e.g. `["UniProt", "HGNC", "Ensembl"]`, case-insensitive), else the earliest added.
  `prefer` may also be a function `group -> id`.
- **Fields**: vector fields are united; empty single-valued fields are filled in;
  differing values are settled by `on_conflict(field, kept, other)` (default: keep the
  canonical's value) and recorded as `Conflict`s.
- **Traceability**: merged ids are added to the canonical resource's `external_ids`, so
  `find_entity(kg, old_id)` still works.
- **Edges**: with `merge_edges = true`, edges that become identical (same class, type, and
  participant sets) are merged too, repeatedly, since merging edges can make edges that
  target them identical.
"""
function resolve!(kg::KnowledgeGraph; prefer = String[], merge_edges::Bool = true, use_xrefs::Bool = false,
                  prefix_aliases = Dict{String,String}(), on_conflict = keep_first)
    mapping, done, conflicts = Dict{String,String}(), Vector{String}[], Conflict[]
    merged = merge_groups!(kg, duplicates(kg; use_xrefs, prefix_aliases), prefer, conflicts, on_conflict, mapping, done)
    rewrite_references!(kg, mapping)
    delete_many!(kg, merged)
    while merge_edges
        groups = duplicate_edges(kg)
        isempty(groups) && break
        round = Dict{String,String}()
        merged = merge_groups!(kg, groups, prefer, conflicts, on_conflict, round, done)
        rewrite_references!(kg, round)
        delete_many!(kg, merged)
        merge!(mapping, round)
        for (k, v) in mapping
            mapping[k] = get(round, v, v)
        end
    end
    # tidy external ids of canonical resources: no self-references, no case/IRI-form repeats
    for c in unique(values(mapping))
        r = get(kg, c, nothing)
        r isa Entity || continue
        seen = Set([normalize_identifier(kg, c; prefix_aliases)])
        filter!(x -> (k = normalize_identifier(kg, x; prefix_aliases); k in seen ? false : (push!(seen, k); true)), r.external_ids)
    end
    reindex!(kg)
    # report each final canonical id once, with everything merged into it, in merge order
    final, order = Dict{String,Vector{String}}(), String[]
    for g in done, x in g[2:end]
        c = mapping[x]
        haskey(final, c) || (push!(order, c); final[c] = String[])
        x in final[c] || push!(final[c], x)
    end
    return Resolution(mapping, [[c; final[c]] for c in order], conflicts)
end

"""
    integrate!(kg, other; kwargs...) -> Resolution

Bring the resources of `other` (e.g. a graph loaded from another database) into `kg`
and resolve duplicates. Resources of `other` are copied; one with an id already in `kg`
is merged field by field into the existing one (same class required). Then `resolve!`
runs with `kwargs`, merging entities linked through `external_ids`. `other` is not modified.
"""
function integrate!(kg::KnowledgeGraph, other::KnowledgeGraph; on_conflict = keep_first, kwargs...)
    conflicts = Conflict[]
    merge!(kg.prefixes, other.prefixes)
    for r in other
        r = deepcopy(r)
        if haskey(kg, r.id)
            merge_into!(kg[r.id], r, conflicts, on_conflict)
            index_remove!(kg, kg[r.id]); index_add!(kg, kg[r.id])
        else
            add!(kg, r)
        end
    end
    res = resolve!(kg; on_conflict, kwargs...)
    return Resolution(res.mapping, res.groups, [conflicts; res.conflicts])
end
