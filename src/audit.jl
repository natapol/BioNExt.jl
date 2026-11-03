# Data-quality audit: checks that go beyond the ontology rules of `validate` — well-formed
# identifiers and sequences, consistent annotations, duplicates, structural problems,
# and annotation gaps.

"""
    Finding(check, severity, id, message)

One audit result. `severity` is `:error` (the data are wrong), `:warning` (probably
wrong or needs attention) or `:info` (worth knowing, often fine).
"""
struct Finding
    check::Symbol
    severity::Symbol
    id::String
    message::String
end

Base.show(io::IO, f::Finding) = print(io, uppercase(string(f.severity)), " [", f.check, "] ", f.id, ": ", f.message)

"""
    AuditReport

Result of `audit`: `findings` (at most `max_per_check` per check), `counts` (all findings,
by check and severity, including those not kept) and `checks` (the checks that ran).
Use `findings(report; check, severity)`, `passed(report)` and `show`.
"""
struct AuditReport
    findings::Vector{Finding}
    counts::Dict{Tuple{Symbol,Symbol},Int}
    checks::Vector{Symbol}
    resources::Int
end

# --- identifier and value formats --------------------------------------------------------

"""
Patterns for the local part of CURIEs, by prefix. `audit` reports ids with one of these
prefixes that do not match. Add entries for other databases as needed.
"""
const ID_PATTERNS = Dict{String,Regex}(
    # accession, optional isoform, optional variant/site suffix as in the ontology's own
    # examples ("UniProt:P04637_p.R273H", "UniProt:P31749_Ser473")
    "UniProt" => r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?(_\S+)?$",
    "CHEBI" => r"^\d+$", "NCBIGene" => r"^\d+$", "NCBITaxon" => r"^\d+$", "HGNC" => r"^\d+$",
    "GO" => r"^\d{7}$", "PMID" => r"^\d+$", "MI" => r"^\d{4}$", "MOD" => r"^\d{5}$", "SO" => r"^\d{7}$",
    "EC" => r"^\d+\.(\d+|-)\.(\d+|-)\.(n?\d+|-)$", "Reactome" => r"^R-[A-Z]{3}-\d+(\.\d+)?$",
    "IntAct" => r"^EBI-\d+$", "SIGNOR" => r"^SIGNOR-[A-Z]*\d+$", "ComplexPortal" => r"^CPX-\d+$",
    "Ensembl" => r"^ENS[A-Z]*[EFGPRT]\d{11}(\.\d+)?$", "RefSeq" => r"^[A-Z]{2}_\d+(\.\d+)?$",
    "DrugBank" => r"^DB\d{5}$", "pubchem.compound" => r"^\d+$", "STRING" => r"^\d+\.\S+$",
    "kegg.compound" => r"^C\d{5}$", "kegg.reaction" => r"^R\d{5}$", "kegg.orthology" => r"^K\d{5}$",
    "kegg.pathway" => r"^[a-z]{2,4}\d{5}$", "kegg.genes" => r"^[a-z]{3,4}:\S+$", "DOI" => r"^10\.\S+/\S+$",
)

const PROTEIN_SEQ = r"^[ACDEFGHIKLMNPQRSTVWYUOXBZJ*]+$"i
const DNA_SEQ = r"^[ACGTNRYSWKMBDHV]+$"i
const RNA_SEQ = r"^[ACGUNRYSWKMBDHV]+$"i
const ONE_LETTER = Dict("Ala" => 'A', "Arg" => 'R', "Asn" => 'N', "Asp" => 'D', "Cys" => 'C', "Glu" => 'E',
    "Gln" => 'Q', "Gly" => 'G', "His" => 'H', "Ile" => 'I', "Leu" => 'L', "Lys" => 'K', "Met" => 'M',
    "Phe" => 'F', "Pro" => 'P', "Ser" => 'S', "Thr" => 'T', "Trp" => 'W', "Tyr" => 'Y', "Val" => 'V',
    "Sec" => 'U', "Pyl" => 'O')

function cas_valid(s)
    m = match(r"^(\d{2,7})-(\d{2})-(\d)$", s)
    m === nothing && return false
    digits = reverse(m[1] * m[2])
    return sum(i * parse(Int, digits[i]) for i in eachindex(digits)) % 10 == parse(Int, m[3])
end

const ACTIVATING = (:Activation, :Increase, :RegulateExpression, :Expression, :Catalysis, :Enhancement)
const INHIBITING = (:Inhibition, :Decrease, :Repression, :Degradation, :Suppression)
effect_sign(t) = any(a -> is_term(t) && is_a(t, a), ACTIVATING) ? 1 : any(a -> is_term(t) && is_a(t, a), INHIBITING) ? -1 : 0

# --- the checks ---------------------------------------------------------------------------
# Each check is (kg, emit, opts) -> nothing, where emit(severity, id, message).

function check_ontology(kg, emit, opts)
    for i in validate(kg)
        emit(i.severity, i.id, i.message)
    end
end

function check_identifiers(kg, emit, opts)
    patterns = opts.patterns
    function test(owner, what, x)
        m = match(CURIE, x)
        m === nothing && return
        pat = get(patterns, m[1], nothing)
        pat === nothing || occursin(pat, m[2]) || emit(:warning, owner, "$what \"$x\" is not a well-formed $(m[1]) identifier")
    end
    for r in kg
        test(r.id, "id", r.id)
        for p in schema(r)
            p.kind in (:iri, :ref) || continue
            for v in field_values(r, p.field)
                test(r.id, string(p.field), v)
            end
        end
    end
end

function check_sequences(kg, emit, opts)
    for r in kg
        hasfield(typeof(r), :sequence) && r.sequence !== nothing || continue
        alphabet = r isa Protein ? PROTEIN_SEQ : r isa Rna ? RNA_SEQ : DNA_SEQ
        s = r.sequence
        isempty(s) && (emit(:error, r.id, "sequence is empty"); continue)
        occursin(alphabet, s) || emit(:error, r.id, "sequence contains characters that are not valid for $(class_name(r))" *
                                     " (e.g. '$(first(filter(c -> !occursin(alphabet, string(c)), s)))')")
        if hasfield(typeof(r), :length) && r.length !== nothing && r.length != length(s)
            emit(:warning, r.id, "length is $(r.length) but the sequence has $(length(s)) residues")
        end
    end
end

function check_coordinates(kg, emit, opts)
    for r in kg
        hasfield(typeof(r), :start_pos) || continue
        s, e, l = r.start_pos, r.end_pos, r.length
        s !== nothing && e !== nothing && e < s && emit(:error, r.id, "end position $e is before start position $s")
        s !== nothing && e !== nothing && l !== nothing && e >= s && l != e - s + 1 &&
            emit(:warning, r.id, "length $l does not match the coordinates $(s)–$(e) ($(e - s + 1))")
    end
end

function check_modification_sites(kg, emit, opts)
    for r in kg
        hasfield(typeof(r), :modification_sites) || continue
        for site in r.modification_sites
            m = match(r"^([A-Z][a-z]{2})(\d+)$", site)
            m === nothing && (emit(:info, r.id, "modification site \"$site\" is not in the form Ser15"); continue)
            pos = parse(Int, m[2])
            seq = r isa Protein ? r.sequence : nothing
            seq === nothing && continue
            if pos > length(seq) || pos < 1
                emit(:error, r.id, "modification site $site is beyond the sequence ($(length(seq)) residues)")
            elseif get(ONE_LETTER, m[1], '?') != uppercase(seq[pos])
                emit(:warning, r.id, "modification site $site: the sequence has $(seq[pos]) at position $pos")
            end
        end
    end
end

function check_values(kg, emit, opts)
    for r in kg
        for f in (:confidence, :weight)
            hasfield(typeof(r), f) && (v = getfield(r, f)) !== nothing || continue
            isfinite(v) || (emit(:error, r.id, "$f is $v"); continue)
            f === :confidence && !(0 <= v <= 1) && emit(:warning, r.id, "confidence $v is outside [0, 1]")
        end
        if r isa SmallMolecule
            r.inchikey === nothing || occursin(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$", r.inchikey) ||
                emit(:warning, r.id, "InChIKey \"$(r.inchikey)\" is malformed")
            r.cas_number === nothing || cas_valid(r.cas_number) ||
                emit(:warning, r.id, "CAS number \"$(r.cas_number)\" is malformed or fails its checksum")
            r.formula === nothing || occursin(r"^(\(?[A-Z][a-z]?\d*\)?\d*|[.·+\-]|\d)+$", r.formula) ||
                emit(:warning, r.id, "formula \"$(r.formula)\" does not look like a chemical formula")
        end
        if r isa NodeEntity && r.hgvs !== nothing
            occursin(r"^([A-Za-z0-9_.\-]+:)?[cgmnopr]\.\S+$", r.hgvs) ||
                emit(:warning, r.id, "HGVS \"$(r.hgvs)\" does not look like an HGVS expression (c./g./p./...)")
        end
    end
end

function check_duplicates(kg, emit, opts)
    for g in duplicates(kg)
        emit(:warning, g[1], "describes the same entity as $(join(g[2:end], ", ")) (run resolve!)")
    end
    for g in duplicate_edges(kg)
        emit(:warning, g[1], "is identical to $(join(g[2:end], ", ")) (same class, type and participants; " *
                             "often one record per experiment; resolve! merges them)")
    end
end

function check_edges(kg, emit, opts)
    effects = Dict{Tuple{Set{String},Set{String}},Vector{Tuple{String,Int}}}()
    for e in edges(kg)
        loops = intersect(e.left, e.right)
        isempty(loops) || emit(:info, e.id, "self-loop: $(join(loops, ", ")) on both sides (e.g. a homodimer)")
        orgs = unique(String[kg[x].organism for x in [e.left; e.right]
                             if get(kg, x, nothing) isa Entity && kg[x].organism !== nothing])
        length(orgs) > 1 && emit(:info, e.id, "participants come from different organisms: $(join(orgs, ", "))")
        s = effect_sign(e.type)
        s == 0 || push!(get!(effects, (Set(e.left), Set(e.right)), Tuple{String,Int}[]), (e.id, s))
    end
    for (_, es) in effects
        pos, neg = [id for (id, s) in es if s > 0], [id for (id, s) in es if s < 0]
        isempty(pos) || isempty(neg) ||
            emit(:info, pos[1], "activating edge(s) $(join(pos, ", ")) and inhibiting edge(s) $(join(neg, ", ")) " *
                                "connect the same participants (context-dependent or contradictory)")
    end
end

function check_structure(kg, emit, opts)
    for c in resources(kg, MolecularComplex)
        isempty(c.subunits) && isempty(c.subcomplexes) && emit(:warning, c.id, "complex has no subunits")
        length(c.subunits) == 1 && isempty(c.subcomplexes) && emit(:info, c.id, "complex has a single subunit (homo-oligomer?)")
    end
    for p in pathways(kg)
        isempty(p.components) && emit(:warning, p.id, "pathway has no components")
    end
    for g in resources(kg, OrthologGroup)
        isempty(g.members) && emit(:warning, g.id, "ortholog group has no members")
    end
    # cycles: a pathway that (indirectly) contains itself, a complex that contains itself
    for (T, field) in ((PathwayModule, :components), (MolecularComplex, :subcomplexes))
        for r in resources(kg, T)
            seen, stack = Set{String}(), copy(getfield(r, field))
            while !isempty(stack)
                x = pop!(stack)
                x == r.id && (emit(:warning, r.id, "$(class_name(r)) contains itself through $field"); break)
                x in seen && continue
                push!(seen, x)
                y = get(kg, x, nothing)
                y isa T && append!(stack, getfield(y, field))
            end
        end
    end
end

function check_connectivity(kg, emit, opts)
    used = Set{String}()
    for r in kg, p in schema(r)
        p.kind === :ref && union!(used, field_values(r, p.field))
    end
    for n in nodes(kg)
        n.id in used || emit(:info, n.id, "not connected to anything (no edge, complex, group or pathway)")
    end
    for r in kg
        (r isa Provenance || r isa Taxon) && !(r.id in used) && emit(:info, r.id, "$(class_name(r)) is not used by any resource")
    end
end

function check_annotation(kg, emit, opts)
    for r in kg
        r isa Entity || continue
        r.name === nothing && !(r isa EdgeEntity) && emit(:info, r.id, "has no name")
        r isa Union{Gene,Protein,Dna,Rna} && r.organism === nothing && emit(:info, r.id, "has no organism")
        isempty(r.provenance) && emit(:info, r.id, "has no provenance")
        r isa EdgeEntity && isempty(r.publications) && emit(:info, r.id, "has no publication evidence")
    end
end

"""
The checks `audit` can run, in order: `check => (description, function)`.
"""
const AUDIT_CHECKS = [
    :ontology => ("ontology rules (`validate`)", check_ontology),
    :identifiers => ("identifiers match their database's format (`ID_PATTERNS`)", check_identifiers),
    :sequences => ("sequences use the right alphabet and match `length`", check_sequences),
    :coordinates => ("start ≤ end, length matches the coordinates", check_coordinates),
    :modification_sites => ("modification sites lie in the sequence and match its residue", check_modification_sites),
    :values => ("scores in range; InChIKey, CAS, formula and HGVS well formed", check_values),
    :duplicates => ("duplicate entities and identical edges", check_duplicates),
    :edges => ("self-loops, cross-species edges, opposite effects on the same participants", check_edges),
    :structure => ("empty complexes, pathways and groups; containment cycles", check_structure),
    :connectivity => ("nodes connected to nothing; unused taxa and provenance", check_connectivity),
    :annotation => ("missing names, organisms, provenance and publications", check_annotation),
]

"""
    audit(kg; checks = all, skip = Symbol[], max_per_check = 1000, patterns = ID_PATTERNS) -> AuditReport

Check the data in `kg` for problems the ontology rules cannot see (see `AUDIT_CHECKS`):

| Check | Looks for | Severities |
|---|---|---|
| `:ontology` | everything `validate` reports | error, warning |
| `:identifiers` | ids and cross-references that do not match their database's format | warning |
| `:sequences` | invalid residues for the molecule type; `length` ≠ sequence length | error, warning |
| `:coordinates` | end before start; length inconsistent with coordinates | error, warning |
| `:modification_sites` | sites beyond the sequence or on a different residue | error, warning |
| `:values` | NaN/Inf scores, confidence outside [0,1], malformed InChIKey/CAS/formula/HGVS | error, warning |
| `:duplicates` | entities that `resolve!` would merge; identical edges | warning |
| `:edges` | self-loops, cross-species edges, activating and inhibiting edges between the same participants | info |
| `:structure` | empty complexes/pathways/ortholog groups; containment cycles | warning, info |
| `:connectivity` | isolated nodes; unused taxa and provenance records | info |
| `:annotation` | entities without name, organism, provenance; edges without publications | info |

`checks` selects checks, `skip` leaves some out. At most `max_per_check` findings are
kept per check; `counts` always has the full numbers.
"""
function audit(kg::KnowledgeGraph; checks = first.(AUDIT_CHECKS), skip = Symbol[], max_per_check::Integer = 1000,
               patterns = ID_PATTERNS)
    selected = [c for c in first.(AUDIT_CHECKS) if c in checks && !(c in skip)]
    unknown = setdiff(collect(checks), first.(AUDIT_CHECKS))
    isempty(unknown) || throw(ArgumentError("unknown audit checks: $(join(unknown, ", "))"))
    out, counts = Finding[], Dict{Tuple{Symbol,Symbol},Int}()
    opts = (; patterns)
    for (name, (_, f)) in AUDIT_CHECKS
        name in selected || continue
        kept = Ref(0)
        emit(sev, id, msg) = begin
            counts[(name, sev)] = get(counts, (name, sev), 0) + 1
            kept[] < max_per_check && (push!(out, Finding(name, sev, id, msg)); kept[] += 1)
        end
        f(kg, emit, opts)
    end
    return AuditReport(out, counts, selected, length(kg))
end

"""    findings(report; check = nothing, severity = nothing) -> Vector{Finding}"""
findings(r::AuditReport; check = nothing, severity = nothing) =
    filter(f -> (check === nothing || f.check === check) && (severity === nothing || f.severity === severity), r.findings)

"""
    passed(report; fail_on = :error) -> Bool

`true` when the report has no findings of severity `fail_on` or worse
(`:error`; `:warning` = errors or warnings; `:info` = anything).
"""
function passed(r::AuditReport; fail_on::Symbol = :error)
    levels = fail_on === :error ? (:error,) : fail_on === :warning ? (:error, :warning) : (:error, :warning, :info)
    return all(p -> !(p.first[2] in levels) || p.second == 0, r.counts)
end

"""    count_findings(report; check = nothing, severity = nothing) -> Int

Number of findings (including those not kept because of `max_per_check`)."""
count_findings(r::AuditReport; check = nothing, severity = nothing) =
    sum((n for ((c, s), n) in r.counts if (check === nothing || c === check) && (severity === nothing || s === severity)); init = 0)

struct QualityError <: Exception
    report::AuditReport
    fail_on::Symbol
end
Base.showerror(io::IO, e::QualityError) =
    (print(io, "QualityError: the graph failed its audit (fail_on = :$(e.fail_on))\n"); show(io, MIME"text/plain"(), e.report))

"""
    assert_quality(kg; fail_on = :error, kwargs...) -> AuditReport

Run `audit` (with `kwargs`) and throw a `QualityError` showing the report if it does not
pass at level `fail_on`; otherwise return the report. Use it as a gate in pipelines and tests.
"""
function assert_quality(kg::KnowledgeGraph; fail_on::Symbol = :error, kwargs...)
    r = audit(kg; kwargs...)
    passed(r; fail_on) || throw(QualityError(r, fail_on))
    return r
end

function Base.show(io::IO, ::MIME"text/plain", r::AuditReport)
    nerr, nwarn, ninfo = (count_findings(r; severity = s) for s in (:error, :warning, :info))
    status = nerr > 0 ? "FAILED" : nwarn > 0 ? "passed with warnings" : "passed"
    println(io, "Audit of $(r.resources) resources: $status — $nerr errors, $nwarn warnings, $ninfo notes")
    for c in r.checks
        n = [(s, get(r.counts, (c, s), 0)) for s in (:error, :warning, :info) if get(r.counts, (c, s), 0) > 0]
        print(io, "  ", rpad(string(c), 20), isempty(n) ? "ok" : join(("$k $s" * (k == 1 ? "" : "s") for (s, k) in n), ", "))
        ex = [f for f in r.findings if f.check === c && f.severity !== :info]
        isempty(ex) && (ex = [f for f in r.findings if f.check === c])
        isempty(ex) || print(io, "\n      e.g. ", ex[1].id, ": ", ex[1].message)
        println(io)
    end
end
Base.show(io::IO, r::AuditReport) =
    print(io, "AuditReport(", count_findings(r; severity = :error), " errors, ", count_findings(r; severity = :warning),
          " warnings, ", count_findings(r; severity = :info), " notes)")
