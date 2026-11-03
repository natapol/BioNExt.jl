# Julia mirror of the BioNExt class hierarchy.
#
#   Resource
#   ├── Provenance                                   bn:Provenance
#   └── Entity                                       bn:Entity
#       ├── Taxon                                    bn:Taxon
#       ├── PathwayModule                            bn:PathwayModule
#       ├── NodeEntity                               bn:NodeEntity
#       │   ├── Gene, Dna, Rna, Protein, SmallMolecule, OrthologGroup
#       │   └── MolecularComplex                     bn:Complex (renamed: Base.Complex)
#       └── EdgeEntity                               bn:EdgeEntity
#           └── Interaction, Influence, Reaction
#
# References to other resources are stored as identifier strings and resolved
# through a KnowledgeGraph. Inverse properties (e.g. bn:hasVariant / bn:isVariantOf)
# are stored in one canonical direction only; the other is derived on query/export.

"""Supertype of everything a `KnowledgeGraph` stores: entities and provenance records."""
abstract type Resource end
"""Abstract `bn:Entity`: any biological or conceptual entity (nodes, edges, pathways, taxa)."""
abstract type Entity <: Resource end
"""Abstract `bn:NodeEntity`: a pool of molecules or complexes acting as a node in pathway models."""
abstract type NodeEntity <: Entity end
"""Abstract `bn:EdgeEntity`: a typed relationship with `left` and `right` participants."""
abstract type EdgeEntity <: Entity end

const OptString = Union{String,Nothing}
const OptFloat = Union{Float64,Nothing}
const OptInt = Union{Int,Nothing}

# Properties whose domain is bn:Entity.
const ENTITY_FIELDS = [
    :(name::OptString = nothing),                  # bn:name
    :(synonyms::Vector{String} = String[]),        # bn:synonyms
    :(organism::OptString = nothing),              # bn:organism -> Taxon
    :(db_xrefs::Vector{String} = String[]),        # bn:hasDbXref
    :(external_ids::Vector{String} = String[]),    # bn:hasExternalId
    :(provenance::Vector{String} = String[]),      # bn:hasProvenance -> Provenance
    :(publications::Vector{String} = String[]),    # bn:publicationEvidence
]

# Properties whose domain is bn:NodeEntity.
const NODE_FIELDS = [
    :(cellular_location::Vector{Symbol} = Symbol[]),  # bn:cellularLocation -> CellularLocation term
    :(hgvs::OptString = nothing),                     # bn:hasHGVS
]

# Properties whose domain is bn:EdgeEntity. `type` has no default: every edge needs one.
const EDGE_FIELDS = [
    :(type::Symbol),                                  # bn:hasType -> RelationshipType term
    :(left::Vector{String} = String[]),               # bn:hasLeft
    :(right::Vector{String} = String[]),              # bn:hasRight
    :(weight::OptFloat = nothing),                    # bn:hasWeight
    :(confidence::OptFloat = nothing),                # bn:hasConfidenceScore
]

# Turns `struct T <: Super ... end` into a mutable keyword struct whose fields are
# `id`, then the fields inherited from the ontology superclass, then its own.
macro bionext(ex)
    ex isa Expr && ex.head === :struct || error("@bionext expects a struct definition")
    ex.args[1] = true
    super = ex.args[2].args[2]
    inherited = super === :NodeEntity ? [ENTITY_FIELDS; NODE_FIELDS] :
                super === :EdgeEntity ? [ENTITY_FIELDS; EDGE_FIELDS] :
                super === :Entity     ? ENTITY_FIELDS : Expr[]
    own = filter(a -> !(a isa LineNumberNode), ex.args[3].args)
    ex.args[3] = Expr(:block, :(id::String), deepcopy(inherited)..., own...)
    return esc(:(Base.@kwdef $ex))
end

"""Organism at species or genus level (`bn:Taxon`). Use `name` for the scientific name."""
@bionext struct Taxon <: Entity end

"""Provenance metadata (`bn:Provenance`); attached to entities via their `provenance` field."""
@bionext struct Provenance <: Resource
    collected_date::Union{DateTime,Nothing} = nothing  # bn:collectedDate
    modified_date::Union{DateTime,Nothing} = nothing   # bn:modifiedDate
    created_by::Vector{String} = String[]              # bn:createdBy
    derived_from::Vector{String} = String[]            # bn:derivedFrom
    sources::Vector{String} = String[]                 # bn:hasSource
end

"""A gene (`bn:Gene`)."""
@bionext struct Gene <: NodeEntity
    functions::Vector{String} = String[]      # bn:hasFunction (e.g. GO term IRIs/CURIEs)
    templates::Vector{String} = String[]      # bn:hasTemplate -> Dna   (inverse: bn:annotates)
    regions::Vector{Symbol} = Symbol[]        # bn:hasRegion -> DnaRegion term
    variant_of::Vector{String} = String[]     # bn:isVariantOf -> Gene  (inverse: bn:hasVariant)
    variant_type::Vector{String} = String[]   # bn:variantType
end

"""A DNA element (`bn:Dna`)."""
@bionext struct Dna <: NodeEntity
    chromosome::OptString = nothing           # bn:chromosome
    sequence::OptString = nothing             # bn:hasSeq
    start_pos::OptInt = nothing               # bn:startPos
    end_pos::OptInt = nothing                 # bn:endPos
    length::OptInt = nothing                  # bn:length
    regions::Vector{Symbol} = Symbol[]        # bn:hasRegion -> DnaRegion term
    variant_of::Vector{String} = String[]     # bn:isVariantOf -> Dna
    variant_type::Vector{String} = String[]   # bn:variantType
end

"""An RNA element (`bn:Rna`)."""
@bionext struct Rna <: NodeEntity
    sequence::OptString = nothing             # bn:hasSeq
    start_pos::OptInt = nothing               # bn:startPos
    end_pos::OptInt = nothing                 # bn:endPos
    length::OptInt = nothing                  # bn:length
    regions::Vector{Symbol} = Symbol[]        # bn:hasRegion -> RnaRegion term
    encoded_by::Vector{String} = String[]     # bn:encodedBy -> Dna     (inverse: bn:encodes)
    modification_sites::Vector{String} = String[]  # bn:hasModificationSitedOn
    modification_type::Vector{String} = String[]   # bn:modificationType
    variant_of::Vector{String} = String[]     # bn:isVariantOf -> Rna
    variant_type::Vector{String} = String[]   # bn:variantType
end

"""A protein (`bn:Protein`)."""
@bionext struct Protein <: NodeEntity
    sequence::OptString = nothing             # bn:hasSeq
    cofactors::Vector{String} = String[]      # bn:hasCofactor -> SmallMolecule
    isoform_of::Vector{String} = String[]     # bn:isIsoformOf -> Protein
    modification_sites::Vector{String} = String[]  # bn:hasModificationSitedOn
    modification_type::Vector{String} = String[]   # bn:modificationType
    variant_of::Vector{String} = String[]     # bn:isVariantOf -> Protein
    variant_type::Vector{String} = String[]   # bn:variantType
end

"""A biomolecular complex (`bn:Complex`)."""
@bionext struct MolecularComplex <: NodeEntity
    subunits::Vector{String} = String[]       # bn:hasSubunit -> NodeEntity (inverse: bn:memberOfComplex)
    subcomplexes::Vector{String} = String[]   # bn:hasSubcomplex -> Complex
end

"""A small molecule (`bn:SmallMolecule`)."""
@bionext struct SmallMolecule <: NodeEntity
    cas_number::OptString = nothing           # bn:hasCasNum
    formula::OptString = nothing              # bn:hasChemFormula
    inchi::OptString = nothing                # bn:hasInChI
    inchikey::OptString = nothing             # bn:hasInChIKey
    net_charge::OptInt = nothing              # bn:hasNetCharge
    smiles::OptString = nothing               # bn:hasSmiles
end

"""A protein ortholog group (`bn:OrthologGroup`)."""
@bionext struct OrthologGroup <: NodeEntity
    clustering_method::OptString = nothing    # bn:clusteringMethod
    confidence::OptFloat = nothing            # bn:hasConfidenceScore
    members::Vector{String} = String[]        # bn:hasMember / bn:hasProteinMember -> Gene|Protein
    species::Vector{String} = String[]        # bn:hasSpeciesCoverage -> Taxon
end

"""A physical/functional interaction between node entities (`bn:Interaction`)."""
@bionext struct Interaction <: EdgeEntity end

"""A directed influence whose target may be a node, an edge or a pathway (`bn:Influence`)."""
@bionext struct Influence <: EdgeEntity end

"""A biochemical reaction transforming left (substrates) into right (products) (`bn:Reaction`)."""
@bionext struct Reaction <: EdgeEntity
    ec_numbers::Vector{String} = String[]     # bn:hasEcNumber
end

"""A pathway module grouping nodes, edges and sub-pathways (`bn:PathwayModule`)."""
@bionext struct PathwayModule <: Entity
    components::Vector{String} = String[]     # bn:pathwayComponent (inverse: bn:participatesInPathway)
end

# rdf:type local names in the bn: namespace.
const CLASS_NAMES = Dict{Type,String}(
    Entity => "Entity", NodeEntity => "NodeEntity", EdgeEntity => "EdgeEntity",
    Taxon => "Taxon", Provenance => "Provenance", Gene => "Gene", Dna => "Dna",
    Rna => "Rna", Protein => "Protein", MolecularComplex => "Complex",
    SmallMolecule => "SmallMolecule", OrthologGroup => "OrthologGroup",
    Interaction => "Interaction", Influence => "Influence", Reaction => "Reaction",
    PathwayModule => "PathwayModule",
)
const CLASS_TYPES = Dict{String,Type}(v => k for (k, v) in CLASS_NAMES)

"""    class_name(x) -> String

Local name of the ontology class of a resource or type (e.g. `"Complex"` for `MolecularComplex`)."""
class_name(::Type{T}) where {T} = CLASS_NAMES[T]
class_name(r::Resource) = class_name(typeof(r))
