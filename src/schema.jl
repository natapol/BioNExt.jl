# One table maps every struct field to its ontology predicate. Validation, RDF
# export and RDF import are all driven from it.

"""
    Prop(field, pred, kind, range)

Mapping of a struct field to a `bn:` predicate. `kind` is one of
- `:ref`     – identifier of another resource in the graph; `range` is a tuple of allowed types
- `:term`    – controlled-vocabulary `Symbol`; `range` is the root term it must descend from
- `:iri`     – external IRI or CURIE (not resolved in the graph); `range` is `nothing`
- `:literal` – data value; `range` is the XSD datatype local name
"""
struct Prop
    field::Symbol
    pred::String
    kind::Symbol
    range::Any
end

const PARTICIPANT = (EdgeEntity, NodeEntity, PathwayModule)

const ENTITY_PROPS = [
    Prop(:name, "name", :literal, "string"),
    Prop(:synonyms, "synonyms", :literal, "string"),
    Prop(:organism, "organism", :ref, (Taxon,)),
    Prop(:db_xrefs, "hasDbXref", :iri, nothing),
    Prop(:external_ids, "hasExternalId", :iri, nothing),
    Prop(:provenance, "hasProvenance", :ref, (Provenance,)),
    Prop(:publications, "publicationEvidence", :iri, nothing),
]

const NODE_PROPS = [
    ENTITY_PROPS;
    Prop(:cellular_location, "cellularLocation", :term, :CellularLocation);
    Prop(:hgvs, "hasHGVS", :literal, "string")
]

edge_props(right_range) = [
    ENTITY_PROPS;
    Prop(:type, "hasType", :term, :RelationshipType);
    Prop(:left, "hasLeft", :ref, PARTICIPANT);
    Prop(:right, "hasRight", :ref, right_range);
    Prop(:weight, "hasWeight", :literal, "float");
    Prop(:confidence, "hasConfidenceScore", :literal, "float")
]

variant_props(T) = [
    Prop(:variant_of, "isVariantOf", :ref, (T,)),
    Prop(:variant_type, "variantType", :iri, nothing),
]

coordinate_props() = [
    Prop(:sequence, "hasSeq", :literal, "string"),
    Prop(:start_pos, "startPos", :literal, "positiveInteger"),
    Prop(:end_pos, "endPos", :literal, "positiveInteger"),
    Prop(:length, "length", :literal, "positiveInteger"),
]

modification_props() = [
    Prop(:modification_sites, "hasModificationSitedOn", :literal, "string"),
    Prop(:modification_type, "modificationType", :iri, nothing),
]

const SCHEMA = Dict{Type,Vector{Prop}}(
    Taxon => ENTITY_PROPS,
    Provenance => [
        Prop(:collected_date, "collectedDate", :literal, "dateTime"),
        Prop(:modified_date, "modifiedDate", :literal, "dateTime"),
        Prop(:created_by, "createdBy", :iri, nothing),
        Prop(:derived_from, "derivedFrom", :iri, nothing),
        Prop(:sources, "hasSource", :iri, nothing),
    ],
    Gene => [
        NODE_PROPS;
        Prop(:functions, "hasFunction", :iri, nothing);
        Prop(:templates, "hasTemplate", :ref, (Dna,));
        Prop(:regions, "hasRegion", :term, :DnaRegion);
        variant_props(Gene)
    ],
    Dna => [
        NODE_PROPS;
        Prop(:chromosome, "chromosome", :literal, "string");
        coordinate_props();
        Prop(:regions, "hasRegion", :term, :DnaRegion);
        variant_props(Dna)
    ],
    Rna => [
        NODE_PROPS;
        coordinate_props();
        Prop(:regions, "hasRegion", :term, :RnaRegion);
        Prop(:encoded_by, "encodedBy", :ref, (Dna,));
        modification_props();
        variant_props(Rna)
    ],
    Protein => [
        NODE_PROPS;
        Prop(:sequence, "hasSeq", :literal, "string");
        Prop(:cofactors, "hasCofactor", :ref, (SmallMolecule,));
        Prop(:isoform_of, "isIsoformOf", :ref, (Protein,));
        modification_props();
        variant_props(Protein)
    ],
    MolecularComplex => [
        NODE_PROPS;
        Prop(:subunits, "hasSubunit", :ref, (NodeEntity,));
        Prop(:subcomplexes, "hasSubcomplex", :ref, (MolecularComplex,))
    ],
    SmallMolecule => [
        NODE_PROPS;
        Prop(:cas_number, "hasCasNum", :literal, "string");
        Prop(:formula, "hasChemFormula", :literal, "string");
        Prop(:inchi, "hasInChI", :literal, "string");
        Prop(:inchikey, "hasInChIKey", :literal, "string");
        Prop(:net_charge, "hasNetCharge", :literal, "integer");
        Prop(:smiles, "hasSmiles", :literal, "string")
    ],
    OrthologGroup => [
        NODE_PROPS;
        Prop(:clustering_method, "clusteringMethod", :literal, "string");
        Prop(:confidence, "hasConfidenceScore", :literal, "float");
        Prop(:members, "hasMember", :ref, (Gene, Protein));
        Prop(:species, "hasSpeciesCoverage", :ref, (Taxon,))
    ],
    Interaction => edge_props((NodeEntity,)),
    Reaction => [edge_props((NodeEntity,)); Prop(:ec_numbers, "hasEcNumber", :iri, nothing)],
    Influence => edge_props(PARTICIPANT),
    PathwayModule => [ENTITY_PROPS; Prop(:components, "pathwayComponent", :ref, PARTICIPANT)],
)

"""    schema(T) -> Vector{Prop}

Field-to-predicate mapping for a resource type."""
schema(::Type{T}) where {T<:Resource} = SCHEMA[T]
schema(r::Resource) = schema(typeof(r))

# Inverse predicate => canonical predicate stored on the struct.
const INVERSES = Dict(
    "hasVariant" => "isVariantOf",
    "encodes" => "encodedBy",
    "annotates" => "hasTemplate",
    "participatesInPathway" => "pathwayComponent",
    "memberOfComplex" => "hasSubunit",
    "hasOrthologGroup" => "hasProteinMember",
)
# Sub-property => the property it is stored under.
const SUBPROPERTIES = Dict("hasProteinMember" => "hasMember")

# Values of a field as a vector (scalars wrapped, `nothing` dropped).
field_values(r, f::Symbol) = (v = getfield(r, f); v isa AbstractVector ? v : v === nothing ? () : (v,))
