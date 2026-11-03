"""
    BioNExt

Build, query, validate and serialize knowledge graphs structured by the
BioNExt ontology (https://w3id.org/bionext/ontology/).
"""
module BioNExt

using Dates

include("vocabulary.jl")
include("types.jl")
include("schema.jl")
include("graph.jl")
include("validate.jl")
include("rdf.jl")
include("resolve.jl")
include("algorithms.jl")
include("io/turtle.jl")
include("io/json.jl")
include("io/rdfxml.jl")
include("io/tables.jl")
include("io/network.jl")
include("io/io.jl")
include("audit.jl")
include("stats.jl")

export
    # vocabulary
    Term, term, terms, is_a, is_term, VOCABULARY,
    # classes
    Resource, Entity, NodeEntity, EdgeEntity, Provenance, Taxon,
    Gene, Dna, Rna, Protein, MolecularComplex, SmallMolecule, OrthologGroup,
    Interaction, Influence, Reaction, PathwayModule, class_name,
    # schema
    Prop, schema,
    # graph
    KnowledgeGraph, add!, resources, nodes, edges, pathways,
    edges_of, neighbors, pathways_of, components, complexes_of, variants,
    ortholog_groups_of, references, referrers,
    # indexing
    index!, reindex!, unindex!, isindexed,
    # algorithms
    degree_of, shortest_path, network_components, subgraph, to_graph,
    # entity resolution
    duplicates, resolve!, integrate!, find_entity, identifier_map, normalize_identifier,
    Resolution, Conflict, keep_first,
    # validation
    Issue, validate,
    # RDF
    IRI, Literal, Triple, triples, expand, compact,
    write_turtle, write_ntriples, parse_ntriples, read_ntriples, from_triples,
    parse_turtle, read_turtle, write_rdfxml,
    write_json, read_json, write_jsonld, read_jsonld, parse_jsonld, parse_json,
    write_table, read_table, assign!,
    Arc, network, write_sif, read_sif, write_graphml, write_dot,
    # any format
    load, load!, save, formats,
    # data quality
    audit, AuditReport, Finding, AUDIT_CHECKS, ID_PATTERNS, findings, passed, count_findings,
    assert_quality, QualityError,
    # statistics
    graph_stats, GraphStats, stats_table, write_stats

end
