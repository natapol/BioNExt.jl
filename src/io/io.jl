# Format-agnostic load/save, dispatching on file extension.

const EXTENSIONS = Dict(
    ".ttl" => :turtle, ".turtle" => :turtle, ".nt" => :ntriples,
    ".jsonld" => :jsonld, ".json" => :json,
    ".rdf" => :rdfxml, ".owl" => :rdfxml, ".xml" => :rdfxml,
    ".csv" => :csv, ".tsv" => :tsv, ".tab" => :tsv, ".txt" => :tsv,
    ".sif" => :sif, ".graphml" => :graphml, ".dot" => :dot, ".gv" => :dot,
)

const READERS = Dict{Symbol,Function}(
    :turtle => read_turtle, :ntriples => read_ntriples, :jsonld => read_jsonld, :json => read_json,
    :csv => (io; kw...) -> read_table(io; delim = ',', kw...),
    :tsv => (io; kw...) -> read_table(io; delim = '\t', kw...),
    :sif => read_sif,
)

const WRITERS = Dict{Symbol,Function}(
    :turtle => write_turtle, :ntriples => write_ntriples, :jsonld => write_jsonld, :json => write_json,
    :rdfxml => write_rdfxml,
    :csv => (io, kg; kw...) -> write_table(io, kg; delim = ',', kw...),
    :tsv => (io, kg; kw...) -> write_table(io, kg; delim = '\t', kw...),
    :sif => write_sif, :graphml => write_graphml, :dot => write_dot,
)

"""    formats() -> Vector of (format, extensions, read, write)

Supported file formats."""
formats() = [(format = f,
              extensions = sort!([e for (e, g) in EXTENSIONS if g === f]),
              read = haskey(READERS, f), write = haskey(WRITERS, f))
             for f in sort!(unique(collect(values(EXTENSIONS))))]

function detect_format(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    haskey(EXTENSIONS, ext) && return EXTENSIONS[ext]
    throw(ArgumentError("cannot infer format from \"$path\"; pass format = one of $(join(sort!(unique(collect(values(EXTENSIONS)))), ", "))"))
end

"""
    load(path; format = <from extension>, kwargs...) -> KnowledgeGraph

Read a graph from any supported format (see `formats()`). Keyword arguments go to the
format's reader, e.g. `load("ppi.sif"; relationships = Dict("pp" => :Binding))` or
`load("genes.csv"; type = Gene, mapping = Dict("symbol" => :name))`.
"""
function load(path::AbstractString; format::Symbol = detect_format(path), kwargs...)
    haskey(READERS, format) || throw(ArgumentError("reading $format is not supported (write-only)"))
    return open(io -> READERS[format](io; kwargs...), path)
end

"""
    load!(kg, path; replace = false, kwargs...) -> kg

Load a file and merge it into `kg`. Tables and SIF are read straight into `kg`, so
their rows can refer to (and update) resources already in it.
"""
function load!(kg::KnowledgeGraph, path::AbstractString; replace::Bool = false,
               format::Symbol = detect_format(path), kwargs...)
    if format in (:csv, :tsv, :sif)
        return load(path; format, into = kg, kwargs...)
    end
    return merge!(kg, load(path; format, kwargs...); replace)
end

"""
    save(path, kg; format = <from extension>, kwargs...) -> path

Write a graph in any supported format (see `formats()`); keyword arguments go to the
writer, e.g. `save("edges.tsv", kg; types = EdgeEntity)` or `save("net.graphml", kg; reify = true)`.
"""
function save(path::AbstractString, kg::KnowledgeGraph; format::Symbol = detect_format(path), kwargs...)
    haskey(WRITERS, format) || throw(ArgumentError("writing $format is not supported"))
    open(io -> WRITERS[format](io, kg; kwargs...), path, "w")
    return path
end

# Path-accepting methods for every reader/writer.
for f in (:write_jsonld, :write_json, :write_rdfxml, :write_table, :write_sif, :write_graphml, :write_dot)
    @eval $f(path::AbstractString, kg::KnowledgeGraph; kw...) = open(io -> $f(io, kg; kw...), path, "w")
end
read_table(path::AbstractString; delim::Char = lowercase(splitext(path)[2]) == ".csv" ? ',' : '\t', kw...) =
    open(io -> read_table(io; delim, kw...), path)
read_sif(path::AbstractString; kw...) = open(io -> read_sif(io; kw...), path)
