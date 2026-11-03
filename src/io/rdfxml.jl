# RDF/XML writer (one rdf:Description per subject).

xml_escape(s::AbstractString) =
    replace(s, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;", '\'' => "&apos;")

const NCNAME = r"^[A-Za-z_][A-Za-z0-9_.\-]*$"

function xml_qname(iri::String, pfx)
    best = ""
    for (p, ns) in pfx
        startswith(iri, ns) && occursin(NCNAME, iri[length(ns)+1:end]) &&
            (isempty(best) || length(ns) > length(pfx[best])) && (best = p)
    end
    isempty(best) && throw(ArgumentError("cannot write predicate <$iri> as an XML qualified name"))
    return best * ":" * iri[length(pfx[best])+1:end]
end

node_attr(iri::String, attr) =
    startswith(iri, "_:") ? "rdf:nodeID=\"$(xml_escape(iri[3:end]))\"" : "rdf:$attr=\"$(xml_escape(iri))\""

"""
    write_rdfxml(io_or_path, kg; inverses = true)

Serialize the graph as RDF/XML (readable by Protégé, rdflib, Jena, ...).
"""
function write_rdfxml(io::IO, kg::KnowledgeGraph; inverses::Bool = true)
    pfx = prefixes(kg)
    println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    print(io, "<rdf:RDF")
    for p in sort!(collect(keys(pfx)))
        print(io, "\n    xmlns:", p, "=\"", xml_escape(pfx[p]), "\"")
    end
    println(io, ">")
    order, groups = IRI[], Dict{IRI,Vector{Triple}}()
    for t in unique(triples(kg; inverses))
        haskey(groups, t.subject) || (push!(order, t.subject); groups[t.subject] = Triple[])
        push!(groups[t.subject], t)
    end
    for s in order
        println(io, "  <rdf:Description ", node_attr(s.value, "about"), ">")
        for t in groups[s]
            q = xml_qname(t.predicate.value, pfx)
            o = t.object
            if o isa IRI
                println(io, "    <", q, " ", node_attr(o.value, "resource"), "/>")
            else
                attr = o.lang !== nothing ? " xml:lang=\"$(o.lang)\"" :
                       o.datatype !== nothing ? " rdf:datatype=\"$(xml_escape(o.datatype))\"" : ""
                println(io, "    <", q, attr, ">", xml_escape(o.value), "</", q, ">")
            end
        end
        println(io, "  </rdf:Description>")
    end
    println(io, "</rdf:RDF>")
end
