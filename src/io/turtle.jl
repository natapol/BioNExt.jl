# Turtle (RDF 1.1) reader: prefixes/base (both @ and SPARQL style), prefixed
# names, `a`, `;`/`,` lists, blank nodes (`_:x`, `[]`, `[ ... ]`), collections,
# short/long string literals with language tags or datatypes, numbers and booleans.

const XSD_INTEGER = XSD * "integer"
const XSD_DECIMAL = XSD * "decimal"
const XSD_DOUBLE = XSD * "double"
const XSD_BOOLEAN = XSD * "boolean"

mutable struct TurtleParser
    s::String
    i::Int
    base::String
    prefixes::Dict{String,String}
    bnodes::Int
    triples::Vector{Triple}
end

turtle_error(p::TurtleParser, msg) =
    throw(ArgumentError("Turtle parse error at line $(count(==('\n'), SubString(p.s, 1, min(p.i, ncodeunits(p.s)))) + 1): $msg"))

eof(p::TurtleParser) = p.i > ncodeunits(p.s)
peekc(p::TurtleParser) = eof(p) ? '\0' : p.s[p.i]
advance!(p::TurtleParser, n = 1) = (p.i = nextind(p.s, p.i, n))
lookahead(p::TurtleParser, str) = startswith(SubString(p.s, p.i), str)

function skipws!(p::TurtleParser)
    while !eof(p)
        c = peekc(p)
        if isspace(c)
            advance!(p)
        elseif c == '#'
            j = findnext('\n', p.s, p.i)
            p.i = j === nothing ? ncodeunits(p.s) + 1 : j
        else
            break
        end
    end
end

function expect!(p::TurtleParser, c::Char)
    skipws!(p)
    peekc(p) == c || turtle_error(p, "expected '$c', found '$(peekc(p))'")
    advance!(p)
end

function grab!(p::TurtleParser, r::Regex)
    m = match(r, p.s, p.i)
    m === nothing && return nothing
    p.i += ncodeunits(m.match)
    return m
end

fresh_bnode!(p::TurtleParser) = IRI("_:b" * string(p.bnodes += 1))

resolve(p::TurtleParser, ref) =
    occursin(r"^[A-Za-z][A-Za-z0-9+.\-]*:", ref) ? String(ref) :
    startswith(ref, '#') ? replace(p.base, r"#.*$" => "") * ref : p.base * ref

function parse_iriref!(p::TurtleParser)
    j = findnext('>', p.s, p.i)
    j === nothing && turtle_error(p, "unterminated IRI")
    ref = unescape_literal(p.s[nextind(p.s, p.i):prevind(p.s, j)])
    p.i = nextind(p.s, j)
    return resolve(p, ref)
end

const PNAME = r"\G([A-Za-z\x80-\xFF][\w\-.\x80-\xFF]*)?:((?:[\w\-.:%\x80-\xFF]|\\[_~.\-!$&'()*+,;=/?#@%])*)"

function parse_pname!(p::TurtleParser)
    start = p.i
    m = grab!(p, PNAME)
    m === nothing && turtle_error(p, "expected IRI or prefixed name")
    prefix = something(m[1], "")
    local_ = m[2]
    while endswith(local_, '.')              # a trailing '.' ends the statement
        local_ = chop(local_); p.i -= 1
    end
    endswith(prefix, '.') && turtle_error(p, "invalid prefix '$prefix'")
    haskey(p.prefixes, prefix) || (p.i = start; turtle_error(p, "undefined prefix '$prefix:'"))
    return p.prefixes[prefix] * replace(local_, r"\\(.)" => s"\1")
end

function parse_iri!(p::TurtleParser)
    skipws!(p)
    return IRI(peekc(p) == '<' ? parse_iriref!(p) : parse_pname!(p))
end

function parse_string!(p::TurtleParser)
    q = peekc(p)
    long = lookahead(p, string(q)^3)
    advance!(p, long ? 3 : 1)
    start = p.i
    while true
        eof(p) && turtle_error(p, "unterminated string")
        c = peekc(p)
        if c == '\\'
            advance!(p, 2)
        elseif c == q && (!long || lookahead(p, string(q)^3))
            body = p.s[start:prevind(p.s, p.i)]
            advance!(p, long ? 3 : 1)
            return unescape_literal(body)
        else
            !long && c == '\n' && turtle_error(p, "newline in short string")
            advance!(p)
        end
    end
end

function parse_literal!(p::TurtleParser)
    value = parse_string!(p)
    if peekc(p) == '@'
        m = grab!(p, r"\G@([A-Za-z]+(?:-[A-Za-z0-9]+)*)")
        return Literal(value, nothing, m[1])
    elseif lookahead(p, "^^")
        advance!(p, 2)
        return Literal(value, parse_iri!(p).value, nothing)
    end
    return Literal(value)
end

function parse_blank!(p::TurtleParser)
    if lookahead(p, "_:")
        m = grab!(p, r"\G_:([\w\-.]*[\w\-])")
        m === nothing && turtle_error(p, "invalid blank node label")
        return IRI("_:" * m[1])
    end
    advance!(p)                                  # '['
    node = fresh_bnode!(p)
    skipws!(p)
    peekc(p) == ']' || parse_predicate_objects!(p, node)
    expect!(p, ']')
    return node
end

function parse_collection!(p::TurtleParser)
    advance!(p)                                  # '('
    items = Union{IRI,Literal}[]
    while (skipws!(p); peekc(p) != ')')
        eof(p) && turtle_error(p, "unterminated collection")
        push!(items, parse_object!(p))
    end
    advance!(p)
    isempty(items) && return IRI(RDF * "nil")
    head = fresh_bnode!(p)
    node = head
    for (k, item) in enumerate(items)
        push!(p.triples, Triple(node, IRI(RDF * "first"), item))
        next = k == length(items) ? IRI(RDF * "nil") : fresh_bnode!(p)
        push!(p.triples, Triple(node, IRI(RDF * "rest"), next))
        node = next
    end
    return head
end

const NAME_END = r"\G(?![\w\-:])"

function parse_object!(p::TurtleParser)
    skipws!(p)
    c = peekc(p)
    c == '"' || c == '\'' ? parse_literal!(p) :
    c == '[' || lookahead(p, "_:") ? parse_blank!(p) :
    c == '(' ? parse_collection!(p) :
    c == '<' ? parse_iri!(p) :
    (m = grab!(p, r"\G[+-]?(?:\d+\.\d*[eE][+-]?\d+|\.\d+[eE][+-]?\d+|\d+[eE][+-]?\d+)")) !== nothing ?
        Literal(m.match, XSD_DOUBLE, nothing) :
    (m = grab!(p, r"\G[+-]?\d*\.\d+")) !== nothing ? Literal(m.match, XSD_DECIMAL, nothing) :
    (m = grab!(p, r"\G[+-]?\d+")) !== nothing ? Literal(m.match, XSD_INTEGER, nothing) :
    (m = grab!(p, r"\G(?:true|false)(?![\w\-:])")) !== nothing ? Literal(m.match, XSD_BOOLEAN, nothing) :
    parse_iri!(p)
end

function parse_verb!(p::TurtleParser)
    skipws!(p)
    if peekc(p) == 'a' && match(r"\Ga(?![\w\-:.])|\Ga(?=\.\s)", p.s, p.i) !== nothing
        advance!(p)
        return IRI(RDF_TYPE)
    end
    return parse_iri!(p)
end

function parse_predicate_objects!(p::TurtleParser, subject::IRI)
    while true
        verb = parse_verb!(p)
        while true
            push!(p.triples, Triple(subject, verb, parse_object!(p)))
            skipws!(p)
            peekc(p) == ',' || break
            advance!(p)
        end
        # one or more ';', optionally followed by another predicate
        peekc(p) == ';' || return
        while (skipws!(p); peekc(p) == ';')
            advance!(p)
        end
        c = peekc(p)
        (c == '.' || c == ']' || c == '\0') && return
    end
end

function parse_directive!(p::TurtleParser)
    sparql = !lookahead(p, "@")
    m = grab!(p, r"\G@?(prefix|base)\b"i)
    if lowercase(m[1]) == "prefix"
        skipws!(p)
        pm = grab!(p, r"\G([A-Za-z][\w\-.]*)?:")
        pm === nothing && turtle_error(p, "expected prefix name")
        skipws!(p)
        peekc(p) == '<' || turtle_error(p, "expected IRI after prefix")
        p.prefixes[something(pm[1], "")] = parse_iriref!(p)
    else
        skipws!(p)
        p.base = parse_iriref!(p)
    end
    sparql || expect!(p, '.')
end

"""
    parse_turtle(io_or_path; base = "") -> (triples, prefixes)

Parse a Turtle document into triples. Also returns the declared prefix map.
"""
function parse_turtle(io::IO; base::AbstractString = "")
    p = TurtleParser(read(io, String), 1, String(base), Dict{String,String}(), 0, Triple[])
    while (skipws!(p); !eof(p))
        if match(r"\G(@prefix|@base)\b|\G(?i:prefix|base)\s"i, p.s, p.i) !== nothing
            parse_directive!(p)
            continue
        end
        c = peekc(p)
        subject = c == '[' || lookahead(p, "_:") ? parse_blank!(p) :
                  c == '(' ? parse_collection!(p) : parse_iri!(p)
        skipws!(p)
        c == '[' && peekc(p) == '.' || parse_predicate_objects!(p, subject)
        expect!(p, '.')
    end
    return p.triples, p.prefixes
end
parse_turtle(path::AbstractString; kw...) = open(io -> parse_turtle(io; kw...), path)

# Use a document's prefixes when compacting ids: its `bnr:` becomes the graph base.
function graph_options(pfx::AbstractDict; base = nothing, prefixes = Dict{String,String}())
    standard = Set(["bn", "rdf", "rdfs", "xsd", "owl", ""])
    b = base === nothing ? get(pfx, "bnr", BNR) : base
    user = Dict{String,String}(k => v for (k, v) in pfx if !(k in standard) && k != "bnr")
    return (base = b, prefixes = merge(user, Dict{String,String}(prefixes)))
end

"""    read_turtle(io_or_path; base = nothing, prefixes = Dict()) -> KnowledgeGraph"""
function read_turtle(src; base = nothing, prefixes = Dict{String,String}())
    ts, pfx = parse_turtle(src)
    return from_triples(ts; graph_options(pfx; base, prefixes)...)
end
