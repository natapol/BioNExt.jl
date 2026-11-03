# Run every ```julia block of the hand-written manual pages, in order, one fresh module
# per page, inside an empty temporary directory. Expected results are checked too: a
# single-line expression ending in `# -> value`, or an expression followed by comment
# lines `# -> value` (continued on further `#` lines). If `value` parses and evaluates as
# Julia, it must equal the result; descriptive text is ignored. A block preceded by the
# line `<!-- norun -->` is skipped.

const MANUAL = joinpath(@__DIR__, "..", "docs", "src")

julia_blocks(md) = [m[2] for m in eachmatch(r"(<!-- norun -->\s*\n\s*)?^```julia[ \t]*\n(.*?)^```"ms, md)
                     if m[1] === nothing]

function expected_value(mod, line)
    m = match(r"#\s*->\s*(.+?)\s*$", line)
    m === nothing && return nothing
    try
        return Some(Core.eval(mod, Meta.parse(m[1])))
    catch
        return nothing          # descriptive text, not a value
    end
end

function check!(mod, path, source, result, expected_line)
    exp = expected_value(mod, expected_line)
    exp === nothing && return 0
    @test isequal(result, something(exp)) ||
          (println(stderr, "Mismatch in $(basename(path)): $(strip(source))\n  got: ", repr(result)); false)
    return 1
end

# Expected value for a chunk (one parsed expression plus the comments that follow it):
# an inline `# -> value` on its only code line, or `# -> value` comment lines after its
# last code line, continued on further `#` lines.
function expectation(chunk)
    lines = split(chunk, '\n')
    code = findall(l -> !isempty(strip(l)) && !startswith(strip(l), "#"), lines)
    isempty(code) && return nothing
    length(code) == 1 && occursin(r"#\s*->", lines[code[1]]) && return lines[code[1]]
    rest = lines[last(code)+1:end]
    i = findfirst(l -> occursin(r"^\s*#\s*->", l), rest)
    i === nothing && return nothing
    all(l -> startswith(strip(l), "#"), rest[1:i]) || return nothing
    text = [strip(rest[i])]
    for l in rest[i+1:end]
        m = match(r"^\s*#(?!\s*->)\s*(.*)$", l)
        m === nothing && break
        push!(text, m[1])
    end
    return join(text, " ")
end

function run_page(path)
    mod = Module(Symbol(splitext(basename(path))[1]))
    Core.eval(mod, :(using Test))
    checked = 0
    for (b, code) in enumerate(julia_blocks(read(path, String)))
        pos = 1
        while pos <= ncodeunits(code)
            ex, next = Meta.parse(code, pos; raise = true)
            chunk = code[pos:prevind(code, next)]
            pos = next
            ex === nothing && continue
            result = try
                Core.eval(mod, ex)
            catch err
                error("$(basename(path)), block $b failed at:\n$(strip(chunk))\n\n$(sprint(showerror, err))")
            end
            e = expectation(chunk)
            e === nothing || (checked += check!(mod, path, chunk, result, e))
        end
    end
    return checked
end

@testset "manual examples" begin
    pages = sort!(vcat([joinpath(MANUAL, "index.md")],
                       [joinpath(MANUAL, "manual", f) for f in readdir(joinpath(MANUAL, "manual")) if endswith(f, ".md")]))
    for page in pages
        @testset "$(basename(page))" begin
            cd(mktempdir()) do
                redirect_stdout(() -> run_page(page), devnull)
            end
        end
    end
end
