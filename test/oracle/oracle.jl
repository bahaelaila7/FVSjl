# =============================================================================
# oracle.jl — the dual-oracle test harness
#
# A modernization rewrite is only safe if every step can be checked against a
# known-good reference. We have two:
#
#   Oracle A — the faithful port at /workspace/FVSjulia. Already matches Fortran
#              to ulp, always present, no rebuild needed. PRIMARY reference.
#   Oracle B — the live Fortran FVSsn. Rebuildable at /tmp/FVSsn_new. Used for
#              final confirmation (see tests/test_fortran_diff.sh in FVSjulia).
#
# This module runs a `.key` through Oracle A (and, once the engine exists, through
# FVSjl) and provides tolerant diff helpers. FVS writes outputs (.sum/.out/.db)
# NEXT TO the keyword file, so each run happens in an isolated temp dir.
# =============================================================================
module Oracle

using Test

const FVSJULIA_DIR = "/workspace/FVSjulia"                       # Oracle A package
const FVSSN_TESTS  = "/workspace/ForestVegetationSimulator/tests/FVSsn"

"Locate a julia binary (mirrors FVSjulia/tests/find_julia.sh)."
function julia_bin()
    f = joinpath(FVSJULIA_DIR, "tests", "find_julia.sh")
    if isfile(f)
        try
            return strip(read(`bash $f`, String))
        catch
        end
    end
    return Base.julia_cmd()[1]
end

"""
    run_in_tempdir(keypath; runner) -> tempdir

Copy `keypath` (+ a sibling `.tre` if present) into a fresh temp dir, invoke
`runner(tmpdir, keyfile)`, and return the temp dir holding the outputs. Caller
is responsible for cleanup.
"""
function run_in_tempdir(keypath::AbstractString, runner)
    tmp = mktempdir()
    key = joinpath(tmp, basename(keypath))
    cp(keypath, key; force=true)
    tre = first(splitext(keypath)) * ".tre"
    isfile(tre) && cp(tre, joinpath(tmp, basename(tre)); force=true)
    runner(tmp, key)
    return tmp
end

"Run Oracle A (FVSjulia) on `keypath`; returns its temp output dir."
function run_oracle_a(keypath::AbstractString)
    jl = julia_bin()
    run_in_tempdir(keypath) do tmp, key
        code = "using FVSjulia; FVSjulia.main([\"--keywordfile=$key\"])"
        run(pipeline(`$jl --project=$FVSJULIA_DIR -e $code`; stdout=devnull, stderr=devnull))
    end
end

"Run FVSjl on `keypath`; returns its temp output dir. (Wired up from C2 onward.)"
function run_fvsjl(keypath::AbstractString)
    jl = julia_bin()
    proj = normpath(joinpath(@__DIR__, "..", ".."))
    run_in_tempdir(keypath) do tmp, key
        code = "using FVSjl; FVSjl.main([\"--keywordfile=$key\"])"
        run(pipeline(`$jl --project=$proj -e $code`; stdout=devnull, stderr=devnull))
    end
end

# ---------------------------------------------------------------------------
# Diff helpers
# ---------------------------------------------------------------------------

"Find a single output file with the given extension inside `dir`."
function find_output(dir::AbstractString, ext::AbstractString)
    hits = filter(f -> endswith(lowercase(f), lowercase(ext)), readdir(dir))
    isempty(hits) && return nothing
    return joinpath(dir, first(hits))
end

const _NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

"""
    diff_text_numeric(path_a, path_b; atol, rtol) -> Vector{String}

Compare two text files line by line. Numbers are compared with tolerance; all
non-numeric text must match exactly. Returns a list of human-readable diffs
(empty == identical within tolerance).
"""
function diff_text_numeric(a::AbstractString, b::AbstractString; atol=1e-3, rtol=1e-4)
    la = readlines(a); lb = readlines(b)
    diffs = String[]
    n = max(length(la), length(lb))
    for i in 1:n
        ra = i <= length(la) ? la[i] : ""
        rb = i <= length(lb) ? lb[i] : ""
        ra == rb && continue
        # tokenize on numbers; compare structure + numeric tolerance
        na = collect(m.match for m in eachmatch(_NUM, ra))
        nb = collect(m.match for m in eachmatch(_NUM, rb))
        skel_a = replace(ra, _NUM => "#")
        skel_b = replace(rb, _NUM => "#")
        if skel_a != skel_b || length(na) != length(nb)
            push!(diffs, "L$i text: |$ra| vs |$rb|")
            continue
        end
        for (x, y) in zip(na, nb)
            fx = tryparse(Float64, x); fy = tryparse(Float64, y)
            (fx === nothing || fy === nothing) && continue
            if !isapprox(fx, fy; atol=atol, rtol=rtol)
                push!(diffs, "L$i num: $fx vs $fy")
            end
        end
    end
    return diffs
end

"Assert FVSjl's `.sum` matches Oracle A's for `keyname` (a basename in FVSsn tests)."
function assert_sum_matches(keyname::AbstractString; kwargs...)
    keypath = joinpath(FVSSN_TESTS, keyname)
    da = run_oracle_a(keypath)
    db = run_fvsjl(keypath)
    try
        sa = find_output(da, ".sum"); sb = find_output(db, ".sum")
        @test sa !== nothing
        @test sb !== nothing
        d = diff_text_numeric(sa, sb; kwargs...)
        isempty(d) || (@info "SUM diffs ($keyname)" first(d, 20))
        @test isempty(d)
    finally
        rm(da; recursive=true, force=true)
        rm(db; recursive=true, force=true)
    end
end

end # module Oracle
