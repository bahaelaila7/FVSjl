# mutation.jl — SENSITIVITY / MUTATION CHECK of the tiered suite.
#
# Applies small, realistic mutations one at a time — each in a fresh temp git worktree of HEAD (so the working tree
# and the shared precompile cache of this checkout are never touched) — and runs the QUICK snapshot tier
# (snapshot.jl check, TIERED=quick) there. A mutation is KILLED when the check reports a difference (exit 1). The kill
# rate measures whether the bit-identity net would catch the kind of subtle change a refactor can introduce.
# Requires committed SNAPSHOT.tsv manifests (snapshot.jl bless) at HEAD.
# Usage: julia --project=. test/harness/tiered/mutation.jl [id ...]     (default: all mutations)
# Env: TIERED_VARIANTS (default: the CORE group — tiered_runner.jl TIERED_GROUPS)
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# (id, kind, file, old, new) — `old` must occur exactly once in the file at HEAD.
const MUTATIONS = [
    ("assoc-ba", "swap a*b*c association", "src/engine/standstats.jl",
     "ba += 0.005454154f0 * (d * (d * t.tpa[i]))", "ba += 0.005454154f0 * ((d * d) * t.tpa[i])"),
    ("tie-rdpsrt", "flip a sort tie-break (RDPSRT >= → >)", "src/engine/cuts.jl",
     "if key[indil] >= key[indiu]; @goto l40; end", "if key[indil] > key[indiu]; @goto l40; end"),
    ("ulp-baptree", "perturb a coefficient by 1 ULP (cuts BA factor)", "src/engine/cuts.jl",
     "const _BA_PER_TREE = 0.005454154f0", "const _BA_PER_TREE = 0.0054541547f0"),
    ("f64-sdimax", "Float32 → Float64 intermediate (stand SDImax weights)", "src/engine/standstats.jl",
     "tb = 0.0054542f0 * t.dbh[i]^2 * t.tpa[i]", "tb = Float32(0.0054542 * Float64(t.dbh[i])^2 * Float64(t.tpa[i]))"),
    ("libm-fpowi", "replace gfortran __powisf2 emulation by generic Float32^Int", "src/core/fmath.jl",
     "    n = unsigned(abs(m)); y = isodd(n) ? x : 1f0", "    return x^m; n = unsigned(abs(m)); y = isodd(n) ? x : 1f0"),
    ("coef-bmdf-bark", "perturb a data coefficient (BM DF bark1 0.903563 → 0.9035631)", "data/bluemountains/species_coefficients.csv",
     "3,DF,202,PSME,50.0,110.0,50.0,0.903563,", "3,DF,202,PSME,50.0,110.0,50.0,0.9035631,"),
]

function run_mutation(m)
    id, kind, file, old, new = m
    wt = mktempdir(); rm(wt)
    run(pipeline(`git -C $ROOT worktree add --detach $wt HEAD`; stdout = devnull, stderr = devnull))
    try
        for f in ("deps/libfvsmath.so",)
            isfile(joinpath(ROOT, f)) && cp(joinpath(ROOT, f), joinpath(wt, f); force = true)
        end
        p = joinpath(wt, file); s = read(p, String)
        n = count(old, s); n == 1 || return (id, kind, "INVALID (target occurs $n×)")
        write(p, replace(s, old => new))
        env = copy(ENV); env["TIERED"] = "quick"
        cmd = `$(Base.julia_cmd()) --project=$wt $(joinpath(wt, "test/harness/tiered/snapshot.jl")) check`
        log = joinpath(mktempdir(), "check.log")
        pr = run(pipeline(ignorestatus(setenv(cmd, env)); stdout = log, stderr = log))
        out = read(log, String)
        verdict = pr.exitcode == 1 && occursin("differing case", out) ? "KILLED" :
                  pr.exitcode == 0 ? "SURVIVED" : "ERROR (exit $(pr.exitcode)): " * last(out, 400)
        firstdiff = something(match(r"DIFF .*", out), (match = "",)).match
        return (id, kind, verdict * (isempty(firstdiff) ? "" : "  — " * first(firstdiff, 160)))
    finally
        run(pipeline(`git -C $ROOT worktree remove --force $wt`; stdout = devnull, stderr = devnull))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    sel = isempty(ARGS) ? MUTATIONS : [m for m in MUTATIONS if m[1] in ARGS]
    res = [run_mutation(m) for m in sel]
    killed = count(r -> startswith(r[3], "KILLED"), res)
    println("\n== mutation check (quick snapshot tier) ==")
    for r in res; println(rpad(r[1], 16), rpad(r[2], 62), r[3]); end
    println("kill rate: $killed / $(length(res))")
    exit(killed == length(res) ? 0 : 1)
end
