# seed_test.jl — the SEED-PERTURBATION adjudicator: does a per-stand jl−oracle divergence FLIP with the RNG
# seed (⇒ a genuine RANSED realization STRADDLE, cornerable) or is it SEED-INVARIANT and one-directional
# (⇒ a REAL deterministic count/height bug, must be FIXED not cornered)?
#
# WHY THIS EXISTS. The ledger's `count_divergence_UNVERIFIED` signature (structure moves only in TPA/QMD while
# BA/SDI/CCF/TopHt are preserved) is AMBIGUOUS: BA converges BY CONSTRUCTION both for a benign self-thin
# realization straddle AND for a real bug that merely redistributes TPA across size-classes. The ONLY way to
# tell them apart is to reseed the main RNG identically in jl AND the oracle across several seeds and watch the
# SIGNED residual: a straddle lands on either side of the oracle depending on the seed; a deterministic bug does
# not move. (2026-09-04 IE `none`-regime re-diagnosis: this test un-cornered a batch of seed-invariant −574…±100
# TPA deterministic bugs that had been falsely tagged "#206 OLDRN/ZRAND straddle".)
#
# ⚠ RANNSEED CAVEATS (READ BEFORE TRUSTING A VERDICT):
#   1. RANNSEED (FVS option 61) reseeds ONLY the MAIN stream (RANSED). The ESTABLISHMENT stream (ESRANN) is
#      FIXED at 55329 regardless. So the ESTABLISHMENT-cohort realization does NOT vary with RANNSEED — an
#      establishment-driven divergence is SEED-INVARIANT BY CONSTRUCTION and this test will (correctly) NOT see
#      it flip. For a bare-establishment / regen-cohort stand, "seed-invariant" therefore does NOT by itself
#      prove a deterministic bug; the verdict there must instead rest on the ONE-DIRECTIONAL population
#      signature (a systematic jl−oracle sign across a POPULATION sample of establishment stands — see the
#      adjudication doc's stratified sign-tally), NOT on this single-stand reseed alone.
#   2. The seed value MUST sit in the keyword RECORD (cols 11-20) — keytext()/run_live() do this via
#      kwrec("RANNSEED", seed). A value on a blank supplemental line is a SILENT FVS no-op (reuses seed 55329).
#      Identical keytext feeds jl and the oracle so both genuinely reseed to the same value.
#
# Usage:
#   JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project=. test/harness/fia/seed_test.jl <VARIANT> <CN> [CN...]
#     [--seeds=55329,1,777,424242,987654] [--bin=/path/FVSxx] [--regime=none]
#   Stands are read from the ledger's MASTER FIA DB (/workspace/SQLite_FIADB_ENTIRE.db, via build_subdb).
#   Oracle binary: --bin, else the ledger BIN dict, else /workspace/.<vl>work/FVS<vl>_g16 (western g16 oracles).
#
# Output: per stand, a per-seed signed-residual (jl − oracle, last common cycle) table for TPA/BA/TopHt/QMD, and
# a PER-METRIC verdict: SEED-VARIANT(straddle) | SEED-INVARIANT(deterministic) | ~0(match). A stand with ANY
# metric SEED-INVARIANT & material is a real-bug candidate (do NOT corner); one whose material metrics all FLIP
# is a genuine straddle (cornerable).

include(joinpath(@__DIR__, "ledger_fia.jl"))   # keytext / run_live / parse_sum10 / build_subdb / kwrec / COLS / VAR / BIN

const DEFAULT_SEEDS = [55329, 1, 777, 424242, 987654]
# metrics adjudicated (name => COLS index) + a per-metric NOISE FLOOR ε below which a residual is "match":
#   a signed residual only counts as a signal beyond this absolute floor (separates a real move from a ±1-2
#   print/AVHT40 jitter). A metric SIGN-FLIPS iff it has one residual < -ε AND one > +ε across the seed sweep.
const METRICS = [("TPA",1,4.0), ("BA",2,3.0), ("TopHt",5,3.0), ("QMD",6,1.0)]
# The STAND verdict keys on the DECISIVE metrics, not "any metric": the classifier flagged a COUNT divergence,
# so TPA is the count metric that adjudicates straddle-vs-bug, with BA as density corroboration. TopHt only
# signals a real HEIGHT bug (the D2/Mode-A small-tree height deficit) when it is LARGE and seed-invariant — a
# ±5-8 ft one-sided TopHt residual is AVHT40 mortality-realization jitter, NOT a height-model signal (per the
# 2026-09-04 adjudication), so a height residual must clear this floor before it flips the stand to a bug.
# QMD is arithmetically derived from TPA/BA (not independent) ⇒ display-only, never drives the stand verdict.
const HEIGHT_DETERM_FLOOR = 10.0

# Oracle binary for ANY variant (eastern dict wins; western → the freshly-relinked g16 oracle).
seed_oracle_bin(v, override) =
    override != "" ? override :
    get(BIN, v) do; "/workspace/.$(lowercase(v))work/FVS$(lowercase(v))_g16"; end
seed_variant(v) = get(VAR, v) do; FVSjl.variant_from_code(v); end

# last common .sum cycle's signed residual (jl − oracle) for every metric col; nothing if no comparable cycle.
function last_cycle_residuals(livetext, jltext)
    L = parse_sum10(livetext); J = parse_sum10(jltext)
    (isempty(L) || isempty(J)) && return nothing
    Jd = Dict(y=>vv for (y,vv) in J)
    ys = sort([y for (y,_) in L if haskey(Jd, y)])
    isempty(ys) && return nothing
    y = last(ys); lv = Dict(L)[y]; jv = Jd[y]
    (year=y, res=Dict(name => jv[k]-lv[k] for (name,k,_) in METRICS))
end

# verdict over the signed residuals of ONE metric across all seeds.
function metric_verdict(vals, eps)
    amax = maximum(abs, vals)
    amax <= eps && return ("~0(match)", amax)
    (minimum(vals) < -eps && maximum(vals) > eps) && return ("SEED-VARIANT(straddle)", amax)
    ("SEED-INVARIANT(determ)", amax)                     # all same sign beyond noise ⇒ one-directional
end

function seedtest_stand(v, cn, seeds, bin, regime)
    dir = mktempdir(); sub = joinpath(dir, "sub.db")
    build_subdb([cn], sub)
    var = seed_variant(v)
    # residual[metric] = Vector aligned with `seeds` (missing seeds recorded as NaN, shown as "·")
    resid = Dict(name => Float64[] for (name,_,_) in METRICS)
    ok_seeds = Int[]; last_year = 0
    for s in seeds
        live, _ = run_live(bin, cn, sub, regime, dir, 0, s)
        keyf = joinpath(dir, "jl.key"); write(keyf, keytext(cn, sub, regime, 0, s))
        jlout = try FVSjl.run_keyfile(keyf; variant=var) catch e; @warn "jl run failed" cn s e; ""; end
        r = last_cycle_residuals(live, jlout)
        if r === nothing
            for (name,_,_) in METRICS; push!(resid[name], NaN); end
        else
            last_year = r.year; push!(ok_seeds, s)
            for (name,_,_) in METRICS; push!(resid[name], r.res[name]); end
        end
    end
    (resid=resid, seeds=seeds, ok=ok_seeds, year=last_year)
end

fmt(x) = isnan(x) ? "    ·" : (x == round(x) ? lpad(string(Int(round(x))),6) : lpad(string(round(x,digits=1)),6))

function main_seedtest(v, cns, seeds, binoverride, regime)
    bin = seed_oracle_bin(v, binoverride)
    isfile(bin) || error("oracle binary not found: $bin  (pass --bin=... )")
    println("SEED-PERTURBATION TEST  variant=$v  oracle=$bin  regime=$regime")
    println("seeds (RANNSEED, cols 11-20 of the record): $(join(seeds, ", "))")
    println("residual = jl − oracle at the last common .sum cycle. ε(TPA/BA/TopHt/QMD)=4/3/3/1.")
    println("⚠ establishment stream (ESRANN=55329) does NOT reseed — an estab-cohort divergence is seed-invariant")
    println("  BY CONSTRUCTION; for bare-estab stands use the population one-directional signature, not this alone.\n")
    ndeterm = 0; nstraddle = 0
    for cn in cns
        t = seedtest_stand(v, cn, seeds, bin, regime)
        println("── CN $cn   (last common cycle year $(t.year); ran $(length(t.ok))/$(length(seeds)) seeds)")
        println("  metric | " * join([lpad(string(s),6) for s in seeds], " ") * "  |  verdict (|max| residual)")
        vd = Dict{String,Tuple{String,Float64}}()
        for (name,_,eps) in METRICS
            vals = t.resid[name]
            good = filter(!isnan, vals)
            verdict, amax = isempty(good) ? ("no-data", 0.0) : metric_verdict(good, eps)
            vd[name] = (verdict, amax)
            println("  " * rpad(name,6) * " | " * join([fmt(x) for x in vals], " ") *
                    "  |  $verdict  (|max|=$(round(amax,digits=1)))")
        end
        isd(m) = occursin("determ", vd[m][1]); iss(m) = occursin("straddle", vd[m][1])
        # decisive-count bug: TPA (or density BA) seed-invariant & material. height bug: TopHt seed-invariant AND
        # ≥ HEIGHT_DETERM_FLOOR (excludes AVHT40 jitter). count straddle: the flagged count metric TPA flips.
        determ_count  = isd("TPA") || isd("BA")
        determ_height = isd("TopHt") && vd["TopHt"][2] >= HEIGHT_DETERM_FLOOR
        stand_determ  = determ_count || determ_height
        stand_straddle = iss("TPA")
        cls = stand_determ ?
                "REAL DETERMINISTIC BUG — do NOT corner; FIX" *
                    (determ_height && !determ_count ? " (height/D2)" : " (count/D1-modeB)") *
                    ".  ⚠ if this is a bare-ESTABLISHMENT stand, seed-invariance is BY CONSTRUCTION (ESRANN fixed)" *
                    " — confirm via the population one-directional signature, not this reseed alone." :
              stand_straddle ?
                "GENUINE STRADDLE — cornerable; the count residual FLIPS sign with the seed." :
                "MATCH / negligible."
        println("  ⇒ $cls\n")
        stand_determ && (ndeterm += 1); (stand_straddle && !stand_determ) && (nstraddle += 1)
    end
    println("SUMMARY: $(length(cns)) stands  |  real-deterministic (any metric seed-invariant & material): $ndeterm" *
            "  |  pure-straddle: $nstraddle")
end

function _cli(args)
    seeds = DEFAULT_SEEDS; binoverride = ""; regime = "none"
    pos = String[]
    for a in args
        if startswith(a, "--seeds=");      seeds = parse.(Int, split(a[9:end], ","))
        elseif startswith(a, "--bin=");    binoverride = a[7:end]
        elseif startswith(a, "--regime="); regime = a[10:end]
        else push!(pos, a); end
    end
    length(pos) >= 2 || error("usage: seed_test.jl <VARIANT> <CN> [CN...] [--seeds=a,b,c] [--bin=...] [--regime=none]")
    main_seedtest(pos[1], pos[2:end], seeds, binoverride, regime)
end

if abspath(PROGRAM_FILE) == @__FILE__
    _cli(ARGS)
end
