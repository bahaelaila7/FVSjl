# ledger_fia.jl — REPRODUCIBLE per-stand divergence ledger for the FIA/FVS differential.
#
# Purpose: a durable, committed record of EVERY sampled stand's status (bit-exact or diverging, by how much,
# with a deterministic signature) so that after a later bug fix you re-run the SAME stands and DIFF the CSV to
# see status flips. No hand-wavey classification — the primary columns are MEASURED FACTS; the `signature`
# column is a deterministic rule over those facts, with UNCLASSIFIED for anything that does not cleanly match.
#
# Self-contained: given a committed stand-CN list, it builds a TEMP indexed sub-DB from the read-only master
# (ATTACH + CREATE TABLE ... WHERE STAND_CN IN (...)), runs live FVS + FVSjl per stand over the full default
# horizon, compares all 10 .sum columns every cycle, and APPENDS one CSV row per stand.
#
# Usage: LEDGER=<out.csv> julia --project=. test/harness/fia/ledger_fia.jl <standlist> <SN|NE|CS|LS> [regime]
#   regime ∈ {none(default), simfire, thinbba, salvage, plant}. Re-run for each (variant,regime) to fill the ledger.
#
# CSV columns:
#   variant, regime, cn, n_cycles, bit_exact, div_cols, worst_col, worst_cycle, max_rel_pct, max_abs_diff,
#   struct_max_rel_pct, vol_max_rel_pct, density_bitexact, converges, signature
#   - div_cols     : '|'-joined names of the columns that EVER diverge (empty if bit-exact)
#   - worst_col/cycle/max_rel_pct/max_abs_diff : the single largest-RELATIVE-divergence cell (max_abs = largest
#                    ABSOLUTE diff across ANY col — may be a different col than worst_col)
#   - struct_max_rel_pct / vol_max_rel_pct : largest relative divergence WITHIN cols 1-6 / cols 7-10 separately
#   - density_bitexact : BA+SDI+CCF+TopHt never diverge by more than 1 unit / 1% (the count-straddle signature)
#   - converges    : last-cycle max-rel < 0.5×peak max-rel, or last cycle bit-exact
#   - signature    : deterministic bucket (see classify()); UNCLASSIFIED ⇒ needs a manual both-sides trace.
# `MATERIAL` threshold below: a divergence < MATERIAL rel AND ≤1 abs unit is a print/ULP straddle, not a
# structural signal — so a ±1 SDI straddle does not masquerade as a dense-phase divergence.
const MATERIAL = 0.01f0

import SQLite, DBInterface
using FVSjl
include(joinpath(@__DIR__, "sweep_db.jl"))   # open_sweepdb / upsert! — durable per-stand coverage record
const MASTER = "/workspace/SQLite_FIADB_ENTIRE.db"
# Oracle binary per variant. Western oracles live in /workspace/.<v>work/ (persistent); a fresh relink lands as
# FVS<v>_g16 (preferred), else the older FVS<v>_clean. ⚠ RESTART⇒RELINK ALL before a sweep (binaries can go stale);
# _oracle_bin picks the g16 if present. Eastern keep their tmp/oracles/*_new paths.
_oracle_bin(v) = let w = "/workspace/." * lowercase(v) * "work"
    isfile("$w/FVS$(lowercase(v))_g16") ? "$w/FVS$(lowercase(v))_g16" : "$w/FVS$(lowercase(v))_clean"
end
const BIN = Dict("SN"=>"/workspace/FVSjl/tmp/oracles/FVSsn_new","NE"=>"/workspace/FVSjl/tmp/oracles/FVSne_new","CS"=>"/workspace/FVSjl/tmp/oracles/FVScs_new","LS"=>"/workspace/FVSjl/tmp/oracles/FVSls_new",
    ("$c"=>_oracle_bin(c) for c in ("CR","KT","IE","EM","TT","UT","BM","CI","BC","NC","OC","AK","WC","PN","EC","CA","SO","WS"))...)
const VAR = Dict("SN"=>FVSjl.Southern(),"NE"=>FVSjl.Northeast(),"CS"=>FVSjl.CentralStates(),"LS"=>FVSjl.LakeStates(),
    ("$c"=>FVSjl.variant_from_code(c) for c in ("CR","KT","IE","EM","TT","UT","BM","CI","BC","NC","OC","AK","WC","PN","EC","CA","SO","WS"))...)
const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD","TCuFt","MCuFt","SCuFt","BdFt"]   # .sum fields 3..12
const DENSITY_COLS = (2,3,4,5)   # BA,SDI,CCF,TopHt — preserved under the self-thinning count-straddle

kwrec(kw, f...) = rpad(kw,10) * join(lpad(string(x),10) for x in f)
# PLANT uses a CALENDAR-year date (the standard form every validated test + real scenario uses), NOT a
# cycle-number date: `PLANT 2.0` hits FVS's cycle-number→age scheduler path where jl has a small residual
# (audit slice 42d-42g), which tanked the PLANT bit-exact rate to ~0% as a HARNESS ARTIFACT. `plantyr` is the
# stand's INV_YEAR+period (cycle 1) passed from main(). The other activities keep cycle-number "2.0" (their
# scheduling is age-independent and already bit-exact at normal rates). plantyr=0 ⇒ fall back to cycle "2.0".
# --- Climate-FVS regime scaffolding: Climate-FVS's growth/mort multipliers fire ONLY when CLIMDATA populates
#     s.climate (LCLIMATE = NATTRS>0 && NYEARS>0, clin.f); a bare GROWMULT card is caught by the generic
#     multiplier handler and DISCARDED without CLIMDATA. FIA stands carry no climate scenario, so — exactly like
#     WSBWE's synthetic synwx.txt (fed byte-identically to jl AND the oracle) — we template an INLINE CLIMDATA
#     block per stand: the real CGCM3_A2 attribute row from the validated fixture (clim_iet.key), re-keyed to
#     THIS stand's CN, at years 1990 & 2100 (bracketing every FIA inventory; algslp interpolates flat). Both
#     sides read the identical inline data, so the regime is fully oracle-comparable.
# ⚠ REVIEW-GATE FIX (2026-09-06): the climate scenario MUST be TIME-VARYING or Climate-FVS is inert. GROWMULT
# (clgmult.f:225) scales only the DEVIATION of the climate score PS from 1.0: TREEMULT=1+((PS−1)·CLGROWMULT). A
# FLAT scenario (identical rows both years) ⇒ every _NOW−_INVYR delta = 0 ⇒ PS=1 ⇒ TREEMULT=1 regardless of
# GROWMULT ⇒ byte-identical to none (tested a WARM 2100 row → it bites). So we template the fixture's REAL
# CGCM3_A2 1990 (cool) row at year 1990 and its 2090 (warm: mat 3.7→8.6, dd5 1095→2087) row at year 2100 — the
# genuine warming trajectory. Both sides read the identical inline data ⇒ oracle-comparable.
const _CLIM_HEADER, _CLIM_TAIL_1990, _CLIM_TAIL_2090 = let hdr = "", t90 = "", t20 = ""
    for ln in eachline(joinpath(dirname(@__DIR__), "..", "fixtures", "climate", "clim_iet.key"))
        startswith(ln, "Stand_ID,Scenario,Year") && (hdr = ln)
        startswith(ln, "S248112,CGCM3_A2,1990,") && (t90 = ln[length("S248112,CGCM3_A2,1990,")+1:end])
        startswith(ln, "S248112,CGCM3_A2,2090,") && (t20 = ln[length("S248112,CGCM3_A2,2090,")+1:end])
    end
    (hdr, t90, t20)
end
# the ClimData sub-block ONLY (the enclosing CLIMATE…End is supplied by regime_block, once): ClimData +
# scenario name + `*` (inline) + header + two DIFFERING year rows re-keyed to `cn` + the -999 terminator.
climate_block(cn) = string(
    "ClimData\nCGCM3_A2\n*\n", _CLIM_HEADER, "\n",
    cn, ",CGCM3_A2,1990,", _CLIM_TAIL_1990, "\n",
    cn, ",CGCM3_A2,2100,", _CLIM_TAIL_2090, "\n-999")

regime_block(r, plantyr=0, cn="") =
    r == "simfire" ? "FMIn\n" * kwrec("SIMFIRE","2.0","10.00","1","50.0") * "\nEnd" :
    r == "thinbba" ? kwrec("THINBBA","2.0","40.0") :
    # SALVAGE is an FFE keyword (fmin.f:47), NOT a base keyword — a bare SALVAGE draws FVS01 INVALID and is
    # DROPPED by the oracle (inert on green stands, but a false-divergence risk on standing-dead stands). Wrap in
    # FMIn…End (review-gate fix 2026-09-06). Removes standing-dead ≥0" up to 999" at 0.9 fraction, cycle 2.
    r == "salvage" ? "FMIn\n" * kwrec("SALVAGE","2.0","0.0","999.0","0.9") * "\nEnd" :
    r == "plant"   ? "ESTAB\n" * kwrec("PLANT", plantyr > 0 ? string(plantyr) : "2.0", "3","400") * "\nEnd" :
    # --- EXTENSION regimes (USER 2026-09-05: enrich the western matrix). Faithful activation keywords per FVSjl's
    #     own dispatch (keyword_dispatch.jl / root_disease.jl).
    # MISTPINF is a SUB-keyword of the MISTOE section (initre.f opt 95 → misin.f keyword loop), NOT a top-level
    # keyword — a bare MISTPINF draws FVS01 INVALID in the oracle (dropped) while jl over-leniently accepts it,
    # manufacturing a false divergence on every host stand. Wrap in MISTOE…END (review-gate fix 2026-09-06).
    r == "mistletoe" ? "MISTOE\n" * kwrec("MISTPINF","1","0","1.0","3.0") * "\nEND" :     # force initial infection all hosts, DMR 3 (misin.f opt10) — bites where host present
    r == "cover"     ? "COVER\nEnd" :                                                   # COVER report SECTION (cvin.f opt12 → activity 900); MUST close with END or it swallows following keywords (FVS01 INVALID). Report-only ⇒ crash/integration coverage
    r == "econ"      ? kwrec("THINBBA","2.0","40.0") * "\nECON\n" * kwrec("ANNUCST","10.0") * "\n" * kwrec("HRVRVN","0","999","30.0") * "\nEnd" : # ECON valuation needs a harvest to value: thin + $10/ac/yr cost + $30/ccf revenue
    r == "climate"   ? "CLIMATE\n" * kwrec("GROWMULT","1","0","1.1") * "\n" * climate_block(cn) * "\nEnd" : # Climate-FVS: inline per-stand CLIMDATA (bites) + growth-mult weight
    # WRD (RDIN block): RRINIT field1<1 ⇒ RANDOM center placement (field1≥1 = DETERMINISTIC placement reads
    # follow-on coordinate cards and DERAILS the reader). fields 2..6 = 1 center, 10 infected + 50 uninfected
    # TPA/ac, 0.5 root-infection prop, 1.0-acre disease area ⇒ actively seeds infection (bites where host present).
    r == "rootdis"   ? "RDIN\n" * kwrec("RRTYPE","1") * "\n" * kwrec("SAREA","1.0") * "\n" * kwrec("RRINIT","0","1","10.0","50.0","0.5","1.0") * "\nEnd" :
    ""

# `seed` (optional): when set, emit a RANNSEED keyword whose value sits in the keyword RECORD columns 11-20
# (via kwrec — field 1) so FVS actually reseeds. ⚠ A RANNSEED value placed on a BLANK supplemental line is
# SILENTLY IGNORED by FVS (a blank field 1 = "restart from the saved seed" ⇒ a no-op that reuses the default
# seed 55329) — the value MUST be in the record itself. RANNSEED (option 61) reseeds ONLY the main RANSED
# stream; the ESTABLISHMENT stream (ESRANN, fixed 55329) is NOT affected — establishment-cohort divergence is
# therefore seed-invariant BY CONSTRUCTION. Identical keytext feeds jl and the oracle, so both reseed the same.
keytext(cn, db, regime, plantyr=0, seed=nothing) = """
STDIDENT
$cn
DATABASE
DSNin
$db
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         5.0
$(seed === nothing ? "" : kwrec("RANNSEED", seed) * "\n")$(regime_block(regime, plantyr, cn))
ECHOSUM
PROCESS
STOP
"""

# The FVS .sum data row is FIXED-WIDTH (sumout.f FORMAT 9014, identical across the eastern
# variants SN/NE/CS/LS): `2I4,I6,I4,I5,2I4,F5.1,9I6,...` ⇒ 1-indexed char columns
#   year[1:4] age[5:8] TPA[9:14] BA[15:18] SDI[19:23] CCF[24:27] TopHt[28:31] QMD[32:36]
#   TCuFt[37:42] MCuFt[43:48] SCuFt[49:54] BdFt[55:60].
# A whitespace split() MISPARSES whenever a field fills its width and abuts its neighbour — notably
# CCF (I4) hits 4 digits at ≥1000 and merges with the preceding SDI (I5), e.g. "  762"+"1019" → "7621019".
# That shift fabricated multi-thousand-percent SDI/CCF/QMD "divergences" (and false structure_densephase
# labels) whenever CCF straddled 1000 differently between live and jl. Parse by fixed columns instead.
const _SUM10_COLS = ((9,14),(15,18),(19,23),(24,27),(28,31),(32,36),(37,42),(43,48),(49,54),(55,60))
function parse_sum10(text)
    rows = Tuple{Int,Vector{Float64}}[]
    for ln in split(text, '\n')
        s = rstrip(ln); length(s) < 60 && continue          # leading cols are position-critical: rstrip only
        y = tryparse(Int, strip(s[1:4])); (y === nothing || y < 1000 || y > 3000) && continue
        vals = Float64[]; ok = true
        for (a,b) in _SUM10_COLS
            v = tryparse(Float64, strip(s[a:b])); v === nothing && (ok=false; break); push!(vals, v)
        end
        ok && push!(rows, (y, vals))
    end
    rows
end

# Returns (sum_text, crashed). crashed=true when the live binary died on a SIGNAL (SIGFPE=8/SIGSEGV=11/SIGABRT=6)
# or exit>128 — distinguishing a live-FVS CRASH (e.g. the FVS40 >1000-TPA-seedling floating-point exception, which
# FVSjl survives) from a clean no-output run. Lets the ledger record `live_crash` instead of silently skipping.
function run_live(bin, cn, db, regime, dir, plantyr=0, seed=nothing)
    key = joinpath(dir,"s.key"); write(key, keytext(cn, db, regime, plantyr, seed))
    for f in ("s.sum","s.out"); fp=joinpath(dir,f); isfile(fp) && rm(fp); end
    crashed = false
    try
        p = run(pipeline(ignorestatus(`$bin --keywordfile=$key`); stdout=devnull, stderr=devnull))
        crashed = (p.termsignal != 0) || (p.exitcode > 128)
    catch
        crashed = true
    end
    sp = joinpath(dir,"s.sum"); (isfile(sp) ? read(sp,String) : "", crashed)
end

# Build a temp indexed sub-DB from the master for exactly `cns` (C-speed; master never modified).
function build_subdb(cns, out)
    inlist = join(["'" * replace(s,"'"=>"''") * "'" for s in cns], ",")
    isfile(out) && rm(out)
    dst = SQLite.DB(out); DBInterface.execute(dst, "ATTACH DATABASE '$MASTER' AS m")
    for tbl in ("FVS_STANDINIT_COND","FVS_TREEINIT_COND")
        DBInterface.execute(dst, "CREATE TABLE $tbl AS SELECT * FROM m.$tbl WHERE STAND_CN IN ($inlist)")
        DBInterface.execute(dst, "CREATE INDEX ix_$(tbl)_cn ON $tbl(STAND_CN)")
    end
    DBInterface.execute(dst, "DETACH DATABASE m"); SQLite.close(dst)
end

# Deterministic signature over the MEASURED facts (first match wins). UNCLASSIFIED = needs a manual trace.
# `struct_mat`/`vol_mat`/`density_be` are MATERIAL flags (divergence beyond a ±1-unit / MATERIAL-rel straddle).
function classify(bit_exact, struct_mat, vol_mat, density_mat, converges, max_rel, max_abs)
    bit_exact && return "bit_exact"
    # every diverging cell is a ±1-unit / sub-MATERIAL straddle ⇒ print/ULP boundary
    (!struct_mat && !vol_mat) && return "print_boundary"
    # structure moves materially but ONLY in TPA/QMD (density BA/SDI/CCF/TopHt preserved).
    # ⚠ This is NOT self-evidently a benign self-thin realization straddle: it is the IDENTICAL signature of a
    # REAL deterministic count/height bug. When a bug redistributes TPA across size-classes (over-/under-kill in
    # density-dependent mortality, or a stunted small-tree height cohort that never crosses the diameter gate),
    # BA/SDI/CCF/TopHt CONVERGE BY CONSTRUCTION (the model agrees on stand density, disagrees on HOW MANY trees
    # carry it), so "density preserved + TPA/QMD move" cannot on its own distinguish RNG straddle from bug. The
    # 2026-09-04 IE re-diagnosis showed this pattern auto-cornered SEED-INVARIANT deterministic bugs (TPA −574,
    # −141, ±100, seed-invariant one-directional) as if they were #206 OLDRN/ZRAND straddles. So we emit a
    # signature that does NOT auto-drop: material count divergences ESCALATE to needs_dig (see sweep_db.jl /
    # filter_digworthy.jl) and MUST be seed-tested (test/harness/fia/seed_test.jl) before any corner — a seed
    # sweep that FLIPS the signed jl−oracle residual is a genuine straddle; a SEED-INVARIANT one-directional
    # residual is a real deterministic bug.
    (struct_mat && !density_mat) && return "count_divergence_UNVERIFIED"
    # structure clean (only ±1 straddles) but volume moves materially ⇒ merch/threshold crossing
    (!struct_mat && vol_mat && converges)  && return "threshold_crossing"
    (!struct_mat && vol_mat && !converges) && return "volume_persistent"        # FLAG: volume-only, no convergence
    # density itself moves materially (DGSCOR record-ordering / LS dense-phase growth-ranking)
    (struct_mat && density_mat) && return "structure_densephase"
    return "UNCLASSIFIED"                                                        # FLAG: manual both-sides trace
end

function main(listfile, v, regime)
    bin = BIN[v]; var = VAR[v]
    cns = [split(strip(l),'\t')[1] for l in eachline(listfile) if !isempty(strip(l))]
    dir = mktempdir(); sub = joinpath(dir, "sub.db")
    print(stderr, "building sub-DB ($(length(cns)) stands)..."); flush(stderr)
    build_subdb(cns, sub); println(stderr, "ok")
    # Per-stand INV_YEAR so the PLANT regime uses a CALENDAR date at a CYCLE BOUNDARY (the faithful form, audit
    # 42d-42i). Offset +10 lands on a boundary for BOTH default cycle lengths: cycle-1 for 10-yr variants
    # (NE/CS/LS) and cycle-2 for 5-yr SN. A MID-CYCLE date (e.g. inv+5 on a 10-yr variant) instead triggers jl's
    # sub-cycle establishment-date age residual on every stand (the 42h "NE/CS/LS 0%" was that harness mistake).
    period = 10
    invyr = Dict{String,Int}()
    let db = SQLite.DB(sub)
        for r in DBInterface.execute(db, "SELECT STAND_CN, INV_YEAR FROM FVS_STANDINIT_COND")
            (r[:STAND_CN] === missing || r[:INV_YEAR] === missing) && continue
            invyr[String(r[:STAND_CN])] = Int(r[:INV_YEAR])
        end
        SQLite.close(db)
    end
    plantyr_of(cn) = (regime == "plant" && haskey(invyr, cn)) ? invyr[cn] + period : 0
    out = get(ENV, "LEDGER", "docs/fia_ledger.csv")
    newfile = !isfile(out)
    io = open(out, "a")
    newfile && println(io, "variant,regime,cn,n_cycles,bit_exact,div_cols,worst_col,worst_cycle,max_rel_pct,max_abs_diff,struct_max_rel_pct,vol_max_rel_pct,density_bitexact,converges,signature,struct_max_abs,vol_max_abs")
    n=0; nbe=0; ndiv=0; nskip=0; ncrash=0
    # optional durable per-stand coverage record: set SWEEP_DB to the local SQLite path (survives sessions /
    # container restart) and every stand's outcome (bit_exact | ulp_class | needs_dig) is upserted as we go.
    sdb = haskey(ENV, "SWEEP_DB") ? open_sweepdb(ENV["SWEEP_DB"]) : nothing
    # SKIP_DONE: skip CNs already recorded (bit_exact/ulp_class/live_crash) so a re-sweep of a covered range only
    # actually runs the UNCOVERED stands (cheap backfill — the historical skips, mostly fast live crashes). Preload
    # the done-set once (indexed lookup).
    ndone = 0; done_cns = Set{String}()
    if haskey(ENV, "SKIP_DONE") && sdb !== nothing
        for r in DBInterface.execute(sdb, "SELECT cn FROM sweep WHERE variant='$v' AND dig_class IN ('bit_exact','ulp_class','live_crash')")
            push!(done_cns, String(r.cn))
        end
    end
    # a per-cell divergence is MATERIAL if rel≥MATERIAL AND abs>1 (beyond a ±1-unit / sub-% straddle)
    ismat(lv,jv) = (ad=abs(lv-jv); ad > 1.0 + 1e-6 && (lv==0 ? true : ad/abs(lv) >= MATERIAL))
    for cn in cns
        n += 1; n % 200 == 0 && (print(stderr, "[$n/$(length(cns))] be=$nbe div=$ndiv\r"); flush(stderr))
        cn in done_cns && (ndone += 1; continue)     # SKIP_DONE: already recorded, don't re-run FVS
        py = plantyr_of(cn)
        live, live_crashed = run_live(bin, cn, sub, regime, dir, py)
        keyf = joinpath(dir,"jl.key"); write(keyf, keytext(cn, sub, regime, py))
        jlout = try FVSjl.run_keyfile(keyf; variant=var) catch; ""; end
        L = parse_sum10(live); J = parse_sum10(jlout)
        if isempty(L)
            # Live emitted no comparable .sum. If it CRASHED (SIGFPE etc.) but FVSjl projected fine, this is an
            # FVS-UB (live-binary bug) that jl SURVIVES — record it as `live_crash` so coverage is honestly
            # accounted (comparable + live_crash + skip), not silently dropped. Otherwise a genuine skip.
            if live_crashed && !isempty(J)
                ncrash += 1
                println(io, join([v, regime, cn, length(J), false, "", "", 0, 0.0, 0.0, 0.0, 0.0,
                                  false, false, "live_crash", 0.0, 0.0], ","))
                flush(io)
                if sdb !== nothing
                    try
                        upsert!(sdb, (variant=v, cn=cn, regime=regime, n_cycles=length(J), bit_exact=false,
                                      div_cols="", worst_col="", worst_cycle=0, max_rel_pct=0.0, max_abs_diff=0.0,
                                      struct_max_rel_pct=0.0, vol_max_rel_pct=0.0, struct_max_abs=0.0, vol_max_abs=0.0,
                                      density_bitexact=false, converges=false, signature="live_crash"))
                    catch e; print(stderr, "sweep_db live_crash upsert failed $cn: $e\n"); end
                end
            else
                nskip += 1
            end
            continue
        end
        isempty(J) && (nskip+=1; continue)
        Jd = Dict(y=>vv for (y,vv) in J)
        div_set = Set{Int}(); worst_rel=0.0; worst_col=0; worst_yr=0; max_abs=0.0
        struct_rel=0.0; vol_rel=0.0; struct_mat=false; vol_mat=false; density_mat=false
        struct_abs=0.0   # largest ABSOLUTE diff among structure cols 1-6 (escalation floor: separates a real
                         # BA/SDI/CCF move of 10s of units from a small-base ±1-5 ULP straddle that inflates to a
                         # big RELATIVE % only because the base is tiny — e.g. young age-3 stand, BA 2 vs 3 = 33%)
        vol_abs=0.0      # largest ABSOLUTE diff of the TCuFt column (7) SPECIFICALLY — the floor for the TCuFt
                         # escalation net. Must be TCuFt-only, NOT all vol cols: BdFt (board feet) magnitudes are
                         # ~10x cubic feet and would dominate, making the floor meaningless. A real volume-equation
                         # bug (FORKOD zero-vol) moves 1000s of cuft; a 35% on a young 310-cuft stand (109 cuft) is
                         # small-base inflation, not a bug. (Mirror of struct_abs.)
        peak_rel=0.0; last_rel=0.0; ncyc=0
        for (y, lv) in sort(collect(L))
            haskey(Jd, y) || continue
            jv = Jd[y]; ncyc += 1; cyc_rel = 0.0
            for k in 1:10
                if lv[k] != jv[k]
                    push!(div_set, k)
                    ad = abs(lv[k]-jv[k]); ad > max_abs && (max_abs = ad)
                    rel = lv[k]==0 ? (jv[k]==0 ? 0.0 : 1.0) : abs(lv[k]-jv[k])/abs(lv[k])
                    mat = ismat(lv[k], jv[k])
                    k <= 6 && ad > struct_abs && (struct_abs = ad)
                    k == 7 && ad > vol_abs && (vol_abs = ad)   # TCuFt column ONLY (see vol_abs comment)
                    if k <= 6; rel > struct_rel && (struct_rel = rel); mat && (struct_mat = true)
                    else;      rel > vol_rel && (vol_rel = rel);       mat && (vol_mat = true); end
                    (k in DENSITY_COLS) && mat && (density_mat = true)
                    rel > cyc_rel && (cyc_rel = rel)
                    if rel > worst_rel; worst_rel=rel; worst_col=k; worst_yr=y; end
                end
            end
            cyc_rel > peak_rel && (peak_rel = cyc_rel); last_rel = cyc_rel
        end
        bit_exact = isempty(div_set)
        converges = bit_exact ? true : (last_rel == 0.0 || last_rel < 0.5*peak_rel)
        density_be = !density_mat
        sig = classify(bit_exact, struct_mat, vol_mat, density_mat, converges, worst_rel, max_abs)
        bit_exact ? (nbe+=1) : (ndiv+=1)
        dcols = join([COLS[k] for k in sort(collect(div_set))], "|")
        wcol = worst_col==0 ? "" : COLS[worst_col]
        println(io, join([v, regime, cn, ncyc, bit_exact, dcols, wcol, worst_yr,
                          round(worst_rel*100,digits=3), round(max_abs,digits=3),
                          round(struct_rel*100,digits=3), round(vol_rel*100,digits=3),
                          density_be, converges, sig, round(struct_abs,digits=3), round(vol_abs,digits=3)], ","))
        flush(io)
        if sdb !== nothing
            try
                upsert!(sdb, (variant=v, cn=cn, regime=regime, n_cycles=ncyc, bit_exact=bit_exact,
                              div_cols=dcols, worst_col=wcol, worst_cycle=worst_yr,
                              max_rel_pct=round(worst_rel*100,digits=3), max_abs_diff=round(max_abs,digits=3),
                              struct_max_rel_pct=round(struct_rel*100,digits=3),
                              vol_max_rel_pct=round(vol_rel*100,digits=3),
                              struct_max_abs=round(struct_abs,digits=3), vol_max_abs=round(vol_abs,digits=3),
                              density_bitexact=density_be, converges=converges, signature=sig))
            catch e
                print(stderr, "sweep_db upsert failed for $cn: $e\n")   # never let the DB write break the sweep
            end
        end
    end
    close(io)
    sdb !== nothing && SQLite.close(sdb)
    println(stderr)
    println("LEDGER $v/$regime → $out : bit_exact=$nbe  diverging=$ndiv  live_crash=$ncrash  skipped(no-both-sum)=$nskip  of $n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: LEDGER=out.csv ledger_fia.jl <standlist> <variant> [regime]")
    main(ARGS[1], ARGS[2], length(ARGS) >= 3 ? ARGS[3] : "none")
end
