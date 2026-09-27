# test_em_autoes_tiered.jl — EM AUTOES establishment / birth-cycle REGENT per record vs the LIVE FVSem_g16 tiered goldens
# (test/fixtures/tiered/em: stands.db + <stand>_none.key + FVS_TreeList.csv written by the live oracle).
#
# Stand 196378260020004 (4 inventory points, PSLO .20/.60/.25/.35) books its first AUTOES ingrowth cohort at 2031
# (records ES020370-426). estb/estab.f:474-479 sets SLO=PSLO(NNID), XCOS/XSIN=COS/SIN(PASP(NNID))·SLO once per
# inventory point, and ESTPP, ESNSPE, ESPADV/ESPSUB/ESPXCS and ESADVH/ESSUBH all read those /ESCOMN/ values. jl fed only
# ESTPP the per-point slope, so the point-4 plots (SLO .35) drew their species from the point-1 (SLO .20) logits:
# live FVSem_g16 books DF on plot 41, LP on 43, PP on 49; jl booked WL, DF, LP.
module EMAutoesTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

# Rows of a table keyed by (Year, TreeId, TreeIndex) => Dict(col => value string).
function _keyed(hdr, rows)
    iy = findfirst(==("Year"), hdr); it = findfirst(==("TreeId"), hdr); ix = findfirst(==("TreeIndex"), hdr)
    Dict((r[iy], strip(r[it]), r[ix]) => Dict(hdr[k] => r[k] for k in eachindex(hdr)) for r in rows)
end
_eq(a, b) = a == b || (let x = tryparse(Float64, a), y = tryparse(Float64, b); x !== nothing && x == y end)

const STAND = "196378260020004"
function _run_stand(cn)
    d = mktempdir()
    _, db, crashed, err = run_case("EM", cn, "none"; dir = d)
    (db = db, crashed = crashed, err = err)
end
const _run = _run_stand(STAND)
# Number of FVS_TreeList cells (every column but the keys) that differ from the live golden in the given years.
function _treelist_diffcells(cn, db, years; cols = nothing)
    gold = _keyed(read_csv(joinpath(fixture_dir("EM"), "$(cn)_none.FVS_TreeList.csv"))...)
    got  = _keyed(db_table_rows(db, "FVS_TreeList")...)
    n = 0
    for (k, g) in gold
        k[1] in years || continue
        haskey(got, k) || (n += 1; continue)
        for (c, v) in g
            c in ("StandID", "CaseID") && continue
            cols === nothing || c in cols || continue
            _eq(v, get(got[k], c, "")) || (n += 1)
        end
    end
    n
end

@testset "EM AUTOES per-point topography (estab.f:474-479) — 196378260020004 cohort 2031 vs FVSem_g16" begin
    @test !_run.crashed
    gold = _keyed(read_csv(joinpath(fixture_dir("EM"), "$(STAND)_none.FVS_TreeList.csv"))...)
    got  = _keyed(db_table_rows(_run.db, "FVS_TreeList")...)
    es = sort([k for k in keys(gold) if k[1] == "2032" && startswith(k[2], "ES")]; by = k -> parse(Int, k[3]))
    @test length(es) == 57                                       # ES020370..426
    @test all(haskey(got, k) for k in es)
    # the booked cohort: species and inventory point per record (live-exact)
    for c in ("SpeciesFVS", "PtIndex")
        @test count(k -> haskey(got, k) && !_eq(gold[k][c], got[k][c]), es) == 0
    end
end

# em/esgent.f → REGENT(LESTB) runs INSIDE ESTAB (estab.f:1493), before gradd.f:244's post-establishment DENSE, so the
# birth-cycle SMHTGF/SMDGF read RELDEN/BA/PCCF from gradd.f:192's DENSE (post-growth, pre-ESNUTR). jl re-DENSEd with the
# new cohort first: live RELDEN 109.1242 / point-1 PCCF 88.5849 vs jl 109.1673 / 88.6849 ⇒ every birth-cycle HTGRR
# ~2e-4 low (DF 370: live 2.2399454, jl 2.23945). The cohort's heights/diameters after ESGENT are live-exact.
@testset "EM ESGENT reads the pre-ESNUTR DENSE (gradd.f:192) — 196378260020004 cohort 2031 vs FVSem_g16" begin
    gold = _keyed(read_csv(joinpath(fixture_dir("EM"), "$(STAND)_none.FVS_TreeList.csv"))...)
    got  = _keyed(db_table_rows(_run.db, "FVS_TreeList")...)
    es = [k for k in keys(gold) if k[1] == "2032" && startswith(k[2], "ES")]
    for c in ("Ht", "HtG", "DBH", "DG", "PctCr")
        @test count(k -> haskey(got, k) && !_eq(gold[k][c], got[k][c]), es) == 0
    end
end

# EM growth kernels at gfortran single precision: em/htgf.f + pothtg.f EXP/ALOG/** are glibc expf/logf/powf (fexp/flog/
# fpow), and the RALPH term folds `-0.1*18.158` to 1.8158001 at compile time; em/dgf.f DGCONS SIN/COS (sinf/cosf), ALOG,
# DGPCCF*(RELDEN**2) and the LM `.01*(-.199592)*RELDEN` second term in Fortran order; dense.f (the IE/BM file) sums BA/
# TSUMD2 as D*(D*P) and RELDEN as species-major RELDSP subtotals in IND1 order. MEASURED FVSem_g16 196378260020004
# cycle 1: RELDEN 42B41C64 (jl record-order 42B41C67) ⇒ BAL/DDS 1-2 ULP ⇒ 16 DG/15 HtG/15 DBH cells off at 2022.
@testset "EM growth kernels at gfortran precision (htgf/pothtg/dgf/dense) vs FVSem_g16" begin
    @test _treelist_diffcells(STAND, _run.db, ("2012", "2022")) == 0
    r2 = _run_stand("3087467010690")
    @test !r2.crashed
    @test _treelist_diffcells("3087467010690", r2.db, ("1988", "1998", "2008")) == 0
end

# em/morts.f label 10: the QMD-convergence loop (≤10 passes; D10=D10N while |D10−D10N|>0.1 and D10N>DIA0), each pass
# re-deriving BA10/RZ and TN10 from the current D10 — jl ran ONE pass. Plus the Fortran shape: stand sums in IND1 order
# with G=(DG/BARK)*(FINT/10), T capped at 35000, AVED over all records, POT DATA literals, EXP/** via gfortran libm, the
# X*0.6 NI-share order, SIZCAP floor. MEASURED FVSem_g16 DEBUG MORTS 2999215010690 cycle 1: DQ10 7.1664 → D10N 7.3036
# (2nd pass), jl killed ~0.1% fewer CW/AS; 3087467010690 cycle 3 whole-stand TPA.
@testset "EM MORTS QMD-convergence passes (em/morts.f label 10) vs FVSem_g16" begin
    r3 = _run_stand("3087467010690")
    @test _treelist_diffcells("3087467010690", r3.db, ("1988", "1998", "2008", "2018", "2028", "2038")) == 0
    r4 = _run_stand("2999215010690")
    @test _treelist_diffcells("2999215010690", r4.db, ("1998",); cols = ("TPA", "MortPA", "DBH", "Ht", "DG", "HtG")) == 0
end
end # module
