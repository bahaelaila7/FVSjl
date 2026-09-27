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
const _run = let d = mktempdir()
    _, db, crashed, err = run_case("EM", STAND, "none"; dir = d)
    (db = db, crashed = crashed, err = err)
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
end # module
