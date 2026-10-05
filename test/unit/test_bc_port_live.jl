# test_bc_port_live.jl — BC (British Columbia, metric) fixes vs the LIVE BC tiered goldens (test/fixtures/tiered/bc: stands.db +
# <stand>_<regime>.key + the .sum / *_Metric DBS tables of the private DB-capable oracle /workspace/.bcwork/dbfix/FVSbc_dbfix,
# ORACLE_SOURCE_AUDIT §7). Each testset names the Fortran it follows and the case it was measured on.
module BcPortLiveTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

function _case(cn, r)
    d = mktempdir()
    txt, db, crashed, err = run_case("BC", cn, r; dir = d)
    return (txt = txt, db = db, crashed = crashed, err = err,
            ms = crashed ? Mismatch[] : compare_case("BC", cn, r, txt, db))
end
_in(ms, files) = [m for m in ms if m.file in files]

# BC FFE (canada/fire/bc fm*.f + the shared fire/base, fire/vbase, fire/ie layers it links, metric/fire I/O): SIMFIRE and
# SALVAGE ran to an "FFE fuel tables not ported" error. Ported: bc/fmvinit.f species/decay/snag parameters, bc/fmcba.f (the
# IE FMCBA on 15-species FULIVE/FULIVI/FUINIE/FUINII), bc/fmbrkt.f, bc/fmcrow.f ISPMAP (+FMCROWE for 11/12/13/15),
# bc/fmsvol.f (Kozak CFVOL + CFTOPK, no cone floor), metric fmin.f SIMFIRE km/h·°C and SALVAGE cm, ICMETRC=1, and the
# metric DBS (FVS_Carbon_Metric, FVS_Mortality_Metric; no BurnReport/PotFire — metric fmfout.f:96 / fmpofl.f:295 windows).
@testset "BC FFE SIMFIRE/SALVAGE (canada/fire/bc, metric/fire)" begin
    for r in ("simfire", "salvage")
        c = _case("YSM029-250", r)
        @test !c.crashed
        @test isempty(_in(c.ms, ("sum", "FVS_Summary_Metric", "FVS_Mortality_Metric", "FVS_Hrv_Carbon")))
        @test !any(m -> m.col == "PRESENCE", c.ms)          # the live table set: no PotFire/BurnReport/InvReference/Error
    end
end

# canada/bc/cfvol.f → min.f → log.f (Kozak 2002 taper), as FMSVOL calls it: LOG's REAL**REAL/ALOG/EXP are glibc
# powf/logf/expf; MIN converts CFVOL's imperial grinit.f thresholds back to metric (SH = 30·CMtoFT·FTtoM, TD =
# 10·CMtoIN·INtoCM — not 0.3 m / 10 cm, so a sub-0.3 m stump takes log.f's partial-first-log branch and DI3 there feeds the
# stump volume); cfvol.f:55 gates merch in INCHES (D ≥ DBHMIN = 12.5·CMtoIN catches the 12.5-cm PL). Replays 499 FMSVOL calls
# dumped from the private oracle (test/fixtures/bc_ffe/fmsvol_calls.txt): VN and VM bit-exact.
@testset "BC Kozak taper volume via FMSVOL (canada/bc/cfvol.f, min.f, log.f; canada/fire/bc/fmsvol.f)" begin
    h2f(x) = reinterpret(Float32, parse(UInt32, x; base = 16))
    n = 0; bad = 0
    for l in eachline(joinpath(@__DIR__, "..", "fixtures", "bc_ffe", "fmsvol_calls.txt"))
        startswith(l, "#") && continue
        f = split(l)
        sp = parse(Int, f[1]); d = h2f(f[2]); h = h2f(f[3]); xht = h2f(f[4])
        ltkil = f[8] == "T"
        vn, vm = FVSjl.bc_fmsvol(nothing, sp, d, h, ltkil ? xht : -1f0)
        n += 1
        (vn == h2f(f[6]) && vm == h2f(f[7])) || (bad += 1)
    end
    @test n == 499
    @test bad == 0
    c = _case("YSM029-250", "salvage")                     # live FVS_Carbon_Metric pools now exact too
    @test isempty(c.ms)
end

# dbstreesin.f reads every INTEGER tree column with fsql3_colint (dbsqlite/fvsqlite3.c:585): NULL ⇒ the default (PLOT_ID 0, not
# 1 — dbstreesin.f:53/93; intree.f:322 forces plot 1 only under IPTINV=1), a REAL cell C-truncated (CrRatio 71.523 ⇒ 71, not
# NINT 72). MEASURED FVSbc_dbfix Fir.20 (NULL Plot_ID, REAL CrRatio): FVS_TreeList_Metric PtIndex 0, inventory PctCr exact.
@testset "BC DB tree reader fsql3_colint semantics (dbstreesin.f, fvsqlite3.c:585)" begin
    c = _case("Fir.20", "none")
    @test !c.crashed
    @test !any(m -> m.year == "2020", c.ms)                # the whole inventory treelist (PtIndex, PctCr, CrWidth, …)
    @test !any(m -> m.col == "PtIndex", c.ms)
end

# canada/newmist/mistoe.f:168 NEWSI = NEWMOD .AND. MISFLG: the NISI spatial model (DMTREG/DMMDMR) runs only under NEWSPRED;
# a MISTOE block alone (with MISTPINF) keeps the base mistoe.f spread/MISINF/MISMRT. jl ran DMTREG on ms.active (any MISTOE)
# and the base spread, then published the DMTREG ratings over it. MEASURED FVSbc_dbfix mistletoe regime: the .sum and
# FVS_Summary_Metric of YSM029-250/-266 exact (SkyRanch-Control 102 → 22 cells; BC mistletoe regime ~990 → ~260).
@testset "BC DMTREG only under NEWSPRED (canada/newmist/mistoe.f:168 NEWSI)" begin
    for cn in ("YSM029-250", "YSM029-266")
        c = _case(cn, "mistletoe")
        @test !c.crashed
        @test isempty(_in(c.ms, ("sum", "FVS_Summary_Metric")))
    end
end

# The metric builds link metric/dbsqlite + metric/base/sstage.f + metric/newmist/misprt.f + rd/rdpr.f (LMTRIC):
# FVS_StrClass_Metric (DBHS·INtoCM, INT(ht·FTtoM)); FVS_DM_Stnd_Sum_Metric / FVS_DM_Spp_Sum_Metric (TPH = /ACRtoHA, BA, Vol
# per ha; the Stnd table keeps only its first row — DBSMIS2's `Commit` with no `Begin` fails ⇒ IDM2=0); FVS_RD_Sum_Metric
# (metric column names; Live_Merch_CuM sums WK1 = canada/bc/vols.f:167 merch cubic, US builds the DG snapshot); an EMPTY
# FVS_Climate_Metric (metric dbsclsum.f INSERTs into FVS_Climate ⇒ ICLIM=0). MEASURED FVSbc_dbfix tiered YSM029-250.
@testset "BC metric extension DBS tables (metric/dbsqlite, sstage.f, misprt.f, rdpr.f LMTRIC)" begin
    for r in ("none", "thinbba", "climate")
        c = _case("YSM029-250", r)
        @test !c.crashed
        @test isempty([m for m in c.ms if m.file != "FVS_TreeList_Metric"])   # (per-record HTG ULPs: open)
    end
    c = _case("YSM029-250", "rootdis")
    @test !any(m -> m.col == "PRESENCE" || (tryparse(Int, m.year) !== nothing && parse(Int, m.year) < 2058), c.ms)
    c = _case("YSM029-250", "mistletoe")
    @test !any(m -> m.col == "PRESENCE" || m.col == "ROWCOUNT" || m.file in ("sum", "FVS_Summary_Metric"), c.ms)
end

# canada/bc/crown.f:510-512 (V3) LSTART dub: a ≥2-cm record with no crown ratio gets ICRI=INT(BACHLO(ICRI,CRSD(ISPC),RANN))
# when DGSD ≥ 1, in the species-major IND1 order of DO 70/DO 60. jl drew nothing, so every later RANN draw (calibration
# DGDRIV, the cycle REGENT ZZRANs) landed on other records. MEASURED FVSbc_dbfix SkyRanch-Control: 207 CRATET draws before
# the calibration DGDRIV (61 records), the REGENT seed at every cycle now equals live's; 2029 HtG 2919 cells → 4.
@testset "BC LSTART crown dub random error (canada/bc/crown.f:510-512)" begin
    c = _case("SkyRanch-Control", "none")
    @test !c.crashed
    @test count(m -> m.year == "2029" && m.col == "HtG", c.ms) <= 4
end

end # module
