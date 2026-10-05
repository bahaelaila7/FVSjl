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

# canada/bc/regent.f:1508/1638 (tripling): D starts as DBH(I) and after a D<3-in pass the crown block re-reads D=DBH(K), so
# the next copy's XWT, D<3 test and DDS use the previous pass's DBH(K) — the record's own reset DBH (HK<4.5) or the DBH
# DGDRIV wrote into the copy slot (dgdriv.f:271/278). MEASURED FVSbc_dbfix SkyRanch-Control: the 6-7 cm / <1 m PL records
# (e.g. 923: copy 2823 XWT 0, keeps DBH 6 cm); 2029 now only one record's DG ULP away.
@testset "BC REGENT tripled copies inherit the previous pass's D (canada/bc/regent.f:1638)" begin
    c = _case("SkyRanch-Control", "none")
    @test count(m -> m.year == "2029" && m.file == "FVS_TreeList_Metric", c.ms) <= 3
end

# canada/bc/initre.f converts the metric keyword fields in place before the shared imperial processing — the same
# statement set as canada/on/initre.f (THINBBA ARRAY(2)·M2pHAtoFT2pACR :3602, DBH limits ·CMtoIN, heights ·MtoFT, …), which
# jl applied for ON only. MEASURED FVSbc_dbfix tiered thinbba (THINBBA 40 m²/ha): live removes nothing at 2028 on
# YSM029-265 (jl thinned to 40 ft²/ac, RTPA 103).
@testset "BC metric keyword fields (canada/bc/initre.f THINBBA etc.)" begin
    c = _case("YSM029-265", "thinbba")
    @test !c.crashed
    @test isempty(_in(c.ms, ("sum", "FVS_CutList_Metric")))
    @test !any(m -> m.col in ("RTpa", "RTPA", "RTph"), c.ms)
end

# canada/bc/crown.f cycling: BARK=BRATIO(ISPC,D,H) (:422, BC's BARK1 constant — jl used the generic linear form, which
# floors at 0.80) backdates D-DG/BARK, and the backdated BAL reads OLDPCT (cratet.f:531 OLDPCT=PCT before the LSTART
# CROWN) and OLDBA (dense.f:239-240, crown.f:344). MEASURED FVSbc_dbfix Fir.20 cycle-1 CROWN dump: backdated D 3.9409
# (= the 2020 DBH) live vs 3.8163 jl, OLDPCT 84.65 vs 0 ⇒ ICR off by 1-4 on 335 records; now none.
@testset "BC crown backdating bark + OLDPCT/OLDBA (canada/bc/crown.f:422-470, cratet.f:531)" begin
    c = _case("Fir.20", "none")
    @test !any(m -> m.year == "2030" && m.col in ("PctCr", "CrWidth", "DBH", "Ht"), c.ms)
end

# MISTOFF (mistoe/misin.f:446, metric/newmist/misin.f:480): MISFLG=.FALSE. ⇒ MISDAM loads no rating from the damage codes
# (misdam.f:64) and BC's NEWSI stays off. Unported in every variant. Live golden: the tiered SkyRanch-Control stand (damage-
# code mistletoe) with a MISTOE/MISTOFF/END block, private FVSbc_dbfix .sum (test/fixtures/bc_port).
@testset "MISTOFF keyword (misin.f MISFLG, misdam.f:64)" begin
    d = mktempdir(); fx = joinpath(@__DIR__, "..", "fixtures", "bc_port")
    cp(joinpath(@__DIR__, "..", "fixtures", "tiered", "bc", "stands.db"), joinpath(d, "stands.db"))
    key = replace(read(joinpath(fx, "skyranch_control_mistoff.key"), String),
                  "\nstands.db\n" => "\n" * joinpath(d, "stands.db") * "\n", "\nout.db\n" => "\n" * joinpath(d, "out.db") * "\n")
    write(joinpath(d, "s.key"), key)
    txt = FVSjl.run_keyfile(joinpath(d, "s.key"); variant = FVSjl.BritishColumbia())
    jl = [l[1:min(end, 120)] for l in split(txt, '\n') if length(l) > 20 && !startswith(l, "-999")]
    live = [l[1:min(end, 120)] for l in eachline(joinpath(fx, "skyranch_control_mistoff.live.sum")) if length(l) > 20]
    @test length(jl) == length(live) == 6
    @test jl == live
end

# The REAL*4 sums feeding BC's small-tree REGENT: RELDEN is base dense.f's species-major RELDSP over IND1 (dense.f:166-223)
# — BC's own flat record-order bc_stand_ccf was 1-2 ULP off — and the subcycle RDNEXT/BANEXT(J+1) accumulations walk
# DO 16 ISPC / DO 15 I3 / I=IND1(I3) (regent.f:1289-1321), not record order. MEASURED FVSbc_dbfix SkyRanch-Control
# cycle-2 REGENT dump: RELDEN 4135C00A live / 4135C073 jl, RDNEXT(2) 4084C978 / 4084C842 ⇒ every HTGRL 1-3 ULP low.
@testset "BC REGENT density sums in IND1 order (dense.f RELDSP, regent.f:1289-1458)" begin
    c = _case("SkyRanch-Control", "none")
    @test length(c.ms) < 100                      # 74,847 cells before
end

# canada/bc/sitset.f:260 BAMAX = BAMAX*M2pHAtoFT2pACR uses the METRIC.F77 PARAMETER 4.3560773; jl used 1/FT2pACRtoM2pHA
# (4.356078, 2 ULP high) ⇒ SDIDEF/XMAX/BAMAX 2 ULP high ⇒ the BAMAX-approach RIPP and every morts.f kill. And morts.f:444-445
# sums SD2SQ = SD2SQ + P*(D*D + CIOBDS), CIOBDS = (2*D*G + G*G) (jl added the three terms left to right). MEASURED FVSbc_dbfix
# Fir.20 MORTS dump: BAMAX 4336F48B live / 4336F48D jl ⇒ cycle-2 MortPH 5e-6 relative on 566 records.
@testset "BC mortality BAMAX constant + SD2SQ order (sitset.f:260, morts.f:444)" begin
    for (cn, cap) in (("Fir.20", 150), ("YSM029-246", 150), ("SkyRanch-0.3m", 10))
        c = _case(cn, "none")
        @test length(c.ms) <= cap
    end
end

# canada/bc/crown.f:471-472 V3 backdated BAL = (1-OLDPCT/100)*OLDBA uses OLDBA itself; the RDM1<100 fallback OBA=BA
# (crown.f:344-348) only feeds the V2 DCRCON/X1. jl backdated with OBA. MEASURED FVSbc_instr YSM029-271 cycle-1 CROWN dump
# (RELDM1<100): backdated BAL 41ED1F47 live / 3x larger jl ⇒ EXPDCR off on 66/66 records ⇒ 2028 PctCr 1 low on 29.
@testset "BC backdated crown BAL uses OLDBA (canada/bc/crown.f:472)" begin
    c = _case("YSM029-271", "none")
    @test !any(m -> m.year == "2028" && m.col == "PctCr", c.ms)
    @test length(c.ms) < 200                      # 13,304 cells before
end

# canada/bc/crown.f:606-628 DO 79: cycle-0 dead records with no crown are dubbed — V3 DUBSCR = CRNMD(…,YSD=0) on
# BAL = (1-PCT/100)*OLDBA, XCRCON from RELDM1, both the dense.f:259-264 FINTH/FINT interpolations of cratet's backdating
# DENSE. jl never ran it for BC. MEASURED FVSbc_instr YSM029-250 DO 79 dump: P 10.95/6.60/25.17, OLDBA 41FF958F,
# RELDM1 4216E1D6 ⇒ ICR 95 on all 3 dead records.
@testset "BC cycle-0 dead-record crown dub (canada/bc/crown.f:606-628)" begin
    c = _case("YSM029-250", "none")
    @test !any(m -> m.year == "2018" && m.col == "PctCr", c.ms)
end

end # module
