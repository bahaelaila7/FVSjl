# test_west_shared_tiered.jl — shared western fixes vs the LIVE western tiered goldens (test/fixtures/tiered/<v>: stands.db +
# <stand>_<regime>.key + the live .sum / DBS tables written by FVS<v>_g16 / FVScr_clean). Each testset names the Fortran
# it follows and the measured case; cells are compared at a 1e-5 relative tolerance (Float32 drift from other open
# residuals is not what these tests pin — the named mechanism is).
module WestSharedTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

function _case(v, cn, r)
    d = mktempdir()
    txt, db, crashed, _ = run_case(v, cn, r; dir = d)
    return (txt = txt, db = db, crashed = crashed, ms = compare_case(v, cn, r, txt, db))
end
_rel(m) = (a = tryparse(Float64, m.gold); b = tryparse(Float64, m.got);
           (a === nothing || b === nothing) ? Inf : abs(a - b) / max(abs(a), 1e-9))
_material(ms) = [m for m in ms if _rel(m) > 1e-5]

# tt/fmcrow.f and ut/fmcrow.f call FMCROWW(SPIE=ISPMAP(SPIW)) for the conifers and FMCROWE only for the hardwoods
# (TT CASE(6,15,16,18), UT CASE(6,18:20,22)); both fmcroww.f/fmcrowe.f are byte-identical to CR's. jl had no TT/UT route,
# so every species took the eastern FMCROWE with `ls_spi` (MEASURED FVStt_g16 2780339010690 2000, the inventory year:
# Aboveground_Total_Live 13.17 live / 14.30 jl, Canopy_Density 0.102 / 0.061 — 196 inventory-year FFE cells on TT).
@testset "TT/UT crown weights route FMCROWW by ISPMAP (tt|ut/fmcrow.f) vs live at the inventory year" begin
    for (v, cn, y) in (("TT", "2780339010690", "2000"), ("TT", "2839796010690", "1999"),
                       ("UT", "42642675010690", "2010"), ("UT", "2390216010690", "2002"))
        c = _case(v, cn, "salvage")
        @test !c.crashed
        @test isempty(_material([m for m in c.ms if m.year == y && m.file in ("FVS_Carbon", "FVS_PotFire")]))
    end
end

# fmcba.f bare-stand cover type: with no basal area the first FFE year takes COVTYP = COVINI(ITYPE) (the habitat's seral
# cover) else the variant's "NO VALID HABITAT" default (tt/ut 7 LP, wc/pn COVINI6→16, ec →3); the CA-FFE top-2 variants
# then load the one cover at weight 1 (nc:359-360 COVCA(1)=COVTYP, COVCAWT(1)=1). jl used DF 3 for TT/UT, 16 for every
# WC/PN habitat, and left NC's COVCAWT at 0 (MEASURED: TT 3333677010690 1992 live COVTYP 7 — DDW 4.05 / jl 2.40; NC
# 450603388489998 2016 DDW 1.40 / jl 0.0; WC 374601116489998 COVTYP 19; PN 26379272010900 19; EC 450507010497 10).
@testset "bare-stand FFE cover type COVINI(ITYPE) + top-2 COVCA (fmcba.f) vs live at the first FFE year" begin
    for (v, cn, y) in (("TT", "3333677010690", "1992"), ("UT", "471768556489998", "2016"), ("NC", "450603388489998", "2016"),
                       ("WC", "374601116489998", "2015"), ("PN", "26379272010900", "2002"), ("EC", "450507010497", "2009"))
        c = _case(v, cn, "salvage")
        @test !c.crashed
        @test isempty(_material([m for m in c.ms if m.year == y && m.file == "FVS_Carbon"]))
    end
end

# ci/habtyp.f:44-90: PVREF4 crosswalks (PV_CODE, PV_REF_CODE) → KODTYP; an unknown pair (FVS34/33/32) or a KODTYP outside
# 10..999 (FVS14) keeps ci/grinit.f ICINDX=21 / ITYPE=4 (habitat 260). jl used PV_CODE mod 1000 and ICINDX 1 (MEASURED
# FVSci_g16: 3369538010690 9999999/491 and 12276084010690 45101/494 "MAPPED TO 260"; jl 999 / 101 ⇒ DGHAB/ITYPE wrong,
# e.g. 3369538010690 2006 HtG 11.67 live / 8.64 jl, BdFt 1246 / 1201).
@testset "CI habitat via PVREF4 + ICINDX 21 default (ci/habtyp.f, ci/pvref4.f) vs FVSci_g16" begin
    for cn in ("3369538010690", "12276084010690")
        c = _case("CI", cn, "none")
        @test !c.crashed
        @test count(m -> m.file == "sum", c.ms) == 0
        @test count(m -> m.file == "FVS_Error", c.ms) == 0          # FVS33 + FVS32 (unknown pair)
    end
    c = _case("CI", "3159852010690", "salvage")                        # bare, no PV code: FVS14 + CIPVG(ICINDX 21)
    @test count(m -> m.file in ("FVS_Error", "FVS_PotFire", "FVS_Carbon"), c.ms) == 0
end

# fortyp.f:1117-1140 California mixed-conifer test (ISTATE 6 or region 5): DF off the north coast (needs ICNTY, read
# from the FIA COUNTY column — dbsstandin.f:392-395), SP/IC, PP/JP with PP < 80%, WF/RF with true fir < 80% → 371
# (MEASURED: NC 44777283020004 2011 live 371 / jl 221, SO 7690240010901 371/261, CA 23742358010900 371/201, WS
# 850400255290487 371/222).
@testset "FORTYP California mixed conifer 371 (fortyp.f:1117-1140) vs live" begin
    for (v, cn) in (("NC", "44777283020004"), ("SO", "7690240010901"), ("CA", "23742358010900"), ("WS", "850400255290487"))
        c = _case(v, cn, "none")
        @test !c.crashed
        @test count(m -> m.col == "ForTyp", c.ms) == 0
    end
end

# dbsreference.f prints VEQNNC/VEQNNB (VOLEQDEF) and CYCLEs a blank JSP. NC/SO/CA/WS/AK keep their equations in the
# volume driver's per-forest tables (species.vol_eq blank ⇒ jl wrote '' — MEASURED FVSnc_g16 450603388489998 PP
# 500WO2W122); AK FORST '04' (Chugach) takes R10_EQN CHUEQN; WC/PN wrote rows for their '__' placeholder slots.
@testset "FVS_InvReference VolEq ids + blank-JSP rows (dbsreference.f, voleqdef.f) vs live" begin
    for (v, cn) in (("NC", "450603388489998"), ("SO", "24079993010900"), ("SO", "374286168489998"), ("CA", "23742358010900"),
                    ("WS", "850400255290487"), ("AK", "10709171010497"), ("AK", "10706662010497"), ("WC", "30194620010497"),
                    ("PN", "26379272010900"))
        c = _case(v, cn, "none")
        @test !c.crashed
        @test count(m -> m.file == "FVS_InvReference", c.ms) == 0
    end
end

# rdpr.f:79 IF (ITRN .EQ. 0) RETURN — a stand with no tree records writes no RDPR row (MEASURED FVStt_g16 3333677010690
# rootdis: no FVS_RD_Sum table; jl wrote 6 rows).
@testset "WRD: no FVS_RD_Sum row without tree records (rdpr.f:79) vs live" begin
    c = _case("TT", "3333677010690", "rootdis")
    @test count(m -> m.file == "FVS_RD_Sum", c.ms) == 0
end

# so/ccfcal.f: SH/WO CCF from R5CRWD on the Region-5 forests (IFOR 4-9) (MEASURED FVSso_g16 15364795010497, forest 514:
# inventory CCF 184 live / 207 jl).
@testset "SO SH/WO CCF by R5CRWD on R5 forests (so/ccfcal.f) vs FVSso_g16" begin
    c = _case("SO", "15364795010497", "none")
    @test count(m -> m.file == "sum", c.ms) == 0
end

# habtyp.f ERRGRO 14/32/33/34 + FORKOD 03 for the western variants (TT/UT PVREF4+CRDECD, WC/PN/EC PVREF6+HBDECD, NC/SO/CA
# R5HABT/PCOML + PVREF5/6, WS R5HABT+PVREF5), each measured against the live FVS_Error rows.
@testset "western habtyp / forkod FVS_Error rows (FVS03/14/32/33/34) vs live" begin
    for (v, cn) in (("TT", "1589567497290487"), ("UT", "434219452489998"), ("WC", "1127545052290487"), ("WC", "25075210010900"),
                    ("PN", "1166897755290487"), ("EC", "30193987010497"), ("SO", "15184869010497"), ("SO", "374286168489998"),
                    ("NC", "30192555010497"), ("CA", "22960323010497"), ("WS", "15353585010497"), ("CI", "3261005010690"),
                    ("KT", "22404917010497"))
        c = _case(v, cn, "none")
        @test count(m -> m.file == "FVS_Error", c.ms) == 0
    end
end

# fmcfmd.f treeless calls reuse the previous cover metagroup(s): CR/TT/UT OLDICT (fmvinit 6), WC/PN OLDICT/OLDICT2/
# OLDICTWT (fmvinit 2/0/1,0) — MEASURED FVStt_g16 3333677010690 1992 FMD 2 (jl 8), FVSpn_g16 26379272010900 2002 FMD 5 (jl 10).
@testset "bare-stand fuel model: FMCFMD OLDICT reuse (cr|tt|ut|wc|pn/fmcfmd.f) vs live" begin
    for (v, cn, y) in (("TT", "3333677010690", "1992"), ("CR", "103573371010661", "2006"), ("PN", "26379272010900", "2002"),
                       ("WC", "1285602241290487", "2022"))
        c = _case(v, cn, "salvage")
        @test isempty(_material([m for m in c.ms if m.year == y && m.file == "FVS_PotFire"]))
    end
end

# so/ and nc/ fmsvol.f & co. are the shared western FMSVOL layer: snag bole = NATCRS TCF, live merch stem = NATCRS MCF
# (MEASURED FVSso_g16 374286168489998 2015: Standing_Dead 2.91 live / 0.05 jl, Aboveground_Merch_Live 34.43 / 35.07).
@testset "SO FFE on the western FMSVOL layer (fmsvol.f) vs FVSso_g16 at the inventory year" begin
    for (cn, y) in (("374286168489998", "2015"), ("15184869010497", "2010"))
        c = _case("SO", cn, "salvage")
        @test isempty(_material([m for m in c.ms if m.year == y && m.file == "FVS_Carbon"]))
    end
end

# cr/regent.f:343-346 HK ≤ 4.5 ⇒ DG = 0 and DBH = D + 0.001·HK at once, so MORTS reads the bumped DBH (MEASURED
# FVScr_clean 46279527020004: SUMDR0 16533.72 live, 16531.86 when jl carried the bump as an UPDATE-time DG ⇒ MortPA drift).
@testset "CR REGENT sub-4.5 ft DBH set before MORTS (cr/regent.f:343-346) vs FVScr_clean" begin
    c = _case("CR", "46279527020004", "none")
    @test isempty(_material([m for m in c.ms if m.file == "FVS_TreeList" && m.year in ("2003", "2013", "2021")]))
end

# tt/cratet.f:99-148 50-yr-base SITEAR conversion reads TEMCCF on the NOTRE-expanded PROB (MEASURED FVStt_g16 335002534489998:
# TEMCCF 152.30, LP 70 → 42.42; at INITRE time TEMCCF was the 125 floor ⇒ LP 43.99 ⇒ FVS_InvReference SiteIndex 44 vs 42).
@testset "TT CRATET site conversion at CRATET time (tt/cratet.f:99-148) vs FVStt_g16" begin
    for cn in ("335002534489998", "2750433010690")
        c = _case("TT", cn, "none")
        @test count(m -> m.file == "FVS_InvReference" && m.col == "SiteIndex", c.ms) == 0
    end
end

# ---------------------------------------------------------------------------------------------------------------------
# west-shared-3 (round 3): Float32-faithfulness clusters. These pin EXACT equality (no tolerance) — the mechanism is a
# 1-ULP order/primitive difference that cascades through record tripling, so a 1e-5 test would not see it.
_cells(c; pred = m -> true) = count(pred, c.ms)

# r4vol.f R4MATTAPER: CF0=.002727*(HT67*STUMPD**2+D67**2*THT) — HT67 times the SQUARE, REAL ** = powf (FVStt_g16 1998 TCuFt).
@testset "R4 Matney CF0 Float32 order (r4vol.f R4MATTAPER) vs FVStt_g16 inventory TCuFt" begin
    for (eq, d, h, live) in (("400MATW108", 10.9f0, 62f0, 18.452423095703125), ("400MATW108", 10.6f0, 59f0, 16.697031021118164),
                             ("400MATW108", 7.7f0, 52f0, 7.950294494628906), ("400MATW108", 13.2f0, 65f0, 27.994131088256836),
                             ("400MATW746", 9.0f0, 64f0, 11.517925262451172))
        @test Float64(FVSjl.r4vol_volumes(eq, d, h, 6f0, 0f0)[1]) == live
    end
end

# dense.f RELDSP species-major RELDEN (all variants) + tt/bratio.f in HTGF/DGFASP + tt/morts.f IND1 sums / PEFF powf:
# the TT stand's whole projection is exact but for the known FVS_InvReference LocationCode source skew.
@testset "TT stand exact: RELDEN RELDSP walk, tt/bratio.f in HTGF, tt/morts.f IND1/powf vs FVStt_g16" begin
    c = _case("TT", "2750433010690", "none")
    @test !c.crashed
    @test _cells(c; pred = m -> m.file != "FVS_InvReference") == 0
end

# so/regent.f:482-512 WB copies read ICR(K) of the unfilled slot (stale_icr) + so/htgf.f glibc ops.
@testset "SO WB tripled-copy ICR(K) + HTGF glibc (so/regent.f, so/htgf.f) vs FVSso_g16" begin
    c = _case("SO", "15184869010497", "none")
    @test _cells(c; pred = m -> m.file == "FVS_TreeList" && m.col in ("DBH", "DG", "HtG") && m.year <= "2040") == 0
    @test _cells(c; pred = m -> m.file == "FVS_Summary" && m.col == "QMD") == 0
end

# Calibration DGF reads cratet.f's dead-inclusive RMSQD (every variant) + ut/morts.f IND1 + UT glibc/Curtis-Arney association.
@testset "UT stand exact: calibration RMSQD, ut/morts.f IND1, UT glibc (FVSut_g16)" begin
    for cn in ("42642675010690", "2402179010690", "286785821489998")
        c = _case("UT", cn, "none")
        @test _cells(c) == 0
    end
end

# ci/dgf.f CASE(15) MC DDS −0.000981·BA + ci/morts.f SD2SQ association.
@testset "CI MC DDS BA term (ci/dgf.f:459-463) vs FVSci_g16" begin
    c = _case("CI", "12276084010690", "none")
    @test _cells(c) == 0
end

# ec/ glibc + __powisf2 htcalc + REAL AB PCTRED, and ec/ccfcal.f D*D*RD3 (RELDEN ⇒ PCTRED).
@testset "EC HTCALC/REGENT Float32 + CCFCAL D*D*RD3 (ec/htcalc.f, regent.f:89, ccfcal.f:221) vs FVSec_g16" begin
    for cn in ("450445010497", "374300286489998")
        c = _case("EC", cn, "none")
        @test _cells(c) == 0
    end
end

# nc/htgf.f AGMAX → label-140 XMOD, nc/bratio.f for every species, nc/morts.f IND1 + gates.
@testset "NC HTGF XMOD / BRATIO / MORTS IND1 (nc/htgf.f, bratio.f, morts.f) vs FVSnc_g16" begin
    c = _case("NC", "7690091010901", "none")
    @test _cells(c; pred = m -> m.year < "2031") == 0
    c = _case("NC", "15303130010497", "none")
    @test _cells(c) == 0
end

# wc|pn/regent.f label-40 small-tree HCOR calibration (was unported ⇒ CON=1) + PN/WC htcalc association.
@testset "WC/PN REGENT small-tree height calibration (regent.f label 40) vs FVSpn_g16/FVSwc_g16" begin
    c = _case("PN", "720634155290487", "none")
    @test _cells(c) == 0
    c = _case("WC", "22404557010497", "none")
    @test _cells(c) == 0
end

# CRATET RELDEN under the read-order IND1 (dead interleaved) ⇒ SO REGENT-calibration PCTRED ⇒ LP HCOR.
@testset "CRATET RELDEN read-order IND1 (dense.f/cratet.f:171) vs FVSso_g16" begin
    c = _case("SO", "449489561489998", "none")
    @test _cells(c; pred = m -> m.file != "FVS_StrClass") == 0
end

# cr/ccfcal.f CCFT·P directly (not 0.001803·CRWDTH²) ⇒ point CCF ⇒ GEMHT CCFTEM.
@testset "CR CCF from CCFT not CRWDTH² (cr/ccfcal.f) vs FVScr_clean" begin
    c = _case("CR", "3026069010690", "none")
    @test _cells(c; pred = m -> !(m.col in ("Ht2TDCF", "Ht2TDBF"))) == 0   # left: report-only NVEL Ht2TD 1-ULP (open)
end

# r6vol.f:105 short-tree Smalian VOL(1)=0.00272708*(DBHIB*DBHIB)*TTH squares first (EC/NC/SO/WC/PN Behre/R6VOL path).
@testset "R6VOL short-tree cylinder Float32 order (r6vol.f:105) vs FVSso_g16" begin
    c = _case("SO", "850566877290487", "none")
    @test _cells(c) == 0
end

# SSTAGE reads WK6=CRWDTH(I) — the TreeList CrWidth (sstage.f:238/276); SO used a stale so_cwcalc.
@testset "SO SSTAGE crown width = CRWDTH (sstage.f:238) vs FVSso_g16" begin
    for cn in ("7690240010901", "449489561489998")
        c = _case("SO", cn, "none")
        @test _cells(c; pred = m -> m.file == "FVS_StrClass") == 0
    end
end

# strp estab.f PLANT RAN∈[0,1.5] for NC + nc/esgent.f birth-cycle REGENT(LESTB).
@testset "NC PLANT RAN draw + ESGENT (estab.f:485-489, nc/regent.f LESTB) vs FVSnc_g16" begin
    for r in ("plant_cyc", "plant_cal")
        c = _case("NC", "450603388489998", r)
        @test _cells(c) == 0
    end
end
# gradd.f:96 MISTOE after MORTS+TRIPLE+REASS: spread, MISINF (MISRAN over the 27 TRIPLED records) and MISMRT on the
# tripled list (cr/mistoe.f:517-522) + cr/mismrt.f REAL rate. Was 254 cells (Mort 2005 16 vs 17, DMR mix off).
@testset "CR MISTOE post-TRIPLE seam: MISINF/MISRAN + MISMRT on tripled records (cr/mistoe.f, misinf.f) vs FVScr_clean" begin
    c = _case("CR", "5278473010690", "mistletoe")
    @test !c.crashed
    @test _cells(c) == 0
end

# f_other.f BRK_OT: Float32-literal BK, REAL DR/DBT roundings (R2 aspen 200FW2W746 Ht2TD via SF_HS/BRK_UP).
@testset "BRK_OT REAL roundings (f_other.f) - CR aspen Ht2TD vs FVScr_clean" begin
    c = _case("CR", "3026069010690", "none")
    @test _cells(c) == 0
end

# NC FFE: FMPOCR LSW white fir (nc/fmvinit.f CASE(4,9)), FMCFMD CWHR reads CRWDTH (R5CRWD, cwcalc.f:382), snag HTX/HTR1
# (nc/fmvinit.f:122,278-283) and the NC TFALL table — the fire stand's PotFire/carbon/snags all within 1e-5 of live.
@testset "NC FFE: LSW, CWHR CRWDTH, snag HTX, TFALL (nc/fmvinit.f, fmcfmd.f) vs FVSnc_g16" begin
    c = _case("NC", "23660512010900", "simfire")
    @test isempty(_material(c.ms))
end

# CR/CI/TT/UT TFALL tables (fmvinit.f, dumped from FMVINIT) for the snag-crown fall.
@testset "CR/CI/TT/UT TFALL (fmvinit.f) snag-crown fall vs live" begin
    c = _case("CI", "3369538010690", "simfire")
    @test isempty(_material(c.ms))
end

# SO: so/fmcons.f (1-3in unburned on the natural path), so/fmcba.f forest-dependent snag/decay parameters + TFALL,
# FMCBA PERCOV from _forest_crwdth, FMCFMD DSTLG from HARVYR/BURNYR.
@testset "SO FFE: fmcons, forest snag/decay params, PERCOV CRWDTH, DSTLG (so/fmcons.f, fmcba.f, fmcfmd.f) vs FVSso_g16" begin
    for (cn, r) in (("15364795010497", "salvage"), ("7690240010901", "salvage"), ("15184869010497", "simfire"),
                    ("850566877290487", "simfire"))
        c = _case("SO", cn, r)
        @test isempty(_material(c.ms))
    end
end

# CR SSTAGE/FMCBA read CRWDTH clamped to [0.5,99.9] (cwcalc.f:2391-2392): seedling stratum species order.
@testset "CR CRWDTH clamp in SSTAGE (sstage.f:276, cwcalc.f:2391) vs FVScr_clean" begin
    c = _case("CR", "46279527020004", "none")
    @test _cells(c) == 0
end
# ut|tt|nc|on/morts.f `IF(ICYC.GT.1 .AND. ABS(T-TPAMRT).GT.1.)` resets CEPMRT/SLPMRT at ICYC=2 already; jl's 0-based
# control.cycle made the test `cycle > 1` (ICYC>2) ⇒ the cycle-2 reset after a cycle-1 trajectory change (here the MISTOE
# kill: live T 4475.32 vs TPAMRT 4568.13) never fired and jl kept the latched line (cycle-2 mortality ~20% low for every
# species). MEASURED FVSut_g16 42642675010690 MISTOE DEBUG MORTS.
@testset "MORTS ICYC>1 TPAMRT reset at cycle 2 (ut/morts.f:234) vs FVSut_g16" begin
    c = _case("UT", "42642675010690", "mistletoe")
    @test !c.crashed
    @test _cells(c) == 0
end
# gradd.f:96 MISTOE precedes gradd.f:118 FMMAIN; on a non-fire tripling cycle jl runs MISTOE post-TRIPLE (mis_post), so the
# deferred R6 FFE annual loop (fmmain.f:228 FMSNAG → FMR6HTLS RANN per snag pool) must wait for the spread's draws too.
# jl drew Y 0.920/0.025 at the pre-TRIPLE seam; live 0.61/0.76 (= 275 spread draws later) ⇒ the 2015 LP snag lost 78% of
# its height in jl, 0 live (MEASURED FVSso_g16 374286168489998 SALVAGE, DEBUG FMR6HTLS FMSNAG).
@testset "R6 FFE annual loop after MISTOE spread on tripling cycles (gradd.f:96/:118, fmsnag.f, fmr6htls.f) vs FVSso_g16" begin
    c = _case("SO", "374286168489998", "salvage")
    @test !c.crashed
    @test _cells(c; pred = m -> m.col in ("Standing_Dead", "Total_Stand_Carbon") || (m.file == "FVS_Carbon" && _rel(m) > 1e-5)) == 0
end
# CI PLANT cohort (ci/esinit.f STOADJ=0 ⇒ esnutr.f catch-all ESTAB, estab.f no-stocking plot chain): per-plot EMSQR/FIRST(2)
# ESSUBH heights, ESGENT HTG·WK4 (HTIMLT 5/5.0001), and the next-cycle MORTS WK1 = the ESGENT DG (ci/dgdriv.f:171).
# MEASURED FVSci_g16 3159852010690 PLANT DF 400 (DEBUG ESTAB/ESGENT/MORTS).
@testset "CI PLANT cohort: estb catch-all plot chain, ESGENT WK4, WK1 (ci/esinit.f, estab.f, esgent.f, morts.f) vs FVSci_g16" begin
    for r in ("plant_cal", "plant_cyc")
        c = _case("CI", "3159852010690", r)
        @test !c.crashed
        @test _cells(c) == 0
    end
end
# tt/dbsreference.f (LOCCODE upgrade, only in the TT build) adds `LocationCode int` = KODFOR to FVS_InvReference.
@testset "TT FVS_InvReference LocationCode = KODFOR (tt/dbsreference.f:36,98-121) vs FVStt_g16" begin
    c = _case("TT", "2780339010690", "none")
    @test _cells(c; pred = m -> m.file == "FVS_InvReference") == 0
end
# EC PLANT cohort: ec/esgent.f HTG·WK4 (HTIMLT 5/5.0001) and the gradd.f:192 DENSE PCCF for the ec/regent.f:285 LESTB crown
# (an INT(CR·100+0.5) tie flipped on one 374300286489998 record when jl read the post-ESTAB point CCF).
@testset "EC PLANT cohort: ESGENT WK4 + pre-ESTAB PCCF crown (ec/esgent.f, ec/regent.f:285) vs FVSec_g16" begin
    for cn in ("450507010497", "374300286489998")
        c = _case("EC", cn, "plant_cal")
        @test !c.crashed
        @test _cells(c) == 0
    end
end
# fvsvol.f:340 HT2TD(IT,2)=MAX(HT1PRD,HT2PRD) with nsvb.f:417 HT2PRD lifted to the bucked-log top HTsaw (NVB equations).
@testset "NVB Ht2TDCF = MAX(HT1PRD, HTsaw) (fvsvol.f:340, nsvb.f:368-417) vs FVScr_clean" begin
    c = _case("CR", "5278473010690", "none")
    @test _cells(c) == 0
end
end # module
