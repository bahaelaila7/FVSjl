# test_ocopon_tiered.jl — OC (Oregon Coast) / OP (Olympic) / ON (Ontario) fixes vs the LIVE tiered goldens
# (test/fixtures/tiered/{oc,op,on}: stands.db + <stand>_<regime>.key + the live .sum / DBS tables of the private
# relinks /workspace/.{oc,op,on}work/tiered/FVS<v>_tiered, ORACLE_SOURCE_AUDIT §8). Each testset names the Fortran it
# follows and the case it was measured on; a testset asserts only what its mechanism moved.
module OcOpOnTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

const _CACHE = Dict{Tuple{String,String,String},Any}()
function _case(v, cn, r)
    get!(_CACHE, (v, cn, r)) do
        d = mktempdir()
        txt, db, crashed, _ = run_case(v, cn, r; dir = d)
        (txt = txt, db = db, crashed = crashed, ms = compare_case(v, cn, r, txt, db))
    end
end
_yr(m) = tryparse(Int, m.year)
_n(ms; file = nothing, col = nothing, year = nothing, before = nothing) =
    count(m -> (file === nothing || m.file == file) && (col === nothing || m.col == col) &&
               (year === nothing || m.year == string(year)) &&
               (before === nothing || ((k = _yr(m)); k !== nothing && k < before)), ms)

# {oc,op}/grinit.f DGSD=0 — FVS has ONE DGSD; jl mirrors it in dg_sd AND dg_stddev_bound (the latter gates the
# dgdriv.f calibration OLDRN seeding / DO 191 BACHLO draws and the ESTAB ±10% random). OC/OP set only dg_sd.
@testset "OC/OP DGSD=0 in both DGSD fields ({oc,op}/grinit.f)" begin
    for v in (FVSjl.OregonCoast(), FVSjl.Olympic())
        s = FVSjl.StandState(v)
        FVSjl.init_blockdata!(s, v)
        @test s.control.dg_sd == 0f0
        @test s.control.dg_stddev_bound == 0f0
    end
end

# OC 645182389126144 (50 measured-DG records, 7 FIA DF dwarf-mistletoe ratings): before this round the stand read 7,688
# tiered cells (cycle-1 Accr 193 vs 307: ORGANON DG half the live value).
@testset "OC LSTART DGDRIV calibration: DO 105 inside-bark DG + COR (oc/dgdriv.f:597-1080, oc/cratet.f:875)" begin
    c = _case("OC", "645182389126144", "none")
    @test !c.crashed
    @test _n(c.ms; file = "FVS_TreeList", col = "DG", year = 2018) == 0     # was 50 (outside-bark measured DG)
end

@testset "ORGANON DGCALIB (organon/diamcal.f) + EXECUTE CALIB decay (execute2.f:231-252) + STOR carry / RVARS MSDI" begin
    c = _case("OC", "645182389126144", "none")
    @test _n(c.ms; file = "sum") == 0                                        # every .sum row, 2018-2043
    @test _n(c.ms; file = "FVS_TreeList", col = "DG", year = 2023) == 0      # DGCALIB CALIB(3)=2.0 ⇒ DG ×2
    @test _n(c.ms; file = "FVS_TreeList", col = "MortPA", year = 2023) == 0  # MORTAL_RUN A1 from RVARS(3) (pre-SDICHK)
end

@testset "OC dwarf mistletoe (misintoc.f; oc/dgdriv.f:562 MISDGF, oc/htgf.f:266 MISHGF, mistoe.f:522 MISMRT)" begin
    c = _case("OC", "645182389126144", "none")
    @test _n(c.ms; file = "FVS_TreeList", col = "MistCD") == 0               # FIA 23023 ratings were dropped
    @test _n(c.ms; file = "FVS_TreeList", col = "TPA", before = 2043) == 0
end

@testset "OC REGENT small trees: HTDBH diameter step, ORGANON HTG floor, DBH+=DG/BRATIO (oc/regent.f:257-370, update.f:115)" begin
    c = _case("OC", "645182389126144", "none")
    @test _n(c.ms; file = "FVS_TreeList", col = "DBH", before = 2038) == 0   # 0.1" seedlings: DBH 0.225 live
end

@testset "NVEL volinit.f:168 DBH<1 ⇒ no volume (OC/OP BLM + R6)" begin
    c = _case("OC", "645182389126144", "none")
    @test _n(c.ms; file = "FVS_TreeList", col = "TCuFt", year = 2018) == 0
end

# OC 15148784010497 (no big-6 ⇒ ORGANON never runs, SMORMT=0): jl killed nothing (.sum Mort 0 vs 2 from 2010).
@testset "OC FVS-native MORTS when SMORMT=0 (oc/morts.f:463-546 == ca/morts.f with FINT/5, oc/varmrt.f == ca/varmrt.f)" begin
    c = _case("OC", "15148784010497", "none")
    @test !c.crashed
    @test isempty(c.ms)
end

# OP 374608587489998 (no big-6 ⇒ every record FVS-native; 48 measured-DG records, 30 of them WH): op/dgf.f skipped the
# ORGANON-valid species (WK2 stale ⇒ WH calibration CSNXX=0 ⇒ COR NaN ⇒ NaN heights), had no red-alder equation, and the
# native kernels ran on Julia exp/log/^ instead of glibc expf/logf/powf / libgcc __powisf2.
@testset "OP op/dgf.f every record + RA equation + glibc kernels (dgf.f, htcalc.f, htgf.f, morts.f, crown.f, smhgdg.f)" begin
    c = _case("OP", "374608587489998", "none")
    @test !c.crashed
    @test _n(c.ms; file = "sum", before = 2025) == 0
    @test _n(c.ms; file = "FVS_TreeList", before = 2020) == 0
end

# ON 245908831010661 (FIA LS stand converted to metric; 16 measured DG, 3 inventory-dead, FIA 356 serviceberry):
# 17,636 cells before this round.
@testset "ON SPCTRN default NC (canada/on/spctrn.f CASE DEFAULT = 49)" begin
    c = _case("ON", "245908831010661", "none")
    @test !c.crashed
    @test _n(c.ms; col = "SpeciesFVS") == 0 && _n(c.ms; file = "FVS_TreeList_East_Metric", col = "ROWCOUNT") == 0
end

@testset "ON CRATET crown dub on the backdated dead-inclusive DENSE (canada/on/cratet.f:159-161, :501)" begin
    c = _case("ON", "245908831010661", "none")
    @test _n(c.ms; file = "FVS_TreeList_East_Metric", col = "PctCr", year = 2009) == 0     # seedling ICR 37/27 vs 39/30
    @test _n(c.ms; file = "FVS_TreeList_East_Metric", col = "BAPctile", year = 2009) == 0  # dead PCT on backdated WK5
end

@testset "ON REGENT small-tree HCOR (canada/on/regent.f:128 CON, :419-563 calibration)" begin
    c = _case("ON", "245908831010661", "none")
    @test _n(c.ms; file = "FVS_TreeList_East_Metric", col = "HtG", year = 2019) == 0       # MH seedlings +4.2% before
    @test _n(c.ms; file = "sum", before = 2029) == 0
end

@testset "metric FVS_StrClass_Metric table name (metric/dbsqlite/dbsstrclass.f:126-128)" begin
    c = _case("ON", "245908831010661", "none")
    @test _n(c.ms; col = "PRESENCE") == 0
end

# ON PLANT (metric/strp/esin.f, metric/strp/estab.f, strp/essubh.f, strp/esgent.f → canada/on/regent.f LESTB): every ON
# plant_cal/plant_cyc case crashed (KeyError :ht_curve_b1 — the SN ESSUBH site curve).
@testset "ON establishment: metric PLANT, strp ESTAB/ESSUBH, ON REGENT(LESTB)" begin
    for (cn, r) in (("92320368010661", "plant_cyc"), ("21769405010661", "plant_cal"), ("720549798290487", "plant_cal"))
        c = _case("ON", cn, r)
        @test !c.crashed
        @test isempty(c.ms)
    end
end

# op/sichg.f:36-37 `DATA REFLOC/10*'B','T',10*'B','T',17*'B'/` — the second 'T' is species 22 (RA), not 21.
@testset "OP SICHG REFLOC (op/sichg.f:36-52): 'T' at LP(11) and RA(22)" begin
    @test findall(==(1), FVSjl.OP_SICHG_REFLOC) == [11, 22]
    sa = FVSjl.op_sichg(22, 60f0)                      # RA site species ('T'): every 'B' species gets DIFF=-1
    @test sa[22] == Float32(FVSjl.OP_SICHG_REFAGE[22]) && sa[11] == Float32(FVSjl.OP_SICHG_REFAGE[11])
    @test sa[21] < Float32(FVSjl.OP_SICHG_REFAGE[21])  # was DIFF=0 (21 mis-tagged 'T')
end

# OP 30194263010497 (2 red alder among ORGANON big-6): RA took the HD_NWO height ratio and never reached RAMORT.
@testset "ORGANON red alder: HTGRO2 RAH40 increment (htgrowth.f:241-255), RAAGE/RASI (execute2.f:297-320), RAMORT (mortality.f:122-127, :355-410)" begin
    r = FVSjl.organon_ramort!(Float32[-3f0, -3f0, 0f0], Int32[351, 351, 202], Float32[10f0, 10f0, 5f0], 3, 55f0, 20f0)
    @test r[1] == r[2] && r[1] > -3f0 && r[3] == 0f0    # alders lifted by one common KR1; DF untouched
    c = _case("OP", "30194263010497", "none")
    @test !c.crashed
    @test _n(c.ms; file = "sum") == 0
    @test _n(c.ms) == _n(c.ms; file = "FVS_TreeList", col = "HtG")   # only the GF HTGRO1 1-ULP HtG corner remains
end

# FVS14 "habitat type/plant association not recognized" (op/habtyp.f:131 == pn/habtyp.f; oc/habtyp.f == ca/habtyp.f):
# OC/OP had no habtyp_errors! method ⇒ every FIA stand lacked the FVS_Error row.
@testset "OC/OP habtyp FVS14 errors (op/habtyp.f:131, oc/habtyp.f)" begin
    @test isempty(_case("OC", "24382439010900", "none").ms)
    @test _n(_case("OP", "30194263010497", "none").ms; file = "FVS_Error") == 0
end

# OP FIA PV_CODE (alpha plant-association code) -> KODTYP (op/habtyp.f HBDECD == pn/habtyp.f, op/pvref6.f == pn/pvref6.f
# byte-identical). The FIA reader decoded PV_CODE only for PN/WC, so OP kept habitat 0 => CHS133 default fan (DF SI 98)
# where live resolves CHS137 (ECOCLS DF 108): FVS_InvReference SiteIndex, TreeList growth, every .sum column.
@testset "OP PV_CODE habitat decode (op/habtyp.f HBDECD)" begin
    c = _case("OP", "12973618010497", "none")
    @test !c.crashed
    @test isempty(c.ms)
end

# OC 645142535126144 thinbba (THINBBA cycle 2 leaves one DF; ESTAB adds 8 ES records afterwards). Two ORGANON mechanisms:
# (1) oc/dgdriv.f:304-327 "MANAGE THE THINNING VARIABLES" — OLDBA(pre-CUTS, grincr.f:283) > BA ⇒ YST(1)=(CYCLG-1)*5,
#     BART(1)=OLDBA-ATBA, BABT=OLDBA, fed to DG_THIN/HG_THIN (diagro.f:578, htgrowth.f:957): the surviving DF grew 11% slow;
# (2) oc/crown.f:281-288 stale slot IORG=1 with CR2=0 for the ES records ⇒ ICR=1 ⇒ ~95% suppression mortality.
# Before: 228 mismatches (.sum TPA 9 vs 132 from 2033); after: only StrClass cover rounding + 1 QMD ULP-class cell.
@testset "OC ORGANON thinning response (dgdriv.f:304-327, diagro.f DG_THIN, htgrowth.f HG_THIN) + stale-IORG CROWN" begin
    c = _case("OC", "645142535126144", "thinbba")
    @test !c.crashed
    @test _n(c.ms; file = "sum") == 0
    @test _n(c.ms; file = "FVS_Summary", before = 2043) == 0
    @test _n(c.ms) <= 4
end

# OP: same two mechanisms through the NWO engine (op/dgdriv.f:304-327 identical block, op/diagro.f DG_THIN == oc, op/htgrowth.f
# HG_THIN == oc; op/crown.f:279-286 stale-IORG + CR2=0). 30194263010497 thinbba: TPA 35 live vs 86 jl at 2023 before.
@testset "OP ORGANON thinning response + stale-IORG CROWN (op/dgdriv.f:304-327, op/crown.f:279)" begin
    for cn in ("24164991010900", "30194263010497")
        c = _case("OP", cn, "thinbba")
        @test !c.crashed
        @test _n(c.ms; file = "sum") == 0
        @test _n(c.ms; file = "FVS_Summary") == 0
    end
end

# REGENT(LESTB) DBH step for a PLANTed record (FINT≤5 ⇒ LSKIPH, HTG=0): oc/regent.f:285-345 (HK≤4.5 ⇒ DBH=D+0.001·HK; else
# HTDBH⁻¹, DIAM floor, +0.001·HK, DG=DBH) and op/regent.f:325-406 (HK<4.5 ⇒ D+0.001·HK; else DBH=DG=0 ⇒ DIAM(sp)), then DGBND and
# ESGENT (esgent.f:20-34, WK4=0.99998 ⇒ DG=DBH). jl left every new record at the strp DBH=0.1 (OC .sum QMD/BA/CCF/SDI off from the
# birth cycle: OC 24382439010900 plant_cyc 62 cells, OP 22964003010497 plant_cyc 88 cells -> 0).
@testset "OC/OP REGENT(LESTB) establishment DBH (oc/regent.f:285-345, op/regent.f:325-406, esgent.f)" begin
    for (v, cn) in (("OC", "24382439010900"), ("OP", "22964003010497"))
        for r in ("plant_cyc", "plant_cal")
            c = _case(v, cn, r)
            @test !c.crashed
            @test isempty(c.ms)
        end
    end
end

# OC FFE runs on the PRE-growth stand: organon_apply_growth! grows DBH/HT/CR inline, but FVS's FMMAIN (FMBURN/FMEFF/FMCBA) runs in
# GRADD before UPDATE (live DEBUG FMEFF 30193600010497 2013: DBH 7.87 / HT 56.8 vs jl 8.65 / 65.5 ⇒ PMORT 0.588 vs 0.329), and its
# crown cover reads the stored CRWDTH (PERCOV 14.603 live). Before: .sum Mort 2013 44 vs 76, 239 mismatches; after: only FVS_Carbon.
@testset "OC FFE sees the cycle-start stand + stored CRWDTH (simfire)" begin
    c = _case("OC", "30193600010497", "simfire")
    @test !c.crashed
    for f in ("sum", "FVS_Summary", "FVS_Mortality", "FVS_BurnReport")   # FVS_PotFire/FVS_Carbon post-fire residuals stay OPEN
        @test _n(c.ms; file = f) == 0
    end
end

# OP FFE crown biomass: fmcrow.f:159-164 SELECT CASE (SPIW) sends only the hardwoods 26,27,34,35,36,37,39 to FMCROWE (Jenkins); every
# other OP species uses FMCROWW (ISPMAP). jl routed ALL OP species through FMCROWE (conifer CROWNW 0.46x live: DF 7.3in/55ft CROWNW(0)
# 14.13 vs 30.033 live, FMDOUT dump FVSop_tiered 24164991010900 2002) ⇒ Aboveground_Total_Live 7.81 vs 8.16 and the whole fire year.
@testset "OP FFE crown biomass = FMCROWW for conifers (op/fmcrow.f:159-164)" begin
    c = _case("OP", "24164991010900", "simfire")
    @test !c.crashed
    @test _n(c.ms; file = "FVS_Carbon", col = "Aboveground_Total_Live", year = 2002) == 0
    @test _n(c.ms; file = "FVS_BurnReport", col = "Flame_length") == 0
    @test FVSjl.op_uses_fmcrowe(26) && FVSjl.op_uses_fmcrowe(39) && !FVSjl.op_uses_fmcrowe(16) && !FVSjl.op_uses_fmcrowe(11)
end

# Stand Carbon Report merch pool (fmcrbout.f:124-130 LMERCH=LVWEST=.TRUE.; fire/base/fmsvol.f:146-153): OC/OP are not CS/LS/NE/SN, so
# VOL2HT=MCF with NO cone floor MAX(X,MCF) — the 0.0055·H floor added 0.2 cuft to every unmerchantable tree (OP 24164991010900 2002
# Aboveground_Merch_Live 1.984 vs 1.4758 live; per-record VT dump FVSop_tiered: tree 8 D6.2 VT 0.0, jl 0.218).
@testset "OC/OP Aboveground_Merch_Live has no cone floor (fmcrbout.f:127, fmsvol.f:146-153)" begin
    c = _case("OP", "24164991010900", "simfire")
    @test _n(c.ms; file = "FVS_Carbon", col = "Aboveground_Merch_Live", year = 2002) == 0
end

# op/fmcba.f carries no COVCA/COVCAWT: FLIVE (:453-457) and the initial dead fuel STFUEL (:503-504) load from the single COVTYP. jl used
# the NC/WS/CA/OC BA-weighted top-2 species pair (OP 24164991010900: 0.899 DF + 0.101 LP) ⇒ Forest_Shrub_Herb 0.19495 vs 0.2000 and
# Forest_Down_Dead_Wood 21.46 vs 22.806 at the inventory year, and the fire's fuel models/flame/scorch.
@testset "OP FMCBA loads live+dead fuel from the single COVTYP (op/fmcba.f:453-457,503-504)" begin
    c = _case("OP", "24164991010900", "simfire")
    @test _n(c.ms; file = "FVS_Carbon", year = 2002) == 0
    @test _n(c.ms; file = "FVS_BurnReport") == 0
    @test _n(c.ms; file = "FVS_Mortality") == 0
    @test _n(c.ms; file = "FVS_PotFire", year = 2002) == 0
end

# FMSVOL/FMSVL2 for OC/OP (fire/base/fmsvol.f:98-142): NATCRS on the stand's volume equation — a Region-6 forest (KODFOR 6xx) runs
# R6_EQN (westside FW2 / INGY FW2 / 616BEHW), a BLM forest BLMVOL — at (D,H) with BARK=BRATIO(JS,D,H), then CFTOPK at INT(XHT*100) on
# every snag report (CTKFLG=T). jl valued FFE snags (and the live-stem carbon) with the BLM total cubic on every forest: OP 24164991010900
# (forest 609 R6) snag D5.2/H28 LP TCF 2.9802 live vs 1.82; DF D7.09/H48.06 5.1 vs 7.29 ⇒ Standing_Dead 12.705 vs 12.473 after the fire.
@testset "OC/OP FFE snag + live-stem volume = FMSVOL on the stand's NATCRS equation (fmsvol.f:98-142)" begin
    c = _case("OP", "24164991010900", "simfire")
    @test _n(c.ms; file = "FVS_Carbon", col = "Standing_Dead", year = 2007) == 0
    @test _n(c.ms; file = "FVS_Carbon", year = 2007) == 0
    @test _n(c.ms; file = "FVS_BurnReport") == 0
    @test _n(c.ms; file = "FVS_Mortality") == 0
end

# OC FFE FMOLDC: OLDHT/OLDCRL are the PRE-growth HT/ICR (fmoldc.f closes FMMAIN, which runs in GRADD before UPDATE). OC's ORGANON hook grows
# the stand inline, so the non-fire cycle's ffe_fmoldc! stored GROWN heights ⇒ next FMSDIT's crown-lift X was 0 for those records
# (OC 30193600010497 2013 FD dump: OLDCRW(1) 0.10202 live vs 0 jl) ⇒ the fire snags' CWD2B2 sizes 1-4 (YRSCYC*OLDCRW term, fmscro.f:147) low.
# oc/fmvinit.f:481-510: snag HTX(I,1:4)=1/0/1/0 + HTR1 0.03406 (fire snags lose height: 48.9 -> 41.1 ft by 2018) and TFALL(I,1)=10, (2)=15,
# (4:5)=TFALL(3) (the generic tfall class rows dropped a fire's crowns in ~4 yr instead of 10-20).
@testset "OC FFE: pre-growth FMOLDC, snag height loss HTX/HTR1, crown TFALL (oc/fmvinit.f, fmoldc.f, fmsdit.f)" begin
    c = _case("OC", "30193600010497", "simfire")
    @test _n(c.ms; file = "FVS_Carbon") == 0
    @test isempty(c.ms)                                  # incl. FVS_PotFire 2028 Fuel_Mod1 5 (CWHR on the stored CRWDTH, oc/fmcfmd.f)
    @test FVSjl._fm_tfall_table(FVSjl.OregonCoast())[18, 1:6] == Float32[3, 10, 10, 10, 10, 10]
end

# OP snag decay: op/fmvinit.f sets DECAYX(I)=1.0 for every species (MEASURED FVSop_tiered FMDOUT dump, all 39); the OP fire_species_props.csv
# carried the OC value 999 ("hard snags never soften") ⇒ FMSNGDK DKTIME = JYRSOFT·DECAYX = 21·999, so TSOFT never bounded a fire snag's
# crown fall (FMSCRO: RLIFE = MIN(TFALL, TSOFT); live LP D5.1 TSOFT 11 ⇒ ILIFE 11, jl 15) ⇒ CWD2B size-3 class 2 29.18 live vs 33.64 lb.
@testset "OP DECAYX = 1.0 (op/fmvinit.f): fire-snag crown fall bounded by TSOFT" begin
    c = _case("OP", "24164991010900", "simfire")
    @test _n(c.ms; file = "FVS_Carbon", col = "Standing_Dead", before = 2022) == 0
    @test _n(c.ms; file = "FVS_Carbon", col = "Forest_Down_Dead_Wood") == 0
    @test _n(c.ms; file = "FVS_PotFire", before = 2017) == 0
end

# OC bare-stand cover type (oc/fmcba.f:529-545): no live BA in the first FFE year ⇒ COVTYP = COVINI6(ITYPE) on an R6/BLM forest
# (IFOR ≥ 6) else COVINI5(ITYPE) (tables byte-identical to ca/fmcba.f), default 7; jl always used 7 ("COVINI waits on habtyp").
# OC 24446507010900 (no trees): inventory Forest_Shrub_Herb 1.2 vs live 1.15, Forest_Down_Dead_Wood 1.65 vs 3.05 (initial FUINI by cover).
@testset "OC bare-stand COVTYP = COVINI(ITYPE) (oc/fmcba.f:529-545)" begin
    c = _case("OC", "24446507010900", "simfire")
    @test !c.crashed
    @test isempty(c.ms)
end

# OC first-year FFE decay rates (oc/fmcba.f:575-737, == ca/fmcba.f): the Oregon base DKRT × DKRADJ(CAHMC(ITYPE),CAWMD(ITYPE),K) with the
# stand's habtyp ITYPE (CDC321 ⇒ ITYPE 3, TEMP 2, MOIST 3 — MEASURED FVSoc_tiered); jl hard-coded the mesic CWC221 row for every stand.
# OC 15148784010497: 2015 Forest_Down_Dead_Wood 0.454 live vs 0.436; salvage/simfire Carbon+PotFire 39 cells each -> 0.
@testset "OC FFE decay rates by habitat ITYPE (oc/fmcba.f:575-737)" begin
    for r in ("simfire", "salvage")
        c = _case("OC", "15148784010497", r)
        @test !c.crashed
        @test isempty(c.ms)
    end
end

# OP CROWN: the `ICRI < 10` bumps (op/crown.f:385 CRMAX, :429 floor 10, stmt 58) are guarded `(.NOT. LORGANON) .AND. (IORG(I) .NE. 1)`;
# LORGANON=.TRUE. (op/grinit.f:340) ⇒ they never fire, for ORGANON and non-ORGANON trees alike (OC already omitted them). jl bumped every
# non-ORGANON tree: OP 24164991010900 fire-scorched LP ICR 1 live / 10 jl at 2017 ⇒ the survivors' DBH, QMD 14.718 vs 14.695, PotFire.
@testset "OP CROWN has no ICRI<10 bump under LORGANON=.TRUE. (op/crown.f:385,429)" begin
    c = _case("OP", "24164991010900", "simfire")
    @test !c.crashed
    @test isempty(c.ms)
end

# ORGANON integer-power precedence: Fortran `B2*DBH**2` is B2*(DBH*DBH) (the power binds first), jl wrote (b2*dbh)*dbh — a 1-ULP
# difference. MCW_NWO/MCW_SWO (crngrow.f:508,567) -> LCW -> CCH crown-closure table -> TCCH (HG_NWO) -> HTGRO1 height growth (OP
# 30194263010497 2013 HtG 7.8015337 live vs 7.8015332, K=29: TCCH 21.724193573 live vs 21.724195; CCH(9,11,12,13,14,18,32,36,40) each 1 ULP off,
# first at one DF D12.8/H71 record's LCW 23.464387894 vs 23.46438980). Also PM_NWO (mortality.f:533-733 `B2*DBH**2`) and the SWO crown
# HCB `B6*OG**2` (start2.f:534-537). OC 645142535126144 none HtG 2038: 3 cells -> 0.
@testset "ORGANON B*X**2 evaluates as B*(X*X) (crngrow.f MCW, mortality.f PM, start2.f HCB)" begin
    for (v, cn, r) in (("OP", "30194263010497", "none"), ("OP", "30194263010497", "thinbba"), ("OC", "645142535126144", "none"))
        c = _case(v, cn, r)
        @test !c.crashed
        @test _n(c.ms; file = "FVS_TreeList", col = "HtG") == 0
        @test _n(c.ms; file = "FVS_CutList") == 0
    end
end

end # module
