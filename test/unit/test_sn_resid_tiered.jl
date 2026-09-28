# test_sn_resid_tiered.jl — SN per-record fixes vs the LIVE FVSsn_g16 tiered goldens (test/fixtures/tiered/sn: stands.db +
# <stand>_<regime>.key + the live .sum / DBS tables). Stand 157577477010854 (LP plantation, NONE + TREELIST, 1972-1997) —
# the SN suite's largest residual block (11,599 cells) — is walked year by year; each testset names the Fortran it follows.
module SNResidTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

const _C = Ref{Any}(nothing)
function _case157()
    _C[] === nothing || return _C[]
    d = mktempdir()
    txt, db, crashed, _ = run_case("SN", "157577477010854", "none"; dir = d)
    _C[] = (crashed = crashed, ms = compare_case("SN", "157577477010854", "none", txt, db))
end
_cells(file, cols, yrs) = count(m -> m.file == file && m.col in cols && m.year in yrs, _case157().ms)

const _C238 = Ref{Any}(nothing)
function _case238()
    _C238[] === nothing || return _C238[]
    d = mktempdir()
    txt, db, crashed, _ = run_case("SN", "238813815010854", "thinbba"; dir = d)
    _C238[] = (crashed = crashed, ms = compare_case("SN", "238813815010854", "thinbba", txt, db))
end
_cells238(file, cols, yrs) = count(m -> m.file == file && m.col in cols && m.year in yrs, _case238().ms)

const _CX = Dict{Tuple{String,String},Any}()
function _casex(st, rg)
    get!(_CX, (st, rg)) do
        d = mktempdir()
        txt, db, crashed, _ = run_case("SN", st, rg; dir = d)
        (crashed = crashed, ms = compare_case("SN", st, rg, txt, db))
    end
end
_cellsx(st, rg, file, cols, yrs) = count(m -> m.file == file && m.col in cols && m.year in yrs, _casex(st, rg).ms)

# sn/crown.f:156-159 ISORT(IND(JJ)) = ITRN−JJ+1 over FVS's IND — at the inventory CRATET's IND1-seeded RDPSRT(.FALSE.)
# (sn/cratet.f:155-157; :261 RDPSRT(.TRUE.) with dead records) — so tied diameters rank in RDPSRT's order, and the
# Weibull X = ISORT/ITRN·SCALE dubs their crowns accordingly. jl ranked by a stable ascending sortperm (MEASURED 1972:
# tied LP pairs 7.1"/49.0', 6.4"/45.4', 5.9"/42.6' had dubbed CRs 34/35, 31/29, 26/27 swapped against live).
@testset "SN crown ISORT from FVS's RDPSRT IND (sn/crown.f:156-159) vs FVSsn_g16" begin
    @test !_case157().crashed
    @test _cells("FVS_TreeList", ("PctCr", "CrWidth"), ("1972",)) == 0
end

# sn/cratet.f:300-305 fits AA as SUMX = SUMX+YY−XX ((SUMX+ALOG(H−4.5))−BX/(D+1)) and dubs H = EXP(AX+BX/(D+1))+4.5, and
# htdbh.f:292-297 is H = 4.5+P2·EXP(−1.·P3·D**P4) — REAL*4 logf/expf/powf. jl used Julia log/exp/^ and SUMX+(YY−XX)
# (MEASURED 1972: the dubbed LP 6.9" height 47.98899460 live vs 47.98899078).
@testset "SN height dub in REAL*4 (sn/cratet.f:300-360, htdbh.f:292-297) vs FVSsn_g16" begin
    @test _cells("FVS_TreeList", ("Ht", "EstHt"), ("1972",)) == 0
end

# htcalc.f:170 AGET = 1./B3*(ALOG(1-((H-HB)/B1/SI**B2)**(1./B4/SI**B5))) divides in turn; jl divided by the products
# (MEASURED FVSsn_g16 private HTGF trace, 1972: AGET 1-2 ULP off on 10 of 35 records ⇒ HTG1 ⇒ 1977 HtG).
@testset "SN HTCALC tree age at htcalc.f's association (htcalc.f:170) vs FVSsn_g16" begin
    @test _cells("FVS_TreeList", ("HtG",), ("1977",)) == 0
end

# htdbh.f:296-310 MODE 1 (REGENT's DK/DKK, sn/regent.f:317-321): HAT3 = 4.5+P2·EXP(−1.·P3·3.0**P4);
# D = (ALOG(MIN(H−4.5,0.9999·P2))−ALOG(P2))/(−1.·P3); D = EXP(ALOG(D)·1./P4) — logf/expf/powf, (ALOG(D)·1.)/P4 not ·(1/P4),
# and MIN on H−4.5 (MEASURED 1977: a 2.7" LP's REGENT DG 0.3353589 live vs 0.3353593).
@testset "SN HTDBH height→diameter inverse in REAL*4 (htdbh.f:296-310) vs FVSsn_g16" begin
    @test _cells("FVS_TreeList", ("DG", "DBH"), ("1977",)) == 0
end

# dense.f (one file in 23 of the 24 builds, not ON) sums TPROB/TSUMD2 over IND1 with WK5 = D·(D·P); RMSQD = QMD. jl did
# that only for BM/EM/IE/AK (MEASURED 1977: QMD 6.26626635 live vs 6.26626682 record order) — the QMD feeds SN's
# SDI/mortality, so every later cycle drifted.
@testset "SN QMD over IND1 with dense.f's WK5 (dense.f:179-188) vs FVSsn_g16" begin
    @test _cells("FVS_Summary", ("QMD", "ATQMD"), ("1977", "1982", "1987", "1992", "1997")) == 0
end

# ptbal.f:148 (one file in every build) XBALT = XBALT + WK5·.005454154·PI/GROSPC, WK5 = D·(D·P): the point BA/BAL that
# DGF's PBAL term reads. jl used P·(0.005454154·D²)·PI/GROSPC outside AK (MEASURED FVSsn_g16 private DGF trace, 1982:
# PBAL 1 ULP off on 263 of 315 records ⇒ DDS). With it the whole stand — .sum, Summary, TreeList, StrClass — equals live.
@testset "SN point BAL from ptbal.f's WK5 form (ptbal.f:148): 157577477010854 NONE all cells vs FVSsn_g16" begin
    @test isempty(_case157().ms)
end

# esnutr.f:119-125 (every build): after ESUCKR, whenever ITRNRM>=1, IREC1=ITRN; CALL SPESRT — IND1 is relisted in
# ascending physical record order, dropping the post-TRIPLE REASS lineage; the sprout records also start with
# WK1=WK2=WK4=0 (esuckr.f:321-328). 238813815010854 THINBBA (2005 thin of the sprouting hardwoods): jl kept the lineage
# key, so cycle-3 DGSCOR walked the WN records 3,1,4,5,2,6 (live 1..6) and every later BACHLO draw landed on another
# tree (MEASURED private DGSCOR trace: cyc-3 RANN draws 105 live vs 108); the reused slots also printed the deleted
# records' MortPA in 2010 (private TREELIST run).
@testset "SN SPESRT after ESUCKR (esnutr.f:119-125, esuckr.f:321-328): 238813815010854 THINBBA vs FVSsn_g16" begin
    @test !_case238().crashed
    @test _cells238("FVS_Summary", ("Tpa", "BA", "QMD", "TCuFt", "Acc", "Mort"), ("2010", "2015", "2020", "2025")) == 0
    @test count(m -> m.file == "sum", _case238().ms) == 0
end

# essprt.f:514-594 CASE('SN') writes PREM = PREM * 1. / (1. + EXP(-(A + B*DSTMP))) — left to right, (PREM*1.)/(1+expf),
# not PREM*(1/(1+e)) (MEASURED 238813815010854 THINBBA: the AE/SU sprouts' TPA 166.51021 vs live 166.51022 ⇒ StrClass
# Stratum_2_DBH 1 ULP). With it the whole THINBBA case equals live.
@testset "SN sprout survival at essprt.f's association (essprt.f:514-594): 238813815010854 THINBBA all cells vs FVSsn_g16" begin
    @test isempty(_case238().ms)
end

# fmsvol.f:98-153 (CS/LS/NE/SN): every report recomputes each snag's bole as NATCRS(DBHS, HTDEAD) cut at the CURRENT height
# (XHT=HTIH ⇒ LTKIL, IHT=INT(XHT·100), CFTOPK), VOL2HT = MAX(0.005454154·HTDEAD, MCF). jl truncated only under SNAGBRK, so
# an inventory snag whose current height sits below its normal height (200267456010854: SO 12.2" HTDEAD 66 / HTIH 51) kept
# its full 20.1 cuft bole (live 19.787945) ⇒ 2002 Standing_Dead 3.3061 vs live 3.2710 (private FVSsn FMDOUT trace).
@testset "SN broken-top snag bole at the current height (fmsvol.f:98-153) vs FVSsn_g16" begin
    for rg in ("salvage", "simfire")
        @test !_casex("200267456010854", rg).crashed
        @test _cellsx("200267456010854", rg, "FVS_Carbon", ("Standing_Dead", "Total_Stand_Carbon"), ("2002",)) == 0
    end
end

# FMSNAG calls FMSNGHT for every standing snag pool even at SN's HTX=0 (sn/fmvinit.f:1089): HTSNEW = HTCURR and
# fmsnght.f:164 `IF (HTSNEW .LT. 1.5) HTSNEW = 0.0`, so a snag under 1.5 ft breaks to fuel (CWD2) and its density is
# zeroed (fmsnag.f:262-270) in its first FMSNAG year. jl skipped FMSNGHT without SNAGBRK, so 200267456010854's 0.1"/1.01-ft
# sp62 snag stood on at 0.239/ac (live gone by 2008, private FMSNAG trace). Also each report now recomputes the eastern
# FMSVOL at (DBHS, HTDEAD) with (SNVIS+SNVIH)·V2T, V2T pre-divided by 2000 (fmdout.f:139-155, fmvinit.f:1094).
@testset "SN FMSNGHT <1.5-ft snag break at HTX=0 (fmsnght.f:164, fmsnag.f:262-270) vs FVSsn_g16" begin
    @test _cellsx("200267456010854", "salvage", "FVS_Carbon", ("Standing_Dead",), ("2002", "2007", "2012", "2017", "2022")) == 0
end

# FMMAIN (FMCBA's TBA/TOTCRA, FMCADD's litterfall/breakage/crown lift) walks DO I=1,ITRN over the list GRADD hands it —
# in a tripling cycle the TRIPLEd list: originals at FMPROB=PROB·.60, then each record's .25/.15 copies (fmcba.f:189-203,
# fmcadd.f). jl ran both on its untripled list, so the sums rounded differently (MEASURED private FMCBA/FMDOUT traces,
# 200267456010854: sp74 TBA 28.677080 live vs 28.677082 ⇒ 6-12" fuel 1.0099999; 2007 litter 3.1324124 vs 3.1324131).
@testset "SN FFE tripled FMPROB walk in FMCBA/FMCADD (fmcba.f:189-203) vs FVSsn_g16" begin
    @test _cellsx("200267456010854", "salvage", "FVS_Carbon", ("Forest_Down_Dead_Wood", "Forest_Floor"), ("2007", "2012")) == 0
end

# fmdout.f:316-343: CWDVOL(I,J,K,L) = CWD(I,J,K,L)·2000/CWDDEN(K,L) per pile/decay class, THEN summed over I and L. jl converted
# the class-summed biomass (MEASURED 200267456010854 2002 6-12" hard 80.929482 vs live 80.929489). Live golden: FVSsn_g16 on the
# tiered SALVAGE key plus DWDVLOUT/DWDVLDB (FVS_Down_Wood_Vol rows 2002, 2007).
@testset "SN FVS_Down_Wood_Vol per-pool volume sum (fmdout.f:316-343) vs FVSsn_g16" begin
    fx = fixture_dir("SN"); d = mktempdir()
    cp(joinpath(fx, "stands.db"), joinpath(d, "stands.db"))
    key = read(joinpath(fx, "200267456010854_salvage.key"), String)
    key = replace(key, "SALVAGE          2.0       0.0     999.0       0.9\n" =>
                       "SALVAGE          2.0       0.0     999.0       0.9\nDWDVLOUT\n", "POTFIRDB\n" => "POTFIRDB\nDWDVLDB\n")
    key = join([l == "out.db" ? joinpath(d, "out.db") : l == "stands.db" ? joinpath(d, "stands.db") : l
                for l in split(key, '\n')], '\n')
    write(joinpath(d, "x.key"), key)
    FVSjl.run_keyfile(joinpath(d, "x.key"); variant = FVSjl.Southern())
    h, rows = db_table_rows(joinpath(d, "out.db"), "FVS_Down_Wood_Vol")
    live = Dict("2002" => (219.55130004882812, 34.45512771606445, 80.92948913574219, 80.92948913574219, 0.0, 0.0, 0.0,
                           415.86541748046875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "2007" => (77.14576721191406, 47.670494079589844, 151.08570861816406, 48.984920501708984, 0.0, 0.0, 0.0,
                           324.88690185546875, 151.7133331298828, 44.99720764160156, 148.01951599121094, 65.89952087402344,
                           0.0, 0.0, 0.0, 410.62957763671875))
    iy = findfirst(==("Year"), h); i0 = findfirst(==("DWD_Volume_0to3_Hard"), h)
    got = Dict(string(r[iy]) => Tuple(parse(Float64, string(r[i0 + k])) for k in 0:15) for r in rows)
    for y in ("2002", "2007")
        @test got[y] == live[y]
    end
end

# fmsadd.f:300-306: FMSADD books each mortality snag's crown with UNFIRE = SNGNEW − FIRKIL, and FMKILL(2) clears FIRKIL only
# after FMSADD — so in a fire cycle, where SNGNEW = WK2 − FIRKIL (fmkill.f:119-123), the crown density is WK2 − 2·FIRKIL.
# jl scheduled WK2 − FIRKIL (MEASURED private FMSCRO trace, 200267456010854 SIMFIRE 2007: sp74 9.66" 0.0623765 vs live
# 0.0369589 = 0.08779417 − 2·0.025417633) ⇒ 2012 Standing_Dead 1.4473 vs live 1.4302.
@testset "SN fire-cycle mortality crowns at SNGNEW−FIRKIL (fmsadd.f:300-306) vs FVSsn_g16" begin
    @test !_casex("200267456010854", "simfire").crashed
    @test _cellsx("200267456010854", "simfire", "FVS_Carbon", ("Standing_Dead",), ("2007", "2012", "2017", "2022")) == 0
end

# fmsnag.f:182-190/200-214: FMSNAG STORES each record's post-burn rates PBFRIS/PBFRIH in the burn year and the year after
# ((IYR−BURNYR) ≤ 1, from that year's DENTTL; a record whose HARD flag has flipped takes PBFRIH = PBFRIS) and reuses them
# through PBTIME on the actual IYR. jl recomputed them every year from the current density, on the cycle-start year, and
# ignored HARD (MEASURED private CWD1 trace, 200267456010854 SIMFIRE 2007: the 1990 SO 8.9" snag fell 1.8242 vs live 2.1741).
@testset "SN post-burn snag fall rates stored per record (fmsnag.f:182-214) vs FVSsn_g16" begin
    @test _cellsx("200267456010854", "simfire", "FVS_Carbon", ("Forest_Down_Dead_Wood", "Forest_Floor"), ("2007", "2012")) == 0
end

# fmsout.f:122-172 / fmssum: the snag reports read FMSNAG's HARD flag, which FMSADD sets TRUE and FMSNAG flips only at the end of
# a year it ran — so the inventory-year report shows every input snag HARD (jl flipped the 1990 snags soft already in 2002),
# and each record's volume is FMSVOL(II,HTIx) on (DBHS, HTDEAD) cut at the current height with the 0.005454154·HTDEAD floor
# (jl used a (DBH, HTIH) tree's merch: SD 5.0" 0 vs live 2.2976; SO 8.9" 96.29 vs 109.53). Live golden: FVSsn_g16 on the
# tiered SALVAGE key plus SNAGOUT/SNAGSUM/SNAGOUDB/SNAGSUDB.
@testset "SN snag reports: inventory-year HARD flag + FMSVOL volume (fmsout.f:122-172) vs FVSsn_g16" begin
    fx = fixture_dir("SN"); d = mktempdir()
    cp(joinpath(fx, "stands.db"), joinpath(d, "stands.db"))
    key = read(joinpath(fx, "200267456010854_salvage.key"), String)
    key = replace(key, "SALVAGE          2.0       0.0     999.0       0.9\n" =>
                       "SALVAGE          2.0       0.0     999.0       0.9\nSNAGOUT\nSNAGSUM\n",
                  "POTFIRDB\n" => "POTFIRDB\nSNAGOUDB\nSNAGSUDB\n")
    key = join([l == "out.db" ? joinpath(d, "out.db") : l == "stands.db" ? joinpath(d, "stands.db") : l
                for l in split(key, '\n')], '\n')
    write(joinpath(d, "x.key"), key)
    FVSjl.run_keyfile(joinpath(d, "x.key"); variant = FVSjl.Southern())
    col(h, r, c) = parse(Float64, string(r[findfirst(==(c), h)]))
    h, rows = db_table_rows(joinpath(d, "out.db"), "FVS_SnagSum")
    live = Dict(2002 => (36.1082763671875, 36.1082763671875, 0.0, 0.0, 36.1082763671875),
                2007 => (42.9259033203125, 42.9259033203125, 8.559402465820312, 8.559402465820312, 51.48530578613281))
    for r in rows
        y = Int(col(h, r, "Year")); haskey(live, y) || continue
        @test Tuple(col(h, r, c) for c in ("Hard_snags_class1", "Hard_snags_total", "Soft_snags_class1", "Soft_snags_total",
                                             "Hard_soft_snags_total")) == live[y]
    end
    h, rows = db_table_rows(joinpath(d, "out.db"), "FVS_SnagDet")
    det = Dict(("SD", 1) => (19.0, 0.0, 2.2976343631744385, 0.0, 12.036091804504395, 0.0),
               ("SO", 1) => (49.999996185302734, 0.0, 109.5284423828125, 0.0, 12.036091804504395, 0.0),
               ("SO", 2) => (51.0, 0.0, 238.16952514648438, 0.0, 12.036091804504395, 0.0))
    got = Dict((strip(string(r[findfirst(==("SpeciesFVS"), h)])), Int(col(h, r, "DBH_Class"))) =>
               Tuple(col(h, r, c) for c in ("Current_Ht_Hard", "Current_Ht_Soft", "Current_Vol_Hard", "Current_Vol_Soft",
                                             "Density_Hard", "Density_Soft")) for r in rows if Int(col(h, r, "Year")) == 2002)
    @test got == det
end

# fmfout.f:326-333 TOTCLS(KSP,ICLS) = TOTCLS(KSP,ICLS) + CURKIL(I) + FMPROB(I): gfortran adds left to right, (TOTCLS+CURKIL)+FMPROB;
# jl grouped CURKIL+FMPROB first (MEASURED 200267456010854 SIMFIRE 2007 Total_class2 8.7718506 vs live 8.7718487).
@testset "SN FVS_Mortality class totals in fmfout.f's order (fmfout.f:326-333) vs FVSsn_g16" begin
    @test count(m -> m.file == "FVS_Mortality", _casex("200267456010854", "simfire").ms) == 0
end

# sn/fmburn.f:586-589 re-runs FMCBA once BURNYR=IYR BEFORE FMEFF reduces FMPROB (and only for VARACD='SN'); jl re-ran it after the
# kills, so FMSNFT typed the fire-thinned stand on the survivors' pine share (MEASURED 216786838010854 SIMFIRE 2009: IFFEFT 3
# pine/hardwood ⇒ FULIV 0.35 vs live hardwood/pine 0.04 ⇒ Shrub_Herb 0.175 vs 0.02).
@testset "SN post-burn FMCBA re-run on the pre-kill list (sn/fmburn.f:586-589) vs FVSsn_g16" begin
    @test _cellsx("216786838010854", "simfire", "FVS_Carbon", ("Forest_Shrub_Herb", "Total_Stand_Carbon"), ("2009",)) == 0
end

# fmtret.f:378-389 builds SMALL one pool at a time (SMALL = SMALL + CWD(I,J1,K,L), J1=1..3, then litter); jl added the pooled
# (1+2+3+10) sum, so the FMDYN fuel-model weights moved (MEASURED private FMFINT trace, 238813815010854 SIMFIRE 2005: FM9 weight
# 0.4597697 vs live 0.45976925 ⇒ flame 3.0070512 vs 3.0070524).
@testset "SN SMALL fuel accumulation order (fmtret.f:378-389) vs FVSsn_g16" begin
    @test _cellsx("238813815010854", "simfire", "FVS_BurnReport", ("Flame_length", "Scorch_height"), ("2005",)) == 0
    @test count(m -> m.file == "FVS_Mortality", _casex("238813815010854", "simfire").ms) == 0
end

# sn/regent.f:284-287: a seedling that stays under breast height (HK = H+HTG ≤ 4.5) gets DG(K)=0 and DBH(K)=D+0.001·HK set
# directly; jl pushed 0.001·HK through UPDATE's DBH += DG/BARK (MEASURED 830602414290487 PLANT 2035: planted PI DBH
# 0.1052496 vs live 0.1043516 ⇒ QMD 0.117378 vs live 0.116691).
@testset "SN REGENT sub-breast-height DBH set directly (sn/regent.f:284-287) vs FVSsn_g16" begin
    for rg in ("plant_cyc", "plant_cal")
        @test !_casex("830602414290487", rg).crashed
        @test isempty(_casex("830602414290487", rg).ms)
    end
end

# dbsstrclass.f:120-128 returns on NTREES=0 before CREATE TABLE, so a stand with no StrClass row leaves no FVS_StrClass table;
# jl created it empty (830602414290487 NONE, a bare stand: live table absent).
@testset "SN FVS_StrClass not created without rows (dbsstrclass.f:120-128) vs FVSsn_g16" begin
    @test isempty(_casex("830602414290487", "none").ms)
end

# sn/regent.f:360 DG(K) = SQRT((D*BARK)**2.0+DDS)-BARK*D — REAL**2.0 compiles to powf, which is not bit-equal to D·D here
# (MEASURED private REGENT trace, 238813815010854 cycle 3 LK record: every input bit-identical, DG 3F0E7B0E live vs
# 3F0E7B0C with (D·BARK)²). With it the whole NONE case equals live.
@testset "SN REGENT DG with REAL**2.0 as powf (sn/regent.f:360): 238813815010854 NONE all cells vs FVSsn_g16" begin
    @test isempty(_casex("238813815010854", "none").ms)
end

# dbsfmmort.f:147 `IF (TOTAL(J,8) .LE. 0) CYCLE` skips the ALL row as well: a SIMFIRE over a treeless stand writes no
# FVS_Mortality row (830602414290487 SIMFIRE 2025: live 0 rows, jl wrote an all-zero ALL row).
@testset "SN FVS_Mortality ALL row only with trees (dbsfmmort.f:147) vs FVSsn_g16" begin
    @test count(m -> m.file == "FVS_Mortality", _casex("830602414290487", "simfire").ms) == 0
end

# morts.f:717-718 (the BAMAX check, the same line in 16 variant builds): BANEW = BANEW + (0.0054542*(D+G)**2.)*P — REAL**2.
# compiles to powf, not (D+G)·(D+G). jl's shared mortality driver squared it, so a stand sitting at the BAMAX knife edge
# iterated a different number of times (MEASURED FVSbm_g16 41136808010497 COVER: Summary/StrClass 4 cells, NONE 543 cells
# ⇒ 127 with powf, all 4 COVER cells gone).
@testset "MORTS BAMAX check (D+G)**2. as powf (morts.f:717-718) vs FVSbm_g16" begin
    fx = fixture_dir("BM")
    d = mktempdir()
    txt, db, crashed, _ = run_case("BM", "41136808010497", "cover"; dir = d)
    @test !crashed
    @test isempty(compare_case("BM", "41136808010497", "cover", txt, db))
end

# FMPOFL for SN/CS (fmpofl.f:83-306 → dbsfmpf.f:121-360, dbsfmpfc.f): FVS_PotFire_East (surface flame, ACTCBH/CBD, INT(POKILL·100),
# INT(POVOLK), PSMOKE·P2T, severe/moderate fuel models) once per FMMAIN — never the final row, after FMBURN in a fire year —
# and FVS_PotFire_Cond at the first call. jl wrote the western FVS_PotFire layout with an approximate mortality/smoke.
@testset "SN FVS_PotFire_East + FVS_PotFire_Cond (fmpofl.f, dbsfmpf.f:121-360) vs FVSsn_g16" begin
    ms = _casex("216786838010854", "salvage").ms
    @test count(m -> startswith(m.file, "FVS_PotFire"), ms) == 0
end

# fmcba.f:228-229 PERCOV = (1.0-EXP(-TOTCRA/43560.))*100 — EXP is expf; jl used Julia exp, a ULP apart often enough to move WMULT
# and so FWIND (MEASURED private FMFINT trace, 238813815010854 2000 potential fire: FWIND 40875A5C vs live 40875A5B ⇒ spread
# rate and flame ULPs; 205045340010854 SIMFIRE 2018 potential flame).
@testset "FFE PERCOV with expf (fmcba.f:228-229) vs FVSsn_g16" begin
    @test count(m -> m.file == "FVS_PotFire_East", _casex("205045340010854", "simfire").ms) == 0
end

# fmoldc.f (fmmain.f:268) closes the fire cycle's FMMAIN after FMBURN, so OLDCRL = HT·FMICR/100 uses the fire-shortened crown;
# the next FMSDIT crown lift (fmsdit.f:103-118) of a scorched survivor starts from that crown base. jl kept the pre-fire ICR
# (MEASURED private FMCADD trace, 238813815010854 SIMFIRE 2010: crown-lift size 1 0.0116946 vs live 0.0131616 ⇒ DDW 2015+).
@testset "FFE FMOLDC with the post-fire FMICR (fmoldc.f, fmsdit.f:103-118) vs FVSsn_g16" begin
    ms = _casex("238813815010854", "simfire").ms
    @test _cellsx("238813815010854", "simfire", "FVS_Carbon", ("Forest_Down_Dead_Wood", "Forest_Floor"), ("2020",)) == 0
    # 2015 was 3.936955 vs live 3.936730; what is left there is a 1-ULP pool residual
    for m in ms
        (m.file == "FVS_Carbon" && m.col in ("Forest_Down_Dead_Wood", "Forest_Floor") && m.year == "2015") || continue
        g = parse(Float64, m.gold); j = parse(Float64, m.got)
        @test abs(j - g) <= 4e-7 * abs(g)
    end
end

# fmcwd.f:176-183 (CWD1) and :230-236 (CWD2): TVOLI = FMSVL2(SP,DBHS,HTDEAD,-1.,…,'D',.FALSE.) recomputed at each fall —
# MAX(0.005454154·HTDEAD, MCF) with no top-kill; jl divided the stored per-stem tons back by V2T/2000, a round trip off by an
# ULP (MEASURED 161035853010854 SALVAGE 2003 DDW 1.3504815 vs live 1.3504814).
@testset "SN snag-fall TVOLI from FMSVL2 (fmcwd.f:176-183) vs FVSsn_g16" begin
    @test _cellsx("161035853010854", "salvage", "FVS_Carbon", ("Forest_Down_Dead_Wood",), ("2003",)) == 0
    @test count(m -> m.file == "FVS_PotFire_East", _casex("216786838010854", "simfire").ms) == 0
end

# cratet.f:150-195 → dense.f:83-87,244: CRATET's backdating DENSE runs over live+dead (ITRN still includes the inventory dead)
# and PCTILE leaves each dead record a PCT (IMC-9 snags carry WK5=0, so they take the cumulative share below them); later
# DENSEs are live-only, so FVS_TreeList reports that PCT for the cycle-0 dead rows. jl printed 0 (200267456010854: live
# 0.635/14.38/29.29).
@testset "SN cycle-0 dead-record PCT from CRATET's DENSE (dense.f:83-87,244) vs FVSsn_g16" begin
    @test isempty(_casex("200267456010854", "none").ms)
end

# fmcba.f TBA per variant: bc/cs/em/ie/kt/sn BA1·FMPROB (BA1 = 3.14159·(DBH/24)²), so FMPROB·5.454153E-03·DBH**2, every other
# build FMPROB·DBH·DBH·0.0054542. PRCL = TBA/TOTBA splits the initial dead fuel into decay classes, so the form shows in the
# first report (MEASURED FVSbm_g16 12827438010497 SALVAGE 2005 Forest_Floor 8.2318602 vs live 8.2318592 with the SN form).
@testset "FMCBA species BA in each build's form (fmcba.f TBA/FMTBA) vs FVSbm_g16" begin
    d = mktempdir()
    txt, db, crashed, _ = run_case("BM", "12827438010497", "salvage"; dir = d)
    @test !crashed
    ms = compare_case("BM", "12827438010497", "salvage", txt, db)
    @test count(m -> m.file == "FVS_Carbon" && m.col == "Forest_Floor" && m.year == "2005", ms) == 0
end

# gradd.f:192 DENSE (post-UPDATE, before ESNUTR) supplies the PCCF that ESGENT→REGENT(LESTB) reads for the new regen crown
# CR = 0.89722−0.0000461·PCCF (regent.f:178); nothing re-DENSEs between ESUCKR and ESTAB. jl read the start-of-cycle point CCF
# (MEASURED private REGENT trace, 161035853010854 PLANT cycle 2: PCCF 105.02 live vs 78.60 ⇒ planted crowns 87/86 vs 88/87).
@testset "SN regen crown from the post-growth PCCF (gradd.f:192, regent.f:178) vs FVSsn_g16" begin
    for (st, rg) in (("161035853010854", "plant_cyc"), ("238813815010854", "plant_cyc"))
        @test isempty(_casex(st, rg).ms)
    end
end

# algslp.f: Y(I)+((Y(I+1)-Y(I))/(X(I+1)-X(I)))*(XX-X(I)) — the slope first, then ×(XX−X(I)); an XX on a breakpoint takes the
# NEXT interval. The western FULIVE/FULIVI and FUINIE/FUINII interpolations by PERCOV (and SN's FULIV2 by shrub age) used
# (Y2−Y1)·(X−X1)/(X2−X1) (MEASURED private FVSie FMCBA trace, 3027007010690 1996: STFUEL now bit-equal to live; Shrub_Herb
# 1996/2006/2026 were 1 ULP off).
@testset "FFE ALGSLP interpolation in algslp.f's association vs FVSie_g16" begin
    d = mktempdir()
    txt, db, crashed, _ = run_case("IE", "3027007010690", "salvage"; dir = d)
    @test !crashed
    ms = compare_case("IE", "3027007010690", "salvage", txt, db)
    @test count(m -> m.file == "FVS_Carbon" && m.col == "Forest_Shrub_Herb", ms) == 0
end

end # module
