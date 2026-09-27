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

# EM DVEW hardwood volume: dvest.f:142 `VOL(2)=ANINT(VOL(2))` rounds every direct-volume-estimator board-foot volume
# (R1KEMP/R2OLDV) to the nearest board foot — live CW 13.4"×47' BdFt 62, jl 62.08; and r1kemp.f:363's polynomial takes
# the integer powers first, C2*DBHOB**2 + C3*DBHOB**3 (jl had (C2*D)*D and ((C3*D)*D)*D ⇒ 1-ULP TCuFt on 5-9.5" CW/AS).
@testset "EM DVE hardwood volume (dvest.f ANINT, r1kemp.f powers) vs FVSem_g16" begin
    r4 = _run_stand("2999215010690")
    @test _treelist_diffcells("2999215010690", r4.db, ("1988", "1998")) == 0
    r5 = _run_stand("888512560290487")
    @test _treelist_diffcells("888512560290487", r5.db, ("2031", "2041", "2051", "2061", "2071"); cols = ("BdFt", "TCuFt", "MCuFt")) == 0
end

# estab.f:579 FTEMP=1/(1+EXP(-(PN+ESB-ESB1(NCOUNT)))) evaluates (PN+ESB)−ESB1 left to right (jl added the precomputed
# ESB−ESB1), and the plot aspect terms are COS/SIN = glibc cosf/sinf: PROB1 was 1 ULP off ⇒ every ingrowth record's
# PROB=(ESPROB*300)/DUPNPT 1 ULP off (196378260020004: 2.0518632 vs live 2.0518634 on 6 cohort records).
# The per-point BAAA/OVER feeding ESTOCK/ESPADV are dense.f's IND1-order BATREE*PI/GROSPC sums, and BAAINV (ESB1) is
# esfltr.f's record-order 0.005454154*D*D*PROB*PIX — both 1-4 ULP off in jl's reassociated forms (684750664126144
# point-3 BAAA 432083F6 vs 432083F8; IE 3356357010690 point-4 ESB1 C0084738 vs C0084734).
@testset "EM/IE AUTOES PROB1 at estab.f precision ((PN+ESB)−ESB1, cosf/sinf, BAAA, BAAINV) vs FVSem/FVSie_g16" begin
    @test _treelist_diffcells(STAND, _run.db, ("2012", "2022", "2032", "2042", "2052", "2062")) == 0
    r6 = _run_stand("684750664126144")
    @test _treelist_diffcells("684750664126144", r6.db, ("2028", "2038", "2048", "2058", "2068")) == 0
    # shared estb code: the IE stand whose point-4 ESB1 was 4 ULP off
    dI = mktempdir(); _, dbI, crI, _ = run_case("IE", "3356357010690", "none"; dir = dI)
    @test !crI
    goldI = _keyed(read_csv(joinpath(fixture_dir("IE"), "3356357010690_none.FVS_TreeList.csv"))...)
    gotI  = _keyed(db_table_rows(dbI, "FVS_TreeList")...)
    @test count(k -> !haskey(gotI, k) || any(!_eq(v, get(gotI[k], c, "")) for (c, v) in goldI[k] if c != "StandID"),
                collect(keys(goldI))) == 0
end

# ICL5 = the INPUT habitat code (dbsstandin.f:590-593 IFIX(PV_CODE); grinit.f:200) — esplt2.f brackets it for the AUTOES
# habitat group, not the translated KODTYP: 3006831010690 PV_CODE 9999999 → KODTYP 260 but ICL5 bracket group 16 (AF
# series, ISER 5): live books AF/ES/LP where jl booked DF/LP/PP; 2999215010690 PV_CODE 356 → 250, group 4 not 3.
@testset "EM AUTOES habitat group from ICL5 (esplt2.f) vs FVSem_g16" begin
    r7 = _run_stand("3006831010690")
    @test _treelist_diffcells("3006831010690", r7.db, ("1999", "2009", "2019", "2029", "2039")) == 0
    r8 = _run_stand("2999215010690")
    @test _treelist_diffcells("2999215010690", r8.db, ("1988", "1998", "2008", "2018", "2028", "2038")) == 0
end

# esnutr.f:264-289 LAUTAL (estb, shared by IE and EM): a thin removing ≥THRES1 of the TPA or cubic volume schedules the
# NTALLY=1 disturbance tally dated at the thin. jl captured the removal fraction for IE only, so an EM thin fell through
# to the ingrowth tally (MEASURED FVSem_g16 196378260020004 THINBBA 2022: XTPA 0.950 ⇒ NTALLY 1; live .sum 2032 TPA 182,
# jl 127). The whole .sum must match live through the tally cycle.
@testset "EM LAUTAL post-thin disturbance tally (esnutr.f) vs FVSem_g16" begin
    d = mktempdir()
    txt, db, crashed, _ = run_case("EM", STAND, "thinbba"; dir = d)
    @test !crashed
    live = sum_rows(read(joinpath(fixture_dir("EM"), "$(STAND)_thinbba.live.sum"), String))
    got = sum_rows(txt)
    row(rows, y) = (i = findfirst(r -> startswith(r, y), rows); i === nothing ? "" : rows[i])
    for y in ("2012", "2022", "2032")
        @test split(row(got, y)) == split(row(live, y))
    end
end

# LOAD (ESHAP): FIA-DB tree records carry plot site data (esplt1.f IPINFO=2 ⇒ esplt2.f:274 LOAD=1), so the first
# disturbance tally takes each plot's site prep from the plot data (1 = none) — no ESPREP default sampling. jl always
# sampled ESPREP ⇒ MECH plots (IPREP 2) where live has none. MEASURED FVSem_g16 196378260020004 thinbba @2031
# "PNONE= 0 PMECH= 0 PBURN= 0", every plot IPREP 1; jl IPREP 2 on plots 2+ (UPRE term ⇒ other heights/species).
@testset "EM/IE AUTOES site prep from plot data when LOAD=1 (esplt2.f/estab.f) vs FVSem_g16" begin
    d = mktempdir()
    txt, db, crashed, _ = run_case("EM", STAND, "thinbba"; dir = d)
    @test !crashed
    ms = compare_case("EM", STAND, "thinbba", txt, db)
    @test count(m -> m.year == "2032" && m.file == "FVS_StrClass", ms) == 0
end

# estab.f:1248-1256 (and the PLANT block :1427-1439): a booked record's DG/HTG, OLDPCT/OLDRN, WK1/WK2 are zeroed
# and MISPUTZ(ITRN,0) clears its mistletoe rating. After a thin TREDEL leaves removed records in the slots the new
# cohort reuses; jl kept their MortPA (WK2) and DMR (MEASURED FVSem_g16 196378260020004 thinbba @2032: ES020123
# MistCD live 0 / jl 3, ES020115-120 MortPA live 0 / jl 0.2-2.5) ⇒ mistletoe loss/spread on the cohort from 2042.
@testset "EM/IE booked records clear the reused slot (estab.f WK1/WK2/MISPUTZ) vs FVSem_g16" begin
    d = mktempdir()
    txt, db, crashed, _ = run_case("EM", STAND, "thinbba"; dir = d)
    @test !crashed
    ms = compare_case("EM", STAND, "thinbba", txt, db)
    @test count(m -> !(m.file == "FVS_Summary" && m.col == "MAI"), ms) == 0
end

# evtstv.f:414 (CASE DEFAULT — every variant; the eastern CASE is commented out): BCYMAI=(TOTREM+CURVOL)/AGE with
# TOTREM accumulating the INTEGER IOSUM(9,ICYC-1)=INT(OMCREM(7)/GROSPC+0.5) (disply.f:341). jl summed the per-acre
# float removal, which agrees at the .sum's F5.1 but not in FVS_Summary's REAL MAI (MEASURED FVSem_g16
# 196378260020004 thinbba 2032: live 10.668750 = 1707/160, jl 10.670339 = 1707.25/160).
@testset "MAI accumulates the integer TOTREM (evtstv.f:414) vs FVSem_g16" begin
    for cn in (STAND, "2999215010690", "684750664126144")
        d = mktempdir()
        txt, db, crashed, _ = run_case("EM", cn, "thinbba"; dir = d)
        @test !crashed
        ms = compare_case("EM", cn, "thinbba", txt, db)
        @test count(m -> m.col == "MAI", ms) == 0
    end
end

# dbstrls.f:215-216: at the inventory (ICYC=0, TEM=0) the TreeList DG is WORK1(I), which dgdriv.f:785-805 set to the
# measured increment when DG>0 .AND. HT>4.5 and to 0 otherwise — not the calibration's −1 "missing" sentinel
# (MEASURED FVSem_g16 3006831010690 1989: every no-DG record DG 0, jl −1 on 11 records).
@testset "Inventory-year TreeList DG is WORK1, not the −1 sentinel (dbstrls.f:215) vs FVSem_g16" begin
    r = _run_stand("3006831010690")
    @test !r.crashed
    @test _treelist_diffcells("3006831010690", r.db, ("1989",); cols = ("DG",)) == 0
    r6 = _run_stand("684750664126144")
    @test _treelist_diffcells("684750664126144", r6.db, ("2018",); cols = ("DG",)) == 0
end

# dbsclsum.f:66-76 builds the FVS_Climate INSERT with a list-directed WRITE: each REAL*4 reaches SQLite as
# 9-significant-digit text (0.775909066), not a bound double — jl stored the exact Float32 (0.7759090662002563), so
# every real cell differed (MEASURED FVSem_g16 196378260020004 climate: 107 of 109 FVS_Climate cells).
@testset "FVS_Climate reals pass through the list-directed text (dbsclsum.f) vs FVSem_g16" begin
    d = mktempdir()
    txt, db, crashed, _ = run_case("EM", STAND, "climate"; dir = d)
    @test !crashed
    ms = compare_case("EM", STAND, "climate", txt, db)
    @test count(m -> m.file == "FVS_Climate" && m.col in ("Viability", "SiteMult", "GrowthMult", "BA", "TPA",
                                                           "dClimMort", "ViabMort", "MxDenMult"), ms) == 0
end
end # module
