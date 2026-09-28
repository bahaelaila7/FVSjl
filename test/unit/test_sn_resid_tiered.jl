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

end # module
