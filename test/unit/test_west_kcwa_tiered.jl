# test_west_kcwa_tiered.jl — KT / WS / CA / AK fixes vs the LIVE western tiered goldens (test/fixtures/tiered/<v>: stands.db +
# <stand>_<regime>.key + the live .sum / DBS tables of FVSkt_clean / FVS<v>_g16) and live keyfile goldens. Each testset names
# the Fortran it follows and the case it was measured on; a testset only asserts what its mechanism moved.
module WestKcwaTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

function _case(v, cn, r)
    d = mktempdir()
    txt, db, crashed, _ = run_case(v, cn, r; dir = d)
    return (txt = txt, db = db, crashed = crashed, ms = compare_case(v, cn, r, txt, db))
end
_yr(m) = tryparse(Int, m.year)
_before(ms, y) = [m for m in ms if (k = _yr(m); k !== nothing && k < y)]

# kt/htgf.f, dgf.f, crown.f, morts.f, regent.f, ktfctr.f: REAL*4 EXP/ALOG/** are glibc expf/logf/powf and 0.985**K is
# libgcc __powisf2 (doctrine §4). Julia's own Float32 exp/log/^ put the 1999 HTG of four records 1-2 ULP off live
# (3021216010690: HtG 2.92662287 live / 2.92662311 jl), compounding through every later cycle.
@testset "KT growth kernels in glibc expf/logf/powf (kt/htgf.f:107-114, dgf.f, crown.f, morts.f)" begin
    c = _case("KT", "3021216010690", "none")
    @test !c.crashed
    @test isempty(_before(c.ms, 2000))
end

# kt/regent.f:560-562: a sub-4.5 ft record gets DBH(K)=0.1+DIAM·0.01+HK·0.001 with DG(K)=0 in its OWN slot — the tripled
# copies too (kt/dgdriv.f:253/261 put the pre-REGENT D in their slots; TRIPLE does not copy DBH). jl gave the copies the
# central DBH plus an equivalent increment (TreeList DG 7.9E-4 vs live 0) and a direct-set central's DBH for the others.
@testset "KT REGENT sub-4.5 ft tripled copies keep their own DBH(K) (kt/regent.f:560-562)" begin
    c = _case("KT", "4718785010690", "none")
    @test !c.crashed
    @test count(m -> m.file == "FVS_TreeList" && m.col == "DG" && m.year == "2014", c.ms) == 0
    c = _case("KT", "1627650292290487", "none")
    @test isempty(_before(c.ms, 2050))          # 8325 cells from 2030 before
end

# dense.f:168-229 sums CCFT·P species-major in IND1 order into RELDSP(ISPC), then RELDT=RELDT+RELDSP(ISPC). The flat
# record-order sum put KT's RELDEN 1 ULP off (3021216010690 cycle 2: DGF CCF2 11904.7773 live / 11904.774 jl), and the
# tripled copies' DG with it (TreeIndex 114 DG 2.35536575 live / 2.3553638 jl).
@testset "KT stand CCF summed per species (dense.f:168-229 RELDSP)" begin
    c = _case("KT", "3021216010690", "none")
    @test count(m -> m.file == "FVS_TreeList" && m.col == "DG" && m.year == "2009", c.ms) == 0
end

# kt/morts.f:201-203 CIOBDS=(2.0*D*G+G*G); SD2SQ=SD2SQ+P*(D*D+CIOBDS) — D²+(2DG+G²), not (D²+2DG)+G².
@testset "KT MORTS SD2SQ association (kt/morts.f:201-203)" begin
    c = _case("KT", "4718785010690", "none")
    @test isempty(_before(c.ms, 2034))          # first divergence 2024 before
end

# kt/morts.f:288-307: X = XMORT inside the MORTMULT DBH window, then the ESTAB best-tree IESTAT guard. jl ran X=1 for every
# record. ktt01 with MORTMULT 1990 all-species ×2.5 on 5-15" and 2010 DF ×0.4: every .sum row equals live FVSkt_clean
# (all 7 rows off without the window).
@testset "KT MORTS MORTMULT window (kt/morts.f:288-292) vs live FVSkt" begin
    dir = joinpath(@__DIR__, "..", "fixtures", "kootenai", "ktt01_mortmult")
    out = mktempdir()
    for f in ("ktmm.key", "ktmm.tre"); cp(joinpath(dir, f), joinpath(out, f)); end
    rows(t) = [l for l in split(t, '\n') if occursin(r"^\s?\d{4} ", l)]
    jl = rows(cd(() -> FVSjl.run_keyfile("ktmm.key"; variant = FVSjl.Kootenai(), output = :sum), out))
    lv = rows(read(joinpath(dir, "ktmm.live.sum"), String))
    @test length(jl) == length(lv)
    for (g, j) in zip(lv, jl); @test rstrip(j) == rstrip(g); end
end

# kt/cratet.f:148-151 IND=IND1; RDPSRT(ITRN,DBH,IND,.FALSE.) is the order the :184 backdating DENSE's PCTILE ranks the
# calibration percentile in (as every other variant's cratet.f, bm_cratet166_ind). KT re-sorted afresh (.TRUE.), which put
# 196396140020004's two 13.4" trees in the other order: DGF BAL 244.979 vs live 239.103 ⇒ the DO-220 dub WK2 2.4734 vs
# 2.4888 ⇒ cycle-1 MORTS WK1 0.4768 vs 0.4840 ⇒ the kill.
@testset "KT calibration PCT in cratet.f's IND1 RDPSRT(.FALSE.) order (kt/cratet.f:148-151)" begin
    c = _case("KT", "196396140020004", "none")
    @test isempty(_before(c.ms, 2032))          # 18712 cells from 2022 before
end

# kt/bratio.f: BRATIO = BKRAT(IS), exactly. jl's linear (0 + BKRAT·D)/D is 1 ULP off for many D (22404917010497 REGENT
# BARK 0.915000021 live / 0.91499996 jl ⇒ small-tree DG/DBH in the 7th digit, compounding).
@testset "KT bark ratio = BKRAT(IS) (kt/bratio.f)" begin
    for cn in ("22404917010497", "31445998010690", "720724089290487", "31454370010690", "251193174489998", "750215552290487")
        c = _case("KT", cn, "none")
        @test !c.crashed
        @test isempty(c.ms)
    end
end

# ca|ws/fmmain.f run FMSDIT and the annual FMSNAG/FMCWD/FMCADD loop like every western FFE build, and ca|ws/fmsvol.f are
# byte-identical to so's (FMSVOL VOL2HT = MAX(0.005454154·H, TCF) on the variant's NATCRS cubic). CA/WS were off both the
# FFE-dynamics gate and the western volume layer: no inventory snags (WS 7689156010901 2006 Standing_Dead 10.61 live / 0
# jl, Belowground_Dead 3.37 / 0), then the eastern R8-Clark bole (0.03).
@testset "CA/WS FFE snags + down wood (ca|ws/fmmain.f, fmsvol.f) vs live at the inventory year" begin
    c = _case("WS", "7689156010901", "simfire")
    @test !c.crashed
    @test count(m -> m.file == "FVS_Carbon" && m.year == "2006", c.ms) == 0
    c = _case("CA", "374401353489998", "salvage")
    @test count(m -> m.file == "FVS_Carbon" && m.year == "2015" && m.col in ("Standing_Dead", "Belowground_Dead"), c.ms) == 0
end
# ak/esnutr.f:282-290: after a removal the LAUTAL automatic tally fires on XTES=MAX(ONTREM/ONTCUR, OCVREM/OCVCUR) — AK's own
# esnutr.f forms the same ratio as estb/esnutr.f. jl only stashed it for IE/EM/KT, so a thinned AK stand ran the ingrowth
# tally (ITPP capped at MAXING=3, not MAXTPP=21): FVSak_g16 10709344010497 thinbba 2016 regen 76.5 TPA live / 33.3 jl.
@testset "AK post-removal LAUTAL tally (ak/esnutr.f:282-290)" begin
    c = _case("AK", "10709344010497", "thinbba")
    @test !c.crashed
    @test count(m -> m.file == "sum", c.ms) == 0
end

# fortyp.f:1166 IF(ICYC.EQ.0) IIFORTP=IFORTP — the inventory forest type, which ak/estab.f:296-302 falls back to when the
# current type is not one of its 17 named types (e.g. 999 nonstocked). jl never set it (0 ⇒ IFT 14 ⇒ TPP 0): FVSak_g16
# 10706662010497 thinbba, the 2037 ingrowth tally on a 999 stand: live IFT 6 (inventory 301), 2047 TPA 571 / jl 52.
@testset "AK ESTAB falls back to the inventory forest type IIFORTP (fortyp.f:1166, ak/estab.f:296-302)" begin
    c = _case("AK", "10706662010497", "thinbba")
    @test !c.crashed
    @test count(m -> m.file == "sum", c.ms) == 0
end

# ca/htgf.f, htcalc.f, dgf.f, crown.f, htdbh.f, regent.f: REAL*4 EXP/ALOG/SIN/COS/** are glibc expf/logf/sinf/cosf/powf
# (doctrine §4), (D*BARK)**2.0 is powf at -O0, and htcalc.f's constant ALOG(50.)/EXP(50.0*(-0.0440853))/50.0**1.51744
# are folded (correctly rounded) by gfortran. 23742358010900: DF HTG 9.54881763 live / 9.54880428 jl at 2011.
@testset "CA growth kernels in glibc expf/logf/powf (ca/htgf.f, htcalc.f, dgf.f)" begin
    c = _case("CA", "23742358010900", "none")
    @test !c.crashed
    @test count(m -> m.file == "FVS_TreeList" && m.col == "HtG" && m.year == "2011", c.ms) == 0
end

# ca/fmcba.f:579-745 sets CA's FFE decay rates at the first FFE year (R5: 0.025/0.0125 woody; R6: the Oregon table ×
# DKRADJ(CAHMC,CAWMD); litter 0.5, duff 0.002; ×DCYMLT from the Dunning code for KODFOR<600); ws/fmvinit.f + fmcba.f:541-581
# the same California table × DCYMLT. jl decayed both with the SN table (woody 0.07-0.11, litter 0.65): bare CA stands'
# initial down wood fell ~4x too fast (23999387010900 2012 Forest_Down_Dead_Wood 0.255 live / 0.296 jl).
@testset "CA/WS FFE decay rates (ca/fmcba.f:579-745, ws/fmvinit.f:120-150, ws/fmcba.f:541-581)" begin
    for cn in ("23999387010900", "647485752126144", "249049837489998"), r in ("simfire", "salvage")
        c = _case("CA", cn, r)
        @test !c.crashed
        @test isempty(c.ms)
    end
    c = _case("WS", "23771657010900", "salvage")          # 1-ULP down-wood residue remains (Float32 decay order)
    rel(m) = (a = tryparse(Float64, m.gold); b = tryparse(Float64, m.got); (a === nothing || b === nothing) ? Inf : abs(a - b) / max(abs(a), 1e-9))
    @test count(m -> m.file == "FVS_Carbon" && rel(m) > 1e-5, c.ms) == 0
end

# ws/dgf.f:618 CR=FLOAT(ICR(I))/100. — a divide; jl multiplied by 0.01 (0.39999998 vs 0.400000006 for ICR 40), putting the
# DGF WK2 2 ULP off (7689398010901 cycle 1 WF record 16: WK2 2.61435175 live / 2.6143513 jl ⇒ DG ⇒ HTG).
@testset "WS DGF crown ratio CR = ICR/100 (ws/dgf.f:618)" begin
    for cn in ("7689398010901", "446847010497", "23762384010900", "850400255290487")
        c = _case("WS", cn, "none")
        @test !c.crashed
        @test isempty(c.ms)
    end
end

# dgdriv.f DO 155/160 (ws/dgdriv.f:459-525): the DG self-calibration sums run over each species' IND1 slice in
# REAL*4 with glibc EXP/ALOG and SNXX=SNXX+P*EDDS*EDDS left-to-right; jl summed in record order with Julia exp/log
# and P*EDDS^2 ⇒ CONSPP 1.48652804 jl vs 1.48652816 live (15353585010497 cycle 1, 1 ULP) ⇒ DG drift.
@testset "DG self-calibration sums in IND1 order, glibc EXP/ALOG (ws/dgdriv.f:459-525)" begin
    for cn in ("15353585010497", "23771657010900", "248626203489998")
        c = _case("WS", cn, "none")
        @test !c.crashed
        @test isempty(c.ms)
    end
end

# kt/blkdat.f:113-128 OCURNF(IFO,sp) (= ie/blkdat.f species 1-11; KT row 11 all 0) multiplies every ESPADV/ESPXCS/ESPSUB
# species probability. KT had no table (1.0 default), so 3021216010690 (IFO 10: OCURNF 1 1 1 1 0 0 1 1 1 0 0) booked PP
# ingrowth live zeroes, and every later regen record's species/height shifted (17,382 cells from 2009).
@testset "KT AUTOES OCURNF national-forest occupancy (kt/blkdat.f:113-128, espadv.f)" begin
    c = _case("KT", "3021216010690", "none")
    @test !c.crashed
    @test count(m -> m.file == "FVS_TreeList" && m.year == "2009" && m.col in ("SpeciesFVS", "Ht", "HtG", "DBH"), c.ms) == 0
    @test count(m -> m.file == "sum", c.ms) == 0
end

end # module
