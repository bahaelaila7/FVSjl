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

# MISCOM DMFLAG: MISTOE resets it each cycle (mistoe.f:193) and sets it for a host species with SMR>0 (:267) or a MISTPINF
# host infection (misinf.f:182), even when that species has no trees; the end-of-cycle MISPRT (fvs.f:400, misprt.f:393)
# then writes FVS_DM_Stnd_Sum/Spp_Sum. AK 10707449010497 MISTPINF on species 1 (absent): live writes the 2011 rows.
@testset "DM report rows follow DMFLAG (misprt.f:393, misinf.f:182)" begin
    for cn in ("10707449010497", "10709171010497", "644809321126144", "666760633126144")
        c = _case("AK", cn, "mistletoe")
        @test count(m -> startswith(m.file, "FVS_DM_"), c.ms) == 0
    end
end
end # module
