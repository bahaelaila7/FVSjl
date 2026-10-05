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

# OP 374608587489998 (no big-6 ⇒ every record FVS-native; 48 measured-DG records, 30 of them WH): op/dgf.f skipped the
# ORGANON-valid species (WK2 stale ⇒ WH calibration CSNXX=0 ⇒ COR NaN ⇒ NaN heights), had no red-alder equation, and the
# native kernels ran on Julia exp/log/^ instead of glibc expf/logf/powf / libgcc __powisf2.
@testset "OP op/dgf.f every record + RA equation + glibc kernels (dgf.f, htcalc.f, htgf.f, morts.f, crown.f, smhgdg.f)" begin
    c = _case("OP", "374608587489998", "none")
    @test !c.crashed
    @test _n(c.ms; file = "sum", before = 2025) == 0
    @test _n(c.ms; file = "FVS_TreeList", before = 2020) == 0
end

end # module
