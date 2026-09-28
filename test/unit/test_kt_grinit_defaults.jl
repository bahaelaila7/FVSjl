# test_kt_grinit_defaults.jl — kt/grinit.f:161,223 defaults ASPECT = 45., SLOPE = 30.0 (every other build: 0., 5.0).
# A KT FIA stand whose FVS_STANDINIT row has no SLOPE/ASPECT (tiered KT 22404917010497, a bare condition) runs on
# them: live FVSkt prints "STDINFO ... ASPECT AZIMUTH IN DEGREES=  45.; SLOPE=  30.%" and its ESTAB tally reads
# "SLO= 0.30 ... ASPECT= 0.785". jl gave KT the 5% / 0° of the other builds.

using Test, FVSjl
const _KG = FVSjl

@testset "KT FIA stand without SLOPE/ASPECT takes kt/grinit.f's 30% / 45°" begin
    db = joinpath(@__DIR__, "..", "fixtures", "tiered", "kt", "stands.db")
    for (v, cn, slope, asp) in ((_KG.Kootenai(), "22404917010497", 0.30f0, 45),
                                (_KG.InlandEmpire(), "22404917010497", 0.05f0, 0))
        s = _KG.StandState(v); s.plot.stand_id = cn
        _KG.load_fia_stand!(s, db, "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "")
        @test (typeof(v), s.plot.slope, s.plot.slope_raw, s.plot.aspect_deg) == (typeof(v), slope, Int32(round(slope * 100)), Int32(asp))
        @test s.plot.aspect == Float32(asp) * 0.0174533f0
    end
end
