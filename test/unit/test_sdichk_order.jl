# test_sdichk_order.jl — SDICHK runs LAST in CRATET (every variant's cratet.f: after the CROWN dub, the DGDRIV/REGENT
# calibration and the final DENSE), so the inventory crown dub reads the UNRESET SDImax. Synthetic over-dense stands
# (3000 TPA of 4.5-6.3" trees, half the crowns missing) trip FVS41 in live FVSbm_g16 / FVSsn_g16; goldens = the live
# 1990 FVS_TreeList PctCr. Resetting SDImax before the dub gave BM 5/10 and SN 4/10 dubbed crowns 1-15 points off.
using FVSjl, Test

@testset "SDICHK after the CRATET crown dub (live FVSbm/FVSsn over-dense inventory)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "sdichk")
    for (v, var) in (("bm", FVSjl.BlueMountains()), ("sn", FVSjl.Southern()))
        live = Dict(parse(Int, split(l, ',')[1]) => parse(Int, split(l, ',')[2])
                    for l in readlines(joinpath(fx, "$(v)_sdichk.live.csv"))[2:end])
        dir = mktempdir()
        cp(joinpath(fx, "$(v)_sdichk.key"), joinpath(dir, "s.key")); cp(joinpath(fx, "$(v)_sdichk.tre"), joinpath(dir, "s.tre"))
        got = cd(dir) do
            FVSjl.run_keyfile("s.key"; variant = var, output = :sum)
            db = FVSjl.SQLite.DB(joinpath(dir, "OUT.db"))
            d = Dict(Int(r[:TreeIndex]) => Int(r[:PctCr])
                     for r in FVSjl.DBInterface.execute(db, "SELECT TreeIndex, PctCr FROM FVS_TreeList WHERE Year = 1990"))
            FVSjl.SQLite.close(db); d
        end
        @test length(live) == 10
        @test got == live
    end
end
