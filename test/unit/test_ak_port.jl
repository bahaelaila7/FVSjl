# test_ak_port.jl — AK (SoutheastAlaska) establishment / small-tree growth / volume vs the live FVSak_g16 oracle.
#
# Fixture test/fixtures/southeastalaska/akplant.key = akt01's bare PLANT stand (NOTREES, ESTAB 1992, PLANT 400 WS +
# 400 SS, NOAUTOES, NOTRIPLE, 10 cycles); akplant.live.sum = FVSak_g16 on it.
#  * establishment: ak/estab.f (the tally + the per-regen-plot PLANT trees at the ESSUBH HTCALC height) and
#    ak/esgent.f → ak/regent.f REGENT(.TRUE.) grow the 100 records in their birth cycle, then the cycling
#    REGENT(.FALSE.) small-tree model (HTCALC Hegyi/Payandeh + Chapman-Richards DBH inverse) — every
#    TPA/BA/SDI/CCF/TopHt/QMD row through 2082.
#  * volume start bark: vols.f:150-151 BARK=BRATIO(DBH_start) before D=D+DG/BARK, so NVEL's DBTBH=D·(1−BARK)
#    scales the F32 profile with the start-of-cycle bark — the small-SS TCuFt 0.1 roundings (2012 231 not 240).
#  * volume start bark: vols.f:150-151 BARK=BRATIO(DBH_start) before D=D+DG/BARK, so NVEL's DBTBH=D·(1−BARK)
#    scales the F32 profile with the start-of-cycle bark — the small-SS TCuFt 0.1 roundings (2012 231 not 240).
using FVSjl, Test

function _ak_sum_rows(path)
    rows = Dict{Int,Vector{Float64}}()
    for l in eachline(path)
        f = split(l)
        (length(f) > 11 && all(isdigit, f[1]) && length(f[1]) == 4) || continue
        rows[parse(Int, f[1])] = parse.(Float64, f[3:12])      # TPA BA SDI CCF TopHt QMD TCuFt MCuFt SCuFt BdFt
    end
    return rows
end

@testset "AK PLANT stand vs live FVSak (establishment + REGENT + start-bark volume)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "southeastalaska")
    live = _ak_sum_rows(joinpath(fx, "akplant.live.sum"))
    dir = mktempdir()
    cp(joinpath(fx, "akplant.key"), joinpath(dir, "akplant.key"))
    jl = cd(dir) do
        write("jl.sum", FVSjl.run_keyfile("akplant.key"; variant = FVSjl.SoutheastAlaska(), output = :sum))
        _ak_sum_rows("jl.sum")
    end
    yrs = sort([y for y in keys(live) if y <= 2082])
    @test length(yrs) == 10
    @testset "establishment + small-tree growth (TPA..QMD)" begin
        for y in yrs
            @test get(jl, y, zeros(10))[1:6] == live[y][1:6]
        end
    end
    @testset "volume with the start-of-cycle bark" begin
        for y in yrs
            @test get(jl, y, zeros(10))[7:10] == live[y][7:10]
        end
    end
end
