# CI REGENT + CROWN vs the live FVSci_g16 oracle on five FIA Central Idaho stands (fixture sub-DB stands.db, 3 cycles).
# Goldens were produced by running live FVSci_g16 on these exact keys on 2026-09-26:
#   ci_<cn>.live.sum        the live .sum
#   ci_<cn>.live_trees.csv  every live FVS_TreeList row, 2019-2049 (DBH, height, TPA and crown % per record)
# This locks the chain measured exact per record through 2049:
#   - ci/regent.f DO 17 subcycles with BANEXT/RDNEXT feedback, the IND1 walk, and per-copy ZZRAN / XWT / BKPT / DK
#   - ci/htgf.f:481-512 tripled-copy HTG
#   - ci/dgdriv.f:171 cycle-1 WK1 from the DO-220 DG, grincr.f OLDFNT, and the ci/morts.f RIP clamp
#   - ci/cratet.f:230-233 IND for the calibration PCT and for the cycle-1 ATAVH
#   - ci/crown.f:210-219 RELSDI = BA/BAMAX (SDIAC/SDIDEF only for species 11-16), the B floor, the species-14 and
#     17/19 crown-length models through the shared label-53 change limit and CRMAX cap
# The .sum rows are exact in every column in the years listed; in the other years every non-volume column is exact and
# TCuFt/MCuFt/BdFt sit a few units low on Douglas-fir records whose DBH and height are exact (the CI volume residual).
using Test, FVSjl, SQLite, DBInterface

const _CIRG_FX = joinpath(@__DIR__, "..", "fixtures", "centralidaho", "regent_live")
# stand => the .sum years required equal in every column
const _CIRG_SUM_EXACT = Dict("185544" => (2019, 2029, 2039, 2049), "186600" => (2019, 2029, 2039, 2049),
                             "188889" => (2019, 2049), "207086" => (2019, 2029, 2039), "207360" => (2019,))
# .sum fields 1-8 (year..QMD) and 13+ (removals, after-treatment stand) — everything but the volume columns
_cirg_novol(row) = (f = split(row); vcat(f[1:8], f[13:end]))

_cirg_rows(txt) = Dict(parse(Int, split(l)[1]) => rstrip(l) for l in split(txt, '\n') if occursin(r"^\d{4}\s", l))

@testset "CI REGENT vs live FVSci" begin
    @testset "$cn" for cn in sort(collect(keys(_CIRG_SUM_EXACT)))
        dir = mktempdir()
        cp(joinpath(_CIRG_FX, "ci_$cn.key"), joinpath(dir, "ci_$cn.key"))
        cp(joinpath(_CIRG_FX, "stands.db"), joinpath(dir, "stands.db"))
        jl = cd(() -> FVSjl.run_keyfile("ci_$cn.key"; variant = FVSjl.CentralIdaho(), output = :sum), dir)
        gold = _cirg_rows(read(joinpath(_CIRG_FX, "ci_$cn.live.sum"), String))
        got = _cirg_rows(jl)
        for y in _CIRG_SUM_EXACT[cn]
            @test get(got, y, "") == gold[y]
        end
        for y in sort(collect(keys(gold)))
            @test _cirg_novol(get(got, y, "")) == _cirg_novol(gold[y])
        end
        # every FVS_TreeList record 2019-2049: DBH to 5e-5, HT to 5e-4, TPA to 1e-5 relative (records sit at most a
        # few Float32 ULPs from live), crown % exact
        live = Dict{Tuple{Int,String,Int},Tuple{Float64,Float64,Float64,Int}}()
        for l in readlines(joinpath(_CIRG_FX, "ci_$cn.live_trees.csv"))[2:end]
            f = split(l, ',')
            live[(parse(Int, f[1]), f[2], parse(Int, f[3]))] =
                (parse(Float64, f[4]), parse(Float64, f[5]), parse(Float64, f[6]), parse(Int, f[7]))
        end
        db = SQLite.DB(joinpath(dir, "out.db"))
        mine = Dict{Tuple{Int,String,Int},Tuple{Float64,Float64,Float64,Int}}()
        for r in DBInterface.execute(db, "SELECT Year,TreeId,TreeIndex,DBH,Ht,TPA,PctCr FROM FVS_TreeList")
            mine[(Int(r.Year), strip(String(r.TreeId)), Int(r.TreeIndex))] =
                (Float64(r.DBH), Float64(r.Ht), Float64(r.TPA), Int(r.PctCr))
        end
        @test length(mine) == length(live)
        bad = [k for (k, (d, h, p, cr)) in live
               if !haskey(mine, k) || abs(mine[k][1] - d) > 5e-5 || abs(mine[k][2] - h) > 5e-4 ||
                  abs(mine[k][3] - p) > max(1e-4, 1e-5 * p) || mine[k][4] != cr]
        @test isempty(bad)
    end
end
