# CI REGENT growth vs the live FVSci_g16 oracle on five FIA Central Idaho stands (fixture sub-DB stands.db, 2 cycles).
# Goldens were produced by running live FVSci_g16 on these exact keys on 2026-09-26:
#   ci_<cn>.live.sum        the live .sum
#   ci_<cn>.live_trees.csv  the live FVS_TreeList rows for 2029 (the end of cycle 1: each record's DBH, height and TPA
#                           after REGENT, the large-tree growth, tripling and morts)
# This locks the chain measured exact per record at cycle 1:
#   - ci/regent.f DO 17 subcycles with BANEXT/RDNEXT feedback, the IND1 walk, and per-copy ZZRAN / XWT / BKPT / DK
#   - ci/htgf.f:481-512 tripled-copy HTG
#   - ci/dgdriv.f:171 cycle-1 WK1 from the DO-220 DG, grincr.f OLDFNT, and the ci/morts.f RIP clamp
#   - ci/cratet.f:230-233 IND for the calibration PCT and for the cycle-1 ATAVH
# The .sum rows are exact through 2029, except where marked. TCuFt/MCuFt on 188889 and 207360 at 2029 are one
# cubic foot below live on records whose DBH and height are exact: the FW2/MATW volume residual, not REGENT.
using Test, FVSjl, SQLite, DBInterface

const _CIRG_FX = joinpath(@__DIR__, "..", "fixtures", "centralidaho", "regent_live")
# stand => the .sum years required equal in every column
const _CIRG_SUM_EXACT = Dict("185544" => (2019, 2029), "186600" => (2019, 2029), "188889" => (2019,),
                             "207086" => (2019, 2029), "207360" => (2019,))

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
        # every FVS_TreeList record at the end of cycle 1 (2029): DBH to 5e-5, HT to 5e-4, TPA to 1e-5 relative (185544's
        # 362-TPA seedling copies sit 8 ULP from live: Float32 exp/log ULPs, not a logic difference)
        live = Dict{Tuple{String,Int},NTuple{3,Float64}}()
        for l in readlines(joinpath(_CIRG_FX, "ci_$cn.live_trees.csv"))[2:end]
            f = split(l, ',')
            live[(f[2], parse(Int, f[3]))] = (parse(Float64, f[4]), parse(Float64, f[5]), parse(Float64, f[6]))
        end
        db = SQLite.DB(joinpath(dir, "out.db"))
        mine = Dict{Tuple{String,Int},NTuple{3,Float64}}()
        for r in DBInterface.execute(db, "SELECT TreeId,TreeIndex,DBH,Ht,TPA FROM FVS_TreeList WHERE Year=2029")
            mine[(strip(String(r.TreeId)), Int(r.TreeIndex))] = (Float64(r.DBH), Float64(r.Ht), Float64(r.TPA))
        end
        @test length(mine) == length(live)
        bad = [k for (k, (d, h, p)) in live
               if !haskey(mine, k) || abs(mine[k][1] - d) > 5e-5 || abs(mine[k][2] - h) > 5e-4 ||
                  abs(mine[k][3] - p) > max(1e-4, 1e-5 * p)]
        @test isempty(bad)
    end
end
