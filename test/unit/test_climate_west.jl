# test_climate_west.jl — Climate-FVS in the western variants that used to run CLIMATE as a no-op (#252).
#
# Every fixture carries the EM tiered CGCM3_A2 ClimData rows re-keyed to the stand (plus declining synthetic
# viability columns for the variant's PLNJSP symbols the header lacks), GROWMULT 1.1 and CLIMREDB. Goldens are the
# live oracle's rows on these exact files (2026-09-26): FVSut_g16, FVSec_g16, FVSkt_clean, FVSop_clean, FVSak_g16.
#
# - UT 39467329010690 and EC 212810748020004: the whole .sum AND every FVS_Climate row equal live. The EC stand
#   guards dgdriv.f:217 DDS=EXP(WK2+XDGROW)·WK4 with WK2 left un-multiplied: DGSCOR damps FRM on WK2>4, so folding
#   ln(WK4) into WK2 flattened its 40"+ DF growth from cycle 3 on.
# - KT ktt01 and OP opt01: the climate columns (Viability/ViabMort/dClimMort/GrowthMult/SiteMult/MxDenMult) equal
#   live in every row; BA/TPA/AutoEstb carry each variant's own no-climate residual (KT cycle-0 volume, OP DF
#   crown change). OP's GrowthMult was 1.0 before its ORGANON DDS took CLGMULT's WK4.
# - AK links exclim.f stubs: CLIMATE raises FVS11 and changes nothing (live .sum identical with and without it,
#   no FVS_Climate table), so jl must not activate it either.
using Test
using FVSjl
using SQLite, DBInterface

const _CW_FX = joinpath(@__DIR__, "..", "fixtures", "climate_west")
const _CW_COLS = ["Viability", "BA", "TPA", "ViabMort", "dClimMort", "GrowthMult", "SiteMult", "MxDenMult",
                  "AutoEstbTPA"]

_cw_sumrows(txt) = [split(strip(l)) for l in split(txt, '\n') if occursin(r"^\d{4}\s", strip(l))]

function _cw_run(name, variant)
    dir = mktempdir()
    for f in readdir(_CW_FX)
        (startswith(f, name) || f == "stands.db") || continue
        cp(joinpath(_CW_FX, f), joinpath(dir, replace(f, name => "k")))
    end
    txt = cd(() -> FVSjl.run_keyfile("k.key"; variant = variant, output = :sum), dir)
    return dir, txt
end

function _cw_climate(db)
    out = Dict{Tuple{Int,String},Vector{Float64}}()
    "FVS_Climate" in [t.name for t in SQLite.tables(db)] || return out
    for r in DBInterface.execute(db, "SELECT Year,SpeciesFVS," * join(_CW_COLS, ",") * " FROM FVS_Climate")
        out[(Int(r.Year), strip(String(r.SpeciesFVS)))] = [Float64(coalesce(r[Symbol(c)], 0.0)) for c in _CW_COLS]
    end
    return out
end

function _cw_golden(name)
    lines = readlines(joinpath(_CW_FX, name * ".FVS_Climate.csv"))
    hdr = split(lines[1], ',')
    g = Dict{Tuple{Int,String},Vector{Float64}}()
    for l in lines[2:end]
        f = split(l, ',')
        g[(parse(Int, f[1]), String(f[2]))] = [parse(Float64, f[findfirst(==(c), hdr)]) for c in _CW_COLS]
    end
    return g
end

_cw_close(a, b) = abs(a - b) <= max(1e-3, 1e-4 * abs(a))

@testset "Climate-FVS western wiring vs live (#252)" begin
    for (name, variant, cols) in (("ut_39467329010690", FVSjl.Utah(), 1:9),
                                  ("ec_212810748020004", FVSjl.EastCascades(), 1:9),
                                  ("kt_ktt01", FVSjl.Kootenai(), (1, 4, 5, 6, 7, 8)),
                                  ("op_opt01", FVSjl.Olympic(), (1, 4, 5, 6, 7, 8)))
        @testset "$name" begin
            dir, txt = _cw_run(name, variant)
            gold = _cw_golden(name)
            got = _cw_climate(SQLite.DB(joinpath(dir, "out.db")))
            @test length(got) == length(gold)
            nbad = count(k -> !haskey(got, k) || !all(_cw_close(gold[k][i], got[k][i]) for i in cols), keys(gold))
            @test nbad == 0
            if cols == 1:9
                live = [split(strip(l)) for l in readlines(joinpath(_CW_FX, name * ".live.sum"))]
                @test _cw_sumrows(txt) == live
            end
        end
    end
    @testset "AK: exclim stubs ⇒ CLIMATE inert" begin
        dir, with = _cw_run("ak_666708911126144", FVSjl.SoutheastAlaska())
        _, without = _cw_run("ak_666708911126144_noclim", FVSjl.SoutheastAlaska())
        @test _cw_sumrows(with) == _cw_sumrows(without)
        @test isempty(_cw_climate(SQLite.DB(joinpath(dir, "out.db"))))
    end
end
