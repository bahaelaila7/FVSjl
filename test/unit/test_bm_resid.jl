# test_bm_resid.jl — BM tiered residual fixes, checked record for record against the live FVSbm_g16 goldens of the
# tiered fixture (test/fixtures/tiered/bm). Float32 bits compared (the golden CSV holds FVS's REAL*4 as text).
using FVSjl, Test, SQLite, DBInterface

function _bm_fixture_treelist(cn::AbstractString, rg::AbstractString)
    fx = joinpath(@__DIR__, "..", "fixtures", "tiered", "bm")
    dir = mktempdir(); cp(joinpath(fx, "stands.db"), joinpath(dir, "stands.db"))
    key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
           for l in readlines(joinpath(fx, "$(cn)_$(rg).key"))]
    write(joinpath(dir, "s.key"), join(key, '\n'))
    FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = FVSjl.BlueMountains())
    gl = readlines(joinpath(fx, "$(cn)_$(rg).FVS_TreeList.csv")); hdr = split(gl[1], ',')
    num = [c for c in hdr if !(c in ("StandID", "TreeId", "SpeciesFVS", "SpeciesPLANTS", "SpeciesFIA"))]
    f32(x) = x isa AbstractString ? Float32(parse(Float64, x)) : Float32(x)
    gold = Dict{Tuple{Int,Int},Vector{Float32}}()
    for l in gl[2:end]
        f = split(l, ','); g = Dict(zip(hdr, f))
        gold[(parse(Int, g["Year"]), parse(Int, g["TreeIndex"]))] = [f32(g[c]) for c in num]
    end
    db = SQLite.DB(joinpath(dir, "out.db"))
    jl = Dict{Tuple{Int,Int},Vector{Float32}}()
    for r in DBInterface.execute(db, "SELECT " * join(num, ",") * " FROM FVS_TreeList")
        jl[(Int(r[:Year]), Int(r[:TreeIndex]))] = [f32(r[Symbol(c)]) for c in num]
    end
    SQLite.close(db)
    return gold, jl
end

@testset "BM tiered 22960873010497 (dense LP/GF/AF/ES regeneration) FVS_TreeList vs live FVSbm_g16" begin
    # 2,443 records over 2007-2057: HTGF/REGENT on glibc expf/logf/powf (bm/htgf.f, regent.f, htdbh.f, findag.f) and
    # RELDEN summed species-major over IND1 (dense.f:95-140). Before: 2027 HtG/BAPctile 1-6 ULP off, every record
    # off from 2037.
    gold, jl = _bm_fixture_treelist("22960873010497", "none")
    @test length(gold) == 2443
    @test length(jl) == length(gold)
    bad = [k for k in keys(gold) if get(jl, k, Float32[]) != gold[k]]
    @test isempty(bad)
end

@testset "BM tiered 645155287126144 simfire FVS_Carbon 2018 + 2028 (fire year) vs live FVSbm_g16" begin
    # BM snags on the western FMSVOL layer (fmsvol.f CFTOPK at HTIH, fmcwd.f CWD2 broken tops): before, the fire-year
    # Standing_Dead was 7.02 vs live 1.99 and Forest_Down_Dead_Wood 3.19 vs 5.95 (snag densities already equal).
    fx = joinpath(@__DIR__, "..", "fixtures", "tiered", "bm")
    dir = mktempdir(); cp(joinpath(fx, "stands.db"), joinpath(dir, "stands.db"))
    key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
           for l in readlines(joinpath(fx, "645155287126144_simfire.key"))]
    write(joinpath(dir, "s.key"), join(key, '\n'))
    FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = FVSjl.BlueMountains())
    gl = readlines(joinpath(fx, "645155287126144_simfire.FVS_Carbon.csv")); hdr = split(gl[1], ',')
    num = [c for c in hdr if !(c in ("StandID",))]
    db = SQLite.DB(joinpath(dir, "out.db"))
    jl = Dict(Int(r[:Year]) => [Float32(r[Symbol(c)]) for c in num]
              for r in DBInterface.execute(db, "SELECT " * join(num, ",") * " FROM FVS_Carbon"))
    SQLite.close(db)
    for l in gl[2:end]
        g = Dict(zip(hdr, split(l, ','))); y = parse(Int, g["Year"])
        y in (2018, 2028) || continue
        gv = [Float32(parse(Float64, g[c])) for c in num]
        for (k, c) in enumerate(num)
            if y == 2028 && c == "Standing_Dead"
                @test isapprox(jl[y][k], gv[k]; rtol = 1f-6)   # 1 ULP (1.9869455 live, the fmdout.f sum order) — open
            else
                @test jl[y][k] == gv[k]
            end
        end
    end
end

@testset "BM tiered 374443645489998 plant_cal FVS_Summary QMD (ESGENT HTG·WK4) vs live FVSbm_g16" begin
    # bm/esgent.f HTG=HTG*WK4 (WK4=HTIMLT 0.99998) and the WK4<1 DBH/DG rescale; before, the planted cohort's QMD was
    # 7e-6 high from 2035 on (1.8450947 vs 1.8450816).
    fx = joinpath(@__DIR__, "..", "fixtures", "tiered", "bm")
    dir = mktempdir(); cp(joinpath(fx, "stands.db"), joinpath(dir, "stands.db"))
    key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
           for l in readlines(joinpath(fx, "374443645489998_plant_cal.key"))]
    write(joinpath(dir, "s.key"), join(key, '\n'))
    FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = FVSjl.BlueMountains())
    gl = readlines(joinpath(fx, "374443645489998_plant_cal.FVS_Summary.csv")); hdr = split(gl[1], ',')
    iy = findfirst(==("Year"), hdr); iq = findfirst(==("QMD"), hdr)
    gold = Dict(parse(Int, split(l, ',')[iy]) => Float32(parse(Float64, split(l, ',')[iq])) for l in gl[2:end])
    db = SQLite.DB(joinpath(dir, "out.db"))
    jl = Dict(Int(r[:Year]) => Float32(r[:QMD]) for r in DBInterface.execute(db, "SELECT Year, QMD FROM FVS_Summary"))
    SQLite.close(db)
    @test length(gold) == 6
    @test all(jl[y] == q for (y, q) in gold)
end
