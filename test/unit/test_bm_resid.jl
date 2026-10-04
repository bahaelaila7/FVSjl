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

# Run one BM tiered fixture case; return (live golden rows, jl rows) of `table` as Year => Dict(column => Float32 or String),
# the jl side `nothing` when jl wrote no such table.
function _bm_case_table(cn::AbstractString, rg::AbstractString, table::AbstractString)
    fx = joinpath(@__DIR__, "..", "fixtures", "tiered", "bm")
    dir = mktempdir(); cp(joinpath(fx, "stands.db"), joinpath(dir, "stands.db"))
    key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
           for l in readlines(joinpath(fx, "$(cn)_$(rg).key"))]
    write(joinpath(dir, "s.key"), join(key, '\n'))
    FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = FVSjl.BlueMountains())
    val(x) = x isa AbstractString ? (y = tryparse(Float64, x); y === nothing ? String(x) : Float32(y)) :
             x isa Real ? Float32(x) : string(x)
    gold = Dict{Int,Dict{String,Any}}()
    gp = joinpath(fx, "$(cn)_$(rg).$(table).csv")
    if isfile(gp)
        gl = readlines(gp); hdr = split(gl[1], ',')
        for l in gl[2:end]
            g = Dict(String(h) => val(v) for (h, v) in zip(hdr, split(l, ',')))
            gold[Int(g["Year"])] = g
        end
    end
    db = SQLite.DB(joinpath(dir, "out.db"))
    has = !isempty(collect(DBInterface.execute(db, "SELECT name FROM sqlite_master WHERE type='table' AND name='$table'")))
    jl = has ? Dict{Int,Dict{String,Any}}() : nothing
    if has
        for r in DBInterface.execute(db, "SELECT * FROM $table")
            nt = NamedTuple(r)
            jl[Int(nt.Year)] = Dict(String(k) => val(v) for (k, v) in pairs(nt))
        end
    end
    SQLite.close(db)
    return gold, jl
end

@testset "BM WRD 177426703020004 rootdis FVS_RD_Sum vs live FVSbm_g16 (WK1, PRINF, NINSIM)" begin
    # bm/dgdriv.f:161/746-769 WK1 keeps the measured DG (Live_Merch_CuFt 2022 was 7.05 vs 14.77); rdcntl.f DO 800 PRINF
    # stored at RDMORT time (Ave_Pct_Root_Inf 2032 was 32.952 vs 33.00258); rdinsd.f:404 /(REAL(NINSIM)+1E-6).
    gold, jl = _bm_case_table("177426703020004", "rootdis", "FVS_RD_Sum")
    @test jl !== nothing
    for y in (2012, 2022), (c, v) in gold[y]
        c in ("StandID", "CaseID") && continue
        @test jl[y][c] == v
    end
    @test isapprox(jl[2032]["Ave_Pct_Root_Inf"], gold[2032]["Ave_Pct_Root_Inf"]; rtol = 1f-5)
    g2, j2 = _bm_case_table("22960873010497", "rootdis", "FVS_RD_Sum")
    # grincr.f:281-285 OLDTPA/ORMSQD at the cycle start feed RDROOT (was 0.487087 / 10.0826 vs 0.484182 / 10.154723)
    for y in (2017, 2027)
        @test isapprox(j2[y]["Ave_Pct_Root_Inf"], g2[y]["Ave_Pct_Root_Inf"]; rtol = 1f-5)
    end
    # rdpr.f:78 ITRN=0 ⇒ no report: a bare stand has no FVS_RD_Sum table in live
    _, j3 = _bm_case_table("722766017290487", "rootdis", "FVS_RD_Sum")
    @test j3 === nothing
end

@testset "BM FFE snag/crown pools 645155287126144 salvage FVS_Carbon vs live FVSbm_g16 (FMSADD slots, TFALL)" begin
    # fmsadd.f:47-62 empty height-class records shift the FMR6HTLS draws; bm/fmvinit.f TFALL (3/10/15/15 yr for DF).
    # Before: Standing_Dead 2048 0.2500 vs live 0.4512.
    gold, jl = _bm_case_table("645155287126144", "salvage", "FVS_Carbon")
    for y in sort(collect(keys(gold)))
        @test isapprox(jl[y]["Standing_Dead"], gold[y]["Standing_Dead"]; rtol = 1f-6)
        @test isapprox(jl[y]["Forest_Down_Dead_Wood"], gold[y]["Forest_Down_Dead_Wood"]; rtol = 1f-6)
    end
end

@testset "BM SIMFIRE 12827438010497 fire-year PotFire + post-fire carbon vs live FVSbm_g16" begin
    # fire-basis FMCFMD on pre-fire PROB/FMTBA, FMPOCR and FMEFF on the scorched FMICR (fmmain.f:188-196), and FMOLDC
    # recording the scorched crown (fmoldc.f:53). Before: Fuel_Wt 79/17 vs 56/42, Canopy_Density 0.016526 vs 0.017175,
    # Mortality_BA_Sev 22 vs 17, Aboveground_Total_Live 2025 22.6185 vs 22.6070, Forest_Down_Dead_Wood 2035 1.915 vs 1.832.
    gold, jl = _bm_case_table("12827438010497", "simfire", "FVS_PotFire")
    for c in ("Surf_Flame_Sev", "Surf_Flame_Mod", "Tot_Flame_Sev", "Tot_Flame_Mod", "Canopy_Density", "Crown_Index",
              "Torch_Index", "Mortality_BA_Sev", "Mortality_BA_Mod", "Mortality_VOL_Sev", "Mortality_VOL_Mod",
              "Fuel_Mod1", "Fuel_Mod2", "Fuel_Mod3", "Fuel_Wt1", "Fuel_Wt2", "Fuel_Wt3")
        @test jl[2015][c] == gold[2015][c]
    end
    gc, jc = _bm_case_table("12827438010497", "simfire", "FVS_Carbon")
    @test jc[2025]["Aboveground_Total_Live"] == gc[2025]["Aboveground_Total_Live"]
    @test isapprox(jc[2035]["Forest_Down_Dead_Wood"], gc[2035]["Forest_Down_Dead_Wood"]; rtol = 1f-6)
    # fmkill.f:129 + fmsadd.f:292-296: fire-cycle mortality crowns on UNFIRE = SNGNEW-FIRKIL (was 5.60227 vs 5.53913)
    g4, j4 = _bm_case_table("41136808010497", "simfire", "FVS_Carbon")
    @test isapprox(j4[2025]["Standing_Dead"], g4[2025]["Standing_Dead"]; rtol = 1f-6)
    # FMCBA's CRWDTH = the last CWIDTH's (gradd.f:254), not REGENT-grown seedlings: PERCOV ⇒ midflame wind exact (was
    # 2.1490381 vs 2.1491208). Flame/scorch carry the year-1 fuel split's ULPs (FMTBA still sums the seedlings' pre-REGENT
    # DBH, fmcba.f:246-270 at FMMAIN reads the grown ones) — open, so compared at 1E-5.
    g5, j5 = _bm_case_table("22960873010497", "simfire", "FVS_BurnReport")
    @test j5[2017]["Midflame_Wind"] == g5[2017]["Midflame_Wind"]
    for c in ("Flame_length", "Scorch_height")
        @test isapprox(j5[2017][c], g5[2017][c]; rtol = 1f-5)
    end
end
