# test_core_resid.jl — CORE (IE/EM/SN) tiered residual fixes, checked against the live-oracle goldens of the tiered
# fixtures (test/fixtures/tiered/<v>, built from FVS{v}_g16). Float32 bits compared (the golden CSVs hold FVS's REAL*4
# as text) unless a test states an upstream ULP residual that is still open.
using FVSjl, Test, SQLite, DBInterface

const _CR_VAR = Dict("ie" => FVSjl.InlandEmpire(), "em" => FVSjl.EasternMontana(), "sn" => FVSjl.Southern())

"Run tiered fixture `<v>/<cn>_<rg>.key` through run_keyfile; returns the output DB path."
function _cr_run(v::AbstractString, cn::AbstractString, rg::AbstractString)
    fx = joinpath(@__DIR__, "..", "fixtures", "tiered", v)
    dir = mktempdir(); cp(joinpath(fx, "stands.db"), joinpath(dir, "stands.db"))
    key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
           for l in readlines(joinpath(fx, "$(cn)_$(rg).key"))]
    write(joinpath(dir, "s.key"), join(key, '\n'))
    FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = _CR_VAR[v])
    return joinpath(dir, "out.db")
end

"Golden rows of table `t` (Vector of Dict col => text) and jl rows (Vector of Dict col => value)."
function _cr_table(v, cn, rg, db, t)
    fx = joinpath(@__DIR__, "..", "fixtures", "tiered", v)
    gl = readlines(joinpath(fx, "$(cn)_$(rg).$(t).csv")); hdr = split(gl[1], ',')
    gold = [Dict(zip(hdr, split(l, ','))) for l in gl[2:end]]
    d = SQLite.DB(db)
    jl = [Dict(String(k) => r[k] for k in propertynames(r)) for r in DBInterface.execute(d, "SELECT * FROM $t")]
    SQLite.close(d)
    return gold, jl
end
_cr_f32(x) = x isa AbstractString ? Float32(parse(Float64, x)) : Float32(x)

@testset "EM 684750664126144 cycle-0 dead aspen R1KEMP KLASS=1 + DVE CFTOPK (FVS_TreeList) vs live FVSem_g16" begin
    # em/vols.f:136-140 LIVEDEAD='D' ⇒ r1kemp.f KLASS=1 (cubic min 1.6, not 2.4) and NATCRS CTKFLG=.TRUE. for DVE
    # (fvsvol.f:531) ⇒ the broken-top trim (vols.f:194-196). Dead AS D6.1 TruncHt 26: live TCuFt 2.2532847; before 2.4.
    db = _cr_run("em", "684750664126144", "none")
    gold, jl = _cr_table("em", "684750664126144", "none", db, "FVS_TreeList")
    num = ["TPA", "MortPA", "DBH", "DG", "Ht", "HtG", "TCuFt", "MCuFt", "BdFt"]
    key(r) = (parse(Int, string(r["Year"])), parse(Int, string(r["TreeIndex"])))
    jd = Dict(key(r) => r for r in jl)
    @test length(jl) == length(gold)
    bad = [key(g) for g in gold if !haskey(jd, key(g)) || any(_cr_f32(jd[key(g)][c]) != _cr_f32(g[c]) for c in num)]
    @test isempty(bad)
    g = only(x for x in gold if x["Year"] == "2018" && x["TreeIndex"] == "2998")
    @test _cr_f32(jd[(2018, 2998)]["TCuFt"]) == _cr_f32(g["TCuFt"]) == 2.2532847f0
end

@testset "EM 684750664126144 salvage FVS_Carbon Standing_Dead (FMSVOL DVE CFTOPK) vs live FVSem_g16" begin
    # em/fmsvol.f:139-140 trims every equation family (CTKFLG from NATCRS): the inventory AS snag VOL2HT 2.3469481,
    # not the untrimmed 2.4 (Standing_Dead 2018 0.7333492 live vs 0.7350840 before).
    db = _cr_run("em", "684750664126144", "salvage")
    gold, jl = _cr_table("em", "684750664126144", "salvage", db, "FVS_Carbon")
    jd = Dict(parse(Int, string(r["Year"])) => r for r in jl)
    for g in gold
        y = parse(Int, g["Year"])
        @test _cr_f32(jd[y]["Standing_Dead"]) == _cr_f32(g["Standing_Dead"])
    end
end

@testset "IE 11855985010690 salvage FVS_PotFire Pot_Smoke_Sev — FMEFF CWD2B head on the FMPOFL-year pools" begin
    # fmmain.f:196 FMPOFL → FMEFF (ICALL=1) burns CRBURN of the waiting snag crowns CWD2B/CWD2B2 (fmeff.f:118-138)
    # BEFORE the year's FMCADD drops them (fmmain.f:241). Before: 2016 0.3109127 vs live 0.3129615 (head 0.0017 vs
    # 0.1941 t/ac). The 2006 cycle-0 Canopy_Density is 1 ULP off (open, upstream of this), so the later years are
    # compared at 2e-6 relative.
    db = _cr_run("ie", "11855985010690", "salvage")
    gold, jl = _cr_table("ie", "11855985010690", "salvage", db, "FVS_PotFire")
    jd = Dict(parse(Int, string(r["Year"])) => r for r in jl)
    for g in gold, c in ("Pot_Smoke_Sev", "Pot_Smoke_Mod")
        y = parse(Int, g["Year"])
        @test isapprox(_cr_f32(jd[y][c]), _cr_f32(g[c]); rtol = 2f-6)
    end
end

@testset "IE 11855985010690 simfire 2016 burn (FMCBA CRWDTH = last CWIDTH dims) vs live FVSie_g16" begin
    # ie/fmcba.f:244 CWIDTH=CRWDTH(I): the end-of-previous-cycle crown width TRIPLE copied, not one recomputed from the
    # seam DBH of this cycle's grown small trees. Before: TOTCRA 30163.06 vs 30151.213 ⇒ Midflame_Wind 1.66760 vs
    # 1.66796 and every FVS_Mortality fire-kill cell of the 2016 burn off (Bakill 132.00223 vs 132.00040).
    db = _cr_run("ie", "11855985010690", "simfire")
    for t in ("FVS_BurnReport", "FVS_Mortality")
        gold, jl = _cr_table("ie", "11855985010690", "simfire", db, t)
        @test length(jl) == length(gold)
        k(r) = (string(r["Year"]), string(get(r, "SpeciesFVS", "")))
        jd = Dict(k(r) => r for r in jl)
        cols = [c for c in keys(gold[1]) if !(c in ("StandID", "Year", "SpeciesFIA")) && tryparse(Float64, gold[1][c]) !== nothing]
        for g in gold, c in cols
            @test _cr_f32(jd[k(g)][c]) == _cr_f32(g[c])
        end
    end
end

@testset "EM 196378260020004 mistletoe stand: MISTOE (gradd.f:96) before FMMAIN (gradd.f:118) vs live FVSem_g16" begin
    # Non-fire tripling cycle: FMPTRH must draw after the post-TRIPLE mistletoe spread (live RANNGET 2036729867 at 2012;
    # jl sampled at 431495394 ⇒ PTorch_Mod 0.13873 vs 0.25494). SIMFIRE tripling cycle: the spread precedes the burn
    # (jl burned first ⇒ 2022 PTorch_Sev 0.23130 vs 0.18498, BA 45 vs 46 at 2042).
    db = _cr_run("em", "196378260020004", "salvage")
    gold, jl = _cr_table("em", "196378260020004", "salvage", db, "FVS_PotFire")
    jd = Dict(string(r["Year"]) => r for r in jl)
    for g in gold, c in ("PTorch_Sev", "PTorch_Mod", "Torch_Index", "Crown_Index", "Mortality_BA_Sev")
        @test _cr_f32(jd[g["Year"]][c]) == _cr_f32(g[c])
    end
    db = _cr_run("em", "196378260020004", "simfire")
    gold, jl = _cr_table("em", "196378260020004", "simfire", db, "FVS_Summary")
    jd = Dict(string(r["Year"]) => r for r in jl)
    for g in gold, c in ("Tpa", "BA", "SDI", "CCF", "TopHt", "QMD", "TCuFt", "MCuFt", "BdFt")
        @test _cr_f32(jd[g["Year"]][c]) == _cr_f32(g[c])
    end
    gold, jl = _cr_table("em", "196378260020004", "simfire", db, "FVS_PotFire")
    jd = Dict(string(r["Year"]) => r for r in jl)
    for g in gold
        g["Year"] in ("2012", "2022") || continue
        @test _cr_f32(jd[g["Year"]]["PTorch_Sev"]) == _cr_f32(g["PTorch_Sev"])
        @test _cr_f32(jd[g["Year"]]["PTorch_Mod"]) == _cr_f32(g["PTorch_Mod"])
    end
end
