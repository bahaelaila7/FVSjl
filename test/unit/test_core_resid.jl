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

@testset "EM 231908428020004 MISTPINF on an absent host: MISPRT DMFLAG from MISINF vs live FVSem_g16" begin
    # misinf.f:182 sets DMFLAG for any host species a MISTPINF card targets, even with no trees of it; MISPRT
    # (fvs.f:400) then writes one zero-infection FVS_DM_Stnd_Sum row and creates an empty FVS_DM_Spp_Sum; the next
    # MISTOE (mistoe.f:193) clears it. jl wrote neither table.
    db = _cr_run("em", "231908428020004", "mistletoe")
    gold, jl = _cr_table("em", "231908428020004", "mistletoe", db, "FVS_DM_Stnd_Sum")
    @test length(jl) == length(gold) == 1
    for c in keys(gold[1])
        c == "StandID" && continue
        @test _cr_f32(jl[1][c]) == _cr_f32(gold[1][c])
    end
    gold, jl = _cr_table("em", "231908428020004", "mistletoe", db, "FVS_DM_Spp_Sum")
    @test isempty(gold) && isempty(jl)
end

@testset "EM 231908428020004 econ: DBSECHARV_open creates FVS_EconHarvestValue on a non-PCT harvest vs live FVSem_g16" begin
    # eccalc.f:745-746 opens (creates) the table for every non-PCT harvest with IDBSECON=2, before any HRVRVN-valued
    # row; this stand's harvest has none, so live carries the empty table (jl omitted it).
    db = _cr_run("em", "231908428020004", "econ")
    d = SQLite.DB(db)
    tabs = [r[:name] for r in DBInterface.execute(d, "SELECT name FROM sqlite_master WHERE type='table'")]
    n = "FVS_EconHarvestValue" in tabs ? only(r[1] for r in DBInterface.execute(d, "SELECT COUNT(*) FROM FVS_EconHarvestValue")) : -1
    SQLite.close(d)
    @test n == 0
end

@testset "IE 3356357010690 / EM 196420598020004 thinbba FVS_StrClass after-thin cover: SSTAGE reads the stored CRWDTH" begin
    # sstage.f:238/276 WK6=CRWDTH(I), filled by CWIDTH at load / gradd.f:254 — the pre-thin stand BA in the Crookston
    # BAREA term, not the residual BA (after-thin Total_Cover 17 vs live 16 before; COVER 16.6095 vs 16.2964).
    for (v, cn) in (("ie", "3356357010690"), ("em", "196420598020004"))
        db = _cr_run(v, cn, "thinbba")
        gold, jl = _cr_table(v, cn, "thinbba", db, "FVS_StrClass")
        k(r) = (string(r["Year"]), string(r["Removal_Code"]))
        jd = Dict(k(r) => r for r in jl)
        @test length(jl) == length(gold)
        for g in gold, c in ("Stratum_1_Crown_Cover", "Stratum_2_Crown_Cover", "Total_Cover", "Structure_Class")
            @test string(jd[k(g)][c]) == g[c] || _cr_f32(jd[k(g)][c]) == _cr_f32(g[c])
        end
    end
end

@testset "IE 1627682513290487 thinbba: no COR/HCOR attenuation for a species with no records (dgdriv.f IF(I1.EQ.0) GO TO 50)" begin
    # The 2031 thin removed every WH; dgdriv.f skips the attenuation for WH that cycle, so ESTAB's REGENT books the new
    # WH regeneration with the cycle-1 HCOR 0.00339259 (jl attenuated to 0.0059646 ⇒ HtG +0.26%, QMD off from 2051).
    db = _cr_run("ie", "1627682513290487", "thinbba")
    gold, jl = _cr_table("ie", "1627682513290487", "thinbba", db, "FVS_Summary")
    jd = Dict(string(r["Year"]) => r for r in jl)
    for g in gold, c in ("Tpa", "BA", "QMD", "ATQMD", "TCuFt", "TopHt")
        @test _cr_f32(jd[g["Year"]][c]) == _cr_f32(g[c])
    end
    gold, jl = _cr_table("ie", "1627682513290487", "thinbba", db, "FVS_StrClass")
    k(r) = (string(r["Year"]), string(r["Removal_Code"]))
    jd2 = Dict(k(r) => r for r in jl)
    for g in gold, c in ("Stratum_1_DBH", "Stratum_2_DBH")
        @test _cr_f32(jd2[k(g)][c]) == _cr_f32(g[c])
    end
end

@testset "IE 3285544010690 thinbba: REGENT(LESTB) seedling crown on the pre-ESNUTR PCCF vs live FVSie_g16" begin
    # regent.f:301-319 CR=0.89722-0.0000461*PCCF reads the gradd.f:192 DENSE PCCF (before ESNUTR); jl's density had
    # the AUTOES cohort in it (point 2 178.16805 vs 177.69788) ⇒ ICR 89 vs 90 ⇒ QMD off at 2052.
    db = _cr_run("ie", "3285544010690", "thinbba")
    gold, jl = _cr_table("ie", "3285544010690", "thinbba", db, "FVS_Summary")
    jd = Dict(string(r["Year"]) => r for r in jl)
    for g in gold, c in ("Tpa", "BA", "QMD", "ATQMD", "TCuFt", "TopHt")
        @test _cr_f32(jd[g["Year"]][c]) == _cr_f32(g[c])
    end
end

@testset "IE 374547584489998 (KODFOR 621) FFE live merch stem on the stand's NVEL region vs live FVSie_g16" begin
    # FMSVL2 → NATCRS uses IREGN=KODFOR/100 (fvsvol.f:90-96): the Colville's region-6 merch rules, as VOLS does.
    # Before: Aboveground_Merch_Live 39.71698 vs live 39.78778 at 2015 (region-1 rules).
    db = _cr_run("ie", "374547584489998", "salvage")
    gold, jl = _cr_table("ie", "374547584489998", "salvage", db, "FVS_Carbon")
    jd = Dict(string(r["Year"]) => r for r in jl)
    for g in gold
        @test _cr_f32(jd[g["Year"]]["Aboveground_Merch_Live"]) == _cr_f32(g["Aboveground_Merch_Live"])
    end
end

@testset "EM 3087467010690 FFE live merch stem (FMSVL2 NATCRS SF_HS merch top) vs live FVSem_g16" begin
    # em_nocut_cuft ran the FW2 profile with the bisection merch top; NATCRS (as in VOLS) solves it with SF_HS
    # (PP D9.2199 H39.4966: MCF 5.5 live vs 5.8 ⇒ 2008 Aboveground_Merch_Live 11.94416 vs 11.94876).
    for rg in ("salvage", "simfire")
        db = _cr_run("em", "3087467010690", rg)
        gold, jl = _cr_table("em", "3087467010690", rg, db, "FVS_Carbon")
        jd = Dict(string(r["Year"]) => r for r in jl)
        for g in gold
            @test _cr_f32(jd[g["Year"]]["Aboveground_Merch_Live"]) == _cr_f32(g["Aboveground_Merch_Live"])
        end
    end
end

@testset "EM 196378260020004 rootdis RD_Sum Live_Merch_CuFt: WK1=DG for last cycle's regeneration too (dgdriv.f:144)" begin
    # RDPR CFVPA sums TCLAS·WK1(I); DGDRIV sets WK1(I)=DG(I) for every record, and a regenerated record's DG is its
    # birth DK (regent.f:941). jl zeroed those (record 393: 0 vs 0.5952508 ⇒ 2042 130.13326 vs 131.93132).
    db = _cr_run("em", "196378260020004", "rootdis")
    gold, jl = _cr_table("em", "196378260020004", "rootdis", db, "FVS_RD_Sum")
    jd = Dict(string(r["Year"]) => r for r in jl)
    for g in gold
        g["Year"] in ("2022", "2032", "2042") || continue
        @test _cr_f32(jd[g["Year"]]["Live_Merch_CuFt"]) == _cr_f32(g["Live_Merch_CuFt"])
    end
end

@testset "IE 11855985010690 salvage FVS_Carbon: the non-fire FMMAIN pass runs after REGENT's direct small-tree DBH (gradd.f:118)" begin
    # FMCBA's TBA (fmcba.f:236-237) reads the seedlings' REGENT DBH (ie/regent.f:881) because FMMAIN follows GRINCR: sp-9 TBA
    # 21.4268 live vs 21.4036 from the start-of-cycle DBH ⇒ PRCL ⇒ the initial dead fuels 1 ULP (2006 Forest_Down_Dead_Wood
    # 14.608339 vs live 14.608341), drifting every later fuel/carbon row.
    db = _cr_run("ie", "11855985010690", "salvage")
    gold, jl = _cr_table("ie", "11855985010690", "salvage", db, "FVS_Carbon")
    jd = Dict(parse(Int, string(r["Year"])) => r for r in jl)
    @test length(jl) == length(gold)
    cols = ["Aboveground_Total_Live", "Aboveground_Merch_Live", "Belowground_Live", "Belowground_Dead", "Standing_Dead",
            "Forest_Down_Dead_Wood", "Forest_Floor", "Forest_Shrub_Herb", "Total_Stand_Carbon"]
    for g in gold, c in cols
        @test _cr_f32(jd[parse(Int, g["Year"])][c]) == _cr_f32(g[c])
    end
end
