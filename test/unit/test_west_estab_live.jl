# test_west_estab_live.jl — western PLANT establishment vs LIVE oracles (tiered FIA stands, PLANT 2.0 DF 400 TPA).
#
# Keys: the tiered-west plant_cyc fixture keys of one stand per variant, plus TREELIST/TREELIDB so the planted cohort
# can be compared per record. Goldens: the live FVS<v>_g16 FVS_TreeList rows of the planted records (TreeIndex >= the
# first new record), every year (<case>_live_tl.csv). stands.db comes from test/fixtures/tiered/<v>/.

using Test, FVSjl, SQLite, DBInterface
const _WE = FVSjl
const _WE_FX = joinpath(@__DIR__, "..", "fixtures", "west_estab")
const _WE_TIER = joinpath(@__DIR__, "..", "fixtures", "tiered")
const _WE_VAR = Dict("so" => _WE.SouthCentralOregon(), "ws" => _WE.WestSierra(), "ca" => _WE.CentralCalifornia(),
                     "nc" => _WE.Klamath(), "kt" => _WE.Kootenai())

function _we_run(case)
    v = split(case, '_')[1]
    dir = mktempdir()
    cp(joinpath(_WE_TIER, v, "stands.db"), joinpath(dir, "stands.db"))
    key = replace(read(joinpath(_WE_FX, case * ".key"), String),
                  "\nstands.db\n" => "\n" * joinpath(dir, "stands.db") * "\n", "\nout.db\n" => "\n" * joinpath(dir, "out.db") * "\n")
    write(joinpath(dir, "k.key"), key)
    err = ""; txt = ""
    try
        txt = _WE.run_keyfile(joinpath(dir, "k.key"); variant = _WE_VAR[v])
    catch e
        err = sprint(showerror, e)
    end
    _WE_SUM[case] = txt
    return err, joinpath(dir, "out.db")
end
const _WE_SUM = Dict{String,String}()

function _we_live(case)
    rows = Dict{Tuple{Int,Int},Vector{Any}}()
    for l in readlines(joinpath(_WE_FX, case * "_live_tl.csv"))[2:end]
        f = split(l, ',')
        rows[(parse(Int, f[1]), parse(Int, f[2]))] = Any[f[3], parse(Float64, f[4]), parse(Float64, f[5]),
                                                         parse(Float64, f[6]), parse(Int, f[7])]
    end
    rows
end

function _we_jl(db, lo)
    rows = Dict{Tuple{Int,Int},Vector{Any}}()
    con = SQLite.DB(db)
    for r in DBInterface.execute(con, "select Year,TreeIndex,SpeciesFVS,TPA,DBH,Ht,PctCr from FVS_TreeList where TreeIndex>=$lo")
        rows[(Int(r.Year), Int(r.TreeIndex))] = Any[string(r.SpeciesFVS), Float64(r.TPA), Float64(r.DBH), Float64(r.Ht),
                                                    Int(r.PctCr)]
    end
    rows
end

const _WE_CASES = ["so_645126898126144_plant_cyc", "ws_15353585010497_plant_cyc", "ca_15320267010497_plant_cyc",
                   "nc_15303130010497_plant_cyc"]
# cases whose first-cycle planted heights/DBH are compared per record (NC: FINT=5 ⇒ LSKIPH, the listed height IS the
# ESTAB height)
const _WE_FIRST = ["nc_15303130010497_plant_cyc"]
const _WE_RUNS = Dict{String,Any}()

@testset "western PLANT establishment vs live: planted cohort per record" begin
    for case in _WE_CASES
        lv = _we_live(case)
        lo = minimum(k[2] for k in keys(lv))
        err, db = _we_run(case)
        _WE_RUNS[case] = (err, db)
        @test (case, err) == (case, "")
        isempty(err) || continue
        jl = _we_jl(db, lo)
        # the ESTAB booking: the planted cohort's records and their species/TPA at the end of the planting cycle
        yr1 = minimum(k[1] for k in keys(lv))
        k1 = sort([k for k in keys(lv) if k[1] == yr1])
        @test (case, sort([k for k in keys(jl) if k[1] == yr1])) == (case, k1)
        for k in k1
            haskey(jl, k) || continue
            @test (case, k, jl[k][1:2]) == (case, k, lv[k][1:2])
        end
    end
end

@testset "western PLANT establishment vs live: first-cycle height and DBH per record" begin
    for case in _WE_FIRST
        err, db = get(() -> _we_run(case), _WE_RUNS, case)
        isempty(err) || continue
        lv = _we_live(case)
        lo = minimum(k[2] for k in keys(lv)); yr1 = minimum(k[1] for k in keys(lv))
        jl = _we_jl(db, lo)
        for k in sort([k for k in keys(lv) if k[1] == yr1])
            haskey(jl, k) || continue
            @test (case, k, jl[k][3], jl[k][4]) == (case, k, lv[k][3], lv[k][4])
        end
    end
end

# ESGENT → REGENT(LESTB): the birth-cycle growth (crown dub from the post-growth PCCF, height/DBH over FINT−5 years,
# the WK4 step). At the end of the planting cycle every planted record equals live exactly. In the later cycles the
# cohort grows with each variant's ordinary growth and mortality, whose Float32 residuals in these stands' `none`
# regime (a few ULP in DBH/HT/TPA; TPA through the stand density) reach the cohort too: bounded at 1e-6 relative,
# crowns exact.
const _WE_FULL = ["so_645126898126144_plant_cyc", "ws_15353585010497_plant_cyc", "ca_15320267010497_plant_cyc",
                  "nc_15303130010497_plant_cyc"]
_we_close(a, b) = a[1] == b[1] && a[5] == b[5] && all(i -> abs(a[i] - b[i]) <= 1e-6 * abs(b[i]), 2:4)

@testset "western PLANT establishment vs live: the planted cohort (ESGENT)" begin
    for case in _WE_FULL
        err, db = get(() -> _we_run(case), _WE_RUNS, case)
        isempty(err) || continue
        lv = _we_live(case)
        lo = minimum(k[2] for k in keys(lv)); yr1 = minimum(k[1] for k in keys(lv))
        jl = _we_jl(db, lo)
        @test (case, sort(collect(keys(jl)))) == (case, sort(collect(keys(lv))))
        for k in sort(collect(keys(lv)))
            haskey(jl, k) || continue
            if k[1] == yr1
                @test (case, k, jl[k]) == (case, k, lv[k])
            else
                @test (case, k, _we_close(jl[k], lv[k])) == (case, k, true)
            end
        end
    end
end

# KT: estb/estab.f (== ie's) with the IE AUTOES tally on KT's 11 species, the estb PLANT path and estb/esgent.f →
# kt/regent.f REGENT(LESTB). Bare FIA stand 22404917010497 (no SLOPE/ASPECT ⇒ kt/grinit.f 30%/45°): natural
# regeneration from the first cycle (live 591 TPA at 2013; jl had none), and PLANT 400 DF on top of it.
const _WE_KT = ["kt_22404917010497_none", "kt_22404917010497_plant_cyc"]

@testset "KT establishment (AUTOES + PLANT + ESGENT) vs live" begin
    for case in _WE_KT
        err, db = _we_run(case)
        @test (case, err) == (case, "")
        isempty(err) || continue
        jl_rows = [l for l in split(_WE_SUM[case], '\n') if occursin(r"^\d{4} ", l)]
        lv_rows = [l for l in readlines(joinpath(_WE_FX, case * "_live.rows"))]
        @test (case, length(jl_rows)) == (case, length(lv_rows))
        for (a, b) in zip(jl_rows, lv_rows)
            @test (case, split(a)) == (case, split(b))
        end
        lv = _we_live(case); jl = _we_jl(db, 1)
        yr1 = 2013
        for k in sort(collect(keys(lv)))
            @test (case, k, haskey(jl, k)) == (case, k, true)
            haskey(jl, k) || continue
            if k[1] <= yr1
                @test (case, k, jl[k]) == (case, k, lv[k])      # the first natural cohort: every record exact
            else
                @test (case, k, jl[k][1] == lv[k][1] && jl[k][5] == lv[k][5] &&
                       all(i -> abs(jl[k][i] - lv[k][i]) <= 2e-6 * abs(lv[k][i]), 2:4)) == (case, k, true)
            end
        end
    end
end

# estb/estab.f's no-stocking branch (STOADJ≈0: NOAUTOES / NATURAL / STOCKADJ 0) still samples the per-plot site preps
# (:333-399, ESPREP defaults on NTALLY=1) and sets each plot's ESSUBH inputs before its PLANT/NATURAL trees. Stock
# ktt01 stand 5 (bare ground, NOAUTOES, PLANT 400 WL + 400 PP in 1992): live plots 1-41 NONE, 42-48 MECH, 49-50 BURN.
@testset "estb no-stocking branch: planted-tree site preps (KT ktt01 bare PLANT) vs live" begin
    case = "kt_ktt01_bareplant"
    err, db = _we_run(case)
    @test (case, err) == (case, "")
    if isempty(err)
        jl_rows = [l for l in split(_WE_SUM[case], '\n') if occursin(r"^\d{4} ", l)]
        lv_rows = readlines(joinpath(_WE_FX, case * "_live.rows"))
        @test length(jl_rows) == length(lv_rows)
        for (a, b) in zip(jl_rows, lv_rows)
            @test split(a) == split(b)
        end
        lv = _we_live(case); jl = _we_jl(db, 1)
        for k in sort(collect(keys(lv)))
            @test (k, get(jl, k, nothing)) == (k, lv[k])
        end
    end
end
