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
                     "nc" => _WE.Klamath())

function _we_run(case)
    v = split(case, '_')[1]
    dir = mktempdir()
    cp(joinpath(_WE_TIER, v, "stands.db"), joinpath(dir, "stands.db"))
    key = replace(read(joinpath(_WE_FX, case * ".key"), String),
                  "\nstands.db\n" => "\n" * joinpath(dir, "stands.db") * "\n", "\nout.db\n" => "\n" * joinpath(dir, "out.db") * "\n")
    write(joinpath(dir, "k.key"), key)
    err = ""
    try
        _WE.run_keyfile(joinpath(dir, "k.key"); variant = _WE_VAR[v])
    catch e
        err = sprint(showerror, e)
    end
    return err, joinpath(dir, "out.db")
end

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

const _WE_CASES = ["so_645126898126144_plant_cyc", "ws_15353585010497_plant_cyc", "ca_15320267010497_plant_cyc"]

@testset "western PLANT establishment vs live: planted cohort per record" begin
    for case in _WE_CASES
        lv = _we_live(case)
        lo = minimum(k[2] for k in keys(lv))
        err, db = _we_run(case)
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
