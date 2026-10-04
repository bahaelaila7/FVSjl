# test_west_sprout_live.jl — ESUCKR stump/root sprouts of PN/WC/CA/WS (strp/essprt.f) and AK (estb/essprt.f) vs LIVE.
#
# These builds had no jl ESSPRT/NSPREC/SPRTHT tables, so a thinning/fire/ECON cut of a sprouting species fell into the
# SN form and crashed (KeyError :essprt_fsp — 40 tiered-west thinbba/simfire/econ cases). Keys: the tiered-west
# fixture keys with TREELIST/TREELIDB added; goldens: the live FVS<v>_g16 FVS_TreeList rows of every sprout record
# (TreeId ES…). The sprouts' first listing is compared exactly; later years to 1e-6 relative (the stands' ordinary
# growth ULPs), crowns exact.

using Test, FVSjl, SQLite, DBInterface
const _SP = FVSjl
const _SP_FX = joinpath(@__DIR__, "..", "fixtures", "west_sprout")
const _SP_VAR = Dict("wc" => _SP.WestCascades(), "pn" => _SP.PacificNorthwest(), "ak" => _SP.SoutheastAlaska(),
                     "ca" => _SP.CentralCalifornia(), "ws" => _SP.WestSierra())
const _SP_CASES = ["wc_33137520020004_thinbba", "pn_248879504489998_thinbba", "ak_644809321126144_thinbba",
                   "ca_23742358010900_econ", "ws_23762384010900_econ"]

function _sp_run(case)
    v = split(case, '_')[1]
    dir = mktempdir()
    cp(joinpath(@__DIR__, "..", "fixtures", "tiered", v, "stands.db"), joinpath(dir, "stands.db"))
    key = replace(read(joinpath(_SP_FX, case * ".key"), String),
                  "\nstands.db\n" => "\n" * joinpath(dir, "stands.db") * "\n", "\nout.db\n" => "\n" * joinpath(dir, "out.db") * "\n")
    write(joinpath(dir, "k.key"), key)
    err = ""
    try
        _SP.run_keyfile(joinpath(dir, "k.key"); variant = _SP_VAR[v])
    catch e
        err = sprint(showerror, e)
    end
    return err, joinpath(dir, "out.db")
end

@testset "western ESUCKR sprouts vs live (PN/WC/CA/WS/AK)" begin
    for case in _SP_CASES
        err, db = _sp_run(case)
        @test (case, err) == (case, "")
        isempty(err) || continue
        lv = Dict{Tuple{Int,String},Vector{Any}}()
        for l in readlines(joinpath(_SP_FX, case * "_live_es.csv"))[2:end]
            f = split(l, ',')
            lv[(parse(Int, f[1]), String(f[2]))] = Any[String(f[3]), parse(Float64, f[4]), parse(Float64, f[5]),
                                                        parse(Float64, f[6]), parse(Int, f[7])]
        end
        jl = Dict{Tuple{Int,String},Vector{Any}}()
        for r in DBInterface.execute(SQLite.DB(db), "select Year,TreeId,SpeciesFVS,TPA,DBH,Ht,PctCr from FVS_TreeList where TreeId like 'ES%'")
            jl[(Int(r.Year), String(r.TreeId))] = Any[string(r.SpeciesFVS), Float64(r.TPA), Float64(r.DBH), Float64(r.Ht), Int(r.PctCr)]
        end
        @test (case, sort(collect(keys(jl)))) == (case, sort(collect(keys(lv))))
        yr1 = minimum(k[1] for k in keys(lv))
        for k in sort(collect(keys(lv)))
            haskey(jl, k) || continue
            a = jl[k]; b = lv[k]
            if k[1] == yr1
                @test (case, k, a) == (case, k, b)
            else
                @test (case, k, a[1] == b[1] && a[5] == b[5] && all(i -> abs(a[i] - b[i]) <= 1e-6 * abs(b[i]), 2:4)) ==
                      (case, k, true)
            end
        end
    end
end
