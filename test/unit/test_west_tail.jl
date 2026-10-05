# test_west_tail.jl — western-tail / shared FFE fixes, each checked against live-oracle goldens in test/fixtures/west_tail
# (keys run on the tiered fixtures' stands.db; goldens dumped from FVS{v}_g16 out.db, CaseID dropped). Float32 bits.
using FVSjl, Test, SQLite, DBInterface

const _WT_FX = joinpath(@__DIR__, "..", "fixtures", "west_tail")
const _WT_VAR = Dict("ws" => FVSjl.WestSierra(), "sn" => FVSjl.Southern())

"Run west_tail/<name>.key for variant `v` (stands.db from the tiered fixture); returns the output DB path."
function _wt_run(v::AbstractString, name::AbstractString)
    dir = mktempdir(); cp(joinpath(@__DIR__, "..", "fixtures", "tiered", v, "stands.db"), joinpath(dir, "stands.db"))
    key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
           for l in readlines(joinpath(_WT_FX, "$(name).key"))]
    write(joinpath(dir, "s.key"), join(key, '\n'))
    FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = _WT_VAR[v])
    return joinpath(dir, "out.db")
end

"Golden rows (Dict col => text) and jl rows (Dict col => value) of table `t`."
function _wt_table(name, db, t)
    gl = readlines(joinpath(_WT_FX, "$(name).$(t).csv")); hdr = split(gl[1], ',')
    gold = [Dict(zip(hdr, split(l, ','))) for l in gl[2:end]]
    d = SQLite.DB(db)
    jl = [Dict(String(k) => r[k] for k in propertynames(r)) for r in DBInterface.execute(d, "SELECT * FROM $t")]
    SQLite.close(d)
    return gold, jl
end
_wt_f32(x) = x isa AbstractString ? Float32(parse(Float64, x)) : Float32(x)
_wt_eq(a, b) = (a isa AbstractString && tryparse(Float64, a) === nothing) ? strip(a) == strip(string(b)) : _wt_f32(a) == _wt_f32(b)

"Compare columns `cols` of every golden row (keyed by Year (+ `extra` key columns)) with jl; returns mismatches."
function _wt_cmp(name, db, t, cols; extra = String[])
    gold, jl = _wt_table(name, db, t)
    key(r) = (string(r["Year"]), (string(r[c]) for c in extra)...)
    jd = Dict(key(r) => r for r in jl)
    bad = String[]
    length(jl) == length(gold) || push!(bad, "rows $(length(jl)) vs $(length(gold))")
    for g in gold, c in cols
        j = get(jd, key(g), nothing)
        (j !== nothing && _wt_eq(g[c], j[c])) || push!(bad, "$(key(g)) $c live $(g[c]) jl $(j === nothing ? "missing" : j[c])")
    end
    return bad
end

@testset "FMSCUT TONRMH/BIOREM(1) + FMDOUT TONREM: Total_Removed_Carbon / Biomass_Removed vs live" begin
    # fmscut.f:60-187 (TONRMH = Σ FMSVL2·V2T·WK3 − LVSNBM + TKCRBM; BIOREM(1) += ABIO·HRVTRE − LVCRWN), fmdout.f:274-291
    # (TONREM = TONRMS+TONRMH+TONRMC), fmcrbout.f:136-141 (V(10) = BIOREM(2) FFE / BIOREM(1) Jenkins). Before: 0.
    # Cut roots: BIOROOT += RBIO·(WK3−SSNG) in record order, then FMSADD's RBIO·SSNG (Belowground_Dead SN 1 ULP before).
    for (v, name) in (("ws", "ws_ffe_thin"), ("sn", "sn_ffe_thin"), ("sn", "sn_jenk_thin"))
        db = _wt_run(v, name)
        @test isempty(_wt_cmp(name, db, "FVS_Carbon", ["Total_Removed_Carbon", "Belowground_Dead"]))
        @test isempty(_wt_cmp(name, db, "FVS_Fuels", ["Biomass_Removed"]))
    end
end

