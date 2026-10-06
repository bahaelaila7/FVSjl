# test_west_tail3.jl — west-tail round 3, each fix checked against LIVE FVS{v}_g16 goldens in test/fixtures/west_tail3
# (keys run on the tiered fixtures' stands.db; goldens dumped from the oracle's out.db, CaseID dropped). Float32 bits.
using FVSjl, Test, SQLite, DBInterface

const _W3_FX = joinpath(@__DIR__, "..", "fixtures", "west_tail3")
const _W3_VAR = Dict("ci" => FVSjl.CentralIdaho(), "nc" => FVSjl.Klamath(), "bm" => FVSjl.BlueMountains(),
                     "ws" => FVSjl.WestSierra(), "wc" => FVSjl.WestCascades(), "ec" => FVSjl.EastCascades())

"Run west_tail3/<name>.key for variant `v` (stands.db from the tiered fixture); returns the output DB path."
function _w3_run(v::AbstractString, name::AbstractString)
    dir = mktempdir(); cp(joinpath(@__DIR__, "..", "fixtures", "tiered", v, "stands.db"), joinpath(dir, "stands.db"))
    key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
           for l in readlines(joinpath(_W3_FX, "$(name).key"))]
    write(joinpath(dir, "s.key"), join(key, '\n'))
    FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = _W3_VAR[v])
    return joinpath(dir, "out.db")
end

_w3_f32(x) = x isa AbstractString ? Float32(parse(Float64, x)) : Float32(x)
_w3_eq(a, b) = (a isa AbstractString && tryparse(Float64, a) === nothing) ? strip(a) == strip(string(b)) :
               (b !== missing && _w3_f32(a) == _w3_f32(b))

"Mismatches of every golden row/column of table `t` (rows keyed by their order) against jl's out.db."
function _w3_cmp(name, db, t; skip = String[])
    gl = readlines(joinpath(_W3_FX, "$(name).$(t).csv")); hdr = split(gl[1], ',')
    gold = [split(l, ','; keepempty = true) for l in gl[2:end]]
    d = SQLite.DB(db)
    jl = [Dict(String(k) => r[k] for k in propertynames(r)) for r in DBInterface.execute(d, "SELECT * FROM $t")]
    SQLite.close(d)
    bad = String[]
    length(jl) == length(gold) || push!(bad, "rows $(length(jl)) vs $(length(gold))")
    for (i, g) in enumerate(gold), (k, c) in enumerate(hdr)
        c in skip && continue
        i <= length(jl) || continue
        _w3_eq(g[k], jl[i][c]) || push!(bad, "row $i $c live $(g[k]) jl $(jl[i][c])")
    end
    return bad
end

# FUELDCAY / FUELMULT + the first-FFE-year decay multiplier. ci/fmcba.f:407-445, nc/fmcba.f:395-414 and ws/fmcba.f:573-580 (DCYMLT) and
# bm|ec|wc|pn|op/fmcba.f (DKRADJ, `(SETDECAY(I,J).LT.0).AND.(ISWTCH.NE.1)`) scale ONLY the elements
# with SETDECAY(I,J) < 0 (the ones FUELDCAY did not set, fmin.f:831-860) of the DKR in force (the variant's fmvinit.f table as
# FUELMULT/FUELDCAY edited it, fmin.f:806-870/1368) by DCYMLT. Before: the keyword path copied the SN table (not the variant's)
# and any FUELDCAY/FUELMULT disabled the multiplier for the whole matrix.
@testset "FUELDCAY/FUELMULT + DCYMLT/DKRADJ per-element SETDECAY semantics (fmcba.f, fmin.f:806-870)" begin
    for (v, cn) in (("ci", "5385337010690"), ("nc", "15303130010497"), ("bm", "22960873010497"), ("ws", "7689398010901"),
                    ("wc", "214084313020004"), ("ec", "248652369489998"))
        name = "fd_$cn"
        db = _w3_run(v, name)
        for t in ("FVS_Fuels", "FVS_Carbon", "FVS_PotFire", "FVS_BurnReport", "FVS_Summary")
            @test isempty(_w3_cmp(name, db, t; skip = String[]))
        end
    end
end
