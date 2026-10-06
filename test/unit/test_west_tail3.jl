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

# crown.f `IF((ITRN.LE.0).AND.(IREC2.LT.MAXTP1)) GO TO 74` (every variant): a stand whose inventory is ALL dead records still dubs
# the dead records' crowns (DO 79; the DUBSCR/BACHLO draws are consumed, the projection's RNG stays in step). jl returned early at
# n == 0 in KT/AK/BC/BM/EM/OC/OP (PctCr 0 on every record) and, for SN/CS/LS/NE, its CRATET init returned at n == 0.
# dead_<v>: the variant's busiest tiered FIA/local stand with every record HISTORY 6, crown blank, no-HT/<1" records dropped
# (an all-dead FVSkt cwcalc SIGFPEs on a zero-height record: ALOG(H)); goldens = live oracle FVS_TreeList (before: PctCr
# mismatches BC 1540/1540, AK 43/43, KT 527/528, BM 55/55, EM 803/814, SN 20/20, CS 59/59, LS 55/55, NE 37/37, OC 73/73, OP 71/71).
const _W3_ALLVAR = merge(_W3_VAR, Dict("bc" => FVSjl.BritishColumbia(), "ak" => FVSjl.SoutheastAlaska(), "kt" => FVSjl.Kootenai(),
    "bm" => FVSjl.BlueMountains(), "em" => FVSjl.EasternMontana(), "sn" => FVSjl.Southern(), "cs" => FVSjl.CentralStates(),
    "ls" => FVSjl.LakeStates(), "ne" => FVSjl.Northeast(), "oc" => FVSjl.OregonCoast(), "op" => FVSjl.Olympic()))

@testset "all-dead inventory dubs the dead crowns (crown.f GO TO 74 / DO 79)" begin
    for v in ("bc", "ak", "kt", "bm", "em", "sn", "cs", "ls", "ne", "oc", "op")
        dir = mktempdir(); cp(joinpath(_W3_FX, "dead_$(v).stands.db"), joinpath(dir, "stands.db"))
        key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
               for l in readlines(joinpath(_W3_FX, "dead_$(v).key"))]
        write(joinpath(dir, "s.key"), join(key, '\n'))
        FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = _W3_ALLVAR[v])
        tbl = v == "bc" ? "FVS_TreeList_Metric" : "FVS_TreeList"
        gl = readlines(joinpath(_W3_FX, "dead_$(v).$(tbl).csv"))[2:end]
        gold = [(parse(Int, f[1]), f[2], f[3], parse(Float64, f[4]), parse(Int, f[5])) for f in (split(l, ',') for l in gl)]
        d = SQLite.DB(joinpath(dir, "out.db"))
        jl = [(Int(r.Year), string(r.TreeId), string(r.SpeciesFVS), Float64(r.DBH), Int(r.PctCr))
              for r in DBInterface.execute(d, "SELECT Year,TreeId,SpeciesFVS,DBH,PctCr FROM $tbl ORDER BY Year,TreeId,SpeciesFVS,DBH,PctCr")]
        SQLite.close(d)
        @test length(jl) == length(gold)
        @test all(a[1:3] == b[1:3] && a[4] ≈ b[4] && a[5] == b[5] for (a, b) in zip(gold, jl))
    end
end
