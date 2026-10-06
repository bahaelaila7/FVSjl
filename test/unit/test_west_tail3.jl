# test_west_tail3.jl — west-tail round 3, each fix checked against LIVE FVS{v}_g16 goldens in test/fixtures/west_tail3
# (keys run on the tiered fixtures' stands.db; goldens dumped from the oracle's out.db, CaseID dropped). Float32 bits.
module WestTail3Test
using FVSjl, Test, SQLite, DBInterface
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

const _W3_FX = joinpath(@__DIR__, "..", "fixtures", "west_tail3")

# a tiered fixture case (goldens = the live oracle's, test/fixtures/tiered/<v>) compared whole
function _w3_case(v, cn, r)
    d = mktempdir()
    txt, db, crashed, err = run_case(v, cn, r; dir = d)
    return (crashed = crashed, ms = crashed ? Mismatch[] : compare_case(v, cn, r, txt, db))
end
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

# AK ESTAB (ak/estab.f:699/816) reads BAAA/PTPAA of the gradd.f:192 DENSE, which walked IND1 in the post-TRIPLE REASS interleave
# (grincr.f:553); esnutr.f:138's SPESRT after ESUCKR's sprouts resets IND1 to ascending record order, which jl used for its
# recompute: BAAA 111.172935 vs live 111.172951 (FVSak_g16 10709344010497 simfire 2016, instrumented dense.f/estab.f relink)
# => ESTOCK PN/PROB1 ULPs => every regen record's PROB and the later QMD/ATQMD 8 cells (2016-2046).
@testset "AK ESTAB BAAA/PTPAA summed in the gradd.f:192 DENSE's IND1 order" begin
    c = _w3_case("AK", "10709344010497", "simfire")
    @test !c.crashed
    @test isempty(c.ms)
end

# The .sum driver compacts the list at CUTS entry (cuts.f:255-275 TREDEL) when activities are due, ahead of grow_cycle!'s DENSE.
# FVS runs no DENSE after a removal-free CUTS and TREMOV does not move PTBALT, so DGDRIV/MORTS read the previous DENSE's
# PTBALT by SLOT. jl recomputed it over the compacted list: AK 666741868126144 plant 2026-2036 species-11 seedlings PTBAL
# 319.1590 vs live 319.1517 (MRT trace of an instrumented morts.f), RIP/WKI ULPs => ESTAB SUM0 (3.0983538 vs 3.0983527) =>
# ESA/PROB1 => 194 records' PROB, QMD/ATQMD 2066.
@testset "PTBALT stays slot-indexed across a removal-free CUTS-entry TREDEL (cuts.f:255, grincr.f:305)" begin
    for r in ("plant_cal", "plant_cyc")
        c = _w3_case("AK", "666741868126144", r)
        @test !c.crashed
        @test isempty(c.ms)
    end
end

# grincr.f:305-320 saves ATAVH/ATBA/ATCCF right after CUTS (and its thin-only DENSE), BEFORE COMCUP (grincr.f:391) whose DENSE
# (NDEL>0) refreshes AVH/BA/RELDEN. ESTAB's REGENT(LESTB) reads that saved trio (TEMAHT/TEMCCF -> DELMAX -> DADJ -> the birth DBH).
# jl captured it after COMCUP: KT 196396140020004 simfire cycle 3 ATAVH 4.3966694 vs live 4.3969569 (56 PROB<=1E-5 stubs deleted),
# DADJ 0.86324865 vs 0.86326247, every birth-cohort DBH 1.38e-5 low => QMD/ATQMD, FVS_Carbon, FVS_PotFire (14 KT cells).
@testset "ATAVH/ATBA/ATCCF are saved before COMCUP (grincr.f:316-320 vs :391)" begin
    c = _w3_case("KT", "196396140020004", "simfire")
    @test !c.crashed
    @test isempty(c.ms)
end

# dgdriv.f:186-205: `IF(I1.EQ.0) GO TO 50` precedes the COR/HCOR attenuation, so a species with no records at a DGDRIV is not
# advanced; its COR keeps the value the last present cycle left and DGF uses that stale value when the species returns.
# jl re-evaluated the elapsed-time curve every cycle for the present species: WS 7689398010901 thinbba species 29 (all thinned at
# the cycle-2 CUTS, back as sprouts): live DGF CONSPP -1.485131 (COR -0.05918458, instrumented dgdriv.f COR0/COR1 relink) vs jl -1.4789,
# 126 sp-29 records' DBH 1e-3 off => QMD/ATQMD/StrClass (WS 12 cells on thinbba + econ).
@testset "COR is advanced only for species present at the previous DGDRIV (dgdriv.f:186-205)" begin
    for r in ("thinbba", "econ")
        c = _w3_case("WS", "7689398010901", r)
        @test !c.crashed
        @test isempty(c.ms)
    end
end

# canada/bc/regent.f:1349 BAL=(1-PCT(I)/100)*BAJ with PCT=0 for ESTAB's new records (the post-ESTAB DENSE sets it later). jl read the
# new records' stale percentile: BAL 0.98792 vs live 0.99426961 (instrumented regent.f SUB trace), HTGR1 3e-6 low => QMD/ATQMD (BC 4 cells).
@testset "BC REGENT(LESTB) BAL uses PCT=0 for the new records (canada/bc/regent.f:1349)" begin
    for r in ("plant_cal", "plant_cyc")
        c = _w3_case("BC", "SkyRanch-0.3m", r)
        @test !c.crashed
        @test isempty(c.ms)
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

end # module WestTail3Test
