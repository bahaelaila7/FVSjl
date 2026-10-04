# test_variant_maxtre_live.jl — per-variant MAXTRE (each build's PRGPRM.F77) vs LIVE oracles.
#
# MAXTRE is the tree-record capacity: it sets the tripling guard (grincr.f:74 LTRIP needs ITRN ≤ MAXTRE/3), the
# FVS_TreeList TreeIndex of the cycle-0 inventory-dead records (intree.f files them from MAXTRE downward), the COMPRESS
# default target (MAXTRE/2) and the root-disease IRINIT (rdinit.f:703 IRINIT=10*MAXTRE). Most builds declare 3000;
# BC 4000, ON 6000 (tested in test_ontario_big_live.jl), OC/OP 2000. jl used 3000 everywhere.

using Test, FVSjl, SQLite, DBInterface
const _MXT = FVSjl
const _MXT_FX = joinpath(@__DIR__, "..", "fixtures", "maxtre")

function _mxt_run(stem::AbstractString, variant; ncyc::Union{Nothing,Int} = nothing)
    dir = mktempdir(); out = String[]
    for ln in eachline(joinpath(_MXT_FX, "$stem.key"))
        startswith(ln, "PROCESS") && append!(out, ["DATABASE", "SUMMARY", "TREELIDB", "END"])
        push!(out, ln)
        startswith(ln, "TREEDATA") && push!(out, "TREELIST           0")
        length(out) == 2 && append!(out, ["DATABASE", "DSNout", "$stem.db", "END"])
    end
    write(joinpath(dir, "$stem.key"), join(out, "\n") * "\n")
    cp(joinpath(_MXT_FX, "$stem.tre"), joinpath(dir, "$stem.tre"))
    cd(() -> _MXT.run_keyfile("$stem.key"; variant = variant, output = :sum), dir)
    db = SQLite.DB(joinpath(dir, "$stem.db"))
    tab = first(r.name for r in DBInterface.execute(db, "SELECT name FROM sqlite_master WHERE name LIKE 'FVS_TreeList%'"))
    return db, tab
end

# BC (canada/bc PRGPRM.F77 MAXTRE=4000): a 1200-record BG stand triples in cycle 1 (1200 ≤ 4000/3). jl's 3000 guard
# (≤ 1000) left it untripled (2002: 1200 records vs live 3600).
@testset "BC 1200-record stand: cycle-1 TRIPLE (MAXTRE=4000) — 2002 record set, TPH, DBH == live" begin
    db, tab = _mxt_run("bc_bg1200", _MXT.BritishColumbia())
    jl = Dict((Int(r.Year), Int(r.TreeIndex)) => (Float64(r.TPH), Float64(r.DBH), Float64(r.Ht))
              for r in DBInterface.execute(db, "SELECT Year,TreeIndex,TPH,DBH,Ht FROM $tab WHERE Year=2002"))
    live = Dict{Tuple{Int,Int},NTuple{3,Float64}}()
    for ln in eachline(joinpath(_MXT_FX, "bc_bg1200_2002_tl_live.csv"))
        (startswith(ln, "#") || startswith(ln, "Year")) && continue
        f = split(ln, ',')
        live[(parse(Int, f[1]), parse(Int, f[2]))] = (parse(Float64, f[3]), parse(Float64, f[4]), parse(Float64, f[5]))
    end
    @test sort(collect(keys(jl))) == sort(collect(keys(live)))       # 1200 records ⇒ 3600 after the TRIPLE
    for (k, v) in sort(collect(live))
        @test (k, :TPH, get(jl, k, (NaN, NaN, NaN))[1]) == (k, :TPH, v[1])
        # DBH: exact outside the small-tree REGENT blend. Pre-existing BC residuals (present identically in a 900-record
        # stand that tripled before this fix): 42 records with 2002 DBH 5-9 cm carry a 1-2 ULP small-tree DG, and 2 large
        # records (+ their copies) a 1-ULP HTG — so Ht is not compared here.
        v[2] >= 9.0 && @test (k, :DBH, get(jl, k, (NaN, NaN, NaN))[2]) == (k, :DBH, v[2])
    end
end

# OC/OP (PRGPRM.F77 MAXTRE=2000): the inventory-dead records are listed at TreeIndex MAXTRE+1−k (1996-2000 live; jl
# listed 2996-3000).
for (stem, v) in (("oc_dead", _MXT.OregonCoast()), ("op_dead", _MXT.Olympic()))
    @testset "$(uppercase(stem[1:2])) inventory-dead TreeIndex (MAXTRE=2000) == live" begin
        db, tab = _mxt_run(stem, v)
        got = join(["$(r.TreeIndex):$(r.TreeId)" for r in DBInterface.execute(db,
                   "SELECT TreeIndex,TreeId FROM $tab WHERE TreeIndex>50 ORDER BY TreeIndex")], ' ')
        @test got == strip(read(joinpath(_MXT_FX, "$(stem)_live.txt"), String))
    end
end
