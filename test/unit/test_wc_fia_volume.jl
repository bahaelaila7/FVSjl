# test_wc_fia_volume.jl — WC FIA volume + crown width vs live FVSwc_g16 (cycle-0 FVS_TreeList, 4 stands).
#
# Guards three WC fixes (branch wc-vol):
#  * VOLEQ by forest (wc/sitset.f:218-243 → voleqdef.f R6_EQN/R7_EQN): jl had only Willamette (618) and sent
#    every other forest's DF/WH/GF/… to region-6 Behre ⇒ +25-45% cubic on 603/605/613/615 stands.
#  * Top-killed trees (wc/vols.f:145-146, 191-193, 390-391): volume at NORMHT, trimmed by CFTOPK/BFTOPK at
#    ITRUNC (FIA HTTOPK). jl volumed the full stem (an NF snapped at 5 ft: 1447 ft³ vs live 83.4).
#  * wc_cwcalc: every WCMAP equation via the national cwcalc library + the R6 forest bias factor by KODFOR
#    (jl had 7 forest-618 species and errored "02206 not yet ported" on the rest, e.g. noble fir).
# Fixture: a self-contained FIA sub-DB of stands 22404962010497 (LOC 613→605, two broken tops),
# 22405142010497 (603), 45080008020004 (615) and 45064714020004 (618); goldens = live FVS_TreeList cycle 0.
# Known 1-tree residual: 45080008020004 RC tree 103 (I00FW2W242 INGY) MCuFt 34.4 vs 34.0, BdFt 170 vs 169.
using FVSjl, Test, SQLite, DBInterface

@testset "WC FIA volume + crown width vs FVSwc_g16 (cycle-0 treelist)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "wc_fia_volume")
    db = abspath(joinpath(fx, "wc_fia_volume.db"))
    # TreeIds repeat within a multi-condition FIA stand, so compare per (stand, TreeId) as sorted multisets.
    gold = Dict{Tuple{String,String},Vector{NTuple{4,Float64}}}()
    for l in eachline(joinpath(fx, "live_cyc0_treelist.csv"))
        startswith(l, "stand") && continue
        f = split(l, ',')
        push!(get!(gold, (String(f[1]), String(f[3])), NTuple{4,Float64}[]),
              (parse(Float64, f[5]), parse(Float64, f[6]), parse(Float64, f[7]), parse(Float64, f[8])))
    end
    for cn in unique(first.(keys(gold)))
        dir = mktempdir()
        key = join(["STDIDENT", cn, "DATABASE", "DSNin", db,
                    "StandSQL", "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
                    "TreeSQL", "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
                    "END", "NUMCYCLE           1", "TREELIST           0",
                    "DATABASE", "TREELIDB", "END", "PROCESS", "STOP"], '\n') * '\n'
        write(joinpath(dir, "k.key"), key)
        cd(dir) do
            FVSjl.run_keyfile("k.key"; variant = FVSjl.WestCascades(), output = :sum)
        end
        out = SQLite.DB(joinpath(dir, "FVSOut.db"))
        y0 = first(DBInterface.execute(out, "SELECT MIN(Year) AS y FROM FVS_TreeList")).y
        got = Dict{String,Vector{NTuple{4,Float64}}}()
        for r in DBInterface.execute(out, "SELECT TreeId, TCuFt, MCuFt, BdFt, CrWidth FROM FVS_TreeList WHERE Year = $y0")
            push!(get!(got, String(strip(r.TreeId)), NTuple{4,Float64}[]), (r.TCuFt, r.MCuFt, r.BdFt, r.CrWidth))
        end
        SQLite.close(out)
        for ((c, id), gs) in gold
            c == cn || continue
            js = get(got, id, NTuple{4,Float64}[])
            @test length(js) == length(gs)
            known = cn == "45080008020004" && id == "103"
            for (g, j) in zip(sort(gs), sort(js))
                @test j[1] ≈ g[1] rtol = 1e-5 atol = 1e-3
                known || (@test j[2] ≈ g[2] rtol = 1e-5 atol = 1e-3)
                known || (@test j[3] ≈ g[3] rtol = 1e-5 atol = 0.51)
                @test j[4] ≈ g[4] rtol = 1e-5
            end
        end
    end
end
