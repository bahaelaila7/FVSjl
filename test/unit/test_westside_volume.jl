# test_westside_volume.jl — NC/CA/SO/WS/EC cycle-0 per-tree volume vs live FVS{v}_g16 (FVS_TreeList).
#
# Guards the westside-vol fixes (branch westside-vol):
#  * R5HARV/DVEST (500DVEW California hardwoods, NC/CA/SO/WS R5 forests): dvest.f:170 VOL(2)=ANINT(VOL(2)) —
#    Scribner rounded to the whole foot — plus r5harv.f's REAL*8-coefficient/REAL-body precision, the red-alder
#    (351) and giant-sequoia board branches, and volinit.f:559 VOL(7)≥0 (big giant chinkapin MCF).
#  * profile.f:293 VOL(1)=NINT(TCVOL*10.0)*1E-1 (×REAL 0.1, not /10) for the WO2W/westside-FW2 total cubic.
#  * CA VOLEQ by forest (ca/sitset.f VOLEQDEF): 518 is R5; 710/711/712 are BLM B00BEHW/B01|B02BEHW202 → BLMVOL
#    (jl ran them through the R6 Behre form-class path — every BLM tree wrong).
#  * SO CRATET: WJ/WB/AS (ISPC 11/16/24) skip HTDBH (so/cratet.f:415/532) ⇒ Wykoff HT1/HT2 NORMHT dub for
#    broken tops (jl dubbed 4.5 ft ⇒ volume at the broken height).
#  * EC BRATIO: ISPC 1:19,31 are the constant BARK1 (ec/bratio.f), not BARK1·D/D (1-ULP off ⇒ BEHPRM-amplified
#    broken-top cubic).
#  * TreeList crown width routed as cwcalc.f does: CA IFOR≤5 → R5CRWD (:385) else CAMAP (+12 national codes);
#    NC IFOR≤3/5 → R5CRWD MAPNC (:382) else NCMAP + Siskiyou BF; SO IFOR 4-9 → R5CRWD MAPSO (:376) else SOMAP + KODFOR BF.
# Fixture: a sub-DB with one synthetic stand per distinct live VOLEQ table (every species × 3 sizes + 2 broken
# tops + 2 missing heights) and two FIA stands per variant; goldens = live FVS_TreeList at the inventory year
# (TCuFt/MCuFt/BdFt/CrWidth, all five variants).
using FVSjl, Test, SQLite, DBInterface

@testset "westside NC/CA/SO/WS/EC cycle-0 volume vs live (FVS_TreeList)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "westside_volume")
    db = abspath(joinpath(fx, "westside_volume.db"))
    ctor = Dict("nc" => FVSjl.Klamath, "ca" => FVSjl.CentralCalifornia, "so" => FVSjl.SouthCentralOregon,
                "ws" => FVSjl.WestSierra, "ec" => FVSjl.EastCascades)
    # TreeIds repeat within a multi-condition FIA stand ⇒ compare per (stand, TreeId) as sorted multisets.
    gold = Dict{String,Dict{Tuple{String,String},Vector{NTuple{4,Float64}}}}()
    for l in eachline(joinpath(fx, "live_cyc0_treelist.csv"))
        startswith(l, "variant") && continue
        f = split(l, ',')
        g = get!(gold, String(f[1]), Dict{Tuple{String,String},Vector{NTuple{4,Float64}}}())
        push!(get!(g, (String(f[2]), String(f[4])), NTuple{4,Float64}[]),
              Tuple(parse.(Float64, f[6:9])))
    end
    for v in ("nc", "ca", "so", "ws", "ec")
        g = gold[v]
        stands = unique(first.(keys(g)))
        dir = mktempdir()
        key = join([join(["STDIDENT", cn, "DATABASE", "DSNin", db,
                          "StandSQL", "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
                          "TreeSQL", "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
                          "END", "NUMCYCLE           1", "TREELIST           0",
                          "DATABASE", "TREELIDB", "END", "PROCESS"], '\n') for cn in stands], '\n') * "\nSTOP\n"
        write(joinpath(dir, "k.key"), key)
        cd(dir) do
            FVSjl.run_keyfile("k.key"; variant = ctor[v](), output = :sum)
        end
        out = SQLite.DB(joinpath(dir, "FVSOut.db"))
        y0 = Dict(String(strip(r.StandID)) => r.y for r in
                  DBInterface.execute(out, "SELECT StandID, MIN(Year) AS y FROM FVS_TreeList GROUP BY StandID"))
        got = Dict{Tuple{String,String},Vector{NTuple{4,Float64}}}()
        for r in DBInterface.execute(out, "SELECT StandID, Year, TreeId, TCuFt, MCuFt, BdFt, CrWidth FROM FVS_TreeList")
            st = String(strip(r.StandID))
            r.Year == y0[st] || continue
            push!(get!(got, (st, String(strip(r.TreeId))), NTuple{4,Float64}[]), (r.TCuFt, r.MCuFt, r.BdFt, r.CrWidth))
        end
        SQLite.close(out)
        @testset "$v" begin
            @test length(got) == length(g)
            for (k, gs) in g
                js = get(got, k, NTuple{4,Float64}[])
                @test length(js) == length(gs)
                for (a, b) in zip(sort(gs), sort(js))
                    @test b[1] ≈ a[1] rtol = 1e-5 atol = 1e-4
                    @test b[2] ≈ a[2] rtol = 1e-5 atol = 1e-4
                    @test b[3] ≈ a[3] rtol = 1e-5 atol = 1e-3
                    @test b[4] ≈ a[4] rtol = 1e-5
                end
            end
        end
    end
end
