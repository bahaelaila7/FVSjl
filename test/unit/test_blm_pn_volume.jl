# test_blm_pn_volume.jl — NVEL BLMVOL (Oregon BLM forests) + PN per-forest VOLEQ/merch/crown width, vs live
# FVSwc_g16 / FVSpn_g16 cycle-0 FVS_TreeList (TCuFt, MCuFt, BdFt, CrWidth).
#
# Guards (branch blm-vol):
#  * BLMVOL (blmvol.f/blmtap.f): VOLEQs starting with 'B' (WC 708-711, PN 708/709/712: B00BEHW<fia>,
#    B01/B02BEHW202) were volumed with region-6 Behre.
#  * PN VOLEQ by forest (pn/sitset.f:213-243): jl knew only Siuslaw (612). Olympic/Quinault (609/800) put SS on
#    F03FW2W263 and DF on F03FW2W202; RA on 612 is NVBM240351 (NSVB); BLM merch specs 5"/7" (pn/sitset.f:189).
#  * pn_cwcalc: pn/cwcalc.f == wc/cwcalc.f, evaluated with PN's KODFOR (was 7 species, error on the rest).
#  * NVB broken tops (fvsvol.f:85-88 BRKHT=ITRNC/100 for LIVE top-killed trees; nsvb.f VOL(1)=Vtotib·Rrem and
#    HT1PRD capped at BRKHT; vols.f:191 no CFTOPK for 'NVB', BFTOPK with BFMAX=truncated TVOL(1)): jl ran CFTOPK
#    on the full NVB volume (PN RA on 612; CR region-3 CB/ES/WF/AS — 89 broken-top CR records were off).
# Fixture: FIA sub-DB with WC stands on 708/709/710/711, PN stands on 609 (SS), 612 (RA), 800, 708, 709, 712, and
# three CR region-3 (forest 302/308) stands with top-killed NVB-equation trees.
using FVSjl, Test, SQLite, DBInterface

@testset "BLMVOL + PN + NVB broken-top volume/crown width vs live (cycle-0 treelist)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "blm_pn_volume")
    db = abspath(joinpath(fx, "blm_pn_volume.db"))
    # TreeIds repeat within a multi-condition FIA stand: compare per (stand, TreeId) as sorted multisets.
    gold = Dict{Tuple{String,String,String},Vector{NTuple{4,Float64}}}()
    for l in eachline(joinpath(fx, "live_cyc0_treelist.csv"))
        startswith(l, "variant") && continue
        f = split(l, ',')
        push!(get!(gold, (String(f[1]), String(f[2]), String(f[4])), NTuple{4,Float64}[]),
              (parse(Float64, f[6]), parse(Float64, f[7]), parse(Float64, f[8]), parse(Float64, f[9])))
    end
    ctor = Dict("WC" => FVSjl.WestCascades(), "PN" => FVSjl.PacificNorthwest(), "CR" => FVSjl.CentralRockies())
    # Measured corner: a dead SS (HISTORY=8) broken 10 ft up a 130-ft NORMHT. CFTOPK's merch ratio is
    # BEHRE(PHT,STUMP) = 4.8e-4 left after cancelling O(1) terms (behre.f), so a 1-3 ULP change in VMAX/BARK swings
    # MCF by ±0.006 (jl 49.6967, live 49.7029, Float64 49.6948). TCF/BF/CrWidth are exact.
    mcf_atol = Dict(("PN", "22961767010497", "324") => 0.01)
    for (v, cn) in unique([(k[1], k[2]) for k in keys(gold)])
        dir = mktempdir()
        key = join(["STDIDENT", cn, "DATABASE", "DSNin", db,
                    "StandSQL", "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
                    "TreeSQL", "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
                    "END", "NUMCYCLE           1", "TREELIST           0",
                    "DATABASE", "TREELIDB", "END", "PROCESS", "STOP"], '\n') * '\n'
        write(joinpath(dir, "k.key"), key)
        cd(dir) do
            FVSjl.run_keyfile("k.key"; variant = ctor[v], output = :sum)
        end
        out = SQLite.DB(joinpath(dir, "FVSOut.db"))
        y0 = first(DBInterface.execute(out, "SELECT MIN(Year) AS y FROM FVS_TreeList")).y
        got = Dict{String,Vector{NTuple{4,Float64}}}()
        for r in DBInterface.execute(out, "SELECT TreeId, TCuFt, MCuFt, BdFt, CrWidth FROM FVS_TreeList WHERE Year = $y0")
            push!(get!(got, String(strip(r.TreeId)), NTuple{4,Float64}[]), (r.TCuFt, r.MCuFt, r.BdFt, r.CrWidth))
        end
        SQLite.close(out)
        @testset "$v $cn" begin
            for ((vv, c, id), gs) in gold
                (vv == v && c == cn) || continue
                js = get(got, id, NTuple{4,Float64}[])
                @test length(js) == length(gs)
                for (g, j) in zip(sort(gs), sort(js))
                    @test j[1] ≈ g[1] rtol = 1e-5 atol = 1e-3
                    @test j[2] ≈ g[2] rtol = 1e-5 atol = get(mcf_atol, (vv, c, id), 1e-3)
                    @test j[3] ≈ g[3] rtol = 1e-5 atol = 0.51
                    @test j[4] ≈ g[4] rtol = 1e-5
                end
            end
        end
    end
end
