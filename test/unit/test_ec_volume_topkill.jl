# =============================================================================
# test_ec_volume_topkill.jl — EC broken/killed-top TOTAL-CUBIC truncation (CFTOPK) regression.
#
# ect01 stand-1 (S248112 UNTHINNED CONTROL) carries two broken-top trees (rec-6 WL D8.0 code 96,
# rec-22 DF D10.4 code 97). FVS (vols.f:145-146,193) computes their full cubic over the dubbed NORMAL
# height (cratet NORMHT), then CFTOPK trims it back to the standing break. jl's compute_volumes_ec!
# was taking the FULL normal-height cubic without the trim → the long-standing "TCuFt +2" corner
# (jl 1617 vs oracle 1615). MEASURED at the tree level + decisively: patching the oracle to DISABLE
# CFTOPK reproduces jl's 1617 exactly; enabling it gives 1615. Fix = call the shared r4_topkill (the
# CFTOPK/BFTOPK port already used by BM/CI/UT/TT) in EC's FW2 branch, EC TOPD=4.5.
#
# cyc0 has NO growth ⇒ deterministic ⇒ this row is a TRUE bit-exact test (not straddle-sensitive).
# Golden = FVSec_g16 (relinked oracle) 1990 row — NOT the stale ect01.sum.save different-build artifact.
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "EC — broken-top total-cubic CFTOPK trim (ect01 cyc0 bit-exact vs FVSec_g16)" begin
    fx  = joinpath(@__DIR__, "..", "fixtures", "eastcascades")
    key = joinpath(fx, "ecvol.key")
    if !isfile(key) || !isfile(joinpath(fx, "ecvol.tre"))
        @test_skip "eastcascades volume fixture not present"
    else
        row = cd(fx) do
            txt = F.run_keyfile("ecvol.key"; variant = F.EastCascades(), output = :sum)
            r = [strip(l) for l in split(txt, '\n') if occursin(r"^1990\s", strip(l))]
            isempty(r) ? String[] : split(r[1])
        end
        gold = split(strip(read(joinpath(fx, "ecvol_oracle.rows"), String)))

        @test !isempty(row)
        @test length(row) == length(gold)
        # Every .sum column bit-identical at cyc0 (Year Age TPA BA SDI CCF TopHt QMD TCuFt MCuFt _ BdFt ...).
        for (j, g) in enumerate(gold)
            @test row[j] == g
        end
        # Spotlight the volume columns the fix targets: TCuFt(9)=1615, MCuFt(10)=986, BdFt(12)=5073.
        @test row[9]  == "1615"   # total cubic — was 1617 before the CFTOPK trim
        @test row[10] == "986"    # merch cubic — unchanged by the trim (break above merch top)
        @test row[12] == "5073"   # board foot  — unchanged
    end
end
