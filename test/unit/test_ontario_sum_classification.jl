# =============================================================================
# test_ontario_sum_classification.jl — ON (Ontario) cyc0 .sum row bit-exact vs FVSon_g16,
# through the FORTYP / size-class / stocking-class classification columns.
#
# The ON variant is metric (cm DBH, m HT, trees/ha internally GROSPC-expanded). The .sum
# classification columns (FORTYP, and the trailing size|stock field) are computed by the
# shared stkval.f/fortyp.f port from the SAME expanded PROB that produces the bit-exact
# TPA/BA/QMD row — so getting them right is a whole-row invariant, not a cosmetic add-on.
#
# ont01 (8 inline metric trees, PW/SW/MH/BE/SB/CE/PJ/BF) classifies to FORTYP=125 (white
# pine), size class 1, stocking class 1 → the row tail "125 11" — byte-identical to the
# oracle. (An earlier port state emitted "999/55" = nonstocked; this locks the fix in.)
# Golden = the FVSon_g16 ont01.sum cyc0 data row.
# =============================================================================

using Test, FVSjl
const F = FVSjl

@testset "ON — cyc0 .sum row bit-exact vs FVSon_g16 (through FORTYP/size/stock)" begin
    fx  = joinpath(@__DIR__, "..", "fixtures", "ontario")
    key = joinpath(fx, "ont01.key")
    if !isfile(key)
        @test_skip "ontario fixture not present"
    else
        # run from the fixture dir so the bare TREEDATA keyword resolves ont01.tre
        row_j = cd(fx) do
            txt = F.run_keyfile("ont01.key"; variant = F.Ontario(), output = :sum)
            only(l for l in split(txt, '\n') if occursin(r"^\s+0\s+0\s+84799", l))
        end
        row_o = strip(read(joinpath(fx, "ont01_cyc0_oracle.row"), String), ['\n'])

        # the trailing "FORTYP size|stock" field is the headline this test guards
        @test split(strip(row_j))[end-1:end] == ["125", "11"]
        # and the entire cyc0 data row is byte-identical to the oracle (whitespace-normalized
        # to be robust to trailing spaces; every numeric field must match)
        @test split(strip(row_j)) == split(strip(row_o))
    end
end
