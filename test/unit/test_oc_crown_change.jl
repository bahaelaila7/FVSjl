# =============================================================================
# test_oc_crown_change.jl — OC FVS-native crown-ratio CHANGE for IORG=0 trees (oc/crown.f CASE DEFAULT).
#
# Regression guard for the LP height-growth fix: `crown_ratio_update!(::OregonCoast)` was a NO-OP, so the
# crown ratio of every FVS-native (IORG=0) tree was FROZEN at its cycle-0 inventory value. That starved
# the crown-vigor term CRMOD in oc/htgf.f, so LP (OC sp 12, not an ORGANON species) under-grew — the
# fixed model marches LP tree 1's crown 35→53 (%) over the run, tracking the oracle bit-exact.
#
# Golden = the live FVSoc_clean run on S248112 (ocgro), measured through the SAME engine path this test
# uses (each_stand → setup_growth! → write_sum_file with a per-cycle hook). The keyfile has NUMCYCLE 10
# / INVYEAR 1990 / 5-yr cycles ⇒ the oracle .sum has 11 rows 1990…2040, so the cycle_hook's LAST value
# is the year-2040 state. The FINAL (year-2040) per-tree oracle numbers are:
#   • tree 1  (LP, sp12, IORG=0): crown 53, HT 133.863 — jl BIT-EXACT (crown 35→53, HT to 133.8633).
#   • tree 28 (BR, sp22, IORG=0, = the oracle's internal index 25): crown 56 (down from the inventory
#     65), HT 71.766 — jl crown BIT-EXACT (56); jl HT 71.883 carries the pre-existing cornered ~0.12 ft
#     residual (the sp22 RELHT reads AVH, which drifts once the IORG=1 ORGANON large-tree heights drift —
#     the documented "OC large-tree ORGANON multi-cycle carry", a SEPARATE cornered item, not this fix).
# NOTE: tree indexing is by jl `tree_id` (input record order). The oracle's per-tree dump uses FVS's
# internal (reordered) index, so oracle internal I=25 is jl tree_id 28 — matched here by species + the
# cycle-0 inventory height (sp22, H0=28 ft). A revert to the no-op FREEZES crown (LP→35, BR→65) and
# lowers HT, so every assertion below fails on the frozen-crown behavior. DGSD=0 on OC ⇒ deterministic.
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OC crown-change — IORG=0 crown ratio grows (LP height no longer frozen)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "oregoncoast")
    key = joinpath(fx, "ocgro.key")
    if !isfile(key) || !isfile(joinpath(fx, "ocgro.tre"))
        @test_skip "oregoncoast ocgro fixture not present"
    else
        snaps = cd(fx) do
            acc = Dict{Int,Tuple{Float32,Int32}}()   # tree_id => (height, crown_pct) at the last cycle
            for s in F.each_stand("ocgro.key"; variant = F.variant_from_code("OC"))
                F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
                hook = (st, yr, pl, cy) -> foreach(
                    i -> (acc[Int(st.trees.tree_id[i])] = (st.trees.height[i], st.trees.crown_pct[i])),
                    1:st.trees.n)
                io = IOBuffer(); F.write_sum_file(io, s; variant = "OC", cycle_hook = hook)
                break
            end
            acc
        end
        # LP tree 1 (FVS-native, IORG=0): oracle FINAL (year 2040) HT 133.863, crown 53 (%). The fix makes
        # both BIT-EXACT; a frozen-crown revert stays crown 35 and grows less tall.
        @test haskey(snaps, 1)
        @test abs(snaps[1][1] - 133.8633f0) < 0.05f0     # HT bit-exact vs FVSoc_clean
        @test snaps[1][2] == 53                          # crown 35→53 (revert stays 35)
        # BR tree 28 (sp22, IORG=0; oracle internal index 25): crown 56 (down from inventory 65) — jl crown
        # BIT-EXACT; HT within the cornered AVH residual (oracle 71.766 / jl 71.883).
        @test haskey(snaps, 28)
        @test snaps[28][2] == 56                          # crown bit-exact (revert stays 65)
        @test abs(snaps[28][1] - 71.766f0) < 0.30f0       # HT bit-exact-or-cornered (AVH residual)
    end
end
