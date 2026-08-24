# =============================================================================
# test_oc_crown_change.jl — OC FVS-native crown-ratio CHANGE for IORG=0 trees (oc/crown.f CASE DEFAULT).
#
# Regression guard for the LP height-growth fix: `crown_ratio_update!(::OregonCoast)` was a NO-OP, so the
# crown ratio of every FVS-native (IORG=0) tree was FROZEN at its cycle-0 inventory value. That starved
# the crown-vigor term CRMOD in oc/htgf.f, so LP (OC sp 12, not an ORGANON species) under-grew ~0.15
# ft/cycle — compounding to −3.25 ft by cycle 9 on the S248112 (ocgro) stand.
#
# Golden = the live FVSoc_clean run (scratchpad/lp_height, instrumented htgf/crown dumps): LP tree 1
# marches crown 0.35→0.53 over 10 cycles and reaches HT 128.81; before the fix it stayed CR 0.35 / HT
# 125.56. Tolerances ACCEPT the fixed (bit-exact) state and FAIL a revert to the frozen-crown behavior.
# DGSD=0 on OC ⇒ deterministic ⇒ this is a true bit-exact target.
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
        # LP tree 1 (FVS-native, IORG=0): oracle final HT 128.81, crown 53 (%). Frozen-crown revert → HT
        # ~125.56, crown 35 — the tolerances below distinguish the two decisively.
        @test haskey(snaps, 1)
        @test abs(snaps[1][1] - 128.81f0) < 0.20f0      # HT bit-exact-or-cornered (revert = −3.25 ft)
        @test snaps[1][2] >= 50                          # crown grew from 35 (revert stays 35)
        # BR tree 25 (also IORG=0) — crown grew from 65 toward 57; height within the cornered AVH residual.
        @test haskey(snaps, 25)
        @test snaps[25][2] <= 60 && snaps[25][2] >= 55   # crown 57 (was frozen 65)
        @test abs(snaps[25][1] - 67.49f0) < 0.30f0
    end
end
