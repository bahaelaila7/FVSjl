# PPE MXHRVP multistand harvest scheduling — HVSEL bit-exact vs gfortran-16 driver-golden.
# Goldens produced by scratchpad/ppe/mxhrvp/driver_hvsel.f (the real pristine hvsel.f linked
# over the historical FVSppe objects, bc6e2377^): each case feeds synthetic policy-1 inputs
# (priority, yield-if-selected, yield-if-not, target, EXACT flag) and dumps IHVSTA + HVPART.
using Test
using FVSjl: hvsel!, hvccut, eval_policy_expr, HarvestVars, EventCtx, eval_event, parse_event_condition, lbmemr

# One golden: (label, lprtct, target, [(priority, yield_sel, yield_notsel)...],
#             expected_status (original stand order), expected_hvpart_hex::UInt32)
const HVSEL_GOLDENS = [
    ("C1 greedy, no-partial",     false, 150.0f0,
        [(10f0,100f0,5f0),(8f0,100f0,5f0),(6f0,100f0,5f0),(4f0,100f0,5f0)],
        [3,2,2,2], 0x00000000),
    ("C2 greedy, EXACT partial",  true,  150.0f0,
        [(10f0,100f0,5f0),(8f0,100f0,5f0),(6f0,100f0,5f0),(4f0,100f0,5f0)],
        [3,4,2,2], 0x3eb33333),                       # HVPART = 0.35
    ("C3 target met by thinning", false, 10.0f0,
        [(10f0,100f0,5f0),(8f0,100f0,5f0),(6f0,100f0,5f0),(4f0,100f0,5f0)],
        [2,2,2,2], 0x00000000),
    ("C4 target 0",               false, 0.0f0,
        [(5f0,50f0,2f0),(4f0,50f0,2f0),(3f0,50f0,2f0)],
        [2,2,2], 0x00000000),
    ("C5 unsorted pri, no-part",  false, 250.0f0,
        [(6f0,100f0,5f0),(10f0,100f0,5f0),(4f0,100f0,5f0),(8f0,100f0,5f0)],
        [2,3,2,3], 0x00000000),
    ("C6 unsorted pri, partial",  true,  250.0f0,
        [(6f0,100f0,5f0),(10f0,100f0,5f0),(4f0,100f0,5f0),(8f0,100f0,5f0)],
        [4,3,2,3], 0x3ecccccd),                       # HVPART = 0.4
    ("C7 round-down at pneed=0.5", false, 160.0f0,
        [(9f0,100f0,5f0),(5f0,100f0,5f0),(1f0,100f0,5f0)],
        [3,2,2], 0x00000000),
    ("C8 tiny-yield stand skip",  false, 50.0f0,
        [(9f0,0.0000001f0,5f0),(5f0,100f0,5f0),(1f0,100f0,5f0)],
        [2,2,2], 0x00000000),
]

@testset "PPE MXHRVP — HVSEL selection bit-exact vs FVSppe driver-golden" begin
    for (label, lprtct, target, stands, exp_status, exp_hvpart) in HVSEL_GOLDENS
        n = length(stands)
        status = fill(1, n)                                   # all candidates (±1 → +1)
        priority     = Float32[s[1] for s in stands]
        yield_sel    = Float32[s[2] for s in stands]
        yield_notsel = Float32[s[3] for s in stands]
        st, hvpart = hvsel!(status, priority, yield_sel, yield_notsel, target; lprtct = lprtct)
        @testset "$label" begin
            @test st == exp_status
            @test reinterpret(UInt32, hvpart) == exp_hvpart
        end
    end
end

@testset "PPE MXHRVP — HVCCUT clearcut-condition test (hvccut.f)" begin
    # -1 iff after-thin TPA < DCTPA AND after-thin top height < DCTOP, else +1.
    @test hvccut(5f0,  20f0, 10f0, 30f0) == -1   # both below thresholds → clearcut
    @test hvccut(15f0, 20f0, 10f0, 30f0) == 1    # TPA not below → not clearcut
    @test hvccut(5f0,  40f0, 10f0, 30f0) == 1    # topht not below → not clearcut
    @test hvccut(10f0, 20f0, 10f0, 30f0) == 1    # TPA == DCTPA (strict <) → not clearcut
    @test hvccut(0.001f0, 0.001f0, 0f0, 0f0) == 1  # DEFCCUT unset (0,0) → never a clearcut
end

@testset "PPE MXHRVP — policy expression eval (ALGCMP/ALGEVL via ported evmon)" begin
    hv = HarvestVars(; avbtpa=250f0, avbtcuft=3000f0, avbmcuft=2600f0, avbbdft=1000f0,
                       avbba=120f0, avbacc=40f0, avbmort=8f0, totalwt=40f0, oldtarg=100f0)
    # TARGET-style expressions over the PTSTV1 landscape aggregates (Float32-exact).
    @test eval_policy_expr("AVBBDFT * 0.5", hv)        == 1000f0 * 0.5f0
    @test eval_policy_expr("AVBBDFT / TOTALWT", hv)    == 1000f0 / 40f0
    @test eval_policy_expr("TOTALWT * OLDTARG - AVBTPA", hv) == 40f0 * 100f0 - 250f0
    @test eval_policy_expr("AVBTCUFT - AVBMCUFT", hv)  == 3000f0 - 2600f0
    # SELECTED (EVSET4 code 9) — 1 selected / 0 not.
    @test eval_policy_expr("SELECTED * 500", HarvestVars(; selected=1f0)) == 500f0
    @test eval_policy_expr("SELECTED * 500", HarvestVars(; selected=0f0)) == 0f0
    # PPE vars read 0.0 when there is no harvest context (matching zeroed PTSTV1).
    @test eval_event(parse_event_condition("AVBBDFT + 7"), EventCtx(1, 1990, nothing)) == 7f0
    # backward-compat: the 3-arg EventCtx still constructs (harvest = nothing).
    @test EventCtx(1, 1990, nothing).harvest === nothing
end

@testset "PPE MXHRVP — LBMEMR label-set membership (lbmemr.f/lb1mem.f) vs Fortran" begin
    # comma-space-delimited label sets; exact, length-gated membership (the candidacy test).
    @test lbmemr("ALL",    "ALL, STAND1, GROUP2") == true
    @test lbmemr("STAND1", "ALL, STAND1, GROUP2") == true
    @test lbmemr("GROUP2", "ALL, STAND1, GROUP2") == true    # last member (no trailing comma)
    @test lbmemr("X",      "ALL, STAND1")         == false
    @test lbmemr("STAND",  "ALL, STAND1")         == false   # prefix, different length → no
    @test lbmemr("ALL",    "ALL")                 == true    # single-member set
    @test lbmemr("B",      "A, B, C, D")          == true
    @test lbmemr("D",      "A, B, C, D")          == true
end
