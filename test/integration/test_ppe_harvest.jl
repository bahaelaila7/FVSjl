# PPE MXHRVP multistand harvest scheduling — HVSEL bit-exact vs gfortran-16 driver-golden.
# Goldens produced by scratchpad/ppe/mxhrvp/driver_hvsel.f (the real pristine hvsel.f linked
# over the historical FVSppe objects, bc6e2377^): each case feeds synthetic policy-1 inputs
# (priority, yield-if-selected, yield-if-not, target, EXACT flag) and dumps IHVSTA + HVPART.
using Test
using FVSjl: hvsel!

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
