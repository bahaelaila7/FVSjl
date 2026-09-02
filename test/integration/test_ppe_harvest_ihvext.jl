# PPE MXHRVP IHVEXT=1 (external harvest selection) — HVSEL's external-selection rule
# bit-exact vs gfortran-16 driver-golden (scratchpad/ppe/mxhrvp/driver_hvsel_ext.f, the real
# hvsel.f with IHVEXT=1 over the historical FVSppe objects). Under external selection the
# priorities are supplied by an external program (staged read); HVSEL never selects a stand
# whose supplied priority is ≤ 0, and the EXACT partial-cut is forced off.
using Test
using FVSjl: hvsel!

# (label, target, [(priority, yield_sel, yield_notsel)...], expected_status)
const IHVEXT_GOLDENS = [
    ("X1 pri 0 and -3 skipped", 150.0f0,
        [(10f0,100f0,5f0),(0f0,100f0,5f0),(-3f0,100f0,5f0),(8f0,100f0,5f0)],
        [3,2,2,2]),
    ("X2 pri -1 skipped, 5&9 selected", 250.0f0,
        [(5f0,100f0,5f0),(-1f0,100f0,5f0),(9f0,100f0,5f0)],
        [3,2,3]),
]

@testset "PPE MXHRVP IHVEXT — external-selection rule bit-exact vs driver-golden" begin
    for (label, target, stands, exp_status) in IHVEXT_GOLDENS
        n = length(stands)
        status = fill(1, n)
        prio = Float32[s[1] for s in stands]
        ysel = Float32[s[2] for s in stands]
        ynot = Float32[s[3] for s in stands]
        st, = hvsel!(status, prio, ysel, ynot, target; ihvext = true)
        @testset "$label" begin
            @test st == exp_status
        end
    end

    # IHVEXT forces the EXACT partial-cut off: even with lprtct=true requested, a marginal
    # stand rounds (>0.5 select / ≤0.5 not) rather than becoming a status-4 partial.
    st = fill(1, 2)
    _, hvpart = hvsel!(st, Float32[9,6], Float32[100,100], Float32[0,0], 40f0;
                       lprtct = true, ihvext = true)
    @test hvpart == 0.0f0   # no partial cut under external selection
end
