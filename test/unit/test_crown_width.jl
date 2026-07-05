# C2 unit tests — crown-width equation library (CSV-driven, engine/crown_width.jl).
# Anchor values were captured from the original baked-in 144-equation library, one
# per formula family, to guard the CSV data + family evaluator against regressions.

using Test
using FVSjl
using FVSjl: crown_width, coefficients, Southern, hopkins_index

@testset "crown_width families (value-identity vs baked library)" begin
    c = coefficients(Southern())
    lat, long, elev = 34.5f0, 86.0f0, 5f0
    hi = hopkins_index(lat, long, elev)

    # All four families evaluate the formula in the SAME Float32 op-order as jl's crown_width, so the
    # identity is EXACT (measured Δ=0 for every case incl. the ^power families — Julia's Float32 `^` is
    # shared). Compare == (was padded atol 1f-4 / 1f-3 that hid a bit-identical evaluation).
    # bechtold M2 (loblolly LP=13101): a + b·D + cr·CR, clamp 55
    @test crown_width(c, "LP", 12f0, 50f0, 50f0, 0, lat, long, elev) ==
          (-0.8277f0 + 1.3946f0*12f0 + 0.0768f0*50f0)
    # bechtold M3 with D² + dcap (red maple RM=31601, dlim=min(D,50))
    @test crown_width(c, "RM", 20f0, 40f0, 60f0, 0, lat, long, elev) ==
          (2.7563f0 + 1.4212f0*20f0 - 0.0143f0*400f0 + 0.0993f0*60f0 - 0.0276f0*hi)
    # ek (power, omind floor, clamp) — American elm open form (97203)
    @test crown_width(c, "AE", 15f0, 50f0, 50f0, 1, lat, long, elev) ==
          (2.8290f0 + 3.4560f0*15f0^0.8575f0)
    # bragg family is value-checked directly (pure power, no scale/clamp)
    @test FVSjl._cw_eval(c.crown_eqs["97202"], 15f0, 50f0, hi) ==
          (-53.239079f0 + 61.327257f0*15f0^0.060166f0)
    # small-tree floor (D<5 ⇒ scaled toward 0) stays finite + clamped ≥ 0.5
    @test crown_width(c, "LP", 0.5f0, 50f0, 50f0, 0, lat, long, elev) ≥ 0.5f0
    # unknown species → 0.5
    @test crown_width(c, "ZZ", 10f0, 50f0, 50f0, 0, lat, long, elev) == 0.5f0
    # data shape
    @test length(c.crown_eqs) == 145
    @test length(c.crown_species) == 167
end
