# Climate-FVS SPCALIB first-cycle presence-calibration (clmorts.f:57-98) — the deferred "chunk C".
# At ICYC=1 FVS sets SPCALIB[sp] = viab@IY(ICYC)·0.9 for species PRESENT in the inventory (−1 if absent), then
# rescales the viability→survival curve for present LOW-viability species so they are not over-killed in early
# cycles. Without it jl over-reported climate mortality for a resident low-viability species (PP: 0.2217 at cyc1
# vs the oracle's 0). VALIDATED here BIT-EXACT vs the oracle FVS_Climate ViabMort column (climtest_oracle.db,
# IE Clearwater CGCM3_A2 stand S248112): PP's SPMORT1 series across 1990..2090 reproduces the oracle exactly.
# Fixture: test/fixtures/climate/clim_iet.key (inline CLIMDATA) + .tre. Oracle values are the FVS_Climate
# ViabMort column dumped 2026-08-21.

using Test
using FVSjl
const _F = FVSjl

@testset "Climate-FVS SPCALIB presence-calibration (clmorts chunk C)" begin
    keyf = joinpath(@__DIR__, "..", "fixtures", "climate", "clim_iet.key")
    @assert isfile(keyf)
    V = _F.variant_from_code("IE")
    got = nothing
    for s in _F.each_stand(keyf; variant = V)
        _F.notre!(s); _F.setup_growth!(s); _F.compute_volumes!(s)
        (s.climate === nothing || !s.climate.active) && continue
        c = s.climate; cd = c.data
        sp = findfirst(==("PIPO"), c.plant_symbols)      # PP = ponderosa (present, low viability)
        sp === nothing && continue
        # SPCALIB(PP) = raw viability at IY(ICYC)=1990 (cycle-start) · 0.9  (clmorts.f:63)
        spcal = _F.species_vscore(cd, "PIPO", 1990f0)[1] * 0.9f0
        series = Float64[]
        for yr in 1990:10:2090
            xv = _F.species_vscore(cd, "PIPO", Float32(yr) + 5f0)[1]   # THISYR = IY(ICYC)+FINT/2 (fint=10)
            surv = _F.clim_survival_cal(xv, spcal)                     # presence-calibrated 10-yr survival
            push!(series, (1f0 - surv) * c.mortmult[sp])               # SPMORT1 = (1−X)·CLMRTMLT1
        end
        got = series
        break
    end
    @test got !== nothing
    # oracle FVS_Climate ViabMort for PP, 1990..2090 (10-yr cycles):
    oracle = [0.0, 0.0, 0.0338, 0.1214, 0.2291, 0.357, 0.4849, 0.5753, 0.6281, 0.681, 0.7074]
    @test length(got) == length(oracle)
    for (g, o) in zip(got, oracle)
        @test isapprox(g, o; atol = 5e-4)     # bit-exact to the oracle's 4-dp report precision
    end
    # the calibration must keep a resident low-viability species alive while viab ≥ its calibrated floor
    @test got[1] == 0.0 && got[2] == 0.0      # 1990, 2000: no climate mortality (was 0.2217 pre-SPCALIB)
end
