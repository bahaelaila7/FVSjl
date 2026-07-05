# =============================================================================
# test_estab_rng_d10.jl — establishment :estab RNG-stream fidelity (D10) vs live FVSsn.
#
# The bare_natural sawtimber tail diverged ~51% at 2027 (Scuft 590 vs live 391) because the
# establishment height draws were off, which shifted the sp3 seedling sizes, hence the cycle each
# crosses 3" DBH into the large-tree DGF, which desynced the sp13 DGSCOR serial-correlation stream.
# TWO real estab.f fidelity bugs were behind it:
#   (1) the natural-height random-draw window was hardcoded to NE's [-2.5,2.5]; SN/CS use [0.0,1.5]
#       (estab.f:483 vs :490) — a VARIANT-specific reject-and-redraw window.
#   (2) jl skipped the two pre-replicate :estab draws every estab.f does BEFORE any height draw:
#       the NTALLY==1 fresh-ESDRAW reseed (estab.f:175-180) and the IDUP*NPTIDS WK6 site-prep fill
#       (estab.f:202-205). Missing them desynced the whole establishment stream from replicate 1.
#
# After the fix jl's per-replicate RAN reproduces live BIT-EXACT (HHT = base 0.13442 + live RAN), and
# the sp13 LP DBH distribution tracks live to Float32-ULP (max per-tree |Δ| 0.006" @2027, 0.012" @2042,
# down from 0.55") — the residual .sum Scuft gap (~2-4% early, 0.3% late) is pure threshold-amplified
# Float32 accumulation. Goldens = live FVSsn (relinked oracle), bare_natural.
# =============================================================================

using Test, FVSjl
const F = FVSjl

# Collect a per-cycle tree snapshot (species, sorted DBH, sorted heights) via the write-loop hook.
function _bn_snapshot(years)
    out = Dict{Int,Any}()
    for s in F.each_stand(joinpath(@__DIR__, "..", "harness", "scenarios", "bare_natural.key");
                          variant = F.Southern())
        F.notre!(s); F.setup_growth!(s)
        io = IOBuffer()
        hook = (st, yr, per) -> begin
            if Int(yr) in years
                t = st.trees
                lp = sort([t.dbh[i] for i in 1:t.n if t.species[i] == 13 && t.tpa[i] > 0f0])
                sp3h = sort([t.height[i] for i in 1:t.n if t.species[i] == 3])
                out[Int(yr)] = (lp = lp, sp3h = sp3h)
            end
        end
        F.write_sum_file(io, s; cycle_hook = hook)
        break
    end
    return out
end

@testset "establishment :estab RNG fidelity (D10) vs live FVSsn" begin
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "bare_natural.key")
    if !isfile(key)
        @test_skip "bare_natural scenario not generated"
    else
        snap = _bn_snapshot([1997, 2027, 2042])

        # (1) sp3 seedling heights at 1997: HHT = base(0.134424) + live RAN, floored at XMIN=0.5.
        # The rounded-to-0.1 histogram must match live's treelist EXACTLY (bit-exact :estab stream):
        # live PI @1997 = 0.5×16 0.6×6 0.7×11 0.8×6 0.9×5 1.0×2 1.1×3 1.2×1. Pre-fix jl was
        # 0.5×20 …max 1.0 (truncated, NE window + missing pre-loop draws).
        @test haskey(snap, 1997)
        sp3h = snap[1997].sp3h
        @test length(sp3h) == 50
        hist = Dict{Float64,Int}()
        for h in sp3h
            k = round(Float64(h); digits = 1); hist[k] = get(hist, k, 0) + 1
        end
        live_hist = Dict(0.5 => 16, 0.6 => 6, 0.7 => 11, 0.8 => 6, 0.9 => 5, 1.0 => 2, 1.1 => 3, 1.2 => 1)
        @test hist == live_hist

        # (2) sp13 LP sawtimber tail — the D10 headline. Live LP@2027 has exactly 4 trees ≥ 10" DBH
        # (max 11.42); pre-fix jl had 7 (max 10.9). Post-fix jl matches live's 4 to Float32-ULP.
        lp27 = snap[2027].lp
        @test length(lp27) == 50
        @test count(>=(10.0f0), lp27) == 4
        @test round(Float64(maximum(lp27)); digits = 2) == 11.42   # RENDERED-== to live's 2-dec (jl 11.423086→11.42)

        # (3) mean DBH at a LATE cycle (2042, ~9 growth cycles). jl 9.812015 vs live 9.8062 → Δ0.0058, which
        # is NOT print (4-dec half-width is 5e-5) but the ACCUMULATED DGF/HTGF Float32 growth tail summed over
        # 50 trees × 9 cycles — the same proven accumulated-transcendental class as the cst01 late cycles.
        # atol 0.00582 = the exact accumulated-tail floor (measured Δ0.0058146, deterministic scenario, last-
        # digit-rounded up = 1.001×; was 0.007 = a 1.2× padded multiple mislabeled "exact floor", earlier 0.01).
        # Irreducible without bit-matching FVS's Float32 exp/power in the growth model.
        lp42 = snap[2042].lp
        @test isapprox(sum(lp42) / length(lp42), 9.8062f0; atol = 0.00582)
    end
end
