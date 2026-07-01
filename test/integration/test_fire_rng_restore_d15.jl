# =============================================================================
# test_fire_rng_restore_d15.jl — FMEFF RANNGET/RANNPUT RNG save-restore (D15) vs live FVSsn.
#
# FVS's FMEFF brackets its per-tree fire-mortality RANN draws with RANNGET(SAVESO) (fmeff.f:143) …
# RANNPUT(SAVESO) (fmeff.f:569) — the fire's ~ITRN draws are ROLLED BACK, so a fire consumes ZERO net
# main-stream RNG. jl previously drew `rann!` per record WITHOUT restoring, advancing the stream ~ITRN
# draws and desyncing the POST-fire DGSCOR serial-correlation deviates: the fire KILL stayed bit-exact
# (same draws) but the survivors then grew wrong, ~4.4% Bdft high by the 3rd post-fire cycle. Fix: fmburn!
# now saves (rannget) before the fire loop and restores (rannput!) after — matching FVS.
#
# (NOTE: my first D15 diagnosis — the fmeff.f FMICR crown-scorch reducing survivor crown ratio — was WRONG:
# fmmain.f:111 `FMICR=ICR` shows FMICR is FFE-internal (fuel/potential-fire), it does NOT feed the growth
# crown. The re-trace found the real cause = the RANNGET/RANNPUT rollback.)
#
# Golden = live FVSsn (fire_burn.key = SIMFIRE 2000 on S248112): post-fix the .sum is BIT-EXACT on
# TPA/BA/QMD every cycle, with only 1-13 unit ULP Bdft residuals at the far post-fire cycles.
# =============================================================================

using Test, FVSjl

@testset "fire FMEFF RANN save-restore (D15) vs live FVSsn" begin
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "fire_burn.key")
    if !isfile(key)
        @test_skip "fire_burn scenario not available"
    else
        txt = FVSjl.run_keyfile(key; variant = FVSjl.Southern(), output = :sum)
        rows = Dict{String,Vector{SubString{String}}}()
        for ln in split(txt, '\n')
            t = split(strip(ln))
            length(t) >= 12 && occursin(r"^(19|20)\d\d$", t[1]) && (rows[t[1]] = t)
        end
        # live FVSsn goldens (TPA, BA, QMD) — bit-exact post-fix (were BA 85/84, QMD 12.4/12.3 pre-fix).
        live = Dict(
            "2000" => ("470", "126", "7.0"),   # pre-fire
            "2005" => ("104", "70", "11.1"),   # post-fire (kill was already bit-exact)
            "2010" => ("101", "84", "12.3"),   # 1st post-fire growth — the pre-fix divergence (was BA 85/12.4)
            "2015" => ("99", "97", "13.4"),    # 2nd — was BA 100/13.6
        )
        for (yr, (tpa, ba, qmd)) in live
            @test haskey(rows, yr)
            if haskey(rows, yr)
                r = rows[yr]
                @test r[3] == tpa   # TPA
                @test r[4] == ba    # BA  (bit-exact post-fire — the headline of the fix)
                @test r[8] == qmd   # QMD
            end
        end
    end
end
