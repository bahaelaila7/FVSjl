# =============================================================================
# test_thinprsc_fragment_d14.jl — THINPRSC residual-fragment deletion (D14) vs live FVSsn.
#
# A proportional THINPRSC (e.g. `THINPRSC 2000 0.999`, cuteff 0.999) leaves a 0.001-scale residual on
# each pre-marked record. FVS cuts.f:1631-1637 DELETES any cut record whose residual TPA (what's left) is
# ≤ 0.0005 — it cuts the ENTIRE tree (PROB→0) and TREDEL compacts it out. jl previously kept those tiny
# fragments as live records, so after the thin jl had 243 tree records vs live's 230 (13 extra fragments;
# BA/normalized-TPA ~bit-exact but the divergent record structure amplified to ~11% Scuft@2010 at the saw
# threshold). Fix: `_thinprsc!` now applies the residual≤0.0005 whole-tree deletion (cuts.jl).
#
# Golden = live FVSsn (cut_thinprsc.key = THINPRSC 2000 0.999 on S248112, with tripling): the fix makes the
# .sum BIT-EXACT through 2030 (TPA/BA/Scuft/Bdft), with only 1-2 unit ULP residuals at 2035/2040.
# =============================================================================

using Test, FVSjl

@testset "THINPRSC residual-fragment deletion (D14) vs live FVSsn" begin
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "cut_thinprsc.key")
    if !isfile(key)
        @test_skip "cut_thinprsc scenario not available"
    else
        txt = FVSjl.run_keyfile(key; variant = FVSjl.Southern(), output = :sum)
        rows = Dict{String,Vector{SubString{String}}}()
        for ln in split(txt, '\n')
            t = split(strip(ln))
            length(t) >= 12 && occursin(r"^(19|20)\d\d$", t[1]) && (rows[t[1]] = t)
        end
        # live FVSsn goldens (TPA, BA, Scuft, Bdft) — bit-exact through 2030 post-fix.
        live = Dict(
            "2000" => ("470", "126", "545", "2484"),
            "2005" => ("192", "124", "1108", "5095"),
            "2010" => ("177", "142", "1977", "9279"),
            "2020" => ("123", "150", "3168", "15919"),
            "2030" => ("86", "152", "3839", "21121"),
        )
        for (yr, (tpa, ba, sc, bd)) in live
            @test haskey(rows, yr)
            if haskey(rows, yr)
                r = rows[yr]
                @test r[3] == tpa      # TPA
                @test r[4] == ba       # BA
                @test r[11] == sc      # Scuft
                @test r[12] == bd      # Bdft
            end
        end
        # 2010 was the worst pre-fix divergence (jl 1759 vs live 1977 Scuft, ~11%) — now exact.
        @test rows["2010"][11] == "1977"
    end
end
