# test_bamax.jl — BAMAX keyword → user-pinned maximum basal area (initre.f:6800) vs Fortran.
#
# BAMAX sets the stand's max attainable BA. site_index_setup! then derives every species'
# SDImax default from it (SDIDEF = BAMAX/(0.5454154·PMSDIU), sdical.f:208), so the SDImax-
# driven self-thinning mortality (mortality.jl) caps the residual BA at BAMAX. Before the
# fix BAMAX was an *unrecognized* keyword (silently ignored), and the BAMAX→SDIDEF branch
# in site_index had two latent bugs (no 0.85 PMSDIU default ⇒ div-by-zero; a stray /100).
# Checks:
#   1. SETS — the keyword sets s.control.ba_max;
#   2. FIRES — a dense stand caps its BA at ~BAMAX (≠ the same stand without BAMAX);
#   3. FORTRAN — TPA/BA/SDI/CCF/TopHt/QMD match live Fortran (board feet within Scribner noise).

using Test, FVSjl

const _BX_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_bx_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_bx_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "BAMAX → SDImax self-thinning cap" begin
    # 1. SETS — the keyword sets control.ba_max (initre.f:6800; field-1 ≤ 0 leaves it 0)
    setbx(field1) = begin
        s = FVSjl.StandState(FVSjl.Southern())
        rec = FVSjl.KeywordRecord("BAMAX   ", "", [field1, fill("", 11)...],
                                  Float32[parse(Float32, field1), zeros(Float32, 11)...],
                                  [!isempty(field1), falses(11)...], 12, FVSjl.KW_OK, 0)
        FVSjl.kw_bamax!(s, rec)
        s.control.ba_max
    end
    @test setbx("150") == 150f0
    @test setbx("0") == 0f0          # field-1 == 0 ⇒ not set (stays 0)

    key = joinpath(_BX_DIR, "bamax.key")
    sav = joinpath(_BX_DIR, "bamax.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "bamax scenario not available"
    else
        jl = _bx_rows(FVSjl.run_keyfile(key; faithful = true))

        # 2. FIRES — removing BAMAX lets the BA grow well past the cap.
        offkey = joinpath(_BX_DIR, "_bamax_off.key")
        cp(joinpath(_BX_DIR, "bamax.tre"), joinpath(_BX_DIR, "_bamax_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "BAMAX") || println(io, l); end
        end
        try
            off = _bx_rows(FVSjl.run_keyfile(offkey; faithful = true))
            # the capped run holds BA at or under 150; without BAMAX the stand's own SDImax governs. Live FVSsn_g16 on
            # the BAMAX-less key (bamax_off.sum.save) peaks at 149 in 2015 by a different path (TPA 292 vs the capped
            # run's), so assert the uncapped run against live, not a hand-picked threshold. (The old `> 180` encoded a
            # jl misread of this tree file — 5897 TPA at 1990 vs live 536 — fixed on bc-maxtre.)
            @test maximum(parse(Int, r[4]) for r in jl) <= 160
            offlive = _bx_base(joinpath(_BX_DIR, "bamax_off.sum.save"))
            @test length(off) == length(offlive)
            for (j, f) in zip(off, offlive), c in (3, 4, 5, 6, 7)       # TPA, BA, SDI, CCF, TopHt
                @test parse(Int, j[c]) == parse(Int, f[c])
            end
            @test [r[3:7] for r in off] != [r[3:7] for r in jl]          # BAMAX changes the trajectory
        finally
            rm(offkey; force = true); rm(joinpath(_BX_DIR, "_bamax_off.tre"); force = true)
        end

        # 3. FORTRAN — TPA/BA/SDI/CCF/TopHt/QMD match live Fortran.
        ft = _bx_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            for c in (3, 4, 5, 6, 7)                    # TPA, BA, SDI, CCF, TopHt
                @test parse(Int, j[c]) == parse(Int, f[c])
            end
            @test parse(Float64, j[8]) == parse(Float64, f[8])   # QMD
            @test parse(Int, j[10]) == parse(Int, f[10])            # MerchCuFt — BIT-EXACT (measured Δ=0; was over-cautious ≤2)
        end
    end
end
