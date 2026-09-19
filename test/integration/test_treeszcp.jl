# test_treeszcp.jl — per-species size cap (TREESZCP / SIZCAP) vs live Fortran.
#
# TREESZCP sets a per-species maximum size: field 1 = species (0 = all), 2 = cap DBH
# (SIZCAP[1]), 3 = annual mortality rate of capped trees (SIZCAP[2]), 4 = no-mortality
# flag IDMFLG (SIZCAP[3]), 5 = height cap (SIZCAP[4]). The cap drives three mechanisms,
# each validated by a scenario whose .sum.save is live-Fortran output for the same key:
#
#   * treeszcp_nomort — DBH cap 10" with IDMFLG=1 (no size-cap mortality): exercises ONLY
#     the diameter-growth bound (dgbnd). Bit-exact every cycle (TPA/BA/TopHt/QMD).
#   * treeszcp_cap    — DBH cap 10" with mortRate 1.0: DG bound + size-cap mortality floor
#     (morts.f:692). QMD is bit-exact every cycle and the endpoint matches; the mid-cycle
#     TPA/BA carry the regen response to the cap-driven mortality (the known regen tail).
#   * treeszcp_htcap  — height cap 30': exercises the htgf.f:286 HT cap. TPA/BA/QMD are
#     bit-exact every cycle; TopHt drifts ≤4' as a declining-stand artifact (the frozen
#     tall trees fall out by mortality slightly faster than Fortran — height→crown→mort).

using Test, FVSjl

const _TSZ_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_tsz_rows(txt) = [split(l) for l in split(txt, "\n")
                  if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_tsz_base(path) = [split(l) for l in eachline(path)
                   if length(split(l)) >= 11 &&
                      (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_col(r, c) = parse(Float64, r[c])

@testset "size cap (TREESZCP / SIZCAP) vs Fortran" begin
    have(nm) = isfile(joinpath(_TSZ_DIR, nm * ".key")) && isfile(joinpath(_TSZ_DIR, nm * ".sum.save"))
    runjl(nm) = (_tsz_rows(FVSjl.run_keyfile(joinpath(_TSZ_DIR, nm * ".key"); faithful = true)),
                 _tsz_base(joinpath(_TSZ_DIR, nm * ".sum.save")))

    # All three scenarios are BIT-EXACT every cycle on TPA / BA / TopHt / QMD (cols 3 4 7 8) vs live FVSsn
    # (re-derived 2026-09-19 on FVSsn_g16 and the corrected FVSsn_g16.new — identical .sum). The earlier
    # "cornered" residuals were real bugs, traced per record with a hex-dumping instrumented FVSsn copy:
    #   * VARMRT `(1-EFFTR)**NPASS` is REAL**INTEGER = libgcc __powisf2 (single-precision square-and-multiply);
    #     Julia's Float32^Int rounds once from Float64 ⇒ 1-ULP power ⇒ ~1e-5 relative kill error through the
    #     1-(…) cancellation (fpowi).
    #   * PCTILE/DENSE: WK5 = D*(D*P) and PCT = cum/(TOT/100.), top record = 100 — jl used (D*D)*P and ×100/TOT.
    #   * The SIZCAP kill P*SIZCAP*FINT/5.0 is 1 ULP BELOW P on some records in FVS (left-to-right Float32), so
    #     a capped record survives with ~2e-7 TPA and is deleted by COMCUP (≤1E-5), not the CUTS-entry TREDEL
    #     (≤1E-10) — the record-order effect that made the endpoint shift once deletion timing was made faithful.
    #   * htgf.f:292-307 caps each TRIPLED copy's uncapped TEMHTG against HT(ITFN) of its future slot BEFORE
    #     TRIPLE writes it ⇒ FVS reads what TREDEL left in that static COMMON slot (deterministic, not UB): jl
    #     keeps that shadow in t.stale_ht.
    for nm in ("treeszcp_nomort", "treeszcp_cap", "treeszcp_htcap")
        if !have(nm); @test_skip "$nm scenario not available"; continue; end
        @testset "$nm" begin
            jl, ft = runjl(nm)
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl), c in (3, 4, 7, 8)
                    @test _col(jl[i], c) == _col(ft[i], c)
                end
            end
        end
    end

    # The cap must visibly act, and the endpoint must land where Fortran does.
    if have("treeszcp_cap")
        jl, ft = runjl("treeszcp_cap")
        @test _col(jl[end], 8) <= 8                 # QMD capped near the 10" DBH limit (base ≈ 15)
        @test _col(jl[end], 3) == _col(ft[end], 3)            # endpoint TPA — bit-exact (was the "DGF ULP×cap" corner)
        @test _col(jl[end], 4) == _col(ft[end], 4)            # endpoint BA — bit-exact
    end
    if have("treeszcp_htcap")
        jl, ft = runjl("treeszcp_htcap")
        @test _col(jl[end], 7) <= 50                # TopHt held well below the uncapped ≈ 79
        # TopHt every cycle — bit-exact via the stale-slot HT(ITFN) cap (htgf.f:292-307; t.stale_ht).
        @test all(_col(jl[i], 7) == _col(ft[i], 7) for i in 1:length(jl))
    end
end
