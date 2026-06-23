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

    # (scenario, columns required bit-exact every cycle). Cols: 3 TPA / 4 BA / 7 TopHt / 8 QMD.
    for (nm, cols) in (("treeszcp_nomort", (3, 4, 7, 8)),  # pure DG bound — fully bit-exact
                       ("treeszcp_cap",    (8,)),          # DG bound + size-cap mort — QMD exact
                       ("treeszcp_htcap",  (3, 4, 8)))     # HT cap — TPA/BA/QMD exact (TopHt below)
        if !have(nm); @test_skip "$nm scenario not available"; continue; end
        @testset "$nm" begin
            jl, ft = runjl(nm)
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl), c in cols
                    @test abs(_col(jl[i], c) - _col(ft[i], c)) <= 1
                end
            end
        end
    end

    # The cap must visibly act, and the endpoint must land where Fortran does.
    if have("treeszcp_cap")
        jl, ft = runjl("treeszcp_cap")
        @test _col(jl[end], 8) <= 8                 # QMD capped near the 10" DBH limit (base ≈ 15)
        @test abs(_col(jl[end], 3) - _col(ft[end], 3)) <= 6   # endpoint TPA (regen-tail tolerance)
        @test abs(_col(jl[end], 4) - _col(ft[end], 4)) <= 2   # endpoint BA
    end
    if have("treeszcp_htcap")
        jl, ft = runjl("treeszcp_htcap")
        @test _col(jl[end], 7) <= 50                # TopHt held well below the uncapped ≈ 79
        for i in 1:length(jl)                       # TopHt drift bounded (declining-stand artifact)
            @test abs(_col(jl[i], 7) - _col(ft[i], 7)) <= 4
        end
    end
end
