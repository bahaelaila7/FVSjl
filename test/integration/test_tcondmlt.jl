# test_tcondmlt.jl — TCONDMLT tree-condition cut weights (cuts.f:1074/1424) vs live Fortran.
#
# TCONDMLT sets two weights added to the RDPSRT cut-priority key:
#   WK2 = ±DBH + IORDER(SPECPREF) + TCWT·IMC + SPCLWT·ISPECL
# TCWT (PRM 1) weights the mortality/condition code IMC (1..3 live); SPCLWT (PRM 2) weights the
# special-status code ISPECL (damage code 55). A positive weight removes the flagged trees first,
# regardless of size. Two scenarios, each marking the sugar maples and thinning from below (THINBBA):
#   * tcondmlt — IMC=3 + TCWT=100 (condition weight);
#   * spclwt   — special-status code 55 (ISPECL=9) + SPCLWT=100 (special-status weight).
# Each must (1) differ from the same thin without the TCONDMLT line, and (2) match live Fortran on
# TPA/BA/cubic columns (board feet within Scribner noise).

using Test, FVSjl

const _TC_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_tc_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_tc_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_tccol(r, c) = parse(Float64, r[c])

@testset "TCONDMLT condition-weighted thin vs Fortran" begin
    for stem in ("tcondmlt", "spclwt")
        key = joinpath(_TC_DIR, stem * ".key"); sav = joinpath(_TC_DIR, stem * ".sum.save")
        tre = joinpath(_TC_DIR, stem * ".tre")
        if !isfile(key) || !isfile(sav)
            @test_skip "$stem scenario not available"; continue
        end
        @testset "$stem" begin
            jl = _tc_rows(FVSjl.run_keyfile(key; faithful = true))

            # 1. NON-VACUOUS: without TCONDMLT the thin ranks by size only → a different cut.
            notc = tempname() * ".key"
            write(notc, join(filter(l -> !startswith(l, "TCONDMLT"), readlines(key)), "\n") * "\n")
            cp(tre, replace(notc, ".key" => ".tre"); force = true)
            @test jl != _tc_rows(FVSjl.run_keyfile(notc; faithful = true))

            # 2. matches live Fortran (cubic columns bit-exact; board feet within Scribner noise).
            ft = _tc_base(sav)
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl)
                    @test _tccol(jl[i], 3) == _tccol(ft[i], 3)   # TPA — BIT-EXACT
                    @test _tccol(jl[i], 4) == _tccol(ft[i], 4)   # BA  — BIT-EXACT
                    for c in (10, 11)                           # MCuFt/SCuFt — BIT-EXACT (measured Δ0 both scenarios)
                        @test _tccol(jl[i], c) == _tccol(ft[i], c)
                    end
                    # TCuFt + BdFt: measured — the tcondmlt stem is BIT-EXACT (Δ0 both cols, all cycles); the
                    # spclwt stem carries a genuine 1-step rendered-integer residual (exposed below the loop).
                    if stem != "spclwt"
                        @test _tccol(jl[i], 9)  == _tccol(ft[i], 9)            # TCuFt — BIT-EXACT (tcondmlt stem, measured Δ0)
                        @test _tccol(jl[i], 12) == _tccol(ft[i], 12)           # BdFt  — BIT-EXACT (tcondmlt stem, measured Δ0)
                    end
                end
                # spclwt stem: the TCuFt/BdFt 1-step render residual is the non-associative Float32 tree-SUM
                # accumulation order (doctrine #9: exposed as @test_broken, not a passing ≤1 hiding in green).
                if stem == "spclwt"
                    @test_broken all(_tccol(jl[i], 9)  == _tccol(ft[i], 9)  for i in 1:length(jl))  # TCuFt — non-associative tree-SUM order
                    @test_broken all(_tccol(jl[i], 12) == _tccol(ft[i], 12) for i in 1:length(jl))  # BdFt  — non-associative tree-SUM order
                end
            end
        end
    end
end

@testset "TCONDMLT point weights (PBAWT/PCCFWT/PTPAWT) are inert — matches live FVS" begin
    # The TCONDMLT point-density terms (cuts.f:1075 +PBAWT·PTBAA+PCCFWT·PCCF+PTPAWT·PTPA) are EMPIRICALLY
    # INERT in live FVSsn: a multi-point (11-point) THINBBA with PTPAWT set produces a .sum byte-identical to
    # the same thin without any TCONDMLT, even at PTPAWT=9999. FVS's cuts.f PTPA(IP)/PTBAA(IP) is uniform/zero
    # across points at thin time (a point-thinning-only path plain TCONDMLT does not arm). jl therefore does NOT
    # add the point term (`_cut_pref_wt`); doing so diverged from live. This guards against re-introducing it:
    # tcond_pw (PTPAWT=1) must equal tcond_base (no TCONDMLT) AND match the committed live oracle.
    base = joinpath(_TC_DIR, "tcond_base.key"); pw = joinpath(_TC_DIR, "tcond_pw.key")
    sav  = joinpath(_TC_DIR, "tcond_pw.sum.save")
    if !isfile(pw) || !isfile(base) || !isfile(sav)
        @test_skip "tcond point-weight scenarios not available"
    else
        jl_pw   = _tc_rows(FVSjl.run_keyfile(pw;   faithful = true))
        jl_base = _tc_rows(FVSjl.run_keyfile(base; faithful = true))
        @test jl_pw == jl_base                      # point weight changes nothing (inert, as in live)
        ft = _tc_base(sav)
        @test length(jl_pw) == length(ft)
        for (j, f) in zip(jl_pw, ft)
            for c in (3, 4)                          # TPA, BA bit-exact vs the live point-weighted run
                @test _tccol(j, c) == _tccol(f, c)
            end
        end
    end
end
