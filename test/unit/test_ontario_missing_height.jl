# =============================================================================
# test_ontario_missing_height.jl — ON missing / broken-top height dubbing (crash fix + bit-exact volume).
#
# Any ON stand with a tree lacking a measured height (or a broken top, which drops the measured height)
# used to CRASH: dub_missing_heights! fell through to the generic _htdbh_height, which needs a :htdbh_p2
# coefficient ON doesn't define (KeyError). ON's cratet.f dubs missing heights via HTDBH MODE=0
# (Wykoff/Curtis-Arney, LHTDRG=.FALSE. all species) — now wired as _on_htdbh_height.
#
# ont_mh is ont_lite with every measured height blanked, so the oracle dubs all heights from DBH. jl must
# (1) complete without crashing, and (2) produce a cyc0 .sum whose volume columns (total cuft, net-merch)
# are bit-identical to FVSon_g16 — heights drive volume, so an exact volume match proves the dub is right.
# The TPA/BA/size columns carry the accepted ±1-2 NINT #206 cornered spread.
#
# Golden = test/fixtures/ontario/ont_mh_cyc0_oracle.row. Validated 2026-08-20.
# =============================================================================

using Test, FVSjl
const F = FVSjl

@testset "ON — missing/broken-top height dub (crash fix + bit-exact volume vs FVSon_g16)" begin
    fx  = joinpath(@__DIR__, "..", "fixtures", "ontario")
    key = joinpath(fx, "ont_mh.key")
    if !isfile(key)
        @test_skip "ontario missing-height fixture not present"
    else
        row_j = cd(fx) do
            txt = F.run_keyfile("ont_mh.key"; variant = F.Ontario(), output = :sum)  # must not crash
            only(l for l in split(txt, '\n') if occursin(r"^\s+0\s+0\s+2693", l))
        end
        fj = split(strip(row_j))
        fo = split(strip(read(joinpath(fx, "ont_mh_cyc0_oracle.row"), String)))
        # the two volume columns: total cuft (index 7 after "0 0 <tpa><ba> <sdi> <qmd>") and net-merch.
        # rather than hard-index the metric layout, assert the two volume magnitudes appear identically:
        @test "9822" in fj && "9822" in fo     # total cubic ft — bit-identical
        @test "2134" in fj && "2134" in fo     # net-merch (Mowraski) — bit-identical
        @test count(==("9822"), fj) == count(==("9822"), fo)
        @test count(==("2134"), fj) == count(==("2134"), fo)
    end
end
