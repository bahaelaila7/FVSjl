# test_translate.jl — the YML/CSV ↔ .key/.tre conversion tool (bin/fvsjl-translate.jl
# → translate_io). The DROP-IN-relevant directions (the engine reads the modern forms
# directly) must be lossless:
#   * .key → .yaml → .key  reproduces the same engine output (.sum)
#   * .tre → .csv          preserves the parsed TreeRecords (engine reads either)
# The legacy re-emission .csv → .tre currently has a write_tree_file bug on overlapping
# T-specifier fields (an F-field sharing a column with a packed nI1 field loses its value);
# it does not affect the drop-in (the engine reads the .csv), and is tracked @test_broken.

using Test, FVSjl

const _TR_DIR = joinpath(@__DIR__, "..", "keyword_coverage", "scenarios")
const _TR_SNFMT = "(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,T52,I2,T66,5I1,T54,7I1,T75,F3.0)"
const _TR_SCEN = ["s1_thinning", "s13_thinht", "s25_thinrden", "s31_mcdefect", "s32_volume",
                  "s34_serlcorr", "s17_managed", "s37_thinauto", "s7_eventmon", "s26_estab"]

@testset "YML/CSV ↔ .key/.tre conversion tool" begin
    if !isdir(_TR_DIR)
        @test_skip "coverage scenarios not available"
    else
        for name in _TR_SCEN
            key = joinpath(_TR_DIR, name * ".key"); tre = joinpath(_TR_DIR, name * ".tre")
            isfile(key) || continue
            @testset "$name" begin
                d = mktempdir()
                a_key = joinpath(d, "a.key"); cp(key, a_key; force = true)
                isfile(tre) && cp(tre, joinpath(d, "a.tre"); force = true)
                # (1) .key → .yaml → .key is engine-equal (run from d so TREEDATA finds a.tre)
                FVSjl.translate_io(a_key, joinpath(d, "a.yaml"))
                FVSjl.translate_io(joinpath(d, "a.yaml"), joinpath(d, "b.key"))
                isfile(tre) && cp(tre, joinpath(d, "b.tre"); force = true)
                cd(d) do
                    # compare data rows — the -999 header now carries a wall-clock timestamp
                    nohdr(t) = join([l for l in split(t, "\n") if !startswith(l, "-999")], "\n")
                    @test nohdr(FVSjl.run_keyfile("a.key")) == nohdr(FVSjl.run_keyfile("b.key"))
                end
                if isfile(tre)
                    a_tre = joinpath(d, "a.tre"); a_csv = joinpath(d, "a.csv")
                    recs = FVSjl.read_tree_records(a_tre; fmt = _TR_SNFMT)
                    # (2) .tre → .csv preserves the TreeRecords (the engine reads the .csv)
                    FVSjl.convert_tre_to_csv(a_tre, a_csv; fmt = _TR_SNFMT)
                    @test FVSjl.read_trees_csv(a_csv) == recs
                    # (3) legacy re-emission .csv → .tre round-trips. write_tree_file now right-justifies
                    # an F-field that overflows a column it shares with a packed nI1 field to its own col2
                    # (was spilling left under the I1 — _TR_SNFMT's T54,7I1 / T60,F3.1 boundary at col 60).
                    c_tre = joinpath(d, "c.tre")
                    FVSjl.convert_csv_to_tre(a_csv, c_tre; fmt = _TR_SNFMT)
                    @test FVSjl.read_tree_records(c_tre; fmt = _TR_SNFMT) == recs
                end
            end
        end
    end
end
