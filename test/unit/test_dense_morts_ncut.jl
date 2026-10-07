# test_dense_morts_ncut.jl — NC/UT morts.f on a synthetic >35,000-TPA seedling stand vs the LIVE oracle (FVSnc_g16 /
# FVSut_g16). Two morts.f mechanisms jl had wrong in klamath/utah mortality.jl (the same ones fia-east fixed for SN/NE/CS/LS):
#   * `IF(T .GT. 35000.) T=35000.` (nc/morts.f:322, ut/morts.f:337) was never applied — it follows the DIA0<0.3 reset;
#   * the IPASS QMD-convergence loop broke out after a background pass; morts.f tests convergence after EVERY pass.
# Stands: one real tiered stand (NC 15303130010497, UT 2390216010690) with its trees replaced by 4 plots × 2 species ×
# {0.4", 0.9", 1.6" seedlings, 9" overstory} (test/fixtures/dense_morts, built by /workspace/.postswap/int34/dense/mk.py).
# Before the fix only the inventory row matched; the rows below match live bit-for-bit. Later cycles still differ
# (growth on this out-of-distribution stand) — tracked as a dig item, not asserted here.
using Test
using FVSjl

const DM_DIR = joinpath(@__DIR__, "..", "fixtures", "dense_morts")
_dm_rows(text) = [split(strip(ln)) for ln in split(text, '\n')
                  if (f = split(strip(ln)); length(f) >= 20 && (y = tryparse(Int, f[1])) !== nothing && y >= 1000)]

@testset "NC/UT dense-stand morts vs live" begin
    for (v, cn, nexact) in (("nc", "15303130010497", 3), ("ut", "2390216010690", 2))
        dir = mktempdir()
        key = joinpath(dir, "s.key")
        write(key, """
STDIDENT
$cn
DATABASE
DSNin
$(abspath(joinpath(DM_DIR, v * "_stands.db")))
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         6.0
ECHOSUM
PROCESS
STOP
""")
        live = _dm_rows(read(joinpath(DM_DIR, v * ".live.sum"), String))
        jl = _dm_rows(FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code(uppercase(v))))
        @testset "$v" begin
            @test length(jl) == length(live)
            @test parse(Int, live[1][3]) > 35000          # the stand really is above the T clamp
            for k in 1:nexact
                @test live[k][1:21] == jl[k][1:21]
            end
        end
    end
end
