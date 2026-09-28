# test_sn_resid_tiered.jl — SN per-record fixes vs the LIVE FVSsn_g16 tiered goldens (test/fixtures/tiered/sn: stands.db +
# <stand>_<regime>.key + the live .sum / DBS tables). Stand 157577477010854 (LP plantation, NONE + TREELIST, 1972-1997) —
# the SN suite's largest residual block (11,599 cells) — is walked year by year; each testset names the Fortran it follows.
module SNResidTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

const _C = Ref{Any}(nothing)
function _case157()
    _C[] === nothing || return _C[]
    d = mktempdir()
    txt, db, crashed, _ = run_case("SN", "157577477010854", "none"; dir = d)
    _C[] = (crashed = crashed, ms = compare_case("SN", "157577477010854", "none", txt, db))
end
_cells(file, cols, yrs) = count(m -> m.file == file && m.col in cols && m.year in yrs, _case157().ms)

# sn/crown.f:156-159 ISORT(IND(JJ)) = ITRN−JJ+1 over FVS's IND — at the inventory CRATET's IND1-seeded RDPSRT(.FALSE.)
# (sn/cratet.f:155-157; :261 RDPSRT(.TRUE.) with dead records) — so tied diameters rank in RDPSRT's order, and the
# Weibull X = ISORT/ITRN·SCALE dubs their crowns accordingly. jl ranked by a stable ascending sortperm (MEASURED 1972:
# tied LP pairs 7.1"/49.0', 6.4"/45.4', 5.9"/42.6' had dubbed CRs 34/35, 31/29, 26/27 swapped against live).
@testset "SN crown ISORT from FVS's RDPSRT IND (sn/crown.f:156-159) vs FVSsn_g16" begin
    @test !_case157().crashed
    @test _cells("FVS_TreeList", ("PctCr", "CrWidth"), ("1972",)) == 0
end

end # module
