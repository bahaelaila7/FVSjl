# LPMPB INVMORT live-inventory damage-code path (dampro.f→MPBDAM→mpsdlp GREINF→colmod).
# A tree flagged with an MPB (agent 2, severity 3) damage pair seeds GREINF (infested TPA per COLIND
# size class), which raises the cycle-1 outbreak mortality. Validated vs FVSie_lpmpb: the with-damage
# minus without-damage MORT DELTA at the outbreak cycle (1990) matches the oracle (+2 TPA). The absolute
# MORT carries the pre-existing IE COL-path multi-cycle straddle (present in the control too), so the
# test asserts the DAMAGE-CODE DELTA, which is what this port adds.

using Test
using FVSjl

const _DFX = joinpath(@__DIR__, "..", "fixtures", "lpmpb_dam")

# field 25 of a .sum data row is MORT; field 1 = year.
function _mort_by_year(sumtext)
    d = Dict{Int,Int}()
    for ln in split(sumtext, '\n')
        f = split(strip(ln))
        (length(f) >= 25 && occursin(r"^(19|20)\d\d$", f[1])) || continue
        d[parse(Int, f[1])] = parse(Int, f[25])
    end
    return d
end
_mort_file(p) = _mort_by_year(read(p, String))

@testset "LPMPB INVMORT damage-code GREINF (delta vs FVSie_lpmpb)" begin
    dam = _mort_by_year(FVSjl.run_keyfile(joinpath(_DFX, "lp_dam.key"); variant = FVSjl.InlandEmpire(), output = :sum))
    ctl = _mort_by_year(FVSjl.run_keyfile(joinpath(_DFX, "lp_ctl.key"); variant = FVSjl.InlandEmpire(), output = :sum))
    ordam = _mort_file(joinpath(_DFX, "oracle_dam.sum"))
    orctl = _mort_file(joinpath(_DFX, "oracle_ctl.sum"))

    # the outbreak fires at 1990 (MPBSTART cycle 1); the damage codes add GREINF mortality there.
    @test haskey(dam, 1990) && haskey(ctl, 1990)
    jl_delta  = dam[1990] - ctl[1990]
    or_delta  = ordam[1990] - orctl[1990]
    @test or_delta == 2                 # oracle: damage codes add +2 TPA at the outbreak cycle
    @test jl_delta == or_delta          # jl reproduces the exact damage-code mortality delta
end
