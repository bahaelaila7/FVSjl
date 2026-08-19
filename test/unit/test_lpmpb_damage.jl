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

# FORTYP/SIZ/STK classification (western): the national FIA stocking CSVs were missing from the
# western variants, so their .sum classification columns defaulted to 999/55. With the CSVs added,
# the IE lp_popdy stand classifies to the oracle's real codes bit-exact across all cycles.
@testset "FORTYP/SIZ/STK classification (IE lp_popdy vs FVSie_lpmpb)" begin
    # trailing two whitespace fields of each data row = FORTYP and SIZ|STK; parse per year.
    fscls(txt) = Dict(parse(Int, f[1]) => (f[end-1], f[end])
        for f in (split(strip(l)) for l in split(txt, '\n'))
        if length(f) >= 25 && occursin(r"^(19|20)\d\d$", f[1]))
    jl = fscls(FVSjl.run_keyfile(joinpath(_DFX, "lp_ctl.key"); variant = FVSjl.InlandEmpire(), output = :sum))
    orc = fscls(read(joinpath(_DFX, "oracle_ctl.sum"), String))
    @test jl[1990] == ("281", "14")     # lodgepole type (was 999/55 before the national CSVs were added)
    for yr in keys(orc)                 # every cycle bit-exact incl the post-outbreak 999/55 nonstocked ones
        @test jl[yr] == orc[yr]
    end
end
