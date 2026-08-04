# BC V2 (LV2ATV) CONSPP resolution validation — the NI branch (dgf.f:2061-2068, LSPPOK=false, e.g. sp14).
# DGCON = DGHAB(MAPHAB(ITYPE,I),I) + DGFOR(MAPLOC(IFOR,I),I) + DGEL·ELEV + DGEL2·ELEV² +
#         (DGSASP·sin(ASP)+DGCASP·cos(ASP)+DGSLOP)·SLOPE + DGSLSQ·SLOPE².  Coeffs = data/britishcolumbia/v2dg_*.csv.
# Reference: instrumented FVSbc DGCONS dump on all_BC_essf (ESSFdk/01), sp14: ITYPE=4 IFOR=4 → DGCON=0.7631423.
const DIR = "/workspace/FVSjl/data/britishcolumbia"
_rd(n) = [split(l, ",") for l in split(strip(read(joinpath(DIR, "v2dg_$n.csv"), String)), "\n")][2:end]
_g(n, sp, c) = parse(Float32, _rd(n)[sp][c + 1])   # c is 1-based column

"""V2 NI-branch DGCON (dgf.f:2061-2068). `itype`,`ifor` from becset (habitat/forest); angles in rad."""
function bc_v2_dgcon_ni(sp::Int, itype::Int, ifor::Int, elev::Float32, asp::Float32, slope::Float32)
    isphab = parse(Int, _rd("MAPHAB")[sp][itype + 1])
    ispfor = parse(Int, _rd("MAPLOC")[sp][ifor + 1])
    _g("DGHAB", sp, isphab) + _g("DGFOR", sp, ispfor) + _g("DGEL", sp, 1)*elev + _g("DGEL2", sp, 1)*elev^2 +
        (_g("DGSASP", sp, 1)*sin(asp) + _g("DGCASP", sp, 1)*cos(asp) + _g("DGSLOP", sp, 1))*slope +
        _g("DGSLSQ", sp, 1)*slope^2
end

function main()
    # oracle reference (BCV2CON dump, sp14, ESSFdk/01)
    dgcon = bc_v2_dgcon_ni(14, 4, 4, 0.2297f0, 5.4978f0, 0.30f0)
    ref = 0.7631423f0
    println("BC V2 CONSPP NI-branch (sp14, ESSF): jl DGCON=", round(dgcon, digits=7),
            " oracle=$ref  Δ=", round(abs(dgcon - ref), digits=7),
            abs(dgcon - ref) < 1f-5 ? "  ✓ MATCH" : "  ✗")
end
main()
