# =============================================================================
# dgcon_validate.jl — PN chunk-3 DGCON (pn/dgf.f ENTRY DGCONS), bit-exact vs FVSpn_g16.
#
# Reference: ref_dgcon_pnt01.txt — per-species {ISPC, JSPC, DGCON, SITEAR} from an instrumented pn/dgf.f
# (unconditional ENTRY-DGCONS WRITE, no DEBUG ⇒ no NATCRS crash), pnt01 (forest 612 → IFOR 2, elev 7.0,
# slope 30%, aspect 315°). Reproduces pn_dgcons!'s DGCON standalone from the shipped PN_* consts + the
# live SITEAR (decoupled from the site-index chunk). This guards the 20-group coefficient transcription +
# the DGCON formula (WO King's-SI at jspc 19, no JFOR remap, DGFOR 20×3, ES ×3.281, group-14 elev cap, RW).
# =============================================================================
using FVSjl
const M = FVSjl

const _IFOR = 2                                # pnt01 forest 612 → forkod IFOR 2 (MEASURED: only ifor giving 0 Δ)
const _ELEV = 7.0f0; const _SLOPE = 0.30f0
const _SINA = sin(deg2rad(315f0)); const _COSA = cos(deg2rad(315f0))

function pn_dgcon_ref(isp::Int, si::Float32)
    jspc = M.PN_MAPSPC[isp]
    isp == 17 && return -3.502444f0 + 0.415435f0 * log(max(si, 1f0))
    isfor = M.PN_MAPLOC[jspc, _IFOR]
    sasp = (M.PN_DGSASP[jspc]*_SINA + M.PN_DGCASP[jspc]*_COSA + M.PN_DGSLOP[jspc])*_SLOPE + M.PN_DGSLSQ[jspc]*_SLOPE*_SLOPE
    xsite = si
    jspc == 10 && (xsite *= 3.281f0)
    jspc == 19 && (xsite = -37.60812f0 * log(1f0 - (xsite/114.24569f0)^0.4444f0))
    temel = _ELEV; (jspc == 14 && temel > 30f0) && (temel = 30f0)
    return M.PN_DGFOR[jspc, isfor] + M.PN_DGEL[jspc]*temel + M.PN_DGEL2[jspc]*temel*temel +
           M.PN_DGSITE[jspc]*log(max(xsite, 1f0)) + sasp
end

function run_validation()
    ref = joinpath(@__DIR__, "ref_dgcon_pnt01.txt")
    worst = 0.0f0; n = 0; nfail = 0
    for line in eachline(ref)
        isempty(strip(line)) && continue
        f = split(strip(line))
        isp = parse(Int, f[1]); dcon = parse(Float32, f[3]); si = parse(Float32, f[4])
        jl = pn_dgcon_ref(isp, si); d = abs(jl - dcon); worst = max(worst, d); n += 1
        d > 1f-3 && (nfail += 1; println("MISMATCH isp=$isp: jl=$jl live=$dcon (Δ$d)"))
    end
    println("PN DGCON: $n species, $(n-nfail) bit-exact, worst |Δ|=$worst")
    nfail == 0 || error("PN DGCON validation FAILED ($nfail)")
end
run_validation()
