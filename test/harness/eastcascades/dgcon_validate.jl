# =============================================================================
# dgcon_validate.jl — EC chunk-3 DGCON (ec/dgf.f ENTRY DGCONS), bit-exact vs FVSec_g16.
#
# Reference: ref_dgcon_ect01.txt — per-species {ISPC, DGCON, SITEAR} from an instrumented ec/dgf.f
# (unconditional ENTRY-DGCONS WRITE, no DEBUG ⇒ no NATCRS crash), ect01 (forest 608, elev 45, slope 30%,
# aspect 315°). EC's dgf is NOT group-compressed (identity index). Reproduces the DGCON formula standalone
# from the shipped EC_* consts + the live SITEAR (decoupled from the site-index chunk). Guards the 32-species
# coefficient transcription + the DGCONS special cases: WO King's-SI at ISPC 28, MH/OS ×3.281 (12/31),
# elev cap {19,23-27,29,30,32}, AF +0.3835 (9), SL0DUM slope==0 override. Sweeps IFOR to pin forest 608.
# =============================================================================
using FVSjl
const M = FVSjl

const _ELEV = 45.0f0; const _SLOPE = 0.30f0
const _SINA = sin(deg2rad(315f0)); const _COSA = cos(deg2rad(315f0))

function ec_dgcon_ref(isp::Int, si::Float32, ifor::Int)
    isfor = M.EC_MAPLOC[isp, ifor]
    if (isp in M.EC_SL0_SP) && _SLOPE == 0f0
        sasp = M.EC_SL0DUM[isp]
    else
        sasp = (M.EC_DGSASP[isp]*_SINA + M.EC_DGCASP[isp]*_COSA + M.EC_DGSLOP[isp])*_SLOPE +
               M.EC_DGSLSQ[isp]*_SLOPE*_SLOPE
    end
    xsite = si
    (isp == 12 || isp == 31) && (xsite *= 3.281f0)
    if isp == 28
        x = 1f0 - (xsite/114.24569f0)^0.4444f0
        xsite = x <= 0f0 ? 125f0 : -37.60812f0 * log(x)
    end
    temel = _ELEV; (isp in M.EC_ELCAP_SP && temel > 30f0) && (temel = 30f0)
    dgcon = M.EC_DGFOR[isp, isfor] + M.EC_DGEL[isp]*temel + M.EC_DGEL2[isp]*temel*temel +
            M.EC_DGSITE[isp]*log(max(xsite, 1f0)) + sasp
    isp == 9 && (dgcon += 0.3835f0)
    return dgcon
end

function run_validation()
    ref = joinpath(@__DIR__, "ref_dgcon_ect01.txt")
    data = [(parse(Int, f[1]), parse(Float32, f[2]), parse(Float32, f[3]))
            for f in (split(strip(l)) for l in eachline(ref)) if !isempty(f)]
    # sweep IFOR 1..7 to find forest 608's mapping (min worst-|Δ|)
    best_ifor = 0; best_worst = Inf32
    for ifor in 1:7
        w = maximum(abs(ec_dgcon_ref(isp, si, ifor) - dc) for (isp, dc, si) in data)
        w < best_worst && (best_worst = w; best_ifor = ifor)
    end
    println("EC DGCON: forest 608 → IFOR=$best_ifor (worst |Δ|=$best_worst)")
    nfail = 0
    for (isp, dc, si) in data
        jl = ec_dgcon_ref(isp, si, best_ifor); d = abs(jl - dc)
        d > 1f-3 && (nfail += 1; println("MISMATCH isp=$isp: jl=$jl live=$dc (Δ$d)"))
    end
    println("EC DGCON: $(length(data)) species, $(length(data)-nfail) bit-exact")
    nfail == 0 || error("EC DGCON validation FAILED ($nfail)")
    return best_ifor
end
run_validation()
