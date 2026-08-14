# =============================================================================
# dgcon_validate.jl — CA chunk-3 DGCON (ca/dgf.f ENTRY DGCONS), bit-exact vs FVSca_dump.
#
# Reference: ref_dgcon_cat01.txt — per-species {ISPC, JSPC, ISPFOR, DGCON, DGDSQ, SITEAR, ELEV, SLOPE,
# ASPECT} from an instrumented ca/dgf.f (unconditional ENTRY-DGCONS WRITE), cat01 (forest 610→IFOR 6,
# elev 35 (hundred-ft), slope 30%, aspect 315°). CA's dgf is GROUP-COMPRESSED (MAPSPC → 13 groups).
# Reproduces the DGCON formula standalone from the shipped CA_* consts + the live SITEAR/ELEV/SLOPE/ASPECT
# (decoupled from the site-index chunk). Sweeps IFOR to pin forest 610 → IFOR (only group 4/DF is
# forest-dependent). GS/RW (ISPC 23/50) use the ln(SITEAR) constant form.
# =============================================================================
using FVSjl
const M = FVSjl

function ca_dgcon_ref(isp::Int, si::Float32, elev::Float32, slope::Float32, asp::Float32, ifor::Int)
    jspc = M.CA_MAPSPC[isp]
    (isp == 23 || isp == 50) && return -3.502444f0 + 0.415435f0 * log(si)
    isfor = M.CA_MAPLOC[jspc, ifor]
    sina = sin(asp); cosa = cos(asp)
    sasp = (M.CA_DGSASP[jspc]*sina + M.CA_DGCASP[jspc]*cosa + M.CA_DGSLOP[jspc])*slope +
           M.CA_DGSLSQ[jspc]*slope*slope
    return M.CA_DGFOR[jspc, isfor] + M.CA_DGEL[jspc]*elev + M.CA_DGELSQ[jspc]*elev*elev +
           M.CA_DGSITE[jspc]*log(si) + sasp
end

function run_validation()
    ref = joinpath(@__DIR__, "ref_dgcon_cat01.txt")
    rows = NTuple{6,Float32}[]   # isp, dgcon_live, sitear, elev, slope, aspect
    for l in eachline(ref)
        f = split(strip(l)); isempty(f) && continue
        push!(rows, (parse(Float32, f[2]), parse(Float32, f[5]), parse(Float32, f[7]),
                     parse(Float32, f[8]), parse(Float32, f[9]), parse(Float32, f[10])))
    end
    # sweep IFOR 1..10 to find forest 610's mapping (min worst-|Δ|)
    best_ifor = 0; best_worst = Inf32
    for ifor in 1:10
        w = maximum(abs(ca_dgcon_ref(Int(isp), si, el, sl, as, ifor) - dc)
                    for (isp, dc, si, el, sl, as) in rows)
        w < best_worst && (best_worst = w; best_ifor = ifor)
    end
    println("CA DGCON: forest 610 → IFOR=$best_ifor (worst |Δ|=$best_worst)")
    nfail = 0; worst = 0f0
    for (isp, dc, si, el, sl, as) in rows
        jl = ca_dgcon_ref(Int(isp), si, el, sl, as, best_ifor); d = abs(jl - dc)
        worst = max(worst, d)
        d > 1f-3 && (nfail += 1; println("MISMATCH isp=$(Int(isp)): jl=$jl live=$dc (Δ$d)"))
    end
    println("CA DGCON: $(length(rows)) species, $(length(rows)-nfail) bit-exact, worst |Δ|=$worst")
    nfail == 0 || error("CA DGCON validation FAILED ($nfail)")
    return best_ifor
end
run_validation()
