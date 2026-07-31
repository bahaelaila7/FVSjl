# =============================================================================
# crown.jl (inlandempire) — IE per-tree CCF (ie/ccfcal.f MODE=1) + crown width (MODE=2).
# stand CCF = Σ CCFT is the RELDEN the DG (dgf!) and height (htgf) chunks read (western polynomial,
# NOT the eastern crown-width→area path). Per-species dispatch (3 classes) + a large-tree/small-tree split.
#   CCFT poly (D large):   RD1(sp) + D*RD2(sp) + D^2*RD3(sp)
#   CCFT power (D small):  RDA(sp) * D^RDB(sp)
# =============================================================================

const IE_RD1 = Float32[.03,.02,.11,.04,.03,.03,.01925,.03,.03,.03,.03, .02,.01925,.03, .01925,.01925,.01925, .03,.03,.03,.03,.03,.03]
const IE_RD2 = Float32[.0167,.0148,.0333,.0270,.0215,.0238,.01676,.0173,.0216,.0180,.0215, .0148,.01676,.0216, .01676,.01676,.01676, .0238,.0215,.0238,.0238,.0215,.0215]
const IE_RD3 = Float32[.00230,.00338,.00259,.00405,.00363,.00490,.00365,.00259,.00405,.00281,.00363, .00338,.00365,.00405, .00365,.00365,.00365, .00490,.00363,.00490,.00490,.00363,.00363]
const IE_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.011109, .007244,.009187,.011402, .009187,.009187,.009187, .008915,.011109,.008915,.008915,.011109,.011109]
const IE_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.7250,1.7800,1.7600,1.7360,1.7560,1.7680,1.7250, 1.8182,1.7600,1.7560, 1.7600,1.7600,1.7600, 1.7800,1.7250,1.7800,1.7800,1.7250,1.7250]

# crown-width B1..B6 (ie/ccfcal.f MODE=2, exp form for sp 1:10,12:22; sp 11/23 special).
const IE_CWB1 = Float32[1.04050,1.02478,1.01685,1.03030,1.02460,1.03597,1.03992,1.02687,1.02886,1.02687,0.0, 1.02478,1.03992,1.02886, 1.03992,1.03992,1.03992, 1.03597,1.02460, 1.03597,1.03597, 1.02460,0.0]
const IE_CWB2 = Float32[1.27990,0.99889,1.48372,1.14079,1.35223,1.46111,1.58777,1.28027,1.01255,1.49085,0.0, 0.99889,1.58777,1.01255, 1.58777,1.58777,1.58777, 1.46111,1.35223, 1.46111,1.46111, 1.35223,0.0]
const IE_CWB3 = Float32[0.11941,0.19422,0.27378,0.20904,0.24844,0.26289,0.30812,0.22490,0.30374,0.18620,0.0, 0.19422,0.30812,0.30374, 0.30812,0.30812,0.30812, 0.26289,0.24844, 0.26289,0.26289, 0.24844,0.0]
const IE_CWB4 = Float32[0.42745,0.59423,0.49646,0.38787,0.41212,0.18779,0.64934,0.47075,0.37093,0.68272,0.0, 0.59423,0.64934,0.37093, 0.64934,0.64934,0.64934, 0.18779,0.41212, 0.18779,0.18779, 0.41212,0.0]
const IE_CWB5 = Float32[0.0,-0.09078,-0.18669,0.0,-0.10436,0.0,-0.38964,-0.15911,-0.13731,-0.28242,0.0, -0.09078,-0.38964,-0.13731, -0.38964,-0.38964,-0.38964, 0.0,-0.10436, 0.0,0.0, -0.10436,0.0]
const IE_CWB6 = Float32[-0.07182,-0.02341,-0.01509,0.0,0.03539,0.0,0.0,0.0,0.0,0.0,0.0, -0.02341,0.0,0.0, 0.0,0.0,0.0, 0.03539,0.0, 0.0,0.0, 0.03539,0.0]

"""Per-tree CCF contribution (ie/ccfcal.f MODE=1), excluding the ×TPA factor."""
@inline function ie_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    D = Float32(d)
    if sp <= 12 || sp == 14 || sp == 23
        return D >= 10f0 ? IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp] : IE_RDA[sp] * D^IE_RDB[sp]
    elseif sp == 19 || sp == 22
        D >= 10f0 && return IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp]
        return D > 0.1f0 ? IE_RDA[sp] * D^IE_RDB[sp] : 0.001f0
    else  # sp 13,15,16,17,18,20,21
        D >= 1f0 && return IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp]
        return D > 0.1f0 ? IE_RDA[sp] * D^IE_RDB[sp] : 0.001f0
    end
end

"""IE crown width (ie/ccfcal.f MODE=2). JCR = crown ratio (%); barea = BA (or OLDBA under thin/fire).
IFOR==5 (Colville) uses R6CRWD (deferred). Returns width clamped [0.1, 99.9]."""
@inline function ie_crown_width(sp::Integer, d::Real, h::Real, jcr::Integer, barea::Real)::Float32
    (jcr <= 0 || h <= 0 || d <= 0 || barea <= 0) && return 0.1f0
    D = Float32(d); H = Float32(h)
    cl = Float32(jcr) * H * 0.01f0
    cl <= 0f0 && return 0.1f0
    cw = if sp == 11 || sp == 23
        if H <= 5f0
            0.8f0 * H * max(0.5f0, jcr * 0.01f0)
        elseif H >= 15f0
            6.90396f0 * D^0.55645f0 * H^(-0.28509f0) * cl^0.20430f0
        else
            c1 = 0.8f0 * H * max(0.5f0, jcr * 0.01f0)
            c2 = 6.90396f0 * D^0.55645f0 * H^(-0.28509f0) * cl^0.20430f0
            w = (H - 5f0) * 0.1f0
            c1 * (1f0 - w) + c2 * w
        end
    else
        ba = barea < 1f0 ? 1f0 : Float32(barea)
        IE_CWB1[sp] * exp(IE_CWB2[sp] + IE_CWB3[sp]*log(cl) + IE_CWB4[sp]*log(D) +
                          IE_CWB5[sp]*log(H) + IE_CWB6[sp]*log(ba))
    end
    cw > 99.9f0 && (cw = 99.9f0)
    cw < 0.1f0 && (cw = 0.1f0)
    return cw
end
