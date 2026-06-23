# =============================================================================
# fire/rothermel.jl — Rothermel surface-fire behavior (FFE chunk F5 core)
#
# Ported from: bin/FVSsn_buildDir/fmfint.f (FMFINT, the single-fuel-model body).
#
# The Rothermel (1972) surface fire-spread model: given a fuel model (loads, surface-
# area-to-volume by size class, bed depth, dead moisture of extinction, heat content)
# and the environment (fuel moistures, wind, slope), it returns the rate of spread,
# reaction intensity, Byram fireline intensity, and flame length. The Byram intensity
# is exactly the input the ported fire-effects chain (F6) consumes:
#   `byram → scorch_height → crown_volume_scorched → fire_tree_mortality`.
#
# Fuel arrays are indexed [category, class] with category 1 = dead (classes 1hr, 10hr,
# 100hr, dead-herb) and category 2 = live (live herb, live woody). This is the per-
# model computation; FMFINT itself weights several models from FMCFMD (F5b).
# =============================================================================

# FMFINT physical constants (fmfint.f:91-93).
const _RG_RHOP   = 32.0f0      # oven-dry particle density (lb/ft³)
const _RG_TMIN   = 0.0555f0    # total mineral content (fraction)
const _RG_SILFRE = 0.01f0      # silica-free mineral content (fraction)

"""
    rothermel_surface_fire(load, sav, depth, mext_dead, mois;
                           lhv=8000, wind=0, slope_tan=0) -> (; byram, flame, spread, sigma, xir)

Rothermel surface fire behavior for one fuel model (FMFINT, fmfint.f). `load[2,4]` are
fuel loads (lb/ft²), `sav[2,4]` the surface-area-to-volume ratios (1/ft), `mois[2,4]`
the fuel moistures (fraction), all indexed [1=dead/2=live, class]. `depth` is bed depth
(ft), `mext_dead` the dead moisture of extinction, `lhv` heat content (BTU/lb), `wind`
midflame wind (mi/h), `slope_tan` slope tangent. Returns the Byram fireline intensity
(BTU/ft/min), flame length (ft), rate of spread (ft/min), characteristic SAV, and
reaction intensity. A fuel bed too moist to carry fire returns zeros.
"""
function rothermel_surface_fire(load::AbstractMatrix{Float32}, sav::AbstractMatrix{Float32},
                                depth::Float32, mext_dead::Float32, mois::AbstractMatrix{Float32};
                                lhv::Float32 = 8000f0, wind::Float32 = 0f0, slope_tan::Float32 = 0f0)
    zero_out = (; byram = 0f0, flame = 0f0, spread = 0f0, sigma = 0f0, xir = 0f0)
    # per-category class counts (a class is present when it carries load)
    noclas = [count(j -> load[i, j] > 0f0, 1:4) for i in 1:2]
    (noclas[1] == 0 && noclas[2] == 0) && return zero_out
    # KMIN per category; with 1-based ranges `1:noclas` is empty when a category has no
    # classes, matching the Fortran's KMAX≠0 guard (so IFINES(2)=0/NL=0 needs no special case).
    ifines = (1, 1)

    # order each category's classes finest-first (descending SAV); drop classes SAV<16
    isize = [sortperm([sav[i, j] for j in 1:4]; rev = true) for i in 1:2]
    mext = Float32[mext_dead, 0f0]
    for i in 1:2
        k = noclas[i]
        for kk in 1:noclas[i]
            if sav[i, isize[i][kk]] < 16f0; k = kk - 1; break; end
        end
        noclas[i] = k
    end

    # weighting factors A (area), F (area fraction), FX (category fraction)
    A  = zeros(Float32, 2, 4); F = zeros(Float32, 2, 4); GS = zeros(Float32, 2, 4)
    WO = zeros(Float32, 2, 4); AI = zeros(Float32, 2)
    for i in 1:2, kk in ifines[i]:noclas[i]
        j = isize[i][kk]
        A[i, j]  = load[i, j] * (sav[i, j] / _RG_RHOP)
        GS[i, j] = exp(-138f0 / (sav[i, j] + 1f-9))
        AI[i]   += A[i, j]
        WO[i, j] = load[i, j] * (1f0 - _RG_TMIN)
    end
    for i in 1:2, kk in ifines[i]:noclas[i]
        j = isize[i][kk]; AI[i] != 0f0 && (F[i, j] = A[i, j] / AI[i])
    end
    at = AI[1] + AI[2]
    fx = (AI[1] / (at + 1f-9), 0f0); fx = (fx[1], 1f0 - fx[1])

    # live moisture of extinction (dead/live fine ratio)
    fined = 0f0; finel = 0f0; wdfmn = 0f0
    for kk in ifines[1]:noclas[1]
        j = isize[1][kk]; ep = sav[1, j] != 0f0 ? exp(-138f0 / sav[1, j]) : 0f0
        wtfac = load[1, j] * ep
        m = j == 4 ? mois[1, 1] : mois[1, j]            # dead herb uses dead 1-hr moisture
        fined += wtfac; wdfmn += wtfac * m
    end
    findm = fined != 0f0 ? wdfmn / fined : 0f0
    for kk in ifines[2]:noclas[2]
        j = isize[2][kk]; ep = sav[2, j] != 0f0 ? exp(-500f0 / sav[2, j]) : 0f0
        finel += load[2, j] * ep
    end
    mext[2] = finel != 0f0 ? max(mext_dead, 2.9f0 * (fined / finel) * (1f0 - findm / mext_dead) - 0.226f0) : 100f0

    # per-category reaction-intensity intermediates
    ir = zeros(Float32, 2); sigma = 0f0
    sum1 = 0f0; sum2 = 0f0; sum3 = 0f0
    mdcsa1 = 0f0
    for i in 1:2
        (noclas[i] == 0 || ifines[i] > noclas[i]) && continue
        mcsa = 0f0; bse = 0f0; sigma1 = 0f0; lhv1 = 0f0
        aa = zeros(Float32, 5)                           # area by SAV bin (for net-load G weighting)
        for kk in ifines[i]:noclas[i]
            j = isize[i][kk]; ax = F[i, j]; s = sav[i, j]
            b = s < 48f0 ? 5 : s < 96f0 ? 4 : s < 192f0 ? 3 : s < 1200f0 ? 2 : 1
            aa[b] += A[i, j]
            m = (i == 1 && j == 4) ? mois[1, 1] : mois[i, j]
            qig = 250f0 + 1116f0 * m
            mcsa += ax * m; bse += ax * _RG_SILFRE; sigma1 += ax * s; lhv1 += ax * lhv
            sum1 += load[i, j]; sum2 += load[i, j] / _RG_RHOP; sum3 += fx[i] * F[i, j] * qig * GS[i, j]
        end
        wo1 = 0f0
        for kk in ifines[i]:noclas[i]
            j = isize[i][kk]; s = sav[i, j]
            b = s < 48f0 ? 5 : s < 96f0 ? 4 : s < 192f0 ? 3 : s < 1200f0 ? 2 : 1
            g = aa[b] / (AI[i] == 0f0 ? 1f-9 : AI[i])
            wo1 += g * WO[i, j]
        end
        beta = mcsa / (mext[i] + 1f-9)
        mdcsa = 1f0 - beta * (2.59f0 - beta * (5.11f0 - beta * 3.52f0))
        mext[i] < mcsa && (mdcsa = 0f0)
        barns = bse != 0f0 ? min(1f0, 0.174f0 / bse^0.19f0) : 0f0
        sigma += fx[i] * sigma1
        ir[i] = wo1 * lhv1 * mdcsa * barns
        i == 1 && (mdcsa1 = mdcsa)
    end
    mdcsa1 <= 0f0 && return zero_out                     # dead fuel too moist to spread

    sigma == 0f0 && (sigma = 1f-9)
    rhop1 = sum1 / depth                                 # bulk density
    beta1 = sum2 / depth                                 # packing ratio
    best  = 3.348f0 / sigma^0.8189f0                     # optimum packing ratio
    rat   = beta1 / (best + 1f-9)
    a1    = 133f0 / sigma^0.7913f0
    v     = sigma^1.5f0
    gamma = (v * rat^a1 * exp(a1 * (1f0 - rat))) / (495f0 + 0.0594f0 * v)  # opt. reaction velocity
    xir   = gamma * ir[1] + gamma * ir[2]                # reaction intensity
    rhobqig = rhop1 * sum3                               # heat sink
    b   = (0.792f0 + 0.681f0 * sqrt(sigma)) * (0.1f0 + beta1)
    xio = (xir * exp(b)) / (192f0 + 0.2595f0 * sigma)    # propagating flux
    phis = beta1 != 0f0 ? 5.275f0 * slope_tan * slope_tan / beta1^0.3f0 : 0f0
    xm1 = 0.02526f0 * sigma^0.54f0
    xn1 = 0.715f0 * exp(-0.000359f0 * sigma)
    c1  = 7.47f0 * exp(-0.133f0 * sigma^0.55f0)
    rat != 0f0 && (c1 = c1 / rat^xn1)
    w = wind * 88f0                                       # mi/h → ft/min
    phiw = c1 * w^xm1
    spread = xio * (1f0 + phis + phiw) / (rhobqig + 1f-9)
    byramt = xir * spread * 384f0 / sigma
    flame = 0.45f0 * (byramt / 60f0)^0.46f0
    return (; byram = byramt, flame = flame, spread = spread, sigma = sigma, xir = xir)
end
