# =============================================================================
# fire/fuel_model.jl — dynamic surface fuel model construction (FFE chunk F5b)
#
# Ported from: bin/FVSsn_buildDir/fmcfmd2.f (FMCFMD3, the SN custom model) +
# fmgfmv.f (the dead-herb moisture split) + the FMINIT USAV/UBD/CANMHT defaults.
#
# The SN FFE does not use a static fuel model — it builds a *custom* one each fire from
# the stand's own fuels: the down-wood pools (FireState.cwd, F3), the live herb/shrub
# load (FireState.flive, F3), and the understory crown biomass (crown_biomass, F2). The
# result (loads, SAV, depth, moisture of extinction) is exactly what the Rothermel model
# (F5) consumes — so this is the keystone that ties F2+F3 into the fire behavior and
# makes the crown-biomass chunk a live input rather than an inert one.
# =============================================================================

const _FM_USAV   = (2000f0, 1800f0, 1500f0)   # USAV: dead-1hr / live-herb / live-woody SAV (fminit.f:826)
const _FM_UBD    = (0.10f0, 0.75f0)           # UBD: fuelbed bulk-density bounds (fminit.f:829)
const _FM_CANMHT = 6.0f0                       # CANMHT: understory height threshold, ft (fminit.f:147)
const _TONS_TO_LBFT2 = 0.04591f0              # tons/acre → lb/ft²

@inline _fm_algslp2(x, x1, x2, y1, y2) =       # 2-point clamped linear interpolation (ALGSLP)
    x < x1 ? y1 : x >= x2 ? y2 : y1 + (y2 - y1) / (x2 - x1) * (x - x1)

"""
    standard_fuel_model(coef, model) -> (load, sav, depth, mext)

Rothermel inputs for one standard fire-behavior fuel model (1–13, the Anderson models,
fminit.f / `fire_fuel_models.csv`). `load[2,4]`/`sav[2,4]` are loads (lb/ft²) and surface-
area-to-volume by [1=dead/2=live, class]; the 10-hr / 100-hr / dead-herb / live-herb SAVs
take the FFE constant defaults (109 / 30 / 1500 / 1500). `depth` is bed depth (ft), `mext`
the dead moisture of extinction.
"""
function standard_fuel_model(coef::SpeciesCoefficients, model::Integer)
    m = @view coef.ffe_fuel_models[model, :]   # [sav_1hr, sav_lwoody, l_1hr,l_10,l_100,l_lwoody,l_lherb, depth, mext]
    load = zeros(Float32, 2, 4); sav = zeros(Float32, 2, 4)
    load[1, 1] = m[3]; load[1, 2] = m[4]; load[1, 3] = m[5]      # dead 1/10/100-hr
    load[2, 1] = m[6]; load[2, 2] = m[7]                          # live woody / herb
    sav[1, 1] = m[1]; sav[1, 2] = 109f0; sav[1, 3] = 30f0; sav[1, 4] = 1500f0
    sav[2, 1] = m[2]; sav[2, 2] = 1500f0
    return (load, sav, m[8], m[9])
end

# ----------------------------------------------------------------------------
# FMCFMD/FMDYN — weighted *standard* fuel-model selection (FFE chunk F4-select).
#
# The SN FFE's real surface-fire path is NOT the custom dynamic model above: FMBURN
# calls FMCFMD, which (a) picks a set of candidate standard fuel models from the FFE
# forest type + fuel moisture, then (b) hands them to FMDYN, which places the stand's
# (SMALL, LARGE) down-wood point in the 2-D fuel space and weights the candidate models
# by inverse perpendicular distance to each model's iso-line. FMFINT then runs Rothermel
# on each weighted model and sums the weighted flame (F4-weight). build_dynamic_fuel_model
# is only used when FUELMODL/IFLOGIC forces a static/custom model — not for snt01 stand 4.
# ----------------------------------------------------------------------------

# XPTS (fmcfmd.f:79): per-model iso-line as (SMALL-intercept, LARGE-intercept), tons/ac.
const _FMD_XPTS = Float32[
    5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  # models 1–7
    5  15;  5  15;                                          # models 8–9
   10  30; 15  30; 30  60; 45 100; 30  60]                  # models 10–14
const _FMD_ICLSS = 14
const _FMD_MXFMOD = 5     # MXFMOD (FMPARM.F77)

"SMALL/LARGE down-wood loads (tons/ac): classes 1–3 + litter(10) are SMALL, 4–9 are LARGE (fmtret.f:382)."
function _small_large_fuel(fs)
    small = 0f0; large = 0f0
    @inbounds for k in 1:2, l in 1:4
        small += fs.cwd[1, k, l] + fs.cwd[2, k, l] + fs.cwd[3, k, l] + fs.cwd[10, k, l]
        for j in 4:9
            large += fs.cwd[j, k, l]
        end
    end
    return small, large
end

"""
    select_fuel_models(s, mois) -> Vector{Tuple{Int,Float32}}

The SN weighted standard fuel models (FMCFMD + FMDYN) for the current stand: candidate
models chosen from the FFE forest type (`ffe_forest_type`/FMSNFT) and dead 100-hr fuel
moisture `mois[1,4]`, then weighted by the (SMALL, LARGE) down-wood point's inverse
distance to each model's iso-line. Returns up to `MXFMOD` (model, weight) pairs whose
weights sum to 1. This is the input FMFINT integrates over for the surface fire.
"""
function select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}; fire_basis::Bool = false)
    eqwt = zeros(Float32, _FMD_ICLSS)
    iffeft = ffe_forest_type(s)
    # An actual SIMFIRE burns on the start-of-cycle + 1-annual-step down wood (FVS interleaves the annual
    # fuel loop with FMBURN); the PotFire report and the no-stash case use the live (period-end) cwd.
    sm, lg = (fire_basis && s.fire.fire_smlg[1] >= 0f0) ? s.fire.fire_smlg : _small_large_fuel(s.fire)
    m14 = mois[1, 4]                                   # dead 100-hr (3+") moisture

    # --- candidate-model selection (fmcfmd.f:131) ---
    if iffeft in (1, 2, 3)                             # hardwood / hwd-pine / pine-hwd
        if sm > 6f0
            eqwt[5] = 1f0
        else
            moiswt8 = 0f0; moiswt9 = 0f0
            if m14 <= 0.15f0
                moiswt9 = 1f0
            elseif m14 > 0.25f0
                moiswt8 = 1f0
            else
                moiswt9 = 1f0 - (m14 - 0.15f0) / 0.1f0
                moiswt8 = 1f0 - (0.25f0 - m14) / 0.1f0
            end
            if sm <= 4f0
                eqwt[8] = moiswt8; eqwt[9] = moiswt9
            else
                eqwt[8] = (1f0 - (sm - 4f0) / 2f0) * moiswt8
                eqwt[9] = (1f0 - (sm - 4f0) / 2f0) * moiswt9
                eqwt[5] = 1f0 - (6f0 - sm) / 2f0
            end
        end
    elseif iffeft in (4, 8)                            # pine / saint francis
        if m14 <= 0.15f0
            eqwt[9] = 1f0
        elseif m14 > 0.25f0
            eqwt[8] = 1f0
        else
            eqwt[9] = 1f0 - (m14 - 0.15f0) / 0.1f0
            eqwt[8] = 1f0 - (0.25f0 - m14) / 0.1f0
        end
    elseif iffeft in (5, 6)                            # pine bluestem / oak savannah
        eqwt[2] = 1f0
    elseif iffeft == 7                                 # eastern redcedar
        rcht = 0f0; rctpa = 0f0
        t = s.trees
        @inbounds for i in 1:t.n
            if Int(t.species[i]) == 2                  # redcedar
                rcht += t.height[i]; rctpa += t.tpa[i]
            end
        end
        rcht = rctpa > 0f0 ? rcht / rctpa : 0f0
        if rcht > 7.5f0
            eqwt[4] = 1f0
        elseif rcht <= 4.5f0
            eqwt[6] = 1f0
        else
            eqwt[6] = 1f0 - (rcht - 4.5f0) / 3f0
            eqwt[4] = 1f0 - (7.5f0 - rcht) / 3f0
        end
    elseif iffeft == 9
        eqwt[6] = 1f0
    end
    # models 10 & 12 are always candidates for natural fuels (fmcfmd.f:202)
    eqwt[10] = 1f0; eqwt[12] = 1f0

    return _fmdyn(sm, lg, eqwt)
end

"""
    _fmdyn(sm, lg, eqwt) -> Vector{Tuple{Int,Float32}}

FMDYN (fmdyn.f): resolve candidate fuel models (`eqwt[i] > 0`) into weighted models by
the inverse perpendicular distance from the point `(sm, lg)` to each model's iso-line.
All SN iso-lines are sloped (ITYP≡0), so only the sloped-line geometry is ported.
Collinear candidates (identical XPTS — e.g. the litter models 1–9) share their bracket's
weight in proportion to `eqwt`. Returns up to `MXFMOD` (model, weight) pairs summing to 1.
"""
function _fmdyn(sm::Float32, lg::Float32, eqwt::Vector{Float32})
    ic = _FMD_ICLSS; mx = _FMD_MXFMOD; xpts = _FMD_XPTS
    out = Tuple{Int,Float32}[]
    (sm < 0f0 || lg < 0f0) && return out

    lok = falses(ic)
    for i in 1:ic
        eqwt[i] > 0f0 && (lok[i] = true)
    end
    # unset candidates with a zero intercept (degenerate line) — none in the SN table
    for i in 1:ic
        if lok[i] && (xpts[i, 1] == 0f0 || xpts[i, 2] == 0f0)
            lok[i] = false; eqwt[i] = 0f0
        end
    end
    # EQMOD: tag each candidate with the first candidate sharing its iso-line (collinear)
    eqmod = zeros(Int, ic)
    for i in 1:ic
        lok[i] || continue
        for j in i:ic
            if eqmod[j] == 0 && lok[j] && xpts[i, 1] == xpts[j, 1] && xpts[i, 2] == xpts[j, 2]
                eqmod[j] = i
            end
        end
    end
    # rescale each collinear group's eqwt to sum to 1
    for i in 1:ic
        lok[i] || continue
        xwt = 0f0
        for j in i:ic
            (lok[j] && eqmod[j] == i) && (xwt += eqwt[j])
        end
        if xwt > 1f-6
            for j in i:ic
                (lok[j] && eqmod[j] == i) && (eqwt[j] /= xwt)
            end
        end
    end
    # XD/YD: signed distance from the point to where each line crosses its LARGE / SMALL
    xd = zeros(Float32, ic); yd = zeros(Float32, ic)
    for i in 1:ic
        lok[i] || continue
        m1 = xpts[i, 2] / (-xpts[i, 1]); b1 = xpts[i, 2]
        xd[i] = (lg - b1) / m1 - sm
        yd[i] = (m1 * sm + b1) - lg
    end
    # nearest left/right (xd) and below/above (yd) candidate lines
    nbr = zeros(Int, 4)
    prv = Float32[-9.99f30, 9.99f30, -9.99f30, 9.99f30]
    for i in 1:ic
        lok[i] || continue
        if xd[i] < 0f0 && xd[i] > prv[1]
            prv[1] = xd[i]; nbr[1] = i
        elseif xd[i] >= 0f0 && xd[i] < prv[2]
            prv[2] = xd[i]; nbr[2] = i
        end
        if yd[i] < 0f0 && yd[i] > prv[3]
            prv[3] = yd[i]; nbr[3] = i
        elseif yd[i] >= 0f0 && yd[i] < prv[4]
            prv[4] = yd[i]; nbr[4] = i
        end
    end
    # perpendicular distance from the point to each neighbor line
    wt = zeros(Float32, 4)
    for k in 1:4
        i = nbr[k]; (i == 0 || !lok[i]) && continue
        m1 = xpts[i, 2] / (-xpts[i, 1]); b1 = xpts[i, 2]
        m2 = -(1f0 / m1); b2 = lg - m2 * sm
        npt1 = (b2 - b1) / (m1 - m2); npt2 = m2 * npt1 + b2
        wt[k] = sqrt((lg - npt2)^2 + (sm - npt1)^2)
    end
    # merge duplicate neighbors, accumulating distance into fmod/fwt
    fmod = zeros(Int, mx); fwt = zeros(Float32, mx); k2 = 0
    for i in 1:4
        nbr[i] == 0 && continue
        k2 += 1; k = k2; found = false
        for j in 1:i
            if fmod[j] == nbr[i]
                k = j; found = true; break
            end
        end
        !found && k <= mx && (fmod[k] = nbr[i])
        k <= mx && (fwt[k] += wt[i])
    end
    # weight by inverse distance, normalize
    xwt = 0f0
    for i in 1:mx
        fmod[i] == 0 && continue
        fwt[i] = 1f0 / (fwt[i] + 1f-6); xwt += fwt[i]
    end
    for i in 1:mx
        fmod[i] != 0 && (fwt[i] /= xwt)
    end
    # compact nonzero entries to the top
    k = 0
    for i in 1:mx
        fmod[i] == 0 && continue
        k += 1
        if i != k && k <= mx
            fmod[k] = fmod[i]; fwt[k] = fwt[i]; fmod[i] = 0; fwt[i] = 0f0
        end
    end
    # split each bracket's weight among its collinear models in proportion to eqwt
    fmod2 = zeros(Int, mx); fwt2 = zeros(Float32, mx); k = 1
    for i in 1:mx
        fmod[i] == 0 && continue
        ii = fmod[i]; xw = fwt[i]
        if eqmod[ii] == 0 && k <= mx
            fwt2[k] = fwt[i]; fmod2[k] = fmod[i]; k += 1
        else
            for j in 1:ic
                if eqmod[j] == eqmod[ii] && eqwt[j] > 0f0 && k <= mx
                    fwt2[k] += eqwt[j] * xw; fmod2[k] = j; k += 1
                end
            end
        end
    end
    # merge any duplicate models that arose from collinear additions
    fmod = zeros(Int, mx); fwt = zeros(Float32, mx); k2 = 0
    for i in 1:mx
        fmod2[i] == 0 && continue
        k2 += 1; k = k2
        for j in 1:i
            if fmod[j] == fmod2[i]
                k = j; break
            end
        end
        k <= mx && fmod[k] == 0 && (fmod[k] = fmod2[i])
        k <= mx && (fwt[k] += fwt2[i])
    end
    for i in 1:mx
        fmod[i] != 0 && push!(out, (fmod[i], fwt[i]))
    end
    return out
end

"""
    build_dynamic_fuel_model(s, mois) -> (load, sav, depth, mext)

Construct the SN dynamic surface fuel model (FMCFMD3, fmcfmd2.f) from the stand's fuel
state. `load[2,4]`/`sav[2,4]` are loads (lb/ft²) and surface-area-to-volume by
[1=dead/2=live, class], `depth` the fuel-bed depth (ft), `mext` the dead moisture of
extinction. Dead loads come from the down-wood pools (`fire.cwd`, with the 1-hr class =
0–.25" + litter), the live-woody load from the understory crown biomass (foliage +
½·fine for trees ≤ CANMHT) plus the live shrub, and the live-herb load from `fire.flive`.
A moisture-dependent share of the live herb is moved to a dead-herb class (fmgfmv.f).
`mois` is the fuel-moisture matrix from `fuel_moisture`.
"""
function build_dynamic_fuel_model(s::StandState, mois::AbstractMatrix{Float32})
    fs = s.fire; t = s.trees
    # down-wood pools → load by size class (tons/ac → lb/ft²)
    currcwd = zeros(Float32, 11)
    @inbounds for j in 1:11, k in 1:2, l in 1:4
        currcwd[j] += fs.cwd[j, k, l] * _TONS_TO_LBFT2
    end
    herb = fs.flive[1] * _TONS_TO_LBFT2
    # understory live-woody load: crown foliage + ½ of the 0–.25" crown for trees ≤ CANMHT
    woody = 0f0
    @inbounds for i in 1:t.n
        (t.tpa[i] > 0f0 && t.height[i] <= _FM_CANMHT) || continue
        xv = crown_biomass(s, t.species[i], t.dbh[i], t.height[i], Int(t.crown_pct[i]))
        woody += (xv[1] + 0.5f0 * xv[2]) * t.tpa[i] * _FM_P2T    # ×P2T undoes crown_biomass's /P2T
    end
    woody = (woody + fs.flive[2]) * _TONS_TO_LBFT2

    load = zeros(Float32, 2, 4); sav = zeros(Float32, 2, 4)
    load[1, 1] = max(0f0, currcwd[1] + currcwd[10])  # 1-hr: 0–.25" + litter
    load[1, 2] = max(0f0, currcwd[2])                # 10-hr
    load[1, 3] = max(0f0, currcwd[3])                # 100-hr
    load[2, 1] = max(0f0, woody)                     # live woody
    load[2, 2] = max(0f0, herb)                      # live herb
    sav[1, 1] = _FM_USAV[1]; sav[1, 2] = 109f0; sav[1, 3] = 30f0
    sav[2, 1] = _FM_USAV[3]; sav[2, 2] = _FM_USAV[2]

    # dead-herb split: when herb is dry (moisture < 1.2) move part of the live herb into a
    # dead-herb class, which takes the live-herb SAV (fmgfmv.f:79/88-97).
    if load[2, 2] > 0f0 && mois[2, 2] < 1.2f0
        wt = _fm_algslp2(mois[2, 2], 0.30f0, 1.2f0, 0f0, 1f0)
        load[1, 4] = (1f0 - wt) * load[2, 2]; sav[1, 4] = sav[2, 2]
        load[2, 2] = wt * load[2, 2]
    end

    # fuel-bed depth from a load-weighted bulk density; moisture of extinction (fmcfmd2.f:582)
    fdfl = currcwd[1] + currcwd[10]
    ffl  = fdfl + herb + woody
    wf   = ffl > 0f0 ? fdfl / ffl : 0f0
    bdavg = _FM_UBD[1] + wf * (_FM_UBD[2] - _FM_UBD[1])
    depth = bdavg > 0f0 ? (ffl + currcwd[2] + currcwd[3]) / bdavg : 0f0
    mext  = (12f0 + 480f0 * bdavg / 32f0) / 100f0
    return (load, sav, depth, mext)
end
