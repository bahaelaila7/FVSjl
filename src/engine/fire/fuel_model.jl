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
