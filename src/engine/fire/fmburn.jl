# =============================================================================
# fire/fmburn.jl — fire event driver (FFE chunk F5b-driver)
#
# Ported from: bin/FVSsn_buildDir/fmburn.f (the SIMFIRE → behavior → effects path) +
# fmeff.f (the per-tree kill application).
#
# `fmburn!` runs one simulated fire: it builds the stand's fuel context (FMCBA), the
# dynamic fuel model (FMCFMD3), and the fire weather (FMMOIS + wind reduction), drives
# the Rothermel surface-fire model (FMFINT) to a Byram intensity and flame length,
# converts that to a scorch height (Van Wagner), and then kills trees (FMEFF): a per-
# tree draw decides whether the record falls in the burned fraction (PSBURN), and the
# burned records lose `PMORT` of their TPA. This composes the already-ported fire
# functions into the actual `.sum`-affecting kill.
# =============================================================================

"Result of a fire: TPA killed per acre and the computed surface fire behavior."
struct FireResult
    killed::Float32      # trees/acre killed
    flame::Float32       # flame length (ft)
    byram::Float32       # Byram fireline intensity (BTU/ft/min)
    scorch::Float32      # scorch height (ft)
    carbon_released::Float32  # carbon released by fuel consumption (tons C/ac)
end

"""
    fmburn!(s; atemp, wind, fmois, psburn, mortcode, burnseas, flmult, crburn) -> FireResult

Run one simulated fire on stand `s` and apply the fire-caused mortality to tree TPA
(FMBURN/FMEFF). `atemp` air temperature (°F), `wind` 20-ft wind (mi/h), `fmois` the
dryness model (1–4), `psburn` percent of the stand burned, `mortcode` 1=FFE mortality
(0=off), `burnseas` burn season (1–4), `flmult` flame-length multiplier (FLAMEADJ),
`crburn` crown-fire fraction. A per-tree draw on the main RNG decides whether each
record is in the burned portion; burned records lose `PMORT·TPA` (plus the crown-fire
share). No-op unless FFE is active.
"""
# FVS_Mortality DBH class lower bounds (LOWDBH, fmvinit.f:53-59): 7 NON-cumulative bins
# [0,5) [5,10) [10,20) [20,30) [30,40) [40,50) [50,∞). A tree falls in the class whose lower bound it meets.
const _FM_LOWDBH = (0f0, 5f0, 10f0, 20f0, 30f0, 40f0, 50f0)
@inline function _fm_mort_class(d::Float32)::Int
    @inbounds for c in 1:7
        d < _FM_LOWDBH[c] && return c - 1
    end
    return 7
end

function fmburn!(s::StandState; atemp::Float32 = 70f0, wind::Float32 = 20f0, fmois::Integer = 1,
                 psburn::Float32 = 100f0, mortcode::Integer = 1, burnseas::Integer = 1,
                 flmult::Float32 = 1f0, crburn::Float32 = 0f0, year::Integer = 0)::FireResult
    fs = s.fire
    (fs === nothing || !fs.active) && return FireResult(0f0, 0f0, 0f0, 0f0, 0f0)
    fmcba!(s)                                            # fuel pools, cover, percent cover
    t = s.trees; coef = s.coef
    mois = fuel_moisture(fmois)
    fwind = wind * fire_wind_reduction(fs.percov)        # 20-ft wind → midflame
    # SN surface fire (FMCFMD + FMDYN + FMFINT): select the weighted standard fuel models
    # for the stand and integrate Rothermel over them, summing the weighted flame & Byram.
    models = select_fuel_models(s, mois; fire_basis = true)   # burn on start-of-cycle + 1-annual-step down wood
    flame_raw = 0f0; byram = 0f0
    for (fm, w) in models
        load, sav, depth, mext = standard_fuel_model(coef, fm)
        r = rothermel_surface_fire(load, sav, depth, mext, mois; wind = fwind)
        flame_raw += r.flame * w
        byram += r.byram * w
    end
    flame = flame_raw * flmult
    # scorch height (Van Wagner) from the (weighted) Byram intensity
    sch = byram > 0f0 ? scorch_height(byram, atemp, fwind) : 0f0

    # pre-fire total live TPA by FVS_Mortality DBH class (LOWDBH bins, 7 non-cumulative classes)
    totcls = zeros(Float32, 7); clskil = zeros(Float32, 7)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        c = _fm_mort_class(t.dbh[i]); c >= 1 && (totcls[c] += t.tpa[i])
    end
    killed = 0f0; killed_ba = 0f0; killed_vol = 0f0
    if mortcode != 0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            (rann!(s.rng) * 100f0 > psburn) && continue  # this record is in the unburned portion
            csv = crown_volume_scorched(sch, t.height[i], Int(t.crown_pct[i]))
            sp = Int(t.species[i]); d = t.dbh[i]
            pmort = fire_tree_mortality(coef, sp, d, flame, csv)
            pmort = fire_mortality_adjust(pmort, sp, d, burnseas)
            (d <= 1f0 && csv > 50f0) && (pmort = 1f0)     # fmeff.f:330
            pmort = clamp(pmort, 0f0, 1f0)
            curkil = pmort * t.tpa[i]
            crburn > 0f0 && (curkil += crburn * (t.tpa[i] - curkil))  # crown-fire share
            t.tpa[i] -= curkil
            t.tpa[i] < 0f0 && (t.tpa[i] = 0f0)
            killed += curkil
            killed_ba += curkil * 0.005454154f0 * d * d   # fire-killed basal area (ft²/ac, fmfout.f:303)
            killed_vol += curkil * t.merch_cuft_vol[i]    # SN: merch cubic volume killed (fmfout.f:306)
            c = _fm_mort_class(d); c >= 1 && (clskil[c] += curkil)
            add_snag!(fs, sp, d, curkil, year; height = t.height[i])  # fire-killed trees become standing snags
        end
    end
    # the fire consumes a share of the surface fuels — releasing carbon, leaving the rest. The CONSUMED
    # loadings (FVS_Consumption) are the before−after difference in the FFE fuel pools.
    fuel_before = ffe_fuel_loadings(s)
    carbon_released = apply_fire_consumption!(fs, mois)
    fuel_after = ffe_fuel_loadings(s)
    consumed = NamedTuple{keys(fuel_before)}(map(-, values(fuel_before), values(fuel_after)))
    # capture the burn-event record for the FVS_BurnReport / Mortality / Consumption DBS tables
    push!(fs.burn_reports, (; year = Int(year), mois = copy(mois), wind = fwind, flame = flame,
          scorch = sch, models = collect(models), killed = killed, killed_ba = killed_ba,
          killed_vol = killed_vol, released = carbon_released,
          clskil = Tuple(clskil), totcls = Tuple(totcls), consumed = consumed))
    return FireResult(killed, flame, byram, sch, carbon_released)
end

# PM2.5 smoke emission factors (fmcons.f:60-70): dead surface fuel by moisture-type (lb/ton consumed) and
# the live/crown classes. SN reports potential smoke = Σ consumed-by-class × factor.
const _FM_SMOKE_DEAD = 19.0f0          # representative dead-fuel PM2.5 factor (22.5/18.3/16.2 by moisture)
const _FM_SMOKE_LIVE = 21.3f0          # live herb/shrub PM2.5 factor (EMFACL)

"""
    potential_fire_report(s) -> NamedTuple

Bundle the FVS_PotFire report row (FMPOFL): the dual-scenario surface fire (`potential_fire`), the canopy
bulk density (`canopy_bulk_density`), and torching probabilities (`torching_probability`). In SN total
flame = surface flame and the crown-fire Torch/Crown indices are −1 (FMCFIR is skipped). Returns `nothing`
without an active fire state.
"""
function potential_fire_report(s::StandState)
    pf = potential_fire(s); pf === nothing && return nothing
    cbd = canopy_bulk_density(s)
    pt = torching_probability(s, pf.severe.flame, pf.moderate.flame)
    return (; surf_flame_sev = pf.severe.flame, surf_flame_mod = pf.moderate.flame,
            tot_flame_sev = pf.severe.flame, tot_flame_mod = pf.moderate.flame,   # = surface (no crown fire in SN)
            ptorch_sev = pt.severe, ptorch_mod = pt.moderate,
            torch_index = -1f0, crown_index = -1f0,                                # FMCFIR skipped in SN
            canopy_ht = cbd.canopy_ht, canopy_density = cbd.cbd,
            mort_ba_sev = pf.severe.ba_kill, mort_ba_mod = pf.moderate.ba_kill,
            mort_vol_sev = pf.severe.vol_kill, mort_vol_mod = pf.moderate.vol_kill,
            smoke_sev = pf.severe.smoke, smoke_mod = pf.moderate.smoke,
            models = pf.severe.models)                                            # severe-case fuel models (fmpofl.f:230)
end

"""
    canopy_bulk_density(s) -> (; cbd, actcbh, canopy_ht, tcload)

Canopy crown-fuel profile for the FVS_PotFire report (FMPOCR, fmpocr.f, SN uniform-distribution path).
Builds a 1-ft-resolution vertical crown-fuel array `CRFILL` (lbs/ac-ft): each live tree spreads its canopy
fuel `(foliage + ½ finest-woody)·TPA` uniformly over its crown (base `HT·(1−ICR/100)` to top `HT`), with
partial top/bottom layers. Returns: `cbd` = the max 13-ft running mean of CRFILL converted to kg/m³ and
capped at 0.35; `actcbh` = the actual crown base height (ft, lowest layer whose 3-ft running mean ≥ 30
lbs/ac-ft, −1 if none); `canopy_ht` = effective canopy top (ft); `tcload` = total canopy fuel (lbs/ft²).
"""
function canopy_bulk_density(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return (cbd = 0f0, actcbh = -1, canopy_ht = 0, tcload = 0f0)
    t = s.trees
    NH = 400
    crfill = zeros(Float32, NH)                         # crown fuel by 1-ft height layer (lbs/ac-ft)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        h = t.height[i]; h > 0f0 || continue
        icr = Float32(t.crown_pct[i])
        crbot = h * (1f0 - icr * 0.01f0); crbot < 0f0 && (crbot = 0f0)
        xv = crown_biomass(s, Int(t.species[i]), t.dbh[i], h, Int(round(icr)))
        crbio = (xv[1] + xv[2] * 0.5f0) * t.tpa[i]      # foliage + ½ finest woody, ×TPA (lbs/ac)
        crbio > 0f0 || continue
        len = h - crbot; len > 0f0 || continue
        adcrwn = crbio / len                            # uniform density over the crown length (lbs/ac-ft)
        i1 = Int(floor(crbot)) + 1; i2 = Int(floor(h)) + 1
        i1 > NH && (i1 = NH); i2 > NH && (i2 = NH)
        i1 <= i2 || continue
        for j in i1:i2
            adj = j == i1 ? clamp(Float32(i1) - crbot, 0f0, 1f0) :
                  j == i2 ? clamp(1f0 - (Float32(i2) - h), 0f0, 1f0) : 1f0  # partial top/bottom layers
            crfill[j] += adcrwn * adj
        end
    end
    tcload = sum(crfill) / 43560f0                       # lbs/ac → lbs/ft²
    # crown start/end = lowest/highest 1-ft layer with > 5 lbs/ac-ft
    j1 = findfirst(>(5f0), crfill); j1 === nothing && return (cbd = 0f0, actcbh = -1, canopy_ht = 0, tcload = tcload)
    j2 = findlast(>(5f0), crfill)
    cbd_lb = 0f0; actcbh = -1; abotmx = 0f0; mxj = -1
    if j1 == j2
        cbd_lb = crfill[j1]; actcbh = j1
    else
        @inbounds for j in j1:j2                         # 13-ft running mean (6 below … 6 above) → max = CBD
            a = 0f0; n = 0
            for k in max(j - 6, j1):min(j + 6, j2); a += crfill[k]; n += 1; end
            a /= n; a > cbd_lb && (cbd_lb = a)
            b = 0f0; nb = 0                              # 3-ft running mean → crown base height (≥ 30)
            for k in max(j - 1, j1):min(j + 1, j2); b += crfill[k]; nb += 1; end
            b /= nb
            b > abotmx + 0.1f0 && (abotmx = b; mxj = j)
            b >= 30f0 && actcbh == -1 && (actcbh = j)
        end
        actcbh == -1 && abotmx > 5f0 && (actcbh = mxj)
    end
    cbd = cbd_lb * 0.45359237f0 / (4046.856422f0 * 0.3048f0)   # lbs/ac-ft → kg/m³
    cbd > 0.35f0 && (cbd = 0.35f0)                              # cap (S. Rebain 2005)
    return (cbd = cbd, actcbh = actcbh, canopy_ht = Int(j2), tcload = tcload)
end

# Standard normal lower-tail CDF (FMPOFL_NPROB) — Abramowitz & Stegun 26.2.17 rational approximation.
@inline function _normal_cdf(z::Float64)::Float64
    s = z < 0 ? -1.0 : 1.0; x = abs(z) / sqrt(2.0)
    tt = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (((((1.061405429tt - 1.453152027) * tt) + 1.421413741) * tt - 0.284496736) * tt + 0.254829592) * tt * exp(-x * x)
    return 0.5 * (1.0 + s * y)
end

"""
    torching_probability(s, flame_sev, flame_mod; reps=30) -> (; severe, moderate)

Probability of crown torching under the severe / moderate flame lengths (FMPOFL_FMPTRH, fmpofl.f). A
`reps`-rep Monte Carlo: each rep draws a virtual plot (trees present with Poisson probability
`1−exp(−TPA·PSIZE)`, PSIZE=0.025), finds the lowest crown base height that must ignite for the plot to
torch given the ladder-fuel rule (a tree carries fire up if the running max height ×1.25 exceeds the next
crown base, until a tree reaches the critical height `CRIT = clamp(0.5·avg-height-of-top-40-TPA, 5, 50)`).
Torching probability for a flame length is the mean over reps of the normal CDF that the required crown
base ≤ the flame's max-needle-torch height `log((FL/0.0775)^1.45 / 30.5)` (log scale, σ=0.25). Uses the
stand RNG (`rann!`), matching FVS's RANN-based stochastic torching. Returns 0 when a flame length is ~0.
"""
function torching_probability(s::StandState, flame_sev::Float32, flame_mod::Float32; reps::Int = 30)
    t = s.trees; psize = 0.025f0
    n = 0; prb = Float32[]; cbh = Float32[]; ht = Float32[]
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 && t.height[i] > 0f0 || continue
        push!(prb, t.tpa[i]); push!(ht, t.height[i])
        push!(cbh, t.height[i] * (1f0 - Float32(t.crown_pct[i]) * 0.01f0)); n += 1
    end
    n == 0 && return (; severe = 0f0, moderate = 0f0)
    ord = sortperm(cbh)                                  # trees by crown base height ascending (RDPSRT)
    # CRIT: half the avg height of the top-40-TPA cohort (in list order), clamped [5, 50]
    avht = 0f0; ssum = 0f0
    @inbounds for i in 1:n
        p = prb[i]; ssum + p > 40f0 && (p = 40f0 - ssum)
        ssum += p; avht += ht[i] * p; ssum >= 40f0 && break
    end
    ssum > 0f0 && (avht /= ssum)
    crit = clamp(0.5f0 * avht, 5f0, 50f0)
    mincb = Float32[]                                    # required ignition crown base per torching rep
    yes = Int[]
    for _ in 1:reps
        empty!(yes); itop = false
        @inbounds for ii in 1:n
            i = ord[ii]
            (prb[i] > 1000f0 || rann!(s.rng) > exp(-prb[i] * psize)) || continue
            push!(yes, i); ht[i] >= crit && (itop = true)
        end
        itop || continue
        mc = -1f0
        @inbounds for jj in length(yes):-1:1            # scan present trees from the top down
            i = yes[jj]
            if ht[i] >= crit
                mc = cbh[i]; break
            elseif jj > 1                                # can this tree ladder fire up to CRIT?
                mxht = ht[i]
                for kk in (jj - 1):-1:1
                    j = yes[kk]
                    if mxht * 1.25f0 > cbh[j]
                        mxht < ht[j] && (mxht = ht[j])
                        if mxht >= crit; mc = cbh[i]; break; end
                    end
                end
                mc > -1f0 && break
            end
        end
        mc > 0f0 && push!(mincb, mc)
    end
    isempty(mincb) && return (; severe = 0f0, moderate = 0f0)
    p = 1.0 / reps
    function ptorch(fl::Float32)
        fl > 1f-4 || return 0f0
        mxnt = log(((Float64(fl) / 0.0775)^1.45) / 30.5)
        acc = 0.0
        # PT1 = the RIGHT tail (fmpofl.f calls NPROB(Z,Q,PT1,…): PT1 receives Q = 1−CDF), so torching is
        # likelier when the flame's reach MXNT exceeds the required crown base log(MINCB) (Z more negative).
        for mc in mincb; acc += (1.0 - _normal_cdf((log(Float64(mc)) - mxnt) / 0.25)) * p; end
        return Float32(acc)
    end
    return (; severe = ptorch(flame_sev), moderate = ptorch(flame_mod))
end

"""
    potential_fire(s) -> (; severe, moderate)

Potential SURFACE-fire behavior under the two FFE fixed weather scenarios (FMPOFL, fmpofl.f:103),
WITHOUT applying mortality — the value-grounded core of the FVS_PotFire report. SEVERE = fmois 1, 20 mph,
70°F; MODERATE = fmois 3, 8 mph, 60°F (fmvinit.f:63-66). Each scenario returns flame length, scorch
height, potential fire-killed basal area / merch volume, an estimated PM2.5 smoke (consumed fuel ×
emission factor), and the weighted standard fuel models. In SN the crown-fire spread model (FMCFIR) is
skipped, so total flame = surface flame, there are no crown indices, and CBD/Canopy are 0.
"""
function potential_fire(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return nothing
    fmcba!(s)
    coef = s.coef; t = s.trees
    function scenario(fmois::Int, wind::Float32, temp::Float32, season::Int)
        mois = fuel_moisture(fmois)
        fwind = wind * fire_wind_reduction(fs.percov)
        models = select_fuel_models(s, mois)
        flame = 0f0; byram = 0f0
        for (fm, w) in models
            load, sav, depth, mext = standard_fuel_model(coef, fm)
            r = rothermel_surface_fire(load, sav, depth, mext, mois; wind = fwind)
            flame += r.flame * w; byram += r.byram * w
        end
        sch = byram > 0f0 ? scorch_height(byram, temp, fwind) : 0f0
        ba_kill = 0f0; vol_kill = 0f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            csv = crown_volume_scorched(sch, t.height[i], Int(t.crown_pct[i]))
            sp = Int(t.species[i]); d = t.dbh[i]
            pm = fire_tree_mortality(coef, sp, d, flame, csv)
            pm = fire_mortality_adjust(pm, sp, d, season)
            (d <= 1f0 && csv > 50f0) && (pm = 1f0)
            pm = clamp(pm, 0f0, 1f0); kil = pm * t.tpa[i]
            ba_kill += kil * 0.005454154f0 * d * d
            vol_kill += kil * t.merch_cuft_vol[i]
        end
        # potential smoke (PM2.5): the surface fuel a fire would consume × emission factor (FMCONS), a
        # NON-mutating estimate — Σ cwd[size]·consumed-fraction × dead factor + live shrub/herb × live factor.
        fr = fire_consumption_fractions(mois)
        consumed = 0f0
        @inbounds for sz in 1:11; consumed += sum(@view fs.cwd[sz, :, :]) * fr[sz]; end
        smoke = consumed * _FM_SMOKE_DEAD + (fs.flive[1] + fs.flive[2]) * _FM_SMOKE_LIVE
        return (; flame, scorch = sch, ba_kill, vol_kill, smoke, models = collect(models))
    end
    return (; severe = scenario(1, 20f0, 70f0, 1), moderate = scenario(3, 8f0, 60f0, 1))
end
