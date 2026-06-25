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
    models = select_fuel_models(s, mois)
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
            add_snag!(fs, sp, d, curkil, year)         # fire-killed trees become standing snags
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
