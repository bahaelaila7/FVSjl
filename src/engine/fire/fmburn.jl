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
function fmburn!(s::StandState; atemp::Float32 = 70f0, wind::Float32 = 20f0, fmois::Integer = 1,
                 psburn::Float32 = 100f0, mortcode::Integer = 1, burnseas::Integer = 1,
                 flmult::Float32 = 1f0, crburn::Float32 = 0f0, year::Integer = 0)::FireResult
    fs = s.fire
    (fs === nothing || !fs.active) && return FireResult(0f0, 0f0, 0f0, 0f0)
    fmcba!(s)                                            # fuel pools, cover, percent cover
    t = s.trees; coef = s.coef
    mois = fuel_moisture(fmois)
    fwind = wind * fire_wind_reduction(fs.percov)        # 20-ft wind → midflame
    load, sav, depth, mext = build_dynamic_fuel_model(s, mois)
    r = rothermel_surface_fire(load, sav, depth, mext, mois; wind = fwind)
    flame = r.flame * flmult
    # scorch height (Van Wagner) from the Byram intensity
    sch = r.byram > 0f0 ? scorch_height(r.byram, atemp, fwind) : 0f0

    killed = 0f0
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
            add_snag!(fs, sp, d, curkil, year)         # fire-killed trees become standing snags
        end
    end
    return FireResult(killed, flame, r.byram, sch)
end
