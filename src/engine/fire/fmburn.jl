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
# MOISTURE keyword (fmin.f opt 5): the (date, 7-%) override active for a fire in calendar `year`.
# Date semantics match the activity schedule (apply_compress!): a date ≥ 1000 is a calendar year, a
# date < 1000 is a 1-based cycle number (cycle index + 1). The most-recently-listed match wins (FVS
# OPGET loops the due activities; the last sets MOIS). Returns the 7-tuple of % or nothing.
function _active_moisture_override(s::StandState, year::Int)
    fs = s.fire
    (fs === nothing || isempty(fs.moisture_ovr)) && return nothing
    fvscyc = Int(s.control.cycle) + 1
    hit = nothing
    @inbounds for (date, pr) in fs.moisture_ovr
        d = Int(date)
        (d == year || (0 < d < 1000 && d == fvscyc)) && (hit = pr)
    end
    return hit
end

# FUELTRET (fmusrfm.f): the fuel-bed DEPTH multiplier active for calendar `year` — applied for 5 years after
# the treatment date (date ≥ 1000 = year; < 1000 = cycle, resolved via the cycle length). 1.0 when none.
function _fueltret_dpmod(s::StandState, year::Int)::Float32
    fs = s.fire
    (fs === nothing || isempty(fs.fueltret)) && return 1f0
    dp = 1f0
    @inbounds for (date, d) in fs.fueltret
        dy = Int(date) >= 1000 ? Int(date) :
             Int(s.control.cycle_year[1]) + (Int(date) - 1) * max(1, round(Int, s.control.year))
        (dy <= year <= dy + 5) && (dp = d)
    end
    return dp
end

# Build the 2×5 fuel-moisture matrix (dead 1hr/10hr/100hr/1000hr/duff in row 1; live woody/herb in row 2)
# from the 7 MOISTURE % values, converting % → fraction (fmburn.f:373-380 MOIS(i,j)=MPRMS·.01).
function _moisture_matrix(pr::NTuple{7,Float32})::Matrix{Float32}
    m = zeros(Float32, 2, 5)
    m[1, 1] = pr[1] * 0.01f0; m[1, 2] = pr[2] * 0.01f0; m[1, 3] = pr[3] * 0.01f0
    m[1, 4] = pr[4] * 0.01f0; m[1, 5] = pr[5] * 0.01f0
    m[2, 1] = pr[6] * 0.01f0; m[2, 2] = pr[7] * 0.01f0
    return m
end

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
                 flmult::Float32 = 1f0, crburn::Float32 = 0f0, year::Integer = 0,
                 cyclen::Real = 5f0)::FireResult
    fs = s.fire
    (fs === nothing || !fs.active) && return FireResult(0f0, 0f0, 0f0, 0f0, 0f0)
    fmcba!(s)                                            # fuel pools, cover, percent cover
    t = s.trees; coef = s.coef
    # MOISTURE keyword (fmburn.f:367): a fire in the cycle a MOISTURE activity is scheduled for uses the
    # user's explicit fuel moistures (FMOIS=0 path) instead of the FMMOIS dryness-model table.
    movr = _active_moisture_override(s, Int(year))
    mois = movr === nothing ? fuel_moisture(fmois, s.variant) : _moisture_matrix(movr)
    fwind = wind * fire_wind_reduction(fs.percov)        # 20-ft wind → midflame
    # SN surface fire (FMCFMD + FMDYN + FMFINT): select the weighted standard fuel models
    # for the stand and integrate Rothermel over them, summing the weighted flame & Byram.
    models = select_fuel_models(s, mois; fire_basis = true)   # burn on start-of-cycle + 1-annual-step down wood
    byram = 0f0
    # FVS's FMFINT loops over the selected fuel models and sets FLAG(1)=1 if ANY model's DEAD-fuel moisture
    # damping MDCSA(1) ≤ 0 (dead moisture ≥ that model's moisture of extinction); fmburn.f:473 then SKIPS the
    # entire FMEFF tree-mortality path (the flame/scorch are still reported). rothermel_surface_fire returns
    # byram=0 exactly when a model's mdcsa1 ≤ 0, so `fire_carries` mirrors NOT(FLAG(1)). Confirmed bit-exact
    # vs live: a wet fire (one model at extinction) reports flame but applies zero mortality.
    fire_carries = true
    dpmod = _fueltret_dpmod(s, Int(year))                # FUELTRET fuel-bed depth multiplier (1.0 if none)
    for (fm, w) in models
        load, sav, depth, mext = fuel_model_resolved(s, fm)
        depth *= dpmod
        r = rothermel_surface_fire(load, sav, depth, mext, mois; wind = fwind, slope_tan = s.plot.slope)
        r.byram <= 0f0 && (fire_carries = false)          # this model does not carry → FVS FLAG(1)=1
        byram += r.byram * w
    end
    # FVS recomputes flame from the WEIGHTED Byram (fmfint.f:541, NLC 2003) — NOT the Σ of per-model
    # flames (x^0.46 is concave, so the per-model sum biases low). fmburn.f:439-464: apply the FLAMEADJ
    # multiplier to flame, then back-compute Byram from the modified flame so scorch tracks it.
    flame = byram > 0f0 ? 0.45f0 * fpow(byram / 60f0, 0.46f0) : 0f0     # ^0.46 via gfortran companion (doctrine #8)
    oldfl = flame
    flmult != 1f0 && (flame = oldfl * flmult)
    flame != oldfl && (byram = 60f0 * fpow(flame / 0.45f0, 1f0 / 0.46f0))
    sch = byram > 0f0 ? scorch_height(byram, atemp, fwind) : 0f0

    # pre-fire total live TPA by FVS_Mortality DBH class (LOWDBH bins, 7 non-cumulative classes), both the
    # stand aggregate (the ALL row) and PER-SPECIES (FVS_Mortality emits one row per species + an ALL row).
    totcls = zeros(Float32, 7); clskil = zeros(Float32, 7)
    sp_tot = Dict{Int,Vector{Float32}}(); sp_kil = Dict{Int,Vector{Float32}}()
    sp_bak = Dict{Int,Float32}(); sp_vol = Dict{Int,Float32}()
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        c = _fm_mort_class(t.dbh[i]); c >= 1 || continue
        totcls[c] += t.tpa[i]
        get!(() -> zeros(Float32, 7), sp_tot, Int(t.species[i]))[c] += t.tpa[i]
    end
    killed = 0f0; killed_ba = 0f0; killed_vol = 0f0
    v2t = coef_col(coef, :v2t)
    is_sprout = coef_col(coef, :is_sprouting)            # ESTUMP sprout-species filter (fmkill.f:80 → estump.f)
    if mortcode != 0 && fire_carries                      # FLAG(1) gate: skip mortality if the fire doesn't carry
        # FMEFF brackets its per-tree RANN draws with RANNGET(SAVESO) (fmeff.f:143) … RANNPUT(SAVESO)
        # (fmeff.f:569): the fire's draws are ROLLED BACK, so the fire consumes ZERO NET main-stream RNG.
        # jl must restore too — else the ~ITRN fire draws ADVANCE the stream and desync the POST-fire DGSCOR
        # serial-correlation deviates, making the survivors grow wrong (the kill stays bit-exact — same draws
        # — but the next cycle's growth drifts, ~4.4% Bdft by the 3rd post-fire cycle). D15.
        _fire_rng_save = rannget(s.rng)                   # RANNGET(SAVESO)
        @inbounds for i in 1:t.n
            # FMEFF draws RANN for EVERY record (DO 100 I=1,ITRN, fmeff.f:144/152), UNCONDITIONALLY
            # before any FMPROB/tpa guard. Draw first so the stream count matches live FVS exactly;
            # the FMPROB>0 guard (fmeff.f:176) applies only after the draw.
            (rann!(s.rng) * 100f0 > psburn) && continue  # unburned portion (fmeff.f:159 GOTO 90)
            t.tpa[i] > 0f0 || continue                   # FMPROB>0 guard (fmeff.f:176), post-draw
            csv = crown_volume_scorched(sch, t.height[i], Int(t.crown_pct[i]))
            sp = Int(t.species[i]); d = t.dbh[i]
            pmort = fire_tree_mortality(coef, sp, d, flame, csv, s.variant)
            pmort = fire_mortality_adjust(pmort, sp, d, burnseas, s.variant)
            (d <= 1f0 && csv > 50f0) && (pmort = 1f0)     # fmeff.f:330
            pmort *= active_fmort_mult(s.control, sp, year, d)   # FMORTMLT per-tree multiplier (fmeff.f:340)
            pmort = clamp(pmort, 0f0, 1f0)
            curkil = pmort * t.tpa[i]
            crburn > 0f0 && (curkil += crburn * (t.tpa[i] - curkil))  # crown-fire share
            t.tpa[i] -= curkil
            t.tpa[i] < 0f0 && (t.tpa[i] = 0f0)
            # Fire-killed sprouting trees feed the ESUCKR stump-sprout pool exactly as cutting does: FVS
            # fmkill.f:80 (ICALL=1, in GRADD after the fire) calls ESTUMP(ISP,DBH,FIRKIL,ITRE,ISHAG) for every
            # record with FIRKIL>0.00001 — the SAME pool cuts.f:1713 fills. ISHAG = IY(ICYC+1)−BURNYR (yrs
            # fire→cycle-end = sprout age); live stamp: SIMFIRE is booked at BURNYR=cycle-start ⇒ ISHAG = cyclen
            # (=fint). Gated by LSPRUT (NOAUTOES/NOSPROUT turn it off) + the per-species is_sprouting flag, so
            # NOAUTOES keeps the fire kill but suppresses the post-fire sprouts (the 456→177 TPA difference).
            # esuckr! (simulate.jl, after the fire) then sprouts them. cut_log was freshly cleared by cuts! this
            # cycle and cuts! runs BEFORE the fire, so these records append cleanly alongside any harvest cuts.
            if s.control.lsprut && is_sprout[sp] == 1f0
                push!(s.control.cut_log,
                      (species = Int32(sp), dstmp = d, prem = curkil,
                       plot = Int32(t.plot_id[i]), ishag = round(Int32, cyclen)))
            end
            killed += curkil
            killed_ba += curkil * 0.005454154f0 * d * d   # fire-killed basal area (ft²/ac, fmfout.f:303)
            killed_vol += curkil * t.merch_cuft_vol[i]    # SN: merch cubic volume killed (fmfout.f:306)
            c = _fm_mort_class(d)
            if c >= 1
                clskil[c] += curkil
                get!(() -> zeros(Float32, 7), sp_kil, sp)[c] += curkil
                sp_bak[sp] = get(sp_bak, sp, 0f0) + curkil * 0.005454154f0 * d * d
                sp_vol[sp] = get(sp_vol, sp, 0f0) + curkil * t.merch_cuft_vol[i]
            end
            # Fire-killed trees become standing snags. Carry the MERCH bole (mcf·v2t/2000) — the same basis
            # as ordinary-mortality snags (mortality.jl) and the carbon_snt-validated StandDead/down-wood
            # bole — so the fall transfers a stem-only bole, NOT the jenkins TOTAL-AGB fallback (which
            # double-counts the crown that belongs in the separate CWD2B path) (fmsvol.f merch MCF).
            mcf = max(0.005454154f0 * t.height[i], t.merch_cuft_vol[i])
            add_snag!(fs, sp, d, curkil, year; bolevol = mcf * v2t[sp] / 2000f0, height = t.height[i])
            # Pool the fire-killed CROWN into the crown-debris pool (CWD2B), as FMEFF does for the dead
            # trees. But FIRST consume the fire-REACHED fine crown the way FMEFF does (fmeff.f:457-460)
            # BEFORE it is booked as snags: in the scorched crown zone the fire burns 100% of the foliage
            # (size 0) and 50% of the 0-0.25" branches (size 1, incl. its OLDCRW crown-lift) — those go to
            # the atmosphere (BCROWN released), NOT to down-wood. PROPCR = the scorched fraction of the
            # crown LENGTH (fmeff.f:435 = sl/CRL; the parabolic `csv` used for mortality is a DIFFERENT,
            # volume measure). Tall trees whose crown sits above the scorch height get PROPCR=0 (crown
            # intact — the prior "above the flame" assumption, correct only for them); small trees get
            # PROPCR=1 (foliage gone, size-1 halved). Live-validated per-tree vs FVSsn CROWNW at the fire:
            # sugar-maple d1.28 size-1 0.2725→0.136 (PROPCR 1), beech d6.9 ×0.822 (PROPCR≈0.36). Sizes 2-5
            # are above the flames / too coarse to burn ⇒ unchanged, so the fine down-wood path is intact.
            xc = crown_biomass(s, sp, d, t.height[i], Int(t.crown_pct[i]))
            ol = crown_lift_at_death(t, i, cyclen)             # YRSCYC·OLDCRW (fmscro.f:147)
            crl = t.height[i] * Float32(t.crown_pct[i]) / 100f0
            sl  = crl > 0f0 ? clamp(sch - (t.height[i] - crl), 0f0, crl) : 0f0
            propcr = crl > 0f0 ? sl / crl : 0f0
            ol2 = 0.5f0 * ol[2]                                # fmeff.f:460 ALWAYS halves OLDCRW(1) for fire-killed
            #   trees (inside IF(ICALL.EQ.0), NOT gated on the scorch zone). The other half is burned (fmeff.f:448
            #   BCROWN += 0.5·YRSCYC·OLDCRW(1)). The old `propcr>0 ? 0.5·ol[2] : ol[2]` over-booked the FULL crown-
            #   lift into CWD2B for propcr=0 trees (bark-killed, crown above the scorch) ⇒ StandDead-high (SN/CS/NE
            #   fire crowns are scorched, propcr>0, so were unaffected; LS bark-driven jack-pine kills expose it).
            xvc = (xc[1] * (1f0 - propcr),                     # foliage burned over the scorched length
                   xc[2] * (1f0 - 0.5f0 * propcr) + ol2,       # half the scorched 0-0.25" branches burned
                   xc[3] + ol[3], xc[4] + ol[4], xc[5] + ol[5], xc[6] + ol[6])
            fmscro!(s, sp, d, xvc, curkil, clamp(ffe_dkr_cls(s, sp), 1, 4))  # FUELPOOL-overridable
            # Fire-killed coarse ROOTS → the dead-root pool (BIOROOT, fmsadd.f:320 BIOROOT+=RBIO·SNGNEW·XDCAY).
            # Freshly killed ⇒ XDCAY=(1−CRDCAY)^0=1, same age-0 basis as ordinary mortality (mortality.jl). The
            # snag-FALL path transfers only the BOLE (not roots), so this is the sole root booking (no double-count).
            _, _, rbio = jenkins_biomass(coef, sp, d)
            fs.bioroot += rbio * curkil
        end
        rannput!(s.rng, _fire_rng_save)                   # RANNPUT(SAVESO): roll back the fire's RANN draws
    end
    # the fire consumes a share of the surface fuels — releasing carbon, leaving the rest. The CONSUMED
    # loadings (FVS_Consumption) are the before−after difference in the FFE fuel pools.
    fuel_before = ffe_fuel_loadings(s)
    carbon_released = apply_fire_consumption!(fs, mois)
    fuel_after = ffe_fuel_loadings(s)
    consumed = NamedTuple{keys(fuel_before)}(map(-, values(fuel_before), values(fuel_after)))
    # per-species mortality rows (FVS_Mortality emits one row per present species + the ALL aggregate),
    # sorted by species index for determinism (FMFOUT/dbsfmmort.f).
    species_mort = NamedTuple[]
    @inbounds for sp in sort!(collect(keys(sp_tot)))
        kil = get(sp_kil, sp, zeros(Float32, 7))
        push!(species_mort, (; fvs = strip(coef.code_alpha[sp]), plants = strip(coef.code_plants[sp]),
              fia = strip(coef.code_fia[sp]), clskil = Tuple(kil), totcls = Tuple(sp_tot[sp]),
              bakill = get(sp_bak, sp, 0f0), volkill = get(sp_vol, sp, 0f0)))
    end
    # capture the burn-event record for the FVS_BurnReport / Mortality / Consumption DBS tables
    push!(fs.burn_reports, (; year = Int(year), mois = copy(mois), wind = fwind, flame = flame,
          slope = s.plot.slope, scorch = sch, models = collect(models), killed = killed, killed_ba = killed_ba,
          killed_vol = killed_vol, released = carbon_released,
          clskil = Tuple(clskil), totcls = Tuple(totcls), species_mort = species_mort, consumed = consumed))
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
    # FMPOFL_FMPTRH (fmpofl.f:506/649) does RANNGET(SAVES0) … RANN draws … RANNPUT(SAVES0): the POTENTIAL-fire
    # REPORT is a hypothetical and must NOT consume the simulation RNG (its stochastic torching draws would
    # otherwise shift the crown-ratio stream — a 1-CCF perturbation seen on a diverse FFE stand). Save/restore.
    saved_s0 = s.rng.s0
    pt = torching_probability(s, pf.severe.flame, pf.moderate.flame)
    rannput!(s.rng, saved_s0)
    return (; surf_flame_sev = pf.severe.flame, surf_flame_mod = pf.moderate.flame,
            tot_flame_sev = pf.severe.flame, tot_flame_mod = pf.moderate.flame,   # = surface (no crown fire in SN)
            ptorch_sev = pt.severe, ptorch_mod = pt.moderate,
            torch_index = torching_index(s, cbd.cbd, cbd.actcbh, 1, s.variant),   # OINIT1 (torching, severe fmois=1); −1 for SN
            crown_index = crowning_index(s, cbd.cbd, 1, s.variant),               # OACT1 (crowning, severe fmois=1); −1 for SN
            canopy_ht = cbd.actcbh, canopy_density = cbd.cbd,   # FVS_PotFire Canopy_Ht = ACTCBH (fmpofl.f:302), the crown base
            mort_ba_sev = pf.severe.ba_kill, mort_ba_mod = pf.moderate.ba_kill,
            mort_vol_sev = pf.severe.vol_kill, mort_vol_mod = pf.moderate.vol_kill,
            smoke_sev = pf.severe.smoke, smoke_mod = pf.moderate.smoke,
            models = pf.severe.models)                                            # severe-case fuel models (fmpofl.f:230)
end

# Canopy minimum height (fminit.f:147 CANMHT=6.0): a tree must be taller than this to enter the canopy
# crown-fuel profile (fmpocr.f:80).
const _FM_CANMHT = 6f0

# LSW (fmvinit.f): the "canopy softwood" species that contribute crown fuel to the canopy bulk-density
# profile — HARDWOODS do NOT (fmpocr.f:19,78). The classification is per-variant BLOCK DATA: NE = species
# 1:25 (ne/fmvinit.f:1151-1156); SN = species 1:17 + 88 (sn/fmvinit.f:1011-1014).
fm_canopy_lsw(sp::Integer, ::Southern)  = sp <= 17 || sp == 88
fm_canopy_lsw(sp::Integer, ::Northeast) = sp <= 25
fm_canopy_lsw(sp::Integer, ::AbstractVariant) = sp <= 25

# PotFire severe/moderate scenario wind (mi/h) + temperature (°F): (sev_wind, sev_temp, mod_wind, mod_temp),
# from each variant's fmvinit.f PREWND/POTEMP BLOCK DATA (SN fmvinit.f:63-66 vs NE fmvinit.f:63-66).
potfire_env(::Southern)  = (20f0, 70f0, 8f0, 60f0)
potfire_env(::Northeast) = (25f0, 80f0, 15f0, 50f0)
potfire_env(::AbstractVariant) = (20f0, 70f0, 8f0, 60f0)

# Fuel model 10 (timber litter + understory) — the fixed crown fuel model FMCFIR overlays for the crown-fire
# indices (fmcfir.f:122-133): 3 dead classes + 1 live, loads (lb/ft²) / SAV (1/ft), depth 1, dead MEXT .25.
const _FM10_LOAD = Float32[0.138 0.092 0.23 0.0; 0.092 0.0 0.0 0.0]
const _FM10_SAV  = Float32[2000.0 109.0 30.0 0.0; 1500.0 0.0 0.0 0.0]

"""
    crowning_index(s, cbd, fmois, variant) -> Float32

The Scott & Reinhardt crowning index O'active — the 20-ft wind (mi/h) at which an active crown fire is
sustained (FMCFIR, fmcfir.f:162-168). NE runs FMCFIR (fmpofl.f:167 ELSE branch); SN/CS skip it ⇒ −1.
Computed from the FM10 crown-fuel-model intermediates at the scenario moisture (propagating flux SIRXI =
`xio`, heat sink SRHOBQ = `rhobqig`, slope factor SPHIS = `phis`) and the canopy bulk density `cbd`.
"""
crowning_index(::StandState, ::Float32, ::Int, ::AbstractVariant) = -1f0
function crowning_index(s::StandState, cbd::Float32, fmois::Int, ::Northeast)::Float32
    cbd > 0f0 || return -1f0
    r = rothermel_surface_fire(_FM10_LOAD, _FM10_SAV, 1f0, 0.25f0,
                               fuel_moisture(fmois, s.variant); slope_tan = s.plot.slope)
    r.xio < 1f-5 && return -1f0
    o = ((2.95f0 * r.rhobqig / (r.xio * cbd)) - r.phis - 1f0) / 0.001612f0
    return o > 0f0 ? o^0.7f0 * 0.01137f0 / 0.4f0 : 0f0
end

"""
    torching_index(s, cbd, actcbh, fmois, variant) -> Float32

The Scott & Reinhardt torching index O'init — the 20-ft wind (mi/h) at which crown fire INITIATES
(FMCFIR, fmcfir.f:197-271). NE only; SN/CS ⇒ −1. The critical surface spread rate for torching
RINIT1 = 60·INIT1/HPA (INIT1 from the foliar-moisture/crown-base-height ladder rule, HPA = the stand's
heat-per-area = `Σxir·w·384/Σsigma·w`); then BISECT the 20-ft wind (× the canopy reduction WMULT) until
the stand's WEIGHTED surface-fuel-model spread = RINIT1. NB the torching bisection uses the FMCFMD weighted
STAND models (fmfint.f:120-134, the ICALL=2 ELSE branch) — NOT the fixed FM10 the crowning index uses.
"""
torching_index(::StandState, ::Float32, ::Integer, ::Int, ::AbstractVariant) = -1f0
function torching_index(s::StandState, cbd::Float32, actcbh::Integer, fmois::Int, ::Northeast)::Float32
    (cbd > 0f0 && actcbh >= 0) || return -1f0
    mois = fuel_moisture(fmois, s.variant)
    models = select_fuel_models(s, mois)
    # HPA = stand heat-per-area = Σxir·w·384/Σsigma·w (fmfint.f:550, wind-independent intermediates)
    sxir = 0f0; ssig = 0f0
    for (fm, w) in models
        r = rothermel_surface_fire(fuel_model_resolved(s, fm)..., mois; slope_tan = s.plot.slope)
        sxir += r.xir * w; ssig += r.sigma * w
    end
    (ssig > 0f0 && sxir > 0f0) || return -1f0
    hpa = sxir * 384f0 / ssig
    folmc = 100f0                                     # foliar moisture content (fminit.f:150 default)
    init1 = ((460f0 + 25.9f0 * folmc) * 0.001333f0 * Float32(actcbh))^1.5f0
    rinit1 = 60f0 * init1 / hpa
    wmult = fire_wind_reduction(s.fire.percov)
    # weighted-model surface spread (ft/min) at a 20-ft wind `oi` (canopy-reduced to midflame `oi·wmult`)
    spr(oi) = sum(rothermel_surface_fire(fuel_model_resolved(s, fm)..., mois;
                  wind = oi * wmult, slope_tan = s.plot.slope).spread * w for (fm, w) in models)
    spr(999f0) < rinit1 && return 999f0               # never reaches the critical rate ⇒ cap at 999
    lo = 0f0; hi = 999f0; o = 0f0
    for _ in 1:1000                                   # bisection (fmcfir.f:237 DO 200)
        o = (lo + hi) / 2f0; d = spr(o) - rinit1
        abs(d) <= 0.001f0 && break
        d > 0.001f0 ? (hi = o) : (lo = o)
    end
    return min(o, 999f0)
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
        # Tree-inclusion filter (fmpocr.f:78-80): canopy-softwood species (LSW; hardwoods excluded), crown
        # ratio > 0 (FMICR), and height > CANMHT. Without it the profile picks up hardwood + understory crown
        # ⇒ CBD too high / crown base too low (the FVS_PotFire Canopy_Density + crown-fire index error).
        h = t.height[i]; h > _FM_CANMHT || continue
        fm_canopy_lsw(Int(t.species[i]), s.variant) || continue
        icr = Float32(t.crown_pct[i]); icr > 0f0 || continue
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
    function scenario(sev::Int, fmois::Int, wind::Float32, temp::Float32, season::Int)
        # POTF* keyword overrides for this severity (sev 1=SEVERE, 2=MODERATE); −1/0 ⇒ scenario default.
        pc = fs.params.potf[sev]
        mois = pc.mois[1] >= 0f0 ? _moisture_matrix(pc.mois) : fuel_moisture(fmois, s.variant)
        pc.wind >= 0f0   && (wind = pc.wind)
        pc.temp >= 0f0   && (temp = pc.temp)
        pc.season > 0    && (season = Int(pc.season))
        fwind = wind * fire_wind_reduction(fs.percov)
        models = select_fuel_models(s, mois)
        dpmod = _fueltret_dpmod(s, Int(current_cycle_year(s)))   # FUELTRET depth multiplier
        byram = 0f0
        for (fm, w) in models
            load, sav, depth, mext = fuel_model_resolved(s, fm)
            depth *= dpmod
            r = rothermel_surface_fire(load, sav, depth, mext, mois; wind = fwind, slope_tan = s.plot.slope)
            byram += r.byram * w
        end
        flame = byram > 0f0 ? 0.45f0 * (byram / 60f0)^0.46f0 : 0f0   # fmfint.f:541 — flame from weighted Byram
        sch = byram > 0f0 ? scorch_height(byram, temp, fwind) : 0f0
        ba_kill = 0f0; vol_kill = 0f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            csv = crown_volume_scorched(sch, t.height[i], Int(t.crown_pct[i]))
            sp = Int(t.species[i]); d = t.dbh[i]
            pm = fire_tree_mortality(coef, sp, d, flame, csv, s.variant)
            pm = fire_mortality_adjust(pm, sp, d, season, s.variant)
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
        # POTFPAB: % area burned scales the potential kill + smoke (FMEFF/FMCONS take POTPAB, default 100 =
        # full, which equals jl's unscaled values; an override < 100 reduces them).
        if pc.pab >= 0f0
            f = pc.pab * 0.01f0; ba_kill *= f; vol_kill *= f; smoke *= f
        end
        return (; flame, scorch = sch, ba_kill, vol_kill, smoke, models = collect(models))
    end
    # PotFire scenario wind (mi/h) + temperature (°F) are per-variant BLOCK DATA (fmvinit.f:63-66
    # PREWND/POTEMP): SN severe 20/70°F + moderate 8/60°F; NE severe 25/80°F + moderate 15/50°F.
    sw, st, mw, mt = potfire_env(s.variant)
    return (; severe = scenario(1, 1, sw, st, 1), moderate = scenario(2, 3, mw, mt, 1))
end
