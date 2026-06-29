# =============================================================================
# fire/fmcba.jl — FFE per-cycle fuel & cover-type update (FFE chunk F3-state)
#
# Ported from: bin/FVSsn_buildDir/fmcba.f (FMCBA, the deterministic body).
#
# Each cycle, FMCBA establishes the stand's fire/fuels context:
#   - the cover type (the species carrying the most basal area),
#   - percent canopy cover (from per-tree crown areas),
#   - the live herb/shrub surface fuel (by FFE forest type), and
#   - in the first FFE year, the dead surface fuel pools, split into decay classes
#     by each species' share of stand basal area.
# Results live in the per-stand `FireState` (no globals). Gated on `fire.active`, so
# it is a no-op for non-FFE stands; FMCBA itself changes nothing the `.sum` reports
# (it feeds the fire-behavior / consumption chunks F5–F8).
# =============================================================================

"""
    fmcba!(s) -> StandState

Update the stand's `FireState` cover type, percent cover, big DBH, live fuels, and
(first FFE year) dead-fuel pools (FMCBA, fmcba.f). No-op unless FFE is active.
"""
function fmcba!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    t = s.trees; coef = s.coef
    nsp = length(coef_col(coef, :dbh_min))            # MAXSP

    # live herb/shrub fuels are re-set every year. NE (ne/fmcba.f:68) uses a flat constant
    # FULIV=(0.31,0.31) for all forest types; SN uses the FFE-forest-type table with a
    # coastal-plain/piedmont/mountain rough-age × site-index override (FULIV2) where it applies.
    if s.variant isa Northeast
        fs.flive = (0.31f0, 0.31f0)
    else
        ovr = ffe_live_fuel_override(s)
        fs.flive = ovr === nothing ? ffe_live_fuel_loading(coef, ffe_forest_type(s)) : ovr
    end

    # per-species basal area, total crown area (for percent cover), and the big DBH
    tba = zeros(Float32, nsp)
    totcra = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); d = t.dbh[i]
        tba[sp] += 3.14159f0 * (d / 24f0) * (d / 24f0) * t.tpa[i]
        d > fs.bigdbh && (fs.bigdbh = d)
        sp2 = s.species.class_codes[sp, 1][1:2]       # forest-grown crown width (CWCALC, iwho=0)
        cw = crown_width(coef, sp2, d, t.height[i], Float32(t.crown_pct[i]), 0,
                         s.plot.latitude, s.plot.longitude, s.plot.elevation)
        totcra += 3.1415927f0 * cw * cw / 4f0 * t.tpa[i]
    end

    # cover type = the species with the most basal area; total BA for the decay split
    bamost = 0f0; covtyp = Int32(0); totba = 0f0
    @inbounds for ksp in 1:nsp
        tba[ksp] > bamost && (bamost = tba[ksp]; covtyp = Int32(ksp))
        totba += tba[ksp]
    end
    # no basal area: default to red oak (75) the first year, else keep the previous cover type
    covtyp == 0 && (covtyp = fs.covtyp == 0 ? Int32(75) : fs.covtyp)
    fs.covtyp = covtyp
    fs.percov = (1f0 - exp(-totcra / 43560f0)) * 100f0

    # dead fuels: loaded once (first FFE year), distributed into decay classes by the species BA share
    # (fmcba.f:375-393). The "hard" (J=2) column comes from ffe_dead_fuel_loading; the "soft" (J=1) column
    # is 0 by default. FUELINIT (hard) / FUELSOFT (soft) override per-size-class values (STFUEL, fmcba.f:320-371).
    # IDC = each species' decay-rate class (DKRCLS).
    if !fs.fuels_init
        deffuel = s.variant isa Northeast ? ne_dead_fuel_loading(s) :
                  ffe_dead_fuel_loading(coef, Int(s.plot.forest_type))
        ovh = fs.params.stfuel_hard; ovs = fs.params.stfuel_soft
        fill!(fs.cwd, 0f0)
        @inbounds for isz in 1:11
            sh = (length(ovh) >= isz && ovh[isz] >= 0f0) ? ovh[isz] : deffuel[isz]   # FUELINIT hard
            ss = (length(ovs) >= isz && ovs[isz] >= 0f0) ? ovs[isz] : 0f0            # FUELSOFT soft
            if totba > 0f0
                for ksp in 1:nsp
                    tba[ksp] > 0f0 || continue
                    idc = ffe_dkr_cls(s, ksp)               # FUELPOOL-overridable decay-rate class
                    prcl = tba[ksp] / totba
                    fs.cwd[isz, 2, idc] += prcl * sh
                    ss > 0f0 && (fs.cwd[isz, 1, idc] += prcl * ss)
                end
            else
                idc = ffe_dkr_cls(s, covtyp)               # FUELPOOL-overridable decay-rate class
                fs.cwd[isz, 2, idc] += sh
                ss > 0f0 && (fs.cwd[isz, 1, idc] += ss)
            end
        end
        fs.fuels_init = true
    end
    return s
end
