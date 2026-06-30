# =============================================================================
# fire/fuel_loading.jl — initial surface fuel loading (FFE chunk F3-data)
#
# Ported from: bin/FVSsn_buildDir/fmcba.f (FMCBA fuel initialization).
#
# FMCBA initializes the stand's surface fuel pools from its forest type:
#   - dead fuels (1hr…1000hr down wood + litter + duff) by FIA forest type, and
#   - live herb/shrub fuels by a coarser live-fuel category.
# The loading tables (FUINI 11×9, FULIV 2×4) live in data/southern/fire_fuel_{dead,live}.csv;
# the two pure classification maps from forest type → fuel-table column are here.
#
# The 11 dead-fuel size classes (MXFLCL, FMPARM.F77): <0.25", 0.25–1", 1–3", 3–6",
# 6–12", 12–20", 20–35", 35–50", >50", litter, duff.
# =============================================================================

"FFE dead-fuel forest-type group (1–9) from the FIA forest-type code (FMCBA, fmcba.f:260)."
@inline function ffe_dead_fuel_type(ifortp::Integer)::Int
    if     101 <= ifortp <= 105                       1   # eastern white pine
    elseif ifortp == 141 || ifortp == 142             2   # longleaf–slash pine
    elseif 161 <= ifortp <= 168                        3   # loblolly–shortleaf pine
    elseif ifortp == 181 || ifortp == 402              4   # eastern redcedar
    elseif ifortp == 401 || (403 <= ifortp <= 409)     5   # oak–pine
    elseif 501 <= ifortp <= 520                        6   # oak–hickory
    elseif 601 <= ifortp <= 608                        7   # oak–gum–cypress
    elseif 701 <= ifortp <= 709                        8   # elm–ash–cottonwood
    elseif 801 <= ifortp <= 809                        9   # maple–beech–birch
    else                                               6   # default → oak–hickory
    end
end

"""
    ffe_forest_type(s) -> Int

The FFE categorical forest type IFFEFT (1–9) for the stand (FMSNFT, fmsnft.f),
used to set default live surface fuels (F3) and in the fuel-model logic (F4):
1 hardwood · 2 hardwood/pine · 3 pine/hardwood · 4 pine · 5 pine-bluestem ·
6 oak-savannah · 7 redcedar · 8 St-Francis · 9 nonstocked. The mixed pine/hardwood
classes split on the stand's pine basal-area fraction (SN species 4–14 are pines);
the savannah/bluestem classes also use top height (ATAVH = AVHT40) and stocking class.
"""
function ffe_forest_type(s::StandState)::Int
    ifortp = Int(s.plot.forest_type)
    t = s.trees
    if ifortp in (997, 504, 505, 510, 512, 515, 519, 520)
        return 1                                       # hardwood
    elseif ifortp in (103, 104, 141, 142, 996, 401, 409) ||
           (161 <= ifortp <= 168) || (403 <= ifortp <= 407)
        pineba = 0f0; npineba = 0f0                     # pine vs non-pine basal area
        @inbounds for i in 1:t.n
            x = t.tpa[i] * t.dbh[i] * t.dbh[i] * 0.0054542f0
            (4 <= t.species[i] <= 14) ? (pineba += x) : (npineba += x)
        end
        tot = pineba + npineba
        iffeft = tot > 0f0 ? (f = pineba / tot; f <= 0.50f0 ? 2 : f <= 0.70f0 ? 3 : 4) : 0
        (ifortp == 162 && stand_top_height(s) > 50f0 && Int(s.plot.stocking_class) >= 3) && (iffeft = 5)
        return iffeft
    elseif ifortp == 501 || ifortp == 503
        return (stand_top_height(s) > 30f0 && Int(s.plot.stocking_class) >= 3) ? 6 : 1
    elseif ifortp == 181 || ifortp == 402
        return 7                                       # eastern redcedar
    elseif ifortp in (602, 605, 701, 706, 708, 807)
        return 8                                       # St. Francis types
    elseif ifortp == 999
        return 9                                       # nonstocked
    else
        return 1                                       # default → hardwood
    end
end

"FFE live-fuel category (1–4) from the FFE forest type IFFEFT (FMCBA, fmcba.f:161)."
@inline function ffe_live_fuel_type(iffeft::Integer)::Int
    if     3 <= iffeft <= 5    1   # pines
    elseif iffeft == 7         3   # redcedar
    elseif iffeft == 6         4   # oak savannah
    else                       2   # hardwoods
    end
end

# FMNEFT (ne/fmneft.f): NE FFE forest type (IFFEFT, 1–11) from the FIA forest type.
@inline function _ne_iffeft(ifortp::Integer)::Int
    if     101 <= ifortp <= 105   1   # white/red/jack pine (+ hemlock)
    elseif 121 <= ifortp <= 127   2   # spruce-fir
    elseif 161 <= ifortp <= 168   3   # loblolly-shortleaf
    elseif 381 <= ifortp <= 383   4   # exotic softwoods
    elseif 401 <= ifortp <= 409   5   # oak-pine
    elseif 501 <= ifortp <= 520   6   # oak-hickory
    elseif 601 <= ifortp <= 608   7   # oak-gum-cypress
    elseif 701 <= ifortp <= 709   8   # elm-ash-cottonwood
    elseif 801 <= ifortp <= 809   9   # maple-beech-birch
    elseif 901 <= ifortp <= 904   10  # aspen-birch
    elseif ifortp == 999          11  # nonstocked
    else                          2   # default → spruce-fir
    end
end

"""
    ne_dead_fuel_loading(s) -> NTuple{11,Float32}

NE initial dead surface fuel loading (FUINI, ne/fmcba.f:214): row `FTDEADFU = (IFFEFT−1)·3 + ISZCL`
of the 31-group `fire_fuel_dead.csv` (IFFEFT=11 nonstocked → 31). IFFEFT from FMNEFT; ISZCL is the
STKVAL stand size class (`plot.size_class`, 1=saw/2=pole/3=seedling). NE live fuel is constant (0.31,0.31).
"""
@inline function ne_dead_fuel_loading(s::StandState)
    iffeft = _ne_iffeft(Int(s.plot.forest_type))
    iszcl = clamp(Int(s.plot.size_class), 1, 3)
    ft = iffeft == 11 ? 31 : (iffeft - 1) * 3 + iszcl
    row = @view s.coef.ffe_fuel_dead[ft, :]
    return ntuple(i -> row[i], 11)
end

"""
    ffe_dead_fuel_loading(coef, ifortp) -> NTuple{11,Float32}

Initial dead surface fuel loading (tons/acre) for a stand of FIA forest type `ifortp`,
by the 11 FFE size classes (FUINI, fmcba.f). The first nine are down-wood size classes,
then litter and duff.
"""
@inline function ffe_dead_fuel_loading(coef::SpeciesCoefficients, ifortp::Integer)
    ft = ffe_dead_fuel_type(ifortp)
    row = @view coef.ffe_fuel_dead[ft, :]
    return ntuple(i -> row[i], 11)
end

"""
    ffe_live_fuel_loading(coef, iffeft) -> (herb, shrub)

Initial live herb and shrub fuel loading (tons/acre) for FFE forest type `iffeft`
(FULIV, fmcba.f). Coastal-plain/piedmont stands override the shrub load from an
understory-age/site-index curve (FULIV2) — see `ffe_live_fuel_override` (ported).
"""
@inline function ffe_live_fuel_loading(coef::SpeciesCoefficients, iffeft::Integer)
    ft = ffe_live_fuel_type(iffeft)
    return (coef.ffe_fuel_live[ft, 1], coef.ffe_fuel_live[ft, 2])
end

# FULIV2 (sn/fmcba.f:82-89): the coastal-plain / piedmont / mountain SHRUB-load override. For
# ecological units 232* / 231* / M221* the herb/shrub live fuel is read from a understory-rough-age
# × site-index-class curve (Southern Forestry Smoke Management Guidebook, GTR-SE-10 p118) instead of
# the flat FULIV table — all biomass goes to SHRUB, herbs are 0. Rows = the 8 rough ages in `_FULIV2_Y`;
# columns = 6 site-index classes (<50, 50-65, 65-80, 80-95, 95-110, ≥110). (Small fixed table — a
# distinct override, kept inline rather than in the live-fuel CSV.)
const _FULIV2_Y = (1f0, 2f0, 3f0, 5f0, 7f0, 10f0, 15f0, 20f0)          # rough age (years)
const _FULIV2 = Float32[                                                # [age 1:8, si-class 1:6]
    0.4 1.2 2.6 4.5 7.0 10.0
    0.4 1.3 2.6 4.5 7.0 10.0
    0.5 1.3 2.7 4.6 7.0 10.0
    0.6 1.5 2.8 4.7 7.2 10.2
    0.9 1.7 3.1 5.0 7.4 10.4
    1.4 2.2 3.5 5.5 7.9 10.9
    2.6 3.4 4.7 6.6 9.1 12.1
    4.2 5.1 6.4 8.3 10.8 13.8]

"ALGSLP (algslp.f): segmented-linear interpolation of `y` at `x` over ascending knots `xs`, clamped flat beyond the ends."
@inline function _ffe_algslp(x::Float32, xs::NTuple{8,Float32}, y::AbstractVector{Float32})
    x <= xs[1] && return y[1]
    x >= xs[8] && return y[8]
    @inbounds for k in 2:8
        if x <= xs[k]
            f = (x - xs[k-1]) / (xs[k] - xs[k-1])
            return y[k-1] + f * (y[k] - y[k-1])
        end
    end
    return y[8]
end

"""
    ffe_live_fuel_override(s) -> (herb, shrub) | nothing

The FULIV2 coastal-plain/piedmont/mountain shrub-load override (sn/fmcba.f:122-154). Returns the
(herb, shrub) live-fuel loading (tons/ac) for ecological units 232*/231*/M221* — herb is always 0,
shrub from the rough-age × site-index curve, ×0.40 for piedmont/mountain (231*/M221*). `nothing` for
other units (use the flat `ffe_live_fuel_loading`). Rough age = years since the last burn, or
(current − inventory year + 5) when unburned, clamped to [1, 20].
"""
function ffe_live_fuel_override(s::StandState)
    eu = s.plot.eco_unit
    (startswith(eu, "232") || startswith(eu, "231") || startswith(eu, "M221")) || return nothing
    si = s.plot.sp_site_index[s.plot.site_species]
    j = si < 50f0 ? 1 : si < 65f0 ? 2 : si < 80f0 ? 3 : si < 95f0 ? 4 : si < 110f0 ? 5 : 6
    # Rough age = years since the LAST ACTUAL burn (FULIV2 resets the understory rough age after a fire),
    # or inventory-based when nothing has burned yet. `fire.fire_year` is the SCHEDULED fire year — using
    # it pre-burn makes (iyr − fire_year) negative → clamps the age to 1 → wrong (too-young) shrub load
    # (231Dd: jl 0.48 vs live 0.60/0.88). Derive the actual last-burn from the accumulated burn_reports
    # (0 before any fire), the same fix as the fire-basis (sm,lg) timing.
    burnyr = (s.fire !== nothing && !isempty(s.fire.burn_reports)) ?
             maximum(Int(br.year) for br in s.fire.burn_reports) : 0
    iyr = current_cycle_year(s); iy1 = Int(s.control.cycle_year[1])
    age = burnyr > 0 ? iyr - burnyr : iyr - iy1 + 5
    age = clamp(Float32(trunc(age)), 1f0, 20f0)
    shrub = _ffe_algslp(age, _FULIV2_Y, @view _FULIV2[:, j])
    (startswith(eu, "231") || startswith(eu, "M221")) && (shrub *= 0.40f0)
    return (0f0, shrub)
end
