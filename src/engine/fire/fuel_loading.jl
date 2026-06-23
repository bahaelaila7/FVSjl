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

"FFE live-fuel category (1–4) from the FFE forest type IFFEFT (FMCBA, fmcba.f:161)."
@inline function ffe_live_fuel_type(iffeft::Integer)::Int
    if     3 <= iffeft <= 5    1   # pines
    elseif iffeft == 7         3   # redcedar
    elseif iffeft == 6         4   # oak savannah
    else                       2   # hardwoods
    end
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
(FULIV, fmcba.f). Coastal-plain/piedmont stands later override the shrub load from
an understory-age/site-index curve (FULIV2) — not yet ported.
"""
@inline function ffe_live_fuel_loading(coef::SpeciesCoefficients, iffeft::Integer)
    ft = ffe_live_fuel_type(iffeft)
    return (coef.ffe_fuel_live[ft, 1], coef.ffe_fuel_live[ft, 2])
end
