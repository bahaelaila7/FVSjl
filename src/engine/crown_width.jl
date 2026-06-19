# =============================================================================
# crown_width.jl — open-grown / forest-grown crown width (CWCALC)
#
# Ported from: base/cwcalc.jl (the eastern-US crown-width equation library).
#
# The ~144 species/equation variants reduce to FOUR formula families; the
# per-equation coefficients and the species→equation map now live in
# data/<variant>/crown_width_{equations,species}.csv (loaded into
# `coef.crown_eqs` / `coef.crown_species`). `iwho==1` selects the open-grown
# equation; `hi` is the Hopkins bioclimatic index (Bechtold 2003), which needs
# stand lat/long/elevation. Used by the crown competition factor (CCF).
#
# Families (Dc = min(D, dbh_cap); small trees scaled toward 0 below a threshold):
#   bechtold : a + b·Dc + c·Dc² + cr_coef·CR + hi_coef·HI   (threshold 5", clamp)
#   bragg    : a + b·D^power                                  (no scaling, no clamp)
#   ek       : a + b·D^power                                  (threshold 3", clamp)
#   smith    : (a + b·Dcm + c·Dcm²)·3.28084,  Dcm = D·2.54    (threshold 3", clamp)
# =============================================================================

"Hopkins bioclimatic index from stand latitude/longitude/elevation (Bechtold 2003)."
@inline function hopkins_index(lat::Real, long::Real, elev::Real)::Float32
    hilong = -abs(Float32(long))
    hielev = Float32(elev) * 100f0
    return ((hielev - 887f0) / 100f0) * 1f0 +
           (Float32(lat) - 39.54f0) * 4f0 +
           (-82.52f0 - hilong) * 1.25f0
end

"Evaluate one crown-width equation family for DBH `d`, crown ratio `cr`, Hopkins `hi`."
@inline function _cw_eval(e::CrownWidthEq, d::Float32, cr::Float32, hi::Float32)::Float32
    fam = e.family
    if fam === :bragg
        return e.a + e.b * d ^ e.power
    end
    # bechtold + braggm use a 5" small-tree floor; ek/smith use 3"
    thr = (fam === :bechtold || fam === :braggm) ? 5f0 : 3f0
    x = d >= thr ? d : thr
    local v::Float32
    if fam === :bechtold
        xc = min(x, e.dbh_cap)
        v = e.a + e.b * xc + e.c * xc * xc + e.cr_coef * cr + e.hi_coef * hi
    elseif fam === :ek || fam === :braggm
        v = e.a + e.b * x ^ e.power
    else  # :smith
        dcm = x * 2.54f0
        v = (e.a + e.b * dcm + e.c * dcm * dcm) * 3.28084f0
    end
    d < thr && (v *= d / thr)
    (e.max_cw > 0f0 && v > e.max_cw) && (v = e.max_cw)
    return v
end

"""
    crown_width(coef, sp2, d, h, cr, iwho, lat, long, elev) -> Float32

Crown width (ft) for a tree of species 2-char code `sp2`, DBH `d`, height `h`,
crown ratio `cr`. `iwho==1` → open-grown crown. Unknown species → 0.5. Clamped to
[0.5, 99.9]. Equation data comes from `coef.crown_species` / `coef.crown_eqs`.
"""
function crown_width(coef::SpeciesCoefficients, sp2::AbstractString, d::Real, h::Real,
                     cr::Real, iwho::Integer, lat::Real, long::Real, elev::Real)::Float32
    cw = 0f0
    pair = get(coef.crown_species, rstrip(String(sp2)), nothing)
    if pair !== nothing
        eqnum = iwho == 1 ? pair[2] : pair[1]
        e = get(coef.crown_eqs, eqnum, nothing)
        e !== nothing && (cw = _cw_eval(e, Float32(d), Float32(cr),
                                        hopkins_index(lat, long, elev)))
    end
    cw < 0.5f0 && (cw = 0.5f0)
    cw > 99.9f0 && (cw = 99.9f0)
    return cw
end
