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
    # sp2 is a pre-stripped key (species.code2 is rstripped at load; structure_stage passes strip()) ⇒
    # no per-call rstrip/String copy. `haskey`+`getindex` instead of `get(…, nothing)` because the latter
    # returns a `Union{Tuple{String,String},Nothing}` whose non-isbits Tuple gets BOXED (32 B/call, per-tree
    # per-cycle); the two-lookup form allocates nothing and is value-identical ⇒ bit-exact.
    if haskey(coef.crown_species, sp2)
        pair = coef.crown_species[sp2]
        eqnum = iwho == 1 ? pair[2] : pair[1]
        if haskey(coef.crown_eqs, eqnum)
            cw = _cw_eval(coef.crown_eqs[eqnum], Float32(d), Float32(cr),
                          hopkins_index(lat, long, elev))
        end
    end
    cw < 0.5f0 && (cw = 0.5f0)
    cw > 99.9f0 && (cw = 99.9f0)
    return cw
end

# -----------------------------------------------------------------------------
# CWIDTH (base/cwidth.f): the per-record CRWDTH(I) array. FVS fills it only at load (fvs.f:207, after CRATET) and at the
# end of every cycle (gradd.f:254, after DENSE/CROWN); TRIPLE copies it to a record's copies and COMPRESS averages it.
# Its consumers in between — SSTAGE (sstage.f:238/276 WK6=CRWDTH(I)) and FFE FMCBA (CWIDTH=CRWDTH(I)) — therefore see the
# dims AND the stand BA (the Crookston BAREA term) of that last CWIDTH call, not a thin's residual BA or a SIMFIRE seam's
# grown small trees. Stored in `t.crown_width` for the variants whose consumers read it.
# -----------------------------------------------------------------------------
# AK/KT/WS/CA run the same base cwidth.f at fvs.f:207 / gradd.f:254: a THINBBA's post-thin StrClass row read the
# residual-BA crown widths in jl (FVSak_g16 10709344010497 2006 Removal_Code 1: cover 23 live / 20 jl, top species YC/MH
# live / RA/YC jl; KT 4718785010690, CA 23742358010900 likewise), and FMCBA's crown biomass the same.
"Variants whose SSTAGE/FMCBA read the stored CRWDTH(I) that `cwidth!` fills."
# TT/UT/CI/CR/NC/SO/EC/WC/PN: the same base/cwidth.f CRWDTH(I) array feeds their sstage.f WK6 and fmcba.f CWIDTH (MEASURED
# FVSpn_g16 504512112126144 THINBBA 2027: the after-thin StrClass cover 21 live from the pre-thin-BA CRWDTH, 22 jl from a
# recomputed post-thin-BA width; FVScr_clean 46279527020004 SIMFIRE 68 -> 24 cells from FMCBA's PERCOV).
# IE/EM (core-resid), AK/KT/WS/CA (west-kcwa-2) and the nine above (west-shared-4) all read the stored array.
_stored_crwdth(v) = v isa InlandEmpire || v isa EasternMontana || v isa SoutheastAlaska || v isa Kootenai ||
                    v isa WestSierra || v isa CentralCalifornia || v isa PacificNorthwest || v isa WestCascades ||
                    v isa SouthCentralOregon || v isa CentralRockies || v isa CentralIdaho || v isa EastCascades ||
                    v isa Utah || v isa Teton || v isa Klamath || v isa BlueMountains ||
                    v isa BritishColumbia   # metric/base/cwidth.f, called at the same fvs.f:207 / gradd.f:254 points

"""
    cwidth!(s) -> s

CWIDTH (cwidth.f): CRWDTH(I) for every live record from its current DBH/HT/ICR and the current stand BA (the forest-grown
`tree_crwdth`, the value FVS_TreeList reports). Called at the two FVS CWIDTH points; no-op outside `_stored_crwdth`.
"""
function cwidth!(s::StandState)
    _stored_crwdth(s.variant) || return s
    t = s.trees
    @inbounds for i in 1:t.n
        t.crown_width[i] = tree_crwdth(s, Int(t.species[i]), t.dbh[i], t.height[i], t.crown_pct[i])
    end
    return s
end

"The CRWDTH(I) a consumer reads: the stored CWIDTH value when set, else computed from the record's current dims."
@inline stored_crwdth(s::StandState, i::Integer) =
    s.trees.crown_width[i] > 0f0 ? s.trees.crown_width[i] :
    tree_crwdth(s, Int(s.trees.species[i]), s.trees.dbh[i], s.trees.height[i], s.trees.crown_pct[i])
