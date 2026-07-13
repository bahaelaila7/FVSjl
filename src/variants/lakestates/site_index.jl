# =============================================================================
# site_index.jl (lakestates) — LS SITSET site-index fan-out (ls/sitset.f)
#
# LS is NE-family here: it converts the input site index (for the site species
# ISISP) to a site index for EVERY species through a full MAXSP×MAXSP linear model
#   SITEAR(I) = SICOEF1(ISISP,I) + SICOEF2(ISISP,I)·SITEAR(ISISP)              (:230)
# then a SECOND pass keyed on aspen (species 41) for any species still unset      (:233)
#   SITEAR(I) = SICOEF1(41,I) + SICOEF2(41,I)·SITEAR(41)
# and a final fallback SITEAR(I) = SITEAR(ISISP) for anything still ≤ .0001       (:239).
#   default site species = 3 (RN, red pine); default SI = 60                     (:215-221)
# SDIDEF (:242-245): per-species max SDI = BAMAX-derived if a BAMAX/SDIMAX keyword
# set it, else SDICON(I) (ls/sitset.f DATA SDICON → species_coefficients.sdi_max_default).
#
# ls/forkod.f maps the user forest code to an LS forest index IFOR (JFOR = 902,903,904,
# 906,907,909,910,913,924) and sets default TLAT/TLONG/ELEV per IFOR (grinit default IFOR=5).
# =============================================================================

# SICOEF1/SICOEF2 68×68 matrices, loaded once from data/lakestates/site_sicoef.csv.
# The CSV holds the Fortran DATA flat order ((C(I,J),I=1,68),J=1,68) — I-fastest — so
# reshape(v,68,68) reproduces M[I,J] == SICOEF(I,J) (Julia is column-major, matching Fortran).
const _LS_SICOEF = let
    path = joinpath(LS_DATADIR, "site_sicoef.csv")
    c1 = Float32[]; c2 = Float32[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        f = split(strip(line), ',')
        push!(c1, parse(Float32, f[1])); push!(c2, parse(Float32, f[2]))
    end
    (reshape(c1, 68, 68), reshape(c2, 68, 68))   # M[I,J] = SICOEF(I,J)
end

"""
    ls_site_index_setup!(s)

SITSET (ls/sitset.f): fill `plot.sp_site_index[sp]` for every species from the site
species' index via the per-species SICOEF1/SICOEF2 linear model (+ aspen second pass +
fallback), then set `plot.sp_sdi_def` from SDICON (unless a BAMAX/SDIMAX keyword set it).
"""
function ls_site_index_setup!(s::StandState)
    p = s.plot; v = s.variant
    C1, C2 = _LS_SICOEF
    sea = p.sp_site_index
    isisp = Int(p.site_species)
    if isisp <= 0                                   # ls/sitset.f:215 default site species = 3 (RN)
        isisp = 3; p.site_species = Int32(3)
    end
    sea[isisp] <= 0f0 && (sea[isisp] = 60f0)        # :219 default SI = 60
    sisp = sea[isisp]
    @inbounds for i in 1:nspecies(v)                # DO 5 (:229-231)
        sea[i] <= 0.0001f0 && (sea[i] = C1[isisp, i] + C2[isisp, i] * sisp)
    end
    if sea[41] > 0.0001f0                            # DO 10 (:233-237) — aspen second pass
        s41 = sea[41]
        @inbounds for i in 1:nspecies(v)
            sea[i] <= 0.0001f0 && (sea[i] = C1[41, i] + C2[41, i] * s41)
        end
    end
    @inbounds for i in 1:nspecies(v)                # DO 15 (:239-241) — fallback
        sea[i] < 0.0001f0 && (sea[i] = sisp)
    end

    # SDIDEF (:242-245). LS uses PMSDIU as a PERCENT (÷100); default SDICON when no BAMAX.
    bamax = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 85f0
    sdicon = s.coef.species[:sdi_max_default]
    @inbounds for i in 1:nspecies(v)
        if p.sp_sdi_def[i] <= 0f0
            p.sp_sdi_def[i] = bamax > 0f0 ? bamax / (0.5454154f0 * (pmsdiu / 100f0)) : sdicon[i]
        end
    end
    return s
end

# ls/forkod.f — user forest code → LS forest index IFOR; the standard national forests are
# JFOR(1..9). Alias/old-region codes (7109/7502/…) also fold into these; ported as needed.
const _LS_JFOR = Int32[902, 903, 904, 906, 907, 909, 910, 913, 924]
# Per-IFOR default lat/long/elev (ls/forkod.f:295-320 SELECT CASE(IFOR); 1/4/8 share, others distinct).
const _LS_FOR_DEFAULTS = Dict{Int,NTuple{3,Float32}}(
    1 => (45.93f0, 90.44f0, 15f0), 4 => (45.93f0, 90.44f0, 15f0), 8 => (45.93f0, 90.44f0, 15f0),
    2 => (47.38f0, 94.60f0, 13f0),
    3 => (44.25f0, 85.40f0,  9f0),
    5 => (46.45f0, 90.17f0, 14f0),
    6 => (46.78f0, 92.11f0, 16f0),
    7 => (45.75f0, 87.06f0,  8f0),
    # ls/forkod.f:297-322 SELECT CASE(IFOR) has NO CASE(9): forest 924 (Manistee, JFOR index 9) gets NO geo
    # default — TLAT/TLONG/ELEV fall through unset (grinit ELEV=0). The prior `9 => IFOR-5 values` was a
    # misread of the fall-through and set ELEV=14 ⇒ inflated Hopkins index ⇒ open-grown crown width ⇒ CCF
    # (audit 43ec: 18447951010661 CCF ~15% high vs live). Omit key 9 ⇒ the fallback below leaves geo unset.
)

"""
    ls_forkod_defaults!(s)

FORKOD (ls/forkod.f): map the user forest code to an LS forest index (`IFOR`) and set the
default TLAT/TLONG/ELEV for that forest when the stand has none. A code matching JFOR(i)
maps to IFOR=i; an unrecognized/absent code keeps the grinit default IFOR=5.
"""
function ls_forkod_defaults!(s::StandState)
    p = s.plot
    code = Int32(p.user_forest_code)
    ifor = 5                                          # ls/grinit.f:186 default
    idx = findfirst(==(code), _LS_JFOR)
    idx !== nothing && (ifor = idx)
    d = get(_LS_FOR_DEFAULTS, ifor, (0f0, 0f0, 0f0))   # unmapped IFOR (e.g. 9) ⇒ no default, matches FVS fall-through
    lat0, long0, elev0 = d
    p.latitude  == 0f0 && (p.latitude  = lat0)
    p.longitude == 0f0 && (p.longitude = long0)
    p.elevation == 0f0 && (p.elevation = elev0)
    p.forest_idx = Int32(ifor)
    return s
end

function site_setup!(s::StandState, ::LakeStates)
    ls_forkod_defaults!(s)
    ls_site_index_setup!(s)
end
