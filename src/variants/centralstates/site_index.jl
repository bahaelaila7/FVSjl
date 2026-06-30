# =============================================================================
# site_index.jl (centralstates) — CS SITSET site-index fan-out (cs/sitset.f)
#
# CS converts the input site index (for the site species ISISP) to a site index
# for EVERY species through a per-species LINEAR model SITEAR(I) = ASITE(I) +
# BSITE(I)·SITEAR(47), where species 47 (white oak) is the reference. The site
# species' index is first mapped to the white-oak index by inverting its own
# (ASITE,BSITE). This is structurally different from NE's 28×28 SICOEF matrix and
# SN's master-group SITSET, so CS gets its own method.
#
#   default site species = 47 (white oak), default SI = 65   (cs/sitset.f:101-102)
#   SITEAR(47) = -ASITE(ISISP)/BSITE(ISISP) + (1/BSITE(ISISP))·SITEAR(ISISP)  (:106)
#   SITEAR(I)  = ASITE(I) + BSITE(I)·SITEAR(47)   for every unset species      (:113)
# =============================================================================

# ASITE/BSITE per-species coefficients, loaded once from data/centralstates/site_coef.csv.
const _CS_SITE_COEF = let
    path = joinpath(CS_DATADIR, "site_coef.csv")
    asite = zeros(Float32, 96); bsite = zeros(Float32, 96)
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        f = split(strip(line), ',')
        i = parse(Int, f[1])
        asite[i] = parse(Float32, f[3]); bsite[i] = parse(Float32, f[4])
    end
    (asite, bsite)
end

"""
    cs_site_index_setup!(s)

SITSET (cs/sitset.f): fill `plot.sp_site_index[sp]` for every species from the site
species' index via the per-species linear ASITE/BSITE model (white-oak-referenced).
"""
function cs_site_index_setup!(s::StandState)
    p = s.plot; v = s.variant
    asite, bsite = _CS_SITE_COEF
    sea = p.sp_site_index
    isisp = Int(p.site_species)
    isisp <= 0 && (isisp = 47; p.site_species = Int32(47))
    sea[isisp] <= 0f0 && (sea[isisp] = 65f0)
    # Map the site species' index to the white-oak (47) reference index.
    if isisp != 47
        sea[47] = (-asite[isisp] / bsite[isisp]) + (1f0 / bsite[isisp]) * sea[isisp]
    end
    wo = sea[47]
    @inbounds for i in 1:nspecies(v)
        sea[i] < 0.0001f0 && (sea[i] = asite[i] + bsite[i] * wo)
    end

    # SDIDEF (sitset.f): per-species max SDI = SDICON unless BAMAX (or a keyword) set it.
    bamax = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    sdicon = s.coef.species[:sdi_max_default]
    @inbounds for i in 1:nspecies(v)
        if p.sp_sdi_def[i] <= 0f0
            p.sp_sdi_def[i] = bamax > 0f0 ? bamax / (0.5454154f0 * pmsdiu) : sdicon[i]
        end
    end
    return s
end

# CS forest codes (cs/forkod.f JFOR) → IFOR index, and per-IFOR default lat/long/elev.
const _CS_JFOR = Int32[905, 908, 912, 911]              # Mark Twain, Shawnee, Hoosier, Wayne-Hoosier
const _CS_FOR_DEFAULTS = Dict(                           # forkod.f:173-186 (IFOR 4 shares IFOR 3)
    1 => (37.95f0, 91.77f0, 10f0),
    2 => (37.74f0, 88.54f0,  4f0),
    3 => (38.86f0, 86.49f0,  6f0),
)

"""
    cs_forkod_defaults!(s)

FORKOD (cs/forkod.f): map the user forest code to a CS forest index (`IFOR`) via `JFOR`,
remap Wayne-Hoosier (IFOR 4 → Hoosier 3), then set the default TLAT/TLONG/ELEV for that
forest when the stand has none (needed by the crown-width CCF Hopkins index). An
unrecognized/absent code keeps the grinit default `IFOR=1` (Mark Twain — 37.95/91.77/10).
Tribal pseudo-codes also map to IFOR=1, identical to the default fallthrough.
"""
function cs_forkod_defaults!(s::StandState)
    p = s.plot
    code = Int32(p.user_forest_code)
    ifor = 1                                              # cs/grinit.f:186 default
    idx = findfirst(==(code), _CS_JFOR)
    idx !== nothing && (ifor = idx)
    ifor == 4 && (ifor = 3)                               # Wayne-Hoosier (911) → Hoosier (912)
    lat0, long0, elev0 = _CS_FOR_DEFAULTS[ifor]
    p.latitude  == 0f0 && (p.latitude  = lat0)
    p.longitude == 0f0 && (p.longitude = long0)
    p.elevation == 0f0 && (p.elevation = elev0)
    p.forest_idx = Int32(ifor)
    return s
end

function site_setup!(s::StandState, ::CentralStates)
    cs_forkod_defaults!(s)
    cs_site_index_setup!(s)
end
