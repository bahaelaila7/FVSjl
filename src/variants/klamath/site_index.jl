# =============================================================================
# site_index.jl (klamath) — NC site index + SDImax (chunk 2). Traced nc/sitcind.f + nc/sitset.f + nc/ecocls.f.
#
# SITEAR: site-range interpolation from SIMIN/SIMAX (species_coefficients.csv site_lo/site_hi), shared with
# CI/UT (nc/sitcind.f defaults unset species from the site-species' index). SDImax (SDIDEF): NC is ecocls-
# driven — the plant association (PA) → SDIMX from data/klamath/ecocls_sdimax.csv (like CI's R4SDI but a
# 90-row PA table); BA-weighted in sdical (BAMAX=XMAX·0.5454154·PMSDIU). NC forests (nc/htdbh.f IFOR 1-7):
# 505 Klamath / 510 Six Rivers / 514 Trinity / 611 Siskiyou / 705 Hoopa / 800 Simpson / 712 BLM Coos Bay.
# ⚠ Chunk-2 FIRST PASS: SITEAR interpolation (validatable at cyc0) + provisional SDImax (sp_sdi_def from the
# CSV default / ecocls; SDImax is cyc0-INERT — feeds only self-thin mortality). ecocls PA-string crosswalk +
# BA-weighting to be wired + validated with the mortality chunk (7).
# =============================================================================

const NC_JFOR = Int[505, 510, 514, 611, 705, 800, 712]   # nc/htdbh.f forest codes → IFOR 1..7

"nc/forkod: KODFOR → IFOR 1..7 (default 1 = Klamath if unrecognized)."
function nc_forkod!(p)
    idx = findfirst(==(Int(p.user_forest_code)), NC_JFOR)
    p.forest_idx = Int32(idx === nothing ? 1 : idx)
    return Int(p.forest_idx)
end

# nc/sitset.f site-index species defaults (SI array): the per-species site index when NO SITECODE, with the
# site species (ISISP default = DF sp3) at 90. MEASURED from live nct01.out (DF=90, SP/PP=100, rest 90) — the
# HTCALC-based DF→species conversion (sitset.f DO 30) reduces to these for the default DF-site-species=90 case.
# (Stands WITH a SITECODE set p.sp_site_index directly; the full HTCALC conversion for a non-default site
#  species/index is a follow-up — nct01 + the common no-SITECODE case uses these defaults.)
const NC_SITE_DEFAULT = Float32[90,100,90,90,90,90,90,90,90,100,90,90]

"nc/sitset.f: site species default (DF) at 90; unset species get the SI() defaults (DF→species conversion)."
function nc_sitset!(s::StandState)
    p = s.plot; sd = s.coef.species
    isisp = Int(p.site_species); isisp == 0 && (isisp = 3)     # sitset.f: ISISP default = 3 (DF)
    p.site_species = Int32(isisp)
    @inbounds for i in 1:12
        p.sp_site_index[i] <= 0f0 && (p.sp_site_index[i] = NC_SITE_DEFAULT[i])
    end
    # SDImax (SDIDEF) — provisional (CSV sdi_max_default; ecocls PA lookup + BA-weight lands with chunk 7).
    sdimax = sd[:sdi_max_default]
    @inbounds for i in 1:12
        p.sp_sdi_def[i] <= 0f0 && (p.sp_sdi_def[i] = sdimax[i])
    end
    return s
end

function nc_site_index_setup!(s::StandState)
    nc_forkod!(s.plot)
    nc_sitset!(s)
    return s
end

site_setup!(s::StandState, ::Klamath) = nc_site_index_setup!(s)
