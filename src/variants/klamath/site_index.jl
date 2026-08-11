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

"nc/sitcind.f: SITEAR site-range default (fill unset species from the site-species' index) + provisional SDImax."
function nc_sitset!(s::StandState)
    p = s.plot; sd = s.coef.species
    slo_a = sd[:site_lo]; shi_a = sd[:site_hi]
    isisp = Int(p.site_species)
    tem = 50.0f0
    (isisp > 0 && p.sp_site_index[isisp] > 0f0) && (tem = p.sp_site_index[isisp])
    isisp == 0 && (isisp = 3)                             # default DF (sp3) like CI
    slossp = slo_a[isisp]; shissp = shi_a[isisp]
    @inbounds for i in 1:12
        tem < slossp && (tem = slossp)
        slo = slo_a[i]; shi = shi_a[i]
        if p.sp_site_index[i] <= 0f0
            p.sp_site_index[i] = slo + (tem - slossp) / (shissp - slossp) * (shi - slo)
        end
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
