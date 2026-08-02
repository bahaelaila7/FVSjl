# =============================================================================
# site_index.jl (bluemountains) — BM site index + SDImax + forest code. Chunk 2.
# (bm/sitset.f, bm/ecocls.f, bm/habtyp.f, bm/sichg.f, bm/htcalc.f, bm/forkod.f)
#
# BM's site chain is UNLIKE the R4 variants (UT/TT interpolate a [SITELO,SITEHI] range).
# BM is Region-6: the habitat code (KODTYP) indexes a 92-entry plant-association list
# (PCOML, bm/habtyp.f) → a PA code → the 283-row ECOCLS eco-class DB (bm/ecocls.f), which
# supplies each eco-class species' (SDImax, SI, site-species flag). The IFLAG=1 species is
# the SITE SPECIES (ISISP); its SI passes through directly (SITEAR fills only unset entries),
# and SICHG+HTCALC translate an equivalent SI onto the OTHER species via age-curve inversion.
#
# ★ VERIFIED end-to-end vs live bmt01: STDINFO forest 614 (Umatilla→IFOR 3) + habitat 12 →
# PCOML[12]='CDS722' → ECOCLS → site species DF(3), SITEAR(DF)=64, SDIDEF(DF)=346.
# =============================================================================

# --- bm/habtyp.f PCOML (92 KODTYP→PA-code) + bm/ecocls.f ECOCLS (283 rows), loaded from CSV. ---
function _bm_load_site_tables()
    pcoml = String[]
    for l in readlines(joinpath(BM_DATADIR, "pcoml.csv"))[2:end]
        push!(pcoml, String(strip(split(l, ',')[2])))
    end
    # ecocls rows in file order: (pa, spc, fvsseq, sdimx, site, numbr, iflag)
    rows = NamedTuple{(:pa,:spc,:fvsseq,:sdimx,:site,:numbr,:iflag),
                      Tuple{String,String,Int,Float32,Float32,Int,Int}}[]
    for l in readlines(joinpath(BM_DATADIR, "ecocls.csv"))[2:end]
        f = split(strip(l), ',')
        push!(rows, (pa=String(f[1]), spc=String(f[2]), fvsseq=parse(Int, f[3]),
                     sdimx=parse(Float32, f[4]), site=parse(Float32, f[5]),
                     numbr=parse(Int, f[6]), iflag=parse(Int, f[7])))
    end
    return pcoml, rows
end
const BM_PCOML, BM_ECOCLS = _bm_load_site_tables()

# bm/sichg.f — linear age-at-breast-height coefficients (SIAGE per species) + reference-age/type.
# Species order WP WL DF GF MH WJ LP ES AF PP WB LM PY YC AS CW OS OH.
const BM_SICHG_A = Float32[18.43316, 8.11668, 8.50000, 9.01840, 45.24969, 34.0, 10.65724,
                           33.72545, 13.88200, 8.0, 34.0, 34.0, 11.56252, 11.56252, 34.0,
                           11.56252, 8.0, 11.56252]
const BM_SICHG_B = Float32[-0.13399, -0.05661, -0.05000, -0.05700, -1.28885, -0.20, -0.10667,
                           -0.274509, -0.06588, -0.04286, -0.20, -0.20, -0.05586, -0.05586,
                           -0.20, -0.05586, -0.04286, -0.05586]
const BM_SICHG_REFLOC = ['T','B','B','B','B','T','T','B','B','B','T','T','B','B','B','B','B','B']
const BM_SICHG_REFAGE = Float32[50,50,50,50,100,100,50,100,100,100,100,100,100,100,80,100,100,100]

# bm/forkod.f — KODFOR → IFOR (location subscript) + IGL. JFOR=[604,607,614,616,619], all KFOR=1.
const BM_JFOR = Int[604, 607, 614, 616, 619]
function bm_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 1; useigl = true
    if kodfor == 8117           # Umatilla Reservation → Umatilla (bm/forkod.f CASE 8117)
        ifor = 3
    else
        idx = findfirst(==(kodfor), BM_JFOR)
        idx === nothing ? (useigl = false; ifor = 1) : (ifor = idx)
    end
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(1))    # IGL = KFOR(IFOR) = 1
    return ifor
end

# bm/habtyp.f — KODTYP (numeric habitat code) → PA code via PCOML. KODTYP≤0 or >NPA → "" (default).
function bm_habtyp(kodtyp::Integer)::String
    (kodtyp <= 0 || kodtyp > length(BM_PCOML)) && return ""
    return BM_PCOML[kodtyp]
end

# bm/ecocls.f — return the ECOCLS rows for plant-association code `pa` (an eco-class = consecutive rows).
bm_ecocls(pa::AbstractString) = filter(r -> r.pa == pa, BM_ECOCLS)

# bm/sichg.f — SIAGE(I): age each species' reference reaches, given the site species (ISISP) SI = SSITE.
function bm_sichg(isisp::Integer, ssite::Float32)
    isiloc = BM_SICHG_REFLOC[isisp]
    siage = Vector{Float32}(undef, length(BM_SICHG_A))
    @inbounds for i in eachindex(siage)
        diff = 0
        (isiloc == 'T' && BM_SICHG_REFLOC[i] == 'B') && (diff = -1)
        (isiloc == 'B' && BM_SICHG_REFLOC[i] == 'T') && (diff = 1)
        age2bh = 0f0
        if diff != 0
            age2bh = BM_SICHG_A[i] + BM_SICHG_B[i] * ssite
            (isisp != 5 && i == 5) && (age2bh = BM_SICHG_A[i] + BM_SICHG_B[i] * (ssite / 3.281f0))
            (isisp == 5 && i != 5) && (age2bh = BM_SICHG_A[i] + BM_SICHG_B[i] * (ssite * 3.281f0))
        end
        siage[i] = BM_SICHG_REFAGE[i] + age2bh * diff
    end
    return siage
end

# bm/htcalc.f — height (HGUESS) at age AG on species ISPC's site curve, given site index SINDX.
# 18 curves verbatim. ★ CASE(3) DF carries a Fortran ** precedence bug in the LAST term:
# term B = (1-EXP(..))**0.966998  (correct), term C = 1-EXP(..)**0.966998  (EXP raised first) —
# replicated faithfully (bm/htcalc.f:69-70).
function bm_htcalc(sindx::Float32, ispc::Integer, ag::Float32)::Float32
    lag = log(ag)
    if ispc == 1
        return sindx / (0.37504453f0 * (1f0 - 0.92503f0 * exp(-0.0207959f0 * ag))^(-2.4881068f0))
    elseif ispc == 2
        q = -0.12528f0 + 0.039636f0*ag - 0.0004278f0*ag*ag + 1.7039f-6*ag^3
        return 4.5f0 + 1.46897f0*ag + 0.0092466f0*ag*ag - 0.00023957f0*ag^3 +
               1.1122f-6*ag^4 + (sindx-4.5f0)*q - 73.57f0*q
    elseif ispc == 3
        # DF — faithful precedence bug: term C uses EXP(..)**p (not (1-EXP)**p).
        termB = -0.2828f0 + 1.87947f0*(1f0 - exp(-0.022399f0*ag))^0.966998f0
        termC = -0.2828f0 + 1.87947f0*(1f0 - exp(-0.022399f0*ag)^0.966998f0)
        return 4.5f0 + exp(-0.37496f0 + 1.36164f0*lag - 0.00243434f0*lag^4) -
               79.97f0*termB + (sindx-4.5f0)*termC
    elseif ispc == 4
        x2 = -0.30935f0 + 1.2383f0*lag + 0.001762f0*lag^4 - 5.4f-6*lag^9 +
             2.046f-7*lag^11 - 4.04f-13*lag^18
        x3 = -6.2056f0 + 2.097f0*lag - 0.09411f0*lag^2 - 0.00004382f0*lag^7 +
             2.007f-11*lag^16 - 2.054f-17*lag^24
        return exp(x2) - 84.73f0*exp(x3) + (sindx-4.5f0)*exp(x3) + 4.5f0
    elseif ispc == 5
        h = (22.8741f0 + 0.950234f0*sindx)*(1f0 - exp(-0.00206465f0*sqrt(sindx)*ag))^
            (1.365566f0 + 2.045963f0/sindx)
        return (h + 1.37f0)*3.281f0
    elseif ispc == 6
        return 0f0
    elseif ispc == 7
        return sindx*(-0.0968f0 + 0.02679f0*ag - 0.00009309f0*ag*ag)
    elseif ispc == 8
        return 4.5f0 + (2.75780f0*sindx^0.83312f0)*(1f0 - exp(-0.015701f0*ag))^
               (22.71944f0*sindx^(-0.63557f0))
    elseif ispc == 9
        return sindx*(-0.07831f0 + 0.0149f0*ag - 4.0818f-5*ag*ag)
    elseif ispc == 10 || ispc == 17
        b = -0.7864f0 + 2.49717f0*(1f0 - exp(-0.0045042f0*ag))^0.33022f0
        return 128.8952205f0*(1f0 - exp(-0.016959f0*ag))^1.23114f0 - b*100.43f0 +
               b*(sindx-4.5f0) + 4.5f0
    elseif ispc == 11 || ispc == 12 || ispc == 15
        return 0f0
    elseif ispc == 13 || ispc == 14 || ispc == 16 || ispc == 18
        return (sindx-4.5f0)/(0.6192f0 - 5.3394f0/(sindx-4.5f0) + 240.29f0*ag^(-1.4f0) +
               (3368.9f0/(sindx-4.5f0))*ag^(-1.4f0)) + 4.5f0
    else
        return 0f0
    end
end

# bm/sitset.f — set site species + SITEAR + SDIDEF from the eco-class (ECOCLS), then translate SI.
const BM_FORMAX = 850f0
function bm_sitset!(s::StandState)
    p = s.plot
    maxsp = nspecies(s.variant)
    nsiset = count(>(0f0), @view p.sp_site_index[1:maxsp])   # #species SI-set by keyword
    pcom = bm_habtyp(Int(p.habitat_code))
    rows = bm_ecocls(pcom)
    isisp = 0; jsisp = 0
    for r in rows
        iseq = r.fvsseq
        iseq == 0 && continue
        rsdi = min(r.sdimx, BM_FORMAX)
        (jsisp == 0 && r.iflag == 1) && (jsisp = iseq)
        (isisp <= 0 && r.iflag == 1) && (isisp = iseq)
        (p.sp_site_index[iseq] <= 0f0 && nsiset == 0) && (p.sp_site_index[iseq] = r.site)
        p.sp_sdi_def[iseq] <= 0f0 && (p.sp_sdi_def[iseq] = rsdi)
        if isisp > 0 && r.iflag == 1 && rsdi <= 0f0 && p.sp_sdi_def[isisp] <= 0f0
            p.sp_sdi_def[isisp] = rsdi
        end
    end
    isisp <= 0 && (isisp = 10)                                # default site species = PP (bm/sitset.f:142)
    p.sp_site_index[isisp] <= 0f0 && (p.sp_site_index[isisp] = 70f0)
    p.site_species = Int32(isisp)

    # SITERANGE for the site species (site_lo/site_hi loaded in species_coefficients.csv).
    site_lo = s.coef.species[:site_lo]; site_hi = s.coef.species[:site_hi]
    slossp = site_lo[isisp]; shissp = site_hi[isisp]

    si = Vector{Float32}(undef, maxsp)
    if isisp == 6 || isisp == 11 || isisp == 12 || isisp == 15
        # UT/TT site species — alternate [SLO,SHI]-range interpolation.
        tem = p.sp_site_index[isisp] > 0f0 ? p.sp_site_index[isisp] : 30f0
        tem < slossp && (tem = slossp)
        @inbounds for i in 1:maxsp
            slo = site_lo[i]; shi = site_hi[i]
            si[i] = slo + (tem - slossp)/(shissp - slossp)*(shi - slo)
            i == 5 && (si[i] = min(si[i]/3.281f0, 28f0))
        end
    else
        siage = bm_sichg(isisp, p.sp_site_index[isisp])
        sindx = p.sp_site_index[isisp]
        @inbounds for ispc in 1:maxsp
            slo = site_lo[ispc]; shi = site_hi[ispc]
            si[ispc] = bm_htcalc(sindx, isisp, siage[ispc])
            ispc == 5 && (si[ispc] = min(si[ispc]/3.281f0, 28f0))
            if ispc == 6 || ispc == 11 || ispc == 12 || ispc == 15
                tem = max(sindx, slossp)
                p.sp_site_index[ispc] <= 0f0 &&
                    (si[ispc] = slo + (tem - slossp)/(shissp - slossp)*(shi - slo))
            end
        end
    end
    @inbounds for i in 1:maxsp
        p.sp_site_index[i] == 0f0 && (p.sp_site_index[i] = si[i])
    end
    # SDIDEF fill (bm/sitset.f:236-244): default to site-species SDIDEF, capped at FORMAX.
    k = p.sp_sdi_def[isisp] > 0f0 ? isisp : jsisp
    @inbounds for i in 1:maxsp
        p.sp_sdi_def[i] > 0f0 && continue
        k > 0 && (p.sp_sdi_def[i] = min(p.sp_sdi_def[k], BM_FORMAX))
    end
    return s
end

function bm_site_index_setup!(s::StandState)
    bm_forkod!(s.plot)
    bm_sitset!(s)
    return s
end

site_setup!(s::StandState, ::BlueMountains) = bm_site_index_setup!(s)
