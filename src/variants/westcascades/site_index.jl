# =============================================================================
# site_index.jl (westcascades) — WC site index + Reineke SDImax + forest code. Chunk 2.
# (wc/forkod.f, wc/habtyp.f, wc/ecocls.f, wc/sichg.f, wc/htcalc.f, wc/sitset.f)
#
# WC's site chain is Region-6 (same family as BM, template): the habitat code (KODTYP) indexes the
# 139-entry plant-association list (PCOML, wc/habtyp.f) → a PA code → the 139-row ECOCLS eco-class
# table (wc/ecocls.f, all single-species NUMBR=1), which supplies the site species (IFLAG=1), its
# SDImax and site index. SICHG converts an equivalent reference age per species; wc_htcalc (chunk 4,
# ALREADY validated 24/24 bit-exact) evaluates the SITE-SPECIES height curve at each species' ref age
# to fan SITEAR to every species; then wc/sitset.f's misc-hardwood reduction factors apply.
#
# ★ VALIDATED end-to-end vs live FVSwc_clean wct01 (DEBUG SITSET): STDINFO forest 618 (Willamette→
# IFOR 6) + habitat 52 → PCOML[52]='CFS551' → ECOCLS → site species DF(16), SITEAR(DF)=73,
# SDIDEF=815. All 39 SITEAR bit-exact vs the live "AFTER SITE ADJUSTMENT FACTORS" dump; all SDI MAX
# = 815. See docs/WC_VARIANT_PORT_AUDIT.md chunk 2.
# =============================================================================

# --- wc/habtyp.f PCOML (139 KODTYP→PA code) + wc/ecocls.f ECOCLS (139 single-species rows), from CSV. ---
function _wc_load_site_tables()
    pcoml = String[]
    for l in readlines(joinpath(WC_DATADIR, "pcoml.csv"))[2:end]
        isempty(strip(l)) && continue
        push!(pcoml, String(strip(split(l, ',')[2])))
    end
    rows = NamedTuple{(:pa,:spc,:fvsseq,:sdimx,:site,:numbr,:iflag),
                      Tuple{String,String,Int,Float32,Float32,Int,Int}}[]
    for l in readlines(joinpath(WC_DATADIR, "ecocls.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(strip(l), ',')
        push!(rows, (pa=String(f[1]), spc=String(f[2]), fvsseq=parse(Int, f[3]),
                     sdimx=parse(Float32, f[4]), site=parse(Float32, f[5]),
                     numbr=parse(Int, f[6]), iflag=parse(Int, f[7])))
    end
    return pcoml, rows
end
const WC_PCOML, WC_ECOCLS = _wc_load_site_tables()

# wc/forkod.f — accepted FVS location codes (KODFOR) → IFOR 1..11 subscript. Reservation pseudo-codes
# 8124→MtBaker(2), 8130→GiffordPinchot(1); IFOR 11 (613) remaps to MtBaker-Snoqualmie(2).
const WC_JFOR = Int[603, 605, 606, 610, 615, 618, 708, 709, 710, 711, 613]

function wc_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 1; useigl = true
    if kodfor == 8124
        ifor = 2                                  # Sauk-Suiattle Reservation → Mt Baker (605)
    elseif kodfor == 8130
        ifor = 1                                  # Yakama Nation Reservation → Gifford Pinchot (603)
    else
        idx = findfirst(==(kodfor), WC_JFOR)
        idx === nothing ? (useigl = false; ifor = 1) : (ifor = idx)
    end
    ifor == 11 && (ifor = 2)                      # 613 Mt Baker-Snoqualmie → Mt Baker-Snoqualmie (605)
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(1))         # IGL = KFOR(IFOR) = 1 (all locations)
    p.user_forest_code = Int32(WC_JFOR[ifor])     # wc/forkod.f: KODFOR = JFOR(IFOR)
    return ifor
end

# wc/habtyp.f — KODTYP (numeric habitat code) → PA code via PCOML. Out-of-range → "" (default CFS551).
function wc_habtyp(kodtyp::Integer)::String
    (kodtyp <= 0 || kodtyp > length(WC_PCOML)) && return ""
    return WC_PCOML[kodtyp]
end

# wc/ecocls.f — return the ECOCLS rows for plant-association code `pa` (single species here).
wc_ecocls(pa::AbstractString) = filter(r -> r.pa == pa, WC_ECOCLS)

# wc/sichg.f — SIAGE(I): reference age each species reaches, given the site species (ISISP) SI = SSITE.
# Reproduces the live "CALLING HTCALC … AG" per species (e.g. LP(11) 52.87 at DF-site 73).
function wc_sichg(s::StandState, isisp::Integer, ssite::Float32)
    sd = s.coef.species
    a = sd[:sichg_a]; b = sd[:sichg_b]; refage = sd[:sichg_refage]; refloc = sd[:sichg_refloc]
    maxsp = nspecies(s.variant)
    isiloc = refloc[isisp]                        # 0='B' (breast-hgt), 1='T' (total age)
    siage = Vector{Float32}(undef, maxsp)
    @inbounds for i in 1:maxsp
        diff = 0
        (isiloc == 1f0 && refloc[i] == 0f0) && (diff = -1)
        (isiloc == 0f0 && refloc[i] == 1f0) && (diff = 1)
        age2bh = 0f0
        if diff != 0
            age2bh = a[i] + b[i] * ssite
            (isisp != 20 && i == 20) && (age2bh = a[i] + b[i] * (ssite / 3.281f0))   # MH ref units
            (isisp == 20 && i != 20) && (age2bh = a[i] + b[i] * (ssite * 3.281f0))
        end
        siage[i] = refage[i] + age2bh * Float32(diff)
    end
    return siage
end

# wc/sitset.f — set site species + SITEAR + SDIDEF from the eco-class, then translate SI to every species.
const WC_FMSDI = Float32[950,950,900,850,825,870,885,870,825,850]   # wc/sitset.f FMSDI(IFOR) forest SDI cap

function wc_sitset!(s::StandState)
    p = s.plot; sd = s.coef.species
    maxsp = nspecies(s.variant)
    formax = WC_FMSDI[clamp(Int(p.forest_idx), 1, length(WC_FMSDI))]
    nsiset = count(>(0f0), @view p.sp_site_index[1:maxsp])
    pcom = wc_habtyp(Int(p.habitat_code))
    isempty(pcom) && (pcom = "CFS551")            # wc/habtyp.f default PA (ITYPE=52 → PCOML[52])
    rows = wc_ecocls(pcom)

    isisp = (1 <= Int(p.site_species) <= maxsp) ? Int(p.site_species) : 0
    @inbounds for r in rows
        iseq = r.fvsseq
        iseq == 0 && continue
        rsdi = min(r.sdimx, formax)                            # wc/sitset.f: cap default SDI at FORMAX
        (isisp <= 0 && r.iflag == 1) && (isisp = iseq)
        (p.sp_site_index[iseq] <= 0f0 && nsiset == 0) && (p.sp_site_index[iseq] = r.site)
        p.sp_sdi_def[iseq] <= 0f0 && (p.sp_sdi_def[iseq] = rsdi)
    end
    isisp <= 0 && (isisp = 10)                                 # wc/sitset.f:111 Region-6 default site sp = ES(10)
    p.sp_site_index[isisp] <= 0f0 && (p.sp_site_index[isisp] = 70f0)   # wc/sitset.f:112
    p.site_species = Int32(isisp)

    # SICHG → wc_htcalc fan-out: evaluate the SITE-SPECIES curve at each species' reference age, then apply
    # the wc/sitset.f misc-hardwood reduction factors (site_redux) + the MH(20)/WO(28) special transforms.
    sindx = p.sp_site_index[isisp]
    siage = wc_sichg(s, isisp, sindx)
    redux = sd[:site_redux]
    si = Vector{Float32}(undef, maxsp)
    @inbounds for ispc in 1:maxsp
        if p.sp_site_index[ispc] > 0f0               # SITECODE / site species already set — keep it
            si[ispc] = p.sp_site_index[ispc]
            continue
        end
        v = wc_htcalc(sindx, isisp, siage[ispc])
        if ispc == 20 && isisp != 20                          # MH: metric→ft, cap 28
            v = v / 3.281f0
            v > 28f0 && (v = 28f0)
        elseif ispc == 28 && isisp != 28                      # WO: Gould max-height transform (King DF SI)
            v = 114.2f0 * (1f0 - exp(-0.0266f0 * v))^2.26f0
        elseif ispc != isisp && redux[ispc] != 1f0
            v = v * redux[ispc]                               # misc-hardwood SI reductions (sp 21,23-27,29,31,33-37)
        end
        si[ispc] = v
    end
    @inbounds for i in 1:maxsp
        p.sp_site_index[i] <= 0f0 && (p.sp_site_index[i] = si[i])
    end
    # SDIDEF fill (wc/sitset.f:168-170): unset species default to the site-species SDIDEF.
    @inbounds for i in 1:maxsp
        p.sp_sdi_def[i] <= 0f0 && (p.sp_sdi_def[i] = p.sp_sdi_def[isisp])
    end
    return s
end

function wc_site_index_setup!(s::StandState)
    wc_forkod!(s.plot)
    wc_sitset!(s)
    return s
end

site_setup!(s::StandState, ::WestCascades) = wc_site_index_setup!(s)
