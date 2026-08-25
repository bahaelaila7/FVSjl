# =============================================================================
# site_index.jl (southcentraloregon) — SO site index + forest code + SICHG fan. Chunk 2.
# (so/forkod.f, so/sichg.f, so/htcalc.f, so/sitset.f)
#
# SO's R6 site chain mirrors EC: forkod maps the FVS location code to IFOR (and, for R6 forests IFOR≤3 or
# =10, resets LZEIDE=.FALSE. ⇒ Reineke); the site species' SITEAR is fanned to every unset species via
# SICHG (per-species reference age SIAGE) + so_htcalc (the site-species height curve inverted to SIAGE).
# so_htcalc is a per-species SELECT CASE of published site curves (Brickell WP / Cochran DF-WF-WL /
# Dolph-red-fir + Dunning-Levitan SH-WO / Dahms LP / Alexander ES / Herman NF / Wiley WH / Curtis misc /
# Harrington RA / Barrett SP-PP). sot01 gives SI=60 via STDINFO (NSISET>0) so the ECOCLS default path is
# inert; the ECOCLS/SDIDEF default-lookup (for no-SI stands) is deferred (sot01 SDImax rides the species
# sdimax column via stand_sdimax, Reineke).
# =============================================================================

# so/forkod.f — location code → IFOR via the JFOR table (11 R6-Oregon/NE-Cal forests).
const SO_JFOR = Int[601, 602, 620, 505, 506, 509, 511, 701, 514, 799, 702]

function so_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 0
    # reservation pseudo-code crosswalk (so/forkod.f): 7710→602, others → R6 forests
    if kodfor == 7711 || kodfor == 7710
        ifor = 2                       # Fort McDermitt → Fremont(602)=IFOR 2
    else
        for (i, f) in enumerate(SO_JFOR)
            kodfor == f && (ifor = i; break)
        end
        ifor == 0 && (ifor = 1)        # default Deschutes(601)=IFOR 1
    end
    p.forest_idx = Int32(ifor)
    p.user_forest_code = Int32(SO_JFOR[ifor])
    return ifor
end

# so/htcalc.f — predicted total height (ft) at age `ag` for species `ispc`, site index `sindx`, forest `jfor`.
const SO_DUNL1 = Float32[-88.9, -82.2, -78.3, -82.1, -56.0, -33.8]
const SO_DUNL2 = Float32[49.7067, 44.1147, 39.1441, 35.4160, 26.7173, 18.6400]
const SO_DUNL3 = Float32[2.375, 2.025, 1.650, 1.225, 1.075, 0.875]

@inline function _so_dunlevitan(sindx::Float32, ag::Float32)::Float32
    indx = sindx <= 44f0 ? 6 : sindx <= 52f0 ? 5 : sindx <= 65f0 ? 4 :
           sindx <= 82f0 ? 3 : sindx <= 98f0 ? 2 : 1
    return ag <= 40f0 ? SO_DUNL3[indx] * ag : SO_DUNL1[indx] + SO_DUNL2[indx] * log(ag)
end

function so_htcalc(jfor::Int, sindx::Float32, ispc::Int, ag::Float32)::Float32
    if ispc == 1                                                  # WP Brickell
        return sindx / (0.37504453f0 * (1f0 - 0.92503f0 * exp(-0.0207959f0 * ag))^(-2.4881068f0))
    elseif ispc == 3 || ispc == 32                                # DF/OS Cochran-251
        return 4.5f0 + exp(-0.37496f0 + 1.36164f0 * log(ag) - 0.00243434f0 * log(ag)^4) -
               79.97f0 * (-0.2828f0 + 1.87947f0 * (1f0 - exp(-0.022399f0 * ag))^0.966998f0) +
               (sindx - 4.5f0) * (-0.2828f0 + 1.87947f0 * (1f0 - exp(-0.022399f0 * ag)^0.966998f0))
    elseif ispc == 4 || ispc == 6 || ispc == 12 || ispc == 14     # WF/IC/GF/SF Cochran-252
        la = log(ag)
        x2 = -0.30935f0 + 1.2383f0 * la + 0.001762f0 * la^4 - 5.4f-6 * la^9 + 2.046f-7 * la^11 - 4.04f-13 * la^18
        x3 = -6.2056f0 + 2.097f0 * la - 0.09411f0 * la^2 - 0.00004382f0 * la^7 + 2.007f-11 * la^16 - 2.054f-17 * la^24
        return exp(x2) - 84.73f0 * exp(x3) + (sindx - 4.5f0) * exp(x3) + 4.5f0
    elseif ispc == 5                                              # MH interim-means (metric)
        h = (22.8741f0 + 0.950234f0 * sindx) * (1f0 - exp(-0.00206465f0 * sqrt(sindx) * ag))^(1.365566f0 + 2.045963f0 / sindx)
        return (h + 1.37f0) * 3.281f0
    elseif ispc == 7                                              # LP Dahms
        return sindx * (-0.0968f0 + 0.02679f0 * ag - 0.00009309f0 * ag * ag)
    elseif ispc == 8                                              # ES Alexander
        return 4.5f0 + ((2.75780f0 * sindx^0.83312f0) * (1f0 - exp(-0.015701f0 * ag))^(22.71944f0 * sindx^(-0.63557f0)))
    elseif ispc == 9                                              # SH: R6 Dolph red fir else Dunning-Levitan
        if jfor <= 3 || jfor == 10
            term = ag * exp(ag * (-0.0440853f0)) * 1.41512f-6
            b = sindx * term - 3.04951f6 * term * term + 5.72474f-4
            term2 = 50f0 * exp(50f0 * (-0.0440853f0)) * 1.41512f-6
            b50 = sindx * term2 - 3.04951f6 * term2 * term2 + 5.72474f-4
            return (sindx - 4.5f0) * (1f0 - exp(-b * ag^1.51744f0)) / (1f0 - exp(-b50 * 50f0^1.51744f0)) + 4.5f0
        else
            return _so_dunlevitan(sindx, ag)
        end
    elseif ispc == 2 || ispc == 10                                # SP/PP Barrett
        t = -0.7864f0 + 2.49717f0 * (1f0 - exp(-0.0045042f0 * ag))^0.33022f0
        return (128.89522f0 * (1f0 - exp(-0.016959f0 * ag))^1.23114f0) - (t * 100.43f0) + (t * (sindx - 4.5f0)) + 4.5f0
    elseif ispc == 13                                             # AF Johnson/DeMars
        return sindx * (-0.07831f0 + 0.0149f0 * ag - 4.0818f-5 * ag * ag)
    elseif ispc == 15                                             # NF Herman
        x1 = -564.38f0 + 22.250f0 * (sindx - 4.5f0) - 0.04995f0 * (sindx - 4.5f0)^2
        x2 = 6.8f0 + 2843.21f0 / (sindx - 4.5f0) + 34735.54f0 / (sindx - 4.5f0)^2
        return 4.5f0 + (sindx - 4.5f0) / (x1 * (1f0 / ag)^2 + x2 * (1f0 / ag) + 1f0 - 0.0001f0 * x1 - 0.01f0 * x2)
    elseif ispc == 17                                             # WL Cochran-424
        t = -0.12528f0 + 0.039636f0 * ag - 0.0004278f0 * ag * ag + 1.7039f-6 * ag^3
        return 4.5f0 + 1.46897f0 * ag + 0.0092466f0 * ag * ag - 0.00023957f0 * ag^3 + 1.1122f-6 * ag^4 +
               (sindx - 4.5f0) * t - 73.57f0 * t
    elseif ispc == 18                                             # RC Hegyi
        return 1.3283f0 * sindx * ((1f0 - exp(-0.0174f0 * ag))^1.4711f0)
    elseif ispc == 19                                             # WH Wiley
        z = 2500f0 / (sindx - 4.5f0)
        return ag * ag / (-1.73070f0 + 0.1394f0 * z + (-0.0616f0 + 0.0137f0 * z) * ag +
               (0.00192f0 + 0.00007f0 * z) * (ag * ag)) + 4.5f0
    elseif ispc == 20 || ispc == 21 || ispc == 23 || ispc == 25 || ispc == 26 ||
           (28 <= ispc <= 31) || ispc == 33                       # misc Curtis
        return (sindx - 4.5f0) / (0.6192f0 - 5.3394f0 / (sindx - 4.5f0) + 240.29f0 * ag^(-1.4f0) +
               (3368.9f0 / (sindx - 4.5f0)) * ag^(-1.4f0)) + 4.5f0
    elseif ispc == 22                                             # RA Harrington
        c = 59.5864f0 + 0.79530f0 * sindx
        return sindx + c * (1f0 - exp((0.00194f0 - 0.000740f0 * sindx) * ag))^0.9198f0 -
               c * (1f0 - exp((0.00194f0 - 0.000740f0 * sindx) * 20f0))^0.9198f0
    elseif ispc == 27                                             # WO: R6 Powers else Dunning-Levitan
        if jfor <= 3 || jfor == 10
            term = sqrt(ag) - sqrt(50f0)
            return (sindx * (1f0 + 0.322f0 * term)) - 6.413f0 * term
        else
            return _so_dunlevitan(sindx, ag)
        end
    else                                                          # WJ(11)/WB(16)/AS(24) + future = no curve
        return 0f0
    end
end

# so/siterange.f — per-species site-index range (SITELO/SITEHI); used to interpolate SITEAR for the
# no-curve species WJ(11)/WB(16)/AS(24).
const SO_SITELO = Float32[13,27,21,5,5,5,5,12,10,7, 5,9,6,4,7,20,60,29,6,5, 5,56,108,30,10,10,21,20,5,5, 5,5,5]
const SO_SITEHI = Float32[137,178,148,195,133,169,140,227,134,176, 40,173,127,221,210,65,147,152,203,75,
                          100,192,142,66,191,104,85,93,100,75, 75,175,125]

# so/sitset.f SDIDEF fan coefficients (chunk 4c): C6 = R6 per-species SDImax-ratio basis (IFOR≤3/10);
# C5 = the non-R6 fallback SDImax. FORMAX = SDImax cap. PMSDIU (so/grinit.f:240) = the BAMAX→SDImax %.
const SO_SDIDEF_C6 = Float32[447,447,767,659,758,447,541,659,750,429, 500,659,659,659,659,250,250,650,650,650,
                             200,300,300,250,250,200,200,150,250,100, 100,447,250]
const SO_SDIDEF_C5 = Float32[272,561,570,800,687,576,679,620,1000,365, 272,800,602,790,1000,621,423,762,682,576,
                             441,441,629,562,452,441,440,447,785,501, 501,365,441]
const SO_FORMAX = 850f0
const SO_PMSDIU = 85f0

# so/habtyp.f + so/ecocls.f + so/pvref6.f — plant-association (PA) → per-species SDImax (SDIDEF). SO's R6
# forests (IFOR≤3 or 10) take the ECOCLS PA-specific SDImax; the FIA DB delivers PV_CODE as the ALPHA PA
# code (e.g. "CWS313"). Without decoding it, so_sitset! defaulted every stand's site-species SDImax to the
# CPS111 default RSDI (285) instead of the stand's own ecoclass (CWS313 → 810), so on a non-default PA the
# BA-weighted SDIMAX ran wrong (~2.8× low on CWS313) and the density self-thin fired at the wrong level.
# PCOML (92 R6 PA codes; KODTYP indexes it; default PA = CPS111 = index 49), ECOCLS (92 PA rows → site
# species FVSSEQ / RSDI / RSI), PVREF6 ((PV_CODE,PV_REF_CODE)→HABPVR crosswalk).
function _so_load_site_tables()
    pcoml = String[]
    for l in readlines(joinpath(SO_DATADIR, "pcoml.csv"))[2:end]
        isempty(strip(l)) && continue
        push!(pcoml, String(strip(split(l, ',')[2])))
    end
    rows = NamedTuple{(:pa,:spc,:fvsseq,:sdimx,:site,:numbr,:iflag),
                      Tuple{String,String,Int,Float32,Float32,Int,Int}}[]
    for l in readlines(joinpath(SO_DATADIR, "ecocls.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(strip(l), ',')
        push!(rows, (pa=String(f[1]), spc=String(f[2]), fvsseq=parse(Int, f[3]),
                     sdimx=parse(Float32, f[4]), site=parse(Float32, f[5]),
                     numbr=parse(Int, f[6]), iflag=parse(Int, f[7])))
    end
    return pcoml, rows
end
const SO_PCOML, SO_ECOCLS = _so_load_site_tables()

# so/pvref6.f — (PV_CODE, PV_REF_CODE) → HABPVR crosswalk. FVS EXITs on the FIRST full match; a partial/no
# match leaves KARD2 blank ⇒ default. Only non-blank HABPVR rows stored (missing key ⇒ blank ⇒ default).
const SO_PVREF6 = let d = Dict{Tuple{String,String},String}()
    for l in readlines(joinpath(SO_DATADIR, "pvref6.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(l, ','; limit = 3)
        c = String(strip(f[1])); r = String(strip(f[2]))
        h = length(f) >= 3 ? String(strip(f[3])) : ""
        (isempty(c) || isempty(h)) && continue
        key = (c, r); haskey(d, key) || (d[key] = h)
    end
    d
end

const SO_HAB_DEFAULT_PA = "CPS111"                        # so/habtyp.f R6 default (PCOML[49])
# so/habtyp.f (R6) — KODTYP → PCOM plant-association code; anything out of 1..92 ⇒ the CPS111 default.
so_habtyp(kodtyp::Integer)::String =
    (1 <= kodtyp <= length(SO_PCOML)) ? SO_PCOML[kodtyp] : SO_HAB_DEFAULT_PA
so_ecocls(pa::AbstractString) = filter(r -> r.pa == pa, SO_ECOCLS)

# so/habtyp.f (R6 path) — decode the FIA alpha PV_CODE (+ optional PV_REF_CODE) into the KODTYP index into
# SO_PCOML. When a reference code is present PVREF6 crosswalks (pv,ref)→HABPVR (blank on partial/no match ⇒
# 0 ⇒ CPS111 default); then HBDECD string-matches the (crosswalked) code against PCOML. No-ref ⇒ the raw
# PV_CODE is matched directly.
function so_habitat_kodtyp(pv::AbstractString, pvref::AbstractString)
    pvs = String(strip(pv)); refs = String(strip(pvref))
    kard2 = pvs
    if !isempty(refs)
        kard2 = get(SO_PVREF6, (pvs, refs), "")           # blank on partial/no match ⇒ default
    end
    isempty(kard2) && return 0
    idx = findfirst(==(kard2), SO_PCOML)
    idx === nothing ? 0 : Int(idx)
end

# so/sichg.f — SIAGE(i) per species (reference age for the site-species curve). RF(5) is metric.
function so_sichg(s::StandState, isisp::Integer, ssite::Float32)
    sd = s.coef.species
    a = sd[:sichg_a]; b = sd[:sichg_b]; refage = sd[:sichg_refage]; refloc = sd[:sichg_refloc]
    maxsp = nspecies(s.variant)
    isiloc = refloc[isisp]
    siage = Vector{Float32}(undef, maxsp)
    @inbounds for i in 1:maxsp
        diff = 0
        (isiloc == 1f0 && refloc[i] == 0f0) && (diff = -1)        # site 'T', species 'B'
        (isiloc == 0f0 && refloc[i] == 1f0) && (diff = 1)         # site 'B', species 'T'
        age2bh = 0f0
        if diff != 0
            age2bh = a[i] + b[i] * ssite
            (isisp != 5 && i == 5) && (age2bh = a[i] + b[i] * (ssite / 3.281f0))
            (isisp == 5 && i != 5) && (age2bh = a[i] + b[i] * (ssite * 3.281f0))
        end
        i == 9 && (age2bh = 18f0)                                 # SH fixed
        i == 27 && (age2bh = 2f0)                                 # WO fixed
        siage[i] = refage[i] + age2bh * Float32(diff)
    end
    return siage
end

function so_sitset!(s::StandState)
    p = s.plot; maxsp = nspecies(s.variant)
    ifor = Int(p.forest_idx)
    # so/sitset.f:81 — R6 forests (IFOR≤3 or =10) use Reineke (LZEIDE=.FALSE.); others keep Zeide.
    (ifor <= 3 || ifor == 10) && (s.control.zeide_sdi = false)

    isisp = (1 <= Int(p.site_species) <= maxsp) ? Int(p.site_species) : 0
    isisp <= 0 && (isisp = 10)                                    # R6 default site species = PP(10)
    p.sp_site_index[isisp] <= 0f0 && (p.sp_site_index[isisp] = 70f0)   # so/sitset.f default SITEAR

    sindx = p.sp_site_index[isisp]
    siage = so_sichg(s, isisp, sindx)
    slossp = SO_SITELO[isisp]; shissp = SO_SITEHI[isisp]
    @inbounds for ispc in 1:maxsp
        ispc == isisp && continue
        if ispc == 11 || ispc == 16 || ispc == 24       # WJ/WB/AS: no site curve → SITERANGE interpolation
            tem = sindx < slossp ? slossp : sindx
            v = SO_SITELO[ispc] + (tem - slossp) / (shissp - slossp) * (SO_SITEHI[ispc] - SO_SITELO[ispc])
        else
            v = so_htcalc(ifor, sindx, isisp, siage[ispc])
            ispc == 5 && (v = v / 3.281f0; v > 28f0 && (v = 28f0))   # MH metric conversion, cap 28
        end
        v < 0f0 && (v = 0f0)
        p.sp_site_index[ispc] = v                        # SO sitset is authoritative (overrides the generic reader)
    end

    # SDIDEF (per-species SDImax) — so/sitset.f:79-125 (R6 ECOCLS PA seed) + :231-247 fan (chunk 4c; prereq
    # for crown+mort RELSDI). For R6 forests (IFOR≤3 or 10) so/sitset.f ECOCLS the stand's plant association
    # (habitat_code → PCOML → PA; the CPS111 default when unresolved) and seeds the site species' SDImax from
    # that PA's ecoclass RSDI (so/sitset.f:112-119: SDIDEF(ISEQ)=RSDI, and SDIDEF(ISISP)=RSDI on the IFLAG=1
    # row). Without this seed every stand rode the CPS111 default RSDI 285 regardless of PA ⇒ on a non-default
    # PA (e.g. CWS313 → 810) the SDIMAX was ~2.8× wrong (measured vs FVSso_clean). The FIA reader now decodes
    # PV_CODE → habitat_code so this picks the stand's real PA (a no-habitat / unresolved stand still defaults
    # to CPS111 285, so sot01 + the default case stay bit-exact vs the FVSso_g16 SDIDEF dump).
    # Fan (IFOR≤3/10 R6, BAMAX unset): SDIDEF[i]=SDIDEF[ISISP]·C6[i]/C6[ISISP] (cap FORMAX); else C5[i].
    # BAMAX-keyword branch (SDIDEF=BAMAX/(0.5454154·PMSDIU/100)) is a follow-on — no BAMAX (=0) here, so the
    # R6 C6-ratio fan below is the exercised path; keep the branch for when a BAMAX keyword lands.
    bamax = 0f0
    if ifor <= 3 || ifor == 10
        # so/sitset.f R6 ECOCLS: seed each of the PA's ecoclass species' SDImax; ISISP (from DB) keeps its
        # PA RSDI on the IFLAG=1 row. Unresolved habitat ⇒ so_habtyp default CPS111 (RSDI 285, site sp PP).
        pa = so_habtyp(Int(p.habitat_code))
        rows = so_ecocls(pa)
        isempty(rows) && (rows = so_ecocls(SO_HAB_DEFAULT_PA))
        @inbounds for r in rows
            iseq = r.fvsseq; (iseq < 1 || iseq > maxsp) && continue
            rsdi = min(r.sdimx, SO_FORMAX)
            (isisp <= 0 && r.iflag == 1) && (isisp = iseq)
            p.sp_sdi_def[iseq] <= 0f0 && (p.sp_sdi_def[iseq] = rsdi)
            (isisp > 0 && r.iflag == 1 && p.sp_sdi_def[isisp] <= 0f0) && (p.sp_sdi_def[isisp] = rsdi)
        end
    end
    p.sp_sdi_def[isisp] <= 0f0 && (p.sp_sdi_def[isisp] = 285f0)   # global fallback (ECOCLS CPS111 site sp)
    k = isisp
    @inbounds for i in 1:maxsp
        p.sp_sdi_def[i] > 0f0 && continue
        if bamax > 0f0
            p.sp_sdi_def[i] = bamax / (0.5454154f0 * (SO_PMSDIU / 100f0))
        elseif ifor <= 3 || ifor == 10
            v = p.sp_sdi_def[k] * (SO_SDIDEF_C6[i] / SO_SDIDEF_C6[k])
            v > SO_FORMAX && (v = SO_FORMAX)
            p.sp_sdi_def[i] = v
        else
            p.sp_sdi_def[i] = SO_SDIDEF_C5[i]
        end
    end

    p.site_species = Int32(isisp)
    return s
end

function so_site_index_setup!(s::StandState)
    so_forkod!(s.plot)
    so_sitset!(s)
    return s
end

site_setup!(s::StandState, ::SouthCentralOregon) = so_site_index_setup!(s)
