# =============================================================================
# site_index.jl (eastcascades) — EC site index + Reineke SDImax + forest code. Chunk 2.
# (ec/forkod.f, ec/habtyp.f, ec/ecocls.f, ec/sichg.f, ec/htcalc.f, ec/sitset.f)
#
# Region-6 chain (habitat → PCOML → PA → ECOCLS → site species/SI/SDImax; SICHG ref age; ec_htcalc fans
# SITEAR to every species; site reductions). EC specifics:
#   • forkod: JFOR=[606,608,617,699,603,613,621]; IFOR corrections (603/613→617, 621→699); reservation
#     crosswalk (8106→2, 8117→1, 8130→6, 8131→3). Forest 608 → IFOR 2 (Okanogan), IGL 1.
#   • sitset: FORMAX=900 (scalar), Region-6 default site species = PP(10), default SITEAR=70; SDICON
#     per-species SDI-ratio expansion; MH/OS(12,31) metric SI, WO(28) King's-DF transform, per-species
#     site_redux for the WC-form hardwoods.
# ect01 forest 608 habitat 12 → PA CDG132 → site species DF(3), SI 69, SDImax 550 (matches live SITEAR).
# =============================================================================

const EC_FORMAX = 900.0f0                          # ec/sitset.f DATA FORMAX/900./
# ec/sitset.f SDICON — per-species Reineke SDImax reference (for the non-site-species SDIDEF expansion).
const EC_SDICON = Float32[645,648,766,766,766,766,674,766,700,645,900,766,900,900,900,766,900,900,900,900,900,900,900,900,900,900,900,900,900,900,766,900]
const EC_PMSDIU = 85.0f0                            # ec/grinit.f PMSDIU

function _ec_load_site_tables()
    pcoml = String[]
    for l in readlines(joinpath(EC_DATADIR, "pcoml.csv"))[2:end]
        isempty(strip(l)) && continue
        push!(pcoml, String(strip(split(l, ',')[2])))
    end
    rows = NamedTuple{(:pa,:spc,:fvsseq,:sdimx,:site,:numbr,:iflag),
                      Tuple{String,String,Int,Float32,Float32,Int,Int}}[]
    for l in readlines(joinpath(EC_DATADIR, "ecocls.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(strip(l), ',')
        push!(rows, (pa=String(f[1]), spc=String(f[2]), fvsseq=parse(Int, f[3]),
                     sdimx=parse(Float32, f[4]), site=parse(Float32, f[5]),
                     numbr=parse(Int, f[6]), iflag=parse(Int, f[7])))
    end
    return pcoml, rows
end
const EC_PCOML, EC_ECOCLS = _ec_load_site_tables()

# ec/forkod.f — KODFOR → IFOR 1..7 with corrections + reservation crosswalk.
const EC_JFOR = Int[606, 608, 617, 699, 603, 613, 621]
const EC_KFOR = Int[1, 1, 1, 1, 2, 3, 1]

function ec_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 1; igl = 1; useigl = true; forfound = false
    if kodfor == 8106
        ifor = 2                                   # Colville Res → Okanogan
    elseif kodfor == 8117
        ifor = 1                                   # Umatilla Res → Mt Hood
    elseif kodfor == 8130
        ifor = 6                                   # Yakama Res → Baker-Snoqualmie
    elseif kodfor == 8131
        ifor = 3                                   # Spokane Res → Wenatchee
    else
        idx = findfirst(==(kodfor), EC_JFOR)
        idx === nothing ? (useigl = false; ifor = 1) : (ifor = idx; forfound = true)
    end
    # forest mapping corrections
    if ifor == 5
        ifor = 3; igl = 2; useigl = false          # Gifford Pinchot → Wenatchee
    elseif ifor == 6
        ifor = 3; igl = 3; useigl = false          # Baker-Snoqualmie → Wenatchee
    elseif ifor == 7
        ifor = 4; igl = 1; useigl = false          # Colville → Okanogan (Tonasket)
    end
    useigl && (igl = EC_KFOR[ifor])
    p.forest_idx = Int32(ifor)
    p.geo_location = Int32(igl)
    p.user_forest_code = Int32(EC_JFOR[ifor])
    return ifor
end

ec_habtyp(kodtyp::Integer)::String =
    (kodtyp <= 0 || kodtyp > length(EC_PCOML)) ? "" : EC_PCOML[kodtyp]
ec_ecocls(pa::AbstractString) = filter(r -> r.pa == pa, EC_ECOCLS)

# ec/sichg.f — SIAGE(I) per species; metric species = MH(12) & OS(31).
function ec_sichg(s::StandState, isisp::Integer, ssite::Float32)
    sd = s.coef.species
    a = sd[:sichg_a]; b = sd[:sichg_b]; refage = sd[:sichg_refage]; refloc = sd[:sichg_refloc]
    maxsp = nspecies(s.variant)
    isiloc = refloc[isisp]
    siage = Vector{Float32}(undef, maxsp)
    @inbounds for i in 1:maxsp
        diff = 0
        (isiloc == 1f0 && refloc[i] == 0f0) && (diff = -1)      # site 'T', species 'B'
        (isiloc == 0f0 && refloc[i] == 1f0) && (diff = 1)       # site 'B', species 'T'
        age2bh = 0f0
        if diff != 0
            age2bh = a[i] + b[i] * ssite
            (isisp != 31 && i == 31) && (age2bh = a[i] + b[i] * (ssite / 3.281f0))
            (isisp == 31 && i != 31) && (age2bh = a[i] + b[i] * (ssite * 3.281f0))
            (isisp != 12 && i == 12) && (age2bh = a[i] + b[i] * (ssite / 3.281f0))
            (isisp == 12 && i != 12) && (age2bh = a[i] + b[i] * (ssite * 3.281f0))
        end
        siage[i] = refage[i] + age2bh * Float32(diff)
    end
    return siage
end

function ec_sitset!(s::StandState)
    p = s.plot; sd = s.coef.species
    maxsp = nspecies(s.variant)
    formax = EC_FORMAX
    nsiset = count(>(0f0), @view p.sp_site_index[1:maxsp])
    pcom = ec_habtyp(Int(p.habitat_code))
    isempty(pcom) && (pcom = "CPS241")             # ec/habtyp.f default (ITYPE=114 → PCOML[114])
    rows = ec_ecocls(pcom)

    isisp = (1 <= Int(p.site_species) <= maxsp) ? Int(p.site_species) : 0
    @inbounds for r in rows
        iseq = r.fvsseq; iseq == 0 && continue
        rsdi = min(r.sdimx, formax)
        (isisp <= 0 && r.iflag == 1) && (isisp = iseq)
        (p.sp_site_index[iseq] <= 0f0 && nsiset == 0) && (p.sp_site_index[iseq] = r.site)
        (isisp > 0 && r.iflag == 1 && p.sp_sdi_def[isisp] <= 0f0) && (p.sp_sdi_def[isisp] = rsdi)
    end
    isisp <= 0 && (isisp = 10)                      # ec/sitset.f Region-6 default site sp = PP(10)
    p.sp_site_index[isisp] <= 0f0 && (p.sp_site_index[isisp] = 70f0)   # ec/sitset.f default SITEAR

    sindx = p.sp_site_index[isisp]
    siage = ec_sichg(s, isisp, sindx)
    redux = sd[:site_redux]
    si = Vector{Float32}(undef, maxsp)
    @inbounds for ispc in 1:maxsp
        if p.sp_site_index[ispc] > 0f0
            si[ispc] = p.sp_site_index[ispc]; continue
        end
        v = ec_htcalc(sindx, isisp, siage[ispc])   # ec/sitset.f: HTCALC(SINDX, ISISP, SIAGE(ISPC)) — site-species curve
        if (ispc == 12 && isisp != 12) || (ispc == 31 && isisp != 31)
            v = v / 3.281f0; v > 28f0 && (v = 28f0)                   # MH/OS metric, cap 28
        elseif ispc == 28 && isisp != 28
            v = 114.2f0 * (1f0 - exp(-0.0266f0 * v))^2.26f0          # WO King's DF
        elseif ispc != isisp && redux[ispc] != 1f0
            v = v * redux[ispc]                                       # WC-form site reduction
        end
        si[ispc] = v
    end
    @inbounds for i in 1:maxsp
        p.sp_site_index[i] <= 0f0 && (p.sp_site_index[i] = si[i])
    end

    # SDIDEF expansion (ec/sitset.f DO 80): SDICON ratio for EC-form, site-species value for WC-form.
    k = isisp
    @inbounds for i in 1:maxsp
        p.sp_sdi_def[i] > 0f0 && continue
        if i in EC_NATIVE_SP                          # BAMAX starts 0 ⇒ SDICON-ratio branch
            v = p.sp_sdi_def[k] * (EC_SDICON[i] / EC_SDICON[k])
            v > formax && (v = formax)
            p.sp_sdi_def[i] = v
        else
            p.sp_sdi_def[i] = p.sp_sdi_def[k]
        end
    end
    p.site_species = Int32(isisp)
    return s
end

function ec_site_index_setup!(s::StandState)
    ec_forkod!(s.plot)
    ec_sitset!(s)
    return s
end

site_setup!(s::StandState, ::EastCascades) = ec_site_index_setup!(s)
