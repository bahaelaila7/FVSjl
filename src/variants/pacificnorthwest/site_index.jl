# =============================================================================
# site_index.jl (pacificnorthwest) — PN site index + Reineke SDImax + forest code. Chunk 2.
# (pn/forkod.f, pn/habtyp.f, pn/ecocls.f, pn/sichg.f, pn/htcalc.f, pn/sitset.f)
#
# Same Region-6 chain as WC (habitat → PCOML → PA → ECOCLS → site species/SI/SDImax; SICHG ref age;
# pn_htcalc fans SITEAR to every species; misc-hardwood reductions). Differences from WC (agent-verified):
#   • PCOML/ECOCLS are PN tables (75 rows; default PA CHS133 vs WC CFS551).
#   • PN_JFOR = [609,612,800,708,709,712] (612→IFOR 2); reservation pseudo-codes 81xx; NO IFOR remap.
#   • FMSDI is a SCALAR 950 (all forests); Region-6 default site species = DF(16) (WC ES/10); default
#     SITEAR = 100 (WC 70). MH(20)/WO(28) SI transforms + misc-hardwood site_redux IDENTICAL to WC.
# pnt01 forest 612 → IFOR 2, habitat 40 → CHS133 → site species DF(16), SI 98, SDIDEF (capped at 950).
# =============================================================================

function _pn_load_site_tables()
    pcoml = String[]
    for l in readlines(joinpath(PN_DATADIR, "pcoml.csv"))[2:end]
        isempty(strip(l)) && continue
        push!(pcoml, String(strip(split(l, ',')[2])))
    end
    rows = NamedTuple{(:pa,:spc,:fvsseq,:sdimx,:site,:numbr,:iflag),
                      Tuple{String,String,Int,Float32,Float32,Int,Int}}[]
    for l in readlines(joinpath(PN_DATADIR, "ecocls.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(strip(l), ',')
        push!(rows, (pa=String(f[1]), spc=String(f[2]), fvsseq=parse(Int, f[3]),
                     sdimx=parse(Float32, f[4]), site=parse(Float32, f[5]),
                     numbr=parse(Int, f[6]), iflag=parse(Int, f[7])))
    end
    return pcoml, rows
end
const PN_PCOML, PN_ECOCLS = _pn_load_site_tables()

# pn/forkod.f — KODFOR → IFOR 1..6. NO post-lookup remap. Reservation pseudo-codes map straight to IFOR.
const PN_JFOR = Int[609, 612, 800, 708, 709, 712]
const _PN_RES_IFOR1 = Set(Int[8110,8111,8113,8114,8115,8116,8119,8120,8121,8122,8123,8125,8126,8127,8128,8129])

function pn_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 1; useigl = true
    if kodfor == 8101 || kodfor == 8102 || kodfor == 8103
        ifor = 2                                  # → 612 Siuslaw
    elseif kodfor == 8104 || kodfor == 8105
        ifor = 6                                  # → 712 BLM Coos Bay
    elseif kodfor in _PN_RES_IFOR1
        ifor = 1                                  # → 609 Olympic
    else
        idx = findfirst(==(kodfor), PN_JFOR)
        idx === nothing ? (useigl = false; ifor = 1) : (ifor = idx)
    end
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(1))
    p.user_forest_code = Int32(PN_JFOR[ifor])
    return ifor
end

pn_habtyp(kodtyp::Integer)::String =
    (kodtyp <= 0 || kodtyp > length(PN_PCOML)) ? "" : PN_PCOML[kodtyp]
pn_ecocls(pa::AbstractString) = filter(r -> r.pa == pa, PN_ECOCLS)

# pn/sichg.f — SIAGE(I) per species (from the CSV sichg_a/b/refage/refloc; PN values differ at 6/16/18).
function pn_sichg(s::StandState, isisp::Integer, ssite::Float32)
    sd = s.coef.species
    a = sd[:sichg_a]; b = sd[:sichg_b]; refage = sd[:sichg_refage]; refloc = sd[:sichg_refloc]
    maxsp = nspecies(s.variant)
    isiloc = refloc[isisp]
    siage = Vector{Float32}(undef, maxsp)
    @inbounds for i in 1:maxsp
        diff = 0
        (isiloc == 1f0 && refloc[i] == 0f0) && (diff = -1)
        (isiloc == 0f0 && refloc[i] == 1f0) && (diff = 1)
        age2bh = 0f0
        if diff != 0
            age2bh = a[i] + b[i] * ssite
            (isisp != 20 && i == 20) && (age2bh = a[i] + b[i] * (ssite / 3.281f0))
            (isisp == 20 && i != 20) && (age2bh = a[i] + b[i] * (ssite * 3.281f0))
        end
        siage[i] = refage[i] + age2bh * Float32(diff)
    end
    return siage
end

const PN_FMSDI = 950.0f0                           # pn/sitset.f DATA FORMAX/950./ (scalar, all forests)

function pn_sitset!(s::StandState)
    p = s.plot; sd = s.coef.species
    maxsp = nspecies(s.variant)
    formax = PN_FMSDI
    nsiset = count(>(0f0), @view p.sp_site_index[1:maxsp])
    pcom = pn_habtyp(Int(p.habitat_code))
    isempty(pcom) && (pcom = "CHS133")            # pn/habtyp.f default PA (ITYPE=40 → PCOML[40])
    rows = pn_ecocls(pcom)

    isisp = (1 <= Int(p.site_species) <= maxsp) ? Int(p.site_species) : 0
    @inbounds for r in rows
        iseq = r.fvsseq; iseq == 0 && continue
        rsdi = min(r.sdimx, formax)
        (isisp <= 0 && r.iflag == 1) && (isisp = iseq)
        (p.sp_site_index[iseq] <= 0f0 && nsiset == 0) && (p.sp_site_index[iseq] = r.site)
        p.sp_sdi_def[iseq] <= 0f0 && (p.sp_sdi_def[iseq] = rsdi)
    end
    isisp <= 0 && (isisp = 16)                     # pn/sitset.f Region-6 default site sp = DF(16)
    p.sp_site_index[isisp] <= 0f0 && (p.sp_site_index[isisp] = 100f0)   # pn/sitset.f default SITEAR
    p.site_species = Int32(isisp)

    sindx = p.sp_site_index[isisp]
    siage = pn_sichg(s, isisp, sindx)
    redux = sd[:site_redux]
    si = Vector{Float32}(undef, maxsp)
    @inbounds for ispc in 1:maxsp
        if p.sp_site_index[ispc] > 0f0
            si[ispc] = p.sp_site_index[ispc]; continue
        end
        v = pn_htcalc(sindx, isisp, siage[ispc])
        if ispc == 20 && isisp != 20
            v = v / 3.281f0; v > 28f0 && (v = 28f0)
        elseif ispc == 28 && isisp != 28
            v = 114.2f0 * (1f0 - exp(-0.0266f0 * v))^2.26f0
        elseif ispc != isisp && redux[ispc] != 1f0
            v = v * redux[ispc]
        end
        si[ispc] = v
    end
    @inbounds for i in 1:maxsp
        p.sp_site_index[i] <= 0f0 && (p.sp_site_index[i] = si[i])
    end
    @inbounds for i in 1:maxsp
        p.sp_sdi_def[i] <= 0f0 && (p.sp_sdi_def[i] = p.sp_sdi_def[isisp])
    end
    return s
end

function pn_site_index_setup!(s::StandState)
    pn_forkod!(s.plot)
    pn_sitset!(s)
    return s
end

site_setup!(s::StandState, ::PacificNorthwest) = pn_site_index_setup!(s)
