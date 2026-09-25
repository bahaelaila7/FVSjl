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

# nc/forkod.f DATA JFOR(11): IFOR 1..7 are the modelled forests; 8..11 (Trinity 518, Los Padres 507, Mendocino
# 508, Simpson 715) exist only to be remapped by the FOREST MAPPING CORRECTION below.
const NC_JFOR = Int[505, 510, 514, 611, 705, 800, 712, 518, 507, 508, 715]
# nc/forkod.f reservation pseudo-codes → pre-correction IFOR (13 tribal lands → 10 ⇒ Six Rivers after mapping).
const NC_FOR_RESERV = Dict{Int,Int}(
    7806=>10, 7807=>10, 7810=>10, 7813=>10, 7815=>10, 7816=>10, 7820=>10, 7821=>10, 7824=>10,
    7830=>10, 7831=>10, 7833=>10, 7834=>10, 7839=>2, 7841=>2, 7843=>2, 7845=>1, 8103=>4, 8105=>4)

"""nc/forkod.f: KODFOR → IFOR 1..7 and KODFOR = JFOR(IFOR). Reservation codes first, else a JFOR match (not found
⇒ ERRGRO 3, IFOR stays 1 = Klamath, IGL untouched), then the mapping correction 8→3 (518→Shasta-Trinity),
9,10→2 (507/508→Six Rivers), 11→6 (715→800). Used to match only the first 7 codes, so Trinity/Mendocino/Los
Padres/Simpson and every reservation stand silently ran as Klamath (IFOR 1)."""
function nc_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = get(NC_FOR_RESERV, kodfor, 0); useigl = true
    if ifor == 0
        idx = findfirst(==(kodfor), NC_JFOR)
        idx === nothing ? (ifor = 1; useigl = false) : (ifor = idx)
    end
    ifor == 8 ? (ifor = 3) : (ifor == 9 || ifor == 10) ? (ifor = 2) : ifor == 11 && (ifor = 6)
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(1))     # IGL = KFOR(IFOR) = 1
    p.user_forest_code = Int32(NC_JFOR[ifor])
    return ifor
end

# nc/sitset.f site-index species defaults (SI array): the per-species site index when NO SITECODE, with the
# nc/sitset.f DO 30/35 site-index conversion (SICHG + HTCALC): the site species' index (ISISP, default DF sp3
# at 90) is converted to a per-species site index via a height-at-age curve. SICHG gives each species' site
# age SIAGE; HTCALC evaluates the SITE SPECIES' potential-height curve at that age. This REPLACES the old
# hardcoded NC_SITE_DEFAULT=[90,100,90,…] (measured/rounded from the DF=90 default only) — that was wrong for
# any stand whose DB site index ≠ 90 (the CA redwood/hardwood stands), giving redwood SITEAR 90 vs the oracle's
# 109.976 at DF site 110 ⇒ a redwood htgr5 (small-tree height) + sp12 DGCON (ln SITEAR) diameter error.

# SICHG (nc/sichg.f): per-species site AGE from the reference species' site index. REFLOC = total-age ('T':
# SP sp2, PP sp10) vs breast-height-age ('B': all others) basis; the age slides to align the reference.
const NC_SICHG_B      = Float32[-0.08,-0.05,-0.08,-0.07,-0.02,-0.05,-0.05,-0.03,-0.06,-0.05,-0.03,-0.08]
const NC_SICHG_A      = Float32[10,12,10,10,3,10,6,4,10,12,4,10]
const NC_SICHG_SIMIN  = Float32[50,40,50,30,50,30,30,50,30,40,50,50]
const NC_SICHG_SIMAX  = Float32[150,120,150,130,100,130,70,90,130,120,90,150]
const NC_SICHG_IREFAG = 50f0

# REFLOC: 'T' (total-age) basis for SP (sp2) and PP (sp10); 'B' (breast-height age) for all others (sichg.f).
@inline _nc_refloc_t(sp::Int) = sp == 2 || sp == 10

"nc/sichg.f SICHG — per-species site age SIAGE from the site species (isisp) index ssite."
function nc_sichg(isisp::Int, ssite::Float32)::NTuple{12,Float32}
    ref_t = _nc_refloc_t(isisp)
    temsi = ssite
    smin = NC_SICHG_SIMIN[isisp]; smax = NC_SICHG_SIMAX[isisp]
    spread = smax - smin
    ntuple(12) do i
        tgt_t = _nc_refloc_t(i)
        # IDIFF: +1 if ref='B' & tgt='T'; -1 if ref='T' & tgt='B'; 0 if same basis (sichg.f:38-40)
        idiff = (!ref_t && tgt_t) ? 1 : (ref_t && !tgt_t) ? -1 : 0
        age2bh = 0f0
        if idiff != 0
            ts = temsi < smin ? smin : temsi > smax ? smax : temsi
            relsi = 100f0 * (ts - smin) / spread
            age2bh = NC_SICHG_A[i] + NC_SICHG_B[i] * relsi
        end
        NC_SICHG_IREFAG + age2bh * idiff
    end
end

# HTCALC is nc_htcalc (defined in height_growth.jl, the htcalc.f port) — the DO 30 loop evaluates the
# SITE SPECIES' (ISISP) curve at each target species' site age.

# nc/sitset.f SDIDEF fan (chunk 7): C6 = R6 per-species SDImax-ratio basis (IFOR 4/7 = Siskiyou 611 / BLM
# Coos Bay 712); C5 = the non-R6 fallback SDImax; FORMAX = SDImax cap; PMSDIU (nc/grinit.f:239) = BAMAX→%.
const NC_C6 = Float32[624,647,547,759,588,706,382,759,800,571,759,1052]
const NC_C5 = Float32[365,561,570,800,515,576,406,785,1000,365,785,1052]
const NC_FORMAX = 850f0
const NC_PMSDIU = 85f0

# nc/habtyp.f + nc/ecocls.f + nc/pvref6.f — plant-association (PA) → per-species SDImax. NC's R6 forests
# (IFOR 4/7) take the ECOCLS PA-specific SDImax; the FIA DB delivers PV_CODE as the ALPHA PA code (e.g.
# "HTS121"). Without decoding it, nc_sitset! rode a provisional per-species default (uniform ~720) instead
# of the stand's ecoclass (CWC221 default → DF 815 fanned; a real PA → its own RSDI), so the SDIMAX (and the
# self-thin it drives) was wrong on EVERY R6 NC stand — even the default. PCOML (90 R6 PA codes; KODTYP
# indexes it; default PA = CWC221 = index 46), ECOCLS (90 PA rows → site species FVSSEQ / RSDI / RSI),
# PVREF6 ((PV_CODE,PV_REF_CODE)→HABPVR crosswalk).
function _nc_load_site_tables()
    pcoml = String[]
    for l in readlines(joinpath(NC_DATADIR, "pcoml.csv"))[2:end]
        isempty(strip(l)) && continue
        push!(pcoml, String(strip(split(l, ',')[2])))
    end
    rows = NamedTuple{(:pa,:spc,:fvsseq,:sdimx,:site,:numbr,:iflag),
                      Tuple{String,String,Int,Float32,Float32,Int,Int}}[]
    for l in readlines(joinpath(NC_DATADIR, "ecocls.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(strip(l), ',')
        push!(rows, (pa=String(f[1]), spc=String(f[2]), fvsseq=parse(Int, f[3]),
                     sdimx=parse(Float32, f[4]), site=parse(Float32, f[5]),
                     numbr=parse(Int, f[6]), iflag=parse(Int, f[7])))
    end
    return pcoml, rows
end
const NC_PCOML, NC_ECOCLS = _nc_load_site_tables()

const NC_PVREF6 = let d = Dict{Tuple{String,String},String}()
    for l in readlines(joinpath(NC_DATADIR, "pvref6.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(l, ','; limit = 3)
        c = String(strip(f[1])); r = String(strip(f[2]))
        h = length(f) >= 3 ? String(strip(f[3])) : ""
        (isempty(c) || isempty(h)) && continue
        key = (c, r); haskey(d, key) || (d[key] = h)
    end
    d
end

const NC_HAB_DEFAULT_PA = "CWC221"                       # nc/habtyp.f R6 default (PCOML[46])
nc_habtyp(kodtyp::Integer)::String =
    (1 <= kodtyp <= length(NC_PCOML)) ? NC_PCOML[kodtyp] : NC_HAB_DEFAULT_PA
nc_ecocls(pa::AbstractString) = filter(r -> r.pa == pa, NC_ECOCLS)

# nc/habtyp.f (R6 path) — decode the FIA alpha PV_CODE (+ optional PV_REF_CODE) → KODTYP index into NC_PCOML.
# Ref present ⇒ PVREF6 crosswalk (blank on partial/no match ⇒ 0 ⇒ CWC221 default); no ref ⇒ the raw PV_CODE
# is HBDECD-matched directly against PCOML.
function nc_habitat_kodtyp(pv::AbstractString, pvref::AbstractString)
    pvs = String(strip(pv)); refs = String(strip(pvref))
    kard2 = pvs
    if !isempty(refs)
        kard2 = get(NC_PVREF6, (pvs, refs), "")
    end
    isempty(kard2) && return 0
    idx = findfirst(==(kard2), NC_PCOML)
    idx === nothing ? 0 : Int(idx)
end

"nc/sitset.f: site species default (DF) at 90; unset species get the SI() defaults (DF→species conversion).
SDImax (SDIDEF) — nc/sitset.f:70-172: R6 forests (IFOR 4/7) seed the site species' SDImax from the stand's
ECOCLS plant association (habitat_code → PCOML → PA; the CWC221 default when unresolved) then C6-ratio fan
capped at FORMAX; non-R6 forests use C5(i). The FIA reader now decodes the alpha PV_CODE into habitat_code."
function nc_sitset!(s::StandState)
    p = s.plot; ifor = Int(p.forest_idx)
    # nc/sitset.f:72-73 — Region-6 forests (IFOR 4 = Siskiyou 611 / IFOR 7 = BLM Coos Bay 712) reset
    # LZEIDE=.FALSE. (Reineke SDI) UNLESS a SDICALC keyword set CALCSDI (sdi_method non-blank). The reported
    # .sum SDI (disply.f branches on LZEIDE) and the SDImax self-thin both follow this, so R6 stands report
    # Reineke (QMD form) while interior forests keep the grinit Zeide default. process_keywords! runs before
    # site_setup!, so a blank sdi_method here == CALCSDI blank at the Fortran SITSET call.
    if (ifor == 4 || ifor == 7) && all(isspace, s.control.sdi_method)
        s.control.zeide_sdi = false
    end
    isisp = Int(p.site_species); isisp == 0 && (isisp = 3)     # sitset.f:119 ISISP default = 3 (DF)
    sref = p.sp_site_index[isisp]
    sref <= 0f0 && (sref = 90f0; p.sp_site_index[isisp] = 90f0) # sitset.f:120 reference SI default = 90
    # sitset.f DO 30/35 — fill each UNSET species' site index from the site species' index via the
    # SICHG (site age) + HTCALC (site-species height-at-age curve) conversion. A species whose SITEAR was
    # already set (site species / keyword / DB) keeps it (sitset.f:152 `IF(SITEAR(I).EQ.0.)`).
    siage = nc_sichg(isisp, sref)
    @inbounds for i in 1:12
        p.sp_site_index[i] <= 0f0 && (p.sp_site_index[i] = nc_htcalc(sref, isisp, siage[i]))
    end
    # R6 ECOCLS PA seed (nc/sitset.f:82-113): seed each PA species' SDImax; SDIDEF(ISISP)=RSDI on IFLAG=1.
    jsisp = 0
    if ifor == 4 || ifor == 7
        pa = nc_habtyp(Int(p.habitat_code))
        rows = nc_ecocls(pa)
        isempty(rows) && (rows = nc_ecocls(NC_HAB_DEFAULT_PA))
        @inbounds for r in rows
            iseq = r.fvsseq; (iseq < 1 || iseq > 12) && continue
            rsdi = min(r.sdimx, NC_FORMAX)
            (jsisp == 0 && r.iflag == 1) && (jsisp = iseq)
            (isisp <= 0 && r.iflag == 1) && (isisp = iseq)
            p.sp_sdi_def[iseq] <= 0f0 && (p.sp_sdi_def[iseq] = rsdi)
            (isisp > 0 && r.iflag == 1 && p.sp_sdi_def[isisp] <= 0f0) && (p.sp_sdi_def[isisp] = rsdi)
        end
    end
    p.site_species = Int32(isisp)
    # DO 40 fan (nc/sitset.f:155-172): unset species → SDIDEF(K)·C6/C6(K) capped (R6) else C5(i). (BAMAX=0.)
    k = isisp
    p.sp_sdi_def[k] <= 0f0 && (k = jsisp > 0 ? jsisp : isisp)
    @inbounds for i in 1:12
        p.sp_sdi_def[i] > 0f0 && continue
        if ifor == 4 || ifor == 7
            v = p.sp_sdi_def[k] * (NC_C6[i] / NC_C6[k])
            v > NC_FORMAX && (v = NC_FORMAX)
            p.sp_sdi_def[i] = v
        else
            p.sp_sdi_def[i] = NC_C5[i]
        end
    end
    return s
end

function nc_site_index_setup!(s::StandState)
    nc_forkod!(s.plot)
    nc_sitset!(s)
    return s
end

site_setup!(s::StandState, ::Klamath) = nc_site_index_setup!(s)
