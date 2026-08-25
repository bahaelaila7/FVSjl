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
    isisp = Int(p.site_species); isisp == 0 && (isisp = 3)     # sitset.f: ISISP default = 3 (DF)
    @inbounds for i in 1:12
        p.sp_site_index[i] <= 0f0 && (p.sp_site_index[i] = NC_SITE_DEFAULT[i])
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
