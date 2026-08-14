# =============================================================================
# site_index.jl (centralcalifornia) — CA site index + Zeide SDImax + forest code. Chunk 2.
# (ca/forkod.f, ca/habtyp.f, ca/ecocls.f, ca/sitset.f)
#
# CA's Region-6 site fan is the R6ADJ LINEAR ratio (Hann-Scrivani DF SI): HGUESS = SITEAR(ISISP)/R6ADJ(ISISP),
# then SITEAR(I)=HGUESS·R6ADJ(I) — simpler than EC's ec_htcalc curve (no sichg/htcalc needed). ECOCLS gives the
# per-plant-association site species / SI / SDImax; SDIDEF fans via BAMAX or the site-species default.
# cat01: forest 610 → IFOR 6; habitat 452 undecodable → default ecoclass CWC221 → DF(7)/SI 92/SDImx 815.
# =============================================================================

const CA_FORMAX = 850.0f0                       # ca/sitset.f DATA FORMAX/850./
# ca/sitset.f R6ADJ — per-species Hann-Scrivani DF-SI adjustment (the R6 site-fan multiplier).
const CA_R6ADJ = Float32[
    0.90, 0.70, 0.80, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.90,
    0.90, 0.90, 0.90, 0.90, 0.94, 1.00, 0.94, 0.94, 0.90, 0.90,
    0.76, 0.76, 1.00, 0.40, 0.76, 0.28, 0.42, 0.34, 0.28, 0.40,
    0.56, 0.76, 0.28, 0.76, 0.56, 0.76, 0.76, 0.76, 0.40, 0.70,
    0.40, 0.76, 0.76, 0.40, 0.76, 0.25, 0.25, 0.25, 0.56, 1.00]

# ca/forkod.f — KODFOR → IFOR via the JFOR table (region 5+6 CA forests); KFOR = all 1 (IGL=1).
const CA_JFOR = Int[505, 506, 508, 511, 514, 610, 611, 710, 711, 712, 518]

# ---- data tables (data/centralcalifornia/pcoml.csv + ecocls.csv) ----
function _ca_load_site_tables()
    pcoml = String[]
    for l in readlines(joinpath(CA_DATADIR, "pcoml.csv"))[2:end]
        isempty(strip(l)) && continue
        push!(pcoml, String(strip(split(l, ',')[2])))
    end
    rows = NamedTuple{(:pa,:spc,:fvsseq,:sdimx,:site,:numbr,:iflag),
                      Tuple{String,String,Int,Float32,Float32,Int,Int}}[]
    for l in readlines(joinpath(CA_DATADIR, "ecocls.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(strip(l), ',')
        push!(rows, (pa=String(f[1]), spc=String(f[2]), fvsseq=parse(Int, f[3]),
                     sdimx=parse(Float32, f[4]), site=parse(Float32, f[5]),
                     numbr=parse(Int, f[6]), iflag=parse(Int, f[7])))
    end
    return pcoml, rows
end
const CA_PCOML, CA_ECOCLS = _ca_load_site_tables()

# ca/forkod.f — set IFOR from KODFOR (cat01 610 → IFOR 6); default IFOR=1 if the code isn't in JFOR.
function ca_forkod!(p)
    kodfor = Int(p.user_forest_code)
    idx = findfirst(==(kodfor), CA_JFOR)
    ifor = idx === nothing ? 1 : idx
    p.forest_idx = Int32(ifor)
    p.geo_location = Int32(1)                     # ca/forkod.f KFOR = all 1
    p.user_forest_code = Int32(CA_JFOR[ifor])
    return ifor
end

# ca/habtyp.f (R6, KODFOR≥600) — KODTYP → PCOM plant-association code. A 1..90 sequence indexes PCOML;
# anything else (e.g. cat01's 452) is undecodable → habtyp default PA = CWC221 (matches ca/sitset.f).
ca_habtyp(kodtyp::Integer)::String =
    (1 <= kodtyp <= length(CA_PCOML)) ? CA_PCOML[kodtyp] : "CWC221"

ca_ecocls(pa::AbstractString) = filter(r -> r.pa == pa, CA_ECOCLS)

# ca/sitset.f — Region-6 path (IFOR≥6): ECOCLS site-index fan via the R6ADJ ratio; SDIDEF via BAMAX or
# the site-species default. (R5 IFOR<6 path with R5ADJ/R5SDI is not exercised by the CA ref stands.)
function ca_sitset!(s::StandState)
    p = s.plot; maxsp = nspecies(s.variant)
    formax = CA_FORMAX
    nsiset = count(>(0f0), @view p.sp_site_index[1:maxsp])
    pcom = ca_habtyp(Int(p.habitat_code))         # R6 ecoclass code (default CWC221 when undecodable)
    isempty(pcom) && (pcom = "CWC221")
    rows = ca_ecocls(pcom)
    isempty(rows) && (rows = ca_ecocls("CWC221")) # ca/sitset.f ICL5==0 default when PA not found

    isisp = (1 <= Int(p.site_species) <= maxsp) ? Int(p.site_species) : 0
    jsisp = 0
    @inbounds for r in rows
        iseq = r.fvsseq; iseq == 0 && continue
        rsdi = min(r.sdimx, formax)
        (jsisp == 0 && r.iflag == 1) && (jsisp = iseq)
        (isisp <= 0 && r.iflag == 1) && (isisp = iseq)
        (p.sp_site_index[iseq] <= 0f0 && nsiset == 0) && (p.sp_site_index[iseq] = r.site)
        p.sp_sdi_def[iseq] <= 0f0 && (p.sp_sdi_def[iseq] = rsdi)
        (isisp > 0 && r.iflag == 1 && p.sp_sdi_def[isisp] <= 0f0) && (p.sp_sdi_def[isisp] = rsdi)
    end
    isisp <= 0 && (isisp = 7)                      # ca/sitset.f R6 global default site sp = DF(7)
    p.sp_site_index[isisp] <= 0f0 && (p.sp_site_index[isisp] = 80f0)   # ca/sitset.f default SI

    # ca/sitset.f: HGUESS = Hann-Scrivani DF SI = SITEAR(ISISP)/R6ADJ(ISISP); fan every unset species.
    hguess = p.sp_site_index[isisp] / CA_R6ADJ[isisp]
    @inbounds for i in 1:maxsp
        p.sp_site_index[i] == 0f0 && (p.sp_site_index[i] = hguess * CA_R6ADJ[i])
    end

    # ca/sitset.f DO 40 — SDIDEF fan: BAMAX branch else the site-species default (K), capped at FORMAX.
    k = isisp
    p.sp_sdi_def[k] <= 0f0 && (k = jsisp > 0 ? jsisp : isisp)
    bamax = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi                 # ca/grinit.f PMSDIU (only used when BAMAX>0)
    @inbounds for i in 1:maxsp
        p.sp_sdi_def[i] > 0f0 && continue
        v = bamax > 0f0 ? bamax / (0.5454154f0 * (pmsdiu / 100f0)) : p.sp_sdi_def[k]
        v > formax && (v = formax)
        p.sp_sdi_def[i] = v
    end
    p.site_species = Int32(isisp)
    return s
end

function ca_site_index_setup!(s::StandState)
    ca_forkod!(s.plot)
    ca_sitset!(s)
    return s
end

site_setup!(s::StandState, ::CentralCalifornia) = ca_site_index_setup!(s)
