# cover.jl — FVS COVER understory-vegetation / canopy-cover extension (covr/ vcovr/ pg/).
#
# SCOPE (see scratchpad/cover/HANDOFF.md for full evidence):
#   COVER is a REPORT-ONLY extension. It reads the tree list + stand density and
#   emits a text report (CANOPY COVER STATISTICS / SHRUB STATISTICS / SUMMARY /
#   SHRUB-SMALL CONIFER COMPETITION) to output unit JOSHRB at end of projection.
#   It does NOT modify tree DBH/HT/PROB/growth/mortality — grep of covr/*.f
#   vcovr/*.f shows ZERO writes to any tree array (verified). Therefore:
#     - INERT  = a run without the COVER keyword is byte-identical (LCVATV=.FALSE.).
#     - VALID  = the emitted report text matches the live oracle FVSem_g16.
#   COVER draws NO random numbers.  There is NO RNG to port.
#
# LINKAGE: real cvin.o (active) is compiled into ALL 18 western variant binaries.
#   Oracle needs NO relink — the shipped western binary fires COVER on keyword 46.
#
# THIS FILE PORTS THE CANOPY-COVER HALF:
#   * CoverState      — CVCOM subset (canopy path)
#   * cover_init!      (CVINIT)         — per-stand reset
#   * kw_coverin!      (CVIN)           — COVER…END keyword-block reader (report-only)
#   * cover_cvcw!      (CVCW)           — crown-area, DUMP-REPLAY BIT-EXACT vs FVSem_g16
#   * cover_accumulate! (CVCNOP driver: CVCW+CVSHAP+CVCBMS+CVSUM) per cycle
#       - _cover_cvshap  (CVSHAP)       — DUMP-REPLAY 351/351 ISHAPE bit-exact
#       - _cover_cvcbms2 (CVCBMS opt 2) — DUMP-REPLAY 337/351 hex, 14 @1-ULP (libm exp/log)
#       - _cover_cvsum!  (CVSUM)        — 10-ft height-class canopy geometry
#   * cover_report      (CVOUT)         — "CANOPY COVER STATISTICS" table emit
#   The shrub half (CVBROW/CVSCON/CVBCAL/CVCLAS) and pg/ serialization are DEFERRED.
#
# VALIDATION: the emitted CANOPY COVER STATISTICS integer table is byte-identical to
#   FVSem_g16 (EM emt01, 3 cycles) — see test/unit/test_cover.jl and HANDOFF.md.

# =====================================================================================
# CVSHAP crown-shape discriminant coefficients (Moeur 1983; original 11 NI species).
# Fortran DATA CONST(5,11) etc. are column-major → each source line = one species'
# 5 shape coefficients.  Stored here species-major: COEF[species, shape].
# =====================================================================================
const _CV_CONST = Float32[
 -13.943 -37.395 -99.000 -32.104 -32.042
 -99.000 -38.278 -99.000 -44.875 -32.742
 -25.380 -31.308 -35.906 -31.435 -28.185
 -19.633 -23.785 -20.062 -22.012 -22.962
 -19.633 -23.785 -20.062 -22.012 -22.962
 -99.000 -27.573 -99.000 -22.434 -26.275
 -22.322 -37.190 -99.000 -38.316 -31.988
 -99.000 -46.116 -99.000 -49.990 -43.456
 -24.002 -37.611 -99.000 -30.234 -36.183
 -20.296 -24.195 -36.526 -26.622 -22.397
 -24.002 -37.611 -99.000 -30.234 -36.183 ]
const _CV_BCR = Float32[
 49.794 100.240 -99.00 86.383 94.554
 -99.00 123.805 -99.00 134.780 116.554
 52.715 60.761 63.156 61.744 59.922
 45.992 54.464 38.042 51.755 53.598
 45.992 54.464 38.042 51.755 53.598
 -99.00 59.560 -99.00 53.273 55.799
 73.920 104.413 -99.00 102.380 99.877
 -99.00 106.784 -99.00 110.115 100.844
 56.463 80.909 -99.00 67.587 78.258
 45.808 57.540 60.760 59.567 57.457
 56.463 80.909 -99.00 67.587 78.258 ]
const _CV_BHT = Float32[
 0.48348 0.92054 -99.000 0.81786 0.85462
 -99.000 1.03836 -99.000 1.00712 1.02891
 1.02304 1.13332 1.28061 1.13364 1.16673
 0.72015 0.88508 0.56218 0.84169 0.88962
 0.72015 0.88508 0.56218 0.84169 0.88962
 -99.000 1.08626 -99.000 0.98588 0.98314
 0.83497 1.16782 -99.000 1.17028 1.15645
 -99.000 1.60533 -99.000 1.62982 1.59203
 1.25278 1.59728 -99.000 1.38587 1.58670
 0.71361 0.91473 0.98645 0.89231 0.91047
 1.25278 1.59728 -99.000 1.38587 1.58670 ]
const _CV_BRAD = Float32[
 0.56507 1.11080 -99.000 1.82564 1.32317
 -99.000 1.08869 -99.000 1.81120 0.79969
 -0.10457 0.28897 -0.49614 0.07177 -0.04803
 0.13463 0.35420 0.24447 0.32306 0.18187
 0.13463 0.35420 0.24447 0.32306 0.18187
 -99.000 0.92215 -99.000 0.91886 0.48536
 0.62301 1.69840 -99.000 1.57036 1.08601
 -99.000 1.12257 -99.000 0.86087 0.89170
 1.64765 1.72301 -99.000 1.83248 1.53816
 0.98573 0.60433 2.22885 0.89341 0.58180
 1.64765 1.72301 -99.000 1.83248 1.53816 ]
const _CV_BCL = Float32[
 -0.8373 -1.5604 -99.000 -1.3933 -1.4959
 -99.000 -1.9316 -99.000 -2.1712 -1.8977
 -1.3082 -1.4288 -1.5055 -1.4393 -1.4367
 -1.1754 -1.3304 -0.9574 -1.3074 -1.3250
 -1.1754 -1.3304 -0.9574 -1.3074 -1.3250
 -99.000 -1.2529 -99.000 -1.1295 -1.1502
 -1.6455 -2.2164 -99.000 -2.1892 -2.1148
 -99.000 -2.2233 -99.000 -2.3029 -2.1475
 -1.6686 -2.0343 -99.000 -1.8193 -2.0154
 -1.0722 -1.3127 -1.3903 -1.3079 -1.2973
 -1.6686 -2.0343 -99.000 -1.8193 -2.0154 ]
const _CV_BTPA = Float32[
 0.00437 0.00911 -99.000 0.00603 0.00899
 -99.000 0.01548 -99.000 0.02342 0.02134
 0.00634 0.00706 0.00673 0.00683 0.00673
 0.00782 0.00841 0.00813 0.00778 0.00816
 0.00782 0.00841 0.00813 0.00778 0.00816
 -99.000 0.00896 -99.000 0.00736 0.00810
 0.00600 0.00406 -99.000 0.00542 0.00613
 -99.000 0.04348 -99.000 0.04908 0.05410
 0.01644 0.01891 -99.000 0.01809 0.01923
 0.01233 0.01068 0.01136 0.01190 0.01070
 0.01644 0.01891 -99.000 0.01809 0.01923 ]
const _CV_BDBH = Float32[
 0.05236 0.08377 -99.0000 -0.17327 0.22765
 -99.0000 -0.14770 -99.0000 0.25077 0.01210
 -0.40495 -0.66563 -0.42355 -0.50661 -0.54780
 0.14819 -0.02081 0.14503 0.14760 0.00202
 0.14819 -0.02081 0.14503 0.14760 0.00202
 -99.0000 -0.81991 -99.0000 -0.99031 -0.33308
 0.21973 -0.17842 -99.0000 0.00182 -0.04777
 -99.0000 0.06301 -99.0000 0.41333 0.02430
 -0.61101 -0.80165 -99.0000 -0.74084 -0.72773
 -0.39061 -0.33475 -0.99571 -0.40740 -0.37839
 -0.61101 -0.80165 -99.0000 -0.74084 -0.72773 ]

# CVCBMS foliage-biomass coefficients (Moeur 1981; original 11 species)
const _CV_BINT11 = Float32[-2.15894,-5.02156,-2.30430,-2.78090,-4.22701,-2.64034,-3.38394,-3.30673,-2.03919,-3.02050,-2.81317]
const _CV_BCL11  = Float32[1.48969,2.31835,1.52896,1.90272,2.22534,1.69973,1.96060,2.27613,1.64942,1.88712,1.47513]
const _CV_BINT12 = Float32[-1.94951,-4.73762,-2.05828,-2.43200,-4.17456,-2.24876,-3.13488,-2.93508,-1.60998,-2.74410,-2.63387]
const _CV_BCL12  = Float32[1.22023,1.98479,1.25837,1.60270,2.00749,1.37600,1.62368,1.96125,1.32649,1.58171,1.35092]
const _CV_BINT2  = Float32[2.666072,1.756537,2.705866,3.115084,2.654572,3.059351,2.622505,3.300852,3.060169,2.452492,2.622505]

# Variant → original-11-species maps (cvshap.f / cvcbms.f), indexed by variant species
# number (1-based).  Only variants needed for the current beachhead are populated;
# add the remaining DATA MAP** blocks from vcovr/cvshap.f as later variants are wired.
const _CV_MAP_NI = Int[1,2,3,4,5,6,7,8,9,10,11,6]                       # KT, NC (original 11)
const _CV_MAP_EM = Int[1,2,3,1,2,6,7,8,9,10, 2,2,2,2,2,2,2,11,2]        # EM (19 sp)

# Return the variant→11-species map for the given variant code, or `nothing`
# (→ default map-to-1, matching cvshap.f CASE DEFAULT) if not yet transcribed.
function _cover_spmap(varcode::AbstractString)
    varcode == "EM"            && return _CV_MAP_EM
    (varcode == "KT" || varcode == "NC") && return _CV_MAP_NI
    return nothing
end

# CRWDTH source (base/cwidth.f → CWCALC).  FVS fills the CRWDTH array — read verbatim by
# CVCW — from the *variant's* forest-grown crown-width equation (IWHO=0, actual CR, stand
# BA).  For EM this is the western CWCALC crown-width library `em_cwcalc`
# (variants/easternmontana/crown.jl), dump-replay bit-exact vs FVSem_g16.

# Dispatch the CRWDTH source by variant.  EM = the western CWCALC crown-width library
# (em_cwcalc, base/cwidth.f→cwcalc.f IWHO=0, all 15 EM CWEQN forms) — dump-replay bit-exact
# vs FVSem_g16.  Other variants need their own CWCALC map (later chunks).
@inline function _cover_crwdth(s::StandState, ispc::Int, d::Float32, h::Float32,
                               jcr::Int, ba::Float32)
    if s.variant isa EasternMontana
        el = s.plot.elevation                                   # EL = ELEV (100's of ft)
        hi = _cr_hopkins(s.plot.latitude, s.plot.longitude, s.plot.elevation)  # western Hopkins
        return em_cwcalc(ispc, d, h, Float32(jcr), ba, el, hi)  # CR = FLOAT(ICR)
    end
    return 0.0f0
end

# =====================================================================================
# CoverState — CVCOM subset for the canopy-cover report path.
# =====================================================================================
"""
    CoverRow

One projection-cycle's canopy-cover statistics (CVSUM output for slot IP1, ITHN=1),
in FVS units. Emitted by `cover_report` as one "CANOPY COVER STATISTICS" block.
"""
struct CoverRow
    year::Int
    txht::Vector{Float32}    # TXHT   : trees/ac by 10-ft height class (16)
    pcxht::Vector{Float32}   # PCXHT  : % canopy cover by class (overlap-adjusted)
    proxht::Vector{Float32}  # PROXHT : crown profile area by class
    volxht::Vector{Float32}  # VOLXHT : crown volume by class
    cfbxht::Vector{Float32}  # CFBXHT : foliage biomass by class
    tretot::Float32          # TRETOT : total trees/ac
    tpctcv::Float32          # TPCTCV : total % canopy cover
    tproar::Float32          # TPROAR : total profile area
    tcvol::Float32           # TCVOL  : total crown volume
    totbms::Float32          # TOTBMS : total foliage biomass
    stdht::Float32           # STDHT  : top/avg height (=AVH)
    sdiam::Float32           # SDIAM  : sum of stem diameters (ft)
    icvage::Int              # ICVAGE : overstory age
end

# Shrub-calibration state (covr/cvbcal.f + the SHRBLAYR/SHRUBHT/SHRUBPC keyword inputs).
# Populated at keyword-parse time; BHTCF/BPCCF (+ report intermediates + carried RESIDC/
# RESIDH) computed once by cvbcal! at cycle 0 and applied to the shrub predictions each
# cycle while calibration is in effect (LCALIB).  Defaults per covr/cvinit.f.
mutable struct ShrubCalib
    lcal1::Bool                    # SHRBLAYR active (calibration by layer)
    lcal2::Bool                    # SHRUBHT/SHRUBPC active (calibration by species)
    lcalib::Bool                   # calibration still in effect (off after a thin)
    avgbht::Vector{Float32}        # (3) observed layer heights, sorted tallest→shortest
    avgbpc::Vector{Float32}        # (3) observed layer covers
    nklass::Int                    # number of observed layers
    sumcvr::Float32                # total observed cover (Σ AVGBPC over input fields)
    shrbht::Vector{Float32}        # (31) observed species heights (-99999 = not input)
    shrbpc::Vector{Float32}        # (31) observed species covers  (-99999 = not input)
    bhtcf::Vector{Float32}         # (31) height correction factors (default 1)
    bpccf::Vector{Float32}         # (31) cover  correction factors (default 1)
    residc::Float32                # RESIDC (SUMCVR-TCOV, method 1) carried from cyc0
    residh::Vector{Float32}        # (31) RESIDH (SHRBHT-SH, method 2) carried from cyc0
    computed::Bool                 # cvbcal! has run (factors fixed at cyc0)
    # report intermediates (cvout.f calibration tables)
    xsh::Vector{Float32}           # (31) uncalibrated SH at cyc0
    xcv::Vector{Float32}           # (31) uncalibrated PBCV at cyc0
    xpb::Vector{Float32}           # (31) uncalibrated PB at cyc0
    htavg::Vector{Float32}         # (3) predicted mean layer height
    cvavg::Vector{Float32}         # (3) predicted layer cover
    htfrac::Vector{Float32}        # (3) height scaling factor by layer
    cvfrac::Vector{Float32}        # (3) cover  scaling factor by layer
    ilayr::Vector{Int}             # (31) layer assigned to each species
end

ShrubCalib() = ShrubCalib(false, false, false,
    fill(-99999.0f0, 3), fill(-99999.0f0, 3), 0, 0.0f0,
    fill(-99999.0f0, 31), fill(-99999.0f0, 31), ones(Float32, 31), ones(Float32, 31),
    0.0f0, zeros(Float32, 31), false,
    zeros(Float32, 31), zeros(Float32, 31), zeros(Float32, 31),
    zeros(Float32, 3), zeros(Float32, 3), zeros(Float32, 3), zeros(Float32, 3),
    zeros(Int, 31))

mutable struct CoverState <: AbstractCoverState
    active::Bool          # LCVATV: COVER activity scheduled → routinely called
    lcov::Bool            # LCOV:   activity has fired this/earlier cycle
    lcnop::Bool           # LCNOP:  CANOPY sub-keyword (crown shape+biomass) on
    lbrow::Bool           # LBROW:  SHRUBS sub-keyword on (shrub half — deferred)
    lcover::Bool          # LCOVER: emit CANOPY COVER STATISTICS table  (NOCOVOUT→false)
    lshrub::Bool          # LSHRUB: emit SHRUB STATISTICS table         (NOSHBOUT→false)
    lcvsum::Bool          # LCVSUM: emit SUMMARY table                  (NOSUMOUT→false)
    covopt::Int           # COVOPT: foliage-biomass eqn option (default 2)
    idt::Int              # IDT:    cycle/date the COVER activity turns on
    joshrb::Int           # JOSHRB: output unit for the report (default = main out)
    crarea::Float32       # CRAREA: sum of crown projection areas (sq ft/ac), CVCW out
    trecw::Vector{Float32}# TRECW(MAXTRE): per-tree crown width used by cover (=CRWDTH)
    rows::Vector{CoverRow}# per-cycle accumulated canopy rows (CVSUM output)
    # -- SHRUBS half (report-only) --
    sage0::Float32        # SAGE:  time since disturbance from SHRUBS card (-1 ⇒ use IAGE)
    ihtype::Int           # IHTYPE: habitat code from SHRUBS card (0 ⇒ stand habitat ICL5)
    iphys::Int            # IPHYS:  physiographic position (default 2 = lower slope)
    idist::Int            # IDIST:  disturbance type (default 1 = none)
    shrub_const::Any      # ShrubConst (CVSCON output) computed at cycle 0, else nothing
    shrub_rows::Vector    # per-cycle ShrubRow (CVBROW/CVSUM/CVCLAS output)
    calib::ShrubCalib     # shrub-calibration inputs/outputs (covr/cvbcal.f)
end

CoverState() = CoverState(false, false, false, false, true, true, true,
                          2, 1, 0, 0.0f0, Float32[], CoverRow[],
                          -1.0f0, 0, 2, 1, nothing, Any[], ShrubCalib())

# CVINIT (covr/cvinit.f): initialise per-stand cover flags/arrays.
function cover_init!(cv::CoverState)
    cv.lcov   = false
    cv.lcnop  = false
    cv.lbrow  = false
    cv.lcover = true
    cv.lshrub = true
    cv.lcvsum = true
    cv.crarea = 0.0f0
    empty!(cv.trecw)
    empty!(cv.rows)
    cv.sage0  = -1.0f0
    cv.ihtype = 0
    cv.iphys  = 2
    cv.idist  = 1
    cv.shrub_const = nothing
    empty!(cv.shrub_rows)
    cv.calib = ShrubCalib()
    return cv
end

# =====================================================================================
# CVIN keyword-block reader.
# =====================================================================================
"""
    kw_coverin!(s, rec, kr)

CVIN (covr/cvin.f): parse the `COVER … END` keyword block.  Report-only — sets the
CoverState flags and schedules the report; never perturbs the tree record, growth,
or mortality (INERT).  Sub-keywords handled: CANOPY, SHRUBS, NOCOVOUT, NOSHBOUT,
NOSUMOUT, COVER (activity-900 schedule → `active`), SHOWSHRB, CVNOHEAD, DEBUG, END.
The shrub-calibration sub-keywords (SHRBLAYR/SHRUBHT/SHRUBPC) are consumed but their
calibration effect is deferred with the shrub half.
"""
function kw_coverin!(s::StandState, rec, kr)
    s.cover === nothing && (s.cover = CoverState())
    cv = s.cover
    cover_init!(cv)
    # The base card that entered CVIN is the keyword 'COVER' itself, which cvin.f finds
    # at TABLE(12) and processes as option 12 (schedule activity 900). Mirror that on
    # `rec` so the report fires even with no inner COVER sub-card.
    cv.active = true
    cv.lcov   = true
    cv.idt    = (length(rec.present) >= 1 && rec.present[1] && rec.values[1] > 0.0f0) ?
                Int(trunc(rec.values[1])) : 1
    (length(rec.present) >= 2 && rec.present[2]) && (cv.joshrb = Int(trunc(rec.values[2])))
    # The block body follows as sub-keyword records read from `kr` until END/EOF/STOP.
    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            break
        elseif k == "CANOPY"
            cv.lcnop  = true
            cv.covopt = 2
            r.present[1] && (cv.covopt = Int(trunc(r.values[1])))
        elseif k == "SHRUBS"          # SHRUBS: field1 SAGE, 2 IHTYPE, 3 IPHYS, 4 IDIST
            cv.lbrow  = true
            cv.sage0  = (r.present[1] ? r.values[1] : -1.0f0)
            cv.ihtype = (r.present[2] ? Int(trunc(r.values[2])) : 0)
            cv.iphys  = (r.present[3] ? Int(trunc(r.values[3])) : 2)
            cv.idist  = (r.present[4] ? Int(trunc(r.values[4])) : 1)
        elseif k == "SHRBLAYR"
            _cvin_shrblayr!(cv.calib, r)          # calibration by shrub layer (cvin.f opt 4)
        elseif k == "SHRUBHT"
            _cvin_shrub_species!(cv.calib, kr, true)   # observed heights (cvin.f opt 5)
        elseif k == "SHRUBPC"
            _cvin_shrub_species!(cv.calib, kr, false)  # observed covers  (cvin.f opt 6)
        elseif k == "NOCOVOUT"
            cv.lcover = false
        elseif k == "NOSHBOUT"
            cv.lshrub = false
        elseif k == "NOSUMOUT"
            cv.lcvsum = false
        elseif k == "COVER"           # schedule activity 900 → report is emitted
            cv.active = true
            cv.lcov   = true
            cv.idt    = r.present[1] && r.values[1] > 0.0f0 ? Int(trunc(r.values[1])) : 1
            r.present[2] && (cv.joshrb = Int(trunc(r.values[2])))
        elseif k == "SHOWSHRB"
            read_keyword!(kr)          # one species-selection data record (deferred)
        elseif k == "CVNOHEAD" || k == "DEBUG"
            # no-op for the canopy report path
        else
            (!isempty(k) && isletter(first(k))) &&
                push!(s.control.unrecognized_keywords, k)
        end
    end
    return nothing
end

# SHRBLAYR (cvin.f opt 4): observed shrub-layer heights/covers → AVGBHT/AVGBPC/NKLASS/
# SUMCVR, bubble-sorted by decreasing height.  Fields 1,3,5 = heights, 2,4,6 = covers.
function _cvin_shrblayr!(cal::ShrubCalib, r)
    cal.lcal1 = true
    cal.lcalib = true
    getv(i) = (i <= length(r.values) && r.present[i]) ? r.values[i] : 0.0f0
    pres(i) = i <= length(r.present) && r.present[i]
    cal.avgbht[1] = getv(1); cal.avgbpc[1] = getv(2)
    cal.avgbht[2] = getv(3); cal.avgbpc[2] = getv(4)
    cal.avgbht[3] = getv(5); cal.avgbpc[3] = getv(6)
    # count observed layers (non-blank COVER fields) and total cover (blanks = 0)
    cal.nklass = 0
    cal.sumcvr = 0.0f0
    for i in (2, 4, 6)
        pres(i) && (cal.nklass += 1)
        cal.sumcvr += getv(i)
    end
    # bubble-sort AVGBHT/AVGBPC by decreasing height (cvin.f DO 120/130)
    while true
        lsort = false
        for i in 1:2
            j = i + 1
            if cal.avgbht[i] < cal.avgbht[j]
                lsort = true
                cal.avgbht[i], cal.avgbht[j] = cal.avgbht[j], cal.avgbht[i]
                cal.avgbpc[i], cal.avgbpc[j] = cal.avgbpc[j], cal.avgbpc[i]
            end
        end
        lsort || break
    end
    return cal
end

# Shrub species abbreviation → index (CH4BSR over the 31 SNAME + blank(32) + "-999"(33)).
function _cvb_species_index(abbrev::AbstractString)
    a = uppercase(strip(abbrev))
    isempty(a) && return 32
    a == "-999" && return 33
    idx = findfirst(==(a), _CV_SNAME)      # _CV_SNAME[1:31] are the species; 32 = OTHR
    (idx === nothing || idx > 31) && return 0
    return idx
end

# SHRUBHT / SHRUBPC (cvin.f opt 5/6): read up to 4 data records, format 8(A4,F6.1);
# set SHRBHT (isheight) or SHRBPC per species.  "-999" terminates.  Blank abbrev keeps
# the -99999 dummy.  cvin.f leaves the -99999 for species not listed.
function _cvin_shrub_species!(cal::ShrubCalib, kr, isheight::Bool)
    cal.lcal2 = true
    cal.lcalib = true
    target = isheight ? cal.shrbht : cal.shrbpc
    for _ in 1:4
        rr = read_keyword!(kr)
        (rr.status == KW_EOF || rr.status == KW_STOP) && return cal
        line = rpad(rr.raw, 80)
        done = false
        for p in 0:7                       # 8 pairs of (A4, F6.1) = 10 cols each
            c = p * 10
            abbrev = line[c+1:c+4]
            num = _cvb_species_index(abbrev)
            if num == 33                   # "-999" → end of data
                done = true
                break
            end
            num == 0 && continue           # unknown abbrev (ERRGRO) → skip
            num == 32 && continue          # blank abbrev → keep dummy
            valstr = strip(line[c+5:c+10])
            v = isempty(valstr) ? 0.0f0 : (tryparse(Float32, valstr))
            v === nothing && (v = 0.0f0)
            target[num] = v
        end
        done && return cal
    end
    return cal
end

# =====================================================================================
# CVCW (covr/cvcw.f) — per-tree crown width & stand crown area.  DUMP-REPLAY BIT-EXACT.
# =====================================================================================
"""
    cover_cvcw!(cv, crown_width, prob, itrn) -> Float32

CRAREA = (Σ TRECW(i)²·PROB(i))·0.785398, TRECW(I)=CRWDTH(I).  Float32 running sum in
tree order; literal 0.785398f0 — both load-bearing.  Validated bit-exact (Float32 hex)
vs FVSem_g16 over 3 cycles.
"""
function cover_cvcw!(cv::CoverState, crown_width::AbstractVector{Float32},
                     prob::AbstractVector{Float32}, itrn::Integer)
    resize!(cv.trecw, itrn)
    acc = 0.0f0
    @inbounds for i in 1:itrn
        cw = crown_width[i]
        cv.trecw[i] = cw
        acc = acc + cw * cw * prob[i]
    end
    cv.crarea = acc * 0.785398f0
    return cv.crarea
end

# CVSHAP (vcovr/cvshap.f): discriminant argmax → crown shape 1..5.  ispi = mapped index.
@inline function _cover_cvshap(ispi::Int, dbh::Float32, ht::Float32, cr::Float32,
                               rad::Float32, tpa::Float32)
    cl = cr * ht
    scorm1 = -999.0f0
    ishape = 0
    @inbounds for j in 1:5
        score = _CV_CONST[ispi,j] + _CV_BCR[ispi,j]*cr + _CV_BHT[ispi,j]*ht +
                _CV_BCL[ispi,j]*cl + _CV_BRAD[ispi,j]*rad +
                _CV_BDBH[ispi,j]*dbh + _CV_BTPA[ispi,j]*tpa
        if score > scorm1
            ishape = j
            scorm1 = score
        end
    end
    return ishape
end

# CVCBMS (vcovr/cvcbms.f) COVOPT=2 foliage biomass (lbs).
@inline function _cover_cvcbms2(ispi::Int, d::Float32, h::Float32, cl::Float32,
                                rd::Float32, dds::Float32, tpa::Float32)
    alntpa = log(tpa)
    if d < 3.5f0
        trf = exp(_CV_BINT12[ispi] + _CV_BCL12[ispi]*log(cl)
                  - 0.12975f0*alntpa + 0.40350f0*log(h))
        trf = trf * 1.13178f0
    else
        trf = exp(_CV_BINT2[ispi] + 1.468547f0*log(d) + 0.308847f0*log(dds)
                  - 1.077047f0*log(h) + 0.690825f0*log(cl)
                  - 0.142096f0*alntpa + 0.399244f0*log(rd))
    end
    return trf
end

# CVSUM (covr/cvsum.f) per-tree canopy binning: crown profile area / volume / foliage
# biomass by 10-ft height class, computed as a sum of crown frustums, plus the crown
# projection-area (CRXHT) and tree-count (TXHT) tallies by top-height class.  Bit-exact
# frustum geometry — DUMP-REPLAY 239/240 vs FVSem_g16, 1 value @1-ULP (see test_cover.jl).
function _cover_cvsum_tree!(txht, crxht, proxht, volxht, cfbxht,
                            ishap::Int, cw::Float32, dbh::Float32, ht::Float32,
                            icr::Int, prob::Float32, trfbms::Float32)
    rad = cw * 0.5f0
    itop = trunc(Int, ht / 10.0001f0) + 1
    itop > 16 && (itop = 16)
    @inbounds txht[itop]  += prob
    @inbounds crxht[itop] += 0.785398f0*cw*cw*prob    # CRXHT → PCXHT
    bot  = ht - (ht * Float32(icr) / 100.0f0)
    ibot = trunc(Int, bot / 10.0f0) + 1
    hc   = Float32(icr) * ht / 100.0f0
    base = bot
    (ishap == 1 || ishap == 5) && (base = bot + hc/2.0f0)
    cnop1 = 10.0f0*ibot - 20.0f0
    cnop2 = cnop1 + 10.0f0
    j1 = itop + 1
    @inbounds for j in ibot:itop
        if ishap != 3
            j1 = j
        else
            j1 = j1 - 1
        end
        cnop1 += 10.0f0
        cnop2 += 10.0f0
        uplim = min(ht, cnop2)
        j == itop && (uplim = ht)
        lowlim = max(bot, cnop1)
        h1 = lowlim - base
        h2 = uplim - lowlim
        y1 = lowlim - base
        y2 = uplim - base
        y1 >= hc && (y1 = hc); y2 >= hc && (y2 = hc)
        y1 < -hc && (y1 = -hc); y2 < -hc && (y2 = -hc)
        r1 = (1.0f0 - h1/hc)*rad
        r2 = (1.0f0 - (h1+h2)/hc)*rad
        local frust::Float32, parea::Float32, cvolum::Float32
        if ishap == 1 || ishap == 5           # sphere / ellipsoid
            cvolum = 2.09439f0*rad*rad*hc
            base >= lowlim && (h1 = base - uplim)
            b1 = hc/2.0f0
            y1 > b1 && (y1 = b1); y1 < -b1 && (y1 = -b1)
            y2 > b1 && (y2 = b1); y2 < -b1 && (y2 = -b1)
            z1 = b1*b1 - y1*y1; z2 = b1*b1 - y2*y2
            z1 < 0.0f0 && (z1 = 0.0f0); z2 < 0.0f0 && (z2 = 0.0f0)
            cst = 1.04720f0*h2*rad*rad/(hc*hc)
            frust = cst*(3f0*hc*hc - 12f0*h1*h1 - 12f0*h1*h2 - 4f0*h2*h2)
            parea = rad/b1*(y2*sqrt(z2)+b1^2*asin(y2/b1)) -
                    rad/b1*(y1*sqrt(z1)+b1^2*asin(y1/b1))
        elseif ishap == 2                      # cone
            cvolum = 1.04719f0*rad*rad*hc
            cst = 1.04720f0*h2
            frust = cst*(r1*r1 + r2*r2 + r1*r2)
            parea = (r1+r2)*h2
        elseif ishap == 3                      # neiloid
            cvolum = 1.57079f0*rad*rad*hc
            frust = (3.14159f0*rad*rad*h2) -
                    (1.57079f0*h2*rad^2/hc)*(2f0*hc - 2f0*h1 - h2)
            parea = rad*((y2 + 0.66667f0*hc*(1f0-y2/hc)^1.5f0) -
                         (y1 + 0.66667f0*hc*(1f0-y1/hc)^1.5f0))
        else                                   # ishap == 4  paraboloid
            cvolum = 1.57079f0*rad*rad*hc
            frust = (1.57079f0*h2*rad^2/hc)*(2f0*hc - 2f0*h1 - h2)
            parea = -rad*1.33333f0*hc*((1f0-(y2/hc))^1.5f0 - (1f0-(y1/hc))^1.5f0)
        end
        prop = frust/cvolum
        proxht[j1] += parea*prob
        volxht[j1] += frust*prob
        cfbxht[j1] += trfbms*prop*prob
    end
    return nothing
end

# =====================================================================================
# cover_accumulate! — the CVCNOP driver (CVCW → CVSHAP → CVCBMS → CVSUM) for one cycle.
# Report-only: reads the tree list + stand stats, appends one CoverRow.  ITHN=1
# (pre-thin) only; the post-thin (LTHIN/ITHN=2) path is deferred with the thin seam.
# =====================================================================================
function cover_accumulate!(cv::CoverState, s::StandState, year::Integer, year0::Integer,
                           fint::Real)
    t = s.trees
    itrn = t.n
    spmap = _cover_spmap(variant_code(s.variant))
    ba    = s.plot.basal_area           # BA (for CRWDTH ccfcal MODE=2)
    fint  = Float32(fint)               # FINT (cycle length; CVCBMS DDS denominator)
    avh   = s.plot.avg_height           # AVH → STDHT
    iage  = Int(s.plot.stand_age)       # IAGE
    # TPROB = Σ PROB over the live records; RMSQD = the stand quadratic mean diameter.
    tpa = 0.0f0
    @inbounds for i in 1:itrn
        tpa += t.tpa[i]
    end
    rmsqd = stand_qmd(s)                 # RMSQD

    txht   = zeros(Float32, 16); pcxht  = zeros(Float32, 16)
    proxht = zeros(Float32, 16); volxht = zeros(Float32, 16)
    cfbxht = zeros(Float32, 16); crxht  = zeros(Float32, 16)
    sd2xht = zeros(Float32, 16)          # SD2XHT: Σ DBH² by top-height class (for CVCLAS)
    trsh   = zeros(Float32, 11)          # TRSH: trees/ac by SHTRHT height threshold (shrub-conifer table)
    htmax  = 0.0f0; htmin = 999.0f0
    sdiam = 0.0f0

    @inbounds for i in 1:itrn
        dbh  = t.dbh[i]
        ht   = t.height[i]
        icr  = Int(t.crown_pct[i])       # ICR (integer %)
        prob = t.tpa[i]
        dg   = t.diam_growth[i]
        itopc = trunc(Int, ht / 10.0001f0) + 1; itopc > 16 && (itopc = 16)
        sd2xht[itopc] += dbh*dbh*prob
        ht > htmax && (htmax = ht); ht < htmin && (htmin = ht)
        # --- CVCW: TRECW(I) = CRWDTH(I) (variant forest-grown crown width) ---
        cw   = _cover_crwdth(s, Int(t.species[i]), dbh, ht, icr, ba)

        ispi = spmap === nothing ? 1 : spmap[Int(t.species[i])]
        cr   = Float32(icr) / 100.0f0
        rad  = cw * 0.5f0

        # --- CVSHAP: crown shape ---
        ishap = _cover_cvshap(ispi, dbh, ht, cr, rad, tpa)

        # --- CVCBMS (COVOPT=2): foliage biomass ---
        cl_b = Float32(icr) * ht / 100.0f0
        rd   = dbh / rmsqd
        dds  = (2f0*dbh*dg + dg*dg) / fint
        dds < 0.0001f0 && (dds = 0.0001f0)
        trfbms = _cover_cvcbms2(ispi, dbh, ht, cl_b, rd, dds, tpa)

        # --- SHRUB-SMALL CONIFER COMPETITION tree tally (cvsum.f DO 300, LBROW): bin HT
        #     into the 11 SHTRHT thresholds (cumulated after the loop). Report-only. ---
        if cv.lbrow
            hm1 = -1.0f0
            @inbounds for ih in 1:11
                hc = _CVSUM_SHTRHT[ih]
                (ht > hm1 && ht <= hc) && (trsh[ih] += prob)
                hm1 = hc
            end
        end

        # --- CVSUM: bin this tree into 10-ft height classes (canopy path) ---
        _cover_cvsum_tree!(txht, crxht, proxht, volxht, cfbxht,
                           ishap, cw, dbh, ht, icr, prob, trfbms)
        sdiam += dbh*prob/12.0f0
    end

    # totals & overlap adjustment (CVSUM DO 350)
    tretot = 0.0f0; tpctcv = 0.0f0; tproar = 0.0f0; tcvol = 0.0f0; totbms = 0.0f0
    @inbounds for j in 1:16
        pcxht[j] = crxht[j]/435.6f0
        tpctcv += pcxht[j]
        pcxht[j] = 100.0f0*(1.0f0 - exp(-0.01f0*pcxht[j]))
        tretot += txht[j]
        tproar += proxht[j]
        tcvol  += volxht[j]
        totbms += cfbxht[j]
    end
    tpctcv = 100.0f0*(1.0f0 - exp(-0.01f0*tpctcv))
    icvage = iage + (Int(year) - Int(year0))

    push!(cv.rows, CoverRow(Int(year), txht, pcxht, proxht, volxht, cfbxht,
                            tretot, tpctcv, tproar, tcvol, totbms, avh, sdiam, icvage))

    # ---- SHRUB half (CVBROW→CVSUM shrub→CVCLAS); report-only -------------------------
    if cv.lbrow
        # cumulate trsh from the top (cvsum.f DO 399: TRSH[j] += TRSH[j+1]) so each column
        # J = trees/ac with HT greater than SHTRHT[J-1] (display heights 0.0,0.5,…,20.0).
        @inbounds for j in 10:-1:1
            trsh[j] += trsh[j+1]
        end
        _cover_shrub_cycle!(cv, s, Int(year), Int(year0), ba, avh, iage, rmsqd,
                            crxht, sd2xht, htmax, htmin, trsh)
    end
    return cv
end

# Drive one cycle of the shrub half.  Computes CVSCON constants at cycle 0, then the
# per-cycle CVBROW predictions, the CVSUM shrub aggregation, and the CVCLAS stage.
function _cover_shrub_cycle!(cv::CoverState, s::StandState, year::Int, year0::Int,
                             ba::Float32, avh::Float32, iage::Int, rmsqd::Float32,
                             crxht::Vector{Float32}, sd2xht::Vector{Float32},
                             htmax::Float32, htmin::Float32, trsh::Vector{Float32})
    idist = cv.idist
    # SAGE: time since disturbance.  Start value from card (−1 ⇒ stand age), floored at 3.
    sstart = cv.sage0 < 0.0f0 ? Float32(iage) : cv.sage0
    sstart < 3.0f0 && (sstart = 3.0f0)          # LSAGE clamp (subtract-3 quirk deferred)
    sage = sstart + Float32(year - year0)
    sage > 40.0f0 && return cv                   # shrubs not computed past 40 yr

    # LSTART = first shrub cycle (== cycle 0): CVSCON + CVBCAL calibration run once here.
    lstart = cv.shrub_const === nothing
    # CVSCON constants (once, at cycle 0 / first shrub cycle).
    if cv.shrub_const === nothing
        ihtype = cv.ihtype == 0 ? _cover_stand_habitat(s) : cv.ihtype
        iov, iun, itun, ok = _cvb_habitat(ihtype)
        ok || (cv.lbrow = false; return cv)      # invalid habitat → shrubs off
        cv.ihtype = ihtype
        inf = _cvb_inf(variant_code(s.variant), Int(s.plot.habitat_input) == 0 ? 1 : 1)
        cv.shrub_const = cvscon(s.plot.slope, s.plot.elevation, s.plot.aspect,
                                iov, iun, cv.iphys, inf, idist, itun)
    end
    sc = cv.shrub_const::ShrubConst

    cal = cv.calib.lcalib ? cv.calib : nothing
    cyc = cvbrow_cycle(sc, sage, ba, idist; cal=cal, lstart=lstart)
    # CVBCAL: compute correction factors at cycle 0, then apply each cycle (cvbrow DO 63/65).
    if cal !== nothing
        lstart && !cal.computed &&
            cvbcal!(cal, cyc.sh, cyc.pbcv, cyc.pb, cyc.htindx, cyc.totlcv)
        cyc = _cvbcal_apply(cal, cyc)
    end
    ss  = cvsum_shrub(cyc)
    relden = s.plot.relative_density
    crarea = sum(crxht)
    istage = cvclas(sage, relden, avh, ss.tallsh, cyc.totlcv,
                    crxht, sd2xht, _cover_txht(cv), crarea, rmsqd, htmax, htmin)

    # Irwin & Peek shrub biomass / twig production (GF/CLUN, WRC/CLUN, WH/CLUN only).
    sbmass = 0.0f0; twigs = 0.0f0
    ih = cv.ihtype
    if ih == 520 || ih == 530 || ih == 570
        igf = ih == 520 ? 1.0f0 : 0.0f0
        ice = ih == 530 ? 1.0f0 : 0.0f0
        sbmass = _f32exp(5.52f0 + 1.29f0*_f32log(sage) - 0.1f0*sage - 0.01f0*relden) * 0.8921791f0
        twigs  = _f32exp(2.67f0 + 0.98f0*_f32log(sage) - 0.08f0*sage - 0.09f0*relden +
                         0.69f0*igf + 0.52f0*ice) * 10.76387f0
    end

    push!(cv.shrub_rows, ShrubRow(year, sage, cyc.pgt0, ss.clow, ss.cmed, ss.ctall,
                                  cyc.totlcv, ss.asht, ss.tallsh, sbmass, twigs, istage,
                                  ss.issp, ss.scv, ss.sht, ss.spb, ss.scov, trsh))
    return cv
end

# TXHT for the current cycle = the tree-count-by-class of the row just pushed.
_cover_txht(cv::CoverState) = isempty(cv.rows) ? zeros(Float32,16) : cv.rows[end].txht

# =====================================================================================
# CVOUT — emit the "CANOPY COVER STATISTICS" table (canopy path; report-only).
# =====================================================================================
_cv_ifix(x) = trunc(Int, 0.5f0 + x)     # FVS IFIX(.5+x)

function _cover_stand_header(io, stand_id, mgmt_id, title)
    println(io)
    println(io, "STAND ID: ", rpad(stand_id, 26), "    MGMT CODE: ",
            rpad(mgmt_id, 4), "  ", title)
    println(io)
end

function cover_report(cv::CoverState, io::IO, stand_id::AbstractString,
                      mgmt_id::AbstractString, title::AbstractString)
    has_canopy = cv.lcover && cv.lcnop && !isempty(cv.rows)
    has_shrub  = cv.lshrub && cv.lbrow && !isempty(cv.shrub_rows)
    has_sum    = cv.lcvsum && (cv.lcnop || cv.lbrow) &&
                 (!isempty(cv.rows) || !isempty(cv.shrub_rows))
    (has_canopy || has_shrub || has_sum) || return io
    if has_canopy
    _cover_stand_header(io, stand_id, mgmt_id, title)
    println(io)
    println(io, "-"^51, "  CANOPY COVER STATISTICS  ", "-"^52)
    println(io, " "^52, "(BASED ON STOCKABLE AREA)")
    println(io)
    println(io, " "^48, "-"^31)
    println(io, " "^49, "ATTRIBUTE BY 10' HEIGHT CLASS")
    println(io, " "^48, "-"^31)
    println(io)
    println(io, " "^25, "  TREES -- TREES PER ACRE")
    println(io, " "^25, "  COVER -- PERCENTAGE OF CANOPY COVER CONTRIBUTED BY TREES IN HEIGHT CLASS")
    println(io, " "^25, "   AREA -- CROWN PROFILE AREA (SQ.FT. PER ACRE)")
    println(io, " "^25, " VOLUME -- CROWN VOLUME (CU.FT. PER ACRE X 100)")
    println(io, " "^25, "BIOMASS -- FOLIAGE BIOMASS (LBS. PER ACRE)")
    println(io); println(io); println(io)
    println(io, "-"^130)
    println(io, " "^53, "STAND HEIGHT CLASS")
    println(io, "            0.0-  10.1-  20.1-  30.1-  40.1-  50.1-  60.1-  ",
            "70.1-  80.1-  90.1- 100.1- 110.1- 120.1- 130.1- 140.1- 150.1+")
    println(io, "YEAR       10.0'  20.0'  30.0'  40.0'  50.0'  60.0'  70.0'  ",
            "80.0'  90.0' 100.0' 110.0' 120.0' 130.0' 140.0' 150.0' ",
            " "^7, "TOTAL")
    println(io, "-"^130)
    for row in cv.rows
        println(io, lpad(row.year, 4))
        _cover_emit_line(io, "    TREES", row.txht, row.tretot)
        _cover_emit_line(io, "    COVER", row.pcxht, row.tpctcv)
        _cover_emit_line(io, "     AREA", row.proxht, row.tproar)
        _cover_emit_line(io, "   VOLUME", row.volxht ./ 100.0f0, row.tcvol/100.0f0)
        _cover_emit_line(io, "  BIOMASS", row.cfbxht, row.totbms)
        println(io)
    end
    end  # has_canopy
    (has_shrub && cv.calib.lcalib && cv.calib.computed) &&
        _cover_calib_stats(cv, io, stand_id, mgmt_id, title)
    has_shrub && _cover_shrub_stats(cv, io, stand_id, mgmt_id, title)
    has_sum   && _cover_summary(cv, io, stand_id, mgmt_id, title)
    # SHRUB-SMALL CONIFER COMPETITION follows the summary on the same page (cvout.f 9090).
    (has_sum && cv.lbrow && !isempty(cv.shrub_rows)) && _cover_shrub_conifer(cv, io)
    return io
end

function _cover_emit_line(io::IO, label::AbstractString,
                          vals::AbstractVector{Float32}, tot::Float32)
    print(io, label)
    for j in 1:16
        print(io, lpad(_cv_ifix(vals[j]), 7))
    end
    print(io, lpad(_cv_ifix(tot), 8))
    println(io)
end


# ####################################################################################
# SHRUB HALF (staged port; validated vs FVSem_g16 DEBUG goldens — see scratchpad/cover)
# ####################################################################################

# ---- from cover_shrub.jl ----
# cover_shrub.jl — SHRUB half of the FVS COVER extension (vcovr/cvbrow.f driver →
# covr/cvscon.f shrub-cover model → covr/cvbcal.f calibration → covr/cvclas.f
# successional stage → the SHRUB portion of covr/cvsum.f aggregation).
#
# REPORT-ONLY.  No tree-record writes, no RNG.  Float32 throughout (Fortran REAL*4);
# transcendentals routed to glibc libm (logf/expf) for Float32 bit-exactness per doctrine.
#
# This file is a STAGED addition to src/engine/cover.jl.  It is written standalone so the
# replay harness (replay_shrub.jl) can validate it against the FVSem_g16 DEBUG goldens
# BEFORE it is merged into cover.jl.  The merge recipe is in HANDOFF.md.

# ---- Float32 libm transcendentals (1-ULP-faithful to Fortran REAL EXP/ALOG) -----------
@inline _f32log(x::Float32) = ccall((:logf, "libm.so.6"), Float32, (Float32,), x)
@inline _f32exp(x::Float32) = ccall((:expf, "libm.so.6"), Float32, (Float32,), x)

# =====================================================================================
# CVSCON coefficient tables (covr/cvscon.f).  Fortran DATA fills each (N,31) array
# column-major, so the flat DATA order is species-major: the first N values are
# species-1's N category coefficients.  reshape(vec, N, 31) therefore gives
# TAB[category, species] — index TAB[cat, isp].
# =====================================================================================

# PHABOV(5,31): overstory habitat series [DF, GF, WRC, WH, AF/MH]
const _CVS_PHABOV = reshape(Float32[
  9.088568,1.112479,0.082179,-6.242122,-4.041104, 0.097835,-0.976795,-1.443615,3.962374,-1.639797,
 -2.281724,-0.300165,0.222035,1.248187,1.111665, -1.559909,-0.066248,0.221251,1.464268,-0.059362,
  0.890738,0.560309,-0.510439,-0.700194,-0.240415, 0.000000,-0.648945,-1.715471,2.257159,0.107257,
 -0.049016,-0.660303,0.000000,0.000000,0.709319, -1.692400,0.206668,0.572810,0.336420,0.576504,
  0.000000,-0.194440,0.358941,-0.738026,0.573525, 1.210348,0.601579,-0.530107,-0.277865,-1.003957,
 -0.404204,0.101982,0.600343,0.049325,-0.347447, 0.119231,0.590724,0.157729,-0.681847,-0.185837,
 -0.991672,0.382548,0.766728,0.175734,-0.333337, 0.558700,5.350603,5.162087,-5.814926,-5.256464,
  1.247055,0.793011,0.220189,-1.791092,-0.469164, -1.702553,0.084306,-0.184811,0.671725,1.131332,
 -1.558622,-0.443300,-1.343065,0.972536,2.372452, -0.198937,-0.085908,0.533011,0.225364,-0.473530,
  0.594133,-0.773278,0.232943,-0.399616,0.345819, -0.647701,0.341590,0.140994,0.334487,-0.169370,
 -1.261175,0.259181,0.512034,-0.013122,0.503083, 0.780336,0.072750,-0.863347,-0.234732,0.244993,
 -0.453181,0.289715,1.287433,0.291517,-1.415481, 0.490208,1.054449,0.220875,0.523104,-2.288634,
 -6.071812,1.450901,2.011634,1.663452,0.945825, 0.247127,0.981953,0.557725,-0.342670,-1.444135,
  0.787814,-0.427235,-0.264309,0.381721,-0.477991, 1.386547,-0.155526,0.346276,-0.395677,-1.181620,
 -0.207683,0.099576,0.616547,0.239518,-0.747958, -1.092577,0.062480,0.373197,-0.055660,0.712560,
 -0.689655,0.289342,-0.797073,-0.017767,1.215152], 5, 31)

# PHABUN(5,31): understory union [CLUN1, TALL, LOW, CLUN2, GRASS]
const _CVS_PHABUN = reshape(Float32[
  2.710610,-6.075198,-4.848361,8.212950,0.000000, 1.043959,1.666373,1.724007,-5.630460,1.196120,
  0.413300,-0.557303,0.000000,0.144004,0.000000, 0.226423,-0.582335,0.000000,0.355912,0.000000,
  0.107452,-0.010810,1.084070,0.119038,-1.299750, 0.892790,1.386127,0.128510,-2.026148,-0.381281,
 -1.352566,-0.236322,0.523485,0.000000,1.065404, 0.636272,0.819635,-0.032862,0.992910,-2.415953,
 -1.609339,-0.226909,0.000000,-0.369997,2.206244, 0.673904,2.001927,0.074934,-0.883838,-1.866927,
  0.117480,-0.055445,-0.400352,0.421247,-0.082930, -0.251528,-0.006640,0.059342,0.260100,-0.061274,
  0.864850,-0.439995,-1.869472,1.444616,0.000000, -4.353611,-0.281722,-0.174776,4.810109,0.000000,
 -0.092514,0.035195,0.718816,0.023768,-0.685265, 0.855702,0.399076,-1.171601,0.586127,-0.669304,
  2.038632,0.568528,-2.315288,0.697174,-0.989048, 0.696647,-0.214418,-1.248584,0.766356,0.000000,
 -0.186215,-0.048786,0.171811,-0.465541,0.528731, 0.536920,0.522479,-0.470882,0.321000,-0.909518,
  0.615945,0.749974,-1.848295,1.039323,-0.556947, 0.874847,0.185377,-0.014645,-0.053301,-0.992279,
 -0.192900,0.160415,0.586891,-0.554406,0.000000, -0.280322,-0.253554,0.098498,0.454597,-0.019220,
  0.000000,0.000000,0.000000,0.000000,0.000000, -1.034255,1.088705,0.245059,-0.876263,0.576753,
  0.000000,0.000000,0.000000,0.000000,0.000000, 0.325657,-1.100824,-0.022447,0.932375,-0.134761,
  0.066676,0.103292,0.081037,0.230391,-0.481397, 0.472095,-0.117253,-0.228589,0.307475,-0.433726,
  0.925955,0.452005,-1.248314,0.659459,-0.789106], 5, 31)

# PPHYS(5,31): physiography [BOTTOM, LOWER, MID, UPPER, RIDGE]
const _CVS_PPHYS = reshape(Float32[
 -1.115882,0.237709,0.397100,-0.029464,0.510538, -0.758886,-0.158827,-0.084437,-0.298196,1.300344,
  0.472246,0.026784,-0.037302,0.132755,-0.594484, -0.330405,0.182488,0.079716,0.357364,-0.289164,
  0.154180,-0.136077,0.109557,-0.109143,-0.018517, 0.000000,-0.121054,-0.200353,0.178785,0.142623,
 -0.291258,0.239395,0.076529,0.185398,-0.210063, -0.189881,0.085693,0.216330,0.267704,-0.379845,
 -0.025618,0.338617,0.651530,0.064811,-1.029338, -1.221598,-0.192430,0.105655,0.604654,0.703721,
 -0.464616,0.029192,0.281638,0.335585,-0.181799, -0.166037,-0.111284,0.309331,0.295109,-0.327118,
  0.507743,0.180563,0.133105,-0.037654,-0.783756, -0.785894,-0.758659,-0.603993,0.285434,1.863111,
 -0.070391,0.165321,0.045400,0.190815,-0.331146, -0.525546,0.089144,0.251261,0.399115,-0.213974,
 -0.145504,0.126120,-0.383385,0.246056,0.156713, 0.558692,0.542812,0.311799,-0.621233,-0.792070,
 -0.463954,0.159999,-0.141970,0.375748,0.070176, 0.247838,0.164348,0.082047,-0.222117,-0.272116,
  0.801413,-0.039444,-0.056379,-0.378473,-0.327117, 0.096562,0.052876,0.133738,-0.059111,-0.224065,
 -0.980035,-0.786041,0.137313,0.447255,1.181507, 0.372962,0.178652,-0.060715,-0.392138,-0.098762,
  0.270737,0.335691,-0.281113,0.194329,-0.519645, -1.438125,0.389894,0.203401,0.163025,0.681806,
 -0.878815,0.338468,0.002664,-0.110483,0.648166, 0.316784,0.097615,-0.272526,-0.140761,-0.001112,
  0.240351,-0.496234,-0.210591,-0.029912,0.496386, 0.036801,-0.009421,-0.478119,-0.399715,0.850454,
 -0.311908,-0.021160,0.028282,0.118770,0.186016], 5, 31)

# PFLOC(4,31): national forest location [SOUTH, NEZPERCE, CENTRAL, NORTH]
const _CVS_PFLOC = reshape(Float32[
  3.600125,-1.499356,-2.554641,0.453873, 0.812938,-0.914264,-0.526915,0.628240,
  0.131632,0.449116,-0.208468,-0.372281, 0.458707,-1.408397,0.403550,0.546141,
  0.460303,-0.398612,-0.059656,-0.002034, 0.682410,-1.902216,0.451361,0.768446,
  0.000000,0.000000,0.000000,0.000000, 0.713447,-0.608892,0.107954,-0.212510,
 -3.879320,0.636311,1.681632,1.561377, 0.123896,0.349324,-0.461448,-0.011772,
 -0.199094,0.561158,-0.093965,-0.268099, 0.025643,0.080601,-0.130028,0.023784,
 -0.225770,0.390307,-0.038916,-0.125622, -2.749750,0.000000,-0.621770,3.371520,
  0.454416,0.867404,-0.068041,-1.253778, -0.327891,-0.368999,0.475865,0.221025,
 -0.656554,1.036110,0.009707,-0.389264, 0.333801,-0.292903,0.169438,-0.210336,
  0.483550,-0.081656,-0.081960,-0.319933, -0.063613,-0.162506,0.073234,0.152886,
 -0.667632,-0.649358,0.299813,1.017177, 0.549982,-0.024194,-0.185210,-0.340578,
 -1.713942,0.890978,0.389472,0.433492, 0.942476,-1.046767,0.600527,-0.496235,
 -0.084171,-0.108050,0.086456,0.105765, -1.034402,1.177656,-0.231173,0.087919,
  2.427337,-0.389151,-0.526386,-1.511799, 1.514727,0.208822,-0.338610,-1.384938,
  0.554700,0.230839,-0.473056,-0.312484, -0.350419,0.620746,-0.012547,-0.257780,
  0.638508,-1.180981,0.142076,0.400397], 4, 31)

# PDIST(4,31): disturbance [NONE, MECH, BURN, ROAD]
const _CVS_PDIST = reshape(Float32[
  0.702636,-1.154640,0.452005,0.0, 0.724008,-0.392686,-0.331322,0.0,
  0.013091,-0.495381,0.482290,0.0, 0.312114,-0.471186,0.159072,0.0,
 -0.018909,0.302284,-0.283375,0.0, -0.273111,0.297179,-0.024068,0.0,
  0.173677,0.085770,-0.259447,0.0, 0.217620,0.106280,-0.323900,0.0,
  0.350763,-0.348382,-0.002381,0.0, 0.432237,0.054443,-0.486681,0.0,
 -0.553423,0.418804,0.134619,0.0, 0.276317,-0.015560,-0.260756,0.0,
 -0.162309,0.306713,-0.144405,0.0, 0.868899,-2.708511,1.839612,0.0,
  0.170876,0.460758,-0.631634,0.0, 0.331513,-0.237219,-0.094294,0.0,
  0.067213,-0.376034,0.308821,0.0, 0.501561,0.161611,-0.663172,0.0,
 -0.647399,0.395331,0.252068,0.0, 0.302101,0.224373,-0.526473,0.0,
 -0.509263,0.222324,0.286939,0.0, 0.376367,-0.110195,-0.266171,0.0,
 -0.424833,-0.571652,0.996485,0.0, -0.483495,0.293256,0.190239,0.0,
  0.686777,1.240241,-1.927018,0.0, -0.023516,-0.488626,0.512142,0.0,
 -0.303771,0.335156,-0.031385,0.0, -0.215707,0.293401,-0.077694,0.0,
 -0.279459,-0.316944,0.596403,0.0, -0.210975,0.038058,0.172917,0.0,
  0.576745,-0.015996,-0.560749,0.0], 4, 31)

# POTHER(6,31): site [INTERCEPT, SLOPE, ELEV, ELEV*ELEV, SXCOS, SXSIN]
const _CVS_POTHER = reshape(Float32[
 -6.234583,-5.926120,0.199323,-0.005927,-0.891665,1.722985, 3.877884,-0.109995,-0.360401,0.003286,-3.621449,0.016028,
 -7.762071,-0.102417,0.255819,-0.004110,1.429693,0.073954, -14.085910,0.607160,0.556593,-0.006863,-0.684819,0.216532,
 -3.715241,1.825835,0.006793,-0.000081,-0.715995,-0.070842, -22.552940,-3.503802,0.689370,-0.006252,-0.835270,0.986894,
 -20.695860,1.276123,0.594190,-0.004380,-0.928040,1.255746, -8.949919,0.616951,0.215076,-0.002198,0.481796,0.441640,
 -16.076040,-2.297629,0.402144,-0.003145,4.851739,0.256853, -3.758699,2.789399,-0.047727,0.000038,-0.295032,0.113204,
 -9.333570,0.065296,0.249074,-0.002068,0.764599,0.114451, -5.426313,0.387720,0.123033,-0.001946,-0.629058,-0.329992,
 -6.268844,2.151100,0.103793,-0.001099,0.542544,0.633462, -0.037685,-4.960156,-0.232049,0.001834,0.612515,-2.289439,
  4.429537,-0.764605,-0.233709,0.001609,-1.650102,-0.442272, -11.674530,0.327773,0.311276,-0.002692,0.538878,-0.088429,
 -20.323050,2.743548,0.469166,-0.003569,-1.730066,-0.102061, -5.463860,0.765169,0.074884,-0.000837,0.163987,-0.040309,
  3.540272,-0.589173,-0.296972,0.002683,0.030100,-0.019474, -8.650060,2.651302,0.161251,-0.001873,0.267403,0.375045,
 -9.934808,-0.615187,0.162745,-0.000919,2.224786,0.678724, -5.050717,1.425095,0.097565,-0.001720,-0.898583,0.151992,
  0.045530,0.764632,-0.154438,0.000545,-1.677844,-2.045689, -11.839120,-0.505556,0.356478,-0.003468,-1.537973,0.596131,
 -9.238182,-0.336484,0.221618,-0.003561,1.363716,1.150854, 4.458409,0.634152,-0.319046,0.002291,-0.636745,0.488674,
 -12.080660,0.121266,0.206421,-0.001605,-1.438883,0.383933, -8.600142,1.689407,0.177508,-0.002386,-2.394606,-1.039733,
 -3.127204,-1.266641,0.049740,-0.000591,0.125282,0.043020, -14.135830,0.401931,0.358190,-0.002911,1.170672,1.220915,
 -18.351950,1.030739,0.467154,-0.004475,0.025292,0.358210], 6, 31)

# TOTAL SHRUB COVER habitat/disturbance constants
const _CVS_THABU1 = Float32[0.65403,0.267,-0.00503,0.454,-0.180,-1.19]      # by ITUN
const _CVS_THABU2 = Float32[0.0,-0.031625,-0.091442,-0.27636,-0.002983,-0.771579]
const _CVS_TDIST  = Float32[0.420077,-0.138963,0.112421,0.0]               # by IDIST

# ---- HEIGHT coefficients (per-species variable-length arrays) -------------------------
const _CVS_IHI = Int[0,1,0,10,10,0,0,9,0,0,1,5,5,0,10,10,0,1,0,0,1,1,5,5,0,0,0,0,1,1,0]
const _CVS_IPI = Int[0,0,0,1,5,0,0,0,0,1,0,0,10,0,5,5,1,6,0,5,10,0,0,14,5,0,0,0,6,0,0]
const _CVS_ILI = Int[0,0,0,6,1,0,0,5,1,6,0,0,15,1,1,1,0,0,1,0,6,6,10,10,0,5,0,1,11,0,1]
const _CVS_IDI = Int[0,0,0,0,0,0,0,1,0,0,6,1,1,0,0,15,0,0,0,1,0,0,1,1,1,1,0,0,0,0,0]

const _CVS_HARUV = Float32[0.500000]
const _CVS_HBERB = Float32[0.635000,0.024375,0.100294,0.000000,0.127500,-0.235000]
const _CVS_HLIBO = Float32[0.500000]
const _CVS_HPAMY = Float32[-2.192232,-0.036425,0.225308,0.075128,0.000000,0.161033,-0.208727,-0.651338,-0.044806,0.000000,0.146765,-0.128785,0.000000,0.007612,-0.069397,0.617216]
const _CVS_HSPBE = Float32[0.087057,0.000000,-0.034163,0.052158,0.000000,0.016369,0.128729,0.109886,0.000000,0.092615,-0.113355,-0.027844,0.000000,0.138574,-0.086942]
const _CVS_HVASC = Float32[-1.584374,1.215803,0.460397,1.560992]
const _CVS_HCARX = Float32[0.500000]
const _CVS_HLONI = Float32[0.447867,0.044723,-0.014173,0.003483,0.000000,-0.054145,-0.288548,-0.059246,0.000000,-0.059526,0.037735,0.000000,0.095682,0.053252,0.175588,0.172200,0.139718]
const _CVS_HMEFE = Float32[2.702792,-0.031333,-0.850066,-0.173478,0.000000,-0.016392,0.116897]
const _CVS_HPHMA = Float32[7.113002,-1.834477,0.170875,-0.030095,0.000000,-0.370640,0.122902,0.590908,0.079706,0.000000,-0.505541,-1.169115]
const _CVS_HRIBE = Float32[0.922317,0.185217,-0.069884,0.000000,-0.017954,-0.213053,-0.073536,-0.128280,-0.121669,0.000000,0.051824,-0.100895,0.065351]
const _CVS_HROSA = Float32[0.751472,-0.023310,-0.114210,-0.128391,0.000000,-0.174293,-0.118963,0.000000,0.111218,-0.165629,0.108633,-0.000216,0.264153]
const _CVS_HRUPA = Float32[0.902969,-0.064163,-0.107400,-0.200936,0.000000,-0.323602,-0.177336,0.000000,-0.170269,-0.160183,0.204840,0.166577,0.122142,0.000000,0.154868,0.255838,0.169348,-0.005363,0.000000,-0.112133]
const _CVS_HSHCA = Float32[0.816644,-0.273829,0.000000,0.000000,0.000000,-0.177801,-0.450658,0.300302]
const _CVS_HSYMP = Float32[-0.431044,-0.026261,-0.150352,-0.261149,0.000000,0.301546,0.135749,0.062181,0.000000,0.003662,0.045338,0.083974,0.000000,0.097050,0.271221,0.011762]
const _CVS_HVAME = Float32[0.643919,-0.274397,-0.217928,-0.077408,0.000000,-0.049784,0.094958,0.083723,0.000000,0.060136,-0.352248,-0.148628,0.000000,-0.067128,-0.159786,-0.118378,-0.378361,-0.653911,0.000000,0.220559,0.097746,0.103688]
const _CVS_HXETE = Float32[3.022810,0.210587,-0.140420,-0.083689,0.000000,0.137891,-0.687047,0.286603]
const _CVS_HFERN = Float32[0.479118,-0.257565,0.108291,0.000000,-0.072578,-0.215387,0.416443,0.484604,0.245805,0.000000,0.305190]
const _CVS_HCOMB = Float32[-2.142289,-0.212620,-0.135207,-0.490167,0.000000,0.525595,0.962807,0.851434,0.662427]
const _CVS_HACGL = Float32[2.606310,-0.030179,0.096289,-0.130003,0.000000,0.118123,-0.089183,0.145533,0.000000,-0.195650,-0.339980,-0.223145]
const _CVS_HALSI = Float32[0.130874,0.096784,0.104609,0.000000,-0.050438,0.065727,-0.061039,-0.385028,-0.102994,0.000000,0.556503,0.203410,0.225201,0.000000,0.001873]
const _CVS_HAMAL = Float32[0.191243,-0.155838,-0.059025,0.000000,0.102939,-0.178142,0.385587,0.366095,0.148307,0.000000,-0.530085]
const _CVS_HCESA = Float32[1.721365,-0.081934,-0.260050,-0.044915,0.000000,-0.269311,-0.006042,0.000000,0.182349,0.019450,0.203430,0.010725,0.111365,0.000000,-0.011168,0.571874,0.162710,-0.246587]
const _CVS_HCEVE = Float32[1.482635,0.545563,0.172759,0.297625,0.000000,-0.634357,-0.093287,0.000000,0.530326,-0.053496,0.628226,-0.290862,0.806147,0.000000,0.054440,0.210558,0.455416,0.000000,0.127645]
const _CVS_HCOST = Float32[1.306008,0.005804,0.048403,-0.213566,0.000000,0.178659,0.131490,0.190049,0.000000,0.000000]
const _CVS_HHODI = Float32[2.860934,0.123307,0.020195,-0.009312,0.000000,0.119827,0.192505,0.052546,0.000000,-0.647201]
const _CVS_HPREM = Float32[1.097208]
const _CVS_HPRVI = Float32[3.260277,-0.469992,-0.130373,-0.173093,0.000000,-0.024267,0.350357,0.311887,0.079295]
const _CVS_HSALX = Float32[0.980963,-0.135857,-0.162329,0.000000,-0.000261,-0.390787,0.062137,-0.153430,-0.024562,0.000000,0.118595,0.290822,0.054484,-0.023476,0.000000,0.240669,-0.254219,-0.235183]
const _CVS_HSAMB = Float32[1.547656,-0.062514,-0.197323,0.000000,-0.281832,-0.310178,-0.061625,0.127680,-0.010855,0.899952]
const _CVS_HSORB = Float32[1.348299,-0.147872,-0.420400,-0.188858,0.000000]

# ---- PERCENT COVER coefficients -------------------------------------------------------
const _CVS_ICH = Int[0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,0,0,1,0,0,1,5,5,5,0,0,1,1,0]
const _CVS_ICP = Int[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0]
const _CVS_ICL = Int[0,1,0,0,5,0,0,6,0,10,6,6,0,0,0,10,0,0,0,0,0,0,6,10,0,0,0,0,0,0,0]
const _CVS_ICD = Int[0,0,0,0,1,0,0,0,6,6,0,0,6,0,6,6,0,0,0,6,0,0,10,1,1,1,0,0,6,0,0]

const _CVS_CARUV = Float32[7.410459,-1.277750]
const _CVS_CBERB = Float32[3.814140,-0.159095,0.267293,0.126369,0.000000,-0.435941,-0.442553]
const _CVS_CLIBO = Float32[2.905366,-1.440789,-0.326447,0.000000,0.095137,-0.383999]
const _CVS_CPAMY = Float32[2.018329,-0.354632,0.154106,0.000000,0.203464,-0.286453,-0.717539]
const _CVS_CSPBE = Float32[3.521268,0.009649,-0.074220,-0.034684,0.000000,0.517227,-0.069129,-0.101223,0.000000,-0.070098,-0.166171,0.124020]
const _CVS_CVASC = Float32[2.669771]
const _CVS_CCARX = Float32[2.735197,0.187976,-0.249295,0.000000,0.000000,0.000000,-0.473890]
const _CVS_CLONI = Float32[0.956830,0.157181,-0.054025,0.000000,-0.076409,-0.259281,-0.095843,0.172093,0.161978,0.000000,-0.766977]
const _CVS_CMEFE = Float32[7.616001,0.000000,-0.573019,0.000000,-0.090717,-0.212948,-0.104867,-0.228262,0.004891,0.000000,-0.874342,0.126315,-0.022442,1.295553]
const _CVS_CPHMA = Float32[14.537630,-0.110563,0.061265,0.000000,-0.022906,0.472747,-0.043693,-0.262411,-0.241079,0.000000,0.153486,0.896645,0.168982,0.000000,-1.879735,-0.744306]
const _CVS_CRIBE = Float32[-0.718346,-0.401155,0.075869,0.000000,0.023182,0.163116,-0.323542,-0.029047,0.002164,0.000000,-0.574932,0.006134]
const _CVS_CROSA = Float32[1.346910,0.089665,0.146983,0.000000,-0.119553,-0.071583,-0.387973,-0.235797,0.002709,0.000000,-0.289654,-0.017397,-0.752834]
const _CVS_CRUPA = Float32[-0.066672,0.070099,0.124684,0.000000,0.329402,0.150745,0.077226,0.181005,0.235033,0.000000,-0.191659]
const _CVS_CSHCA = Float32[-12.254140]
const _CVS_CSYMP = Float32[4.432121,-0.223466,-0.228996,0.000000,-0.589140,-0.228178,0.020206,-0.109190,-0.024860,0.000000,-0.265650,-0.135566,-0.034793,-0.411302]
const _CVS_CVAME = Float32[4.941563,-1.709876,-0.665779,0.000000,-0.302910,-0.646100,0.204189,-0.313452,-0.514592,0.000000,-1.274516,-0.737516,-0.337436,0.000000,0.862965,0.830979,0.532588]
const _CVS_CXETE = Float32[0.573361,0.101778,0.000000,0.000000,0.000000,0.270742,0.046568,-0.157317,-0.505288,0.010342]
const _CVS_CFERN = Float32[3.192776,-1.153183,0.379981,0.019392]
const _CVS_CCOMB = Float32[2.522206,0.054718,-0.505550,-0.586245]
const _CVS_CACGL = Float32[5.452397,0.164089,0.107552,0.000000,-0.231238,0.127618,-0.019414,0.297514,-0.001329,0.000000]
const _CVS_CALSI = Float32[6.588399,1.044745]
const _CVS_CAMAL = Float32[3.053453,0.015746,0.195635,0.271514,0.167479]
const _CVS_CCESA = Float32[11.957993,-0.941885,-0.780114,-0.730843,-0.354437,0.000000,0.000000,0.007548,-0.131303,-0.521329,-0.507620,-1.100015,-0.544426,0.000000,1.298341,0.813414,-0.602230,-1.422588]
const _CVS_CCEVE = Float32[5.657328,0.773755,0.311426,0.141040,0.000000,-0.470481,-0.057275,0.000000,0.483218,-0.262381,0.420124,-0.480422,0.681898,0.000000,0.842243]
const _CVS_CCOST = Float32[1.842154,0.701067,0.661132,-0.391999,0.000000,0.000000,0.743327,0.000000,0.564952,-1.138883,-0.868155,-1.125811,0.764758]
const _CVS_CHODI = Float32[0.488701,-0.229604,-0.218397,-0.221773,0.000000,0.150929,-0.097314,0.000000,-0.108503,-0.655863,-0.433806,0.223717,-0.234565]
const _CVS_CPREM = Float32[-3.412103,0.053313]
const _CVS_CPRVI = Float32[1.581330]
const _CVS_CSALX = Float32[3.484206,-0.433030,-0.084864,0.000000,0.038140,-0.021767,0.319310,0.346655,0.058027,0.000000,-0.356754,-0.160517,-0.036213,0.000000,-0.102810,0.226843]
const _CVS_CSAMB = Float32[1.372619,-0.916968,0.064605,0.000000,0.286853,0.039102,-1.040168]
const _CVS_CSORB = Float32[-0.731129,-0.729834]

# =====================================================================================
# CVSCON — compute PCON(31), TCON(2), HCON(31), CCON(31) from site + habitat indices.
# All inputs Float32 (SLOPE fraction, ELEV in 100's ft, ASPECT radians).
# =====================================================================================
struct ShrubConst
    pcon::Vector{Float32}
    tcon::NTuple{2,Float32}
    hcon::Vector{Float32}
    ccon::Vector{Float32}
end

function cvscon(slope::Float32, elev::Float32, aspect::Float32,
                iov::Int, iun::Int, iphys::Int, inf::Int, idist::Int, itun::Int)
    slopec = slope - 0.266f0
    sxcos  = slope * cos(aspect)
    sxsin  = slope * sin(aspect)
    alnel  = _f32log(elev)

    pcon = Vector{Float32}(undef, 31)
    @inbounds for i in 1:31
        pscal = _CVS_POTHER[1,i] + _CVS_POTHER[2,i]*slope + _CVS_POTHER[3,i]*elev +
                _CVS_POTHER[4,i]*elev*elev + _CVS_POTHER[5,i]*sxcos + _CVS_POTHER[6,i]*sxsin
        pcon[i] = _CVS_PHABOV[iov,i] + _CVS_PHABUN[iun,i] + _CVS_PPHYS[iphys,i] +
                  _CVS_PFLOC[inf,i] + _CVS_PDIST[idist,i] + pscal
    end
    tcon1 = -0.513f0 + 1.35f0*slope + 0.0957f0*elev - 0.00124f0*elev*elev + _CVS_THABU1[itun]
    tcon2 = 3.747585f0 + 0.852539f0*slopec - 0.009437f0*elev -
            0.975758f0*slopec*slopec + _CVS_TDIST[idist] + _CVS_THABU2[itun]

    # Height/cover index functions (per-species offset + category index)
    ih = Vector{Int}(undef,31); ip = Vector{Int}(undef,31)
    il = Vector{Int}(undef,31); id = Vector{Int}(undef,31)
    @inbounds for i in 1:31
        ih[i] = iov   + _CVS_IHI[i]
        ip[i] = iphys + _CVS_IPI[i]
        il[i] = inf   + _CVS_ILI[i]
        id[i] = idist + _CVS_IDI[i]
    end
    H = _CVS_HARUV; # aliases for readability handled inline below
    hcon = Vector{Float32}(undef, 31)
    hcon[1]  = _CVS_HARUV[1]
    hcon[2]  = _CVS_HBERB[1]+_CVS_HBERB[ih[2]]
    hcon[3]  = _CVS_HLIBO[1]
    hcon[4]  = _CVS_HPAMY[1]+_CVS_HPAMY[ip[4]]+_CVS_HPAMY[il[4]]+_CVS_HPAMY[ih[4]]+_CVS_HPAMY[16]*alnel
    hcon[5]  = _CVS_HSPBE[1]+_CVS_HSPBE[il[5]]+_CVS_HSPBE[ip[5]]+_CVS_HSPBE[ih[5]]
    hcon[6]  = _CVS_HVASC[1]+_CVS_HVASC[2]*slope+_CVS_HVASC[3]*sxcos+_CVS_HVASC[4]*sxsin
    hcon[7]  = _CVS_HCARX[1]
    hcon[8]  = _CVS_HLONI[1]+_CVS_HLONI[id[8]]+_CVS_HLONI[il[8]]+_CVS_HLONI[ih[8]]+_CVS_HLONI[15]*slope+_CVS_HLONI[16]*sxcos+_CVS_HLONI[17]*sxsin
    hcon[9]  = _CVS_HMEFE[1]+_CVS_HMEFE[il[9]]+_CVS_HMEFE[6]*elev+_CVS_HMEFE[7]*slope
    hcon[10] = _CVS_HPHMA[1]+_CVS_HPHMA[ip[10]]+_CVS_HPHMA[il[10]]+_CVS_HPHMA[11]*slope+_CVS_HPHMA[12]*alnel
    hcon[11] = _CVS_HRIBE[1]+_CVS_HRIBE[ih[11]]+_CVS_HRIBE[id[11]]+_CVS_HRIBE[11]*sxcos+_CVS_HRIBE[12]*sxsin+_CVS_HRIBE[13]*slope
    hcon[12] = _CVS_HROSA[1]+_CVS_HROSA[id[12]]+_CVS_HROSA[ih[12]]+_CVS_HROSA[11]*sxcos+_CVS_HROSA[12]*sxsin+_CVS_HROSA[13]*slope
    hcon[13] = _CVS_HRUPA[1]+_CVS_HRUPA[id[13]]+_CVS_HRUPA[ih[13]]+_CVS_HRUPA[ip[13]]+_CVS_HRUPA[il[13]]+_CVS_HRUPA[20]*alnel
    hcon[14] = _CVS_HSHCA[1]+_CVS_HSHCA[il[14]]+_CVS_HSHCA[6]*sxcos+_CVS_HSHCA[7]*sxsin+_CVS_HSHCA[8]*slope
    hcon[15] = _CVS_HSYMP[1]+_CVS_HSYMP[il[15]]+_CVS_HSYMP[ip[15]]+_CVS_HSYMP[ih[15]]+_CVS_HSYMP[16]*elev
    hcon[16] = _CVS_HVAME[1]+_CVS_HVAME[il[16]]+_CVS_HVAME[ip[16]]+_CVS_HVAME[ih[16]]+_CVS_HVAME[id[16]]+_CVS_HVAME[20]*sxcos+_CVS_HVAME[21]*sxsin+_CVS_HVAME[22]*slope
    hcon[17] = _CVS_HXETE[1]+_CVS_HXETE[ip[17]]+_CVS_HXETE[7]*alnel+_CVS_HXETE[8]*slope
    hcon[18] = _CVS_HFERN[1]+_CVS_HFERN[ih[18]]+_CVS_HFERN[ip[18]]
    hcon[19] = _CVS_HCOMB[1]+_CVS_HCOMB[il[19]]+_CVS_HCOMB[6]*alnel+_CVS_HCOMB[7]*slope+_CVS_HCOMB[8]*sxcos+_CVS_HCOMB[9]*sxsin
    hcon[20] = _CVS_HACGL[1]+_CVS_HACGL[id[20]]+_CVS_HACGL[ip[20]]+_CVS_HACGL[11]*alnel+_CVS_HACGL[12]*slope
    hcon[21] = _CVS_HALSI[1]+_CVS_HALSI[ih[21]]+_CVS_HALSI[il[21]]+_CVS_HALSI[ip[21]]
    hcon[22] = _CVS_HAMAL[1]+_CVS_HAMAL[ih[22]]+_CVS_HAMAL[il[22]]+_CVS_HAMAL[11]*slope
    hcon[23] = _CVS_HCESA[1]+_CVS_HCESA[id[23]]+_CVS_HCESA[ih[23]]+_CVS_HCESA[il[23]]+_CVS_HCESA[15]*elev+_CVS_HCESA[16]*slope+_CVS_HCESA[17]*sxcos+_CVS_HCESA[18]*sxsin
    hcon[24] = _CVS_HCEVE[1]+_CVS_HCEVE[id[24]]+_CVS_HCEVE[ih[24]]+_CVS_HCEVE[il[24]]+_CVS_HCEVE[ip[24]]
    hcon[25] = _CVS_HCOST[1]+_CVS_HCOST[id[25]]+_CVS_HCOST[ip[25]]
    hcon[26] = _CVS_HHODI[1]+_CVS_HHODI[id[26]]+_CVS_HHODI[il[26]]+_CVS_HHODI[10]*alnel
    hcon[27] = _CVS_HPREM[1]
    hcon[28] = _CVS_HPRVI[1]+_CVS_HPRVI[il[28]]+_CVS_HPRVI[6]*elev+_CVS_HPRVI[7]*sxcos+_CVS_HPRVI[8]*sxsin+_CVS_HPRVI[9]*slope
    hcon[29] = _CVS_HSALX[1]+_CVS_HSALX[ih[29]]+_CVS_HSALX[ip[29]]+_CVS_HSALX[il[29]]+_CVS_HSALX[16]*slope+_CVS_HSALX[17]*sxcos+_CVS_HSALX[18]*sxsin
    hcon[30] = _CVS_HSAMB[1]+_CVS_HSAMB[ih[30]]+_CVS_HSAMB[7]*sxcos+_CVS_HSAMB[8]*sxsin+_CVS_HSAMB[9]*elev+_CVS_HSAMB[10]*slope
    hcon[31] = _CVS_HSORB[1]+_CVS_HSORB[il[31]]

    ihc = Vector{Int}(undef,31); ipc = Vector{Int}(undef,31)
    ilc = Vector{Int}(undef,31); idc = Vector{Int}(undef,31)
    @inbounds for i in 1:31
        ihc[i] = iov   + _CVS_ICH[i]
        ipc[i] = iphys + _CVS_ICP[i]
        ilc[i] = inf   + _CVS_ICL[i]
        idc[i] = idist + _CVS_ICD[i]
    end
    ccon = Vector{Float32}(undef, 31)
    ccon[1]  = _CVS_CARUV[1]+_CVS_CARUV[2]*alnel
    ccon[2]  = _CVS_CBERB[1]+_CVS_CBERB[ilc[2]]+_CVS_CBERB[6]*alnel+_CVS_CBERB[7]*slope
    ccon[3]  = _CVS_CLIBO[1]+_CVS_CLIBO[ihc[3]]
    ccon[4]  = _CVS_CPAMY[1]+_CVS_CPAMY[ihc[4]]+_CVS_CPAMY[7]*slope
    ccon[5]  = _CVS_CSPBE[1]+_CVS_CSPBE[idc[5]]+_CVS_CSPBE[ilc[5]]+_CVS_CSPBE[10]*sxcos+_CVS_CSPBE[11]*sxsin+_CVS_CSPBE[12]*slope
    ccon[6]  = _CVS_CVASC[1]
    ccon[7]  = _CVS_CCARX[1]+_CVS_CCARX[ihc[7]]+_CVS_CCARX[7]*slope
    ccon[8]  = _CVS_CLONI[1]+_CVS_CLONI[ihc[8]]+_CVS_CLONI[ilc[8]]+_CVS_CLONI[11]*slope
    ccon[9]  = _CVS_CMEFE[1]+_CVS_CMEFE[ihc[9]]+_CVS_CMEFE[idc[9]]+_CVS_CMEFE[11]*sxcos+_CVS_CMEFE[12]*sxsin+_CVS_CMEFE[13]*elev+_CVS_CMEFE[14]*slope
    ccon[10] = _CVS_CPHMA[1]+_CVS_CPHMA[ihc[10]]+_CVS_CPHMA[idc[10]]+_CVS_CPHMA[ilc[10]]+_CVS_CPHMA[15]*alnel+_CVS_CPHMA[16]*slope
    ccon[11] = _CVS_CRIBE[1]+_CVS_CRIBE[ihc[11]]+_CVS_CRIBE[ilc[11]]+_CVS_CRIBE[11]*slope+_CVS_CRIBE[12]*elev
    ccon[12] = _CVS_CROSA[1]+_CVS_CROSA[ihc[12]]+_CVS_CROSA[ilc[12]]+_CVS_CROSA[11]*sxcos+_CVS_CROSA[12]*sxsin+_CVS_CROSA[13]*slope
    ccon[13] = _CVS_CRUPA[1]+_CVS_CRUPA[ihc[13]]+_CVS_CRUPA[idc[13]]+_CVS_CRUPA[11]*slope
    ccon[14] = _CVS_CSHCA[1]
    ccon[15] = _CVS_CSYMP[1]+_CVS_CSYMP[ihc[15]]+_CVS_CSYMP[idc[15]]+_CVS_CSYMP[11]*sxcos+_CVS_CSYMP[12]*sxsin+_CVS_CSYMP[13]*elev+_CVS_CSYMP[14]*slope
    ccon[16] = _CVS_CVAME[1]+_CVS_CVAME[ihc[16]]+_CVS_CVAME[idc[16]]+_CVS_CVAME[ilc[16]]+_CVS_CVAME[15]*slope+_CVS_CVAME[16]*sxcos+_CVS_CVAME[17]*sxsin
    ccon[17] = _CVS_CXETE[1]+_CVS_CXETE[ihc[17]]+_CVS_CXETE[7]*elev+_CVS_CXETE[8]*slope+_CVS_CXETE[9]*sxcos+_CVS_CXETE[10]*sxsin
    ccon[18] = _CVS_CFERN[1]+_CVS_CFERN[2]*slope+_CVS_CFERN[3]*sxcos+_CVS_CFERN[4]*sxsin
    ccon[19] = _CVS_CCOMB[1]+_CVS_CCOMB[2]*sxcos+_CVS_CCOMB[3]*sxsin+_CVS_CCOMB[4]*slope
    ccon[20] = _CVS_CACGL[1]+_CVS_CACGL[ihc[20]]+_CVS_CACGL[idc[20]]
    ccon[21] = _CVS_CALSI[1]+_CVS_CALSI[2]*slope
    ccon[22] = _CVS_CAMAL[1]+_CVS_CAMAL[2]*elev+_CVS_CAMAL[3]*slope+_CVS_CAMAL[4]*sxcos+_CVS_CAMAL[5]*sxsin
    ccon[23] = _CVS_CCESA[1]+_CVS_CCESA[ihc[23]]+_CVS_CCESA[ilc[23]]+_CVS_CCESA[idc[23]]+_CVS_CCESA[15]*slope+_CVS_CCESA[16]*sxcos+_CVS_CCESA[17]*sxsin+_CVS_CCESA[18]*alnel
    ccon[24] = _CVS_CCEVE[1]+_CVS_CCEVE[idc[24]]+_CVS_CCEVE[ihc[24]]+_CVS_CCEVE[ilc[24]]+_CVS_CCEVE[15]*slope
    ccon[25] = _CVS_CCOST[1]+_CVS_CCOST[idc[25]]+_CVS_CCOST[ihc[25]]+_CVS_CCOST[11]*sxcos+_CVS_CCOST[12]*sxsin+_CVS_CCOST[13]*slope
    ccon[26] = _CVS_CHODI[1]+_CVS_CHODI[idc[26]]+_CVS_CHODI[ihc[26]]+_CVS_CHODI[11]*sxcos+_CVS_CHODI[12]*sxsin+_CVS_CHODI[13]*slope
    ccon[27] = _CVS_CPREM[1]+_CVS_CPREM[2]*elev
    ccon[28] = _CVS_CPRVI[1]
    ccon[29] = _CVS_CSALX[1]+_CVS_CSALX[ihc[29]]+_CVS_CSALX[idc[29]]+_CVS_CSALX[ipc[29]]+_CVS_CSALX[16]*alnel
    ccon[30] = _CVS_CSAMB[1]+_CVS_CSAMB[ihc[30]]+_CVS_CSAMB[7]*slope
    ccon[31] = _CVS_CSORB[1]+_CVS_CSORB[2]*slope

    return ShrubConst(pcon, (tcon1, tcon2), hcon, ccon)
end

# ---- from cover_shrub_cvbrow.jl ----
# cover_shrub_cvbrow.jl — CVBROW per-cycle shrub predictions (vcovr/cvbrow.f).
# Appended to cover_shrub.jl at merge time.  RESIDC/RESIDH = 0 (no-calibration path;
# the SHRBLAYR/SHRUBHT/SHRUBPC calibration terms are deferred with CVBCAL).

# POTHER(4,31): [TCOV, BA, 1/SAGE, BA/SAGE]  (POTH1/POTH2)
const _CVB_POTHER = reshape(Float32[
  0.011558,-0.003476,-4.717070,-0.015309, 0.013073,-0.006486,-4.885640,0.058535,
  0.017331,0.008484,-2.697824,0.018135, 0.016048,-0.001808,-2.666703,-0.081739,
  0.017789,0.003092,0.755702,-0.015269, 0.006343,0.016910,-4.969932,-0.077932,
 -0.003413,0.005090,-7.844033,-0.020012, 0.012374,-0.004095,-0.422440,0.006656,
  0.015961,-0.001766,2.348577,0.023978, 0.019857,0.003355,2.063106,-0.033853,
  0.008987,-0.012976,-3.379406,0.042276, 0.013767,0.005637,4.700924,-0.024066,
  0.015314,-0.006543,0.647053,0.001973, 0.014120,0.001780,-20.147840,0.034757,
  0.017693,0.001188,0.381343,0.000390, 0.013852,0.000796,0.872434,0.029311,
  0.019001,-0.026023,0.629758,0.127159, -0.002198,-0.011778,1.334858,0.052637,
  0.012473,-0.000119,2.681394,-0.001130, 0.020034,-0.000620,3.167380,0.014952,
  0.018084,-0.003382,-5.644269,0.007870, 0.014920,-0.001082,-1.860893,0.020735,
  0.020936,-0.013398,-0.968741,0.013742, 0.018872,-0.021482,-5.591056,0.016878,
  0.007739,-0.014455,-2.976664,0.046610, 0.017762,0.003649,-1.424454,-0.013778,
  0.021698,-0.005112,-8.778732,0.024828, 0.017873,-0.004751,-1.536393,-0.009653,
  0.014568,0.002879,-5.250633,-0.062526, 0.003294,-0.014991,1.209158,0.038918,
  0.019810,-0.004930,3.130198,-0.007978], 4, 31)

# PDTSD(4,31): 1/SAGE*disturbance-type [NONE, MECH, BURN, ROAD]  (PDTSD1/PDTSD2)
const _CVB_PDTSD = reshape(Float32[
 -6.107311,5.902798,0.204514,0.0, -4.723021,3.865583,0.857437,0.0,
  0.945456,5.703588,-6.649044,0.0, -2.480139,1.423352,1.056787,0.0,
  0.512271,-0.751467,0.239196,0.0, 1.471227,0.730630,-2.201857,0.0,
 -2.008472,1.021768,0.986705,0.0, -0.259681,-0.052514,0.312195,0.0,
  0.961347,2.235664,-3.197011,0.0, -1.837033,0.455432,1.381601,0.0,
  0.725868,-2.567226,1.841358,0.0, -0.905012,-0.585221,1.490232,0.0,
 -0.094418,-1.813987,1.908404,0.0, -5.870269,17.217840,-11.347580,0.0,
 -1.099482,-2.018118,3.117599,0.0, -0.533537,0.883687,-0.350149,0.0,
 -1.197717,3.193539,-1.995821,0.0, -0.306601,-1.556538,1.863138,0.0,
  3.543876,-2.077322,-1.466554,0.0, 1.083944,-1.509908,0.425963,0.0,
  5.573135,-0.359563,-5.213573,0.0, -1.058366,2.130129,-1.071762,0.0,
 -2.008768,1.897653,0.111116,0.0, -2.802399,-1.026632,3.829032,0.0,
  0.971046,-7.029215,6.058168,0.0, 0.462000,4.296114,-4.758114,0.0,
  7.517502,-2.886484,-4.631018,0.0, 1.465091,-2.525816,1.060724,0.0,
 -0.984921,2.809616,-1.824695,0.0, 0.372300,-1.068950,0.696649,0.0,
 -1.836557,0.267258,1.569299,0.0], 4, 31)

const _CVB_TDTSD  = Float32[-0.02001682,0.01211542,0.00767395,0.00]  # by IDIST
const _CVB_H11DTS = Float32[0.001316,0.002771,0.001187,-0.005310]
const _CVB_H16DTS = Float32[-0.013230,0.049648,0.143607,-0.102255]

struct ShrubCycle
    pgt0::Float32
    tcov::Float32
    pb::Vector{Float32}
    sh::Vector{Float32}
    cv::Vector{Float32}
    pbcv::Vector{Float32}
    cabht::Vector{Float32}
    totlcv::Float32
    htindx::Vector{Int}
end

const _CVB_ZERO31 = zeros(Float32, 31)

function cvbrow_cycle(sc::ShrubConst, sage::Float32, ba::Float32, idist::Int;
                      cal::Union{ShrubCalib,Nothing}=nothing, lstart::Bool=false)
    pcon = sc.pcon; hcon = sc.hcon; ccon = sc.ccon
    tcon1, tcon2 = sc.tcon
    E = _f32exp

    pgt0 = 1.0f0/(1.0f0 + E(-(tcon1 + 0.0138f0*sage - 0.00412f0*ba - 0.000213f0*ba*sage)))
    tcov = E(tcon2 + 0.02885483f0*sage - 0.00020194f0*ba*sage + sage*_CVB_TDTSD[idist] + 0.3238953f0)

    # RESIDC (calibration method 1): SUMCVR - TCOV, set at cycle 0, carried through.
    if cal !== nothing && cal.lcal1 && lstart
        cal.residc = cal.sumcvr - tcov
    end
    residc = cal === nothing ? 0.0f0 : cal.residc

    pb = Vector{Float32}(undef, 31)
    @inbounds for i in 1:31
        cscal = _CVB_POTHER[1,i]*tcov + _CVB_POTHER[2,i]*ba + _CVB_POTHER[3,i]/sage +
                _CVB_POTHER[4,i]*ba/sage + _CVB_PDTSD[idist,i]/sage
        pb[i] = 1.0f0/(1.0f0 + E(-(pcon[i] + cscal)))
    end

    barea = ba <= 1.0f0 ? 1.0f0 : ba
    alnba = _f32log(barea)
    aldts = _f32log(sage)

    sh = Vector{Float32}(undef, 31)
    sh[1]  = hcon[1]
    sh[2]  = hcon[2]
    sh[3]  = hcon[3]
    sh[4]  = E(hcon[4]-0.582791f0/sage+0.001952f0*residc+0.007388f0*tcov+0.07339f0)
    sh[5]  = E(hcon[5]+0.002937f0*tcov+0.002149f0*residc-0.000621f0*ba+0.06798f0)
    sh[6]  = E(hcon[6]+0.068334f0*alnba+0.174239f0*aldts+0.007149f0*tcov+0.007753f0*residc+0.16350f0)
    sh[7]  = hcon[7]
    sh[8]  = E(hcon[8]+0.001631f0*residc+0.001889f0*tcov+0.178257f0*aldts-0.012620f0*sage+0.023587f0*alnba-0.000341f0*ba+0.06231f0)
    sh[9]  = hcon[9]+0.070653f0*alnba+0.267706f0*aldts-0.005634f0*sage+0.016909f0*tcov+0.008538f0*residc
    sh[10] = hcon[10]+0.010595f0*alnba-1.313471f0/sage+0.009158f0*residc+0.015028f0*tcov
    sh[11] = E(hcon[11]+0.001716f0*residc+0.003271f0*tcov-0.014447f0*ba/sage+ba*_CVB_H11DTS[idist]+0.05623f0)
    sh[12] = E(hcon[12]-0.009299f0*alnba+0.001613f0*tcov+0.001325f0*residc+0.07756f0)
    sh[13] = E(hcon[13]+0.039276f0*alnba-0.001362f0*ba+0.024601f0*aldts+0.004278f0*tcov+0.003161f0*residc+0.07116f0)
    sh[14] = E(hcon[14]-0.692003f0/sage+0.036054f0*alnba+0.003684f0*residc+0.004300f0*tcov+0.02645f0)
    sh[15] = E(hcon[15]-0.000522f0*ba+0.094560f0*aldts+0.005903f0*tcov+0.002175f0*residc+0.10828f0)
    sh[16] = E(hcon[16]+0.000797f0*ba+_CVB_H16DTS[idist]*aldts+0.004205f0*tcov+0.003230f0*residc+0.06272f0)
    sh[17] = E(hcon[17]-0.941824f0/sage-0.002898f0*tcov+0.000816f0*residc+0.03773f0)
    sh[18] = E(hcon[18]+0.08319f0)
    sh[19] = E(hcon[19]+0.004244f0*tcov+0.003993f0*residc+0.073122f0*alnba+0.334423f0*aldts+0.29502f0)
    sh[20] = E(hcon[20]+0.010526f0*tcov+0.003868f0*residc-0.391740f0/sage-0.000384f0*ba+0.11944f0)
    sh[21] = E(hcon[21]+0.003185f0*residc+0.008998f0*tcov+0.103704f0*aldts+0.002357f0*ba+0.09106f0)
    sh[22] = E(hcon[22]+0.228213f0*aldts-0.013028f0*sage+0.003357f0*ba-0.000003f0*ba*ba+0.013626f0*tcov+0.002758f0*residc+0.12963f0)
    sh[23] = E(hcon[23]+0.002668f0*tcov+0.003559f0*residc-1.708042f0/sage+0.002716f0*(ba-11.945f0)-0.000027f0*(ba-11.945f0)*(ba-11.945f0)+0.05472f0)
    sh[24] = hcon[24]+0.025771f0*tcov+0.017823f0*residc-1.349370f0/sage
    sh[25] = E(hcon[25]+0.002177f0*tcov+0.002429f0*residc+0.06477f0)
    sh[26] = E(hcon[26]+0.038519f0*alnba+0.274757f0*aldts-0.012626f0*sage+0.004241f0*tcov+0.002252f0*residc+0.07819f0)
    sh[27] = E(hcon[27]+0.008813f0*tcov+0.005733f0*residc+0.08837f0)
    sh[28] = E(hcon[28]+0.004620f0*residc-3.798911f0/sage+0.128902f0)
    sh[29] = E(hcon[29]+0.005564f0*residc+0.010411f0*tcov+0.115970f0*aldts+0.004220f0*(ba-20.119f0)-0.000023f0*(ba-20.119f0)*(ba-20.119f0)+0.11118f0)
    sh[30] = E(hcon[30]+0.163050f0*aldts-0.002211f0*tcov+0.001493f0*residc+0.04912f0)
    sh[31] = E(hcon[31]+0.004520f0*tcov+0.001510f0*residc-0.173325f0/sage+0.08724f0)

    # RESIDH (calibration method 2): SHRBHT - SH at cycle 0, carried through (cvbrow DO 48).
    if cal !== nothing && lstart
        @inbounds for i in 1:31
            cal.shrbht[i] >= 0.0f0 && (cal.residh[i] = cal.shrbht[i] - sh[i])
        end
    end
    residh = cal === nothing ? _CVB_ZERO31 : cal.residh

    alnht = Vector{Float32}(undef, 31)
    @inbounds for i in 1:31
        shh = sh[i] <= 0.0f0 ? 0.1f0 : sh[i]
        alnht[i] = _f32log(shh)
    end

    htindx = Vector{Int}(undef, 31)
    rdpsrt!(31, sh, htindx, true)

    cv    = zeros(Float32, 31)
    pbcv  = zeros(Float32, 31)
    cabht = zeros(Float32, 31)
    totlcv = 0.0f0
    pbcum = 0.0f0
    @inbounds for j in 1:31
        k = htindx[j]
        if j == 1
            cabht[k] = 0.0f0
        else
            pbcum += pbcv[htindx[j-1]]
            cabht[k] = pbcum
        end
        cv[k] = _cvb_cv(k, ccon, sh, cabht[k], alnht, alnba, aldts, sage, ba, residh)
        pbcv[k] = pb[k]*cv[k]
        totlcv += pbcv[k]
    end
    return ShrubCycle(pgt0, tcov, pb, sh, cv, pbcv, cabht, totlcv, htindx)
end

# per-species % cover (cvbrow.f labels 51..531).  RESIDH(k) (calibration method 2)
# enters at its exact Fortran position; residh = 0 for the no-calibration path (x+0f0 = x
# exactly, so the uncalibrated result is byte-identical).
@inline function _cvb_cv(k::Int, ccon, sh, cbh::Float32, alnht, alnba::Float32,
                         aldts::Float32, sage::Float32, ba::Float32, residh)
    E = _f32exp
    c = ccon[k]; r = residh
    if k == 1;      return E(c-0.007130f0*cbh+0.23176f0)
    elseif k == 2;  return E(c+0.078842f0)
    elseif k == 3;  return E(c+0.006036f0*sage+0.28364f0)
    elseif k == 4;  return E(c+1.128531f0*alnht[4]+0.208575f0*r[4]+0.003494f0*sage+0.001445f0*ba+0.123307f0)
    elseif k == 5;  return E(c+5.617948f0*alnht[5]-2.102114f0*sh[5]+0.083065f0*r[5]-0.00488f0*cbh+0.002444f0*ba+0.195726f0)
    elseif k == 6;  return E(c+0.112634f0*alnba-0.010360f0*cbh+0.620216f0*alnht[6]+0.009096f0*r[6]+0.253561f0)
    elseif k == 7;  return E(c-0.004303f0*cbh+0.306118f0)
    elseif k == 8;  return E(c-0.042881f0*alnba+0.745749f0*sh[8]-0.039454f0*sh[8]*sh[8]+0.102498f0*r[8]+0.147753f0)
    elseif k == 9;  return 200.0f0/(E(c-3.612059f0*alnht[9]-0.199540f0*r[9]+0.009250f0*cbh+0.010110f0*sage)+1.0f0)
    elseif k == 10; return 200.0f0/(E(c-0.002451f0*(ba-30.573f0)+0.000035f0*(ba-30.573f0)*(ba-30.573f0)+0.248107f0*aldts+0.013535f0*cbh-4.936937f0*alnht[10]-0.123863f0*r[10])+1.0f0)
    elseif k == 11; return E(c-0.005229f0*cbh+2.786906f0*alnht[11]+0.145728f0*r[11]+0.418976f0*aldts-0.057897f0*sage+0.153522f0)
    elseif k == 12; return E(c-0.0034056f0*sage+1.493548f0*alnht[12]+0.083533f0*r[12]-0.002010f0*cbh+0.09473f0)
    elseif k == 13; return E(c-0.069935f0*alnba-0.007007f0*cbh+1.429518f0*sh[13]-0.084391f0*sh[13]*sh[13]+0.202652f0*r[13]+0.156056f0)
    elseif k == 14; return c-0.151413f0*cbh+25.163600f0*alnht[14]+0.603551f0*r[14]
    elseif k == 15; return E(c+1.511617f0*alnht[15]-0.251851f0*sh[15]+0.087134f0*r[15]-0.003701f0*cbh-0.091528f0*aldts-0.052584f0*alnba+0.249182f0)
    elseif k == 16; return 200.0f0/(E(c+0.003013f0*ba-0.147043f0*r[16]-5.284045f0*alnht[16]+0.468977f0*sh[16]+0.014973f0*cbh+0.119072f0*aldts)+1.0f0)
    elseif k == 17; return E(c-0.004972f0*cbh+1.385094f0*alnht[17]+0.089854f0*r[17]-0.226805f0*aldts+0.283218f0)
    elseif k == 18; return E(c-0.005018f0*cbh+0.245277f0)
    elseif k == 19; return E(c+0.040967f0*sh[19]+0.044127f0*r[19]+0.267191f0)
    elseif k == 20; return 200.0f0/(E(c+0.007244f0*sage-0.037231f0*alnba+0.005470f0*cbh-1.580698f0*alnht[20]-0.066696f0*r[20])+1.0f0)
    elseif k == 21; return 200.0f0/(E(c+0.011024f0*cbh-2.590103f0*alnht[21]-0.076644f0*r[21]+0.005074f0*ba)+1.0f0)
    elseif k == 22; return 200.0f0/(E(c-0.642331f0*alnht[22]-0.065580f0*r[22]+0.002359f0*cbh)+1.0f0)
    elseif k == 23; return 200.0f0/(E(c-3.761727f0*alnht[23]-0.120113f0*r[23]+0.009849f0*cbh+0.817357f0*aldts+0.118734f0*alnba)+1.0f0)
    elseif k == 24; return 200.0f0/(E(c+0.002371f0*ba+0.016261f0*sage-3.633245f0*alnht[24]-0.184807f0*r[24]+0.017470f0*cbh)+1.0f0)
    elseif k == 25; return E(c+0.160346f0)
    elseif k == 26; return E(c-0.003541f0*ba-0.094480f0*aldts-0.004642f0*cbh+0.814695f0*sh[26]-0.044342f0*sh[26]*sh[26]+0.090125f0*r[26]+0.172664f0)
    elseif k == 27; return E(c-0.011936f0*cbh+0.718063f0*sh[27]-0.024411f0*sh[27]*sh[27]+0.020491f0*r[27]+0.252399f0)
    elseif k == 28; return E(c+0.172636f0*sh[28]-0.004438f0*sh[28]*sh[28]+0.068645f0*r[28]+0.185518f0)
    elseif k == 29; return 200.0f0/(E(c+0.181524f0*aldts-1.260551f0*alnht[29]-0.039281f0*r[29]+0.006094f0*cbh)+1.0f0)
    elseif k == 30; return E(c+0.300947f0*sh[30]+0.082204f0*r[30]-0.004439f0*ba+0.172504f0)
    else;           return E(c-0.004196f0*cbh+2.298975f0*alnht[31]+0.129692f0*r[31]-0.059705f0*alnba+0.196114f0)
    end
end

# =====================================================================================
# CVBCAL (covr/cvbcal.f) — compute shrub height/cover correction factors BHTCF/BPCCF
# from user-observed shrub data, at cycle 0.  Two methods: by layer (SHRBLAYR, LCAL1)
# or by species (SHRUBHT/SHRUBPC, LCAL2).  Report-only downstream (DGSD=0 path).
# =====================================================================================
function cvbcal!(cal::ShrubCalib, sh::Vector{Float32}, pbcv::Vector{Float32},
                 pb::Vector{Float32}, htindx::Vector{Int}, totlcv::Float32)
    @inbounds for i in 1:31
        cal.xcv[i] = pbcv[i]; cal.xpb[i] = pb[i]; cal.xsh[i] = sh[i]
    end
    cal.lcal2 ? _cvbcal_m2!(cal, sh, pbcv) : _cvbcal_m1!(cal, sh, pbcv, htindx, totlcv)
    cal.computed = true
    return cal
end

# Method 1: calibration by shrub layer (cvbcal.f lines 98-238).
function _cvbcal_m1!(cal::ShrubCalib, sh, pbcv, htindx::Vector{Int}, totlcv::Float32)
    nklass = cal.nklass
    relpc = Vector{Float32}(undef, 3)
    @inbounds for i in 1:nklass
        relpc[i] = 0.00001f0 + cal.avgbpc[i]/cal.sumcvr
    end
    iht = zeros(Int, 3, 2)
    nextsp = 1
    j = 0; k = 0
    reached130 = false
    @inbounds for i in 1:nklass
        if nextsp > 31                       # GO TO 130
            reached130 = true
            break
        end
        iht[i, 1] = nextsp
        cumuli = 0.0f0; rcumpc = 0.0f0
        i == nklass && break                 # last group → GO TO 110 (handled below)
        # DO 70: find the class-1 boundary
        dropflag = false; keepflag = false
        j = nextsp
        while j <= 31
            k = htindx[j]
            diff1 = abs(rcumpc - relpc[i])
            cumuli += pbcv[k]
            rcumpc = cumuli/totlcv
            diff2 = abs(rcumpc - relpc[i])
            if relpc[i] == 0.0f0             # GO TO 80 (drop)
                dropflag = true; break
            end
            if rcumpc <= relpc[i]            # GO TO 70 (continue)
                j += 1; continue
            end
            if diff1 <= diff2                # GO TO 80 (drop)
                dropflag = true; break
            else                             # GO TO 90 (keep)
                keepflag = true; break
            end
        end
        if !keepflag                         # DO 70 fell through OR drop → 80
            j > 31 && (j = 32; k = htindx[31])   # Fortran J=32, K=HTINDX(31) after loop
            cumuli -= pbcv[k]
            nextsp = j
            iht[i, 2] = j - 1                # LASTSP
        else                                 # 90 keep
            nextsp = j + 1
            iht[i, 2] = j                    # LASTSP
        end
    end
    reached130 || (iht[nklass, 2] = 31)      # 110: last group = the rest
    # DO 150: per-class predicted mean height and cover.
    sumcv = zeros(Float32, 3); sumht = zeros(Float32, 3)
    @inbounds for i in 1:nklass
        cal.avgbpc[i] <= 0.0f0 && continue
        cal.cvavg[i] = 0.0f0; cal.htavg[i] = 0.0f0
        cal.cvfrac[i] = 0.0f0; cal.htfrac[i] = 0.0f0
        i1 = iht[i, 1]; i2 = iht[i, 2]
        i2 < i1 && continue
        for jj in i1:i2
            kk = htindx[jj]
            cal.ilayr[kk] = i
            sumcv[i] += pbcv[kk]
            sumht[i] += sh[kk]
        end
        cal.cvavg[i] = sumcv[i]
        cal.htavg[i] = sumht[i]/Float32((i2 - i1) + 1)
    end
    # DO 180: ratios observed/predicted → correction factors per species.
    @inbounds for i in 1:nklass
        cal.avgbpc[i] <= 0.0f0 && continue
        cal.cvavg[i] > 0.0f0 && (cal.cvfrac[i] = cal.avgbpc[i]/cal.cvavg[i])
        cal.htavg[i] > 0.0f0 && (cal.htfrac[i] = cal.avgbht[i]/cal.htavg[i])
        i1 = iht[i, 1]; i2 = iht[i, 2]
        i2 < i1 && continue
        for jj in i1:i2
            kk = htindx[jj]
            cal.bhtcf[kk] = cal.htfrac[i]
            cal.bpccf[kk] = cal.cvfrac[i]
        end
    end
    return cal
end

# Method 2: calibration by individual species (cvbcal.f lines 239-295).
function _cvbcal_m2!(cal::ShrubCalib, sh, pbcv)
    @inbounds for ispi in 1:31
        hh = cal.shrbht[ispi]; pc = cal.shrbpc[ispi]
        if pc >= 0.0f0
            if pc <= 0.0f0                   # pc == 0 → both factors 0
                cal.bpccf[ispi] = 0.0f0; cal.bhtcf[ispi] = 0.0f0
                continue
            end
            cal.bpccf[ispi] = pc/pbcv[ispi]  # 240: pc > 0
            hh < 0.0f0 && continue           # height not observed → BHTCF stays 1
            cal.bhtcf[ispi] = hh > 0.0f0 ? hh/sh[ispi] : 0.0f0
        else                                 # pc < 0 (cover not observed)
            hh < 0.0f0 && continue           # both not observed → unchanged
            if hh > 0.0f0
                cal.bhtcf[ispi] = hh/sh[ispi]  # 210: BPCCF stays 1
            else                             # hh == 0
                cal.bpccf[ispi] = 0.0f0; cal.bhtcf[ispi] = 0.0f0
            end
        end
    end
    return cal
end

# Apply BHTCF/BPCCF to the cycle predictions (cvbrow.f DO 63 / DO 65) → new totlcv.
# Returns a calibrated ShrubCycle (report/aggregation consume the scaled arrays).
function _cvbcal_apply(cal::ShrubCalib, cyc::ShrubCycle)
    sh = copy(cyc.sh); pb = copy(cyc.pb); pbcv = copy(cyc.pbcv)
    totlcv = 0.0f0
    if cal.lcal1
        @inbounds for i in 1:31
            sh[i] = sh[i]*cal.bhtcf[i]
            pbcv[i] = pbcv[i]*cal.bpccf[i]
            totlcv += pbcv[i]
        end
    else  # lcal2
        @inbounds for i in 1:31
            if !(cal.shrbht[i] != 0.0f0 && cal.shrbpc[i] != 0.0f0)
                pb[i] = 0.0f0; sh[i] = 0.0f0; pbcv[i] = 0.0f0
            end
            sh[i] = sh[i]*cal.bhtcf[i]
            pbcv[i] = pbcv[i]*cal.bpccf[i]
            totlcv += pbcv[i]
        end
    end
    return ShrubCycle(cyc.pgt0, cyc.tcov, pb, sh, cyc.cv, pbcv, cyc.cabht, totlcv, cyc.htindx)
end

# ---- from cover_shrub_cvsum.jl ----
# cover_shrub_cvsum.jl — SHRUB aggregation half of covr/cvsum.f + covr/cvclas.f.
# Groups the 31 species into LOW(1-7)/MED(8-19)/TALL(20-31), picks the top-3 by cover
# (RDPSRT) with the rest combined into "OTHR", and computes CLOW/CMED/CTALL, ASHT,
# TALLSH, and the shrub-cover-by-height SCOV vector.  Report-only.

const _CVSUM_NUM  = (7,12,12)
const _CVSUM_INC  = (0,7,19)
const _CVSUM_IGSP = (0,4,8)
const _CVSUM_SHTRHT = Float32[0.5,1.0,2.0,3.0,4.0,5.0,7.5,10.0,15.0,20.0,400.0]

# Species abbreviations (32 = OTHR)
const _CV_SNAME = ["ARUV","BERB","LIBO","PAMY","SPBE","VASC","CARX","LONI","MEFE","PHMA",
                   "RIBE","ROSA","RUPA","SHCA","SYMP","VAME","XETE","FERN","COMB","ACGL",
                   "ALSI","AMAL","CESA","CEVE","COST","HODI","PREM","PRVI","SALX","SAMB",
                   "SORB","OTHR"]

# Aggregated shrub-summary values for one cycle (report inputs).
struct ShrubSummary
    asht::Float32          # weighted avg shrub height (all species)
    tallsh::Float32        # avg "tall shrub" (sp 20-31) height, weighted by PBCV
    clow::Float32          # Σ PBCV low group
    cmed::Float32
    ctall::Float32
    issp::Vector{Int}      # 12 species codes (top-3+OTHR × 3 groups)
    scv::Vector{Float32}   # 12 covers
    sht::Vector{Float32}   # 12 heights
    spb::Vector{Float32}   # 12 probs (×100)
    scov::Vector{Float32}  # 11 cumulative shrub cover by height
end

function cvsum_shrub(cyc::ShrubCycle)
    sh = cyc.sh; pb = cyc.pb; pbcv = cyc.pbcv
    height = 0.0f0; sumpb = 0.0f0; tsht = 0.0f0; tsumpb = 0.0f0
    scov = zeros(Float32, 11)
    @inbounds for isp in 1:31
        sumpb  += pbcv[isp]
        height += sh[isp]*pbcv[isp]
        if isp >= 20
            tsumpb += pbcv[isp]
            tsht   += sh[isp]*pbcv[isp]
        end
        hm1 = -1.0f0
        for ih in 1:11
            hc = _CVSUM_SHTRHT[ih]
            (sh[isp] > hm1 && sh[isp] <= hc) && (scov[ih] += pbcv[isp])
            hm1 = hc
        end
    end
    asht   = sumpb  <= 1f-4 ? 0.0f0 : height/sumpb
    tallsh = tsumpb <= 1f-4 ? 0.0f0 : tsht/tsumpb
    # cumulate scov from top (SCOV[j] += SCOV[j+1], j=10..1)
    @inbounds for ih in 1:10
        j = 11 - ih
        scov[j] += scov[j+1]
    end

    issp = Vector{Int}(undef, 12)
    scv  = Vector{Float32}(undef, 12)
    sht  = Vector{Float32}(undef, 12)
    spb  = Vector{Float32}(undef, 12)
    clow = cmed = ctall = 0.0f0
    for ig in 1:3
        num = _CVSUM_NUM[ig]; inc = _CVSUM_INC[ig]; igsp = _CVSUM_IGSP[ig]
        othcov = 0.0f0; othht = 0.0f0; othpb = 0.0f0
        pass = Vector{Float32}(undef, num)
        @inbounds for i in 1:num
            k = i + inc
            pass[i] = pbcv[k]
            othcov += pbcv[k]
            othht  += sh[k]*pbcv[k]
            othpb  += pb[k]
        end
        ig == 1 && (clow = othcov); ig == 2 && (cmed = othcov); ig == 3 && (ctall = othcov)
        indx = Vector{Int}(undef, num)
        rdpsrt!(num, pass, indx, true)
        local j = 0
        @inbounds for i in 1:3
            k = indx[i] + inc
            othcov -= pbcv[k]; othht -= sh[k]*pbcv[k]; othpb -= pb[k]
            j = i + igsp
            scv[j] = pbcv[k]; spb[j] = 100.0f0*pb[k]; sht[j] = sh[k]; issp[j] = k
        end
        jp1 = j + 1
        scv[jp1] = othcov
        sht[jp1] = othcov <= 1f-4 ? 0.0f0 : othht/othcov
        spb[jp1] = 100.0f0*othpb
        issp[jp1] = 32
    end
    return ShrubSummary(asht, tallsh, clow, cmed, ctall, issp, scv, sht, spb, scov)
end

# CVCLAS successional-stage classifier (covr/cvclas.f).  Shrub-dominated stages 1-6 use
# SAGE/RELDEN/AVH/TALLSH/TOTLCV; tree-dominated stages 7-10 use a selective stand DBH
# (DIAM) derived from the crown-area-by-height-class tallies (CRXHT/SD2XHT/TXHT) built by
# the canopy CVSUM tree loop.  Returns ISTAGE 1..10.
function cvclas(sage::Float32, relden::Float32, avh::Float32, tallsh::Float32,
                totlcv::Float32, crxht::Vector{Float32}, sd2xht::Vector{Float32},
                txht::Vector{Float32}, crarea::Float32, rmsqd::Float32,
                htmax::Float32, htmin::Float32)
    sage > 40.0f0 && return 0
    # --- shrub-dominated stages ---
    if sage <= 5.0f0
        if !((relden > 30.0f0) && (avh > 1.0f0))
            if !((tallsh > 1.0f0) || (totlcv > 25.0f0))
                return 1
            end
        end
    end
    shrub_to_tree = false
    if (relden > 30.0f0) && (avh > 1.0f0)
        shrub_to_tree = true
    else
        if tallsh <= 2.5f0; return 2; end
        if tallsh <= 5.0f0; return 3; end
        if totlcv >= 70.0f0; return 4; end
        shrub_to_tree = true
    end
    if shrub_to_tree
        # stage 5 / 6 gates before tree-dominated
        stage5 = true
        if (tallsh < 5.0f0 && totlcv < 50.0f0); stage5 = false; end
        if avh > 5.0f0; stage5 = false; end
        if relden > 50.0f0; stage5 = false; end
        stage5 && return 5
        stage6 = true
        if (tallsh < 5.0f0 && totlcv < 30.0f0); stage6 = false; end
        if relden > 100.0f0; stage6 = false; end
        stage6 && return 6
    end
    # --- tree-dominated: selective DBH (DIAM) then stage 7-10 ---
    diam = _cvclas_diam(crxht, sd2xht, txht, crarea, rmsqd, htmax, htmin)
    diam > 4.0f0  || return 7
    diam > 11.0f0 || return 8
    diam > 24.0f0 || return 9
    return 10
end

function _cvclas_diam(crxht, sd2xht, txht, crarea::Float32, rmsqd::Float32,
                      htmax::Float32, htmin::Float32)
    htrnge = htmax - htmin
    if htrnge <= 30.0f0
        return rmsqd
    end
    indx = Vector{Int}(undef, 16)
    rdpsrt!(16, crxht, indx, true)
    one = indx[1]; two = indx[2]; three = indx[3]
    test = crxht[one]/crarea
    if test >= 0.50f0
        return sqrt(sd2xht[one]/txht[one])
    end
    top2 = crxht[one]+crxht[two]; test = top2/crarea
    itest = max(one,two) - min(one,two)
    if test >= 0.70f0 && itest <= 2
        return sqrt((sd2xht[one]+sd2xht[two])/(txht[one]+txht[two]))
    end
    top3 = crxht[one]+crxht[two]+crxht[three]; test = top3/crarea
    itest = max(one,two,three) - min(one,two,three)
    if test >= 0.90f0 && itest <= 2
        return sqrt((sd2xht[one]+sd2xht[two]+sd2xht[three])/(txht[one]+txht[two]+txht[three]))
    end
    # two-storied test
    crax20_1 = 0.0f0; mtag1 = 1
    @inbounds for i in 1:15
        d = crxht[i]+crxht[i+1]
        if d >= crax20_1; crax20_1 = d; mtag1 = i; end
    end
    crax20_2 = 0.0f0; mtag2 = 1; it = mtag1+1
    @inbounds for i in 1:15
        (i+1 == mtag1 || i == mtag1 || i == it) && continue
        d = crxht[i]+crxht[i+1]
        if d >= crax20_2; crax20_2 = d; mtag2 = i; end
    end
    if abs(mtag1-mtag2) >= 4
        mp1 = mtag1+1
        return sqrt((sd2xht[mtag1]+sd2xht[mp1])/(txht[mtag1]+txht[mp1]))
    end
    # all-aged
    return sqrt((sd2xht[one]+sd2xht[two]+sd2xht[three])/(txht[one]+txht[two]+txht[three]))
end

# ---- from cover_shrub_wire.jl ----
# cover_shrub_wire.jl — habitat mapping + per-cycle shrub row + report emit.  Appended
# to cover.jl.  Habitat/union/physiography index tables from vcovr/cvbrow.f.

# Valid habitat codes and their overstory/understory/total-cover union subscripts.
const _CVB_IHCODE = Int[210,220,260,310,320,330,340,380,390,395,
                        505,510,511,515,520,525,530,540,550,570,
                        590,620,635,645,650,670,690,705,710,720,
                        721,730,790,830]
const _CVB_IPHABO = Int[1,1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2, 3,3,3, 4, 2, 5,5,5,5,5,5,5,5,5,5,5,5,5]
const _CVB_IPHABU = Int[5,5,2,3,5,5,3,3,2,3, 3,1,1,2,1,2,1,1,1,4, 1,4,1,2,5,2,1,3,1,2, 2,2,5,1]
const _CVB_ITHABU = Int[6,6,2,4,6,6,4,4,2,4, 4,1,1,2,1,2,3,3,3,5, 1,5,1,2,6,2,1,4,1,2, 2,1,6,1]
const _CVB_INFOR  = Int[3,4,3,3,4,4,4,4,3,2,3]   # IFOR→INF (IE); other variants fixed

# INF (national-forest grouping) by variant, per cvbrow.f SELECT CASE(VARACD).
function _cvb_inf(varcode::AbstractString, ifor::Int)
    varcode == "EM" && return 3
    varcode == "CI" && return (ifor == 1 ? 2 : 1)
    if varcode == "IE"
        (1 <= ifor <= 11) && return _CVB_INFOR[ifor]
        return 1
    end
    return 1                                     # CASE DEFAULT → Boise/Payette
end

# One projection cycle's shrub statistics (ITHN=1 pre-thin only).
struct ShrubRow
    year::Int
    sage::Float32          # TIMESD
    pgt0::Float32
    clow::Float32; cmed::Float32; ctall::Float32; totlcv::Float32
    asht::Float32; tallsh::Float32
    sbmass::Float32; twigs::Float32
    istage::Int
    issp::Vector{Int}; scv::Vector{Float32}; sht::Vector{Float32}; spb::Vector{Float32}
    scov::Vector{Float32}
    trsh::Vector{Float32}   # trees/ac by SHTRHT height threshold (SHRUB-SMALL CONIFER COMPETITION)
end

# Recover the raw ICL5 habitat CODE (e.g. 260) used by cvbrow.f as the default IHTYPE.
# EM overwrites s.plot.habitat_code with IEMTYP (the 1..118 bucket index into EM_JTYPE)
# during em_sitset!, so map it back through EM_JTYPE; other variants keep the raw code.
function _cover_stand_habitat(s::StandState)
    hc = Int(s.plot.habitat_code)
    if s.variant isa EasternMontana && 1 <= hc <= length(EM_JTYPE)
        return Int(EM_JTYPE[hc])
    end
    return hc
end

# Resolve habitat code → (iov, iun, itun, ok).  ok=false ⇒ invalid habitat (LBROW off).
function _cvb_habitat(ihtype::Int)
    idx = findfirst(==(ihtype), _CVB_IHCODE)
    idx === nothing && return (0, 0, 0, false)
    return (_CVB_IPHABO[idx], _CVB_IPHABU[idx], _CVB_ITHABU[idx], true)
end

# ---- from cover_shrub_report.jl ----
# cover_shrub_report.jl — CVOUT SHRUB STATISTICS + CANOPY AND SHRUBS SUMMARY emit.
# Header text is verbatim from cvout.f FORMATs (byte-identical to FVSem_g16 except the
# GROHED page-header/timestamp line, which is excluded from the parity check).

_f61(x::Float32) = @sprintf("%6.1f", x)     # Fortran F6.1
_f71(x::Float32) = @sprintf("%7.1f", x)     # Fortran F7.1
_f81(x::Float32) = @sprintf("%8.1f", x)     # Fortran F8.1

# place string `s` into buf ending at 1-based column `endcol` (right-justified field).
function _place!(buf::Vector{Char}, endcol::Int, s::AbstractString)
    n = length(s)
    @inbounds for k in 1:n
        c = endcol - n + k
        1 <= c <= length(buf) && (buf[c] = s[k])
    end
    return buf
end

const _CV_SHRUB_HDR = [
"                   LOW SPECIES (0-1.7 FT)        MEDIUM SPECIES (1.7-7 FT)         TALL SPECIES (7+ FT)",
"                -----------------------------   ---------------------------     --------------------------",
"                ARUV:ARCTOSTAPHYLOS UVA-URSI    LONI:LONICERA SPP.              ACGL:ACER GLABRUM         ",
"                BERB:BERBERIS SPP.              MEFE:MENZIESIA FERRUGINEA       ALSI:ALNUS SINUATA        ",
"                LIBO:LINNAEA BOREALIS           PHMA:PHYSOCARPUS MALVACEUS      AMAL:AMELANCHIER ALNIFOLIA",
"                PAMY:PACHISTIMA MYRSINITES      RIBE:RIBES SPP.                 CESA:CEANOTHUS SANGUINEUS ",
"                SPBE:SPIRAEA BETULIFOLIA        ROSA:ROSA SPP.                  CEVE:CEANOTHUS VELUTINUS  ",
"                VASC:VACCINIUM SCOPARIUM        RUPA:RUBUS PARVIFLORUS          COST:CORNUS STOLONIFERA   ",
"                CARX:CAREX SPP.                 SHCA:SHEPHERDIA CANADENSIS      HODI:HOLODISCUS DISCOLOR  ",
"                                                SYMP:SYMPHORICARPOS SPP.        PREM:PRUNUS EMARGINATA    ",
"                                                VAME:VACCINIUM MEMBRANACEUM     PRVI:PRUNUS VIRGINIANA    ",
"                                                XETE:XEROPHYLLUM TENAX          SALX:SALIX SPP.           ",
"                                                FERN:FERNS                      SAMB:SAMBUCUS SPP.        ",
"                                                COMB:OTHER SHRUBS COMBINED      SORB:SORBUS SPP.          ",
]
const _CV_SHRUB_HDR2 = [
"                        --------------------------------------------------------------------------------",
"                         ATTRIBUTES OF THE FIRST THREE SPECIES WITH GREATEST COVER IN EACH HEIGHT GROUP",
"                                    (ALL OTHERS WITHIN GROUP COMBINED INTO CATEGORY \"OTHR\")",
"                        --------------------------------------------------------------------------------",
"                                 COVER -- SPECIES COVER ",
"                                HEIGHT -- AVERAGE SPECIES HEIGHT (FEET)  ",
"                                  PROB -- SPECIES PROBABILITY OF OCCURRENCE ",
]
const _CV_SUM_STAGE = [
"                         -------- DEFINITIONS OF SUCCESSIONAL STAGE CODES USED IN OUTPUT --------",
"                         1: RECENT DISTURBANCE                 6: TALL SHRUB WITH MOSTLY CONIFERS",
"                         2: LOW SHRUB                          7: SAPLING TIMBER",
"                         3: MEDIUM SHRUB                       8: POLE TIMBER",
"                         4: TALL SHRUB WITH NO CONIFERS        9: MATURE TIMBER",
"                         5: TALL SHRUB WITH FEW CONIFERS      10: OLD-GROWTH TIMBER",
]
const _CV_SUM_HDR = [
"         ------------------- UNDERSTORY ATTRIBUTES -------------------         ------------ OVERSTORY ATTRIBUTES -------------",
"         TIME     PROB.   -----SHRUB COVER----  AVG.    DORMANT                                                 SUM OF        ",
"         SINCE    (SHRUB                        SHRUB   SHRUB    TWIGS  SUCC.  STAND  TOP     CANOPY   FOLIAGE  STEM    NUMBER",
"         DISTURB. COV>0)  LOW  MED TALL TOTAL   HEIGHT  BIOMASS  (NO./  STAGE  AGE    HEIGHT   COVER   BIOMASS  DIAMS.   OF   ",
"DATE     (YEARS)   (%)    (%)  (%)  (%)  (%)    (FEET)  (LB/AC)  SQFT)  CODE   (YRS)  (FEET)    (%)    (LB/AC)  (FEET)  STEMS ",
]

# Place a preformatted string at an absolute column (Fortran T-descriptor).
function _cvcol!(buf::Vector{Char}, col::Int, s::AbstractString)
    n = length(s)
    while length(buf) < col - 1 + n; push!(buf, ' '); end
    @inbounds for (i, c) in enumerate(s); buf[col + i - 1] = c; end
    return buf
end

# SHRUB MODEL CALIBRATION STATISTICS table (covr/cvout.f 1000/1010/1020/1030/1040 for
# LCAL1; 1050/1060-1063 for LCAL2).  Emitted once, before the SHRUB STATISTICS table,
# when calibration is in effect.
function _cover_calib_stats(cv::CoverState, io::IO, stand_id, mgmt_id, title)
    cal = cv.calib
    _cover_stand_header(io, stand_id, mgmt_id, title)
    println(io); println(io)
    println(io, "-"^46, "  SHRUB MODEL CALIBRATION STATISTICS  ", "-"^47)
    println(io); println(io); println(io)
    if cal.lcal1
        println(io, "CALIBRATION BY SHRUB LAYER (SHRBLAYR KEYWORD CARD):")
        println(io)
        println(io, "          AVERAGE HEIGHT (FEET)                    AVERAGE PERCENT COVER       ")
        println(io, "    -----------------------------------     -----------------------------------")
        println(io, "    SHRUB  OBSERVED  PREDICTED  SCALING     SHRUB  OBSERVED  PREDICTED  SCALING")
        println(io, "    LAYER  VALUES    VALUES     FACTORS     LAYER  VALUES    VALUES     FACTORS")
        println(io, "    -----  --------  ---------  -------     -----  --------  ---------  -------")
        println(io)
        for i in 1:cal.nklass
            buf = Char[]
            _cvcol!(buf, 5,  @sprintf("%4d", i))
            _cvcol!(buf, 9,  @sprintf("%10.1f", cal.avgbht[i]))
            _cvcol!(buf, 19, @sprintf("%11.1f", cal.htavg[i]))
            _cvcol!(buf, 30, @sprintf("%9.2f", cal.htfrac[i]))
            _cvcol!(buf, 46, @sprintf("%4d", i))
            _cvcol!(buf, 50, @sprintf("%10.1f", cal.avgbpc[i]))
            _cvcol!(buf, 60, @sprintf("%11.1f", cal.cvavg[i]))
            _cvcol!(buf, 71, @sprintf("%9.2f", cal.cvfrac[i]))
            println(io, rstrip(String(buf)))
        end
        println(io); println(io)
        println(io, "                               HEIGHT         % COVER")
        println(io, "      SHRUB     ASSIGNED       SCALING        SCALING")
        println(io, "     SPECIES     LAYER         FACTOR         FACTOR")
        println(io, "     -------    --------      ---------      ---------")
        println(io)
        for i in 1:31
            buf = Char[]
            _cvcol!(buf, 8,  _CV_SNAME[i])
            _cvcol!(buf, 20, @sprintf("%1d", cal.ilayr[i]))
            _cvcol!(buf, 33, @sprintf("%5.2f", cal.bhtcf[i]))
            _cvcol!(buf, 48, @sprintf("%5.2f", cal.bpccf[i]))
            println(io, rstrip(String(buf)))
        end
    else
        println(io, "CALIBRATION BY INDIVIDUAL SPECIES (SHRUBHT AND/OR SHRUBPC KEYWORD CARDS):")
        println(io)
        println(io, "                        SHRUB HEIGHT (FEET)                        PERCENT COVER")
        println(io, "                  ----------------------------------     ----------------------------------")
        println(io, "     SHRUB         OBSERVED     PREDICTED    SCALING      OBSERVED     PREDICTED    SCALING")
        println(io, "    SPECIES          VALUE        VALUE      FACTOR         VALUE        VALUE      FACTOR")
        println(io, "    -------       -----------   ---------   --------     -----------   ---------   --------")
        for i in 1:31
            hh = cal.shrbht[i]; pc = cal.shrbpc[i]
            buf = Char[]
            _cvcol!(buf, 7, _CV_SNAME[i])
            hh != -99999.0f0 && _cvcol!(buf, 20, @sprintf("%7.1f", hh))
            _cvcol!(buf, 33, @sprintf("%7.1f", cal.xsh[i]))
            _cvcol!(buf, 45, @sprintf("%7.2f", cal.bhtcf[i]))
            pc != -99999.0f0 && _cvcol!(buf, 59, @sprintf("%7.1f", pc))
            _cvcol!(buf, 72, @sprintf("%7.1f", cal.xcv[i]))
            _cvcol!(buf, 84, @sprintf("%7.2f", cal.bpccf[i]))
            println(io, rstrip(String(buf)))
        end
    end
    return io
end

function _cover_shrub_stats(cv::CoverState, io::IO, stand_id, mgmt_id, title)
    _cover_stand_header(io, stand_id, mgmt_id, title)
    println(io); println(io)
    println(io, "-"^54, "  SHRUB STATISTICS  ", "-"^55)
    println(io, " "^52, "(BASED ON STOCKABLE AREA)")
    println(io)
    for l in _CV_SHRUB_HDR; println(io, l); end
    println(io)
    for l in _CV_SHRUB_HDR2; println(io, l); end
    println(io)
    println(io, "YEAR                     LOW SPECIES             MEDIUM SPECIES             TALL SPECIES           USER-SELECTED INDICATOR SPECIES")
    println(io, "----               ----------------------    ----------------------    ----------------------    ----------------------------------")
    println(io)
    for r in cv.shrub_rows
        r.sage > 40.0f0 && continue
        println(io, lpad(r.year, 4))
        sp = "SPECIES" * " "^10; cvl = "  COVER" * " "^10
        htl = " HEIGHT" * " "^10; pbl = "   PROB" * " "^10
        for ig in 1:3
            for kk in 1:4
                j = (ig-1)*4 + kk
                sp  *= "  " * rpad(_CV_SNAME[r.issp[j]], 4)
                cvl *= _f61(r.scv[j]); htl *= _f61(r.sht[j]); pbl *= _f61(r.spb[j])
            end
            sp *= "  "; cvl *= "  "; htl *= "  "; pbl *= "  "
        end
        println(io, rstrip(sp)); println(io, rstrip(cvl))
        println(io, rstrip(htl)); println(io, rstrip(pbl))
        println(io)
    end
    return io
end

function _cover_summary(cv::CoverState, io::IO, stand_id, mgmt_id, title)
    _cover_stand_header(io, stand_id, mgmt_id, title)
    println(io); println(io)
    println(io, "-"^48, "  CANOPY AND SHRUBS SUMMARY  ", "-"^49)
    println(io, " "^50, "(BASED ON STOCKABLE AREA)")
    println(io); println(io); println(io)
    if cv.lshrub
        for l in _CV_SUM_STAGE; println(io, l); end
        println(io); println(io)
    end
    for l in _CV_SUM_HDR; println(io, l); end
    println(io, "-"^126)
    for (i, sr) in enumerate(cv.shrub_rows)
        cr = i <= length(cv.rows) ? cv.rows[i] : nothing
        buf = fill(' ', 126)
        _place!(buf, 4,  string(sr.year))
        if sr.sage <= 40.0f0
            _place!(buf, 15, string(_cv_ifix(sr.sage)))
            _place!(buf, 22, string(_cv_ifix(sr.pgt0*100f0)))
            _place!(buf, 29, string(_cv_ifix(sr.clow)))
            _place!(buf, 34, string(_cv_ifix(sr.cmed)))
            _place!(buf, 39, string(_cv_ifix(sr.ctall)))
            _place!(buf, 44, string(_cv_ifix(sr.totlcv)))
            _place!(buf, 53, _f71(sr.asht))
            _place!(buf, 62, string(_cv_ifix(sr.sbmass)))
            _place!(buf, 70, _f81(sr.twigs))
            _place!(buf, 76, string(sr.istage))
        end
        if cr !== nothing
            _place!(buf, 83,  string(cr.icvage))
            _place!(buf, 90,  string(_cv_ifix(cr.stdht)))
            _place!(buf, 98,  string(_cv_ifix(cr.tpctcv)))
            _place!(buf, 109, string(_cv_ifix(cr.totbms)))
            _place!(buf, 117, string(_cv_ifix(cr.sdiam)))
            _place!(buf, 125, string(_cv_ifix(cr.tretot)))
        end
        println(io, rstrip(String(buf)))
    end
    println(io, "-"^126)
    return io
end

# CVOUT SHRUB-SMALL CONIFER COMPETITION table (cvout.f 9090/9092/9094).  Emitted after
# the CANOPY AND SHRUBS SUMMARY (same page, no header), gated on the shrub model (LBROW)
# and KODE≠3.  Per-cycle rows with TIMESD (sage) ≤ 40 only (ITHN=1 pre-thin).  Report-only.
function _cover_shrub_conifer(cv::CoverState, io::IO)
    println(io); println(io)
    println(io, "-"^45, "  SHRUB-SMALL CONIFER COMPETITION  ", "-"^46)
    println(io)
    println(io, " "^30, "SHRUB COVER -- TOTAL COVER OF SHRUBS GREATER THAN HEIGHT")
    println(io, " "^30, "TREES/ACRE  -- TOTAL NUMBER OF TREES PER ACRE GREATER THAN HEIGHT")
    println(io)
    println(io, "-"^126)
    println(io, " "^56, "  HEIGHT (FEET)")
    println(io, "YEAR", " "^24, "0.0    0.5    1.0    2.0    3.0    4.0    5.0    ",
            "7.5   10.0   15.0   20.0")
    println(io, "-"^126)
    for r in cv.shrub_rows
        r.sage > 40.0f0 && continue
        l1 = lpad(r.year, 4) * "        SHRUB COVER "
        l2 = "            TREES/ACRE  "
        @inbounds for j in 1:11
            l1 *= lpad(_cv_ifix(r.scov[j]), 7)
            l2 *= lpad(_cv_ifix(r.trsh[j]), 7)
        end
        println(io, l1); println(io, l2); println(io)
    end
    println(io, "-"^126)
    return io
end
