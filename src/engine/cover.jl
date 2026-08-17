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
end

CoverState() = CoverState(false, false, false, false, true, true, true,
                          2, 1, 0, 0.0f0, Float32[], CoverRow[])

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
        elseif k == "SHRUBS"          # shrub half deferred; flag only
            cv.lbrow = true
        elseif k == "SHRBLAYR" || k == "SHRUBHT" || k == "SHRUBPC"
            # shrub calibration cards — deferred (consumed; SHRUBHT/PC read 4 data recs)
            if k != "SHRBLAYR"
                for _ in 1:4
                    rr = read_keyword!(kr)
                    (rr.status == KW_EOF || rr.status == KW_STOP) && break
                end
            end
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
    sdiam = 0.0f0

    @inbounds for i in 1:itrn
        dbh  = t.dbh[i]
        ht   = t.height[i]
        icr  = Int(t.crown_pct[i])       # ICR (integer %)
        prob = t.tpa[i]
        dg   = t.diam_growth[i]
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
    return cv
end

# =====================================================================================
# CVOUT — emit the "CANOPY COVER STATISTICS" table (canopy path; report-only).
# =====================================================================================
_cv_ifix(x) = trunc(Int, 0.5f0 + x)     # FVS IFIX(.5+x)

function cover_report(cv::CoverState, io::IO, stand_id::AbstractString,
                      mgmt_id::AbstractString, title::AbstractString)
    (cv.lcover && cv.lcnop && !isempty(cv.rows)) || return io
    println(io)
    println(io, "STAND ID: ", rpad(stand_id, 26), "    MGMT CODE: ",
            rpad(mgmt_id, 4), "  ", title)
    println(io)
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
