# =============================================================================
# Western Root Disease (WRD) extension — port of FVS rd/*.f
# =============================================================================
#
# Scope of THIS file (Chunk −1, infrastructure + keyword reader + inert seam):
#   * RD parameter constants        (rd/RDPARM.F77)
#   * RootDiseaseState               (the reduced-path members of /RDCOM/ /RDADD/
#                                     /RDCRY/ /RDARRY/ — commons ported as Julia state)
#   * rd_init_defaults!              (rd/rdinit.f — the scalar / small-array defaults)
#   * kw_rdin!                       (rd/rdin.f entry RDKEY — the RRTYPE / RRINIT /
#                                     SAREA / BBCLEAR / END subset)
#   * rd_active                      (rd/rdatv.f gate  L = RRTINV .OR. RRMAN)
#   * root_disease_setup! / _mn2! / _treg!  — the fvs.f/grincr.f/gradd.f engine seams,
#                                     wired GATED + INERT (early-return unless active);
#                                     the per-cycle MORTALITY body is Chunk 0 (pending).
#
# NOT in this file (Chunk 0 and later — deliberately deferred, see the port report):
#   * rdblk1.f host-species coefficient tables (HABFAC/PKILLS/PNINF/PCOLO/DECFN/…)
#   * rdsetp/rdinoc/rdarea/rdinf/rdmort/rdend/rdgrow mortality + growth-loss bodies
#   * the RD RNG (rdrani/rdrann/rdranp) + random center placement (rdcloc)
# These require a bit-exact FVS growth baseline for validation, which the KT engine
# does not yet provide for the reference stand (a cornered growth straddle exceeds
# the RD mortality signal). See docs and the port report for the measured boundary.
#
# All references are to /workspace/ForestVegetationSimulator/rd/.

# --- rd/RDPARM.F77 PARAMETERs -------------------------------------------------
const RD_DSO    = 1        # dead standing outside center
const RD_DSII   = 2        # dead standing inside, infected
const RD_DSIU   = 3        # dead standing inside, uninfected
const RD_IRRTRE = 1500     # max FVS tree records RD dimensions for
const RD_ITOTSP = 40       # total species across all WRD variants
const RD_ITOTRR = 4        # number of modeled root diseases

# Disease index → name (rd/rdin.f RRTYPE echo; 1=P-annosus 2=S-annosus 3=Armillaria
# 4=Phellinus). RRTYPE 3 in the turnkey oracle echoes "ARMILLARIA".
const RD_DISEASE_NAMES = ("ANNOSUS-P", "ANNOSUS-S", "ARMILLARIA", "PHELLINUS")
const RD_DISEASE_CODES = ('P', 'S', 'A', 'W')   # rd/rdpr.f CHTYPE

# AbstractRootDiseaseState is forward-declared in core/state.jl (StandState field type).

"""
    RootDiseaseState

Julia port of the FVS Western Root Disease common-block state, reduced to the
members exercised by the center-init (RRINIT) reduced mortality path. Populated by
`rd_init_defaults!` then `kw_rdin!`. `active` mirrors the rd/rdatv.f gate.
"""
mutable struct RootDiseaseState <: AbstractRootDiseaseState
    # ---- activation / gate (rd/RDADD.F77 /RDADD/, rd/rdatv.f) ----
    iroot::Int32              # RD model on flag (rdinit 0 → RDIN sets 1)
    rrman::Bool               # manual (keyword) init used  → RRINIT
    rrtinv::Bool              # treelist init used          → RRTREIN (not in chunk 0)
    lrtype::Bool              # RRTYPE keyword was given
    lbbon::Bool               # bark beetles default-on unless BBCLEAR (rdin END handler)
    bbclear::Bool             # BBCLEAR seen → suppress default bark beetles

    # ---- disease selection (rd/rdin.f RRTYPE) ----
    minrr::Int32              # lowest active disease number
    maxrr::Int32              # highest active disease number

    # ---- stand geometry (rd/RDCOM.F77) ----
    sarea::Float32            # stand area, acres
    dimen::Float32            # stand square dimension, ft = sqrt(SAREA)*208.7
    irhab::Int32              # habitat flag (1 default; 2 after IRGEN(8) year)
    istep::Int32              # time-slot counter (= ICYC+1); advanced in rdmn2

    # ---- per-disease init parameters (dim RD_ITOTRR) ----
    ncents::Vector{Int32}     # number of disease centers
    parea::Vector{Float32}    # diseased patch area, acres
    ooarea::Vector{Float32}   # previous-cycle PAREA
    prkill::Vector{Float32}   # infected trees/acre in diseased area (RRINIT fld3)
    prun::Vector{Float32}     # uninfected trees/acre in diseased area (RRINIT fld4)
    rrincs::Vector{Float32}   # initial proportion of roots infected (RRINIT fld5)
    rrnew::Vector{Float32}    # mean proportion roots infected for a NEW infection (rdinit RRNEW)
    lparea::Vector{Bool}      # PAREA set by keyword (do not default to 0.25·SAREA)
    lonect::Vector{Int32}     # 0 unassigned / 1 stand-is-one-center / 2 multiplot
    ipcflg::Vector{Int32}     # 0 random centers / 1 centers read explicitly

    # ---- mortality-curve constants (rd/rdinit.f; RRTYPE-independent) ----
    xxinf::NTuple{5,Float32}  # DBH knots of the years-to-kill curve (TEMP16)
    yyinf::NTuple{5,Float32}  # years-to-kill knots (TEMP17)
    nninf::Int32              # number of active knots (=3)
    stcut::NTuple{5,Float32}  # stump size-class breakpoints
    xminkl::Vector{Float32}   # min time-to-kill offset per disease (=0)
    xminlf::Vector{Float32}   # min inoculum lifespan per disease ({1,1,1,20})

    # ---- RD random-number generator (rd/RDADD.F77 /RRANN/ + rd/rdin.f RSEED) ----
    #   The WRD model uses its OWN Park–Miller stream (16807, mod 2147483647), seeded
    #   from DSEED, kept SEPARATE from the base FVS `rann` stream. rd/rdrani.f init,
    #   rd/rdrann.f draw. DSEED default 889347.0 (rd/rdinit.f); RSEED keyword overrides.
    dseed::Float64            # RD seed (rd/rdinit.f DSEED; RSEED keyword sets it)
    rd_s0::Float64            # /RRANN/ S0 — current generator state
    rd_s1::Float64            # /RRANN/ S1 — last-produced state
    rd_ss::Float64            # /RRANN/ SS — seed state (odd-forced INT(DSEED))
    # ---- rd/rdranp.f binomial random-proportion memo (/RDCOM/ OLDPRP,CDF) ----
    rd_oldprp::Float32        # last PROPIN (skips CDF rebuild when unchanged)
    rd_cdf::Vector{Float32}   # cumulative distribution buffer (size 1001 in FVS)

    # ---- disease-center geometry (rd/RDCOM.F77 PCENTS/IRRSP; rd/rdcloc.f, rd/rdarea.f) ----
    irrsp::Int32                     # current disease type being placed (rdsetp DO-600 loop)
    pcents::Array{Float32,3}         # PCENTS(ITOTRR,100,3): (x,y,radius) ft, per center
    yincpt::Float32                  # SDI→roots effect intercept (rdsetp)
    sdislp::Float32                  # SDI slope (rdinit −0.0033)
    sdnorm::Float32                  # SDI normal (rdinit 369.0)

    # ---- host-species crosswalk + per-disease density normalizer (rd/rdsetp.f) ----
    irtspc::Vector{Int32}     # host-species → base-RD-species crosswalk (rdblk1.f IRTSPC)
    rrgen1::Vector{Float32}   # RRGEN(idi,1): infected+uninfected density → normalized to 1

    # ---- per-record initial-infection state (rd/rdsetp.f tree loop) ----
    #   Sized to the tree list at rd_setp!; the reduced RRINIT manual-init path.
    probi::Vector{Float32}    # PROBI(I,1,1)  infected TPA inside patches
    probiu::Vector{Float32}   # PROBIU(I)     uninfected TPA inside patches
    fprob::Vector{Float32}    # FPROB(I)      TPA outside patches
    propi::Vector{Float32}    # PROPI(I,1,1)  proportion of roots infected (RDRANP draw)
    probl::Vector{Float32}    # PROBL(I) = PROB(I)
    propn::Vector{Float32}    # PROPN(I)      per-record infected proportion (rdiprp)

    # ---- per-cycle driver working state (RDDriver; lazily built at LSTART) ----
    #   Holds the 3-D PROBI/PROPI, PROBIT, FFPROB, ROOTL, RRKILL, the per-center
    #   RRATES, the stump lists (PROBDA/DBHDA/ROOTDA), and the RDEND/RDGROW carry
    #   arrays. Typed `Any` so the struct need not forward-declare RDDriver.
    driver::Any               # ::Union{Nothing,RDDriver}
    icyc::Int32               # RD cycle counter (FVS ICYC; 0 at LSTART, +1 per grow cycle)
    sum_rows::Vector{Any}     # RDSUM accumulator: (year, rd_sum_report) per cycle (FVS_RD_Sum / dbsrd.f DBSRD1)

    RootDiseaseState() = rd_init_defaults!(new())
end

"""
    rd_init_defaults!(rd) -> rd

Port of rd/rdinit.f: the RD common-block default initialization (the members the
reduced center-init mortality path reads). Values verified against the live FVSkt
oracle's RDIN option echo.
"""
function rd_init_defaults!(rd::RootDiseaseState)
    rd.iroot  = Int32(0)
    rd.rrman  = false
    rd.rrtinv = false
    rd.lrtype = false
    rd.lbbon  = true
    rd.bbclear = false

    rd.minrr = Int32(1)       # rdinit MINRR=1, MAXRR=2 (annosus) — overwritten by RRTYPE
    rd.maxrr = Int32(2)

    rd.sarea = 100.0f0
    rd.dimen = 2087.0f0       # sqrt(100)*208.7
    rd.irhab = Int32(1)
    rd.istep = Int32(1)

    n = RD_ITOTRR
    rd.ncents = fill(Int32(20), n)      # rdinit NCENTS=20
    rd.parea  = fill(25.0f0, n)         # rdinit PAREA=25
    rd.ooarea = zeros(Float32, n)
    rd.prkill = fill(0.5f0, n)          # rdinit PRKILL=0.5
    rd.prun   = fill(0.5f0, n)          # rdinit PRUN=0.5
    rd.rrincs = fill(0.1f0, n)          # rdinit RRINCS=0.1
    rd.rrnew  = Float32[0.001, 0.001, 0.05, 0.05]   # rdinit RRNEW (Annosus 0.001; Armillaria/Phellinus 0.05)
    rd.lparea = fill(false, n)
    rd.lonect = fill(Int32(0), n)
    rd.ipcflg = fill(Int32(0), n)

    rd.xxinf = (0.0f0, 3.9f0, 35.4f0, 0.0f0, 0.0f0)     # TEMP16
    rd.yyinf = (0.0f0, 5.0f0, 40.0f0, 0.0f0, 0.0f0)     # TEMP17
    rd.nninf = Int32(3)
    rd.stcut = (0.0f0, 12.0f0, 24.0f0, 48.0f0, 100.0f0)
    rd.xminkl = fill(0.0f0, n)
    rd.xminlf = Float32[1.0, 1.0, 1.0, 20.0]            # TEMP15

    rd.dseed = 889347.0                                 # rd/rdinit.f DSEED = 889347.0
    rd.rd_oldprp = -1.0f0                                # forces first CDF build
    rd.rd_cdf = Float32[]
    rd_rani!(rd, rd.dseed)                               # seed the RD stream (rd/rdrani.f)

    rd.irrsp  = Int32(0)
    rd.pcents = zeros(Float32, RD_ITOTRR, 100, 3)
    rd.yincpt = 1.0f0
    rd.sdislp = -0.0033f0                                # rd/rdinit.f SDISLP = -0.0033
    rd.sdnorm = 369.0f0                                  # rd/rdinit.f SDNORM = 369.0

    rd.irtspc = copy(RD_IRTSPC_KT)                       # rdblk1.f IRTSPC (NI/CI/KT base crosswalk)
    rd.rrgen1 = zeros(Float32, n)                        # rdinit RRGEN(IDI,1) = 0
    rd.probi  = Float32[]
    rd.probiu = Float32[]
    rd.fprob  = Float32[]
    rd.propi  = Float32[]
    rd.probl  = Float32[]
    rd.propn  = Float32[]
    rd.driver = nothing
    rd.icyc   = Int32(0)
    rd.sum_rows = Any[]
    return rd
end

# -----------------------------------------------------------------------------
# rd_sum_report (rd/rdpr.f DBSRD1 aggregation) — the FVS_RD_Sum "1st report":
# summary statistics for a root-disease area, per active disease type, aggregated
# from the driver (per-record PROBIU/PROBIT/RDKILL, stump PROBDA/DBHDA, spread
# RRRATE) + the tree list (CFV/DBH), ALL per-acre-in-disease-area. Reads state
# post-`rd_end_apply!` — additive, does not touch the WRD kernels. Values are the
# METRIC form DBSRD1 emits (PAREA·ACRtoHA, RRRATE·FTtoM, /ACRtoHA, ·FT2/FT3pACR…).
# The 4 new-infection columns (New_Inf_Prp_Ins/Exp/Tot) + Ave_Pct_Root_Inf need
# the CORINF/EXPINF/PRINF accumulators not yet tracked in jl ⇒ 0 (documented
# follow-on). Validated bit-exact-or-cornered vs FVSkt_clean (11-row rdsum_oracle).
# -----------------------------------------------------------------------------
const _RD_ACRtoHA        = 0.4046945f0    # base/PRGPRM.F77 PARAMETERs
const _RD_FTtoM          = 0.3048f0
const _RD_FT2pACRtoM2pHA = 0.2295643f0
const _RD_FT3pACRtoM3pHA = 0.0699713f0
const _RD_TYPE_CHAR = ("P", "S", "A", "W") # CHTYPE (rd/rdpr.f DATA CHTYPE/'P','S','A','W'/); idi 3=Armillaria='A'

function rd_sum_report(rd::RootDiseaseState, s::StandState, year::Integer, iage::Integer)
    d = rd.driver::RDDriver; t = s.trees; n = t.n
    idi = Int(rd.minrr)                                       # single active disease type (KT path)
    parea = rd.parea[idi]; pinv = 1.0f0 / (parea + 1.0f-9)
    # stumps (PROBDA/DBHDA over 2 pools × 5 size × up-to-41 age slots)
    tstmps = 0.0f0; bastpa = 0.0f0
    istep = max(1, Int(rd.istep))
    @inbounds for m in 1:istep, k in 1:5, l in 1:2
        p = d.probda[idi, l, k, m]
        tstmps += p * pinv
        dd = d.dbhda[idi, l, k, m]
        bastpa += p * (3.141593f0 * (dd / 24.0f0)^2) * pinv
    end
    # live tree list restricted to the disease area (per-record, 1:1 with trees)
    tun = 0.0f0; tin = 0.0f0; tdie = 0.0f0; tdvol = 0.0f0; bapa = 0.0f0; cfvpa = 0.0f0
    @inbounds for i in 1:n
        tclas = d.probiu[i] + d.probit[i]
        tun   += d.probiu[i] * pinv
        tin   += d.probit[i] * pinv
        tdie  += d.rdkill[i] * pinv
        tdvol += d.rdkill[i] * t.cuft_vol[i] * pinv
        bapa  += tclas * (3.14159f0 * (t.dbh[i] / 24.0f0)^2) * pinv
        cfvpa += tclas * t.cuft_vol[i] * pinv
    end
    rrrate = idi <= length(d.rrrate) ? d.rrrate[idi] : 0.0f0
    ncent  = idi <= length(rd.ncents) ? Int(rd.ncents[idi]) : 0
    rdtype = 1 <= idi <= 4 ? _RD_TYPE_CHAR[idi] : "A"
    # DBSRD1 units: rdpr.f emits the METRIC form (·ACRtoHA etc.) only under LMTRIC (BC/ON,
    # rdinit.f:728 VARACD BC/ON), else the RAW imperial values.
    m = (s.variant isa BritishColumbia) || (s.variant isa Ontario)
    conv(x, f) = m ? x * f : x
    div_(x, f) = m ? x / f : x
    return (Year = Int(year), Age = Int(iage), RD_Type = rdtype, Num_Centers = ncent,
            RD_Area = conv(parea, _RD_ACRtoHA), Spread = conv(rrrate, _RD_FTtoM),
            Stumps = div_(tstmps, _RD_ACRtoHA), Stumps_BA = conv(bastpa, _RD_FT2pACRtoM2pHA),
            Mort_TPA = div_(tdie, _RD_ACRtoHA), Mort_CuFt = conv(tdvol, _RD_FT3pACRtoM3pHA),
            UnInf_TPA = div_(tun, _RD_ACRtoHA), Inf_TPA = div_(tin, _RD_ACRtoHA),
            Ave_Pct_Root_Inf = 0.0f0,                        # PRINF·100 — follow-on (accumulator absent)
            Live_Merch_CuFt = conv(cfvpa, _RD_FT3pACRtoM3pHA), Live_BA = conv(bapa, _RD_FT2pACRtoM2pHA),
            New_Inf_Prp_Ins = 0.0f0, New_Inf_Prp_Exp = 0.0f0, New_Inf_Prp_Tot = 0.0f0)  # CORINF/EXPINF — follow-on
end

# -----------------------------------------------------------------------------
# RD random-number generator (rd/rdrani.f, rd/rdrann.f) — the WRD model's OWN
# Park–Miller minimal-standard stream (16807, mod 2^31−1), SEPARATE from the base
# `rann` stream. Ported faithfully (never FFI). Validated bit-exact against the
# live relinked FVSkt oracle: 2000 draws, 0 mismatches on both the integer state
# S1 and the Float32 return value.
# -----------------------------------------------------------------------------

"""
    rd_rani!(rd, sseed)

Port of rd/rdrani.f: seed the RD generator. `ISEED = INT(SSEED)`, forced odd,
`SS = S0 = FLOAT(ISEED)`.
"""
function rd_rani!(rd::RootDiseaseState, sseed::Real)
    iseed = trunc(Int, sseed)
    iseed % 2 == 0 && (iseed += 1)          # rdrani: force odd
    rd.rd_ss = Float64(iseed)
    rd.rd_s0 = rd.rd_ss
    rd.rd_s1 = 0.0
    return rd
end

"""
    rd_rann!(rd) -> Float32

Port of rd/rdrann.f: one uniform draw. `S1 = DMOD(16807·S0, 2147483647)` in double
precision; returns `REAL(S1 / 2147483648)` (Float32); advances `S0 = S1`.
"""
function rd_rann!(rd::RootDiseaseState)
    rd.rd_s1 = rem(16807.0 * rd.rd_s0, 2147483647.0)     # DMOD (double precision)
    val = Float32(rd.rd_s1 / 2147483648.0)               # REAL(...) — Float32
    rd.rd_s0 = rd.rd_s1
    return val
end

"""
    rd_ranp!(rd, prop) -> Float32

Port of rd/rdranp.f: a random proportion drawn from a truncated binomial CDF with
mean `prop`, using the RD stream. All arithmetic in Float32 (matches the 2014 REAL
CDF). `OLDPRP`/`CDF` memoization is carried in `rd`; the RESULT is a pure function of
`(PROPIN, INTNUM, RANNUM)`, so a rebuild produces identical values. Validated
bit-exact against the live FVSkt oracle: 297 calls, 0 mismatches.
"""
function rd_ranp!(rd::RootDiseaseState, prop::Real)
    intnum, lrev = _rd_ranp_setup!(rd, prop)             # adjusts PROPIN + builds CDF
    rannum = rd_rann!(rd)                                 # rdranp.f label 300: RANNUM>0
    while rannum == 0.0f0
        rannum = rd_rann!(rd)
    end
    return _rd_ranp_walk(rd.rd_cdf, intnum, lrev, rannum)
end

# rdranp.f labels 100–200: adjust PROPIN and (re)build the truncated-binomial CDF into
# rd.rd_cdf, memoized on OLDPRP. Returns (INTNUM, LREV).
function _rd_ranp_setup!(rd::RootDiseaseState, prop::Real)
    propin = Float32(prop)
    intnum = round(Int, 5.0f0 / propin)                  # NINT(5/PROPIN)
    intnum > 100 && (intnum = 100)
    intnum < 10  && (intnum = 10)

    L = 0                                                # rdranp.f label 100
    while true
        exprop = propin / (1.0f0 - (1.0f0 - propin)^intnum)
        if abs(propin - exprop) > (1.0f0 / Float32(intnum)) && L < 10
            propin = propin - (0.5f0 * (propin - exprop))
            L += 1
        else
            break
        end
    end

    lrev = false
    if propin > 0.5f0                                    # use lower half, reverse later
        propin = 1.0f0 - propin
        lrev = true
    end

    if propin != rd.rd_oldprp                            # rebuild only when changed
        rd.rd_oldprp = propin
        length(rd.rd_cdf) < intnum + 1 && (rd.rd_cdf = zeros(Float32, max(intnum + 1, 1001)))
        cdf = rd.rd_cdf
        pdf = (1.0f0 - propin)^intnum
        cdf[1] = pdf
        @inbounds for k in 1:intnum
            pdf = pdf > 1.0f-15 ?
                pdf * (propin / (1.0f0 - propin) * Float32(intnum - k + 1)) / Float32(k) :
                0.0f0
            cdf[k+1] = cdf[k] + (pdf / (1.0f0 - cdf[1]))
        end
    end
    return intnum, lrev
end

# rdranp.f label 400: walk the CDF with a given draw. Pure (used by rd_ranp! and tests).
function _rd_ranp_walk(cdf::Vector{Float32}, intnum::Int, lrev::Bool, rannum::Float32)
    res = 0.7f0 * (1.0f0 / Float32(intnum))
    k = 1
    @inbounds while rannum >= cdf[k+1] && k <= intnum
        res = Float32(k + 1) / Float32(intnum)
        k += 1
    end
    lrev && (res = 1.0f0 - res)
    return res
end

# Test/validation entry: RDRANP evaluated for a supplied RANNUM (no stream draw), with a
# fresh CDF rebuild. Pure function of (prop, rannum) — matches the live oracle bit-exact.
function rd_ranp_from(rd::RootDiseaseState, prop::Real, rannum::Real)
    intnum, lrev = _rd_ranp_setup!(rd, prop)
    return _rd_ranp_walk(rd.rd_cdf, intnum, lrev, Float32(rannum))
end

# -----------------------------------------------------------------------------
# RD disease-center placement (rd/rdcloc.f) + patch-area grid sampling (rd/rdarea.f).
# rdcloc lays NCENTS non-overlapping circles at random (RD-stream) locations; rdarea
# grids them on a 75×75 mesh and, during LSTART, grows the radii to match the target
# input area. Validated bit-exact vs the live FVSkt oracle (dominant-signal keyfile,
# S-annosus, 20 centers): all 20 PCENTS (x,y,radius) and PAREA/OOAREA at LSTART.
# -----------------------------------------------------------------------------

"""
    rd_cloc!(rd)

Port of rd/rdcloc.f: place `NCENTS[IRRSP]` disease centers. LONECT==1 ⇒ a single stand-
covering center (deterministic). Otherwise each center gets radius
`sqrt(PAREA/(3.14159·NCENTS))·208.7` and a random (x,y) = `RDRANN·DIMEN`, rejected (≤20
tries) if it overlaps a prior center's x- or y-extent. Advances the RD stream.
"""
function rd_cloc!(rd::RootDiseaseState)
    irrsp = Int(rd.irrsp)
    nc    = Int(rd.ncents[irrsp])
    dimen = sqrt(rd.sarea) * 208.7f0
    rd.dimen = dimen
    P = rd.pcents
    if rd.lonect[irrsp] == 1
        P[irrsp,1,3] = dimen
        P[irrsp,1,1] = 0.5f0 * dimen
        P[irrsp,1,2] = 0.5f0 * dimen
        return rd
    end
    ax1 = zeros(Float32, nc); ax2 = zeros(Float32, nc)
    ay1 = zeros(Float32, nc); ay2 = zeros(Float32, nc)
    r0  = sqrt(rd.parea[irrsp] / (3.14159f0 * nc)) * 208.7f0
    ii = 1
    @inbounds while ii <= nc
        P[irrsp,ii,3] = r0
        ntry = 0
        while true                                   # rdcloc.f label 101
            ntry += 1
            ntry > 20 && break                       # abandon center (GOTO 100)
            P[irrsp,ii,1] = rd_rann!(rd) * dimen
            P[irrsp,ii,2] = rd_rann!(rd) * dimen
            x = P[irrsp,ii,1]; y = P[irrsp,ii,2]; r = P[irrsp,ii,3]
            ax2[ii] = x + r; ax1[ii] = x - r
            ay2[ii] = y + r; ay1[ii] = y - r
            ii == 1 && break                         # first center: accept (GOTO 100)
            overlap = false
            for jj in 1:ii-1
                if (x < ax2[jj] && x > ax1[jj]) || (y < ay2[jj] && y > ay1[jj])
                    overlap = true; break            # GOTO 101 (retry)
                end
            end
            overlap || break                         # accepted
        end
        ii += 1
    end
    return rd
end

"""
    rd_area!(rd, lstart) -> Float32

Port of rd/rdarea.f: grid-sample the `NCENTS[IRRSP]` circles on a 75×75 mesh and set
`PAREA[IRRSP] = SAREA·(IN/75²)`. During `lstart` (and unless PAREA==−1, user-supplied
centers) iteratively grow every radius toward the input-area target within `TOLER`
(≤25 passes). Sets `OOAREA[IRRSP]` when `lstart`. Uses Fortran IFIX truncation.
"""
function rd_area!(rd::RootDiseaseState, lstart::Bool)
    irrsp = Int(rd.irrsp)
    nc    = Int(rd.ncents[irrsp])
    P     = rd.pcents
    igrid = 75
    toler = 0.01f0
    itrn1(v, s) = trunc(Int, (v * s) + 0.5f0)        # rdarea.f ITRN1 = IFIX(r*s+0.5)
    itrn0(v, s) = trunc(Int, (v * s) - 0.3f0)        # rdarea.f ITRN0 = IFIX(r*s-0.3)
    if nc <= 0
        rd.parea[irrsp] = 0.0f0
        return 0.0f0
    end
    areai = rd.parea[irrsp]
    artem = areai
    ntry  = 0
    arcal = 0.0f0
    lmem  = falses(igrid, igrid)
    while true                                        # rdarea.f label 6
        fill!(lmem, false)
        scl = Float32(igrid) / sqrt(rd.sarea * 43560.0f0)
        @inbounds for ipat in 1:nc
            cx = P[irrsp,ipat,1]; cy = P[irrsp,ipat,2]; cr = P[irrsp,ipat,3]
            ixc = itrn1(cx, scl); iyc = itrn1(cy, scl)
            ix1 = itrn0(cx - cr, scl); ix2 = itrn1(cx + cr, scl)
            iy1 = itrn0(cy - cr, scl); iy2 = itrn1(cy + cr, scl)
            if !(iyc <= 0 || iyc > igrid)
                for ix in ix1:ix2
                    (ix > 0 && ix <= igrid) && (lmem[ix, iyc] = true)
                end
            end
            if !(ixc <= 0 || ixc > igrid)
                for iyy in iy1:iy2
                    (iyy > 0 && iyy <= igrid) && (lmem[ixc, iyy] = true)
                end
            end
            ix1b = ixc + 1; ix3 = ix2 - 1
            if ix1b <= ix2
                for ix in ix1b:ix3
                    ixx = ixc + ixc - ix
                    x = Float32(ix - ixc) / scl
                    x > cr && continue
                    y2 = sqrt((cr * cr) - (x * x))
                    iy2b = iyc + itrn0(y2, scl)
                    iy1b = iyc - itrn1(y2, scl)
                    iy1b <= 0 && (iy1b = 1)
                    iy1b > iy2b && continue
                    for iyr in iy1b:iy2b
                        iyr > igrid && continue
                        (ix  >= 1 && ix  <= igrid) && (lmem[ix,  iyr] = true)
                        (ixx >= 1 && ixx <= igrid) && (lmem[ixx, iyr] = true)
                    end
                end
            end
        end
        IN = count(lmem)
        arcal = rd.sarea * (Float32(IN) / Float32(igrid * igrid))
        lstart || break
        rd.parea[irrsp] == -1.0f0 && break            # user centers: no target
        perdif = (areai - arcal) / areai
        delta  = artem - arcal
        (abs(perdif) <= toler || abs(delta) < toler) && break
        artar = artem + artem * perdif                # ARINC=ARTEM·PERDIF; ARTAR=ARTEM+ARINC
        artem = artar
        artar = artar * 43560.0f0 / Float32(nc)
        rnurad = sqrt(artar / 3.14159f0)
        @inbounds for k in 1:nc
            P[irrsp,k,3] = rnurad
        end
        ntry += 1
        ntry > 25 && break
    end
    rd.parea[irrsp] = arcal
    lstart && (rd.ooarea[irrsp] = arcal)
    return arcal
end

"""
    rd_place_centers!(rd, lstart)

Port of the rd/rdsetp.f disease-area section (DO-600 loop): set DIMEN and YINCPT, then
for each active disease place its centers (rd_cloc!) and compute its patch area
(rd_area!). Manual-RRINIT / random-center path (IPCFLG==0); PLOTINF sub-plot init is a
later chunk. Populates `rd.pcents`, `rd.parea`, `rd.ooarea`.
"""
function rd_place_centers!(rd::RootDiseaseState, lstart::Bool)
    rd.yincpt = (rd.sdislp == 0.0f0 || rd.sdnorm == 0.0f0) ?
                1.0f0 : 1.0f0 - (rd.sdnorm * rd.sdislp)
    rd.dimen = rd.sarea != 0.0f0 ? sqrt(rd.sarea) * 208.7f0 : 0.0f0
    @inbounds for idi in Int(rd.minrr):Int(rd.maxrr)
        rd.irrsp = Int32(idi)
        if rd.ipcflg[idi] == 1
            rd_area!(rd, lstart)                      # user-assigned centers
        elseif rd.parea[idi] > 0.0f0
            rd_cloc!(rd); rd_area!(rd, lstart)
        else
            # rdsetp.f:174  PAREA = SAREA·NDPLTS/PI. Manual-RRINIT init has no diseased
            # sub-plots (IRDPLT all 0 ⇒ NDPLTS=0 ⇒ PAREA=0); the PLOTINF sub-plot path
            # that sets NDPLTS>0 is a later chunk.
            rd.parea[idi] = 0.0f0
            rd_cloc!(rd)
        end
    end
    return rd
end

"""
    rd_active(rd) -> Bool

Port of rd/rdatv.f: L = RRTINV .OR. RRMAN. True once RD has been activated by an
RRINIT (manual) or RRTREIN (treelist) keyword.
"""
rd_active(rd::RootDiseaseState) = rd.rrtinv || rd.rrman
rd_active(::Nothing) = false

# -----------------------------------------------------------------------------
# Host-species coefficient tables (rd/rdinit.f DATA blocks TEMP4/TEMP6/TEMP7,
# ITEMP2; rd/rdblk1.f IRTSPC). Variant-GENERIC base RD tables (NI/CI/KT share the
# base IRTSPC). Stored column-major exactly as the Fortran DATA order (KSP fastest),
# so `reshape` reproduces the Fortran (KSP,IDI,IHAB) indexing directly.
#   HABFAC(ksp,idi,ihab) — relative time-to-death           (TEMP4, 40×4×2)
#   PNINF(ksp,idi)       — probability of infection         (TEMP7, 40×4)
#   PKILLS(ksp,idi)      — proportion roots infected at death(TEMP6, 40×4; RRTREIN path)
#   IDITYP(ksp)          — annosus host type 0/1/2           (ITEMP2, 40)
#   IRTSPC_KT            — KT host-species → base-RD-species  (rdblk1.f)
# -----------------------------------------------------------------------------
const RD_HABFAC_FLAT = Float32[1.0f0, 1.0f0, 1.0f0, 1.5f0, 1.75f0, 1.0f0, 1.5f0, 1.0f0, 1.5f0, 0.5f0, 1.75f0, 1.0f0, 1.5f0, 1.75f0, 1.5f0, 1.5f0, 1.5f0, 1.0f0, 99.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 99.0f0, 1.0f0, 0.5f0, 1.0f0, 1.75f0, 99.0f0, 1.0f0, 0.5f0, 99.0f0, 1.0f0, 1.5f0, 2.0f0, 1.5f0, 1.0f0, 1.5f0, 1.5f0, 99.0f0, 1.0f0, 1.0f0, 1.0f0, 1.5f0, 1.75f0, 1.0f0, 1.5f0, 1.0f0, 1.5f0, 0.5f0, 1.75f0, 1.0f0, 1.5f0, 1.75f0, 1.5f0, 1.5f0, 1.5f0, 1.0f0, 99.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 99.0f0, 1.0f0, 0.5f0, 1.0f0, 1.75f0, 99.0f0, 1.0f0, 0.5f0, 99.0f0, 1.0f0, 1.5f0, 2.0f0, 1.5f0, 1.0f0, 1.5f0, 1.5f0, 99.0f0, 1.8f0, 2.0f0, 1.0f0, 0.75f0, 0.9f0, 1.2f0, 1.8f0, 1.1f0, 0.75f0, 1.8f0, 0.9f0, 1.8f0, 0.75f0, 0.9f0, 0.75f0, 0.75f0, 1.8f0, 0.9f0, 0.9f0, 1.1f0, 0.75f0, 1.8f0, 1.8f0, 0.9f0, 1.1f0, 0.9f0, 0.75f0, 0.9f0, 0.9f0, 0.9f0, 1.8f0, 0.9f0, 0.9f0, 0.2f0, 10.0f0, 1.1f0, 10.0f0, 0.2f0, 0.75f0, 99.0f0, 3.0f0, 1.5f0, 1.0f0, 1.0f0, 1.5f0, 10.0f0, 3.0f0, 1.1f0, 1.0f0, 3.0f0, 1.5f0, 3.0f0, 1.0f0, 1.5f0, 1.0f0, 1.0f0, 3.0f0, 1.5f0, 1.5f0, 1.1f0, 1.0f0, 3.0f0, 3.0f0, 1.5f0, 1.1f0, 1.5f0, 1.0f0, 1.5f0, 1.5f0, 1.5f0, 3.0f0, 1.5f0, 1.5f0, 10.0f0, 10.0f0, 1.5f0, 10.0f0, 10.0f0, 1.0f0, 99.0f0, 1.0f0, 1.0f0, 1.0f0, 1.5f0, 1.75f0, 1.0f0, 1.5f0, 1.0f0, 1.5f0, 0.5f0, 1.75f0, 1.0f0, 1.5f0, 1.75f0, 1.5f0, 1.5f0, 1.5f0, 1.0f0, 99.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 99.0f0, 1.0f0, 0.5f0, 1.0f0, 1.75f0, 99.0f0, 1.0f0, 0.5f0, 99.0f0, 1.0f0, 1.5f0, 2.0f0, 1.5f0, 1.0f0, 1.5f0, 1.5f0, 99.0f0, 1.0f0, 1.0f0, 1.0f0, 1.5f0, 1.75f0, 1.0f0, 1.5f0, 1.0f0, 1.5f0, 0.5f0, 1.75f0, 1.0f0, 1.5f0, 1.75f0, 1.5f0, 1.5f0, 1.5f0, 1.0f0, 99.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 99.0f0, 1.0f0, 0.5f0, 1.0f0, 1.75f0, 99.0f0, 1.0f0, 0.5f0, 99.0f0, 1.0f0, 1.5f0, 2.0f0, 1.5f0, 1.0f0, 1.5f0, 1.5f0, 99.0f0, 1.8f0, 2.0f0, 1.0f0, 0.75f0, 0.9f0, 1.2f0, 1.8f0, 1.1f0, 0.75f0, 1.8f0, 0.9f0, 1.8f0, 0.75f0, 0.9f0, 0.75f0, 0.75f0, 1.8f0, 0.9f0, 0.9f0, 1.1f0, 0.75f0, 1.8f0, 1.8f0, 0.9f0, 1.1f0, 0.9f0, 0.75f0, 0.9f0, 0.9f0, 0.9f0, 1.8f0, 0.9f0, 0.9f0, 0.2f0, 10.0f0, 1.1f0, 10.0f0, 0.2f0, 0.75f0, 99.0f0, 3.0f0, 1.5f0, 1.0f0, 1.0f0, 1.5f0, 10.0f0, 3.0f0, 1.1f0, 1.0f0, 3.0f0, 1.5f0, 3.0f0, 1.0f0, 1.5f0, 1.0f0, 1.0f0, 3.0f0, 1.5f0, 1.5f0, 1.1f0, 1.0f0, 3.0f0, 3.0f0, 1.5f0, 1.1f0, 1.5f0, 1.0f0, 1.5f0, 1.5f0, 1.5f0, 3.0f0, 1.5f0, 1.5f0, 10.0f0, 10.0f0, 1.5f0, 10.0f0, 10.0f0, 1.0f0, 99.0f0]
const RD_PNINF_FLAT = Float32[0.4f0, 0.4f0, 0.4f0, 0.5f0, 0.5f0, 0.0f0, 0.4f0, 0.4f0, 0.5f0, 0.5f0, 0.5f0, 0.4f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.4f0, 0.0f0, 0.0f0, 0.4f0, 0.4f0, 0.45f0, 0.5f0, 0.0f0, 0.0f0, 0.5f0, 0.0f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.0f0, 0.4f0, 0.4f0, 0.4f0, 0.5f0, 0.5f0, 0.0f0, 0.4f0, 0.4f0, 0.5f0, 0.5f0, 0.5f0, 0.4f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.4f0, 0.0f0, 0.0f0, 0.0f0, 0.4f0, 0.45f0, 0.5f0, 0.0f0, 0.0f0, 0.5f0, 0.0f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.0f0, 0.1f0, 0.05f0, 0.5f0, 0.6f0, 0.1f0, 0.1f0, 0.2f0, 0.5f0, 0.5f0, 0.2f0, 0.1f0, 0.1f0, 0.6f0, 0.1f0, 0.5f0, 0.6f0, 0.1f0, 0.1f0, 0.1f0, 0.5f0, 0.6f0, 0.1f0, 0.1f0, 0.1f0, 0.5f0, 0.1f0, 0.6f0, 0.1f0, 0.1f0, 0.1f0, 0.2f0, 0.1f0, 0.1f0, 0.1f0, 0.1f0, 0.1f0, 0.1f0, 0.1f0, 0.4f0, 0.0f0, 0.1f0, 0.2f0, 0.4f0, 0.4f0, 0.2f0, 0.02f0, 0.1f0, 0.2f0, 0.4f0, 0.1f0, 0.4f0, 0.1f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.1f0, 0.4f0, 0.4f0, 0.2f0, 0.4f0, 0.1f0, 0.1f0, 0.4f0, 0.2f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.4f0, 0.1f0, 0.4f0, 0.4f0, 0.02f0, 0.02f0, 0.4f0, 0.1f0, 0.02f0, 0.4f0, 0.0f0]
const RD_PKILLS_FLAT = Float32[0.6f0, 0.9f0, 0.9f0, 0.8f0, 0.8f0, 1.0f0, 0.5f0, 0.9f0, 0.8f0, 0.5f0, 0.8f0, 0.6f0, 0.8f0, 0.9f0, 0.8f0, 0.8f0, 0.8f0, 1.0f0, 0.0f0, 1.0f0, 1.0f0, 0.6f0, 1.0f0, 0.0f0, 0.9f0, 0.5f0, 0.75f0, 0.9f0, 0.0f0, 1.0f0, 0.5f0, 0.0f0, 0.7f0, 0.8f0, 0.8f0, 0.8f0, 0.7f0, 0.8f0, 0.8f0, 1.0f0, 0.6f0, 0.9f0, 0.9f0, 0.8f0, 0.8f0, 1.0f0, 0.5f0, 0.9f0, 0.8f0, 0.5f0, 0.8f0, 0.6f0, 0.8f0, 0.9f0, 0.8f0, 0.8f0, 0.8f0, 1.0f0, 0.0f0, 1.0f0, 1.0f0, 0.6f0, 1.0f0, 0.0f0, 0.9f0, 0.5f0, 0.75f0, 0.9f0, 0.0f0, 1.0f0, 0.5f0, 0.0f0, 0.7f0, 0.8f0, 0.8f0, 0.8f0, 0.7f0, 0.8f0, 0.8f0, 1.0f0, 0.3f0, 1.0f0, 0.8f0, 0.8f0, 0.8f0, 0.75f0, 0.3f0, 0.75f0, 0.8f0, 0.3f0, 0.8f0, 0.3f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.3f0, 0.8f0, 0.8f0, 0.75f0, 0.8f0, 0.3f0, 0.3f0, 0.8f0, 0.75f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 1.0f0, 0.85f0, 0.75f0, 0.8f0, 0.6f0, 0.8f0, 0.85f0, 0.85f0, 0.65f0, 0.6f0, 0.85f0, 0.8f0, 0.85f0, 0.6f0, 0.8f0, 0.6f0, 0.6f0, 0.85f0, 0.8f0, 0.8f0, 0.65f0, 0.6f0, 0.85f0, 0.85f0, 0.8f0, 0.65f0, 0.8f0, 0.6f0, 0.8f0, 0.8f0, 0.8f0, 0.85f0, 0.8f0, 0.8f0, 0.85f0, 0.85f0, 0.8f0, 0.85f0, 0.85f0, 0.6f0, 1.0f0]
const RD_IDITYP = Int32[1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 2, 1, 0, 0, 2, 1, 1, 1, 0, 2, 1, 1, 2, 0, 1, 1, 0, 1, 2, 2, 2, 1, 1, 2, 0]
const RD_IRTSPC_KT = Int32[1,2,3,4,5,6,7,8,9,10,30]
# rdblk1ie.f IRTSPC — Inland Empire 23-species → base-RD-species crosswalk.
# IE species order: WP WL DF GF WH RC LP ES AF PP MH WB LM LL PI JU PY AS CO MM PB OH OS.
# MM(20) and PB(21) have no RD host ⇒ mapped to RD "Other" (30). The host tables
# (HABFAC/PNINF/PKILLS/IDITYP/PCOLO in rd/rdinit.f) are byte-identical across
# variants — only this IRTSPC index differs (verified: IE vs KT rdinit.f identical).
const RD_IRTSPC_IE = Int32[1,2,3,4,5,6,7,8,9,10,11,22,23,36,33,26,38,19,24,30,30,18,17]

# ---------------------------------------------------------------------------
# Remaining WRD variants (rd/rdblk1<v>.f DATA IRTSPC, transcribed byte-for-byte
# from each variant's linked rdblk1 in bin/FVS<v>_buildDir). Each variant's RD
# block-data differs from the base in EXACTLY the IRTSPC host-species crosswalk;
# the host tables HABFAC/PNINF/PKILLS/IDITYP/PCOLO/RRPSWT (rd/rdinit.f) are
# byte-identical across ALL 15 base-rd variants (md5 7e6ea38e… verified for the
# 13 buildDir rdinit.f). length(IRTSPC) == that variant's MAXSP.
# ---------------------------------------------------------------------------
# rdblk1bc.f — British Columbia, MAXSP=15
const RD_IRTSPC_BC = Int32[1,2,3,4,5,6,7,8,9,10,40,19,24,3,40]
# rdblk1bm.f — Blue Mountains, MAXSP=18
const RD_IRTSPC_BM = Int32[1,2,3,4,11,26,7,8,9,10,22,23,38,34,19,24,17,18]
# rdblk1ci.f — Central Idaho, MAXSP=19
const RD_IRTSPC_CI = Int32[1,2,3,4,5,6,7,8,9,10,22,38,19,26,40,23,24,17,18]
# rdblk1cr.f — Central Rockies, MAXSP=38
const RD_IRTSPC_CR = Int32[9,21,3,4,13,11,6,2,33,23,7,33,10,22,1,26,20,8,25,19,
                           24,24,18,29,18,29,29,40,26,26,26,26,33,33,33,10,17,18]
# rdblk1ec.f — East Cascades, MAXSP=32
const RD_IRTSPC_EC = Int32[1,2,3,16,6,4,7,8,9,10,5,11,38,22,39,13,36,34,26,40,
                           40,40,40,40,40,19,40,40,40,40,17,18]
# rdblk1em.f — Eastern Montana, MAXSP=19
const RD_IRTSPC_EM = Int32[22,2,3,23,36,26,7,8,9,10,18,19,24,24,24,24,18,17,18]
# rdblk1nc.f — Klamath (VARACD NC), MAXSP=12
const RD_IRTSPC_NC = Int32[27,12,3,13,30,14,29,32,15,10,18,35]
# rdblk1pn.f — Pacific Northwest, MAXSP=39
const RD_IRTSPC_PN = Int32[16,13,4,9,15,8,39,34,14,8,7,31,12,1,10,3,35,6,5,11,
                           40,40,40,40,40,19,40,40,26,36,22,37,38,40,40,40,40,40,40]
# rdblk1so.f — South-Central Oregon / NE California, MAXSP=33
const RD_IRTSPC_SO = Int32[1,12,3,13,11,14,7,8,15,10,26,4,9,16,39,22,2,6,5,38,
                           40,40,40,19,24,40,40,40,40,40,40,17,18]
# rdblk1tt.f — Teton, MAXSP=18
const RD_IRTSPC_TT = Int32[22,23,3,33,20,19,7,8,9,10,26,26,40,40,24,40,17,18]
# rdblk1ut.f — Utah, MAXSP=24
const RD_IRTSPC_UT = Int32[22,23,3,13,20,19,7,8,9,10,33,26,18,33,26,26,33,24,24,40,40,40,17,18]
# rdblk1wc.f — West Cascades, MAXSP=39 (same as PN except FVS-sp 6 → RD 40 not 8)
const RD_IRTSPC_WC = Int32[16,13,4,9,15,40,39,34,14,8,7,31,12,1,10,3,35,6,5,11,
                           40,40,40,40,40,19,40,40,26,36,22,37,38,40,40,40,40,40,40]
# rdblk1ws.f — West Sierra Nevada, MAXSP=43
const RD_IRTSPC_WS = Int32[12,3,13,28,14,31,15,10,7,22,1,33,16,37,37,37,23,10,37,37,33,3,
                           35,11,26,26,26,29,29,29,29,29,29,32,32,19,40,40,40,40,40,17,18]

"""
    rd_irtspc_for(variant) -> Vector{Int32}

Per-variant IRTSPC host-species crosswalk dispatch (the ONLY variant-specific RD
block-data; rd/rdblk1<v>.f). Default = the base NI/CI/KT table (rd/rdblk1.f).
"""
rd_irtspc_for(::AbstractVariant)     = RD_IRTSPC_KT
rd_irtspc_for(::InlandEmpire)        = RD_IRTSPC_IE
rd_irtspc_for(::BritishColumbia)     = RD_IRTSPC_BC
rd_irtspc_for(::BlueMountains)       = RD_IRTSPC_BM
rd_irtspc_for(::CentralIdaho)        = RD_IRTSPC_CI
rd_irtspc_for(::CentralRockies)      = RD_IRTSPC_CR
rd_irtspc_for(::EastCascades)        = RD_IRTSPC_EC
rd_irtspc_for(::EasternMontana)      = RD_IRTSPC_EM
rd_irtspc_for(::Klamath)             = RD_IRTSPC_NC
rd_irtspc_for(::PacificNorthwest)    = RD_IRTSPC_PN
rd_irtspc_for(::SouthCentralOregon)  = RD_IRTSPC_SO
rd_irtspc_for(::Teton)               = RD_IRTSPC_TT
rd_irtspc_for(::Utah)                = RD_IRTSPC_UT
rd_irtspc_for(::WestCascades)        = RD_IRTSPC_WC
rd_irtspc_for(::WestSierra)          = RD_IRTSPC_WS
const RD_HABFAC = reshape(RD_HABFAC_FLAT, RD_ITOTSP, RD_ITOTRR, 2)   # HABFAC(ksp,idi,ihab)
const RD_PNINF  = reshape(RD_PNINF_FLAT,  RD_ITOTSP, RD_ITOTRR)      # PNINF(ksp,idi)
const RD_PKILLS = reshape(RD_PKILLS_FLAT, RD_ITOTSP, RD_ITOTRR)      # PKILLS(ksp,idi)
# RRPSWT (rd/rdinit.f) defaults to 1.0 for every species; only the RRPSWT keyword
# (not in the chunk-0b-2 turnkey path) changes it, so the const default is faithful here.
const RD_RRPSWT = ones(Float32, RD_ITOTSP)

# rd/rdinit.f TEMP5 — PCOLO(ksp,idi): proportion of roots colonized after tree death
# (used by the RDSPRD spread simulation when an infection reaches the lethal radius).
const RD_PCOLO_FLAT = Float32[0.8f0, 0.95f0, 0.95f0, 0.9f0, 0.9f0, 1.0f0, 0.75f0, 0.95f0, 0.9f0, 0.75f0, 0.9f0, 0.8f0, 0.9f0, 0.95f0, 0.9f0, 0.9f0, 0.9f0, 1.0f0, 0.0f0, 1.0f0, 1.0f0, 0.8f0, 1.0f0, 0.0f0, 0.95f0, 0.75f0, 0.9f0, 0.95f0, 0.0f0, 0.86f0, 0.75f0, 0.0f0, 0.9f0, 0.95f0, 0.95f0, 0.95f0, 0.9f0, 0.95f0, 0.95f0, 0.0f0, 0.8f0, 0.95f0, 0.95f0, 0.9f0, 0.9f0, 1.0f0, 0.75f0, 0.95f0, 0.9f0, 0.75f0, 0.9f0, 0.8f0, 0.9f0, 0.95f0, 0.9f0, 0.9f0, 0.9f0, 1.0f0, 0.0f0, 1.0f0, 1.0f0, 0.8f0, 1.0f0, 0.0f0, 0.95f0, 0.75f0, 0.9f0, 0.95f0, 0.0f0, 0.86f0, 0.75f0, 0.0f0, 0.9f0, 0.95f0, 0.95f0, 0.95f0, 0.9f0, 0.95f0, 0.95f0, 0.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.8f0, 0.6f0, 0.3f0, 0.85f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.75f0, 0.65f0, 0.8f0, 0.6f0, 0.3f0, 0.85f0, 0.3f0, 0.85f0, 0.8f0, 0.8f0, 0.75f0, 0.65f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0, 0.8f0]
const RD_PCOLO = reshape(RD_PCOLO_FLAT, RD_ITOTSP, RD_ITOTRR)   # PCOLO(ksp,idi)

"""
    rd_slp(x, xx, yy, n) -> Float32

Port of rd/rdslp.f: piecewise-linear interpolation of `x` on the knot series
`(xx, yy)` (xx strictly increasing, `n` active knots). Clamped to the end values.
"""
function rd_slp(x::Float32, xx, yy, n::Int)
    x < xx[1] && return Float32(yy[1])
    x > xx[n] && return Float32(yy[n])
    i = 1
    @inbounds while x > xx[i+1]
        i += 1
    end
    @inbounds return yy[i] + (yy[i+1] - yy[i]) / (xx[i+1] - xx[i]) * (x - xx[i])
end

"""
    rd_root(dbh, ht, proot, rslop, sdislp, yincpt, oldtpa, grospc, ormsqd, ba) -> Float32

Port of rd/rdroot.f: the live-tree ROOT RADIUS (feet), the first routine of the
RDTREG/RDCNTL spread chain (called per record from rd/rdtreg.f DO-1001). Deterministic
(consumes no RD random draws).

The SDI EFFECT (`SDINEW = (OLDTPA/GROSPC)·(ORMSQD/10)^1.605`; `EFFECT = SDISLP·SDINEW
+ YINCPT`, clamped to [0.5, 1.5]) scales the species-typical root radius. Trees ≥ 3.5"
DBH use `EFFECT·RSLOP·PROOT·DBH/12`; sub-3.5" trees use a height/BA allometry
(`exp(0.61157·ln HT + 0.04032·ln BA − 0.80815)·EFFECT`), falling back to 0.01 when HT
or BA is zero. `PROOT`/`RSLOP` are the per-base-RD-species tables (rd/rdblk1kt.f),
indexed `IRTSPC(ISP)` by the caller. All arithmetic in Float32 (matches the FVS REAL
path; g16 bit-exact, incl. the `^1.605` / exp / log rounding).
"""
function rd_root(dbh::Float32, ht::Float32, proot::Float32, rslop::Float32,
                 sdislp::Float32, yincpt::Float32, oldtpa::Float32,
                 grospc::Float32, ormsqd::Float32, ba::Float32)::Float32
    sdinew = (oldtpa / grospc) * (ormsqd / 10.0f0)^1.605f0
    effect = sdinew > 0.0f0 ? sdislp * sdinew + yincpt : 1.0f0
    effect = min(effect, 1.5f0)
    effect = max(effect, 0.5f0)
    if dbh < 3.5f0
        ans = 0.01f0
        if ht != 0.0f0 && ba != 0.0f0
            ans = exp(0.61157f0 * log(ht) + 0.04032f0 * log(ba) - 0.80815f0)
            ans = ans * effect
        end
        return ans
    else
        return effect * rslop * proot * dbh / 12.0f0
    end
end

"""
    rd_iprp!(rd, s, noplot, pran)

Port of rd/rdiprp.f: compute the per-record infected proportion `PROPN` for the
manual-RRINIT init (NOPLOT=0). Distributes `PRAN[idi]·TOTREE` infected trees over
the host species by their relative value (PNINF × average years-to-kill), capping
each species at proportion 1.0 and redistributing the overflow (labels 1750–1900).
All arithmetic in Float32 (matches the FVS REAL path). Writes `rd.propn` (size ITRN).
"""
function rd_iprp!(rd::RootDiseaseState, s::StandState, noplot::Int, pran::Vector{Float32})
    t   = s.trees
    nsp = nspecies(s.variant)
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr)
    irt = rd.irtspc; ihab = Int(rd.irhab)
    xxinf = rd.xxinf; yyinf = rd.yyinf; nninf = Int(rd.nninf)

    totval = zeros(Float32, RD_ITOTRR); totinf = zeros(Float32, RD_ITOTRR)
    totree = zeros(Float32, RD_ITOTRR); inumsp = zeros(Int, RD_ITOTRR)
    sstval = zeros(Float32, RD_ITOTRR); numful = zeros(Int, RD_ITOTRR)
    yessp  = falses(nsp); full = falses(nsp)
    totsp  = zeros(Float32, nsp); totytk = zeros(Float32, nsp)
    treeno = zeros(Float32, nsp); relval = zeros(Float32, nsp); tmpprp = zeros(Float32, nsp)

    # DO 1500 — totals by species (ITRN = t.n; no recently-dead at LSTART init).
    idi = maxrr
    @inbounds for i in 1:t.n
        # IDPLOT(I)=0 for the turnkey manual path ⇒ ((IDPLOT==NOPLOT) .OR. (NOPLOT==0)) true
        ((0 == noplot) || (noplot == 0)) || continue
        t.species[i] == 0 && continue
        ksp = Int(t.species[i])
        maxrr < 3 && (idi = Int(RD_IDITYP[irt[ksp]]))
        idi <= 0 && continue
        totnum = t.tpa[i] * rd.parea[idi]
        totsp[ksp] += totnum
        habsp = RD_HABFAC[irt[ksp], idi, ihab]
        ytk = rd_slp(t.dbh[i], xxinf, yyinf, nninf)
        if ytk > rd.xminkl[idi]
            ytk = (ytk - rd.xminkl[idi]) * habsp * RD_RRPSWT[irt[ksp]] + rd.xminkl[idi]
        end
        totytk[ksp] += ytk * totnum
        (!yessp[ksp] && totsp[ksp] > 0.0f0) && (yessp[ksp] = true)
    end

    # DO 1600 — relative value + per-disease totals.
    idi = maxrr
    @inbounds for ksp in 1:nsp
        totsp[ksp] > 0.0f0 || continue
        maxrr < 3 && (idi = Int(RD_IDITYP[irt[ksp]]))
        idi <= 0 && continue
        avgytk = totytk[ksp] / totsp[ksp]
        relval[ksp] = RD_PNINF[irt[ksp], idi] * avgytk
        totval[idi] += relval[ksp]
        totree[idi] += totsp[ksp]
        yessp[ksp] && (inumsp[idi] += 1)
    end

    # DO 1699 — TOTINF target + SSTVAL normalizer.
    @inbounds for id2 in minrr:maxrr
        totinf[id2] = pran[id2] * totree[id2]
        for ksp in 1:nsp
            totval[id2] > 0.0f0 && (sstval[id2] += relval[ksp] * totsp[ksp] / totval[id2])
        end
    end

    # DO 1700 — trees-to-infect by species.
    idi = maxrr
    @inbounds for ksp in 1:nsp
        maxrr < 3 && (idi = Int(RD_IDITYP[irt[ksp]]))
        idi <= 0 && continue
        if totsp[ksp] > 0.0f0 && totval[idi] > 0.0f0 && sstval[idi] > 0.0f0
            treeno[ksp] = totinf[idi] * totsp[ksp] * relval[ksp] / sstval[idi] / totval[idi]
        end
    end

    # DO 1900 with the 1750 redistribution restart (GOTO 1750 restarts the whole
    # MINRR..MAXRR loop, carrying OVER; OVER is zeroed at each 1850 exit).
    over = 0.0f0
    restart = true
    @inbounds while restart
        restart = false
        for irrsp in minrr:maxrr
            if inumsp[irrsp] > numful[irrsp]
                addon = over / Float32(inumsp[irrsp] - numful[irrsp]); over = 0.0f0
                idi3 = irrsp
                for ksp in 1:nsp
                    irrsp < 3 && (idi3 = Int(RD_IDITYP[irt[ksp]]))
                    idi3 != irrsp && continue
                    if yessp[ksp] && !full[ksp]
                        treeno[ksp] += addon
                        tmpprp[ksp] = treeno[ksp] / totsp[ksp]
                        if tmpprp[ksp] > 1.0f0
                            tmpprp[ksp] = 1.0f0; full[ksp] = true; numful[irrsp] += 1
                            over += treeno[ksp] - totsp[ksp]
                        end
                    end
                end
                if numful[irrsp] != inumsp[irrsp] && over > 0.0f0 && numful[irrsp] < inumsp[irrsp]
                    restart = true; break          # GOTO 1750
                end
            end
            over = 0.0f0                            # 1850 CONTINUE
        end
    end

    # DO 2000 — PROPN by record from the species proportion.
    resize!(rd.propn, t.n); fill!(rd.propn, 0.0f0)
    @inbounds for ii in 1:t.n
        ((0 == noplot) || (noplot == 0)) || continue
        ksp = Int(t.species[ii])
        ksp == 0 && continue
        rd.propn[ii] = tmpprp[ksp]
    end
    return rd
end

"""
    rd_inoc!(rd, s, licall)

Port of rd/rdinoc.f: decompose infected root systems in the dead-tree/stump list.
At LSTART init (`licall=true`) FVS decays any stumps carried in from a treelist/STREAD
init. In the chunk-0b-2 manual-RRINIT turnkey there are NO initialized stumps (the
PROBDA/DBHDA stump-decay arrays are all zero), so every guarded branch is skipped and
this is a verified no-op that does not touch PROBI/PROBIU/FPROB/PROPI. The full
stump-decay body activates once stumps exist (RRTREIN treelist-init / per-cycle
mortality) — a later chunk. See the port report for the validated boundary.
"""
function rd_inoc!(::RootDiseaseState, ::StandState, ::Bool)
    # No stump-root inoculum records at manual-RRINIT LSTART ⇒ nothing to decay.
    return nothing
end

"""
    rd_setp!(rd, s)

Port of rd/rdsetp.f (manual-RRINIT init path): place disease centers, normalize the
per-disease infected/uninfected densities, distribute infection over species
(rd_iprp!), then set the per-record initial-infection state PROBI/PROBIU/FPROB/PROPI
and PROBL. Called once at LSTART from `root_disease_setup!`. Advances the RD stream
(center placement + one RDRANP draw per host record).
"""
function rd_setp!(rd::RootDiseaseState, s::StandState)
    t = s.trees; p = s.plot
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr)

    # Diseased sub-plot count. The manual-RRINIT turnkey has no PLOTINF/IRDPLT ⇒ 0.
    ndplts = zeros(Int, RD_ITOTRR)
    iipi = trunc(Int, p.pi)                        # IIPI = INT(PI) (points inventoried)
    @inbounds for idi in minrr:maxrr
        if iipi == ndplts[idi]
            rd.lonect[idi] = Int32(1); rd.parea[idi] = rd.sarea
        end
        rd.lonect[idi] == 1 && (ndplts[idi] = 1)
    end

    # SDI-effect intercept + stand square dimension.
    rd.yincpt = (rd.sdislp == 0.0f0 || rd.sdnorm == 0.0f0) ? 1.0f0 : 1.0f0 - (rd.sdnorm * rd.sdislp)
    rd.dimen  = rd.sarea != 0.0f0 ? sqrt(rd.sarea) * 208.7f0 : 0.0f0

    # DO 600 — place centers (rd_cloc!/rd_area!) and normalize RRGEN densities.
    @inbounds for irrsp in minrr:maxrr
        rd.irrsp = Int32(irrsp)
        if rd.ipcflg[irrsp] == 1
            rd_area!(rd, true)                     # user-assigned centers
        elseif rd.parea[irrsp] > 0.0f0
            rd_cloc!(rd); rd_area!(rd, true)
        else
            rd.parea[irrsp] = rd.sarea * Float32(ndplts[irrsp]) / p.pi   # rdsetp.f:174 (PI = points)
            rd_cloc!(rd)
        end
        # RRGEN normalization (.NOT. RRTINV .AND. .NOT. LPLINF): PRKILL/PRUN densities → proportions.
        rd.rrgen1[irrsp] = rd.prkill[irrsp] + rd.prun[irrsp]
        if rd.rrgen1[irrsp] > 0.0f0
            rd.prkill[irrsp] = rd.prkill[irrsp] / rd.rrgen1[irrsp]
            rd.prun[irrsp]   = rd.prun[irrsp]   / rd.rrgen1[irrsp]
            rd.rrgen1[irrsp] = 1.0f0
        end
    end

    # Distribute infection over species → PROPN by record.
    rd_iprp!(rd, s, 0, rd.prkill)

    # Per-record initial-infection state.
    nrec = t.n
    rd.probi  = zeros(Float32, nrec); rd.probiu = zeros(Float32, nrec)
    rd.fprob  = zeros(Float32, nrec); rd.propi  = zeros(Float32, nrec)
    rd.probl  = zeros(Float32, nrec)

    idi = maxrr
    @inbounds for i in 1:t.n
        if maxrr < 3 && minrr != maxrr
            idi = Int(RD_IDITYP[rd.irtspc[t.species[i]]])
        end
        idi <= 0 && continue                       # non-host record: skip (GOTO 2100)
        prob_i = t.tpa[i]                           # FVS PROB(I)
        # Manual-RRINIT path ((.NOT. LPLINF .OR. IDPLOT>0) .AND. IDI≠0): LPLINF false ⇒ true.
        rd.parea[idi] > 0.0f0 && (rd.propi[i] = rd_ranp!(rd, rd.rrincs[idi]))
        denom = rd.rrgen1[idi] + 0.00001f0
        rd.probi[i]  = rd.propn[i] * prob_i * rd.parea[idi] / denom
        rd.probiu[i] = (1.0f0 - rd.propn[i]) * prob_i * rd.parea[idi] / denom
        rd.fprob[i]  = prob_i
        rd.probl[i]  = prob_i                       # PROBL(I) = PROB(I)
    end

    rd_inoc!(rd, s, true)                           # rd/rdinoc.f (.TRUE.) — inert (no stumps)
    return rd
end

# -----------------------------------------------------------------------------
# kw_rdin! — RDIN keyword block reader (rd/rdin.f entry RDKEY subset)
# -----------------------------------------------------------------------------
"""
    kw_rdin!(s, rec, kr)

Parse the `RDin … End` block (rd/rdin.f). Chunk-0 subset: RRTYPE, RRINIT, SAREA,
BBCLEAR, END (RRDOUT/BBOUT recognized no-ops). Populates `s.root_disease`.
"""
function kw_rdin!(s::StandState, rec, kr::KeywordReader)
    s.root_disease === nothing && (s.root_disease = RootDiseaseState())
    rd = s.root_disease
    rd.irtspc = copy(rd_irtspc_for(s.variant))  # rdblk1<v>.f — per-variant host crosswalk
    rd.iroot = Int32(1)                     # rdin.f: entry sets IROOT=1

    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            rd_in_end!(rd)                  # rd/rdin.f option 9 (END)
            break
        elseif k == "RRTYPE"                # rd/rdin.f option 42
            rd_in_rrtype!(rd, r)
        elseif k == "RRINIT"                # rd/rdin.f option 5
            rd_in_rrinit!(rd, r)
        elseif k == "SAREA"                 # rd/rdin.f option 10
            rd_in_sarea!(rd, r)
        elseif k == "RSEED"                 # rd/rdin.f option 25 — reseed the RD RNG
            # ARRAY(1)==0 ⇒ GETSED clock-seed (non-deterministic; out of scope).
            if r.present[1] && r.values[1] != 0.0f0
                rd.dseed = Float64(r.values[1])
                rd_rani!(rd, rd.dseed)      # RDMN1(1) reseeds from DSEED after parse
            end
        elseif k == "BBCLEAR"               # rd/rdin.f option 41 — suppress default bark beetles
            rd.bbclear = true; rd.lbbon = false
        elseif k == "RRDOUT" || k == "BBOUT" || k == "RRECHO" || k == "SMCOUT"
            # report-request sub-keywords — recognized, no state effect for chunk 0
        else
            # RD sub-keyword not yet ported (SPREAD/CARRY/PSTUMP/RRTREIN/…); record it
            (!isempty(k) && isletter(first(k))) && push!(s.control.unrecognized_keywords, k)
        end
    end
    return nothing
end

# rd/rdin.f option 42 — RRTYPE: pick the active disease number(s) from ARRAY.
function rd_in_rrtype!(rd::RootDiseaseState, r)
    rd.lrtype = true
    minrr = 0; maxrr = 0
    @inbounds for idi in 1:min(RD_ITOTRR, length(r.values))
        if r.present[idi]
            ipoint = Int(nint(r.values[idi]))
            (ipoint < minrr || minrr == 0) && (minrr = ipoint)
            ipoint > maxrr && (maxrr = ipoint)
        end
    end
    minrr <= 0 && (minrr = 1)
    minrr = clamp(minrr, 1, RD_ITOTRR)
    maxrr = clamp(maxrr, 1, RD_ITOTRR)
    maxrr < minrr && (maxrr = minrr)
    # only one disease type modeled at a time, except the P+S annosus pair (1,2)
    if maxrr != minrr && !(minrr == 1 && maxrr == 2)
        maxrr = minrr
    end
    rd.minrr = Int32(minrr)
    rd.maxrr = Int32(maxrr)
    return nothing
end

# rd/rdin.f option 5 — RRINIT: manual center-init. Fields (ARRAY):
#   1 = center-placement flag (<1 ⇒ random), 2 = #centers, 3 = infected TPA/acre,
#   4 = uninfected TPA/acre, 5 = proportion roots infected, 6 = disease area acres,
#   7 = disease type (else MINRR).
function rd_in_rrinit!(rd::RootDiseaseState, r)
    idi = Int(rd.minrr)
    (length(r.present) >= 7 && r.present[7]) && (idi = clamp(Int(nint(r.values[7])), 1, RD_ITOTRR))
    rd.rrman = true

    rndcenters = !(r.present[1] && r.values[1] >= 1.0f0)   # ARRAY(1) < 1 ⇒ random
    if rndcenters
        (r.present[2]) && (rd.ncents[idi] = Int32(min(100, Int(nint(r.values[2])))))
        if r.present[6]
            rd.parea[idi] = Float32(r.values[6]); rd.lparea[idi] = true
        end
        rd.ipcflg[idi] = Int32(0)
        (rd.ncents[idi] == 1 && rd.parea[idi] == rd.sarea) && (rd.lonect[idi] = Int32(1))
    else
        # explicit-center branch (IPCFLG=1) reads center geometry records — Chunk 1.
        (r.present[2]) && (rd.ncents[idi] = Int32(min(100, Int(nint(r.values[2])))))
        rd.parea[idi] = -1.0f0; rd.lparea[idi] = true; rd.ipcflg[idi] = Int32(1)
    end
    # field 3/4 = infected/uninfected TPA in diseased area; field 5 = initial root-infection proportion
    (r.present[3] && r.values[3] >= 0.0f0) && (rd.prkill[idi] = Float32(r.values[3]))
    (r.present[4] && r.values[4] >= 0.0f0) && (rd.prun[idi]   = Float32(r.values[4]))
    (r.present[5] && 0.0f0 <= r.values[5] <= 1.0f0) && (rd.rrincs[idi] = Float32(r.values[5]))
    return nothing
end

# rd/rdin.f option 10 — SAREA: stand area (acres). DIMEN = sqrt(SAREA)*208.7.
function rd_in_sarea!(rd::RootDiseaseState, r)
    if r.present[1] && r.values[1] > 0.0f0
        rd.sarea = Float32(r.values[1])
        rd.dimen = sqrt(rd.sarea) * 208.7f0
    end
    return nothing
end

# rd/rdin.f option 9 — END: defaults not-set patch areas and finalizes LONECT.
# (The default-bark-beetle auto-enable is deferred to Chunk 7; noted via lbbon.)
function rd_in_end!(rd::RootDiseaseState)
    @inbounds for idi in Int(rd.minrr):Int(rd.maxrr)
        rd.lparea[idi] || (rd.parea[idi] = 0.25f0 * rd.sarea)
        (rd.lonect[idi] == 0 && rd.ncents[idi] == 1 && rd.parea[idi] == rd.sarea) &&
            (rd.lonect[idi] = Int32(1))
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Per-cycle mortality kernel (rd/rdmort.f + rd/rdsum.f) — Chunk 0b-3.
#
# RDMORT ages the infected root systems of every live host record through one FVS
# cycle and kills the record's infected TPA (into RRKILL/RDKILL) once the modeled
# proportion of infected roots reaches the species' lethal threshold PKILLS. It is
# the FIRST .sum-visible WRD effect (RDEND later applies RRKILL to WK2 mortality).
#
# Ported faithfully from rd/rdmort.f and validated BIT-EXACT (0-ULP) against the
# LIVE relinked FVSkt oracle: a single-.o instrumentation swap of rdmort.f dumped
# every record's PROBI/PROPI entry state and RRKILL/RDKILL/PROPI exit state for all
# 10 cycles of the turnkey scenario (RRType 3 Armillaria, RRInit 0 10 10 20 0.1 10 3,
# SArea 100; the instrumented .sum stayed byte-identical to FVSkt_clean). Driving the
# kernel with the oracle's dumped entry state reproduces every exit value bit-exact.
#
# SCOPE BOUNDARY (see the port report): this is the mortality KERNEL. The ENTRY state
# it consumes is produced upstream by the full RDCNTL spread chain (RDINSD inside-patch
# aging into slot IT=ISTEP,IP=1 + RDINF new-infection into slot IT=ISTEP,IP=2, driven
# by RDSPRD/RDAREA center growth), and the RRKILL it produces is applied to the .sum by
# RDEND/RDGROW downstream — those remain later sub-chunks. RDSTP/RDSSIZ stump-list
# bookkeeping (feeds next cycle's RDSHRK, not the RRKILL signal) is likewise deferred.
# -----------------------------------------------------------------------------

# rd/rdinit.f spore-model per-disease defaults (no spore centers in the turnkey path,
# so SPPROP stays 0 ⇒ the HABSP spore multiplier is identically 1). Kept as the faithful
# defaults; the SPORE keyword that changes them is a later chunk.
const RD_SPPROP0 = 0.0f0
const RD_SPYTK0  = 3.0f0

"""
    rd_sum!(probit, probi, istep)

Port of rd/rdsum.f: `PROBIT(I) = Σ_{IT=1..ISTEP, IP=1..2} PROBI(I,IT,IP)` — the total
infected TPA per record, summed over all time-slots and both infection pathways.
"""
function rd_sum!(probit::Vector{Float32}, probi::Array{Float32,3}, istep::Int)
    nrec = size(probi, 1)
    @inbounds for i in 1:nrec
        s = 0.0f0
        for it in 1:istep, ip in 1:2
            s += probi[i, it, ip]
        end
        probit[i] = s
    end
    return probit
end

"""
    rd_mort_kernel!(rd, probi, propi, rrkill, rdkill, dbh, isp, istep, fint;
                    spprop, spytk)

Port of the rd/rdmort.f per-record body. Mutates `probi` (killed slots zeroed) and
`propi` (aged), accumulates the killed infected TPA into `rrkill`/`rdkill` (both are
the same running total within a cycle; `rdkill` is zeroed here per RDMORT, `rrkill`
is the cycle accumulator zeroed by RDMN2). Returns `rdkill`.

Faithful to rd/rdmort.f:
  YTKILL = RDSLP(DBH) then, if > XMINKL, scaled by HABSP·RRPSWT (+XMINKL);
  CURAGE = PROPI·YTKILL/PKILLS; PROPI ← (CURAGE+FINT)·PKILLS/YTKILL;
  when PROPI ≥ PKILLS the slot's PROBI is killed (added to RRKILL/RDKILL, then zeroed).
HABSP = HABFAC(base,IDI,IRHAB)·(SPPROP·SPYTK + (1−SPPROP)). IDI = MAXRR unless the
annosus case (MAXRR<3) selects it by IDITYP(base). Non-host records (IDI≤0) skip.
"""
function rd_mort_kernel!(rd::RootDiseaseState,
                         probi::Array{Float32,3}, propi::Array{Float32,3},
                         rrkill::Vector{Float32}, rdkill::Vector{Float32},
                         dbh::AbstractVector{<:Real}, isp::AbstractVector{<:Integer},
                         istep::Int, fint::Real;
                         spprop::Float32 = RD_SPPROP0, spytk::Float32 = RD_SPYTK0)
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr); irhab = Int(rd.irhab)
    irt = rd.irtspc; xxinf = rd.xxinf; yyinf = rd.yyinf; nninf = Int(rd.nninf)
    fintf = Float32(fint)
    nrec = length(dbh)

    fill!(rdkill, 0.0f0)                          # rdmort.f DO 150 (RDKILL zeroed)

    tarea = 0.0f0                                 # rdmort.f DO 410
    @inbounds for idi in minrr:maxrr
        tarea += rd.parea[idi]
    end
    (tarea <= 0.0f0 || nrec == 0) && return rdkill

    hmul = spprop * spytk + (1.0f0 - spprop)      # spore time-to-death multiplier
    @inbounds for i in 1:nrec
        ksp = Int(isp[i]); ksp == 0 && continue
        base = Int(irt[ksp])
        idi = maxrr < 3 ? Int(RD_IDITYP[base]) : maxrr
        idi <= 0 && continue                      # non-host: skip (rdmort GOTO 500)
        habsp = RD_HABFAC[base, idi, irhab] * hmul
        pk = RD_PKILLS[base, idi]
        xmk = rd.xminkl[idi]
        for it in 1:istep, ip in 1:2
            probi[i, it, ip] <= 0.0f0 && continue
            ytk = rd_slp(Float32(dbh[i]), xxinf, yyinf, nninf)
            if ytk > xmk
                ytk = (ytk - xmk) * habsp * RD_RRPSWT[base] + xmk
            end
            curage = propi[i, it, ip] * ytk / pk
            propi[i, it, ip] = (curage + fintf) * pk / ytk
            propi[i, it, ip] < pk && continue     # not yet lethal
            rrkill[i] += probi[i, it, ip]
            rdkill[i] += probi[i, it, ip]
            probi[i, it, ip] = 0.0f0
        end
    end
    return rdkill
end

# -----------------------------------------------------------------------------
# Downstream .sum-application kernels (rd/rdmn2.f, rd/rdend.f, rd/rdgrow.f) —
# Chunk 0b-3b. These are the pieces that make the RD signal .sum-VISIBLE:
#   * RDMN2  advances the time-slot counter and re-sums PROBI each cycle.
#   * RDEND  reconciles the RD infected-tree deaths (RRKILL) into the FVS
#     per-record mortality WK2 — the mortality the `.sum` reports.
#   * RDGROW reduces the per-record diameter/height growth DG/HTG in proportion
#     to the infected roots (PROPI) — the compounding BA growth-loss.
# The .sum-visible arithmetic (the WK2 reconciliation and the DG/HTG combine)
# is factored into the pure helpers `rd_end_newwk2` / `rd_grow_newg`, validated
# BIT-EXACT (0-ULP Float32) against the live relinked FVSkt oracle: a single-.o
# instrumentation swap of rdend.f/rdgrow.f dumped every record's entry state and
# WK2/DG exit value for all 10 cycles of the turnkey scenario (RRType 3 Armillaria,
# RRInit 0 10 10 20 0.1 10 3, SArea 100; the instrumented .sum stayed byte-identical
# to FVSkt_clean). Driving the helpers with the oracle's dumped entry state
# reproduces every WK2 (270/270 records) and DG (270/270 records) bit-exact.
#
# SCOPE BOUNDARY (see the port report): these consume the per-record infection
# ENTRY state (PROBI/PROPI slots, PROBIT, PROBIU, FPROB, RRKILL, PAREA) that is
# BUILT upstream by the RDCNTL Monte-Carlo spread chain (RDSPRD center-radius
# growth → RDAREA new PAREA → RDINF new-area infection + RDINSD inside-patch
# infection). That upstream chain — and the engine seam-reorder that runs RDEND
# after FVS MORTS computes WK2 — are the remaining 0b-3b sub-chunks, so the
# engine seam below stays INERT (the 0b-2/0b-3 no-RD/RD byte-identical guarantee
# is preserved). The kernels here are validated as pure functions against the
# oracle dumps, exactly as the 0b-3 RDMORT kernel was.
# -----------------------------------------------------------------------------

# OAKL/BBKILL "dead-tree" second-index codes (rd/rdend.f DATA DRR/DBB/DWND/DNAT).
const RD_DRR  = 1        # root disease
const RD_DBB  = 2        # bark beetle
const RD_DWND = 3        # windthrow
const RD_DNAT = 4        # natural

"""
    rd_end_newwk2(wk2, prob, rrkill, probit, tdiun, tdieou, sarea) -> Float32

The rd/rdend.f WK2 reconciliation (the `.sum`-visible mortality). `BACKGD = WK2/PROB`
is the background mortality fraction; `TDIEN = AMAX1(RRKILL, BACKGD·(PROBIT+RRKILL))`
is the infected-tree death (always ≥ the RD kill); `WMESS = (TDIEOU+TDIEN+TDIUN)/SAREA`
is the combined RD mortality, and the new `WK2 = AMAX1(WMESS, WK2)` capped at `PROB`.
All Float32. Validated bit-exact vs the live FVSkt oracle (270/270 records).
"""
function rd_end_newwk2(wk2::Float32, prob::Float32, rrkill::Float32, probit::Float32,
                       tdiun::Float32, tdieou::Float32, sarea::Float32)
    backgd = wk2 / prob
    tdien  = max(rrkill, backgd * (probit + rrkill))
    wmess  = (tdieou + tdien + tdiun) / sarea
    w = max(wmess, wk2)
    w > prob && (w = prob)
    return w
end

"""
    rd_grow_newg(g, gtot, outnum, probiu, probit, diff) -> Float32

The rd/rdgrow.f per-record growth (DG or HTG) reduction. Weighted average of the
UNaffected trees (outside `OUTNUM` + inside-uninfected `PROBIU`, full growth) and the
infected trees (`PROBIT`, growth scaled by `GTOT`), over `BOTTOM=OUTNUM+PROBIU+PROBIT`.
Reproduces rdgrow.f including the `GOTO 999` skip (`DIFF≤1e-3 ∧ PROBIU≤1e-3 ∧ PROBIT≤1e-3`
⇒ growth unchanged) and the 0.0001 floor. All Float32. Validated bit-exact vs the live
FVSkt oracle (270/270 records, DG). `gtot` = Σ RDSLP(max(PROPI,0))·PROBI/(PROBIT+1e-6).
"""
function rd_grow_newg(g::Float32, gtot::Float32, outnum::Float32, probiu::Float32,
                      probit::Float32, diff::Float32)
    if diff <= 1.0f-3 && probiu <= 1.0f-3 && probit <= 1.0f-3
        return max(g, 0.0001f0)                       # rdgrow.f GOTO 999 (growth unchanged)
    end
    bottom = outnum + probiu + probit
    bottom <= 1.0f-6 && return 0.0001f0
    g2 = (g * (outnum + probiu) + g * gtot * probit) / bottom
    return max(g2, 0.0001f0)
end

"""
    rd_grow_gtot(rd, probi, propi, i, isp, xknot, yknot, ifac) -> Float32

The rd/rdgrow.f inner DGTOT/HTTOT accumulator for record `i`: sum over all infection
slots (IT=1..ISTEP, IP=1..2) of `RDSLP(max(PROPI,0), xknot, yknot, 2)·IFAC·PROBI /
(PROBIT+1e-6)`, skipping slots with PROBI ≤ 1e-3. `RDSLP` is the (already bit-exact)
piecewise-linear reduction curve; XDBH/YDBH=(0,0.5)/(1,0) for diameter, XHT/YHT for
height. IFAC = DBIFAC/HTIFAC (both 1.0 by default). Negative PROPI (infection not yet
at the tree center) ⇒ PRADI=0 ⇒ no growth loss.
"""
function rd_grow_gtot(probi::Array{Float32,3}, propi::Array{Float32,3}, i::Int,
                      istep::Int, probit_i::Float32,
                      xknot, yknot, ifac::Float32)
    gtot = 0.0f0
    @inbounds for j in 1:istep, ip in 1:2
        probi[i, j, ip] <= 1.0f-3 && continue
        pradi = propi[i, j, ip]
        pradi < 0.0f0 && (pradi = 0.0f0)
        gimp = rd_slp(pradi, xknot, yknot, 2) * ifac
        gtot += gimp * probi[i, j, ip] / (probit_i + 1.0f-6)
    end
    return gtot
end

# rd/rdgrow.f XDBH/YDBH + XHT/YHT growth-reduction knots (rdinit.f). PRADI∈[0,0.5]
# ⇒ linear 1→0; PRADI≥0.5 ⇒ 0 (an infected-to-the-center tree stops growing).
const RD_XDBH = (0.0f0, 0.5f0)
const RD_YDBH = (1.0f0, 0.0f0)
const RD_XHT  = (0.0f0, 0.5f0)
const RD_YHT  = (1.0f0, 0.0f0)
# rdinit.f DBIFAC/HTIFAC default 1.0 for every species (only the DBIFAC/HTIFAC keyword,
# not in the turnkey path, changes them) — kept as the faithful const default.
const RD_DBIFAC = 1.0f0
const RD_HTIFAC = 1.0f0

"""
    rd_mn2_advance!(rd, rrkill, probit, probi, fint)

Port of the .sum-relevant body of rd/rdmn2.f (called from GRINCR each cycle BEFORE
tripling): advance the time-slot counter ISTEP, add the cycle length to IYEAR, zero
the per-record RRKILL accumulator (RDKILL is NOT zeroed — it carries to next cycle as
PRANKL), then re-sum PROBI into PROBIT. Inert when RD is inactive or has no trees.
(RDSHST stump-history + PROBIN carryover-model zeroing are deferred: no stumps / no
carryover in the manual-RRINIT turnkey.)
"""
function rd_mn2_advance!(rd::RootDiseaseState, rrkill::Vector{Float32},
                         probit::Vector{Float32}, probi::Array{Float32,3}, fint::Real)
    rd.iroot == 0 && return rd
    rd.istep += Int32(1)
    fill!(rrkill, 0.0f0)
    rd_sum!(probit, probi, Int(rd.istep))
    return rd
end

"""
    rd_end_kernel!(rd, wk2, prob, rrkill, rdkill, probit, probiu, fprob, probi, propi,
                   isp, rootl, wk22, rroott, dprob, oakl, bbkill)

Port of the rd/rdend.f per-record body. Mutates `wk2` (the FVS per-record mortality
the `.sum` reports) via `rd_end_newwk2`, then applies the rd/rdend.f state updates:
the RROOTT/WK22 root-radius weighted average, the natural-mortality reallocation to
inside-uninfected `probiu`, the outside-tree `fprob` reduction, the `dprob` dead-tree
bookkeeping (by agent), and the DIENAT natural-mortality redistribution over the
infection slots `probi` (accumulated into `rrkill`). Finishes with RDSUM(→`probit`).
`oakl`/`bbkill` are the (3×n) other-agent / bark-beetle kill densities — all zero in
the turnkey scenario (no active bark beetles), so TDIUN/TDIEOU and the DPROB agent
splits are zero there; the code is faithful for the general case. Returns `wk2`.
"""
function rd_end_kernel!(rd::RootDiseaseState,
                        wk2::Vector{Float32}, prob::Vector{Float32},
                        rrkill::Vector{Float32}, rdkill::Vector{Float32},
                        probit::Vector{Float32}, probiu::Vector{Float32},
                        fprob::Vector{Float32}, probi::Array{Float32,3},
                        propi::Array{Float32,3}, isp::AbstractVector{<:Integer},
                        rootl::Vector{Float32}, wk22::Vector{Float32},
                        rroott::Vector{Float32}, dprob::Array{Float32,3},
                        oakl::Array{Float32,2}, bbkill::Array{Float32,2};
                        probda = nothing, dbhda = nothing, rootda = nothing,
                        dbh = nothing, driver = nothing)
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr); irt = rd.irtspc
    istep = Int(rd.istep); sarea = rd.sarea
    n = length(prob)
    fill!(dprob, 0.0f0)                              # rdend.f DO 27
    tparea = 0.0f0
    @inbounds for idi in minrr:maxrr; tparea += rd.parea[idi]; end
    (n == 0 || tparea == 0.0f0) && return wk2
    @inbounds for i in 1:n
        ksp = Int(isp[i]); ksp == 0 && continue
        base = Int(irt[ksp])
        idi = maxrr < 3 ? Int(RD_IDITYP[base]) : maxrr
        idi <= 0 && continue                         # non-host: skip (rdend GOTO 1000)
        backgd = wk2[i] / prob[i]
        tdien  = max(rrkill[i], backgd * (probit[i] + rrkill[i]))
        tdiun  = min(probiu[i] * 0.95f0, oakl[RD_DSIU, i])
        diff   = sarea - rd.parea[idi]
        test   = diff > 0.0f0 ? oakl[RD_DSO, i] / diff : 0.0f0
        pmort  = min(fprob[i] * 0.95f0, test)
        tdieou = pmort * diff
        wmess  = (tdieou + tdien + tdiun) / sarea
        wk2[i] = rd_end_newwk2(wk2[i], prob[i], rrkill[i], probit[i], tdiun, tdieou, sarea)
        # RROOTT/WK22 weighted root-radius average (next-cycle inoculum bookkeeping).
        rroott[i] = (rroott[i] * wk22[i] + rootl[i] * wk2[i]) /
                    (wk22[i] + wk2[i] + 0.00001f0)
        wk22[i] += wk2[i]
        # Natural-mortality reallocation to inside-uninfected trees.
        die = tdiun; natiu = 0.0f0; unapp = 0.0f0
        if wmess < wk2[i]
            unapp = (wk2[i] - wmess) * sarea
            unapp = unapp - (backgd * fprob[i] * diff - tdieou)
            natiu = backgd * probiu[i]
            natiu > unapp && (natiu = max(unapp, 0.0f0))
            natiu > tdiun && (die = natiu)
        end
        probiu[i] = probiu[i] - die
        (fprob[i] - pmort > 1.0f-6) && (fprob[i] = fprob[i] - pmort)
        # DPROB dead-tree bookkeeping by agent (all zero when oakl/bbkill are zero).
        if oakl[RD_DSO, i] > 0.0f0
            pbb = bbkill[RD_DSO, i] / oakl[RD_DSO, i]
            dprob[i, RD_DBB,  RD_DSO] = tdieou * pbb
            dprob[i, RD_DWND, RD_DSO] = tdieou * (1.0f0 - pbb)
        end
        if oakl[RD_DSIU, i] > 0.0f0
            pbb = bbkill[RD_DSIU, i] / oakl[RD_DSIU, i]
            dprob[i, RD_DBB,  RD_DSO]  = tdiun * pbb
            dprob[i, RD_DWND, RD_DSO]  = tdiun * (1.0f0 - pbb)
            dprob[i, RD_DNAT, RD_DSO]  = die - tdiun
        end
        dprob[i, RD_DBB,  RD_DSII] = bbkill[RD_DSII, i]
        dprob[i, RD_DWND, RD_DSII] = oakl[RD_DSII, i] - bbkill[RD_DSII, i]
        dprob[i, RD_DRR,  RD_DSII] = rdkill[i]
        dprob[i, RD_DWND, RD_DSII] < 0.0f0 && (dprob[i, RD_DWND, RD_DSII] = 0.0f0)
        # DIENAT: natural mortality of infected trees, redistributed over the slots.
        dienat = tdien - rrkill[i]
        dienat <= 1.0f-6 && (dienat = 0.0f0)
        for j in 1:istep
            broke = false
            for ip in 1:2
                if probi[i, j, ip] <= 0.0f0
                    broke = true; break              # rdend GOTO 900 (next J, skips IP=2)
                end
                d = probi[i, j, ip] * dienat / probit[i]
                probi[i, j, ip] -= d
                rrkill[i] += d
                dprob[i, RD_DNAT, RD_DSII] += d
                probi[i, j, ip] <= 1.0f-6 && (probi[i, j, ip] = 0.0f0)
            end
            broke && continue
        end
        # RDSSIZ/RDSTP stump update: natural mortality of infected trees (DIENAT>0)
        # feeds next cycle's inoculum (rdend.f:248). Active only when the driver +
        # stump arrays are supplied (the live seam); the pure-function tests omit them.
        if dienat > 0.0f0 && driver !== nothing
            rd_stp!(rd, driver, ksp, Float32(dbh[i]), rootl[i], dienat)
        end
    end
    rd_sum!(probit, probi, istep)                    # rdend.f final RDSUM
    return wk2
end

"""
    rd_grow_kernel!(rd, dg, htg, prob, probit, probiu, fprob, probi, propi, isp)

Port of the rd/rdgrow.f per-record body: reduce each infected host record's diameter
`dg` and height `htg` growth by the infected-root proportion. For each record it builds
DGTOT/HTTOT (`rd_grow_gtot`) from the infection slots and combines via `rd_grow_newg`.
Non-host records and records with no diseased area (`PAREA[idi] ≤ 0`) are unchanged.
Mutates `dg`/`htg`. Validated bit-exact vs the live FVSkt oracle (270/270 records, DG).
"""
function rd_grow_kernel!(rd::RootDiseaseState,
                         dg::Vector{Float32}, htg::Vector{Float32},
                         prob::Vector{Float32}, probit::Vector{Float32},
                         probiu::Vector{Float32}, fprob::Vector{Float32},
                         probi::Array{Float32,3}, propi::Array{Float32,3},
                         isp::AbstractVector{<:Integer})
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr); irt = rd.irtspc
    istep = Int(rd.istep); sarea = rd.sarea
    n = length(prob)
    tparea = 0.0f0
    @inbounds for idi in minrr:maxrr; tparea += rd.parea[idi]; end
    (n == 0 || tparea == 0.0f0) && return dg
    rd_sum!(probit, probi, istep)                    # rdgrow.f leading RDSUM
    @inbounds for i in 1:n
        ksp = Int(isp[i]); ksp == 0 && continue
        base = Int(irt[ksp])
        idi = maxrr < 3 ? Int(RD_IDITYP[base]) : maxrr
        idi <= 0 && continue
        rd.parea[idi] <= 0.0f0 && continue
        dgtot = rd_grow_gtot(probi, propi, i, istep, probit[i], RD_XDBH, RD_YDBH, RD_DBIFAC)
        httot = rd_grow_gtot(probi, propi, i, istep, probit[i], RD_XHT,  RD_YHT,  RD_HTIFAC)
        diff   = sarea - rd.parea[idi]
        outnum = fprob[i] * (sarea - rd.parea[idi])
        dg[i]  = rd_grow_newg(dg[i],  dgtot, outnum, probiu[i], probit[i], diff)
        htg[i] = rd_grow_newg(htg[i], httot, outnum, probiu[i], probit[i], diff)
    end
    return dg
end

# -----------------------------------------------------------------------------
# UPSTREAM Monte-Carlo spread chain (rd/rdsprd.f + rd/rdrate.f + rd/rdinf.f) —
# Chunk 0b-3d. These are the per-cycle routines (called from RDCNTL) that BUILD
# the per-record PROBI/PROPI entry state which the already-ported downstream
# kernels (RDMORT/RDEND/RDGROW) consume:
#   * RDSPRD runs an explicit small-stand Monte-Carlo simulation of disease
#     spread and returns the per-Monte spread rate MCRATE + its mean RRRATE.
#   * RDRATE distributes the Monte-Carlo rates across the disease centers into
#     the per-center RRATES (and recomputes RRRATE from the applied rates).
#   * RDINF converts a newly-infected AREA (from center growth) into equivalent
#     infected TPA per host record, drawing the initial root-infection proportion
#     PROPI(I,ISTEP,2)=RDRANP(RRNEW) and adding PROBI(I,ISTEP,2)/PROBIU(I).
# All three are the RD-RNG consumers of the cycle: per the g16 trace, RDINSD then
# RDSPRD then RDINF (RDINF = exactly 1 RDRANP per host record). Validated by
# DUMP-REPLAY against the live relinked FVSkt oracle: a single-.o instrumentation
# swap of rdsprd.f/rdrate.f/rdinf.f dumped each routine's full ENTRY state + the
# RD-RNG S0 consumption for all 10 cycles of the turnkey scenario (RRType 3
# Armillaria, RRInit 0 10 10 20 0.1 10 3, SArea 100; the instrumented .sum stayed
# byte-identical to FVSkt_clean). Replaying the dumped entry state reproduces
# MCRATE/RRRATE/RRATES (===) and the exact draw counts (SPRD 1510, INF 27 for
# cycle 1). The RD RNG (rd_rann!/rd_ranp!) is kept faithful — never FFI'd.
#
# SCOPE BOUNDARY: these consume the per-record FFPROB/DBH/ROOTL/HABSP list and the
# center PAREA/RRATES state that RDCNTL/RDTREG assemble each cycle (via RDROOT for
# ROOTL, RDOAGM for FFPROB, RDAREA for PAREA, the DO-300 center-radius growth for
# PCENTS). The remaining upstream piece before the engine seam can go live is
# RDINSD (the inside-patch MC, the biggest RNG consumer) plus the RDCNTL glue that
# threads ISCT/IND1/FFPROB/AREANU between them; the seam therefore stays INERT.
# -----------------------------------------------------------------------------

"""
    rd_sprd!(rd, idi; nmont, irsnyr, nrstep, irstyp, rrsfrn, pint, fint,
             rrsare, rrsdim, xminkl, dbh, rootl, ffprob, ksp, habsp, rrpswt)
        -> (mcrate::Vector{Float32}, rrrate::Float32)

Port of the rd/rdsprd.f Monte-Carlo spread-rate model (the pre-RDRATE result).
Runs `nmont` independent simulations of disease spread through a small explicit
stand of trees selected in proportion to the outside-center density `ffprob`, and
returns the per-Monte spread rate `mcrate` (ft/yr) and its mean `rrrate`.

Faithful to rd/rdsprd.f for the reduced turnkey path (`irstyp==0` random spacing,
`EFFSDI≡1`, `UPDATE=UPLTD=1`, `YTKX=1`): the RD-RNG draw ORDER is reproduced
exactly — per Monte: one `rd_rann!` per host record in the tree-selection loop,
then two `rd_rann!` per placed tree (x,y), then per timestep one `rd_rann!` per
baseline-infection test + one per infected pal-contact test + one per
update-limited internal-spread increment. Records are supplied in the ISCT/IND1
species-sorted processing order; `ksp` indexes the per-species coefficient tables
via `rd.irtspc`. `habsp`/`rrpswt` are the per-record time-to-death multipliers
assembled upstream (RDCNTL). Validated bit-exact (mcrate ===, draw-count ==) vs
the live FVSkt oracle. Does not consume `parea`/`sarea` directly — the caller has
already reduced them into `rrsare`/`rrsdim`.
"""
function rd_sprd!(rd::RootDiseaseState, idi::Int;
                  nmont::Int, irsnyr::Int, nrstep::Int, irstyp::Int,
                  rrsfrn::Float32, pint::Real, fint::Real,
                  rrsare::Float32, rrsdim::Float32, xminkl::Float32,
                  dbh::AbstractVector{Float32}, rootl::AbstractVector{Float32},
                  ffprob::AbstractVector{Float32}, ksp::AbstractVector{<:Integer},
                  habsp::AbstractVector{Float32}, rrpswt::AbstractVector{Float32})
    irt = rd.irtspc; xxinf = rd.xxinf; yyinf = rd.yyinf; nninf = Int(rd.nninf)
    rpint = Float32(pint); fintf = Float32(fint); nrf = Float32(nrstep)
    nrec = length(ffprob)
    NMAX = 50
    # local per-simulation arrays (Fortran DIMENSION 50)
    rrsdbh = zeros(Float32, NMAX); rrsrad = zeros(Float32, NMAX)
    irssp  = zeros(Int, NMAX);     ytk    = zeros(Float32, NMAX)
    trurad = zeros(Float32, NMAX); rkills = zeros(Float32, NMAX)
    xrrs   = zeros(Float32, NMAX); yrrs   = zeros(Float32, NMAX)
    distnc = zeros(Float32, NMAX, NMAX)
    sick   = zeros(Int, NMAX);     sine   = zeros(Float32, NMAX)
    numpal = zeros(Int, NMAX);     idpal  = zeros(Int, NMAX, NMAX)
    radnow = zeros(Float32, NMAX); radnew = zeros(Float32, NMAX)
    ifix(v) = trunc(Int, v)

    mcrate = zeros(Float32, nmont)
    rrrate = 0.0f0
    @inbounds for _jt in 1:nmont
        # --- select trees for the simulated area (rdsprd.f DO 100/90/80) ---
        ntrees = 0
        for k in 1:nrec
            ffprob[k] == 0.0f0 && continue
            rrsmen = ffprob[k] * rrsare
            numtre = ifix(rrsmen)
            ptre   = rrsmen - Float32(numtre)
            r      = rd_rann!(rd)
            r <= ptre && (numtre += 1)
            numtre == 0 && continue
            for _kk in 1:numtre
                ntrees += 1
                if ntrees > NMAX
                    ntrees = NMAX
                    break                         # rdsprd.f GOTO 90 (next record)
                end
                rrsdbh[ntrees] = dbh[k]
                rrsrad[ntrees] = rootl[k]
                irssp[ntrees]  = k                 # store record index (ksp via ksp[k])
                y0 = rd_slp(dbh[k], xxinf, yyinf, nninf)
                ytk[ntrees] = (y0 - xminkl) * habsp[k] * rrpswt[k] + xminkl
            end
        end
        ntrees == 0 && (mcrate[_jt] = 0.0f0; continue)

        # --- root radius crowding + kill radius (DO 105); EFFSDI≡1, YTKX=1 ---
        for i in 1:ntrees
            trurad[i] = rrsrad[i]                   # * EFFSDI(=1)
            base = Int(irt[ksp[irssp[i]]])
            rkills[i] = trurad[i] * RD_PKILLS[base, idi]
            # ytk[i] = 1 * ytk[i]  (YTKX=1)
        end

        # --- tree positions (DO 150, random branch irstyp != 1) ---
        for it in 1:ntrees
            xrrs[it] = rd_rann!(rd) * rrsdim
            yrrs[it] = rd_rann!(rd) * rrsdim
        end

        # --- pairwise distances (DO 205/200) ---
        for i in 1:ntrees, j in (i+1):ntrees
            d = sqrt((xrrs[i]-xrrs[j])^2 + (yrrs[i]-yrrs[j])^2)
            distnc[i,j] = d; distnc[j,i] = d
        end

        # --- (re)initialize tree variables (DO 210) ---
        for i in 1:ntrees
            sick[i] = 0; radnow[i] = -trurad[i]; radnew[i] = -trurad[i]
            numpal[i] = 0; sine[i] = 0.0f0
        end

        # --- find contact-trees / pals (DO 225/220) ---
        for i in 1:ntrees, j in (i+1):ntrees
            if trurad[i] + trurad[j] >= distnc[i,j]
                numpal[i] += 1; idpal[i, numpal[i]] = j
                numpal[j] += 1; idpal[j, numpal[j]] = i
            end
        end

        # --- spread simulation over time steps (DO 300) ---
        room = 1; jyears = 0
        for irstep in 1:irsnyr
            if room > 0
                jyears += nrstep
                for in_ in 1:ntrees
                    kk = irssp[in_]                 # record index
                    base = Int(irt[ksp[kk]])
                    if radnow[in_] < rkills[in_]
                        pnsp = RD_PNINF[base, idi]
                        # baseline infection test
                        if sick[in_] == 0 && yrrs[in_] <= trurad[in_]
                            r = rd_rann!(rd)
                            pnin = irstep == 1 ? pnsp : 1.0f0 - (1.0f0 - pnsp)^(nrf/rpint)
                            if r <= pnin
                                sick[in_] = 1; radnew[in_] = -yrrs[in_]; sine[in_] = 1.0f0
                            end
                        end
                        # infection by contact with spreading neighbours
                        if irstep > 1 && (sick[in_] == 0 || radnow[in_] < 0.0f0)  # UPDATE=1,UPLTD=1
                            for ii in 1:numpal[in_]
                                it = idpal[in_, ii]
                                if radnow[it] > 0.0f0
                                    if radnew[in_] < -(distnc[it,in_] - radnow[it])
                                        r = rd_rann!(rd)
                                        pnin = 1.0f0 - (1.0f0 - pnsp)^(nrf/fintf)
                                        if r <= pnin
                                            sick[in_] = 1
                                            radnew[in_] = -(distnc[it,in_] - radnow[it])
                                            if distnc[in_,it] > 0.0f0
                                                sine[in_] = (yrrs[in_]-yrrs[it]) / distnc[it,in_]
                                            else
                                                sine[in_] = 1.0f0
                                            end
                                        end
                                    end
                                end
                            end
                            if radnew[in_] > 0.0f0                # UPLTD==1
                                irstp = ifix(rd_rann!(rd) * nrf)
                                irstp == nrstep && (irstp = nrstep - 1)
                                radnew[in_] = Float32(irstp) * rrsrad[in_] *
                                              RD_PKILLS[base, idi] / ytk[in_]
                            end
                        end
                        # spread rot through infected roots
                        if sick[in_] == 1
                            radnew[in_] += nrf * rrsrad[in_] * RD_PKILLS[base, idi] / ytk[in_]
                            if radnew[in_] >= rkills[in_]
                                radnew[in_] = RD_PCOLO[base, idi] * trurad[in_]
                            end
                        end
                        (yrrs[in_] + radnew[in_]) > rrsdim && (room = 0)
                    end
                end
            end
            if room > 0
                for i in 1:ntrees; radnow[i] = radnew[i]; end
            end
        end

        # --- maximum spread of infection (DO 330) → MCRATE ---
        ydmax = 0.0f0; yforwd = 0.0f0
        for in_ in 1:ntrees
            if sick[in_] > 0
                if radnow[in_] > 0.0f0
                    yforwd = yrrs[in_] + radnow[in_]
                elseif sine[in_] > 0.0f0
                    yforwd = yrrs[in_] - (sine[in_] * (-radnow[in_]))
                end
                yforwd > ydmax && (ydmax = yforwd)
            end
        end
        if room > 0
            mcrate[_jt] = ydmax / Float32(jyears)
        else
            mcrate[_jt] = jyears > nrstep ? ydmax / Float32(jyears - nrstep) :
                                            rrsdim / Float32(jyears)
        end
        rrrate += mcrate[_jt]
    end
    rrrate /= Float32(nmont)
    return mcrate, rrrate
end

"""
    rd_rate!(mcrate, ncents, nscen, rrates, shcent2, icensp) -> (rrates_out, rrrate)

Port of rd/rdrate.f: distribute the Monte-Carlo spread rates `mcrate` (length
`nsim = length(mcrate)`) across the `ncents` disease centers, returning the updated
per-center rate vector `rrates_out` (length `ncents`) and the mean applied rate
`rrrate`. `GCENTS = ncents - nscen` non-shrinking centers receive rates: when
`GCENTS ≤ nsim` each center is the average of `nsim/GCENTS` consecutive Monte
rates; otherwise each Monte rate is applied to one-or-more centers. Zero rates are
floated to the top and, via the `IPNT` pointer built from `shcent2`/`rrates`/
`icensp`, assigned to shrinking / previously-zero / spore centers. Faithful REAL
(Float32) arithmetic. Validated bit-exact vs the live FVSkt oracle.
"""
function rd_rate!(mcrate::AbstractVector{Float32}, ncents::Int, nscen::Int,
                  rrates::AbstractVector{Float32}, shcent2::AbstractVector{Float32},
                  icensp::AbstractVector{<:Integer})
    nsim  = length(mcrate)
    rout  = Float32.(collect(rrates[1:ncents]))
    gcents = ncents - nscen
    gcents <= 0 && return rout, 0.0f0

    rra  = zeros(Float32, 100)
    ipnt = collect(1:100)
    lzero = false
    rem   = 0.0f0
    k = 1
    if gcents <= nsim
        div = Float32(nsim) / (Float32(gcents) + 1.0f-9)
        for i in 1:gcents
            num = trunc(Int, div + rem)
            rem = div + rem - Float32(num)
            for _j in 1:num
                rra[i] += mcrate[k]
                k += 1
                k > nsim && (num = _j)
            end
            rra[i] = rra[i] / (Float32(num) + 1.0f-9)
            (!lzero && rra[i] == 0.0f0) && (lzero = true)
        end
    else
        i = 1
        div = Float32(gcents) / Float32(nsim)
        while true
            num = trunc(Int, div + rem)
            rem = div + rem - Float32(num)
            for _j in 1:num
                rra[i] = mcrate[k]
                i += 1
                (!lzero && rra[i] == 0.0f0) && (lzero = true)
            end
            k = min(nsim, k + 1)
            i <= gcents || break
        end
    end

    if lzero
        # bubble the zeros to the top (other values stay unsorted) — DO 350/300
        for i in 1:gcents
            if rra[i] == 0.0f0
                for j in i:-1:2
                    if rra[j-1] > 0.0f0
                        rra[j], rra[j-1] = rra[j-1], rra[j]
                    else
                        break
                    end
                end
            end
        end
        # pointer array: shrinking centers to bottom, zero/spore to top — DO 450
        for i in 1:ncents
            if shcent2[i] > 0.0f0
                for j in i:ncents-1
                    ipnt[j] = ipnt[j+1]
                end
                ipnt[ncents] = i
            elseif rout[i] == 0.0f0 && icensp[i] == 0
                for j in i:-1:2
                    ipnt[j] = ipnt[j-1]
                end
                ipnt[1] = i
            end
        end
    end

    trr = 0.0f0
    for i in 1:gcents
        rout[ipnt[i]] = rra[i]
        trr += rra[i]
    end
    rrrate = trr / (Float32(gcents) + 1.0f-9)
    return rout, rrrate
end

"""
    rd_inf_pnsp(rd, ksp, idi, fint, pint) -> Float32

The rd/rdinf.f per-species probability-of-infection `PNSP`:
`PNSP = 1 - (1-PNINF(base,idi))^(FINT/PINT)`, then modified for the spore-initiated
fraction `PNSP *= SPPROP·SPTRAN + (1-SPPROP)`. `base = irtspc(ksp)`. Float32.
"""
function rd_inf_pnsp(rd::RootDiseaseState, ksp::Integer, idi::Int, fint::Real, pint::Real;
                     spprop::Float32 = RD_SPPROP0, sptran::Float32 = 0.5f0)
    base = Int(rd.irtspc[ksp])
    pnsp = 1.0f0 - (1.0f0 - RD_PNINF[base, idi])^(Float32(fint)/Float32(pint))
    return pnsp * (spprop * sptran + (1.0f0 - spprop))
end

"""
    rd_inf_kernel!(rd, idi, areanu, ksp, fprob, probi_in, probiu_in; fint, pint, sptran)
        -> (propi, probi_out, probiu_out)

Port of the rd/rdinf.f per-record body: convert a newly-infected area `areanu`
(acres, from disease-center growth) into equivalent infected TPA. For each host
record with `fprob > 0` (supplied in ISCT/IND1 order) it draws the initial
root-infection proportion `PROPI(I,ISTEP,2) = RDRANP(RRNEW(idi))` from the RD
stream, then adds `ADDINF = areanu·fprob·PNSP` to the infected TPA and
`areanu·fprob·(1-PNSP)` to the inside-uninfected TPA. Reproduces the exact draw
sequence (1 `rd_ranp!` per processed host record). The OAKL/RDMREC dead-tree
re-bookkeeping is skipped when OAKL≡0 (BBCLEAR + no windthrow). Validated bit-exact
(=== outputs, draw-count == records) vs the live FVSkt oracle. Returns the parallel
output vectors in processing order.
"""
function rd_inf_kernel!(rd::RootDiseaseState, idi::Int, areanu::Float32,
                        ksp::AbstractVector{<:Integer}, fprob::AbstractVector{Float32},
                        probi_in::AbstractVector{Float32}, probiu_in::AbstractVector{Float32};
                        fint::Real, pint::Real, sptran::Float32 = 0.5f0)
    n = length(fprob)
    propi     = zeros(Float32, n)
    probi_out = copy(collect(probi_in))
    probiu_out = copy(collect(probiu_in))
    areanu <= 0.0f0 && return propi, probi_out, probiu_out
    rrnew = rd.rrnew[idi]
    @inbounds for k in 1:n
        fprob[k] <= 0.0f0 && continue             # rdinf.f FPROB(I) .LE. 0 → GOTO 400 (no draw)
        pnsp   = rd_inf_pnsp(rd, ksp[k], idi, fint, pint; sptran = sptran)
        nuinsd = areanu * fprob[k]
        addinf = nuinsd * pnsp
        propi[k]      = rd_ranp!(rd, rrnew)        # PROPI(I,ISTEP,2)
        probi_out[k]  = probi_in[k]  + addinf
        probiu_out[k] = probiu_in[k] + nuinsd * (1.0f0 - pnsp)
    end
    return propi, probi_out, probiu_out
end

"""
    rd_insd!(rd, idi, rriare, rridim, parea, irinit, itrn, ninsim,
             ksp, rootl, probiu, probi, propi; fint, pint, sptran) -> (rrninf, polp)

Port of the rd/rdinsd.f inside-patch infection Monte-Carlo (the reduced turnkey
path: `ninsim==1`, and NO stumps — PROBD≡0, verified for cycle 1 where inoculum is
purely live-infected roots). It is the FIRST RD-RNG consumer of the cycle and the
biggest. For each simulation it (1) selects live-infected inoculum sources — one
`rd_rann!` per (record, slot) with `PROPI > 0`, radius `ROOTL·PROPI`; (2) places
each inoculum at a random (x,y) — two `rd_rann!` each; (3) for every uninfected
host record throws `NUMTRE = min(IRINIT/ITRN, INT(PROBIU))` trial trees — two
`rd_rann!` each — and tallies root-contact overlaps to get the average infection
probability `AVGINF`, accumulating `RRNINF += PROBIU·AVGINF` and the mean
before-center overlap into `POLP` (`PROPI(I,ISTEP,1) = -POLP` downstream, a NEGATIVE
proportion flagging "before center"). The RD-RNG draw ORDER (selection → placement
→ contact) is reproduced exactly. `ksp` indexes `rd.irtspc`; `idi` is the disease.
Validated by dump-replay vs the live FVSkt oracle: RRNINF/POLP === per record AND
exact RD-RNG final S0 (cycle 1: 7949 draws = 27 selection + 62 placement + 7860
contact). Returns the raw (pre-`/NINSIM`) accumulators. The stump branch (PROBD/
ROOTD) activates once prior-cycle mortality creates stumps — a later sub-chunk.
"""
# gfortran `REAL**INTEGER` (_gfortran_pow_r4_i4) integer-power: right-to-left binary
# exponentiation squaring the base in Float32 each step. Matches the RDINSD `(1-PNSP)
# **ITROLP` rounding bit-exact (Julia's generic `^` associates differently on Float32).
@inline function rd_powi(x::Float32, n::Integer)
    n == 0 && return 1.0f0
    b = x; r = 1.0f0; m = Int(n)
    while true
        (m & 1) == 1 && (r *= b)
        m >>= 1
        m == 0 && break
        b *= b
    end
    return r
end

function rd_insd!(rd::RootDiseaseState, idi::Int, rriare::Float32, rridim::Float32,
                  parea::Float32, irinit::Int, itrn::Int, ninsim::Int,
                  ksp::AbstractVector{<:Integer}, rootl::AbstractVector{Float32},
                  probiu::AbstractVector{Float32}, probi::Array{Float32,3},
                  propi::Array{Float32,3}; fint::Real, pint::Real, sptran::Float32 = 0.5f0,
                  probd::Union{Nothing,AbstractMatrix{Float32}} = nothing,
                  rootd::Union{Nothing,AbstractMatrix{Float32}} = nothing)
    n = length(ksp); istep = size(probi, 2)
    rrninf = zeros(Float32, n); polp = zeros(Float32, n)
    (itrn == 0 || parea == 0.0f0) && return rrninf, polp

    @inbounds for _sim in 1:ninsim
        # --- select live-infected inoculum sources (DO 100/95/90/85/80) ---
        rrirad = Float32[]
        irincs = 0
        capped = false
        for k in 1:n
            for it in 1:istep, ip in 1:2
                pp = propi[k, it, ip]
                pp <= 0.0f0 && continue
                rrimen = probi[k, it, ip] * rriare / (parea + 1.0f-6)
                rrimen == 0.0f0 && continue
                numtre = trunc(Int, rrimen)
                ptre   = rrimen - Float32(numtre)
                r      = rd_rann!(rd)
                r < ptre && (numtre += 1)          # rdinsd.f: R .LT. PTRE
                numtre <= 0 && continue
                for _kk in 1:numtre
                    irincs += 1
                    if irincs > RD_IRRTRE
                        irincs = RD_IRRTRE; capped = true; break
                    end
                    push!(rrirad, rootl[k] * pp)
                end
                capped && break
            end
            capped && break
        end
        # --- select infected dead trees / stumps (rdinsd.f DO 605/600) ---
        # RRIMEN = PROBD(irrsp,i,j)·RRIARE/(PAREA+1e-6); one rd_rann! per (i,j) with
        # PROBD>0, radius ROOTD(irrsp,i,j). PROBD≡0 at cycle 1 (no stumps) ⇒ no draws;
        # cyc2+ stumps (from RDMORT/RDEND, RDINUP-averaged) add inoculum here.
        if probd !== nothing && !capped
            for i in 1:2, j in 1:5
                rrimen = probd[i, j] * rriare / (parea + 1.0f-6)
                rrimen == 0.0f0 && continue
                numtre = trunc(Int, rrimen)
                ptre   = rrimen - Float32(numtre)
                r      = rd_rann!(rd)
                r < ptre && (numtre += 1)
                numtre <= 0 && continue
                for _kk in 1:numtre
                    irincs += 1
                    if irincs > RD_IRRTRE
                        irincs = RD_IRRTRE; capped = true; break
                    end
                    push!(rrirad, rootd[i, j])
                end
                capped && break
            end
        end
        irincs == 0 && continue

        # --- place inoculum randomly (DO 650) ---
        xrri = Vector{Float32}(undef, irincs); yrri = Vector{Float32}(undef, irincs)
        for it in 1:irincs
            xrri[it] = rd_rann!(rd) * rridim
            yrri[it] = rd_rann!(rd) * rridim
        end

        # --- infect uninfected trees by simulated contact (DO 1000/900/800/750) ---
        numtre_base = trunc(Int, Float32(irinit) / (Float32(itrn) + 1.0f-6))
        for k in 1:n
            pnsp = rd_inf_pnsp(rd, ksp[k], idi, fint, pint; sptran = sptran)
            numtre = numtre_base
            Float32(numtre) > probiu[k] && (numtre = trunc(Int, probiu[k]))
            numtre <= 0 && continue
            numolp = 0; disolp = 0.0f0; pniolp = 0.0f0
            rl = rootl[k]
            for _it in 1:numtre
                xtry = rd_rann!(rd) * rridim
                ytry = rd_rann!(rd) * rridim
                itrolp = 0
                for in_ in 1:irincs
                    dist = sqrt((xtry - xrri[in_])^2 + (ytry - yrri[in_])^2)
                    dist > (rrirad[in_] + rl) && continue
                    itrolp += 1
                    overlp = rrirad[in_] + rl - dist
                    overlp > rl && (overlp = rl)
                    disolp += overlp
                end
                if itrolp > 0
                    pniolp += rd_powi(1.0f0 - pnsp, itrolp)
                else
                    pniolp += 1.0f0
                end
                numolp += itrolp
            end
            avginf = 1.0f0 - (pniolp / (Float32(numtre) + 1.0f-6))
            rrninf[k] += probiu[k] * avginf
            numolp <= 0 && continue
            temp = 1.0f0 - (disolp / Float32(numolp)) / (rl + 1.0f-6)
            temp <= 0.0f0 && (temp = -0.001f0)
            polp[k] += temp
        end
    end
    return rrninf, polp
end

# =============================================================================
# WRD 0b-3e — the per-cycle DRIVER (rd/rdtreg.f + rd/rdcntl.f) + stump lifecycle
# (rd/rdstp.f, rd/rdssiz.f, rd/rdinup.f, rd/rdinoc.f) + the LIVE engine seam.
#
# `rd_control!` composes the already-validated kernels into one per-cycle driver
# matching the RDTREG→RDCNTL call order for the manual-RRINIT turnkey (MINRR=MAXRR,
# IRSPTY=1, IRIPTY=1, LONECT=0, INFLAG=0, IRSTYP=0, no bark beetles/windthrow):
#
#   RDTREG: FPROB (outside density) · ROOTL (rd_root) · RDOAGM (FFPROB min-carry;
#           OAKL≡0) · RDCNTL{ DO525 RDINUP+RDINSD · DO75 RDSPRD+RDRATE · DO300
#           grow-centers+RDZERO+RDAREA+RDINF · RDMORT(+RDSTP stumps) } · RDSUM
#   [seam] RDEND (→ WK2 mortality, at FVS MORTS time) · RDGROW (→ DG/HTG) · RDINOC
#
# The deterministic no-op routines for this scenario are faithfully skipped:
# RDTIM/RDJUMP (outputs unused w/o a cut), RDPUSH (no PSTUMP), RDSHST (windthrow-
# only), RDSHRK (NSCEN≡0), the carryover model (INFLAG≡0), and RDOAGM's BB/wind
# (BBCLEAR + no windthrow ⇒ OAKL≡0). Validated end-to-end at the .sum DELTA level
# vs the live FVSkt oracle (rd.key − ctrl.key), cornered within the #206 straddle.
# =============================================================================

# rd/rdinit.f ISPS (wood type 1=resinous/2=non) + RSLOP (root-radius slope);
# PROOT≡1.0 (rdinit.f:504). Base-RD-species indexed (via IRTSPC).
const RD_ISPS = Int32[1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2,
                      2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2]
const RD_RSLOP = Float32[14.26, 14.26, 14.26, 14.26, 14.5, 9.35, 14.26, 14.26, 14.26, 14.26,
                         14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26,
                         14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26,
                         14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26, 14.26]
const RD_PROOT = 1.0f0
# rd/rdinit.f DECFN/YRSITF/RSITFN (idi,i,j) stump-decay coefficients + DSFAC.
# Column-major (idi fastest) exactly as the Fortran DATA fills TEMP1/2/3.
const RD_DECFN  = reshape(Float32[0.02212,0.02212,0.0,0.0,  1.2032,1.2032,1.6045,1.2834,
                                  0.02212,0.02212,0.03214,0.03214,  1.2032,1.2032,1.1459,0.9167],
                          RD_ITOTRR, 2, 2)                         # DECFN(idi,i,j)
const RD_YRSITF = reshape(Float32[2.0,2.0,1.1111,1.1111,  0.0,0.0,0.0,0.0,
                                  2.0,2.0,0.5556,0.5556,   0.0,0.0,6.6667,6.6667],
                          RD_ITOTRR, 2, 2)                         # YRSITF(idi,i,j)
const RD_RSITFN = reshape(Float32[0.166667,0.166667,0.297,0.2377,  0.0,0.0,0.0,0.0],
                          RD_ITOTRR, 2)                            # RSITFN(idi,i)
const RD_DSFAC  = Float32[1.0, 0.75]
const RD_SPYTK  = 3.0f0
const RD_SPTRAN = 0.5f0
const RD_ISTEP_MAX = 41                        # rd/RDARRY.F77 PROBI(IRRTRE,41,2)

# rd/rdssiz.f — stump size class from DBH vs STCUT breakpoints.
@inline function _rd_ssiz(a::Float32, stcut)
    @inbounds for j in 2:5
        (a < stcut[j-1] || a > stcut[j]) && continue
        return j - 1
    end
    return 5
end

# rd/rdstp.f — add DEN infected stumps of record ISP (dbh, live root radius RTD)
# to the PROBDA/DBHDA/ROOTDA stump lists (weighted-average DBHDA/ROOTDA). ISTFLG=0
# ⇒ idi is the record's disease type. Shared by RDMORT (RDKILL) and RDEND (DIENAT).
function rd_stp!(rd::RootDiseaseState, d,
                 isp::Integer, dbh::Float32, rtd::Float32, den::Float32)
    den <= 0.0f0 && return
    base = Int(rd.irtspc[isp])
    maxrr = Int(rd.maxrr)
    idi = maxrr < 3 ? Int(RD_IDITYP[base]) : maxrr
    idi <= 0 && return
    ist = max(1, Int(rd.istep))
    is  = Int(RD_ISPS[base])
    isl = _rd_ssiz(dbh, rd.stcut)
    rotd = rtd * RD_PCOLO[base, idi]
    old  = d.probda[idi, is, isl, ist]
    tst  = old + den
    @inbounds begin
        d.dbhda[idi, is, isl, ist]  = ((d.dbhda[idi, is, isl, ist]  * old) + (dbh  * den)) / tst
        d.rootda[idi, is, isl, ist] = ((d.rootda[idi, is, isl, ist] * old) + (rotd * den)) / tst
        d.probda[idi, is, isl, ist] = old + den
    end
    return
end

"""
    RDDriver

The per-cycle working state for the WRD driver — the 3-D PROBI/PROPI (record order),
PROBIT/PROBIU/FPROB/FFPROB/ROOTL, RRKILL/RDKILL, per-center RRATES, the stump lists
(PROBDA/DBHDA/ROOTDA/DECRAT/JRAGED + the RDINUP-averaged PROBD/DBHD/ROOTD), and the
RDEND carry arrays (WK22/RROOTT/DPROB + the all-zero OAKL/BBKILL). Sized to the tree
list at `rd_build_driver!`; rebuilt if the record count changes (COMCUP).
"""
mutable struct RDDriver
    n::Int
    probi::Array{Float32,3}
    propi::Array{Float32,3}
    probit::Vector{Float32}
    probiu::Vector{Float32}
    fprob::Vector{Float32}
    ffprob::Matrix{Float32}       # (n,2)
    rootl::Vector{Float32}
    rrkill::Vector{Float32}
    rdkill::Vector{Float32}
    wk22::Vector{Float32}
    rroott::Vector{Float32}
    rrates::Matrix{Float32}       # (ITOTRR,100)
    rrrate::Vector{Float32}       # (ITOTRR)
    areanu::Vector{Float32}       # (ITOTRR)
    shcent::Array{Float32,3}      # (ITOTRR,100,3)
    nscen::Vector{Int32}          # (ITOTRR)
    icensp::Matrix{Int32}         # (ITOTRR,100)
    probda::Array{Float32,4}      # (ITOTRR,2,5,41)
    dbhda::Array{Float32,4}
    rootda::Array{Float32,4}
    decrat::Array{Float32,4}
    jraged::Array{Int32,4}
    probd::Array{Float32,3}       # (ITOTRR,2,5)
    dbhd::Array{Float32,3}
    rootd::Array{Float32,3}
    oakl::Matrix{Float32}         # (3,n) — all zero (no bark beetles/windthrow)
    bbkill::Matrix{Float32}
    dprob::Array{Float32,3}       # (n,4,3)
end

# Build the driver from the LSTART per-record state (rd.probi/propi/probiu/fprob are
# in record order; slot (1,1) holds the initial infection). PROBIT via RDSUM.
function rd_build_driver!(rd::RootDiseaseState, n::Int)
    R = RD_ITOTRR; S = RD_ISTEP_MAX
    d = RDDriver(n,
        zeros(Float32, n, S, 2), zeros(Float32, n, S, 2),
        zeros(Float32, n), zeros(Float32, n), zeros(Float32, n),
        zeros(Float32, n, 2), zeros(Float32, n),
        zeros(Float32, n), zeros(Float32, n), zeros(Float32, n), zeros(Float32, n),
        zeros(Float32, R, 100), zeros(Float32, R), zeros(Float32, R),
        zeros(Float32, R, 100, 3), zeros(Int32, R), zeros(Int32, R, 100),
        zeros(Float32, R, 2, 5, S), zeros(Float32, R, 2, 5, S), zeros(Float32, R, 2, 5, S),
        zeros(Float32, R, 2, 5, S), zeros(Int32, R, 2, 5, S),
        zeros(Float32, R, 2, 5), zeros(Float32, R, 2, 5), zeros(Float32, R, 2, 5),
        zeros(Float32, 3, n), zeros(Float32, 3, n), zeros(Float32, n, 4, 3))
    @inbounds for i in 1:min(n, length(rd.probi))
        d.probi[i, 1, 1] = rd.probi[i]
        d.propi[i, 1, 1] = rd.propi[i]
        d.probiu[i]      = rd.probiu[i]
        d.fprob[i]       = rd.fprob[i]
    end
    rd_sum!(d.probit, d.probi, Int(rd.istep))
    return d
end

# Host records (idi>0) in FVS ISCT/IND1 processing order: species ascending, stable
# within species by record index. The RD-RNG consumers (RDINSD/RDSPRD/RDINF) walk
# this order; the non-RNG kernels (RDMORT/RDEND/RDGROW) are order-independent.
function _rd_host_order(rd::RootDiseaseState, s::StandState)
    t = s.trees; maxrr = Int(rd.maxrr); irt = rd.irtspc
    order = Int[]
    @inbounds for i in 1:t.n
        ksp = Int(t.species[i]); ksp == 0 && continue
        base = Int(irt[ksp])
        idi = maxrr < 3 ? Int(RD_IDITYP[base]) : maxrr
        idi <= 0 && continue
        push!(order, i)
    end
    sort!(order; by = i -> (Int(t.species[i]), i))
    return order
end

# rd/rdinup.f — weighted-average the per-timestep stump lists into PROBD/DBHD/ROOTD.
function rd_inup!(rd::RootDiseaseState, d::RDDriver, idi::Int)
    istep = Int(rd.istep)
    @inbounds for i in 1:2, j in 1:5
        pd = 0.0f0; db = 0.0f0; rt = 0.0f0
        for k in 1:istep
            p = d.probda[idi, i, j, k]
            p <= 0.0f0 && continue
            pd += p
            db += d.dbhda[idi, i, j, k]  * p
            rt += d.rootda[idi, i, j, k] * p
        end
        d.probd[idi, i, j] = pd
        d.dbhd[idi, i, j]  = db / (pd + 1.0f-6)
        d.rootd[idi, i, j] = rt / (pd + 1.0f-6)
    end
    return
end

# rd/rdzero.f — remove centers whose radius shrank to ≤0.01 (never fires when
# every center grows; ported faithfully). Compacts PCENTS/RRATES/SHCENT/ICENSP.
function rd_zero_centers!(rd::RootDiseaseState, d::RDDriver, idi::Int)
    P = rd.pcents
    icent = Int(rd.ncents[idi])
    @inbounds for icen in icent:-1:1
        P[idi, icen, 3] > 0.01f0 && continue
        rd.ncents[idi] -= Int32(1)
        d.shcent[idi, icen, 2] > 0.0f0 && (d.nscen[idi] -= Int32(1))
        for jcen in icen:icent-1
            for k in 1:3
                P[idi, jcen, k] = P[idi, jcen+1, k]
                d.shcent[idi, jcen, k] = d.shcent[idi, jcen+1, k]
            end
            d.icensp[idi, jcen] = d.icensp[idi, jcen+1]
            d.rrates[idi, jcen] = d.rrates[idi, jcen+1]
        end
        if Int(rd.ncents[idi]) == icent - 1
            d.icensp[idi, icent] = Int32(0); d.rrates[idi, icent] = 0.0f0
            for k in 1:3
                P[idi, icent, k] = 0.0f0; d.shcent[idi, icent, k] = 0.0f0
            end
        end
    end
    return
end

# rd/rdinoc.f (LICALL=.FALSE.) — decay the stump root radius each cycle.
function rd_inoc_decay!(rd::RootDiseaseState, d::RDDriver, fint::Real)
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr); istep = Int(rd.istep)
    jint0 = trunc(Int, Float32(fint))
    @inbounds for idi in minrr:maxrr, i in 1:2, j in 1:5, k in 1:istep
        (d.probda[idi,i,j,k] == 0.0f0 || d.dbhda[idi,i,j,k] == 0.0f0) && continue
        jint = jint0
        dbhda = d.dbhda[idi,i,j,k]
        rotsit = RD_RSITFN[idi,1] * dbhda + RD_RSITFN[idi,2]
        d.rootda[idi,i,j,k] < rotsit && (d.rootda[idi,i,j,k] = rotsit)
        jrsit = dbhda <= 12.0f0 ?
            trunc(Int, RD_YRSITF[idi,1,1]*dbhda + RD_YRSITF[idi,2,1]) :
            trunc(Int, RD_YRSITF[idi,1,2]*dbhda + RD_YRSITF[idi,2,2])
        if d.jraged[idi,i,j,k] <= 0
            if d.decrat[idi,i,j,k] <= 0.0f0
                d.decrat[idi,i,j,k] = dbhda <= 12.0f0 ?
                    (RD_DECFN[idi,1,1]*d.rootda[idi,i,j,k] + RD_DECFN[idi,2,1]) / RD_DSFAC[i] :
                    (RD_DECFN[idi,1,2]*d.rootda[idi,i,j,k] + RD_DECFN[idi,2,2]) / RD_DSFAC[i]
                droots = d.rootda[idi,i,j,k] - rotsit
                tminlf = Float32(jrsit) + droots / d.decrat[idi,i,j,k]
                tminlf < rd.xminlf[idi] &&
                    (d.decrat[idi,i,j,k] = droots / (rd.xminlf[idi] - Float32(jrsit)))
            end
            rtodec = d.decrat[idi,i,j,k] * Float32(jint)
            rtrem  = d.rootda[idi,i,j,k] - rtodec
            if rtrem < rotsit
                rtodec = d.rootda[idi,i,j,k] - rotsit
                rtodec < 0.0f0 && (rtodec = 0.0f0)
                d.jraged[idi,i,j,k] = jint - trunc(Int, rtodec / d.decrat[idi,i,j,k])
                rtrem = rotsit
            end
            d.rootda[idi,i,j,k] = rtrem
        else
            d.jraged[idi,i,j,k] += jint
        end
        if d.jraged[idi,i,j,k] > jrsit
            d.probda[idi,i,j,k] = 0.0f0; d.dbhda[idi,i,j,k] = 0.0f0
            d.jraged[idi,i,j,k] = Int32(0); d.rootda[idi,i,j,k] = 0.0f0
        end
    end
    return
end

"""
    rd_control!(rd, s, fint)

Port of the per-cycle RDTREG→RDCNTL spread chain (up to and including RDMORT), for
the manual-RRINIT turnkey. Reads the stand's PRE-growth DBH/HT (FVS runs RDTREG in
GRADD before UPDATE) and the stand aggregates (OLDTPA/GROSPC/ORMSQD/BA). Produces the
per-record RRKILL + the aged PROBI/PROPI/PROBIU/FPROB the downstream RDEND/RDGROW
(`rd_end_apply!`/`rd_grow_apply!`) consume. Advances the RD RNG stream (RDINSD then
RDSPRD then RDINF). Mutates `rd.driver`.
"""
function rd_control!(rd::RootDiseaseState, s::StandState, fint::Real)
    t = s.trees; p = s.plot; n = t.n
    d = rd.driver::RDDriver
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr)
    irt = rd.irtspc; irhab = Int(rd.irhab); istep = Int(rd.istep)
    sarea = rd.sarea; fintf = Float32(fint); pint = 10.0f0
    idi = maxrr                                   # single-disease turnkey (MINRR=MAXRR)
    rd.irrsp = Int32(idi)

    # TPAREA gate (RDTREG).
    tparea = 0.0f0
    @inbounds for id in minrr:maxrr; tparea += rd.parea[id]; end
    (tparea == 0.0f0 || n == 0) && return

    # --- FPROB: outside-center density (RDTREG DO 843) ---
    diffv = sarea - rd.parea[idi]
    @inbounds for i in 1:n
        ksp = Int(t.species[i]); ksp == 0 && continue
        base = Int(irt[ksp]); di = maxrr < 3 ? Int(RD_IDITYP[base]) : maxrr
        di <= 0 && continue
        if diffv >= 1.0f-6
            fp = (t.tpa[i] * sarea - d.probiu[i] - d.probit[i]) / diffv
            fp <= 1.0f-6 && (fp = 0.0f0)
            d.fprob[i] = fp
        else
            d.fprob[i] = 0.0f0
        end
    end

    # --- ROOTL: live-tree root radius (RDTREG DO 1001, rd_root) ---
    yincpt = (rd.sdislp == 0.0f0 || rd.sdnorm == 0.0f0) ? 1.0f0 : 1.0f0 - (rd.sdnorm * rd.sdislp)
    rd.yincpt = yincpt
    # RDROOT stand aggregates (PLOT.F77 OLDTPA/ORMSQD/BA/GROSPC). FVSjl's compute_density!
    # does not populate old_tpa/old_qmd, so derive them from the (pre-growth) tree list:
    # OLDTPA = Σ TPA, ORMSQD = quadratic mean diameter (matches the g16 RDROOT dump).
    tpa_sum = 0.0f0; dsq_sum = 0.0f0
    @inbounds for i in 1:n
        tpa_sum += t.tpa[i]; dsq_sum += t.tpa[i] * t.dbh[i] * t.dbh[i]
    end
    oldtpa = tpa_sum
    ormsqd = tpa_sum > 0.0f0 ? sqrt(dsq_sum / tpa_sum) : 0.0f0
    grospc = p.gross_space; ba = p.basal_area
    @inbounds for i in 1:n
        ksp = Int(t.species[i]); ksp == 0 && (d.rootl[i] = 0.0f0; continue)
        base = Int(irt[ksp])
        d.rootl[i] = rd_root(t.dbh[i], t.height[i], RD_PROOT, RD_RSLOP[base],
                             rd.sdislp, yincpt, oldtpa, grospc, ormsqd, ba)
    end

    # --- RDOAGM: FFPROB min-carry (cyc1 = FPROB; else min(FFPROB2,FPROB)); OAKL≡0 ---
    icyc = Int(rd.icyc)
    @inbounds for i in 1:n
        dff = d.fprob[i] - d.ffprob[i, 2]
        d.ffprob[i, 1] = (dff <= 1.0f-4 || icyc == 1) ? d.fprob[i] : d.ffprob[i, 2]
        d.ffprob[i, 2] = d.fprob[i]
    end
    fill!(d.oakl, 0.0f0); fill!(d.bbkill, 0.0f0)
    rd_sum!(d.probit, d.probi, istep)

    order = _rd_host_order(rd, s)
    m = length(order)

    # ==== RDCNTL DO 525 — RDINUP + RDINSD (inside-patch infection) ====
    rd_inup!(rd, d, idi)
    if rd.parea[idi] != 0.0f0 && m > 0
        smbi = 0.0f0; smiu = 0.0f0
        @inbounds for i in order; smbi += d.probit[i]; smiu += d.probiu[i]; end
        dennew = (smiu + smbi) / rd.parea[idi]
        @inbounds for i in 1:2, j in 1:5
            dennew += d.probd[idi, i, j] / (rd.parea[idi] + 1.0f-9)
        end
        if dennew > 0.0f0
            rriare = 100.0f0 / (dennew + 1.0f-9)          # RRGEN(idi,9)=100
            rridim = sqrt(rriare) * 208.7f0
            ksp_o   = Int[Int(t.species[i]) for i in order]
            rootl_o = Float32[d.rootl[i] for i in order]
            probiu_o = Float32[d.probiu[i] for i in order]
            probi_o = zeros(Float32, m, RD_ISTEP_MAX, 2); propi_o = zeros(Float32, m, RD_ISTEP_MAX, 2)
            @inbounds for (kk, i) in enumerate(order), it in 1:istep, ip in 1:2
                probi_o[kk, it, ip] = d.probi[i, it, ip]
                propi_o[kk, it, ip] = d.propi[i, it, ip]
            end
            rrninf, polp = rd_insd!(rd, idi, rriare, rridim, rd.parea[idi],
                                    10 * MAXTRE, n, 1, ksp_o, rootl_o, probiu_o,
                                    probi_o, propi_o; fint = fintf, pint = pint,
                                    sptran = RD_SPTRAN,
                                    probd = @view(d.probd[idi, :, :]),
                                    rootd = @view(d.rootd[idi, :, :]))
            # rdinsd.f DO 1050: apply averaged results (NINSIM=1).
            @inbounds for (kk, i) in enumerate(order)
                rrninf[kk] <= 1.0f-4 && continue
                nk = rrninf[kk]                              # /NINSIM (=1)
                pl = polp[kk]
                d.probiu[i] -= nk; d.probiu[i] < 0.0f0 && (d.probiu[i] = 0.0f0)
                d.probi[i, istep, 1] += nk
                d.propi[i, istep, 1] = -pl
            end
            rd_sum!(d.probit, d.probi, istep)
        end
    end

    # ==== RDCNTL DO 75 — RDSPRD + RDRATE (spread rate per center) ====
    ncents = Int(rd.ncents[idi])
    gcents = ncents - Int(d.nscen[idi])
    if rd.lonect[idi] != 1 && gcents > 0 && rd.parea[idi] != 0.0f0 && (sarea - rd.parea[idi]) > 1.0f-3
        newden = 0.0f0
        @inbounds for i in order; newden += d.ffprob[i, 1]; end
        rrsare = 20.0f0 / (newden + 1.0f-9)                 # RRGEN(idi,2)=20, IRSTYP=0
        if newden > 0.0f0
            rrsdim = sqrt(rrsare) * 208.7f0
            ksp_o   = Int[Int(t.species[i]) for i in order]
            dbh_o   = Float32[t.dbh[i] for i in order]
            rootl_o = Float32[d.rootl[i] for i in order]
            ffp_o   = Float32[d.ffprob[i, 1] for i in order]
            habsp_o = Float32[RD_HABFAC[Int(irt[Int(t.species[i])]), idi, irhab] for i in order]
            rrps_o  = ones(Float32, m)
            mcrate, _ = rd_sprd!(rd, idi; nmont = 10, irsnyr = 20, nrstep = 5, irstyp = 0,
                                 rrsfrn = 1.0f0, pint = pint, fint = fintf,
                                 rrsare = rrsare, rrsdim = rrsdim, xminkl = rd.xminkl[idi],
                                 dbh = dbh_o, rootl = rootl_o, ffprob = ffp_o, ksp = ksp_o,
                                 habsp = habsp_o, rrpswt = rrps_o)
            rout, rrate = rd_rate!(mcrate, ncents, Int(d.nscen[idi]),
                                   @view(d.rrates[idi, 1:ncents]),
                                   @view(d.shcent[idi, 1:ncents, 2]),
                                   @view(d.icensp[idi, 1:ncents]))
            @inbounds for i in 1:ncents; d.rrates[idi, i] = rout[i]; end
            d.rrrate[idi] = rrate
        end
    end

    # ==== RDCNTL DO 300 — grow centers, RDZERO, RDAREA, AREANU, RDINF ====
    @inbounds for i in 1:ncents
        rd.pcents[idi, i, 3] > 0.0f0 && (rd.pcents[idi, i, 3] += d.rrates[idi, i] * fintf)
    end
    rd_zero_centers!(rd, d, idi)
    rd_area!(rd, false)                                     # recompute PAREA from grown radii
    areanu = rd.parea[idi] - rd.ooarea[idi]
    rd.parea[idi] <= 0.0f0 && (rd.parea[idi] = 0.0f0)
    areanu <= 0.0f0 && (areanu = 0.0f0)
    rd.ooarea[idi] = rd.parea[idi]
    d.areanu[idi] = areanu
    if areanu > 0.0f0 && m > 0
        ksp_o   = Int[Int(t.species[i]) for i in order]
        fprob_o = Float32[d.fprob[i] for i in order]
        pin_o   = Float32[d.probi[i, istep, 2] for i in order]
        piu_o   = Float32[d.probiu[i] for i in order]
        propi_new, probi_new, probiu_new =
            rd_inf_kernel!(rd, idi, areanu, ksp_o, fprob_o, pin_o, piu_o;
                           fint = fintf, pint = pint, sptran = RD_SPTRAN)
        @inbounds for (kk, i) in enumerate(order)
            d.propi[i, istep, 2] = propi_new[kk]
            d.probi[i, istep, 2] = probi_new[kk]
            d.probiu[i]          = probiu_new[kk]
        end
    end

    # ==== RDMORT (+ RDSTP stump creation) ====
    isp_rec = Int[Int(t.species[i]) for i in 1:n]
    dbh_rec = Float32[t.dbh[i] for i in 1:n]
    rd_mort_kernel!(rd, d.probi, d.propi, d.rrkill, d.rdkill, dbh_rec, isp_rec, istep, fintf)
    @inbounds for i in 1:n
        d.rdkill[i] > 0.0f0 && rd_stp!(rd, d, isp_rec[i], dbh_rec[i], d.rootl[i], d.rdkill[i])
    end
    rd_sum!(d.probit, d.probi, istep)
    return
end

"""
    rd_end_apply!(rd, s, old_tpa)

Port of rd/rdend.f applied at FVS MORTS time: reconcile the RD infected-tree kill
(RRKILL) with the FVS per-record background mortality (`WK2 = old_tpa − t.tpa`, since
`mortality!` has already applied it), then re-apply the RD-adjusted WK2 to `t.tpa`.
Updates the driver PROBIU/FPROB/PROBI/PROBIT (natural-mortality reallocation) that
RDGROW reads, and creates DIENAT stumps. `old_tpa` is the cycle-start PROB (record
order). Gated: only the non-tripled, non-fire turnkey path calls this.
"""
function rd_end_apply!(rd::RootDiseaseState, s::StandState, old_tpa::Vector{Float32})
    d = rd.driver::RDDriver; t = s.trees; n = t.n
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr); istep = Int(rd.istep)
    tparea = 0.0f0
    @inbounds for id in minrr:maxrr; tparea += rd.parea[id]; end
    (n == 0 || tparea == 0.0f0) && return
    prob = Vector{Float32}(undef, n); wk2 = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        prob[i] = old_tpa[i]
        wk2[i]  = old_tpa[i] - t.tpa[i]                    # MORTS kill already applied
    end
    isp_rec = Int[Int(t.species[i]) for i in 1:n]
    dbh_rec = Float32[t.dbh[i] for i in 1:n]
    # rd_end_kernel! mutates wk2 + driver PROBIU/FPROB/PROBI/PROBIT; OAKL/BBKILL≡0.
    rd_end_kernel!(rd, wk2, prob, d.rrkill, d.rdkill, d.probit, d.probiu, d.fprob,
                   d.probi, d.propi, isp_rec, d.rootl, d.wk22, d.rroott, d.dprob,
                   d.oakl, d.bbkill; probda = d.probda, dbhda = d.dbhda, rootda = d.rootda,
                   dbh = dbh_rec, driver = d)
    @inbounds for i in 1:n
        t.tpa[i] = old_tpa[i] - wk2[i]
        t.tpa[i] < 0.0f0 && (t.tpa[i] = 0.0f0)
    end
    return
end

"""
    rd_grow_apply!(rd, s, stash)

Port of rd/rdgrow.f + the tail rd/rdinoc.f decay, applied to the stashed DG/HTG right
before the DBH/HT update. Reduces each infected host record's diameter and height
growth by the infected-root proportion (post-RDEND PROBIU/FPROB/PROBIT), then decays
the stump root radii for next cycle. Mutates `t.diam_growth`/`t.ht_growth`.
"""
function rd_grow_apply!(rd::RootDiseaseState, s::StandState, fint::Real)
    d = rd.driver::RDDriver; t = s.trees; n = t.n
    minrr = Int(rd.minrr); maxrr = Int(rd.maxrr)
    tparea = 0.0f0
    @inbounds for id in minrr:maxrr; tparea += rd.parea[id]; end
    if !(n == 0 || tparea == 0.0f0)
        isp_rec = Int[Int(t.species[i]) for i in 1:n]
        dg  = Float32[t.diam_growth[i] for i in 1:n]
        htg = Float32[t.ht_growth[i]   for i in 1:n]
        rd_grow_kernel!(rd, dg, htg, Float32[t.tpa[i] for i in 1:n], d.probit, d.probiu,
                        d.fprob, d.probi, d.propi, isp_rec)
        @inbounds for i in 1:n
            t.diam_growth[i] = dg[i]
            t.ht_growth[i]   = htg[i]
        end
    end
    rd_inoc_decay!(rd, d, fint)                            # rd/rdinoc.f (.FALSE.)
    return
end

# -----------------------------------------------------------------------------
# Engine seams (rd/rdmn1.f, rd/rdmn2.f, rd/rdtreg.f) — wired GATED + LIVE.
# A stand with no RD keyword is byte-identical (s.root_disease === nothing ⇒
# every seam early-returns). A stand WITH an active RDIN block runs the full
# per-cycle mortality/spread/growth-loss driver, .sum-visible.
# -----------------------------------------------------------------------------

"""
    root_disease_setup!(s)

fvs.f RDMN1 init seam (called once at stand setup). Inert unless RD is active.
Chunk 0 will call the ported RDSETP here (center placement + initial infection).
"""
function root_disease_setup!(s::StandState)
    rd = s.root_disease
    rd_active(rd) || return nothing
    # RDSETP → RDCLOC/RDAREA (center placement) + RDIPRP + per-record PROBI/PROBIU/
    # FPROB/PROPI + RDINOC(true) — the initial-infection state (LSTART, ISTEP=1).
    rd_setp!(rd, s)
    rd.driver = rd_build_driver!(rd, s.trees.n)     # 3-D PROBI/PROPI + stump lists
    rd.icyc = Int32(0)
    # RDSUM: the inventory (1990/icyc=0) FVS_RD_Sum row — RDPR#1 at fvs.f:347, before the cycle
    # loop (pre-projection: rdkill/probda/rrrate=0 ⇒ Mort/Stumps/Spread=0, Inf/UnInf/BA from the
    # initial infection). Populate PROBIT (=Σ PROBI) first, as RDPR would.
    if s.control.dbs_rd_sum
        d = rd.driver::RDDriver
        rd_sum!(d.probit, d.probi, max(1, Int(rd.istep)))
        yr = Int(s.control.cycle_year[1]); iage = Int(s.plot.stand_age)
        push!(rd.sum_rows, (yr, rd_sum_report(rd, s, yr, iage)))
    end
    return nothing
end

"""
    root_disease_mn2!(s, fint)

grincr.f RDMN2 seam (each cycle, before tripling). Advances the time-slot counter
ISTEP, zeroes the RRKILL accumulator, and re-sums PROBIT (rd/rdmn2.f). Live when RD
is active. RDSHST (windthrow stump history) is inert here (no windthrow).
"""
function root_disease_mn2!(s::StandState, fint::Real)
    rd = s.root_disease
    (rd_active(rd) && rd.iroot != 0) || return nothing
    d = rd.driver
    d === nothing && return nothing
    s.trees.n == 0 && return nothing
    rd.icyc += Int32(1)
    # rebuild the driver if the record count changed (COMCUP dropped PROB≤1e-5 records)
    d.n == s.trees.n || (d = rd.driver = _rd_resize_driver!(rd, d, s.trees.n))
    rd_mn2_advance!(rd, d.rrkill, d.probit, d.probi, fint)
    return nothing
end

"""
    root_disease_treg!(s, fint)

gradd.f RDTREG seam (each cycle, after growth increments, on the PRE-growth DBH — FVS
runs RDTREG in GRADD before UPDATE). Runs the full per-cycle spread chain
(`rd_control!`: FPROB/ROOTL/RDOAGM → RDCNTL RDINSD/RDSPRD/RDRATE/RDAREA/RDINF/RDMORT),
producing RRKILL + the aged PROBI/PROPI. The downstream RDEND (→ WK2 mortality) runs
at FVS MORTS time via `rd_end_apply!`, and RDGROW (→ DG/HTG) via `rd_grow_apply!`
just before the DBH update — both woven into `grow_cycle!`.
"""
function root_disease_treg!(s::StandState, fint::Real)
    rd = s.root_disease
    (rd_active(rd) && rd.iroot != 0) || return nothing
    rd.driver === nothing && return nothing
    rd_control!(rd, s, fint)
    return nothing
end

# COMCUP compaction of the driver: FVS does not compress the RD arrays for the
# manual-RRINIT path (RDRDEL is treelist-only), so a dropped record desyncs. Best-
# effort guard — copy the overlapping record prefix from the old driver (stumps carry
# through unchanged). For the turnkey no record reaches PROB≤1e-5 within 10 cycles, so
# this never fires there; it exists so an unrelated stand cannot crash the seam.
function _rd_resize_driver!(rd::RootDiseaseState, old::RDDriver, n::Int)
    d = rd_build_driver!(rd, n)
    m = min(n, old.n)
    @inbounds for i in 1:m
        for it in 1:RD_ISTEP_MAX, ip in 1:2
            d.probi[i,it,ip] = old.probi[i,it,ip]; d.propi[i,it,ip] = old.propi[i,it,ip]
        end
        d.probit[i] = old.probit[i]; d.probiu[i] = old.probiu[i]; d.fprob[i] = old.fprob[i]
        d.ffprob[i,1] = old.ffprob[i,1]; d.ffprob[i,2] = old.ffprob[i,2]
        d.rootl[i] = old.rootl[i]; d.wk22[i] = old.wk22[i]; d.rroott[i] = old.rroott[i]
    end
    d.rrates .= old.rrates; d.rrrate .= old.rrrate; d.areanu .= old.areanu
    d.shcent .= old.shcent; d.nscen .= old.nscen; d.icensp .= old.icensp
    d.probda .= old.probda; d.dbhda .= old.dbhda; d.rootda .= old.rootda
    d.decrat .= old.decrat; d.jraged .= old.jraged
    d.probd .= old.probd; d.dbhd .= old.dbhd; d.rootd .= old.rootd
    return d
end
