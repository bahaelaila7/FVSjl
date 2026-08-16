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
    return rd
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
const RD_HABFAC = reshape(RD_HABFAC_FLAT, RD_ITOTSP, RD_ITOTRR, 2)   # HABFAC(ksp,idi,ihab)
const RD_PNINF  = reshape(RD_PNINF_FLAT,  RD_ITOTSP, RD_ITOTRR)      # PNINF(ksp,idi)
const RD_PKILLS = reshape(RD_PKILLS_FLAT, RD_ITOTSP, RD_ITOTRR)      # PKILLS(ksp,idi)
# RRPSWT (rd/rdinit.f) defaults to 1.0 for every species; only the RRPSWT keyword
# (not in the chunk-0b-2 turnkey path) changes it, so the const default is faithful here.
const RD_RRPSWT = ones(Float32, RD_ITOTSP)

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
# Engine seams (rd/rdmn1.f, rd/rdmn2.f, rd/rdtreg.f) — wired GATED + INERT.
# The per-cycle mortality/spread/growth-loss bodies are Chunk 0 (pending); until
# then these are no-ops, so a stand with no RD keyword is byte-identical and a
# stand WITH an RD keyword parses/initializes but applies no mortality yet.
# -----------------------------------------------------------------------------

"""
    root_disease_setup!(s)

fvs.f RDMN1 init seam (called once at stand setup). Inert unless RD is active.
Chunk 0 will call the ported RDSETP here (center placement + initial infection).
"""
function root_disease_setup!(s::StandState)
    rd = s.root_disease
    rd_active(rd) || return nothing
    # Chunk 0b-2: RDSETP → RDCLOC/RDAREA (center placement) + RDIPRP + per-record
    # PROBI/PROBIU/FPROB/PROPI + RDINOC(true). Sets the initial-infection state; the
    # per-cycle mortality driver (RDMORT/RDGROW) that makes this .sum-visible is 0b-3.
    rd_setp!(rd, s)
    return nothing
end

"""
    root_disease_mn2!(s, fint)

grincr.f RDMN2 seam (each cycle, before tripling). Advances the time-slot counter.
Inert (no mortality) until Chunk 0.
"""
function root_disease_mn2!(s::StandState, fint::Real)
    rd = s.root_disease
    (rd_active(rd) && rd.iroot != 0) || return nothing
    # Chunk 0 (pending): ISTEP+=1; IYEAR+=fint; zero RRKILL; RDSUM; RDSHST.
    return nothing
end

"""
    root_disease_treg!(s, fint)

gradd.f RDTREG seam (each cycle, after growth). The main per-cycle RD driver:
Chunk 0 will run RDCNTL (RDSETP→spread→RDMORT) + RDEND (→WK2 mortality) + RDGROW
(growth loss) + RDINOC. Inert until then.
"""
function root_disease_treg!(s::StandState, fint::Real)
    rd = s.root_disease
    (rd_active(rd) && rd.iroot != 0) || return nothing
    # Chunk 0 (pending): the reduced RDCNTL→RDMORT mortality path.
    return nothing
end
