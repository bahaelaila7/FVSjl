# =============================================================================
# Mountain Pine Beetle (LPMPB) impact model — lpmpb/*.f
# =============================================================================
# Fifth insect/pathogen event-extension port (after DFB/DFTM/WPBR/WWPB), and the
# first STAND-LEVEL beetle model that is end-to-end .sum-validatable (its own
# per-cycle MPBGO gate → MPBCUP → COLDRV/MPBDRV, uses stand ACRES, no PPE harness).
#
# FVS structure (lpmpb/, 54 Fortran files) — the DEFAULT "rate of loss" path:
#   * MPBINT  (mpbint.f)  — one-time init of MPBCOM/COLCOM defaults + RNG seed
#                           55329 (MPRANN). Ported as `mpb_defaults!`.
#   * MPBIN   (mpbin.f)   — the `MPB … END` keyword-block reader (keywds.f 'MPB',
#                           40 sub-keywords). Ported as `kw_mpbin!`.
#   * MPBGO   (mpbgo.f)   — per-cycle outbreak gate (grincr → MPBCUP). Fires when
#                           MPBER minimum conditions are met AND an outbreak is
#                           scheduled this cycle (OPFIND activity 555) — the
#                           MANUAL/MPBSTART DETERMINISTIC path; RANSTART draws
#                           MPRANN < PROTBK (stochastic).
#   * MPBER   (mpber.f)   — DETERMINISTIC stand stats (CNTLP/BALPP/A45DBH/TLP45/
#                           PBALPP) + minimum-condition flag NOER. `mpb_er`.
#   * MPBCUP  (mpbcup.f)  — dispatch: LPOPDY ? MPBDRV : COLDRV. DEFAULT LPOPDY=F
#                           ⇒ COLDRV (Cole 1983 rate-of-loss model).
#   * COLDBH  (coldbh.f)  — DETERMINISTIC START[1..10] = LP trees/ac by DBH class.
#                           `mpb_coldbh_start`.
#   * COLIND  (colind.f)  — DBH → size class 1..10 (IFIX/2 + MOD). `mpb_colind`.
#   * COLMOD  (colmod.f)  — DETERMINISTIC epidemic loop over NUMYRS: initial
#                           DEAD=START·ZINMOR; then GREEN/DEAD recurse via the
#                           Q-value survival PRNOIN^DEAD. `mpb_colmod`. No RNG.
#   * COLMRT  (colmrt.f)  — PRKILL[i]=(START-GREEN(NUMYRS))/START; per LP record
#                           XT=PRKILL(colind(DBH))·PROB into WK2 (max-combine,
#                           PROB−1e-6 cap). `mpb_prkill` + `mpb_colmrt!`.
#   * MPRANN  (mprann.f)  — the model's OWN double-precision MINSTD LCG, seed
#                           55329 (identical algorithm to DFBRAN but seeded直接,
#                           NO +1128 reset): S1=DMOD(16807·S0,2147483647);
#                           SEL=S1/2^31. Entries MPRNSD/MPRNGT/MPRNPT. Only drawn
#                           in MPBGO's RANSTART branch. NEVER FFI'd. `mpb_rand!`.
#   * mpblkd<v>.f          — per-variant BLOCK DATA. IDXLP (FVS species # of
#                           lodgepole pine) = 7 for every lpmpb variant EXCEPT
#                           CR (=11). Variants linking real lpmpb: BM CI CR EC EM
#                           IE SO TT UT (+ generic). Every other variant links the
#                           base/exmpb.f NO-OP stub (LPMPB inert there).
#
# VALIDATION (scratchpad/lpmpb, relinked FVSie_lpmpb oracle, g16 dump-replay):
#   The DEFAULT deterministic COL path is BIT-EXACT (Float32) on the IE lodgepole
#   host stand lp_on.key (MPBSTART cycle 1, POPDYN off): COLDBH START 10/10,
#   COLMOD GREEN(NUMYRS=10) 10/10 (incl. the full 9-iteration PRNOIN^DEAD pow
#   chain), COLMRT PRKILL 10/10, COLIND 18/18, per-record XT 18/18. Instrumented
#   .sum data rows byte-identical to the clean relink. The end-to-end .sum-DELTA
#   (on−off) is CORNERED by the IE #206 growth straddle: the COL math is exact on
#   equal inputs; the absolute START/DBH it bins scale with the straddling cyc-1
#   growth (mirrors the DFB verdict).
#
# DETERMINISTIC vs STOCHASTIC split:
#   DETERMINISTIC (dump-replay bit-exact): COLIND, COLDBH, COLMOD (default
#       IBOUSE=0), COLMRT, MPBER, and the MANUAL/MPBSTART gate (OPFIND 555).
#   STOCHASTIC (MPRANN realization): only the RANSTART stand-inclusion draw in
#       MPBGO (MPRANN < PROTBK, MPOTPR probability). The COL mortality itself is
#       fully deterministic once the outbreak fires.
#
# Float32 discipline: LPMPB's Fortran is REAL. Every result here is Float32 so the
# port is bit-identical to the oracle (the 0.005454154 BAF, accumulation order,
# and IFIX truncation reproduced exactly).
# =============================================================================

const MPB_BAF = 0.005454154f0   # basal-area factor (mpber.f)

# MPBINT (mpbint.f) Cole rate-of-loss defaults — model constants, overridable by
# the INITMORT / QVALUES keywords.
const MPB_ZINMOR0 = Float32[0.0,0.0,0.0038,0.0128,0.0206,0.0353,0.0500,0.1429,0.1500,0.1500]  # XZIN
const MPB_PRNOIN0 = Float32[1.0,1.0,0.9935,0.982,0.965,0.909,0.743,0.309,0.285,0.285]         # XPRN
# IBOUSE=1 (Bousfield) bounds — only used on the NOPOPDYN <n≥1> branch (non-default).
const MPB_PMIN = Float32[.970,.950,.900,.800,.600,.500,.400,.150,.100,.050]
const MPB_PMAX = Float32[1.0,1.0,1.0,1.0,.990,.980,.970,.900,.800,.700]
const MPB_WAYNE = Float32[150.,100.,95.,66.,35.5,16.5,5.5,3.,1.5,.80]

"""
    MpbState

Mountain Pine Beetle model state (FVS `MPBCOM` + `COLCOM`). `active` mirrors
`LMPB1` (set true as soon as any MPB sub-keyword is read). The per-cycle engine
seam (`mpb_apply!`) is gated on `active`, so a stand without an MPB block projects
byte-identically.
"""
mutable struct MpbState <: AbstractMpbState
    active::Bool            # LMPB1
    debug::Bool             # DEBUIN
    lpopdy::Bool            # LPOPDY — false=Cole rate-of-loss (default), true=population dynamics (deferred)
    ibouse::Int32           # IBOUSE — 0 (default) simple Q, 1 Bousfield-bounded Q
    nclass::Int32           # NCLASS — requested classes (default 10)
    mpmxyr::Int32           # MPMXYR — max years in projection (default 10)
    lranst::Bool            # LRANST — RANSTART stochastic inclusion
    istdt::Int32            # ISTDT  — RANSTART start cycle/year (default 1)
    lepi::Bool              # LEPI   — EPIPROB user probability set
    epiprb::Float32         # EPIPRB — specified outbreak probability (default 0.5)
    prbscl::Float32         # PRBSCL — probability scaling (default 1.0)
    forlat::Float32         # FORLAT — forest latitude (default 44.0)
    lcurmr::Bool            # LCURMR — CURRMORT current-mortality data present
    linvmr::Bool            # LINVMR — INVMORT use inventoried attack data
    mpbon::Int32            # MPBON  — a MANSTART/MPBSTART outbreak was requested
    orseed::Float32         # MPRANN seed (default 55329)
    zinmor::Vector{Float32} # ZINMOR — initial mortality by class (INITMORT)
    prnoin::Vector{Float32} # PRNOIN — Q survival by class (QVALUES)
    currmr::Vector{Float32} # CURRMR — current mortality by class (CURRMORT)
    # scheduled MANUAL/MPBSTART outbreak cycles/dates (OPNEW activity 555) — matched by mpb_outbreak_due.
    outbreak_years::Vector{Int32}
    rng_s0::Float64         # MPRANN COMMON S0 — NaN until first draw, then seeded to ORSEED.
end

"""
    mpb_defaults!() → MpbState

FVS `MPBINT` (mpbint.f): the LPMPB run-time defaults, set when the first MPB
keyword is seen.
"""
function mpb_defaults!()
    return MpbState(
        false,          # active (LMPB1)
        false,          # debug
        false,          # lpopdy (COL rate-of-loss default)
        Int32(0),       # ibouse
        Int32(10),      # nclass
        Int32(10),      # mpmxyr
        false,          # lranst
        Int32(1),       # istdt
        false,          # lepi
        0.5f0,          # epiprb
        1.0f0,          # prbscl
        44.0f0,         # forlat
        false,          # lcurmr
        false,          # linvmr
        Int32(0),       # mpbon
        55329.0f0,      # orseed
        copy(MPB_ZINMOR0),
        copy(MPB_PRNOIN0),
        zeros(Float32, 10),   # currmr
        Int32[],        # outbreak_years
        NaN,            # rng_s0 (MPRANN S0) — lazy-seeded to ORSEED on first draw
    )
end

# -----------------------------------------------------------------------------
# COLIND (colind.f) — DBH → LP size class 1..10 (deterministic). VALIDATED 18/18.
# -----------------------------------------------------------------------------
"""
    mpb_colind(dbh) -> Int

FVS `COLIND`: `IVAL=IFIX(DBH)`; `INDEX=IVAL/2 + MOD(IVAL,2)`, clamped [1,10].
"""
@inline function mpb_colind(dbh::Float32)::Int
    ival = trunc(Int, dbh)              # IFIX
    idx = div(ival, 2) + mod(ival, 2)
    idx > 10 && (idx = 10)
    idx < 1  && (idx = 1)
    return idx
end

# -----------------------------------------------------------------------------
# COLDBH (coldbh.f) — START[1..10] = LP trees/ac by size class. VALIDATED 10/10.
# -----------------------------------------------------------------------------
"""
    mpb_coldbh_start(lp_dbh, lp_tpa) -> Vector{Float32}(10)

FVS `COLDBH`: for every lodgepole record (in the stand's IND1 order), add its TPA
(PROB) to `START[mpb_colind(DBH)]`. Float32 accumulation, in order.
"""
function mpb_coldbh_start(lp_dbh::AbstractVector{Float32}, lp_tpa::AbstractVector{Float32})
    start = zeros(Float32, 10)
    @inbounds for i in eachindex(lp_dbh)
        start[mpb_colind(lp_dbh[i])] += lp_tpa[i]
    end
    return start
end

# -----------------------------------------------------------------------------
# COLMOD (colmod.f) — deterministic Cole epidemic loop. VALIDATED GREEN 10/10.
# -----------------------------------------------------------------------------
"""
    mpb_colmod(start, numyrs; ibouse, zinmor, prnoin, icyc, lcurmr, linvmr, currmr, greinf) -> GREEN::Matrix{Float32}(numyrs,10)

FVS `COLMOD` (default IBOUSE=0, single-cycle ELSE branch): initial
`DEAD(1,i)=START·ZINMOR(i)`, `GREEN(1,i)=START−DEAD(1,i)`. For K=2..NUMYRS (while
TDEAD(K-1)>5e-4 and TGREEN(K-1)>5e-4): `DEAD(K,J)=GREEN(K-1,J)·(1−PNEW^DEAD(K-1,J))`
with `PNEW=PRNOIN(J)` (IBOUSE=0). Returns the full GREEN(year,class) matrix.

The CURRMORT/INVMORT ICYC=1 branch (colmod.f:93-99, active when
`icyc==1 && (lcurmr||linvmr)`): DEAD(1,i)=GREINF(i)+CURRMR(i), GREEN(1,i)=START−GREINF(i)
where GREINF is bumped up to `START·ZINMOR−CURRMR` when the natural initial mortality
exceeds the observed infested+current, then capped at START. `greinf` is treated as
read-only (a local copy is mutated, mirroring the COLCOM array). CURRMR comes from the
CURRMORT keyword; GREINF comes from MPBDAM inventory damage/severity codes (DAMCDS) which
the FVSjl treelist does NOT carry, so on loadable stands `greinf` is all-zero — and with
`currmr` also zero the branch is byte-identical to the default (verified: INVMORT with no
damage == MPBSTART). Supplying nonzero `currmr` (CURRMORT keyword) or synthetic `greinf`
exercises the divergent path (unit-tested against the FVSie_lpmpb oracle golden).
"""
function mpb_colmod(start::Vector{Float32}, numyrs::Int; ibouse::Int=0,
                    zinmor::Vector{Float32}=MPB_ZINMOR0, prnoin::Vector{Float32}=MPB_PRNOIN0,
                    icyc::Int=0, lcurmr::Bool=false, linvmr::Bool=false,
                    currmr::Vector{Float32}=zeros(Float32,10),
                    greinf::Vector{Float32}=zeros(Float32,10))
    numyrs < 1 && (numyrs = 1)
    DEAD  = zeros(Float32, numyrs, 10)
    GREEN = zeros(Float32, numyrs, 10)
    TDEAD  = zeros(Float32, numyrs)
    TGREEN = zeros(Float32, numyrs)
    use_curr = (icyc == 1 && (lcurmr || linvmr))
    gi = use_curr ? copy(greinf) : greinf     # COLCOM GREINF is mutated in place
    @inbounds for i in 1:10
        if use_curr
            if start[i] * zinmor[i] > gi[i] + currmr[i]
                gi[i] = start[i] * zinmor[i] - currmr[i]
            end
            gi[i] > start[i] && (gi[i] = start[i])
            DEAD[1,i]  = gi[i] + currmr[i]
            GREEN[1,i] = start[i] - gi[i]
        else
            DEAD[1,i]  = start[i] * zinmor[i]
            GREEN[1,i] = start[i] - DEAD[1,i]
        end
        TDEAD[1]  += DEAD[1,i]
        TGREEN[1] += GREEN[1,i]
    end
    @inbounds for k in 2:numyrs
        (TDEAD[k-1] <= 0.0005f0 || TGREEN[k-1] <= 0.0005f0) && break
        for j in 1:10
            pnew = prnoin[j]
            if ibouse == 1
                if GREEN[1,j] < 0.03f0
                    pnew = 0.0f0
                else
                    pnew = prnoin[j] ^ (MPB_WAYNE[j] / GREEN[1,j])
                end
                pnew < MPB_PMIN[j] && (pnew = MPB_PMIN[j])
                pnew > MPB_PMAX[j] && (pnew = MPB_PMAX[j])
            end
            DEAD[k,j]  = GREEN[k-1,j] * (1.0f0 - (pnew ^ DEAD[k-1,j]))
            TDEAD[k]  += DEAD[k,j]
            GREEN[k,j] = GREEN[k-1,j] - DEAD[k,j]
            TGREEN[k] += GREEN[k,j]
        end
        TDEAD[k] += TDEAD[k-1]
    end
    return GREEN
end

# -----------------------------------------------------------------------------
# COLMRT (colmrt.f) — PRKILL + per-record WK2 mortality. VALIDATED PRKILL/XT.
# -----------------------------------------------------------------------------
"""
    mpb_prkill(start, green_final) -> Vector{Float32}(10)

FVS `COLMRT` proportion killed by class: `PRKILL(i)=(START-GREEN(NUMYRS))/START`
(0 when START≤0). `green_final` = GREEN(NUMYRS,·). VALIDATED 10/10.
"""
function mpb_prkill(start::Vector{Float32}, green_final::AbstractVector{Float32})
    prkill = zeros(Float32, 10)
    @inbounds for i in 1:10
        prkill[i] = start[i] <= 0.0f0 ? 0.0f0 : (start[i] - green_final[i]) / start[i]
    end
    return prkill
end

"""
    mpb_colmrt!(t, old_tpa, lpidx, lp_dbh, lp_tpa, prkill)

FVS `COLMRT` per-record apply (mirror of DFBMRT's WK2 max-combine): for each
lodgepole record `XT = PRKILL(colind(DBH))·PROB`; `WK2 = max(WK2, XT)`, capped at
`PROB−1e-6`; the record's TPA is written back `old_tpa − WK2`. `WK2` is the current
periodic mortality (`old_tpa − t.tpa`). XT VALIDATED 18/18. Mutates `t.tpa`.
"""
function mpb_colmrt!(t, old_tpa::AbstractVector{Float32}, lpidx::AbstractVector{<:Integer},
                     lp_dbh::AbstractVector{Float32}, lp_tpa::AbstractVector{Float32},
                     prkill::Vector{Float32})
    @inbounds for k in eachindex(lpidx)
        j = lpidx[k]
        idx = mpb_colind(lp_dbh[k])
        xt = prkill[idx] * lp_tpa[k]                 # XT = PRKILL(INDEX)*PROB
        wk2 = old_tpa[j] - t.tpa[j]                  # current MORTS periodic mortality
        xt > wk2 && (wk2 = xt)                        # IF (XT .GT. WK2) WK2 = XT
        (old_tpa[j] - wk2 < 1.0f-6) && (wk2 = old_tpa[j] - 1.0f-6)   # PROB−WK2<1e-6 cap
        t.tpa[j] = old_tpa[j] - wk2
        t.tpa[j] < 0.0f0 && (t.tpa[j] = 0.0f0)
    end
    return nothing
end

# -----------------------------------------------------------------------------
# MPBER (mpber.f) — LP stand stats + minimum-outbreak-condition flag. Deterministic.
# -----------------------------------------------------------------------------
"""
    mpb_er(lp_dbh, lp_tpa, ba) -> (cntlp, balpp, a45dbh, tlp45, pbalpp, noer)

FVS `MPBER`. `lp_*` = lodgepole records in IND1 order; `ba` = stand basal area
(PLOT BA). NOER (minimum conditions met) = there is ≥1 LP record, CNTLP≥0.01, and
TLP45≥1 (≥1 TPA of LP ≥4.5″). All Float32, accumulation in order.
"""
function mpb_er(lp_dbh::AbstractVector{Float32}, lp_tpa::AbstractVector{Float32}, ba::Float32)
    cntlp = 0.0f0; balpp = 0.0f0; a45dbh = 0.0f0; tlp45 = 0.0f0
    isempty(lp_dbh) && return (cntlp=cntlp, balpp=balpp, a45dbh=a45dbh, tlp45=tlp45, pbalpp=0.0f0, noer=false)
    @inbounds for ii in eachindex(lp_dbh)
        p = lp_tpa[ii]; d = lp_dbh[ii]
        cntlp += p
        balpp += MPB_BAF * d * d * p
        if d >= 4.5f0
            tlp45  += p
            a45dbh += d * p
        end
    end
    (cntlp < 0.01f0 || tlp45 < 1.0f0) && return (cntlp=cntlp, balpp=balpp, a45dbh=a45dbh, tlp45=tlp45, pbalpp=0.0f0, noer=false)
    pbalpp = balpp / ba
    a45dbh = a45dbh / tlp45
    return (cntlp=cntlp, balpp=balpp, a45dbh=a45dbh, tlp45=tlp45, pbalpp=pbalpp, noer=true)
end

# -----------------------------------------------------------------------------
# MPRANN (mprann.f) — the LPMPB model's own double-precision MINSTD LCG. NEVER FFI.
# -----------------------------------------------------------------------------
"""
    mpb_rand!(m) -> Float32

FVS `MPRANN`: `S1=DMOD(16807·S0,2147483647)`; `SEL=REAL(S1/2147483648)`; `S0=S1`.
Double-precision modular step (exact); only the returned uniform is Float32.
Lazy-seeds `S0=ORSEED` (default 55329) the first time — MPBINT seeds it directly
(no DFBSCH +1128). Only consumed by the RANSTART branch of MPBGO.
"""
@inline function mpb_rand!(m::MpbState)::Float32
    isnan(m.rng_s0) && (m.rng_s0 = Float64(m.orseed))
    s1 = rem(16807.0 * m.rng_s0, 2147483647.0)
    m.rng_s0 = s1
    return Float32(s1 / 2147483648.0)
end

# -----------------------------------------------------------------------------
# MPOTPR (mpotpr.f) — deterministic MPB outbreak PROBABILITY (RANSTART path).
# -----------------------------------------------------------------------------
# glibc single-precision expf — matches gfortran REAL EXP() bit-exact. The Fortran
# MPOTPR expression is all-REAL, so the transcendental is expf, not exp(::Float64).
@inline _mpb_expf(x::Float32)::Float32 = ccall(:expf, Float32, (Float32,), x)

"""
    mpb_mpotpr(pbalpp, relden, reldsp_lp, a45dbh, cntlp, istdt, icyc, iy) -> Float32

FVS `MPOTPR`: the probability of a mountain-pine-beetle outbreak, used only by the
RANSTART branch of MPBGO. `PROTBK=0` unless ALL minimum conditions hold (at/after the
RANSTART start date; average DBH of ≥4.5″ LP ≥6.0; ≥25% of stand BA is LP; ≥20% of
stand relative-density (CCF) is LP; ≥40 LP TPA). Otherwise the logistic
`1/(1+exp(9.583 − 0.08967·(X·RELDEN)))` with `X=min(PBALPP,0.8)` and `RELDEN` the
stand relative density (total CCF). `icyc` is the 1-based FVS cycle, `iy` its calendar
year, `istdt` the RANSTART start date (≤40 ⇒ cycle number, >40 ⇒ calendar year).
VALIDATED bit-exact (Float32) 5/5 vs FVSie_lpmpb dump-replay.
"""
function mpb_mpotpr(pbalpp::Float32, relden::Float32, reldsp_lp::Float32,
                    a45dbh::Float32, cntlp::Float32, istdt::Int, icyc::Int, iy::Int)::Float32
    protbk = 0.0f0
    if (istdt <= 40 && icyc < istdt) || (istdt > 40 && iy < istdt)
        return protbk
    end
    a45dbh < 6.0f0 && return protbk
    pbalpp < 0.25f0 && return protbk
    relden <= 0.0f0 && return protbk                 # guard the RELDSP/RELDEN divide
    (reldsp_lp / relden) < 0.20f0 && return protbk
    cntlp < 40.0f0 && return protbk
    x = pbalpp
    x > 0.8f0 && (x = 0.8f0)
    protbk = 1.0f0 / (1.0f0 + _mpb_expf(9.583f0 - 0.08967f0 * (x * relden)))
    return protbk
end

"""
    mpb_lp_ccf(s, idxlp) -> Float32

Lodgepole share of the stand relative density (FVS `RELDSP(IDXLP)`): the same per-tree
crown-competition factor `stand_ccf` sums for `RELDEN`, summed over lodgepole records
only. Used for the MPOTPR ≥20%-LP minimum-condition gate. Dispatches to the variant's
per-tree CCF exactly as `stand_ccf` does (CR via the crown-width→area path).
"""
function mpb_lp_ccf(s::StandState, idxlp::Int)::Float32
    t = s.trees; ccf = 0.0f0
    @inbounds for i in 1:t.n
        Int(t.species[i]) == idxlp || continue
        d = t.dbh[i]; p = t.tpa[i]
        c = if s.variant isa InlandEmpire;          ie_tree_ccf(idxlp, d)
            elseif s.variant isa Kootenai;           kt_tree_ccf(idxlp, d)
            elseif s.variant isa EasternMontana;     em_tree_ccf(idxlp, d)
            elseif s.variant isa Utah;               ut_tree_ccf(idxlp, d)
            elseif s.variant isa BlueMountains;      bm_tree_ccf(idxlp, d)
            elseif s.variant isa CentralIdaho;       ci_tree_ccf(idxlp, d)
            elseif s.variant isa Teton;              tt_tree_ccf(idxlp, d)
            elseif s.variant isa EastCascades;       ec_tree_ccf(idxlp, d)
            elseif s.variant isa SouthCentralOregon; so_tree_ccf(idxlp, d, t.height[i])
            elseif s.variant isa CentralRockies
                cw = cr_crown_width(idxlp, d, Int(s.plot.model_type))
                d > 0.1f0 ? 0.001803f0 * cw * cw : 0.001f0
            else
                0.0f0
            end
        ccf += c * p
    end
    return ccf
end

# -----------------------------------------------------------------------------
# IDXLP — the FVS species number of lodgepole pine for a lpmpb-linked variant.
# -----------------------------------------------------------------------------
"""
    mpb_idxlp(variant) -> Int

`mpblkd<v>.f` `IDXLP`: 7 for every lpmpb variant (BM/CI/EC/EM/IE/SO/TT/UT +
generic) EXCEPT CentralRockies (=11). `0` (lpmpb not linked / no host) elsewhere,
which keeps the LPMPB seam inert.
"""
function mpb_idxlp(variant)::Int
    variant isa CentralRockies && return 11
    (variant isa BlueMountains || variant isa CentralIdaho || variant isa EastCascades ||
     variant isa EasternMontana || variant isa InlandEmpire || variant isa SouthCentralOregon ||
     variant isa Teton || variant isa Utah) && return 7
    return 0
end

# -----------------------------------------------------------------------------
# MPBGO outbreak gate (mpbgo.f) — MANUAL/MPBSTART deterministic path (OPFIND 555).
# -----------------------------------------------------------------------------
"""
    mpb_outbreak_due(m, s) -> Bool

Whether a MANUAL/MPBSTART regional outbreak (OPNEW activity 555 at date IDT,
mpbin.f opt 5/6) is due this cycle — the OPFIND(555) test in MPBGO. Mirrors the
OPNEW/OPCYCL date matching: `IDT==0` every cycle; `0<IDT<1000` a 1-based cycle;
`IDT≥1000` a calendar year in this cycle's window. In the COL path MPBYR stays 0,
so the outbreak fires purely on this scheduled test (RANSTART would additionally
draw MPRANN < PROTBK — the only stochastic branch, deferred).
"""
function mpb_outbreak_due(m::MpbState, s::StandState)::Bool
    isempty(m.outbreak_years) && return false
    fvscyc = Int(s.control.cycle) + 1
    cyc0 = Int(s.control.cycle)
    cs = cycle_year_at(s.control, cyc0)
    ce = cycle_year_at(s.control, cyc0 + 1); ce <= cs && (ce = cs + 1)
    @inbounds for m0 in m.outbreak_years
        mi = Int(m0)
        mi == 0 && return true
        (0 < mi < 1000) && mi == fvscyc && return true
        (mi >= 1000) && (cs <= mi < ce) && return true
    end
    return false
end

# -----------------------------------------------------------------------------
# MPBCUP→COLDRV seam (mpbcup.f/coldrv.f: COLDBH→COLMOD→COLMRT), gated by MPBGO —
# called from the grow/mortality path (gradd.f:63, IF LMPBGO CALL MPBCUP).
# -----------------------------------------------------------------------------
"""
    mpb_apply!(s, old_tpa, fint)

The LPMPB cycle driver + gate (FVS `MPBGO` → `MPBCUP` → `COLDRV`), wired into the
FVSjl mortality path. No-op unless an MPB block is active, an outbreak is due this
cycle (MANUAL/MPBSTART), the variant has a lodgepole host (IDXLP≠0), and the stand
meets the MPBER minimum condition. When it fires it bins the cycle-start lodgepole
(COLDBH), runs the deterministic Cole epidemic (COLMOD, NUMYRS=min(MPMXYR,IFINT)≤10),
and distributes the class PRKILL to the per-record mortality (COLMRT). Byte-identical
when no outbreak fires.

Outbreak-decision paths handled (MPBGO): the deterministic MANUAL/MPBSTART schedule
(OPFIND 555 ⇒ fire), the CURRMORT/INVMORT auto-cycle-1 schedule (ICYC=1 GREINF branch in
COLMOD), and RANSTART (LRANST) — which draws `mpb_rand! < PROTBK` (PROTBK from MPOTPR,
scaled by PRBSCL, or EPIPRB when LEPI), drawing once per eligible cycle so the RNG stream
stays in sync with the oracle. The LPOPDY population-dynamics path (MPBDRV/MPBMOD) stays
DEFERRED (early-return). IBOUSE from the keyword. GREINF (live inventory-attack) is 0 on
loadable stands (no treelist damage codes); CURRMR is the CURRMORT keyword.
"""
function mpb_apply!(s::StandState, old_tpa::Vector{Float32}, fint::Real)
    m = s.mpb
    (m === nothing || !m.active) && return nothing
    m.lpopdy && return nothing                             # population-dynamics path deferred
    idxlp = mpb_idxlp(s.variant)
    idxlp == 0 && return nothing
    t = s.trees; n = t.n
    n == 0 && return nothing
    lpidx = Int[i for i in 1:n if Int(t.species[i]) == idxlp]
    isempty(lpidx) && return nothing
    lp_dbh = Float32[t.dbh[i] for i in lpidx]
    lp_tpa = Float32[old_tpa[i] for i in lpidx]
    # PLOT BA (stand basal area) from cycle-start PROB (mpber.f uses PLOT BA).
    ba = 0.0f0
    @inbounds for i in 1:n
        ba += MPB_BAF * t.dbh[i] * t.dbh[i] * old_tpa[i]
    end
    # MPBGO: MPBER minimum conditions (L=NOERR). If not met, no outbreak AND no draw.
    r = mpb_er(lp_dbh, lp_tpa, ba)
    r.noer || return nothing
    # MPBGO label-100 decision. MPBYR==0 in the COL path, so the choice is: a user-scheduled
    # outbreak this cycle (OPFIND 555, MANSTART/MPBSTART/CURRMORT) fires deterministically;
    # otherwise the RANSTART branch draws MPRANN once and fires when X < PROTBK*PRBSCL.
    fire = false
    if mpb_outbreak_due(m, s)                              # OPFIND(555): NTODO>0
        fire = true
    elseif m.lranst
        icyc = Int(s.control.cycle) + 1
        iy   = Int(cycle_year_at(s.control, Int(s.control.cycle)))
        relden    = s.plot.relative_density               # RELDEN (stand CCF)
        reldsp_lp = mpb_lp_ccf(s, idxlp)                   # RELDSP(IDXLP)
        protbk = m.lepi ? m.epiprb :
                 mpb_mpotpr(r.pbalpp, relden, reldsp_lp, r.a45dbh, r.cntlp,
                            Int(m.istdt), icyc, iy) * m.prbscl
        x = mpb_rand!(m)                                   # MPRANN draw — advance RNG each eligible cycle
        fire = x < protbk
    end
    fire || return nothing
    start = mpb_coldbh_start(lp_dbh, lp_tpa)               # COLDBH
    numyrs = min(Int(m.mpmxyr), round(Int, fint)); numyrs > 10 && (numyrs = 10)
    numyrs < 1 && (numyrs = 1)
    icyc = Int(s.control.cycle) + 1                        # 1-based FVS ICYC (CURRMORT branch fires at ICYC==1)
    green = mpb_colmod(start, numyrs; ibouse=Int(m.ibouse), zinmor=m.zinmor, prnoin=m.prnoin,
                       icyc=icyc, lcurmr=m.lcurmr, linvmr=m.linvmr, currmr=m.currmr)  # COLMOD
    prkill = mpb_prkill(start, @view green[numyrs, :])     # COLMRT PRKILL
    mpb_colmrt!(t, old_tpa, lpidx, lp_dbh, lp_tpa, prkill) # COLMRT apply
    return nothing
end

# -----------------------------------------------------------------------------
# kw_mpbin! (mpbin.f) — MPB keyword-block reader (keywds.f 'MPB', 40 sub-keywords)
# -----------------------------------------------------------------------------
"""
    kw_mpbin!(s, rec, kr)

Parse the `MPB … END` block (lpmpb/mpbin.f). Sets the MPBCOM/COLCOM-equivalent
state on `s.mpb` and consumes sub-keyword records up to `END`. `LMPB1` (active) is
set true as soon as the block is entered. The gated `mpb_apply!` seam realizes the
deterministic COL mortality; the reader is otherwise faithful to MPBIN.

Sub-keywords ported (state-setting): END, INVMORT, MANSTART, MPBSTART, DEBUG,
PRBSCALE, RANNSEED, RANSTART, NUMCLASS, MAXYEARS, POPDYN, NOPOPDYN, CURRMORT,
EPIPROB, NODEBUG, QVALUES, INITMORT, LATITUDE. The remaining pheromone/spray/
population-dynamics-only keywords (PSFOUND/DCFOUND/AGGPHERM/REPPHERM/BETTER/AMP/
ACTSRF/CRITAD/MPBECHO/BEETLES/HABSUIT/STNDSIZE/GENOTYPE/STRONG/EMERINC/AGGTHRES/
MPBGRF/NOMPBGRF/PSDBHLIM/PSPKILL/DCPKILL/PARTIAL) are recognized (no effect on the
COL mortality path); several read a trailing data record which the block reader
must skip — see the handoff notes.
"""
function kw_mpbin!(s::StandState, rec, kr::KeywordReader)
    s.mpb === nothing && (s.mpb = mpb_defaults!())
    m = s.mpb
    m.active = true                        # mpbin.f: LMPB1 = .TRUE.
    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            break
        elseif k == "INVMORT"              # opt 4
            m.linvmr = true; m.lcurmr = false
        elseif k == "MANSTART" || k == "MPBSTART"   # opt 5/6 — schedule outbreak
            m.mpbon = Int32(1)
            idt = r.present[1] ? Int32(trunc(Int, r.values[1])) : Int32(1)
            push!(m.outbreak_years, idt)
        elseif k == "DEBUG"                # opt 9
            m.debug = true
        elseif k == "PRBSCALE"             # opt 10
            r.present[1] && (m.prbscl = Float32(r.values[1]))
        elseif k == "RANNSEED"             # opt 11
            if r.present[1]
                seed = Float32(r.values[1])
                seed % 2.0f0 == 0.0f0 && (seed += 1.0f0)   # MPRNSD forces odd
                m.orseed = seed
            end
        elseif k == "RANSTART"             # opt 13
            m.lranst = true
            r.present[1] && (m.istdt = Int32(trunc(Int, r.values[1])))
        elseif k == "NUMCLASS"            # opt 14
            (r.present[1] && 1.0f0 <= r.values[1] <= 30.0f0) && (m.nclass = Int32(trunc(Int, r.values[1])))
        elseif k == "MAXYEARS"           # opt 21
            (r.present[1] && 1.0f0 <= r.values[1] <= 30.0f0) && (m.mpmxyr = Int32(trunc(Int, r.values[1])))
        elseif k == "POPDYN"             # opt 26
            m.lpopdy = true
        elseif k == "NOPOPDYN"           # opt 27
            m.lpopdy = false
            r.present[1] && (m.ibouse = r.values[1] >= 1.0f0 ? Int32(1) : Int32(0))
        elseif k == "CURRMORT"           # opt 29 — reads a trailing record (see handoff)
            @inbounds for i in 1:7
                r.present[i] && (m.currmr[i] = Float32(r.values[i]))
            end
            m.lcurmr = true; m.linvmr = false
            push!(m.outbreak_years, Int32(1))   # OPNEW(555) cycle 1
        elseif k == "EPIPROB"            # opt 30
            if r.present[1] && 0.0f0 <= r.values[1] <= 1.0f0
                m.lepi = true; m.epiprb = Float32(r.values[1])
            end
        elseif k == "NODEBUG"            # opt 32
            m.debug = false
        elseif k == "QVALUES"           # opt 33 — reads a trailing record (see handoff)
            @inbounds for i in 1:7
                r.present[i] && (m.prnoin[i] = Float32(r.values[i]))
            end
        elseif k == "INITMORT"          # opt 34 — reads a trailing record (see handoff)
            @inbounds for i in 1:7
                r.present[i] && (m.zinmor[i] = Float32(r.values[i]))
            end
        elseif k == "LATITUDE"          # opt 39
            r.present[1] && (m.forlat = Float32(r.values[1]))
        elseif k in ("PSFOUND","DCFOUND","AGGPHERM","REPPHERM","DEBUG","PRBSCALE",
                     "BETTER","AMP","ACTSRF","CRITAD","MPBECHO","BEETLES","HABSUIT",
                     "STNDSIZE","GENOTYPE","STRONG","EMERINC","AGGTHRES","MPBGRF",
                     "NOMPBGRF","PSDBHLIM","PSPKILL","DCPKILL","PARTIAL")
            # recognized; no effect on the COL rate-of-loss mortality path.
        else
            (!isempty(k) && isletter(first(k))) && push!(s.control.unrecognized_keywords, k)
        end
    end
    return nothing
end
