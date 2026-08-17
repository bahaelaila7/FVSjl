# =============================================================================
# Douglas-fir Beetle (DFB) impact model — dfb/*.f
# =============================================================================
# PILOT port of the FVS insect/pathogen event-extension family (roadmap item 2).
#
# FVS structure (dfb/, 29 Fortran files):
#   * DFBINT  (dfbint.f)  — one-time init of the DFBCOM defaults + RNG seed
#                           (called from INITRE); ported here as `dfb_defaults!`.
#   * DFBIN   (dfbin.f)   — the DFBEETLE keyword-block reader (initre option 100,
#                           analogous to RDIN option 101); ported as `kw_dfbin!`.
#   * DFBGO   (dfbgo.f)   — per-cycle "does an outbreak fire in this stand this
#                           cycle" gate (grincr.f). STOCHASTIC unless MANSTART.
#   * DFBDRV  (dfbdrv.f)  — cycle driver: DFBDBH → DFBMOD → DFBMRT (gradd.f).
#   * DFBDBH  (dfbdbh.f)  — DETERMINISTIC: build START[] = DF trees/ac per size
#                           class. Ported as `dfb_dbh_start`.
#   * DFBER   (dfber.f)   — DETERMINISTIC: stand statistics BA9/BADF9/A45DBH/
#                           PBADF4 and the min-condition flag. Ported as `dfb_er`.
#   * DFBPRB  (dfbprb.f)  — DETERMINISTIC: stand-outbreak probability PROTBK.
#                           Ported as `dfb_prb`.
#   * DFBIND  (dfbind.f)  — DETERMINISTIC: DBH → size class 1..20. `dfb_ind`.
#   * DFBMOD  (dfbmod.f)  — STOCHASTIC: DFKILL = BACHLO(EXPCTD,EXSTDV,DFBRAN)*…
#                           (a normal draw) + the CUROUTBK/LINPRG in-progress branch. `dfb_mod`.
#   * DFBMRT  (dfbmrt.f)  — applies DFKILL to WK2 (mortality). Depends on DFBMOD.
#   * DFBSCH  (dfbsch.f)  — RANSCHED regional-outbreak auto-scheduler. `dfb_schedule!`.
#   * DFBWIN  (dfbwin.f)  — WINDTHR windthrow mortality + OKILL feed. `dfb_win!`/`dfb_win_kernel!`.
#   * DFBINV  (dfbinv.f)  — CUROUTBK pre-killed-DF count from the treelist. `dfb_inv!`.
#   * DFBRAN  (dfbran.f)  — the model's OWN Lehmer/MINSTD LCG (NOT the FVS ZRAND
#                           stream): S1 = mod(16807·S0, 2147483647); u = S1/2^31.
#                           seed default ORSEED = 55329. Entries DFBNSD/DFBGSD/
#                           DFBCSD. Reimplemented in Julia when DFBMOD is ported;
#                           NEVER FFI'd.
#   * dfblkd<v>.f          — per-variant BLOCK DATA (WINSUC windthrow suscept.,
#                           IFVSSP crosswalk, PERDD, ROWDOM, IDFSPC = FVS species
#                           number of Douglas-fir). Variants with real dfb: BC,
#                           BM, CI, CR, EC, IE, PN, TT, UT (+ generic dfblkd.f).
#                           Every other variant links the base/exdfb.f NO-OP stub.
#
# SCOPE (complete): the keyword reader + the deterministic routines + the full mortality path
# (DFBRAN/BACHLO/DFBMOD/DFBMRT) + all outbreak-activation modes — MANSTART/MANSCHED, RANSTART
# (stochastic stand inclusion), RANSCHED (DFBSCH auto-scheduler), CUROUTBK (LINPRG in-progress
# outbreak + DFBINV), and DFBWIN (windthrow). Every numeric path is validated BIT-EXACT (Float32)
# by dump-replay against a relinked FVSie_dfb g16 oracle (scratchpad/dfb; instrumented .sum
# byte-identical to the clean relink first), plus a PN variant-sweep smoke A/B. The end-to-end
# .sum-DELTA is CORNERED by the documented IE growth straddle (the DFB math is exact on equal
# inputs; the absolute BADF9/BA9/PCT/HT the modes read scale with the straddling cyc-1 growth).
#
# DETERMINISTIC vs STOCHASTIC split:
#   DETERMINISTIC (dump-replay bit-exact): DFBIND, DFBDBH(START), DFBER(BA9, BADF9, A45DBH, PBADF4,
#       LMIN), DFBPRB(PROTBK), the MANSTART/MANSCHED gate, DFBSCH's schedule, DFBWIN's windthrow
#       kill, and the CUROUTBK user-PREKLL branch.
#   STOCHASTIC (DFBRAN realization, .sum-DELTA straddle): the RANSTART inclusion draw (DFBGO,
#       precedes BACHLO), the DFKILL magnitude (DFBMOD via BACHLO→DFBRAN), the CUROUTBK BACHLO
#       branch, and RANSCHED's regional scheduling draws (seeded at ORSEED, not +1128).
#
# Float32 discipline: DFB's Fortran is REAL (Float32). Every arithmetic result
# here is Float32 so the port is bit-identical to the oracle. The 0.005454154
# basal-area-factor constant, the accumulation order, and IFIX truncation are
# reproduced exactly.
# =============================================================================

const DFB_BAF = 0.005454154f0   # basal-area factor: BA = DFB_BAF·DBH²·TPA (dfber.f/dfbmrt.f)
# PERDD (dfblkd*.f) — cumulative fraction of the 4-year outbreak's kill realized by outbreak-year
# IYOUT. Identical across every DFB-linked variant (a model constant, not variant data). Used by the
# CUROUTBK/LINPRG in-progress-outbreak branch of DFBMOD (1-based; IYOUT ∈ 1..4).
const DFB_PERDD = (0.30f0, 0.80f0, 0.95f0, 1.00f0)

"""
    DfbState

Douglas-fir Beetle model state (FVS `DFBCOM`). Populated by `kw_dfbin!` and
seeded with the `DFBINT` defaults by `dfb_defaults!`. `active` mirrors `LDFBON`
(set true as soon as any DFB sub-keyword is read). The per-cycle engine seams
(`dfb_win!`/`dfb_apply!`) and the pre-loop `dfb_setup!` are gated on `active`, so
a stand without a DFB block projects byte-identically.
"""
mutable struct DfbState <: AbstractDfbState
    active::Bool        # LDFBON  — DFB model is in use
    debug::Bool         # DEBUIN  — DEBUG keyword
    lbamod::Bool        # LBAMOD  — mortality distribution: false=DBH (default), true=BA (MORTDIS)
    ismeth::Int32       # ISMETH  — 1=MANSTART, 2=RANSTART (default 2)
    idbsch::Int32       # IDBSCH  — 1=MANSCHED, 2=RANSCHED (default 2)
    ilenth::Int32       # ILENTH  — outbreak length, years (default 4)
    iwait::Int32        # IWAIT   — min wait between regional outbreaks (default 10)
    ipast::Int32        # IPAST   — year of last recorded outbreak (default 1950)
    iyout::Int32        # IYOUT   — CUROUTBK: years an in-progress outbreak has run
    expctd::Float32     # EXPCTD  — expected DF trees killed per outbreak year (default 6.0)
    exstdv::Float32     # EXSTDV  — std-dev for the DFKILL normal draw (default 2.0)
    epiprb::Float32     # EPIPRB  — STOPROB user-specified inclusion probability (default 0.5)
    lepi::Bool          # LEPI    — STOPROB set
    prpwin::Float32     # PRPWIN  — WINDTHR proportion (default 0.8)
    minden::Float32     # MINDEN  — WINDTHR min stems/ac (default 0.0)
    mwinht::Float32     # MWINHT  — WINDTHR min tree height (default 20.0)
    okill::Float32      # OKILL   — extra windthrow mortality carried into DFB (default 0.0)
    orseed::Float32     # ORSEED  — DFBRAN seed (default 55329)
    dbevnt::Float32     # DBEVNT  — annual regional-outbreak probability (default 0.05)
    prekll::Float32     # PREKLL  — TPA already killed before model start (CUROUTBK/treelist)
    linprg::Bool        # LINPRG  — CUROUTBK: outbreak in progress at start
    linv::Bool          # LINV    — mortality-in-progress data comes from the treelist
    # scheduled manual outbreaks (MANSCHED IDT, OPNEW activity 2209) — matched by dfb_outbreak_due.
    mansched_years::Vector{Int32}
    windthr_years::Vector{Int32}   # WINDTHR IDT (activity 2210) — matched by dfb_activity_due (DFBWIN)
    # RANSCHED (IDBSCH=2): the auto-scheduled regional-outbreak cycles (1-based FVS ICYC),
    # computed once before the cycle loop by DFBSCH (dfbsch.f) from the DFBRAN stream seeded at
    # ORSEED (NOT +1128). Empty on the MANSCHED path (mansched_years is used instead).
    scheduled_cycles::Vector{Int}
    scheduled::Bool                # DFBSCH has run for this stand (idempotent guard)
    # DFBRAN Lehmer/MINSTD state (dfbran.f COMMON /DFRCOM/ S0). NaN until the first draw, when it is
    # seeded to ORSEED+1128 — the DFBSCH end-of-routine reset (dfbsch.f:180, TSEED=ORSEED+1128D0;
    # DFBCSD) that FVS applies once before the cycle loop. On the MANSTART/MANSCHED path DFBMOD's
    # BACHLO is the only DFBRAN consumer, so lazy-seeding at first draw is stream-identical.
    rng_s0::Float64
end

"""
    dfb_defaults!(...) → DfbState

FVS `DFBINT` (dfbint.f): the DFB run-time defaults. Called when the first
DFBEETLE keyword is seen.
"""
function dfb_defaults!()
    return DfbState(
        false,          # active (LDFBON)
        false,          # debug
        false,          # lbamod
        Int32(2),       # ismeth  (RANSTART)
        Int32(2),       # idbsch  (RANSCHED)
        Int32(4),       # ilenth
        Int32(10),      # iwait
        Int32(1950),    # ipast
        Int32(0),       # iyout
        6.0f0,          # expctd
        2.0f0,          # exstdv
        0.5f0,          # epiprb
        false,          # lepi
        0.8f0,          # prpwin
        0.0f0,          # minden
        20.0f0,         # mwinht
        0.0f0,          # okill
        55329.0f0,      # orseed
        0.05f0,         # dbevnt
        0.0f0,          # prekll
        false,          # linprg
        false,          # linv
        Int32[],        # mansched_years
        Int32[],        # windthr_years
        Int[],          # scheduled_cycles (RANSCHED)
        false,          # scheduled
        NaN,            # rng_s0 (DFBRAN S0) — lazy-seeded to ORSEED+1128 on first draw
    )
end

# -----------------------------------------------------------------------------
# DFBIND (dfbind.f) — map a tree DBH to a DFB size class 1..20 (deterministic)
# -----------------------------------------------------------------------------
"""
    dfb_ind(dbh) -> Int

FVS `DFBIND`: `IVAL = IFIX(DBH)`; `SZNDX = IVAL/2 + MOD(IVAL,2)`, clamped to
[1,20]. `IFIX` truncates toward zero (matches Fortran on the non-negative DBH
this model ever sees).
"""
@inline function dfb_ind(dbh::Float32)::Int
    ival = trunc(Int, dbh)          # IFIX — truncate toward zero
    szndx = div(ival, 2) + mod(ival, 2)
    szndx > 20 && (szndx = 20)
    szndx < 1  && (szndx = 1)
    return szndx
end

# -----------------------------------------------------------------------------
# DFBDBH (dfbdbh.f) — START[] = DF trees/ac binned by DFB size class
# -----------------------------------------------------------------------------
"""
    dfb_dbh_start(df_dbh, df_tpa) -> Vector{Float32}(20)

FVS `DFBDBH`: for every Douglas-fir record (in the stand's species-sorted IND1
order), add its TPA (PROB) to `START[dfb_ind(dbh)]`. `df_dbh`/`df_tpa` are the
DF records already in IND1 order. Float32 accumulation, in order.
"""
function dfb_dbh_start(df_dbh::AbstractVector{Float32}, df_tpa::AbstractVector{Float32})
    start = zeros(Float32, 20)
    @inbounds for i in eachindex(df_dbh)
        i3 = dfb_ind(df_dbh[i])
        start[i3] += df_tpa[i]
    end
    return start
end

# -----------------------------------------------------------------------------
# DFBER (dfber.f) — stand statistics + minimum-outbreak-condition flag
# -----------------------------------------------------------------------------
"""
    dfb_er(all_dbh, all_tpa, df_dbh, df_tpa, ba) -> (ba9, badf9, a45dbh, pbadf4, lmin)

FVS `DFBER`. `all_*` are every ITRN tree record (for BA9 = stand basal area in
stems ≥9″); `df_*` are the Douglas-fir records in IND1 order; `ba` is the stand
basal area (PLOT `BA`). Returns the stand statistics and `lmin` (NOER: minimum
conditions met = at least 1 TPA of DF ≥4.5″). All Float32, accumulation in order.
"""
function dfb_er(all_dbh::AbstractVector{Float32}, all_tpa::AbstractVector{Float32},
                df_dbh::AbstractVector{Float32}, df_tpa::AbstractVector{Float32},
                ba::Float32)
    ba9    = 0.0f0
    badf45 = 0.0f0
    badf9  = 0.0f0
    a45dbh = 0.0f0
    tdf45  = 0.0f0
    pbadf4 = 0.0f0

    # BA9: basal area of ALL stand trees with DBH ≥ 9.0
    @inbounds for ii in eachindex(all_dbh)
        if all_dbh[ii] >= 9.0f0
            ba9 += DFB_BAF * all_dbh[ii] * all_dbh[ii] * all_tpa[ii]
        end
    end

    lmin = true
    if isempty(df_dbh)
        return (ba9 = ba9, badf9 = badf9, a45dbh = a45dbh, pbadf4 = pbadf4, lmin = false)
    end

    @inbounds for ii in eachindex(df_dbh)
        p = df_tpa[ii]
        d = df_dbh[ii]
        badf = DFB_BAF * d * d * p
        if d >= 4.5f0
            tdf45  += p
            a45dbh += d * p
            badf45 += badf
        end
        if d >= 9.0f0
            badf9 += badf
        end
    end

    if tdf45 < 1.0f0
        # not enough DF for an outbreak
        return (ba9 = ba9, badf9 = badf9, a45dbh = a45dbh, pbadf4 = pbadf4, lmin = false)
    end

    pbadf4 = badf45 / ba
    a45dbh = a45dbh / tdf45
    return (ba9 = ba9, badf9 = badf9, a45dbh = a45dbh, pbadf4 = pbadf4, lmin = lmin)
end

# -----------------------------------------------------------------------------
# DFBPRB (dfbprb.f) — stand-outbreak probability
# -----------------------------------------------------------------------------
"""
    dfb_prb(a45dbh, pbadf4, badf9, ba9) -> Float32

FVS `DFBPRB`: 0 unless A45DBH ≥ 9 and PBADF4 ≥ 0.25; then BADF9/BA9, capped 0.9.
"""
function dfb_prb(a45dbh::Float32, pbadf4::Float32, badf9::Float32, ba9::Float32)::Float32
    protbk = 0.0f0
    a45dbh < 9.0f0  && return protbk
    pbadf4 < 0.25f0 && return protbk
    protbk = badf9 / ba9
    protbk > 0.9f0 && (protbk = 0.9f0)
    return protbk
end

# -----------------------------------------------------------------------------
# DFBRAN (dfbran.f) — the DFB model's own double-precision Lehmer/MINSTD LCG.
# -----------------------------------------------------------------------------
"""
    dfb_rand!(d) -> Float32

FVS `DFBRAN`: `S1 = DMOD(16807·S0, 2147483647)`; `SEL = REAL(S1 / 2147483648)`; `S0 = S1`.
`S0`/`S1` are DOUBLE (COMMON /DFRCOM/), so the modular step is exact integer arithmetic; only
the returned uniform `SEL` is truncated to Float32. Lazy-seeds `S0 = ORSEED + 1128` the first
time (the DFBSCH end reset, dfbsch.f:180). NEVER FFI'd — a faithful reimplementation.
"""
@inline function dfb_rand!(d::DfbState)::Float32
    isnan(d.rng_s0) && (d.rng_s0 = Float64(d.orseed) + 1128.0)   # DFBSCH: S0 = ORSEED + 1128
    s1 = rem(16807.0 * d.rng_s0, 2147483647.0)                   # DMOD (exact, values < 2^31)
    d.rng_s0 = s1
    return Float32(s1 / 2147483648.0)                            # SEL = REAL(S1 / 2^31)
end

# -----------------------------------------------------------------------------
# BACHLO (base/bachlo.f) — normal draw via Batchelor composite rejection, over DFBRAN.
# -----------------------------------------------------------------------------
"""
    dfb_bachlo(xbar, stdev, d) -> Float32

FVS `BACHLO(XBAR,STDEV,DFBRAN)`: a Batchelor composite-rejection normal draw consuming three
DFBRAN uniforms per attempt. Returns `XBAR` when `STDEV ≤ 0`. In the common `U ≤ 2/3` (uniform)
branch the RETURNED variate `X = 1.5·U` carries NO transcendental, so the draw is bit-exact
regardless of the platform `log`; `log` only enters the accept/reject test `Y = −ln(R1) ≤ Z`
and the `U > 2/3` exponential branch. Float32 throughout to mirror the REAL Fortran.
"""
function dfb_bachlo(xbar::Float32, stdev::Float32, d::DfbState)::Float32
    stdev <= 0.0f0 && return xbar
    @inbounds while true
        u  = dfb_rand!(d)
        r1 = dfb_rand!(d)
        r2 = dfb_rand!(d)
        local x::Float32, z::Float32
        if u <= (2.0f0 / 3.0f0)
            x = 1.5f0 * u
            z = 0.5f0 * x * x
        else
            zz = 3.0f0 * u - 2.0f0
            zz < 0.001f0 && continue                # GOTO 10: redraw all three
            x = 1.0f0 - 0.5f0 * log(zz)             # ALOG
            z = 0.5f0 * (x - 2.0f0)^2
        end
        y = -log(r1)                                # −ALOG(R1)
        y <= z && continue                          # reject → GOTO 10
        r2 >= 0.5f0 && (x = -x)
        return x * stdev + xbar
    end
end

# -----------------------------------------------------------------------------
# DFBMOD (dfbmod.f) — projected DF trees/ac killed this outbreak (the NEW-outbreak path).
# -----------------------------------------------------------------------------
"""
    dfb_mod(d, badf9, ba9, numyrs, icyc) -> Float32 (DFKILL)

FVS `DFBMOD`. New-outbreak (`.ELSE.`) branch:
`DFKILL = BACHLO(EXPCTD,EXSTDV,DFBRAN)·(BADF9/BA9)·NUMYRS + OKILL`, redrawn while `< 0`.

CUROUTBK in-progress-outbreak branch (`LINPRG .AND. ICYC == 1`, dfbmod.f:74):
`NUMYRS = NUMYRS − IYOUT` (recomputed but unused in the kill), then
* `IYOUT > 4` → `DFKILL = 0` (outbreak assumed over);
* `PREKLL ≤ 0` (mortality unknown / from the treelist) → redraw-while-<0
  `DFKILL = BACHLO(EXPCTD,EXSTDV)·(BADF9/BA9)·4.0·(1 − PERDD[IYOUT])`;
* `PREKLL > 0` (user-entered kill) → deterministic `DFKILL = PREKLL/PERDD[IYOUT] − PREKLL`.

`NUMYRS` is INTEGER promoted to REAL; Float32 with left-to-right associativity as the Fortran.
"""
function dfb_mod(d::DfbState, badf9::Float32, ba9::Float32, numyrs::Int, icyc::Int)::Float32
    if d.linprg && icyc == 1
        iyout = Int(d.iyout)
        iyout > 4 && return 0.0f0                     # outbreak too long: assumed over
        p = @inbounds DFB_PERDD[iyout]                # IYOUT ∈ 1..4
        if d.prekll <= 0.0f0
            @inbounds while true
                dfkill = dfb_bachlo(d.expctd, d.exstdv, d) * (badf9 / ba9) * 4.0f0 * (1.0f0 - p)
                dfkill >= 0.0f0 && return dfkill
            end
        else
            return d.prekll / p - d.prekll           # deterministic — no DFBRAN draw
        end
    end
    @inbounds while true
        dfkill = dfb_bachlo(d.expctd, d.exstdv, d) * (badf9 / ba9) * Float32(numyrs) + d.okill
        dfkill >= 0.0f0 && return dfkill
    end
end

# -----------------------------------------------------------------------------
# DFBMRT (dfbmrt.f) — distribute DFKILL to WK2 over the DF records with DBH ≥ 9 (DCLAS ≥ 5).
# -----------------------------------------------------------------------------
"""
    dfb_mrt!(t, old_tpa, dfidx, df_dbh, df_tpa, dfkill, badf9, lbamod, okill)

FVS `DFBMRT` (no-windthrow, OKILL≤0 path is the max-combine; OKILL>0 adds). For each Douglas-fir
record (`dfidx` in IND1 order, `df_dbh`/`df_tpa` the cycle-start DBH/PROB) with `dfb_ind(DBH) ≥ 5`
(DBH ≥ 9″): `TAMORT = DFKILL·(D·P/SUMDBH)` (DBH method, default) or `DFKILL·BA/BADF9` (BA method,
MORTDIS), capped at `P`; `WK2 = OKILL>0 ? WK2+TAMORT : max(WK2,TAMORT)`, then capped at `PROB`.
`WK2` is the current periodic mortality (`old_tpa − t.tpa`); the record's TPA is written back as
`old_tpa − WK2`. Records with DBH < 9″ are untouched (byte-identical). SUMDBH accumulates in IND1
order (Float32). Mutates `t.tpa` (and `t.mort_pa`).
"""
function dfb_mrt!(t, old_tpa::AbstractVector{Float32}, dfidx::AbstractVector{<:Integer},
                  df_dbh::AbstractVector{Float32}, df_tpa::AbstractVector{Float32},
                  dfkill::Float32, badf9::Float32, lbamod::Bool, okill::Float32)
    sumdbh = 0.0f0
    if !lbamod
        @inbounds for k in eachindex(dfidx)
            df_dbh[k] >= 9.0f0 && (sumdbh += df_tpa[k] * df_dbh[k])
        end
    end
    @inbounds for k in eachindex(dfidx)
        dd = df_dbh[k]
        dfb_ind(dd) < 5 && continue                 # DCLAS < 5 (DBH < 9″): no DFB kill
        j = dfidx[k]; p = df_tpa[k]
        tamort = lbamod ? dfkill * (DFB_BAF * dd * dd * p) / badf9 :
                          dfkill * (dd * p / sumdbh)
        tamort > p && (tamort = p)
        wk2 = old_tpa[j] - t.tpa[j]                 # current MORTS periodic mortality
        if okill > 0.0f0
            wk2 += tamort                           # windthrow: add DFB to windthrow kill
        elseif tamort > wk2
            wk2 = tamort                            # else take the greater of DFB / background
        end
        wk2 > old_tpa[j] && (wk2 = old_tpa[j])      # WK2 > PROB → PROB
        t.tpa[j] = old_tpa[j] - wk2
        t.tpa[j] < 0.0f0 && (t.tpa[j] = 0.0f0)
    end
    return nothing
end

# -----------------------------------------------------------------------------
# DFB variant gate — the FVS species number of Douglas-fir (dfblkd<v>.f IDFSPC).
# -----------------------------------------------------------------------------
"""
    dfb_idfspc(variant) -> Int

The FVS species index of Douglas-fir for a DFB-linked variant (dfblkd<v>.f `IDFSPC`): 3 for
BC/BM/CI/CR/EC/IE/TT/UT, 16 for PN. `0` (no DF species / DFB not linked) for every other variant,
which keeps the DFB seam inert there.
"""
function dfb_idfspc(variant)::Int
    (variant isa BritishColumbia || variant isa BlueMountains || variant isa CentralIdaho ||
     variant isa CentralRockies || variant isa EastCascades || variant isa InlandEmpire ||
     variant isa Teton || variant isa Utah) && return 3
    variant isa PacificNorthwest && return 16
    return 0
end

# -----------------------------------------------------------------------------
# DFBWIN windthrow tables (dfblkd<v>.f) — shared WINSUC + per-variant IFVSSP crosswalk.
# -----------------------------------------------------------------------------
# WINSUC: the DFB model's 39 species-susceptibility-to-windthrow values (identical in every
# dfblkd<v>.f — a model constant, indexed by the DFB species number via IFVSSP).
const DFB_WINSUC = Float32[
    0.028, 0.083, 0.056, 0.139, 0.111, 0.111, 0.028, 0.139, 0.139, 0.056,
    0.111, 0.042, 0.139, 0.111, 0.056, 0.098, 0.028, 0.056, 0.056, 0.139,
    0.139, 0.0,   0.0,   0.056, 0.139, 0.0,   0.028, 0.056, 0.0,   0.042,
    0.056, 0.0,   0.0,   0.111, 0.056, 0.056, 0.0,   0.0,   0.139,
]
# IFVSSP: FVS-species → DFB-species crosswalk (into WINSUC). Per-variant (dfblkd<v>.f). ROWDOM is
# 80.0 for every species in every variant, MWINHT the DFBINT default 20.0 — carried as scalars.
const _DFB_IFVSSP_IE = Int[1,2,3,4,5,6,7,8,9,10,11,22,23,36,33,26,38,19,24,30,30,18,17]
const _DFB_IFVSSP_PN = Int[16,13,4,9,15,8,39,34,14,8,7,31,12,1,10,3,35,6,5,11,18,18,
                           18,18,32,19,24,29,26,36,22,37,38,18,30,18,30,30,30]
"""
    dfb_ifvssp(variant) -> Vector{Int} | nothing

The FVS-species → DFB-species (WINSUC) crosswalk for a DFB-linked variant (dfblkd<v>.f `IFVSSP`).
Only IE and PN are tabled here (the two windthrow-validated variants); returns `nothing` elsewhere,
which keeps DFBWIN inert for a variant whose table has not been ported.
"""
function dfb_ifvssp(variant)
    variant isa InlandEmpire     && return _DFB_IFVSSP_IE
    variant isa PacificNorthwest && return _DFB_IFVSSP_PN
    return nothing
end

# -----------------------------------------------------------------------------
# DFBWIN (dfbwin.f) — Douglas-fir Beetle windthrow: WK2 windthrow mortality + OKILL feed.
# -----------------------------------------------------------------------------
"""
    dfb_win_kernel!(species, dbh, ht, prob, pct, wk2, ifvssp, rowdom, mwinht, crash, thresh, idfspc, maxsp)
        -> (df9kil, telig)

FVS `DFBWIN` (dfbwin.f) numeric core, on the scheduled-windthrow cycle. A tree record is *eligible*
when `PCT ≥ ROWDOM(=80)` and `HT > MWINHT(=20)`. Accumulating over the FVS species blocks
(`ISPI = 1..MAXSP`, in species-major / ascending-record order = IND1): `SPCNUM` = present species,
`TOTSUC` = Σ WINSUC(IFVSSP(ISPI)) over present species, `ELIGBL(ISPI)` = Σ eligible PROB. If total
eligible `TELIG ≥ THRESH(=MINDEN)` a windthrow fires: per species the eligible proportion killed is
`PRPMRT = CRASH·WINSUC/(TOTSUC/SPCNUM)` capped at 0.95 (the `ELIGBL` factor cancels), each eligible
record's `WK2 += PROB·PRPMRT` (capped at `PROB−1e-6`), and `DF9KIL` sums the Douglas-fir (`IDFSPC`)
`WK2` for `DBH ≥ 9`. Returns `(DF9KIL, TELIG)`; mutates `wk2`. All Float32, accumulation in FVS order.
"""
function dfb_win_kernel!(species::AbstractVector{<:Integer}, dbh::AbstractVector{Float32},
                         ht::AbstractVector{Float32}, prob::AbstractVector{Float32},
                         pct::AbstractVector{Float32}, wk2::AbstractVector{Float32},
                         ifvssp::AbstractVector{<:Integer}, rowdom::Float32, mwinht::Float32,
                         crash::Float32, thresh::Float32, idfspc::Int, maxsp::Int)
    n = length(species)
    eligbl = zeros(Float32, maxsp)
    totsuc = 0.0f0
    spcnum = 0.0f0
    elig(i) = (pct[i] >= rowdom) & (ht[i] > mwinht)
    @inbounds for sp in 1:maxsp                       # DO 300 ISPI=1,MAXSP (species-major, IND1 order)
        present = false
        for i in 1:n
            species[i] == sp || continue
            present = true
            elig(i) && (eligbl[sp] += prob[i])
        end
        if present                                    # ISCT(ISPI,1) > 0
            spcnum += 1.0f0
            totsuc += DFB_WINSUC[ifvssp[sp]]
        end
    end
    telig = 0.0f0
    @inbounds for sp in 1:maxsp                        # DO 400: TELIG = Σ ELIGBL
        telig += eligbl[sp]
    end
    df9kil = 0.0f0
    if telig >= thresh                                 # windthrow occurs
        avg = totsuc / spcnum
        @inbounds for sp in 1:maxsp                     # DO 600 ISPI=1,MAXSP
            prp = 0.0f0
            if eligbl[sp] > 0.0f0
                prp = crash * DFB_WINSUC[ifvssp[sp]] / avg   # ELIGBL cancels in PRPMRT
                prp > 0.95f0 && (prp = 0.95f0)
            end
            for i in 1:n                                # DO 500 J=I1,I2 (IND1 order)
                species[i] == sp || continue
                elig(i) || continue
                w = wk2[i] + prob[i] * prp
                (prob[i] - w < 1.0f-6) && (w = prob[i] - 1.0f-6)
                wk2[i] = w
                (sp == idfspc && dbh[i] >= 9.0f0) && (df9kil += w)   # DF ≥9″ windthrown
            end
        end
    end
    return (df9kil, telig)
end

"""
    dfb_win!(s, old_tpa)

FVS `DFBWIN` engine seam (gradd.f:72, BEFORE `DFBDRV`), wired into the FVSjl mortality path just
ahead of `dfb_apply!`. No-op unless a DFB block is active, a `WINDTHR` event is scheduled this cycle,
and the variant has an IFVSSP table. When it fires it computes the windthrow kill (`dfb_win_kernel!`)
over the current stand — `WK2` seeded from the background mortality (`old_tpa − t.tpa`), `PCT` from
`t.crown_ratio`, `HT`/`DBH`/`PROB` current — writes the raised `WK2` back to `t.tpa`, and (when a
regional outbreak is active this cycle and ≥1 large DF blew down) sets `OKILL = DF9KIL` so the
following `dfb_apply!` runs the windthrow-add DFBMOD/DFBMRT path. Byte-identical when no windthrow
is due. The end-to-end `.sum`-DELTA is CORNERED by the IE growth straddle (the cyc-2 PCT/HT/DBH the
windthrow reads scale with cyc-1 growth); the windthrow math itself is bit-exact on equal inputs.
"""
function dfb_win!(s::StandState, old_tpa::Vector{Float32})
    d = s.dfb
    (d === nothing || !d.active) && return nothing
    dfb_activity_due(d.windthr_years, s) || return nothing    # OPFIND(2210): windthrow scheduled?
    idfspc = dfb_idfspc(s.variant)
    idfspc == 0 && return nothing
    ifvssp = dfb_ifvssp(s.variant)
    ifvssp === nothing && return nothing                      # variant table not ported ⇒ inert
    t = s.trees; n = t.n
    n == 0 && return nothing
    maxsp = length(ifvssp)
    species = Int[Int(t.species[i]) for i in 1:n]
    dbh  = Float32[t.dbh[i] for i in 1:n]
    ht   = Float32[t.height[i] for i in 1:n]
    prob = Float32[old_tpa[i] for i in 1:n]                    # PROB
    pct  = Float32[t.crown_ratio[i] for i in 1:n]             # PCT (stand BA percentile)
    wk2  = Float32[old_tpa[i] - t.tpa[i] for i in 1:n]        # WK2 = background periodic mortality
    crash  = d.prpwin                                          # PRMS(1) = PRPWIN
    thresh = d.minden                                          # PRMS(2) = MINDEN
    df9kil, telig = dfb_win_kernel!(species, dbh, ht, prob, pct, wk2, ifvssp,
                                    80.0f0, d.mwinht, crash, thresh, idfspc, maxsp)
    telig >= thresh || return nothing                          # no windthrow (reschedule; MINDEN>0 edge)
    @inbounds for i in 1:n                                     # apply the raised WK2 to the surviving TPA
        nt = old_tpa[i] - wk2[i]
        nt < 0.0f0 && (nt = 0.0f0)
        t.tpa[i] = nt
    end
    # OKILL: only when a regional outbreak is active this cycle (DFBGO) and ≥1 large DF blew down.
    (df9kil >= 1.0f0 && dfb_outbreak_due(d, s)) && (d.okill = df9kil)
    return nothing
end

"""
    dfb_activity_due(years, s) -> Bool

OPFIND date-match for a DFB activity list (`years` = the scheduled IDT dates): `IDT == 0` every
cycle; `0 < IDT < 1000` a 1-based cycle; `IDT ≥ 1000` a calendar year in this cycle's window. Shared
by the MANSCHED outbreak gate and the WINDTHR event.
"""
function dfb_activity_due(years::AbstractVector{<:Integer}, s::StandState)::Bool
    isempty(years) && return false
    fvscyc = Int(s.control.cycle) + 1
    cyc0 = Int(s.control.cycle)
    cs = cycle_year_at(s.control, cyc0)
    ce = cycle_year_at(s.control, cyc0 + 1); ce <= cs && (ce = cs + 1)
    @inbounds for m in years
        mi = Int(m)
        mi == 0 && return true
        (0 < mi < 1000) && mi == fvscyc && return true
        (mi >= 1000) && (cs <= mi < ce) && return true
    end
    return false
end

# -----------------------------------------------------------------------------
# DFBSCH (dfbsch.f) — RANSCHED auto-scheduler: draw the regional-outbreak cycles.
# -----------------------------------------------------------------------------
"""
    dfb_schedule!(d, c)

FVS `DFBSCH` (dfbsch.f), called once from MAIN (fvs.f:143) BEFORE the cycle loop. On the RANSCHED
path (`IDBSCH == 2`) it walks the DFBRAN stream — seeded at `ORSEED` (the DFBINT/blockdata seed,
NOT the `+1128` cycle-loop seed) — treating each draw as one calendar year and scheduling a regional
outbreak (OPNEW activity 2208) whenever `RAND ≤ DBEVNT`. Between outbreaks it waits `IWAIT` years
plus the `NYR` draws it took to hit the probability, snaps the outbreak to the enclosing FVS cycle
(the cycle year ≤ the drawn year), and continues until the next outbreak would fall past
`IY(NCYC)+10` or the schedule reaches the final cycle. The resulting 1-based cycle numbers are the
regional-outbreak cycles `DFBGO`'s OPFIND then matches. `DFBSCH` ends by resetting the seed to
`ORSEED+1128` (dfbsch.f:180) — modelled here by leaving `d.rng_s0` lazy-seeded to `ORSEED+1128`, so
the cycle-loop draws are stream-identical whether or not scheduling consumed draws. On MANSCHED the
loop is skipped (only the seed reset applies) so this is a no-op there.
"""
function dfb_schedule!(d::DfbState, c)
    d.scheduled = true
    d.idbsch == Int32(2) || return nothing         # only RANSCHED walks the stream (else: seed reset only)
    empty!(d.scheduled_cycles)
    ncyc = Int(c.ncycle); ncyc < 1 && (ncyc = 1)
    iy(k) = Int(cycle_year_at(c, k - 1))            # IY(k), 1-based Fortran subscript
    invyr  = iy(1)
    ifinyr = iy(ncyc) + 10
    iwait  = Int(d.iwait)
    dbevnt = d.dbevnt
    iprev  = Int(d.ipast)
    i1     = 1
    s0     = Float64(d.orseed)                       # DFBSCH start seed = ORSEED (no +1128)
    @inbounds while true                             # label 100
        nyr = 0
        found = false
        while true                                   # label 200
            s1 = rem(16807.0 * s0, 2147483647.0); s0 = s1
            rand = Float32(s1 / 2147483648.0)        # DFBRAN
            if rand <= dbevnt                        # GOTO 400: outbreak year found
                found = true; break
            end
            nyr += 1
            (iprev + iwait + nyr) >= ifinyr && return nothing   # past the run end (GOTO 1000)
            nyr < 100 || return nothing              # DBEVNT too small: warning + exit
        end
        found || return nothing
        next = iprev + iwait + nyr                   # label 400
        next < invyr && (next = invyr)
        iprev = next
        i = i1                                        # find enclosing cycle: first IY(I) > NEXT
        while i <= ncyc && iy(i) <= next
            i += 1
        end
        iyi = i - 1
        push!(d.scheduled_cycles, iyi)               # OPNEW(2208) at IY(IYI) ⇒ regional outbreak in cycle IYI
        i1 = iyi
        iyi >= ncyc && return nothing
    end
end

"""
    dfb_inv!(d, killed_tpa)

FVS `DFBINV` (dfbinv.f), called once from MAIN (fvs.f:174) before the cycle loop. When CUROUTBK
requested the in-progress kill be taken from the tree list (`LINV`, no user PREKLL), it sums the
TPA (PROB) of the Douglas-fir records read in as recently DFB-killed — the trees `DAMCDS`/`DFBDAM`
flagged (damage code 3, severity ≥ 3, tree-history 6/7 ⇒ IMC 7) — into `PREKLL`. `killed_tpa` is
that list of per-record TPA. FVSjl's tree reader does not carry DFB damage/tree-history codes, so on
every currently-loadable stand this list is empty and `PREKLL` stays 0 (the DFBMOD `PREKLL ≤ 0`
BACHLO branch) — bit-exact with the oracle on iet01 (no DFB damage codes). The summation itself is
ported and unit-tested against synthetic records so the arithmetic is faithful when a reader ever
supplies the codes. No-op unless `LINV`.
"""
function dfb_inv!(d::DfbState, killed_tpa::AbstractVector{Float32})
    d.linv || return nothing
    isempty(killed_tpa) && return nothing
    acc = 0.0f0
    @inbounds for p in killed_tpa                    # DFBINV: PREKLL = Σ PROB(IPT(II)), in order
        acc += p
    end
    d.prekll = acc
    return nothing
end

"""
    dfb_setup!(s)

FVS `DFBSCH`/`DFBINV` seam (fvs.f:143/174), wired into `setup_growth!` after the cycle schedule
(`IY`) is built. Inert unless a DFB block is active; runs the RANSCHED auto-scheduler and the
CUROUTBK pre-killed-DF inventory count once per stand. Idempotent.
"""
function dfb_setup!(s::StandState)
    d = s.dfb
    (d === nothing || !d.active || d.scheduled) && return nothing
    dfb_schedule!(d, s.control)
    # DFBINV: FVSjl's treelist carries no DFB damage/tree-history codes, so the killed-DF list is
    # empty and PREKLL stays 0 (matches the oracle on stands without DFB damage codes).
    dfb_inv!(d, Float32[])
    return nothing
end

# -----------------------------------------------------------------------------
# DFBGO outbreak gate (dfbgo.f) — MANSTART/MANSCHED deterministic path.
# -----------------------------------------------------------------------------
"""
    dfb_outbreak_due(d, s) -> Bool

Whether a regional DFB outbreak activity (MANSCHED's code-2209, scheduled by dfbin.f option 9 at
date `IDT`) is due this cycle — the OPFIND(2208/2209) test in `DFBGO`. Mirrors OPNEW/OPCYCL date
matching (`_compute_due`): `IDT == 0` = every cycle; `0 < IDT < 1000` = a 1-based CYCLE number
(fires when the FVS cycle equals it); `IDT ≥ 1000` = a calendar year (fires in the cycle whose
[start,end) window contains it). A single MANSCHED thus fires in exactly one cycle (OPDONE).
"""
function dfb_outbreak_due(d::DfbState, s::StandState)::Bool
    # RANSCHED (IDBSCH=2): the DFBSCH-computed regional-outbreak cycles.
    (Int(s.control.cycle) + 1) in d.scheduled_cycles && return true
    # MANSCHED (IDBSCH=1): the user-scheduled OPNEW(2209) dates.
    return dfb_activity_due(d.mansched_years, s)
end

# -----------------------------------------------------------------------------
# DFBDRV seam (dfbdrv.f: DFBDBH→DFBMOD→DFBMRT), gated by DFBGO — called from the grow/mortality
# path (gradd.f:74) right after MORTS has set WK2, on the non-tripled cycle stand.
# -----------------------------------------------------------------------------
"""
    dfb_apply!(s, old_tpa, fint)

The DFB cycle driver + gate (FVS `DFBGO` → `GRADD` `IF (LDFBGO) CALL DFBDRV`), wired into the
FVSjl mortality path. No-op unless a `DFB` block is active AND a regional outbreak is due this
cycle (MANSTART deterministic gate) AND the stand meets the minimum DF condition (`DFBER` LMIN).
When it fires it computes `DFKILL` (`dfb_mod`) from the cycle-start DF statistics (`dfb_er` on
`old_tpa`/cycle-start DBH) and distributes it to the per-record mortality (`dfb_mrt!`), raising
`t.tpa` reductions on the large Douglas-fir. Byte-identical when no outbreak fires.

MANSTART (ISMETH=1) is the validated deterministic path; RANSTART (ISMETH=2) draws a DFBRAN
inclusion test against the stand-outbreak probability (stochastic) before DFKILL.
"""
function dfb_apply!(s::StandState, old_tpa::Vector{Float32}, fint::Real)
    d = s.dfb
    (d === nothing || !d.active) && return nothing
    dfb_outbreak_due(d, s) || return nothing                 # OPFIND(2208/2209): outbreak scheduled?
    idfspc = dfb_idfspc(s.variant)
    idfspc == 0 && return nothing                            # DFB not linked for this variant
    t = s.trees; n = t.n
    n == 0 && return nothing
    # DF records in ascending record index (IND1 species-block order); cycle-start DBH/PROB.
    dfidx = Int[i for i in 1:n if Int(t.species[i]) == idfspc]
    isempty(dfidx) && return nothing
    all_dbh = Float32[t.dbh[i] for i in 1:n]
    df_dbh  = Float32[t.dbh[i] for i in dfidx]
    df_tpa  = Float32[old_tpa[i] for i in dfidx]
    # PLOT BA (stand basal area) from the cycle-start PROB.
    ba = 0.0f0
    @inbounds for i in 1:n
        ba += DFB_BAF * all_dbh[i] * all_dbh[i] * old_tpa[i]
    end
    r = dfb_er(all_dbh, old_tpa, df_dbh, df_tpa, ba)         # DFBER: BA9/BADF9/A45DBH/PBADF4/LMIN
    r.lmin || return nothing                                 # minimum outbreak conditions not met
    r.ba9 <= 0.0f0 && return nothing                         # guard the BADF9/BA9 division
    # ISMETH=2 (RANSTART) stochastic stand-inclusion test, mirroring DFBGO's DFBRAN draw < PROTBK.
    if d.ismeth == 2
        protbk = d.lepi ? d.epiprb : dfb_prb(r.a45dbh, r.pbadf4, r.badf9, r.ba9)
        dfb_rand!(d) < protbk || return nothing
    end
    # NUMYRS = min(ILENTH, IFINT) capped at 10; IFINT is the cycle length in years (fint).
    numyrs = min(Int(d.ilenth), round(Int, fint))
    numyrs > 10 && (numyrs = 10)
    icyc = Int(s.control.cycle) + 1                          # FVS ICYC (1-based) — CUROUTBK fires only ICYC==1
    dfkill = dfb_mod(d, r.badf9, r.ba9, numyrs, icyc)        # DFBMOD (DFBDBH's START not needed here)
    dfb_mrt!(t, old_tpa, dfidx, df_dbh, df_tpa, dfkill, r.badf9, d.lbamod, d.okill)  # DFBMRT
    d.okill = 0.0f0                                          # DFBMRT clears OKILL
    return nothing
end

# -----------------------------------------------------------------------------
# kw_dfbin! (dfbin.f) — DFB keyword-block reader (initre option 100)
# -----------------------------------------------------------------------------
"""
    kw_dfbin!(s, rec, kr)

Parse the `DFBEETLE … END` block (dfb/dfbin.f). Faithfully sets the DFBCOM-
equivalent state on `s.dfb` and consumes sub-keyword records up to `END`,
exactly as FVS `DFBIN` does. `LDFBON` (active) is set true as soon as the block
is entered. No per-cycle engine seam is wired yet, so this is INERT: a stand
carrying a DfbState still projects byte-identically to one without DFB.

Sub-keywords ported (state-setting): END, CUROUTBK, DEBUG, NODEBUG, MANSTART,
RANSTART, RANNSEED, RANSCHED, MANSCHED, WINDTHR, STOPROB, MORTDIS, OLENGTH,
EXYRMORT. (DFBECHO opens a post-processor report file — recognized, no state.)
The MANSCHED/WINDTHR OPNEW activity-scheduler seam is deferred (years recorded).
"""
function kw_dfbin!(s::StandState, rec, kr::KeywordReader)
    s.dfb === nothing && (s.dfb = dfb_defaults!())
    d = s.dfb
    d.active = true                       # dfbin.f: LDFBON = .TRUE.

    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"                     # option 1
            break
        elseif k == "CUROUTBK"            # option 2
            if r.present[1]
                d.iyout = Int32(trunc(Int, r.values[1]))
                d.linprg = true
                if r.present[2]
                    d.prekll = Float32(r.values[2])
                else
                    d.linv = true
                end
            end
        elseif k == "DEBUG"               # option 3
            d.debug = true
        elseif k == "DFBECHO"             # option 4 — post-processor report file (no state effect)
            # DFBIN reads a supplemental file-name record here; the KeywordReader
            # supplies the block records, so nothing to consume in the parse model.
        elseif k == "MANSTART"            # option 5
            d.ismeth = Int32(1)
        elseif k == "RANSTART"            # option 6
            d.ismeth = Int32(2)
        elseif k == "RANNSEED"            # option 7
            if r.present[1]
                seed = Float32(r.values[1])
                seed % 2.0f0 == 0.0f0 && (seed += 1.0f0)   # DFBNSD: force an ODD seed
                d.orseed = seed
            end
        elseif k == "RANSCHED"            # option 8
            d.idbsch = Int32(2)
            r.present[1] && (d.iwait  = Int32(trunc(Int, r.values[1])))
            r.present[2] && (d.dbevnt = Float32(r.values[2]))
            r.present[3] && (d.ipast  = Int32(trunc(Int, r.values[3])))
        elseif k == "MANSCHED"            # option 9 — schedule regional outbreak at IDT (default 1)
            d.idbsch = Int32(1)
            idt = r.present[1] ? Int32(trunc(Int, r.values[1])) : Int32(1)
            push!(d.mansched_years, idt)
        elseif k == "WINDTHR"             # option 10
            idt = r.present[1] ? Int32(trunc(Int, r.values[1])) : Int32(1)
            r.present[2] && (d.prpwin = Float32(r.values[2]))
            r.present[3] && (d.minden = Float32(r.values[3]))
            push!(d.windthr_years, idt)
        elseif k == "STOPROB"             # option 11
            if 0.0f0 <= r.values[1] <= 1.0f0 && r.present[1]
                d.lepi = true
                d.epiprb = Float32(r.values[1])
            end
        elseif k == "NODEBUG"             # option 12
            d.debug = false
        elseif k == "MORTDIS"             # option 13
            d.lbamod = true
        elseif k == "OLENGTH"             # option 15
            if 1.0f0 <= r.values[1] <= 10.0f0
                d.ilenth = Int32(trunc(Int, r.values[1]))
            end
        elseif k == "EXYRMORT"            # option 16
            r.present[1] && (d.expctd = Float32(r.values[1]))
            r.present[2] && (d.exstdv = Float32(r.values[2]))
        else
            # DFB sub-keyword not recognized — record it so it can't hide as a
            # silent gap (mirrors the RD/keyword-dispatch policy).
            (!isempty(k) && isletter(first(k))) && push!(s.control.unrecognized_keywords, k)
        end
    end
    return nothing
end
