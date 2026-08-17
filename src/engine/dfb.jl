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
#                           (a normal draw). NOT yet ported.
#   * DFBMRT  (dfbmrt.f)  — applies DFKILL to WK2 (mortality). Depends on DFBMOD.
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
# BEACHHEAD SCOPE (this chunk): the keyword reader (gated, inert — no engine seam)
# + the four DETERMINISTIC routines (DFBIND/DFBDBH/DFBER/DFBPRB), validated
# bit-exact (Float32) against an instrumented g16 build of the pristine Fortran
# (test/engine/test_dfb.jl replays scratchpad/dfb/dfb_golden.txt).
#
# DETERMINISTIC vs STOCHASTIC split:
#   DETERMINISTIC (dump-replay bit-exact): DFBIND, DFBDBH(START), DFBER(BA9,
#       BADF9, A45DBH, PBADF4, LMIN), DFBPRB(PROTBK), and the MANSTART/MANSCHED
#       outbreak gate.
#   STOCHASTIC (RNG realization, .sum-DELTA straddle): the RANSTART inclusion
#       draw (DFBGO), the DFKILL magnitude (DFBMOD via BACHLO→DFBRAN), and
#       RANSCHED regional scheduling. These are the NEXT chunks.
#
# Float32 discipline: DFB's Fortran is REAL (Float32). Every arithmetic result
# here is Float32 so the port is bit-identical to the oracle. The 0.005454154
# basal-area-factor constant, the accumulation order, and IFIX truncation are
# reproduced exactly.
# =============================================================================

const DFB_BAF = 0.005454154f0   # basal-area factor: BA = DFB_BAF·DBH²·TPA (dfber.f/dfbmrt.f)

"""
    DfbState

Douglas-fir Beetle model state (FVS `DFBCOM`). Populated by `kw_dfbin!` and
seeded with the `DFBINT` defaults by `dfb_defaults!`. `active` mirrors `LDFBON`
(set true as soon as any DFBEETLE sub-keyword is read). No per-cycle engine seam
is wired yet, so a DfbState-carrying stand still projects byte-identically.
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
    # scheduled manual outbreaks (MANSCHED IDT) — recorded here; the OPNEW/OPFIND
    # activity-scheduler seam (activity code 2209) is a later chunk.
    mansched_years::Vector{Int32}
    windthr_years::Vector{Int32}   # WINDTHR IDT (activity 2210), likewise deferred
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
# kw_dfbin! (dfbin.f) — DFBEETLE keyword-block reader (initre option 100)
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
            d.orseed = Float32(r.values[1])
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
