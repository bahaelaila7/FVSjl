# =============================================================================
# Western Spruce Budworm (WSBWE) defoliation model — wsbwe/*.f
# =============================================================================
# Sixth insect/pathogen event-extension port (after DFB/DFTM/WPBR/WWPB/LPMPB).
# STAND-LEVEL (audit-confirmed: no PPE dependency in the reachable call graph —
# BWEGO is called from base GRINCR/GRADD exactly like DFB's DFBGO, and BWEINT
# from INITRE). Unlike PPE-gated WWPB, a stand oracle IS relinkable and runs
# (FVSem_wsbwe = FVSem .o − base/exbudl.o + wsbwe/*.o; off ≡ stock byte-identical).
#
# BEACHHEAD SCOPE (this file): the RNG (bit-exact) + the `WSBW … END` keyword
# block reader (gated, INERT), mirroring the WWPB/WPBR beachhead. There is NO
# active engine seam yet: the defoliation → growth-loss/mortality effect path
# (BWEGO→BWEDR→BWEDAM/BWEDIE) is DEFERRED to the handoff (it was not yet
# dump-replay validated). A stand with no WSBW block — and, with this reader,
# a stand WITH one — projects byte-identically because no seam is wired.
#
# FVS structure (wsbwe/, 57 files, internal prefix `bwe`):
#   * BWEINT (bweint.f)  — one-time init of BWECOM/BWECM2/BWEBOX defaults
#                          (called unconditionally from INITRE:186). `wsbwe_defaults!`.
#   * BWEIN  (bwein.f)   — the `WSBW … END` block reader (base keywds.f option 8
#                          'WSBW' → initre.f:680 CALL BWEIN). 25 sub-keywords in
#                          the extension's own TABLE (loaded in bwebk<v>.f block
#                          data). Ported as `kw_wsbwe!`.
#   * BWEGO  (bwego.f)   — per-cycle outbreak gate (grincr.f:414 CALL BWEGO;
#                          gate consumed at gradd.f:108 IF (IPMODI==1 .AND.
#                          LBWEGO)). Fires when LDEFOL (manual DEFOL) OR
#                          LCALBW/LBUDL (BUDLITE regional outbreak scheduled this
#                          cycle). DEFERRED (no Julia seam yet).
#   * BWEDR  (bwedr.f)   — the annual within-cycle driver (defoliation dynamics,
#                          BWEDEF/BWEDAM/BWEDIE growth-loss & top-kill & mortality).
#                          DEFERRED — the deep effect path.
#   * BWERAN (bweran.f)  — the model's OWN double-precision MINSTD LCG, seed 55329
#                          set DIRECTLY in every bwebk<v>.f block data (like LPMPB's
#                          MPRANN — NO +1128 offset): S1=DMOD(16807·S0,2147483647);
#                          SEL=REAL(S1/2^31); S0=S1. Entries BWERSD (reseed),
#                          BWERGT/BWERPT (get/put seed). `wsbwe_rand!`/`wsbwe_seed!`.
#                          Consumers (BUDLITE only): BWEBET (beta draws) + BWERNP
#                          (random defoliation given a mean). The MANUAL DEFOL path
#                          draws NO random numbers.
#   * bwebk<v>.f          — per-variant BLOCK DATA (BLOCK DATA BWEBK): host-species
#                          map IBWSPM, TABLE, seed, NEMULT, OBTABL, etc. bwebms<v>.f
#                          (SUBROUTINE BWEBMS): per-variant retained-biomass coeffs.
#                          Host variants with own coeffs: BC BM CI EC EM SO TT
#                          (+ generic bwebk.f/bwebms.f for others). Every FVS*_
#                          buildDir currently links the base/exbudl.f NO-OP stub —
#                          NO shipped binary has WSBWE active (mirrors WWPB).
#
# RNG identity — VALIDATED bit-exact (scratchpad/wsbwe/driver_bweran.f, gfortran-16,
# pristine wsbwe/bweran.f). Seed 55329, first 8 draws (Float32 hex):
#   3EDDB57A 3F5AB21B 3F630A07 3F2734C9 3EF4FB2C 3F4AF646 3F6E7B68 3F67EA59
# IDENTICAL to LPMPB(MPRANN)/DFTM(TMRANN)/WWPB(BMRANN) — same LCG + seed. BWERGT→
# BWERPT save/restore reproduces the stream exactly (validated: draw2 after
# get/put == draw2). Ported in Float64 modular arithmetic, Float32 uniform. NEVER FFI.
#
# ACTIVATION vs OUTPUT-ONLY (bwein.f) — the base keyword is 'WSBW'; the outbreak is
# turned on by a SUB-keyword:
#   ACTIVATION:
#     * DEFOL   (opt 8)  → LDEFOL=.TRUE.  — user-supplied annual defoliation, the
#                          DETERMINISTIC path (NO RNG, NO weather). OPNEW act 2151.
#     * GENDEFOL(opt 14) → LBUDL=.TRUE., IBWCHK=1 — the BUDLITE population model
#                          (STOCHASTIC: draws BWERAN; needs a weather source, else
#                          bwein.f:312 EOF on fort.40). OPNEW act 2150.
#     * SETPRBIO(opt 9)  → OPNEW act 2153 (recompute retained-biomass proportions).
#     * BWSPRAY (opt 24) → OPNEW act 2157/2159 (insecticide).
#     * OBSCHED (opt 25) → outbreak start/end schedule for BUDLITE.
#   OUTPUT-ONLY / CONFIG (no growth effect): OPEN CLOSE MGMTID RANNSEED STARTYR
#     COMMENT DAMAGE NODAMAGE PERDAM NOPERDAM OUTBRLOC RECOVERY WEATHER NEMULT
#     BWOUTPUT PARASITE FQUALDEV FQUALWT TITLE.
#
# HOST species (IBWSPM, per-variant) — budworm defoliation classes:
#   1=WF 2=DF 3=GF 4=AF(subalpine fir) 5=ES(Engelmann spruce) 6=WL(western larch)
#   7=NOT-A-HOST. EM (19 spp) IBWSPM = [7,6,2,7,7,7,7,5,4,7,7,7,7,7,7,7,7,7,7]
#   (WL→6, DF→2, ES→5, AF→4; NHOSTS=6). Other variants' IBWSPM live in their
#   bwebk<v>.f and must be extracted for the seam (see HANDOFF.md).
# =============================================================================

const WSBWE_SEED0 = 55329.0f0   # bweran.f DATA S0/55329D0/ (every bwebk<v>.f)

# -----------------------------------------------------------------------------
# Platform single-precision transcendentals (glibc). gfortran's ALOG/EXP/`**`
# (real exponent) resolve to glibc logf/expf/powf, which differ from Julia's
# native/openlibm Float32 math by ~1 ULP on ~20-30% of inputs. For BIT-EXACT
# parity with the relinked FVS oracle, WSBWE routes every transcendental through
# glibc via ccall. (Measured: 0/60000 mismatch vs gfortran; native Julia log/exp
# mismatched 16k/60k.) This is NOT an RNG FFI — the BWERAN LCG stays pure Julia.
# -----------------------------------------------------------------------------
# wsbwe_apply! LIVE switch. true = the manual-DEFOL effect path is active (EM host). The feeder is
# end-to-end bit-exact vs FVSem_wsbwe (per-year AVPRBO/PEDDS/PEHTG/AVYRMX match the DAM/PDMPE golden)
# after two faithful glue fixes (see wsbwe_apply!): (1) OLDTPA/ORMSQD = DENSE's TPROB/RMSQD
# (grincr.f:281/285), NOT the stale s.plot.old_tpa/old_qmd (0 at cycle 1 ⇒ ALOG(0)/÷0 ⇒ NaN foliage);
# (2) the DEFOL species field is SPDECD-decoded (bwein.f:519) — an alpha "DF"/"ES"/"AF" → species
# index, not read as a raw 0 (which defoliated ALL hosts once per record ⇒ triple defoliation).
# End-to-end host_defol .sum DELTA is bit-exact-or-cornered by the EM #206 growth straddle.
const WSBWE_APPLY_LIVE = true

const WSBWE_LIBM = "libm.so.6"
@inline wsbwe_log(x::Float32) = ccall((:logf, WSBWE_LIBM), Float32, (Float32,), x)
@inline wsbwe_exp(x::Float32) = ccall((:expf, WSBWE_LIBM), Float32, (Float32,), x)
@inline wsbwe_pow(x::Float32, y::Float32) = ccall((:powf, WSBWE_LIBM), Float32, (Float32,Float32), x, y)
# gfortran real4 ** int4 (libgfortran pow_r4_i4): right-to-left binary exponentiation.
# Julia's `^`(Int) uses a different multiply order ⇒ 1-ULP drift; replicate exactly.
@inline function wsbwe_powi(base::Float32, n::Integer)::Float32
    x = base; e = Int(n); pw = 1.0f0
    while true
        (e & 1) != 0 && (pw *= x)
        e >>= 1
        e != 0 ? (x *= x) : break
    end
    return pw
end
@inline wsbwe_ifix(x::Float32) = trunc(Int, x)   # Fortran IFIX (truncate toward zero)

"""
    WsbweState

Western Spruce Budworm model state (FVS `BWECOM`/`BWECM2`/`BWEBOX`). `active`
mirrors "a WSBW sub-keyword was read". No engine seam is wired yet, so `active`
alone changes nothing in the projection — the reader is faithfully inert. The
activation flags/schedules are captured so the deferred BWEGO→BWEDR seam can be
wired without re-reading keywords.
"""
mutable struct WsbweState <: AbstractWsbweState
    active::Bool                 # a WSBW block was entered
    ldefol::Bool                 # LDEFOL — user DEFOL manual defoliation (deterministic path)
    lbudl::Bool                  # LBUDL — GENDEFOL BUDLITE outbreak requested (stochastic)
    lbwdam::Bool                 # LBWDAM — DAMAGE output flag
    lbwpdm::Bool                 # LBWPDM — PERDAM periodic-damage output flag
    iobloc::Int32                # IOBLOC — outbreak geographic location (1=SW,2=NW,3=MT); default 2
    iwopt::Int32                 # IWOPT — weather option (1=incl. variation, 2=means only)
    iobopt::Int32                # IOBOPT — OBSCHED option (1..4)
    ipara::Int32                 # IPARA — PARASITE option
    dseed::Float32               # damage-model RNG seed (RANNSEED); default 55329
    obseed::Float32              # OUTBRLOC outbreak-timing RNG seed
    wseed::Float32               # WEATHER RNG seed
    mgmidb::String               # MGMTID — budworm management id tag
    # scheduled DEFOL activities (OPNEW act 2151): (idt, bwspecies, crown, new,1yr,2yr,remaining)
    defol_sched::Vector{NTuple{7,Float32}}
    # OBSCHED user outbreak windows (start,end) years (IOBOPT==3)
    obsched::Vector{Tuple{Int32,Int32}}
    rng_s0::Float64              # BWERAN COMMON S0 — NaN until first draw, then seeded to `dseed`.
end

"""
    wsbwe_defaults!() → WsbweState

FVS `BWEINT` + bwebk<v>.f block-data run-time defaults, set when the first WSBW
sub-keyword is seen.
"""
function wsbwe_defaults!()
    return WsbweState(
        false,          # active
        false,          # ldefol
        false,          # lbudl
        false,          # lbwdam
        false,          # lbwpdm
        Int32(2),       # iobloc (NW default)
        Int32(1),       # iwopt
        Int32(1),       # iobopt
        Int32(1),       # ipara
        WSBWE_SEED0,    # dseed (55329)
        WSBWE_SEED0,    # obseed
        WSBWE_SEED0,    # wseed
        "",             # mgmidb
        NTuple{7,Float32}[],
        Tuple{Int32,Int32}[],
        NaN,            # rng_s0
    )
end

# -----------------------------------------------------------------------------
# BWERAN (bweran.f) — the WSBWE model's own double-precision MINSTD LCG. NEVER FFI.
# -----------------------------------------------------------------------------
"""
    wsbwe_rand!(w) -> Float32

FVS `BWERAN`: `S1=DMOD(16807·S0,2147483647)`; `SEL=REAL(S1/2147483648)`; `S0=S1`.
Double-precision modular step (exact); only the returned uniform is Float32.
Lazy-seeds `S0=dseed` (default 55329) the first time. VALIDATED bit-exact (8-draw
golden, seed 55329).
"""
@inline function wsbwe_rand!(w::WsbweState)::Float32
    isnan(w.rng_s0) && (w.rng_s0 = Float64(w.dseed))
    s1 = rem(16807.0 * w.rng_s0, 2147483647.0)
    w.rng_s0 = s1
    return Float32(s1 / 2147483648.0)
end

"""
    wsbwe_seed!(w, seed; lset=true)

FVS `BWERSD(LSET,SEED)`. `lset=true`: force the seed odd (`AMOD(SEED,2)==0 ⇒
SEED+1`) then store it as the stream start and current state. `lset=false`:
restart the stream from the stored start seed (`w.dseed`). Mirrors the RANNSEED
keyword's `BWERSD(LNOTBK(1),DSEEDR)` call.
"""
function wsbwe_seed!(w::WsbweState, seed::Real; lset::Bool=true)
    if lset
        s = Float32(seed)
        (s % 2.0f0 == 0.0f0) && (s += 1.0f0)
        w.dseed = s
        w.rng_s0 = Float64(s)
    else
        w.rng_s0 = Float64(w.dseed)
    end
    return nothing
end

# -----------------------------------------------------------------------------
# kw_wsbwe! (bwein.f) — the `WSBW … END` block reader (keywds.f option 8 'WSBW').
# -----------------------------------------------------------------------------
"""
    kw_wsbwe!(s, rec, kr)

Parse the `WSBW … END` block (wsbwe/bwein.f). Sets the BWECOM/BWECM2-equivalent
state on `s.wsbwe` and consumes sub-keyword records up to `END`. `active` is set
true as soon as the block is entered. INERT: no engine seam is wired, so the
projection is byte-identical whether or not this block is present (validated A/B).

Sub-keywords (extension TABLE, bwebk<v>.f): END OPEN CLOSE MGMTID RANNSEED STARTYR
COMMENT DEFOL SETPRBIO DAMAGE NODAMAGE PERDAM NOPERDAM GENDEFOL OUTBRLOC RECOVERY
WEATHER NEMULT BWOUTPUT PARASITE FQUALDEV FQUALWT TITLE BWSPRAY OBSCHED. The reader
captures the activation-relevant state (LDEFOL/LBUDL/schedules/seeds); output-file
opening and report formatting are intentionally omitted (no .sum effect). Several
keywords read a trailing data record in Fortran (COMMENT to 'END', WEATHER file
name, TITLE) — see HANDOFF.md before wiring the effect seam.
"""
function kw_wsbwe!(s::StandState, rec, kr::KeywordReader)
    s.wsbwe === nothing && (s.wsbwe = wsbwe_defaults!())
    w = s.wsbwe
    w.active = true
    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            break
        elseif k == "MGMTID"                 # opt 4 (reads a trailing A4 record in Fortran)
            # management id tag; no .sum effect.
        elseif k == "RANNSEED"               # opt 5 — damage-model seed
            if r.present[1] && r.values[1] != 0.0f0 && r.values[1] > 10000.0f0
                wsbwe_seed!(w, Float32(r.values[1]); lset=true)
            end
        elseif k == "STARTYR"                # opt 6
            # stand-alone first projection year; no linked-mode .sum effect.
        elseif k == "DEFOL"                  # opt 8 — MANUAL defoliation (deterministic activation)
            w.ldefol = true
            idt  = r.present[1] ? Float32(trunc(r.values[1])) : 1.0f0
            # SPDECD (bwein.f:519 `CALL SPDECD(2,IS,...)`): the species field is an ALPHA code
            # ("DF"/"ES"/"AF") or a numeric species index — decode it to the FVS species SEQUENCE
            # index IS, NOT read as a raw number (a raw `r.values[2]` on alpha "DF" is 0, which
            # `spc<=0` treats as ALL hosts, so each of DF/ES/AF defoliates every host → the host
            # gets defoliated 3× ⇒ 0.15³ retained-biomass instead of 0.15). IS==0/-999 ⇒ 0 (ALL),
            # matching `IF (IS.EQ.-999.OR.IS.EQ.0) ARRAY(2)=0.0`. A nonhost IS (IBWSPM==7) is left
            # for the feeder's `is1>6` skip (bwein.f GOTO 10), same defoliation outcome.
            spp  = Float32(species_selector(s, length(r.fields) >= 2 ? r.fields[2] : ""))
            crn  = (r.present[3] && r.values[3] <= 15.0f0) ? Float32(r.values[3]) : 0.0f0
            new  = r.present[4] ? Float32(r.values[4]) : 0.0f0
            y1   = r.present[5] ? Float32(r.values[5]) : 0.0f0
            y2   = r.present[6] ? Float32(r.values[6]) : 0.0f0
            rem_ = r.present[7] ? Float32(r.values[7]) : 0.0f0
            push!(w.defol_sched, (idt, spp, crn, new, y1, y2, rem_))
        elseif k == "SETPRBIO"               # opt 9 — schedule biomass-proportion recompute
            # OPNEW act 2153; effect path deferred.
        elseif k == "DAMAGE"                 # opt 10
            w.lbwdam = true
        elseif k == "NODAMAGE"               # opt 11
            w.lbwdam = false
        elseif k == "PERDAM"                 # opt 12
            w.lbwpdm = true
        elseif k == "NOPERDAM"               # opt 13
            w.lbwpdm = false
        elseif k == "GENDEFOL"               # opt 14 — BUDLITE outbreak (stochastic, needs weather)
            w.lbudl = true
        elseif k == "OUTBRLOC"               # opt 15
            (r.present[1]) && (w.iobloc = Int32(trunc(r.values[1])))
            if r.present[2]
                wsbwe_seed!(w, Float32(r.values[2]); lset=true); w.obseed = w.dseed
            end
        elseif k == "RECOVERY"               # opt 16 — under development, no effect
        elseif k == "WEATHER"                # opt 17 (reads trailing file-name record[s])
            (r.present[1]) && (w.iwopt = Int32(trunc(r.values[1])))
            if r.present[3]
                wsbwe_seed!(w, Float32(r.values[3]); lset=true); w.wseed = w.dseed
            end
        elseif k == "NEMULT"                 # opt 18 — natural-enemy multipliers (BUDLITE)
        elseif k == "BWOUTPUT"               # opt 19 — output-table flags
        elseif k == "PARASITE"               # opt 20
            (r.present[1]) && (w.ipara = Int32(trunc(r.values[1])))
        elseif k == "FQUALDEV"               # opt 21
        elseif k == "FQUALWT"                # opt 22
        elseif k == "TITLE"                  # opt 23 (reads a trailing A72 record)
        elseif k == "BWSPRAY"                # opt 24 — insecticide (OPNEW 2157/2159)
        elseif k == "OBSCHED"                # opt 25 — outbreak schedule
            (r.present[1]) && (w.iobopt = Int32(trunc(r.values[1])))
            if w.iobopt == 3
                # up to 3 (start,end) pairs in fields (2,3),(4,5),(6,7)
                kk = 0
                for _ in 1:3
                    kk += 2
                    st = r.present[kk]   ? Int32(trunc(r.values[kk]))   : Int32(0)
                    en = r.present[kk+1] ? Int32(trunc(r.values[kk+1])) : Int32(0)
                    (st != 0 && en != 0) && push!(w.obsched, (st, en))
                end
            end
        elseif k == "OPEN" || k == "CLOSE" || k == "COMMENT"
            # file I/O / comment block; no .sum effect (COMMENT reads to its own 'END').
        else
            (!isempty(k) && isletter(first(k))) && push!(s.control.unrecognized_keywords, k)
        end
    end
    return nothing
end

# =============================================================================
# DEFOLIATION EFFECT KERNEL (deterministic payload math) — DUMP-REPLAY VALIDATED
# bit-exact vs the instrumented live oracle FVSem_wsbwe (single-.o swap of bwedam/
# bwepdm, Float32-hex TRANSFER dump; instrumented .sum byte-identical to the clean
# relink). Host stand = DF/ES/AF-relabelled emt01, sustained 85/80/70/60% DEFOL
# 1990-1999 (crashes TPA 536→52 by 2000). See scratchpad/wsbwe/GOLDENS.md.
#
# These two kernels are the arithmetic core of the DEFOL path's growth-loss and
# mortality. They are pure Float32 and reproduce gfortran's `**`/EXP bit-exactly.
# The SURROUNDING chain that feeds them — BWESIT foliage biomass (BWEBMS/BWEADV),
# BWEAGE aging + PRBIO, BWEDAM AVPRBO/CUMDEF accumulation, and BWEPDM's per-tree
# application (which draws the damage RNG via BWERNP/BWEBET + topkill) — is the
# NEXT chunk (see HANDOFF.md "NEXT CHUNKS"); until it is ported the seam below is
# gated INERT.
# -----------------------------------------------------------------------------

# BWEDAM host-indexed small-tree height-growth coeff (STHTGR, bwedam.f:52); classes
# 1=WF 2=DF 3=GF 4=AF 5=ES 6=WL. (Used by the small-tree Ferguson HTG branch.)
const WSBWE_STHTGR = Float32[-2.3661, -2.4757, -2.0008, -2.3661, -2.9171, 0.0]

"""
    wsbwe_rdds(avprbo_whole, rddsm1) -> Float32

FVS `BWEDAM` proportional diameter-growth (DDS) multiplier — Nichols (1984/88)
model (bwedam.f:184). `avprbo_whole` = average whole-tree proportion of retained
biomass; `rddsm1` = last-period RDDS carry (1.0 at outbreak start). Capped at 1.0.
Pure Float32. VALIDATED 90/90 bit-exact vs FVSem_wsbwe (bwedam dump-replay).
"""
@inline function wsbwe_rdds(avprbo_whole::Float32, rddsm1::Float32)::Float32
    r = 0.083861f0 * wsbwe_pow(avprbo_whole * 100.0f0, 0.4725f0 + 0.07f0 * rddsm1) *
        wsbwe_pow(rddsm1, 0.3241f0)
    return r > 1.0f0 ? 1.0f0 : r
end

"""
    wsbwe_rhtg(avprbo_top, rhtgm1) -> Float32

FVS `BWEDAM` proportional height-growth multiplier for MEDIUM/LARGE trees — Nichols
model (bwedam.f:188). `avprbo_top` = average top-third proportion of retained
biomass; `rhtgm1` = last-period RHTG carry (1.0 at outbreak start). Capped at 1.0.
Pure Float32. VALIDATED 90/90 bit-exact vs FVSem_wsbwe. (Small trees use the
Ferguson `exp(STHTGR*(1-AVPRBO_top))` branch instead — see WSBWE_STHTGR.)
"""
@inline function wsbwe_rhtg(avprbo_top::Float32, rhtgm1::Float32)::Float32
    r = 0.193013f0 * wsbwe_pow(avprbo_top * 100.0f0, 0.3814f0 - 0.0212f0 * rhtgm1) *
        wsbwe_pow(rhtgm1, 0.5509f0)
    return r > 1.0f0 ? 1.0f0 : r
end

# BWEPDM Marsden logistic MORTALITY coefficients (bwepdm.f B0..B7), host-indexed.
const WSBWE_B0 = Float32[46.27900, 57.75010, 46.27900, 46.27900, 57.75010, 57.75010]  # intercept
const WSBWE_B1 = Float32[-1.80810, -2.28210, -1.80810, -1.80810, -2.28210, -2.28210]  # elev
const WSBWE_B2 = Float32[0.01930, 0.02430, 0.01930, 0.01930, 0.02430, 0.02430]        # elev^2
const WSBWE_B3 = Float32[0.00550, 0.0, 0.00550, 0.00550, 0.0, 0.0]                     # point BA
const WSBWE_B4 = Float32[-0.00808, -0.00793, -0.00808, -0.00808, -0.00793, -0.00793]  # host BA
const WSBWE_B5 = Float32[0.57450, 0.92870, 0.57450, 0.57450, 0.92870, 0.92870]        # missing-fol top
const WSBWE_B6 = Float32[-0.24050, -0.22330, -0.24050, -0.24050, -0.22330, -0.22330]  # topkill cat
const WSBWE_B7 = Float32[-0.09640, -0.13180, -0.09640, -0.09640, -0.13180, -0.13180]  # MFT*MFM

"""
    wsbwe_mort_pr(ih, elev, pntba, pnthba, mft, mfm, ktk) -> Float32

FVS `BWEPDM` probability-of-mortality logistic (bwepdm.f:640) — Marsden analysis of
Hostetler's data. `ih`=host class 1..6; `elev`=stand elevation (hundred-ft, FVS
`ELEV`); `pntba`/`pnthba`=point total / host basal area; `mft`/`mfm`=missing-foliage
top/middle category (0..10, from PRBIO); `ktk`=topkill category (0..10). Returns 0
when either missing-foliage category is 0 (bwepdm.f:645 guard). Pure Float32.
VALIDATED 81/81 bit-exact vs FVSem_wsbwe (bwepdm dump-replay). This is the SURVIVAL/
mortality probability BEFORE the period-scaling `PR**IBWYR·(1-BASE)**FA` and the
background-rate `max(BASE,PR)` composition (both deterministic; ported in the next
chunk with the per-tree WK2 application).
"""
@inline function wsbwe_mort_pr(ih::Integer, elev::Float32, pntba::Float32,
                               pnthba::Float32, mft::Float32, mfm::Float32,
                               ktk::Float32)::Float32
    (mft > 0.0f0 && mfm > 0.0f0) || return 0.0f0
    e = WSBWE_B0[ih] + (WSBWE_B1[ih] * elev) + (WSBWE_B2[ih] * elev * elev) +
        (WSBWE_B3[ih] * pntba) + (WSBWE_B4[ih] * pnthba) + (WSBWE_B5[ih] * mft) +
        (WSBWE_B6[ih] * ktk) + (WSBWE_B7[ih] * mft * mfm)
    return 1.0f0 / (1.0f0 + wsbwe_exp(e))
end

# =============================================================================
# DEFOLIATION EFFECT PATH — the live BWECUP chain (feeder + per-tree apply).
# Dump-replay VALIDATED bit-exact vs instrumented FVSem_wsbwe on host_defol.key:
#   * BWEBMS  WK4 per-tree foliage biomass ...... 81/81 bit-exact
#   * BWESIT  FOLPOT/FOLADJ/BWTPHA/FNEW/FOLD1/FOLD2/FREM/HOSTST .... bit-exact
#   * BWEAGE/BWEDAM per-year PRBIO/TOTR/AVPRBO/PEDDS/PEHTG/AVYRMX/CDEF . bit-exact
#   * BWEPDM  per-tree topkill+DG+HTG+WK2 + BWERNP/BWEBET RNG ... 81/81 trees,
#             519/519 damage draws bit-exact (draw ORDER exact).
# EM host coefficients (bwebkem.f/bwebmsem.f). Other host variants (BC BM CI EC SO
# TT + generic) need their own IBWSPM/IBIOMP/BWEBMS coeffs — see HANDOFF.md.
# -----------------------------------------------------------------------------

# --- EM block data (bwebkem.f) ---
const WSBWE_IBWSPM_EM = Int[7,6,2,7,7,7,7,5,4,7,7,7,7,7,7,7,7,7,7]     # FVS EM sp → host class (7=nonhost)
const WSBWE_PRCRN3 = Float32[0.05,0.3,0.65,0.15,0.45,0.40,0.15,0.45,0.40]
# THEOFL(4,9), RELFX(4,9), RELFY(4,9): Fortran column-major (age fastest) ⇒ [age,crown]
const WSBWE_THEOFL = reshape(Float32[
 0.70,0.24,0.04,0.02, 0.45,0.30,0.20,0.05, 0.35,0.30,0.25,0.10,
 0.50,0.35,0.12,0.03, 0.30,0.25,0.20,0.25, 0.10,0.20,0.20,0.50,
 0.50,0.30,0.15,0.05, 0.30,0.25,0.20,0.25, 0.10,0.20,0.20,0.50], (4,9))
const WSBWE_RELFX = reshape(Float32[
 0.0,0.04,0.95,1.0, 0.0,0.04,0.95,1.0, 0.0,0.04,0.95,1.0,
 0.0,0.04,0.70,1.0, 0.0,0.12,0.80,1.0, 0.0,0.15,0.90,1.0,
 0.0,0.04,0.70,1.0, 0.0,0.12,0.80,1.0, 0.0,0.15,0.90,1.0], (4,9))
const WSBWE_RELFY = reshape(Float32[
 0.05,0.05,1.0,1.0, 0.05,0.05,1.0,1.0, 0.05,0.05,1.0,1.0,
 0.05,0.05,1.0,1.0, 0.15,0.15,1.0,1.0, 0.20,0.20,1.0,1.0,
 0.05,0.05,1.0,1.0, 0.15,0.15,1.0,1.0, 0.20,0.20,1.0,1.0], (4,9))
# BWEBMS EM (bwebmsem.f) — ICVOPT=2 (DDS model)
const WSBWE_IBIOMP_EM = Int[1,2,3,4,2,11,7,8,9,10,11,11,11,11,11,11,11,11,11]

# --- TT (Teton) block data (bwebktt.f / bwebmstt.f). TT MAXSP=18. Hosts (IBWSPM<7):
# sp3=DF(class2), sp8=ES(class5), sp9=AF(class4) — same host classes as EM. The host-
# class-indexed tables (PRCRN3/THEOFL/RELFX/RELFY) and the ICVOPT=2 biomass coeff arrays
# (BINT2/BINT12/BCL12/…) are IDENTICAL to EM (verified vs bwebmstt.f); only IBWSPM (species
# → host class) and IBIOMP (species → biomass eqn) differ per variant.
const WSBWE_IBWSPM_TT = Int[7,7,2,7,7,7,7,5,4,7,7,7,7,7,7,7,7,7]           # bwebktt.f DATA IBWSPM
const WSBWE_IBIOMP_TT = Int[1,4,3,11,8,11,7,8,9,10,11,11,11,11,11,11,11,11] # bwebmstt.f DATA IBIOMP

# --- BM (Blue Mountains) block data (bwebkbm.f / bwebmsbm.f). MAXSP=18. Hosts: sp3=DF(2),
# sp4=GF(3), sp8=ES(5), sp9=AF(4). Host-class tables (PRCRN3/THEOFL/RELFX/RELFY) AND the
# ICVOPT=2 biomass coeff arrays are IDENTICAL to EM (verified) — only IBWSPM/IBIOMP differ.
const WSBWE_IBWSPM_BM = Int[7,7,2,3,7,7,7,5,4,7,7,7,7,7,7,7,7,7]           # bwebkbm.f DATA IBWSPM
const WSBWE_IBIOMP_BM = Int[1,2,3,4,5,11,7,8,9,10,11,11,11,6,11,11,11,11]  # bwebmsbm.f DATA IBIOMP

# --- SO (SouthCentralOregon) block data (bwebkso.f / bwebmsso.f). MAXSP=33. Hosts:
# sp3=DF(2), sp4=WF(1), sp8=ES(5), sp12=GF(3), sp13=AF(4), sp17=WL(6) — adds WF (class 1)
# and WL (class 6). Host-class tables + biomass coeffs identical to EM (verified).
const WSBWE_IBWSPM_SO = Int[7,7,2,1,7,7,7,5,7,7,7,3,4,7,7,7,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]  # bwebkso.f
const WSBWE_IBIOMP_SO = Int[1,1,3,4,5,6,7,8,3,10,11,4,9,9,3,1,2,6,5,11,11,11,11,11,11,11,11,11,11,11,11,11,11] # bwebmsso.f

# CI (Central Idaho, MAXSP=19): sp3=DF(2), sp4=WF/GF(3), sp8=ES(5), sp9=AF(4). bwebkci.f/bwebmsci.f
# — the biomass coeffs (BINT11/BCL11/BINT12/BCL12/BINT2) + host-defol (PRCRN3) are byte-identical to
# EM (verified); only IBWSPM/IBIOMP differ ⇒ a clean host-class swap.
const WSBWE_IBWSPM_CI = Int[7,7,2,3,7,7,7,5,4,7,7,7,7,7,7,7,7,7,7]           # bwebkci.f DATA IBWSPM
const WSBWE_IBIOMP_CI = Int[1,2,3,4,5,6,7,8,9,10,1,11,11,11,11,4,11,11,11]   # bwebmsci.f DATA IBIOMP

# EC (East Cascades, MAXSP=32): sp3=DF(2), sp6=GF(3), sp8=ES(5), sp9=AF(4). biomass +
# host-defol coeffs byte-identical to EM (verified) ⇒ clean host-class swap.
const WSBWE_IBWSPM_EC = Int[7,7,2,7,7,3,7,5,4,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]  # bwebkec.f
const WSBWE_IBIOMP_EC = Int[1,2,3,4,6,4,7,8,9,10,5,5,11,1,3,4,2,6,11,11,11,11,11,11,11,11,11,11,11,11,11,11] # bwebmsec.f

const WSBWE_IBWSPM_BC = Int[7,6,2,3,7,7,7,5,4,7,7,7,7,2,7]
const WSBWE_IBIOMP_BC = Int[1,2,3,4,5,6,7,8,9,10,11,11,11,3,11]

# Per-variant host/biomass dispatch (mirrors the insect-model `mpb_idxlp`/`dfb_idfspc`
# per-variant dispatch). A variant is WSBWE-host-supported iff it returns non-nothing.
wsbwe_ibwspm_for(v) = v isa EasternMontana    ? WSBWE_IBWSPM_EM :
                      v isa Teton             ? WSBWE_IBWSPM_TT :
                      v isa BlueMountains      ? WSBWE_IBWSPM_BM :
                      v isa CentralIdaho       ? WSBWE_IBWSPM_CI :
                      v isa EastCascades       ? WSBWE_IBWSPM_EC :
                      v isa BritishColumbia    ? WSBWE_IBWSPM_BC :
                      v isa SouthCentralOregon ? WSBWE_IBWSPM_SO : nothing
wsbwe_ibiomp_for(v) = v isa EasternMontana    ? WSBWE_IBIOMP_EM :
                      v isa Teton             ? WSBWE_IBIOMP_TT :
                      v isa BlueMountains      ? WSBWE_IBIOMP_BM :
                      v isa CentralIdaho       ? WSBWE_IBIOMP_CI :
                      v isa EastCascades       ? WSBWE_IBIOMP_EC :
                      v isa BritishColumbia    ? WSBWE_IBIOMP_BC :
                      v isa SouthCentralOregon ? WSBWE_IBIOMP_SO : nothing
const WSBWE_BINT2  = Float32[2.666072,1.756537,2.705866,3.115084,2.654572,3.059351,2.622505,3.300852,3.060169,2.452492,2.622505]
const WSBWE_BINT12 = Float32[-1.94951,-4.73762,-2.05828,-2.43200,-4.17456,-2.24876,-3.13488,-2.93508,-1.60998,-2.74410,-2.63387]
const WSBWE_BCL12  = Float32[1.22023,1.98479,1.25837,1.60270,2.00749,1.37600,1.62368,1.96125,1.32649,1.58171,1.35092]
# BWEAGE recovery curves + BWEPDM topkill coeffs
const WSBWE_RECRX = Float32[0.,.2,.8,1.]
const WSBWE_RECRY = Float32[.6,.6,0.,0.]
const WSBWE_RECHST = Float32[0.5,1.0,0.5,0.5,1.0,0.5]
const WSBWE_IFIR = Int[0,1,0,0,1,1]                                    # bwepdm.f IFIR (fir flag)
const WSBWE_TK0=Float32[21.89880,17.76920,21.89880,21.89880,17.76920,17.76920]
const WSBWE_TK1=Float32[-0.00714,-0.00535,-0.00714,-0.00714,-0.00535,-0.00535]
const WSBWE_TK2=Float32[6.99f-6,5.22f-6,6.99f-6,6.99f-6,5.22f-6,5.22f-6]
const WSBWE_TK3=Float32[0.0,-0.01270,0.0,0.0,-0.01270,-0.01270]
const WSBWE_TK4=Float32[-5.26520,-7.94630,-5.26520,-5.26520,-7.94630,-7.94630]
const WSBWE_TK5=Float32[-0.13630,-0.20360,-0.13630,-0.13630,-0.20360,-0.20360]
const WSBWE_TK6=Float32[0.28030,0.12690,0.28030,0.28030,0.12690,0.12690]
const WSBWE_TK7=Float32[0.14030,0.38010,0.14030,0.14030,0.38010,0.38010]
const WSBWE_TK8=Float32[-0.03770,-0.04120,-0.03770,-0.03770,-0.04120,-0.04120]

# --- BWEBET (Cheng 1978 rejection sampler): 2 BWERAN draws per pass. Returns (x,kode). ---
function wsbwe_bebet(w::WsbweState, a0::Float32, b0::Float32)
    (a0 <= 0.00001f0 || b0 <= 0.00001f0) && return (0.0f0, -1)
    kode = 0
    if min(a0,b0) > 1.0f0
        a = min(a0,b0); b = max(a0,b0)
        alf = a + b; bet = sqrt((alf-2.0f0)/(2.0f0*a*b-alf)); gam = a + (1.0f0/b)
        local w_
        while true
            kode += 1
            u1 = wsbwe_rand!(w); u2 = wsbwe_rand!(w)
            v = bet*wsbwe_log(u1/(1.0f0-u1)); v > 20.0f0 && (v = 20.0f0)
            w_ = a*wsbwe_exp(v); z = u1*u1*u2; rr = gam*v - 1.3862944f0; s = a + rr - w_
            s + 2.609438f0 >= 5.0f0*z && break
            t = wsbwe_log(z); s >= t && break
            rr + alf*wsbwe_log(alf/(b+w_)) < t && continue
            break
        end
        return ((a == a0) ? w_/(b+w_) : b/(b+w_), kode)
    else
        a = max(a0,b0); b = min(a0,b0)
        alf = a + b; bet = 1.0f0/b; del = 1.0f0 + a - b
        c1 = del*(0.0138889f0 + 0.0416667f0*b)/(a*bet - 0.777778f0)
        c2 = 0.25f0 + (0.5f0 + 0.25f0/del)*b
        local w_
        while true
            kode += 1
            u1 = wsbwe_rand!(w); u2 = wsbwe_rand!(w)
            if u1 >= 0.5f0
                z = u1*u1*u2
                if z <= 0.25f0
                    v = bet*wsbwe_log(u1/(1.0f0-u1)); v > 20.0f0 && (v = 20.0f0)
                    w_ = a*wsbwe_exp(v); break
                end
                z >= c2 && continue
                v = bet*wsbwe_log(u1/(1.0f0-u1)); v > 20.0f0 && (v = 20.0f0)
                w_ = a*wsbwe_exp(v)
                alf*(wsbwe_log(alf/(b+w_))+v) - 1.3862944f0 < wsbwe_log(z) && continue
                break
            else
                y = u1*u2; z = u1*y
                if 0.25f0*u2 + z - y >= c1
                    continue
                else
                    v = bet*wsbwe_log(u1/(1.0f0-u1)); v > 20.0f0 && (v = 20.0f0)
                    w_ = a*wsbwe_exp(v)
                    alf*(wsbwe_log(alf/(b+w_))+v) - 1.3862944f0 < wsbwe_log(z) && continue
                    break
                end
            end
        end
        return ((a == a0) ? w_/(b+w_) : b/(b+w_), kode)
    end
end

# --- BWERNP (bwernp.f): beta variate whose variance is a function of the mean. ---
const WSBWE_BC = 2.7f0
const WSBWE_BSC = 0.9304f0
function wsbwe_bernp(w::WsbweState, xmean::Float32, xmxvar::Float32)::Float32
    (xmean < 0.001f0 || xmean > 0.999f0) && return xmean
    xm = xmean > 0.5f0 ? (1.0f0 - xmean) : xmean
    var = (WSBWE_BSC*xmxvar)*(WSBWE_BC*wsbwe_pow(xm*(WSBWE_BC-1.0f0), WSBWE_BC-1.0f0)*
          wsbwe_exp(-wsbwe_pow((WSBWE_BC-1.0f0)*xm, WSBWE_BC)))
    ww = (xmean*(1.0f0-xmean)/var) - 1.0f0
    v = xmean*ww; ww = (1.0f0-xmean)*ww
    (x, kode) = wsbwe_bebet(w, v, ww)
    return kode < 0 ? xmean : x
end

# --- BWESLP linear interpolation (bweslp.f), BWECRC crown-class (bwecrc.f) ---
function wsbwe_bweslp(xx::Float32, x, y, n::Int)::Float32
    @inbounds for i in 1:n-1
        (xx < x[i] || xx > x[i+1]) && continue
        return y[i] + ((y[i+1]-y[i])/(x[i+1]-x[i]))*(xx-x[i])
    end
    r = y[n]; xx < x[1] && (r = y[1]); return r
end
@inline function wsbwe_bwecrc(xht::Float32)
    ihtc = 1
    if !(xht < 23.0f0); ihtc = 2; !(xht < 46.0f0) && (ihtc = 3); end
    (ihtc, ihtc*3-2, ihtc*3)
end

# -----------------------------------------------------------------------------
# BWEGO gate + wsbwe_apply! seam (mirrors the DFB LDFBGO / MPB seam).
# -----------------------------------------------------------------------------
"""
    wsbwe_go(w) -> Bool

FVS `BWEGO` per-cycle gate (bwego.f) for the MANUAL-DEFOL branch: fires when a WSBW
block is active with `LDEFOL` and at least one DEFOL activity is scheduled. The
`LCALBW/LBUDL` regional-outbreak (BUDLITE) branch is deferred (returns false).
"""
@inline wsbwe_go(w::WsbweState)::Bool =
    w.active && w.ldefol && !isempty(w.defol_sched)

# DEFOL crown-index range (bwesin.f) for a crown code (0..15) → (icrc1,icrc2,step).
@inline function wsbwe_crown_range(crn::Int)
    crn <= 0 && return (1,9,1)
    crn <= 9 && return (crn,crn,1)
    crn == 10 && return (1,3,1)
    crn == 11 && return (4,6,1)
    crn == 12 && return (7,9,1)
    crn == 13 && return (1,7,3)
    crn == 14 && return (2,8,3)
    return (3,9,3)   # 15
end

"""
    wsbwe_feeder(w, ns, sp, ht, dbh, dg, icr, prob, ibwspm, ibiomp,
                 elev, oldtpa, ormsqd, fint, ifint, iy_start) -> (PRBIO, PEDDS, PEHTG, AVYRMX, IFHOST)

The full BWESIT→BWEBMS→BWEADV→[year loop: BWEDR/BWEDEF/BWEAGE/BWEDAM] feeder for the
manual-DEFOL single-outbreak path. Pure — takes the FVS tree list + stand scalars,
returns the four arrays that feed BWEPDM (all validated bit-exact by dump-replay).
Multi-cycle POFPOT carryover (BWEPRB) is not exercised on the validated fixture and
is DEFERRED; the BWESIT LSKBIO reset (IY(ICYC) > IPRBYR) applies for a fresh outbreak.
"""
function wsbwe_feeder(w::WsbweState, ns::Int, sp, ht, dbh, dg, icr, prob,
                      ibwspm, ibiomp, elev::Float32, oldtpa::Float32, ormsqd::Float32,
                      fint::Float32, ifint::Int, iy_start::Int)
    # --- BWEBMS (ICVOPT=2) foliage biomass WK4 per tree ---
    alntpa = wsbwe_log(oldtpa)
    WK4 = zeros(Float32, ns)
    @inbounds for i in 1:ns
        D = dbh[i]; ispi = ibiomp[sp[i]]; H = ht[i]
        CL = (Float32(icr[i])*H)/100.0f0
        if D < 3.5f0
            WK4[i] = wsbwe_exp(WSBWE_BINT12[ispi] + WSBWE_BCL12[ispi]*wsbwe_log(CL) -
                     0.12975f0*alntpa + 0.40350f0*wsbwe_log(H)) * 1.13178f0
        else
            RD = D/ormsqd
            DDS = (2.0f0*D*dg[i] + dg[i]*dg[i])/fint
            WK4[i] = wsbwe_exp(WSBWE_BINT2[ispi] + 1.468547f0*wsbwe_log(D) +
                     0.308847f0*wsbwe_log(DDS) - 1.077047f0*wsbwe_log(H) +
                     0.690825f0*wsbwe_log(CL) - 0.142096f0*alntpa + 0.399244f0*wsbwe_log(RD))
        end
    end
    # --- BWESIT accumulation ---
    FOLPOT = zeros(Float32,6,9,4); FOLADJ = zeros(Float32,6,9,4)
    BWTPHA = zeros(Float32,7,3);  IFHOST = zeros(Int,7)
    @inbounds for i in 1:ns; IFHOST[ibwspm[sp[i]]] = 1; end
    @inbounds for i in 1:ns
        (ihtc,icrc1,icrc2) = wsbwe_bwecrc(ht[i]); iszi = ihtc
        PROBI = prob[i]*2.47103f0
        ihost = ibwspm[sp[i]]
        BWTPHA[ihost,iszi] += PROBI
        BIO = WK4[i]*453.6f0*PROBI
        if ihost < 7
            for ic in icrc1:icrc2, ia in 1:4
                FOLPOT[ihost,ic,ia] += WSBWE_THEOFL[ia,ic]*WSBWE_PRCRN3[ic]*BIO
            end
        end   # nonhost FOLNH not needed downstream for the DEFOL apply
    end
    for ihost in 1:6
        IFHOST[ihost]==0 && continue
        for ic in 1:9
            iszi = div(ic+2,3); DIV = BWTPHA[ihost,iszi]
            for ia in 1:4
                FOLPOT[ihost,ic,ia] = DIV > 0.0f0 ? FOLPOT[ihost,ic,ia]/DIV : 0.0f0
            end
        end
    end
    for ihost in 1:6, ic in 1:9, ia in 1:4
        IFHOST[ihost]==0 && continue
        FOLADJ[ihost,ic,ia] = FOLPOT[ihost,ic,ia]*1.0f0   # ×POFPOT(=1 on reset)
    end
    IFHOST[6] = 0   # BWEADV folds larch into nonhost then drops it as a host
    FNEW=zeros(Float32,9,6); FOLD1=zeros(Float32,9,6); FOLD2=zeros(Float32,9,6); FREM=zeros(Float32,9,6)
    PRBIO=ones(Float32,6,9,4)
    for ihost in 1:6
        IFHOST[ihost]==0 && continue
        for ic in 1:9
            FNEW[ic,ihost]=FOLADJ[ihost,ic,1]; FOLD1[ic,ihost]=FOLADJ[ihost,ic,2]
            FOLD2[ic,ihost]=FOLADJ[ihost,ic,3]; FREM[ic,ihost]=FOLADJ[ihost,ic,4]
        end
    end
    # --- year loop ---
    RDDSM1=ones(Float32,6,3); RHTGM1=ones(Float32,6,3)
    PEDDS=zeros(Float32,6,3); PEHTG=zeros(Float32,6,3)
    AVYRMX=zeros(Float32,6,3); BWMXCD=zeros(Float32,6,3)
    CUMDEF=zeros(Float32,6,3,5); APRBYR=zeros(Float32,6,3,2,5); CDEF=zeros(Float32,6,3)
    NCUMYR=0; ICUMYR=0
    BWFINT = fint < 1.0f0 ? 1.0f0 : fint
    IBWYR1=iy_start; IBWYR2=iy_start+ifint-1
    for IYRCUR in iy_start:(iy_start+ifint-1)
        # BWEDR TOTP
        TOTP=zeros(Float32,6,3)
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for ic in 1:9
                iszi=div(ic+2,3)
                for ia in 1:3
                    TOTP[ihost,iszi]=(TOTP[ihost,iszi]+FOLPOT[ihost,ic,ia]*0.6f0)+FOLADJ[ihost,ic,ia]*0.4f0
                end
                TOTP[ihost,iszi]+=FREM[ic,ihost]
            end
        end
        # BWEDEF (scheduled DEFOLs firing this year)
        for d in w.defol_sched
            Int(trunc(d[1])) == IYRCUR || continue
            spc = Int(trunc(d[2]))
            is1 = spc <= 0 ? 1 : ibwspm[spc]
            is2 = spc <= 0 ? 5 : is1
            is1 > 6 && continue
            (c1,c2,c3) = wsbwe_crown_range(Int(trunc(d[3])))
            pnew=d[4]; p1=d[5]; p2=d[6]; prem=d[7]
            for ih in is1:is2
                IFHOST[ih]==0 && continue
                for ic in c1:c3:c2
                    FNEW[ic,ih]  *= (1.0f0-(pnew/100.0f0))
                    FOLD1[ic,ih] *= (1.0f0-(p1/100.0f0))
                    FOLD2[ic,ih] *= (1.0f0-(p2/100.0f0))
                    FREM[ic,ih]  *= (1.0f0-(prem/100.0f0))
                end
            end
        end
        # BWEAGE
        TOTR=zeros(Float32,6,3)
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for ic in 1:9
                iszi=div(ic+2,3)
                DIV=FOLADJ[ihost,ic,1]; PRBIO[ihost,ic,1]=0.0f0
                FNEW[ic,ihost]<0.00001f0 && (FNEW[ic,ihost]=0.0f0)
                DIV>0.00001f0 && (PRBIO[ihost,ic,1]=FNEW[ic,ihost]/DIV)
                DIV=FOLADJ[ihost,ic,2]; PRBIO[ihost,ic,2]=0.0f0
                FOLD1[ic,ihost]<0.00001f0 && (FOLD1[ic,ihost]=0.0f0)
                DIV>0.00001f0 && (PRBIO[ihost,ic,2]=FOLD1[ic,ihost]/DIV)
                DIV=FOLADJ[ihost,ic,3]; PRBIO[ihost,ic,3]=0.0f0
                FOLD2[ic,ihost]<0.00001f0 && (FOLD2[ic,ihost]=0.0f0)
                DIV>0.00001f0 && (PRBIO[ihost,ic,3]=FOLD2[ic,ihost]/DIV)
                DIV=FOLADJ[ihost,ic,4]; PRBIO[ihost,ic,4]=0.0f0
                FREM[ic,ihost]<0.00001f0 && (FREM[ic,ihost]=0.0f0)
                DIV>0.00001f0 && (PRBIO[ihost,ic,4]=FREM[ic,ihost]/DIV)
                if BWTPHA[ihost,iszi]>0.0f0
                    TOTR[ihost,iszi]=(((TOTR[ihost,iszi]+FNEW[ic,ihost])+FOLD1[ic,ihost])+FOLD2[ic,ihost])+FREM[ic,ihost]
                end
            end
        end
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for ic in 1:9
                IREM=1; PRB=PRBIO[ihost,ic,1]
                if PRB < 0.80f0
                    IREM=2; iszi=div(ic+2,3)
                    if NCUMYR>0
                        ICMYR=NCUMYR; ICUMYR>0 && (ICMYR=ICUMYR)
                        CUM = ICMYR>0 ? CUMDEF[ihost,iszi,ICMYR] : 0.0f0
                        CUM>=40.0f0 && (IREM=3)
                    end
                end
                FR=FREM[ic,ihost]; PR=FOLPOT[ihost,ic,4]; PRALL=PR
                for ia in 1:3; PRALL+=FOLPOT[ihost,ic,ia]; end
                XMULT=0.0f0; DIV=FOLPOT[ihost,ic,3]
                DIV>0.000001f0 && (XMULT=FOLPOT[ihost,ic,4]/DIV)
                FA=FOLD2[ic,ihost]*XMULT
                if IREM==1
                    DIF=FR-PR
                    if DIF>=0.0f0
                        FA=PR; (DIF/PR>0.10f0) && (FA=FR-(DIF*0.40f0)); FR=FA
                    else
                        FR=FR+FA; FR>PR && (FR=PR)
                    end
                elseif IREM==2
                    FR=FR+FA*wsbwe_bweslp(PRB,WSBWE_RECRX,WSBWE_RECRY,4)*WSBWE_RECHST[ihost]
                    FR>PRALL && (FR=PRALL)
                else
                    FR=(FR+FA*wsbwe_bweslp(PRB,WSBWE_RECRX,WSBWE_RECRY,4))*0.85f0
                    FR>PRALL && (FR=PRALL)
                end
                FREM[ic,ihost]=FR
                XMULT=0.0f0; DIV=FOLPOT[ihost,ic,2]
                DIV>0.00001f0 && (XMULT=FOLPOT[ihost,ic,3]/DIV)
                FOLD2[ic,ihost]=FOLD1[ic,ihost]*XMULT
                XMULT=0.0f0; DIV=FOLPOT[ihost,ic,1]
                DIV>0.00001f0 && (XMULT=FOLPOT[ihost,ic,2]/DIV)
                FOLD1[ic,ihost]=FNEW[ic,ihost]*XMULT
                XMULT=0.0f0; DIV=FOLPOT[ihost,ic,2]
                DIV>0.00001f0 && (XMULT=FOLD1[ic,ihost]/DIV)
                FNEW[ic,ihost]=FOLPOT[ihost,ic,1]*wsbwe_bweslp((PRBIO[ihost,ic,2]*0.6f0)+(XMULT*0.4f0),
                    view(WSBWE_RELFX,:,ic), view(WSBWE_RELFY,:,ic), 4)
                FAa=FOLADJ[ihost,ic,4]+FOLADJ[ihost,ic,3]
                PRp=FOLPOT[ihost,ic,4]; FAa>PRp && (FAa=PRp)
                FRp=FREM[ic,ihost]; FAa<FRp && (FAa=FRp)
                FOLADJ[ihost,ic,4]=FAa
                for ii in 1:2
                    ia=4-ii; XMULT=0.0f0; DIV=FOLPOT[ihost,ic,ia-1]
                    DIV>0.00001f0 && (XMULT=FOLPOT[ihost,ic,ia]/DIV)
                    FOLADJ[ihost,ic,ia]=FOLADJ[ihost,ic,ia-1]*XMULT
                end
                FOLADJ[ihost,ic,1]=FNEW[ic,ihost]
            end
        end
        # BWEDAM
        AVPRBO=zeros(Float32,6,3,2); CDEF.=0.0f0
        STREES=0.0f0
        for ihost in 1:6
            IFHOST[ihost]!=1 && continue
            for iszi in 1:3; STREES+=BWTPHA[ihost,iszi]; end
        end
        STREES==0.0f0 && continue
        if NCUMYR<5; NCUMYR+=1; ICUMYR=NCUMYR else; ICUMYR+=1; ICUMYR>5 && (ICUMYR=1) end
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for iszi in 1:3
                BWTPHA[ihost,iszi]<=0.0f0 && continue
                DIV=TOTP[ihost,iszi]; CU=0.0f0
                DIV>0.00001f0 && (CU=100.0f0-((TOTR[ihost,iszi]/DIV)*100.0f0))
                CU<0.0f0 && (CU=0.0f0)
                CUMDEF[ihost,iszi,ICUMYR]=CU
                APRBYR[ihost,iszi,1,ICUMYR]=PRBIO[ihost,3*iszi-2,1]
                APRBYR[ihost,iszi,2,ICUMYR]=0.0f0
            end
            for ic in 1:9
                iszi=div(ic+2,3); BWTPHA[ihost,iszi]<=0.0f0 && continue
                APRBYR[ihost,iszi,2,ICUMYR]+=PRBIO[ihost,ic,1]*0.3333333f0
            end
        end
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for iszi in 1:3
                BWTPHA[ihost,iszi]<=0.0f0 && continue
                XMULT=1.0f0/Float32(NCUMYR)
                for i in 1:NCUMYR
                    CDEF[ihost,iszi]+=CUMDEF[ihost,iszi,i]
                    for it in 1:2; AVPRBO[ihost,iszi,it]+=APRBYR[ihost,iszi,it,i]*XMULT; end
                end
            end
        end
        if NCUMYR>=4
            for ihost in 1:6
                IFHOST[ihost]==0 && continue
                for iszi in 1:3
                    BWTPHA[ihost,iszi]<=0.0f0 && continue
                    AVYRMX[ihost,iszi] < 1.0f0-AVPRBO[ihost,iszi,1] && (AVYRMX[ihost,iszi]=1.0f0-AVPRBO[ihost,iszi,1])
                end
            end
        end
        XFIFTH=0.0f0; INRUN=IYRCUR-IBWYR1+1; IFIFTH=INRUN%5
        IFIFTH==0 && (XFIFTH=5.0f0/BWFINT)
        (IFIFTH!=0 && IYRCUR==IBWYR2) && (XFIFTH=Float32(IFIFTH)/BWFINT)
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for iszi in 1:3
                BWTPHA[ihost,iszi]<=0.0f0 && continue
                RDDS=wsbwe_rdds(AVPRBO[ihost,iszi,2], RDDSM1[ihost,iszi]); RDDSM1[ihost,iszi]=RDDS
                RHTG=wsbwe_rhtg(AVPRBO[ihost,iszi,1], RHTGM1[ihost,iszi]); RHTGM1[ihost,iszi]=RHTG
                PEDDS[ihost,iszi]+=RDDS/BWFINT
                if iszi>1
                    PEHTG[ihost,iszi]+=RHTG/BWFINT
                else
                    RHTGs=1.0f0
                    AVPRBO[ihost,iszi,2]<0.98f0 && (RHTGs=wsbwe_exp(WSBWE_STHTGR[ihost]*(1.0f0-AVPRBO[ihost,iszi,1])))
                    PEHTG[ihost,iszi]+=RHTGs*XFIFTH
                end
                if INRUN>=5 || IYRCUR==IBWYR2
                    BWMXCD[ihost,iszi]<CDEF[ihost,iszi] && (BWMXCD[ihost,iszi]=CDEF[ihost,iszi])
                end
            end
        end
    end
    return (PRBIO, PEDDS, PEHTG, AVYRMX, IFHOST)
end

"""
    wsbwe_apply!(s, old_tpa, fint)

FVS `BWECUP` seam (gradd.f:108 `IF (IPMODI==1 .AND. LBWEGO) CALL BWECUP`). LIVE for the
manual-DEFOL path (EM host coefficients): runs the validated feeder (`wsbwe_feeder`)
then the BWEPDM per-tree application — topkill (BWERNP AVDEF gate + Marsden logistic +
HTGSTP truncation), DG reduction (`DDS·BWERNP(PEDDS,.03)`), HTG reduction, and WK2
mortality (`max(BASE, PROB·(1−survival))`) — drawing the damage RNG via BWERNP/BWEBET
in FVS ISCT/IND1 order (the damage seed DSEEDD=55329 lives in `w.rng_s0`, saved/restored
across cycles by the persistent state, matching BWERPT/BWERGT). Every routine here is
dump-replay bit-exact vs FVSem_wsbwe (host_defol.key). Early-returns leave the stand
byte-identical: inactive state, the BUDLITE/GENDEFOL branch (`lbudl`, deferred), a
non-EM variant (host block-data not yet ported), or a stand with no budworm host.
"""
function wsbwe_apply!(s::StandState, old_tpa, fint)
    w = s.wsbwe
    (w === nothing || !(w::WsbweState).active) && return nothing
    ww = w::WsbweState
    ww.lbudl && return nothing                 # BUDLITE/GENDEFOL deferred
    ibwspm = wsbwe_ibwspm_for(s.variant)             # per-variant host-class map (EM, TT ported+validated)
    ibwspm === nothing && return nothing             # unsupported variant → inert
    ibiomp = wsbwe_ibiomp_for(s.variant)
    # LIVE gate (const, normally true). The feeder is end-to-end bit-exact vs FVS<v>_wsbwe after
    # the two faithful fixes below (OLDTPA/ORMSQD=TPROB/RMSQD, and SPDECD species decode). Kept as a
    # switch for A/B byte-identity proofs. See scratchpad/wsbwe/HANDOFF.md.
    WSBWE_APPLY_LIVE || return nothing
    t = s.trees; ns = t.n
    ns <= 0 && return nothing
    # host present?
    hashost = false
    @inbounds for i in 1:ns
        (ibwspm[t.species[i]] < 7) && (hashost = true; break)
    end
    hashost || return nothing
    # stand scalars
    elev   = s.plot.elevation
    # OLDTPA/ORMSQD — grincr.f:281/285 set these to the CURRENT cycle-start stand density
    # (OLDTPA=TPROB, ORMSQD=RMSQD) via DENSE (dense.f:182/250): TPROB=ΣP over all trees,
    # RMSQD=sqrt(ΣD²·P/TPROB). These are what BWEBMS (bwebmsem.f:124/130) reads. They are NOT
    # `s.plot.old_tpa`/`old_qmd` (the PREVIOUS-cycle stored scalars — 0 at cycle 1, which drove
    # ALOG(0)=-Inf and D/0=+Inf ⇒ NaN foliage biomass). P = cycle-start tpa = the `old_tpa` arg
    # (t.tpa at cycle start, pre-MORTS), matching DENSE's PROB(I). Bit-exact to the oracle
    # (589.6528/5.1449676 vs FVSem_wsbwe OLDTPA/ORMSQD).
    tprob = 0.0f0; tsumd2 = 0.0f0
    @inbounds for i in 1:ns
        p = Float32(old_tpa[i]); tprob += p; tsumd2 += t.dbh[i]*t.dbh[i]*p
    end
    oldtpaS = tprob
    ormsqd  = tprob > 0.0f0 ? sqrt(tsumd2/tprob) : 0.0f0
    fintf  = Float32(fint)
    ifint  = round(Int, Float64(fint))                  # IFINT — cycle length in yr (the fint arg; s.plot.forecast_interval is 0 here)
    iy_st  = Int(cycle_year_at(s.control, Int(s.control.cycle))) # IY(ICYC) — s.control.cycle is 0-based
    # BWEGO (bwego.f:88): `CALL OPFIND(1,2151,I); LDEFOL=I.GT.0` — the manual-DEFOL gate fires ONLY
    # in a cycle that actually has a DEFOL activity scheduled within its year window [IY(ICYC),
    # IY(ICYC)+IFINT-1]. Without this, apply! runs (and draws the damage RNG) every cycle after the
    # outbreak, where the oracle never calls BWECUP — non-faithful and it consumes the RNG stream.
    # A byte-identical no-op on any cycle with no scheduled DEFOL year (must precede any RNG draw).
    lastyr_cyc = iy_st + ifint - 1
    any(d -> (iy_st <= Int(trunc(d[1])) <= lastyr_cyc), ww.defol_sched) || return nothing
    # --- FEEDER ---
    (PRBIO, PEDDS, PEHTG, AVYRMX, IFHOST) = wsbwe_feeder(ww, ns,
        t.species, t.height, t.dbh, t.diam_growth, t.crown_pct, t.tpa,
        ibwspm, ibiomp, elev, oldtpaS, ormsqd, fintf, ifint, iy_st)
    # NOBWYR/IBWYR/FA — the outbreak spans the whole cycle on the DEFOL path (FA=0)
    IBWYR = ifint; NOBWYR = 0; FA = 0.0f0
    # --- point basal areas (PNTBA all, PNTHBA host≥? bwepdm sums all; host uses IBWSPM<6) ---
    npt = Int(s.plot.points_inv)                        # IPTINV
    PNTBA = zeros(Float32, max(npt,1)); PNTHBA = zeros(Float32, max(npt,1))
    @inbounds for i in 1:ns
        p = Int(t.plot_id[i]); d = t.dbh[i]
        PNTBA[p] += d*d*0.005454154f0
        ibwspm[t.species[i]] < 6 && (PNTHBA[p] += d*d*0.005454154f0)
    end
    # --- BWEPDM per-tree, in FVS ISCT/IND1 order (damage seed restored to w.rng_s0) ---
    isnan(ww.rng_s0) && (ww.rng_s0 = Float64(ww.dseed))   # BWERPT(DSEEDD)
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    MAXSP = length(ibwspm)
    @inbounds for ispi in 1:MAXSP
        isct[ispi,1] == 0 && continue
        ihost = ibwspm[ispi]
        ihost >= 6 && continue                 # nonhost or larch bypass
        for ii in isct[ispi,1]:isct[ispi,2]
            i = Int(ind1[ii])
            H = t.height[i]; DBH = t.dbh[i]; DGI = t.diam_growth[i]; HTGI = t.ht_growth[i]
            BARK = bark_ratio(ba_a, ba_b, t.species[i], DBH)
            (ihtc, ICRC1, _) = wsbwe_bwecrc(H); ISZI = ihtc
            IC2 = ICRC1 + 1
            MFT = 10.0f0 - Float32(wsbwe_ifix((((PRBIO[ihost,ICRC1,1]*0.25f0)+(PRBIO[ihost,ICRC1,2]*0.25f0)+
                    (PRBIO[ihost,ICRC1,3]*0.25f0)+(PRBIO[ihost,ICRC1,4]*0.25f0))*10.0f0)+0.5f0))
            MFM = 10.0f0 - Float32(wsbwe_ifix((((PRBIO[ihost,IC2,1]*0.25f0)+(PRBIO[ihost,IC2,2]*0.25f0)+
                    (PRBIO[ihost,IC2,3]*0.25f0)+(PRBIO[ihost,IC2,4]*0.25f0))*10.0f0)+0.5f0))
            ITR = Int(t.trunc[i]); NRM = Int(t.norm_ht[i])
            PCTK = ITR > 0 ? Float32(ITR)/Float32(NRM) : 0.0f0
            IMC0 = Int(t.mort_code[i]); ICR0 = Int(t.crown_pct[i])
            HTG = HTGI
            # --- topkill (LTOPK is always on) ---
            AVDEF = wsbwe_bernp(ww, AVYRMX[ihost,ISZI], 0.06f0)
            if ISZI == 1
                X = H > 10.0f0 ? 10.0f0 : H
                PRTOPK = 1.0f0/(1.0f0+wsbwe_exp(-(-2.5817f0 - 0.027635f0*Float32(ICR0) +
                          3.709f0*sqrt(AVDEF) + 0.0488f0*X)))
            else
                PRTOPK = WSBWE_IFIR[ihost]==0 ?
                    0.96f0*(1.0f0-wsbwe_exp(-(wsbwe_pow(0.65f0*(AVDEF+1.0f0),14.0f0)))) :
                    0.90f0*(1.0f0-wsbwe_exp(-(wsbwe_pow(0.60f0*(AVDEF+1.0f0),11.0f0))))
            end
            Xr = wsbwe_rand!(ww)
            if Xr < PRTOPK
                PART = (MFT>0.0f0 && MFM>0.0f0) ?
                    1.0f0/(1.0f0+wsbwe_exp(WSBWE_TK0[ihost]+WSBWE_TK1[ihost]*elev+WSBWE_TK2[ihost]*elev*elev+
                        WSBWE_TK3[ihost]*Float32(t.dmr[i])+WSBWE_TK4[ihost]*PCTK+WSBWE_TK5[ihost]*MFT+
                        WSBWE_TK6[ihost]*MFM+WSBWE_TK7[ihost]*MFM*PCTK+WSBWE_TK8[ihost]*MFT*MFM)) : 0.0f0
                PART > 0.9f0 && (PART = 0.9f0)
                if PART > 0.0f0
                    FTKILL = H*PART; TOPH = H - FTKILL
                    ITRC2 = wsbwe_ifix(TOPH*100.0f0+0.5f0)
                    if ITR > 0
                        ITR > ITRC2 && (ITR = ITRC2)
                        H = TOPH
                    else
                        D = DBH*BARK
                        if H >= 25.0f0 && D >= 6.0f0
                            AFv = t.cuft_vol[i]/(0.00545415f0*D*D*H)
                            AFv = 0.44244f0 - (0.99167f0/AFv) - (1.43237f0*wsbwe_log(AFv)) +
                                  (1.68581f0*sqrt(AFv)) - (0.13611f0*AFv*AFv)
                            DTK = FTKILL/H; DTK = (DTK/((AFv*DTK)+(1.0f0-AFv)))*D
                            if DTK > 4.0f0
                                ITR = ITRC2; NRM = wsbwe_ifix(H*100.0f0+0.5f0); IMC0 = 3
                            elseif DTK > 2.0f0 && IMC0 < 2
                                IMC0 = 2
                            end
                        end
                        H = TOPH
                        if ICR0 >= 0
                            CN = (Float32(ICR0)/100.0f0*H) - H + TOPH
                            NEW = wsbwe_ifix(CN/TOPH*100.0f0+0.5f0); NEW < 5 && (NEW = 5)
                            ICR0 = -NEW
                        end
                    end
                end
                HTG = 0.0f0
            end
            # --- DG reduction ---
            DDS = DGI*(2.0f0*BARK*DBH+DGI)
            XD = wsbwe_bernp(ww, PEDDS[ihost,ISZI], 0.03f0)
            DDS *= XD
            DGn = sqrt((DBH*BARK)^2 + DDS) - BARK*DBH
            DGn < 0.0f0 && (DGn = 0.0f0)
            # --- HTG reduction ---
            if HTG > 0.0f0
                if PEDDS[ihost,ISZI] < 0.99f0
                    XH = PEHTG[ihost,ISZI]+((1.0f0-PEHTG[ihost,ISZI])/(1.0f0-PEDDS[ihost,ISZI])*(XD-PEDDS[ihost,ISZI]))
                else
                    XH = wsbwe_bernp(ww, PEHTG[ihost,ISZI], 0.03f0)
                end
                HTG *= XH; HTG < 0.0f0 && (HTG = 0.0f0)
            end
            # --- mortality (Marsden) ---
            KTK = (ITR==0 || NRM==0) ? 0.0f0 :
                  10.0f0 - Float32(wsbwe_ifix(Float32(ITR)/Float32(NRM)*10.0f0 + 0.5f0))
            PR = wsbwe_mort_pr(ihost, elev, PNTBA[Int(t.plot_id[i])], PNTHBA[Int(t.plot_id[i])], MFT, MFM, KTK)
            PROBv = old_tpa[i]; WK20 = old_tpa[i] - t.tpa[i]
            BASE = WK20/PROBv
            PR = 1.0f0 - PR
            PR = NOBWYR > 0 ? wsbwe_powi(PR, IBWYR)*wsbwe_pow(1.0f0-BASE, FA) : wsbwe_powi(PR, IBWYR)
            PR = 1.0f0 - PR
            if BASE > PR
                # keep background mortality (tpa unchanged)
            else
                PR > 0.98f0 && (PR = 0.98f0)
                WK2n = PROBv*PR
                t.tpa[i] = PROBv - WK2n
            end
            # write back growth / topkill state
            t.diam_growth[i] = DGn
            t.ht_growth[i]   = HTG
            t.height[i]      = H
            t.crown_pct[i]   = Int32(ICR0)
            t.trunc[i]       = Int32(ITR)
            t.norm_ht[i]     = Int32(NRM)
            t.mort_code[i]   = Int32(IMC0)
        end
    end
    return nothing
end
