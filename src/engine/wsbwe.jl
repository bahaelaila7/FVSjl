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
            spp  = r.present[2] ? Float32(r.values[2]) : 0.0f0
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
