# =============================================================================
# Westwide Pine Beetle (WWPB) model — wwpb/*.f  (internal prefix "bm")
# =============================================================================
# Fourth port of the FVS insect/pathogen event-extension family (after DFB,
# DFTM, WPBR). WWPB is the "Westwide Pine Beetle" bark-beetle model — Mountain
# Pine Beetle (MPB) / Western Pine Beetle (WPB) / Ips on lodgepole (sp 7) and
# ponderosa (sp 10) pine.
#
# ARCHITECTURE — CRITICAL, and UNLIKE the prior three:
# ----------------------------------------------------------------------------
# DFB/DFTM/WPBR are stand/tree-level extensions: a single keyword block
# activates them and a per-cycle seam (called from grincr/gradd) applies their
# mortality to the single-stand treelist. WWPB is NOT stand-level — it is a
# *landscape* Parallel-Processing-Extension (PPE) model. Its outbreak runs over
# MANY stands with spatial location/area data and beetle dispersal *between*
# stands. Its structure (per the pristine wwpb/*.f):
#
#   * Two independent keyword entry points:
#       - BMIN   (wwpb/bmin.f)   — the STAND-level keyword block, keywds.f
#                 OPTION 126 (base/keywds.f TABLE(126)='BMIN'; initre.f
#                 `CALL BMIN(LKECHO)`). Its TABLE has ONLY output-scheduling
#                 keywords: MAINOUT/TREEOUT/BKPOUT/VOLOUT/END. Each opens a
#                 .bmm/.bmt/.bmb/.bmv report file and schedules an OPNEW output
#                 activity (MYACT 2701..2704). It does NOT activate the outbreak.
#       - BMPPIN (wwpb/bmppin.f) — the PPE-LANDSCAPE keyword block, read from
#                 PPIN (the parallel-processing input driver). 40 keywords; the
#                 activation keyword is DISPERSE (GPNEW activity 301 → outbreak
#                 start year IBMYR1 + duration → IBMYR2). Beetle species host
#                 designations (HOST/PBSPEC), pheromone/salvage/sanitation mgmt,
#                 RANNSEED reseed, etc. all live here.
#   * Driver BMDRV (wwpb/bmdrv.f) — called from ALSTD2 (the PPE all-stands
#                 annual driver). Runs only `IF (LBMSPR .AND. IBMYR1 .GT. 0)`.
#                 LBMSPR ("spread mode on") is set true in bmsetp.f only when
#                 SPLAEX/SPLAAR supply per-stand spatial location + area.
#   * Mortality hand-back BMKILL (wwpb/bmkill.f) — from PPMAIN (not BMDRV):
#                 converts the accumulated TPBK ledger to per-record FVS
#                 mortality WK2(I). This is the single quantity fed to growth.
#   * RNG BMRANN (wwpb/bmrann.f) — the model's own MINSTD LCG (see below).
#
# BUILD REALITY (measured): every one of the 24 FVS*_buildDir variants links the
# NO-OP stub base/exbm.f (its entries BMIN/BMDRV/BMSETP/BMPPIN/BMKILL/… all
# RETURN or emit "**NO BM"). NO sourceList references wwpb/*.f, and the PPE
# framework (PPIN/PPMAIN/ALSTD1/ALSTD2/SPLAEX/GPGET/GPNEW; PPEPRM.F77 lives only
# under archive/PPEcommons/) is ABSENT from the tree. Consequently the beetle
# *outbreak* is unreachable via any shipped or straightforwardly-relinkable
# single-stand run — a full WWPB port would first have to port/stub the entire
# PPE spatial multi-stand harness. This is a materially larger undertaking than
# DFB/DFTM/WPBR and is deferred; see the handoff in scratchpad/wwpb/GOLDENS.md.
#
# THIS PORT (beachhead — the reachable BMIN block + the RNG, INERT seam):
#   * `kw_wwpbin!`  — faithful BMIN block reader (MAINOUT/TREEOUT/BKPOUT/VOLOUT/
#                     END) → `s.wwpb` (WwpbState). FVSjl runs no PPE loop and
#                     emits no .bm* files, so scheduling these output activities
#                     is faithfully INERT: a stand with a BMIN block projects
#                     BYTE-IDENTICALLY to one without (matching the fact that even
#                     the relinked model needs spatial data + the PPE driver to
#                     do anything). No per-cycle engine seam is wired.
#   * `wwpb_rand!` / `wwpb_seed!` — the BMRANN LCG + BMRNSD reseed (bit-exact;
#                     golden-validated against pristine wwpb/bmrann.f via
#                     scratchpad/wwpb/driver_bmrann.f). NEVER FFI'd.
#   * `wwpb_defaults!` — BMINIT/BLOCK-DATA defaults (seed 55329, host size-class
#                     breakpoints UPSIZ, beetle-species selector PBSPEC).
#
# Float32 discipline: WWPB's Fortran is REAL (Float32). Every ported arithmetic
# result and coefficient array is Float32 to stay bit-identical to the model.
# =============================================================================

# -----------------------------------------------------------------------------
# BLOCK DATA (bmblkd*.f) / BMINIT defaults.
# -----------------------------------------------------------------------------
# UPSIZ(NSCL=10) — the DBH size-class upper breakpoints shared by every variant
# BLOCK DATA (bmblkdbm.f:73 …). Transcribed in DATA order.
const WWPB_UPSIZ_DEFAULT = Float32[3, 6, 9, 12, 15, 18, 21, 25, 30, 50]

# WPSIZ(MXDWSZ=3) — dead-wood pool DBH classes (bmblkd*.f).
const WWPB_WPSIZ_DEFAULT = Float32[10, 20, 60]

# ISCMIN(NPBSP) — smallest attractive size class per beetle {MPB, WPB, Ips}
# (bmblkd*.f:96): MPB/WPB won't enter stands < class 3 (6"); Ips needs class ≥ 2.
const WWPB_ISCMIN_DEFAULT = Int32[3, 3, 2]

const WWPB_DEFAULT_SEED = 55329.0f0   # bmblkd*.f DATA BMS0/55329D0/, BMSS/55329./ (odd)

# BMIN output-activity codes (bmin.f OPNEW MYACT), for reference / reproduction.
const WWPB_MYACT_MAINOUT = Int32(2701)
const WWPB_MYACT_TREEOUT = Int32(2702)
const WWPB_MYACT_BKPOUT  = Int32(2703)
const WWPB_MYACT_VOLOUT  = Int32(2704)

"""
    WwpbOutReq

One scheduled WWPB stand-level output request (bmin.f MAINOUT/TREEOUT/BKPOUT/
VOLOUT). `myact` is the OPNEW activity code (2701..2704); `idt` the start
date/cycle (field 1, default 1); `nyears` the number of years to print (PRMS(1),
default 100); `incr` the year increment (PRMS(2), default 5). In FVSjl these are
recorded faithfully but never fire — no PPE output loop exists — so they are
inert w.r.t. the projection and the .sum.
"""
struct WwpbOutReq
    myact::Int32
    idt::Int32
    nyears::Int32
    incr::Int32
end

"""
    WwpbState

Westwide Pine Beetle model state. Populated by `kw_wwpbin!` (the reachable BMIN
stand-level output block) and seeded with the BLOCK-DATA/BMINIT defaults by
`wwpb_defaults!`. `active` mirrors "a BMIN keyword block was read". No per-cycle
engine seam is wired (the outbreak lives in the absent PPE landscape driver), so
a stand carrying a WwpbState projects byte-identically to one without.

Fields beyond the BMIN outputs (`pbspec`, `upsiz`, `iscmin`, the RNG) are the
deterministic setup a future full port builds on; they are defaults here and are
not yet consumed by any seam.
"""
mutable struct WwpbState <: AbstractWwpbState
    active::Bool                       # a BMIN block was read (engine seam gate)
    # --- BMIN stand-level output requests (bmin.f) ---
    lbmain::Bool                       # LBMAIN — MAINOUT requested (.bmm)
    lbmtre::Bool                       # LBMTRE — TREEOUT requested (.bmt)
    lbmbkp::Bool                       # LBMBKP — BKPOUT  requested (.bmb)
    lbmvol::Bool                       # LBMVOL — VOLOUT  requested (.bmv)
    outreqs::Vector{WwpbOutReq}        # the OPNEW-scheduled output activities, in read order
    # --- BLOCK-DATA / BMINIT deterministic defaults (not yet seam-consumed) ---
    pbspec::Int32                      # PBSPEC — beetle selector 1=MPB 2=WPB 3=Ips 4=MPB+WPB (default 1)
    upsiz::Vector{Float32}             # UPSIZ(NSCL) — DBH size-class breakpoints
    iscmin::Vector{Int32}             # ISCMIN(NPBSP) — min attractive class per beetle
    # --- BMRANN LCG state (COMMON /BMRNCM/) ---
    rng_s0::Float64                    # BMS0 — current generator state (double, exact mod)
    rng_ss::Float32                    # BMSS — the reseed default (BMRNSD LSET=false resets to it)
    # --- outbreak activation (DISPERSE, the reconstructed BMPPIN seam) ---
    outbreak::Bool                     # DISPERSE seen → run the per-cycle outbreak (bmdrv/PPMAIN)
    iyr1::Int32                        # IBMYR1 — outbreak start year (0 = inactive, bmdrv.f:35 gate)
    iyr2::Int32                        # IBMYR2 — outbreak end year (= IBMYR1 + duration − 1)
    seed_class::Int32                  # synthetic inventory-damage seed: size class (BKP starts at 0)
    seed_tpa::Float32                  # synthetic inventory-damage seed: TPA beetle-killed (kick-off)
end

"""
    wwpb_defaults!(variant) -> WwpbState

FVS BMINIT (bminit.f) + BLOCK DATA (bmblkd*.f): the WWPB run-time defaults.
Called when the first BMIN keyword is seen. Seed 55329 (odd), PBSPEC=1 (MPB),
the shared UPSIZ/ISCMIN tables. `variant` is accepted for parity with the other
extensions (the per-variant BLOCK DATA differs only in the species-list length
and ISPFLL/HSPEC shapes, none of which the BMIN beachhead consumes).
"""
function wwpb_defaults!(variant)
    return WwpbState(
        false,                                  # active
        false, false, false, false,             # lbmain, lbmtre, lbmbkp, lbmvol
        WwpbOutReq[],                            # outreqs
        Int32(1),                               # pbspec (MPB)
        copy(WWPB_UPSIZ_DEFAULT), copy(WWPB_ISCMIN_DEFAULT),   # upsiz, iscmin
        Float64(WWPB_DEFAULT_SEED), WWPB_DEFAULT_SEED,         # rng_s0, rng_ss
        false, Int32(0), Int32(0), Int32(3), 0.0f0,           # outbreak, iyr1, iyr2, seed_class, seed_tpa
    )
end

# -----------------------------------------------------------------------------
# BMRANN (bmrann.f) — the WWPB model's own double-precision Lehmer/MINSTD LCG.
# -----------------------------------------------------------------------------
"""
    wwpb_rand!(w) -> Float32

FVS `BMRANN`: `BMS1 = DMOD(16807·BMS0, 2147483647)`; `SEL = REAL(BMS1/2147483648)`;
`BMS0 = BMS1`. `BMS0`/`BMS1` are DOUBLE, so the modular step is exact integer
arithmetic (values < 2^31); only the returned uniform `SEL` truncates to Float32.
IDENTICAL algorithm to DFTM's TMRANN and DFB's DFBRAN. Default seed 55329.
NEVER FFI'd. Golden-validated bit-exact (scratchpad/wwpb/driver_bmrann.f).
"""
@inline function wwpb_rand!(w::WwpbState)::Float32
    s1 = rem(16807.0 * w.rng_s0, 2147483647.0)   # DMOD (exact, values < 2^31)
    w.rng_s0 = s1
    return Float32(s1 / 2147483648.0)            # SEL = REAL(BMS1 / 2^31)
end

"""
    wwpb_seed!(w, seed, present)

FVS `BMRNSD(LSET, SEED)` (bmrann.f ENTRY). RANNSEED with a value (`present`):
force an ODD seed, store it as both the working state `BMS0` and the reseed
default `BMSS`. RANNSEED with no value (`present=false`): reset `BMS0` to `BMSS`.
"""
function wwpb_seed!(w::WwpbState, seed::Float32, present::Bool)
    if !present
        w.rng_s0 = Float64(w.rng_ss)             # LSET=false: SEED=BMSS; BMS0=BMSS
        return nothing
    end
    (seed % 2.0f0 == 0.0f0) && (seed += 1.0f0)   # AMOD(SEED,2)==0 → SEED+1 (odd)
    w.rng_ss = seed
    w.rng_s0 = Float64(seed)
    return nothing
end

# -----------------------------------------------------------------------------
# kw_wwpbin! (bmin.f) — WWPB stand-level keyword-block reader (keywds.f opt 126).
# -----------------------------------------------------------------------------
"""
    kw_wwpbin!(s, rec, kr)

Parse the `BMIN … END` block (wwpb/bmin.f). The block is entered by the
top-level BMIN keyword (keywds.f option 126, initre.f `CALL BMIN`). Its TABLE has
only the four output-scheduling keywords plus END:

* `MAINOUT` — stand main-summary report (.bmm), OPNEW activity 2701.
* `TREEOUT` — detailed by-size-class tree report (.bmt), activity 2702.
* `BKPOUT`  — detailed beetle-killing-potential/brood report (.bmb), activity 2703.
* `VOLOUT`  — detailed volume report (.bmv), activity 2704.
* `END`     — end of the single-stand WWPB option block.

Each output keyword takes 3 optional fields: field1 = start date/cycle `IDT`
(default 1); field2 = number of years `PRMS(1)` (default 100); field3 = year
increment `PRMS(2)` (default 5) — exactly as bmin.f reads them (INT-truncated).

FVSjl runs no PPE landscape loop and opens no .bm* files, so recording these
requests is faithfully INERT: a stand carrying a WwpbState projects
byte-identically to one without. (Even the relinked Fortran model produces no
output for these unless the absent PPE spatial driver runs — see the module
header.) No per-cycle engine seam is wired.
"""
function kw_wwpbin!(s::StandState, rec, kr::KeywordReader)
    s.wwpb === nothing && (s.wwpb = wwpb_defaults!(s.variant))
    w = s.wwpb
    w.active = true

    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"                       # option 5 — end of block
            break
        elseif k == "MAINOUT"               # option 1
            w.lbmain = true
            push!(w.outreqs, _wwpb_outreq(WWPB_MYACT_MAINOUT, r))
        elseif k == "TREEOUT"               # option 2
            w.lbmtre = true
            push!(w.outreqs, _wwpb_outreq(WWPB_MYACT_TREEOUT, r))
        elseif k == "BKPOUT"                # option 3
            w.lbmbkp = true
            push!(w.outreqs, _wwpb_outreq(WWPB_MYACT_BKPOUT, r))
        elseif k == "VOLOUT"                # option 4
            w.lbmvol = true
            push!(w.outreqs, _wwpb_outreq(WWPB_MYACT_VOLOUT, r))
        elseif k == "DISPERSE"              # activate the outbreak (reconstructed BMPPIN DISPERSE, GPNEW 301)
            # field1 = start year IBMYR1 (default 0 = this cycle); field2 = duration years
            # (→ IBMYR2 = IBMYR1 + dur − 1); field3 = seed size class; field4 = seed TPA
            # (the synthetic inventory-damage kick-off, since BKP starts at 0 and FVSjl
            # has no PPE Outside-World immigration).
            w.outbreak = true
            iyr1 = r.present[1] ? Int32(trunc(Int, r.values[1])) : Int32(1)
            dur  = r.present[2] ? Int32(trunc(Int, r.values[2])) : Int32(1)
            w.iyr1 = iyr1
            w.iyr2 = iyr1 + max(dur, Int32(1)) - Int32(1)
            r.present[3] && (w.seed_class = Int32(clamp(trunc(Int, r.values[3]), 1, WWPB_NSCL)))
            r.present[4] && (w.seed_tpa = Float32(r.values[4]))
        else
            # options 6..10 are blank placeholders in bmin.f (GOTO 10); any other
            # token is an unrecognized WWPB sub-keyword (FVS ERRGRO warns, skips).
        end
    end
    return nothing
end

# bmin.f field decode: IDT=INT(ARRAY(1)) def 1; PRMS(1)=INT(ARRAY(2)) def 100;
# PRMS(2)=INT(ARRAY(3)) def 5. INT truncates toward zero (trunc).
@inline function _wwpb_outreq(myact::Int32, r)::WwpbOutReq
    idt    = r.present[1] ? Int32(trunc(Int, r.values[1])) : Int32(1)
    nyears = r.present[2] ? Int32(trunc(Int, r.values[2])) : Int32(100)
    incr   = r.present[3] ? Int32(trunc(Int, r.values[3])) : Int32(5)
    return WwpbOutReq(myact, idt, nyears, incr)
end
