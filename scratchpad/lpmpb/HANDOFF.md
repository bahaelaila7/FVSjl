# LPMPB port — HANDOFF (commit + wiring + dependency-ordered next steps)

## ⚠ WHY THIS IS STAGED, NOT COMMITTED
This agent was launched isolated in a **FVS Fortran** worktree
(`/workspace/ForestVegetationSimulator/.claude/worktrees/agent-af98293a390d7ba6e`).
The harness git-isolation guard HARD-REFUSES every git operation (and every
`cd`/`-C`) against the FVSjl shared checkout (`/workspace/FVSjl`) — including from
FVSjl worktree subpaths — and cwd does not persist between bash calls, so cd+git
cannot be split. I therefore **could not create the `lpmpb-port` branch or commit**.
All deliverables are staged in `scratchpad/lpmpb/`, validated, ready to apply.
The port math is proven BIT-EXACT; only the git packaging is blocked.

## APPLY (from a FVSjl checkout on kt-variant-port @ b926d894)
1. `git switch -c lpmpb-port b926d894`
2. `cp scratchpad/lpmpb/lpmpb.jl  src/engine/lpmpb.jl`
   `cp scratchpad/lpmpb/test_lpmpb.jl  test/unit/test_lpmpb.jl`
3. Wire (4 one-line edits, mirroring dfb/dftm/wpbr):
   - `src/core/state.jl`  (~line 933, after AbstractWpbrState):
        `abstract type AbstractMpbState end`
     and in `mutable struct StandState` (~line 963, after `wpbr::…`):
        `mpb::Union{AbstractMpbState,Nothing}`
     …and add `nothing` to the StandState constructor(s) at the mpb position
     (grep the dfb/dftm/wpbr `nothing` initializers — same spots).
   - `src/FVSjl.jl`: `include("engine/lpmpb.jl")` next to the dfb/dftm/wpbr
     includes; export `mpb_colind, mpb_coldbh_start, mpb_colmod, mpb_prkill,
     mpb_er, mpb_idxlp, mpb_defaults!, MpbState, kw_mpbin!` if the test needs them
     (the test uses `FVSjl.mpb_*` — exports optional).
   - `src/engine/keyword_dispatch.jl` (~line 2445, after the BRUST elseif):
        `elseif kw == "MPB";      kw_mpbin!(s, rec, kr)   # Mountain Pine Beetle (keywds.f 'MPB'; MPBSTART/POPDYN/INITMORT/QVALUES/… → s.mpb + gated COL mortality seam)`
   - `src/engine/simulate.jl` (~line 654, next to the wpbr seam):
        ```
        if !tripled && s.mpb !== nothing && (s.mpb::MpbState).active
            mpb_apply!(s, old_tpa, fint)   # MPBGO→MPBCUP→COLDRV (gradd.f:63)
        end
        ```
     Place it with the other `!tripled` insect seams; `old_tpa`/`fint` are already
     in scope there (the dfb seam uses them).
4. Test: `JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project -e 'using Pkg;
   Pkg.test()'` (or run test/unit/test_lpmpb.jl). Provision gitignored *.tre first:
   `cp -n /workspace/FVSjl/test/harness/scenarios/*.tre <worktree>/test/harness/scenarios/`.
5. Commit trailers (each on its own line):
   `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`
   `Claude-Session: https://claude.ai/code/session_01X241b5qt3DpTtEtZcjUwRH`

## WHAT IS VALIDATED (commit-ready)
- The DEFAULT deterministic Cole rate-of-loss mortality path: COLIND, COLDBH,
  COLMOD (IBOUSE=0), COLMRT (PRKILL + XT), MPBER, mpb_idxlp — BIT-EXACT vs the
  relinked FVSie_lpmpb oracle (see GOLDENS.md; test_lpmpb.jl embeds the goldens).
- Keyword reader kw_mpbin! (state-setting for all COL-path options) + inert seam.

## 2026-08-17 UPDATE — RANSTART + CURRMORT/INVMORT MODES DONE (this diff = wiring.diff)
Base HEAD `74d4f5c5` already has lpmpb.jl merged; this change edits **only**
`src/engine/lpmpb.jl` + `test/unit/test_lpmpb.jl` (git-apply `wiring.diff`; APPLY-CHECK OK).
No new engine seam is added — the existing gated `mpb_apply!` seam (simulate.jl ~L672)
is reused; `mpb_apply!` was internally reordered (MPBER first, then the MPBGO label-100
decision) but its external contract is unchanged (still `!tripled && s.mpb.active`).

- **RANSTART** (`mpb_mpotpr` + the LRANST draw in `mpb_apply!`):
  * `mpb_mpotpr(pbalpp,relden,reldsp_lp,a45dbh,cntlp,istdt,icyc,iy)` — mpotpr.f logistic
    `1/(1+expf(9.583−0.08967·min(PBALPP,0.8)·RELDEN))` + the 5 min-condition gates.
    **PROTBK dump-replay BIT-EXACT 5/5** vs FVSie_lpmpb (lp_ran.key, DBGOT dump). Uses a
    glibc `ccall(:expf,...)` for the Float32 transcendental (matches gfortran REAL EXP;
    Julia native `exp(::Float32)` also matched 5/5 here, no ULP straddle on this stand).
  * Draw ORDER (mpbgo.f 100→129): MPBER→MPOTPR(read stand stats)→**one MPRANN draw/cycle**
    whenever NOERR holds and the outbreak is NOT user-scheduled (NTODO=0, MPBYR=0 in COL).
    **MPRANN sequence BIT-EXACT 5/5** (reused `mpb_rand!`, never FFI'd). RNG advances every
    eligible cycle (even when PROTBK=0 / no fire) to stay in sync with the oracle.
  * LEPI: PROTBK=EPIPRB (un-scaled); else PROTBK=MPOTPR·PRBSCL. RELDEN=`s.plot.relative_density`
    (stand CCF); RELDSP(IDXLP) via new `mpb_lp_ccf` (variant per-tree CCF over LP records only).
  * END-TO-END (FVSjl run_keyfile, IE lodgepole): `lp_ran` (unscaled) → **no fire →
    BYTE-IDENTICAL to OFF** (matches oracle RAN==OFF). `lp_ranfire` (PRBSCALE 100) →
    **fires cyc1 → BYTE-IDENTICAL to MPBSTART** in FVSjl; vs oracle cyc1 mort **152/151**
    (±1 = the IE #206 OLDRN straddle, same cornered verdict the deterministic path carries).
    The stochastic inclusion decision (fire-cycle) matches the oracle exactly in both cases.
- **CURRMORT / INVMORT** ICYC=1 GREINF branch (colmod.f:93-99), added to `mpb_colmod`
  (`icyc/lcurmr/linvmr/currmr/greinf` kwargs): DEAD(1,i)=GREINF+CURRMR, GREEN(1,i)=START−GREINF,
  GREINF bumped to START·ZINMOR−CURRMR then capped at START.
  * **GREEN/PRKILL dump-replay BIT-EXACT 12/12** vs FVSie_lpmpb (lp_curr.key, CURRMORT
    classes 3-7 = 2,5,10,3,1 → colmrt DBGCM dump). END-TO-END cyc1 mort **171/170** (±1 #206).
  * GREINF stays **0 on loadable stands** (no treelist damage codes) — so INVMORT with no
    inventory damage is **BYTE-IDENTICAL to MPBSTART** (proven: oracle lp_inv==lp_on; unit
    test g0==default). Faithful & validated; the live-inventory GREINF only differs when
    damage codes exist, which FVSjl cannot yet carry (see deferred #2).

## WHAT IS DEFERRED (dependency-ordered NEXT)
1. **LPOPDY population-dynamics path** (POPDYN keyword → MPBDRV/MPBMOD, ~1050 lines
   incl. GARBEL classifier, SURFCE, phloem model, beetle brood dynamics). Large,
   DEFERRED; mpb_apply! early-returns when `lpopdy` (byte-identical to not wiring).
   This is the last major LPMPB payload.
2. **Live-inventory GREINF damage-code plumbing** (mpbdam.f + mpsdlp.f): MPBDAM scans
   each tree's DBH damage/severity code pairs (code 2 = MPB, severity 3 = successful
   attack) → IPT list; MPSDLP bins IPT into GREINF (live infested, PROB) or CURRMR (dead,
   PROB/GROSPC/(FINT/FINTM)) by size class. FVSjl's treelist carries **no** damage/severity
   codes, so GREINF/CURRMR-from-inventory stay 0 (matches oracle on loadable stands — the
   validated case). To exercise live INVMORT, add damage/severity columns to the treelist
   reader + port MPBDAM/MPSDLP binning. The GREINF **math** is already ported & unit-tested
   synthetically; only the treelist column feed is missing.
3. **The parser's trailing-record sub-keywords** (CURRMORT/QVALUES/INITMORT read a
   SECOND `(3F10.0,T1,3A10)` record for classes 8-10; PSFOUND/DCFOUND/AGGPHERM/
   REPPHERM/ACTSRF/PSDBHLIM/PSPKILL/DCPKILL read `IFIX(ARRAY(1))` extra 8F10.0
   rows). FVSjl's KeywordReader yields one record per keyword — the classes-8-10
   overrides are NOT consumed (kw_mpbin! reads classes 1-7). The stray numeric record is
   harmlessly skipped (starts with a digit ⇒ not treated as a keyword; verified no crash).
   For the COL path only classes 1-7 matter. LOW priority; extend the reader to pull the
   continuation record if a keyfile sets classes 8-10.

## SCOPE FACTS (for the next agent)
- keywds.f main keyword = **'MPB'** (in the extension keyword list, dispatched to
  MPBIN). 40 sub-keywords (see kw_mpbin!). Activation = MANSTART/MPBSTART
  (OPNEW activity **555** at cycle/year) — deterministic; RANSTART = stochastic;
  CURRMORT auto-schedules a cycle-1 outbreak.
- Variants linking real lpmpb: BM CI CR EC EM IE SO TT UT (+ generic mpblkd.f).
  Every other variant links base/exmpb.f (LPMPB inert). IDXLP=7 except CR=11.
- Seam site: gradd.f:63 `IF (IPMODI.EQ.1 .AND. LMPBGO) CALL MPBCUP` (LMPBGO set by
  MPBGO in grincr). MPBCUP → LPOPDY ? MPBDRV : COLDRV.
- WK2 mortality is the DFB-style max-combine (`WK2=max(WK2_background, XT)`, capped
  PROB−1e-6) — mpb_colmrt! mirrors dfb_mrt!.
