# WSBWE apply / wire / commit recipe (from a REAL FVSjl checkout)

This agent ran isolated in an FVS-Fortran worktree and CANNOT commit to the FVSjl
checkout (git-isolation guard). Everything is staged in `scratchpad/wsbwe/`. The
main loop applies + wires + commits from a real `kt-variant-port` checkout.

Beachhead = RNG (bit-exact) + `WSBW … END` keyword reader (INERT). No engine seam
(the defoliation effect path is deferred — see "NEXT CHUNKS"). This mirrors the
WWPB/WPBR beachhead exactly, but WSBWE's oracle RUNS (stand-level), so the effect
seam is later dump-replay-validatable (unlike WWPB).

## APPLY
1. `cp scratchpad/wsbwe/wsbwe.jl        src/engine/wsbwe.jl`
2. `cp scratchpad/wsbwe/test_wsbwe.jl   test/unit/test_wsbwe.jl`

## WIRE — 5 mechanical seams (mirror LPMPB's mpb wiring)

### (1) src/core/state.jl
- After `abstract type AbstractMpbState end` (~line 946) add the forward decl:
  ```julia
  # Forward declaration for the Western Spruce Budworm (WSBWE) defoliation model
  # state; concrete `WsbweState` is in engine/wsbwe.jl (included after this file).
  # `nothing` until a WSBW keyword block activates it — inert for every non-WSBWE
  # run. Stand-level; the defoliation effect seam is deferred, so no engine seam.
  abstract type AbstractWsbweState end
  ```
- In `mutable struct StandState`, immediately AFTER the `mpb::Union{...}` field
  (line 971) add:
  ```julia
      wsbwe::Union{AbstractWsbweState,Nothing}           # Western Spruce Budworm (WSBWE) defoliation model; nothing until WSBW (stand-level; inert reader, effect seam deferred)
  ```
- In `function StandState(variant; faithful)` the constructor call (line 985-989)
  currently ends the field list with **10** trailing `nothing`s (fire … mpb). Add
  ONE more `nothing` (the 11th, for `wsbwe`) — e.g. change the tail to
  `..., nothing, nothing,\n    )` so there are 11 nothings after `DbsState()`.

### (2) src/FVSjl.jl
- After `include("engine/lpmpb.jl")` (line 336) add:
  ```julia
  include("engine/wsbwe.jl")                # Western Spruce Budworm (WSBWE) — beachhead: BWERAN RNG (bit-exact) + WSBW keyword reader (inert; effect seam deferred)
  ```

### (3) src/engine/keyword_dispatch.jl
- After the `elseif kw == "MPB"; kw_mpbin!(...)` line (line 2446) add:
  ```julia
          elseif kw == "WSBW";     kw_wsbwe!(s, rec, kr)     # Western Spruce Budworm block (keywds.f opt 8 'WSBW'; DEFOL/GENDEFOL/OBSCHED/DAMAGE/… → s.wsbwe; INERT — effect seam deferred)
  ```

### (4) src/engine/simulate.jl
- The WSBWE seam is now WIRED (but INERT — apply! early-returns; see below). Replace
  the existing WSBWE doc-note block (the paragraph after the `if !tripled && s.mpb …`
  block, ~line 674-679, ending "...so no per-cycle apply is wired.") with the gated
  apply, mirroring the mpb seam immediately above it:
  ```julia
      # WSBWE (Western Spruce Budworm, wsbwe/*.f): stand-level DEFOL defoliation
      # (BWEGO→BWECUP→BWEDR→BWEDAM/BWEPDM, grincr.f:414 / gradd.f:108). The reader is
      # ported; the deterministic payload kernels (wsbwe_rdds/wsbwe_rhtg/wsbwe_mort_pr)
      # are dump-replay validated bit-exact. INERT: wsbwe_apply! early-returns until the
      # BWESIT→BWEAGE→BWEDAM foliage/PRBIO feeder + BWEPDM per-tree apply are ported
      # (scratchpad/wsbwe/HANDOFF.md "NEXT CHUNKS"), so this projects byte-identically.
      if !tripled && s.wsbwe !== nothing && (s.wsbwe::WsbweState).active &&
         wsbwe_go(s.wsbwe::WsbweState)
          wsbwe_apply!(s, old_tpa, fint)
      end
  ```
  NOTE the gate `wsbwe_go` fires only on the manual-DEFOL branch with a scheduled
  activity; a WSBW block with only GENDEFOL (`lbudl`) or no DEFOL never enters apply!.
  Because apply! is currently a no-op, the seam is byte-identical even when it fires
  (proven: live oracle off ≡ stock, and DEFOL-on-non-host ≡ off, both byte-identical).

### (5) test/runtests.jl
- After `include("unit/test_lpmpb.jl")` (line 27) add:
  ```julia
      include("unit/test_wsbwe.jl")          # WSBWE beachhead: BWERAN RNG seed 55329 bit-exact vs pristine wsbwe/bweran.f + WSBW keyword reader + INERT seam
  ```

## VERIFY (from the real checkout)
```
JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project -e 'using Pkg; Pkg.test()'
```
The test_wsbwe.jl set must pass (RNG 8-draw golden, reseed, defaults, reader,
INERT A/B). Because the seam is inert, the full suite must remain green and every
existing .sum byte-identical (WSBWE touches nothing without a WSBW block, and is
inert even with one).

## ORACLE (regenerate + instrument)
```
# edit build_em_wsbwe.sh WORK= to YOUR worktree scratchpad/wsbwe first!
bash build_em_wsbwe.sh                                  # → FVSem_wsbwe (+ emobj/ cache)
```
- **Run convention**: the binary reads the keyfile name from stdin (unit 15):
  `echo host_defol.key | ./FVSem_wsbwe` (exit 10 = normal FVS stop; .sum written).
- **Faithfulness**: off-run `.sum` == stock `bin/FVSem` off-run `.sum` byte-identical
  except the line-1 timestamp (compare `tail -n +2`). GOLDENS.md.
- **Instrumented dump-replay** (single-.o swap): put an instrumented `bwe*.f` in
  `ovr/` and pass its ABSOLUTE path + an ABSOLUTE out path:
  ```
  W=$PWD; bash build_em_wsbwe.sh "$W/ovr" "$W/FVSem_wsbwe_ins"
  ```
  TWO GOTCHAS (both bit me): the override loop uses `"$OVR"/*.f` AFTER the script has
  `cd`'d into $BWE, so a RELATIVE ovr dir silently matches nothing (no "OVERRIDE"
  line printed — watch for it). And the link runs from inside `$OBJ`, so a RELATIVE
  out name lands the binary in `emobj/` — pass an absolute out path.
- **ALWAYS A/B the instrumented `.sum` against the clean relink `.sum`** (must be
  byte-identical modulo timestamp) BEFORE trusting a dump. Dump at Float32-hex via
  `TRANSFER(x,IHX)` + `Z8.8` to a private unit (78/79) opened APPEND to a side file —
  NEVER to JOSTND/the .sum. `ovr/bwedam.f` + `ovr/bwepdm.f` are the working exemplars.

## FIXTURE (host stand that fires the effect)
- `host.tre` (DF/S/AF) + `host_defol.key` / `host_off.key` (OPEN-3 treedata recipe).
  Regenerate host.tre from emt01.tre relabelling species cols 34-36 to DF/S /AF.
  `host_defol` DEFOL 85/80/70/60% on DF/ES/AF for 1990-1999 ⇒ TPA 536→52 by 2000.
- `nonhost.tre` / `nonhost_*.key` (all-LP): the non-host INERT proof (WSBW-on ≡ off).

## DONE THIS CHUNK (effect kernel — dump-replay validated bit-exact)

- **Host stand + oracle-fires CONFIRMED.** `host.tre` = emt01 relabelled to DF/S(ES)/
  AF (the EM budworm hosts — note EM has NO WF/GF), read via `OPEN 3 / host.tre /
  TREEDATA 3` (the plain `TREEDATA\n<file>` form yields a 0-TPA degenerate stand —
  that was the prior "no effect" red herring). `host_defol.key` = sustained heavy
  DEFOL (85/80/70/60%) on DF/ES/AF for years 1990-1999. **Effect is dramatic**:
  TPA 536→52 and BA 77→8 by 2000 (vs off 528/102), sustained through 2090. GOLDENS.md.
- **Two deterministic payload kernels ported + VALIDATED bit-exact** (Float32-hex
  dump-replay, instrumented `.sum` byte-identical to the clean relink first):
  - `wsbwe_rdds` / `wsbwe_rhtg` — BWEDAM Nichols growth-loss (bwedam.f:184/188):
    **90/90 bit-exact** (bwedam.f single-.o swap → `bwedam_dump.txt`).
  - `wsbwe_mort_pr` — BWEPDM Marsden mortality logistic (bwepdm.f:640): **81/81
    bit-exact** (bwepdm.f single-.o swap → `bwepdm_dump.txt`).
  gfortran `**`/EXP reproduce bit-exactly in Julia Float32 — no ULP drift.
- **Seam wired + inert-proven.** `wsbwe_go` (BWEGO manual-DEFOL gate) + `wsbwe_apply!`
  (BWECUP, early-returns). Live-oracle inert proof: off ≡ stock byte-identical AND
  DEFOL-on-**non-host** stand ≡ off byte-identical (`nonhost_*.key`).

## NEXT CHUNKS (dependency-ordered) — complete the DEFOL apply, then BUDLITE

The two payload kernels are done; what remains to make `wsbwe_apply!` NON-inert is the
FEEDER chain that produces their inputs (FOLPOT/PRBIO → AVPRBO/MFT/MFM/PEDDS/PEHTG),
then the per-tree application. Instrument as before (single-.o swap, hex TRANSFER,
verify instrumented `.sum` == clean relink) on `host_defol.key`.

1. **Foliage setup — BWESIT + BWEBMS + BWECRC + BWEADV** (bwesit.f, bwebmsem.f,
   bwecrc.f [ported trivially: <23ft→sz1, <46→sz2, else sz3; ICRC=3·sz-2..3·sz],
   bwesin.f [ported trivially], bweslp.f [piecewise-linear interp]). Builds FOLPOT/
   FOLADJ/BWTPHA/FOLNH per host×crown×age from the FVS tree list. BWEBMS (Moeur
   foliage biomass, ICVOPT=2, EM IBIOMP map) is the one nontrivial numeric piece —
   dump-replay WK4 (per-tree biomass) then FOLPOT. IBWSPM(EM)=[7,6,2,7,7,7,7,5,4,7,…].
2. **Aging + PRBIO — BWEAGE** (bweage.f): ages foliage, computes PRBIO (retained-
   biomass proportion) and TOTR. Needs RELFX/RELFY/THEOFL/PRCRN3 (all in bwebkem.f,
   already read). Dump-replay PRBIO + TOTR.
3. **Manual defoliation apply — BWEDEF** (bwedef.f): FNEW/FOLD1/FOLD2/FREM ×
   (1−PARMS/100) per scheduled DEFOL. Deterministic. This is what makes PRBIO drop.
   The DEFOL activity is scheduled via OPNEW act 2151 (`defol_sched` already captured
   in the reader; you need OPGET-style per-year matching in the BWECUP year loop).
4. **BWEDAM accumulation** (bwedam.f full): CUMDEF/APRBYR/AVPRBO/CDEF/AVYRMX +
   PEDDS/PEHTG accumulate (PEDDS+=RDDS/BWFINT). Feeds the two validated kernels.
5. **BWEPDM per-tree application** (bwepdm.f full): topkill (BWERNP AVDEF gate +
   Marsden TK0..TK8 logistic + HTGSTP-style truncation), DG reduction
   (`DDS*BWERNP(PEDDS,.03)`), HTG reduction, and WK2 mortality
   (`max(BASE, PROB·(1−survival))`). **DRAWS THE DAMAGE RNG** — port BWERNP (beta
   variate, calls BWEBET) + BWEBET (Cheng 1978 rejection sampler, 2 BWERAN draws/pass)
   in EXACT draw order, seeded DSEEDD (55329, saved/restored via BWERPT/BWERGT around
   the tree loop). The RNG stream is already bit-exact (`wsbwe_rand!`); the challenge
   is the rejection-loop draw ORDER + Float32 ALOG/EXP in BWEBET matching gfortran.
6. **BUDLITE population model + weather** (GENDEFOL) — STOCHASTIC, needs a weather
   source (WEATHER / block-data), NEMULT/PARASITE/OBSCHED. Largest; last. Its own
   RNG consumers are BWEBET + BWERNP (same as step 5) plus BWELIT→BWEDIE larvae
   survival. `wsbwe_apply!` already early-returns on the `lbudl` branch.

Expected end-to-end .sum verdict once steps 1-5 are wired: bit-exact-or-CORNERED by
the EM #206 growth straddle (the defoliation math is exact on equal inputs; absolute
mortality scales with the straddling cycle growth — same verdict family as DFB/LPMPB).

# =============================================================================
# STEPS 1–5 DONE — apply! is LIVE (feeder + BWEPDM per-tree + BWERNP/BWEBET RNG),
# ALL routines dump-replay BIT-EXACT vs FVSem_wsbwe (host_defol.key, cycle 1).
# =============================================================================

## What is validated (measure-don't-infer; instrumented .sum ≡ clean relink first)
- **BWEBMS** WK4 per-tree foliage biomass ............. 81/81 bit-exact
- **BWESIT** FOLPOT/FOLADJ/BWTPHA/FNEW/FOLD1/FOLD2/FREM/HOSTST ... bit-exact
- **BWEAGE** per-year PRBIO + TOTR ..................... bit-exact (10 yrs)
- **BWEDAM** per-year AVPRBO/PEDDS/PEHTG/AVYRMX/CDEF .... bit-exact (10 yrs)
- **final PEDDS/PEHTG/AVYRMX/PRBIO** (feed BWEPDM) ...... bit-exact
- **BWEPDM** per-tree topkill+DG+HTG+WK2+ICR/ITRUNC/NORMHT/IMC/KTK/PR ... 81/81
- **BWERNP/BWEBET** damage RNG stream .................. 519/519 draws bit-exact
  (draw ORDER exact: per host tree = BWERNP(AVYRMX,.06)[beta draws] → BWERAN(topkill)
   → BWERNP(PEDDS,.03)[beta draws]; XH is the deterministic branch, no draw, since
   PEDDS<0.99). Damage seed DSEEDD=55329 lives in `w.rng_s0`.
Reproduce: `bash build_em_wsbwe.sh "$PWD/ins" "$PWD/FVSem_wsbwe_ins"` then
`echo host_defol.key | ./FVSem_wsbwe_ins` (writes bwe_{sit,age,dam,pdm}.txt +
bwe_rng.txt); then `julia replay_bwepdm.jl` (BWEPDM+RNG), `julia -e
'include("replay_feeder.jl");include("replay_yearloop.jl")'` (feeder), and
`julia test_staged.jl` (validates the STAGED wsbwe.jl functions themselves).

## ⚠ KEY FINDING — glibc transcendentals (load-bearing for ALL FVSjl Float32 math)
gfortran's `ALOG`/`EXP`/`**`(real exp) resolve to **glibc** `logf`/`expf`/`powf`,
which differ from Julia's native/openlibm Float32 math by ~1 ULP on ~20-30% of
inputs (measured: native log 16412/60000 mismatch, exp 3006/60000; glibc-ccall
0/60000). WSBWE therefore routes EVERY transcendental through glibc via
`ccall((:logf,"libm.so.6"),Float32,(Float32,),x)` (`wsbwe_log/exp/pow`), and uses
`wsbwe_powi` (libgfortran pow_r4_i4 right-to-left binary) for `PR**IBWYR` (integer
power — Julia's `^Int` also drifts 1 ULP). The two PRE-EXISTING kernels
(`wsbwe_rdds/rhtg/mort_pr`) were native-Julia and only *lucky-exact* on their tested
rows; this chunk switched them to glibc (still pass their goldens). NB: `:logf`
without the `"libm.so.6"` path binds Julia's openlibm — must name the lib.
NOT an RNG FFI: the BWERAN LCG stays pure Julia (`wsbwe_rand!`).

## APPLY (staged files → checkout)
1. `cp scratchpad/wsbwe/wsbwe.jl      src/engine/wsbwe.jl`    (feeder+apply added; kernels→glibc)
2. `cp scratchpad/wsbwe/test_wsbwe.jl test/unit/test_wsbwe.jl` (new kernel+RNG golden testset)

## WIRE — the 5 seams from the beachhead are ALREADY merged (HEAD 868e95cc). The
gated call `if !tripled && s.wsbwe!==nothing && (s.wsbwe).active && wsbwe_go(...)`
already invokes `wsbwe_apply!(s, old_tpa, fint)` in simulate.jl — apply! is now
non-inert, so NO wiring change is needed; just replace wsbwe.jl + test_wsbwe.jl.

## StandState field mappings used by wsbwe_apply! (verified vs src/core/state.jl)
- trees: `t.species`(ISP) `t.height`(HT) `t.dbh`(DBH) `t.diam_growth`(DG)
  `t.ht_growth`(HTG) `t.crown_pct`(ICR) `t.tpa`(PROB) `t.cuft_vol`(CFV)
  `t.trunc`(ITRUNC) `t.norm_ht`(NORMHT) `t.mort_code`(IMC) `t.plot_id`(ITRE) `t.dmr`(IMIST)
- scalars: `s.plot.elevation`(ELEV) `s.plot.old_tpa`(OLDTPA) `s.plot.old_qmd`(ORMSQD)
  `s.plot.forecast_interval`(IFINT) `s.plot.points_inv`(IPTINV)
  `s.control.cycle`(ICYC) `s.control.cycle_year[cycle]`(IY(ICYC))
- FVS order: `s.control.sp_count_tab`(ISCT) + `s.scratch.idx1`(IND1) — species-major,
  the SAME index dfb/forest_type use; gives the exact BWEPDM draw order.
- bark: `bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, dbh)` = EM BRATIO (H unused, Wykoff).
- WK2 timing (verified simulate.jl:449-636): MORTS runs BEFORE the GRADD seam on the
  NON-tripled stand, so at apply! `PROB=old_tpa[i]`, `WK2=old_tpa[i]-t.tpa[i]`
  (background mort). apply! writes `t.tpa[i]=PROB*(1-PR_final)` only when PR>BASE.
- variant gate: `s.variant isa EasternMontana` (only EM host block-data ported).

# =============================================================================
# 2026-08-17 — apply! is LIVE + END-TO-END VALIDATED. Feeder root-caused & fixed.
# =============================================================================

## WHAT WAS DONE (this chunk — flip LIVE + root-cause the feeder NaN/under-kill)
`WSBWE_APPLY_LIVE=true`. The isolated dump-replay validated the feeder MATH on the
oracle's cycle-1 inputs, but the LIVE FVSjl tree-list inputs differed. Root-caused by
instrumenting BOTH sides (oracle bwebmsem.f WRITE(77) OLDTPA/ORMSQD dump A/B'd vs the
clean relink; FVSjl feeder per-year AVPRBO/PEDDS dumps) — TWO faithful bugs:

1. **OLDTPA/ORMSQD were `s.plot.old_tpa`/`old_qmd`** (the PREVIOUS-cycle stored scalars,
   = 0 at cycle 1) instead of the CURRENT cycle-start density. FVS grincr.f:281/285 sets
   `OLDTPA=TPROB`, `ORMSQD=RMSQD` every cycle via DENSE (dense.f:182/250: `TPROB=ΣP`,
   `RMSQD=√(ΣD²P/TPROB)`), and bwebmsem.f:124/130 reads those. With 0 ⇒ `ALOG(0)=-Inf`
   and `D/ORMSQD=+Inf` ⇒ WK4 foliage biomass = NaN ⇒ PRBIO/AVPRBO/PEDDS NaN ⇒ the DG
   `sqrt((D·BARK)²+DDS)` = NaN (the crash) and 0.15-vs-NaN under-kill. FIX: compute
   `oldtpaS=TPROB`, `ormsqd=RMSQD` from the cycle-start tree list (`old_tpa` arg + dbh)
   in wsbwe_apply!. Bit-exact: FVSjl 589.6528/5.1449676 vs oracle 589.6527/5.144968.

2. **DEFOL species field read as raw numeric.** The keyfile species is an ALPHA code
   ("DF"/"ES"/"AF"); `r.values[2]` on alpha = 0, and `spc<=0` in BWEDEF means ALL hosts —
   so each of the 3 records (DF,ES,AF) defoliated EVERY host, once per record ⇒ the host
   was defoliated 3× ⇒ retained biomass 0.15³=0.003375 instead of 0.15 (44× over-kill →
   AVYRMX pinned to 1.0, PEDDS NaN). FVS bwein.f:519 `CALL SPDECD(2,IS,...)` decodes it
   to the species index. FIX: `spp = species_selector(s, r.fields[2])` (FVSjl's SPDECD)
   ⇒ DF→3, ES→8, AF→9, each defoliating only its own host once. Bit-exact: year-1990
   AVPRBO=0.15, final PEDDS=0.20520946/PEHTG=0.12192549/AVYRMX=0.85 == PDMPE golden.

3. **BWEGO per-cycle gate.** Added `any(scheduled DEFOL year ∈ [IY(ICYC),IY(ICYC)+IFINT-1])
   || return` (bwego.f:88 `OPFIND(1,2151,I); LDEFOL=I>0`) so apply! fires ONLY in the
   DEFOL cycle (cycle 0 here), matching the oracle (which calls BWECUP once). Placed
   before any RNG draw ⇒ post-outbreak cycles stay byte-identical.

No NaN-masking guard was added — every fix corrects an upstream input to the value the
Fortran itself uses.

## END-TO-END .sum-DELTA (host_defol_jl − host_off_jl vs FVSem_wsbwe host_defol − host_off)
Oracle (81 tripled trees) vs FVSjl (27 un-tripled). ΔTPA/ΔBA = defol − off:
```
 YR | ORACLE off  defol  Δ        | JL off   defol  Δ
2000|  528/102    52/8  -476/-94  | 528/102  52/9  -476/-93   ΔTPA EXACT
2010|  520/126    51/12 -469/-114 | 520/127  52/13 -468/-114
2050|  459/206    48/31 -411/-175 | 464/212  48/31 -416/-181
2090|  394/264    45/54 -349/-210 | 380/267  45/52 -335/-215
```
BIT-EXACT-or-CORNERED by the EM #206 straddle (the off-run itself straddles: 2090 off
BA jl267/or264). Damage RNG first draw = 3EDDB57A (unit golden). Cycle-0 draw count =
183 on 27 UN-tripled host trees vs oracle 519 on 81 TRIPLED (~6.5 draws/tree both) — the
WSBWE seam runs PRE-tripling (`if !tripled`, like DFB/DFTM/MPB/WPBR) while FVS BWECUP
runs POST-tripling (gradd.f:111, after grincr.f:537). Per the memory doctrine ("per-record
treelist diff INVALID post-tripling — use .sum") this is the CORNERED residual, not a bug;
the .sum DELTA is the acceptance oracle and passes. Regressions: multicycle 339/11
byte-identical, unit 33/33, non-host WSBW-DEFOL ≡ off byte-identical.

## APPLY (this chunk)
1. `cp scratchpad/wsbwe/wsbwe.jl src/engine/wsbwe.jl`  (4 code changes: LIVE flag=true,
   SPDECD species decode, OLDTPA/ORMSQD=TPROB/RMSQD, BWEGO per-cycle gate; +comments).
   No wiring change — the seam call in simulate.jl already invokes wsbwe_apply!.
   test_wsbwe.jl unchanged (still 33/33).
VERIFY: `JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project -e 'using Pkg; Pkg.test()'`
and the multicycle set must stay 339/11.

## REMAINING (next handoffs, dependency-ordered)
0b. **Tripled-seam draw parity (OPTIONAL).** To hit the exact 519-draw / per-tree bit-exact
   mortality, run wsbwe_apply! on the TRIPLED tree list (move the seam after tripling, à la
   FVS gradd.f:111). Not required — the .sum DELTA already validates cornered; deferred.
1. ~~END-TO-END .sum-DELTA~~ **DONE** (above).
2. **Other host variants** (BC BM CI EC SO TT + generic): add their `bwebk<v>.f`
   IBWSPM/IBIOMP + `bwebms<v>.f` coeffs (only EM `WSBWE_*_EM`/`WSBWE_BINT*` here);
   the arithmetic + RNG are variant-independent. Gate currently `isa EasternMontana`.
3. **BWEPRB multi-cycle POFPOT carryover** (bweprb.f) — not exercised by the single
   1990-outbreak fixture (later cycles schedule no DEFOL ⇒ don't fire); needed only
   for OVERLAPPING multi-cycle outbreaks. DEFERRED.
4. **BUDLITE/GENDEFOL** stochastic population path (step 6) — still early-returns on
   `lbudl`. Its RNG consumers are BWEBET+BWERNP (done) plus BWELIT→BWEDIE; needs a
   weather source. Largest; last.
