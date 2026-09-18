# ⚑ SUPERSEDED 2026-09-18 — IE re-sweep DRIVEN TO THE CORNER FLOOR (master `9aa565c4`)

The 2026-09-04 paused state below is HISTORICAL. The one-directional establishment/mortality
under-build it named was closed over the subsequent merges, and the post-restart re-sweep
(N=400, master `9aa565c4`) resolved the 8 remaining `needs_dig` via 4 faithful IE-gated merges
+ corners (all gate 339/11):
- `d9db7c78` broken-top FW2 volume · `2e0db16f` Colville forest-keyed VOLEQ ·
  `cb14152e` AUTOES BAA[1,400] clamp · `b927f4a3` per-point AUTOES species-selection
  (the +117-TPA mono-redcedar over-retention; population sign-tally 8-over/0-under → 0/0).
- #206 straddles (3076266010690, 474165204489998) + establishment realization outlier
  (1856050912290487) CORNERED (seed-proven / sibling-straddle, mean-unbiased).
Re-sweep bit_exact 72→82; residual = dense-phase ULP + establishment realization + cornered
straddles + a deferred region-1 volume tail = the CORNER FLOOR (regimes plant/thinbba/simfire/
salvage show the same signature, no new cluster). Full `Pkg.test` on main-at-master:
**55329 pass / 25 fail / 75 broken, ZERO IE regressions** — all 25 fails pre-existing + IE-
independent (22 ON per-tree-volume ULP, 1 SO-WRD corner, 1 LPMPB oracle-dump, 1 CI simfire
306-vs-307 stale golden). Current state in `RESTART_CHECKPOINT.md` + memory `fvsjl-ie-none-capfix.md`.
Deferred (not IE-none blockers): region-1 full per-forest VOLEQDEF port, wider western
tripling-cycle mistletoe sweep, ON per-tree-volume ULP, CI simfire stale golden.

────────────────────────────────────────────────────────────────────────────

# RESUME NOTE — IE honest re-sweep, PAUSED on a REAL one-directional bug (2026-09-04)

**State:** ie_west.cursor = 2000 / 17808. DIG_PAUSED set (IE real bug). Engine = master
(`git diff master HEAD -- src/` EMPTY; D2 fix 37c66546 IN). Oracle /workspace/.iework/FVSie_g16.
DB: 587 IE rows (0 bit_exact, 580 needs_dig, 7 ulp). Registry UNTOUCHED at 722 (no ie-estab corners).

**Result (cursor 0→2000, FIXED harness):** the reopened establishment cluster ESCALATES (580 needs_dig, not
auto-cornered). Seed test (12 stands): 9 real-determ / 1 straddle / 2 match. **POPULATION signature (40 unbiased):
BA one-directional UNDER 23:1 (mean −5.5), TPA 10over/15under (seed-invariant ±hundreds tail), TopHt 12:2 under
(D2 largely fixed), QMD converges.** Five truly-identical habitat-221 plots give IDENTICAL residual (BA −8) =
determinism proof. ⇒ **REAL one-directional deterministic bug, NOT a straddle. NOT cornerable. PAUSED.**
The FALSE 3360-stand #206 corner is NOT re-created.

**Named real bugs (coordinator to fix on a fix branch, then re-sweep IE from 0):**
1. Mode-B natural-mortality count-partition + density under-build (BA 23:1 under; TPA count NON-conserved by
   ±hundreds seed-invariant ⇒ a mortality-RATE/partition bug, NOT an RDPSRT tie-break which conserves count).
2. Residual established/small-tree height+BA under-build (TopHt 12:2 under; habitat-221 identical-plot BA −8).
3. bitter-cherry dense OH-monoculture (22404926010497): jl over-thins (TPA −574) + under-builds BA (13 vs 75).

Full measurement + fix plan: docs/fia_dig_session_ie_honest_resweep_2026-09-04.md.

**Do NOT:** raise DIGCAP (100) · corner this cluster (one-directional, cornering would MASK the bug) ·
advance the IE cursor past 2000 · continue the sweep · clear DIG_PAUSED until the Mode-B/estab bug is fixed.

**EM:** the 203 pre-existing EM dig-queue rows were moved to .sweep_work/em_digqueue_saved_pre_ie_resweep.csv;
EM is separately paused on the SAME family of real bug (EM_RESUME_NOTE.md, ADJUDICATION_EM_IE_ESTAB doc).

**Cluster order after IE:** EM BM UT CI EC WC PN SO CA WS CR AK (NC done, TT done).
