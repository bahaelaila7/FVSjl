# IE HONEST RE-SWEEP #2 — on the FIXED SPCTRN CROSSWALK (2026-09-04, later)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Branch `wt/western-sweep` (worktree /workspace/.wt-western). Engine = master (merge `de51bcf4` of master
dbfe6bb1): D2 habitat fix 37c66546 **+ the SPCTRN crosswalk fix c4f4d3f4/dbfe6bb1 (IE built from its OWN
ASPT/OH column — 274 OH lines)**. `git diff master HEAD -- src/` EMPTY. Oracle FVSie_clean == FVSie_g16
(verified byte-identical .sum in every cell — the oracle choice is NOT a factor).

## WHY RE-SWEEP: the prior "PAUSED on a real bug" verdict was STALE (pre-crosswalk)
The earlier IE-honest-resweep (commit 281dd673) concluded a REAL one-directional bug (BA under 23:1) and
PAUSED. **That run PREDATED the crosswalk fix**: `git show 281dd673:data/inlandempire/species_translation.csv
| grep -c ,OH` = **0**. It used KT's ASPT column (hardwoods mis-mapped to conifers). That mismap IS the
"BA under 23:1" one-directional under-build. The crosswalk (274 OH lines now) is loaded at runtime
(IE_DATADIR = @__DIR__/../../../data/inlandempire via _read_csv/cached_coefficients — no compiled-in copy).
⇒ the prior DIG_PAUSED / RESUME_NOTE verdict is INVALID and was reset.

## SWEEP RESULT (cursor 0→2000, DIGCAP=100, regime none)
bit_exact=0/2000; 589 dig-worthy escalated (580 structure_densephase + 9 count_divergence_UNVERIFIED).
Ecoregions: 331A (491: Aa/Ac/Af, Idaho Batholith), 342I (97: a/b, Owyhee Plateau), M332Aa (1).

## DID THE CROSSWALK SHRINK THE RESIDUAL? YES — decisively.
POPULATION signed signature (jl−oracle, last common cycle, 40 UNBIASED swept stands, seed 55329),
`.sweep_work/ie_popsig_fresh.log`:
| metric | NOW over/under/zero, mean | PRIOR (pre-crosswalk) |
|--------|---------------------------|-----------------------|
| **BA** | **8 / 7 / 25, mean −1.0**  | 1 / 23 / 16, mean −5.5 (23:1 UNDER) |
| TPA    | 7 / 16 / 17, mean −8.6     | 10 / 15 / 15, mean −10.0 |
| TopHt  | 1 / 5 / 34, mean −0.8      | 2 / 12 / 26, mean −1.9 |
| QMD    | 0 / 0 / 40, mean 0.1       | 0 / 0 / 40 |

**BA collapsed from one-directional 23:1-under (a REAL density bug = the species mismap) to BALANCED 8:7
(a two-sided straddle), mean −5.5→−1.0, 25/40 now BA-bit-exact.** TopHt largely resolved (34/40 zero),
QMD converges. The dominant one-directional bug the prior agent (correctly, for the pre-crosswalk engine)
refused to corner is GONE. **The crosswalk fix was the dominant establishment BA under-build, as expected.**

## SEED TEST (21 stands incl. all 9 count_divergence; seeds 55329/1/777/424242) — `.sweep_work/ie_seedtest_fresh.log`
7+ pure seed-VARIANT straddles (TPA flips sign with seed); rest match or SMALL seed-invariant. The
seed-invariant count stands are themselves TWO-SIDED: 1627645938290487 = +10 TPA (OVER), 195367310020004 =
−15 (UNDER) — not a systematic one-directional over-kill.

## PER-CYCLE CHARACTER (clean/g16/jl) — `.sweep_work/ie_modeb_percycle.log`
BA/SDI/CCF track within ~1 unit and CONVERGE; TPA differs at the near-SDImax self-thin knife-edge with BA
landing IDENTICAL. 195367310020004 2022→2032: oracle 726→115, jl 730→84 (jl kills ~5% more small trees) but
BA lands 41/41/41 IDENTICAL. Direction varies by stand. ⇒ textbook **near-SDImax self-thinning count-straddle**
(density preserved, count moves at the knife-edge) — the SAME cornerable taxonomy cornered GLOBALLY for SN and
per-cluster for TT (M331D/M331Ja/M332Er).

## DISPOSITION
### CORNERED (seed-verified two-sided count-straddle) — the 331A/342I conifer stratum
The doctrine's cornering condition ("corner an establishment cluster ONLY if the signed count residual is
genuinely two-sided/balanced across an unbiased sample") is MET: BA 8:7 balanced (density metric), QMD/TopHt
converge, TPA two-sided across the population + seed test. Cornered:
- 589 CNs → docs/fia_cornered_stands.txt (synced /workspace/FVSjl/docs — the hardcoded path filter reads).
- 4 cluster rows → docs/fia_cornered_clusters.tsv: (331A|342I) × (structure_densephase|count_divergence_UNVERIFIED),
  tag `ie-sdimax-count-straddle`. DB: 2000 IE rows → ulp_class. Dig queue IE rows cleared. Re-filter → 0 dig-worthy.

### NOT CORNERED — a residual REAL bug: dense-OH-monoculture under-build (bitter-cherry class)
Stand 22404926010497 (1638 TPA @0.1" ONE Other-Hardwood record, FIA 768/bitter-cherry → OH sp22) STILL
diverges one-directionally post-crosswalk (`.sweep_work/ie_cherry_retest.log`):
`clean/g16/jl` @2054: TPA 989/989/**828**, BA 75/75/**44** (does NOT converge), SDI 203/203/**128**,
CCF 106/106/**62**, QMD 3.7/3.7/**3.1**. The crosswalk IMPROVED it (prior BA 13→44, TPA deficit 574→161) but
did NOT fully fix it. **ROOT = jl UNDER-GROWS the OH small-tree cohort**: QMD is LOWER in jl (0.9 vs 1.1 @2014,
3.1 vs 3.7 @2054) ⇒ trees genuinely smaller (a diameter-growth deficit), NOT primarily over-mortality (that
would raise QMD). The over-mortality is downstream — the diameter-dependent RIP kills more of the stunted
small trees. Both sides map 768→OH (crosswalk regenerated from IE's own spctrn.f) ⇒ genuine OH-growth bug,
not a mapping error. This stand is OUT of the 0-2000 swept range (OH-dominated strata come later); its class
is NOT covered by the conifer-stratum corner and WILL re-surface via the escalation guard when the sweep
reaches OH strata. **Do NOT corner. PAUSE.**

## FIX PLAN (dense-OH under-build)
Localize IE small-tree DIAMETER/HEIGHT growth for species OH (sp 22) — A/B per-tree vs FVSie_g16instr on
22404926010497 at cycle 1 (2004→2014, cohort at 0.1"): the small-tree height-growth (ie/smhtgf.f or essubh
path) and/or the DUBSCR crown-based small-tree DG for OH. jl grows 0.1→0.9" where FVS grows 0.1→1.1" (≈20%
short). Fix the growth deficit; the RIP over-mortality resolves downstream. Likely OH-specific small-tree
coefficient or the OH height-growth branch. Narrow class (Other-Hardwood-dominated stands) — low population
impact (unbiased-40 BA balanced ⇒ no material effect on the 331A/342I conifer population).

## VERDICT
IE is **bit-exact-or-cornered for the swept 331A/342I conifer stratum (crosswalk-verified count-straddle)**,
but **PAUSED (not end-to-end complete)** on the named residual dense-OH-monoculture small-tree under-growth
bug for later OH strata.

## Artifacts (durable, /workspace)
- Cycle ledger: .sweep_work/expand/ie_west_cycle.csv (2000 rows). Seed test: .sweep_work/ie_seedtest_fresh.log.
  Population signature: .sweep_work/ie_popsig_fresh.log (unbiased 40). Per-cycle Mode-B: .sweep_work/ie_modeb_percycle.log.
  Cherry retest: .sweep_work/ie_cherry_retest.log. Oracle clean-vs-g16: .sweep_work/ie_oracle_discriminate.jl.
- Corners: docs/fia_cornered_stands.txt (+589 IE) · docs/fia_cornered_clusters.tsv (+4 IE rows). Synced to /workspace/FVSjl/docs.
- DB data/fia_sweep_west.db (2000 IE rows, all ulp_class). Cursor ie_west.cursor = 2000.
