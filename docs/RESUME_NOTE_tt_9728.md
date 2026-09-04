# RESUME NOTE — TT COMPLETE at 9728 (2026-09-04)

## State
- Branch `wt/western-sweep`, engine synced to master **ab67d239** (`git diff master HEAD -- src/` EMPTY).
- **TT cursor = 9728 = population → TT COMPLETE.** Sweep DB `data/fia_sweep_west.db` holds all 9728 TT rows.
- Forward batch 8025→9728 (1703 stands): **bit_exact 915/1703**, 788 diverging, live_crash 0.
  Filter escalated 25 dig-worthy → all root-caused + **CORNERED** this session.
- DIG_PAUSED flag is CLEARED. Dig-queue reset to header-only (archived under
  `data/provisional_archive/fia_dig_queue_west_TT_sweep9728_cornered_*.csv`).

## Verdict for the 25 escalated (docs/fia_dig_session_tt_9728_smhtgf_countstraddle_2026-09-04.md)
All 25 are `structure_densephase`, worst_col CCF/TPA (never BA/SDI), 15.1–87.5%. Ecoregions:
11 in M331D-prefix (already cornered, surfaced only by the escalation guard) + 7 M331Ja + 7 M332Er.
**CORNERED**: same tt small-tree growth ZRAND draw-order straddle as the M331D corner (#206 class),
extended to M331Ja + M332Er — an ecoregion-independent RNG-stream primitive that seeds a
near-SDImax self-thinning count-straddle (BA/SDI carrying capacity preserved/reconverges; TPA/QMD
trade; direction balanced across the cluster). 5 stands per-cycle A/B traced (both new ecoregions;
seedling/dense/mature morphologies). NO engine fix — high-regression RNG rework, no net improvement.

## For the fix-coordinator (master sync)
Corner recorded on the sweep branch ONLY. `filter_digworthy.jl` reads
`/workspace/FVSjl/docs/fia_cornered_stands.txt` (main repo) — sync the 25 CNs there (as the M331D 32
were) so a re-sweep drops them. Cluster rows M331Ja/M332Er also added to `docs/fia_cornered_clusters.tsv`.

## Next variant in the sweep order
Order = NC TT IE EM BM UT CI EC WC PN SO CA WS CR AK. NC + TT now at population; **IE is next**
(cursor check: `cat test/harness/fia/expand/ie_west.cursor`). Relaunch the watchdog to continue
(see test/harness/fia/SWEEP_DURABILITY.md LAUNCH line) — only if directed to sweep beyond TT.
