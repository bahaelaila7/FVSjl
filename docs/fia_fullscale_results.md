# Full-population sweep results — dig_class breakdown (Pillar-2 at max scale)

Snapshot of `data/fia_sweep.db` (read from a filesystem copy — no contention with the live sweep). Each stand
projected the full default horizon through live FVS + FVSjl, all 10 `.sum` cols per cycle, bucketed into one
dig_class. SN/NE/CS are FULLY swept; LS ~51% at snapshot time.

| Variant | total | bit_exact | ulp_class (cornered) | needs_dig | live_crash | bit-exact-or-cornered |
|---------|-------|-----------|----------------------|-----------|------------|-----------------------|
| SN | 633,628 | 360,162 (56.8%) | 273,466 (43.2%) | **0** | 0 | **100.000%** |
| NE | 178,148 | 140,088 (78.6%) | 38,031 (21.3%) | 19 | 10 | 99.984% |
| CS | 255,951 | 228,323 (89.2%) | 27,541 (10.8%) | 58 | 29 | 99.966% |
| LS (~51%) | 204,770 | 143,920 (70.3%) | 60,311 (29.5%) | 527 | 12 | 99.737% |

`bit_exact` = FVSjl == live on all 10 cols every cycle. `ulp_class` = diverges but cornered by the escalation
guard to a named primitive (DGSCOR/RDPSRT/volume-ULP). `live_crash` = FVS itself crashes (the D38 r9clark /
essprt-family FVS bugs — FVSjl runs clean; cornered as FVS-bugs in FVS_SOURCE_BUGS.md). `needs_dig` = the guard
could NOT auto-corner ⇒ the dig backlog.

## The needs_dig backlog (Pillar-4 open item — corrects an earlier overstatement)
604 stands (SN 0 / NE 19 / CS 58 / LS 527), by stored signature:
- **~601 `structure_densephase`** (NE 17, CS 58, LS 526) — the dense-stand class. This session established
  (43dn/43ds, docs/fia_divergence_taxonomy.md Class B) that structure_densephase on dense stands is the
  compounded-ULP self-thinning RDPSRT tie-break primitive (BA bit-exact, TPA-only diverges), verified on samples.
  So the backlog is cornered-BY-CLASS. **BUT** this session ALSO proved structure_densephase can HIDE a real
  growth bug (the LS aspen HCOR gap, FIX #7) — so blanket-cornering is not safe; each needs per-stand
  confirmation that BA is bit-exact (growth right) with only TPA diverging (self-thinning). The HCOR fix would
  reclassify the aspen subset on re-evaluation, but these records are from the pre-fix pass.
- 1 `volume_persistent` (NE) = stand 207147469020004, the NVEL VOLINIT extreme-height FVS-bug (already cornered, 43do).
- 3 `threshold_crossing` (NE 1, LS 2) — not yet triaged.

HONEST STATUS: full scale is **99.7–100% bit-exact-or-cornered**; the ~604 needs_dig (0–0.26% per variant) are
flagged-for-dig, mostly the self-thinning primitive class (both-sides-traced) but NOT each individually
reclassified post-HCOR-fix. Pillar-4's "no unexplained divergence" holds for the CORNERED classes and processed
batches; the needs_dig backlog is the remaining per-stand-verification work (dig at DIGCAP + a post-fix
re-sweep of the needs_dig CNs to reclassify the aspen subset). This is the campaign's genuine remaining frontier.
