# FIA/FVS behaviour-compatibility campaign — status capstone

Single-view status of the four pillars (charter: `FIA_FVS_COMPAT_GOAL.md`). Floor held at suite **38587 pass /
0 fail / 75 broken** (≥ the 38527/143/0 floor). The **off-switch (`touch docs/FIA_FVS_COMPAT_COMPLETE`) is the
USER's call** — this doc records readiness, it does not decide it.

## Pillar 1 — Scale & Stratification ✅ artifact landed
Full-population coverage of all 4 campaign variants: **SN 637,641 / NE 178,149 / CS 255,952 / LS 400,649 =
1,472,391 stands** (~9000× the 162-stand modernization baseline), spanning 95/85/69/74 distinct forest types,
227/117/103/99 ecoregions, all age classes (seedling→old-growth) and site classes (lo→vhi).
Artifact: `docs/fia_pillar1_coverage.md` + reproducible generator `test/harness/fia/pillar1_coverage.py`.
Sweep cursors: SN/NE/CS == population (FULLY swept); LS ~43% and progressing.

## Pillar 2 — Multi-cycle projection compatibility ✅ sample done, full-scale in progress
Every cycle's 10 `.sum` cols bit-exact-or-cornered vs freshly-relinked live FVS, over the full default horizon.
Stratified sample (`.sweep_work/pillar2_sample_results.txt`): **NE 12/12, CS 12/12 bit-exact; SN 7/12, LS 10/12**
bit-exact, all residuals = named ULP/volume primitives (see Pillar 4). Fidelity is stratum-UNIFORM (no forest-type/
site/age/geography model gap). Full-scale: the running sweep validates the entire population per variant; growing
stands are 100% bit-exact-or-cornered (~28% pure-bit-exact + DGSCOR/RDPSRT-cornered remainder). SN/NE/CS complete;
LS in progress → completion at ALL_VARIANTS_EXHAUSTED.

## Pillar 3 — Management-scenario compatibility ✅ done
Real plots under standard silvicultural regimes match live across the projection, bit-exact-or-cornered: 4 variants
× thinning (BA/TPA/DBH — THINBTA/THINDBH/THINBBA) × salvage / plant / regen / SIMFIRE-fire × 3 volume/board paths
(audit slices ~43dc–43dh + the keyword-coverage/regime sweeps). No-op regime rate == the Pillar-2 growth-only rate
(consistency check); fired cuts bit-exact and traced. Surfaced + fixed real bugs along the way (e.g. mid-cycle
SIMFIRE, MORTMSB, fire tripling-order, FFE phasing).

## Pillar 4 — Divergence taxonomy & cornering ✅ artifact landed
Consolidated taxonomy `docs/fia_divergence_taxonomy.md`: every non-bit-exact class both-sides-traced and FIXED or
CORNERED — **7 FVSjl bugs fixed** (floor held), **4 ULP-class named primitives** (RDPSRT self-thinning tie-break,
direct DGSCOR/volume-ULP, non-native cycle drift, COMPRESS eigensolver), **4 FVS bugs** FVSjl is correct on and
doesn't replicate (D38 r9clark SIGFPE, CS essprt SIGFPE, NE VOLINIT extreme-height zeroing, shared SDI overflow).
Full scale (docs/fia_fullscale_results.md): 99.7-100% bit-exact-or-cornered (SN 100% / NE 99.98% / CS 99.97% /
LS 99.74% partial). CORRECTION: a needs_dig BACKLOG remains (SN 0 / NE 19 / CS 58 / LS 527 = 604), ~601 of them
structure_densephase = the self-thinning RDPSRT primitive class (cornered-by-class, verified on samples) but NOT
each individually reclassified post-HCOR-fix — and this session proved that class CAN hide a real growth bug
(aspen HCOR). So Pillar-4 is SUBSTANTIALLY met (cornered classes + processed batches) with the needs_dig backlog
as the genuine remaining per-stand-verification frontier.

## Remaining to full completion
- Let the sweep finish LS (the last ~57%) → ALL_VARIANTS_EXHAUSTED = full-population Pillar-1/2 done at max scale.
- Process any dig batches the sweep raises before then (pause → both-sides-trace/corner → resume), extending the
  Pillar-4 taxonomy as needed.
- The 2 currently-queued LS densephase candidates are pending triage (expected: self-thinning-tie-break primitive
  post-HCOR-fix, but to be verified at the next dig pause, not assumed).
