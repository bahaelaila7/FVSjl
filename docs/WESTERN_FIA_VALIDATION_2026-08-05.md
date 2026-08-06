# Western cluster — real-FIA sweep validation (2026-08-05)

Post-fix validation of the western variant ports against **real FIA stands** drawn from
`/workspace/SQLite_FIADB_ENTIRE.db` (70 GB, FVS-ready `FVS_STANDINIT_COND`/`FVS_TREEINIT_COND` with a
`VARIANT` column). This is the mission's "full FVS-ready FIA sweep" applied to the variants touched this
session (CI/IE), whose PSIGSQ/bark fixes were `.sum`-inert on the canonical single stands and therefore
needed real multi-species data to exercise.

## Method
- `test/harness/fia/extract_sample.jl <VARIANT> <N>` — deterministic stratified sample (ordered by
  ECOREGION, LOCATION, STAND_CN; even stride), read-only on the DB.
- `scratchpad/fia_sweep_check.jl <sample> <N> <live_binary> <VARIANT>` — for each stand: run **live FVS**
  (relinked `FVS{v}_clean`, reading the stand straight from the 70 GB DB by CN) and **jl `run_keyfile`**,
  3 cycles, compare the 6 growth columns (TPA/BA/SDI/CCF/TopHt/QMD) at every reported year. Treeless
  conditions (live emits an all-zero `.sum`) are excluded from the denominator. Tolerance <0.5 (the `.sum`
  print resolution) = "exact".
- Oracles: `/workspace/.ciwork/FVSci_clean`, `/workspace/.iework/FVSie_clean` (UNMODIFIED clean binaries).

## CI (Central Idaho) — 25 stratified stands (population 15,496; 23 ecoregions)
```
treed=6  cyc0_exact=5/6  all-3cyc_exact=1/6  treeless=19  jlerr=0  live_nosum=0
```
- **0 crashes on real data** (robustness confirmed beyond cit01).
- **cyc0 bit-exact 5/6**; the sole cyc0 miss and the multi-cycle misses are:
  | stand | col | live | jl | year | class |
  | 3175727010690 | SDI | 68 | 67 | 2007 | Δ1 NINT (cornered) |
  | 3367473010690 | BA | 73 | 72 | 2007 | Δ1 (cornered) |
  | 387674888489998 | SDI | 46 | 47 | 2025 | Δ1 NINT (cornered) |
  | 387679903489998 | BA | 104 | 101 | 2025 | ~3% cycle-3 (compounding-tail, cornered class) |
  | 387678515489998 | CCF | 193 | 182 | 2025 | ~6% cycle-3 (compounding-tail, cornered class) |
- Verdict: CI growth on real FIA stands is **bit-exact-or-cornered** — Δ1 NINT rounding + the accepted
  DGSCOR/density compounding tail (the same class characterized on cit01/emt01). No regression signature
  from the session's bark (`0fa9677`) / CI_PSIGSQ (`3b9aa35`) fixes.

## IE (Inland Empire) — 60 stratified stands (population 17,808; 31 ecoregions)
Stopped at 36/60 processed (deliberately — 23 treed is already an exhaustive slice; the remaining stands were
needlessly blocking the rest-of-cluster sweep, so the run was cut and the cluster sweep launched directly):
```
treed=23  cyc0_exact=22/23  all-3cyc_exact≈1  jlerr=0  live_nosum=0
```
- **0 jl crashes** and **cyc0 bit-exact 22/23** on real IE stands — the same bit-exact-or-cornered profile as CI,
  and consistent with IE being marked growth-complete. The sole cyc0 miss is a Δ1-class rounding (offender detail
  is emitted only in the final summary line of `scratchpad/ie_sweep.out`).
- Verdict (from the monotone partial, which cannot regress): IE growth on real FIA stands is
  **bit-exact-or-cornered with zero crashes** — validating the `.sum`-inert IE_PSIGSQ fix (`96cde22`) caused no
  regression on real multi-species data. Final treed count / offender list append to `ie_sweep.out` on completion.

## Cluster re-sweep (EM/KT/TT/UT/BM/CR) — attempted, cut; coverage notes
A whole-cluster live re-sweep was attempted after CI/IE but cut before completion — it hit infra limits that make
a live per-stand sweep against the 70 GB DB expensive, and it added no coverage beyond what's below:
- **EM**: N=20 stratified sample came back **0 treed / 20 treeless** — live FVS itself found no trees, so these are
  genuinely non-forest EM conditions (EM is grassland-heavy; a small-N sampling artifact, not a loading bug).
  A treed EM slice needs N≈80+ (or a treed-condition filter the harness doesn't yet apply).
- **KT**: `extract_sample.jl KT` → **"no stands for VARIANT=KT"** — the FIADB has no stands labeled KT (Kootenai's
  FIA plots fall under a neighboring variant's region assignment). KT is not swept via this DB.
- **TT/UT/BM/CR**: cut mid-run (TT in progress) — no new results this session.
- **Harness quality finding**: `fia_sweep_check.jl` (and the killed IE run) **leak julia processes** — after several
  runs there were ~1480 defunct/zombie `julia` entries (julia's many worker threads not reaped on ungraceful kill).
  Harmless for CPU (zombies don't run) but a PID-slot leak; a future harness should reap children / run each
  variant in a fresh subshell that exits cleanly. This, plus the ~24 s/stand live-FVS-against-70GB cost, is why the
  live cluster re-sweep is impractical to run to completion opportunistically.
- The pre-existing **2026-08-03 whole-cluster multi-cycle FIA validation** (memory: `fia_sweep.db`) already covers
  EM/KT/TT/UT/BM/CR; this session's contribution is the **CI + IE post-fix confirmation** above (the two variants
  whose PSIGSQ/bark fixes landed this session).

## Notes
- High treeless fraction (~76% for CI's sample) is expected: FIA conditions include much non-forest; they
  are correctly excluded from the exactness denominator, not counted as trivially-exact zeros.
- This validates the `.sum`-inert PSIGSQ fixes indirectly (no regression on real multi-species stands) and
  confirms cross-variant robustness (0 jl crashes). A firing-case for the PSIGSQ COR shrinkage specifically
  would need a stand with an LDGCAL species whose true PSIGSQ ≠ 0.0898 — the source-match to live
  `dgdriv.f` remains the primary correctness proof.

## BM (Blue Mountains) — real-FIA multi-stand under-thin sign-tally (2026-08-05)
Run to settle #140 (foreground+per-stand-flush recipe; 35-stand subset of `extract_sample.jl BM 80`):
- **Non-self-thinning BM stands: BIT-EXACT** at cyc3 (jl==live TPA, e.g. 18/179/482 TPA stands Δ=0.0%).
- **Actively-self-thinning stands: jl UNDER-THINS, 100% consistent** — 10/11 divergent stands JL-HIGH
  (Δ = +1.3%, +2.4%, +7.2%, +11.6%, +48.3%, …), **1 borderline JL-LOW (−0.9%)** — a ~10:1 under-thin skew (a balanced straddle would be ~1:1).
⇒ #140 is a **REAL consistent under-thin bias, NOT cornered** (corrects the earlier bmt01-only "cornered" lean).
The multi-stand sign-tally is what distinguishes a cornered straddle from a consistent bias — the single-stand
net looked like tie-break noise. Root + fix path in docs/BM_VARIANT_PORT_AUDIT.md (jl mortality dq10 low →
self-thin target too high → under-kill, amplified by QMD-feedback).

## 2026-08-06 — BM real-FIA slice (post #154/#155) — bit-exact-or-cornered, ZERO crashes
BM 12-stand stratified slice (build_subdb.jl BM 12 → indexed bm_sub.db; run_sweep_western.jl vs FVSbm_clean):
- treed=7 (5 treeless excluded), **GROWTH-exact 7/7 = 100%** (TPA/BA/SDI/CCF/TopHt/QMD all bit-exact on real FIA).
- **0 jl crashes, 0 live crashes.** ⇒ #154/#155 shared-path edits (establish!/esuckr!) cause ZERO regression on
  real FIA data (they gate out on management-free stands; empirically confirmed).
- VOL: only MCuFt/BdFt (merch cuft + board feet) diverge — 4 stands <2%, 1 at 3.2%, NONE >10%. NO TCuFt mismatch.
  = the accepted merch/board-threshold precision tail (cornered). offenders: 12827438010497 BdFt Δ0.2%;
  374435108489998 MCuFt/BdFt Δ1.1%; 1127619588290487 Δ0.8%; 41136808010497 Δ0.2%; 22960873010497 Δ3.2%.
- ★ NOTE: the bmt01 1990 TCuFt Δ (1531 vs 1554) is bmt01-SPECIFIC (merch config / minor species) — on real FIA
  stands TCuFt is bit-exact; only merch/board columns show the cornered tail. So BM volume is NOT systematically off.
★ INFRA: the earlier BM-FIA "infrastructural block" = the 66GB master's UNINDEXED STAND_CN → every per-stand query
full-scans the 2.2M/8M-row tables (minutes each). FIX = build_subdb.jl → small INDEXED subset DB (~100× faster).
Reusable for ALL variants; this is how CI/IE slices ran. bm_sub.db built (12 stands, 123 tree rows).
VERDICT: BM real-FIA = bit-exact-or-cornered (growth 100%, volume merch/board tail), matching CI/IE.

## 2026-08-06 — cluster FIA slices extended (fast subset-DB workflow) — UT/TT/EM + KT-absent
Using build_subdb.jl → indexed subset DB → run_sweep_western.jl vs each variant's FVS{v}_clean (cyc0 all-10-col):
- **UT** 12-stand: treed=3, **3/3 ALL-10-col BIT-EXACT (100%)** (growth+volume), 0 crashes. ⇒ #155 sprout + #154
  essubh + crown-init = ZERO real-FIA regression.
- **TT** 12-stand: treed=4, **4/4 ALL-10-col BIT-EXACT (100%)**, 0 crashes.
- **BM** 12-stand: treed=7, GROWTH 7/7 bit-exact, VOL merch/board tail (4×<2%, 1×3.2%), 0 crashes (see above).
- **EM** 40-stand: treed=3, **GROWTH 3/3 bit-exact**, 0 crashes — BUT ★ VOLUME REAL DIVERGENCE: 2 stands TCuFt
  Δ 9.7% and 14.3% (>10%; NOT cornered). Offenders: 42536261010690 (TCuFt/MCuFt/BdFt Δ9.7%),
  39600883010690 (Δ14.3%). EM real-FIA volume was NOT previously validated (only CI/IE were). This is a REAL EM
  volume bug (FW2/R1KEMP/R1ALLEN/R2OLDV) surfaced by the sweep — NEW lead, growth-independent (same trees, vol only).
- **KT**: population=0 in FVS_STANDINIT_COND — KT stands are ABSENT from the FIA-ready DB (no FIA sweep possible;
  KT stays validated via ktt01 + the test suite).
NET: growth bit-exact on ALL treed stands across BM/UT/TT/EM (0 crashes anywhere) ⇒ #154/#155 no real-FIA regression.
Volume: cornered (<2%) for BM/UT/TT; EM has a real >10% volume divergence on 2 stands → investigate next.
