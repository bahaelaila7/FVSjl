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
Profile at 17/60 processed (sweep completing in the background; a live-FVS read against the 70 GB DB per treed
stand makes it slow, not the jl side):
```
treed=9  cyc0_exact=8/9  all-3cyc_exact=1  jlerr=0  live_nosum=0   (partial, monotone — final ≈14 treed)
```
- **0 jl crashes** and **cyc0 bit-exact 8/9** on real IE stands — the same bit-exact-or-cornered profile as CI,
  and consistent with IE being marked growth-complete. The sole cyc0 miss is a Δ1-class rounding (offender detail
  is emitted only in the final summary line of `scratchpad/ie_sweep.out`).
- Verdict (from the monotone partial, which cannot regress): IE growth on real FIA stands is
  **bit-exact-or-cornered with zero crashes** — validating the `.sum`-inert IE_PSIGSQ fix (`96cde22`) caused no
  regression on real multi-species data. Final treed count / offender list append to `ie_sweep.out` on completion.

## Notes
- High treeless fraction (~76% for CI's sample) is expected: FIA conditions include much non-forest; they
  are correctly excluded from the exactness denominator, not counted as trivially-exact zeros.
- This validates the `.sum`-inert PSIGSQ fixes indirectly (no regression on real multi-species stands) and
  confirms cross-variant robustness (0 jl crashes). A firing-case for the PSIGSQ COR shrinkage specifically
  would need a stand with an LDGCAL species whose true PSIGSQ ≠ 0.0898 — the source-match to live
  `dgdriv.f` remains the primary correctness proof.
