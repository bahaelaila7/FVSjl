# Western cluster — real-FIA MORTALITY population sweep (2026-08-06)

Cluster-wide validation of jl self-thin/mortality against **real FIA stands** vs the live relinked oracles,
extending the BM #140 investigation to the whole western cluster. Method (per variant): stratified N-stand sample
from the 70 GB `FVS_STANDINIT_COND` (VARIANT column) → indexed sub-DB (`{v}_sub.db`, ~100× faster than the master)
→ run live `FVS{v}_clean` + jl `run_keyfile(variant)` 3 cycles → compare TPA/BA/QMD at the final year. mean|Δ% vs
live| over the genuinely-treed stands (jl>0 & live>0); "AUTOES-0" = jl 0 TPA where live auto-established (#143).
Harness: `test/harness/fia/extract_sample.jl` + the `bm_*`/`em_*`/… scratchpad sweep scripts (reusable).

## Results (this session)

| variant | N | treed | AUTOES-0 | TPA mean\|Δ%\| | TPA sign (>2%) | verdict |
|---------|---|-------|----------|----------------|-----------------|---------|
| **BM** | 80 | 56 | 0 | 3.45 → **0.96** | 31→24 HIGH (fixed) | ★ FIXED (BMTMRT+IPASS routing, commit 56d0931) |
| **UT** | 150 | 43 | 0 | **1.16** | 3 HIGH : 7 LOW | bit-exact-or-cornered — no fix warranted |
| **CI** | 150 | 45 | 0 | **2.28** | ~3 HIGH : 10 LOW | mild OVER-KILL lean (at bar; small-tree stochastic on dense) |
| **TT** | 150 | 56 | 0 | **0.58** | 1 HIGH : 3 LOW | bit-exact-or-cornered — excellent |
| **EM** | 200 | ~18 | **~180** | (AUTOES-dominated) | treed: 5 HIGH : 2 LOW | ★ AUTOES gap dominates (#143); treed mortality under-thins (hybrid #137) |

- **KT** — 0 stands in the FIADB (not DB-sweepable; the goal-file note confirmed).
- **IE / CR** — prior work (IE cyc0 bit-exact + AUTOES residual #143; CR DB-sweep 39/40 done).
- **jl robustness: 0 crashes** across all 730 stands swept (after the BM regent-HTDBH crash fix e278edd).

## Verdicts
1. **BM** self-thin under-thin was a REAL systematic bug — FIXED (routed through the shared BMTMRT+IPASS driver;
   TPA 3.45%→0.96%, BA 8.65%→6.45%, QMD 6.12%→5.29%, 0 aggregate regressions).
2. **UT / TT** mortality is genuinely bit-exact-or-cornered on real FIA (1.16% / 0.58% TPA) — no systematic bias,
   no fix. The BM fix does NOT generalize (measured, doctrine #2 — refuted the "UT likely under-thins" hypothesis).
3. **CI** — the goal-file "~2% straddle" is at population scale a MILD OVER-KILL LEAN (TPA 2.28%; ~3 HIGH:10 LOW >2%),
   not a balanced straddle — the same characterization trap BM's "cornered straddle" verdict fell into. Magnitude is
   at the bar; outliers are the open SMHTGF/DGSCOR small-tree stochastic on dense stands. Verdict refined, not a new
   bug. To close the lean: instrument CI mortality vs ci/morts.f on JL-LOW dense stand 1856022712290487.
4. **EM** — the DOMINANT real-FIA divergence is the **AUTOES establishment gap (#143)**: ~180/200 stands jl=0 where
   live auto-establishes 100-260 (corrects the "EM validated cyc0-bit-exact" claim = both treeless at cyc0). The
   genuinely-treed EM stands under-thin like BM, but EM's mortality! is a HYBRID (SDI + added-species Hamilton) so
   the BM routing needs adaptation (#137).

## Cluster mortality state after this sweep
BM fixed · UT/TT/CR bit-exact-or-cornered · CI at-bar (mild lean, refined) · EM = AUTOES-gap-dominated + hybrid
under-thin (#137). The remaining OPEN mortality/establishment items, in priority order: **#143 AUTOES** (the single
biggest real-FIA fidelity gap cluster-wide — EM/grassland-heavy variants) > **#137 EM hybrid self-thin** > CI lean
closure (optional, at bar). AUTOES = shared estb/ logic + variant data; jl has IE only; extending to EM = a fresh
focused chunk (extract EM estb/ data tables + generalize ie_autoes).
