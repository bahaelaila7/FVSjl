# ADJUDICATION — IE "#206 establishment straddle" corner vs EM "real bug" (2026-09-04)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Branch: `diag-em-ie-estab` (off master 6eec330d). Measurement only — NO src changes.
Oracles: `/workspace/.iework/FVSie_g16` (04:07), `/workspace/.emwork/FVSem_g16` (15:33), both fresh.
Regime: `keytext(cn, db, "none")` = the harness's standard sweep keyfile (NUMCYCLE 5 + THINBBA 2.0 40.0),
identical keyfile to jl and oracle. Metrics: TPA / QMD / TopHt (NOT BA — BA converges by construction).

## VERDICT: (A) REAL one-directional DETERMINISTIC count/height bug. The IE #206 corner is FALSE for the
## material residual. A small genuine two-sided straddle also exists but is swamped by the deterministic bug.

The dispute crux — RNG realization straddle vs deterministic bug — is settled by the SEED-PERTURBATION test
(reseed RANNSEED identically in jl AND oracle across 5 seeds spanning the full range). A straddle varies/flips
with the seed; a deterministic bug persists. Result: **the large-magnitude residuals are SEED-INVARIANT.**

### Seed-perturbation test (last-cycle jl−oracle gap; seeds 55329/1/777/424242/987654)
| CN | variant | mode | TPA gap across seeds | TopHt gap across seeds | verdict |
|----|---------|------|----------------------|------------------------|---------|
| 2978686010690  | IE (CORNERED) | dense/post-thin | −1272,−1273,−1302,−1303,−1281 | ~0 | **DETERMINISTIC** |
| 22404926010497 | IE (CORNERED) | ultra-dense seedling | −574,−598,−621,−605,−593 | −1,+3,−5,−5,−1 | **DETERMINISTIC** |
| 3115334010690  | IE (CORNERED) | dense/post-thin | +88,+90,+92,+89,+90 | +4,+2,+3,+3,+4 | **DETERMINISTIC** |
| 3109970010690  | IE (CORNERED) | estab cohort height | −4,+3,+1,0,+63 | **−28,−24,−26,−25,−25** | **DETERMINISTIC (TopHt)** |
| 66672618010661 | EM | dense/post-thin | −73,−72,−72,−73,−73 | −1,−2,−3,−3,−2 | **DETERMINISTIC** |
| 343659700489998| EM | dense/post-thin | −270,−271,−269,−269,−270 | −4,−3,−3,−6,−2 | **DETERMINISTIC** |
| 1143317126290487| IE (CORNERED)| sparse bare-estab | −2,0,−3,+2,+3 | −5,−2,+8,0,+6 | **RNG straddle (small)** |
| 753179741290487 | IE (CORNERED)| sparse bare-estab | −68,+58,−47,+15,+2 | +2,+5,0,−1,−4 | **RNG straddle** |

Reseeding leaves 2978686 at −1300, 22404926 at −600, 3115334 at +90, 66672618 at −72, 343659700 at −270, and
3109970's TopHt at −25 — UNCHANGED. These cannot be OLDRN/ZRAND realization straddles. Only the sparse
bare-establishment TPA (1143317126, 753179741) flips sign with the seed = a genuine (but small) straddle.

## Stratified count-metric tallies (jl − oracle)

### IE — 26 cornered establishment stands, stratified across all 8 CN-suffix sections
- LAST-CYCLE:  TPA over=7 under=18 zero=1 (mean **−75.9**);  TopHt over=11 under=13 zero=2 (mean **−1.19**);  QMD over=9 under=15
- ALL-CYCLE cell:  TPA over=42 under=68 zero=46;  **TopHt over=30 under=76 zero=50** (2.5:1 UNDER)

### EM — 24 establishment stands (named Mode-1/2 + stratified DB sample)
- LAST-CYCLE:  TPA over=6 under=9 zero=9 (mean **−106.9**);  **TopHt over=1 under=20 zero=3** (mean **−2.92**);  **QMD over=17 under=0 zero=7**
- ALL-CYCLE cell:  TPA over=25 under=36 zero=83;  **TopHt over=11 under=85 zero=48** (7.7:1 UNDER)

EM TopHt is one-directionally LOW 20:1 (last-cycle) / 7.7:1 (all-cycle) and QMD one-directionally HIGH 17:0.
A two-sided RNG straddle is ~50/50 BY DEFINITION. This is categorically deterministic. IE shows the SAME
direction (TopHt under 2.5:1, count-driven QMD skew) — the SHARED signature across two variants with
independent RNG streams rules out a per-variant RNG realization.

## The three deterministic failure modes + one real minority straddle
- **D1 — post-disturbance dense tree-COUNT divergence (dominant magnitude).** Seed-invariant ±72…±1300 TPA.
  BA converges, QMD diverges (fewer trees ⇒ bigger QMD, or more ⇒ smaller). Direction is FIXED per stand but
  input-dependent (both signs across stands): jl over-kills some (2978686 −1272, 22404926 −574, 66672618 −73,
  343659700 −270), over-retains others (3115334 +90, 3304746 +134, 1856139665 +444). **NO-THIN test: the
  catastrophe is HARVEST-TRIGGERED** — removing the cyc-2 THINBBA collapses the gap (66672618 −174→+18/+37;
  2978686 −1272→small & reconverges; 3115334 +174→−6). ⇒ localizes to POST-THIN / post-canopy-opening record
  management: density-dependent mortality partition (which/how-many trees die: morts.f + VARMRT + RDPSRT
  which-trees-die) and/or post-harvest AUTOES ingrowth count.
- **D2 — established/seedling cohort HEIGHT-growth deficit.** TopHt one-directional LOW, seed-invariant
  (3109970: jl 25 vs oracle 53 = −28, identical across all seeds; the established cohort grows in height at
  ~half rate). On ultra-dense seedling stands the cohort fails to build BA at all (22404926 2054 BA 13 vs 75,
  the one case where BA does NOT converge). Kernel: established/small-tree height — ESSUBH/ESADVH AA-vs-HT1
  intercept + subsequent small-tree height growth (essubh.f / esadvh.f). Matches the EM doc's "AA-vs-HT1
  established-height class."
- **D3 — QMD one-directional HIGH (EM 17:0).** Arithmetic dual of D1 (fewer trees carrying the converged BA).
- **S1 — genuine #206 straddle (real but small).** Sparse bare-establishment TPA jitter, seed-SENSITIVE,
  sign-flips with the seed (1143317126 ±3, 753179741 ±60). This is the ONLY component that is a true OLDRN/
  ZRAND realization straddle. It is small and swamped by D1/D2 in the material stands.

## Why the IE cornering evidence was invalid (all three "proofs" fail the count metrics)
1. "BA sign-tally balanced (38/44, mean −0.85)." BA converges by construction (growth + BA-mortality correct);
   a balanced BA tally is EXPECTED under BOTH hypotheses and distinguishes nothing. The corner measured the
   ONE metric that cannot see this bug. Measuring TPA/QMD/TopHt instead exposes the deterministic residual.
2. "Adjacent same-plot pairs flip sign." These are DIFFERENT stand CNs with different inputs. D1's sign is
   input-dependent (over-kill vs over-retain per stand), so different stands flipping sign is exactly what a
   deterministic input-dependent bug produces — not RNG.
3. "Truly-identical plots give identical residual (3315450 +22 == 3315451 +22)." Identical inputs → identical
   output is what DETERMINISM predicts. This is evidence FOR a deterministic bug, mis-read as an RNG proof.
4. The IE agent never ran the seed-perturbation test. It shows the material residual is seed-INVARIANT.

## Confirmation that cornered stands carry the deterministic bug
2978686010690 (−1272, seed-invariant), 22404926010497 (−574), 31702493010690 (−229…−635), 3115334010690
(+90) are all present in `docs/fia_cornered_stands.txt` tagged "IE establishment #206 OLDRN/ZRAND straddle"
— yet are seed-invariant deterministic count bugs, several not even establishment-driven (they start fully
stocked at hundreds–thousands TPA; the divergence is in post-thin mortality/regen, not the AUTOES cohort).

## Bottom line
The EM sweep was correct: this is a REAL one-directional deterministic tree-COUNT / cohort-HEIGHT bug shared
by EM and IE (Wykoff-variant establishment/mortality path), NOT a #206 realization straddle. The IE
"COMPLETE, cornered" status is a FALSE corner for the material residual (D1+D2) and must be un-cornered; the
2773 (+587) IE establishment corners should be reopened. A genuine small sparse-establishment straddle (S1)
does exist but is not what dominates the cornered population. Recommended fix targets: (D1) post-disturbance
density-dependent mortality-partition + post-harvest ingrowth record count; (D2) established/small-tree
height (ESSUBH/ESADVH AA intercept + small-tree height growth). Fix is EM+IE and likely Wykoff-variant-general.
