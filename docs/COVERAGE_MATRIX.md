# FVSjl per-variant coverage matrix

Pillar-1 done-state deliverable: a documented per-variant coverage matrix showing that the
**100% drop-in** claim is defensible for **SN + NE + CS + LS**, every exercised path validated
**bit-exact (`==` vs freshly-relinked live FVS) or cornered to a named primitive**.

Correctness floor (never regress): **38515 pass / 143 broken / 0 fail** (`julia --project=. test/runtests.jl`).
(38479/141 → 38511 via S96 NE snag height-loss + S97 NE carbon-report; → 38513 via S100 Pillar-2
allocation-floor guard; → 38515 via S101 Pillar-4 hot-path type-stability guard.)
The `broken` are the documented cornered set (ULP-class / FVS-UB / eigensolver-tie / accepted primitives).

## Coverage by variant

| Axis | SN | NE | CS | LS |
|------|----|----|----|----|
| **Cycle-0 inventory** (TPA/BA/SDI/CCF/TopHt/QMD/volumes) | bit-exact | bit-exact | bit-exact | bit-exact |
| **Canonical multi-cycle** (`*t01`) | snt01 (+ FFE/sprout/estab/econ) | net01 | cst01 | lst01 (+estab/ffe/fire_sprout) |
| **Growth spine** (DG/HTG/mortality/crown) | bit-exact¹ | bit-exact¹ | bit-exact¹ | bit-exact¹ |
| **Isolated-keyword scenarios** (KCV harness) | via snt01/harness | **65** | **65** | **65** |
| **Volume** (cuft/merch/saw + board) | bit-exact³ | bit-exact | bit-exact | bit-exact² |
| **FFE fire/fuel/snag/carbon** | bit-exact-or-cornered | SIMFIRE bit-exact (KCV) | SIMFIRE cornered* | SIMFIRE cornered* + bit-exact (lst01_ffe) |
| **FIA-DB reader** (real inventory) | bit-exact | bit-exact | bit-exact | bit-exact |
| **FIA-DB, all 10 `.sum` cols, 162 real stands** | — | part of the 159/162 (7/10 cols perfect on all 162) | — | — |

¹ bit-exact for uncalibrated species at the variant's native cycle length; the accepted residual is the
  non-native-cycle-length DGSCOR drift (documented) and the cross-cutting WK3 DGSCOR tail.
² LS board is Scribner OR International ¼" per national forest (S43 `_R9_INTL_BDFT_FORESTS`, `volinit.f`).
³ SN saw/board is IFOR-aware incl. the North Carolina (IFOR=11) merch overrides (S69 `_sn_merch`,
  setcubicdflts.f region-8) — bit-exact vs live on both default and NC forests.

## KCV isolated-keyword coverage (NE/CS/LS — 65 each, `test/fixtures/kwcov/`, `test_kwcov_variants.jl`)

Each = the variant's canonical tree data + ONE isolated keyword, validated bit-exact-or-cornered vs live:

baimult · bamax · bfdefect · bfvolume · compress · compute · crnmult · cuteff · cycleat · dgstdev · estab · eventmon · fertiliz · fixhtg ·
fixdg · fixmort · htgmult · htgstop · leavesp · managed · mcdefect · mcfdln · bffdln · minharv · mortmsb · mortmult · mult · nocalib · notriple · numtrip ·
rannseed · readcord · readcorh · readcorr · resetage · salvage · sdicalc · sdimax · serlcorr · simfire · setsite · specpref · spgroup · spleave · strclass · tcondmlt ·
tfixarea · thinaba · thinata · thinbba · thinbta · thincc · thindbh · thinht · thinmult · thinprsc · thinpt · thinqfa · thinrden · thinsdi · timeint · topkill · treeszcp · volume · yardloss

★ Keyword INTERACTIONS: {ne,cs,ls}_thinfire (thin-then-fire, cut-residue→fire-fuel) are FULL-ROW bit-exact vs live (S92).

(SN is exercised through the broader snt01 multi-stand harness + the SN keyword-coverage suite rather than
the KCV isolated set.)

## FIA-database mass validation (audit S44/S46/S47/S48)

Native `DATABASE`/`DSNIN` reader validated on **162 real cross-variant stands** (SN/NE/CS/LS) extracted
from `SQLITE_FIADB_ENTIRE.db` (read-only, never modified): **159/162 bit-exact on all 10 `.sum` columns**;
TPA/BA/SDI/QMD/MCuFt/SCuFt/BdFt (**7/10**) bit-exact on **every** stand. Fixes landed this campaign:
S43 (LS board-type), S45 (SN missing-elevation→forkod Hopkins index), S47 (SN seedling height ≤0.1→1.01),
S48 (AVH top-height RDPSRT tie-break). The 4 remaining residuals are cornered: 1 TopHt (tree-storage-order
tie), 1 CCF (Δ1 ULP), 2 TCuFt (Δ1-8 ULP).

## Named cornered residuals / accepted deferrals (not bugs)

- **COMPRESS eigensolver/RNG-tie** (s22) — IBM EIGEN + RDPSRT sub-ULP partition tie; faithful port.
- **Non-native-cycle-length DGSCOR** drift — Float32 serial-correlation, variant-native bit-exact.
- **CS SETSITE** — sub-ULP height amplified by the ill-conditioned NC-128 anamorphic inversion at raised
  site; re-converges bit-exact @2040.
- **THINAUTO on NE/CS/LS canonical stands** — live FVS is ill-posed (FPE / no data — full-stocking UB);
  jl handles gracefully. SN THINAUTO validated.
- ~~**Soft-snag soft-DDW residual** (#73)~~ — **RESOLVED (S78, bit-exact vs live).** The soft CWD1 cone-split
  now uses the FVS per-hardness LOHT (fmcwd.f: soft LOHT(1)=1.0 in both the R1 cone-base radius :347 and the
  LOCUT floor; hard LOHT(2)=0.10). `_cwd_cone_fractions` returns `(frac_soft, frac_hard)`, soft normalized by
  the invariant pat_hard(0.10) base. Live-validated on carbon_snagpsft.key: DDW 5.8/5.2/7.9/10.7 (was 8.0/10.8)
  == live on every carbon column. Floor safe — ordinary/fire snags are hard (DFIS=0 ⇒ soft split ×0), so
  carbon_snt + fire scenarios are bit-identical. Regression fixture: `carbon_snagpsft.*` + test_carbon.jl.
- **FIA-reader ULP cells** — TopHt tree-storage-order tie, CCF/TCuFt integer-rounding boundaries.
- **LS PERCOV** (~3.4) — forest-grown crown CR-timing; feeds only the coarse cover-class-binned fire WMULT.
- **ls_simfire FFE fire-mortality — FIXED (S87).** The LS ~5-TPA under-kill was jl applying the SN fuel
  decay-rate table (litter DKR 0.65) to LS (which uses 0.31, ls/fmvinit.f) ⇒ litter ~2× low ⇒ SMALL down-wood
  low ⇒ FMDYN under-weighted the hot fuel model ⇒ under-scorch. `_fm_dkr_default(::LakeStates)` now uses the LS
  DKR table (variant-gated; SN/NE/CS byte-identical). Live-validated: ls_simfire 2020/2040 full-row bit-exact,
  2030 QMD 11.1/11.2 sub-print straddle only (ULP-class @test_broken). Suite 37633/140/0.
- **cs_simfire FFE fire-mortality — FIXED (S89).** After the CS variant DKR table (`_fm_dkr_default(::CentralStates)`
  = cs/fmvinit.f), cs_simfire TPA + BA are now BIT-EXACT vs live all cycles; only a post-fire volume ULP (<0.1%,
  MCuFt/SCuFt/BdFt Δ1-3) remains — an accepted per-tree Clark straddle, @test_broken. NE SIMFIRE is bit-exact.
- **ls_thinpt point-thin** — THINPT/SETPTHIN is BIT-EXACT vs live (2010 removal + 1990-2020); only 2 downstream
  ULP straddles on the heavily-thinned residual stand (2030 MAI 35.2/35.1, 2040 BdFt 7464/7463). NE/CS thinpt
  full-row bit-exact (S91).

## Known parsed-but-unported keywords (honest scope of the drop-in claim)

These FVS keywords are recognized but **not** acted on by FVSjl. Each was investigated and found to be
either `.sum`-invisible (so an isolated KCV test would be vacuous — violates "test must exercise the
semantic") or a specialty/extension path outside the exercised eastern scope. Documented so the
"100% drop-in" claim is honest rather than silently incomplete.

**Coverage-probe conclusion (S99):** the two findable common unimplemented keywords (FIXCW, PRUNE) were
both *measure-verified* `.sum`-invisible (FIXCW at 10×; PRUNE fired-but-identical) — their effects live in
crown-width/fuel/fire pools, never the standard summary. This is the strongest completeness evidence for
the `.sum` drop-in claim: the remaining gaps do not move any `.sum` column on the exercised scenarios.

- **FIXCW** (crown-width multiplier, `cwidth.f` activity 90) — investigated S99. Multiplies per-tree
  crown WIDTH within a (species/group, DBH-window). Effect is confined to PERCOV/crown-width `.out`
  reporting and the coarse cover-class-binned fire WMULT — **fully `.sum`-invisible** (proven: base vs
  FIXCW-10× on the NE canonical stand → byte-identical `.sum`, every column, every cycle). The `.sum`
  CCF column is the DBH-based crown-competition factor, not crown width. Porting is a per-tree multiplier
  threaded through the ~9 `crown_width` consumers, inert by default; deferred because it is `.sum`-
  invisible and its only observable (fire-PERCOV) is coarse-binned — no non-vacuous isolated test exists.
- **PRUNE** (`cuts.f` activity 249 / `fmprun.f`) — investigated S99, **measured `.sum`-invisible**. Prunes
  lower crown to a height (method/feet/maxCRprop/species/DBH-window) and books the pruned crown material
  into the CWD fuel pools (fmprun.f). On the NE canonical stand the activity FIRED ("PRUNE … DONE IN 2000"
  in the `.out`) yet the `.sum` was byte-identical to baseline every column/cycle — its effect is confined
  to fuels→fire and crown/fuel reporting. Non-vacuous isolated `.sum` test does not exist without a
  fire-coupled scenario (and even then via the coarse-binned fuel path). Specialty FFE-fuel keyword.
- **Extension keywords** (`DFB`/`DFTM`/`MPB`/`BRUST`/`WSBW`/`MISTOE`/`CLIMATE`/`ORGANON`/`FIAVBC`) —
  insect/disease/climate/growth-model extensions outside the SN+NE+CS+LS core scope.

Full slice-by-slice detail: `docs/MODERNIZATION_AUDIT.md`.
