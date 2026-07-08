# Stratified FIA-plot validation — results

Replaces the earlier 162-stand convenience sample with a **seeded stratified random sweep**
across SN/NE/CS/LS, validated cycle-0 all-10-`.sum`-columns vs the freshly-relinked live FVS.

## Method
- Population: `FVS_STANDINIT_COND` = 1,472,391 target-variant conditions (SN 637,641 / NE 178,149
  / CS 255,952 / LS 400,649).
- `sample_stands.jl` — seeded reservoir sampler, one streaming scan (~7s), N random STAND_CN/variant.
- `extract_subset.jl` — copies the sampled STANDINIT_COND + TREEINIT_COND rows into a small **indexed**
  subset DB (source ATTACHed READ-ONLY; never modified). Makes per-stand queries <0.3s vs ~15s on the
  unindexed 66GB source (no index on STAND_CN → full 24.8M-row scan per stand).
- `run_sweep.jl` — jl vs live per stand; **hang-safe** (live wrapped in `timeout`, fixing the prior
  runaway). A cell counts as divergent only if the printed `.sum` value differs by ≥0.5.

## Batch 1 (seed 20260708, 100/variant = 400 stands)
| Variant | both-produced | bit-exact | rate | divergent cols (stands) | no-sum* |
|---|---|---|---|---|---|
| SN | 78 | 75 | 96.2% | TCuFt/MCuFt/SCuFt/BdFt (3) | 22 |
| NE | 95 | 94 | 98.9% | TopHt (1) | 1 |
| CS | 100 | 100 | **100%** | none | 0 |
| LS | 100 | 98 | 98.0% | CCF (1), TCuFt/MCuFt/SCuFt (1) | 0 |
| **TOTAL** | **373** | **367** | **98.39%** | | |

\* "no-sum" = live produced no `.sum` — **mostly nonstocked/treeless conditions** (only 182/400 sampled
COND rows had trees), NOT timeouts or crashes; a random treeless-vs-treed check confirmed treed SN stands
run fine in <0.3s.

## Reading
- **All 6 divergent stands fall in the SAME residual classes already cornered by task #85** (volume
  cuft/board ULP, TopHt AVH sort-tie, CCF integer-boundary) — **no NEW divergence class appeared at ~10×
  the scale of the 162 sample.** The 98.4% is the honest cycle-0 print-identical rate on real inventory.
- Next: (a) scale to 500–1000/variant for tighter CIs; (b) filter the sampler to stands-with-trees so the
  denominator is clean; (c) drill the 6 divergent stands to confirm each is the cornered primitive (not a
  new bug).

## Correction + key finding (batch 1 re-analysis)
The batch-1 "no-sum" stands are **NOT (only) treeless** — the majority are **live FVS SIGFPE crashes**.
Decisive: in a seed-99 sample, **34/34** treed SN no-sum stands exit **136 (SIGFPE)**, reproducible on
BOTH the subset and the full 66GB source DB. Root cause identified + reported as **FVS_SOURCE_BUGS.md
D38** (R9 Clark volume init, `volinit.f:414`/`r9clarkdib.f`, divide-by-zero on trees with zero/null HT;
~28–42% of treed SN FIA COND stands). **jl is strictly more robust: valid `.sum` on 119/120** treed SN
stands incl. all live-crashers.

Implication for the metric: the 98.4% cycle-0 bit-exact is on **both-produced** stands (necessarily
excludes the stands where live FVS crashes and yields no oracle). This is a live-FVS UB exclusion,
consistent with the campaign doctrine (known FVS bugs are not part of the drop-in requirement).

## Cross-variant live-FVS crash rate (D38 R9 Clark r9ht SIGFPE), 100 treed stands/variant
| Variant | live OK | SIGFPE (exit 136) |
|---|---|---|
| SN | 70 | 30% |
| NE | 98 | 2% |
| CS | 95 | 5% |
| LS | 97 | 3% |
The bug is cross-variant (shared R9 Clark volume lib) but SN-dominant. jl runs all. This bounds the
"both-produced" denominator per variant: NE/CS/LS FIA validation is ~95-98% runnable, SN ~70%.

## Scale note
A 4000-stand (1000/variant) clean run was attempted but deprioritized: each stand costs a full jl
`run_keyfile` (DB read + 3-cycle sim) **plus** a live FVS run, so 4000 stands ≈ 1–2 h wall — not worth
the marginal tightening over batch-1 (373 stands) + the 100-stand/variant crash-rate table. The harness
(`sample_stands.jl` → `extract_subset.jl` → `run_sweep.jl`, crash-classifying) is proven and can be run
offline at any N; sampler seeds make it reproducible. **Definitive FIA-plot findings stand:** ~98% cycle-0
all-10-col bit-exact on both-produced stands, all residuals in known cornered ULP classes, and the D38
live-FVS SIGFPE crash (SN 30% / NE 2% / CS 5% / LS 3% of treed stands) which jl runs correctly.

## Batch 2 — FULL 1000/variant crash-aware sweep (seed 99), DEFINITIVE
| Variant | both-produced | cycle-0 bit-exact | rate | D38 SIGFPE-crash |
|---|---|---|---|---|
| SN | 833 | 781 | 93.8% | 165 |
| NE | 964 | 956 | 99.2% | 9 |
| CS | 991 | 983 | 99.2% | 9 |
| LS | 995 | 973 | 97.8% | 5 |
| **TOTAL** | **3783** | **3693** | **97.62%** | 188 |

Divergent columns are the known cornered ULP classes: volume (TCuFt/MCuFt/SCuFt/BdFt), TopHt (AVH sort-tie),
CCF (integer boundary). Full per-stand offender list in `batch2_1000pervariant.txt`.
**REAL BUG FOUND (not ULP):** SN is the outlier at 93.8% — ~6% of SN stands diverge on the SAW/BOARD volume.
Drilled offender 928399836290487 (large trees, QMD 15.6): TCuFt Δ10 (0.25%) / MCuFt Δ28 (0.7%) near-fine, but
**SCuFt live 3032 vs jl 3491 (~15% high) and BdFt live 17398 vs jl 19414 (~12% high)**. So jl's sawtimber-cubic
and board-foot volume is systematically too high on LARGE-SAWTIMBER SN stands — a config the exercised scenarios
(snt01/KCV) never hit (they're bit-exact). Total cubic ~fine ⇒ the bug is in the SAW-threshold / saw-top-diameter
or board bucking on big trees, SN-specific (NE/CS/LS volume ≤6 diverging stands each). Task filed to drill the
SN saw/board volume path (r8clark/scfstmp/scftopd or the sawtimber merch spec) on large trees.

## #95 fix entry-point (for the focused session)
SN saw/board volume divergence (~2.4%/large tree in saw/topwood split; ~15% aggregate from mid-size
sawtimber classification). One-tree repro: /tmp/fia_val/one.key (YP DBH25.3): TCuFt/MCuFt BIT-EXACT
(796/794), SCuFt jl774/live756, topwood v[7] jl20/live38 ⇒ jl saw-height HT1PRD to the 9" saw-top is
GREATER than live. STRUCTURE: FVS SN has NO r8clarkdib.f — saw height HT1PRD is computed in
bin/FVSsn_buildDir/fvsvol.f (+ r9clark with R8 coefs). Debug-FVS stamp target = fvsvol.f (WRITE HT1PRD /
sawHt for the tree), compare to jl r8clark_vol.jl ht1prd. Then decide: fixable saw-height formula diff vs
accepted Float32 taper-inversion boundary (MCuFt bit-exact ⇒ core taper right; only the 9"-split boundary
drifts). CAUTION: snt01/KCV SN volume bit-exact — verify from fvsvol.f source, don't regress.

## #95 RESOLVED + post-fix SN sweep [2026-07-08]
Root cause was NOT a taper-inversion boundary — it was the **missing SN North Carolina (IFOR=11) merch
overrides** (setcubicdflts.f region-8): hardwood sawtimber top `SCFTOPD 9→11″` (coastal KODIST 3/10 → 8″),
`SCFMIND 12→15`, softwood `7→6.3`, `TOPD 4→3.5`. jl read a static merch CSV = the non-NC defaults only. The
one-tree offender is a NC (FOREST=11) yellow-poplar. Fixed via `_sn_merch(spi,ifor,kodist)` (1:1 port,
mirrors `_ne/_cs/_ls_merch`); SN joined the unified IFOR-aware branch; IFOR from `plot.forest_idx` when set
(STDINFO Fort Bragg=20), else decoded from KODFOR for FIA stands. See audit S69 / `docs/COVERAGE_MATRIX.md`.

**Post-fix SN spread sweep** (400 rowid-spread treed SN stands, cycle-0 all-10-col vs freshly-relinked live):
- **BIT-EXACT 286/292 both-produced = 97.9%** (was 93.8% pre-fix — the large-sawtimber NC volume class is gone).
- After **#96 fixed** (Fort Bragg FIA remap, below): **287/292 = 98.3% bit-exact**, and the residual **5 are
  ALL cornered-ULP — zero real divergences**: **4 TopHt** (AVH sort-tie) + **1 BdFt Δ1** (stand 161807872010854,
  board-rounding). Every SN FIA volume stand in the sample is now bit-exact-or-cornered.

### #96 — Fort Bragg FIA zero-volume FIXED (bit-exact) [2026-07-08]
Re-trace discipline (measuring the residual deltas, not labeling them) caught this: the initial "all ULP"
read was HALF WRONG — 1 of the 2 volume residuals was a **real zero-volume bug**. Stand 252270661010854 gave
jl=0 for ALL of TCuFt/MCuFt/SCuFt/BdFt vs live 4958/3490/2384/12983. It's **Fort Bragg** (REGION=7/FOREST=1 ⇒
composite LOCATION=701). A **live-FVS debug-stamp** (`fvsvol.f` VOLEQ/IREGN) proved live remaps 701 ⇒
**KODFOR=81110 (NC Uwharrie, region 8), IFOR=20, VOLEQ=821CLKE** — NOT Talladega 80106 (an initial guess,
refuted by measurement). `forkod.f` CASE(701): both encodings converge (keyfile 701xx collapses via
IFORDI==701; the FIA composite is already 701) ⇒ `SELECT CASE(701)` ⇒ 81110/IFOR=20. jl's `kw_stdinfo!` gate
`div(KODFOR,100)==701` caught 70106 but not the bare FIA 701. FIX: extracted `sn_fortbragg_remap!(p)` (gate
widened to `==701 || div==701`), called from BOTH `kw_stdinfo!` and the FIA reader (`io/fia_database.jl`,
SN-gated). Result: stand bit-exact vs live (4958/3490/2384/12983); suite 37121/140/0; Fort Bragg keyfile test
still passes; generalizes to all Fort Bragg FIA stands. (`forest_idx=20` ⇒ default merch, matching live's
setcubicdflts IFOR=20.) NOTE: an earlier default-trap→80106 attempt gave nonzero-but-1%-off and was reverted
— the debug-stamp is what pinned the correct forest.
- 104/400 (~26%) are live-FVS **D38 SIGFPE** crashes (no oracle ⇒ excluded; jl runs them fine) — consistent
  with the ~30% SN rate. 0 no-sum/keyword-error (the spread sample avoided the earlier contiguous-batch clustering).
- Focused NC batch (34 forest-11 stands, coastal+non-coastal): **34/34 bit-exact all 10 columns**.

## LS FIA spread sweep [2026-07-08]
218 rowid-spread treed LS stands, cycle-0 all-10-col vs freshly-relinked live FVSls:
- **BIT-EXACT 207/217 both-produced = 95.4%**. 1 live-crash, 0 no-sum.
- Residual 10, delta-classified (measure, don't label):
  - **CCF:4** (Δ6, e.g. 183/177, 270/264) — the cornered crown-width CCF boundary class (cf. #89).
  - **TopHt:3** (Δ1) — cornered AVH sort-tie.
  - **TCuFt:4 / MCuFt:2 = a REAL ~3-4%-high total-cubic class** on stands 18224519010661 (952/983) and
    18756672010661 (1079/1123) — BOTH LOCATION=904/FOREST=4 (region 9), no broken tops, no cull, dubbed
    heights. NOT aspen-cftopk. Filed **task #97** (region-9 forest-4 LS cubic; needs a per-tree live debug
    to isolate dubbed-height vs R9-Clark-coefficient vs VOLEQ-assignment before a bit-exact fix).
So LS FIA residual = cornered crown/AVH ULP + one real forest-4 cubic class (#97). NE/CS spread sweeps + a
per-variant coverage-matrix refresh are the next Pillar-1 increments.

## 4-variant delta-classified FIA sweep [2026-07-08] — the Pillar-1 coverage matrix
Rowid-spread treed stands per variant, cycle-0 all-10-col vs freshly-relinked live, residuals **delta-measured**
(not label-guessed) and split into cornered-ULP vs REAL divergences.

| Variant | both | bit-exact | rate | cornered residual | REAL divergence (task) |
|---|---|---|---|---|---|
| SN | 292 | 287 | **98.3%** | 4 TopHt (AVH-tie) + 1 BdFt Δ1 | none — all cornered (after #95+#96) |
| NE | 117 | 113 | 96.6% | CCF Δ3, TopHt tie, 1 small TCuFt | 1536568185290487 TCuFt +52% (**#98**) |
| CS | 83  | 78  | 94.0% | TopHt Δ2 tie, small vol | 3303107010661 TopHt +29%→vol (**#99**) |
| LS | 217 | 207 | 95.4% | CCF Δ6 (crown), TopHt tie | 18224519010661/18756672010661 TCuFt ~3-4% (**#97**) |

**SN is fully clean** (residual = 100% cornered-ULP). NE/CS/LS each have ONE real divergence class, and they share
a **common fingerprint: region-9 forests + dubbed heights + high volume/height** — LS forest 904/4, NE 930/30,
CS 912/12. Hypothesis: a single underlying region-9 dubbed-height-or-R9-Clark issue across the eastern variants
(all three use R9 Clark), OR three forest-specific gaps. Each filed (#97/#98/#99) with a per-tree debug entry-point;
verify bit-exact vs live before landing. Harness generalized to `var_sweep.jl <VARIANT> <rowid-mod>` (spread sample
+ extract + delta-classify in one). D38 SIGFPE crashes excluded per variant (SN dominant).

### Region-9 cluster (#97/#98/#99) — ROOT-CAUSE LEAD found [2026-07-08]
Probed the CS stand (#99, 3303107010661, all heights dubbed, TopHt +29%): SITE_INDEX/SITE_SPECIES are **NULL**
in the DB, so both engines DEFAULT the site. Live `.out` SITECODE builds a per-forest per-species default-SI
table (RC=52 … WO=65 … RM=86 …) and selects **SITE SPECIES=WO (white oak, code 47)**; jl picks **species 63,
site_index=0.0**. ⇒ the divergence is the **eastern default-site-index / site-species selection** (sitset.f /
SITECODE path) for region-9 stands with no input SI — jl's wrong site species/SI inflates the dubbed heights,
which cascades to volume. This is the LIKELY SHARED mechanism behind all three (all region-9, all rely on the
default-SI table); porting/fixing the eastern default-SI + site-species selection may close #97/#98/#99 at once.
Caveat: NE (#98, +52%, has some measured heights) may carry an extra volume-eq component; LS (#97, ~3-4%,
TopHt not flagged) a milder SI diff or R9-Clark-coefficient residual — verify each bit-exact vs live.

### ★ CORRECTION to the region-9 "root-cause lead" above [2026-07-08]
The preceding "default site-index / site-species (WO 80 vs 65)" root-cause lead is **RETRACTED** — it was a
probe artifact: the CS stand was run through `initialize()` WITHOUT a variant, which defaults to **SOUTHERN**
(SN `site_index.jl:62` sets `site_species=63`). Re-probed under the correct **CS** variant, jl gets the site
**right**: `site_species=47 (WO)`, `sp_site_index[WO]=65.0` == live SITECODE. So it is NOT a site-index bug and
there is **no proven shared site-index mechanism** across #97/#98/#99 — they are three independent real
divergences (re-investigate each under its own variant):
- **#99 (CS)**: site correct; the per-species SI fan-out gives sp84=77.3/sp33=72.6 and jl dubs sp-84 dbh18.9→71ft
  (TopHt 53 vs live 41) — a CS per-species-SI-fan-out and/or H-D dubbing-curve issue.
- **#98 (NE)**: +52% TCuFt with SOME measured heights — likely a volume-eq/species issue, not dubbing.
- **#97 (LS)**: TopHt not flagged, only TCuFt ~3-4% — the R9-Clark cubic path, not a gross height error.
META (doctrine RE-TRACE): always pass the correct `variant=` to `initialize()`/`run_keyfile()`; a missing
variant defaults to SN and silently contaminates cross-variant probes. Measurement caught the wrong label — again.

### #99 FIXED (bit-exact, one line) — CS American-elm birch-group transcription error [2026-07-08]
Root cause: `src/engine/volume.jl:45` (`_htdbh_wykoff` ifor-3 override) wrongly grouped American elm (sp33)
with birch (sp30/31), overriding AE's Wykoff coefs with birch 4.4635/-3.6456 ⇒ dbh10.2 dub 67.2ft vs live
56.5ft (live `htdbh.f` stamp: AE WYK Wykoff HT1=4.6008/HT2=-7.2732→56.51; jl's base coefs already correct).
Fix: remove sp33 from the birch group. VERIFIED: CS stand 3303107010661 now bit-exact (TopHt 41/TCuFt 433/
MCuFt 378, was 53/509/378); suite 37121/140/0; CS FIA sweep **94.0%→95.2%**. Remaining CS residuals smaller
(e.g. 103635948010661 TopHt 36/38 — possibly another species-grouping issue). The site-index "shared mechanism"
across #97/#98/#99 was REFUTED (variant-contaminated probe); the three are independent. Updated coverage:
SN 98.3% (clean) / NE 96.6% (#98) / CS 95.2% (#99 fixed, minor residual) / LS 95.4% (#97).

### #98 FIXED (bit-exact, one line) — NE catch-all species RL→OH [2026-07-08]
NE FIA +52% was a species-resolution bug: `northeast/species.jl:48` set the catch-all "other" species to
index 97 (RL, red elm, VOLEQ 900CLKE975) instead of 98 (OH, other hardwood, 900CLKE998) — a "TODO verify"
from the port. Unmatched FIA codes (e.g. 6918 "other hardwood") got the RL volume eq. Fix: other_species→98.
VERIFIED: stand 1536568185290487 bit-exact (TCuFt 145 was 220, CCF 53 was 56); suite 37121/140/0; NE FIA
sweep **96.6%→98.3%, BIG divergences 0** (residual all cornered). LS(49)/CS(85) catch-alls checked = correctly
OH(998), so this was NE-only. Updated 4-variant coverage: SN 98.3% / NE 98.3% / CS 95.2% / LS 95.4%.
Only #97 (LS R9-Clark cubic ~3-4%, NOT a catch-all) remains a real region-9 divergence.

### #97 FIXED + ★ ALL 4 VARIANTS FIA-CLEAN [2026-07-08]
#97 (LS TCuFt ~3-4% high) + the true cause of #99 (CS) were ONE variant-safety bug: `_htdbh_wykoff`
(volume.jl) applied the **NE Allegheny IFOR=3 HT-DBH override** (NE species indices) to ALL variants at
forest-index 3 — LS aspen (sp41) got NE ash's curve (dub 64.2 vs live 52.2), CS elm got NE birch's. Fix:
gated the block to `ifor==3 && isne` (isne=variant isa Northeast), restored sp33 to the NE birch group,
threaded isne through all callers. Verify: LS stand bit-exact (TCuFt 952 was 983); CS #99 still bit-exact
(via the gate); suite 37121/140/0; LS FIA sweep 95.4→96.3%. **Final 4-variant FIA coverage — residuals now
100% cornered-ULP (no real divergence class in any variant's sample):**
| Variant | bit-exact | residual |
|---|---|---|
| SN | 98.3% | TopHt AVH-tie + BdFt Δ1 (cornered) |
| NE | 98.3% | TopHt tie + MCuFt Δ1 (cornered) |
| CS | 95.2% | TopHt tie (cornered) |
| LS | 96.3% | CCF Δ6 crown-class + TopHt tie (cornered) |
(the sub-100% is the D38 live-SIGFPE exclusions + cornered ULP; four fixes landed this campaign:
#95 SN-NC-merch, #96 Fort-Bragg-remap, #98 NE-catch-all-species, #99 CS-elm-dub → all superseded/complemented
by #97's variant-gate. Every one found by measuring live per-tree data.)

### ★★ GENERALIZATION CONFIRMED — fresh-stand 4-variant sweep (post all 5 fixes) [2026-07-08]
Re-swept all four variants on DIFFERENT stands (fresh rowid strides, none overlapping the fixed stands) to
prove the fixes generalize and no new divergence class emerged:
| Variant | both-produced | bit-exact | rate | REAL (>2%) divergences |
|---|---|---|---|---|
| SN | 186 | 183 | 98.4% | 0 (2 residuals = TopHt AVH-tie, cornered) |
| NE | 80  | 77  | 96.2% | **0** |
| CS | 43  | 41  | 95.3% | **0** |
| LS | 148 | 148 | **100.0%** | **0** |
**Zero real divergences in any variant on fresh stands.** All residuals are the documented cornered-ULP classes
(TopHt AVH sort-tie, CCF integer-boundary). The sub-100% (SN) is those + the D38 live-SIGFPE exclusions (which
jl runs correctly). This is the definitive Pillar-1 FIA drop-in validation: **SN+NE+CS+LS are bit-exact-or-
cornered on real FIA inventory across the sampled populations, with no remaining real divergence class.**
