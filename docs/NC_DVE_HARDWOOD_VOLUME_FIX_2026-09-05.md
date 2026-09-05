# NC (Klamath) DVE California-hardwood VOLUME fix — 2026-09-05

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Engine base: master d415384f (NC Siskiyou R6 vol + redwood broken-top + hardwood DG fixes). This fixes a
SEPARATE, VOLUME-only crosswalk-exposed bug the master-engine NC re-sweep exposed at cursor 0→4000 (the one
`volume_persistent` dig entry 23950554010900 + a large fraction of the BdFt-worst caps).

## THE BUG (real, one-directional, deterministic, cycle-0)
NC's **DVE California-hardwood** volume path (r5harv.f — species MA/BO/TO/OH, and every hardwood crosswalked
into those 4 slots: canyon/coast live oak, chinkapin, etc.) had TWO un-oracle-validated gaps, both present at
CYCLE 0 (pure equation, no growth), structure BIT-EXACT:
1. **Total/merch cubic OVER by up to +44% on BROKEN-TOP hardwoods.** jl computed r5harv with the dubbed
   unbroken NORMHT (correct — the oracle does too, V1=33.55 matches) but then SKIPPED the CFTOPK broken-top
   trim. The jl comment "the DVE California-hardwood path skips CFTOPK (fvsvol.f)" was WRONG: **vols.f:193
   applies CFTOPK for METHC=6 whenever `CTKFLG .AND. TKILL .AND. VMAX>0`, and BOTH the WO2W conifers AND the
   DVE hardwoods are NVEL method 6** (fvsvol.f NATCRS sets VMAX=TCF & CTKFLG=.TRUE., line 510/531). So a
   broken-top hardwood's full-NORMHT cubic was never trimmed back to the standing break.
2. **Board feet DEFERRED (jl BdFt=0 vs oracle >0).** The r5harv Scribner branch (VOL(2), r5harv.f:342-408)
   was never ported — `t.bdft_vol[i]=0f0  # DVE board deferred (no hardwood on nct01)`.
Also the merch cubic used only CV6 (VOL(4)) instead of MCF = VOL(4)+VOL(7) = cuftgros+topwood (fvsvol.f:512).

Isolated via instrumented FVSnc_g16 (dvest.f R5HARVDBG print) on stand 23950554010900 (pure BO-slot,
coast-live-oak 801): oracle passes H=NORMHT≈45 → V1=33.55 (== jl), then trims to 21.7 in the .sum; jl stayed
at 33.5. Broken-top trees (TRC HT>0) over-predicted; sound trees (TRC HT=0) were already bit-exact.

Invisible on the nct01 conifer validation set (100% conifers, WO2W path, no hardwood) ⇒ classic
crosswalk-exposed unvalidated path, same taxonomy as the redwood broken-top (cc24e470) and Siskiyou-R6
(d415384f) fixes.

## THE FIX (src/variants/klamath/volume.jl)
- `nc_r5harv_vol` now returns `(tcuft VOL1, MERCH = VOL4+VOL7, scribner-BF VOL2)` — added the faithful
  r5harv.f misc-hardwood Scribner board branch (D<11 ⇒ cubic×4 ratio; D≥11 ⇒ TARIF→RS616/SV616(5"top) /
  SV816(7"top)), and MCF = cuftgros+topwood.
- `compute_volumes_nc!` DVE branch computes the FULL (NORMHT) cubic/merch/board, then applies
  `r4_topkill` (CFTOPK/BFTOPK, same call the WO2W path uses) for broken tops. DVE is NO LONGER skipped.

## VALIDATION (live FVSnc_clean A/B, cyc0)
Stand 23950554010900 (BO/coast-live-oak): TOT 1289==1289, MERCH 939==939, BDFT 3321 vs 3316 (−0.15%) —
was TOT 1861 (+44%), BDFT 0. cyc0 TCuFt now BIT-EXACT across every hardwood species:
- canyon live oak 1123837615290487: TOT 4791==4791 (MCuFt −0.6%, BdFt −0.6%)
- tanoak         1123828024290487: TOT 4785==4785 (MCuFt −0.02%, BdFt −0.06%)
- other-hardwood 1123829821290487: TOT 4384==4384 (dense 3300-TPA MIXED DF+MA+OH: MCuFt −2.2%, BdFt −4.7%)
- madrone        1123831101290487: TOT 4636==4636 (MCuFt −0.1%, BdFt +0.07%)
Residual MCuFt/BdFt = small Float32 accumulation in the Scribner ratio + bftopk trim over dense stems
(secondary merch metric; primary TCuFt + all COUNT metrics bit-exact). No conifer regression (redwood
302012523489998 53625==53625, Siskiyou 301903713489998 8057==8057 unchanged).
Gate: `test/integration/test_multicycle.jl` **339 pass / 11 broken** (default bounds), byte-identical.
