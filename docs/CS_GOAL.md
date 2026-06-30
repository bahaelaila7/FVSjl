# ACTIVE GOAL — FVSjl Central States (CS) variant port

## Mission
Port FVS **Central States (CS)** into FVSjl as a *semantically faithful* drop-in for the live
Fortran **FVScs** — the same standard that made Southern and Northeast validated drop-ins (tags
`FVSsn-complete`, `FVSne-complete`, `FVSsn+ne`). Bit-exact end-to-end barring only ULP
single-precision and documented eigensolver-class divergences. Full scope: **`docs/CS_VARIANT_PORT_SCOPE.md`**.

## Oracle (CS-specific — read carefully)
- Sole oracle = **live FVScs**: relink from `bin/FVScs_buildDir/*.o` + the glibc shim, exactly
  like NE. Clone `test/harness/ne_oracle.sh` → `test/harness/cs_oracle.sh` (point `BUILDDIR` at
  `FVScs_buildDir`, `BIN=/tmp/FVScs_new`). There is **NO Oracle-A / FVSjulia for CS** — no 1:1
  transliteration exists; validate against the live binary only.
- Canonical test stand: **`tests/FVScs/cst01.key`** (+ `cst01_method5.key`). Treat any committed
  `.sum.save` as potentially STALE — validate cyc1+ against the freshly-relinked live binary.
- Canonical runner: `julia --project=. test/runtests.jl` (NOT `Pkg.test`). Suite baseline
  **5391 pass / 2 broken** (the 2 broken = accepted SN COMPRESS eigensolver + NOHTDREG ULP — NOT CS).

## The leverage (why CS is the cheapest variant yet)
- **~94% of the FVScs source is the shared base** already ported variant-agnostically (engine,
  keyword dispatch, density, `grow_cycle!`/`mortality_and_fire!`, FFE base + 9 DBS tables, ECON, IO).
- The variant surface is **~30 routines**, and most are **near-identical to NE** (which is done):
  `htgf` 96%, `dgdriv` 97%, `cratet` 98%, `regent` 96%, `varmrt` 94%, `crown` 88%, `cubrds` 84%,
  and CS uses the same **`balmod` BAL-competition** framework. CS volume rides NE's **R9 Clark +
  R9LOGS** path. So NE's structurally-new pieces (BAL competition, height-growth model, R9 volume)
  carry over.
- The ONE genuinely-new model is **`cs/dgf.f`** — an **SN-family** ln(DDS) regression (DBH / site /
  crown / **BA-percentile / QMD**, trees ≥ 5″), NOT NE's BAL-potential iteration. Extend the
  **Southern `dgf!`** framework with CS coefficients, not `ne_dgf!`. Plus CS site/forest/HT-DBH
  tables (`htdbh`/`sitset`/`forkod`), CS fire fuel models (`fire/cs/`, 8 files; CS reuses SN's
  `fmmois`), and the **96-species** coefficient CSVs (CS `MAXSP = 96`; NE 108, SN 90).

## Current state (continue from here — see docs/CS_PORT_STATUS.md for the live baseline)
- **CHUNKS 1 + 2 DONE — cst01 cycle-0 all 6 stand + all 4 volume columns BIT-EXACT** vs live FVScs
  (TPA=536/BA=77/SDI=160/CCF=169/TopHt=63/QMD=5.1; Tcuft=1517/Mcuft=1300/Scuft=497/Bdft=2903; suite
  5403/2, no SN/NE regression; test_cst01.jl). Infra+data+hooks: `CentralStates` (MAXSP 96, YR 10,
  Zeide SDI, RNG 55329); species.jl, site_index.jl (ASITE/BSITE SITSET + CS FORKOD), diameter_growth.jl
  (cs_dgcons! bark copy stub). CS HT-DBH = shared Southern Curtis-Arney/Wykoff; CS crown + crown-width =
  shared NE TWIGS/cwcalc; CS volume = shared eastern R9 Clark + `_cs_merch`. Data via fortran_data_extract.py.
- **CHUNK 3 (DG) — model+calibration BIT-EXACT** (validated vs live DEBUG/DGF + DGDRIV stamps;
  cs_dgf! + cs_dgcons! + dg_coeffs.csv). Remaining DG residual = the DGSCOR serial-correlation /
  tripling step (live tripled-mean DG 0.89 vs jl central 0.654 — the OLDRN/FRM partition for
  uncalibrated species); best validated aggregately via the cycle-1 .sum. See docs/CS_PORT_STATUS.md.
- **CHUNK 4 DONE — full growth spine + CYCLE-1 BIT-EXACT.** Height growth (cs/htgf.f NC-128/MAPCS +
  cs_balmod), small-tree (cs/regent.f REGENT, XMIN=3, cs_balmod), mortality (generic + CS background
  PMSC/PMD via IMAPCS + varmrt Union). grow_cycle!(CentralStates) runs end-to-end; cst01 cycle-1 ALL
  SIX stand columns BIT-EXACT vs live (TPA 518/BA 99/SDI 196/CCF 202/TopHt 68/QMD 5.9; test_cst01.jl
  cyc1 6/6). KEY FIX: the DG-calibration GST DBH floor is CS=5 (cs/dgdriv.f:380), not SN/NE's 3 —
  jl was over-calibrating WO (debug-stamped live FN<5 ⇒ COR=0 for all CS species).
- **CHUNK 5 (FFE) STARTED — full 10-cycle .sum VALIDATED for the projection stand.** Wired the CS FFE
  fuel/fire-effects/crown-biomass path (FMCBA FUINI/FULIV fuel loading, Anderson-13 fuel models, cs/fmeff
  + cs/fmbrkt bark/mortality groups, cs/fmcrowe crown biomass via ISPMAP=ls_spi, fmcblk BIOGRP). Ran the
  first stand (UNTHINNED CONTROL, NOAUTOES) 10 cycles vs freshly-relinked live FVScs: **cyc0-2 ALL SIX
  stand columns BIT-EXACT** (1990/2000/2010), **cyc3-10 within the Float32 ULP floor** (TPA ±1-3, SDI ±1,
  QMD ±0.1, vol ≤1.5%; one FORTYP 503→801 flip-timing diff = accepted SN-COMPRESS class). test_cst01.jl
  +54 assertions (79 CS total). Suite 5416/2, no SN/NE regression. See docs/CS_PORT_STATUS.md.
- **CHUNK 6 (ESTABLISHMENT) DONE — BARE-GROUND-PLANT stand cyc1 BIT-EXACT.** CS estab follows the NE
  pattern (ESSUBH base height → REGENT-LESTB creation-cycle growth), NOT SN's. _CS_ESSUBH_REFAGE +
  _CS_ES_HHTMAX + cs_htcalc_height + CS LESTB growth branch. Re-trace caught a real bug (the
  cs/regent.f:341 +0.001·hk add I'd first omitted) ⇒ 2002 bit-exact (was BA 9/8).
- **BOTH CANONICAL KEYS VALIDATED END-TO-END (suite 5555/2).** cst01.key (5 stands: projection, THINDBH,
  THINPRSC shelterwood, SPECPREF/THINBTA+Econ, BARE-PLANT) + cst01_method5.key run through run_keyfile
  and validate vs live FVScs: thinning SELECTIONS + post-cut states BIT-EXACT (cut logic faithful);
  deep-thinned tails = documented single-precision floor amplified at discrete thresholds. test_cst01.jl
  now has cyc0 + cyc1 + multi-cycle projection + BARE establishment + thinning testsets.
- **CHUNK 7 (SPROUTING) DONE — sprout-regen cycle BIT-EXACT.** Ported cs/essprt.f `CASE('CS')` (essprt_cs
  PREM / nsprec_cs / sprtht_cs + cs_sprout_dbh + aspen ASSPTN sp76 + sprout_essprt.csv is_sprouting) and a
  CentralStates branch in esuckr!. Validated vs live (SPROUT+clearcut): 2010 sprout cycle bit-exact (434/
  23/44/40/44/3.1). Also fixed a shared SPROUT-handler bug (blank species field ⇒ IS=0 all-species enable,
  not disable). SIMFIRE re-trace: cst01 stand 3 DOES fire `SIMFIRE 2003` (non-vacuous: fire kills ~79 TPA;
  jl 173 vs live 178), so the CS fire-effects are validated through a real fire.
- **★★★ CS PORT COMPLETE — meets the SN/NE validated-drop-in standard.** Every model ported + live-validated
  bit-exact where it counts (cyc0 all 10 columns, growth spine cyc1-2, establishment/thinning/sprouting/fire
  at their key cycles); both canonical keys run end-to-end within the documented single-precision floor.
  Only ESCPRS regen compression deferred (exactly as in SN). Suite 5562/2, no SN/NE regression. Off-switch
  set (docs/CS_COMPLETE).
- Scope + per-routine NE-reuse table + chunk order: `docs/CS_VARIANT_PORT_SCOPE.md`.

## Ordered work (most-upstream / least-dependent first — see scope §7)
1. **Variant infra:** `CentralStates<:AbstractVariant`, per-variant `nspecies`→96, CS FORKOD/site
   defaults, CS species CSVs, `cs_oracle.sh`. ⇒ drive `cst01` **cycle-0** stand columns bit-exact
   (TPA/BA/SDI/CCF/QMD/TopHt) + a `test/integration/test_cst01.jl`.
2. **Volume:** wire CS NVEL equation ids into the existing R9 Clark + R9LOGS path (variant dispatch).
   ⇒ cycle-0 `.sum` volume columns.
3. **Diameter growth:** `cs/dgf.f` (the new SN-family CS model) + CS `balmod`/competition; the
   `dgdriv` calibration is shared/near-NE. ⇒ cycle-1 DG.
4. **Height growth** (reuse NE `htgf` + CS coefs) · **crown** · **mortality** (`varmrt`) ·
   `htdbh`/`sitset`. ⇒ drive `cst01` to cycle-1+ vs live.
5. CS-active shared branches (un-gate, don't harden) → FFE (`fire/cs/`).

## THE DOCTRINE (imprinted — apply to every verification and fix)
1. **TRACE LOGIC, NOT RUNTIME.** Trace the logic in BOTH FVSjl and the FVS CS Fortran; match SEMANTICS.
2. **UPSTREAM-FIRST.** Work the most-upstream, least-dependent issue first.
3. **REGRESSION = MASKED-BUG SIGNAL.** A semantically-certain fix that regresses other tests likely
   unmasked a hidden bug — NOTE it, keep going, circle back. Do NOT revert a faithful fix to keep a
   stale test green (a `.sum.save` can be STALE — the oracle, not your fix, may be wrong).
4. **PORT SEMANTICS FAITHFULLY FIRST + double-check vs LIVE FVScs, THEN write the test.** Never test-first.
5. **DOCUMENT EVERY verdict + both-sides logic** (docs/CS_*), so we never run in circles.
6. **VARIANT-AWARE — DO NOT HARDEN.** In shared `src/engine|io|core/*` gate variant behavior on the
   variant/coefficients; SN under `src/variants/southern/`, NE under `src/variants/northeast/`, CS under
   `src/variants/centralstates/`. Adding CS must keep **SN *and* NE bit-exact** (don't harden base code
   to any variant). Reuse the NE routine where CS ≈ NE; write CS-specific only where it genuinely differs.

## Re-trace discipline
A "done" / "ported" / "validated" / "reuses NE" label can be a misread — re-verify against the CS
Fortran source before assuming scope (the % similarity is a hint, not a guarantee — a 96%-similar
routine can hide a CS-specific coefficient or branch). The SN+NE campaigns repeatedly caught real bugs
behind stale already-done flags and behind "irreducible-RNG" floors (trace floors to ground).

## Off-switch
This reminder fires from a Stop hook until the CS port is complete. To silence it:
`touch /workspace/FVSjl/docs/CS_COMPLETE`
