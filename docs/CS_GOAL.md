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
- **CHUNK 1 DONE — cst01 cycle-0 all 6 stand columns BIT-EXACT** vs live FVScs (TPA=536, BA=77,
  SDI=160, CCF=169, TopHt=63, QMD=5.1; suite 5399/2, no SN/NE regression; test_cst01.jl).
  Variant infra + data + hooks landed: `CentralStates` (CS, MAXSP 96, YR 10, Zeide SDI, RNG 55329);
  `species.jl` (blkdat init), `site_index.jl` (ASITE/BSITE SITSET + CS FORKOD). CS HT-DBH reuses the
  Southern Curtis-Arney/Wykoff dub; CS crown + crown-width reuse the NE TWIGS/cwcalc paths (cwcalc.f
  byte-identical, keys on 2-char alpha). Real CS data extracted via tools/fortran_data_extract.py.
- **NEXT — chunk 2: CS volume.** Wire CS NVEL equation ids into the existing R9 Clark + R9LOGS path
  (src/engine/r9clark_vol.jl + volume_equations.jl) so the cyc0 .sum volume columns (Tcuft 1517 /
  Mcuft 1300 / Scuft 497 / Bdft 2903) + forest-type (503) come in. Then chunk 3: `cs/dgf.f`.
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
