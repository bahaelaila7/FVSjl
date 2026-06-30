# ACTIVE GOAL — FVSjl Non-ULP Divergence Fix Campaign

## Mission
Drive **every** place FVSjl differs from the **live** Fortran FVS (SN/NE/CS) that is *neither*
Float32 ULP *nor* the accepted COMPRESS eigensolver to **ULP-class**, OR prove it irreducible and
document exactly why. The end state: SN, NE, CS are faithful bit-exact drop-ins barring only ULP +
eigensolver. Ledger + inventory: **`docs/DIVERGENCE_FIX_CAMPAIGN.md`**.

## Oracle & tools (use these every item)
- **Live binaries** are the only oracle: `test/harness/{sn,ne,cs}_oracle.sh <key> <outdir>` relink
  `/tmp/FVS{sn,ne,cs}_new` from `bin/FVS{sn,ne,cs}_buildDir/*.o` + the glibc shim. NOT Oracle A / FVSjulia
  (it is a thing to FIX not match — proven wrong on DG-COR & s12_phys by ~7% / ~1 TPA).
- **Debug-stamp** the live `.f`: edit in the buildDir, `gfortran -c -I. <f>.f -o <f>.o`, `rm -f
  /tmp/FVS*_new`, run via the oracle, then RESTORE the `.f` + recompile clean. (Used to nail live COR.)
- **Discovery / regression sweep:** `julia --project=. test/harness/divergence_sweep.jl {sn|ne|cs} [keys…]`
  — runs many stands through live+jl, aligns by (stand,year), ranks by max non-ULP rel diff. This is the
  "FIA-plots" backbone: push many diverse plots through BOTH engines and rank where they differ.
- Canonical runner: `julia --project=. test/runtests.jl` (suite baseline **5970 pass / 2 broken** = the
  accepted SN COMPRESS eigensolver + NOHTDREG ULP; NOT campaign targets).

## Current state (continue from here — DIVERGENCE_FIX_CAMPAIGN.md is the live ledger)
- **Sweeps run.** SN = 260 stands (the big inventory). NE all-species + net01: all-species BIT-EXACT
  incl. volume; only net01 BARE-regen ~4% Mcuft late. CS all-species + cst01: at the documented ULP
  floor (cs_allsp 1.52% TPA late) + cst01 deep-thinned tail (accepted). ⇒ **the real work is SN-heavy**,
  with a regen-volume residual also visible in NE.
- **D1 (LP-calibration tail) — ✅ DISPROVEN** (was MY measurement artifact: a probe loop omitting the
  per-cycle `compute_forest_type!`, which feeds diameter growth; `run_keyfile` is bit-exact vs live).
- **Open targets (upstream→down):**
  - **D7** per-species merch/saw/board volume (all_GA/PC/BY cyc0: Tcuft bit-exact but Mcuft/Scuft/Bdft
    ~28% off ⇒ a merchandising-standard, top-dia/min-DBH, gap). DETERMINISTIC at cyc0 → start here.
  - **D2** FINT≠5 calibration volume (~0.4% cuft; FINT=5 bit-exact).
  - **D8** multiplier keywords (REGDMULT/MORTMULT/REGHMULT/BAIMULT) — large diffs.
  - **D10** regen small-tree volume (bare_* Scuft ~50%; net01 BARE Mcuft 4%).
  - **D9** mid-cycle SIMFIRE timing (s10_fire/fire_repeat huge TPA — verify real vs structure).
  - **D4** crown-biomass FMCROWE carbon ~0.9 ton · **D5** #28 snag-fall-timing ~0.2-0.4 ton (downstream).
  - **D6** CS ESCPRS regen-compression (feature gap).
  - **triage:** carbon_* `Scuft=0@2005` (probably a .sum/Volume-keyword artifact — the carbon REPORT is
    validated bit-exact; confirm not a real model diff), compress 50% (recheck vs accepted ~1%).

## THE DOCTRINE (imprinted — apply to every verification and fix)
1. **TRACE LOGIC, NOT RUNTIME.** Trace the logic in BOTH FVSjl and the FVS Fortran; match SEMANTICS.
2. **UPSTREAM-FIRST.** Work the most-upstream, least-dependent divergence first (a growth/volume fix
   removes many downstream report diffs).
3. **REGRESSION = MASKED-BUG SIGNAL.** A semantically-certain fix that regresses other tests likely
   unmasked a hidden bug — NOTE it, keep going, circle back. Do NOT revert a faithful fix to keep a
   stale test/golden green (a `.sum.save` or Oracle-A golden can be WRONG — re-ground vs the live binary).
4. **VALIDATE vs LIVE FVS FIRST + double-check, THEN tighten/write the test.** Never trust Oracle A.
   A "divergence" can be a MEASUREMENT artifact — reproduce it through `run_keyfile` (the production
   path) before believing a probe-loop number (the D1 lesson).
5. **DOCUMENT EVERY verdict** (DIVERGENCE_FIX_CAMPAIGN.md): both-sides logic, the live debug-stamp
   evidence, and the fix OR the irreducibility proof.
6. **VARIANT-AWARE — DO NOT HARDEN.** Gate variant behaviour on the variant/coefficients in shared
   `src/engine|io|core/*`; SN/NE/CS specifics under `src/variants/*`. A fix must keep all THREE variants
   bit-exact — don't harden base code to one variant.

## Re-trace discipline
"documented residual" / "known tail" / "accepted" can be a misread — re-verify against the live binary
(current, freshly relinked) before assuming it's irreducible. Loose tolerances and Oracle-A goldens have
repeatedly HIDDEN real divergences; the sweep + live debug-stamp are how we trace each floor to ground.

## Off-switch
This reminder fires from a Stop hook until the campaign is complete (every ledger item ✅ or 📌 with a
documented reason). To silence it: `touch /workspace/FVSjl/docs/DIVERGENCE_COMPLETE`.
