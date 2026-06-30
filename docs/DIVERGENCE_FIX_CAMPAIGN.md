# FVSjl — Non-ULP Divergence Fix Campaign (2026)

Drive every live-FVS divergence that is **neither Float32 ULP nor the accepted COMPRESS
eigensolver** to ULP-class, or prove it irreducible and document why. Oracle = live FVS
(sn/ne/cs_oracle.sh; debug-stamp the .f, relink, restore). Doctrine: trace logic both
sides → upstream-first → validate vs LIVE before writing the test → document every
verdict → variant-aware (gate, don't harden; keep all three variants bit-exact).

Status: ⬜ open · 🔬 investigating · ✅ fixed-to-ULP · 📌 irreducible/deferred (why documented)

| # | Divergence | Layer (upstream→down) | Magnitude vs live | Status |
|---|---|---|---|---|
| D1 | ~~LP-growth-calibration tail~~ | growth | — | ✅ NOT REAL (artifact) |
| D2 | FINT≠5 calibration volume | growth | ~0.4% cuft | 🔬 localized (deferred, low-impact) |
| D3 | Multi-point density (PCCF/TCONDMLT/structure-stage) | density | per-point approx | ⬜ |
| D4 | Crown-biomass FMCROWE carbon residual | carbon report | ~0.9 ton AGL | ⬜ |
| D5 | #28 carbon snag-fall-timing residual | carbon report | ~0.2-0.4 ton | ⬜ |
| D6 | CS ESCPRS regen-compression not ported | regen | feature gap | ⬜ |
| D7 | Per-species merch/saw/board volume (GA/PC/BY) | volume | cyc0 ~28% Bdft | ✅ FIXED-to-bit-exact |
| D8 | Multiplier keywords (REGDMULT/MORTMULT/REGHMULT/BAIMULT) | growth | large | ⬜ NEW |
| D9 | Mid-cycle SIMFIRE timing (s10_fire, fire_repeat) | fire | TPA huge | ⬜ NEW (verify) |
| D10 | bare_* regen volume (Scuft) | regen | ~50% Scuft | ⬜ NEW |

## Discovery tool — `test/harness/divergence_sweep.jl`
The campaign's plot-based differential (the user's "FIA-plots" principle). Runs many stands through the
live binary ({sn,ne,cs}_oracle.sh) + jl `run_keyfile`, aligns by (stand, year), and ranks scenarios by
max NON-ULP relative diff (skips ≤1 print unit AND ≤0.2%). `julia --project=. test/harness/
divergence_sweep.jl sn`. SN run = 260 stands; the live-vs-jl inventory below is its output.

### SN sweep inventory (2026, ranked) — triaged
- **Real, cycle-0 (deterministic) → D7:** all_PC/GA/BY/GA Bdft@1990 10-35% — Tcuft bit-exact but
  Merch/Saw/Board off ⇒ per-species merchandising standard (top-dia / min-DBH) wrong for these species.
  (all-species test gap: it asserts stand cols but NOT volume — extend it.)
- **Real, growth → D2/D8:** growth_fint10 1.24% (FINT), timeint10 1.96% (non-native cycle), mult_*
  (REGDMULT/MORTMULT/REGHMULT/BAIMULT) large — multiplier-keyword application.
- **Real, regen → D10:** bare_natural/plant/multipoint/mp3 Scuft ~50% — regen small-tree volume.
- **Fire — verify D9:** s10_fire 789% / fire_repeat 288% TPA (mid-cycle SIMFIRE timing?); fire_burn/early
  4.38% Bdft (documented post-fire DG residual); fuelmodl/defulmod/salvage few-%.
- **Carbon scenarios:** carbon_* Scuft jl=0.0 @2005 — likely a .sum-structure/Volume-keyword artifact
  (the CARBON REPORT itself is validated bit-exact); verify not a real model diff.
- **Known/accepted:** compress (s22 eigensolver — but 50% needs a recheck vs the accepted ~1%),
  treeszcp_cap/htcap (declining-stand), dense_long/s09_cyc20 0.76% (long-run ULP).

### D7 — per-species merch/saw/board volume — 🔬 NARROWED to the R9 Clark merch EXTRACTION
all_GA (homogeneous green ash) cyc0: TPA/BA/SDI/**Tcuft BIT-EXACT**, but Mcuft live 900/jl 977, Scuft 47/60,
Bdft 174/223 (~28%, jl HIGHER). Ruled out:
- merch STANDARDS: GA(37) has the SAME top_dib=4/dbh_min=4/scf_top_dib=9/scf_min_dbh=12/bf_top_dib=9 as
  the bit-exact snt01 species (HI 27, SO 64) ⇒ NOT a standards-data gap.
- gross Clark equation: GA uses its own Clark eq `CLKE544` (FIA 544); Total cubic is bit-exact ⇒ the
  profile coefficients are right for TOTAL volume.
⇒ The divergence is in the **R9 Clark merch/saw EXTRACTION** — the DIB (diameter-inside-bark) profile
integrated from stump to the merch-top-diameter height (vol[4]+vol[7], r9clarkdib.f). jl over-extracts
merch (higher Mcuft/Scuft) for Clark eq 544 (and the PC/BY eqs) while matching total. NEXT: debug-stamp
live r9clark/r9clarkdib for a single GA tree (dump DIB-at-height + the merch-cut height + vol[4]/vol[7])
vs jl's `compute_volumes!` for the same tree; the merch-cut height or a profile-segment term differs for
this Clark-equation family. (Note: this is volume-extraction, downstream of growth — but a real cyc0
divergence, so high-value: deterministic, no RNG/timing confound.)

**✅ FIXED (bit-exact).** Root cause = `COEFFSO%DIB17` (the secondary-coefficient inside-bark diameter at
17.3 ft). Live r8prep.f gates the whole fcmin block on `IF(SPEC.NE.221.AND..NE.222.AND..NE.544)`: for
those three species the `(FCLSS−AFI)/BFI` line (r8prep.f:366) is SKIPPED, COEFFSO%DIB17 stays 0, and the
unconditional `:507` floor `IF(COEFFSO%DIB17 < COEFFS%DIB17) COEFFSO%DIB17 = COEFFS%DIB17` then sets it =
COEFFS%DIB17 (= the raw dib17). jl's `_r8_clark` computed `dob17 = (dib17−AFI)/BFI` for ALL species
(missing both the special-case and the :507 floor) ⇒ a too-large dob17 (BFI<1) ⇒ over-extracted
merch/saw/board. Fix (r8clark_vol.jl): `dob17 = (spec∈221/222/544) ? dib17 : (dib17−AFI)/BFI; dob17 =
max(dob17, dib17)`. The :507 floor is a no-op for every other species (proven: all_WO/LP + snt01 stands
1-4 stay bit-exact) and yields dib17 for the three. all_GA/PC/BY cyc0 now BIT-EXACT (1253/900/47/174 ==
live). Suite 6234/2. (snt01 stand-5 BARE residual that remains = D10 regen volume, separate.)

## Verdict log

### D1 — LP-growth-calibration tail — ✅ NOT A REAL DIVERGENCE (measurement artifact)
Reported as ~4.8 TPA / 0.8″ QMD on mix_lp_hi. **Disproven**: `run_keyfile` on mix_lp_hi is BIT-EXACT vs
live FVSsn every cycle (only 1995 Tcuft Δ1 = ULP). The "drift" came from my tolerance-probe loop OMITTING
the per-cycle `compute_forest_type!` — FORTYP (520, ported) feeds diameter growth, so a stale forest type
shifted DG. With FORTYP recomputed each cycle (as the real test does), all 10 multicycle scenarios match
live to print-rounding (TPA ≤0.57, cuft ≤1.0). multicycle re-tightened to uniform atol=1 (bogus mix_lp_hi
carve-out removed). LESSON: re-trace a "tail" through the actual production path before believing a probe.

### D2 — FINT≠5 calibration volume residual — 🔬 REAL, ~0.4% cuft
growth_fint10 (GROWTH diameter-measurement FINT=10, SCALE=YR/FINT=0.5, dgdriv.f:325): TPA/SDI/TopHt
bit-exact, BA ±1, but Tcuft 1995 live 2848 / jl 2835 (Δ13, 0.46%), 2000 live 3308 / jl 3295. Committed
growth_fint10.sum.save MATCHES live ⇒ not stale; genuine. growth_idg1 (FINT=5) is fully bit-exact ⇒
FINT-specific.

**Localized (live debug-stamp).** growth_fint10 = 6 loblolly (sp 13), measured DG=1.5"/10yr. Per-tree:
jl central tree-1 DBH 8.824 vs live 8.9 (jl grows ~0.08" LESS). Stamped live dgdriv COR(13)=0.547359 vs
jl 0.552651 (~1% high). The calibration term flow is BIT-EXACT-identical both sides: live dgdriv.f:423
TERM=DG*(2*BARK*WK3+DG)*SCALE == jl term=dg*(2*bark*wk3+dg)*scale, then RESLOG=log(TERM)-WK2 ->
DEV/DEVSQ/SNX/SNY/SNXX/SNXY -> the regcor/WC*cornew formula (jl matches dgdriv:520-590). Both use
SCALE=0.5. SCALE only shifts cornew by log(scale); WC depends nonlinearly on cornew^2 (temp). Since the
FORMULA matches and FINT=5 is bit-exact, the residual ~1% must be a SCALE-DEPENDENT INPUT differing: most
likely the DENSE density-backdating (BAL/PCCF at start-of-period, feeding WK2=DGF) using a 5-yr vs the
FINT=10 period. NEXT: debug-stamp live WK2/EDDS + DEV/DEVSQ/SNY for sp13 vs jl (one stamp = input-vs-
formula). LOW IMPACT (0.4% volume, non-default GROWTH FINT!=5; default FINT=5 bit-exact) — characterized,
deferred behind higher-impact items.
