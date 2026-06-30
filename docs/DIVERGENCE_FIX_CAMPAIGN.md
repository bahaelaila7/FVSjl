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
| D2 | FINT≠5 / non-native cycle-length calibration (volume) | growth | ~0.4% cuft | 🔬 |
| D3 | Multi-point density (PCCF/TCONDMLT/structure-stage) | density | per-point approx | ⬜ |
| D4 | Crown-biomass FMCROWE carbon residual | carbon report | ~0.9 ton AGL | ⬜ |
| D5 | #28 carbon snag-fall-timing residual | carbon report | ~0.2-0.4 ton | ⬜ |
| D6 | CS ESCPRS regen-compression not ported | regen | feature gap | ⬜ |
| D7 | Per-species merch/saw/board volume standards (GA/PC/BY…) | volume | cyc0 ~28% Bdft | 🔬 NEW |
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

### D7 — per-species merch/saw/board volume standard — 🔬 evidence
all_GA (homogeneous green ash) cyc0: TPA/BA/SDI/Tcuft BIT-EXACT, but Mcuft live 900/jl 977, Scuft 47/60,
Bdft 174/223 (~28%). Total cubic exact ⇒ the merchandising cut (merch top diameter / min DBH / stump)
for GA (and PC, BY) differs from live. Next: compare jl `init_merch_standards!` per-species topd/dbhmin
vs the SN volume-standard source for these species.

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
FINT-specific. Next: debug-stamp the live dgdriv FINT-scale path; BA-exact-but-cuft-off ⇒ sub-integer DBH
diff, likely the DDS→DBH back-calc under SCALE=0.5.
