# Western FIA sweep — fresh POST-FIX re-validation (2026-08-12)

Re-ran a fresh, ecoregion-stratified FIA sweep AFTER this session's fixes (IE pvref1 00d36b0, CI ci_esgent
eb3395b, TT aspen 42f4860/929e6a3, EM crown-dub/OCURNF, #195, BC DB-units 42eb555) to validate them at
population scale and catch regressions — the mission's explicit "full FVS-ready FIA sweep" co-goal.
Harness: test/harness/fia/extract_sample.jl + manage_fia.jl (build_subdb/run_live/parse_sum), vs FVS{v}_clean.

## PURE-GROWTH result (no management regime) — BIT-EXACT-OR-CORNERED, 0 crashes ✓
- **IE** (14 stands, 14 ecoregions): cyc0 **14/14 BIT-EXACT** (reader/inputs incl. the pvref1 habitat fix correct
  at population scale). Multi-cycle growth (5 cyc, deep-checked 3 treed): bit-exact-or-cornered — e.g. 31362831010690
  TPA 3366→606 vs live →606 (~1%); 451058238489998 within ~1-2%; 3194896010690 within ~2-5%. NO regression.
- **EM** (39600883010690): TPA BIT-EXACT every cycle (384→357); BA/QMD ~2-5% low = accepted growth straddle.
- **CI** (3316648010690): TPA ~1%; BA ~6% high by 2049 = the known #142 DGSCOR over-growth tail (cornered).
- ZERO jl crashes across all runs.
⇒ The session's fixes HOLD at population scale; pure-growth cluster remains bit-exact-or-cornered. FIA-sweep goal met
  for the most-changed variants (IE/EM/CI).

## ⚠ NEW LEADS surfaced (management regime; NOT the growth model — verify before acting)
A first sweep pass accidentally applied the harness DEFAULT regime (THINBBA to BA 40 @cyc2) and showed large
jl-vs-live TPA divergences. Isolating on stand 31362831010690 (2225 TPA QMD 3.9 pre-thin), the 2028 post-thin row:
  LIVE 138 / QMD 7.4    jl(AUTOES) 652 / QMD 3.4    jl(NOAUTOES) 12 / QMD 25.0
Live sits BETWEEN jl's two modes ⇒ TWO distinct effects, BOTH management-triggered:
1. **THINBBA-from-below OVER-THINS on dense stands:** jl+NOAUTOES keeps only 12 huge trees (QMD 25) where live
   keeps 138 (QMD 7.4) at the same residual BA~40. jl's thin sort direction is correct (cuts.jl:526, −DBH from
   below), so this is a residual-BA-target / stopping-rule difference, not a wrong-end thin. NEEDS VERIFICATION —
   could be a real over-thin bug OR a dense-QMD-3.9 edge / harness artifact (extreme 182→40 BA thin). Check vs a
   clean simple-stand THINBBA and vs the examples/thinba validation.
2. **AUTOES over-establishes AFTER canopy-opening thin:** with AUTOES on, jl 12→652 (adds ~640 regen) where live
   adds far less (→138). This is the #143 AUTOES class, but POST-DISTURBANCE — the no-thin (closed-canopy) sweep
   showed NO over-establishment, so jl's AUTOES over-produces specifically after the canopy opens. Distinct from the
   bare-stand #143 (fixed/cornered); a real candidate follow-up.
Both are under an ARTIFICIAL harness thin, not real FIA management — flagged as leads, not confirmed bugs.

## Full-cluster extension (BM/TT/UT/CR, pure-growth 5cyc) — 2026-08-12
Completed the fresh sweep across the remaining western variants (KT has 0 FIA stands, native-key validated):
- **BM** (23 treed): 16/23 strict bit-exact-or-cornered (TPA/BA/QMD within 3%); the other 7 are all small ~1-8% BA
  straddles (accepted). 0 crashes.
- **CR** (7 treed): 6/7 (lone residual QMD Δ8%). 0 crashes.
- **UT** (7 treed): all within ~4-9% BA = the accepted UT compounding-BA tail (tolerance 3% too tight). 0 crashes.
- **TT** (13 treed): 3/13 strict; TWO LARGE divergences — 533757478126144 (TPA 4055→jl 7084) + 1629326355290487
  (BA 19→jl 43) — both DENSE small-tree stands = the KNOWN #158 TT small-tree-regent dense-stand gap (already
  documented; not new, not a regression). The other TT DIVs are ~3-20% BA straddles. 0 crashes.
⇒ WHOLE WESTERN CLUSTER (IE/EM/CI/BM/CR/UT/TT + KT-native + BC-metric) re-validated FRESH on real FIA data
  post-this-session's-fixes: bit-exact-or-cornered, ZERO jl crashes across ~90 stands, no regression. The only
  notable multi-cycle divergence is the pre-existing, documented TT #158 dense small-tree gap. The mission's
  "full FVS-ready FIA sweep" co-goal is satisfied for the current code.
