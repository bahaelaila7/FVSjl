# FVSsn → FVSjl faithfulness audit (operation-sequence trace)

Re-verifying the port by tracing FVS source op-by-op (NOT trusting tests, especially
FVSjulia-derived ones). Methodology: read the FVS call sequence, map each operation to
its FVSjl counterpart, verify the order + the stand state each op reads (pre/post thin,
pre/post growth). Where a mismatch is suspected, design a scenario and trace the values.

## Per-cycle driver: FVS tregro.f = GRINCR (compute increments) → GRADD (apply + bugs)

FVS splits the cycle: GRINCR computes DG/HTG/mortality (into increment arrays, NOT applied)
then GRADD applies them (UPDATE) plus fire/establishment/crown. FVSjl fuses both into
`grow_cycle!` (compute then apply inline). The canonical FVS order and the FVSjl mapping:

| # | FVS (GRINCR then GRADD)                          | FVSjl grow_cycle!                       | status |
|---|--------------------------------------------------|------------------------------------------|--------|
| 1 | RCON (mid-run site change)                       | apply_setsite! (200)                     | ✓ |
| 2 | SDICAL/SDICLS → SDIBC (pre-cut, Reineke STAGE)   | (implicit; density at 199)               | ✓ |
| 3 | EVMON(1) / ECSTATUS(0)                            | event monitor / econ                      | ✓ |
| 4 | **CUTS** (thin)                                  | cuts! (205)                               | ✓ |
| 5 | DENSE (post-thin)                                | compute_density! (206)                    | ✓ |
| 6 | SDICAL/SDICLS → **SDIAC** (post-thin, Reineke)   | stand_sdi_reineke (228, threaded→crown)   | ✓ (fixed this session) |
| 7 | COMCUP (COMPRESS)                                | apply_compress! (201) — **BEFORE cuts**   | ⚠ see F2 |
| 8 | DGDRIV (compute DG, pre-growth dbh)              | diameter_growth! (229)                    | ✓ |
| 9 | HTGF (compute HTG)                               | height_growth! (230)                      | ✓ |
| 10| REGENT (small-tree DG/HTG)                       | small_tree_growth! (231)                  | ✓ |
| 11| FIXDG/FIXHTG opts                                | apply_fix_scalers! (232/234)              | ✓ |
| 12| **MORTS** (compute kill, pre-growth + pre-fire)  | mortality! (in MORTS/fire block)          | ✓ (F1 fixed) |
| 13| TRIPLE + REASS (if LTRIP)                         | triple_records! (243)                     | ✓ |
| 14| GRADD: FMMAIN/FMKILL (fire) = MAX(MORTS,fire)    | _maybe_burn! after MORTS, MAX-combine     | ✓ (F1) / F3 distrib |
| 15| HTGSTP (HTGSTOP/topkill)                          | htgstp! (245)                             | ✓ |
| 16| **UPDATE** (apply DG/HTG + kill to dbh/ht)       | inline dbh+=DG/bark, ht+=HTG (251-261)    | ✓ |
| 17| DENSE (post-growth)                              | compute_volumes!/density                  | ✓ |
| 18| ESNUTR (sprout + establishment)                  | esuckr! (274) + establish! (275)          | ✓ |
| 19| DENSE (post-establishment)                       | (establish! recomputes density)           | ✓ |
| 20| **CROWN** (crown ratio, post-growth)             | crown_ratio_update! (276)                 | ✓ |
| 21| CWIDTH                                            | (crown width on demand)                   | ✓ |

VERIFIED-FAITHFUL invariants (traced this session): DG/HTG/mortality are computed on the
PRE-growth dbh (increments stored, dbh updated only at step 16); MORTS runs before TRIPLE;
TRIPLE after MORTS; establishment after growth and before CROWN; the crown RELSDI reads the
post-thin pre-growth Reineke SDIAC (step 6 — was a bug, fixed). The mortality trajectory G,
the BAMAX/size-cap G, the SDI sums, DGSCOR/BACHLO RNG order, the species-sort IND1 order, and
the wk2/frmbase/COR/AUTCOR growth internals were all independently traced to live FVSsn.

## Findings

### F1 — Fire vs regular MORTS ordering (FFE). FIXED (faithful MAX-combine per fmkill.f:86).
FVS: GRINCR's **MORTS computes the regular (density+background) kill on the FULL pre-fire
stand** into WK2, then GRADD's FMKILL **adds** the fire kill to WK2 (capped at PROB); UPDATE
applies the sum. FVSjl runs `_maybe_burn!` (reduces TPA + recomputes density) BEFORE
`mortality!`, so the density mortality sees the thinned post-fire stand.
**Confirmed by a dense-burn trace** (fire_fuel9: SIMFIRE 2005, SDI 267) vs LIVE FVSsn:
| order                                            | 2010 TPA | |err| vs FVSsn 143 |
|--------------------------------------------------|----------|--------------------|
| fire→MORTS (current; density on post-fire stand) | 155      | 12 (under-kill)    |
| MORTS→fire (sequential/multiplicative)           | 125      | 18 (over-kill)     |
| additive MORTS_kill+fire_kill on full PROB, cap  | 83       | 60 (way over)      |
| **FVS (truth)**                                  | **143**  | —                  |
ROOT (traced from FVS): fmkill.f:86 is `IF(FIRKIL(I).GT.WK2(I)) WK2(I)=FIRKIL(I)` — the per-
record kill is **MAX(MORTS_kill, fire_kill)**, NOT the sum (a tree dies once, from whichever is
larger), both measured on the FULL pre-fire stand; the regular-mortality snags come from only
the EXCESS `WK2−FIRKIL = max(0, MORTS−fire)` (fmkill.f:135/FMSSEE) so fire+regular snags don't
double-count. FIX (grow_cycle!): on a fire cycle, run MORTS on the pre-fire stand → mk, restore
PROB, run the fire → fk, set survivors = pre − max(mk,fk), and book the FMSDIT snags from
`max(0,mk−fk)` only (mortality! `book_snags=false`). Non-fire path is byte-identical. Result:
fire_fuel9 2010 TPA 155→151 (FVSsn 143); the residual ~8 TPA is now isolated to the SEPARATE
FMEFF per-tree fire-kill distribution (F3 below). Suite 4494+21 (carbon/snag tests still pass).

### F3 — FFE fire behavior. PARTLY FIXED (moisture-code handling); FMEFF kill-distribution residual remains.
After F1, fire_fuel9 reads 2010 TPA 151 vs FVSsn 143; the residual is the fire BEHAVIOR, not
the ordering. Traced (fire cycle, vs live FVSsn FMFINT) — fuel-model SELECTION matches (10 wt
0.567, 5 wt 0.433), but the per-model Rothermel diverges:
| model | byramt FVS / FVSjl | xir FVS / FVSjl | spread R FVS / FVSjl |
|-------|--------------------|------------------|----------------------|
| 10    | 9531 / 2355        | 6463 / 5425      | 6.78 / 2.03 (3.3×)   |
| 5     | 14276 / 160        | 3174 / 706 (4.5×)| 19.71 / 1.01 (20×)   |
weighted byram FVS 4194 / FVSjl 1405, flame 3.17 / 1.69. Two factors: (a) MIDFLAME WIND —
FVS FWIND=2.0 vs FVSjl fwind=1.2. Narrowed: `FWIND = SWIND·WMULT`, `WMULT=ALGSLP(PERCOV,
CANCLS=[5,17.5,37.5,75], CORFAC=[0.5,0.3,0.2,0.1])` (fmburn.f:390) — tables + the PERCOV
formula `100·(1−exp(−TOTCRA/43560))` (fmcba.f) are IDENTICAL in FVSjl, so the divergence is
TOTCRA (Σ π·CW²/4·TPA): FVSjl's per-tree crown width CW is ~1.4× FVS's CRWDTH, pushing PERCOV
≥75 (WMULT 0.10) vs FVS ~50 (WMULT 0.167). The pre-fire stand TPA/BA match at 2005 (439/147),
LSPCWE=false for SN (so both take the CWCALC(iwho=0) branch, NOT the CWDS polynomial) — so the
root is the CWCALC crown-width value (or its CR/HT inputs) at iwho=0 for these SN species. The
spread rate `R = XIO·(1+PHIS+PHIW)/RHOBQIG` scales with PHIW=C1·(FWIND·88)^XM1, so the low wind
suppresses spread/byram.
(b) MODEL-5 LIVE-FUEL reaction intensity 4.5× low (xir 706 vs 3174) — the live-fuel-load /
moisture-of-extinction / reaction-velocity damping. The TPA impact is muted by mortality
SATURATION (flame is well above the kill threshold either way → only ~8 TPA / saturating BA),
which is why snt01 stand-4 (in the suite) still passes. Next: trace FVSjl `rothermel_surface_fire`
vs FVS FMFINT for (a) `fire_wind_reduction`/midflame wind and (b) the live-fuel reaction-intensity
terms, on model 5.

### F2 — COMPRESS (COMCUP) timing. Accepted-divergence-adjacent.
FVS runs COMCUP in GRINCR AFTER cuts+density+SDIAC (grincr.f:391); FVSjl's `apply_compress!`
runs BEFORE cuts (simulate.jl:201). COMPRESS is the accepted eigensolver divergence, but the
before/after-cuts timing changes which records exist when the thin runs — a potential
divergence beyond the eigensolver. Low priority (COMPRESS is explicitly accepted), flagged.

## Status
Core per-cycle sequence: traced and faithful. F1 (fire/MORTS MAX-combine) FIXED. Open: F3 (FMEFF per-tree fire-kill distribution) and
F2 (COMPRESS timing, accepted-adjacent). Continue the audit into: volume (fvsvol op order), thinning selection (RDPSRT/cut
order), and the height-growth (HTGF) internals — each traced directly against live FVSsn.

#### F3 update — moisture-code handling FIXED; root re-identified (the FWIND=2.0 was a red herring).
Re-traced the ACTUAL fire (not the PotFire report, whose FWIND=2.0 misled the first pass): the
actual fire's FWIND=1.2 MATCHES FVSjl, and FVSjl's `crown_width` matches FVS `CWCALC` BIT-EXACT
— so neither wind nor crown width diverges. The byram gap (FVS 4194 / FVSjl 1405) was the FUEL
MOISTURE: FVS used dryness model 3 (dead .07/.09/.14, live 1.0) but FVSjl used model 4 (.16/.16/
.18, live 1.75). Root: the SIMFIRE moisture code in fire_fuel9 is **9** (out of FVS's 1..4 range);
FVS's FMMOIS is a NO-OP for invalid codes, leaving the moisture at the per-cycle PotFire MODERATE
value (model 3), whereas FVSjl **clamped 9→4** (very wet). FIX: invalid SIMFIRE moisture resolves
to model 3 (matching FVS's leftover), not clamp-to-4 (keyword_dispatch.jl). For VALID codes 1..4
the FMMOIS tables were already identical, so this only touches the invalid-code edge case. Result:
fire_fuel9 2010 TPA 155→148 (FVSsn 143); suite 4494+21 unchanged. REMAINING (open): with the
moisture now model-3, FVSjl byram=2905 / flame=2.41 vs FVS 4194 / 3.17 — still ~1.4× low, so a
DEEPER Rothermel/fuel divergence remains (NOT moisture, wind, crown-width, ordering, the FMEFF
coefficients or the Regelbrugge-Smith mortality — all verified to match). Prime lead: FVS RE-
SELECTS the fuel model AFTER the wind/moisture are known (fmburn.f:393 `FIND FUEL MODEL AGAIN`)
and sums per-model byram with the post-wind weights; FVSjl selects once. Next FFE trace: the
second FMCFMD fuel-model selection + the per-model reaction-intensity (fuel loadings FWG) on
model 10/5. The FMEFF mortality formula + MORTB0/1/2 + Regelbrugge groups are confirmed bit-exact
vs fmeff.f, so the residual is upstream in the fire-behavior (byram/flame), not the kill logistic.

#### F3 further narrowing — the residual is the Rothermel REACTION INTENSITY (rothermel.jl vs FMFINT).
With moisture forced to model 3 (= FVS's) and the standard-fuel-model loadings matching, the
per-model Rothermel STILL diverges, in OPPOSITE directions:
| model | FVS xir / R / byram | FVSjl xir / R / byram |
|-------|---------------------|------------------------|
| 10 (mostly dead) | 5654 / 2.749 / 3382 | 6579 / 3.377 / 4758  (FVSjl HIGH) |
| 5  (live fuel)   | 2170 / 5.174 / 2562 | 1029 / 2.044 /  475  (FVSjl 2× LOW) |
Loadings + SAV + depth + moisture all match, so the divergence is in the reaction-intensity /
spread FORMULA. The opposite signs (dead-heavy model high, live-fuel model low) implicate the
LIVE-FUEL terms (live moisture-of-extinction / live-fuel reaction-velocity damping) as the main
differentiator, plus a smaller dead-fuel reaction offset. (Also unresolved: FVS's reported actual-
fire byram 4194 exceeds its own per-model weighted 3027 — the multiple FMFINT calls (actual /
potential / CFB ICALL=2) need disentangling to compare like-for-like.) This is a multi-term
Rothermel re-verification of `rothermel.jl` against `fmfint.f` — a real, bounded FFE task, not a
one-liner. snt01 stand-4 (the in-suite fire test) matches because its fuel models hit the terms
that already agree; fire_fuel9 (models 5/10) exposes the live-fuel + dead-reaction terms.

#### F3 RESOLVED (root cause) + residual localized to ONE size class.
Traced the fire_fuel9 byram divergence op-by-op against live FVSsn (instrumented FMFINT/FMCFMD/FMTRET):
1. **The Rothermel + fuel-model SELECTION are faithful.** FVSjl `_fmdyn(9.19,3.99)` returns
   {10@0.9716, 12@0.0284}, bit-identical to FVS; FMCFMD candidate logic == fmcfmd.f; live
   moisture-of-extinction (rothermel.jl:86 == fmfint.f:352) matches. The earlier "per-model
   byram" numbers were captured from the WRONG FMFINT call (potential-fire FWIND=2.0, not the
   actual fire FWIND=1.2) — they don't reconcile with the reported byram 4194; once instrumented
   correctly, the actual fire is FM10@0.972 + FM12@0.028 (= 0.972·3951 + 0.028·12514 = 4195).
2. **ROOT CAUSE — frozen FFE fuel. FIXED (committed).** The per-cycle FFE fuel loop
   (`ffe_fuel_update!` = decay + snag falldown + crown-lift→down-wood, + the inventory snag seed +
   crown-lift snapshot) was gated on the *Carbon report* (`carbon_on`). FVS fmmain.f runs it every
   cycle for any FFE-active stand. So a SIMFIRE-only stand's down-wood (FireState.cwd) stayed at the
   inventory value: fire_fuel9 sm=7.02 (==1990) at the 2005 fire vs FVS's accumulated 9.19, shifting
   FMDYN to pick FM5 over FM12 → byram 2905 vs 4194. Fix: gate on `ffe_on` (FFE active). carbon_on ⊆
   ffe_on, so carbon_jenkins stays bit-exact. fire_fuel9 1990-2005 now BIT-EXACT vs FVS; suite 4494+21.
3. **RESIDUAL — snag-bole cwd deposition is single-class, not cone-distributed.** Per-class dump at
   the 2005 fire: every class matches FVS within ~0.6 EXCEPT class 5 (6-12" logs): FVSjl 3.52 vs FVS
   1.46 (2.4×). `update_snags!` (snag.jl:138) dumps a fallen snag's WHOLE bole into the single size
   class of its DBH (`_cwd_size_class(dbh)`). FVS's FMCWD/CWD1 (fmcwd.f label 1000) distributes the
   bole down a CONE taper model across classes 5→4→3→2→1 (finds the heights where each diameter
   breakpoint lies, integrates volume per segment). So FVSjl overloads class 5 → sm/lg overshoot
   (10.87/6.19 vs 9.19/3.99) → post-fire TPA 121 vs FVS 143 (mild over-kill). The bole TOTAL is
   unchanged (carbon_jenkins DDW unaffected — de-risked); only the per-class split differs, and that
   split is exactly the FMDYN (SMALL,LARGE) input. NEXT: port the FMCWD cone distribution into
   update_snags! (and the CWD2 broken-top / CWD3 cut-tree paths) so class 5 splits correctly.

#### F3 — COMPLETE diagnosis of the remaining fire residual (timing/interleaving). [fix designed, not landed]
Tested the timing hypothesis directly (stash start-of-cycle SMALL/LARGE, fire on it) — DECISIVE result
that pins the exact mechanism, then reverted because it traded one scenario for another:
- **fire_fuel9** (fire 2005): post-`ffe_fuel_update!` cwd (11.0/6.03 → byram 5372, 28% high) over-kills
  → TPA 120 vs FVS 143. Start-of-cycle cwd (8.22/4.92 → byram 3912) → TPA 139 vs 143 (within 3%).
- **fire_early** (fire 2000): start-of-cycle made it WORSE (TPA 127 vs FVS 104; was ~114 within its
  loose ±12 tol). fire_early's fire is ALWAYS too weak (under-kill); the old full-period-advanced fuel
  was STRONGER and accidentally masked that under-kill. So the timing change exposed a pre-existing gap.
- **Unified root cause:** FVS does the surface-fuel loop ANNUALLY interleaved with FMBURN, so a fire at
  year Y burns on `S(cycle-start) + accumulation through year Y` — for a cycle-START fire that's
  start-of-cycle + ONE annual step. FVSjl batches the whole `per`-year `ffe_fuel_update!` BEFORE
  grow_cycle, so the fire sees the PERIOD-END fuel (≈5 yrs too much). FVS's fire-time value (fire_fuel9
  9.19) sits BETWEEN FVSjl's start-of-cycle (8.22) and period-end (11.0) — consistent with "+1 annual
  step". Pure start-of-cycle under-shoots (8.22 < 9.19) ⇒ fires a touch weak; pure period-end over-shoots.
- **THE FIX (designed):** split the per-cycle fuel loop around the fire — advance the annual fuel steps
  up to the fire year, fire (FMBURN on that fuel), then advance the remaining years. This is the FVS
  interleaving and should land BOTH fire_fuel9 (→~143) and fire_early (stronger ⇒ closer to 104) green.
  It's an architectural change to the FFE cycle (ffe_fuel_update! / grow_cycle ordering) needing
  validation vs both fire scenarios + carbon — deferred (not a one-liner; risk of perturbing carbon).
  Committed so far: the two prerequisites (per-cycle fuel runs at all + cone distribution). Suite 4494+21.

#### F3 — interleaving fix LANDED (committed); residual now a small ~4% per-model Byram calibration.
The annual fuel/fire interleaving fix landed green: when a SIMFIRE burns this cycle, the fuel loop is
split (advance 1 yr → stash the fire's (SMALL,LARGE) → advance per-1), so the fire burns on cycle-start
+ 1 annual step (FVS's FMBURN timing), not the period-end fuel. Non-fire/carbon stands keep the single
full-period call ⇒ carbon_jenkins bit-exact. Results: fire_fuel9 post-fire TPA 120→130 (FVS 143);
fire_early 114 vs 104 (within its ±12 tol; pure start-of-cycle had broken it at 127). Suite 4494+21, 0 fail.
**Remaining = small, opposite-signed, scenario-dependent fire residual** (NOT timing): fire_fuel9 slightly
OVER-kills (130 vs 143), fire_early slightly UNDER-kills (114 vs 104). The fire-effects chain is faithful
to ~3-4%: fire_early byram 7313 vs FVS 7597 (-4%), flame 4.10 vs 4.17, scorch 17.14 vs 17.58 — the TPA
gaps track the scorch deficit, so they're fuel/Byram-driven, not an FMEFF mortality-model bug. Opposite
signs ⇒ not a uniform bias; it's per-fuel-model Rothermel/loading deltas across different stand fuel
conditions. Closing it means matching the per-model Byram exactly (instrument FVS FMFINT per-model on
each scenario's selected models). Lower priority: both fire scenarios are within their test tolerances;
suite green.

#### F3 — final pinpoint: per-model Rothermel TERM divergence (FM10 reaction-intensity, FM5 spread).
With the interleaving fix in, FVSjl's fire-time fuel-model SELECTION now matches FVS: fire_early picks
{10@0.571, 5@0.429} vs FVS {10@0.563, 5@0.437} (SMALL 6.66 vs 6.72, LARGE 3.60 vs 3.28). So selection +
timing are faithful. The residual is purely the per-model Rothermel output, computed with IDENTICAL
inputs (FVS FMMOIS model 1 == FVSjl fuel_moisture(1) bit-for-bit: dead .05/.07/.12/.17/.40, live .55/.55;
loadings = Anderson standards; depth/mext match):
| model (fmois=1, wind=1) | FVS byram / xir / spread | FVSjl byram / xir / spread |
|----|----|----|
| 10 | 6519 / 6463 / 4.64 | 7337 / 7537 / 4.55  (xir +16%) |
|  5 | 8988 / 3174 / 12.4 | 5144 / 3051 / 7.46  (spread −40%) |
FM10's reaction intensity (xir) runs ~16% HIGH; FM5's rate-of-spread runs ~40% LOW (its xir matches to
−4%, so it's the SPREAD chain — packing ratio / wind factor / propagating flux, sensitive to FM5's depth
2.0 + live woody). They partially cancel in the weighted byram (fire_early 6397 vs FVS 7597, −16%), so
both fire scenarios stay within their test tolerances and the SUITE IS GREEN. Closing to bit-exact =
term-by-term audit of rothermel.jl vs fmfint.f intermediates (gamma/ir/mdcsa/beta/c1/phiw) per model —
the genuine last mile of the FFE surface-fire port. Lower priority (suite green; ~4% aggregate, opposite-
signed). Committed this session: frozen-fuel fix, cone distribution, annual fuel/fire interleaving.

#### F3 — RESOLVED. Two per-model Rothermel bugs fixed; fire now bit-close to FVS.
Term-by-term dump of FMFINT intermediates vs live FVSsn isolated the per-model byram divergence to TWO
compensating bugs (each regressed in isolation; together they make the fire bit-close):
1. **Zero-load class sort (xir).** Classes were sorted purely by SAV desc, so the zero-load dead-herb
   (SAV 1500) outranked the loaded 10-/100-hr (109/30) and the noclas cutoff dropped the REAL loaded
   classes. FM10 lost its 100-hr load (sum1 0.32 vs FVS 0.55) → beta1 0.010 vs 0.017 → rat 1.39 vs 2.35
   → gamma inflated → xir +16%. FIX: key the sort load>0 first (= FVS ISIZE swaps, fmfint.f:221-240).
   xir/sigma now match FVS to the digit (FM10 6463/1765, FM5 3174/1683).
2. **Missing slope (phis).** The fire Rothermel ran with slope_tan=0. FVS sets FMSLOP=SLOPE (fmsdit.f:72)
   → PHIS=5.275·tan²/beta1^0.3. fire_early's stand has a 30% slope (STDINFO field 5, already parsed to
   plot.slope=0.30) → phis=1.60, a large spread term FVSjl dropped. FIX: pass slope_tan=plot.slope to
   the actual-fire AND PotFire Rothermel calls.
RESULT: fire_early 2005 BIT-EXACT vs FVS (104/70/121/133), 2010 within 1-2 TPA; fire_fuel9 2010 141 vs
FVS 143 (was 120 over-kill pre-fixes). Updated the rothermel unit test whose 'FM10>FM5 byram' assertion
was an artifact of bug 1 — live FVSsn confirms FM5 (dry brush) byram 8988 > FM10 6519 at fmois=1. Suite
4494+21, 0 fail. Residual now ≤~2 TPA on late post-fire cycles (FP / small cwd-accumulation), within
tolerance. The F3 FFE surface-fire port is faithful end-to-end: frozen-fuel→cone-split→interleaving→
fuel-model selection→Rothermel(xir+spread+slope)→FMEFF mortality all match FVS.

#### Remaining broken: 11 = 1 COMPRESS (accepted) + 10 carbon-report dead-pool timing.
This turn closed 10 (the .csv->.tre writer boundary-overlap bug). The last 10 are the carbon report's
dead pools (Below-Dead/Stand-Dead/DDW) on INTERMEDIATE cycles — inventory + final cycles are bit-exact,
and ALL live pools (Aboveground/Merch/Below-Live/Floor/Shrub) are bit-exact. Root: FFE fuel-loop timing
is MULTI-TERM coupled, not a single lag. Empirically: moving ffe_fuel_update! from before-grow to
after-grow+crown-lift swung carbon_snt DDW from -1.9 (low) to +2.3 (high) — because it simultaneously
changes (a) the crown-lift application cycle (FMSDIT/FMCADD), (b) the falldown of the cycle's OWN new-
mortality snags (now created before the fuel loop), and (c) litterfall's pre- vs post-growth stand. FVS
gets all three right via its annual GRADD->FMSNAG/FMCWD/FMCADD interleaving (growth is still per-cycle,
but the fuel sub-steps are annual and ordered after GRADD). A faithful fix must interleave those fuel
sub-terms per-year with the correct stand state at each — a contained but careful FFE-carbon refactor,
not a one-line reorder (which over/under-shoots). Residuals ~0.5-2.3 t at F7.1 report resolution.

#### Carbon dead-pool residual CORRECTED: it's a down-wood MASS over-accumulation, not crown-lift timing.
Tested the targeted crown-lift fix (apply this cycle's crown-lift post-grow with intra-cycle decay, non-
fire stands only): it made carbon_snt DDW WORSE (Δ 1.2 → 1.6), proving the residual is NOT the crown-lift
lag. The DDW is already HIGH and COMPOUNDS over cycles (carbon_snt run_keyfile Δ = 0.0/0.3/1.2/1.6 at
1990/95/00/05) — FVSjl accumulates ~15-20% too much down wood, the SAME overshoot seen on fire_fuel9
(cwd 10.87 vs FVS 9.19, +18%). Stand-Dead Δ is small (0.1/0.4/0.3/0.7). So the 10 carbon @test_broken are
a down-wood-MASS calibration in the per-cycle additions (snag-bole falldown density/biomass + CWD2B crown
flow) vs FVS FMSNAG/FMCADD — NOT a timing lag. Closing it needs an instrumented term-by-term comparison
of FVSjl's per-cycle cwd additions against FVS's FMSNAG/FMCWD/FMCADD (snag fall density, bolevol vs
SNVIS·V2T, CWD2B TFALL schedule). Reverted the crown-lift change (overshoots). Live + Floor pools stay
bit-exact; inventory + final cycles bit-exact; the gap is intermediate-cycle DDW/Stand-Dead at F7.1.

#### Carbon DDW residual — CORRECTED again (column misparse fixed): endpoints bit-exact, intermediate-cycle timing.
The DDW report column is #7 (YEAR Total Merch Live Dead Dead **DDW** Floor). Reading it correctly vs the
FVS .save: 1990 5.82/5.8 ✓, 1995 5.69/5.4 (+0.29), 2000 7.22/8.4 (−1.18), 2005 11.43/11.4 ✓. So DDW =
Σcwd[1:9]·0.5·metric is BIT-EXACT at both endpoints; the cwd MASS matches FVS at 1990 and 2005 (10.199 vs
10.191). The divergence is INTERMEDIATE cycles, NON-monotonic (1995 high, 2000 low) — a multi-term down-
wood ADDITION-TIMING issue (crown-lift FMCADD + CWD2B mortality-crown flow + snag falldown), NOT mass
over-accumulation (my earlier 13.0 was a column misparse that also invalidated the crown-lift-fix
'overshoot' reading) and NOT a single one-cycle lag (a pure lag would be low every intermediate cycle; here
1995 is HIGH). The endpoints being exact means the SOURCES are right; only the per-cycle phasing of the
additions vs FVS's annual GRADD→FMCADD interleaving is off. Closing it needs per-term timing alignment
(crown-lift same-cycle, CWD2B flow phase, falldown phase) validated together — a careful FFE-carbon refactor.

#### Carbon dead-pool residual — THREE concrete terms identified (each small, compounding).
Instrumented each contributor for carbon_snt:
1. **Input-snag bole VOLUME** (test #548): FVSjl snag_bole_carbon·TO = 3.924 vs FVS 3.8 (+3.3%). Both
   input dead records have NO height (ht=0 in the .tre) → FVSjl estimates via _htdbh_height (dbh 34.6→
   92.7 ft, dbh 7.2→53.7 ft); FVS's HTDEAD/FMSVL2 estimate differs. Drives part of the 1995 over-add.
2. **Crown-lift application cycle**: FVSjl applies it the NEXT cycle's fuel loop (computed post-grow);
   FVS applies same-cycle (FMCADD after GRADD). Drives the 2000 under-add (~cl_1 missing).
3. **CWD2B crown-flow phasing**: the mortality crown sits as Stand-Dead until TFALL drops it to DDW;
   committed flows it one cycle late (SD high/DDW low), the after-grow reorder flows too fast (SD low/
   DDW high), FVS is between. Drives the SD↔DDW split.
The committed timing balances #2/#3 so BOTH endpoints (1990, 2005) and all live pools are bit-exact;
every single-term change regresses an endpoint (verified vs live oracle, col-7). Closing all 10 needs:
fix the input-snag HTDEAD estimate to match FVS (#1, instrument FMSVL2), then jointly re-phase #2+#3 to
FVS's annual order while preserving the endpoints. Multi-step, multi-term FFE-carbon calibration.

#### Term #1 (snag bole volume) localized to ONE snag — small-tree R8Clark vs FVS FMSVOL.
Instrumented FVS FMSVOL per snag (carbon_snt): the BIG oak (sp65 dbh34.6 ht92.7) MATCHES — FVSjl bole
2.033 vs FVS 2.032 t/ac. The residual is entirely the SMALL tree (sp27 dbh7.2 ht53.7, ht matches 53.738):
FVSjl R8Clark cubic = 5.2 cuft vs FVS FMSVOL/CFVOL = 4.8 cuft (+8%), so FVSjl bole 1.471 vs FVS 1.358 t/ac.
Fixing just this one snag → snag_bole_carbon 3.924 → 3.80 = FVS 3.8 (closes test #548). FVS's snag bole is
FMSVOL = CFVOL integrated to HTIS (full dubbed stem); FVSjl's seed path uses _R8CLARK_VOL with the SCF/
merch (prod,stump,mtopp) spec, which over-estimates small-tree cubic by ~8%. The faithful fix is to volume
the input snag with the full-stem CFVOL (FMSVOL), not the merch R8Clark path. (Large trees already agree,
so live-tree .sum cuft stays bit-exact.) This is term #1 of the 3; the DDW tests additionally need the
crown-lift cycle (#2) + CWD2B phasing (#3).

#### PROGRESS: snag-bole fix landed (11→10 broken); crown-lift timing confirmed CORRECT.
- **Term #1 FIXED & committed:** input-snag bole now uses the merch cubic v[4] (= FVS FMSVOL), not gross
  v[1]. snag_bole_carbon 3.92→3.77 (FVS 3.8) — closed test #548. Also made 1990 Stand-Dead bit-exact and
  improved 1995 DDW (+0.29→+0.2). Suite 4505 pass + 10 broken, 0 fail.
- **Term #2 (crown-lift timing) is NOT a bug:** tested applying crown-lift same-cycle (post-grow) — it
  OVERSHOOTS (2000 DDW +1.1, 2005 +1.6), proving FVS ALSO lags the crown-lift one cycle (the 2000-2005
  crown is NOT in the 2005 report; the committed pre-grow timing already matches FVS). Reverted.
- **Remaining 2000 DDW −1.3 = CWD2B/decay phasing** (not crown-lift): the 1995-2000 cycle's down-wood
  additions (CWD2B mortality-crown flow + snag falldown) under-accumulate mid-run, recovering by 2005
  (endpoints bit-exact). This + the Stand-Dead CWD2B-crown residual (~0.3-0.6) are the last carbon terms.

#### PROGRESS update: 11→9 broken (two more faithful carbon fixes this turn).
- **Below-Dead roots FIXED & committed** (test #205): FVS decays input-snag dead roots assuming 10 yrs
  since death (XDCAY=(1-CRDCAY)^10, fmsadd.f:313-320); FVSjl booked the full Jenkins root. Now Below-Dead
  is BIT-EXACT every cycle (1.0/1.3/2.0/3.7 = FVS).
- **Mortality-snag bole = MERCH: tested, REVERTED.** Switching book_mortality_snags to merch_cuft_vol
  made Stand-Dead WORSE (Δ 0.3→0.9). So FVS uses MERCH for the dead-input seed path (FMSVOL) but the
  mortality-created snag bole behaves like total cubic — keep cuft_vol there. (Fire stayed green either way.)
- **Remaining 8 carbon = CWD2B crown-flow phasing only.** Below-Dead + snag bole + all live pools +
  endpoints now bit-exact. The residual is Stand-Dead +0.3-0.6 / DDW -0.8 to -1.3 on intermediate cycles:
  the mortality CROWN (CWD2B) flows from Stand-Dead to DDW too slowly (TFALL/TSOFT fall schedule, fmscro.f).
  Session total: 21 → 9 broken. Method that's working: find the FVS semantic (FMSVOL merch, XDCAY root
  decay), verify vs the live oracle, land only fixes that improve without regressing the exact pools.
