# FVSsn → FVSjl faithfulness audit (operation-sequence trace)

> ## ⚠ STATUS RECONCILIATION (2026) — all open items below are now CLOSED; this doc was stale
> The open items flagged in this op-sequence audit were resolved by subsequent work but the doc was not
> updated (which is misleading — fixed here). Each re-verified vs live FVS:
> - **KNOWN GAPS (IDG==1, G==0)** — RESOLVED. `_backdate_dbh!` gates on IDG; growth_idg1 init crown +
>   .sum bit-exact vs live. (See the corrected KNOWN-GAPS entry below.)
> - **F1 (fire/MORTS MAX-combine)** — FIXED (in this doc).
> - **F3 (fire_fuel9 per-model Rothermel byram / fire mortality)** — RESOLVED by the 2026 fire-basis
>   (sm,lg) start-of-cycle fix (docs/audit/INDEX.md "FIRE UNDER-KILL"): fire_fuel9 is now BIT-EXACT vs
>   live every cycle (2010 TPA 143 = live; was 148-155).
> - **DDW size-5 / carbon_jenkins** — RESOLVED by the CRATET init-crown backdating (4519/0/1 spec met).
> - **F2 (COMPRESS timing)** — accepted-adjacent (the accepted COMPRESS eigensolver divergence).
> The live current state: suite 4530 pass / 1 broken (= accepted COMPRESS); the campaign ledger is
> docs/audit/INDEX.md. No open op-sequence divergence remains beyond the accepted COMPRESS + the
> documented FFE fuel-phasing intermediate-snapshot residual (INDEX.md #22).

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
[2026-07-07] RE-TRACE (doctrine #9) — a LS "1.4× forest CW" claim from FIA validation was a MISREAD,
now RETRACTED: it compared FVS's FOREST-grown text-.trl CW column against FVSjl's DBS treelist
CrWidth column, which BY DESIGN uses the OPEN-grown CW (iwho=1, CR=90; see dbs_output.jl:532). Direct
check proves FVSjl's `crown_width` is BIT-EXACT for LS both ways: GA dbh7.5 forest(iwho=0,cr35)=
14.8142 == FVS bechtold 54401 (2.9672+1.3066·D+0.0585·CR); open(iwho=1,cr90)=21.0394 == FVS ek 54403
(4.755·D^0.7381). So there is NO LS forest-grown crown-width bug — jl's iwho=0 path is correct. The
small LS-hardwood `.sum` residual (CCF ~4%, BdFt ~8%) is a SEPARATE, still-open small item (open-grown
CCF aggregation / hardwood board), NOT this crown-width equation. The SN PERCOV finding below is
SN-specific (LSPCWE/CWCALC) and should be re-verified the same careful way before any shared-code change.
The
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

#### Carbon Stand-Dead/DDW — final localization (2 coupled deep terms remain).
Decomposed Stand-Dead into bole + crown vs the live oracle:
- **CROWN (CWD2B) MATCHES FVS** — FVSjl 0.44/0.77/1.56 vs FVS 0.42/0.76/1.54 t-C/ac. So the crown
  amount + TFALL flow are CORRECT. Not the bug.
- **Mortality snag BOLE +0.13-0.25** (FVSjl 2.97 vs FVS 2.72 @2005). BUT switching it to the merch cubic
  (merch_cuft_vol OR fresh R8Clark v[4]) makes it WORSE (1.67→1.93) — because the dying record is a TRIPLED
  record whose volume ≠ the full-tree volume; cuft_vol is the closest available. So the merch fix that
  worked for INPUT snags does NOT transfer to mortality snags (tripling interaction). Reverted.
- **DDW deficit (−0.8 to −1.3 @intermediate) = new-mortality falldown LAGS one report-cycle.** update_snags!
  runs at the cycle start (before that cycle's mortality is created) and ages by (current−deathyr), so a
  snag created this cycle (deathyr=cycle yr) falls 0 this cycle and a full cycle next — it stands through
  the current report (Stand-Dead high) and reaches DDW a cycle late (DDW low). FVS spreads deaths over the
  cycle so they fall ~per/2 in the death cycle. The fix is a partial (~half-cycle) falldown of fresh
  mortality, but the prior FULL-falldown attempt collapsed Stand-Dead (memory), so it needs the spread.
Both remaining terms are genuinely coupled to TRIPLING and the death-spread — deeper than the input-snag
fixes that closed #548/#205. Session: 21 → 9 broken. Growth + cuts audited & bit-exact (not in the 9).

#### Verified-before-revert (per the regression-may-unmask guidance):
- **Mortality snag bole = TOTAL cubic (cuft_vol), NOT merch.** Re-applied the merch-v[4] fix (faithful per
  fmdout→FMSVOL) to mortality snags; it broke carbon_jenkins's snag_bole @test (3.72/3.28) which were
  PASSING and bit-correct with cuft_vol, AND it over-counted carbon_snt DDW. Per-record v[4] IS lower than
  cuft for big trees yet the aggregate rose — confirming the dying-record volume state interacts with
  tripling. The passing tests were CORRECT (not masking), so the merch fix is wrong for MORTALITY snags
  (only the INPUT-snag seed path uses merch, #548). Reverted with that verification.
- **Death-spread falldown (force_yrs≈per/2) breaks the bit-exact 2005 ENDPOINT** (DDW Δ0→+1.4). The
  committed one-cycle lag is correct at the endpoints (the lag washes out); the death-spread, while
  conceptually faithful (FVS dates deaths across the cycle), over-adds to DDW cumulatively. So it trades
  the exact endpoint for the intermediate cycle — not a net win. Reverted with verification.
Net: the 8 carbon intermediate-cycle terms are coupled to tripling + the death-spread in ways where the
faithful single-term fixes break CORRECT (passing/endpoint) values. They need a joint death-time-volume +
death-spread treatment that preserves the endpoints — genuinely deep. Session: 21 → 9 broken; growth+cuts
audited bit-exact.

#### RESOLVED — snag BOLE is MERCH cubic (merch_cuft_vol), correcting the prior "total" conclusion
The earlier "mortality snag bole = TOTAL cubic" conclusion above was WRONG — it rested on two flawed
checks: (a) it tried merch = v[4] ALONE, which undershoots loblolly pine (carbon_jenkins) where v[7]≠0;
(b) it validated against the 1-decimal `.report.save`, whose print-rounding masked the real signal.

Verifying the LOGIC (not just the oracle) settled it: FVS's snag Stand-Dead bole is `FMSVOL` →
`VOL2HT = MAX(X, MCF)` (fmsvol.f, SN branch), and for SN every species uses NATCRS, whose MCF is exactly
jl's per-tree `merch_cuft_vol` = `v[4]+v[7]` (with the DBHMIN gate + the Region-8 <10ft-product rule).
That is NOT the gross `cuft_vol` (= v[1]). Differential vs an instrumented FMDOUT (per-snag BOLE/CRWN +
standing DENSITY) on the bit-exact-growth carbon_snt:
  - standing snag DENSITY is BIT-EXACT every cycle (14.76/42.70/48.03/71.04) — falldown was never wrong;
  - gross v[1] ran the per-snag bole 2-8% high on mid/large snags → Stand-Dead +0.6 over four cycles;
  - merch `merch_cuft_vol` matches FVS per-snag (8.1/4.8/16.5/18.7/9.8 …) → Stand-Dead now ≤0.03.
Fix: `book_mortality_snags!` books `bolevol = max(0.005454154·H, merch_cuft_vol) · V2T/2000` (the X-floor
is FMSVOL's tiny-tree cone volume, so sub-merch snags keep a small positive bole instead of hitting the
Jenkins fallback). Input-snag seed path unchanged. crown (CWD2B) was already exact.

Test consequence: Stand-Dead validation moved off the 1-decimal save onto the HIGH-PRECISION instrumented
oracle (BOLE+CRWN, consistent with the bole/crown component tests that already used 3.72/1.46) — four
Stand-Dead @test_broken/@test flipped to PASSING (carbon_snt LIVE pools, age-aware falldown, carbon_jenkins
standing_dead). Full suite 9 → **7 broken, 0 fail**. The remaining 6 non-COMPRESS broken are ALL the one
DDW (down-wood) gap: jl DDW 6.84 vs FVS 8.41 at 2000 — the snag→down-wood falldown/decay phasing, the sole
remaining FFE dead-pool divergence. StandDead, live pools, growth, cuts: bit-exact.

#### DDW (down-wood) — the sole remaining FFE dead-pool gap: crown-lift PHASING
After the snag-bole fix, the 6 remaining non-COMPRESS broken tests are ALL the DDW column. Localized via
an instrumented FMDOUT dump of `CWD(3,k,1,5)+CWD(3,k,2,5)` per size class on carbon_snt:
```
size:     1      2      3   |  4      5      6     7     8    9     (biomass t/ac @2000)
FVS:    0.368  1.426  2.358 | 1.245  1.219  0.527 ...               Σ1-3 = 4.15
FVSjl:  0.194  0.762  1.622 | 1.138  1.409  0.540 ...               Σ1-3 = 2.58
```
The gap is ENTIRELY in the FINE classes (1-3, <3"); the bole/large classes (4-9) match. At 1995 the small
classes match exactly — the gap opens at 2000. Per-cycle additions to sizes 1-3: 1995→2000 FVS +1.94 vs
FVSjl +0.34; 2000→2005 FVSjl +2.56 vs FVS +1.88 — FVSjl's fine-material additions are shifted ~ONE CYCLE
LATE. Source = the live-tree crown-lift (FMCADD fmcadd.f:95-101: shed lower crown OLDCRW→down-wood, sizes
1-5), the dominant fine-fuel term.

LOGIC (FMMAIN): FVS grows (base GRINCR) → FMDOUT report → year-loop {FMSNAG, FMCWD decay, FMCADD
(litterfall+LIMBRK+crown-lift from post-grow crown vs OLDCRW + CWD2B fall)} → FMOLDC. The crown-lift is
applied WITHIN the cycle's year-loop, per-year, with decay interleaved. FVSjl's loop is report →
ffe_fuel_update (applies the PREVIOUS cycle's crown_lift_annual) → grow → compute_crown_lift → snapshot:
the crown-lift computed in cycle c is applied in c+1's fuel_update → an extra cycle of lag.

Two experiments (carbon_snt, scripts only — no code changed):
- Apply crown-lift 1 cycle earlier as a LUMP at end-of-cycle: StandD preserved but DDW OVERSHOOTS
  (9.77 vs 8.41 @2000) — the lump skips the per-year decay the lagged path applies.
- Move the whole ffe_fuel_update after grow: DDW overshoots AND StandD collapses (snag/CWD2B double-proc).
So the fix is NOT a simple reorder: the crown-lift must be applied in its OWN cycle's year-loop WITH the
per-year decay interleaved — i.e. align FVSjl's loop to FVS's FMMAIN phasing (grow→report→year-loop→oldcrw)
so compute_crown_lift feeds the same cycle's fuel loop. Risk: the main path (snt01) has REGEN that grows
the tree list each cycle, so the OLDCRW per-tree state must survive record tripling/compaction (FVS does
this in FMTDEL/FMTRIP/FMCMPR) — a stable record id, not an index snapshot. That is the focused next step
for the DDW column (the last FFE dead-pool item). StandD/live/growth/cuts are bit-exact.

#### DDW progress — three verified fixes (crown-lift snapshot + SCNV); residual now ~0.2 in size 5
The DDW gap (jl 6.84 vs FVS 8.41 @2000, i.e. -1.57) was closed to +0.2 by two further logic-verified fixes
after the snag-bole merch fix:
1. **Inventory crown snapshot (FMOLDC).** FVS calls FMOLDC in the inventory FMMAIN before the first grow,
   so cycle 1's crown-lift has a valid OLDCRW. FVSjl only snapshotted at cycle END → the 1st cycle's
   crown-lift was skipped (ffe_oldht=0), losing ~1.9 t/ac of fine down-wood. Fix: snapshot at inventory in
   write_carbon_report + write_sum_file. Result: down-wood size classes 1-3 now BIT-MATCH FVS (jl
   0.38/1.46/2.36 vs 0.37/1.43/2.36 @2000). Verified jl crown_biomass == FVS CROWNW bit-exactly (so the
   crown weight model is right; the gap was purely the missing 1st-cycle term).
2. **FMCWD SCNV soft/hard factor.** A falling snag's bole biomass is DIF·V2T·SCNV(K), SCNV=(0.80 soft,
   1.00 hard) (fmcwd.f:61) — a SOFT (decayed) snag contributes only 0.80× its volume. update_snags! used
   `a·dfall` with no soft/hard split, over-counting fallen soft boles 1.25× (the size-4/5 overshoot that
   grew as snags softened). Fix: `add = a·(dfis·0.80 + dfih)`. Overshoot +0.5 → +0.2.
Remaining: DDW @2000 8.6 vs 8.4 (+0.2), concentrated in size class 5 (6-12": jl 1.27 vs 1.22 @2000, 1.67
vs 1.55 @2005) — a small, distributed fallen-bole-taper / decay residual in the medium-bole class, growing
with bole falls. Sizes 1-4 + 6-9 match within ~0.03. Suite: 4512 pass / 7 broken / 0 fail (no regression).
The 6 DDW tests are still @test_broken (the +0.2 print-rounds past 0.05) but the column is now ~5x closer.

#### DDW size-5 residual — ruled out normalization + hard/soft pool; still open
Two more candidates for the size-5 (+0.05..0.12) residual were tested and RULED OUT:
- Cone-taper normalization: removing jl's `f[j]/=total` had NO effect (the raw conic fractions already
  sum to ~1.0 for these snags), so the taper distribution is not the cause.
- Hard/soft pool assignment: FVS splits the fallen bole into the soft pool (K=1, DIS, decays ×1.1
  faster) vs the hard pool (K=2, DIH); jl had dumped all fallen bole into the hard pool (index 2).
  Splitting it (addS=a·dfis·0.80→cwd[:,1,:], addH=a·dfih→cwd[:,2,:]) is MORE faithful (kept) but had no
  measurable effect on carbon_snt — its snags stay mostly HARD over the 4 cycles, so dfis≈0.
Confirmed: ALL species' crown weight in size 5 = 0 (jl crown_biomass == FVS CROWNW), so size 5 is fed
ONLY by fallen boles; input-decay matches (size-6 input 1.01→0.73 bit-matches FVS). The residual is jl's
medium-snag (d 6-12") boles contributing ~10%/bole more to size 5 than FVS, growing with bole falls.
Next probe: instrument FVS CWD1/CWD2 (fmcwd.f) per-snag ADD to size 5 vs jl, to find whether it's the
medium-snag merch bolevol or a taper edge. DDW now +0.2 (was -1.57); 6 DDW tests still @test_broken.

#### DDW residual ROOT CAUSE — snag CLASS-AGGREGATION (FMSADD) vs jl per-record snags
After six verified fixes took DDW from -1.57 to -0.1/-0.2, the residual is now traced to a STRUCTURAL
difference, not a tunable sub-term: FVS aggregates snags into discrete CLASSES and jl keeps per-record.
FMSADD (fmsadd.f:255-340) bins every new dead-tree record into a (SPCL, DBHCL, HTCL) class — DBHCL =
INT(DBH/2 + 1) (2-inch DBH classes, 19 for ≥36"), HTCL = 1/2 by MIDHT, SPCL = species class — and
DENSITY-WEIGHT-AVERAGES the new deaths into that class's aggregated record: DBHS(X) = (DBHS(X)·DEND(X) +
DBH(I)·SNGNEW(I))/TOTDEN, likewise HTDEAD(X). So FVS carries ~43 aggregated class-records with averaged
DBH/HT; FVSjl carries ~109 per-tree records (incl. tripling).

Why it shows ONLY in bole-fall: the bole-fall biomass = Σ vol(DBH)·fallrate(DBH)·density is NON-LINEAR in
DBH (vol≈DBH², fallrate linear-decreasing), so falling N per-record snags ≠ falling their density-averaged
class. Stand-Dead (a near-linear Σ vol·density of STANDING snags) is insensitive to the averaging → it
matches bit-exact; the DBH-rate-weighted FALL is not → bole-fall diverges ~13% on the new-mortality snags
(cycles 2-3; cycle-1 input snags, already one class each, match exactly). Verified every other sub-term
(cone fractions ≤1%, crown-lift per-tree EXACT, decay rates, soft fraction 6%, merch bolevol cycle-start,
fall density linear-in-origden) matches FVS — none is the cause.

To CLOSE bit-exact: replicate FMSADD's snag classification (bin into SPCL/DBHCL/HTCL, density-weight-average
DBH/HTDEAD/bolevol per class) so the snag list IS the aggregated class set, then fall/decay the classes.
This is a structural rewrite of the snag representation (FireState.snags), touching StandD/SnagSum/fire/
carbon — a scoped multi-session port, NOT a quick fix. DDW currently 8.3 vs 8.4 (was 6.84); the 6 DDW tests
remain @test_broken pending this. Everything else (StandD, live pools, floor, growth, cuts, DDW sizes
1-3/6-9) is bit-exact.

#### CORRECTION to the snag-aggregation claim above
Tested empirically: the FMSADD class-aggregation effect is only ~0.3-1%, NOT the ~13% bole-fall gap. A
2-inch DBH class averages to nearly the same vol·rate (over [8,10]": jl per-record Σvol·rate ≈ 7.96 vs
class-averaged 7.94, ~0.3%; within-class vol convexity ~0.6%). carbon_snt has 109 records → 25 classes,
20 multi-DBH, but the spread is ≤2" per class so averaging barely moves it. So aggregation is a real FVS/
FVSjl difference but does NOT explain the residual — the earlier "ROOT CAUSE" claim is RETRACTED.

HONEST state: the ~13% bole-fall under on new-mortality snags (cycles 2-3) is NOT root-caused. The
measurements are mutually inconsistent under conservation: standing snag DENSITY matches FVS bit-exact
(SDEN 14.76/42.70/48.03/71.04), Stand-Dead BOLE matches (so per-snag merch bolevol matches), yet the
measured bole-FALL is 13% under (jl 1.27/1.43 vs FVS-instrumented CWD1 1.46/1.63). Fallen density jl
7.3/38.4/41.5. If standing density + bolevol both match, conservation says bole-fall should match — so
either a measurement artifact remains or a subtle interaction (e.g. WHICH snags fall by size, the
hard/soft split timing vs CWD1's DFIH/DFIS) not yet isolated. DDW is 8.3 vs 8.4. This is the honest
open question; the 6 verified fixes (merch bole, FMOLDC snapshot, SCNV, hard/soft pools, decay-timing,
step-softening) all stand and took DDW from -1.57 to -0.1/-0.2.

#### DDW BREAKTHROUGH — soft-transition was the bole-fall bug; bole now BIT-EXACT
The retracted "aggregation" guess was wrong; the real bug was found by instrumenting FVS's per-cycle
hard/soft FALLEN density (FMSNAG DFIH/DFIS): FVS's fallen snags are 100% HARD (DFIS=0.000 EVERY cycle on
carbon_snt). FVS's DENIH/DENIS are the snag's INITIAL hard/soft state at CREATION (mortality snags are
created hard → DENIH); the per-snag HARD flag that flips at DKTIME is a separate decay/reporting state and
does NOT move the fall density. FVSjl wrongly moved den_hard→den_soft at DKTIME, so it applied the SCNV
0.80 soft factor to ~13% of fallen bole that FVS treats as hard → the bole-fall under-count.

Fix: removed the hard→soft density transition (snags fall hard; SCNV 0.80 applies only to snags SEEDED
soft). Result: down-wood size classes 4-9 (the BOLE-fed classes) are now BIT-EXACT vs FVS (jl 1.248/1.222/
0.526/0.305/0.045 vs FVS 1.245/1.219/0.527/0.307/0.045 @2000, within 0.003). SnagSum still passes (it
reports the initial hard/soft state). DDW -1.57 → +0.1/+0.2.

REMAINING (small, localized): DDW now +0.1/+0.2 OVER, entirely in the FINE crown classes 1-2 (size 2 +0.05/
+0.08). Instrumented FVS FMCADD per-source: the snag-crown (CWD2B) fall MATCHES at the report (jl year-1
pool 0.438 vs FVS 0.429); the over is the live-tree CROWN-LIFT to size 2 (~+5%: jl 0.174 vs FVS 0.166/yr)
+ LIMBRK (~+5%). Crown-lift matches PER-TREE (tree-1 size-2 lift 0.0942 = FVS exact), so the aggregate +5%
is likely the tripling × crown-lift interaction (ffe_old* per-record state across record tripling — the
known OLDCRW-survives-tripling concern). DDW 8.5 vs 8.4. The 6 DDW tests remain @test_broken on this last
~+5% crown over; bole, Stand-Dead, live pools, floor all bit-exact.

#### DDW CLOSED on bit-exact growth — crown small-tree bole = MERCH (FMSVL2), not gross
The +5% crown-lift-to-size-2 over (the last DDW residual) was the crown_biomass SMALL-TREE BOLE. For
D<DBHMIN, FVS FMCROWE adds an estimate of the whole-tree bole to the crown weight (TTOPW) using FMSVL2
(= MAX(X, MCF), the MERCH cubic with the tiny-tree cone floor X=0.005454154·H; fmcrowe.f:256-285). FVSjl
used `_fm_cuft` = GROSS cubic (v[1]) → the small-tree bole over by 1.5-2× → crown size-2 over (sp33 d1.5:
jl 2.643 vs FVS 1.759; d2.2: 7.82 vs 3.953) → fed the crown-lift/LIMBRK size-2 → DDW +0.1/+0.2. (Found by
checking crown_biomass per-SPECIES vs instrumented FVS CROWNW: sp22/27 matched, sp33 d≥1 was 1.5-2× over.)
Fix (crown_biomass.jl): `vt = max(0.005454154·hmin, _fm_cuft(...; merch=true))`. Result: sp33 crown
bit-exact (1.761 vs 1.759), and **carbon_snt DDW BIT-EXACT: Δ 0.002/0.007/0.001** (was -1.57).

Eight FFE-carbon fixes this session took DDW from -1.57 to bit-exact on the bit-exact-growth fixture. The
carbon_snt DDW + the live run_keyfile DDW @test_broken FLIPPED to passing @test. Suite 4515 pass / 5 broken
/ 0 fail. Remaining 5 broken = 1 COMPRESS + 4 carbon_jenkins DDW (2000/2005, in the driver + .out-writer
tests): carbon_jenkins is the SYNTHETIC non-bit-exact-GROWTH LP fixture (diameter-growth calibration tail),
so its mortality diverges (jl DDW DROPS 2.5→2.0 @2000 where FVS grows to 3.8 — fewer deaths) → its DDW
can't match FVS's save. That is the LP growth/mortality subsystem (WK3 past-dbh calib), NOT the carbon
model, which is now bit-exact. The carbon_jenkins DDW @test_broken honestly tracks that growth tail.

#### carbon_jenkins DDW — CRATET init crown closed most of it; residual is init-crown ~5% precision
carbon_jenkins's DDW was NOT a growth tail per se — its StandD + snag bole MATCH FVS (5.13 vs 5.2; bole
3.67 vs 3.72), so the snags are right; the DDW dropped because the down-wood ADDITIONS (crown-lift) were 0.
Cause: carbon_jenkins's .tre has NO input crown (crown_pct=0), and FVSjl never estimated it at init (FVS
does, in CRATET/INITRE, which calls CROWN for missing crowns). Fix: init_crown_ratios! runs the CROWN model
on inventory trees with crown_pct=0 (input crowns untouched), wired into setup_growth!. carbon_jenkins DDW
Δ 1.75→0.1 @2000; carbon_snt byte-identical. The carbon_snt + run_keyfile DDW tests PASS; suite 7→5 broken.

Residual: jl's init crown is ~5% low vs FVS's CRATET ICR (jl [21,26,30,34,38,42] vs FVS [22,27,32,36,40,47]
for the 6 inventory trees). Both call CROWN, so it's the CROWN model's INIT context (the Weibull crnew at
init, before the ±1%/yr limit which only applies once a prior crown exists — so the per-cycle crown matches
snt01 bit-exact via the limit, but the raw init Weibull is ~5% low). This ~5% init-crown feeds cycle-1 DGF →
slightly non-bit-exact LP growth → carbon_jenkins DDW Δ0.1 (still > 0.05). The 4 carbon_jenkins DDW
@test_broken track this init-crown precision. Likely the init CCF/SDI (computed with crown_pct=0) feeding
the Weibull `scale`/`relsdi`. The DDW CARBON MODEL itself is bit-exact (carbon_snt).

#### carbon_jenkins residual — init-crown ~5% = CCFCAL formula + DENSE calibration BACKDATING (deep)
The carbon_jenkins init-crown ~5% (jl ICR 21 vs FVS 22) traces to RELDEN (the CCF feeding CROWN's percentile
SCALE): jl stand_ccf=241 vs FVS RELDEN=178. Instrumenting FVS DENSE/CCFCAL showed TWO causes, NOT the crown
ratio (sp13's CW is CR-independent): (1) FVS's CCFCAL returns CCFT directly with CW=0 — a different CCF
formula than jl's 0.001803·crown_width²·tpa (which matches FVS for carbon_snt's hardwoods but not sp13/LP);
(2) at CRATET, DENSE runs with calibration BACKDATING — FVS's dbh is the PAST dbh (6.24 vs the current 8.0),
so the CCF is computed on smaller trees → RELDEN 178 < 241 → SCALE 0.869 > 0.764 → higher init crown. Both
need porting CCFCAL's per-species CCF + DENSE's calibration-time dbh backdating to bit-match the CRATET init
crown. That ~5% init crown feeds cycle-1 DGF → slightly non-bit-exact LP growth → carbon_jenkins DDW Δ0.1.
A deep, multi-routine calibration on ONE synthetic LP fixture; the DDW CARBON MODEL is bit-exact (carbon_snt).

#### carbon_jenkins — CONFIRMED: jl crown_width matches FVS exactly; residual is ONLY the backdated dbh
Definitive: computing jl's stand_ccf on FVS's calibration-backdated dbh (6.24/8.27/… from the CCFDBG dump)
gives 178.32 — matching FVS RELDEN=178.39 to ULP. So jl's crown_width + the 0.001803·CW²·TPA CCF formula
are CORRECT (the earlier "CCFCAL formula" note was wrong; CW=0 in the dump was just an unused output arg).
The ENTIRE carbon_jenkins init-crown error is that FVS computes the CRATET crown's CCF on the BACKDATED dbh
(DENSE's calibration backdating, ~1.7" younger), while jl's init_crown_ratios! uses the current dbh →
CCF 241 vs 178 → SCALE 0.764 vs 0.869 → init crown ~5% low → cycle-1 LP growth slightly off → DDW Δ0.1.
To close: init_crown_ratios! must compute its CCF on the backdated dbh — a circular dependency
(crown→CCF→backdated dbh→projected DG→crown) FVS resolves in CRATET; jl HAS the backdating pass
(diameter_growth.jl:295-313) but it only recomputes the DG's PCT, not the crown ICR. Wiring the crown-ICR
estimate into the backdated-dbh pass is the close. The DDW CARBON MODEL is bit-exact (carbon_snt).

#### ★ SPEC MET — carbon_jenkins DDW closed; only broken test is the accepted COMPRESS
The carbon_jenkins init-crown was closed: init_crown_ratios! now computes the CCF on the DENSE-backdated dbh
(past dbh = sqrt(d²·r), r from the measured DG via diam_growth + bark_ratio; stand-average for unmeasured),
and passes it to crown_ratio_update! as relden_override (current dbh/SDI/rank for everything else, exactly as
FVS CROWN). Init crown is now [22,27,32,36,40,47], BIT-EXACT vs FVS, so cycle-1 LP DGF is bit-exact and
carbon_jenkins DDW matches FVS to ULP every cycle (3.8/2.5/3.8/8.0). The carbon_jenkins driver + .out-writer
DDW @test_broken flipped to passing.

FINAL SUITE: 4519 pass / 0 fail / 1 broken — the one broken is s22_compress (the accepted different-
eigensolver divergence). The drop-in spec is met: the only divergences are ULP FP and the COMPRESS
eigensolver. The full FFE carbon down-wood model + the CRATET init-crown are bit-exact drop-ins.

## init_crown_ratios! — post-hoc source verification (honest re-audit)

Verified the CRATET init-crown backdating against FVS source (not just test match):
- dense.f:257-264: LBKDEN=TRUE DENSE ends with RELDEN=RELDM1 = first-pass BACKDATED RELDT (dense.f:239).
- cratet.f:186 backdated DENSE → :532 CALL CROWN → :613 non-backdated DENSE. CROWN reads backdated RELDEN. ✓
- Backdating algo (dense.f:89-132: BAGR + WK3=SQRT(D²·R), R=BAGR fallback) matches the jl replica. ✓ FAITHFUL.

KNOWN GAPS — ✅ BOTH RESOLVED (this entry was STALE; written before the #14 cleanup extracted the shared
`_backdate_dbh!`, which gates on IDG). Verified 2026 vs live FVS:
1. IDG==1 (past-dbh input): RESOLVED. `init_crown_ratios!` now delegates to `_backdate_dbh!`
   (diameter_growth.jl), which does `gadj = idg==1 ? g : g/bark_ratio(...)` — skipping the bark division for
   IDG==1 exactly as dense.f:101,122. VALIDATED vs live on growth_idg1 (IDG=1, 6 trees all with missing
   crowns ⇒ init_crown fires, all 6 with measured DG ⇒ the idg==1 backdating runs): init crown
   [22,27,31,36,40] and the .sum TPA/BA/TopHt are BIT-EXACT vs live (a wrong init crown would diverge the
   crown→growth path; it doesn't).
2. G==0 exactly: NON-ISSUE. CRATET sets input-zero DG to −1.0 BEFORE DENSE (dense.f:95 comment), so G==0
   never reaches the BAGR loop; jl's `g<=0` exclusion gives the SAME outcome as FVS's convert-to-−1-then-
   `g<0`. (`_backdate_dbh!` also uses `g<0` for IDG==1/3 to match dense.f:100 exactly.)

## CONFIRMED BANDAID: fix #6 snapshot_ffe_oldcrown! at inventory (instrumented)

Source: fmsdit.f:93 `IF (ICYC.GT.1) THEN` guards the ENTIRE crown-lift OLDCRW scaling — cycle 1 is
explicitly skipped ("If we are on the first cycle, then the old crown is not known"). fminit.f:971
sets OLDCRL=0. fvs.f never calls TREGRO at ICYC=0.
Instrumented live FVSsn (carbon_snt), per-cycle annual crown-lift material per acre:
  ICYC=1 (1990): 0.00   ICYC=2 (1995): 947.2   ICYC=3 (2000): 1103.7   ICYC=4 (2005): 1166.2
=> FVS crown-lift down-wood is ZERO in cycle 1. FVSjl's summary.jl:101 inventory snapshot makes the
first compute_crown_lift! (end of loop c=0) nonzero (inventory→cycle1 rise), which FVS does NOT do.
MECHANISM UNFAITHFUL. The code comment "else the 1st cycle's fine down-wood is lost" describes
down-wood FVS deliberately never creates. PENDING DECISION: replicate the ICYC>1 guard (gate the first
crown-lift to cycle≥2 / drop the inventory snapshot) and find what real gap, if any, it was masking.

## LEDGER (this session's DDW fixes) — 6/7 cited-faithful, 1 confirmed bandaid
1 fall keyed on initial state — fmcwd.f:149-150 ✓
2 snag bole=MAX(X,MCF) merch — fmsvol.f:146-150 (SN) ✓
3 SCNV 0.80/1.00 by pool — fmcwd.f:61,404-414 ✓
4 decay order FMSNAG→FMCWD→FMCADD — fmmain.f:228-241 ✓
5 crown small-tree bole=merch — fmsvol.f:80,146-150 ✓
6 inventory crown snapshot — CONTRADICTED by fmsdit.f:93 + instrumented cyc1=0 ✗ BANDAID
7 init crown backdated CCF + IDG gating — cratet.f:186/532, dense.f:100-127,257-264 ✓

## RETRACTION: fix #6 (inventory crown snapshot) is FAITHFUL, not a bandaid

The earlier "CONFIRMED BANDAID: fix #6" section above is WRONG and is retracted. Corrected trace:
FVS applies growth in UPDATE (update.f:65 HT=HT+HTG, :115 DBH=DBH+DG/BRATIO), called at gradd.f:180 AFTER
FMMAIN (gradd.f:118). So FMOLDC (fmmain.f:268, every cycle, ungated) snapshots the PRE-growth crown; cyc1's
FMOLDC captures the INVENTORY crown. jl's inventory snapshot_ffe_oldcrown! is the faithful analog of FMOLDC(cyc1).
The fmsdit.f:93 IF(ICYC.GT.1) gate only suppresses cyc1's FMSDIT (lift computation), NOT cyc1's FMOLDC (state
capture). jl mirrors this exactly: cyc1 fuel loop adds the zero-init array (no cyc1 lift), and the inventory-based
lift(inv→postcyc1) is computed at end of loop c=0 and ADDED during cyc2 — matching FVS's ICYC=2 FMSDIT lift(947).
My earlier instrumentation correctly showed FVS FMSDIT cyc1=0 but I mis-mapped it onto the jl loop. Corroboration:
carbon_snt DDW is bit-exact vs the LIVE-FVS golden WITH #6; removing #6 would zero the cyc2 lift and break that
live-FVS match. CONCLUSION: keep #6. Lesson: bit-exact-vs-LIVE-FVS + a correct both-sides trace = faithful; a
static timing argument alone misled both me and the confirm agent.
