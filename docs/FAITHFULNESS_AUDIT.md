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
