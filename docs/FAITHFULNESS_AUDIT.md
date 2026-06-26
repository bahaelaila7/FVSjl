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

### F3 — Rothermel fire-behavior (byram/flame/spread) divergence (FFE). OPEN, characterized.
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
