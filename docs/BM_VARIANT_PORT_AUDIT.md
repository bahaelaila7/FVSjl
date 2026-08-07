# BM (Blue Mountains) Variant Port — Audit Log

FVS Region-6 Blue Mountains variant → FVSjl. Validated bit-exact-or-cornered per chunk vs the live
relinked oracle `/workspace/.bmwork/FVSbm_clean` (gfortran-16 + isoc23 shim). MAXSP=18, Stage SDI,
DGSD=1.5, 10-yr cycle, seed 55329. Test stand `bmt01` (S248112, forest 614 Umatilla, habitat 12).
Growth validation uses `bmt01_growth.key` (FFE keywords stripped — bmt01.key is a full FMIN/SIMFIRE demo).

## Chunk verdicts

| # | Chunk | Files | Verdict |
|---|-------|-------|---------|
| 0-1 | Scaffold + species | bluemountains.jl, species.jl | ✓ loads, 18 species, seed 55329, POWER bark |
| 2 | Site / habitat | site_index.jl | ✓ **bit-exact** — ISISP=DF(3), SITEAR(DF)=64, SDIDEF=346, all 18 SITEAR match live SITECODE echo. Novel eco-class DB chain: habitat 12 → PCOML[12]=CDS722 → ECOCLS → SICHG/HTCALC (incl. the faithful CASE(3) `**` precedence bug) |
| 3 | Large-tree DG | diameter_growth.jl, dg_coefficients.jl | ✓ **bit-exact-or-cornered** — DGCON bit-exact (DF/WL/LP 1.5531/1.3957/1.2295), WK2 24/27 bit-exact (MSS spline + all 3 species groups), DGSCOR COR bit-exact (GF 0.01963=live). Misses = 3 D=0.08 seedlings |
| 4 | Height | height_growth.jl | ✓ **21/24 bit-exact** — WL/DF/LP/ES exact; bm_findag + POTHTG=HTCALC(SITAGE+10)−SITHT + HGMDCR/HGMDRH modifiers |
| 5 | Crown | crown.jl | ✓ Weibull (= validated TT/UT form) + bm_tree_ccf (RELDEN, 369/369 per-tree CCF bit-exact) |
| 6 | Regent (small-tree) | regent.jl | ✓ bm_smhtgf height matches live (POTHTG 11.42, HTG 5.72); ht-dbh uses bm_htdbh (Curtis-Arney) for LHTDRG=false species |
| 7 | Mortality | mortality.jl | ✓ **faithful** — Hamilton RI + SDI self-thin; RN bit-close for identical input (525.769/0.0114 vs 525.759/0.011), SDIMAX=346 both |
| 8 | Volume | volume.jl | ✓ FW2W Flewelling conifers (total cubic cornered, bark-fixed); **616BEHW COMPLETE** — total cubic (bm_r6vol3 + bm_formcl) + merch cubic + board (bm_r6dibs + bm_r6vol1) all bit-exact vs live (per-tree fort.84 6/6; pure-WP .sum MCuFt 1068/1068 + BdFt 5341/5341 exact, TCuFt ±1) |

## Real bugs found + fixed (all via measurement vs the live binary)

1. **Species crosswalk** — `species_translation.csv` was an EM placeholder; bm/spctrn.f selects ASPT col 5 (BM), not 10 (EM). Grand-fir/white-fir trees (18%+15% of bmt01) were mis-mapped to AF/OS. Fixed → RELDEN 109.10 bit-exact.
2. **Calibration bark** — shared `_backdate_dbh!`/`calibrate_diameter_growth!` had CR/TT bark branches but no BM (POWER model); fell to the linear bark_ratio. Added _bm_bd/_bm_cal.
3. **DIB/D bark ratio** — `bm_bratio` computed DIB (BARK1·D^BARK2), not DIB/D → bark ≈ 0.999 instead of ~0.9. **One fix unmasked FOUR "cornered" residuals**: volume +18%→−4.5%, growth BA 96→98, DG WK2 18/27→24/27, and the GF DGSCOR COR (now bit-exact).
4. **HTDBH dispatch** — regent used Wykoff AX/BX ht-dbh for all conifers, but bm/regent.f:522 uses bm/htdbh.f (Curtis-Arney, forest-dependent) when LHTDRG=false (all except WJ/WB/LM/AS). Ported bm_htdbh + dispatch. Faithful but .sum-neutral (seedling DBH is BA-negligible).

## Hypotheses rejected by measurement (would have been cargo-culted fixes)

- FW2 geosub-substitution as the volume cause (was the DIB/D bark bug)
- GFSUB F-coefficient substitution (implemented → regressed → reverted per doctrine #4)
- Mortality kill-distribution bug (RN math matches live for identical input)
- Regent seedling height over-growth (height matches live)
- Regent seedling DBH as the divergence driver (HTDBH fix is .sum-neutral)
- GF DGSCOR COR as "accepted cornered" (was the bark bug — now bit-exact)

## Remaining residuals (measured, precisely attributed)

- **Multi-cycle divergence** (2090 TPA 163/96): the DGSCOR-FRM serial-correlation (dgdriv.f:271, WK2 bit-exact but DG-from-WK2 is ZZRAN-affected) + tripling, amplified by the faithful mortality feedback. The accepted cornered stochastic class doctrine #3 forbids chasing per-record.
- **Volume total-cubic** ~cornered: INGY SF_SHP form precision (jsp-13 coefficients bit-identical to NVEL source; GFSUB substitution ruled out).
- **Deferred leaves** (scoped, both substantial — not bounded polish):
  - **616BEHW minor species** (WP/MH/WJ/WB/LM/PY/YC/AS/CW/OS/OH — not in bmt01) — ✅ **COMPLETE**.
    **MERCH + BOARD LAYER DONE**: ported `bm_r6dibs` (r6vol3-companion log-bucking, ZONE-1/A=0.62/
    total-height path, Behre taper bucked to MTOPP=TOPD·bark) + `bm_r6vol1` (per-log Scribner from the
    132-entry IFTR table + merch-cubic butt-log/Smalian; INTL14 International skipped = VOL(10) not in
    .sum) + accumulation VOL(2)=Σ ANINT(Scribner) [BdFt], VOL(4)=Σ round(merch·10)/10 [MCuFt]. Gated:
    merch/board=0 for h≤17.3 (R6VOL returns cylinder, skips R6DIBS/R6VOL1 — the bug that first made
    merch>total). Wired MCF=VOL(4) (D≥DBHMIN), BdFt=VOL(2) (D≥BFMIND, METHB≠9). **Validated bit-exact**:
    live-fvsvol.f instrument-replay (fort.84) 6/6 trees VOL(2)+VOL(4) exact; pure-WP .sum 1990
    MCuFt 1068/1068 + BdFt 5341/5341 exact (TCuFt ±1); projected cornered = growth-tail (= total cubic).
    ★ Multi-species hardening (distinct bark + form-class paths): pure-MH (sp5, power bark, FIAJSP FC 264)
    = ALL THREE bit-exact 1507/1001/5048; pure-WB (sp11, constant bark 0.969, default FC 80) = BdFt 5479
    exact, MCuFt ±1, TCuFt ±4 (INGY cornered). Confirms the leaf across bark groups + in-list/default FC.
    ── (history) **TOTAL CUBIC LAYER**: ported `bm_r6vol3` (r6vol3.f Behre Smalian taper
    DR=HRATIO/(0.62·HRATIO+0.38), D17=FCLASS/100·DBHOB, H17=17.3, all 3 branches: DBHIB<TOPD cylinder /
    D17<TOPD two-log / full taper) + `bm_formcl` (formclas.f FORMCL_BM form-class lookup: FIAJSP binary
    search, IFCDBH=(D−1)/10+1, 4 forest tables MALH/OCHO/UMAT/WLWH by IFORST=KODFOR%100) + the R6VOL
    short-tree guard (TTH≤17.3→cylinder) + glue (spec=VEQNNC(8:10), DBTBH=D·(1−bm_bratio) per fvsvol.f:153).
    **Validated**: bm_r6vol3 **bit-exact vs live R6VOL3 (96/96** grid pts, all 3 branches, drv_r6vol3.f);
    bm_formcl **bit-exact vs live FORMCL_BM (120/120**, drv_formcl.f); synthetic pure-WP stand end-to-end
    **TCuFt jl 1556 vs live 1557 (±1**, INGY-cornered class), growth cols bit-exact.
    **REMAINING sub-leaf** — merch cubic + board foot (data flow now fully traced, r6vol.f):
    live 616BEHW yields MCuFt/BdFt (pure-WP 1990: 1068/5341) via the full R6VOL path:
    - **R6DIBS**(IAPZ=ZONE,DBHOB,BTR,FCLASS,MTOPP,TLH,TH, → XLOGS,LOGDIA,SL,XL,A): 297-line log-bucking.
      For our case (ZONE 1, HTTYPE='F' ⇒ TH>0/TLH=0) the live path is label 70→80→130 (16.3-ft logs,
      Behre taper DR=HR/(A·HR+B), MTOPP top). Fills NOLOGP(=XLOGS), LOGDIA(21,3) small-end diams,
      SL/XL(20) scaling+actual lengths. ⚠ heavy computed-GOTO + an IRET-dispatched label-1000 taper-A
      setup (unread) + ZONE-2 (32-ft) paths (labels 200-500, not needed for 616).
    - **R6VOL1**(ZONE,DBHOB,FCLASS,NOLOGP,LOGDIA, → LOGVOL,INTBF): 90 lines. Per-log Scribner board
      LOGVOL(1,·), merch cubic LOGVOL(4,·), International INTBF(·).
    - **Accumulate** (r6vol.f:170-186): VOL(2)=Σ ANINT(LOGVOL(1,I)) [Scribner board=.sum BdFt col],
      VOL(4)=Σ round(LOGVOL(4,I)·10)/10 [merch cubic; MCF=VOL(4)+VOL(7), VOL(7)=0], VOL(10)=Σ INTBF.
    - **Wire**: compute_volumes_bm! BEHW branch sets merch_cuft_vol=VOL(4), bdft_vol=VOL(2) (currently 0).
    ~390 NVEL lines + accumulation. bm_r6vol3/bm_formcl already supply the total-cubic core.
    Validate: standalone R6DIBS/R6VOL1 drivers (as with drv_r6vol3/drv_formcl) then pure-WP end-to-end
    MCuFt/BdFt vs FVSbm_clean. Low-impact (species absent from bmt01), substantial — a fresh-session chunk.
  - **FFE** — bmt01.key is a full FMIN/SIMFIRE/PotFIRE/FuelOut demo; the fire/fuel/snag/carbon
    subsystem (needs BM biomass + fuel coefficients). Large; validate vs the full bmt01.key .sum.

## Verdict

BM growth core (site + DG + COR + height + crown + regent + mortality) is **bit-exact/faithful chunk-by-chunk**
against the live Fortran; volume bark-fixed with a cornered form-precision tail. Every measurable deterministic
layer is validated; the sole remaining residual is the DGSCOR-FRM/tripling/ZZRAN stochastic class. BM joins
EM and UT as a ported, validated core western variant.

## 2026-08-05 — #140 BM growth-only under-thinning: MECHANISM localized (self-thin target too high, feedback-amplified)
bmt01 first stand (control, NOAUTOES, 10 cyc) vs live FVSbm: jl UNDER-thins badly by 2090 — TPA 163 vs live 96
(+70%), BA 165 vs 146, QMD 13.6 vs 16.7. jl holds **SDI constant at 267** (2030-2090) while live's SDI **declines
268→219** (constant BA=146). MEASURED + RULED OUT three hypotheses:
- **MSB (mature-stand breakup)**: NOT it — bm/grinit.f:284-286 defaults QMDMSB=999, SLPMSB=0 ⇒ MSB is OFF by
  default (fires only via keyword; QMD never reaches 999). This also explains why task #144's SLPMRT/CEPMRT attempt
  was falsified — the MSB/MRT alternate-mortality path is simply inert for the default BM stand.
- **RIP blend**: NOT it — despite the "weighted average of RI and RN" comment (morts.f:81-82), the actual code
  (morts.f:508-517) is EITHER/OR (RIP=RN unless T≤TEM or RN≤0 → RIP=RI), which jl matches exactly.
- **SDImax decline**: NOT it — instrumented live morts.f SDIMAX dump: constant **346.00** every cycle, == jl.
ROOT (localized): jl's mortality self-thin TARGET tn10 (=T85D10=CONST·D10^−1.605·PMSDIU) is consistently HIGHER
than live's because **jl's projected dq10 (mortality QMD) is LOWER** — icyc1, SAME trees (T=589.65 both): jl
dq10=6.02 vs live d10=6.10 (−1.3%). Lower QMD ⇒ higher TMD10 ⇒ higher T85D10 ⇒ fewer trees killed. Then BM's
self-thin **QMD-feedback loop** amplifies it: under-thin → more (smaller) trees → lower QMD → higher target →
under-thin, compounding the ~1% icyc1 seed into +70% TPA by 2090 (jl tn10 vs live: 530.8/525.8, 490.9/471.1,
433.1/409.9, 367.8/334.5 — diverging).
NEXT (definitive settle, not yet done): the icyc1 dq10 difference is on the SAME trees ⇒ it is a
projected-diameter (d+g) difference. dq10=sqrt(Σpr·(d+g)²/tt). Dump jl vs live per-tree (d, g=diam_growth/bark) at
icyc1 and diff — if g matches and only the SUM differs, it's an aggregation/order (tie-break) artifact (cornered-
but-amplified, like the EM over-growth tail); if jl's g is systematically ~1% low, it's a real BM DG bias feeding
the mortality (fixable). Either way BM's self-thin QMD-feedback is a strong amplifier that turns a ~1% seed into a
large visible divergence — that sensitivity itself may warrant a damping review. (Live oracle /workspace/.bmwork/
FVSbm_clean; both sides restored clean.)

### SETTLE (same day) — #140 seed is the tie-break class (like EM), amplified by BM self-thin feedback
Ran the cyc0-DG test on BM (NOTRIPLE bmt01, jl projected d2000 at simulate.jl:533 vs live BM dgf D@ICYC=2).
sp2 projected 2000 DBH: jl {8.82, 8.93, 9.55, 9.65} mean 9.24 vs live {8.90, 9.36, 9.39, 9.44} mean 9.27. ⇒
per-tree DG diffs are ~0.2-0.4″ MIXED-SIGN (jl more spread), the SAME crown-ratio/BA-percentile tie-break
signature as the EM cyc0 test, netting to a small −0.3% mean. That −0.3% lower projected mean ⇒ lower dq10 ⇒
higher T85D10 self-thin target ⇒ under-thin, which BM's self-thin QMD-feedback (faithful to live's identical
self-thin form) amplifies into the +70% TPA by 2090. VERDICT: #140's SEED is the accepted tie-break precision
class (irreducible, same as EM/IE); the LARGE visible magnitude is BM's self-thin feedback amplifying it — the
feedback itself is a faithful port of live's model, so there is no wrong equation to fix. CAVEAT (why not a clean
"cornered"): unlike EM (net ~0, aggregate BA bit-exact), BM's cyc0 net is −0.3% and gets amplified 200× — so
whether the −0.3% net is pure tie-break realization (cornered) or carries a tiny systematic component needs a
MULTI-STAND check (does the sign/magnitude of the net vary stand-to-stand, or is jl consistently low?). If multi-
stand shows jl consistently ~0.3% low ⇒ a small real projected-DBH bias worth hunting (crown/PCT); if it straddles
⇒ fully cornered. That multi-stand net-bias check is the one remaining measurement for #140.

### #140 multi-stand check — harness note for next session (2026-08-05)
The remaining net-bias check (does jl's cyc0 −0.3% projected-DBH net straddle across BM stands → cornered, or is it
consistently jl-low → fixable) needs MULTIPLE treed BM stands. tests/FVSbm has only bmt01 (one stand, S248112), so
this requires FIA stands. Friction hit this session (documented so it isn't repeated): (1) fia_sweep_check-style
scripts need `Dict(y=>v for (y,v) in J)` not `Dict(J...)`; (2) background julia via `nohup … &` inside a Bash tool
call dies silently — use Bash `run_in_background: true` on the julia command directly; (3) a stratified N=20 BM
sample risks being treeless-heavy (EM's N=20 cluster sample was 0-treed) → extract N≥60 and/or filter to treed
conditions first. RELIABLE recipe: `extract_sample.jl BM 80` → run each treed stand 3 cycles jl-vs-live, tally the
sign of the final-cycle TPA delta. If JL-HIGH dominates ⇒ small real projected-DBH bias (hunt in crown/PCT feeding
dq10); if it straddles ⇒ #140 fully cornered like EM/IE. The mechanism + tie-break seed are already established
(this file, prior entries); only this sign-tally remains.

### #140 multi-stand check — INFRASTRUCTURALLY BLOCKED in this environment (2026-08-05)
Attempted the multi-stand sign-tally three ways (N=20 FIA, N=80 FIA, via nohup and via proper background) — all
failed to complete: backgrounded `julia` processes die SILENTLY here (main process gone with 0 output), Julia
BUFFERS stdout when redirected to a file (so any partial per-stand results are lost on the silent death, not
flushed), and each killed attempt LEAKS orphaned `timeout`/FVS children that then contend for CPU and slow the next
attempt. Net: the live-FVS-per-stand multi-stand sweep is not runnable-to-completion in this session's environment.
The single-stand measurements that DID complete (bmt01 cyc0-DG test, full-precision) are what established the #140
mechanism + tie-break seed; only the cross-stand sign-tally (cornered-vs-tiny-bias) is blocked. To finish it in a
stable environment: run `bm_dir.jl` (fixed) on `extract_sample.jl BM 80` synchronously (foreground, no redirect
buffering issue) or with per-stand `flush(stdout)` after each println so partial results survive; tally the
final-cycle TPA-delta sign across treed stands. All BM state restored clean (oracle 2090 96/146; repo clean).

### #140 CORRECTED VERDICT (2026-08-05) — REAL consistent under-thin bias, NOT cornered
The multi-stand sign-tally FINALLY ran (foreground-auto-bg + per-stand flush + no grep pipe — the recipe that
survives this env's background quirks). Result on real BM FIA stands (35-stand subset, 3-cycle jl-vs-live TPA):
- **Non-self-thinning stands: bit-exact** (≈EQ, Δ=0.0% — e.g. 18/179/482 TPA stands jl==live exactly).
- **Actively-self-thinning stands: jl UNDER-THINS, 100% consistent direction** — 10/11 divergent stands JL-HIGH
  (Δ = +1.3%, +2.4%, +7.2%, +11.6%, **+48.3%**), just **1 borderline JL-LOW (−0.9%)** — a ~10:1 under-thin skew, not a balanced straddle.
⇒ This **CORRECTS the earlier "cornered tie-break" lean** (SETTLE entry above), which over-generalized from bmt01's
single −0.3% cyc0 net. The consistent JL-HIGH direction proves the −0.3% projected-DBH net is **consistently
signed, NOT mixed-sign tie-break** — i.e. #140 is a **real, systematic under-thin bias** that manifests whenever
BM self-thinning is active (RIP=RN), amplified by the QMD-feedback (small on most stands, large on dense/rapid-
self-thin stands like bmt01 +70% by 2090 and FIA stand 248913820489998 +48% by cyc3). META (doctrine #2, again):
a single-stand net can look like tie-break noise; only the MULTI-STAND sign-tally distinguishes cornered-straddle
from a consistent bias — and here it flipped the verdict. This is why the measurement was worth the harness fight.
ROOT (for the fix): jl's self-thin kills fewer trees than live on active-self-thin stands. From the cyc0-DG test,
jl's mortality-input dq10 runs slightly low ⇒ higher T85D10 target (∝ dq10^−1.605) ⇒ under-kill. Next: on one
JL-HIGH FIA stand (e.g. 374430545489998), instrument jl BM mortality! vs live bm/morts.f at the FIRST divergent
cycle — dump dq10, T85D0/T85D10, TN10, RN, and the per-tree kill — to localize whether the gap is the projected
dq10 (DG/bark into g), the TN10 target formula, or the RN→WKI kill application. The self-thin QMD-feedback then
compounds it; fixing the per-cycle under-kill closes #140.

### #140 ROOT NAILED (2026-08-05) — jl mortality reads a WRONG (partial) diam_growth ⇒ dq10 too low ⇒ under-thin
Instrument-replay on repro FIA stand **374430545489998** (a JL-HIGH +7.2% stand), cyc1 mortality, jl vs live
bm/morts.f, SAME trees (DQ0 bit-exact 5.821):
- Projection formula IDENTICAL: live `G=(DG(I)/BARK)*(FINT/10)`, `CIOBDS=2*D*G+G²`, `SD2SQ+=P*(D²+CIOBDS)` (morts.f
  :222-224); jl `g=diam_growth/bark`, `sd2sq+=pr*(d²+2dg+g²)` (bluemountains/mortality.jl:24-26). FINT=10 ⇒ the
  `(FINT/10)` factor is 1. Both use bm_bratio/BRATIO — **bark MATCHES** (jl 0.8621 vs live-implied 0.862).
- PER-TREE growth term: **live DG(I) = 0.575–0.589″** (G≈0.667″ outside-bark) vs **jl t.diam_growth = 0.156–0.161″**
  (g≈0.182″) — jl's mortality growth is **~1/3.7 of live's**. (Small tree i=6 d=5.1: live G 0.667 vs jl g 0.391.)
- CONSEQUENCE: jl dq10 = sqrt(SD2SQ/T) comes out **5.908 vs live 6.078**; the lower projected QMD keeps
  T=414.6 **below the self-thin threshold** (t55d10≈505) ⇒ jl RN=0 (background only) ⇒ TPA 403; live's higher
  dq10 crosses the threshold ⇒ self-thins RN=0.0067→0.0096 ⇒ TPA 376. That is the +7.2% under-thin, and (compounded
  by the QMD-feedback) the +48% on dense stands and +70% on bmt01.
★ THE BUG: jl's **applied** DG is correct (BA bit-exact ⇒ dbh+=diam_growth/bark applies ~0.575″), but the value in
`t.diam_growth[i]` AT THE MORTALITY READ POINT is only ~0.156″ — a PARTIAL/pre-final value. So `t.diam_growth` is
being read by mortality! before it holds the full cycle DG (or BM stores a pre-scaled/pre-converted increment that
GRADD later finalizes). FIX: make jl's BM mortality use the SAME full DG live uses — either (a) reorder so
mortality reads the finalized diam_growth, or (b) have mortality apply the same scaling/conversion GRADD does, or
(c) BM should use dg_prev/WK1 semantics like KT/IE/TT (check bm/morts.f DG(I) provenance vs jl's t.diam_growth).
NEXT SESSION: dump jl t.diam_growth at the END of the cycle (post-GRADD) for these same trees — confirm it becomes
~0.575″ — then trace where the ~1/3.7 partial value at mortality-time comes from (dgf output? subcycle? FINT). This
is a REAL, high-value fix: it closes #140 and likely tightens BM self-thinning cluster-wide.

### #140 ROOT REFINED (2026-08-05, post-restart) — it's the MORTALITY growth-term projection, NOT applied DG
Correcting the prior "jl mortality reads partial diam_growth" framing with the full trajectory on repro stand
374430545489998 (jl vs live FVSbm):
| year | jl TPA/BA/SDI/QMD | live TPA/BA/SDI/QMD |
| 2015 | 415/77/174/5.8 | 415/77/174/5.8  → **BIT-EXACT** |
| 2025 | 403/77/173/5.9 | 376/80/177/6.2  → jl under-thins (+7.2% TPA) AND under-grows (QMD 5.9 vs 6.2) |
Since 2015 is bit-exact, the **APPLIED DG matches** — this is NOT a DG-application bug. The divergence enters at
the 2015→2025 cycle. KEY per-tree measurement (cyc0 mortality, jl `t.diam_growth` == same at mortality AND apply):
jl mortality growth-term **g ≈ 0.156–0.161″ (constant across d=8–27″)**; live mortality **G = 0.6667″ (also constant
across sizes)**; jl g ≈ 1/3.7 of live G. Both being SIZE-INDEPENDENT constants ⇒ each side hits a bound/floor/default
in the mortality's QMD projection, and they differ. Consequence: jl dq10 = 5.908 vs live 6.078 ⇒ jl's projected QMD
keeps T=414.6 BELOW the self-thin threshold (t55d10≈505) ⇒ jl RN=0 (background only) while live crosses it and
self-thins ⇒ the +7.2% (and, feedback-amplified, +48% / +70%). At the FIRST cycle both are below threshold ⇒ the
g/G difference is inert ⇒ 2015 bit-exact; it only bites at the threshold-crossing cycle.
⇒ **THE BUG IS IN THE MORTALITY's dq10 GROWTH-TERM**, not the applied DG. live morts.f:222 `G=(DG(I)/BARK)*(FINT/10)`
uses a DG(I) that is ~3.7× jl's `g=diam_growth/bark` — and since applied DG matches, live's mortality DG(I) is a
DIFFERENT (un-reduced/projected) quantity than the applied (DGBND-reduced) DG. NEXT: instrument live to dump BOTH
DG(I)@MORTS and the applied per-tree DBH increment for the same tree — confirm live's MORTS DG(I) is the
pre-DGBND-cap DG while the applied is post-cap; then jl's fix = have BM mortality project dq10 from the un-capped
DG (store it alongside the capped diam_growth, or recompute), NOT from the DGBND-reduced t.diam_growth. This closes
#140 (a real mortality-projection bug) and should tighten BM self-thinning cluster-wide.

### #140 ROOT — FINAL, VERIFIED (2026-08-05): mortality dq10 uses REGENT-reduced DG; must use large-tree POTENTIAL DG
Step-by-step trace on repro stand 374430545489998, cyc0, tree i=1 (d=8.1"):
- `diameter_growth!` raw `dgc` = **0.6729** (large-tree DG, un-capped: bounded==raw, so _bound_scale/DGBND is NOT
  the culprit; and it varies with size: d16.5→1.03, d27.8→1.15 — a real large-tree DG).
- AFTER `small_tree_growth!` (simulate.jl:470): `t.diam_growth[1]` = **0.1567** ⇐ REGENT reduces it.
- AFTER `apply_fix_scalers!`: 0.1567 (unchanged). mortality reads 0.1567; apply uses 0.1567.
- Live `bm/morts.f` MORTS uses DG(I) ⇒ G=0.667″ (≈ the large-tree potential DG 0.575″, NOT the reduced value).
⇒ **THE FIX (verified mechanism):** jl's BM mortality `g = t.diam_growth/bark` reads the REGENT-reduced DG (0.156),
but FVS's mortality projects `DQ10 = "QMD at end of cycle if TPA held constant"` (morts.f:62) using the LARGE-TREE
POTENTIAL DG (~0.6″), not the REGENT-reduced applied DG. The REGENT reduction is correct for the *applied* growth
(2015 .sum BIT-EXACT ⇒ jl and live apply the same reduced DG), but it must NOT feed the mortality's dq10. Because
jl's dq10 is built from the reduced 0.156, it comes out low (5.908 vs 6.078), keeping the stand below the self-thin
threshold ⇒ jl RN=0 (background only) while live self-thins ⇒ the +7.2%/+48%/+70% under-thin (13:1 JL-HIGH multi-
stand). Inert at cyc0 (both below threshold ⇒ 2015 bit-exact); bites at the threshold-crossing cycle.
IMPLEMENTATION: in bluemountains/mortality.jl, compute the mortality `g` from the PRE-REGENT (large-tree) DG — either
(a) snapshot t.diam_growth into a scratch field right after diameter_growth! (before small_tree_growth!) and have
mortality read that, or (b) confirm live's exact DG(I) (0.575 vs jl raw 0.67 — a ~14% gap to reconcile, possibly the
FINT/10 factor or a subcycle detail) and match it. Then re-run the multi-stand sign-tally (foreground+flush recipe,
extract_sample.jl BM 80) and confirm the JL-HIGH skew collapses toward ≈EQ. This closes #140 and likely tightens
BM (and possibly EM/UT which share the western small-tree+mortality structure) self-thinning.

### #140 COMPLETE ROOT (2026-08-05) — TWO COUPLED components; explains why #144 was falsely reverted
Implemented + tested the mortality-DG fix (snapshot pre-REGENT DG into dg_prev, BM mortality reads it). VERIFIED it
works AS FAR AS dq10: jl dq10 5.908 → **6.249** (now ≥ live's 6.078), dg_prev = 0.67/1.03/1.0 (correct large-tree
DGs). BUT the repro stand's TPA was UNCHANGED (still 403, not live's 376) — because jl STILL doesn't self-thin:
branch dump `tt=414.6, t55d0=517.7, t55d10=462.0, t85d10=714.1 → tn10=tt=414.6, rn=0` (tt < t55d10 ⇒
bluemountains/mortality.jl line 44 `elseif tt<=t55d10: tn10=tt` ⇒ NO self-thin). LIVE self-thins to **tn10=387.8,
rn=0.0067** even though tt(414.6) < t55d10 — because live reaches it via the **MRT path** bm/morts.f:434
`TN10=EXP(CEPMRT+SLPMRT*TEM)` (SLPMRT defaults to SLP at morts.f:430, so it's ACTIVE), which jl DOES NOT HAVE.
⇒ ★★ **#140 = TWO COUPLED FIXES:** (a) mortality dq10 must use the pre-REGENT large-tree DG (verified above), AND
(b) the CEPMRT/SLPMRT MRT self-thin path (jl mortality lines 39-48 lack it; it lets a stand below t55d10 still
self-thin). This EXPLAINS TASK #144 ("SLPMRT/CEPMRT implemented+measured+REVERTED — FALSIFIED as root"): #144 tested
the MRT path ALONE, WITHOUT the DG fix — so dq10 was still built from the REGENT-reduced 0.156 ⇒ TEM (=const·dq10^
−1.605·pmsdil, the MRT path's input) was wrong ⇒ the MRT path couldn't produce live's target ⇒ looked falsified.
The two are COUPLED: the MRT path needs the correct (pre-REGENT-DG) dq10/TEM. NEXT SESSION: re-implement BOTH
together — (a) the dg_prev snapshot (simulate.jl before small_tree_growth! + BM mortality reads dg_prev), (b) the
CEPMRT/SLPMRT branch (port bm/morts.f:420-450, SLPMRT=SLP default, TN10=exp(CEPMRT+SLPMRT·TEM) capped at T85D10) —
then validate with the multi-stand sign-tally (13:1 JL-HIGH should collapse). The DG fix was REVERTED here (inert on
this stand alone, untested for regressions without the MRT half). Repo clean, oracle clean.

### #140 REAL ROOT — CORRECTED & FIXED (2026-08-06): PVREF6 habitat crosswalk → SDIMAX (NOT DG/MRT-path)
★★ The 2026-08-05 "TWO COUPLED FIXES (DG + CEPMRT/SLPMRT MRT-path)" conclusion (commit ca37c27) was **WRONG** —
a 5th mis-diagnosis from not measuring SDIMAX. Decisive dual-side instrument-replay on repro stand 374430545489998
(FVSbm_trc: WRITE(16,…) in morts.f RN-point + sdical.f + sitset.f + habtyp.f), cyc0 (2015→2025):

| quantity | jl (before fix) | live | note |
|----------|-----------------|------|------|
| **sdimax** | **395.0** | **205.17** | ← THE ROOT (1.93×) |
| const | 15907 | 8262 | = sdimax/0.02483133 |
| dq10 | 5.908 | 6.078 | tiny diff — NOT the driver |
| t85d0 / t55d0 | 800 / 518 | 415.6 / 268.9 | all ∝ sdimax |
| t55d10 | 505.6 | 250.9 | (the prior "462" was WRONG — measured 505.6) |
| branch | `tt<=t55d10 → HOLD, rn=0` | `T55D0<T≤T85D0`, special-case `|t85d0−tt|=1.0≤5 → tn10=t85d10=387.8, rn=0.0067` | |
| **slpmrt / cepmrt** | — | **0 / 0** | ★ the MRT path was NEVER used by live |

⇒ Live self-thins via the ORDINARY branch (tn10=t85d10), NOT the CEPMRT/SLPMRT MRT path (both 0). The DG-projection
and MRT-path hypotheses were BOTH artifacts of jl's inflated SDIMAX: with sdimax=395, jl's thresholds are ~2× too
high, so tt=414.6 falls below t55d10=505 ⇒ HOLD. The prior session's "branch dump t55d10=462, tn10=tt" was measured
with the SAME wrong SDIMAX (they never dumped SDICAL). MEASURE-don't-infer: 5th wrong root-cause on this bug.

ROOT of the SDIMAX error (traced sdical.f → sitset.f → habtyp.f → pvref6.f):
- SDICAL weights per-species SDIDEF by BA (jl's stand_sdimax matches this formula exactly). The divergence is the
  SDIDEF VALUES: live SDIDEF=**166 uniform** (all 18 species; site-species PP value propagated), raised to 205.17 by
  SDICHK's high-stocking bump. jl SDIDEF averaged **395** because it fell into the **CWG113 default** (SDIMAX PP=395).
- WHY: this real-FIA stand's DB `PV_CODE="CJG111"` is NOT a canonical PCOML code. Live's `bm/habtyp.f` calls
  **PVREF6** (bm/pvref6.f, a 3144-row `(PV_CODE, PV_REF_CODE)→HABPVR` crosswalk) whenever PV_REF_CODE is present:
  `("CJG111","622") → "CPG111"` (row 1438), which HBDECD then matches to KODTYP=49 → ECOCLS SDIDEF(PP)=166. jl's FIA
  reader did an EXACT `findfirst(==(pv), BM_PCOML)` on the raw "CJG111" ⇒ miss ⇒ habitat_code=0 ⇒ bm_sitset! CWG113
  fallback ⇒ SDIDEF=395. (Verified: live HABTYP received KARD2="CPG111" already — the J→P translation is PVREF6.)

THE FIX (commit pending): extracted pvref6.f's three 3144-entry arrays → `data/bluemountains/pvref6.csv`
(1591 non-blank rows), loaded as `BM_PVREF6`/`bm_pvref6` in bluemountains/site_index.jl; fia_database.jl BM branch
now crosswalks (PV_CODE, PV_REF_CODE)→HABPVR before the PCOML match (no PV_REF_CODE ⇒ raw PV code, as live). RESULT
on repro: habitat_code 0→49, pcom CWG113→CPG111, **sdimax 395→205.168 (bit-exact vs live)**, jl self-thins the
NORMAL branch (tn10=387.95 vs live 387.80, rn=0.00662 vs 0.00666), **2025 TPA 403→388** (live 376). No CEPMRT/SLPMRT
path, no DG-snapshot needed — both prior hypotheses RETRACTED.

### #140 SCOPE — corrected by MULTI-STAND before/after (2026-08-06): PVREF6 is MINOR; the systematic root is BMTMRT+IPASS
★ Ran the 80-stand real-FIA sign-tally BEFORE vs AFTER the PVREF6 fix (sub-DB `bm_sub.db`, live oracle FVSbm_clean,
final-year TPA). Result RETRACTS the "PVREF6 = dominant root / closed 56%" claim above:
```
BEFORE (CWG113 fallback): HIGH=30 LOW=3 EQ=23   mean|Δ%|=3.53
AFTER  (PVREF6 fix)     : HIGH=30 LOW=3 EQ=23   mean|Δ%|=3.41
```
- PVREF6 changed only **2 of 56 treed stands** (374430545489998: 403→388 = 7.2%→3.2%; 504368405126144: 3.3%→0.4%).
  The other 78 were byte-identical — they either PCOML-match their PV code directly or lack a PV_REF_CODE, so PVREF6
  never fires. The 56% figure was REPRO-STAND-ONLY; across the population PVREF6 is a real but MINOR fix (the small
  non-PCOML-PV subset), ZERO regressions (both changed stands moved toward live).
- ⇒ The DOMINANT, SYSTEMATIC #140 under-thin is the **30:3 skew that PVREF6 does NOT touch** = the missing MORTS
  **IPASS QMD-convergence loop** (morts.f:578-618: apply kill → recompute surviving DQ10N, small trees die ⇒ QMD
  rises 6.08→6.19, if `|D10−D10N|>0.1 & D10N>DIA0` set D10=D10N & re-derive tn10/rn → 2nd-pass tn10=376.4) AND
  **BMTMRT** (bm/bmtmrt.f, 249 lines: distributes the self-thin kill by PERCENTILE + species shade-tolerance VARADJ,
  not the uniform per-tree RN jl applies). jl does a single uniform-RN pass ⇒ stops at the 1st-pass (higher) tn10.
- Two big outliers unchanged by PVREF6 (449441010497 CDS711/no-ref +60.9%; 248913820489998 CDS624/622 +48.3%) — too
  large for IPASS convergence alone; need per-stand study (SDIDEF-value or dense-regen small-tree, not the crosswalk).
- One NEW crash exposed: 374361232489998 (CES411/622) jl **DomainError** — a separate robustness bug to trace.
⇒ NEXT (the real #140 fix): port the BMTMRT percentile/tolerance self-thin distribution + the MORTS IPASS
convergence loop into bluemountains/mortality.jl (a shared-with-EM/UT self-thin subsystem). PVREF6 lands first as a
correct, self-contained, zero-regression prerequisite (correct SDIMAX is needed for the convergence to target right).

### #140 part 2/2 — BMTMRT + IPASS convergence: full implementation recipe (2026-08-06, scoped not-yet-landed)
The DOMINANT #140 under-thin (30:3 skew, mean 3.4%) = jl's BM mortality! does a SINGLE UNIFORM-RN pass; FVS does
(a) the **MORTS IPASS QMD-convergence loop** (morts.f:578-618) and (b) **BMTMRT** percentile/shade-tolerance kill
distribution (bm/bmtmrt.f). KEY DISCOVERY: jl **already implements both** in the shared
`mortality!(::AbstractVariant)` (southern/mortality.jl:261-389) + `_varmrt!` (the BMTMRT geometric-progression
distribution, byte-identical algorithm) + `_pretzsch_tn10` (which is BM's `_em_tn10_iter` PLUS the CEPMRT/SLPMRT
persistence BM lacks — so the shared tn10 is MORE correct than BM's). **CR (a western variant) already routes
through this shared path** with only a `_varmrt_efftr!(::CentralRockies)` — the exact template for BM.

TWO implementation options:
- **Option B (route through shared, à la CR — preferred long-term):** add `mort_ri_scale(::BlueMountains)=0.5`,
  `_varmrt_efftr!(::BlueMountains)`, a BM branch in the shared `_mbark` (bm_bratio POWER bark — the ONE shared-fn
  change), and set BM's species-CSV `mort_bkgd_intercept`/`mort_bkgd_dbh`/`varmrt_varadj` columns to the correct
  values, then DELETE BM's custom mortality!. ⚠ BLOCKER MEASURED: BM's CSV `mort_bkgd_*` and `varmrt_varadj`
  columns are WRONG (placeholder order/values — e.g. mort_bkgd_intercept[1]=5.9617 not BM_PMSC[1]=6.5112;
  varmrt_varadj=[0.7,0.9,0.5,…] not bmtmrt.f VARADJ=[1,1,1,1,1,1.1,1,1,1,1,.8,.8,.5,.5,1.3,.85,1,1]). Fix the CSV
  first (verify the row↔species-index order). Also verify BM sets zeide_sdi=false (STAGE) + mort_dbh_threshold=0.
- **Option A (BM-local rewrite — lower blast radius):** keep BM's own mortality!, but restructure it to mirror the
  shared body (IPASS loop + `_varmrt!` + `_pretzsch_tn10`) using BM-local `BM_PMSC`/`BM_PMD` + a BM-local
  `BM_VARADJ` const from bmtmrt.f + `bm_bratio`. Reuses only `_varmrt!`/`_pretzsch_tn10` (both variant-generic).

BMTMRT PEFF (for `_varmrt_efftr!(::BlueMountains)`), species order WP..OH: original species (1-5,7-10,17)
`PEFF=(14.94435−0.69929·DBH+0.00868·DBH²)·0.1, DBH>40 ⇒ 0.086`; added species (6,11-16,18)
`PEFF=0.84525−0.01074·PCT+0.0000002·PCT³` (PCT = t.crown_ratio after stand_pct!); clamp[0.01,1]; `EFFTR=PEFF·VARADJ·0.01`
(★ ·0.01, not the SN ·0.1). Shared MORTS flow: SUMTRE = (self-thin) T−TN10 else Σbackground-WK2 → `_varmrt!`.

GOTCHAS to verify before landing: (1) BM must run species_sort! so `s.scratch.idx1`/`sp_count_tab` (the morts.f DO-20
IND1 order) are current — the shared SDI sums rely on it for Float32-order bit-exactness; (2) the persisted self-thin
line reset needs `s.density.tpa_mort` updated each cycle (confirm BM participates); (3) PCT lifecycle — confirm
`t.crown_ratio` holds the BA percentile at BM mortality time (as it does for SN/CR). VALIDATION GATE: bmt01 (needs the
separate run_keyfile establishment DomainError resolved first — see below) + the 80-stand sign-tally (scratchpad
`bm_sub.db`+`bm_live.txt`, expect 30:3 skew → ≈EQ) + SN/NE/CR/CS/LS no-regression (Option B only).

Also OPEN (exposed by the 80-stand sweep, separate from the mortality port):
- **jl DomainError crash** on stand 374361232489998 (CES411/622) — a robustness bug (likely a neg sqrt/log in the
  BM growth or mortality path) to trace + fix (doctrine: crash ⇒ fix).
- **Two big outliers** 449441010497 (CDS711, no PV_REF ⇒ resolves directly, +60.9%) and 248913820489998
  (CDS624/622, +48.3%) — too large for IPASS convergence alone; per-stand study needed (SDIDEF value or dense-regen
  small-tree), likely resolved once BMTMRT+IPASS lands but verify individually.

### #140 part 2/2 — LANDED + VALIDATED (2026-08-06, commit 56d0931): BM routes through shared BMTMRT+IPASS
Implemented Option B (route BM through the shared MORTS driver, à la CR): added `_varmrt_efftr!(::BlueMountains)`
(bmtmrt.f efficiency), `mort_ri_scale(::BlueMountains)=0.5`, a BM branch in the shared `_mbark` (bm_bratio POWER
bark — the only shared-code change, additive/inert for SN/NE/CR), fixed the BM species-CSV mort_bkgd/varmrt_varadj
columns (were placeholder values; the row order was already correct = species-index), and DELETED BM's custom
single-pass mortality!. BM now gets the full machinery it was missing: the BMTMRT percentile/tolerance kill
distribution (via _varmrt!), the IPASS QMD-convergence loop, the CEPMRT/SLPMRT persistence (via _pretzsch_tn10),
and BAMAX/SIZCAP caps.

VALIDATION — 80-stand real-FIA sweep vs FVSbm_clean (sub-DB bm_sub.db), mean|Δ% vs live| at final year, OLD (naive
uniform) → NEW (shared), ALL THREE metrics improved, ZERO aggregate regressions:
```
TPA  3.45 → 0.96      BA  8.65 → 6.45      QMD  6.12 → 5.29
TPA sign-tally: HIGH 31→24, LOW 3→3, EQ 23→30
```
Repro stand 374430545489998: 2025 TPA 388→**377** (live 376, Δ1), BA/SDI/CCF **bit-exact**. Only 1/56 treed stands
regressed on TPA (22404092010497, the pre-existing worst over-thin outlier, −6.1%→−8.2%). The residual BA/QMD
divergence is PRE-EXISTING (identical old vs new — the small-tree-DBH/DGSCOR seedling tail; two BA=0 stands are
1554/170-TPA seedling plots, TPA bit-exact vs live but DBH≈0). ⇒ #140's DOMINANT root is RESOLVED; BM self-thin
mortality is now bit-exact-or-cornered on real-FIA data. The same shared driver applies to EM #137 (verify next).
REMAINING BM tail (all cornered/pre-existing, not mortality): small-tree seedling DBH (regent) + the DGSCOR
compounding tie-break class; the 2 big pre-fix outliers (CDS711/CDS624) now fall under the ≤~1% TPA band.

### Cluster mortality-structure audit (2026-08-06) — generalizing the BM #140 shared-driver fix
The BM #140 fix (route the custom naive single-pass mortality! through the shared BMTMRT+IPASS driver) applies to
any western variant whose mortality! is SDI-self-thin-only. Audited all western custom mortality! for the added-
species-Hamilton hybrid branch (which the shared mortality!(::AbstractVariant) does NOT implement):
- **BM** — SDI-only, POWER bark → FIXED (56d0931).
- **CR** — already routes through the shared driver (the template).
- **UT** — SDI-only (utah/mortality.jl, 74 lines, identical structure to BM pre-fix), ZEIDE SDI (already set),
  ut/utmrt.f distribution, LINEAR bark (⇒ NO _mbark branch needed — simpler than BM). ★ CLEAN BM-ROUTING CANDIDATE.
  Recipe: measure UT real-FIA under-thin (near-certain — same single-pass structure), then _varmrt_efftr!(::Utah)
  from ut/utmrt.f + mort_ri_scale(::Utah)=0.5 + fix UT-CSV mort_bkgd/varmrt_varadj + delete UT custom mortality!.
- **EM** (170L), **TT** (244L), **KT** (152L) — HYBRID: original species SDI self-thin + added species a separate
  KT/IE Hamilton potential-mortality regression (GMULT/REIN/RZ/BAMAX) the shared driver lacks. Need the harder
  hybrid fix (BMTMRT+IPASS on the original-species path, preserve the added-species branch).
⇒ NEXT tractable mortality win: UT (clean, like BM). The hybrids (EM #137/TT/KT) are a separate, more involved
effort. All were validated on their single canonical stand (utt01/emt01/ttt01) which — like bmt01 — does NOT exercise
the BMTMRT+IPASS convergence, so they likely ALL carry the same latent real-FIA self-thin under-thin BM had.

### UT MEASURED (2026-08-06) — already bit-exact-or-cornered; BM routing NOT warranted (corrects the audit hypothesis)
Ran the UT real-FIA sweep (150 stands → ut_sub.db vs FVSut_clean, 0 jl crashes). ★ UT does NOT share BM's under-thin
bias — the "UT likely carries the same latent under-thin" hypothesis above is REFUTED by measurement (doctrine #2):
```
UT treed-stand mortality: 43 compared, JL-HIGH(under-thin)=3, JL-LOW(over-kill)=7, EQ=33, mean|Δ%|=1.16
```
Already bit-exact-or-cornered (33/43 EQ, 1.16% mean). The divergence is a slight OVER-kill lean (the 2 biggest are
−14.5%/−12.5% on dense-regen stands = the accepted small-tree/DGSCOR cornered class), NOT a systematic under-thin.
⇒ Applying the BM shared-driver routing to UT would ADD IPASS convergence (kills MORE) and WORSEN the over-kill
outliers — a REGRESSION. Do NOT route UT through the shared driver. (Why UT differs from BM despite the same
single-pass structure: UT uses ZEIDE SDI vs BM's STAGE, and its sample self-thins less; the naive-single-pass gap
that bit BM is simply within the cornered band for UT's stands.) META: the BM #140 fix is NOT a universal western
mortality fix — it must be MEASURED per variant. UT is DONE on real-FIA mortality; no fix needed. (0 AUTOES-0 stands
in the UT sample — unlike EM, UT's real-FIA slice is genuinely treed.)

## 2026-08-06 — #140 DQ10 seed measurement: confirms the tripling/DGSCOR RNG class (doctrine #3)
Instrumented bm/morts.f per-tree at the SD2SQ accumulation (DQ10=SQRT(SD2SQ/T), SD2SQ=Σ P·(D+G)², G=(DG/BARK)·
FINT/10). Captured live per-tree D/DG/BARK/G at icyc1 (bmt01). KEY REALIZATION: DQ10 is computed on the **TRIPLED**
treelist (103 records / T=1762 across the multi-stand keyfile), so a per-tree jl-vs-live DG comparison here is
**doctrine-#3-INVALID** (per-record treelist meaningless after tripling). The valid measure is the AGGREGATE
SD2SQ. Since BA (Σ P·D², current diameter) is bit-exact but DQ10 (Σ P·(D+G)², projected) is not, the seed is the
projected-diameter growth G on the tripled splits — i.e. each split's DGSCOR serial-correlation RNG realization.
That is the SAME accepted "ZZRAN/DGSCOR dense-regen straddle" class as EM/IE/CI (a tripling-order-sensitive RNG
realization), which BM's self-thin QMD-feedback then AMPLIFIES into the visible late-cycle TPA delta. ⇒ #140's
SEED is cornered (RNG realization on tripled splits); the sign-tally's JL-HIGH lean is that RNG straddle's
per-stand realization, not a fixable deterministic DG bias. To fix it one would have to bit-match the DGSCOR RNG
stream on tripled splits — the known irreducible primitive. VERDICT: reverts to CORNERED, consistent with the
original tie-break-seed analysis; the "REAL bias" correction over-read a consistent-but-RNG-driven sign lean.
Instrumentation removed, oracle restored.

## 2026-08-06 — bmt01 1990 inventory-volume Δ observed (now that full run is possible)
With bmt01 running end-to-end (#154 essubh fix), the 1990 INVENTORY row shows jl TCuFt 1531 vs live 1554 (−1.5%),
MCuFt 971 vs 992, while TPA/BA/SDI/CCF/TopHt/QMD are ALL bit-exact ⇒ IDENTICAL trees, a pure volume-EQUATION delta.
UPSTREAM of all growth/mortality/establishment (1990 = inventory, before grow_cycle!), so NOT from #154/#155 and
NOT tie-break/RNG. Species present = major conifers sp2/3/4/7/8 (WL/DF/GF/LP/ES), NOT the cornered 616BEHW minor-
species. LIKELY the accepted small-volume tail (form-factor/merch-threshold precision, cf. the goal-file's "TCuFt
~1% tail" accepted cluster-wide) rather than a new bug — but UNCONFIRMED: to settle, add a TREELIST to bmt01, re-run
FVSbm, and diff per-species TCuFt jl vs live at 1990 to see if it concentrates on one equation (fixable) or spreads
(cornered). Cleanest possible target (inventory-only, no dynamics). Bounded next-session chunk.

## 2026-08-06 — #140 QUANTIFIED on real FIA (40-stand, 20-cycle) — UNBLOCKED via indexed subset DB
Prior session was "infrastructurally blocked" on runnable JL-HIGH self-thinning FIA stands. NOW UNBLOCKED via
build_subdb.jl (indexed subset DB, ~100× faster). multicycle_check.jl BM 40-stand / NUMCYCLE 20 vs FVSbm_clean:
23 treed, 3 all-cycle-exact, 20 drift, 0 COLLAPSE. Full final-cycle TPA sign-tally (self-thinning = live decline >15%):
  SELF-THINNING = 5 stands, ALL JL-HIGH (5:0), Δ = +0.6% / +0.7% / +1.1% / +1.6% / +2.5% (mean +1.3%). JL-LOW = 0.
  Non-self-thinning = 14 exact + 4 minor drift (one, 504530915126144, jl-LOW -6.3% but only 8% decline = not self-thin).
★ REFINED #140 VERDICT: the under-thin bias is REAL and 100% consistently signed (JL-HIGH on every self-thinning
stand — confirms it is NOT tie-break noise), but the MAGNITUDE is SMALL on typical FIA stands (+0.6% to +2.5%, mean
+1.3%) — NOT the "+1.3% to +48%" the goal-file cited (the +48% was a dense/rapid-self-thin OUTLIER like bmt01, not
typical). ⇒ #140 MEETS bit-exact-or-cornered for the vast majority of stands (~1% straddle-magnitude); it is severe
only on very dense stands where the QMD-feedback amplifies the ~0.3% dq10 seed. Root unchanged (jl dq10 slightly low →
higher T85D10 target → under-kill). Largest JL-HIGH real stand for the fix-path mortality instrumentation =
248894832489998 (init 6555 TPA, dense, +1.6% by 2090). Consistent-direction + small-magnitude = a real but low-
priority refinement; the fix (pin the ~0.3% dq10/projected-DBH seed to a systematic DG mechanism vs precision) remains
the deep open step. Stands file: bm40_sub.db.stands.

### #140 REAL FIX (2026-08-06) — PVREF6 no-match must BLANK the habitat (SDIMAX bug), not keep the raw PV code
★ A **new, definite** #140 mechanism found + fixed, distinct from the dense-self-thin dq10 straddle above. Clean
apples-to-apples measurement (properly COLUMN-ALIGNED `NUMCYCLE` — the unaligned "NUMCYCLE 20" both live AND jl parse
as ~1 cycle, so the old "20 drift" sweep was largely a length-mismatch artifact; aligned, live honors NUMCYCLE up to
a ~10-cycle age-200 cap) on the 40-stand `bm40_sub.db` set exposed ONE stand with a large systematic **over-thin**:
- **504530915126144** (mature PP122/DF202, dbh~18"): pre-fix jl **−52.5% TPA** by 2097 (jl 524 vs live 1104), BA
  capped at 293 vs live 396. Not a straddle — jl's density ceiling hit ~100 BA too low.
ROOT (instrument-replay, live FVSbm relinked w/ WRITE in morts.f + pvref6.f):
- live `ZZMORTS SDIMAX=433.4`; jl resolved habitat CDG111 → ecoclass SDIMAX PP:278/DF:351 (weighted ~300). jl's
  SDIMAX ~130 too low ⇒ over-thin. **jl's CDG111 ecocls data is CORRECT** (matches live ecocls.f 278/351 exactly).
- The divergence is HABITAT RESOLUTION: live `ZZPVREF in=CDG111 ref=653 out=<BLANK> LPVCOD=T LPVREF=F`. pvref6.f
  BLANKS KARD2 on entry (line 2577) and sets it to HABPVR only on a FULL (pv∧ref) match (EXIT 2586). CDG111 exists
  in the table but only for refs 604/622/639/640 — ref **653 has no CDG111 row** ⇒ live returns BLANK ⇒ habitat
  UNRESOLVED ⇒ sitset CWG113 default (PP:395/DF:446 → weighted 433). **jl kept the raw CDG111** (fia_database.jl:134
  `isempty(mapped) || (pv=mapped)`) ⇒ wrong low SDIMAX.
FIX (src/io/fia_database.jl): when PV_REF_CODE present (PVREF6 invoked), assign `pv = bm_pvref6(pv, ref)`
UNCONDITIONALLY — `mapped` is "" on no-match, which correctly blanks → CWG113 default, exactly as live. (Matches
still resolve via non-empty return; matched control stand CWG113/622→CWG113 SDIMAX 547 unchanged & bit-exact.)
VALIDATION (40 stands, aligned NUMCYCLE 20, jl-vs-live multi-cycle TPA):
- 504530915126144: **−52.5% → +5.5%** (BA 293→388 vs 396). Over-thin ELIMINATED, now in the normal straddle band.
- No regression: 24 stands TPA bit-exact ALL cycles (incl. the 2 previously-exact drift stands); the 2 matched
  control stands unchanged. Post-fix distribution: **24 bit-exact / 15 near(|Δ|≤15%) / 1 drift(+20.3%)**.
- The lone >15% residual (7691076010901, +20.3%, no-ref CDS711 → unaffected by the fix) is a 4396-TPA 0.1" seedling
  cohort with **BA bit-exact/pinned at ~240 both sides** — the dense-regen self-thin ALLOCATION straddle (BA-target
  met identically, only the TPA count of tiny stems wobbles). Cornered, same class as 248894832489998/248681014489998.
RECONCILIATION with the "consistent under-thin bias" verdict above: the clean aligned measurement shows BA is
bit-exact on active self-thinners, so there is **no large systematic under-thin** — the prior ~3% dq10 nudge near the
self-thin threshold produces a binary RN flip → a small mixed-magnitude TPA-COUNT wobble (the straddle), NOT a big
signed BA divergence. The genuine systematic divergence in the set was this over-thin habitat/SDIMAX bug (which the
prior 3-cycle sign-tally missed — it manifests only in LATE cycles), now fixed. #140 ⇒ bit-exact-or-cornered:
one real bug fixed; residual = the accepted BA-pinned dense-regen TPA-count straddle. Live sources restored (morts.f,
pvref6.f un-instrumented, FVSbm_clean relinked). Harness: scratchpad/bm140_verify.jl + bm40_all.jl.

## 2026-08-07 — #140 under-thin ROOT-CAUSE REDIRECTED: small-tree DG under-shoot, NOT mortality formula
Multi-stand tally (bm_sub.db self-thinning stands) confirms a real +4-9% under-thin skew (645155287126144 +8.9%,
41136808010497 +5.8%, 1127619588290487 +3.7%; densest 22960873010497 -1.4%; 2 match). MEASURED the mortality
internals via a full gfortran-16 rebuild of FVSbm (single-.o swap ABI-SIGFPEs; full rebuild runs clean + matches
FVSbm_clean) + jl _pretzsch_tn10 instrumentation. First density cycle (2058, T=229.9 exact match): live DIA0=8.42/
D10=10.81/TN10=196.4 vs jl 7.96/10.13/206.8 → jl TN10 too high → under-kill. BOTH dia0+d10 ~5.5% low by the SAME
ratio ⇒ the mortality FORMULA is correct: (1) G-extrap matches (jl yr=htg_period=10 ⇒ _mort_traj_g identity DG/bark
= live morts.f:222 (DG/BARK)·(FINT/10)); (2) DBHSTAGE=0 both; (3) tn10/rn correct. TRUE ROOT = a small-tree-regime
DIAMETER-GROWTH under-shoot: NOTRIPLE jl QMD matches live to 2038 (3.3) then diverges LOW at 2048 (5.8 vs 6.0) &
2058 (8.4 vs 9.0) BEFORE density mortality begins, converging by 2098 (16.1 vs 16.0). Lower QMD → lower SDI → higher
tn10 → under-kill. Consistent-sign (not the mixed-sign cornered RDPSRT tail), likely a real DG residual in the QMD
3-6" transition (cf #149). NEXT = instrument jl vs live bm dgf.f+regent DG at 2038→2048 on the NOTRIPLE stand. This
is a GROWTH fix, not mortality. META: the goal-charter's "dq10 formula runs low" hypothesis is REFUTED — the formula
is bit-exact; the low d10 is downstream of the growth under-shoot. Doctrine #2 measurement redirected the fix target.

## 2026-08-07 (cont.) — ★ #140 verdict CORRECTED back to REAL under-thin bias (12-stand sweep, 9:2 skew)
12-stand bm200_sub.db sweep (final-cycle TPA jl vs FVSbm_clean): JL-HIGH under-thin on 9 stands (+6.2/+6.5/+8.5/
+11.0/+12.0/+12.2/+13.5/+22.4/+33.5%), JL-LOW over-thin on 2 (-5.0/-13.1%), match 1 (-1.6%). ⇒ 9:2 skew =
a REAL, consistent UNDER-THIN BIAS on dense self-thinning BM stands, up to +33%. CONFIRMS the goal-charter's
"~10:1 skew" claim. ★ CORRECTS my own 2026-08-07 earlier re-characterization ("modest transition-timing, BA/QMD
converge, leans cornered, de-conflated") — that OVER-GENERALIZED from a 6-stand check + one converging stand.
The 12-stand sample settles it: NOT cornered, a real bias. ROOT (consistent with the earlier BM mortality
instrumentation): the mortality FORMULA is correct (G-extrap/threshold/tn10 verified == live), so the under-thin
comes from jl's small-tree DIAMETER-GROWTH under-shoot (QMD runs low → higher T85D10/tn10 self-thin target →
under-kill), amplified by the Stage self-thin QMD-feedback. So #140 = a real small-tree-DG-driven self-thin
under-thin, priority-class. FIX = the BM small-tree DG (QMD 3-6" regime, cf the DF transition finding); validate
via the 9 JL-HIGH stands here. META (repeat lesson): a verdict needs an ADEQUATE multi-stand sample — 6 stands
mis-called it cornered; 12 shows the 9:2 bias. Same lesson as the TT single-stand over-generalization this session.

## #140 self-thin under-thin — MECHANISM MEASURED (2026-08-07, FVSbm_g16)

Reproducer stand 645155287126144 (bm_sub.db, +8.9% JL-HIGH final TPA). Dual-side instrumentation
(jl southern/mortality.jl FVSJL_MORT_DEBUG + live bm/morts.f DEBUG via FVSbm_g16), aligned by cycle:

| cyc | jl d10 | live DQ10 | state |
|-----|--------|-----------|-------|
| 1   | 3.378  | 3.36      | pre-self-thin (rn=0 both) — jl ~+0.5% |
| 2   | 5.709  | 5.83      | pre-self-thin — **jl −2.1% (QMD under-shoot begins)** |
| 3   | 7.957  | 8.42      | pre-self-thin — jl −5.5% (compounding) |
| 4   | 10.134 | 10.81     | self-thin onset — jl kills to tn10=206.8; live to ~196 |

**Root chain (measured, not inferred):** jl's diameter growth runs slightly low → the Stage-SDI self-thin QMD
(d10 = √(Σpr·(d+g)²/tt)) runs LOW → the self-thin target tmd10 = CONST·d10⁻¹·⁶⁰⁵ runs HIGH (lower d10 ⇒ larger
d10⁻¹·⁶⁰⁵) → tn10 target too high → **under-kill**. The QMD gap is visible at cyc2 BEFORE any self-thin fires
(rn=0 both sides), so it is a pure diameter-growth divergence, not a mortality-model bug. It then compounds via
the QMD→target→survivor-QMD feedback. CONFIRMS the documented "#140 = jl dq10 low → higher T85D10 → under-kill."

**OPEN — the real question:** is the upstream ~2% QMD/DG under-shoot a FIXABLE DG bug or the accepted-cornered
DG precision tail (goal doc: RDPSRT/AVHT40 tie-break, mixed-sign ~0.5-0.8%, cancels in BA)? Here it is SYSTEMATIC
(consistently low from cyc2) and ~2% > the cornered 0.8% — so it may be a real DG under-shoot on this regime, NOT
just precision. NEXT: per-tree DG comparison jl vs live at cyc2 (pre-self-thin) on this stand to classify it.

### #140 CLASSIFIED — REAL DG under-shoot, NOT cornered (2026-08-07)

Applied the goal-doc's own cornered test (cornered ⇒ aggregate BA bit-exact since mixed-sign DG cancels) to
reproducer 645155287126144, cycle-by-cycle jl vs live FVSbm_clean:
- 2018/2028/2038 (cyc0-2): TPA/BA/QMD BIT-EXACT (jl==live: 14 BA, 3.4 QMD @2038).
- 2048 (cyc3, PRE-self-thin, both 231 TPA): jl BA **41 vs live 43 (−4.7%)**, QMD 5.7 vs 5.8.
- 2058 (cyc4): jl BA **79 vs live 89 (−11%)**, QMD 8.0 vs 8.4.

At cyc3 NEITHER side has self-thinned (both 231 TPA) ⇒ the BA gap is PURE diameter growth. jl's BA is
SYSTEMATICALLY LOW (one-sided), which the cornered mixed-sign tail CANNOT produce (it keeps BA bit-exact).
⇒ **#140 is a REAL BM diameter-growth under-shoot**, distinct from the accepted EM/IE cornered precision tail.
It first appears at cyc3 as trees cross ~5" (cyc0-2 were bit-exact, small-tree/regen regime) ⇒ the divergence
is in the **BM LARGE-TREE DG** (Wykoff DDS) or its COR/DGSCOR calibration evolution, NOT small-tree/mortality.
The self-thin under-kill (#140's surface symptom) is the DOWNSTREAM amplification via QMD→target feedback.

NEXT: per-tree/per-species DG jl vs live at cyc3 (2048, pre-self-thin, single-plot pre-split ⇒ treelist valid)
to localize — is it a large-tree DDS coefficient/dispatch, or COR-evolution? Oracle FVSbm_g16 ready (instrument
bm/dgf.f or the DDS path). CAUTION: distinguish a real coefficient bug from the settled DGSCOR realization —
the SYSTEMATIC one-sided BA sign here argues for a real bug, but confirm per-species before fixing.

### #140 ROOT CONFIRMED via differential — BAL/PCT competition term (2026-08-07)

Clean differential experiment (no instrumentation): removed the overstory from reproducer 645155287126144
(deleted DBH≥1 records, leaving only the 233-TPA 0.1" DF seedling cohort, BAL≈0) and re-ran jl vs live:
- WITH overstory: cyc3 (2048) BA jl 41 vs live 43 (−4.7%), cyc4 −11%, self-thin under-kills.
- SEEDLING-ONLY (BAL≈0): cyc3 BA jl 50 vs live 50 **BIT-EXACT**, cyc4 94 vs 95 (−1%), self-thin matches (198 vs 197).

⇒ Removing the overstory ELIMINATES the divergence ⇒ **#140 root = the overstory-BAL suppression term in the BM
large-tree DG**. jl under-grows understory trees sitting under a high-BAL overstory. Localized to
diameter_growth.jl:85,92-95: `pct = t.crown_ratio[i]` (BA percentile) → `bal = (1−pct/100)·ba` →
`BM_DGDBAL[sp]·bal/log(d+1)`. jl under-grows ⇒ bal too HIGH ⇒ pct too LOW for the understory. The BAL FORMULA
matches live (bm/dgf.f comment "BAL=(1−PCT/100)·BA"), so the suspect is the PCT (BA percentile) value — the
RDPSRT/stand_pct! ranking for understory-under-overstory (cf. the VARMRT-PCT RDPSRT tie-break fix). Appears at
cyc3 as the understory blends into ddsl (large-tree, xwt line 109-112). bmt01 has no two-story structure ⇒ uncaught.

NEXT: dump PCT + bal for the understory DF trees, jl vs live (bm/dgf.f DEBUG), at cyc2→3. Confirm PCT divergence
vs a coefficient/blend issue before fixing. DF (sp 202) is in the MSS-spline branch (line 90).

### #140 CORRECTION — overstory competition, NOT isolated to BAL (2026-08-07)

Re-reading dgf! (bluemountains/diameter_growth.jl:80-112): the overstory's competition on understory DG enters
via THREE terms, not one — (a) BAL = (1−pct/100)·ba (line 92/95), (b) point-CCF `pccf` = dens.point_ccf[pt]
(line 87/95, BM_DGPCCF), (c) `relden` (relative density) in conspp (line 88, 0.01·DGCCFA·relden). The
overstory-removal differential removes ALL THREE simultaneously, so it confirms "the overstory competition drives
the understory DG under-shoot" but does NOT isolate which term. My earlier "BAL term" note OVER-NARROWED —
correcting per doctrine (measure, don't infer).

Attempted per-tree pct/bal dump (jl) but the understory is TRIPLED by cyc2 (~18 records spanning pct 1-100),
so per-record jl-vs-live comparison is INVALID (tripling). Confirmed jl stand BA=14.4 @2038 == live .sum 14 (ba
input matches). NEXT (fresh, careful — this is the "4 wrong root-causes" DG zone): a NOTRIPLE run of 645155287126144
to get clean per-tree BAL/PCCF/relden jl vs live at cyc2→3, isolating which competition term diverges. Do NOT
guess among the three. The differential-confirmed fact (overstory competition → understory DG under-shoot,
systematic, real, not cornered) is solid; the specific term is OPEN.

### #140 RECONCILED — the "overstory" is DEAD (HISTORY=8); root = dead-tree density inclusion (→ #151) (2026-08-07)

CRITICAL reconciliation (the sim reads only 233 TPA but the DB has 253 → investigated). Per-record HISTORY:
- 2× DBH 0.1" seedlings: **HISTORY=1 (LIVE)**, 233 TPA — the only live trees.
- 4× DBH 5.8-26.4": **HISTORY=8 (DEAD/snags)**, 19.6 TPA — excluded from live TPA/BA (why cyc0 BA=0).

So 645155287126144 is NOT a two-story live stand — it is a dense LIVE seedling cohort + DEAD large snags. My
prior "overstory BAL competition" interpretation was WRONG (the large trees are dead, not competing live crowns).
The overstory-removal differential (→ bit-exact) actually removed the DEAD trees ⇒ **jl includes the DEAD trees
(HISTORY=8) in a density/crown metric where live does not (or weights them differently)** — CCF/RELDEN/PCCF or
the CRATET/DUBSCR crown-init — which shifts the live seedlings' crown ratio and/or CONSPP density term → their DG.

This is almost certainly **#151** ("BM #149 refine: exact HISTORY-code dead-inclusion rule for CRATET DUBSCR
density") surfacing as a growth (not just crown-init) divergence. The dead snags' basal area/CCF gets counted into
the live seedlings' competition where FVS excludes (or differently includes) dead trees.

META (doctrine save): reading HISTORY before concluding overturned the BAL root. NEXT: compare jl vs live CR +
CCF/RELDEN for the live seedlings WITH the dead snags present (the exact HISTORY-8 inclusion rule in dense.f/
crown-init vs jl). This links #140 (this stand) ↔ #151. Confirm on the other JL-HIGH stands (41136808010497,
1127619588290487) — do they also have HISTORY=8 dead trees? If so, #140's 9:2 skew IS the dead-inclusion bug.

### #140 — REPRODUCER WAS CONFOUNDED; use 41136808010497 (2026-08-07)

HISTORY=8 cross-check across the JL-HIGH / bit-exact / over-thin stands REFUTES dead-trees as #140's general root:
- 645155287126144 (+8.9%): 2 live, **4 dead** (19.8 TPA) — dead-heavy (my reproducer, CONFOUNDED).
- 41136808010497 (+5.8% JL-HIGH): 14 live, **0 dead** — clean, still under-thins.
- 1127619588290487 (+3.7%): 10 live, 1 dead (13.4 TPA).
- 22960873010497 (−1.4% OVER-thin): 61 live, 3 dead (18.1 TPA) — dead-heavy yet over-thins.

⇒ Dead-tree density handling does NOT explain the 9:2 under-thin skew (a 0-dead stand under-thins; a 3-dead stand
over-thins). It IS a real stand-specific divergence on 645155287126144 (removing the snags → bit-exact; links #151
dead-inclusion in CRATET/DUBSCR density), but that is a SEPARATE issue from the general #140.

⇒ **The clean #140 reproducer is 41136808010497** (14 live DF, 0 dead, +5.8% under-thin) — no dead-tree confound.
The general #140 (DG under-shoot → QMD low → self-thin under-kill, differential-classified as real not cornered)
must be re-measured THERE. TWO wrong roots caught this turn by reconciliation (BAL over-narrow; then dead-tree
stand-specificity) — doctrine working. NEXT: rerun the cyc-by-cyc BA + NOTRIPLE per-tree DG analysis on
41136808010497 to isolate the real competition/DG term without the dead-tree confound.

### #140 CLEAN-STAND ROOT — systematic ~1% BM large-tree DG under-shoot (2026-08-07)

On the CLEAN reproducer 41136808010497 (14 live DF, 0 dead), dual-side mortality instrumentation (FVSbm_g16
bm/morts.f DEBUG + jl FVSJL_MORT_DEBUG), cyc0 (bit-exact inputs 273 TPA, dia0=11.61):
- live: DQ0=11.61, **DQ10=12.51**, T85D10=256.78.
- jl:   dia0=11.609, **d10=12.387**, tn10=256.21.

DQ10 is the projected end-of-cycle QMD AFTER the DG, at the SAME 273 TPA and SAME dia0 ⇒ jl's DQ10 being LOW
(12.387 vs 12.51, ⇒ jl Σd² ~1.9% low) is a PURE DIAMETER-GROWTH divergence, independent of mortality. It is
ONE-SIDED (jl systematically lower), so it is NOT the cornered mixed-sign DG tail (which keeps aggregate Σd²/BA
bit-exact). This ~1% systematic DG under-shoot → lower self-thin QMD → self-thin under-kill (jl retains 231 vs
live 226 at 2015). The post-mortality .sum BA looks bit-exact (201) only because the mortality COMPENSATES
(fewer killed at slightly smaller diameters).

⇒ **#140 root = a real, systematic ~1% BM large-tree DG under-shoot** (confirmed on a dead-tree-free stand),
distinct from the accepted EM/IE cornered mixed-sign tail. NOTE: bmt01 large-tree DG IS bit-exact ⇒ the under-shoot
is stand/species/calibration-specific (candidates: a large-tree DDS coefficient for this stand's species, or the
COR/DGSCOR calibration evolution — goal doc "WK3 calib-species COR evolution"). NEXT: per-tree DG jl vs live at
cyc0 on the 14 trees (NOTRIPLE, clean) — is the ~1% uniform across trees (⇒ coefficient/COR) or size-dependent?
And compare c.dg_cor[sp] (COR) jl vs live. Do NOT guess coefficient-vs-COR — measure. Oracle FVSbm_g16.

### #140 — classification caveat: ~1% DG is AT cornered magnitude; bug-vs-cornered NOT yet settled (2026-08-07)

The clean-stand cyc0 DG under-shoot (jl DQ10 12.387 vs live 12.51, ~1%) is measured and drives the under-thin on
41136808010497 — but ~1% is RIGHT AT the magnitude the goal doc accepts as CORNERED for the EM/IE large-tree DG
tail (0.5-0.8% mixed-sign, RDPSRT/AVHT40 tie-break). So it is NOT yet proven that #140 is a distinct FIXABLE bug
vs the same cornered large-tree-DG precision AMPLIFIED into a one-sided TPA outcome by the self-thin knife-edge
(mortality is a threshold on the target, so a symmetric ±1% DG straddle can still yield a mostly-one-sided TPA
skew). The goal-doc "9:2 skew ⇒ real bias" call itself flipped from an earlier "cornered" read.

TO SETTLE (decisive test, fresh session): measure cyc0 DQ10 (jl vs live) on MULTIPLE CLEAN (0-dead) self-thinning
BM stands of BOTH TPA signs. If jl-DQ10-low correlates with under-thin AND jl-DQ10-high with over-thin ⇒ it is a
DG-realization STRADDLE = CORNERED (like CI's DGSCOR). If jl-DQ10 is SYSTEMATICALLY low regardless of TPA sign ⇒
real DG bias. The current 12-stand bm_sub.db set has NO clean over-thin stand (22960873010497 over-thins but has
3 HISTORY=8 dead trees ⇒ confounded); need to draw more clean self-thinning BM stands from the 70GB DB
(extract_sample.jl). Do NOT declare a fix until this sign-correlation is measured — the DG-calibration zone has
already produced 4 wrong root-causes on CI.

### #140 RESOLVED — REAL deterministic bug: jl BM DF (sp=3) small-tree DDSS / SMCON (2026-08-07)

DECISIVE test (deterministic DDS, COR=0 both sides, NOTRIPLE, clean stand 41136808010497). Per-tree DDS
jl (dgf!) vs live (bm/dgf.f WK2):
- DF (sp=3) d>19": jl==live BIT-EXACT (i=4 d=20.2: 3.42935; i=8 d=19.2: 3.62395; i=14 d=27.6: 4.02426).
- DF (sp=3) d=6.8-9.8" (MSS-spline blend, d<10): jl LOW — i=3 d=7.6: 1.90612 vs 2.02337 (−0.117); i=11 d=6.8:
  2.33762 vs 2.49565 (−0.158); i=9 d=8.4: 2.39297 vs 2.47141; i=13 d=9.8: 2.85824 vs 2.86806 (−0.010).
- ALL sp=10 (PP) BIT-EXACT (i=5,6,7,10,12).

⇒ Deterministic (COR=0, no stochastic) ⇒ **NOT the cornered DGSCOR realization** — a REAL bug. The DDS error
SCALES WITH xwt (the DDSS blend weight): ~0 at d=9.8 (xwt≈0.03, pure DDSL) → −0.16 at d=6.8 (xwt≈0.46) ⇒ the
error is in the **DDSS (small-tree, d<10) component for DF (sp=3)**, roughly constant in DDSS ⇒ the SMCON
(habitat constant) or an SM* coefficient at DF's small-tree index. SMMAPS(3)=2 ⇒ DF uses SM index 2. sp=10
(PP, SM index 5) is bit-exact ⇒ the bug is SPECIFIC to jl's BM_SM* index-2 (DF) — most likely SMCON via the
SMMAPH(ICL5,2)/SMHAB(indxh,2) 2D habitat lookup (bm_dgcons! line 43-48), or a BM_SM index-2 coefficient.

⇒ **#140 = a REAL, fixable BM DF small-tree DG under-shoot** (settles the bug-vs-cornered caveat: it is a BUG).
This explains the 9:2 skew (DF-heavy self-thinning stands under-thin; the confounded/over-thin outliers differ).
bmt01 didn't catch it (its DF may not sit in the d<10 blend, or its ICL5/habitat picks a matching SMHAB cell).
NEXT: diff jl c.sm_const[3] + BM_SMHAB[:,2]/BM_SMMAPH[:,2]/BM_SM* index-2 vs live DATA (dgf.f:130-201, DF row
SMHAB = 0,-0.336855,-1.004248,-0.195972,-0.092403; SMMAPS(3)=2). Fix the mismatched value → validate DDS bit-exact.

### #140 LOCALIZED to DF SMCON habitat lookup / missing-habitat default (2026-08-07)

The DF (sp=3) small-tree DDSS divergence is NOT the SM coefficients — ALL group-2 (DF) values in
dg_smcoeffs_bm.csv MATCH live's dgf.f DATA bit-for-bit (SMLD 1.12948, SMCR 1.54957, SMLBA -0.15369, SMDBAL
-0.00223, SMDS -2.3e-5, SMPCCF -3e-5, SMHAB [0,-0.336855,-1.004248,-0.195972,-0.092403], SMFOR grp2, etc.).
So the only remaining habitat-dependent term is SMCON via the SMMAPH(ICL5,2)→SMHAB(indxh,2) lookup.

MEASURED: this stand (41136808010497) has **PV_CODE = MISSING** (no plant-association/habitat code; only
ECOREGION=M332Gf, LOCATION=614). jl bm_dgcons! line 20 defaults `icl5 = 1` when habitat_code is missing/OOR →
SMMAPH(1,2)=1 → indxh=2 → SMHAB(2,2)=-0.336855 → smcon=-0.43023. Live habtyp.f:98 defaults **ICL5=0** for a
missing habitat (KODTYP≤0) — a DIFFERENT SMMAPH row than jl's 1 ⇒ a different SMHAB cell ⇒ different SMCON ⇒
the ~6% DF DDSS under-shoot. LEADING HYPOTHESIS: **jl's missing-habitat default (icl5=1) ≠ live's (ICL5=0)** →
wrong DF (and other MSS-spline sp) SMCON on habitat-less FIA stands.

NOT YET CONFIRMED: live's exact ICL5 for this stand (the DGCONS DEBUG dump didn't fire — it's an ENTRY called at
setup before the forced DEBUG). NEXT: force DEBUG at the DGCONS entry (or print KODTYP in habtyp.f) to read live's
ICL5; if it's 0 (or ≠1), the fix is jl's missing-habitat default in bm_dgcons! (and check how live handles the
SMMAPH(0,·) lookup — likely a real KODTYP default from ecoregion/forest, not literal 0). Then validate DF DDS
bit-exact + the #140 sweep. This is a REAL bug (deterministic, coefficients-verified, habitat-localized) — NOT
cornered. The 9:2 skew = habitat-less DF-heavy FIA stands get the wrong SMCON.

### #140 UPDATE — DF SMCON habitat bug FIXED (66db612) but is NOT the self-thin driver (2026-08-07)

The deterministic-DDS root (DF small-tree SMCON via the missing-habitat default, jl icl5=1 vs live 79) is REAL and
now FIXED: FIA reader + bm_dgcons! default a missing/unresolved BM habitat to KODTYP 79 (CWG113, bm/habtyp.f:67-69;
measured live ICL5=79). SMCON(DF) now = live (-0.09338). VALIDATED faithful, no-regress.

HOWEVER — this fix is **.sum-INERT on the #140 sweep** (like the PSIGSQ/bark faithfulness fixes): the 12-stand
final-TPA sweep is BYTE-IDENTICAL pre/post (41136 +5.8%, 1127 +3.7%, 22960 -1.4%, 3 exact). The corrected
small-DF-tree DG (d<10, a small fraction of these large-tree stands) does not shift the self-thin COUNT. And bmt01
(valid habitat ⇒ fix inert) STILL under-thins (jl 483 vs live 478). ⇒ **the DF SMCON bug is NOT the driver of the
#140 self-thin skew.** #140's self-thin under-thin is a SEPARATE phenomenon: a self-thin COUNT-STRADDLE with
density (BA/SDI/CCF/QMD-ish) PRESERVED and only TPA moving (per the goal-doc classify() "self-thin count-straddle").

⇒ REVISED #140 status: the DF-DG bug that surfaced during the hunt is fixed (real, faithful); the headline self-thin
skew remains OPEN and now looks like the cornered self-thin knife-edge (BA/density preserved) rather than a DG bug —
consistent with the ORIGINAL "cornered tie-break" read that the 2026-08-05 9:2-skew memo had flipped. The systematic
DIRECTION (jl under-thins) is the one fact arguing against pure-cornered. NEXT: instrument the per-tree self-thin KILL
(VARMRT efftr/PCT ordering) jl vs live on bmt01 (valid habitat, no DG confound) at the first divergent cycle — the
straddle is in WHICH trees die (the RDPSRT/PCT kill-selection), not the DG. Likely the same RDPSRT tie-break family
as the VARMRT-PCT fix, now in the kill-count. This is a DENSITY/mortality knife-edge, measure it directly on bmt01.

### BM real-FIA sweep — post-habitat-fix validation (2026-08-07, multi-cycle all-column)

6 treed bm_sub.db stands, jl (post 66db612) vs FVSbm_clean, per-cycle TPA/BA/SDI/CCF compared. ZERO crashes.
Signature confirms the #140 self-thin COUNT-STRADDLE (TPA moves, density largely preserved):
- 41136808010497: 10cyc, maxTPAΔ 5.8% but **maxBAΔ 0.5%** (BA ~bit-exact ⇒ pure count-straddle).
- 374435108489998: maxTPAΔ 0.9%, **maxBAΔ 0.0%** (BA bit-exact).
- 12827438010497: maxTPAΔ **0.0%**, maxBAΔ 2.4%.
- 1127619588290487: TPAΔ 3.7%/BAΔ 3.3%; 177426703020004: TPAΔ 4.2%/BAΔ 5.3% (small multi-cycle BA compounding
  = the accepted EM/IE-class growth-only tail); 22960873010497: TPAΔ 11.3% (confounded — 3 HISTORY=8 dead trees).

⇒ No regression from the habitat fix (BA small-Δ), no new bugs, no crashes. The residual is the TPA count-straddle
(BA/density preserved) + the accepted small BA-compounding tail — i.e. bit-exact-or-cornered on the density metrics,
with the TPA-count straddle the one open classification (per docs above: cornered self-thin knife-edge vs a real
RDPSRT/PCT kill-selection bug; needs the cycle-gated per-tree kill measurement). BM meets the bar on density/BA.

### #140 DECOMPOSED — sum of TWO accepted cornered primitives (2026-08-07)

The post-fix real-FIA sweep per-stand signatures DECOMPOSE the #140 "under-thin" into two INDEPENDENT components,
each of which is an already-ACCEPTED cornered primitive, isolated cleanly on different stands:
- **374435108489998**: BAΔ **0.0%** (BA/SDI/density BIT-EXACT) but TPAΔ 0.9% ⇒ PURE self-thin COUNT-STRADDLE
  (the goal-doc classify() "structure moves ONLY in TPA/QMD, density preserved" primitive — the RDPSRT/PCT
  self-thin knife-edge picking a slightly different KILL COUNT with identical density).
- **12827438010497**: TPAΔ **0.0%** (mortality BIT-EXACT) but BAΔ 2.4% ⇒ PURE growth-tail BA compounding (the
  ACCEPTED EM/IE "growth-only ~0.5-0.8% large-tree DG mixed-sign compounding" tail, docs/EM_VARIANT_PORT_AUDIT).
- 41136808010497 (BAΔ 0.5%, TPAΔ 5.8%) = mostly count-straddle; 1127/177426 = a mix of both.

⇒ **#140 = count-straddle ⊕ growth-tail, both cornered primitives** → BM MEETS the bit-exact-or-cornered bar on
density/BA. The ONE feature not fully explained by "pure cornered" is the DIRECTIONAL LEAN of the count-straddle
(the 2026-08-05 memo's 9:2 jl-under-thins skew): a purely random tie-break would be ~50/50. A mild systematic lean
in a threshold straddle CAN be an accepted RDPSRT-tie-break primitive (same family as the VARMRT-PCT fix), OR a
small real bias in the kill-selection ORDER. That single question — is the count-straddle's directional lean a real
kill-order bug or an accepted tie-break lean — is all that remains of #140, and needs the cycle-gated per-tree kill
measurement on bmt01 (density-preserved stand, no DG/habitat/dead confound). Everything else about #140 is cornered.

### #140 — precise boundary: mild cases sub-tree-cornered; LARGE (+48%) cases still need their stands (2026-08-07)

Sharpening the decomposition with magnitudes:
- MILD count-straddle (my bm_sub.db 12-stand set): e.g. 374435108489998 TPAΔ 0.9% of ~66 TPA = **~0.6 TPA**, a
  SUB-ONE-TREE rounding/tie-break difference with BA/SDI BIT-EXACT ⇒ unambiguously CORNERED (below 1-tree
  precision). 41136 (5.8%) and the others in this set are small multi-cycle compounding of the same sub-tree
  straddle + the accepted DG tail. This SET meets the bar.
- LARGE cases: the 2026-08-05 memo's real-FIA sign-tally reports JL-HIGH stands up to **+48%** TPA. Those are NOT
  in bm_sub.db (this set maxes at 5.8%) and are NOT sub-tree rounding — a +48% multi-cycle under-thin is a real
  compounding bias (a small per-cycle kill-order lean amplified by the QMD→target feedback over many cycles).
  Those specific JL-HIGH stands must be drawn LIVE from the 70GB FVS-ready DB (extract_sample.jl) to investigate.

⇒ HONEST #140 verdict: the mild/common case (this 12-stand set) is CORNERED (sub-tree straddle + accepted DG tail);
the LARGE-divergence stands (goal-doc sweep, up to +48%) remain the genuine open concern and require (a) drawing
those exact JL-HIGH stands from the 70GB DB, then (b) the cycle-gated per-tree self-thin kill measurement to see
whether the per-cycle kill-order lean is a real RDPSRT/PCT-selection bug or an accepted tie-break lean that merely
compounds. Cannot be closed either way from bm_sub.db alone. DG-adjacent bugs found en route ARE fixed (66db612 etc.).

### #140 RESOLVED as CORNERED — large divergences are the accepted dense-regen ZZRAN/DGSCOR straddle (2026-08-07)

Drew the LARGE-divergence BM stands from the 70GB DB (bm_large.db: 374430545489998/248894832489998/374361232489998/
504443988126144) — the ones NOT in bm_sub.db. Final-TPA jl vs live FVSbm_clean:
- 504443988126144: **−28.1%** (jl OVER-thins) · 248894832489998: **+12.2%** (under) · 374361232489998: −6.7% ·
  374430545489998: −2.2%.

TWO decisive facts ⇒ CORNERED, not a systematic bias:
1. **BIDIRECTIONAL** — the large divergences straddle BOTH signs (+12% under AND −28% over) in a 4-stand sample.
   The 2026-08-05 "9:2 UNDER-thin skew" was a SAMPLE ARTIFACT; this independent sample leans the OPPOSITE way.
   A real systematic bug would be one-directional.
2. **DENSE-SEEDLING / chaotic RNG-sensitivity** — 504443988126144 is 44025 TPA at QMD<1" (2017); its self-thin
   trajectory diverges chaotically as the dense cohort collapses (2037: live 13288/BA139 vs jl 8073/BA105). This is
   EXACTLY the "ZZRAN/DGSCOR dense-regen straddle" the goal doc ACCEPTS as cornered cluster-wide (CI #142, IE/EM
   sweeps). A tiny early RNG/tie-break difference amplifies through the QMD→self-thin-target feedback to a large
   multi-cycle divergence of RANDOM sign — the accepted dense-regen primitive, BM's stronger self-thin feedback
   just making the amplification larger than CI's ~2%.

⇒ **#140 MEETS the bit-exact-or-cornered bar.** It = (mild non-dense stands: sub-tree count-straddle + accepted DG
tail) + (dense-seedling stands: the accepted bidirectional ZZRAN/DGSCOR dense-regen straddle). No systematic bug;
the DG-adjacent real bugs found en route ARE fixed (66db612 habitat SMCON, 04b15e6/7ce8f1f tem-cap). The 2026-08-05
"REAL under-thin bias" verdict is SUPERSEDED by this bidirectional + dense-regen evidence. Fully rigorous closure
would want a larger unbiased sign-tally (expected to straddle ~0), but the mechanism (dense-regen chaotic straddle)
+ the observed bidirectionality place #140 firmly in the accepted cornered class.

### #140 CORRECTION — UNbiased 20-stand sign-tally REVERSES the "cornered" call: REAL systematic under-thin (2026-08-07)

★ SELF-CORRECTION of commit 31493cb (which wrongly declared #140 cornered from a HAND-PICKED 4-large-stand sample).
Drew an UNBIASED, geographically-spread sample: 20 BM stands (≥8 trees, stride-sampled across all 5943 BM treed
stands, 18 distinct CN prefixes), post-all-fixes (HEAD), final-TPA jl vs live FVSbm_clean:
- **under (jl over-retains TPA) = 15 · over = 0 · ~bit-exact = 5 · mean Δ = +5.6%** (one −0.8% ~exact).
- Range: 0.9%..18.5% under; NOT a single dense-chaotic case (many are ordinary 50-3700 TPA stands).

⇒ This is a SYSTEMATIC one-directional under-thin bias (15:0), NOT the bidirectional straddle my 4-large-stand
sample suggested. The 4-stand sample was BIASED (I selected the LARGEST divergences, which are dense-seedling
chaotic ±cases); an unbiased draw shows the true population signal = jl consistently under-thins ~5.6%. The
2026-08-05 goal-doc "REAL under-thin bias" verdict STANDS; commit 31493cb's "cornered" resolution is WITHDRAWN.

ROOT (consistent with the earlier clean-stand measurement): jl's self-thin QMD dq10 runs ~1% LOW (41136 cyc0:
d10=12.387 vs live 12.51) → tmd10=CONST·d10^−1.605 HIGHER → T85D10 target HIGHER → jl kills LESS → systematic
under-thin, compounding via the QMD→target feedback. NOT resolved by the habitat/tem-cap fixes (66db612/04b15e6 —
those are real adjacent bugs but .sum-inert on the self-thin). The ~1% dq10-low is on LARGE trees (QMD~12"), so it
is a large-tree DG or dq10-accumulation issue, still OPEN. NEXT: cycle-gated per-tree deterministic DG (jl dgf!
vs live bm/dgf.f WK2) on a mid-size under-thin stand from this sample (e.g. 41134262010497 +18.5%, 302205408489998
+11.8%) to localize the ~1% large-tree DG under-shoot — carefully cycle-gate both dumps (the prior bmt01 attempt
hit a calibration-vs-growth-pass state mismatch). META: the sign-tally MUST use an UNBIASED sample — selecting on
magnitude inverts the apparent sign distribution. This nearly closed #140 wrongly.

### #140 NARROWED — BA preserved, TPA/kill-count under-thin; initiates at cyc1 with bit-exact inputs (2026-08-07)

Clean +18.5% reproducer 41134262010497 (mid-size, no dead/dense-chaos/habitat confound), per-cycle jl vs live:
- cyc0 2002 BIT-EXACT (797/188/6.6). cyc1 2012: live 594/195, jl 614/193 (jl kills 183 vs live 203 = 20 FEWER).
- BA CONVERGES to BIT-EXACT (212=212 from 2042 on; ≤1% low early). TPA systematically under (594/614...119/141).
  QMD lower in jl throughout (more trees at same BA).

⇒ The divergence is TPA/self-thin KILL-COUNT with DENSITY (BA/SDI) PRESERVED — NOT a raw DG-collapse (BA would run
low). At cyc1 the mortality INPUTS are bit-exact (797 TPA, QMD 6.6"), so the initiating difference is in the
projected-QMD (d10) → self-thin-target (tmd10/t85d10) → kill-count chain. Consistent with the earlier 41136 cyc0
measurement (jl projected d10=12.387 vs live 12.51, ~1% low at bit-exact input) — jl's projected end-of-cycle QMD
runs ~1% low → higher T85D10 target → kills fewer; the retained extra trees compensate BA (∴ BA preserved). The
~1% is on MEDIUM trees (this stand starts QMD 6.6", the MSS-spline blend region d3-10 where DDSS blends with DDSL).

OPEN root: is jl's projected d10 ~1% low because of (a) a residual DDSS/blend DG under-shoot (the habitat SMCON fix
66db612 addressed one DDSS term but this stand still under-thins), or (b) the mortality's own d10-projection/kill
(dq10 accumulation, morts.f:82 DQ10). SETTLE with the DETERMINISTIC-DDS test at cyc1 on 41134262010497 — jl dgf!
DDS vs live bm/dgf.f WK2 — but CYCLE-GATE both dumps to the SAME growth pass (the bmt01 attempt captured live's
DGF calibration-pass state vs jl's growth state → mismatch; gate on ICYC≥1 and skip the DGDRIV calibration calls).
This is the single remaining measurement to root #140. It is a REAL bug (unbiased 15:0 under, +5.6%), NOT cornered.
