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
