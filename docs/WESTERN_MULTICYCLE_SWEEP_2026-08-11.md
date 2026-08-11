# Western cluster — GENUINE multi-cycle FIA differential (2026-08-11)

## Why this is new
The multi-cycle validation harness (`test/harness/fia/sweep/multicycle_check.jl`) and every
scratchpad sweep this campaign emitted `NUMCYCLE $nc.0` — a free-form decimal whose value lands at
**column 10**, outside FVS's fixed-column keyword field 2 (cols 11-20). FVS's reader mis-parses it and
**both live and jl silently fall back to 1 projection cycle**. So the "multi-cycle" differential had
been validating only cycle-0 + cycle-1 for every western variant — the doctrine "FIA validation MUST
be multi-cycle" was violated by an invisible harness bug.

Fix (commit 93d7a1a): emit the value column-aligned (`"NUMCYCLE" * lpad(string(Int(nc)),12)`). Verified
on BM 645068643126144: `NUMCYCLE 5.0` → 2 summary rows; aligned → 6 rows (5 projections), **both engines
identical** (jl reproduces FVS's column-parse quirk exactly → refutes the suspected jl INVYEAR/multi-cycle
bug #168). Harness also now self-classifies non-bit-exact stands: `cornered<2%` (RNG straddle) vs
`DRIFT≥2%`. Reproducer: `scratchpad/bm_numcyc.jl`, per-cycle dumper `scratchpad/percycle.jl <VAR> <cn>`.

## Genuine 5-projection results (NUMCYCLE aligned, 6-stand stratified sample per variant)

| Variant | Drift stand | maxrel | Mechanism (per-cycle live vs jl) | Verdict |
|---------|-------------|--------|----------------------------------|---------|
| **EM** | 11864108010690 | 30.9% | TPA jumps at 2046 (live 96→ jl 125, +30%), **BA matched** (27 vs 28), QMD depressed (6.4 vs 7.2) | **REAL — AUTOES over-establishment (#143 EM-side)** |
| **IE** | 373781950489998 | 13.3% | TPA matches live, **BA/QMD/SDI run high** (BA 136 vs 120, QMD 5.1 vs 4.8) | **REAL — large-tree DG OVER-growth** |
| **UT** | 286792451489998 | 8.5% | TPA matches exactly, **jl BA FROZEN ~119 while live grows to 130** (QMD 26.8 vs 28.0) | **REAL — large-tree DG UNDER-growth** |
| TT | 753177780290487 | 5.0% | TPA matches; BA +1.7% at end; the 5% is a 2-ft TopHt integer straddle at cyc-2 | cornered / mild-DG |
| CI | 51046399020004 | 3.5% | bit-exact except final-cycle TPA 166 vs 172 (mortality straddle) | cornered |
| BM | 645068643126144 / 658202425126144 | 2.4% / 1.2% | BA bit-exact-or-Δ1; 2.4% is a QMD 4.1-vs-4.0 print-rounding artifact; TPA ±0.25% | cornered (RNG straddle) |
| CR | 3307550010690 | 0.5% | — | cornered ✓ |

## ⚠ UPDATE 2026-08-11 (measured): UT is NOT a DG bug — hypothesis refuted
Deep per-tree instrumentation of UT 286792451489998 (pure Utah juniper) REFUTED the "large-tree DG
under-growth" reading, and nearly prevented a wrong fix (doctrine save):
- jl DGF LN(DDS) is bit-exact vs live FVSut_clean DEBUG (small trees) / +0.02 on 2 big trees (bark).
- BUT the APPLIED jl DG is a constant ~0.5455"/tree — UT REGENT (`small_tree_growth!`, regent.jl:129)
  overwrites `diam_growth` with a height-derived value for all `d < UT_RG_XMAX` (woodland XMAX>44, so
  the 43.6" juniper is regent-grown), discarding the DGF DDS. Looked like the bug.
- Live TREELIST proves live does the SAME: 43.6→44.15 (+0.55) == jl +0.5455; live grows ALL junipers
  ~+0.5" (11.2→11.7, 30.1→30.6, 43.6→44.1), NOT the DGF ~1.0". So jl woodland DG MATCHES live.
⇒ UT drift is NOT large-tree DG. With TPA matched (36→35 both), it's MORTALITY SELECTION (which juniper
dies each cycle → different BA removed) or a slow straddle. Reclassified — see task #170.
NOTE: IE(#171)/EM(#172) mechanisms below are inferred from .sum column signatures, NOT yet per-tree
verified — the UT lesson (regent-overwrite mimics DG divergence; live may match) means they need the
same per-tree measurement before any fix.

## The headline: multi-causal, mechanism differs by variant/regime
Confirms the TOP-PRIORITY note ("LARGE cross-variant large-tree-DG issue on MATURE/low-density stands")
AND refutes it being one root or one sign:
- **EM** — AUTOES **over-establishment** (TPA-driven, BA matched). jl adds ~30% too many small stems at a
  regen cycle (2046). Matches memory: "jl AUTOES OVER-establishes +63-86% TPA; jl+NOAUTOES bit-exact w/
  live; ie_autoes fires unconditionally simulate.jl:565." IE's #143 NSTORE fix (d089b78) did NOT close
  the EM-side. → **AUTOES investigation**.
- **IE** — large-tree DG **over-growth** on a DENSE young stand (2868 TPA, QMD 2.2). TPA/mortality match;
  diameter grows too fast (BA +13%). → **DGF/DGDRIV over-growth**.
- **UT** — large-tree DG **under-growth** on a SPARSE MATURE stand (36 TPA, QMD 24.6", near-deterministic,
  minimal RNG/regen confound — the CLEANEST reproducer). jl BA barely grows (119 vs live 130). → **DGF/DGDRIV
  under-growth**; opposite sign from IE ⇒ not a global DG scale error.

## Next (per TOP-PRIORITY plan: instrument dgf/dgdriv vs jl at first divergent cycle)
1. **UT 286792451489998** (cleanest, deterministic): FVSut_g16-instrument `dgf`/`dgdriv` per-tree LN(DDS)
   vs jl at 2013→2023 (first cycle jl BA already lags: live 121 vs jl 119). Isolate which DG term
   under-shoots on sparse mature large trees (candidates from memory: DGSCOR COR mis-calib on sparse-GST
   mature / low-competition DG term / DG size-cap).
2. **IE 373781950489998**: FVSie_g16 dgf per-tree at 2015→2025 — over-growth, opposite sign.
3. **EM 11864108010690**: FVSem_g16 AUTOES NSTORE dump at 2036→2046 — confirm over-establishment count
   vs live; port the EM-side of the #143 ie_autoes NSTORE fix if the mechanism matches.

## EM AUTOES deep-dive (2026-08-11) — root narrowed to per-point BAAA/ESTOCK attribution
Measured the EM reproducer (11864108010690) end-to-end (scratchpad/em_autoes.jl NOAUTOES toggle + jl
FVSJL_AUTOES_DEBUG + live estb/estab.f + base/dense.f source read):
1. **jl-NOAUTOES == live-NOAUTOES BIT-EXACT** ⇒ 100% of the divergence is AUTOES; growth+mortality faithful.
2. TPA/cyc: live 30,29,78,76,96,94 | jl 30,29,77,75,125,123. **1st ingrowth (2026) MATCHES** (jl+47.9 ≈
   live+49); the **2nd ingrowth (2046) over-establishes** — jl +51.9 vs live +20.
3. jl AUTOES_IN: icyc=2 baaa=7.35 total=47.9 | icyc=4 baaa=13.99 total=51.9 — jl ingrowth does NOT taper
   as the plot fills.
4. NSTORE suppression RULED OUT: live PLPROB counts only DBH<REGNBK small trees (estab.f:177) == jl point_small.
5. **ROOT = ESTOCK PROB1 via BAAA (per-inventory-point BA).** Source `base/dense.f:206-213`: live
   `BAAA(IP)` counts ONLY overstory (`D≥REGNBK`) trees, `BATREE·PI/GROSPC` per point. At 2046 the 2026
   cohort has grown into the overstory ⇒ live BAAA rises ⇒ ESTOCK suppresses to +20; jl `point_ba[1]`
   (all-tree, no PI/GROSPC) = 13.99 ⇒ under-suppresses ⇒ +52.
   ⚠ Contradiction to resolve by measurement: jl comment establishment.jl:1165-1169 says a prior
   overstory-only BAAA attempt was refuted by live measurement on an all-small-tree stand (live BAAA=20.66
   with no overstory) — so PI/GROSPC per-point expansion and/or backdated WK3 diameters matter.
   NEXT: FVSem_g16-dump BAAA(NNID)+PI/GROSPC+per-point overstory count at the icyc=4 tally; fix jl's
   per-point BA attribution as measured. Shared `ie_autoes_establish!` (EM+IE) ⇒ also closes IE #143 ingrowth.
   Task #172.

### ROOT PROVEN (instrumented FVSem_g16, 2026-08-11)
Instrumented `estab.f:430` (unguarded `BAAA(NNID)` dump) on stand 11864108010690 → **NPTIDS=4 inventory
points**. Live BAAA per point: 2026 `[9.61, 71.85, 0, 0]`; 2046 `[14.37, 92.69, …]`. jl's `s.density.point_ba`
holds the SAME per-point data (2026 `[7.35, 61.3, 0, 0]`; 2046 `[13.99, 82.59, …]`) — **but jl uses only
`point_ba[1]`** (establishment.jl:1170) for the whole-stand ESTOCK. Live computes ESTOCK **per inventory point**
(estab.f `DO 245`/`DO 2451` over NPTIDS: `SUM1=Σmax(BAAA,1)`, `SUM2=Σ SUM1/max(BAAA,1)`, then per-plot ITPP), so
its high-BAAA point 2 (61→82) suppresses ingrowth; jl ignores points 2-4 ⇒ over-establishes, worsening as point
2 fills (2026 ≈matches; 2046 jl+52 vs live+20). **FIX**: restructure `ie_autoes_run`/`ie_autoes_establish!` to
loop ESTOCK over all `nptids` with each point's `point_ba[pt]` and combine per estab.f, instead of the single
`point_ba[1]` scalar. Shared EM+IE (simulate.jl:586) ⇒ also closes the IE #143 ingrowth residual. Instrumentation
recipe durable in FVSem_buildDir/estab.f + build_g16.sh; oracle restored clean. Task #172.

### FIXED (f598868) — but the driver was plot_id, not PROB1
Continued measurement REFUTED the per-point-PROB1 fix (ie_estock CLAMPS BAA above ~13 ⇒ per-point prob1 all-equal,
inert) and re-examined the g16 tally: live's 2nd ingrowth tally **NEWTPP sum = 20** (vs tally-1's 69) — suppressed
by the carried per-point stock (a tally-2 plot showed `PLPROB=1.003, NSTORE=3, NEWTPP=2`). jl at tally-2 DID see
the stock (`tpacre_ingro=44.67`) but **all on one point** (`nz_points=1`). **REAL ROOT**: `ie_autoes_establish!`
set `t.plot_id[n]=Int32(1)` for every established tree ⇒ the whole cohort lands on point 1 ⇒ only point-1 plots get
NSTORE-suppressed next tally; the other points re-establish full. **FIX**: `ie_autoes_tally` accumulates a
per-inventory-point tally (plot n → point `div(n-1,idup)+1`), `ie_autoes_run` returns it, `ie_autoes_establish!`
creates a seedling per `(species,point)` with the true `plot_id` + that point's `point_ccf`. Summed tally
byte-unchanged ⇒ single-point stands & every first tally identical (canonical + single-point FIA inert by
construction); unit suite passes. **EM 11864108010690 drift 30.9% → 10.4%** (2046 TPA 125→86 vs live 96), 0
collapse. Residual ~10% = per-point cohort-distribution / BAAA magnitude (jl point_ba ~1.3× low; PROB1-magnitude
secondary since ie_estock clamps). Doctrine win: two measured refutations (per-point PROB1, NSTORE-absent) before
the plot_id root. Task #172.

Doctrine reminder: per-record treelist INVALID after tripling — use .sum aggregates per cycle (as here).
All reproducers durable in `test/harness/fia/{BM,ci,cr,em,ie,tt,ut}_sample.txt` + `scratchpad/percycle.jl`.
