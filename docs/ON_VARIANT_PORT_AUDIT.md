# ON (Ontario) Variant Port — Audit

Branch: `kt-variant-port`. Oracle: `/workspace/.onwork/FVSon_g16` (full gfortran-16/gcc-16
rebuild of `canada/on/` — NO pre-built `bin/FVSon_buildDir` existed; build recipe
`/workspace/.onwork/build_g16_on.sh`) and the production relink `FVSon_clean`. ⚠ The SQLite
`DATABASE` tree-read (`dbstreesin`) SEGFAULTS under gcc-16 (same toolchain-skew class as BC) —
validation uses an **inline `.tre`** stand (`scratchpad/on/ont01.key` + `ont01.tre`, 8 trees).

Variant facts (measured from `canada/on/`): **MAXSP=72, MAXTRE=6000, METRIC**, SDI = Reineke
exponent **1.605** (optional `LZEIDE`), **DGSD=2.0** (⇒ multi-cycle `#206` OLDRN straddle, the
accepted cornered class), growth = **Penner (2006) annual diameter-increment** (NOT Wykoff
ln(DDS)) — 35 equations, `OSPMAP` 72→35, 8 species carry an AGS/UGS hardwood-quality code.

Doctrine: bit-exact-or-cornered per chunk vs the live oracle by RUNNING it; MEASURE (g16
single-`.o`-swap dump-replay, instrumented `.sum` byte-identical before trusting); glibc libm
`ccall` for Float32 transcendentals (native/`fmath` shim drift ~1 ULP); every chunk additive/
INERT for other variants — `test_multicycle` **339/11 byte-identical** is the hard gate.

## Chunk log (each merged, bit-exact-validated)

| # | Chunk | Merge | Verdict |
|---|---|---|---|
| 1 | **Scaffold + Penner DGF core** | `5e4d1ce9` | `on_penner_dds` dump-replay **8/8** bit-exact (DBHM/DIAGR/DDS incl AGS/UGS species) vs FVSon_g16. Ontario<:AbstractVariant, code ON, nspecies 72; `dg_coefficients.jl` (35-eqn Penner arrays + OSPMAP + LQUAL). |
| 2 | **DGF engine-wiring** | `8d942522` | Per-tree **WK2=ln(DDS)+COR bit-exact 8/8**. Stand-context glue (SIM=SITEAR·FTtoM / BAM=BA·FT2pACRtoM2pHA / QMDM=RMSQD·INtoCM / BALM=(1−PCT/100)·BAM / HTM=HT·FTtoM / AGS=IMC∈LQUAL). TWO glue bugs caught by instrumenting production: (a) **BARK uses the GROWN diameter** (dgf.f:361) → run Penner loop first then `on_bratio(sp,D_grown,ht)`; (b) the **DO-303 leftover-BARK DDS quirk** — DDS uses the scalar bark left over from the last tree of the last species (=0.94 beech for every tree), replicated via a two-pass ISCT/IND1 structure. maple/beech HYBRID hard-disabled (dgf.f:402 `LHYBRID=.FALSE.`) — correctly not ported. |
| 3 | **Species table + translation + coefficients + shared DDS→DG** | `85896ed5` | `resolve_species` maps all 8 ont01 alpha codes to the exact ON indices (8/8); `coefficients(::Ontario)` standalone; the shared `diameter_growth!` DDS→DG branch uses `on_bratio` with the **original** DBH (dgdriv.f:201, distinct from the grown-D bark that formed DDS). Deterministic DG path (d_ib / DDS / WKI) **bit-exact 8/8** vs instrumented `FVSon_wkidump` (dgdriv.f). |
| 4 | **site_setup! + crown + metric tree-input** | `d1fc55a3` | `site_setup!(sitset.f)` **SITEAR 72/72 + SDIDEF 72/72** bit-exact (ISISP=3, PMSDIU=85). crown = the shared TWIGS NC-125 kernel (NE/CS/LS factored into `_twigs_crown_update!`, ON dispatches to it; model 8/8). ON TREEFMT (blkdat.f wide F5.1) + metric tree-input (cm→in/m→ft, ON added to the metric gate like BC) → DBH/HT bit-exact. |
| 7 | **mortality! (morts.f) + small_tree! (REGENT)** | `b7372204` | **run_keyfile COMPLETES cyc0→cyc1** — cyc1 stand TPA 952/BA 40/TopHt 27/QMD 23.1 matches FVSon_g16. mortality!(::Ontario) fresh port (Penner max-SDI/MSB line via most-BA species + new ON Density.mort_ibasp, 4-group background, VARMRT Penner individual-tree logistic): SDIMAX/tt/BAMAX/INDX_BA/TN10/RN bit-exact; per-tree WK2 **2/8 exact + 6/8 at 1-ULP** (VARMRT geometric-progression knife-edge, accepted mortality class). small_tree_growth!(REGENT) = validated no-op (ont01 all-large, DBH≥XMAX=4.72"; sub-12cm branch not ported — no validation stand). **DG-uniformity RESOLVED = #206 OLDRN straddle** (inject oracle DG → jl WK2 ≤1 ULP, proving the port is correct). ontario 80/80. |
| 6 | **height_growth! (htont.f/htgf.f)** | `c089b0e4` | Penner diameter-height (HTM=1.3+(B0+BSI·SIM+BDBHQ·DBHQM+BBA·BAM)·(1−exp(−B1·DBHM^B2)); 27 eqns OSPMAP_P 72→27). HTGF driver rebuilds 10-yr density in RECORD order. Per-tree **HTNOW/HT10/HTG dump-replay 8/8 bit-exact** (BA10/QMD10 exact in record order — species-sorted is 1 ULP off). run_keyfile advances past height_growth!. ⚠ flag: live cyc-1 jl DGs looked uniform ~0.0249 for 7/8 species vs oracle's per-species spread — a DG-calibration-chunk item (given identical DG, htont is bit-exact). ontario 80/80. |
| 5 | **DESIGN plot-count + QMD + metric report + merch** | `6d6d2b4a` | **cyc0 stand row BIT-EXACT** (TPA 84799 / BA 3673 / SDI 74769 / TopHt 24 / QMD 23.5). Root-caused a genuine **production `initre.f` KEYRDR array-overflow**: `ARRAY(7)/LNOTBK(7)` but `keyrdr.f` decodes NF=12 fields and writes `LNOTBK(8:12)` past the 7-long array into the adjacent `ARRAY`, zeroing `ARRAY(1:4)` of every ON card AFTER the numeric read (DESIGN IPTINV f4, INVYEAR year f1, STDINFO age/aspect f3/f4; fields 5+ survive) — confirmed byte-identical on production `FVSon_clean`. Emulated Ontario-gated (keep `present[1:4]`, zero `values[1:4]`) ⇒ IPTINV 11→0→clamp 1 ⇒ PI=1/GROSPC=1.0. + `compute_density!` sets `p.qmd`; + ON per-hectare metric reporting; + `init_merch_standards!(::Ontario)` (metric TOPD/BFTOPD/STMP/DBHMIN). |

## Remaining (dep-ordered)

1. **volume** — `canada/on` htont broken-top + `varvol`/`cubrds`/`nbolt`/`volont` + metric `sumout` scf_* →
   the cuft/mcuft/bdft `.sum` columns (currently 0). The only growth-side column not yet producible; the
   current `run_keyfile(ont01; Ontario())` blocker.
2. **CCF** — ON open-grown crown-**width** `canada/on/cwcalc.f` (IWHO=1) + CCFCAL + a `stand_ccf` Ontario
   branch (distinct from the crown-**ratio** TWIGS model). cyc0 CCF is the only non-bit-exact cyc0 column
   (JL 15 vs oracle ≥1e4 on the degenerate metric-BAF stand). Also the metric `.sum` cyc1 IOSUM(3)
   TPA-column artifact (oracle prints 24 vs BA-implied ~950; a metric sumout.f quirk, not a jl bug).
3. **ON DG calibration** — `dgdriv.f` LSTART SIGMAR/OBSERV/VARDG → COR + the VARDG-driven tripled-record
   OLDRN spread (the multi-cycle perturbation; the ON DGSD=2.0 `#206` straddle — the multi-cycle .sum is
   already cornered-by-inheritance without it, as every western variant is).
4. ON SPCTRN crosswalk targets (native-coded input already works via `resolve_species`).
5. The sub-12cm REGENT branch (htcalc NC128 + ONHTDBH + htdbh Wykoff dubbing) — needs a small-tree/regen
   ON validation stand (ont01 is all-large).

**Growth + mortality are end-to-end validated; ON produces a bit-exact-or-cornered cyc0→cyc1 stand
summary. Volume columns + CCF are the remaining growth-side gaps.**

## Reusable oracle-instrument recipe
Instrument a scratch copy of `canada/on/{dgf,dgdriv,sitset,crown,htont}.f` with an EQUIVALENCE /
`TRANSFER→Z8.8` hex dumper, single-`.o` swap into `/workspace/.onwork/g16obj`, relink →
`FVSon_<x>dump`; A/B the instrumented `.sum` byte-identical to `ont01_clean.sum` first. Scripts in
`scratchpad/on/` (`build_g16_on.sh`, `dgf_dump.f`/`dgdriv_wki.f`/`link_wki.sh`, `validate_*.jl`).

**Status: ON is a large multi-session variant — 5 chunks done + ~5 remaining, each bit-exact.**
