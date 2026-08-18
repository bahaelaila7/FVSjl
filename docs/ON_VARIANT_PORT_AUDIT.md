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
| 8 | **compute_volumes! (volont.f ZAK/HONER + Mowraski)** | `740898ed` | **cyc0 .sum ALL columns bit-exact ⇒ ON growth+mort+vol END-TO-END.** ON default vol = METHC/METHB=8 → varvol OCFVOL/OBFVOL → volont.f: GTV (ZAKVOL/HONER total cubic), GMV (merch cubic; =0 = faithful BC-style metric quirk), NMV (Mowraski age-cull board). Age dubbed via findag→htcalc Carmean NC-128. Per-tree GTV/GMV/NMV **29/29 bit-exact**; cyc0 .sum cuft 28068/mcuft 0/bdft 16514/accr 5/mort 2777 all match. Fixes: gfortran `Z0**4`(=(Z0²)²) vs Julia z0^4 1-ULP→344-ULP in the Zak taper (→ `_on_p4`, ZAK 13/29→29/29); MAPLS/LTBHEC SI-curve 2-row transcription drop. Metric two-stage rounding Ontario-gated (BC one-stage path verified unbroken). ontario 178/178. |
| 7 | **mortality! (morts.f) + small_tree! (REGENT)** | `b7372204` | **run_keyfile COMPLETES cyc0→cyc1** — cyc1 stand TPA 952/BA 40/TopHt 27/QMD 23.1 matches FVSon_g16. mortality!(::Ontario) fresh port (Penner max-SDI/MSB line via most-BA species + new ON Density.mort_ibasp, 4-group background, VARMRT Penner individual-tree logistic): SDIMAX/tt/BAMAX/INDX_BA/TN10/RN bit-exact; per-tree WK2 **2/8 exact + 6/8 at 1-ULP** (VARMRT geometric-progression knife-edge, accepted mortality class). small_tree_growth!(REGENT) = validated no-op (ont01 all-large, DBH≥XMAX=4.72"; sub-12cm branch not ported — no validation stand). **DG-uniformity RESOLVED = #206 OLDRN straddle** (inject oracle DG → jl WK2 ≤1 ULP, proving the port is correct). ontario 80/80. |
| 6 | **height_growth! (htont.f/htgf.f)** | `c089b0e4` | Penner diameter-height (HTM=1.3+(B0+BSI·SIM+BDBHQ·DBHQM+BBA·BAM)·(1−exp(−B1·DBHM^B2)); 27 eqns OSPMAP_P 72→27). HTGF driver rebuilds 10-yr density in RECORD order. Per-tree **HTNOW/HT10/HTG dump-replay 8/8 bit-exact** (BA10/QMD10 exact in record order — species-sorted is 1 ULP off). run_keyfile advances past height_growth!. ⚠ flag: live cyc-1 jl DGs looked uniform ~0.0249 for 7/8 species vs oracle's per-species spread — a DG-calibration-chunk item (given identical DG, htont is bit-exact). ontario 80/80. |
| 10 | **FORTYP + size/stock classification (stkval.f/fortyp.f VBASE)** | _(this chunk)_ | **cyc0 `.sum` classification columns now bit-exact (`125 1 1` = red pine / sawtimber / overstocked) ⇒ the WHOLE cyc0 `.sum` row is byte-identical to FVSon_g16.** Root cause: ON never loaded the Arner-2001 VBASE stocking tables — `stock_b0`/`fia_stock_eq`/`fia_group`/`forest_type_codes` were all empty. ON is in the `stkval.f` no-western-redefine `VARACD` CASE (`CS/LS/NE/SN/ON`), so its TAB2/TAB3 are byte-identical to the eastern variants; copied `stocking_coeffs.csv` + `fia_stocking_map.csv` + `forest_type_codes.csv` into `data/ontario/`. jl's existing `_stkval_stocking` + `compute_forest_type!` (already VARACD-generic) then produce `125/1/1`. Report-only (ON's Penner `dgf!` never reads `forest_type`). ontario 193/193; multicycle 339/11 (ON-scoped data, no other variant touched). |
| 9 | **CCF — open-grown crown width (cwcalc.f IWHO=1) + ccfcal.f** | `28af93ef` | **cyc0 `.sum` CCF + AT_CCF now bit-exact (`****`).** Per-tree open-grown crown width **8/8 bit-exact** vs an instrumented FVSon_cwdump (cwcalc.f IWHO=1 dump). ON's cwcalc remaps ISPC→a US 2-char code (`ON_JSP2`) then evaluates the eastern Bechtold/Ek families `_cw_eval` already has; the generic `stand_ccf` `else` fed ON's raw code2 ⇒ miss ⇒ cw=0.5 ⇒ CCF≈15. Beech (ISPC 28) HI nails the oracle's **TLAT=46.78/TLONG=92.11** (forkod.f US-Superior 915/916 default, location lost to the initre overflow) + **metric STDINFO elevation 300 m→984.25 ft** (`kw_stdinfo!` ON branch, mirrors BC). stand_ccf=20844 ⇒ overflows the sumout.f **I4** field; added `_fi` (Fortran `Iw` integer edit → `*×w` on overflow, byte-identical to `%wd` when it fits) so the `.sum` renders `****`. ontario 190/190; multicycle **339/11 byte-identical** (summary format change is inert for non-overflow). |
| 5 | **DESIGN plot-count + QMD + metric report + merch** | `6d6d2b4a` | **cyc0 stand row BIT-EXACT** (TPA 84799 / BA 3673 / SDI 74769 / TopHt 24 / QMD 23.5). Root-caused a genuine **production `initre.f` KEYRDR array-overflow**: `ARRAY(7)/LNOTBK(7)` but `keyrdr.f` decodes NF=12 fields and writes `LNOTBK(8:12)` past the 7-long array into the adjacent `ARRAY`, zeroing `ARRAY(1:4)` of every ON card AFTER the numeric read (DESIGN IPTINV f4, INVYEAR year f1, STDINFO age/aspect f3/f4; fields 5+ survive) — confirmed byte-identical on production `FVSon_clean`. Emulated Ontario-gated (keep `present[1:4]`, zero `values[1:4]`) ⇒ IPTINV 11→0→clamp 1 ⇒ PI=1/GROSPC=1.0. + `compute_density!` sets `p.qmd`; + ON per-hectare metric reporting; + `init_merch_standards!(::Ontario)` (metric TOPD/BFTOPD/STMP/DBHMIN). |

## ✅ ON CORE COMPLETE + the ENTIRE cyc0 `.sum` row is byte-identical to FVSon_g16 (all numeric + classification columns).

## Remaining (dep-ordered)

1. ~~**FORTYP / size / stock classification**~~ — ✅ DONE (chunk 10): `125/1/1` bit-exact via the VBASE
   stocking-table wiring. cyc0 `.sum` row now fully byte-identical.
2. ~~**CCF**~~ — ✅ DONE (chunk 9): open-grown crown-**width** `cwcalc.f` (IWHO=1) + `ccfcal.f` + a
   `stand_ccf` Ontario branch, per-tree 8/8 bit-exact, `.sum` CCF/AT_CCF now `****`. Remaining CCF-adjacent
   metric artifact (NOT this chunk): the cyc1 IOSUM(3) TPA-column (oracle prints 24 vs BA-implied ~950; a
   metric sumout.f quirk). Also: `point_density!`'s PCCF for ON still uses the generic crown-width `else`
   (un-exercised — ON regen/small-tree not wired; add an `on_tree_ccf` branch there when regen lands).
3. **broken-top CFTOPK/BFTOPK + method-6 R9CLARK/method-5 Gevorkiantz** volume paths — no ont01 exercise;
   scoped, port when a keyword/stand exercises them.
4. **ON DG calibration** — `dgdriv.f` LSTART SIGMAR/OBSERV/VARDG → COR + the VARDG-driven tripled-record
   OLDRN spread (the multi-cycle #206 straddle; the multi-cycle `.sum` is already cornered-by-inheritance
   without it, as every western variant is).
5. The sub-12cm REGENT branch (htcalc NC128 + ONHTDBH + htdbh Wykoff dubbing) — needs a small-tree/regen
   ON validation stand (ont01 is all-large).
6. ON SPCTRN crosswalk targets (native-coded input already works via `resolve_species`).

**ON's core FVS projection (growth, mortality, volume) is bit-exact-or-cornered end-to-end vs FVSon_g16
on ont01 — 8 bit-exact chunks. The remaining items are stand classification (FORTYP/CCF) and
un-exercised alternate paths, not the core model.**

## Reusable oracle-instrument recipe
Instrument a scratch copy of `canada/on/{dgf,dgdriv,sitset,crown,htont}.f` with an EQUIVALENCE /
`TRANSFER→Z8.8` hex dumper, single-`.o` swap into `/workspace/.onwork/g16obj`, relink →
`FVSon_<x>dump`; A/B the instrumented `.sum` byte-identical to `ont01_clean.sum` first. Scripts in
`scratchpad/on/` (`build_g16_on.sh`, `dgf_dump.f`/`dgdriv_wki.f`/`link_wki.sh`, `validate_*.jl`).

**Status: ON is a large multi-session variant — 5 chunks done + ~5 remaining, each bit-exact.**
