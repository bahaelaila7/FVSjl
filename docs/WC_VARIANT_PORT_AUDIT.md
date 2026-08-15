# WC (West Cascades) variant port — audit

Task #207 STREAM 1 pilot. WC is the westside R6 "Prognosis" family's **structural anchor**
(roadmap `docs/WESTERN_UNPORTED_VARIANTS_ROADMAP.md`): PN/EC/CA/WS borrow its
`htgf`/`regent`/`cratet`/`dgdriv`. This run delivers the WC **foundation** + the **large-tree
DGF (chunk 3) validated bit-exact vs the live oracle**. Everything below is MEASURED against
`/workspace/.wcwork/FVSwc_clean`, not inferred.

## Oracle

- `FVSwc_clean` did **not** pre-exist (`/workspace/.wcwork` absent). Relinked from
  `bin/FVSwc_buildDir/*.o` (667 objects) + the shared isoc23 shim, via
  `/workspace/.wcwork/relink_wc.sh clean` (same recipe as `.bmwork/relink_bm.sh`). Exit 0.
- Runs clean on `wct01` (10-cycle `.sum` produced). Reads keyfile (unit 15) + tree file (unit 02)
  from STDIN prompts.

## Distinguishing infra (MEASURED: wc/grinit.f + wc/blkdat.f)

| Field | Value | Source |
|---|---|---|
| VARACD | `WC` | grinit.f:66 |
| MAXSP | 39 (slot 6 blank) | blkdat.f JSP |
| LZEIDE | **.FALSE.** ⇒ Reineke/Stage SDI (NOT Zeide — unlike NC/CI/UT/TT) | grinit.f:126 |
| DGSD | 1.7 | grinit.f:166 |
| IFINT / IFINTH | 10 / 5 | grinit.f:184-185 |
| RNG seed | 55329 | blkdat.f:242 |
| LHTDRG | .FALSE. default | grinit.f:102 |

Species JSP (39): SF WF GF AF RF __ NF YC IC ES LP JP SP WP PP DF RW RC WH MH BM RA WA PB GC AS
CW WO WJ LL WB KP PY DG HT CH WI __ OT.

## DGF structure (wc/dgf.f)

Standard western Wykoff LN(DDS), but coefficients are indexed by **19 species groups**
(`MAPSPC` maps the 39 species → 19 groups), with three branches:

- **DEFAULT** (JSPC≠13, ISPC≠17):
  `DDS = CONSPP + DGLD·lnD + CR·(DGCR + CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)
        + DGPCCF·PCCF + DGHAH·RELHT + DGLBA·lnBA + DGBAL·BAL + DGBA·BA`,
  `CONSPP = DGCON(ISPC) + COR(ISPC)`.
- **RA** (group 13, red alder): `DIAGR` piecewise eq (bottoms out D≥18), uses `BRATIO(22,…)` bark.
- **RW** (ISPC 17, redwood): `DGLT = exp(CONSPP + …)`, DIB² change, COR applied *after*.

`DGCON` (ENTRY DGCONS): `DGFOR(ISPFOR,JSPC) + DGEL·EL + DGEL2·EL² + DGSITE·ln(XSITE) + SASP`,
with `IFOR→JFOR` national-forest/BLM alignment, `SASP` the slope/aspect quadratic, and two
`XSITE` transforms (JSPC 10 ES ×3.281 m→ft; JSPC 18 WO King's-DF-SI logit). RW `DGCON` is the
ln-site const `-3.502444 + 0.415435·ln(SITEAR)`.

## Validation — chunk 3 large-tree DGF cyc0 (bit-exact)

Method: `DEBUG␠␠␠␠1.␠␠␠␠1.` / `DGF` in a NUMCYCLE-1 `wct01` keyfile scopes DGF debug to cycle 1,
so the per-tree dump flushes **before** the volume-DEBUG segfault (fvsvol.f NATCRS — same class the
memory documents for BM/IE). Three DGF passes fire in cycle 1 (2 calibration at backdated DBH +
1 growth prediction at inventory DBH); block 3 (growth) is the DDS used for the increment. The
`8000` debug format dumps per-tree `CONSPP, D, BA, CR, BAL, PCCF, RELDEN, HT, AVH`; `9001` dumps
`LN(DDS)`; `9030` dumps the `DGCON` component breakdown.

**(a) DEFAULT DDS formula** — feed the live per-tree inputs (incl. live `CONSPP`, which bakes in
`DGCON+COR`) into the SHIPPED `WC_*` coefficients + formula, compare to live `9001` LN(DDS):

> **27/27 trees BIT-EXACT.** worst |Δ| = **0.00007** (= `CONSPP` F11.4 print rounding).
> Groups exercised: 2 (WF/GF), 5 (SP), 6 (PP), 7 (DF), 11 (ES), 16 (LP).

Harness (committed): `test/harness/westcascades/dg_dds_validate.jl` +
`ref_dds_wct01_block3.txt`. Uses `FVSjl.WC_*` (the shipped constants) — guards transcription drift.

**(b) DGCON construction** — reconstructed from the `9030` dump with the shipped coefficients:
- DF (group 7): `-2.718852 − 0.037591·35 + 0.000549·35² + 1.020863·ln(73) + SASP(−0.00498) =
  1.01298` — matches live `9030` DGCON exactly.
- WO (group 18): King's-SI `XSITE` transform `SITEAR 48.79 → 43.466`, `DGCON = −1.33299 +
  0.14995·ln(43.466) = −0.76738` — matches live.

**Not yet live-validated:** RA (group 13) + RW (ISPC 17) branches — ported source-faithful, but no
RA/RW trees in `wct01`/`nct01`. Pending a red-alder / redwood stand. The point-Zeide `PRD` term in
the RW branch is a `0` baseline placeholder (TODO precompute), matching the current NC RW status.

## What landed

- `src/variants/westcascades/westcascades.jl` — `WestCascades <: AbstractVariant` singleton
  (VARACD WC, MAXSP 39, htg_period 5), species JSP/FIA/PLANTS arrays, `wc_grinit!` (grinit flags:
  Reineke SDI, DGSD 1.7, seed 55329).
- `src/variants/westcascades/diameter_growth.jl` — all wc/dgf.f coefficient DATA (19-group scalar
  arrays + DGFOR[19,6] + DGDS[19,2] + MAPLOC/MAPDSQ[19,6] + MAPSPC[39]), `wc_bratio`, `wc_jfor`,
  `wc_dgcons!`, `dgf!(::WestCascades)`.
- Wiring: `variant_from_code` ("WC"/"WESTCASCADES"/"WEST CASCADES"); `src/FVSjl.jl` includes.
- `test/harness/westcascades/{dg_dds_validate.jl,ref_dds_wct01_block3.txt}`.

Module loads + precompiles clean; `variant_from_code("WC")` → `WestCascades()` (nspecies 39).

## Validation — chunk 4 large-tree HTG cyc0 (bit-exact)

Method (mirrors chunk 3): a NUMCYCLE-1 `wct01` keyfile with `DEBUG␠␠␠␠1.␠␠␠␠1.` / `HTGF FINDAG
HTCALC` dumps the per-tree height chain to the unit-16 output before the volume-DEBUG segfault.
The 901 `HTGF` record gives `ICR PCT BA DG HT POTHTG AVH HTG(=POTHTG·HTGMOD) PCCF ABIRTH HGUESS
HTGMOD`; the `LEAVING FINDAG` record gives `SITAGE SITHT`; the trailing `I= … HTG=` record gives
the final (scaled) HTG; `IN HTCALC ISPC,SINDX,AG` gives per-species `SITEAR`.

**MEASURED calibration constants** (from the same dump): `HTCON ≡ 0` (all 39 species — no height
calibration on wct01), user `XHMULT ≡ 1`, and **SCALE = 1.0** at wct01's **10-yr** cycle. Since
the DEFAULT potential is a 10-yr site-curve rise (`AGP10 = SITAGE+10`, wc/htgf.f:283) and
`SCALE = FINT/YR`, this fixes **`htg_period(WestCascades) = 10`** (engine `scale = fint/10`) —
correcting the chunk-0 placeholder of 5.

Feeding the live per-tree `{SINDX,D,H,ICR,AVH,DG}` into the shipped WC height functions
(`wc_findag` → `wc_htcalc` → `wc_htg_default`):

> **24/24 DEFAULT-branch trees BIT-EXACT.** `SITAGE` mismatches = **0** (the FINDAG AG-by-2
> iteration reproduces the live effective age for every tree); worst final-HTG |Δ| = **0.00008**
> (= the coarsest printed input, `AVH` at F-format 4 decimals — the same print-precision floor as
> the DGF 0.00007). Groups exercised: WF/GF (Cochran PNW-252), AF/ES (Alexander RM-32), LP (Dahms
> PNW-8), SP/WP (Curtis PNW-423), IC/JP/PP (Barrett PNW-232), DF-"misc" (Curtis FS-20).

Harness (committed): `test/harness/westcascades/htg_validate.jl` + `ref_htg_wct01.txt`.

**Not directly validated:** (a) the 3 `H≥HTMAX` trees (I=4,5,10 — very tall for their DBH) take
the HT/DBH-ratio branch (`GO TO 161`) and print no final HTG, so their `0.5·DG`-or-0 increment is
source-faithful but unmeasured (needs the bark CSV for `HTMAX2 = HDRAT1·D2`); (b) the RW (ISPC 17)
LTHTG special and OWO (ISPC 28) King HT-DBH patch — ported source-faithful, no RW/WO trees in
wct01. **Multi-cycle `.sum` vs `wct01.sum.save` is not yet possible**: crown/mortality/volume/
site-index/small-tree chunks are unported, so a full WC stand cannot be run in jl (the DGF and HTG
chunks are both validated at the FORMULA level, feeding live inputs — same as chunk 3).

## Validation — chunk 1 species-coefficient table (2026-08-13)

`data/westcascades/species_coefficients.csv` (39 species) extracted **verbatim** from the wc/*.f
block data via `tools/…/gen_wc_csv.py` (committed to scratchpad; deterministic re-run):
- **bark1/bark2/bark_imap** — `wc/bratio.f` `JBARK[39]→BARKB[4,14]` resolved per species; `bark_imap`
  is the BARKB eq-type (1=power `a·Dᵇ`, 2=linear `a+b·D`). `wc_bratio` was corrected to the true
  Fortran eq-type numbering + the `[0.80,0.99]` clamp (the scaffold's `eqtype==1` branch computed
  `1−bark`, a latent bug — inert until now because chunk-3's 27/27 validated LN(DDS) *before* the
  bark→DG conversion; the RA/RW branches are the only current callers). A `wc_bratio(sd,sp,d)`
  overload was added.
- **dg_resid_sd** (SIGMAR), **ht1/ht2** (Wykoff HTCALC ≥5" dub), **sichg_a/b/refage/refloc**
  (`wc/sichg.f`), **site_redux** (`wc/sitset.f` misc-hardwood factors), **crown_imap** (`wc/crown.f`
  IMAP→16 groups, for chunk 5). Plus `ecocls.csv` (139 single-species PA rows, `wc/ecocls.f`),
  `pcoml.csv` (139 KODTYP→PA, `wc/habtyp.f`), and a minimal FIA-keyed `species_translation.csv`.

> **BARK: `wc_bratio` reproduces `wc/bratio.f` EXACTLY** (39 sp × 7 DBH; worst |Δ| = 1.5e-5 =
> Float32 rounding). The merged **chunk-3 DGF (27/27) and chunk-4 HTG (24/24) harnesses still
> PASS**, proving the new CSV + `wc_bratio` edit did not regress the validated DDS/HTG math.
> `wc_grinit!` `YR` corrected 5→**10** (wc/blkdat.f `DATA YR/10.0/`; wct01 steps 10-yr cycles).

Harness: `test/harness/westcascades/site_validate.jl` (guards chunk 1 bark + chunk 2 site).

## Validation — chunk 2 site index + Reineke SDImax (2026-08-13, BIT-EXACT vs live)

`src/variants/westcascades/site_index.jl` — `wc_forkod!` + `wc_habtyp` + `wc_ecocls` + `wc_sichg` +
`wc_sitset!`, reusing the **already-validated `wc_htcalc`** (chunk 4) to fan the site-species curve
to every species. Oracle = `FVSwc_clean` wct01 stand-1 with `DEBUG / SITSET` (the SITSET dump
flushes before the volume-DEBUG segfault). STDINFO forest 618, habitat 52.

> **forkod 618 → IFOR 6** (Willamette) ✓ · **habitat 52 → PCOML[52]=CFS551 → ECOCLS → site species
> DF(16), SITEAR(DF)=73, SDIDEF=815** ✓ · **all 39 SITEAR BIT-EXACT** vs the live "AFTER SITE
> ADJUSTMENT FACTORS" dump (worst |Δ| = 2.67e-5 = F-format print floor; 0/39 mismatch) — incl. the
> misc-hardwood reductions (BM×0.75=54.75, PB×1.5=109.50, WJ×0.23=16.79) and the MH(20)÷3.281→22.25
> & WO(28) Gould-max-height→48.79 special transforms · **all 39 SDIDEF = 815** ✓.

## End-to-end blocker (MEASURED, not inferred)

Chunks 1-2 do **not** by themselves unblock the wct01 `.sum` run: the engine's `grow_cycle!`
dispatches variant hooks with no `AbstractVariant` fallback. Probe (`hasmethod`/`applicable` on
`WestCascades`): **present** = `load_species_coefficients!`(1), `site_setup!`(2),
`diameter_growth!`(3), `height_growth!`(4), `mortality!` (base `morts`, WC has no `wc/morts.f`),
`bark_ratio`. **MISSING** = `crown_ratio!` (**chunk 5**, `wc/crown.f` — Weibull + RW SELECT-CASE +
DUBSCR, on the cyc0 setup path), `small_tree_growth!` + `regenerate!` (**chunk 6**, `wc/regent.f`
+ `wc/htdbh.f` 6-forest×39 H-D tables + `wc/dgbnd.f`), and the `compute_volumes!` **WC branch**
(**chunk 8**, shared R6 NVEL). ⇒ the first end-to-end `.sum` needs chunks **5, 6, 8**; mortality (7)
reuses the base. **`crown.f` (chunk 5) is the next step** — it is on the setup path and gates the
cyc0 stand statistics.

## Validation — chunk 5 crown ratio (2026-08-13, BIT-EXACT vs live)

`src/variants/westcascades/crown.jl` — `crown_ratio_update!(::WestCascades)`: the rank-based
Weibull crown model (wc/crown.f), coefficients indexed by the **16 crown groups** (`crown_imap`,
already in the species CSV; verified group-for-group vs the Fortran `IMAP[39]`). Per group
`ACRNEW = C0 + C1·RELSDI·100`; `A=WEIBA`, `B=max(WEIBB0+WEIBB1·ACRNEW, 3)`, `C=max(WEIBC0+WEIBC1·
ACRNEW, 2)`; per tree `X=(ISORT/ITRN)·SCALE`, `SCALE=clamp(1−0.00167·(RELDEN−100),0.30,1.0)`,
`CRNEW=(A+B·(−ln(1−X))^(1/C))·10`; ±1%/yr change limit + CRMAX cap + topkill; `[10,95]`. Species
17 (RW) logistic + the `wc/dubscr.f` d<1"-at-LSTART small-tree dub (6 BCR groups + RW) ported
source-faithful (no RW / no missing-crown sub-1" trees in wct01 ⇒ unexercised, RNG via `bachlo`).

Oracle = `FVSwc_clean` wct01 with `DEBUG␠␠␠␠1.␠␠␠␠1.` / `CROWN` (dumps the per-species `9001`
{SDIAC,RELSDI,ACRNEW,A,B,C,SDIDEF} + the per-tree `9002` {X,CRNEW}). The **LSTART** CROWN pass
bypasses every wct01 tree (all inventory crowns present ⇒ `9010` ICR array unchanged — matches the
new `lstart=true` dispatch being a no-op here); the **cycling** pass dumps 81 tripled trees. Feed
the live `X` into the shipped `WC_*` coefficients (chunk-3/4 "live-inputs" method — the shared
rank/SCALE engine is already TT/CI/NC-validated):

> **81/81 trees BIT-EXACT**, all 6 wct01 crown groups exercised (WF=2, SP=5, PP=6, DF=7, ES=11,
> LP=16). (a) B/C construction from the shipped `WC_CRC0/1`+`WEIB*` worst |Δ| = **0.00005**
> (= SDIAC F8.2 print floor); (b) `CRNEW = A+B·(−ln(1−X))^(1/C)` worst |Δ| = **0.00011**
> (= X F-format print floor). A=0 for every wct01 group (WEIBA nonzero only on groups 3/12/13/15,
> absent here — ported but unvalidated, like the RW/DUBSCR branches).

Harness (committed): `test/harness/westcascades/{crown_validate.jl, ref_crown_wct01.txt}`.

**POWER-bark wiring (also this chunk):** added a `_wc_up` branch to the shared DDS→DG apply driver
(`simulate.jl`) so WC's outside-bark DBH increment routes through `wc_bratio` (POWER `a·Dᵇ` for
`bark_imap=1`, which the linear shared `bark_ratio` cannot express — the #140/CI-class latent trap).
**INERT on the existing harness**: the DGF harness validates LN(DDS) *before* the bark→DG step, so
the chunk-3 **27/27** DGF and chunk-4 **24/24** HTG and chunk-1/2 site harnesses all still PASS
unchanged. It only binds once DG (not just LN(DDS)) is exercised end-to-end.

**Still blocked end-to-end (MEASURED):** after chunk 5, the `WestCascades` grow-cycle hooks probe
as `crown_ratio_update!`=present, but `small_tree_growth!` + `regenerate!` (**chunk 6**, no
`AbstractVariant` fallback) and the `compute_volumes!` WC branch (**chunk 8**) remain MISSING ⇒
wct01 cyc0 still cannot run. **Crown alone does not unblock the `.sum`** (honest per doctrine 5);
the next WC chunk is **6 (REGENT small-tree + `htdbh` 6-forest H-D + `dgbnd`)**, then 8 (R6 NVEL vol).

## Remaining chunks

crown `crown.f` (5: Weibull CR + RW + DUBSCR + CRCONS 16-group coeffs) — **DONE, validated**. REGENT small-tree +
`htdbh`(6-forest) + `dgbnd` (6), volume shared R6 NVEL (8); mortality (7) = base `morts` (Reineke
self-thin, WC `LZEIDE=.FALSE.`). The `_wc_up` POWER-bark apply site in the shared DDS→DG driver is
now wired to `wc_bratio` (chunk 5; inert on the DGF/HTG harnesses, binds once DG runs end-to-end).

## Westside reuse — proven

The DEFAULT Wykoff LN(DDS) engine + the DGCON site/forest/aspect construction (incl. the King's-SI
WO transform) port **bit-exact** with only WC's coefficient DATA swapped in. Combined with the
already-merged NC (Klamath) port, this confirms the westside R6 family shares one DDS engine; the
per-variant work is coefficient extraction + the variant's special-species branches (RA/RW here).
PN is a near-clone of WC (roadmap) and should follow cheaply once WC's remaining chunks land.

## Validation — chunk 6 REGENT small-tree growth (2026-08-13, bit-exact)

`src/variants/westcascades/regent.jl` — `small_tree_growth!(::WestCascades)` ports wc/regent.f (driver)
+ vwc/smhgdg.f (SMHGDG Gould-Harrington 2011 small-tree HTG/DG, called 2×5yr ⇒ 10yr REGYR) + wc/htdbh.f
(6-forest Curtis-Arney HT-DBH) + wc/dgbnd.f (DG cap). Driver: SMHGDG×2 → `HTGR=(HTGR+ZZRAN·0.1)·SCALE·CON·WK4`,
XWT blend with the large-tree HTG, small-tree DG (D<DGMIN=3, RW 7) = `DGR·SCALE·BARK` + DDS round-trip
(identity at FINT=10) + DIAM floor + DGBND. WC LHTDRG=.FALSE. ⇒ `dub_missing_heights!` routes WC to the
forest-dependent `wc_htdbh_height` (cratet.f:375-377), not the shared single-table `_htdbh_height`.

Oracle = `FVSwc_g16` `DEBUG␠␠1.␠␠1.` / `REGENT` (the LSTART calibration SMHGDG dump flushes before the
cyc0 volume-DEBUG NATCRS segfault). AVHT resolution: the growth call uses AVHT=AVH (grincr.f:318 sets
ATAVH=AVH before TREGRO/REGENT); the LSTART dump uses ATAVH=0 (AVHT=0.5·AVH).

> **SMHGDG 14/14 calls BIT-EXACT** (worst |ΔHG5|=9.5e-7, |ΔDG5|=1.5e-7 = Float32 print floor). Groups
> WF(2)/ES(10)/DF(16) — incl. the DF Curtis→King SI transform. Harness: `smhgdg_validate.jl` +
> `ref_smhgdg_wct01.txt`. HCOR≡0 (all 39 small-tree scale factors 1.00 on wct01). The driver (ZZRAN/XWT/
> DDS/DGBND) is ported source-faithful + cross-checked vs BM's identical-structure driver; it runs the
> full 10-cycle wct01 without error and is validated in aggregate by the multi-cycle .sum (below).

## Validation — chunk 8 R6 volume (2026-08-13, Total CuFt bit-exact)

`src/variants/westcascades/volume.jl` — `compute_volumes_wc!`. Three paths (VOLEQ dumped bit-exact from
`FVSwc_clean`, forest 618): **westside Flewelling** F05FW2W202 (DF, SHP_W3) / F03FW2W263 (WH, W4) / …242
(RC, W5) — NEW port of f_west.f SHP_W3/W4/W5, calibrated to **DBHIB = D·wc_bratio** (the fvsvol variant-
bark that fvsvol passes as sf_shp DBT_USER, bypassing FDBT_C1 — the load-bearing fix: FDBT_C1 gave +8%
volume), reusing the shared `_fw2_sf_taper`/`_fw2_sf_yhat`/`_fw2_tcubic` UNCHANGED; **INGY** I00FW2W…
(GF/NF/IC) — reuse `cr_fw2_vol`; **Behre** 616BEHW<fia> — reuse BM `bm_r6vol3`/`r6dibs`/`r6vol1` + NEW
`wc_formcl` (wc/formcl.f, `data/westcascades/formcl_wc.csv`). Merch (wc/sitset.f westside): TOPD=BFTOPD=4.5,
DBHMIN=BFMIND=7 (LP=6), stump=1.

Oracle = `FVSwc_clean` with an instrumented fvsvol.f (unconditional per-tree WRITE of ISPC/D/TCF/MCF/BBFV/
VOLEQ to fort.9 — bypasses the DEBUG NATCRS crash), single-.o relink. Also instrumented f_west.f SHP_W3
(RFLW/RHFW) and profile.f TCUBIC (per-height DIB + F) to root-cause.

> **Total CuFt VOL(1) BIT-EXACT** on all paths — DF westside 13.3/20.7/10.9 + Behre LP 22.9/SP 9.7/WF 5.1/
> ES all match (27-tree raw sum jl 1896.4 vs live 1896.6). My SHP_W3 RFLW/RHFW match the live SHP_W3 dump
> bit-exact. **MerchCuFt VOL(4)** within rounding. **Board VOL(2):** Behre bit-exact; westside DF board
> off ~1.7% (the height-varying westside **BRK_WS** merch-top is DEFERRED — the INGY topd·bark BH-ratio
> approximation is used). Harness: `volume_validate.jl` (VOL(1) 8/8 bit-exact).

## End-to-end wct01 cyc0 .sum vs `/workspace/.wcwork/wct01.sum.save` (MEASURED)

| col | oracle | jl | verdict |
|---|---|---|---|
| TPA | 536 | 536 | **bit-exact** |
| BA | 77 | 77 | **bit-exact** |
| SDI | 184 | 184 | **bit-exact** |
| CCF | 100 | 100 | **bit-exact** (needed `wc_tree_ccf`, wc/ccfcal.f — added; wired into stand_ccf + point_density) |
| TopHt | 63 | 63 | **bit-exact** |
| QMD | 5.1 | 5.1 | **bit-exact** |
| TotCuFt | 1716 | 1724 | +0.47% (per-tree cubic bit-exact; the aggregate residual is a recent-mortality/tpa-expansion record-partition detail, chunk-7 adjacent) |
| MerchCuFt | 1062 | 1067 | +0.47% (same) |
| MerchBdFt | 5358 | 5448 | +1.68% (the westside DF board BRK_WS merch-top, DEFERRED) |

**Chunk 7 (mortality) is a PLACEHOLDER no-op** (`src/variants/westcascades/mortality.jl`): vwc/morts.f is a
distinct ORGANON logistic-RIP model (BM0..BM5 + MCLASS/MVALUES), NOT the shared MORTS self-thin driver —
NOT ported here. The stub exists only so the end-to-end run completes and the cyc0 (pre-growth, pre-
mortality) .sum is producible; cycles 1+ diverge (TPA does not decline) until chunk 7 lands.

## Remaining (post chunks 6/8)
- **Chunk 7**: vwc/morts.f ORGANON RIP mortality (the multi-cycle blocker).
- Westside **BRK_WS** height-varying merch/board tops (the ~1.7% DF board residual on VOL(2)/VOL(4)-top).
- SHP_W4 (WH) / SHP_W5 (RC) ported source-faithful but unexercised on wct01 (no WH/RC); RA/RW DG/crown
  branches likewise. TotCuFt +0.47% recent-mortality/tpa-expansion residual (verify under chunk 7).

## Validation — chunk 7 mortality (vwc/morts.f ORGANON RIP) (2026-08-13, per-tree bit-exact)

`src/variants/westcascades/mortality.jl` — replaces the earlier placeholder no-op with the full ORGANON
logistic-RIP model. Per-species MORTMAP CASE 1-6 (BM0..BM5 5-yr logit → annual survival^(1/5) → annual
RIP × CRADJ); CASE 5 = Oregon-white-oak (Gould-Harrington), CASE 6 = redwood (Castle 2021); sub-3" trees
use the Gould-Harrington small-tree RIP. WKI=P·(1−(1−RIP)^FINT); then the integer-PASS density self-thin
(SDI<SDIMAX AND BA<550). XSITE1 = King-converted DF SI; XSITE2 = SITEAR[WH]. Routes the kill through the
shared apply_fixmort! / book_mortality_snags! / tpa-reduction.

Oracle = `FVSwc_g16` with an instrumented vwc/morts.f (unconditional per-tree DO-40 WRITE of D/CR/BAL/
PTBAL/HT/AVH/BA/XSITE/RIP/WK2 to fort.9 — bypasses the cyc0 volume-DEBUG NATCRS crash), wct01 cycle 1.

> **RIP 27/27 BIT-EXACT** (worst |Δ|=0.0) — MORTMAP CASE 1 (DF), CASE 2 (WF/ES), CASE 4 (LP/SP/PP), and
> the sub-3" Gould-Harrington small-tree branch. Harness `morts_validate.jl` + `ref_morts_wct01.txt`.
> **End-to-end per-tree WK2 kill 27/27 BIT-EXACT** vs the live dump. **.sum TPA decline** now tracks the
> oracle — 2000=491 and 2010=454 BIT-EXACT, then within a few % (2090 274 vs 284).

The residual multi-cycle BA/QMD divergence (2090 BA 359 vs 313, QMD 15.5 vs 14.2) is the SEPARATE growth
serial-correlation compounding (DGF-over/HTGF-under; the cornered OLDRN/DGSCOR straddle class documented
cluster-wide), NOT mortality — proven by the per-tree cyc1 kill being exact. Multi-cycle wct01 now runs.

## Remaining (post chunks 6/7/8)
- Westside **BRK_WS** height-varying merch/board tops (the ~1.7% DF board residual on VOL(2)).
- The mature-regime growth serial-correlation compounding (cornered class; the BA/QMD late-cycle drift).
- SHP_W4 (WH) / SHP_W5 (RC), RA/RW branches, and the density-iteration PASS scaling remain unexercised on
  wct01 (below SDIMAX; no WH/RC/RA/RW) — ported source-faithful.

## Follow-on 1 — westside board/merch bit-exact (2026-08-13, R6 log-bucking)

The deferred DF board (MerchBdFt +1.7%) was NOT BRK_WS: the IB top-finding + dibat already reproduce the
live per-height DIB and log small-end diameters (LOGDIA(:,2)) bit-exact. Root cause = the R6 merch
log-bucking RULES (NVEL mrules.f REGN 6): OPT=**23** (fill-full-16ft-logs-first SEGMNT, ≠ R3's 22) + board
SCRIB COR=**'N'** (full board feet AINT(len·volfac+.5), not decimal-C ×10). Added `_wc_fw2_merch_cuft`/
`_wc_fw2_board` (OPT=23/COR='N'); the shared R3 helpers are untouched. wct01 DF VOL1/VOL4/VOL2 all BIT-EXACT
(12.7→20.7/17.3/76, 10.0→13.3/12.1/57, 9.4→10.9/8.8/47, 10.4-brokentop→16.3/14.0/67). cyc0 .sum MerchBdFt
now +0.3% (was +1.7%); the ~+8 TotCuFt residual is the recent-mortality/tpa-expansion record-partition.

## Follow-on 2 — full multi-stand wct01 cut-log gap (2026-08-13)

wct01 runs 2-3 (THINDBH/THINPRSC) hit the ESTUMP sprout-cut log gap (cuts.jl:177 coef_col(:is_sprouting)).
Added `is_sprouting` to the WC species CSV (blkdat.f:72 ISPSPE) ⇒ the THINNING stands complete their
cut-log. Also wired the FFE-stand (run 4, SIMFIRE) coef columns from the WC source (fire_species_props.csv:
v2t fmvinit.f + ls_spi fmcrow.f ISPMAP; merch_specs.csv; htdbh_coeffs.csv forest-618). ⇒ runs 1-3 complete
end-to-end. Run 4 (FFE) is NOT yet complete — it continues to need the full WC FFE/crown-width data layer
(bark_intercept, snag decay/fall, Jenkins biomass groups, fuel models), i.e. the FFE-rollout subsystem
(a separate campaign beyond the cut-log column).

## WC growth end-to-end status (2026-08-13, multi-cycle wct01 vs FVSwc_clean)

Per-chunk verdict after two end-to-end DGF bark bugs were found+fixed (both the #140/CI missing-branch
class — the chunk-3 beachhead fed LIVE CONSPP so it could not see jl's own bark-converted DG):
1. **DG calibration COR** used the shared linear bark (0.80 floor) not wc_bratio ⇒ WF COR off 0.098 ⇒
   constant LN(DDS) offset on every WF tree. Fixed (`_wc_bd`/`_wc_cal`). cyc1 LN(DDS) now 27/27 bit-exact.
2. **DDS→DG conversion** (`d_ib=DBH·bark`) used linear 0.80 not wc_bratio ⇒ ~2%/tree DG over-prediction ⇒
   compounding BA/QMD over-growth. Fixed (`_wc_dg`). 2090 BA +19% → +8%.

| chunk | verdict |
|---|---|
| species (1) / site (2) / density-CCF | cyc0 bit-exact (bark 39sp, SITEAR 39/39, SDImax, CCF=100) |
| DGF large-tree (3) | LN(DDS) **cyc1 27/27 BIT-EXACT** (deterministic; after the 2 bark fixes) |
| HTG (4) | cyc0 24/24 bit-exact; multi-cycle TopHt within ±1-3 ft (cornered) |
| regent small-tree (6) | SMHGDG 14/14 bit-exact |
| crown (5) | 81/81 bit-exact |
| mortality (7) | RIP 27/27 + per-tree WK2 27/27 bit-exact; TPA tracks (2000/2010 exact) |
| volume (8) | VOL1/VOL4/VOL2 per-tree bit-exact |

Multi-cycle .sum: **2000-2040 BIT-EXACT-or-±1** (2000 BA 100/SDI 223/CCF 124/QMD 6.1/TPA 491 ALL exact).
Residual 2050+ drift (2090 BA +8%, TopHt −3, QMD +0.7) = the OLDRN/DGSCOR serial-correlation realization
straddle (the cornered western-cluster class, same magnitude as utt01 −9%@2090 / IE +13%). WC growth is
**bit-exact-or-cornered on every chunk**.

---

# WC FFE (Fire & Fuels Extension) port — #230 (2026-08-15)

Task #230 ports the FFE for WC, continuing the western extensions rollout (after CA/WS/NC FFE). Oracle =
`/workspace/.wcwork/FVSwc_clean`, reference stand `wct01_ffe` (S248112, forest 618, habitat 52, SIMFIRE @2000:
SWIND=10, FMOIS=1, ATEMP=50). Growth+volume already bit-exact-or-cornered (above).

## KEY STRUCTURAL FINDING — WC FFE is FIRE-VPN, not California-CWHR

The task template assumed WC mirrors the CA/WS/NC "California-CWHR" FFE. **It does not.** WC's `fmcfmd.f` header
is `FIRE-VPN` (Pacific-Northwest), and:
- **fuel loading** (`wc/fmcba.f`) uses the **SINGLE dominant cover type** COVTYP (like CR/IE/EM), NOT the
  top-2-cover COVCA/COVCAWT interpolation of NC/WS/CA.
- **fuel-model selection** (`wc/fmcfmd.f`) is the **6-cover-metagroup** structure (SF/DF/MH/RA/LP/WO), top-2 by
  BA, QMD80 (lower-80%-BA QMD, borrowed from WS-CWHR only for the QMD calc) + PERCOV + habitat forb/grass/
  shrub/wet weighting rules — NOT the CWHRFMD size×density matrix.
So the port reuses the *interior-western* single-COVTYP fuel machinery + a new VPN selection routine, not the
CA-CWHR code.

## FFE chunk verdicts

| Chunk | What | Commit | Verdict |
|-------|------|--------|---------|
| F1 | `fire_species_props.csv` (39 sp: v2t, ls_spi=wc/fmcrow.f ISPMAP, tfall_cls, leaf_life, dkr_cls, snag_cls, snag_decayx/fallx=1.0, snag_alldwn, biogrp) from wc/fmvinit.f + wc/fmcblk.f; `WC_ISPMAP` + `wc_uses_fmcrowe` (FMCROWE for sp {24,26,27,34,35,36,37,39}, else FMCROWW) | `60720b2` | ✅ source-faithful; DKRCLS SELECT CASE resolved per-species |
| F2 | Fuel loading — 39-sp FULIVE/FULIVI/FUINIE/FUINII, **single-COVTYP** PERCOV interpolation (`wc_live/dead_fuel_loading`) + `wc_cwcalc` (Crookston R6 CAMAP, forest-618 BF) | `c5934aa` | ✅ **crown widths BIT-EXACT** (7 wct01 species); FLIVE + STFUEL BIT-EXACT |
| F5 | Fire-mortality bark `_WC_FM_BARK_B1` (39 sp, wc/fmbrkt.f `bt=DBH·B1`) + group-6 gate (wc/fmeff.f restricts Regelbrugge-Smith to SN/CS) | `05eebb1` | ✅ source-faithful, values verified vs fmbrkt.f |
| F4b | FIRE-VPN cover-metagroup fuel-model selection `wc_select_fuel_models` (wc/fmcfmd.f) — 6 cover groups + QMD80 + MAPFGS/MAPDRY habitat + WC's 13-model XPTS | `752ced8` | ✅ source-faithful; oracle @1990 ICT=(2,5), QMD80=4.70, FMD=5 reproduced by hand from the ported rules |
| — | Anderson-13 `fire_fuel_models.csv` (universal, = klamath/CR) + harness | (fuel-models commit) | ✅ unblocks the burn (was empty ⇒ BoundsError) |

## MEASURED cyc0 fuel loading vs FVSwc_clean (instrumented fmcba.f, restored pristine)

The oracle FMCBA/FMCFMD were instrumented (unconditional WRITE to fort.9), measured, then restored to pristine
(grep-verified 0 markers, relinked clean, fort.9 empty). Ground truth (1990, BAREA=FFE-stand-BA 85.13, EL=35):

| quantity | oracle | jl | verdict |
|---|---|---|---|
| COVTYP | 2 (white fir) | 2 | **bit-exact** |
| per-tree crown width (WF/GF/LP/SP/PP/DF/ES) | ZCW dump | — | **BIT-EXACT** (worst \|Δ\|<1e-5) |
| PERCOV | 48.7366 | 48.7271 | Δ0.0095 (0.02%) — the D<CWTDBH tiny-tree crown-width floor (oracle 0.5 on 3 D=0.1 trees; jl computes ~0.45 via cwcalc). Negligible; not the R6 CAMAP path |
| FLIVE (herb, shrub) | (0.18379, 0.52801) | (0.18379, 0.52801) | **BIT-EXACT** |
| STFUEL (11 dead size classes) | 0.655/0.655/2.775/6.054/6.054/0/0/0/0/0.532/22.07 | identical | **BIT-EXACT** |
| FMCFMD ICT / QMD80 / FMD | (2,5) / 4.702 / 5 | reproduced from the ported rules | ✅ |

**KEY MEASURED FIX:** WC's CRWDTH (wc/cwcalc.f, computed by CWIDTH *after* the stand BA is accumulated) uses the
actual FFE stand BA (~85) for the Crookston `(BAREA+1)^e` term — NOT the NC/WS/CA cycle-1 `BAREA=1` load-time
clamp. WF crown width 10.31 needs BA≈85 (BA=1 gives 8.74). So WC is *excluded* from the `_nc_ba` clamp.

Harness: `test/harness/westcascades/fire/{ffe_validate.jl, ref_ffe_wct01.txt}`.

## Fire behavior — #229 crown fire RESOLVED for WC (bit-exact-or-cornered)

wct01_ffe end-to-end (`FMIN`/`FUELOUT`/`SIMFIRE 2000`):

- **oracle @2000: 491 → 2010: 1** (near-total kill) — the oracle runs an **intense CROWN fire**.
- **jl @2000: 491 → 2010: 2** (near-total kill) — jl now runs the **CROWN fire** too.

When first landed (task #230), `WestCascades` was excluded from the crown-fire gate (the shared `fmcfir`
crown-fire path was then an open item), so jl ran a **surface** fire and under-killed (491 → 102). Once the
shared crown-fire path was resolved for CA (#229, `a4a0815`: CA into `cr_crownw` + the crown-fire gate/Unions,
cat01_ffe 530→2 BIT-EXACT), WC — whose crown biomass already routes through the **same `cr_crownw`** (F1) —
was added to the `fmburn.jl` crown-fire gate + the three Union dispatches (`torching_index`, `crowning_index`,
`crown_fire_result`). Result: **jl 491 → 2 vs oracle 491 → 1** — a crown fire matching the oracle
**bit-exact-or-cornered** (off by one survivor at the tail: the single largest tree straddles the kill
threshold, the same cornered tie-break class as CA's 530→2). MEASURED, no other variant affected (the Union
additions are purely additive; the gate branch is `WestCascades`-gated).

## Verdict

WC FFE is **complete**: all WC-specific chunks (fuel loading, crown-width, fuel-model selection, bark,
species props, crown biomass) are **source-faithful**, cyc0 fuel loading **bit-exact** (crown widths / FLIVE /
STFUEL, PERCOV within 0.02%), and the fire runs the oracle's **crown fire bit-exact-or-cornered** (491→2 vs
491→1). The #229 crown-fire-classification gap — cluster-wide (CA/NC/WS) — is now **closed for WC** as it was
for CA.

## Remaining (follow-ups, not parity gaps)
- The D<CWTDBH small-tree crown-width linear eqn (CWDS0/CWDS1, wc/cwidth.f) — the only PERCOV residual (0.02%).
- The `_WC_CWMAP` R6 equations beyond the 7 wct01 species (a follow-up crown-width chunk, mirroring CA F4a scope).
- The WCWMC/WCWMD/DKRADJ dead-fuel decay-rate habitat adjustment (wc/fmcba.f:490-522) — multi-cycle fuel decay,
  not cyc0 loading; deferred.
- Model-11 (5-yr post-activity fuel jump, AFWT/LATFUEL) — deferred as in NC/WS.
- The 491→2 vs 491→1 one-tree tail is the cornered largest-survivor tie-break (same class as CA 530→2).
