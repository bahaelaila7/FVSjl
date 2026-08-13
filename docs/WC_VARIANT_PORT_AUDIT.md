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

## Remaining chunks (TODO — out of this bounded run's scope)

Species-coefficient CSV (chunk 1: bark/crown/site/SDImax/volume — needed to RUN a WC stand
end-to-end), site index + Reineke SDImax (2), crown `crown.f` (5), REGENT small-tree + `htdbh` +
`dgbnd` DG-bound (6), mortality (base `morts`, Reineke self-thin) (7), volume (shared R6 NVEL) (8).
Then end-to-end `.sum` vs `wct01.sum.save`. **Chunk 1 (species CSV) is the highest-leverage next
step** — it unblocks the first end-to-end run and the multi-cycle `.sum` trajectory that both the
DGF and HTG formula-level validations currently cannot reach.

## Westside reuse — proven

The DEFAULT Wykoff LN(DDS) engine + the DGCON site/forest/aspect construction (incl. the King's-SI
WO transform) port **bit-exact** with only WC's coefficient DATA swapped in. Combined with the
already-merged NC (Klamath) port, this confirms the westside R6 family shares one DDS engine; the
per-variant work is coefficient extraction + the variant's special-species branches (RA/RW here).
PN is a near-clone of WC (roadmap) and should follow cheaply once WC's remaining chunks land.
