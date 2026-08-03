# CI (Central Idaho) Variant Port — Audit Log

FVS Central Idaho variant → FVSjl. Tier-1 discount of the western Wykoff-DDS cluster.
Validated bit-exact-or-cornered per chunk vs the live relinked oracle
`/workspace/.ciwork/FVSci_clean` (gfortran-16 + isoc23 shim; `relink_ci.sh`).
Canonical test stand `cit01` (tests/FVSci/cit01.key — S248112, forest 412, habitat 520).

## Measured infra (chunk 0)

| Property | Value | Source |
|---|---|---|
| VARACD | CI | ci/blkdat.f |
| MAXSP | 19 | ci/blkdat.f JSP |
| Growth engine | **western Wykoff DDS** (dgf.f reads DGLD/DGCR/DGDBAL/DGHAB; NO GEMDG) | ci/dgf.f |
| SDI | **ZEIDE** (LZEIDE=.TRUE.) — like UT/TT, *not* Stage (key diff from IE) | ci/grinit.f |
| Cycle / seed | IFINT=10 (10-yr) / 55329 | ci/grinit.f, ci/blkdat.f |
| DGSD | 1.7 | ci/grinit.f |
| LHTDRG | .TRUE. all except sp15 (MC) | ci/grinit.f |
| Habitat | ci/habtyp.f ITYPE → NI 30 habitat types → OCURHT(16,MAXSP) | ci/habtyp.f |
| Template | **IE** (shared N-Rockies conifer set; +UT-style Zeide SDI) | source diff |

### Species (JSP / FIAJSP)
```
 1 WP 119   2 WL 073   3 DF 202   4 GF 017   5 WH 263   6 RC 242   7 LP 108
 8 ES 093   9 AF 019  10 PP 122  11 WB 101  12 PY 231  13 AS 746  14 WJ 064
15 MC 475  16 LM 113  17 CW 747  18 OS 299  19 OH 998
```
First 10 identical to IE's Northern-Rockies conifers.

### cit01 live baseline (FVSci_clean, growth cols)
```
1990  60  536  77 160  94 63 5.1  TCuFt 1541 MCuFt 833  BdFt 3912
2000  70  439  99 194 116 71 6.4        2176      1356       6380
2010  80  369 119 222 136 79 7.7        2914      2162      10203
```

## Chunk verdicts

| # | Chunk | Verdict |
|---|-------|---------|
| 0 | Scaffold: `CentralIdaho` singleton + registration (variant.jl) + include/export; oracle relinked; baseline captured | ✓ **DONE** — loads, `variant_from_code("CI")`→CentralIdaho(), nspecies=19; un-ported hooks error loudly (doctrine #5) |
| 1 | Species + grinit defaults (19 sp, FIA map, Zeide SDI, DGSD 1.7, LHTDRG, seed) — data/centralidaho/*.csv | ✓ **DONE** — species_coefficients.csv (19, IE-template non-DG + REAL ci/bratio.f BARK1/BARK2, ci/siterange.f SITELO/HI, ci/blkdat.f SIGMAR) + species_translation.csv (442, ci/spctrn.f ASPT col7); species.jl (Zeide SDI, DGSD 1.7, LHTDRG-15, seed 55329). cit01 loads species; coeffs verified (bark1/SIGMAR match live). Next hook site_setup! errors loudly |
| 2 | Site/habitat: ci/sitset.f + ci/habtyp.f | ✓ **DONE + validated** — site_index.jl (ci_habtyp ICITYP(130) bracket→ICINDX + NIHMAP→ITYPE; ci_forkod! JFOR/KFOR; ci_sitset! site-range SITEAR + forest-dep SDImax). STDINFO field-2→habitat_code wired for CI. cit01: habitat_code 520→**ICINDX 66, IFOR 4, IGL 3, SDIDEF=R4SDI [529,423,570,562,682,762…] bit-exact** (Zeide). Next hook dgf! errors loudly. |

### Chunk 2 — exact port recipe (ci/habtyp.f + ci/sitset.f traced)
CI has TWO habitat indices: **ICINDX** (1..130, CI-type) and **ITYPE** (1..30, NI). habtyp: bracket
`KODTYP` against **ICITYP(130)** → I; ICINDX=I−1, **ITYPE=NIHMAP(I−1)** (I==1 ⇒ both 1). DG uses
`MAPHAB=ICHBCL(ICINDX,ISPC)+1` (dgf.f:546 — so chunk 3 needs ICHBCL(130,19), not ITYPE directly).
sitset SITEAR default (NO MAPSIT — unlike IE): `TEM=SITEAR(ISISP)` (else 50; ISISP←3 if 0);
`SITERANGE(ISISP)→SLOSSP/SHISSP`; per sp `SITERANGE(sp)→SLO/SHI`; if SITEAR(sp)≤0 ⇒
`SITEAR(sp)=SLO+(TEM−SLOSSP)/(SHISSP−SLOSSP)·(SHI−SLO)` — SITERANGE = the SITELO/SITEHI already in
species_coefficients.csv ✓. SDImax: `LZEIDE=.FALSE.` iff CALCSDI blank & IFOR<2 (forest 1=Stage,
else Zeide); if BAMAX>0 ⇒ SDIDEF=BAMAX/(0.5454154·PMSDIU/100); else BAMAX=**BAMAXA(ICINDX)**, then
IFOR<2 ⇒ same formula, **IFOR≥2 ⇒ SDIDEF(sp)=R4SDI(sp)** (Zeide). cit01 forest 412→IFOR 4≥2 ⇒
**R4SDI** = [529,423,570,562,682,762,679,620,602,446,621,576,562,272,501,409,452,409,452].
forkod: JFOR=[117,402,406,412,413,414] NUMFOR=6 KFOR=[1,2,2,3,2,2]. JTYPE brackets = IE's (shared).
**Still to extract:** ICITYP(130), NIHMAP(130), ICHBCL(130,19) [ch3]. Data in hand: R4SDI, BAMAXA(130),
SITELO/HI, JFOR/KFOR. Validate SITEAR/SDIDEF echo vs FVSci_clean (instrument sitset.f or read SITECODE).
| 3 | Large-tree DG: ci/dgf.f (western Wykoff DDS) + DGHAB coeffs | ◐ **equation traced (recipe below) + 1-D coeffs extracted**. `data/centralidaho/dg_coeffs_1d.csv` = 14 per-species DDS arrays ×19 (DGLD/DGCR/DGCRSQ/DGBAL/DGDBAL/DGBA/DGLBA/DGPCCF/DGEL/DGEL2/DGSLOP/DGSLSQ/DGCASP/DGSASP) via tools/centralidaho/extract_dg_coeffs.py (grab, verified len 19; GF/WH share). **Remaining:** 2-D arrays (DGFOR/DGDS/DGDSQ/DGCCFA/DGHAB=OCURHT/ICHBCL(130,19)/MAPLOC/MAPCCF/MAPDSQ/IBSERV/OBSERV — per-column DATA forms), then port dgf!/ci_dgcons! (clone IE) + DGFTRC WK2 validate. |

### Chunk 3 — CI DG equation (ci/dgf.f, traced) — standard western Wykoff DDS (KT/UT family)
DGCONS (once/stand, per sp): `MAPHAB=ICHBCL(ICINDX,sp)+1`; `DHAB=DGHAB(MAPHAB,sp)` [DGHAB=OCURHT,
16-17 hab groups]; `ISPFOR=MAPLOC(IFOR,sp)`; `TEMEL=ELEV` (capped 30 for some sp);
`DGCON = DHAB + DGFOR(ISPFOR,sp) + DGEL·TEMEL + DGEL2·TEMEL² + (DGSASP·sin(asp)+DGCASP·cos(asp)+DGSLOP)·SLOPE`
`+ site term by sp-group {.001766·XSITE | .006460·XSITE | 0.227307·ln(XSITE)}` (check DGSLSQ·SLOPE² too).
Per cycle: `CONSPP = DGCON(sp) + COR(sp) + 0.01·DGCCFA(sp)·RELDEN`.
Per tree (default sp): `BAL=(1−PCT/100)·BA`; `DDS = CONSPP + DGLD·ln(D) + DGBAL·BAL + CR·(DGCR+CR·DGCRSQ)`
`+ DGDSQ·D² + DGBA·BAL/ln(D+1) − 0.000981·BA` (clamp ≥ −9.21). Special branches: **aspen** `DDS=ASPDG+ln(COR2)+COR`
(ci/dgfasp-like); **WB/DF-limited** `DF=0.25897+1.03129·DPP−0.0002025464·BATEM+0.00177·SI; DIAGR=(DF−DPP)·BARK;`
`DDS=ln(DIAGR·(2·DPP·BARK+DIAGR))+CONSPP`. DGLBA/DGPCCF used in a further sp sub-branch (ICLS path).
**Coeffs to extract (×19 sp unless noted):** DGLD,DGCR,DGCRSQ,DGBAL,DGDBAL,DGDSQ(DGDS),DGBA,DGLBA,DGPCCF,
DGCCFA,DGEL,DGEL2,DGSLOP,DGSLSQ,DGCASP,DGSASP; DGFOR(nloc,sp); DGHAB/OCURHT(17,sp); ICHBCL(130,sp);
MAPLOC(6,sp); IBSERV/OBSERV. Clone ie/diameter_growth.jl + ie/dg_coefficients.jl; ICINDX already in
p.habitat_input (ch2). Validate WK2 via instrument-replay (relink_ci.sh + DGFTRC dump) vs FVSci_clean on cit01.
| 3b | DG port (loader + ci_dgcons!/dgf!) | ◐ **PORTED; DGCON bit-exact; WK2 ~98% (2 stand-stat residuals left)**. Validated via DGFTRC instrument-replay on cit01: **DGCON bit-exact** (sp2 1.5005572, sp3 1.0425863, sp7 0.9978017 — ci_dgcons! + all coeffs correct). **BUG FOUND+FIXED**: DEFAULT DDS missing the 4th continuation line `+ CR·(DGCR+CR·DGCRSQ) + DGBAL·BAL` (dgf.f:493) → WK2 residual 0.64→~0.05. **REMAINING (2 stand-stat inputs; DDS formula itself is bit-exact)**:
(a) **RELDEN** live 81.89 vs jl 0.38 — jl `stand_ccf` returns 0.38 (broken: CI crown/CCF coeffs not ported ⇒ **BLOCKED on chunk 5** ci/ccfcal.f; once CCF works RELDEN follows). Only affects species with DGCCFA≠0.
(b) **BA scale** — TPA is FINE (jl loads 589.7 raw = BM identical; .sum normalizes to 536, validated for BM — false alarm cleared). Open: jl dgf uses p.basal_area=85.13 (raw) but live dgf-internal BA=66.80; trace ci/dense.f BA (large-tree/normalized). Feeds DGLBA·ln(BA).
⇒ Both residuals couple to the density machinery ⇒ **port chunk 5 (crown/CCF, ci/ccfcal.f) next** to unblock RELDEN + provide dense stats, then re-validate WK2 bit-exact. DG formula itself is validated (DGCON bit-exact + crown-term fix). |
| 3c | DG WK2 validation strategy | ★ Per-tree WK2 harness (notre!→setup_growth!→dgf!) is **CONFOUNDED by backdating**: live dgf reads DENSE-computed BA/RELDEN on the growth-period-START (backdated) stand (BA 66.80/RELDEN 81.89, ~0.79× the current 85.13/102.87), which the manual harness can't reproduce. Per **doctrine #3** (per-record invalid, use stand .sum), validate DG via the **full-cycle .sum after height (chunk 4) lands**. DG formula itself is validated (DGCON bit-exact + crown-term fix). |
| 4 | Height: ci/htgf.f | ◐ **PORTED + RUNS** — height_growth.jl: NI-conifer exp form (CON=HTCON+H2COF·HTI²+HGLD·lnD+HGLH·lnHTI; HTG=exp(CON+HDGCOF·lnDG)+BIAS) for sp 1-10,18; Weibull curve (COFLM 11/12/16, COFAS 13/17/19); WJ(14)/MC(15) exp with HTCON=0. HGHC/HGLDD/HGH2 by ITYPE=NIHMAP[ICINDX] (30-elem, pre-expanded). Coeffs extracted from ci/htgf.f (COFLM/COFAS = IE identical). cit01 projects past height → next hook small_tree_growth! (chunk 6) errors loudly. Validate via .sum after regent/mort/vol land. |
| 5a | CCF (ci/ccfcal.f MODE=1) | ✓ **PORTED** — crown.jl ci_tree_ccf (D≥10 → RD1+D·RD2+D²·RD3, else RDA·D^RDB; CI_RD1/2/3/A/B ×19 from ci/ccfcal.f, match KT for shared conifers) wired into stand_ccf/point_ccf (standstats.jl). RELDEN 0.38→102.87 (raw CCF ≈ .sum 94×raw-TPA). **Refine**: live dgf RELDEN=81.89 = dense.f RELDM1 (BACKDATED relative density) + TPA-normalization; couples to the BA-scale item. Crown-WIDTH (MODE=2 B1..B6) + crown_ratio_update! = chunk 5b (TODO). |
| 5b | Crown ratio: ci/crown.f | TODO |
| 6 | Small-tree: ci/regent.f | ◐ **coeffs extracted (17×19) + NIVAR height eqn traced**. regent_coeffs_1d.csv = CNST/CR/CRSQ/DGMAX/DHCN/DHCR/DHHT/HDM1/HDM2/HTBA/RLHT/STBA/XMAX/XMIN/BAL/PBAL/PTBA. NIVAR small-tree HEIGHT (ci/regent.f:641): `HTGRL = CNST + HTBA·RELHT·PTBAA + PBAL·PTBALI + STBA·BA + RLHT·RELHT + CRSQ·RCR² + CR·RCR + BAL·TBAL + HDM1·RHDM1 + HDM2·RHDM2 + PTBA·PTBAA` (RELHT=HT/AVH clamp1.5; PTBALI=PTBAA·(1−PCT/100); TBAL=(1−PCT/100)·BA; RCR=crown class 1-9) → CALL SMHTGF. Aspen path: HITE=26.9825·SITAGE^1.1752 (FINDAG). Bounds XMAX=5/XMIN=2. Port small_tree_growth!(::CentralIdaho) (subcycle scaffold from ie/regent.jl + this CI NIVAR eqn + SMHTGF + calib) next → unblocks .sum. |
| 7 | Mortality: ci/morts.f (Hamilton) + Zeide SDImax self-thin | TODO |
| 8 | Volume: ci/sitset.f VOLEQ via shared NVEL driver | TODO |
| 9 | Full-cycle differential vs FVSci_clean on cit01 | TODO |

## Port strategy (discount)
Clone IE's 15 Julia module files → retarget dispatch to `CentralIdaho` → swap in CI coefficient
CSVs extracted from ci/*.f. Two deltas from IE to watch: (a) **Zeide SDI** (reuse UT's Zeide
SDImax/self-thin path, not IE's Stage); (b) 19-vs-23 species + CI-specific DG/height/regent
coefficient DATA. Validate each chunk bit-exact via instrument-replay (relink_ci.sh) + cit01 .sum.
