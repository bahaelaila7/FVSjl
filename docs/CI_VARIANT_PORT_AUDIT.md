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
| 5b | Crown ratio: ci/crown.f | ◐ **PORTED + RUNS** — crown.jl crown_ratio_update!(::CentralIdaho): Weibull (ACRNEW=CRC0+CRC1·RELSDI·100; A=WEIBA; B=WEIBB0+WEIBB1·ACRNEW; C=WEIBC0 [WEIBC1=0]; crnew=(A+B·(−ln(1−x))^(1/C))·10, x=ISORT/N·SCALE) + CRMAX cap. CI coeffs WEIBA/WEIBB0/WEIBB1/WEIBC0/CRC0/CRC1 (from ci/crown.f); CRNMLT=1/DLOW=0/DHI=99 (no-op mods). RELSDI Zeide (sdiac/sp_sdi_def). Clone of BM/UT/TT Weibull form. cit01 projects past crown → next hook volume (:v2t, chunk 8). |
| 6 | Small-tree: ci/regent.f | ◐ **coeffs extracted (17×19) + NIVAR height eqn traced**. regent_coeffs_1d.csv = CNST/CR/CRSQ/DGMAX/DHCN/DHCR/DHHT/HDM1/HDM2/HTBA/RLHT/STBA/XMAX/XMIN/BAL/PBAL/PTBA. NIVAR small-tree HEIGHT (ci/regent.f:641): `HTGRL = CNST + HTBA·RELHT·PTBAA + PBAL·PTBALI + STBA·BA + RLHT·RELHT + CRSQ·RCR² + CR·RCR + BAL·TBAL + HDM1·RHDM1 + HDM2·RHDM2 + PTBA·PTBAA` (RELHT=HT/AVH clamp1.5; PTBALI=PTBAA·(1−PCT/100); TBAL=(1−PCT/100)·BA; RCR=crown class 1-9) → CALL SMHTGF. Aspen path: HITE=26.9825·SITAGE^1.1752 (FINDAG). Bounds XMAX=5/XMIN=2. Port small_tree_growth!(::CentralIdaho) next → unblocks .sum.

**Regent full scope (traced — intricate, multi-model, fresh-session-scale chunk):**
(1) subcycle scaffold (REGYR=5, NPER=FINT/5, KPER split — clone ie/regent.jl); (2) NIVAR small-tree HEIGHT
= HTGRL eqn (above) → **SMHTGF** (smhtgf.f, STOCHASTIC: BETA1=exp(1.17527−0.42124·lnTPCCF),
BETA2=exp(−2.56002−0.58642·lnTPCCF), HTG1=BETA1+BETA2·CR, HTGRTH=HTG1+ZRAND·STDDEV, ZRAND=BACHLO normal;
FINDAG aspen path HITE=26.9825·SITAGE^1.1752) → H2=H1+HTGRL·SCALE·XRHGRO·CON·WK4; (3) small-tree **DG**
DDS=DG·(2·BARK·D+DG)·SCALE2 (ci/regent.f:1138/1248, DG from small-tree model TBD-read); (4) **HCOR/DCOR**
calibration (CON=RHCON·exp(HCOR)); (5) special species (aspen FINDAG, PI/JU potential-height). cit01 is
large-tree-dominated (few D<5), so regent .sum impact is small — but needed to unblock the full-cycle run.
Given intricacy, this is the CI port's largest remaining single chunk.

**REGENT FULLY TRACED (all components; ready for focused multi-turn write-up):**
- Calibration: `CON = RHCON(sp)·exp(HCOR(sp))` (RHCON≈1, HCOR=0 pre-calib); H2=H1+HTGRL·SCALE·XRHGRO·CON·WK4.
- Small-tree final DBH from HTDBH: `DK = exp(DHCN + DHHT·ln(HK) + DHCR·ln(RCR))` (CIVAR: DBH held during subcycles, set at end); DG=DK−D_old; DDS=DG·(2·BARK·D+DG)·SCALE2.
- SMHTGF (smhtgf.f) OVERWRITES HTGRTH with a stochastic model — BETA path HTG1=BETA1+BETA2·CR + ZRAND·STDDEV, or FINDAG (LESTB establishment) HITE=26.9825·SITAGE^1.1752 — gated on ZRAND(≠−999, |·|≤2)/LESTB/D. **Write-time decision:** determine the exact call-condition selecting HTGRL-regression vs SMHTGF per tree (establishment/size).
- Subcycle: REGYR=5, NPER=FINT/5, KPER split, 0.985^k survival, per-subcycle BA/RELDEN backdate (clone ie/regent.jl scaffold). XMAX=5/XMIN=2 window.
The large-tree growth core (DG+height, DG bit-exact) is done; this regent is the last big growth chunk before mortality/volume/full-.sum.

**◐ regent PORTED + RUNS** (regent.jl): NIVAR path — subcycle (REGYR=5) HTGRL regression → H2=H1+HTGRL·SCALE·XRHGRO·exp(HCOR); HTDBH DK=exp(DHCN+DHHT·lnHK+DHCR·lnRCR); DG=DK−D_old. RHDM1/RHDM2 habitat-code flags (cit01 hab520→RHDM1=1). SMHTGF stochastic BETA/FINDAG establishment path deferred (cit01 small trees are conifers). cit01 projects past regent → next hook mortality _varmrt_efftr! (chunk 7). Validation via full-cycle .sum after mortality/volume land. |
| 7 | Mortality: ci/morts.f (Hamilton) | ◐ **PORTED + RUNS** — mortality.jl: full Hamilton RIP (2.76253+0.222310·√D−0.0460508·√BA+11.2007·G−0.554421/DD+PMSC+0.246301·RELDBH+6.07129·G/DD) + POTENT (POT/IPDG) + BA10 self-thin (IDENTICAL to IE form; ci/morts.f:269/329 confirmed BA-based, not Zeide-SDI). Coeffs PMSC/PMD/POT/IPDG(30,11)/IPDG2 from ci/morts.f; ITYPE=NIHMAP[ICINDX]; IFOR 1-6; bark=ci_bratio; bamax=CI_BAMAXA[ICINDX]. cit01 projects past mortality → next hook crown_ratio_update! (chunk 5b). |
| 8 | Volume: ci VEQNNC via shared NVEL | ✓ **PORTED + RUNS** — volume.jl compute_volumes_ci! (= UT: 400MATW r4vol + I15FW2W cr_fw2_vol + 400DVEW r4d2h; TOPD=6/DBHMIN=8/sp7=7). CI_VOL_EQ dumped from live cit01.out; wired in volume_equations.jl + compute_volumes! dispatch. |
| 9 | Full-cycle .sum vs FVSci_clean (cit01) | ✓ **1990 inventory FULLY BIT-EXACT** (TPA/BA/SDI/TopHt/QMD + BdFt 3912; MCuFt 821 vs 833 ~1.4%). Projected divergence root-caused to the DG **backdated-stand-BA**: live dgf uses BA=66.80 (inventory backdated ~1 period), jl uses current 85.13 → DGLBA·ln(BA) makes jl DG ~4.5% low → QMD lag (2000 6.1 vs 6.4) → less self-thin → TPA jl-high (490 vs 439). SHARED-engine (dgdriv CALL DGF(WK3) backdated) = accepted cornered class (present in EM/UT/BM), risky to change broadly. |

## Conclusion (CI)
CI stood up end-to-end this session: **9/9 chunks ported + running, cit01 1990 fully bit-exact** (growth + BdFt),
**4 real bugs caught by instrument-replay** (DG crown-term continuation; mortality self-thin=BA10-not-Zeide;
mortality BAMAX=SDIMAX·0.5454154·0.85; ci_bratio per-species branches). Remaining = the DG backdated-stand-BA
cornered tail (shared-engine, projected QMD/TPA) + SMHTGF-stochastic/special-species regent + FFE/harvest
(:v2t). CI proves the "cheaper ones" are a bounded one-session western-Wykoff discount (IE template + UT Zeide).

## Port strategy (discount)
Clone IE's 15 Julia module files → retarget dispatch to `CentralIdaho` → swap in CI coefficient
CSVs extracted from ci/*.f. Two deltas from IE to watch: (a) **Zeide SDI** (reuse UT's Zeide
SDImax/self-thin path, not IE's Stage); (b) 19-vs-23 species + CI-specific DG/height/regent
coefficient DATA. Validate each chunk bit-exact via instrument-replay (relink_ci.sh) + cit01 .sum.

## 2026-08-05 — CI cit01 over-kill LOCALIZED (goal re-anchored to full western cluster)
Session re-anchored: CR variant COMPLETE (docs/CR_VARIANT_PORT_COMPLETE); stop hook now points at the whole
western cluster + extensions (docs/WESTERN_ROLLOUT_GOAL.md), not CR alone. Resumed CI (#142).
MEASURED cit01 control vs live FVSci_clean: 1990 BIT-EXACT; then jl OVER-KILLS TPA — 2010 347/364 (−5%),
growing to 243/262 (−7%) by 2030, QMD slightly HIGH (7.8/7.7) ⇒ fewer-but-larger (over-thin of small trees).
Instrumented live ci/morts.f (BA/BAMAX/DQ10/TB/RZ per cycle) vs jl: **BA + BAMAX BIT-EXACT** (85.131, 265.146);
self-thin RZ negligible (0.0086); the divergence is **DQ10 jl 5.968 vs live 6.003** (−0.6%) ⇒ per-tree G lower in
the Hamilton RIP ⇒ slightly higher background mortality ⇒ the ~2% early over-kill (concentrated cyc0-1, matches
by 2010+). ROOT (confirms the memory-noted tail): the CI DGF reads CURRENT stand BA (p.basal_area=85.13) in the
DG **calibration**, where live's dgf reads the BACKDATED growth-period-start BA (66.80; dgdriv CALL DGF(WK3)) ⇒
jl's COR is fit against the wrong density ⇒ DG marginally low (sub-.sum-QMD-rounding but amplified by mortality's
G-sensitivity). SAME family as the CR 8th-bug (backdated calibration density). FIX DIRECTION: recompute stand
density from the backdated dbh for the CI DG-calibration DGF call — but this is in the SHARED
calibrate_diameter_growth! path, so it must be scoped/validated to NOT regress the bit-exact CR/KT/IE/EM/BM/TT/UT
(their DGF density-reads are already correct). Magnitude ~2% (DG-precision/self-thin tail = memory's "cornered
class"); real + systematic, not a straddle. NEXT: trace whether the shared calibration recomputes density post-
backdate for the other western variants (if yes, CI's dgf just needs to read that; if no, they compensate elsewhere).

## 2026-08-05 (CORRECTION, doctrine #2) — CI over-kill root-cause RE-MEASURED: NOT "reads current 85.13"
The prior entry's root-cause was WRONG — inferred from a stale memory note without measuring. Instrumented BOTH
jl (CI_DGBA env dump in dgf!) AND live (ci/dgf.f WRITE BA/RELDEN): the CI DGF calibration call reads the BACKDATED
density in BOTH — jl ba=66.881/relden=81.955 vs live ba=66.796/relden=81.890; the PREDICTION call is BIT-EXACT in
both (ba=85.131, relden=102.874). So the calibration backdating WORKS (shared calibrate_diameter_growth! line ~359
compute_density! on the backdated dbh) — jl does NOT read current 85.13. The ACTUAL residual is a tiny BACKDATING-
PRECISION difference: jl's backdated stand BA is 0.13% HIGH (66.881 vs 66.796; relden 81.955 vs 81.890). That fits
the COR against slightly-high density ⇒ DG marginally low ⇒ mortality DQ10 0.6% low ⇒ ~2% early over-kill. LEAD
(precise): the DENSE backdating (_backdate_dbh!) — WK3=√(d²·r), r from measured DG / bark. jl's backdated dbh is
slightly LARGE ⇒ ci_bratio (POWER bark DIB=BARK1·D^BARK2) gadj=g/bark slightly low, OR the r/rounding. Since the
OTHER western variants share _backdate_dbh! and are bit-exact, suspect a CI-specific input to it (ci_bratio value,
or the measured-DG field). NEXT: instrument the per-tree backdated dbh (WK3) jl vs live for a few cit01 trees
(pre-split, doctrine-#3-valid at LSTART) to see if it's bark or the r formula. Magnitude ~2% (fine precision tail,
memory's "cornered class"); real + systematic. META: measuring corrected a wrong inferred root-cause AGAIN — same
lesson as the CR session (don't trust stale notes; instrument).
