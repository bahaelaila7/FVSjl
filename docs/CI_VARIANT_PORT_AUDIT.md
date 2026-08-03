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
| 1 | Species + grinit defaults (19 sp, FIA map, Zeide SDI, DGSD 1.7, LHTDRG, seed) — data/centralidaho/*.csv | TODO |
| 2 | Site/habitat: ci/sitset.f + ci/habtyp.f (NI 30 habitat types → OCURHT groups) | TODO |
| 3 | Large-tree DG: ci/dgf.f (western Wykoff DDS) + DGHAB coeffs | TODO |
| 4 | Height: ci/htgf.f | TODO |
| 5 | Crown: ci/crown.f + ci/ccfcal.f | TODO |
| 6 | Small-tree: ci/regent.f | TODO |
| 7 | Mortality: ci/morts.f (Hamilton) + Zeide SDImax self-thin | TODO |
| 8 | Volume: ci/sitset.f VOLEQ via shared NVEL driver | TODO |
| 9 | Full-cycle differential vs FVSci_clean on cit01 | TODO |

## Port strategy (discount)
Clone IE's 15 Julia module files → retarget dispatch to `CentralIdaho` → swap in CI coefficient
CSVs extracted from ci/*.f. Two deltas from IE to watch: (a) **Zeide SDI** (reuse UT's Zeide
SDImax/self-thin path, not IE's Stage); (b) 19-vs-23 species + CI-specific DG/height/regent
coefficient DATA. Validate each chunk bit-exact via instrument-replay (relink_ci.sh) + cit01 .sum.
