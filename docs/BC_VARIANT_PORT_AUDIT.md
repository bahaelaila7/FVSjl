# BC (British Columbia) Variant Port — Audit Log

FVS British Columbia variant → FVSjl. Second of the "cheaper ones" (after CI). Tier-1 discount of the
western Wykoff cluster — **closest template IE** (Stage SDI, same DG structure; even closer than CI which
is Zeide). Validated bit-exact-or-cornered per chunk vs `/workspace/.bcwork/FVSbc_clean`.

## Measured infra (chunk 0)
| Property | Value |
|---|---|
| Growth engine | **western Wykoff DDS** (dgf.f reads DGLD; NO GEMDG) |
| SDI | **STAGE** (LZEIDE=.FALSE.) — like IE/EM/BM (NOT Zeide) |
| MAXSP | 15 | 
| Cycle / seed | IFINT=10 / 55329 | DGSD=2.0 |
| Template | **IE** (Stage + Wykoff; shared N-Rockies conifers) |

### Species (JSP → FIAJSP) — BC codes
```
 1 PW 119  2 LW 073  3 FD 202  4 BG 017  5 HW 263  6 CW 242  7 PL 108  8 SE 093
 9 BL 019 10 PY 122 11 EP 375 12 AT 746 13 AC 747 14 OC 202 15 OH
```
(PW=w.white pine, FD=Douglas-fir, BG=grand fir, HW=w.hemlock, CW=w.redcedar, PL=lodgepole, SE=Engelmann
spruce, BL=subalpine fir, PY=ponderosa, EP=paper birch, AT=trembling aspen, AC=black cottonwood. First 10
FIA codes = the IE/CI N-Rockies conifers.)

## Oracle (chunk 0 — DONE)
`/workspace/.bcwork/FVSbc_clean` relinked (539 objs + isoc23 shim + **bc_stubs.o**). relink_bc.sh.
⚠ BC buildDir was missing 4 DBS **output** symbols (dbs_fiavbc_trls/atrtls/cutlst, dbsreference) — referenced
in fvs.f but never defined; stubbed empty (output-only, irrelevant to growth) to link. Test stands are
**DB-based** (tests/FVSbc/StandStructure.key + YSM-SkyRanch.key use DATABASE/DSNin FIA SQLite) — a simple
.tre synthetic stand (or the FIA-DB reader path) needed for the per-chunk .sum differential.

## Chunk plan (= CI arc; IE template, Stage SDI so SIMPLER than CI's Zeide)
0. ✓ Scaffold groundwork. 1. ✓ species (15 sp, bark/SIGMAR validated). 2. site/habitat ⚠BC-SPECIFIC.
3. DG (Wykoff). 4. height. 5. crown/CCF. 6. regent. 7. mortality (Hamilton). 8. volume (NVEL). 9. full-cycle .sum.

## Chunk verdicts
- **0** ✓ scaffold (BritishColumbia loads, nsp=15). **1** ✓ species (data/britishcolumbia/*.csv; BARK1 constant
  imap=2, SIGMAR; STAGE SDI; validated .964/.851/.969 + .1907/.2679 vs live).
- **2 site/habitat ◐ SCOPED — BC-SPECIFIC (main IE departure)**: bc/habtyp.f parses a **BEC (Biogeoclimatic
  Ecosystem Classification) STRING** — e.g. "IDFdk1" — via INDEX into RGN(3)/ZN(7)/SZ(19) tables → a
  {Region,Zone,SubZone} record; bc/sitset.f then `SELECT CASE(BEC%Zone)` → `CASE(SubZone)` → `CASE(iSeries)`
  → **BAMAX** (55/58/60/89/53/…). NOT the IE numeric NI-habitat model. So the STDINFO habitat field is a BEC
  string, and site defaults are zone-keyed. Port = the BEC-string parser + the zone/subzone/series BAMAX tables
  + site defaults (bc/sitset.f). Also: BC test stands are FIA-DB-based ⇒ build a synthetic BC .tre + STDINFO
  with a known BEC zone to validate SDIDEF/BAMAX vs FVSbc_clean. This is BC's largest single departure.

## Reuse from CI (proven this session)
Same extraction tooling (tools/centralidaho/*.py adapts by variant), same chunk methodology, and the 4
CI bug classes are the watch-list: DG DDS continuation lines (read ALL of them); mortality self-thin type
(BA10 vs Zeide — BC is Stage ⇒ BA10 like IE); mortality BAMAX source; bratio.f per-species branches.
