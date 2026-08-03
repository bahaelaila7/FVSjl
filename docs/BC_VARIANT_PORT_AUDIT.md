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
  - **BEC → ITYPE needed for DG**: bc/dgf.f uses `MAPHAB(ITYPE,sp)→ISPHAB→DGHAB(6,sp)` + MAPDSQ/MAPCCF —
    the SAME structure as IE. So habtyp must map the BEC {Zone,SubZone} → an **ITYPE (1..30)** (not just BAMAX).
    SDIDEF=BAMAX/(0.5454154·PMSDIU/100) (Stage). BC is metric-I/O but imperial-internal (ACRtoHA only for display).
    Synthetic-stand shortcut: user BAMAX (>0) skips the zone lookup (LBAMX=false), but ITYPE still needed for DG.
  - **↑ silver lining**: once BEC→ITYPE lands, **DG/height/crown/regent/mortality/volume are direct IE clones**
    (BC dgf = IE MAPHAB/DGHAB/MAPDSQ form; Stage SDI; constant bark). The BEC-zone site is the ONLY real new work.

## Reuse from CI (proven this session)
Same extraction tooling (tools/centralidaho/*.py adapts by variant), same chunk methodology, and the 4
CI bug classes are the watch-list: DG DDS continuation lines (read ALL of them); mortality self-thin type
(BA10 vs Zeide — BC is Stage ⇒ BA10 like IE); mortality BAMAX source; bratio.f per-species branches.

### Chunk 2 CORE — PORTED + RUNS (site_index.jl)
bc_kodtyp_itype (becset.f:1116 KODTYP→ITYPE, 130→1…730→27, NI codes=IE MTYPE, DEFAULT 21) +
SDIDEF=BAMAX/(0.5454154·PMSDIU/100) (Stage). STDINFO field-2→habitat_code wired for BC. Runs:
KODTYP 520→ITYPE 21, SDIDEF 129.42 (BAMAX 60 default). DEFERRED: full BEC-STRING parser (becset.f
~1000 lines) — synthetic BC stands supply numeric KODTYP + BAMAX (all DG/mortality need). Next hook = dgf!.

## ⚠⚠ CHUNK 3 (DG) — SCOPE CORRECTION: BC IS **NOT** AN IE CLONE (source-verified, bc/dgf.f)
The earlier "silver lining: DG/height/… are direct IE clones once BEC→ITYPE lands" was **WRONG**. bc/dgf.f
has **TWO complete DG models** selected by `LV2ATV`:
- **LV2ATV=.TRUE.** (becset: BEC zone NOT in {ICH,IDF,SBS,SBPS} — e.g. ESSF/MS/BWBS): imperial, the IE NI
  form `DDS=CONSPP+DGLD·lnD+DGBAL·BAL+CR·(DGCR+CR·DGCRSQ)+DGDSQ·D²+DGDBAL·BAL/ln(D+1)` using the
  `DGLD/DGCR/DGBAL/DGDBAL` arrays I extracted (these DO match IE). Plus SEILTDG(BECADJ,SEICN2,…) additions.
  **This is the MINORITY path** (non-productive high zones).
- **LV2ATV=.FALSE.** (BEC zone ∈ {ICH,IDF,SBS,SBPS} — the productive interior forests — **AND the empty-BEC
  grinit default**): a **DISTINCT metric "version-3 calibration"** in **centimetres**: CASE DEFAULT
  `DDS=exp(CONSPP+DGLD1·D+DGDSQ1·D²+DGBAL1·BAL+DGDBAL1·BAL/D+DGDBAL2·BAL/ln(D+1)+DGCR1·CR)`; then
  `DDS=(DDS²+2·DDS·D·BRATIO)·CMtoIN²`. Birch(11,15) & aspen(12,13) use power forms `DGLD1·D2^(DGDBAL1+DGDBAL2·BAL/D2)·exp(DGDSQ1·D2²)`.
  This is the COMMON path and is a **different model with different coefficients** — not IE.

**Root of BC's cost — the v3 coefficients are BEC-STRING-DRIVEN table lookups** (bc/dgf.f ENTRY DGCONS):
- `DGCON(I)=SSKONST(JP)%CON+ZNKONST(IP)%CON`; `DGLD1=ZNKONST(IP)%LD`, `DGDSQ1=%DSQ`, `DGBAL1=%BAL`,
  `DGDBAL1=%DBAL1`, `DGDBAL2=%DBAL2`, `DGCR1=%CR`, `DGCCF1=%CCFA` (+ aspect/elev via %CASP/%SASP/%EL/%EL2).
- `IPOS(I)`/`JPOS(I)` are found by string-INDEX-matching the stand's parsed `BEC2%{Zone,SubZone,Series,PrettyName}`
  against **28 `ZNKONST` (MD_STR) zonal records** + **25 `SSKONST` (SS_STR) site-series records**. Each MD_STR =
  {SPP(15) list, ZONE(30) 15-char patterns, OBSERV, CON, CASP, SASP, EL, EL2, CCFA, LD, DSQ, DBAL1, DBAL2, CR,
  BAL, SIGMAR}. RELDN2 caps and CR caps are also `INDEX(BEC%Zone,'ICH'/'IDF')`-gated.
- **The SAME P_SS/P_ZN string-table machinery recurs in htgf.f / morts.f / regent.f** — so BC's height, mortality,
  and regen are ALSO BEC-string models, not IE clones.

**Verdict:** BC's real cost = port the BEC-string classifier (becset.f ~1000 lines: parse "IDFdk1"→{Region,Zone,
SubZone,Series} via RGN/ZN/SZ tables) + the ZNKONST/SSKONST DATA tables (28+25 records × ~16 fields) + the
species×zone string-matching selection + the metric v3 DDS/bratio-cm math — repeated for DG/HTG/MORT/REGENT.
This is a **full variant port dominated by the BEC system**, comparable in cost to a fresh western variant, NOT a
"cheap clone". (The Tier-1 "cheap" projection held for **CI** — done, bit-exact — but not BC.) Chunks 0–2 (scaffold/
species/site-core) remain valid; chunk 3 is where the IE-clone assumption breaks. Extracted this session:
`data/britishcolumbia/dg_*.csv` = the LV2ATV(imperial/IE) arrays — reusable for the minority Seiler path only.

**Recommendation:** BC is a genuine multi-chunk variant port (BEC parser first). Given the "cheaper ones + then
extensions across all variants" directive, CI delivered the cheap win; BC should either be taken as a full port in
its own right or sequenced after the extensions pass on the 11 done + CI. Flagged for the user's sequencing call.
