# NE variant — chunk plan & data-extraction map (working tracker)

Goal: port FVS NE into FVSjl with the SN standard + methodology (most-upstream first, port faithfully
from the Fortran, validate bit-exact per chunk). Sole oracle: live `bin/FVSne_buildDir/FVSne` +
`tests/FVSne/net01.sum.save` (FVSjulia/Oracle-A has NO NE). Scope: `docs/NE_VARIANT_PORT_SCOPE.md`.

## Validation target — net01.sum.save, stand 1 (UNTHINNED), cycle-0 (1990) row
`1990  60  536  77  160 176  63  5.1  1558  1347  292  1633 …`
= year, age, **TPA 536, BA 77, SDI 160, CCF 176, TopHt 63, QMD 5.1**, cuft 1558, mcuft 1347, scuft 292,
bdft 1633. The cycle-0 row needs only **tree parse + density + volume** (no growth) — the first gate.

## Architectural findings (uncovered by driving net01 through the load path)
1. **Tolerant loader — DONE.** `load_species_coefficients` now treats per-subsystem/FFE CSVs as optional
   (empty if absent) so a variant-in-progress loads with what it has; SN unaffected (all files present,
   suite 4392+10). **NE now loads its 108 species.**
2. **`MAXSP=90` is Southern-baked (NEXT architectural blocker — ATTEMPTED, harder than it looks).**
   Every per-species array (`plot.sp_*`, the coefficient vectors) is sized `MAXSP=90`; NE has **108**.
   **Attempted the obvious fix** (bump `MAXSP=108` + zero-pad SN's coefficient/code vectors to 108) and
   it **broke SN in 57 tests** — and not subtly: a `.sum` row came out TPA `0` (vs 536) while SDI was
   `189` (vs 160) in the SAME row, i.e. a structural/column/count coupling shifted, not just an
   out-of-bounds. So `MAXSP` is wired into SN behavior more deeply than array capacity (some `1:MAXSP`
   loop's iteration count or a MAXSP-derived quantity feeds an SN result). **Reverted** to keep SN green
   (4392+10). CONCLUSION: the generalization is a genuine, careful refactor — each of the 57 couplings
   must be traced, and the right design is almost certainly **per-variant species count** (a `nspecies(v)`
   the shared loops iterate, with arrays sized to the MAXSP capacity) rather than a blanket bump. This is
   its own validated chunk and the hard gate before any NE growth (`site_setup!`/DG/mortality) can run on
   net01 — they all touch per-species arrays.
3. **NE site model is structurally different** — `ne/sitset.f` converts SI between species via a 28×28
   `SICOEF` matrix + the `IPOINT` group map (extracted), not SN's `site_species`/`master_group`. So
   `site_setup!(::Northeast)` is a genuine port, not a CSV reshape. (NE uses Zeide SDI like SN —
   `ne/grinit.f:129 LZEIDE=.TRUE.`; RNG seed 55329, same as SN.)

## ★ net01 ORACLE VALIDATION — cycle-0 stand state BIT-EXACT (6 columns)
`net01` cycle-0 now matches `net01.sum.save` exactly on **TPA 536, BA 77, SDI 160, CCF 176, QMD 5.1,
TopHt 63** (test/integration/test_net01.jl, 8 asserts). **CCF landed this pass** via two faithful fixes:
(1) the per-species crown-width equation map (`data/northeast/crown_width_species.csv`) was re-extracted
**mechanically from `bin/FVSne_buildDir/cwcalc.f`** (`SELECT CASE(NSP(ISPC,1)(1:2))` → forest/open CWEQ;
165 species codes, zero transcription) — the prior heuristic (default=01/open=03) was wrong (e.g. YB open
is 37101 not 37102; HI was unmapped). (2) The **NE FORKOD default location** (`ne_forkod_defaults!`,
ne/forkod.f): net01 carries no STDINFO/LOCATION, so the variant default `IFOR=2` → lat=43.53/long=71.47/
elev=20 must be applied; without it lat/long=0 drove the Hopkins bioclimatic index to −250 (vs +13.28),
inflating every Bechtold crown (CCF 227). The reused national `crown_width_equations.csv` was verified
correct against the NE source (eqn 01201 etc. match: linear-D Bechtold Model 2, no D²/HI). NOTE: the
`crown_width_equations` are shared-national and confirmed identical for NE — only the species→equation
**selection** is variant-specific. CANONICAL TEST RUNNER: `julia --project=. test/runtests.jl`
(**4400 pass + 10 broken**); `Pkg.test()` reports spurious `NTuple{9} index[0]` errors from a
stale-sandbox-manifest artifact, NOT a code bug (runtests.jl is green).

### (historical) FIRST net01 validation — TPA/BA/QMD/TopHt
`net01` cycle-0 matched `net01.sum.save`: **TPA 536, BA 77, QMD 5.1, TopHt 63**
(test/integration/test_net01.jl). This validates the whole chain landed: CR-tolerant IO + roster +
translation + init + MAXSP generalization + SICOEF site model + the shared density. The breakthrough was
an IO bug: net01.key uses **old-Mac CR-only line endings** (133 CR, 0 LF) → `readline` read it as one
4192-byte line → total desync; fixed `KeywordReader` to normalize CR/CRLF/LF→LF (SN unaffected, suite
4398+10). Remaining cycle-0 columns (SDI 160, CCF 176, cuft 1558/mcuft 1347/scuft 292/bdft 1633) need the
density-SDImax/crown + volume subsystems; cycle-1+ needs growth.

### Next: the VOLUME subsystem (cycle-0 cuft/mcuft/scuft/bdft) — CORRECTED scope (a real NE port)
INVESTIGATED: NE volume is NOT a shared-engine data swap. Two findings from the code:
- `setup_volume_equations!` (volume_equations.jl) is **Region-8-specific**:
  `vol_eq = (iregn==8 && ifia>0) ? _r8_ceqn(...) : blank`. So for NE (Region 9) it assigns NO equation —
  **NE does not use R8 Clark.** R8 Clark is shared *infrastructure* but NE uses the **Region-9 / eastern
  volume model**: `ne/cubrds.f` (cubic) + `ne/nbolt.f`/`ne/logs.f` (board/log bucking) + `ls/gvrvol.f`
  (gross-volume ratio), via the `ie/vols.f` NE branch. This is a genuine subsystem port, not a CSV.
- Merch specs vary by softwood/hardwood (4 distinct rows in SN); NE needs its own grouping from the
  eastern merch defaults (`ls/vols.f` / MRULES).
So the volume chunk = port the NE cubic/board volume equations (cubrds/nbolt/gvrvol) + the eastern merch
specs + the `ie/vols.f` NE branch + a `compute_volumes!`/`setup_volume_equations!` variant-dispatch.
Substantial — the "reuse R8 Clark" shortcut does NOT apply to NE.

## Chunk status
- [x] **C1 — skeleton + roster.** `struct Northeast`, `variant_code="NE"`, `NE_DATADIR`; the 108-species
  roster `data/northeast/species_translation.csv` (alpha/FIA/PLANTS from ne/blkdat.f JSP/FIAJSP/PLNJSP).
- [~] **C2 — NE data foundation (in progress).** `data/northeast/species_coefficients.csv` started with the
  per-species blkdat blocks: `estab_min_ht` (XMIN), `estab_hht_max` (HHTMAX), `dg_resid_sd` (SIGMAR).
  **Reused the 5 genuinely-national CSVs verbatim** (faithful — not variant-specific): `species_translation`
  (the 563-species master crosswalk), `crown_width_equations` (equation FORMS keyed by name),
  `stocking_coeffs`, `fia_stocking_map`, `forest_locations`. **Still NE-specific & pending** (the C3–C7
  per-subsystem extractions): `site_species`, `site_master_group`, `valid_habitat_codes`,
  `forest_type_codes`, `crown_width_species`, `crown_ratio_coeffs`, `htdbh_coeffs`, `merch_specs`, plus the
  growth/mortality/volume coefficient COLUMNS of `species_coefficients`. The loader needs ALL present, so
  net01 won't LOAD until each subsystem's data lands — there is no shortcut to a cycle-0 validation that
  bypasses the crown/site/volume coefficients. This is why C3–C7 are genuine per-subsystem porting, not a
  data-copy.
- [ ] **C3 — site + height + diameter growth (incl. BAL competition badist/balmod).** net01 cyc-1 DG.
- [ ] **C4 — mortality + crown + small-tree/establishment.** net01 density columns.
- [ ] **C5 — volume** (cubrds/nbolt/logs/gvrvol + ie/vols.f NE branch). net01 cuft/mcuft/scuft/bdft.
- [ ] **C6 — NE-active shared branches** (THINRDSL/RDCLS2, FMEFF LS/NE/ON season+maple, fmsvol NE).
- [ ] **C7 — NE FFE fire data** (`data/northeast/fire_*.csv` from fire/ne/* + fmneft).

## Data-extraction map — which CSV/column comes from which NE Fortran routine
The species CSV accumulates columns as each subsystem is ported (exactly how SN's was built). Extract
VERBATIM from the loaded Fortran DATA values (dump → CSV; zero transcription risk), as for SN/FAPROP.

| CSV / column | NE Fortran source | for chunk |
|---|---|---|
| species roster (alpha/fia/plants) | `ne/blkdat.f` JSP / FIAJSP / PLNJSP | C1 ✅ |
| estab_min_ht / estab_hht_max | `ne/blkdat.f` XMIN / HHTMAX | C2 ✅ |
| dg_resid_sd | `ne/blkdat.f` SIGMAR | C2 ✅ |
| DBHMID (10 dbh class midpoints) | `ne/blkdat.f` DBHMID | C2 (shared const) |
| JTYPE / forest-type codes | `ne/blkdat.f` JTYPE, `ne/forkod.f` | C3 |
| site-index curves / mapping | `ne/sitset.f`, `ne/findag.f`, `so/adjmai.f` | C3 |
| height-growth coeffs (Chapman-Richards) | `ne/htgf.f`, `ls/htcalc.f`, `ne/htdbh.f` | C3 |
| diameter-growth coeffs + **BAL term (B3)** | `ne/dgf.f`, `ne/dgdriv.f`, `ie/dgbnd.f`, **`ne/badist.f`/`ne/balmod.f`** | C3 ★ |
| bark ratio | `kt/bratio.f` | C3 |
| background-mortality b0/b1, density (Pretzsch) | `vls/morts.f`, `ne/varmrt.f` | C4 |
| VARMRT shade-adjust | `ne/varmrt.f` | C4 |
| crown-ratio (Weibull) + crown-width | `ne/crown.f`, `ne/cratet.f`, cwidth | C4 |
| SDImax / SDI defaults | `ne/dgf.f` / blkdat | C4 (needed for net01 SDI col) |
| regen min-diam, establishment | `ne/regent.f`, `ne/essubh.f`, `ls/estab.f` | C4 |
| volume eqn numbers + cubic/board coeffs | `ne/cubrds.f`, `ne/nbolt.f`, `ne/logs.f`, `ls/gvrvol.f` | C5 |
| FFE fuel models / bark / snag / moisture / decay | `fire/ne/*` (fmcfmd/fmbrkt/fmsfall/fmmois/fmvinit/fmcblk) | C7 |

## NE-active branches to un-gate in SHARED code (already-ported SN code, add the NE path)
- `cuts!` — **THINRDSL** (RDCLS2 relative-density SDI-line thin); SN marks ⚪ N/A.
- `fire_effects.jl` (FMEFF) — the LS/NE/ON season + maple mortality adjustments SN skips.
- volume dispatch (`ie/vols.f` lines 179/200/288/320/338/365/437) — NE total-vs-merch + Region-9 path.
- `fmsvol.f:149` — NE snag-volume branch.

## Variant-interface methods to implement for `Northeast` (`variants/variant.jl`)
`load_species_coefficients!`, `site_setup!`, `height_growth!`, `height_from_dbh`, `diameter_growth!`
(+ BAL), `mortality!`, `crown_ratio!`, `regenerate!`, `form_class`, `bark_ratio`, `max_sdi`,
`variant_noop_keywords`. Mirror `src/variants/southern/` file-by-file into `src/variants/northeast/`.

## Key structural difference from SN
NE diameter growth adds a **BAL (basal-area-in-larger-trees) distance-independent competition** term
(`badist` builds the BAL by species/size; `balmod` is the growth modifier) — SN has no equivalent. This
is the one genuinely-new growth mechanism; everything else is the same model structure, different coeffs.
