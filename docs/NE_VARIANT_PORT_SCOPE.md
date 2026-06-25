# NE (Northeast) variant — port scope

Scoping for adding the **NE** variant to FVSjl alongside SN. Driven by a diff of the Fortran build
source lists (`bin/FVSsn_sourceList.txt` vs `bin/FVSne_sourceList.txt`), the `AbstractVariant` interface
already in FVSjl, and the shared/variant split.

## TL;DR — the effort is much smaller than SN-from-scratch

- **490 of ~528 Fortran files (~93%) are SHARED** between SN and NE — the entire engine, IO, keyword
  processor, event monitor, ECON, DBS, base routines, and the NVEL volume *library* (incl. the R8 Clark
  coefficients). **FVSjl already implements all of this variant-agnostically** and it is bit-exact for SN.
- FVSjl was built for this: `src/variants/variant.jl` defines an `AbstractVariant` interface and the
  comment literally says *"Adding a variant later (CS/NE/LS) = define a new `<: AbstractVariant`
  singleton."* Adding NE = a `struct Northeast <: AbstractVariant` + its method implementations + NE data.
- So the NE port is the **~38 variant-specific files** (the growth / mortality / volume / site / crown
  *equations* and the NE FFE fuel constants), **NE species data (CSVs)**, and **un-gating a handful of
  VARACD='NE' branches** that were inert in SN.
- Oracle is ready: `bin/FVSne_buildDir/FVSne` (the NE executable) + `tests/FVSne/net01.key` (the NE analog
  of `snt01.key`). So the same bit-exact differential discipline applies.

## The shared 93% (already done — just parametrize)
The engine/IO/keyword/event-monitor/ECON/DBS/base/NVEL is shared and variant-agnostic in FVSjl. Adding NE
touches almost none of it. The few seams that need a variant branch are listed in §4.

## §1 — NE variant equations to port (the bulk)
The `ne/` directory + the eastern-borrowed routines. Same STRUCTURE as the SN equivalents FVSjl already
has (Chapman-Richards height, logistic background mortality, Pretzsch-style density), but **different
coefficients** and a couple of structural additions. Each maps to an `AbstractVariant` interface method:

| Subsystem | Fortran (NE) | FVSjl interface method | Notes / difference from SN |
|---|---|---|---|
| Diameter growth | `ne/dgf.f`, `ne/dgdriv.f`, `ie/dgbnd.f`, **`ne/badist.f`, `ne/balmod.f`** | `diameter_growth!` | ★ **KEY DIFFERENCE: NE uses a BAL (basal-area-in-larger-trees) distance-independent competition term** (badist builds the BAL distribution, balmod is the BAL growth modifier). SN has no BAL term — this is genuinely new logic, not just new coefficients. |
| Height growth | `ne/htgf.f`, `ls/htcalc.f`, `ne/htdbh.f` | `height_growth!`, `height_from_dbh` | Chapman-Richards; NE/LS coefficient set. |
| Mortality | `vls/morts.f`, `ne/varmrt.f` | `mortality!` | NE borrows the `vls/` (variant lib) morts; check the background/density split + VARMRT shade table. |
| Crown | `ne/crown.f`, `ne/cratet.f` | `crown_ratio!`, `regenerate!` | NE crown-ratio change model. |
| Small-tree / establishment | `ne/regent.f`, `ne/essubh.f`, `ls/estab.f` | `regenerate!` + establishment | NE borrows `ls/estab.f` (Lake-States establishment). |
| Site | `ne/sitset.f`, `ne/findag.f`, `so/adjmai.f` | `site_setup!` | `adjmai` = adjusted MAI for NFI plots (borrowed from `so/`). |
| Bark ratio | `kt/bratio.f` | `bark_ratio` | NE borrows the `kt/` bark-ratio. |
| Init / constants / classify | `ne/blkdat.f`, `ne/grinit.f`, `ne/forkod.f`, `ak/habtyp.f`, `ne/grohed.f` | `load_species_coefficients!`, `compute_forest_type!`, habitat | `blkdat` = the NE coefficient DATA blocks → CSVs (§5). `forkod`/`habtyp` = NE forest-type / habitat decode. |
| COMMON serialize | `pg/varget.f`, `pg/varput.f` | (existing PUTSTD/GETSTD) | NE field set for stop/restart; FVSjl's serialization is generic, needs the NE variable list. |

## §2 — Volume (the R8 Clark question, answered)
**R8 Clark is SHARED, not SN-specific** — `r8cfo.inc / r8clkcoef.inc / r8clist.inc / r8dib.inc` are in the
**NVEL national volume library, present in BOTH source lists**. So `src/engine/r8clark_vol.jl` (and the
extracted `data/southern/volume/*.csv`) is reusable infrastructure, not SN-only.

BUT *which* equation a stand uses is selected per species/region by the volume-equation number. NE stands
primarily use the **NE/Region-9 path**, not R8 (Region-8 = the South). The NE-specific volume routines are
`ne/cubrds.f` (cubic), `ne/nbolt.f` / `ne/logs.f` (board/log bucking), and `ls/gvrvol.f` (the gross-volume
ratio). So:
- Keep `r8clark_vol.jl` as shared (it already is, in `engine/`, not `variants/southern/`).
- Port the NE volume selection + `gvrvol`/`cubrds`/`nbolt`/`logs` path.
- The shared dispatcher `ie/vols.f` has NE/CS/LS branches (eastern "total" vs "merch" volume) — §4.

## §3 — NE FFE fire variant (`fire/ne/`, ~9 files)
`fmcba.f, fmcblk.f, fmcfmd.f, fmcrow.f, fmmois.f, fmsfall.f, fmvinit.f, fmbrkt.f, fmneft.f`. These are
mostly DATA tables (fmvinit alone has ~626 numeric lines) → **CSVs**, exactly like the SN fire data
(`data/southern/fire_*.csv`). The fire ENGINE (Rothermel, FMEFF, snag dynamics, the carbon report, the 9
DBS tables) is shared and done — NE needs its **fuel models / bark coefficients / snag fall rates /
fuel-moisture scenarios / decay constants** as `data/northeast/fire_*.csv`, plus `fmneft` (the NE analog
of SN's `fmsnft` snag-fall-after-fire).

## §4 — Branches inert in SN but ACTIVE in NE (un-gate in shared code)
These are the "keywords/logic that were inert in SN" the request anticipated. They live in SHARED code,
gated on `VARACD`:

1. **THINRDSL** (relative-density SDI-line / SILVAH RD thin) — FVSjl marks it ⚪ "N/A in SN (RDCLS2 gated
   VARACD≠NE)". It is **ACTIVE in NE**. Needs the RDCLS2 relative-density path in `cuts!`.
2. **`ie/vols.f` eastern volume branches** (lines 179/200/288/320/338/365/437): `VARACD.EQ.'NE'` selects
   total-vs-merch volume handling and the Region-9 equations. Un-gate for NE.
3. **`fire/base/fmsvol.f`** snag-volume has an `…OR. VARACD.EQ.'NE'…` branch (line 149) — affects the FFE
   snag bole volume. FVSjl's `snag_bole_carbon` path needs the NE branch.
4. **FMEFF season + maple mortality adjustments** — `src/engine/fire/fire_effects.jl:110` notes the
   `IF (VARACD.EQ.'LS'/'ON'/'NE')` fmeff season/maple reductions that FVSjl **skips for SN**. These are
   **active in NE** — re-introduce them in the NE fire-effects path.

## §5 — NE species data (CSVs) — the mechanical bulk
NE has **MAXSP = 108** species (vs SN's 90). Every per-species CSV I built for SN needs an NE counterpart
under `data/northeast/`, extracted verbatim from `ne/blkdat.f` + the NE coefficient includes (the same
dump-from-loaded-values discipline used for the SN tables and FAPROP):
- `species_coefficients.csv` (108 rows): DG/HTG coefficients, background-mortality b0/b1, SDI maxima,
  bark ratios, crown-width, volume-eq numbers, the BAL-term coefficients (B3), shade-adjust, etc.
- the site / forest-type / habitat tables, the establishment min-ht/min-diam columns,
- `fire_*.csv` (biomass, fuel models, fuel-moisture, decay, snag props, bark eqn) from `fire/ne/*`.

## §6 — Validation plan
Identical discipline to SN: bit-exact differential vs the live NE Fortran.
- Oracle: `bin/FVSne_buildDir/FVSne` + `tests/FVSne/net01.key` (NE's `snt01`). Build/run it for the
  `.sum` + DBS baselines, exactly as the SN harness does.
- Gate each chunk on the net01 `.sum` row (cycle-0 start-state bit-exact, then per-cycle), as for SN.

## §7 — Recommended chunk order (most-upstream-first)
1. **Northeast variant skeleton** — `struct Northeast`, `variant_code="NE"`, wire `each_stand`/dispatch;
   stub the interface methods. Get net01 to PARSE + initialize.
2. **NE species data CSVs** (§5) — extract from `ne/blkdat.f`; the most-upstream dependency.
3. **NE site + height + diameter growth** incl. the **BAL competition (badist/balmod)** — the one
   structurally-new piece; validate net01 cycle-1 DG bit-exact.
4. **NE mortality + crown + small-tree/establishment** — validate net01 `.sum` density columns.
5. **NE volume** (cubrds/nbolt/logs/gvrvol + the `ie/vols.f` NE branch) — validate net01 volume columns.
6. **NE-active shared branches** (THINRDSL, fmeff season/maple) — turn on + validate with NE thin/fire scenarios.
7. **NE FFE fire data** (`data/northeast/fire_*.csv` + fmneft) — validate the NE fire/carbon stand.

## What does NOT need porting (shared / done)
The IO (.key/.tre + YAML/CSV), keyword processor + event monitor, ECON, the DBS output framework + all 9
FFE DBS tables, the snag/fuel/carbon ENGINE, the NVEL volume library + R8 Clark, the RNG, the cycle driver,
stop/restart serialization framework, `.sum`/report writers. NE reuses all of it. The genuinely-new
NE-only logic is small: **BAL competition (badist/balmod)** + the NE coefficient sets + the NE-active
VARACD branches.

> Net: the SN port built the whole shared machine; NE is "swap the equation coefficients + add one
> competition term + un-gate a few branches + supply 108-species data," validated against an oracle that
> already exists.
