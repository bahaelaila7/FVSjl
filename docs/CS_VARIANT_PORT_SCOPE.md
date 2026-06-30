# CS (Central States) variant — port scope

Scope for adding the **Central States (`FVScs`)** variant to FVSjl, written after the NE
port (tag `FVSne-complete`). Grounded in a direct comparison of the FVS CMake source lists.

## TL;DR — CS is the *cheapest* next variant, easier than LS

The FVS `FVScs_sourceList.txt` shares **556 / 587 files (94.7%) with NE** and 547/587 (93.2%)
with SN — i.e. ~94% is the base engine + `common` + FFE base + extensions that FVSjl already
ports **variant-agnostically**. The variant surface that actually differs is **~30 routines**,
and most of them are **near-identical to the NE versions** I already ported:

| routine | CS≈NE | reuse |
|---|---|---|
| `dgdriv` (ARMA serial-corr driver) | 97% | ~drop-in (shared with NE+SN) |
| `cratet` (HT dubbing / broken-top) | 98% | ~drop-in |
| `htgf` (height growth) | **96%** | **reuse NE's height-growth model** |
| `regent` (small-tree/estab growth) | 96% | ~drop-in |
| `varmrt` (variable mortality) | 94% | ~drop-in |
| `crown` (crown ratio) | 88% | mostly reuse |
| `cubrds` (cuft/bdft volume helper) | 84% | mostly reuse |
| `balmod` (BAL competition modifier) | 35% | reuse NE's BAL framework, re-fit |
| `htdbh`/`sitset`/`forkod` | 44/48/48% | partial — CS site/forest-type/HT-DBH tables |
| **`dgf`** (diameter growth) | **13%** | **the one genuinely-new model** |

The structurally-new pieces NE introduced — the **BAL competition** and the **height-growth
model** — carry over (CS shares NE's `htgf` at 96% and uses the same `balmod` framework),
whereas in LS `htgf` is a completely different model. So **CS leans on the NE port harder
than LS does**; the new work concentrates in a single routine (`dgf`) plus coefficient swaps.

CS facts: **MAXSP = 96** (NE 108, SN 90); live oracle buildable from `bin/FVScs_buildDir/*.o`
(551 objects — relink exactly like `ne_oracle.sh`); canonical test stand **`tests/FVScs/cst01.key`**
(+ `cst01_method5.key`).

## The shared ~94% (already done — just runs as CS)

Everything outside `src/variants/` is variant-agnostic and already validated for SN+NE: the
stand loop / keyword dispatch / `each_stand` / `run_keyfile`, density, the cycle driver
(`grow_cycle!`, `mortality_and_fire!` with the tripling-before-fire order + `_fire_due`
OPCYCL gate), the FFE base (snags/CWD/carbon/consumption + the 9 DBS tables), ECON, modern
IO. Adding CS means: register `Northeast`-style `CentralStates <: AbstractVariant`, bump the
per-variant `nspecies` to 96, point `site_setup!`/FORKOD at the CS defaults, and supply the
~30 variant methods below under `src/variants/centralstates/`.

## §1 — CS variant equations to port (`cs/*.f`, 19 files)

Classified by how much reuses NE (`src/variants/northeast/`) vs needs fresh work:

- **Drop-in from NE (coefficient swap only):** `dgdriv` (shared), `cratet`, `htgf` (96%),
  `regent`, `varmrt`, `crown`, `cubrds`, `grinit`, `grohed`, `findag`, `essubh`, `cutstk`,
  `nbolt`. → re-use the NE Julia routine, CSV-swap the CS coefficients.
- **NE framework, CS re-fit:** `balmod` (BAL competition — reuse `ne_balmod` structure, CS
  B3/coefficients), `blkdat` (the CS per-species DATA: YR, SIGMAR, FU/FM/FL, ISPGRP, etc.),
  the `cs/common/{ESPARM,PRGPRM}.F77` parameters (MAXSP 96).
- **Genuinely CS-specific (the bulk of new work):**
  - **`dgf` (the one real new model).** CS DG is the **SN-family ln(DDS) regression** —
    predicted from DBH, site index, crown ratio, **BA percentiles and QMD** for trees ≥ 5″
    (per the routine header) — *not* NE's BAL-potential 10-iteration. So it extends the
    **Southern `dgf!`** framework (`src/variants/southern/diameter_growth.jl`) with CS
    coefficients + CS's own competition term (its `balmod`/percentile inputs), rather than
    `ne_dgf!`. This is the single most important routine to trace both-sides.
  - **`htdbh` / `sitset` / `forkod`** (44/48/48% NE) — CS HT-DBH curves, site-index setup,
    and forest-type code mapping (eastern, closer to NE than SN but with CS tables).

## §2 — Volume (NVEL eastern, like NE)

CS volume rides the same **NVEL** path as NE, not SN's R8 Clark in-source: `cs/cubrds.f`
(84% shared), `cs/nbolt.f`, and **`nc/logs.f`** (the North-Central log routine, which CS and
LS both pull). So FVSjl reuses the **R9 Clark cubic + R9LOGS Scribner board-foot** machinery
already ported for NE (`src/engine/r9clark_vol.jl`), parameterized by CS's VOLEQ defaults /
NVEL equation ids. Confirm CS's default equation numbers (`VOLEQDEF`/region) when porting.

## §3 — CS FFE fire variant (`fire/cs/`, 8 files)

`fmbrkt fmcba fmcblk fmcfmd fmcrow fmcsft fmsfall fmvinit` + CS uses **`fire/sn/fmmois.f`**
(shares SN's fuel-moisture). The FVSjl FFE base (`src/engine/fire/`) is variant-agnostic and
done; this is the CS fuel-model / crown-bulk / scenario-weather surface — the same shape as
the NE fire-model path (`FMCFMD`/`FMDYN` weighted standard fuel models), re-coefficiented.
`fmcsft` (softwood-specific) is CS-only — trace it.

## §4 — Branches inert in SN/NE but ACTIVE in CS (un-gate in shared code)

Audit the shared engine for SN/NE-gated branches that CS turns on (the NE port surfaced a few
of these). Likely candidates: FORKOD/forest-type set, the habitat/ecological-unit decode
(CS uses its own, not SN's `habtyp` nor NE's `ak/habtyp`), and any density/site default keyed
on variant. Gate on the variant/coefficients — **do not** harden; SN+NE must stay bit-exact.

## §5 — CS species data (CSVs) — the mechanical bulk

Add the 96-species CS coefficient tables (the CSV-driven approach from SN/NE): DG/HTG/bark/
crown/volume-eq/SDImax/site/forest-type/sprout coefficients, FIA/alpha species codes, and the
`blkdat` per-species arrays (SIGMAR, FU/FM/FL, ISPGRP). Mechanical but voluminous — extract
from `cs/*.f` DATA statements + the CS coefficient files, validated against the live oracle.

## §6 — Validation plan

Sole oracle = **live `FVScs`**, relinked from `bin/FVScs_buildDir/*.o` (clone `ne_oracle.sh`
→ `cs_oracle.sh`). Drive **`cst01`** to bit-exact, cycle-0 first (TPA/BA/SDI/CCF/QMD/TopHt),
then cycle-1+ growth/volume, then thinning/establishment/FFE — the same ladder as NE. There is
**no** faithful-port (FVSjulia) oracle for CS, exactly as for NE. Tighten tests to `==` once a
column is bit-exact; classify any residual as ULP / documented-divergence (the NE doctrine).

## §7 — Recommended chunk order (most-upstream-first)

1. Variant infra: `CentralStates<:AbstractVariant`, MAXSP→96, CS FORKOD/site defaults, CS
   species CSVs; drive `cst01` cycle-0 stand columns bit-exact.
2. Volume: wire CS NVEL eq ids into the existing R9 Clark + R9LOGS path → cycle-0 .sum volume.
3. Diameter growth: **`dgf` (the new SN-family CS model)** + CS `balmod`/competition + the
   `dgdriv` calibration (shared) → cycle-1 DG.
4. Height growth (reuse NE `htgf` + CS coefs), crown, mortality (`varmrt`), `htdbh`/`sitset`.
5. Drive `cst01` to cycle-1+ vs live; then CS-active shared branches, then FFE (`fire/cs/`).

## What does NOT need porting (shared / done)

The ~94% shared base + the whole variant-dispatch architecture, FFE base, ECON, IO, and the
NE routines CS reuses near-verbatim (dgdriv/cratet/htgf/regent/varmrt/crown/cubrds). Net new
surface ≈ **`dgf` + the CS site/forest/HT-DBH tables + the CS fire fuel models + the 96-species
CSVs** — materially smaller than NE was, because NE already proved the eastern engine, BAL
competition, height-growth model, and the R9 volume path that CS shares.
