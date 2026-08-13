# Western Cluster Post-NEWSPRED Readiness Scoping

Date: 2026-08-13
Scope: what remains in the FVSjl western-variant rollout AFTER the in-progress
BC NEWSPRED port (#196) lands, distinguishing PARITY gaps (a model compiled/linked
into the variant's oracle binary but missing from jl) from ADDITIONAL-FEATURE gaps
(an optional model not present in the oracle binary, or present but structurally dormant).

All findings are MEASURED from source. Repos:
- Julia port: `/workspace/FVSjl` (branch kt-variant-port)
- Live FVS Fortran: `/workspace/ForestVegetationSimulator`
- Per-variant oracle build sources: `/workspace/ForestVegetationSimulator/bin/FVS<v>_buildDir/`

Western cluster = CR, KT, IE, EM, BM, TT, UT, CI, BC.

---

## (a) Executive Summary

**After NEWSPRED lands for BC, the western cluster is effectively closed for parity:
there are NO remaining PARITY gaps on the tested regime (growth + volume + the
default-active extensions). The two things left unported are ADDITIONAL-FEATURE gaps,
not parity gaps.**

1. **NEWSPRED is structurally BC/Canada-only.** The spatial dwarf-mistletoe
   (NEWMOD / NEWSPRED) code path exists only in BC's `misin.f` and the 28 `dm*.f`
   files from `canada/newmist`, which are compiled into **only** `FVSbc_buildDir`.
   Every other western variant links the base non-spatial `mistoe.f`, whose `misin.f`
   does not even contain the NEWSPRED keyword. **After BC #196, no other ported variant
   needs NEWSPRED wiring.**

2. **WRD (Western Root Disease, `rd/`) — CAMPAIGN-PREMISE CORRECTION.** The prior
   belief that there are "0 rd*.f in any FVS<v>_buildDir" is **wrong**. WRD **is
   compiled and linked into every western oracle binary** (75 `rd*.f` + matching `rd*.o`
   in each of bc/bm/ci/cr/em/ie/kt/tt/ut). However, WRD is **dormant by default**: its
   activation flag (`RDATV` returns `L = RRTINV .OR. RRMAN`, both default `.FALSE.`,
   set `.TRUE.` only inside the `rdin.f` keyword parser). Reference stands and FIA
   sweeps carry no root-disease keywords, so WRD never activates and therefore never
   perturbs the growth+volume parity the campaign validated. jl has no WRD port.
   Verdict: **cluster-wide ADDITIONAL-FEATURE gap, NOT a parity gap** on any tested stand.

3. **All other default-active extensions are ported.** FFE fire, base non-spatial
   mistletoe, ECON, and Climate-FVS all have jl implementations and are DONE-or-cornered.
   The only other unported thing linked into western binaries is **SVS** (Stand
   Visualization System, ~67 `sv*.f`), which is visualization/reporting output with no
   effect on tree state — an additional-feature gap, never a parity gap.

4. **Insect/pathogen models (DFTM, DFB, MPB, WSBW/BWE, blister rust) are NOT in any
   western binary.** Only their event-monitor keyword recognizers (`exdfb.f`,
   `exdftm.f`, `exmpb.f`, `ESWSBW.F77`) are linked; the model code is absent (0 model
   files). jl matches this exactly (keywords recognized as no-ops in variant dispatch).
   Additional-feature, not a parity gap.

**Genuine remaining gaps (all additional-feature, none block "growth+volume bit-exact-
or-cornered"):**
- **WRD** — cluster-wide, ~23.4 KLOC, dormant without RD keywords, unported in jl.
- **SVS** — cluster-wide, visualization output only, unported in jl.
- **BC NEWSPRED merch/board volume + V2/non-ICH** — pre-existing BC-specific residuals
  noted in memory, independent of #196's spatial-DM growth path.

---

## (b) NEWSPRED Reach Finding

**Question:** Is `canada/newmist` (NEWSPRED / the NEWMOD spatial-DM path in `misin.f`)
reachable by any ported variant other than BC?

**Answer: No. NEWSPRED is structurally BC/Canada-only.**

Evidence:

- `dm*.f` (newmist spatial-DM source) file count per western buildDir:
  ```
  bc: 28    bm: 0   ci: 0   cr: 0   em: 0   ie: 0   kt: 0   tt: 0   ut: 0
  ```
  Only `FVSbc_buildDir` links the spatial model (`dmtreg.f`, `dmntrd.f`, `dminitbc.f`,
  `misintbc.f`, and the `dmcw*.f` crown-width family).

- The NEWSPRED keyword exists only in BC's `misin.f`:
  - `bin/FVSbc_buildDir/misin.f` is the `METRIC-NEWMIST` variant (header line 3) and
    contains: `NEWSPRED` in the keyword table (line 134), `NEWMOD = .FALSE.` init
    (line 141), and Option 12 `NEWSPRED ... NEWMOD = .TRUE.` (lines 595–599).
  - `bin/FVSem_buildDir/misin.f` (and the other non-BC variants) is the base `MISTOE`
    `misin.f`: `grep -niE 'NEWMOD|NEWSPRED'` returns **nothing**. The base file has no
    NEWSPRED keyword and no `NEWMOD` flag, so it cannot activate the spatial model even
    if the `dm*.f` objects were present.

- jl side: the spatial-DM port is BC-scoped —
  `src/variants/britishcolumbia/newspred.jl` (the #196 work in progress). No other
  variant directory references newspred/dmtreg/nisi.

**Conclusion:** After BC #196 completes, no other ported western variant requires
NEWSPRED wiring. The base variants use the non-spatial `mistoe.f` DM model, already
ported (see section (d)).

---

## (c) WRD (Western Root Disease, `rd/`) Finding

**Size:**
- `/workspace/ForestVegetationSimulator/rd/`: **~23,382 lines** across the `.f` sources
  (`wc -l rd/*.f` total), **85 SUBROUTINE/FUNCTION** definitions in the `.f` files,
  plus `.F77` include/common files (`RDCOM.F77`, `RDADD.F77`, `RDARRY.F77`,
  `RDPARM.F77`, `RDCRY.F77`, …) and per-variant block-data `rdblk1<v>.f`.
- Compiled footprint per western variant: **75 `rd*.f` files (+ 75 `rd*.o`)** in each of
  bc/bm/ci/cr/em/ie/kt/tt/ut buildDirs.

**(a) Is it compiled into western buildDirs?**
**YES — this corrects the campaign premise of "0 rd*.f in any buildDir."** Every western
buildDir contains 75 `rd*.f` plus matching `.o` object files (e.g. `rdarea.o`,
`rddisp`-family, `rdinit.o`, `rdgrow.o`, `rdmort.o`, `rdtreg.o`, per-variant
`rdblk1<v>.o`). WRD is compiled and linked into every western oracle binary.

**(b) Is it in jl?**
**NO.** `grep -rilE 'rddisp|rdinit|rootdis|root.?disease'` over `src/` matches only
`src/variants/southern/southern.jl`, which merely lists `RDIN`/`RRIN`/`RDBBMORT`/
`RDSUM`/`RDDETAIL` as keyword recognizers that "recognize, do nothing" (southern.jl:45–50).
There is no WRD model in jl and no engine stub (`grep` of `src/engine`, `src/core` for
`rddisp`/`root_disease`/`wrd` returns nothing).

**(c) Variant-generic vs variant-specific; would the oracles exercise it?**
WRD is **mostly generic with per-variant host-species block data**: the shared engine
(`rdarea`, `rdgrow`, `rdmort`, `rdtreg`, `rdspread`, …) is common, while host-species
susceptibility tables are variant-specific block data — `rdblk1bc.f`, `rdblk1bm.f`,
`rdblk1ci.f`, `rdblk1cr.f`, `rdblk1em.f`, `rdblk1ie.f`, `rdblk1tt.f`, `rdblk1ut.f`
(KT uses the generic `rdblk1.f`). So it is host-supported across the whole cluster.

**Crucially, it is dormant unless activated by root-disease keywords.** The model's
"is-active" test is `RDATV` (`bin/FVSem_buildDir/rdatv.f`):
```
L = RRTINV .OR. RRMAN
LTREE = RRTINV
```
`RRTINV` and `RRMAN` default `.FALSE.` (block data) and are set `.TRUE.` **only** inside
`rdin.f` (the RD keyword parser: `RRMAN = .TRUE.` at rdin.f:584, `RRTINV = .TRUE.` at
rdin.f:2204). `INITRE` calls `RDINIT`/`RDATV` unconditionally, but `IF (RRGO) CALL RDESIN`
and all downstream growth/mortality effects gate on that false flag. With no `RRINIT`/
`RDIN` keyword block in the keyfile — the case for all reference stands and FIA sweeps —
WRD is a no-op and produces byte-identical results to a binary without it.

**Verdict:** WRD is a real *capability* gap (the model is entirely unported in jl,
cluster-wide) but **NOT a parity gap on the validated regime**. The ported-cluster
oracles cannot exercise it under the tested keyfiles, so growth+volume bit-exact-or-
cornered status is unaffected. It is a purely-additional-feature gap.

---

## (d) Other-Extensions Gap Scan

Method: enumerate optional model families linked into a western buildDir
(`FVSem_buildDir` as representative; `.o` files = actually compiled), map each to jl.

| Model / dir | Compiled into western binary? | Default-active? | jl implementation | Gap class |
|---|---|---|---|---|
| **FFE fire** (`fire/`, `fof_*`, `fm*`) | Yes | Yes (when FMIN keywords) | Yes — `src/engine/fire/` (fire_effects.jl, fuel_model.jl) | Ported / DONE-cornered |
| **Base DM** (`mistoe.f`, `misin*.f`) | Yes | Effect only with DM ratings | Yes — CR `dwarf_mistletoe_model.jl`, IE `mistoe_coefficients.jl`, per-variant `mortality.jl` | Ported / DONE-cornered |
| **NEWSPRED spatial DM** (`canada/newmist`, `dm*.f`) | **BC only** | BC, keyword-gated | In progress — `britishcolumbia/newspred.jl` (#196) | BC-only; lands with #196 |
| **ECON** (`econ/`, `ec*.f`) | Yes | Keyword-gated | Yes — `src/engine/econ.jl` | Ported / DONE |
| **Climate-FVS** (`clim/`, `cl*.f`) | Yes | Keyword-gated | Yes — `src/engine/climate.jl` | Ported / DONE-cornered |
| **WRD root disease** (`rd/`, 75 files) | **Yes (all western)** | **No — dormant w/o RD keywords** | **No** | **Additional-feature gap (not parity)** |
| **SVS** (`sv*.f`, ~67 files incl `svs*`) | Yes | Output only, no tree-state effect | No | Additional-feature gap (visualization; not parity) |
| **DBS database output** (`dbs*.f`) | Yes | Output layer | Yes — `src/io/` | Ported (IO) |
| **DFTM / DFB / MPB insects** | **No model files** (only `exdftm.f`/`exdfb.f`/`exmpb.f` keyword stubs) | n/a — not linked | Keyword no-ops (southern.jl style) | Additional-feature; models not even in binary |
| **WSBW / BWE budworm** | **No model files** (only `ESWSBW.F77` keyword stub) | n/a — not linked | Keyword no-op | Additional-feature; not in binary |
| **Blister rust / MPB stand models** (`wpbr/`, `lpmpb/`, `wwpb/`, `dfb/`, `dftm/`) | **No** (0 files in western buildDirs) | n/a | n/a | Additional-feature; not in binary |

Evidence for the "not compiled" rows: in `FVSem_buildDir`, file-family counts are
`dftm:0 dfb:0 mpb:0 wsbw:0 wpbr:0 lpmpb:0 wwpb:0 blister:0 brust:0`; only the
event-monitor recognizer stubs `exdfb.f`, `exdftm.f`, `exmpb.f`, `ESWSBW.F77` are
present. These stubs only parse keywords; without the model objects they cannot run.

**True parity gaps found: none.** Every model that is (a) compiled into a western
binary AND (b) default-active-affecting-tree-state is ported (FFE, base DM, ECON,
Climate). The two linked-but-unported models — WRD and SVS — are either dormant (WRD)
or non-state-affecting output (SVS), so neither is a parity gap on the tested regime.

---

## (e) Per-Variant Readiness Table (after NEWSPRED lands)

Legend — PARITY = model in the variant's oracle binary + affects tested growth/volume +
missing from jl; ADDITIONAL = optional model not exercised by tested keyfiles (dormant
or not linked or output-only).

| Variant | Growth+Vol | Default extensions ported | PARITY gaps remaining | ADDITIONAL-FEATURE gaps |
|---|---|---|---|---|
| **CR** centralrockies | bit-exact-or-cornered | FFE, base DM (own), ECON, Climate | **none** | WRD (dormant), SVS |
| **KT** kootenai | bit-exact-or-cornered (IE-parity) | FFE, base DM, ECON, Climate | **none** | WRD (dormant, generic rdblk1), SVS; FFE/estab minor residual (memory) |
| **IE** inlandempire | bit-exact-or-cornered | FFE, base DM, ECON, Climate, AUTOES | **none** | WRD (dormant), SVS; FFE-fuel residual (memory) |
| **EM** easternmontana | bit-exact-or-cornered | FFE, base DM, ECON, Climate, ESTAB | **none** | WRD (dormant), SVS; FFE residual (memory) |
| **BM** bluemountains | bit-exact-or-cornered | FFE, base DM, ECON, Climate | **none** | WRD (dormant), SVS; minor-species volume residual (memory) |
| **TT** teton | bit-exact-or-cornered | FFE, base DM, ECON, Climate, woodland/R4VOL | **none** | WRD (dormant), SVS; plant-estab sub-inch + DVEW NATCRS residuals (memory) |
| **UT** utah | bit-exact-or-cornered | FFE, base DM, ECON, Climate | **none** | WRD (dormant), SVS |
| **CI** centralidaho | bit-exact-or-cornered | FFE, base DM, ECON, Climate | **none** | WRD (dormant), SVS |
| **BC** britishcolumbia | bit-exact-or-cornered | FFE, **NEWSPRED spatial DM (#196)**, ECON, Climate | **none** (once #196 lands) | WRD (dormant), SVS; NEWSPRED merch/board vol + V2/non-ICH (memory) |

Notes:
- The "residuals (memory)" items are pre-existing small/cornered items already tracked
  in the variant port docs; they are not extension-model parity gaps and are listed only
  for completeness.
- WRD appears in every row because it is linked into every western binary, but it is an
  ADDITIONAL-FEATURE gap everywhere (dormant without RD keywords). If a future goal is to
  run FVS root-disease scenarios through jl, WRD (~23.4 KLOC, 85 routines + per-variant
  host block data) is the single largest unported western subsystem.

---

## Evidence Index (commands run)

- `dm*.f`/`rd*.f` counts per western buildDir (`ls FVS<v>_buildDir/{dm,rd}*.f | wc -l`).
- `grep -niE 'NEWMOD|NEWSPRED|DMTREG'` on bc vs em `misin.f`; `diff` of the two `misin.f`.
- `rdatv.f` full read (activation flag `L = RRTINV .OR. RRMAN`).
- `grep 'RRTINV|RRMAN'` → set `.TRUE.` only in `rdin.f` (keyword parser).
- `initre.f:180–205` (unconditional RDINIT/RDATV call, `IF (RRGO)` gate).
- `wc -l rd/*.f` = 23,382; routine count 85; `wc -l canada/newmist/*.f` = 8,761.
- Insect-model family counts in `FVSem_buildDir` (all 0; only `ex*`/`ES*` keyword stubs).
- jl scan: `grep -rilE` over `src/` for rd / mistletoe / newmist / svs; econ/clim/fire
  modules under `src/engine/`.
