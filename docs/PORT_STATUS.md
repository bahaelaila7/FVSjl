# FVSjl — port & validation status

_Last updated 2026-08-31. `master` tracks the validated state (active work on the
`fia-resweep` branch off master). Validation doctrine: **[DOCTRINE.md](DOCTRINE.md)**.
Hard gate: `test/integration/test_multicycle.jl` = **339 pass / 11 broken,
byte-identical**._

FVSjl is a Julia reimplementation of the USFS Forest Vegetation Simulator — a drop-in
replacement for the live Fortran FVS (same `.key`/`.tre` in, same SQLite/`.sum` out).
The validation doctrine is **bit-exact-or-cornered vs the live relinked Fortran
oracle** (`FVS{v}_g16`, gfortran-16), measured per subsystem/chunk, never inferred.

**"Cornered"** means a divergence that was *measured* and reduced to a named
floating-point primitive — the #206 OLDRN serial-correlation growth straddle (active
when DGSD≥1), the RDPSRT unstable-quicksort self-thin tie-break, a DGSCOR/volume ULP —
**not** an unexplained difference. (Caveat, learned the hard way: a "cornered" label
is a best-effort verdict, and re-measurement has occasionally found a real bug hiding
behind one — so corners are periodically re-audited against the live oracle.)

## Geographic variants (24) — all ported & validated

| Cluster | Variants |
|---|---|
| Eastern | Southern (SN), Northeast (NE), Central States (CS), Lake States (LS) |
| Western Rockies | Central Rockies (CR), Kootenai (KT), Inland Empire (IE), Eastern Montana (EM), Blue Mountains (BM), Teton (TT), Utah (UT), Central Idaho (CI) |
| Pacific / coastal | Central California (CA), East Cascades (EC), West Cascades (WC), West Sierra (WS), South Central Oregon (SO), Klamath/NC, Pacific Northwest (PN), Oregon Coast (OC, ORGANON), Olympic (OP, ORGANON) |
| Other | British Columbia (BC), Ontario (ON), Southeast Alaska (AK) |

Each has growth + volume, and most have FFE/ECON/mistletoe/Climate/establishment,
validated bit-exact-or-cornered vs the live oracle per subsystem.

## FIA behaviour-compat validation

**Eastern four — EXHAUSTIVE (full FVS-ready FIA population).** Every stand projected
the full horizon, all 10 `.sum` columns/cycle vs freshly-relinked live FVS
(`docs/fia_fullscale_results.md`, `data/fia_sweep.db`):

| Variant | Stands | bit-exact-or-cornered |
|---|--:|--:|
| SN | 633,628 | 99.994% |
| NE | 178,148 | 99.991% |
| CS | 255,951 | 99.986% |
| LS | 400,649 | 99.988% |
| **Total** | **1,468,376** | **99.990%** |

The 77 residual `needs_dig` all classify to named cornered primitives; the 60
`live_crash` are cases where **live FVS itself** SIGFPEs on extreme FIA geometry while
FVSjl runs clean.

**Western / non-eastern — clean full-population sweep IN PROGRESS on final code.** 15 of
the 20 non-eastern variants carry an FVS-ready FIA population (688,903 stands); 5 have zero
FIA (KT/BC/OC/OP/ON — Canada/no-FIA/ORGANON-BLM, N/A). The sweep runs like the eastern one
— every cycle's `.sum` vs freshly-relinked live FVS.

An earlier sweep surfaced and fixed **11 real bugs** the prior sampled validation could not
reach: the EM bare-plot AUTOES under-production (EZCRUISE `INADV=1` skips the ESB
inventory-stocking calibration — the dominant western bug, ~98k EM stands), two NWCMRT
density-mortality omissions (NC, UT), six alpha-`PV_CODE` habitat-decode bugs
(PN/WC/CA/SO/NC/EC → correct SDIMAX), a 3-cause BM seedling small-tree-growth bug, and a
CA/SO forkod forest-index crash.

That earlier run's DB mixed code versions across the fixes (and was inflated by a
discipline lapse that let the dig-queue balloon), so it is **archived, not trusted as
final**. The sweep is now being **re-run clean on the final fixed engine** under the
cap-and-fix discipline of [DOCTRINE.md](DOCTRINE.md): it pauses at ~100 unexplained
divergences (`DIGCAP=100`), each is dug to a named cornered primitive or a real bug, real
bugs are fixed upstream-first, and only then does the sweep resume — with a status-flip
ledger re-checking both directions after every fix. Per-variant bit-exact-or-cornered
numbers are re-established as the clean sweep progresses, converging toward the eastern
standard; the harness is durable/resumable (state on `/workspace/.wt-western`).

## Extensions

FFE fire, FVS-Climate, **WWPB beetle + PPE landscape**, Western Root Disease, dwarf
mistletoe, budworm/tussock/beetle insect models, COVER, Event Monitor, ECON, DBS
database output, establishment — all validated bit-exact-or-cornered.

**PPE (Parallel Processing Extension) landscape** — the recovered-source harness
(deleted from the FVS tree in 2014, rebuilt from git history at `bc6e2377^`; a runnable
historical `FVSppe` oracle lives at `/workspace/.ppework/FVSppe`):
- **mode-1** (independent per-stand projection + area-weighted aggregation) — the
  CMADDS/CMPRT2 composite aggregation is **bit-exact vs the FVSppe COMPOSITE table**;
- **mode-2** (interstand beetle dispersal) — the spatial-redistribution kernels
  (`bmatct_multi!`/`bmdrv_multi!`) are **bit-exact vs pristine `bmatct.f` goldens**;
  live in-run cross-stand coupling is **implemented** (`ppe_run_landscape_live!`, a
  Julia-Task/Channel lockstep barrier), equivalence-validated (live in-flight ==
  premade-decisions replay, bit-exact).

**PPE MXHRVP (multistand harvest scheduling)** — the landscape harvest-flow allocator is
now **ported and oracle-validated** (`ppe_run_landscape_harvest!`, the hvaloc.f coordinator).
The kernels are bit-exact vs gfortran-16 driver-goldens over the historical `FVSppe`
(`HVSEL` greedy priority-ranked target-constrained selection; `HVCCUT`; `LBMEMR`/`LBUNIN`
label sets) and the event-monitor evaluator was extended with the PPE policy variables
(PTSTV1 + `SELECTED`). End-to-end vs authored `MSPOLICY` keyfiles run through the historical
`FVSppe`: the HVSEL selection table (per-stand PRIORITY/CREDIT/SELECT + SELECTED RESOURCE) is
bit-exact at cyc0 and across the full/partial/not + differential (HVYLDS≠HVTHIN) branches,
with the cyc1+ numeric residual cornered to the EC-variant before-thin-BA growth straddle
(inherited by CREDIT=BBA) and the equal-priority tie-break cornered to the RDPSRT unstable
sort (#206-class). The `hvreps` composite-**materialization** (the actual harvest effect) is
**oracle-validated**: the selected harvest fires through the base-FVS activity-group-label path
(`AGPLABEL`; hvreps keeps the MSPLABEL for selected stands), and the post-harvest COMPOSITE
(`_ppe_aggregate` over the thinned per-stand rows) is A/B'd vs the `FVSppe` `msp_thin` golden —
post-thin TPA **bit-exact** every cycle (536→36→35→35; the harvest genuinely removes trees to
residual 40 TPA and the landscape regrows identically), cyc0 volumes ~bit-exact, cyc1+ volumes
cornered to the EC post-thin-regrowth growth straddle. (An earlier "THINBTA over-thin" was traced
to a keyfile column-misalignment artifact — FVSjl's base `THINBTA` is correct and matches the
oracle with a column-aligned card.) `IHVEXT=1` external selection, `hvproj`, and `LHVMXC`
max-contiguous-clearcut are deferred-documented.

## Known exceptions / not-yet-closed

- **ADDTREES** (ESTAB opt 28, `estb/esaddt.f`) — **PORTED + oracle-validated** (staged-read A/B vs live
  `FVSie_g16`). The in-tree code is a *bridge* to an external regeneration-model executable: the keyword
  (esin.f opt 28) schedules an activity `432`; when it fires at the ESNUTR establishment seam
  (`addtrees_bridge!`, establishment.jl) FVSjl writes the `<KWDFIL>_<NPLT>_<KDT>_<VARACD>.es1`
  stand-summary, runs `SYSTEM(CMDLN)` (the external model), then reads back the model's `.es2` **activity
  block** (`OPRDAT` → first line IKEEP, then `IACTK IDT NPRMS PRMS…` records keyed by the stand id until
  `End`) and OPADD-schedules each activity — `430`/`431` route into the already-validated PLANT/NATURAL
  regen path (`s.control.schedule`, the `due` filter in `establish!`). Only the external model exe is
  out-of-tree. **A/B (test_ie_addtrees.jl + manual FVSie_g16):** `ADDTREES(es2:430)` is **byte-identical
  to the direct `PLANT` keyword on BOTH the oracle and FVSjl** (the bridge injects exactly a native PLANT
  card), and a **null bridge** (no `.es2`) is **byte-identical to the same packet with no ADDTREES card on
  both sides**. The oracle-vs-jl gap on the *plain PLANT path itself* (a synthetic bare NOTREES stand:
  oracle 911 vs jl 364 TPA @2002) is the **pre-existing IE bare-plot establishment straddle**
  (EZCRUISE/INADV) — not introduced by the bridge, which is bit-exact to the PLANT path. IMET≠1 (in-tree
  data-base variants of ADDTREES) has no in-tree effect (esin.f falls through) and is not wired.
- **ON database-read path** — reader VALIDATED on a real ON DB (`FVSDataHardwood.db`, stand
  `LD3001`, 94 metric trees) via **inline-equivalence**: the DB path and the inline `.tre` path
  share `ingest_tree_records!`, and on LD3001 they produce byte-identical tree state (species /
  cm→in DBH / m→ft height / id / plot-index), with DB TPA == inline TPA × `ACRtoHA` (per-ha→
  per-acre) bit-exact; stand attrs (age/slope/aspect/elevation/latitude) match the DB row
  (`test_ontario_db_path_equiv.jl`). Surfaced+fixed a real reader bug: a **blank ElevFt string**
  (this DB stores unset numerics as `""`, which `_fia_present` counts as present) dropped the real
  metres `Elevation`; now gated on a NUMERIC value so it falls through to `ELEVATION` (moves the ON
  DB cyc0 `.sum` col 269→258). A **live DB oracle remains blocked**: FVSon_g16 SIGSEGVs in
  `dbstreesin_`'s prologue (`__memset_avx2`) on the SQLite tree-read (stand-read of 86 cols
  succeeds) — reproduces with fresh gfortran-16 objects, the shipped Jun-4 `dbstreesin.o`, and a
  static link; gcc-15/gfortran-15 (the workaround for this toolchain-skew class) is unavailable and
  not installable here, and the isoc23 sscanf shim is unrelated to the memset fault. The archived
  `Hardwood.sum.save` is from a DIFFERENT FVS build (MCuFt=178 vs FVSon_g16/FVSjl's faithful 0), so
  not a valid FVSon_g16 oracle. ⇒ measured toolchain block + inline-equivalence, not a live DB A/B.
- **Western full-population FIA sweep** — clean re-run on final code under way
  (cap-and-fix, `DIGCAP=100`); not yet at the eastern exhaustiveness.
