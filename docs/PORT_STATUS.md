# FVSjl — port & validation status

_Last updated 2026-09-21. `master` tracks the validated state (active work on the
`bm-regime-close` branch off master). Validation doctrine: **[DOCTRINE.md](DOCTRINE.md)**.
Gates: the **tiered integration suite** (`test/integration/test_tiered.jl`, see
[test/harness/tiered/README.md](../test/harness/tiered/README.md)) and
`test/integration/test_multicycle.jl` = **350 pass / 0 broken** (TPA and cuft are now compared
at print precision, like BA/SDI/QMD, so a scenario that prints what live prints passes)._

## How this is gated (2026-09-20)

The multicycle gate is a narrow SN-only smoke test — 10 scenarios, 5 columns, 2 thinning
keywords. Everything it never touched used to pass silently: a whole variant's FFE fuel loop
never running, the ECON summary table never being written, blister rust being inert, a
regeneration date landing a cycle late. The **tiered suite** is now the real gate:

* 12 stratified FIA stands × 11 regimes (none / thin / salvage / plant, both date forms /
  simfire + mistletoe / climate / root disease / econ / cover) per variant, with **live-oracle
  goldens** for the `.sum` AND the DBS output tables, compared at print precision;
* a **signed population tally** per variant × regime — a one-directional bias fails even when
  every individual cell is small;
* a **known-residual allowlist** where every entry is either `CORNER` (with a per-record proof)
  or `OPEN` (with the tracked bug). An unlisted mismatch fails; an entry that stops matching
  fails ("unexpected pass — remove it"); a missing fixture fails rather than skipping;
* **aggregate statistics** printed beside the strict per-cell result: match rate per variant and
  per variant/regime, the magnitude distribution, and the columns carrying the most cells.

Baseline on 2026-09-20 (2,506,861 compared values): **BM 94.1% · SN 92.0% · IE 77.0% ·
EM 72.7%**, 80.4% overall. Those rates are per-tree and per-table, far stricter than the
stand-level `.sum` — which is the point: `.sum` aggregates hide per-tree divergence.

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

**Per-variant certification (the honest state, 2026-09-21).** "Ported" is not "certified":

| Variant | State |
|---|---|
| **BM** | At the oracle floor: **398/400** stratified FIA stands bit-exact with no management (the other two differ by 1 in a single printed cell), **zero one-directional bias in all 10 regimes**. ~20 faithful fixes in the 2026-09-19/20 campaign. |
| **IE** | **RE-OPENED.** The 2026-09-19 "closed, no caveat" claim was retracted — 7 items were real bugs (see `IE_VARIANT_PORT_AUDIT.md` → RETRACTION). All fixed; re-certification requires the tiered suite to pass on the master tip. |
| **SN** | Tiered baseline 92.0%; residuals tracked as OPEN, not yet dug. |
| **EM** | Tiered baseline 72.7% — one-directional volume bias in every regime. Its "COMPLETE-HONEST" claim is false; EM campaign owns it. |
| others | Ported and subsystem-validated, but not yet swept through the tiered suite. |

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
oracle with a column-aligned card.) **`LHVMXC` max-contiguous-clearcut** is now ported
(`spcntg!`/`hvcntg!`/`hvsel!` veto, bit-exact vs a gfortran-16 driver-golden) and **`IHVEXT=1`
external selection** is ported (the deterministic selection rule bit-exact, and the external
read-back path validated via a staged `PPE_FFERdAccess.txt`). Only `hvproj` (the external-only
project-ahead) remains deferred.

## Westside shared fixes (2026-09-23, master `b56bfeda`)

Rewriting the root-disease fixture tests from `rd − ctrl` **deltas** to **absolute** live rows
removed the cancellation that had been hiding *control-stand* bugs: PN was +62 BA over live on
the control run, SO +14, BC's long-documented "+36% baseline-BA straddle" was a real bug, and
EC/CI/NC/WS all carried smaller one-directional offsets. **A delta test is not a validation of
either side** — that is the reusable lesson.

Digging those out surfaced six defects that were each a *shared* FVS mechanism implemented as a
hand-kept per-variant whitelist in FVSjl, drifted out of sync with the Fortran:

1. **Dead-record crown dub** — `crown.f` `DO 79 I=IREC2,MAXTRE` dubs crowns on the DEAD records
   too, in all 23 variants; FVSjl ran it for BM/IE/PN only. The missing records skip their
   `DUBSCR` `BACHLO` draws, so the whole downstream DGSCOR stream is offset (PN: 3 draws).
   Now one shared `dub_dead_crowns!` in `src/engine/crown_init.jl`.
2. **ATTEN** (`dgdriv.f:554` SIGMA pooling) was never set for PN/WC/NC/EC/CA.
3. **RELDEN** was set only for whitelisted variants, so PN/WC/NC/CA read 0 — which pins the crown
   `SCALE` at 1.0 and starves the LPMPB/COVER/DFTM/establishment consumers.
4. **Species-major draw order** — the REGENT `ZZRAN` loop walked storage order; Fortran walks
   species-major via `IND1`. Now `species_major_order(s)`.
5. **Point density** — `point_density!` had no SO/WS branch, so those variants fell back to the
   generic crown width and computed `PCCF` ~100× low.
6. **One `variant_bratio`** — DGDRIV, UPDATE, backdating, calibration and the MORTS `DQ10` each
   carried their own bark-ratio dispatch. UPDATE lacked SO/CA/CI/NC-redwood/BC/ON/WS; MORTS was
   CR/BM/EC/CA-only. This was the BC bug.

Plus: the per-variant **RMAI `grinit` default** (50.0 in 16 variants, 0.0 in BC/CI/IE/KT/SN/WS)
was set nowhere but SO; WS `DUBSCR` was a fixed-10 stub; NC `DUBSCR` consumed no draw and the
redwood branch used stand rather than point density; CI `TEMMAI` was unported.

Four goldens moved and were reconciled **honestly against live**, not re-pinned to silence:
EC `ecvol` col-24 ACCRETION is now `@test_broken` (its per-tree cycle-1 diffs went 27 → 12, i.e.
*closer* to live), BC `dmntrd` SDI 1303 → 1241 (live 936) and CI simfire pre-fire TREES 304 → 302
(live 300) are jl self-snapshots re-pinned in the direction of live.

**Still open on the fixtures** (absolute rows, control + rd): BM exact; CR/PN/SO/BC control ±1–2;
OPEN = EC (TPA 4), CI (BA 5), NC (TPA 10), EM (BA 15), TT (TPA 15), UT (TPA 11), WS (TPA 28/BA 10),
KT (BA 6/TPA 8).

## Crown initialisation, three shared fixes (2026-09-23, master `8b31b046`)

The westside pass above left three more defects in the same area, each a *shared* FVS mechanism that
FVSjl had implemented per-variant and let drift:

1. **SDIAC was 0 at every LSTART crown dub.** `base/fvs.f:193-196` calls
   `SDICLS(0,0.,999.,1,SDIAC,SDIAC2,…)` immediately before `CALL CRATET`, with the source comment
   *"SDICLS IS CALLED HERE SO CROWNS WILL DUB CORRECTLY IN VARIANTS USING THE WEIBULL DISTRIBUTION"*.
   FVSjl passed no `crown_sdi` there, so `RELSDI = SDIAC/SDIDEF` was 0 and `ACRNEW` sat at its maximum
   `C0` for every Weibull-dubbed inventory crown, in every variant.
2. **OC never ran the LSTART crown dub at all.** `oc/cratet.f:831-856` calls `CROWN` when a live or a
   cycle-0 dead record still lacks a crown after ORGANON PREPARE, which only reloads the `IORG=1`
   trees — so every FVS-native record with a missing crown, and every dead record, kept `crown_pct=0`.
3. **The crown ISORT keyed the wrong diameter.** `crown.f` ranks `IND`, which is always a sort of the
   *current* `DBH(I)`: cycling is `gradd.f:177-186` (`UPDATE` applies `DBH += DG/BRATIO`, *then*
   `RDPSRT`), and FVSjl's apply-loop likewise precedes the crown call without clearing
   `t.diam_growth`; LSTART is the read diameter. Eleven variants re-added a cycle of growth that had
   already been applied. BM's correct form is now the shared `crown_isort(s; lstart)`.

**Measured** on an OC control stand (S248112) carrying a missing crown on an `IORG=0` record: SDIAC
live 196.15 vs jl 196.14972; that record's crown 0 (never dubbed) → 85 (with SDIAC=0) → **76, live's
exact value**; `.sum` 1995 TPA live 511, jl 485 → 512.

On the absolute WRD fixtures **PN (control+rd), SO (control+rd) and EC (rd) all reach BIT-EXACT**,
joining BM at the oracle floor; NC's BA tolerance goes 7→3, CI's 5→2, UT's 6→3, while TT and UT each
gain one TPA count. FIA ledgers: BM 400 stands unchanged at 398/400 bit-exact (same two diverging
stands, zero signature changes), PN 5 improved / 1 worse, SO **23 → 25/40 bit-exact**.

Two lessons worth keeping:

* A fixture that exercises a mechanism is not the same as a fixture that *discriminates* on it. The
  SDIAC bug was invisible on all 22 WRD fixture rows because every missing-crown record in those
  `.tre` files is a **dead** record, which takes the `DO 79` DUBSCR path and never reads RELSDI.
* `test/harness/fia/extract_sample.jl` now stratifies on an inventory **tree-count class**
  (0 / 1-9 / 10-49 / 50+), balanced rather than proportional. Without it an EM draw came out 37/40
  bare — nonstocked conditions FVS fills by AUTOES — and the committed 12-stand EM *tiered fixture*
  came out 12/12 bare, so both instruments measured only the establishment path and were blind to the
  growth and crown models they exist to certify. `TREE_CLASS_STRATA=0` reproduces the old draw.

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
