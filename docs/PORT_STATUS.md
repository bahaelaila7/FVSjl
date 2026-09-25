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

## Crown initialisation, round two (2026-09-25, master `5930ffb5`)

Five more shared defects in the same path, all found by building instruments the shipped
fixtures could not provide.

**The instruments.** `/workspace/.postswap/wrd/genblank.jl` reruns each WRD control fixture with every
**LIVE** inventory crown BLANKED, and `gencrn.jl` layers the CRNMULT cases on top. This matters
because the shipped fixtures **cannot** exercise the LSTART dub on a live record — every blank-crown
record in them is a *dead* one, which takes the `DO 79` DUBSCR path and never reads RELSDI or PCT.
All twelve fixtures do carry two top-killed live records, so blanking turns each into a full test of
the dub: Weibull/PARM path, DUBSCR path, statement-55 top-kill, SDIAC, PCT.

**The defects.**

1. **CRNMULT was parsed, stored, and inert.** The record mapping was right (`initre.f:3416`, option
   96) but a blank UPPER DBH defaults to **99.0**, not "no limit" (`initre.f:3439`), and `ARRAY(6)` —
   the "DUB FLAG" that makes the multiplier scale the LSTART dub and then revert (`crown.f` stmt 60)
   — was not stored at all. BM already had the `CRNMLT`/`DLOW`/`DHI` structure at all six call sites
   but read the blkdat DEFAULTS; CR passed them as the literals `1.0/0.0/99.0`. The new shared
   `crn_mult_band` returns the **triple**, because `crown.f` gates the `CHG*CRNMLT` scalings on the
   tree's DBH being in band but gates the two `ICRI<10` floor bumps on `CRNMLT(ISPC) == 1.0` itself.
2. **`crown.f` statement 55** — a top-killed inventory record has its crown re-expressed on the
   NORMAL height at LSTART — was missing from nine variants' live loop. Now shared as `topkill_icri`.
   IE and AK needed a `dubbed` flag: label 58 enters *below* 55, so a DUBSCR-dubbed record must not
   be re-expressed.
3. **TT dubbed with RELSDI = 0.** It reaches CROWN through its own `tt_crown_init_lstart!`, which
   passed no `crown_sdi` — bypassing the fix the shared helper had received.
4. **KT and BC had no LSTART dub at all** (`kt/cratet.f:598`, `canada/bc/cratet.f`), the same defect
   OC had.
5. **LSTART PCT was ordered by the BACKDATED diameter.** `dense.f:244` accumulates the percentile
   over `IND`, which `cratet.f` sorted on the REAL `DBH`, while the weight `WK5 = D*D*PROB` uses the
   backdated diameter (`dense.f:184`). FVSjl derived *both* from the backdated diameters. PCT feeds
   the PCR crown model directly (`b13*P + b14*log(P)` in KT/IE/EM).

**Measured**, worst |jl − live| cell (TPA/BA) over 11 cycles:

| | blanked crowns, before → after | |
|---|---|---|
| KT | **31/51 → 4/17** | had no dub at all |
| BC | 8/7 → 3/4 | had no dub at all |
| EM | 12/7 → 7/4 | |
| IE | 1/1 → **0/1** | |
| TT | 13/7 → 11/7 | from SDIAC, not statement 55 |
| NC | 3/3 → 5/1 | BA tightens, TPA +2 |
| CR/BM/PN/SO | 0/0 throughout | |

With CRNMULT present: BM 0/0 on all four cases (validating the harness), CR **4/1, 17/4, 7/1 → 0/0**.
Per-record traces where it mattered: KT's PCT and EXPPCR now reproduce live exactly (45.693924 /
0.28728095 against live 45.694 / 0.287); TT's SDIAC is 202.939 against live's 202.94.

**Two lessons.**

* A crash this work introduced was caught only by the **tiered** suite: `_pctile!` reads `idx[1]`
  unguarded, so a stand with no live records *and* no dead ones threw a `BoundsError` — BM stand
  647500316126144, all 11 regimes. The WRD fixtures and the blanked A/B both have trees, so neither
  could see it. Different instruments fail differently; keep all of them.
* `FVSkt_clean` **SIGFPEs on a valid CRNMULT card**. The `ICRI<5` floor is gated on `CRNMLT == 1`, so
  an active multiplier disables FVS's own guard, the crown reaches 0, and the next cycle's
  `PDIFPY = CHG/REAL(ICR(I))/FINT` divides by it. KT cannot be A/B'd with CRNMULT active — an oracle
  limitation, not a port gap.

**Still open here:** CRNMULT remains inert in NC/CI/EC/PN/SO/EM/TT/UT/WS/IE — all measured, and
large (PN 83 BA, NC 57 BA). ON's statement 55 delegates to the shared eastern `_twigs_crown_update!`
that NE/CS/LS also use, so it needs an eastern top-kill fixture first. BC's sub-2cm route is not a
separate model: `canada/bc/dubscr.f` for V3 calls `CRNMD`, and `CRNMD` itself substitutes D = 2 cm
and a height re-derived from the species height–DBH curve — porting it needs BC's `AA`/`BB` and the
`LMHTDUB` flag, which may be metric-fitted.

## CRNMULT rollout, STDINFO habitat, site setup and western volume (2026-09-25, branch `crnmult-rest`)

Four linked findings. The last two were **hidden by the first two**: once each variant ran on the
habitat the user actually gave it, TPA and BA came to within ±1 of live, and that exactness exposed
volume columns that were still wrong on master.

**1. CRNMULT at every `crown.f` site in 13 more variants.** PN, SO, EC, WS and CA use the PN-family five
sites. NC, CI, TT, UT and WC use variant-specific patterns. EM and IE have six sites, including the
NIVAR PCR path, with a strict upper band (`.LT. DHI`); the westside variants use `.LE.`. AK is also
done. Faithfulness fixes found along the way:
- NC's missing CRMAX cap, and a redwood floor of 5 that has no Fortran basis.
- CI's missing `ICRI<10` CRMAX bump.
- TT and UT statement 59's `[10,95]` bounds on the main path.

Measured with `gencrn.jl` on the blanked-crown fixtures, as worst dTPA/dBA before → after:

| Variant | Before | After | Note |
|---|---|---|---|
| PN | 16/83 | 0/0 | |
| SO, EC, WC | — | 0/0 | all four cases |
| CI | 20/27 | 3/1 | |
| NC | 150/57 | 10/8 | |
| UT | 18/6 | 8/2 | |
| IE | 6/29 | 1/0 | |
| AK | 49/12 | 24/6 | effect exact through cycle 2 |
| TT | — | — | keyword effect exact for 4 cycles |

KT can't be A/B'd: `FVSkt_clean` SIGFPEs on a CRNMULT card.

**2. STDINFO field 2 (habitat) was dropped for WC, PN, SO, CA, NC and UT.** `kw_stdinfo!` filed it as a
southern ecological unit, so every keyword-path stand in these variants ran on the default plant
association, whatever the user gave. The fix follows `habtyp.f`:
- WC and PN are pure R6, decoded by HBDECD against PCOML.
- SO, CA and NC split on the post-FORKOD KODFOR region. Initre runs `CALL FORKOD` before HABTYP.
  - R6 decodes with HBDECD, plus the NR5-offset sequence fallback (SO habitat 460 → PCOML(54)).
  - R5 decodes to 0, because no R5 site path reads it.
- UT stores the raw KODTYP, which `ut_habtyp` CRDECDs.

The WC shipped fixture went from −40 BA to bit-exact.

**3. Site setup did not follow `sitset.f`.**
- **SO and NC:** ECOCLS never picked the site species or seeded SITEAR on R6. Every R6 stand ran as
  PP@70 (SO) or DF@90 (NC), right only for the default PA. Live 601/CPS311 uses PP SI 85; live
  611/CPC511 uses PP SI 52.
- **SO R5:** used the R6 default instead of WF/50 with C5 SDIDEF, and it overwrote keyworded SITEAR.
- **CA:** had no Region-5 branch at all (R5ADJ fan, R5SDI). R5 stands keep Zeide (`ca/grinit` sets
  LZEIDE true; only the R6 sitset resets it).
- **`nc/forkod`:** now covers all 11 JFOR codes, the reservation codes and the mapping correction.
  Trinity, Mendocino, Los Padres, Simpson and every reservation stand used to run as Klamath.
- **PMSDIU:** stored as a fraction. CA divided it by 100 again, and UT defaulted it to 85.

Result on the shipped PN tree list across forests (`habtest.jl`): SO and CA TPA/BA are within ±1 on
every case (they were up to −23 and +91 BA). NC is at its known base residual (7 TPA / 1 BA, was
175/35).

**4. Western volume.** Every finding was measured per tree from `FVS_TreeList`; the tool is
`/tmp/claude-1000/voldiff.sh`.
- The oracle's `fvsvol.f` passes western merch/board tops **inside bark**, as TOPD·BARK and
  BFTOPD·BARK. `vols.f:150` takes BARK at the **start-of-cycle** DBH, which jl stashes as `vol_bark`.
  NC, CA and SO used the grown-DBH bark, and NC/CA R5 a fixed 6″ top: NC R5 MCuFt was −11% and BdFt
  −20%. The old NC audit put this down to "TPA normalisation"; that doesn't hold, because TCuFt is
  exact.
- `nvel_r5_vol` is a new shared Region-5 dispatcher (WO2W / DVEW / INGY FW2). SO and CA R5 forests
  used the R6 INGY/Behre tables; they now use `voleqdef.f R5_EQN`, which depends only on species and
  variant.
- The CFTOPK/BFTOPK broken-top trim was missing in SO and CA. Example: a DF broken at 49 ft, D 15.9,
  was 53.8 cuft in jl vs live 40.2.
- NC had no `init_merch_standards!` branch, so it ran on generic CSV specs. It now follows
  `nc/sitset.f`: top 4.5 on Siskiyou, 5.0 on the BLM forests, 6.0 on R5.
- NC 712 (Coos Bay BLM) now uses BLMVOL (it had been sent down the R5 path), with the `BLM712` form
  classes and DF B02 profile/taper 2.

Result: cycle 0 is per-tree bit-exact on SO 601/505, CA 610/505 and NC 611/712/518. SO and CA stay
exact through 2020, apart from one CA WF board value 10 bf off in 2020. NC's later-cycle differences
are only on trees whose DBH already differs (the known NC DG residual).

**5. The self-thinning line is latched (EM, NC, UT, TT).** In `{em,nc,ut,tt}/morts.f` the 55%–85% line
solved the first time (SLPMRT/CEPMRT, VARCOM) is **kept and reused in every later cycle**. It is reset only
when RMSQD==0, or when ICYC>1 and |T−TPAMRT|>1 (TPAMRT is the post-mortality TNEW). jl's shared southern
driver and ON already did this, but `_em_tn10_iter` and `_tt_tn10_iter` re-solved the line every cycle.

Measured: with DEBUG MORTS, live UT latches 7.72198 / −0.714817 in cycle 5 and reuses it in cycle 6. jl
killed a uniform 0.736× live's mortality on every record from 2050. That was the whole "UT base residual".

Worst TPA/BA on the shipped fixture, versus master:

| Variant | Master | Now |
|---|---|---|
| UT | 13 / 3 | 0 / 0 over 10 cycles (bit-exact) |
| NC | 7 / 1 | 1 / 1 |
| EM | 8 / 6 | 3 / 3 |
| TT | 20 / 7 | 10 / 3 |

The RMSQD test uses `stand_qmd`, because `p.qmd` is only populated for Ontario.

**6. WS.**
- `ws/ccfcal.f` CASE(9,10,12,14:17,19,20,25:27) is R5CRWD crown width² × 0.001803. jl had a 0.001 stub,
  so cycle-0 CCF was 58 vs live 65; it is now 65.
- DVE hardwood board-foot volume had been zeroed as "deferred"; it is now r5harv at BFTOPD·BARK.

**7. NC's first growth cycle used the wrong previous period (branch `nc-dg`).** In the first projection
cycle, `grincr.f` sets OLDFNT = FINT, the DG measurement period from grinit. AUTCOR uses it as the previous
period when it correlates this cycle's random DG error with the calibration residual. jl fell back to
`htg_period` (YR), which equals FINT in every variant except NC (nc/blkdat YR=5, nc/grinit FINT=10).

Measured with a separate instrumented FVSnc (`/workspace/.ncwork/FVSnc_oldrn`, which prints OLDRN and
SSIGMA/RHO; the oracle is untouched):
- The calibration OLDRN seeds were already bit-exact (26 of 27 trees; one WF is 10 ULP off).
- With `DGSTDEV 0` NC already matched live, so the whole residual was the random component.
- Live's cycle-1 CORR is 0.3906 = AUTCOR(5,10), giving RHO 0.4036. jl used AUTCOR(5,5) = 0.3196.

The fix is `dg_measure_period(v)` (grinit FINT: defaults to `htg_period`, Klamath overrides to 10).
Together with the latch in §5, every NC habtest case (611/30, 712/CWC221, 611/440, 8103, 518, 508) is
bit-exact over 10 cycles on TPA/BA/TCuFt/MCuFt. NC's WRD absolute row (control and root disease) is also
bit-exact and has been promoted to a passing `@test`.

**Oracle changes (user-approved 2026-09-25)**, recorded in `/workspace/ORACLE_SOURCE_AUDIT_2026-09-19.md` §6:
- Debug WRITEs removed from the BM, EM, IE, SN and CR buildDirs.
- The CR `varmrt.f` TEMSUM guard is now in `FVScr_clean`.
- Every shipped `.sum` is byte-identical before and after.

**Gate reconciliation.** The tiered fast tier moved in both directions:
- 38 allowlist entries improved and were tightened (25 IE, 11 EM, 2 SN).
- 10 already-OPEN EM entries grew by 1–16 cells: mistletoe CCF 35→36, TreeList BdFt 2111→2127, PtBAL
  9603→9614, plant MAI 32→33.

The EM latch was verified against live before accepting that: jl latches in the same cycle as live
(cycle 5) and holds 7.90189 / −0.76509 in every later cycle, matching live to ~1e-5. The residual
difference comes from EM's already-differing D10 (8.397 vs 8.347). So the spread is trajectory churn
inside open EM residuals, and the allowlist was widened with a note on each entry.

UT's WRD absolute row became bit-exact, so it was promoted from `@test_broken` to a passing test.

**Still open from this round:**
- WS growth residual: cycle 1 BA is 136 vs 141.
- EM and TT residuals: 3 BA.
- R5 WO2W board-foot: single-tree 10-bf Scribner steps on 1–2 trees per year (NC 518/508 BdFt 70; one CA
  WF in 2020). The standalone VOLINITNVB driver could not reproduce live's merch at the treelist D/H, so
  it needs the exact start-of-cycle volume inputs.
- CRNMULT for OP, OC, KT and BC.
- The SDICALC min-DBH filter on morts.f's T (`D < DBHZEIDE/DBHSTAGE`) is not applied in the EM/NC/UT/TT
  kernels. It is inert at the default 0.

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
