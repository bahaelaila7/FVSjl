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

## WS height growth, REGENT and site; per-variant PSIGSQ (2026-09-26, branch `ws-htg`)

**1. WS `htgf.f` surrogate branches.** Four branches were stubs; all are ported and measured on a species-swap
fixture (the PN tree list with RW/GS/MC above and below HTMAX, and GB large and small) against `FVSws_g16`:
- **GB** uses the Alexander curve. BAU and AGERNG are identically 0 in WS, and the ZZRAN draws run
  species-major.
- **MC** uses Curtis potential plus Hoerl/Chapman-Richards modifiers, with no SCALE on the normal path and
  HTMAX 20.
- **The CA surrogates** use FINDAG, `ws_htcalc` and the Ritchie–Hann XMOD.
- **RW/GS** use the Castle LTHTG form.
- `ws/findag.f` AGEMAX/HTMAX are per species.

**2. WS REGENT (`ws/regent.f`).**
- REGYR is 10 only for GB and MC; every other species uses 5.
- BKPT is 99 for GB, 7 for RW/GS and 3 for the rest. At or above it the large-tree DG is left untouched.
- Below 4.5 ft, DBH is set directly and DG is 0.
- The CA/SO species take the HTDBH override when `!LHTDRG` or IABFLG = 1.
- GB DG was the dgf large-tree value (0.549 vs live 0.054). It is now the subtraction form, capped at
  2·SCALE and DDS-scaled.

**3. WS point PRD.** `ws_point_prd` read a Density field that does not exist, so PRD was 0 everywhere. That
affected the RW/GS DG term (−0.42078·PRD) and every WS DUBSCR crown. It now runs `point_zeide!` per point.
Live RW PRD is 0.0665, and DGLT matches.

**4. WS site and volume.**
- `ws/forkod.f`: an unmatched forest code keeps grinit IFOR 6 (517), the reservation codes are mapped, and
  TLAT is set only for 5-digit codes (grinit default 39).
- CFTOPK/BFTOPK broken-top trims are applied.
- The R5 WO2W board pass uses MERLEN (`nc_merlen`), in WS and in NC's `nc_wo2w_vol`.
- SO Fremont (IFOR 2/3) uses its own volume-equation table.

**5. Per-variant PSIGSQ (`dgdriv.f`).** Every variant fell through to SN's 0.089827273. `src/variants/psigsq.jl`
now carries each variant's own DATA:
- SN and AK use 0.089827273.
- CS, LS, NE, ON, CA, WC, PN, OC and OP use 0.0898.
- SO, WS and EC use per-species tables.

PSIGSQ enters the empirical-Bayes COR shrinkage. On SO WF, WC was 0.5701 vs live 0.5783.

Result: WS 516/505, SO ×5, CA ×4, UT ×3 and NC ×6 habtest cases are exact. On the species-swap fixture,
`DGSTDEV 0` is byte-identical over 11 rows, and with the random component on, one year differs by 1 bdft.

**Gate reconciliation.** Three PPE MXHRVP assertions are jl self-snapshots of EC growth. The EC PSIGSQ fix
moved each one toward the FVSppe oracle:

| Assertion | Before | After | Oracle |
|---|---|---|---|
| 2000 selected resource | 326.6 | 327.3 | 330.48 |
| 2010 selected resource | 432.3 | 433.7 | 439.90 |
| 2000 HVPART | 0.378 | 0.375 | 0.362 |

The same stand run alone against live FVSec agrees: 2010 BA went from 144 to 145 (live 147), and 2020 from
180 to 181 (live 184). The snapshots were re-pinned. No other test changed by name.

Three `@test_broken` rows became passing tests: the EC and WS WRD absolute rows (control and root disease
TPA/BA equal live in every cycle; WS was TPA +28 / BA −10) and the EC `ect01` cycle-0 ACCRETION cell.

**Found and still open:**
- EC `ec_cwcalc` has no crown-width equation for species 6 and later, so TREELIST crashes on the PPE stand.
- On that stand, EC 1990 TCuFt is 1640 vs live 1602 while TPA, BA and SDI are exact.
- WS default-branch tripled-copy HTG.

## EM REGENT in the Fortran's shape, em/bratio.f, cycle-1 WK1, LL crown (2026-09-26, branch `em-vol`)

**1. REGENT (`em/regent.f`).** jl ran the five EM small-tree sub-models (EMVAR SMHTGF/SMDGF, NIVAR LL, TTVAR LM,
CRVAR, UTVAR) as separate passes. `small_tree_growth!` is now one subcycle loop (`DO 17 J / DO 16 ISPC / IND1`) and
one DO-30 assembly with the 918 tripling loop-back. That makes three things match live:
- the order of the random draws;
- the running RDNEXT/BANEXT(J+1) density feedback;
- the variables Fortran carries from tree to tree: BARK (the CR/UT DGK), H1 (the NIVAR RELH), D (across tripled
  copies) and HTGR (CR/UT).

Along the way:
- EMVAR ZRAND is persistent: drawn at −999, reset when the increment floors at 0.1.
- TPCCF is PCCF·PPCCF.
- SMDGF runs at every height, floored at DIAM; DKK = SMDGF(HT) even below 4.5 ft.
- The calibration NPER comes from IFINTH.
- The NIVAR ZZRAN bound is [−1.5, 1] on HTGR1·e^(Z·HSIGMA).
- The DO-6 CCF carries P.
- XMAX/XMIN for OH (species 19) are 2.0/0.5.
- **ESTAB** (`em/esgent.f`) runs the same routine in LESTB mode:
  - the DO-13 crown draw in storage order;
  - FINT−5 years of subcycling from TEMBA/TEMCCF/TEMAHT;
  - XWT = 0, the ESTAB diameter, no tripling or DUBSCR;
  - then HTG·WK4 and the HHTMAX cap.

**2. `em/bratio.f`.** EM bark had been mapped onto the generic b + a/d form, which drops TEMD = max(D,1) (and the
≤ 19 cap on the RM curve). A sub-1" GA/CW/BA/PW/NC/OH therefore got 0.80 instead of 0.806. `em_bratio` is now wired
into `variant_bratio` and every EM call site.

**3. Cycle-1 mortality WK1.** jl passed the raw input DG, so every unmeasured tree had WK1 = 0 and the added-species
Hamilton G fell to the DGT floor: an LM of 1.5" killed 14.3 TPA vs live's 5.9. `em_cycle0_wk1!` now follows the
`dgdriv.f` DO-220 precedence:
- a measured DG is kept;
- HT ≤ 4.5 gives 0;
- otherwise √(D² + e^(WK2+OLDRN)·SCALE) − D, where WK2/WK3 come from the first calibration DGF with the current
  RMSQD.

**4. LL crown (`em/crown.f`).** OBA/RDM1 were never threaded (OBA == BA, so DCRCON == XCRCON), and DCR used the
current PCT instead of OLDPCT. Both are fixed, along with the backdated-D < 3 skip. LL crowns now step
55 → 53 → 51 → 49 as in live; jl had held them at 55.

**5. HTGF.** HTGMULT (XHT) is applied, and each LL tripled copy gets its own large-tree HTG (`em_triple_htg!`).

**Measured against FVSem_g16:**
- habtest `em` 102/250 and 114/0 are bit-exact on TPA, BA, TCuFt, MCuFt and BdFt over 10 cycles (before: up to
  3 TPA / 3 BA / 853 BdFt off).
- A species-swap fixture with seedlings of every sub-model is exact per tree in cycle 1 with DGSTDEV 0, with the
  random component on, and with tripling. It is committed as `test/unit/test_em_regent_dk.jl`.
- The EM WRD absolute row is now a passing bit-exact test.
- The EM tiered fast tier dropped from 177,864 to 130,863 mismatch cells; the allowlist was redrafted from
  measurement (751 → 711 entries).

**Still open:**
- **Climate tally.** The climate TALLY is now one-directional (TPA over = 8, BA over = 9). The gap was already
  there: master was +277 TPA on stand 231908428020004, where live kills the stand in 2023. The REGENT fixes
  removed the stands that used to balance it. Tracked as task #247.
- **Two small residuals:**
  - The aspen calibration DGFASP uses RMSQD 2.998 in live vs the current QMD 2.853 in jl (source not yet found).
  - An LM backdated-PCT tie breaks in a different order (0.013 TPA).
- **REGCAL.** The CR/UT EDH uses a static PCTRED, which is 0 at LSTART (#244).

## EC crown width, volume, REGENT and ESGENT; western fuel moisture; CR NORMHT (2026-09-26, branch `ec-cw`)

Driven by the PPE MXHRVP landscape (three copies of Mt Hood stand S248112, forest 606), whose cycle-1/2 CREDIT
was off the FVSppe oracle because of EC growth, not MXHRVP.

**1. Crown width (`ec/cwcalc.f`).** `ec_cwcalc` hard-coded forest 608 for six species and errored on the other
26, so FVS_TreeList crashed on any EC stand with WH, GF or RC. It is now ECMAP → the national dispatcher with the
shared Region-6 forest BF, with KODFOR passed from both call sites (FMCBA PERCOV, TreeList CrWidth). Codes 02206,
63102 and 81505 were added to the national library. S248112 1990 CrWidth: 29/29 trees equal to live.

**2. Volume.** EC is an eastside variant in `voleqdef.f` R6_EQN, so it now reuses the BM forest table
(`_bm_r6_eqn`, DISTNUM 0). jl had known only Okanogan/Wenatchee and sent every other forest to Behre. On Mt Hood,
live runs DF on westside F05FW2W202 and the rest on I11–I13 INGY. The broken-top CFTOPK/BFTOPK trim now runs for
every equation, not just INGY. S248112 cycle 0 TCuFt/MCuFt/BdFt 1640/1103/5572 → 1602/1064/5456 = live (#243).

**3. REGENT (`ec/regent.f`, `ec/smhtgf.f`).**
- SMHTGF reads the tree's own SITEAR, unclamped. jl had clamped it to the site species' range, so an ES at SI 148
  grew 7.9 ft instead of 17.4.
- DK/DKK are formed per species. The HTDBH inventory form applies only when `.NOT.LHTDRG` or IABFLG = 1.
- Tripling gives each copy its own ZZRAN, HTG blend and small-tree DBH increment.
- REGHMULT and REGDMULT are applied.
- S248112, 3 cycles: per-tree exact in all four setups (DGSTDEV 0 and random, with and without tripling).
- PPE MXHRVP end-to-end now pins cycles 1–2 to the oracle (330.4816 / 439.8994).
- The PPE composite-materialization test (`stand_thin.key`) is now pinned to **current live FVSec_g16**
  volumes. Those volumes equal jl's on every row (1602/1064/5456 … 1780/1651/9038). FVSppe's composite volumes
  (1624/1102/5567 at 1990) come from the historical source's volume equations. TPA is still pinned to FVSppe,
  and is exact there.

**4. ESGENT (`ec/esgent.f`).** EC was missing from the birth-cycle ESGENT dispatch, so planted records sat at
their ESSUBH height for the whole establishment cycle (ect01 PLANT stand: 2002 BA 0 vs live 22). `ec_esgent!` is
REGENT(LESTB) for the new records, with the crown and height draws interleaved per record, then the HHTMAX cap.
The ect01 PLANT stand is per-tree exact over 10 cycles (#245).

**5. FFE fuel moisture.** SIMFIRE/POTFIRE moisture presets for EC and SO now use the IE table, and CA and WS use
the NC table. They had fallen through to the SN table. Each was checked value by value against the variant's
`fmmois.f`.

**6. CR broken-top volume (`cr/vols.f` TKILL).** H = NORMHT/100 for top-killed trees, as for EM/KT. The CR
cycle-0 .sum on the shipped PN list is now 1611/1381/3490, equal to live (#241).

**Still open:**
- **EC FFE fire behaviour (#246).** On the ect01 FFE stand, flame length is 3.92 vs 5.37 ft because the
  fuel-model weights differ (FM10 55% vs 46%), and 2003 SIMFIRE under-kills. jl does not yet emit FVS_Fuels under
  DATABASE FUELSOUT.
- **BC/ON/AK NORMHT (#248).**

## LSTART small-tree height calibration (REGCAL) for EM and IE (2026-09-26, branch `em-regcal`)

The small-tree height calibration in REGENT (label 40, called from CRATET) sets each species' HCOR from seedlings
and saplings that carry a measured HTG. Live's `DEBUG REGENT` prints per-species "SUMS FOR SPECIES n: SNP SNX
SNY", which gives a direct oracle; CORNEW = SNY/SNX.

**EM.**
- **RHCON.** Only NIVAR (LL) carries the NI constant REGCH + 1.0667 + RHHAB. Every other sub-model has
  RHCON = 1.0, or RCOR2 under READCORR (`em/regent.f` REGCON). jl gave all 19 species the NI constant (0.47 on the
  S248112 fixture), so CW's CORNEW came out 2.50 instead of 1.004 and CW height growth was about double.
- **Stand values.** REGCAL runs inside CRATET right after the `cratet.f:182` backdating DENSE, with no DENSE in
  between. It therefore reads that DENSE's BA, RELDEN, AVH, point PCCF and PCT, plus RELDM1 interpolated to the
  FINTH-year start (dense.f:259). jl read the current live-only values (PCTRED 0.459 vs live 0.5379).
  `crown_init_lstart_dead_inclusive!` now snapshots them (`Calibration.cratet_*`).
- **NTYR = IFINTH.** The new `Control.growth_ifinth` is 5 by default and is set only by the DB `HTG_MEASURE`
  column (dbsstandin.f:711). The GROWTH keyword never sets it: its assignment is commented out at initre.f:834. For
  IFINTH > 5 the DO-49 subcycle density projection runs; it matches live on three FIA stands with
  HTG_MEASURE = 10 (NPER 2).
- **Tie order (shared with BM).** The backdated PCT, and EM's DG-calibration PCT, now use the `cratet.f:153`
  IND (`IND=IND1; RDPSRT(.FALSE.)`, dead included). Before, jl used the identity `.TRUE.` re-sort or a stable
  `sortperm`, which permuted PCT inside equal-DBH groups. That reached NIVAR BAL and the DGF WK2, hence the DO-220
  WK1 dub and the LM/LL Hamilton G.
- **CCF at D = 0 (shared with IE/UT).** The IMC = 9 dead that the backdating DENSE zeroes take `ccfcal.f`'s
  D ≤ 0.1 branch, giving 0.001·P. jl gave 0.
- **LESTB crown draw.** It sets ICR only; PCT stays estab.f's 0.

**IE.**
- **NTYR = IFINTH.** jl used FINT (10) with REGYR 5, which ran 2 subcycles where live runs 1. Every NIVAR SNX was
  about 2× live, so CORNEW was halved.
- **Stand values.** The same `cratet.f:218` snapshot is now used (RELDEN 204.36 → 205.42; RELDM1 had been 0).
- **CR/UT arms.** These were never ported. They are now: CRVAR and PI/UJ take the UT form on the unclamped SJ,
  the aspen group gets Sheppard, and the result is ×0.5 using the last subcycle's value.
- **Smaller fixes.** IHTG < 2 gates the backdate, and READCORR is honoured.

**Measured:**
- New test `test/unit/test_regcal_em_ie.jl`, with fixtures under `test/fixtures/{easternmontana,inlandempire}/regcal`.
  Each fixture is the PN inventory plus six HTG seedlings for each of 12 species, including 12-way DBH ties and
  two dead records.
- All 15 EM and all 14 IE per-species sums equal live, across every sub-model.
- EM runs 3 cycles per-tree exact (0/99 trees off; em-vol was 89–99/99), and every .sum field matches live.
  That includes the trailing size/stocking class, after the STKVAL fix below.
- IE per-tree mismatches at 2000 fall from 91/99 on master to 29/99. The remainder is IE's growth-side TTVAR,
  which is not yet in Fortran shape (#250).
- Tiered fast tier: EM 130,863 → 130,769 cells (allowlist redrafted); IE and BM unchanged.

**STKVAL (all western variants, #251).** `stkval.f:325-333` redefines TAB3 at run time for every variant except
CS/LS/NE/SN/ON: FIA 299 ("west other softwood") → stocking equation 8, and 998/999 → 26. The per-variant CSVs
carry only the eastern DATA values (299 → 0, 998/999 → 25). A western OS record was therefore stocked on
equation 25: on the EM fixture, OS SS was 9.57 vs live 4.49, so TOTSTK was 102.85 vs 97.77 and the stocking class
was 1 vs 2. This feeds the FORTYP group array and the .sum size/stocking class; western growth does not read the
forest type.

**Still open:**
- **#249** — KT and CI have no REGCAL port at all, and TT needs an audit (in progress on branch
  `regcal-kt-ci`).
- **#250** — the IE REGENT growth rewrite.

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
