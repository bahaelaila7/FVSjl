# FVSjl — port & validation status

_Last updated 2026-08-24. Branch `kt-variant-port` (mainline; `master` tracks the
validated state). Hard gate: `test/integration/test_multicycle.jl` = **339 pass / 11
broken, byte-identical**._

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

**Western / non-eastern — full-population sweep IN PROGRESS** (`docs/WESTERN_FIA_FULLSCALE_2026-08-24.md`,
`data/fia_sweep_west.db`). 15 of the 20 non-eastern variants carry an FVS-ready FIA
population (688,903 stands); 5 have zero FIA (KT/BC/OC/OP/ON — Canada/no-FIA/ORGANON-BLM,
N/A). The sweep runs like the eastern one — every cycle's `.sum` vs freshly-relinked live
FVS — and has already surfaced and fixed **10 real bugs** the prior sampled validation
could not reach: two NWCMRT density-mortality omissions (NC, UT), six alpha-PV_CODE
habitat-decode bugs (PN/WC/CA/SO/NC/EC → correct SDIMAX), a 3-cause BM seedling
small-tree-growth bug, and a CA/SO forkod forest-index crash. On corrected code the BM
re-sweep reaches **99.6% bit-exact-or-cornered** (from a buggy 95.5%); NC's full re-sweep
is **88.9%** — most of its residual is a newly-identified under-mortality in the extreme
`>10k-TPA` mega-record seedling regime (a real-bug candidate one regime past the mature
NWCMRT fix, plus the FVS-flagged mega-record instability the eastern sweep cornered). Full
population coverage is a long-running background job; per-variant numbers accumulate as the
sweep (now crash-hardened) runs, converging toward the eastern standard.

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

## Known exceptions / not-yet-closed

- **ADDTREES** — external program, no in-tree source. Not ported.
- **PPE MXHRVP (multistand harvest scheduling)** — the one remaining un-wired PPE
  seam (mode-2 live beetle coupling is now done; harvest allocation across stands is not).
- **ON database-read path** — no real ON DB data + a characterized gcc-16×sqlite
  SIGSEGV; ON is validated via the inline path, not a live DB oracle.
- **Western full-population FIA sweep** — under way; not yet at the eastern
  exhaustiveness.
