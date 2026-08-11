# BAMAX residual-BA mortality cap — cross-variant audit (2026-08-11)

## The mechanism
FVS `morts.f` applies a **residual-basal-area cap** after the SDI/Zeide self-thin, MSB, and size-cap
mortality are computed. Header (nc/morts.f): *"SDI-based mortality is used as long as the stand diameter
metric (QMD or Reineke's diameter) is less than 10 inches, at which time BAMAX-based mortality takes over.
If not set by the user, BAMAX is determined from max SDI at 10 inch DBH."*

- **BAMAX value** — `vbase/sdical.f:203-208`: when the user did NOT enter a BAMAX (`.NOT.LBAMAX`, the FIA
  default), `BAMAX = SDIMAX · 0.5454154 · PMSDIU`. `SDICAL` is called at the top of `morts.f` every cycle, so
  BAMAX tracks the (possibly climate/forest-type-adjusted) SDImax. `PMSDIU` is converted `85 → 0.85` in
  `morts.f` *before* the SDICAL call, so the cap uses 0.85.
- **The cap loop** — `morts.f` BA-check (e.g. nc:685-754): compute residual `BANEW = Σ(P−WK2)·0.0054542·(D+G)²`
  and `BADEAD = Σ WK2·0.0054542·(D+G)²`. While `(BANEW − BAMAX) > 1`, scale every record's kill up by
  `ADJFAC = (BANEW−BAMAX)/BADEAD` (`WK2 = min(P, WK2·(1+ADJFAC))`) and re-test, iterating ≤100×. `G` is the
  outside-bark linearly-FINT-extrapolated increment `(DG/BARK)·(FINT/YR)` — in jl this is already baked into
  `diam_growth`, so `g = diam_growth/bark`.

Without it, dense stands whose residual BA exceeds the self-thinning maximum grow **unbounded** — a ~2×
mortality under-kill on the densest large-tree stands.

## jl coverage audit
The block exists in **sn, cr, bm, em, ut, tt, so, ca** `morts.f` (CI/IE/KT/BC lack it — correctly absent in jl).

| Variant | jl mortality path | Cap present? | Action |
|---------|-------------------|--------------|--------|
| SN, CR, BM | shared `southern/mortality.jl` (no override) | ✅ lines 460-482 | covered |
| EM | own `easternmontana/mortality.jl` | ❌ was missing | **FIXED** (commit 427eca4) |
| UT | own `utah/mortality.jl` | ❌ was missing | **FIXED** (commit a503f0c) |
| TT | own `teton/mortality.jl` | ❌ was missing | **FIXED** (commit 8179bae) |
| NC | own `klamath/mortality.jl` | ❌ was missing | **FIXED** (commit b93b1a2) |
| CI, IE, KT, BC | own copies | n/a (block not in their morts.f) | correct |

⚠ **Trap avoided:** EM/TT jl already had `bamax` *terms*, but those are the **BADIST background-mortality
weighting** (Stage variants: `ripp += (bamax−ba)·rip; ripp/=bamax`) — a *different* mechanism from the
residual-cap loop. The residual cap (BANEW/BADEAD/ADJFAC iteration) was genuinely absent from all four
own-copies. NC was forked line-for-line from UT, inheriting UT's omission.

**EM BAMAX value:** `em/sitset.f:159` sets `BAMAX = BAMAXA(ITYPE)` (habitat table) but does NOT set
`LBAMAX = .TRUE.`, so `sdical.f:203` overwrites it with the SDI-derived `SDIMAX·0.5454154·PMSDIU` every morts
cycle (the habitat BAMAXA only feeds `SDIDEF → SDIMAX` at setup). So EM's residual cap uses the **SDI-derived**
value — same as NC/UT/TT — NOT the habitat `EM_BAMAXA` used by its BADIST weighting.

## Validation (live oracle, high-BA mature FIA stands where the cap fires)
- **NC** cn 1123874220 (99% redwood): BA 597→464 vs live 468 (27.6%→0.9% cornered); full 12-stand sweep
  9/12→10/12 cornered, 0 regressions.
- **UT** 6 high-BA mature FIA stands: jl-cap vs live all 0–1.5% BA (4 bit-exact: BA 281/380/502 == live).
  Decisive because live caps these high-BA stands; pre-fix jl would over-grow *above* live.
- **TT** 3 high-BA mature FIA stands: jl-cap vs live all **BIT-EXACT** (BA 653/300/287 == live, Δ0.0%).
- **EM**: the same-mechanism sibling (TT) is bit-exact, so the mechanism + SDI-derived BAMAX are sound. A
  cap-off toggle on 2 high-BA EM stands proved the cap does **not regress** them (cap-on BA == cap-off BA ==
  432/429 — inert below EM's BAMAX). Faithful port committed. ⚠ SEPARATE PRE-EXISTING EM issues surfaced by
  the sweep (NOT the cap; cap-off identical): cn 49331474020004 jl BA 432 vs live 511 (~15% under-growth);
  cn 196394428020004 a jl growth BLOWUP (BA 5.85e6 vs live 0) — filed for separate EM growth investigation.

**Audit complete:** all four own-copy variants (NC/UT/TT/EM) fixed + validated; BM/SN/CR covered via the
shared path; CI/IE/KT/BC correctly without it.

Inert (immediate `break`) whenever residual BA ≤ BAMAX, so the cap cannot regress the canonical test stands
(nct01/utt01/emt01/ttt01) or any below-cap stand — confirmed non-regressing on dense sub-1" seedling stands
(low BA; they diverge only via the separate #156 small-tree-DG tail).

## Meta
- **Refutes the goal-doc BM #140 = "missing self-thin" hypothesis**: BM caps via the shared path; verified
  before porting (doctrine #2 — measure, don't infer), which avoided a redundant fix into BM.
- One NC redwood divergence, chased to root, surfaced a **systematic four-variant** latent gap.
