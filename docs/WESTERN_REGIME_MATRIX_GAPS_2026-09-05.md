# Western regime-matrix enrichment — coverage-breadth ranking + gap catalog (2026-09-05)

Coverage harness: `test/harness/fia/coverage_breadth.jl` + `.sh`. Ran each western variant on 2 stratified
stands × regimes {none,plant,thinbba,salvage,simfire} under `julia --code-coverage=user`; counted distinct
FVSjl src lines covered. (Extensions WRD/mistletoe/climate/econ/cover NOT yet added — base+existing regimes.)

## BREADTH RANKING (widest FVSjl code first; population = TIEBREAK only)
IE 5461 (10/10) · BM 4895 (10/10) · SO 4823* (8/10) · EC 4511 (10/10) · EM 4400 (10/10) ·
PN 3867* (7/10) · NC 3732* (8/10) · CA 3696* (5/10) · CR 3684 (10/10) · WC 3604* (7/10) ·
TT 3441 (10/10) · UT 3355 (10/10) · CI 3235 (10/10) · WS 2995* (6/10) · AK 2898* (6/10)
(* = coverage UNDERSTATED — some regimes crashed, see gaps below. cov-files ~113-124 all variants.)
⇒ **IE = widest-breadth, START HERE.** CR (pop #1 = 338645) is only breadth #9 ⇒ breadth-first ≠ pop-first (validated).

## GAP CATALOG (real deferred/unported items the coverage pass EXPOSED — "leave no item deferred")
1. **`KeyError: :estab_min_ht` under `plant`/ESTAB** — SO, CA, PN, WC, WS, AK crash ⇒ ESTABLISHMENT (PLANT) path
   UNWIRED/incomplete for these ~6 variants (IE/EM/BM/EC/TT/UT/CI/CR/NC have it). WIDEST gap.
2. **`KeyError: :essprt_fsp` under `thinbba`** — CA (+ likely others) ⇒ post-thin SPROUTING unwired (sprouting was
   IE/EM-only this session; x-variant sprouting gap confirmed).
3. **Crown-width follow-up chunks under `simfire`/FFE**: CA `cr_crownw` large-tree SPIE grp 21; WC `cr_crownw`
   small-tree TOTWT grp 10; PN `pn_cwcalc` eq 26305 (sp19); **NC `nc_cwcalc` eq 81802 (sp11) — MA/IC/BO/TO/OH/RW
   crown-width NOT ported (NEW NC gap, hardwoods)**.
4. **FFE fully UNPORTED**: WS (`ffe_fuel_live` empty, CR-FFE chunk pending) + AK (ORGANON species no FFE coeffs).
5. **Misc config-key gaps**: AK `:dkr_cls` under simfire.

## PLAN (next)
- Add extension regimes to regime_block (WRD RRTINV/RDATV, mistletoe MISTMULT, CLIMATE, ECON, COVER) → will
  surface MORE gaps (variants not wiring those extensions).
- INDEPENDENT-AGENT REVIEW of the whole rework.
- Sweep IE-first across the full matrix; fix every gap above (estab wiring, sprouting, crown-width chunks, FFE
  WS/AK, config keys) per "leave no item deferred".
