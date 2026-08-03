# FVS Extensions Rollout Plan — FFE / ECON / Mistletoe / Climate across the 12 ported variants

User directive (after CI done, BC re-scoped): "port all of the extensions to those two and the previous
ones" → chosen sequencing "Extensions first (11 done + CI), defer BC's full BEC port." This plan drives
the FFE (+ ECON/mistletoe/climate) rollout across the 12 variants that have a validated growth core:
**sn, ne, cs, ls (eastern) + cr, ie, kt, em, bm, tt, ut, ci (western)**. Validate each per doctrine:
bit-exact-or-cornered vs the LIVE relinked binary, one chunk/variant at a time, instrument don't infer.

## Current extension coverage (measured, src/engine/fire + econ.jl + variants/*)
| Extension | Have | Missing |
|---|---|---|
| FFE (fire/fuels/carbon) | **CR** (western template, 23 dispatch pts) + **SN/CS/LS/NE** (eastern) | **ie, kt, em, bm, ut, tt, ci** (7 western) |
| ECON | econ.jl (variant-general: log-graded HRVRVN, discount) — validated SN/eastern | confirm western wiring (should be ~drop-in; ECON is variant-agnostic) |
| Dwarf mistletoe | CR (variants/centralrockies/dwarf_mistletoe_model.jl); IE MISTOE done (per memory) | per-variant DMR tables for other western (situational) |
| Climate-FVS | none in tree | all (lowest priority — climate is an optional keyword layer) |

## FFE rollout — the dominant phase. Order set by measured table reuse (bin/FVS*_buildDir/fm*.f diffs):
- **fmmois.f (moisture): ALL western identical (ie=kt=em=bm=ci, 0 diff).** Port the western moisture table ONCE.
- **fmcfmd.f (cover-type→fuel-model, the big table ~1250 lines): ie/kt = 0 (IDENTICAL).** em=205, ci=299,
  bm=622 differ from ie. ⇒ **IE port auto-covers KT** for the fuel-model half.
- **fmbrkt.f (bark→fire mortality): per-species DBH·B1, small tables.** ie=23 sp, kt=11 (=first 11 of ie).
  em/bm/ci differ by 36-83 lines (species-count). fmsnag.f: CR=IE identical (shared snag base).

### Sequence
1. **IE** (N-Rockies template) — full western FFE: fmmois(western-shared) + fmcfmd + fmbrkt + crown biomass +
   snag + fuel loading + dispatch `s.variant isa InlandEmpire`. Validate vs live FVSie SIMFIRE on iet01.
2. **KT** — nearly free after IE: fmcfmd identical, fmbrkt = IE[1:11], fmmois shared. Just add dispatch + 11-sp table.
3. **EM, CI, BM, UT, TT** — each: own fmcfmd (fuel-model), own fmbrkt, shared fmmois, crown-biomass species map.
4. Then ECON confirm-western + mistletoe (per-variant DMR, situational) + climate (optional, last).

### Per-variant FFE data needed (the extension points, mirrored from CR in src/engine/fire/*.jl)
- fire_effects.jl: bark B1 (fmbrkt) → mortality.  ✔ IE+KT bark DONE (this commit).
- fuel_moisture.jl: `_FM_MOIS_<v>` (fmmois) — western-shared, port once.
- fuel_model.jl: `_FMD_XPTS_<v>` breakpoints + fmcfmd cover-type→fuel-model map + fuel loading (FUINIE/FUINII).
- fmcba.jl: cover-type default fuel model, live/dead fuel loading, large-fuel photo-series species default.
- crown_biomass.jl: species→crown-biomass-eqn map (ISPMAP) + bratio (already have per-variant bark_ratio).
- fuel_decay.jl: `_FM_DKR_<v>` decay rates.  snag.jl: snag fall/decay (fmsnag shared with CR).
- fmburn.jl: conifer-species mortality gate + flame/byram kill.

## Progress log
- **2026-08-03** Assessed architecture; set order (IE→KT free, western-shared moisture, then EM/CI/BM/UT/TT).
  Ported IE+KT fire bark-thickness (`_IE_FM_BARK_B1`, ie/fmbrkt.f 23 sp; KT reuses [1:11]) into
  fire_effects.jl with dispatch. First FFE increment for the western N-Rockies cluster. Next: IE fmmois
  (western-shared) + fmcfmd fuel-model + iet01 SIMFIRE .sum validation vs live FVSie.
