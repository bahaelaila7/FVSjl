# CI cit01 ESTAB/PLANT path CRASHES in jl (pre-existing) — measurement-only finding

## Measured (cross-variant smoke test, run_keyfile)
- TT ttt01.key: RUNS clean (cyc0 1990 TPA=536 BA=77, 50 rows).
- CI cit01.key: CRASHES — `MethodError: no method matching getindex(::Nothing, ::Int64)` at
  src/variants/southern/height_growth.jl:26 (_htcalc_coef) ← htcalc_height(bc=Nothing,...) ← establish!
  (src/engine/simulate.jl:562) ← grow_cycle!.
- cit01.key has `ESTAB 1992` + `PLANT 1992 sp2 400` + `PLANT 1992 sp10 400`. So the crash is the CI ESTABLISHMENT
  path: establish! calls the shared htcalc_height with bc=Nothing (CI's height-calc coefficient struct not wired for
  the establishment height calc).

## Confirmed PRE-EXISTING (not a session regression)
Reverted my 4 changed code files (bluemountains/site_index.jl, bluemountains/crown.jl, easternmontana/regent.jl,
engine/simulate.jl) to 73ea2b3 (parent of my first code commit 1fbdd92) ⇒ cit01 STILL crashes identically. My
session's changes are BM/EM-gated and do not touch establish!/htcalc_height ⇒ EXONERATED. This is a pre-existing
jl bug in the CI establishment height-calc, same CLASS as the bmt01/utt01 ESTAB-path crashes noted earlier this
session (getindex(::Nothing) in the estab path) and the EM #137 essubh crash (fixed).

## Significance (reconciles the goal-file)
The goal-file's "CI cit01 1990 bit-exact" = only the cyc0/1990 INVENTORY row (produced before the 1992 ESTAB fires
in the first growth cycle). The "CI 25-stand real-FIA slice → 0 jl crashes" used FVS_STANDINIT DB stands, which carry
NO PLANT/NATURAL keyword ⇒ they never enter this establish! crash path. So the CI ESTAB/PLANT projection has been
UNVALIDATED (it crashes), while the DB growth-only path is fine. cit01's MULTI-cycle .sum (2000+) cannot be produced
in jl.

## Handoff (NO fix-hypothesis — per this session's fatigue lesson, 5 refuted roots)
NEXT (fresh): trace why the CI establishment path passes bc=Nothing to htcalc_height — is CI's height-calc coefficient
struct (the `bc` arg) unwired for establishment (vs the growth path where CI height works), or does CI need a variant
branch in the shared establish!/htcalc_height like the EM essubh fix (#137)? Wire/guard it so CI ESTAB/PLANT stands
run. Reproduction: /workspace/.ciwork/cit01.key + run_keyfile(variant=CentralIdaho()). Same class: bmt01/utt01 ESTAB.
