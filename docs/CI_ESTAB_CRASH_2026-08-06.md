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

## ★★ CLUSTER-WIDE ROOT (deterministic trace, supersedes the CI-only framing)
The crash is NOT CI-specific. establishment.jl:132-136 sets `bc = nothing` for ALL western variants (NE/CS/LS/CR/IE/
TT/EM/BM/UT/CI — "western variants use a fixed/XMIN base, not the SN ht-curve"). The establishment PHASE-2 height
dispatch (establishment.jl:258-266) has explicit branches ONLY for: NE (ne_htcalc_height), CS (cs_htcalc_height),
LS (ls_htcalc_height), EM (em_essubh_hht). Everything ELSE → `htcalc_height(bc, ...)` (line 266) = the SN Chapman-
Richards curve, which indexes bc[1][sp] ⇒ getindex(::Nothing) CRASH when bc=nothing. ⇒ CR/IE/TT/BM/UT/CI ALL crash
on any PLANT/ESTAB/NATURAL stand (the phase-2 established-tree height). This is the SINGLE ROOT of the observed
bmt01 / utt01 / cit01 ESTAB crashes noted throughout this session — one dispatch gap, not per-variant bugs.
(EM was the only western with its establishment height wired — em_essubh_hht — likely from the #137-adjacent work.)
FIX (real, per-variant — NOT a guard: a placeholder height would silently produce WRONG establishment heights, worse
than crashing): wire each western variant's establishment/planted BASE HEIGHT branch in establishment.jl:258, mirroring
em_essubh_hht (EM essubh.f) — port ci/bm/ut/tt/kt/ie's essubh-equivalent (or their "fixed/XMIN base" per the comment).
CR is "COMPLETE" but would ALSO crash here on a PLANT/ESTAB .key — its DB-sweep validation used no-PLANT DB stands.
⇒ this is a cluster-wide ESTAB-height gap; only EM (this session) + NE/CS/LS (eastern) are wired. Reproduction:
any western .key with PLANT/ESTAB (cit01/bmt01/utt01) + run_keyfile. Handoff — NO fix-hypothesis on the exact heights
(fatigue lesson); the DISPATCH GAP is the solid measured root.

## EMPIRICAL CONFIRMATION (ran bmt01/utt01 through jl)
- BM bmt01.key: CRASHES at the SAME getindex(::Nothing) (htcalc_height bc=nothing) ⇒ #154 dispatch gap CONFIRMED for BM.
- UT utt01.key: crashes EARLIER at `KeyError: :essprt_fsp` (a SEPARATE UT establishment bug — stump-sprout ESSPRT
  coefficients unwired; noted earlier this session) ⇒ UT has an ADDITIONAL crash that precedes the htcalc_height path.
So: htcalc_height dispatch gap (#154) = CONFIRMED BM+CI, INFERRED CR/IE/TT (same ELSE→htcalc_height(bc=nothing) path,
deterministic); UT would also hit it but crashes first at :essprt_fsp. ⇒ #154 fix (wire western establishment-height
branches) closes CR/IE/TT/BM/CI; UT needs BOTH #154 + the :essprt_fsp stump-sprout wiring. Two distinct western-estab
bugs, both blocking PLANT/ESTAB projection.

## MEASURED FIX DATA — CI establishment height is essubh-driven (~1.5-1.9'), NOT a fixed base
Instrumented live estb/estab.f PLANT ESSUBH (line 1023) on cit01 (FVSci_trc): CI planted base height HHT = 1.57 /
1.63 / 1.76 / 1.90 ... ' for IPNSPE=2 (per-plot variation, TIME=FINT=10). So CI's establishment height IS the essubh
model (site/BAA-dependent + stochastic per-plot, like EM's em_essubh_hht ~1.0-1.18'), NOT the "fixed/XMIN base" the
code comment guessed. ⇒ the #154 fix per western variant = PORT that variant's essubh (mirror em_essubh_hht: the
essubh.f species equations are variant-specific DATA + a per-plot ESRANN perturbation) and wire it into the
establishment.jl:258 dispatch. This is a real per-variant port (CI/BM/UT/TT/KT/IE), NOT a one-line guard. Measured
target for CI sp2: ~1.5-1.9'. (UT additionally blocked by :essprt_fsp before reaching here.) Fresh-session chunk.
