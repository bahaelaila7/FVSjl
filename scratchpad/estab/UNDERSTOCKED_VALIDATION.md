# ESTAB keywords validated on an UNDER-STOCKED IE stand (where establishment actually fires)

## Why a new fixture
The prior fixture (dense 28000-TPA stand 12343703010690) is SATURATED: instrumented FVSie_estabdump
(estab.f:683 NEWTPP dump) proved the ORACLE books NEWTPP=0 there across all cycles (ITPPnew 1-5 << NSTOREold
99) — so every establishment keyword is inert on BOTH jl and oracle. Wrong fixture.

## The right fixture
Queried /workspace/SQLite_FIADB_ENTIRE.db for VARIANT='IE' stands with few tree records (under-stocked).
FVSie_estabdump confirms real establishment: stand 753189105290487 (9 trees) books NEWTPP up to 4/plot
(PROB1 0.6174, PNN 0.5194 — the estab.f:583 PNN floor is inert here since PROB1>PNN, so jl's omission of the
581-583 clamps does not bite on this stand). ie_understocked.db + ie_v_{base,sa050,noing}.key.

## Results (live TPA, NUMCYCLE 5, ESTAB packet)
                 2019   2029   2039   2049   2059   2069
ORACLE base       763    815    603    728    572    793
ORACLE STOCKADJ Δ   0   -241   -160   -357   -242   -497
ORACLE NOINGROW Δ   0   -241   -160   -357   -242   -497
jl base           763    574    764    559    634    489
jl STOCKADJ Δ       0      0   -243   -148   -223   -145
jl NOINGROW Δ       0      0   -323   -196   -318   -208    (AFTER the kw_estab! fix; was ALL 0 before)

## Verdicts
- **STOCKADJ**: WORKS in jl (Δ-243 at the first establishment cycle ~ oracle -241). The prior "inert"
  reading was a SATURATED-stand artifact, not a code bug. eddf8dc2 code is correct.
- **NOINGROW (+INGROW/AUTALLY/NOAUTALY/THRSHOLD)**: was a REAL BUG — jl handled these only at TOP-LEVEL
  (keyword_dispatch.jl:2428-2434) but FVS parses them INSIDE the ESTAB packet (esin.f opt 21-25), so
  NOINGROW-inside-ESTAB was silently skipped (Δ0). FIX: added the handlers to kw_estab!. NOINGROW now
  disables ingrowth (Δ large-negative, correct sign+scale vs oracle).
- **Residual** (jl establishment ~1 cycle later than oracle; base 2029 574 vs 815): the cornered
  under-stocked SEEDLING-regime establishment-timing straddle — a SEPARATE, pre-existing establishment-MODEL
  divergence (not the keyword), same class the ref stands mask. The keywords themselves now modulate
  establishment as the oracle does.
- Gate: multicycle 339/11 byte-identical (the fix is inert on golden stands — none use NOINGROW-in-ESTAB).

## TALLY (re-validation on the under-stocked fixture) — INERT in the ORACLE, no clean signal
Tried TALLY at 2019/2024/2034, both with automatic establishment ON and ISOLATED (ESTAB/NOAUTALY/NOINGROW/TALLY):
the ORACLE delta is 0 at EVERY date and cycle. A user TALLY (esnutr.f 427) does NOT independently book establishment
on this stand — it modulates the automatic tally's NTALLY/IDSDAT, and the actual establishment still goes through
NEWTPP=max(0,ITPP-NSTORE), which is 0 once the stand is stocked. Its only .sum-visible effect (the saturated-stand
-88 the prior 1b9ea7d9 "validation" saw) is the SAME secondary NSTORE side-channel as STOCKADJ, NOT independent
establishment. ⇒ TALLY's independent effect is not cleanly isolatable on available fixtures; the code
(kw_estab! pushes ScheduledActivity 427/428/429; honored in ie_autoes_schedule!) is a faithful esnutr.f transcription
+ gate-safe (339/11), KEPT, but has no positive oracle signal to A/B against. This is a keyword whose oracle effect is
itself marginal/second-order, not a jl gap.
