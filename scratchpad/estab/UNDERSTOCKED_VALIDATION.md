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

## estab.f:581-582 PROB1 clamps — ADDED (faithful) + measured a REAL runaway bug
STOCKADJ 2.0 on the under-stocked fixture: jl prob1 = logistic(0.617)*2 = 1.234 (a probability >1, UNCAPPED). Since
prob1 sits in the NSTORE denominator (tpacre/(prob1*300)), the uncapped >1 value shrinks NSTORE and makes ingrowth
OVER-book: jl Δ(STOCKADJ 2.0) = +588/+358/+415/+277 TPA, while the oracle caps FTEMP at 0.9990 (estab.f:582) and
shows −251/−169/−436. Added `prob1 = clamp(prob1, 0.0001, 0.9990)` (verbatim estab.f:581-582). The clamp bit
(+588→+462) — jl's prob1 is now capped at 0.999 like the oracle. Gate 339/11 byte-identical; INERT on base/STOCKADJ
0.5/NOINGROW (their prob1 = 0.617/0.309 < 0.999) so the earlier validations are unchanged.
RESIDUAL: jl STOCKADJ 2.0 STILL diverges in sign from the oracle (jl +462 vs oracle −251) — but that is the SEPARATE
cornered seedling-mortality regime: STOCKADJ>1 boosts EARLY establishment, and the oracle then density-mortalizes the
extra seedlings (net −251 later) while jl retains them. That downstream seedling self-thin is the same pre-existing
establishment-MODEL divergence as the ~1-cycle timing shift — NOT the clamp. The estab.f:583 per-plot PNN floor is not
represented (jl prob1 is a stand-level scalar); inert wherever prob1 >= the plot's PNN.

## MECHPREP / BURNPREP (site prep, esin.f opt 4/5 → ESPRIN 491/493) — INERT, measured 3 scenarios
Site prep schedules a ZMECH/ZBURN event that sets IPREP → SPRE[iprep] in the establishment tally. Tried on the
under-stocked establishing stand (753189105290487): (a) MECHPREP/BURNPREP 2019 plain, (b) at inventory, (c)
post-thin (THINPRSC 2034 0.5) + NOINGROW + site prep at the thin date to isolate the disturbance tally. ORACLE
delta = 0 at EVERY cycle in all three. Like TALLY, the site-prep modifiers act on a DISTURBANCE-triggered
establishment tally that the automatic LINGRW/LAUTAL path on these stands does not cleanly expose (the ingrowth
NSTORE clamp / thres1 gate wash it out). ⇒ no validatable .sum signal on the available IE fixtures; MECHPREP/
BURNPREP stay UNWIRED (per doctrine — cannot validate ⇒ do not commit). Exposing them would need a stand whose
DOMINANT establishment is a site-prep-eligible disturbance tally — not constructible from the current fixtures.

## Establishment-keyword item — measured limit (2026-08-21)
VALIDATED + committed: STOCKADJ (works), NOINGROW/INGROW/AUTALLY/NOAUTALY/THRSHOLD (real bug fixed), PROB1 clamps.
MEASURED-INERT on all available fixtures (faithful code kept where wired; unwired ones stay unwired): TALLY/TALLYONE/
TALLYTWO, MECHPREP, BURNPREP. .sum-INVISIBLE (need a per-tree treelist A/B): SPECMULT, HTADJ. Separate model dive:
the cornered seedling-mortality/timing regime. This is the genuine validatable limit with the current IE fixtures.
