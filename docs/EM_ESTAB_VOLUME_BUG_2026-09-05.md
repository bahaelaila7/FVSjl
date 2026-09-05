# EM establishment-cycle volume bug — stale-slot cuft_vol leak (FOUND + FIXED 2026-09-05)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Branch: fix-em-estab-volume (off master 72602f09). Found during the EM FIA cap-and-fix sweep at cap-3 (cursor
30000): the supervisor's guard STOPPED on a VOLBUG (worst_col=TCuFt, vol_abs>=300) instead of auto-cornering.

## Symptom
~41 dense EM stands (cursor 26000-30000) flagged worst_col=TCuFt with total-cubic-volume divergences of
300-2601 cuft and huge max_rel% (up to 15337%). Signed TCuFt tally on a sample: mostly downstream of the
benign count-straddle, BUT the extreme-% cases (5440744010690, 373819114489998, ...) showed a DISTINCT
one-directional bug: jl over-reports cubic volume for SEEDLING-stage stands.

## Root cause (instrumented, definitive)
Trajectory of 5440744010690: at 2015 O[TPA=120,BA=0,QMD=0.1,TCuFt=0], J[same TPA/BA/QMD, TCuFt=2601]. Identical
tree list (0.1" seedlings), yet jl TCuFt=2601 vs oracle 0. VTOT_DBG dump: a live seedling sp=3 dbh=0.1012
carried cuft=22.6 (×114.5 TPA = 2588). compute_volumes_em!'s OWN per-tree output was CORRECT and small
(the d<1 gate zeros sub-inch trees) — so the summary was reading a value compute_volumes! never produced.

Mechanism: the stand's inventory trees are ALL dead (HISTORY=6 snags, 5-18" ponderosa); the live stand is built
entirely by AUTOES establishment. `compute_volumes!` computes cuft_vol for BOTH live AND dead trees (loops
`1:(t.n+t.ndead)`), so the dead-block slots (t.n+1…t.n+ndead) hold the snags' real cubic volume. `compute_volumes!`
runs BEFORE establishment in grow_cycle!. Establishment then inserts regen at `n = t.n+1` (reusing dead-block
slots) and sets species/dbh/height/tpa/crown/etc. but NOT the volume fields — so each regen INHERITS the stale
dead-tree cuft_vol. The .sum (built AFTER establishment) sums `1:t.n` of cuft_vol ⇒ includes the stale seedling
volume. One-directional (jl>oracle), intermittent by cycle (depends on regen count/slot reuse), reconverging once
next cycle's compute_volumes! overwrites with the correct (tiny) values. FVS semantics (grow_cycle! comment):
"newly-established trees get NO volume in their birth cycle."

## Fix (faithful; 2 sites — the third already had it)
Zero cuft_vol/merch_cuft_vol/saw_cuft_vol/bdft_vol on regen insertion:
- src/engine/establishment.jl (~622) — shared PLANT/NATURAL/scheduled ESTAB path.
- src/variants/inlandempire/establishment.jl (~2066) — ie_autoes_establish! (AUTOES; used by IE AND EM).
The sprout path src/engine/sprout.jl (esuckr!, ~794) ALREADY zeroed them with the comment "the slot may hold a
previously-deleted record" — confirming the pattern and the two missed sites. Purely corrective: all cuft_vol
reads are REPORTING-only (accretion/mortality .sum columns + snapshots), none drive growth/mortality, and the
accretion loop covers only original trees (1:n captured pre-ESTAB) — no feedback.

## Validation (vs live FVSem_clean, both fixes)
- 5440744010690: worst TCuFt 2601→0 (2015), 472→8 (2035); max abs dTCuFt now -10 (residual = benign count-straddle).
- 373819114489998: 1768→0 (2025), 1235→8 (2045); max abs dTCuFt now -9.
- 9865845020004 / 886213320290487: volume tracks their cornered VARMRT count-straddle (dTCuFt follows dTPA), unchanged.
Gate: Pkg.test() (339/11) — [see commit].

## Scope / who else is affected
General across variants that establish regen INTO reused dead/removed slots that carry volume. Most visible on
all-dead-inventory (HISTORY=6) + heavy-establishment FIA stands (EM/IE Great Plains). Affects TCuFt/MCuFt/BdFt
in the establishment cycle(s) only. IE also benefits (shared AUTOES path).

## Sweep impact
The EM sweep is PAUSED at cursor 30000 (cap-3). This bug produced the TCuFt escalations that tripped the guard.
Once this fix is on master, re-sweep 26000+ (the cap-3 TCuFt escalations should collapse to the benign class).
Cap-1 (0-16000) + cap-2 (16000-26000) cornerings were TPA/BA count-straddle adjudications, largely independent
of this volume fix, but the TCuFt-worst stands among them are worth a re-check post-merge.
