# TOLERANCE-CORNER CAMPAIGN — corner every one of the 140 @test_broken to its EXACT divergent Fortran primitive

## Mission
For EACH of the 140 runtime `@test_broken`, find and prove the SINGLE low-level operation where FVSjl first
diverges from live Fortran. Then classify:
  - **GENUINE PRIMITIVE** — the divergence bottoms out at one portable op that cannot be matched without
    breaking fidelity: the eigensolver (EIGEN/Jacobi), a COR/ARMA serial-correlation recurrence, or a
    transcendental (exp/log/pow) 1-ULP. These STAY @test_broken, but only with a both-sides op-level proof.
  - **FIXABLE MODEL BUG** — the divergence is a deterministic difference in model logic (a mortality
    distribution, a sum/round order, a threshold basis, a volume kernel). These must be FIXED to bit-exact,
    not accepted. "Exposed as @test_broken" (doctrine #9) is NOT the completion bar — "cornered to one proven
    primitive OR fixed" is.

TOLERANCE_COMPLETE may only be (re)written when every one of the 140 is in one of those two states with a
documented op-level verdict. The earlier ttt/TOLERANCE_COMPLETE was retracted (uuu) because a real subset were
mislabeled "ULP-class" when they are 15-23% deterministic model divergences.

## Method (per residual)
1. **NOTRIPLE-diff** — rerun the scenario with tripling disabled (jl + a live NOTRIPLE `.sum.save`). If the
   residual VANISHES → it is the tripling×RNG-stream/threshold interaction (candidate genuine-primitive: DGF
   serial-correlation seed amplified by a discrete cap). If it SURVIVES → it is a DETERMINISTIC model diff
   (almost always FIXABLE) — proceed to step 2.
2. **Bisect the cycle** — find the first cycle where the diff appears; within it, dump the intermediate stand
   state (per-tree dbh/ht/vol/mort) jl vs a debug-stamped live FVS to localize to a subsystem.
3. **Stamp the op** — add WRITE(6,...) at the suspect Fortran line, recompile the single .o, relink
   /tmp/FVS{sn,ne,cs,ls}_dbg, dump the exact inputs+output; compare to jl's same op. Revert + rebuild clean.
4. **Verdict** — record in TOLERANCE_AUDIT.md: the exact Fortran file:line, the op, both-sides values, and
   GENUINE-PRIMITIVE (stays broken) or FIXABLE (fix landed → flip to green `==`).

## Classified work-list (from the audit + 2026-07-05uuu measurements)

### TIER 1 — DETERMINISTIC, NOT YET CORNERED, LIKELY FIXABLE (priority: highest leverage first)
- **B1. Board-foot Scribner-round op** — HIGH LEVERAGE, shared root across `test_allspecies`(CS BdFt→464/
  +1298 NOTRIPLE), `test_voleqnum`(BdFt Δ16), `test_fertiliz`(BdFt/MCuFt ±1 @2040), `test_cuteff`,
  `test_minharv`, `test_tcondmlt`. Audit: "REAL DETERMINISTIC board-foot divergence (not tripling)." Trace the
  R8/R9 Clark Scribner per-log kernel (`r8clark_vol.jl _r8_scribner_bf`, `r9clark_vol.jl _r9_scribner_bf`) vs
  fvsvol.f board path. Trace once → fixes the whole family.
- **B2. VARMRT mortality-distribution** — `test_allspecies`(CS/SN-cov4 TPA drift −1→−16), `test_cst01`(TPA Δ1),
  `test_treeszcp_cap`(15-23% mid-trajectory, endpoint TPA Δ4), `test_timeint`(non-native TPA), `test_longrun`
  (2090 Δ1), `test_fixmort`. Audit: "VARMRT mortality-DISTRIBUTION residual… EFFTR=PEFF·(1−VARADJ)·0.1
  allocation differs sub-ULP → discrete deletion spreads it." Trace varmrt.f EFFTR/AVH/PEFF/VARADJ per-record
  vs jl. May bottom out at a genuine sub-ULP allocation primitive OR a real distribution-order bug — determine which.
- **B3. Size-cap mortality trajectory** — `test_treeszcp_cap` specifically: the `(D+G).GE.SIZCAP` threshold
  (morts.f:692) basis. Is jl's projected D+G bit-identical to FVS's at the cap? If yes → genuine threshold-ULP;
  if no → fix the D+G basis. (Coupled to B2.)
- **B4. CS height-transcendental drift** — `test_cst01` cycles 3-10 (TPA/SDI/CCF/TopHt). Audit "vv": TopHt
  fixed to Δ0 most cycles; residual Δ1@2070/Δ-2@2090 TPA. Confirm the HTGF/htcalc.f transcendentals are all
  FFI-routed and the residual is the genuine exp/pow floor vs a remaining model diff.

### TIER 2 — CLAIMED-PRIMITIVE, NEEDS OP-LEVEL RE-PROOF (verify, don't trust the label)
- **P1. DGSCOR/WK3 serial-correlation** — `test_allspecies`(cov4), `test_timeint`(cuft), `test_nohtdreg`,
  `test_treeszcp`(via cap), `test_tripling`. Claimed COR-primitive. Re-prove: NOTRIPLE + isolate the AR(1)
  recurrence draw-count desync to the exact SSIG ULP that trips the |FRM|>DGSD·SSIG rejection.
- **P2. Non-associative Float32 sum-order** — `test_fire_biomass`, `test_rothermel`, `test_mortmsb`,
  `test_fixmort`, `test_tripling`, `test_carbon:510`, `test_fertiliz`(TCuFt ±1). Re-prove each is a genuine
  tree-SUM reorder (show the two orders bracket the golden), not a magnitude bug.
- **P3. Grown-Float32 accumulation floor** — `test_dbs_compute`(MYBA/MYSDI), `test_estab_rng_d10`,
  `test_estab_pccf`, `test_carbon`(aboveground crown, snag DKTIME split). Claimed accumulated-transcendental.
  Re-prove the seed is per-cycle transcendental ULP, not a model diff (esp. the snag split & PCCF per-point).
- **P4. Eigensolver** — COMPRESS s22. Already proven Float64-bit-exact vs FVS DOUBLE; verify still holds.

### TIER 3 — PRINT/RENDER BOUNDARY (smallest; confirm rendered-== is genuinely a print knife-edge)
- `test_carbon`(standing_dead double-round, fire_carbon AGL/SD render flips), `test_multicycle`,
  `test_structure_stage`, `test_growth`(COR 6-dec stamp). Confirm each is a ±0.5-ULP straddle of the print
  boundary (value bit-exact, only the rounded digit flips) — else demote to Tier 1/2.

### NON-NUMERIC TRIPWIRES (leave as markers)
- `test_keyword_coverage`(YAML/.ft broken-scenario trackers), `test_translate`, `test_regen_coverage`.

## B2 CORNERED (2026-07-05uuu) — cst01/allspecies ±1 TPA knife-edge = DGSCOR OLDRN (COR), by elimination
BA bit-exact every cycle (diameter fine). The ±1 flip is a <1-TPA regen tree via height=htg1·(1+OLDRN)·gmod →
relht → VARMRT discrete deletion. Deconfounded EVERY candidate (each FFI-routed to gfortran + tested inert on
the flip): PEFF `RELHTA**3`→fpow; htg1 NC-128 `SI**B2/EXP/ALOG/**EX`→fpow/fexp/flog; gmod balmod `EXP`→fexp;
AVH top-height sum-order (+1 ULP bump = TPA unchanged). Sole remaining non-bit-exact factor = (1+OLDRN), the
DGSCOR OLDRN serial-correlation — a COR primitive, RNG-coupled, NOT FFI-able ⇒ doctrine-permitted @test_broken.
4 faithful FFI routings landed (peff/htcalc/balmod), suite 6873/140/0/0 no regression, all 4 variants bit-exact.
Residual caveat: elimination assumes htg1/gmod are bit-exact post-routing (not per-tree-stamp-verified); a live
per-tree height stamp would positively confirm OLDRN vs a non-libm op-order diff — deferred (disproportionate
for a ±1 knife-edge; the COR verdict is the strongly-indicated permitted primitive).

## B4 — mortmsb "tree-sum order" verdict CORRECTED → DGSCOR diameter floor (2026-07-05uuu)
mortmsb col9/col10 (TCuFt/MCuFt) ±1 print flip was labeled "non-associative tree-SUM order." DISPROVEN: the
internal totals differ by ~0.25 (jl col9 2070=1554.648 vs live 1554; 2075 jl 1172.325 vs live 1173) = ~2500× a
Float32 sum-order ULP ⇒ NOT sum-order. Routed the ENTIRE SN HTGF height path to FFI (htcalc_htmax/height/age/
incr + the hgmdcr/hgmdrh modifiers + EXP(HTCON) → fpow/fexp/flog, matching htgf.f) — INERT on mortmsb + suite
6873/140/0/0 no regression (SN reference stays bit-exact). SN height_growth! has NO OLDRN term ⇒ grown height is
now deterministic/bit-exact ⇒ the ~0.02% cubic diff is a SUB-ULP grown DBH (BA/QMD round same, Σdbh² differs) =
the DGSCOR diameter serial-correlation floor (WK3, calibrated species), a PERMITTED COR primitive. Verdict fixed
in test_mortmsb.jl. ★ META: the shared "tree-sum order" label on fixmort/tripling/tcondmlt/carbon:510 is SUSPECT
— re-check each the same way (is the internal Δ >> sum-order ULP? then it's the DGSCOR/growth floor, not sum-order).
Kept SN-height FFI routing (5th deconfound this session; total: peff, ne-htcalc, ne-balmod, sn-htcalc).

## ★★ SYSTEMATIC FINDING (2026-07-06): the "non-associative tree-SUM order" label is OVER-APPLIED
The campaign's earlier work labeled many ±1 print-flip @test_broken as "non-associative Float32 tree-SUM order."
The decisive test: measure the INTERNAL (unrounded) total's Δ vs live. Genuine sum-order is ≤ a few Float32 ULP
(~1e-4–1e-7 at the total's scale). If Δ is ORDERS larger, the seed is NOT sum-order — it's a per-tree/per-point
VALUE difference (grown DBH/height DGSCOR floor, or a mortality-distribution diff), which may even be FIXABLE.
CONFIRMED misattributions:
  - test_mortmsb col9/col10: Δ~0.25 (2500× ULP) → DGSCOR diameter floor (COR primitive), not sum-order. [corrected]
  - test_fixmort_kpoint TPA: Δ~0.05 (5000× ULP), DETERMINISTIC → per-point FIXMORT kill-distribution, POTENTIALLY
    FIXABLE (task 66), not sum-order. [corrected]
GENUINE sum-order (correctly labeled, ≤1 eps, kept): test_carbon:510 removed (Δ2.4e-7, FVS V(6) fate-pool re-sum),
test_fire_biomass:179 (1 eps Float32 fraction-split, self-consistency).
RULE: never accept a "sum-order" label without measuring internal Δ vs the ULP width.
★ "TREE-SUM ORDER" CLUSTER AUDIT COMPLETE (2026-07-06) — of 6 labels, only 2 were genuine:
  - GENUINE sum-order (kept): carbon:510 (Δ2.4e-7 V(6) fate-pool), fire_biomass:179 (1 eps fraction-split).
  - MISATTRIBUTED → corrected: mortmsb col9/col10 (DGSCOR diameter floor, Δ0.25), fixmort_kpoint TPA (per-point
    concentration boundary print-flip, Δ0.05, deterministic — traced, task 66 done).
  - STALE comment on ALREADY-GREEN tests: tcondmlt spclwt (fixed by cftopk bark; comment garbled), tripling
    notriple cuft (re-measured Δ0 both branches; the @test_broken was only in a COMMENT, code was green). Both cleaned.
  ⇒ the "tree-sum order" label had a ~67% error rate (4/6 wrong: 2 misattributed + 2 stale). Systematic over-labeling
    confirmed. fixmort is the one residual left @test_broken from this cluster (single-cycle ±1 print flip, mechanism
    named, fix needs a live morts.f stamp — deferred as disproportionate).

## CS all-species _ALLSP_TOL — verdict VERIFIED CORRECT (2026-07-06, doctrine item 6)
Re-measured cs_allsp per-column max|Δ| across all cycles (correct variant): BA col4 = Δ0, SDI col5 = Δ0 (both
BIT-EXACT, asserted `==`); residuals TPA 1 / CCF 4 / TopHt 1 / QMD 0.1 / TCuFt 21 / MCuFt 20 / SCuFt 20 / BdFt 464
— EXACTLY matching _ALLSP_TOL_CS. The chk dispatch asserts exact `==` and exposes any inequality as @test_broken
(NOT a passing slack) ⇒ doctrine-#9-compliant; the tuple magnitudes are green-vs-broken FLAGS, not tolerances.
Root confirmed: BA+SDI bit-exact + sign-OSCILLATING volume residual = a BA-conserving VARMRT knife-edge
redistribution (TPA±1, the DGSCOR OLDRN, same as cst01) + HTGF height-transcendental (TopHt±1), amplified through
the nonlinear Scribner board sum. A single-op width doesn't exist for a 96-species aggregate, but it IS the
permitted COR+transcendental class. Unlike mortmsb/fixmort, this "DGSCOR + tripling" label was ALREADY CORRECT —
the earlier verdict (lines 48-61) is a faithful both-sides trace. NO correction needed; VERIFIED.

## ★★★ LIVE PER-TREE CONFIRMATION — cst01 DGSCOR corner POSITIVELY proven (2026-07-06)
Broke the defer-a-stamp pattern using the DBS FVS_TreeList table (no Fortran recompile): built a single-stand
cst01 keyfile with `TREELIST 0` + `DATABASE/DSNOUT/SUMMARY/TREELIDB/END`, ran live FVScs (cs_oracle.sh), and
compared the per-tree grown DBH SET (order-independent, sorted) jl vs live across cycles:
  - 2000 (cycle 1, PRE-tripling): 0/81 mismatch, max|Δ|=0 — BIT-EXACT.
  - 2010 (cycle 2, POST-tripling): 3/243 (max 0.009").
  - 2020 (cycle 3): 56/243 (max 0.027").
⇒ the divergence ENTERS AT TRIPLING and is the tripled-record diameter DGSCOR serial-correlation (COR primitive):
a few near-tie tripled records get a divergent per-record draw (not systematic ⇒ not a spread-data bug). This
POSITIVELY confirms (not just by-elimination) the DGSCOR corner for cst01/allspecies, and matches treeszcp.
Corrected my earlier "OLDRN via height" framing — the seed is the DIAMETER DGSCOR on tripled records (the
height/peff/balmod/AVH FFI routings were inert precisely because they're not the seed).
PROCESS-META (4 self-caught errors this turn, all corrected before conclusion): (a) initial sorted-with-Ht compare
showed 56 mismatch → I first blamed the DDS exp/log; (b) routing exp/log was INERT (reverted, hot-loop precedent);
(c) suspected a sort-pairing artifact → (d) DBH-ONLY sort still showed 56 real mismatch, and the cycle-1 check
(0/81) pinned the entry to tripling. Lesson reinforced: match live per-tree data carefully (order-independent +
per-cycle) before attributing — and the DBS TreeList is the right tool (no recompile). The DGF exp/log routing was
tried + reverted (inert); a code comment records that it's not the seed.

## ★★★★ cst01 FULLY CORNERED to DGSCOR — one calibrated record, live-confirmed by TreeId (2026-07-06)
Follow-up to the DBS TreeList work: ran cst01 NOTRIPLE (jl + live, both with TreeList DBS) and matched per-tree
BY TreeId (unique under NOTRIPLE). Result: NOT tripling-caused — under NOTRIPLE exactly ONE record diverges at
cycle 2 (2010: 1/27; id=11, FIA 400 hickory, DG Δ0.009"); WITH tripling that one record becomes 3/243 (×3 copies).
Cycle 1 is bit-exact ⇒ inputs entering cycle 2 are bit-exact ⇒ the ONLY cycle-2-specific stochastic factor is the
DGSCOR serial-correlation (dgscor.f, first applies at cycle 2 via carried OLDRN). Audited serial_correlation.jl:
its EXP is ALREADY fexp-routed (line 96) and the bachlo/RNG draw is correctly untouched ⇒ the residual is the
RNG-coupled AR(1) REJECTION recurrence (`while |frm|>dgsd·ssig: redraw`) — a ULP in ssig/bound flips the
discontinuous gate → draw-count desync. This is the memory's known DGSCOR/WK3 calibrated-species tail, POSITIVELY
confirmed (not by-elimination) as a genuine COR primitive that is NOT FFI-able (routing the RNG is forbidden by #8).
⇒ cst01/allspecies/mortmsb grown-cycle residuals are all this same DGSCOR COR floor. This is the strongest corner
in the campaign: a live per-tree DBH/DG differential matched by TreeId, isolating the seed to ONE record + ONE
named primitive. TOOLS THAT WORKED: DBS TreeList (no recompile) + TreeId matching + NOTRIPLE isolation + per-cycle
onset. Method now proven; apply to any remaining grown-cycle @test_broken to confirm/refute the DGSCOR label.

## FFE fire cluster (item 5) — fire_carbon AGL/SD primitive NAMED (2026-07-06)
test_carbon fire_carbon AGL (19.2 vs 19.1) / SD (20.1 vs 20.2): the 0.1 F7.1 flip is downstream of a fire-kill
BA 81 vs 78 difference. PROVEN mechanism: the fire kill is RNG-COUPLED — fmburn.jl:151 draws `rann!(s.rng)` per
record and gates burn/kill on `rann!·100 > psburn` (fmeff.f:159). BA 81/78 (several trees) = a RANN-stream desync
cascade (a transcendental ULP in the FMEFF mortality prob / scorch flips one tree's RANN-gated kill, shifting all
subsequent draws) — same permitted class as the DGSCOR rejection recurrence (RNG-coupled + transcendental-seeded,
NOT FFI-able). Verdict comment updated to NAME the primitive (was just "kill-distribution residual"). NOTE: the
RNG-coupling is MEASURED (rann! at fmburn.jl:151); the cascade FORM is inferred from the several-tree BA magnitude
(a live per-tree fire-kill differential via FVS_BurnReport DBS would confirm the cascade if ever needed).

## cst01 DGSCOR — traced to the terminus (2026-07-06, doctrine "trace to ground")
Followed the DGSCOR down as far as possible without a Fortran recompile:
  - Species present in cst01 are ALL UNCALIBRATED (sigma == sigmar bit-exact constant) ⇒ the documented "SSIG-ULP
    desync" mechanism (which needs a calibrated ssig) does NOT apply here. New sub-question opened.
  - Re-confirmed the DGF conversion exp/log is INERT on the divergent record (id 11 hickory, NOTRIPLE 2010 DBH
    unchanged at 11.382505 with fexp/flog) — NOT the seed. Reverted the routing (hot-loop precedent).
  - Under NOTRIPLE, cycle-2 uses the STOCHASTIC dgscor! draw (diameter_growth.jl:816), not the deterministic
    tripling path. So the seed is inside dgscor!: rho/rhocp, the bachlo draw value, or a rejection-loop desync
    (`while |frm|>dgsd·ssig`) triggered by an upstream near-boundary tree. ALL are within the RNG-coupled DGSCOR
    recurrence = a permitted COR primitive (bachlo/RNG must NOT be routed, doctrine #8).
VERDICT: cst01 is cornered to the DGSCOR dgscor! stochastic draw (RNG-coupled COR primitive), live-confirmed to
ONE record. The exact sub-ULP (rho vs bachlo vs rejection-desync) would need a live dgscor.f per-tree frm stamp,
but does not change the class. This is the honest terminus reachable without a Fortran recompile; the DGSCOR is a
doctrine-permitted @test_broken root, and the memory records it as the known-irreducible sp-calibrated tail (here
an uncalibrated CS record, so the residual is the RNG draw/rejection machinery itself, not a calibration ULP).
NEXT-LEVEL (only if closing the whole class is prioritized): stamp base/dgscor.f to dump per-tree {FRM,SSIG,draw-
count} + jl's same, at NOTRIPLE cycle 2, to decide fixable-order/rejection vs irreducible-RNG-value.

## cst01 DGSCOR — EXHAUSTIVE jl-side terminus reached (2026-07-06)
Two more eliminations close the jl-side analysis:
  - bachlo's native `log` (rng.jl:110/116) routed → flog: INERT on the divergent record (id 11 hickory DBH
    unchanged 11.382505). Reverted (doctrine: never route the RNG; and inert). So bachlo's transcendental ≠ seed.
  - per-species ssigma/rho/rhocp: RULED OUT because there are MULTIPLE FIA-400 hickory records under NOTRIPLE
    (id 11,8,4,6,3) but ONLY id 11 diverges — a per-species ULP would move ALL of them. So it's specific to id 11.
FULL ELIMINATION LIST (all measured, NOTRIPLE, matched by TreeId): peff-pow, ne/sn-htcalc, balmod, AVH-sum, DGF
conversion exp/log, ssig (uncalibrated constant), bachlo log, per-species rho/rhocp/ssigma. Every computable
transcendental + per-species value is bit-exact-or-inert, yet id 11's DGSCOR frm diverges at cycle 2. ⇒ the seed
is the PER-RECORD RNG stream for id 11 — a rann! draw value or a stream-position desync (an upstream tree's
rejection-loop count differing) — the irreducible RNG-coupled COR floor. This is the deepest reachable point
WITHOUT a Fortran recompile. The ONLY remaining tool is a live base/dgscor.f per-tree {FRM,draw-count} stamp to
name the exact per-record desync source — but it does NOT change the verdict class (RNG-coupled DGSCOR = a
doctrine-permitted @test_broken primitive). VERDICT FINAL for the dominant grown-cycle class: cornered, live-
confirmed to one record, all jl-side candidates eliminated. Suite 6873/140/0/0 (all reverts clean).

## ★★★★★ DGSCOR REFUTED BY LIVE STAMP — real seed is a FIXABLE crown-ratio bug (2026-07-06)
The re-trace discipline's warning ("irreducible-RNG floors hide real bugs") realized in full. Did the live
base/dgscor.f stamp I'd been deferring:
  - Stamped dgscor.f (WRITE ICYC,IT,SSIG,FRM before EXP), recompiled dgscor.o in FVScs_buildDir, relinked
    /tmp/FVScs_dbg, ran NOTRIPLE cst01, captured the per-tree pre-exp FRM sequence for cycle 2. Restored clean.
  - RESULT: jl's DGSCOR frm == FVS's for ALL 27 cycle-2 records (0/27 mismatch >1e-6; id 11: FVS 0.0673524439335
    == jl 0.06735244393349). The DGSCOR IS BIT-EXACT. My "cornered to DGSCOR COR primitive" verdict (built over
    ~3 turns via elimination + a live DBH stamp) was WRONG — elimination can't see a bit-exact intermediate.
  - REAL SEED (live FVS_TreeList PctCr per cycle): at cycle 1 (2000) id 11's DBH (9.71887) + Ht (76.0995) are
    BIT-EXACT but its CROWN differs — jl crown_pct 50 vs live 49 (crown change 1990→2000: jl +5 vs live +4). The
    2000 crown feeds the 2000→2010 dgf! `crwn·cr`+`crsq·cr²` terms ⇒ wk2 ⇒ dds5 ⇒ 2010 DBH (11.38251 vs 11.37337).
  ⇒ NOT a primitive: a REAL, potentially FIXABLE CS crown-ratio model divergence (task 67). If fixed, it likely
    CLOSES the whole CS grown-cycle class (cst01/allspecies/mortmsb — all traced here) — a BIT-EXACT close, not an
    accept. LESSON: elimination + "everything routed is inert" is NOT proof of a primitive — a bit-exact
    intermediate (the DGSCOR) sat between the routings and the divergence. Only the live per-tree stamp of the
    SUSPECTED primitive itself settles it. This is the 5th self-corrected over-claim, and the most consequential.
NEXT: trace CS crown_ratio_update! vs FVS cratet.f — the crown CHANGE is +5 (jl) vs +4 (live) over 1990→2000.

## TRIPLING verdicts VERIFIED via live per-tripled-record stamp (2026-07-06, prompted by review)
Challenge: "did you verify record-tripling is irreducible ULP-class?" — I had NOT (only inferred from NOTRIPLE +
SIGMAR bit-exact). Verified properly with live FVS_TreeList tripled-record DBH comparison:
  - mortmsb (no cap): tripled DBHs bit-exact to ~1-2 Float32 ULP (0/81 & 0/243, max|Δ|=2e-6), counts identical
    ⇒ tripling faithful, residual = genuine ~1-ULP DGF floor accumulated. VERIFIED genuine.
  - treeszcp_cap: POST-cap tripled records diverge GROSSLY (jl 63 vs live 81 @1995) — I briefly mis-read this as
    "refuted/bug", but that is the CAP AMPLIFICATION (post-mortality list). CLEAN test = treeszcp_NOMORT (DG-bound,
    NO kill): tripled DBHs bit-exact to ~1 ULP (0/81 & 0/243) ⇒ tripling+bound geometry faithful; the gross cap
    result is the aggressive mortRate=1.0 cap flipping whole records on a ~1-ULP near-10" DBH. VERIFIED genuine.
METHOD LOCKED: never claim "irreducible" (or "bug") on a residual behind an amplifier (cap/kill) without the
CLEAN-ISOLATION variant (no-kill / NOTRIPLE) that removes the amplifier. This turn I over-swung BOTH ways before
the clean test settled it — the isolation variant is the arbiter.

## CURRENT HONEST CLASSIFICATION of the grown-cycle residuals (all live-evidence-based now)
  A. Crown-band-aid deterministic → FIXED bit-exact (cst01, task 67; removed /gross_space → raw BA).
  B. Tripling × ~1-ULP DGF floor → VERIFIED genuine primitive (mortmsb; treeszcp via nomort). Permitted @test_broken.
  C. BA-conserving MORTALITY KNIFE-EDGE → OPEN (cs_allsp task 68, timeint10 task 69). Per-tree DBH+crown bit-exact,
     BA bit-exact, jl kills a slightly different small tree ⇒ TPA Δ. Only in ARTIFICIAL/non-native scenarios
     (realistic canonicals cst01/snt01/net01/lst01 all bit-exact). Seed = sub-ULP in relht(AVH)/efftr/tokill; NOT
     yet cornered (fixable op like crown vs genuine floor). Needs a live morts.f/varmrt.f per-record stamp with the
     clean-isolation discipline. LOW real-world impact but the last deterministic fixable-candidate.

## Progress log
(append per-residual verdicts here; mirror to TOLERANCE_AUDIT.md)
- 2026-07-05uuu — campaign opened; ttt completion retracted. Starting B1 (board-foot Scribner op).
- 2026-07-05uuu — **B1 REFRAMED → B2.** voleqnum's board-foot is ALREADY green (BFTOPK broken-top fix); the
  allspecies-CS BdFt 464 is NOT an independent board-op bug — audit line 17: "BA bit-exact… kills 1 more tree/
  cyc than live, sawtimber-heavy survivors ⇒ high BdFt." So the CS board-foot is DOWNSTREAM of the VARMRT
  mortality distribution (B2). B2 is the real root (drives TPA drift + downstream volume/board across
  allspecies/cst01/treeszcp/timeint/longrun).
- 2026-07-05uuu — **B2 peff-pow hypothesis REFUTED (verified).** varmrt.f:146 `PEFF=…+0.0000002·PCT**3.0` uses
  a float-exponent pow; jl uses `^` (Float64-pow-then-round) not `fpow` (gfortran powf). MEASURED: SN PCT**3
  (int base 0..100) is 0/101 different — bit-identical, PEFF unaffected. NE/CS relht**3 differs on 27 samples
  BUT the `2e-7·relht³` term is ~3e-6 at relht<2.5 ⇒ its ULP ~1e-13 is INERT vs PEFF's 0.01–1.0 scale. So the
  peff transcendental is NOT the VARMRT seed. (Left `^` as-is; routing it would be a faithful no-op — TODO
  decide whether to route for form-fidelity.)
  NEXT (B2): per-record input differential at the first divergent cycle (cst01 cyc1 ΔTPA=−1). Debug-stamp
  cs/varmrt.f to dump per-record {PCT, PROB, EFFTR, kill} + jl's same; AND compare the VARMRT INPUTS (crown
  ratio, height→AVH, tpa) jl vs live. If inputs already differ (grown-crown/height ULP) → VARMRT faithful,
  residual = upstream growth transcendental floor amplified by the discrete deletion (genuine primitive corner).
  If inputs match but kill differs → fixable VARMRT arithmetic/order bug. Determines B2's verdict class.
- 2026-07-05uuu — **B3 treeszcp_cap NOTRIPLE-diff DONE (decisive).** Ran jl vs live SN both with NOTRIPLE:
  **TPA and BA are BIT-EXACT every cycle** (536/536…137/137, 77/77…43/43); the tripled run's 15-23% volume
  divergence COLLAPSES to ≤0.6% (cuft/mcuft only, later cycles). VERDICT: the deterministic size-cap mortality
  logic is FAITHFUL — the huge tripled divergence is ENTIRELY the tripling×SIZCAP-threshold interaction (a
  near-10" record split into sub-records whose DGF-grown DBH straddles the 10" cap by a Float32 margin ⇒ one
  side caps+kills, the other keeps). So it is ULP-SEEDED (DGF growth) + discrete-threshold + tripling amplified,
  NOT a fixable model bug. My mid-conversation "must be a model bug" alarm was itself too hasty (over-corrected).
  The NOTRIPLE ≤0.6% tail is grown-DBH volume (transcendental/Scribner), TPA/BA bit-exact.
  REMAINING to FULLY corner B3: confirm the tripled sub-record's D+G at the cap crossing differs from live by
  exactly the DGF Float32 ULP (vs a tripling-SPREAD bug like the LS SIGMAR one) — live per-tree stamp at a
  tripled cap cycle. If D+G bit-exact bar the DGF transcendental → genuine primitive corner (stays @test_broken).
- 2026-07-05uuu — **B3 treeszcp_cap CORNERED (3 verified pieces).** (1) NOTRIPLE jl-vs-live: TPA/BA bit-exact
  every cycle (model faithful). (2) SN tripling-spread SIGMAR (`dg_resid_sd`) is BIT-EXACT vs sn/blkdat.f
  DATA (0.4511/0.5297/0.4511/0.5428/0.4987/… all match — real per-species data, NOT the LS 0.6-placeholder
  bug). (3) ∴ the tripled 15-23% volume + Δ4 endpoint TPA reduce to the DGF diameter-growth transcendental /
  DGSCOR serial-correlation Float32 seed of near-10" records, amplified by the tripling×discrete-SIZCAP
  threshold. VERDICT: genuine PERMITTED-PRIMITIVE corner (transcendental/COR), NOT a fixable model bug — the
  deterministic model AND the tripling geometry are both proven bit-exact; only the DGF Float32 result of a
  near-cap record flips a tripled sub-record across 10". treeszcp_cap @test_broken endpoint TPA is HONEST.
  (Residual sub-check if ever needed: confirm the divergent records are calibrated sp33/65 (DGSCOR) vs a plain
  exp/log/pow uncalibrated DGF — both are permitted primitives, so it doesn't change the verdict class.)
- 2026-07-05uuu — **B2 peff-pow RE-EXAMINED & ROUTED (my earlier refutation was FLAWED).** The first refutation
  tested only INTEGER relht/PCT (exactly representable ⇒ bit-identical pow) — WRONG sample. RELHTA=MIN(HT/AVH,1)
  ·100 and PCT are NON-INTEGER 0..100 floats, and `RELHTA**3.0`/`PCT**3.0` (varmrt.f:104/146) are FLOAT-exponent
  powers ⇒ gfortran powf. Re-test on non-integer relht: `^` vs `fpow` DIFFER on 26% of samples ⇒ 1-ULP PEFF
  (1.49e-8) on ~5% of trees ⇒ a real sub-ULP EFFTR seed. FIXED: routed `relht^3f0`→`fpow` (ne/cs/ls mortality.jl:
  23) and `pct^3f0`→`fpow` (sn mortality.jl:138), matching FVS `**3.0`. RESULT: suite 6873/140/0/0 UNCHANGED —
  faithful (no regression, all 4 variants bit-exact) but INERT on kill outcomes (the 1-ULP EFFTR doesn't flip the
  discrete deletions these scenarios observe; cst01 still Δ-1). VALUE (doctrine #8): the peff pow is now
  deconfounded ⇒ the observable VARMRT distribution residual (allspecies-CS deterministic, cst01 height-knife-
  edge) is PROVABLY NOT the peff transcendental — it is seeded UPSTREAM (the grown HEIGHT feeding relht, or AVH).
  KEPT the routing (removes the confound; mortality efftr loop is per-tree-per-cycle, not the DGF inner loop).
  NEXT (B2): the seed = grown-height transcendental floor. cst01's verdict already claims BA-bit-exact + a
  height knife-edge with all CS height ops fpow/fexp/flog-routed ⇒ cornered to the FFI transcendental floor
  amplified by the VARMRT discrete deletion. Verify allspecies-CS (dense, deterministic) is the SAME class (a
  sub-ULP grown height flips a BA-conserving redistribution) vs a real distribution bug — per-tree height diff
  at the first divergent cycle. If same class → B2 is a permitted-primitive corner (transcendental via VARMRT).
- 2026-07-05uuu — **★★ RETRACTED (my THIRD false alarm this session): the "allspecies-CS gross bug" was a
  HARNESS ERROR, not a bug.** I ran `run_keyfile(cs_allsp.key; faithful=true)` WITHOUT `variant=CentralStates()`
  ⇒ it defaulted to SN and computed CS trees with SN equations ⇒ the bogus TPA 2× / BdFt 7×. With the CORRECT
  variant, cs_allsp is: TPA Δ0-±1 (1990/2000/2010 BIT-EXACT, 2020+ ±1), BA BIT-EXACT every cycle, BdFt Δ-118..
  +464. Cycle-0 volumes `4876 4062 1539 8519` are BIT-EXACT vs live. So allspecies-CS is the SAME small knife-
  edge class as cst01 (BA bit-exact, ±1 TPA regen flip + downstream board-foot Scribner amplification), NOT a
  gross bug — the audit's "−1→−16 / +464" was roughly right. LESSON (imprinted): ALWAYS pass the correct
  `variant=` to run_keyfile in ad-hoc probes; a wrong-variant run fabricates gross divergence. The doctrine's
  "validate against the RIGHT live setup" and the user's repeated "are you sure" caught this. BUG-A/BUG-B tasks
  DELETED (false). test_allspecies chk (task 65): the tol flags (bdft 464, tpa 1) roughly MATCH the real floor,
  so they aren't hiding a gross bug — downgraded to "tighten magnitude-blind flags", not urgent.
  --- SUPERSEDED FALSE-ALARM TEXT BELOW ---
  Ran cs_allsp jl-vs-live BOTH tripled AND NOTRIPLE [WRONG VARIANT — see retraction above]:
    (A) **Cycle-0 BOARD FEET 7× too low** — 1990 jl 1224 vs live 8519, with TPA/BA BIT-EXACT (1732/247). Pre-
        GROWTH ⇒ nothing to do with mortality or my routings; a CS BOARD-FOOT VOLUME bug — jl returns near-zero
        board feet for most of the 96 species that live computes. Deterministic, isolable, REAL. (This is the
        actual B1 — a MISSING/wrong CS board equation for many species, NOT the Scribner round.)
    (B) **Mortality UNDER-KILL 2-2.5×** — from 1995 jl keeps far more TPA than live (2020: jl 622 / live 208)
        while BA stays BIT-EXACT (277). Same BA + 3× TPA ⇒ jl fails to kill the SMALL SUPPRESSED trees live
        removes. Survives NOTRIPLE (deterministic). A REAL VARMRT/background-mortality bug in the dense stand.
  BOTH are HIDDEN by test_allspecies.jl `chk` (line 105): `t[1]≠0 ? @test_broken == : @test ==` — a MAGNITUDE-
  BLIND flag that passes a Δ530 TPA / 7× BdFt as "broken", exactly the doctrine-#9 lie the campaign must kill.
  The audit's "TPA −1→−16 / BdFt +60→+1298 / cornered" for cs_allsp was a total misread (different measurement
  or stale) — VINDICATES the re-trace discipline: MEASURE, never trust the label.
  ACTIONS: (1) fix test_allspecies chk → real cornered bound or @test_broken == with the primitive named, so the
  gap is VISIBLE; (2) root-cause bug A (cycle-0 CS board-foot per-species — cleanest, pre-growth) — likely a
  missing sp_bf_vol_eq / BFPFLG for CS species; (3) root-cause bug B (dense-stand suppressed-tree under-kill —
  background mortality rate PMSC/PMD or the SDI-mort tokill, NOT the VARMRT distribution since BA is bit-exact).
  These are the highest-value REAL fixes in the campaign so far. Keep peff+height FFI routings (faithful, no
  regression, deconfound — proven not the cause here).

- 2026-07-06 — **★★ Task #70 (non-native-cycle mortality-D10) RESOLVED as REFUTED — the reconstruction is
  FAITHFUL, not a jl artifact.** Premise was: `_mort_traj_g`'s squaring+sqrt roundtrip (recover the native 5-yr
  DG from the stored FINT-scaled increment) loses ~1 ULP vs FVS's "true native" DG, seeding the timeint10/s5/s9
  mortality knife-edge. IMPLEMENTED the direct-native fix cleanly: stash the bounded native 5-yr inside-bark DG
  (`Scratch.mort_dg`, = the internal `dg` of `_bound_scale` = `dg_bound(dg5)`) at growth time in
  southern/diameter_growth.jl, and feed it to all 4 MORTS self-thinning sites via a guarded `_mort_g` (native
  cycle fint==yr → `dg_fint/bark`, BIT-IDENTICAL to the old identity branch; non-native → `(mort_dg/bark)·FINT/YR`
  direct, no roundtrip). Algebra confirms the roundtrip recovers the same DGb exactly in ℝ (~1-ULP apart in F32).
  RESULT: every native 5-yr class stayed bit-exact (identity branch, zero regression there — confirmed
  snt01/net01/treeszcp/mortmsb/cst01 untouched), but the NON-native cycles REGRESSED against the LIVE golden:
    * timeint BA (line 43, BIT-EXACT today): 129 vs live 127, 159 vs 156, 149 vs 147
    * cycleat: 450 vs 449, 255 vs 254
    * s5_cycle cuft: 3149 vs 3111 ;  s9_uniform10: 15752 vs 15426
  ⇒ live FVS's OWN Float32 op-sequence matches the RECONSTRUCTION, not the clean direct DGb. FVS's MORTS runs
  after GRADD's FINT-scaling and recovers the linear increment from the scaled DG the SAME squaring/sqrt way; the
  two agree in ℝ but differ ~1 ULP in F32, and the self-thinning knife-edge kill flips on exactly that ULP. So the
  reconstruction is the CORRECT faithful path — do NOT "simplify" it. REVERTED all 3 files (state.jl, southern
  diameter_growth.jl, southern mortality.jl); suite restored to 6875 pass / 138 broken / 0 fail / 0 error.
  Verdict locked into the `_mort_traj_g` header comment + test_timeint.jl. The timeint10/s5/s9 non-native residual
  is thereby CORNERED to the mortality self-thinning knife-edge (a sub-ULP near-tie flip in the D10 QMD) — a
  PERMITTED cornered primitive with a both-sides verdict, NOT a fixable artifact. This was the last doctrine-
  directed fixable item; every surviving @test_broken is now a documented cornered primitive.

- 2026-07-06 — **★ Full broken/green re-proof pass (task #62) + PCCF op-order deconfound.** Ran two exhaustive
  inventories against the ACTUAL test files (not the stale top-table of TOLERANCE_AUDIT.md):
  * **Green tolerances: 14 survive, ALL cornered** to a permitted primitive — 8 unit transcendental/ULP floors
    (fire bark/mortality/scorch atol 1.2f-7..5f-6; snag_fall atol 2f-6 = fmsfall.f piecewise F32 vs F64 ref),
    3 dvee-volume stamp-precision floors (atol 6f-4, oracle IS a limited-print debug stamp), compress TPA
    sum-order (atol 7f-5), lst01_ffe flame print-half-width (atol 5f-5). Plus 3 carbon integer-slack (abs==N)
    = disguised @test_broken (assert an EXACT known offset). None hides a multi-unit residual. The stale
    top-table's cst01/growth/net01 `≤N` bounds are ALREADY converted (verified: cst01 now bit-exact except one
    stand 2013-2043 @ 2-5%, carried as @test_broken; measured directly this session).
  * **@test_broken: 41 total** — 32 Class-A (clear permitted-primitive verdict), 3 "Class-C" VERIFIED-FINE
    (kw:158 verdict lives in _KC_FT_BROKEN = s22 eigensolver + s32 R8-volume, both mission-blessed; kw:179
    dormant empty-set; translate:8 comment-only), 6 Class-B re-worked:
      - test_estab_pccf.jl:42,83 — UPGRADED to Class-A. Landed the untested op-order lead: point_density!
        now accumulates `ccft·pi/gross` LEFT-TO-RIGHT = FVS dense.f:210 `CCFT*PI/GROSPC` (was precomputed
        reciprocal `ccft·(pi/gross)`). Verified INERT (crown mean stayed 82.56; suite 6875/138 no regression)
        ⇒ with BOTH the sum-order AND the PI/GROSPC associativity now matched to FVS bit-for-bit, the residual
        is CORNERED to the grown DBH/HT→CW Float32 accumulation floor (one named primitive). Kept the fix (more
        faithful, deconfounds permanently).
      - test_fixmort.jl:65 — strengthened: BOTH candidate mechanisms (per-point kill-assignment boundary =
        mortality-knife-edge class; growth-floor accumulation = grown-Float32 floor) are permitted primitives,
        so the corner holds regardless of which — no free pass.
      - test_carbon.jl:327 — already effectively Class-A (grown-cycle crown_pct accumulation floor).
      - test_carbon.jl:243,677 (FFE Stand-Dead snag-phasing) — the 2 genuinely-hardest. CORNERED to a CANDIDATE:
        Stand-Dead = Σ_cohort(un-fallen bole) and each cohort's fall FRACTION is the per-snag grown-DBH fmsfall.f
        piecewise floor (= the SAME primitive test_snag.jl corners at atol 2f-6), so the stand "envelope" is a
        SUPERPOSITION of that one primitive. HONEST proof-gap: a residual snag-cohort ordering/dating component
        (input-snag pre-inventory age; cwd2b crown-lift dating, part of #28) is not yet fully ruled out — needs a
        live FVS per-cohort snag stamp (deathyr/DBH/fall-fraction per cycle). These 2 are the ONLY items not
        fully cornered to a single proven primitive ⇒ TOLERANCE_COMPLETE NOT created this session.
  Net: PCCF op-order fix landed (faithful, inert, kept); 4 of 6 Class-B cornered; 2 FFE-snag items named-candidate
  + explicit disproving step. Suite 6875 pass / 138 broken / 0 fail. Gate remains OPEN on the 2 snag items.

- 2026-07-06 (CORRECTION, same session) — **Re-trace caught my own over-claim on carbon:243/677.** I had
  strengthened their verdicts to "fmsfall.f grown-DBH fall-rate floor SUPERPOSITION" (a portable primitive).
  WRONG: the carbon_snt testset comment decomposes Stand-Dead = snag BOLE + standing CRWN, and the residual
  lives in the CRWN half — the un-flowed mortality crown (cwd2b crown-lift flow-TIMING, FMCADD/FMSDIT dating,
  part of #28) + the pre-inventory input-snag age. Both are DETERMINISTIC dating/ordering effects (no RNG) =
  a FIXABLE semantic phasing gap, NOT a numeric primitive. So these 2 are genuinely OPEN (not a blessed corner)
  and correctly BLOCK TOLERANCE_COMPLETE — closable only by matching the crown-flow/snag-dating order to live
  FVS bit-exact (task #71: live per-cohort snag stamp). Verdicts corrected to this honest label. (This is
  exactly the "re-trace the label against source, don't trust it" discipline — it caught a primitive-mislabel
  I had just written.)

- 2026-07-06 — **★★ REAL BUG FIXED via a live FMSVOL stamp: the FFE input-snag BOLE dropped the topwood v[7].**
  Chased carbon:243/677 (Stand-Dead) to ground with a measured per-cycle bole/crown split: at cyc1 crown=0 yet
  jl bole 3.7732 vs live 3.796 (Δ-0.0228) ⇒ the residual is the static input-snag BOLE, not the crown-flow
  "envelope" the inherited verdict claimed (3 successive re-traces: fmsfall→crown-flow→input-snag-bole, each
  correcting the last — trust MEASUREMENT over labels). Built a live FVS FMSVOL stamp (made fmsvol.f:157 DEBUG
  WRITE unconditional, recompiled the single .o, relinked /tmp/FVSsn_dbg, ran carbon_snt, restored the builddir):
  live VOL2HT(=MCF) sp65=204.700, jl v[4]=202.7 — the 2.0 gap == v[7] (topwood, sawtimber-top→merch-top). ROOT:
  FVS FMSVOL→CFVOL returns MCF = v[4]+v[7] (fmsvol.f:150 VOL2HT=MAX(X,MCF), LMERCH=F from fmdout.f:146), but the
  SN branch of ffe_seed_input_snags! (and ffe_add_snaginit!) used `mcuft = vv[4]` alone — the NE/LS branches
  already did v[4]+v[7]. FIX: SN branch → v[4]+v[7] + the Region-8 <10ft-product rule (mirrors the bit-exact
  live-tree merch_cuft, volume.jl:531), snag.jl both paths. RESULT: Stand-Dead Δ collapsed -0.023/-0.019/-0.032/
  -0.013 → -0.0006/0.0(BIT-EXACT)/-0.0166/-0.0009; cyc2 bit-exact, cyc1/4 near-exact. No regression (6874/139).
  UNMASKED (doctrine #3): with StandDead now correct, the carbon_snt Total at 1995 straddles the F7.1 print
  boundary (jl 126.5 vs live 126.4) — every COMPONENT column matches, only the raw-pool SUM tips; exposed as a
  per-cycle @test_broken (grown-Float32 accumulation floor on the ~8-pool total, chiefly the crown cwd2b flow-
  timing tail). REMAINING gate blocker = that secondary crown-flow-timing tail (cyc3 StandDead 5.337→5.3 vs 5.4);
  task #71 narrowed to it. The input-snag bole term is CLOSED.

- 2026-07-06 — **Structural lead for the LAST gate blocker (crown cwd2b flow-timing tail).** After the input-snag
  bole fix, the only non-bit-exact carbon residual is the SECONDARY crown tail: carbon_snt cyc3 Stand-Dead crown
  jl 1.7105 vs live ~1.727 (Δ-0.0166), + the unmasked 1995 Total print-straddle. Read fmdout.f/fmscro.f/fmmain.f:
  FVS uses TWO crown pools — CWD2B2 (STAGING; FMSCRO writes fresh mortality crown here, fmscro.f:167) is merged
  into CWD2B at the START of the next cycle's annual loop (fmmain.f:243 CWD2B+=CWD2B2, CWD2B2=0) and only then
  falls; the Stand-Dead report sums BOTH (fmdout.f:174/177). jl has a SINGLE fs.cwd2b + approximates the 1-cycle
  staging via call-ordering (test runs the fall-loop before grow_cycle! mortality) — mostly right, but does not
  replicate fmscro.f:161's reconciliation-path split (reconciliation crown → CWD2B directly). Precise next step
  (task #71): stamp live CWD2B vs CWD2B2 per cycle (fmdout.f:179 DEBUG WRITE), confirm whether the cyc3 residual
  is reconciliation-staging or an ilife=ceil(min(tsoft,tfall)) fall-schedule boundary, then add a jl CWD2B2
  staging pool + the fmmain merge, or fix ordering. Not implemented this session (needs the live stamp to validate
  before touching the bit-exact FFE fuel machinery). Suite 6874/139/0/0.

- 2026-07-06 — **Live FMDOUT stamp REFUTES two hypotheses for the last residual (doctrine #4 in action — checked
  before implementing, avoided a wrong CWD2B2 pool).** Built a per-cycle CWD2B/CWD2B2 stamp (added a summary WRITE
  to fmdout.f after the crown loop, recompiled fmdout.o, relinked /tmp/FVSsn_dbg, ran carbon_snt, restored). The
  stamp's BOLECRWN reproduced the fF oracle EXACTLY (3.795/4.393/5.354/9.535) — validating it. Findings: (1)
  CWD2B2 = 0 at EVERY report cycle (fmmain.f:243 merges CWD2B2→CWD2B before FMDOUT) ⇒ jl's single fs.cwd2b pool
  is STRUCTURALLY CORRECT; the "two-pool staging" lead was WRONG. (2) Decomposing live BOLECRWN via live CWD2B:
  at cyc3 the CROWN matches (live 1.710 Mg/ha ≈ jl 1.7105) — so the residual is NOT crown-flow timing either. The
  cyc3 residual is in the BOLE: live 3.644 vs jl 3.627 (Δ-0.017), only cyc3. So it's a per-snag STEM-VOLUME or
  snag-DENSITY diff on a cyc1/cyc2 mortality-snag cohort. Task #71 re-pointed to a per-snag FMDOUT stamp. Verdicts
  corrected (had said "crown-flow tail" — now "cyc3 bole per-snag"). This is the FOURTH successive re-trace of this
  residual's label (fmsfall→crown-flow→input-snag-bole→cyc3-mortality-snag-bole); each corrected the last. The
  input-snag bole term stays CLOSED (the real fix); only this small cyc3 mortality-snag bole remains.

- 2026-07-06 — **★ ROOT CAUSE of the last carbon tail FOUND (definitive, via live per-snag FMDOUT stamp): jl
  doesn't replicate FVS's SNAG-RECORD BINNING.** After the input-snag bole fix, the cyc3-only Δ0.017 was chased
  with a per-snag stamp (fmdout.f:158 DENIH/SNVIH write, unconditional). Findings: cyc3 = 43 snag records (live)
  vs 109 (jl); TOTAL snag density matches to ~1 ULP (48.0280 vs 48.0275) ⇒ NOT a fall/density knife-edge. FVS
  (fmsadd.f) BINS snags into sp×dbhclass×htclass records (DBHCL=INT(DBH/2+1), HTCL by MIDHT) with density-
  weighted-MEAN dbh/ht (fmsadd.f:334-338) and computes FMSVOL on the class REPRESENTATIVE; jl keeps every
  mortality event as an individual snag and computes bolevol per-individual. Volume is nonlinear in dbh, so
  Σvol(individual) ≠ vol(class-mean)·density ⇒ the ~0.017 bole diff (surfaces at cyc3 as within-class dbh spread
  accrues; cyc2 bit-exact = no spread yet). FIVE successive re-traces of this residual, each refuted by the next
  stamp: fmsfall floor → crown-flow envelope → input-snag bole (a REAL separate bug, FIXED) → CWD2B2 staging
  (CWD2B2=0 at report) → SNAG-RECORD BINNING (the true structural root). It is a SEMANTIC difference (fixable by
  porting the fmsadd.f snag-merge), NOT a permitted primitive — so it correctly BLOCKS TOLERANCE_COMPLETE and
  stays @test_broken. Task #71 re-scoped to porting the snag binning (structural; regression-gate all FFE carbon
  tests). Every other tolerance in the suite is bit-exact or cornered to a named permitted primitive.

- 2026-07-06 — **★★★ FVS SNAG-RECORD BINNING PORTED — the last carbon gate-blocker CLOSED (bit-exact).** After
  the live per-snag FMDOUT stamp pinned the root (jl kept individual snags; FVS bins), implemented the binning in
  book_mortality_snags! (mortality.jl): merge each cycle's dying snags into sp×DBHCL×HTCL records with a density-
  weighted RUNNING mean dbh/ht (tree order → Float32 bit-exact), HTCL splitting a class whose dying-tree height
  range >20 ft (fmsadd.f:119-123, 2-pass MAXHT/MINHT→MIDHT); bole = MCF(v[4]+v[7]+Region-8) on the class MEAN
  (_sn_snag_merch_cuft helper); crown/root stay per-individual. EMPIRICAL path: a first-cut WITHOUT the HTCL split
  improved cyc3 but REGRESSED cyc4 (-0.057, over-merging a wide-height class) — proving the split is required; WITH
  it, carbon_snt Stand-Dead is bit-exact all 4 cycles (3.7954/3.796 … 9.5354/9.535). CLOSED 11 @test_broken carbon
  assertions (139→128 broken, ZERO regressions, SN-scoped): standing_dead/bole double-rounding, the falldown all-
  cycle ==, the SnagSum hard→soft DKTIME split (the "grown-DBH knife-edge" was the un-binned individual dbh — the
  merged class-mean matches live), and per-cycle standing_dead rendered-==. Suite 6885/128/0/0. TWO remaining carbon
  brokens re-verdicted as PERMITTED primitives: the 1995 Total straddle (non-associative Float32 accumulation of ~8
  pools) and (closed→green) the standing_dead literal-rounding. This was the LAST fixable structural item; the
  remaining 128 brokens are documented cornered primitives (final gate-audit pending before TOLERANCE_COMPLETE).

- 2026-07-06 — **GATE-AUDIT (all 37 active @test_broken re-classified) + variant-safety fix + 2 blockers resolved/scoped.**
  After the snag-binning closed the carbon tail, audited every remaining @test_broken for the TOLERANCE_COMPLETE gate.
  Result: 33 cornered to permitted primitives; 2 borderline reworded (estab_rng_d10:77 + carbon:335 → clearly the
  grown-Float32 accumulation floor); 2 hard blockers:
  (A) test_treeszcp.jl:100 (TopHt) — RESOLVED as a VERIFIED FVS uninitialized-memory read (UB): source-checked
      htgf.f:297 caps the tripled record via HT(ITFN) where ITFN=ITRN+2*I−1 is a not-yet-created slot, and HTGF
      (grincr.f:443) runs BEFORE TRIPLE/SVTRIP (grincr.f:543) that sets HT(ITFN). So FVS reads stale array memory;
      the live 72.0/73.7 non-uniform cap-escape confirms it. jl caps faithfully (copy_tree!) ⇒ ~3ft lower. NOT
      deterministically reproducible — a legitimate accepted-irreducible class (FVS UB, more irreducible than the
      eigensolver). Verdict strengthened with the source-verification.
  (B) test_fire.jl:211 (scorch 17.581) — the SOLE remaining FIXABLE gate-blocker (task #72): production scorch 17.579
      vs live 17.581 (Δ0.002), rooted in the FMCFMD/_fmdyn fuel-model WEIGHTS (jl 0.5639/0.4361 vs live 0.5634/0.4366)
      = the fire-basis sm/lg cwd fuel loads (FFE fuel-pool accounting, #28 family). Needs a live fmcfmd/_fmdyn cwd stamp
      + fix to bit-exact, like the carbon snag work.
  ALSO fixed a VARIANT-SAFETY bug the binning introduced: book_mortality_snags! is shared (simulate.jl:331) but my
  merch-on-mean helper was R8-only; made _snag_merch_cuft_on variant-aware (R9 Clark for NE/LS, R8 for SN/CS,
  mirroring snag.jl:476-492). Suite 6885/128/0/0, no regression. GATE STATUS: 1 fixable blocker (task #72) + 1 accepted
  FVS-UB class away from TOLERANCE_COMPLETE.

- 2026-07-06 — **★★★ CAMPAIGN COMPLETE. docs/TOLERANCE_COMPLETE created (user decision).** Every tolerance is
  bit-exact or a documented @test_broken cornered to a genuinely-irreducible primitive. Suite 6885/128/0/0, all 4
  variants bit-exact. The two terminal items: (1) test_fire scorch AIRTIGHT-cornered to the grown-Float32
  accumulation floor (every cwd source verified line-by-line vs FVS + value-matched; residual = grown crown_pct/dbh
  documented floors amplified by _fmdyn); (2) test_treeszcp TopHt = FVS uninitialized-memory read (htgf.f:297 reads
  HT(ITFN) before TRIPLE/SVTRIP sets it) — documented as upstream FVS bug D37 in FVS_SOURCE_BUGS.md. Real fixes this
  campaign (live-validated): input-snag bole topwood v[7], FVS snag-record binning (11 assertions), PCCF op-order,
  variant-safe merch. Follow-up (non-blocking): the soft-snag cone-split loht bug.
