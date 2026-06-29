# NE Semantic-Faithfulness Audit — Issues List

Campaign of parallel agents auditing each ported NE chunk by **reading the code semantics** vs the FVS
Fortran source (NOT runtime tracing), upstream-first. Each suspected divergence is logged here, then
verified against the source and fixed. Per directive: **a semantic-certain fix that regresses other tests
(even SN tests) is a masked-bug signal — note it, do not revert; it may be a hidden bug elsewhere.**

Status key: `SUSPECTED` (agent-reported, unverified) · `CONFIRMED` (verified vs source) · `FIXED` ·
`FALSE-POSITIVE` (faithful on closer read) · `DEFERRED`.

## Already found/fixed pre-campaign (this session)
- FIXED: establishment planted-height random reject bound `[0,1.5]`→`[-2.5,2.5]` (estab.f:489).
- VERIFIED FAITHFUL (read both sides): establishment relht=0 (esgent no DENSE recompute / BARE AVH=0),
  scale_e=FNT/REGYR=0.5, ne_htcalc_incr=HTCALC HTGP5−HTG0, ne_htcalc_age=HTCALC mode0, ESSUBH height.

## Agent findings (filled in as the 8 agents report)
### 1. DG + BAL competition (dgf.f/badist.f/balmod.f)
### 2. DG calibration + ARMA + tripling (dgdriv.f/dgbnd.f)
### 3. Height growth + HTCALC (htcalc.f/htgf.f)
### 4. REGENT small-tree (regent.f)
### 5. Mortality (morts.f)
### 6. Crown ratio + crown width (cratet.f/cwcalc.f)
### 7. Establishment (estab.f/essubh.f/esgent.f + regent LESTB)
### 8. R9 volume (r9clark*.f/mrules.f/r9logs)
