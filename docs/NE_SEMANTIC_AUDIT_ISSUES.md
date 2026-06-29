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

## AGENT FINDINGS (7 of 8 reported; volume pending)

### 1. DG + BAL (dgf.f/badist.f/balmod.f)
- FIXED: **ATTEN=0→1000** (DGCONS dgf.f:195) — feeds calib σ→tripling spread. ⇒ stand-1 BA 193→**194 = live**! Real WP-tail contributor.
- FIXED: **DIAGR≥.0001 floor** (dgf.f:159, always-active) + **removed the SN −9.21 clamp** (NE dgf.f:169 has none; floor guarantees DDS>0) + **COR2 multiply** (dgf.f:158, LDCOR2-gated).
- FALSE-POSITIVE (mine): SDIDEF PMSDIU — keyword_dispatch.jl:469 already /100s it; jl fraction = FVS PMSDIU/100. Faithful.
- LOW/inert: SMCON not zeroed (jl has no smcon field), BADIST d≤0 skip (DBH=0 edge), bark ULP. Faithful: B1/B2/B3 tables, POTBAG, BADIST/BALMOD, calib_dbh data-flow.

### 2. DG calibration + ARMA (dgdriv.f)
- **CONFIRMED — DG_PSIGSQ wrong for NE**: jl hardcodes SN 0.089827273; ne/dgdriv.f:95 is **0.0898** (108*). Affects the empirical-Bayes COR shrink (WC/TEMP) for CALIBRATED species (WP). NET01-AFFECTING. TO FIX (variant-aware). [HIGH]
- gated/edge: DGSD≥1 guard on calib OLDRN seeding (inert at DGSD=2); first-cycle AUTCOR old=YR vs measurement-FINT (GROWTH-keyword only); no-room tripling deterministic-vs-stochastic (capacity edge); cormlt_h split (assumption, likely faithful); fnmin kwarg. Verified faithful: AUTCOR/ssigma/corr/rhocp/FU·FL·FM/BJRHO/COR-shrink/VARDG/DGBND/gradd/backdating/OLDRN order.

### 3. Height growth (htgf.f/htcalc.f) — HIGHLY FAITHFUL
- LTBHEC (127×6) + MAPNE bit-exact; site curve/age/incr/HTMAX/GMOD·0.8/OLDRN all match. Only gated: HCOR2 not ported (height calib), base≤0 guard (benign).

### 4. REGENT small-tree (regent.f)
- **CONFIRMED — tripling diameter-blend uses CENTRAL dgk for all 3 records** (small_tree_growth.jl:71,100) instead of per-triple U/L (regent.f:373 uses DG(ITRIPU/ITRIPL)). stash.dgU/dgL available. Affects tripled 1.5–5" trees (net01 has them). NET01-AFFECTING. TO FIX. [MED-HIGH]
- gated (FINT≠10): diameter-blend double-conversion at FINT≠10; DK<0 fallback extra blend; hk≤4.5 recorded-DG; RHCON=1. Faithful: SCALE/SCALE2/XMIN/XMAX/no-0.8/xwt/±10% RNG/size-cap.

### 5. Mortality (morts.f) — NE hooks FAITHFUL
- All NE-specific verified (PMSC/PMD, RI 0.5 halving, YR=10, DBHSDI, VARMRT EFFTR, BAMAX). Shared-driver gated: sdimax<5 whole-stand-kill passes bg_tokill not full (climate-only; SN-shared); BAMAX-keyword ignored (gated); BADEAD>0 guard (pathological). NOTE: these are SHARED (SN too) — masked-bug candidates.

### 6. Crown ratio (crown.f/cwcalc.f) — FAITHFUL
- CRNEW + BCR1-4 (4×108) bit-exact; BA/gross_space faithful; crown-width library faithful. Only gated CRNMULT corner cases (window-vs-raw guard, ICFLG one-shot, DHI default 99 vs 99999) — inert without a windowed CRNMULT keyword.

### 7. Establishment (estab.f/essubh.f/esgent.f)
- **CONFIRMED-MISSING — WK4=HTIMLT per-period scaling** (esgent.f:52 HTG=HTG*WK4; HTIMLT=FTEMP/(GENTIM+0.0001), estab.f:520). jl omits it entirely. The agent's value (0.40 from min(TRAGE,GENTIM)) is a MISREAD — FTEMP=BAA (estab.f:451, floored 1.0), GENTIM=FINT−5. The BARE-stand WK4 value needs the BAA basis at estab.f:520 determined before applying (likely the BARE first-cycle deficit cause). [HIGH, value-TBD] + the coupled esgent.f:54-62 WK4<1 DBH re-derivation.
- CONFIRMED-MISSING: establishment.jl phase-2 lacks the **hk≤4.5 → DBH=0.1+0.001·hk** guard (regent.f:290-293) — the Wykoff inverse log(hk−4.5) NaNs for hk≤4.5 (reachable: NE base heights 2.5–4 ft). [MED-HIGH]
- gated/RNG: pre-replicate ESRANN/WK6 draw seed-shift (masked by XMIN floor); FINT≤5 LSKIPH extra draw; CON·HGADJ·XRHGRO omitted (HCOR/keyword); HTADJ add (=0). Verified faithful: ESSUBH formula+REFAGE table, reject bounds (my fix), crown draw, age, TRAGE clamp.
