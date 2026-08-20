# OC/OP REGENT small-tree port — turnkey spec (2026-08-20)

Root cause (measured, docs/OC_VARIANT_PORT_AUDIT.md + memory fvsjl-oc-organon-blm-volume): OC/OP route ALL
IORG=0 trees to the large-tree `oc_htgf_native` (organon_hook.jl:86), but FVS routes trees with **DBH < XMAX
(4.0", sp23/50=10.0)** to `oc/regent.f` (small-tree height-age model). `small_tree_growth!(::OregonCoast/
::Olympic)` is a NO-OP (organon_hook.jl:132). Symptom: sub-4.5ft seedlings over-grow ~7ft/cycle (tn2 DF 2.0→
jl 10.7 vs oracle 3.8) — because oc_htgf_native forces `relht=1` when PCCF<100 (htgf.f:230) removing the
competition suppression, while smhtgf KEEPS it. Over-tall seedlings perturb CCH → large trees ~0.15ft short by
cyc2 → the ~1% TCF/CCF .sum drift (BA/TPA bit-exact).

## KEY SIMPLIFICATION: OC/OP are DGSD=0 ⇒ the whole small-tree path is DETERMINISTIC
regent.f:234 `IF(DGSD.GE.1.0) ZZRAN=BACHLO(...)` — DGSD=0 ⇒ ZZRAN=0, NO RNG in the height increment. And the
BACHLO crown draw (regent.f:~193) is ONLY in the LESTB (new-regen) branch, not for inventory small trees. So
inventory small-tree height+DBH is fully deterministic ⇒ **bit-exact achievable** (not just cornered).

## The model (inventory small tree, NOT LESTB), per regent.f + smhtgf.f
For each tree with D < XMAX(ISPC):
1. `HTGRR = SMHTGF(sp, D, H, CR, BA, BAL, SI, RELHT)` — smhtgf.f, 5 equations by MAPSP(sp):
   - MAPSP: 1=pines 2=firs(+DF) 3=blackoak 4=tanoak 5=redwood. FACTOR = SPADJF(sp)·(0.80+0.004·(SI−50)).
     TEMBAL=max(BAL,5).
   - pines(1):    HTGR = exp(0.7452 −0.003271·BAL −0.1632·CR +0.0217·CR² +0.00536·SI)·FACTOR·1.75  [CR here = ICR/10, i.e. 0..10 scale — NOT divided again]
   - firs(2):     DOMHTGR = 5·(2.2227+0.4314·SI)/(29.0−0.05·SI); CRf=CR/10 (→fraction); CRMOD=1−exp(−4.26558·CRf);
                  RHMOD=exp(2.54119·(RELHT^0.250537−1)); SMHMOD=1.016605·CRMOD·RHMOD; HTGR=DOMHTGR·SMHMOD
   - blackoak(3): HTGR = exp(3.817 −0.7829·ln(TEMBAL))·FACTOR
   - tanoak(4):   HTGR = exp(3.385 −0.5898·ln(TEMBAL))·FACTOR
   - redwood(5):  height-age inversion (HTMAX=2.242202·SI; if HTMAX−H≤1 →0 else AGE1 from H, AGE2=AGE1+5, HTGR=H2−H1)
   - floor HTGR≥0.1.
   ⚠ CR CONVENTION: regent.f sets `CR=REAL(ICR(I))/10.` BEFORE calling smhtgf (so CR∈0..10). smhtgf firs branch
     divides AGAIN (`CR=CR/10`) → fraction. pines/oak branches use CR∈0..10 directly. MATCH exactly.
   ⚠ RELHT: regent.f `RELHT=H/AVH; IF(RELHT>1.05) RELHT=1.05` (NO pccf<100→1 override — that's the whole bug).
2. `CON = RHCON(sp)·exp(HCOR(sp))`. RHCON default 1.0 (regent.f:573; =RCOR2 if small-tree-height calibrated).
   **HCOR = the height-calibration correction (regent LSTART, ~line 40+/560+) — THE TRAP (cf IE/LS HCOR fixes).**
   Must MEASURE whether ocmin/oct01 calibrates: if no measured-height small trees, HCOR=0 ⇒ CON=1.0 (simple).
   VERIFY via scoped DEBUG REGENT (⚠ generic DEBUG segfaults FVSoc_clean — use scoped, or g16 dump-replay).
3. `HTGR = (HTGR + ZZRAN·0.1)·XRHGRO·SCALE`. ZZRAN=0 (DGSD=0). XRHGRO=XRHMLT(sp) height mult (1.0 default, no
   HTGMULT kw). SCALE=FNT/REGYR, REGYR=5. FNT=FINT (=5 for these stands; LESTB adjusts, N/A for inventory).
4. XWT blend: `XWT=(D−XMN)/(XMX−XMN); IF(D≤XMN.OR.LESTB) XWT=0`. XMN=2.0, XMX=4.0 (sp23/50: 10.0).
   `HTG = HTGR·(1−XWT) + XWT·HTG_largetree`. (sp23/50 have a special (HTGR+LTHG)/2 pre-average — see regent.f:255.)
   So D≤2" ⇒ pure small-tree; 2<D<4 ⇒ blend with oc_htgf_native; D≥4 ⇒ skipped (large-tree only).
5. Size cap: `IF(H+HTG > SIZCAP(sp,4)) HTG=SIZCAP(sp,4)−H` (floor 0.1).
6. DBH assignment (only for D < DGMIN(sp)=3.0, sp23/50=7.0): HK=H+HTG.
   - HK≤4.5: DG=0, DBH=D+0.001·HK  (TRIVIAL — covers tn2; no HTDBH needed).
   - HK>4.5: DBH from H-D function via HTDBH(IFOR,sp,DK,HK,1) [htdbh.f, 179 ln] — DK=(BX/(ln(HK−4.5)−AX))−1,
     BX=HT2(sp), AX=HT1(sp) if IABFLG=1 else AA(sp); DKK likewise from H (DKK=D if H≤4.5). Then DG via XDWT blend
     ((D−1.5)/1.5 clamped, or (D−XMN)/(DGMIN−XMN) for sp23/50). Uses inventory eqns if .NOT.LHTDRG or IABFLG=1.

## Coefficients to extract (block data / rcon.f / grinit.f)
SPADJF(50),MAPSP(50) — DONE (in smhtgf.f above). XMAX/XMIN/DGMIN/DIAM(50) — regent.f DATA lines 102-108 (DONE).
RHCON(50)=1.0 default. HCOR(50) — MEASURE (calibration). HT1/HT2/AA(50),IABFLG(50),LHTDRG(50) — for HTDBH
(htdbh.f + rcon.f loaders). SITEAR — jl already has (site setup). SIZCAP(sp,4) — max height, jl has. XRHMLT=1.

## Wiring
Replace organon_hook.jl:86 branch: for IORG=0 trees, if D < XMAX(sp) → new `oc_regent_smtree!` (steps 1-6),
else → oc_htgf_native (D≥4). Keep `small_tree_growth!(::OregonCoast/::Olympic)` a no-op (work done in the hook,
matching the ORGANON-drives-everything pattern). Same for Olympic (organon_hook_op.jl). OP MAXSP=39 has its own
species→CA-50 map for smhtgf — verify op/smhtgf.f MAPSP (may differ from oc).

## Validation (vehicle READY, scratchpad/oc/)
oracle: oc_or_trl.pkl (per-tree DBH/HT per cycle from FVSoc_clean TREELIST). jl: capture via write_sum_file
cycle_hook (see the harness that made oc_jl_trl.txt). TARGET: tn2/tn15 HT bit-exact at cyc1, then the whole
oct01 stand-1 .sum bit-exact-or-cornered multi-cycle (currently ~1% drift). Then re-confirm OP opt01 still
bit-exact-or-cornered (test_op_multicycle_sum) + multicycle 339/11. MIND THE TRAPS:
[[fvsjl-ls-regent-stalehtgr-fix]] [[fvsjl-ie-regent-hcor-calibration]] [[fvsjl-ci-regent-xwt-height-blend]].

## ⚠ HCOR CALIBRATION IS LIVE FOR THIS STAND — the trap must be ported (2026-08-20)
HCOR(sp) defaults 0.0 (regent.f:409) ⇒ CON=RHCON=1.0 — BUT regent.f:493 sets `HCOR(sp)=ALOG(CORNEW)` in the
LSTART calibration branch (regent.f ~400-506) whenever the stand has measured-height small trees. oct01
stand-1 HAS them (tn10 h20, tn12 h11, tn13 h13 are measured sub-4"-ish trees), so **HCOR is almost certainly
non-zero here** ⇒ a port that assumes CON=1.0 will NOT be bit-exact on tn2/tn15. The full port therefore
MUST include the regent LSTART small-tree height calibration (observed-vs-predicted small-tree height → CORNEW
→ HCOR=ln(CORNEW)), exactly the [[fvsjl-ie-regent-hcor-calibration]] / [[fvsjl-ls-regent-hcor-fix]] trap that
bit the IE/LS regent ports. MEASURE it: build a scoped-DEBUG REGENT (or g16 single-.o regent.f dump) to read
the live HCOR(sp) per species, then port the calibration to reproduce it before validating the increment.
This is why the port is a focused unit, not a quick height-equation drop-in. Order for the port session:
(1) instrument regent.f → dump live HCOR(sp) + per-tree HTGRR/HTGR for oct01 stand-1 cyc0;
(2) port smhtgf + the LSTART HCOR calibration to match HCOR bit-exact;
(3) port the increment wrapper (CON/XWT/SCALE) + DBH-from-height (HK≤4.5 trivial; HK>4.5 via HTDBH);
(4) wire into organon_hook.jl (route D<XMAX IORG=0), A/B tn2/tn15 vs oc_or_trl.pkl, then whole .sum.
