# ON sub-12cm REGENT small-tree growth — port handoff (measurement foundation READY)

## Measurement foundation (BUILT this turn)
- Validation stand: `scratchpad/on/ont_sm.{key,tre}` — 8 trees DBH 4-11 cm (all < XMAX=12 cm ⇒ every
  record hits the REGENT branch). Oracle `.sum` = `scratchpad/on/ont_sm.sum` (cyc0 TPA 296519/QMD 7.8/FORTYP 122).
- **Per-tree REGENT dump oracle**: `/workspace/.onwork/FVSon_regdump` (regent.f patched to dump, per sub-XMX
  tree at stmt 23: `I ISPC K  HT(I) D HTG(K) DG(K)` as Z8.8 hex → fort.773). Source: `scratchpad/on/regent_regdump.f`.
  Instrumented `.sum` VERIFIED byte-identical to clean (dump inert). Run:
  `cd /workspace/.onwork/run_sm && printf "ont_sm.key\nont_sm.tre\n" | ../FVSon_regdump` → fort.773.
  Each tree tripled: base K=I + 2 stash records K=ITRN+2*I-2+L (L=1,2). Most tiny trees floor DG at 0.1
  (=3DCCCCCD); HTG varies via the NC128/ONSTHG blend.

## Algorithm (canada/on/regent.f, D < XMX=4.72" branch; ON LSTART=false live path)
Per species ISPC (constants: CON=RHCON(ISPC)*exp(HCOR(ISPC)); SCALE=FNT/REGYR; SCALE2=YR/FNT;
DGMX=DGMAX(ISPC)*SCALE; XMX=4.72; XMN=3.15), per tree D<XMX, per triple K:
1. **HTMAX gate**: HTCALC MODE0=0 → HTMAX (=_on_htcalc_htmax, ALREADY PORTED). If HTMAX-H ≤ 1 → HTGR=0.10.
2. **else HTGR**: HTCALC MODE0=9 = `BAL=BA*(100-PCT)*0.01; HTG1=ONSTHG(ISPC,D,H,SI,BAL)*YRS(=10)`.
   Then `HTGR = HTG1 * CON * SCALE * HGADJ(ISPC) * XRHGRO`.  (XRHGRO=XRDMLT user mult, default 1.)
3. **BALMOD**: GMOD ≡ 1.0 for ON (balmod.f). RELHTA=min(HT/AVH,1); GMOD=1-((1-GMOD)*(1-RELHTA))=1 ⇒ HTGR unchanged. Floor 0.1.
4. **XWT blend**: XWT=(D-XMN)/(XMX-XMN); if D≤XMN → XWT=0. `HTGR = HTGR*(1-XWT) + XWT*HTG(K)` (HTG(K)=large-tree htont HTG). Floor 0.1.
5. **random**: if DGSD≥1 (ON=2.0): RAN=BACHLO(0,1,RANN) in [-1,1]; `HTGR += RAN*0.1*HTGR`. (★ port BACHLO, never FFI; RNG stash = same as htont/dgf.)
6. **SIZCAP**: if H+HTG(K) > SIZCAP(ISPC,4) → HTG(K)=SIZCAP-H (floor 0.1). Else HTG(K)=HTGR (floor 0.1).
7. **DBH incr** (stmt 4): HK=H+HTG(K). If HK≤4.5 → DG=0, DBH=D+0.001*HK. Else Wykoff:
   - BX=HT2(ISPC); AX = IABFLG==1 ? HT1(ISPC) : AA(ISPC). DKK=(BX/(log(HK-4.5)-AX))-1. DK = H≤4.5 ? D : (BX/(log(H-4.5)-AX))-1.
   - ON LHTDRG=.FALSE. all ⇒ CALL HTDBH(IFOR,ISPC,DKK,HK,1) then (H>4.5) HTDBH(IFOR,ISPC,DK,H,1). (HTDBH already shared — feed ON HT1/HT2/SNDBAL.)
   - BARK=on_bratio(ISPC,D,H). If DK<0 or DKK<0: DG=HTG(K)*0.2*BARK*XRDGRO; else DGSM=(DKK-DK)*BARK*XRDGRO (≥0); DDS=DGSM*(2*BARK*D+DGSM)*SCALE2; DGSM=sqrt((D*BARK)^2+DDS)-BARK*D (≥0); DGGR=DGSM*(1-XWT)+XWT*DG(K); DG=max(DGGR,0.1).
   - DG floor 0.1, cap DGMX. If DBH+DG < DIAM(ISPC) → DG=DIAM(ISPC)-DBH.
8. **DGBND**(ISPC,DBH(K),DG(K)) size cap (port dgbnd.f — small).
9. Tripling: L=0,1,2 → K=I then ITRN+2*I-2+L.

## ONSTHG (canada/on/onsthg.f — PORT verbatim, Penner small-tree annual HT growth, ft):
OSPMAP(72)→KSP(1-28); B00/B01/B02/BSI/BBAL/B95(28). HM=max(0.05,HT*0.3048); LHM=log(HM); SIM=SI*0.3048;
BALM=BAL*0.2295643(=ON_FT2pACRtoM2pHA). HTG=B00+B01*LHM+B02*HM+BSI*SIM+BBAL*BALM; clamp[-5,5]; exp;
clamp[0.0001,B95(KSP)]; *3.28084(MtoFT). Transcendentals via glibc logf/expf ccall (Float32 bit-exact).
(Full coeff arrays transcribed in onsthg.f — copy the 6 DATA blocks + OSPMAP.)

## Coefficient sources to wire (canada/on block data / grinit.f):
RHCON, HGADJ, DGMAX, HT1, HT2, AA, IABFLG, HCOR(calib, =0 non-calibrated), SIZCAP(have), DIAM(min dbh),
XRHMLT/XRDMLT (MULTS, default 1). Grep `DATA RHCON|HGADJ|DGMAX|HT1|HT2|IABFLG` in canada/on/*.f + grinit.f.

## Validation plan
Dump-replay: port full chain → run jl on ont_sm → compare per-tree HTG(K)/DG(K) hex to fort.773 base records
(the RNG ±0.1 makes it stochastic; inject the oracle's BACHLO draws or match the RNG stream to get bit-exact,
same as htont/mortality VARMRT). Then end-to-end ont_sm .sum cyc0 (TPA/QMD/TopHt). Multicycle 339/11 must hold
(ON-gated). This is a LARGE multi-routine chunk — ONSTHG + HTCALC-mode9 + REGENT blend + DGBND + coeff tables.

## FINAL de-risk (all coefficients sourced — implementation is a mechanical SN-template clone)
The SN `src/variants/southern/small_tree_growth.jl` is the template — ON's regent.f is the SAME
base/regent.f. Clone it into `small_tree_growth!(::Ontario)` swapping:
- XMIN=3.15, XMAX=4.72 (have ON_REG_XMIN/XMAX), REGYR=10.0.
- **Simplifications (all sourced): con=1 (RHCON=1 regent.f:624 non-calibrated, HCOR=0), HGADJ=1
  (regent.f:75 MAXSP*1), XRHGRO=XRDGRO=1 (MULTS default), SCALE=FNT/REGYR=10/10=1, SCALE2=YR/FNT=1,
  DGMX=DGMAX*SCALE=5.0 (regent.f:76), GMOD=1 (balmod.f ON).** So HTGR = ONSTHG*10, then XWT-blend, random, SIZCAP.
- Height increment: **ONSTHG** (onsthg.f — port the 6 DATA arrays B00/B01/B02/BSI/BBAL/B95[28] + OSPMAP[72];
  HM=max(.05,HT*.3048),LHM=log,SIM=SI*.3048,BALM=BAL*.2295643; HTG=Σ; clamp[-5,5];exp;clamp[.0001,B95];*3.28084;
  BAL=BA*(100-PCT)*.01; *YRS=10; glibc logf/expf ccall). Replaces SN's htcalc_incr.
- htmax gate: `_on_htcalc_htmax` (HAVE, volume.jl). If htmax-h≤1 → HTGR=0.1.
- HTDBH (H→D) `_on_htdbh_dbh(sp,h)`: IWYKCA(sp)==0 → `(HT2/(log(H-4.5)-HT1))-1`; else Curtis-Arney
  (SNALL P2/P3/P4, htdbh.f:342-353); then `D=max(D,SNDBAL(sp))`. Coeffs: HT1 blkdat.f:276, HT2 blkdat.f:290,
  IWYKCA htdbh.f:306, SNALL htdbh.f (P2/P3/P4), SNDBAL htdbh.f:293 (all 72; SNALL only for IWYKCA=1 spp).
- bark: on_bratio(sp,d,h). DGBND: `dg=min(dg, 6*exp(-0.03*min(dbh,150)))`, floor 0 (dgbnd.f).
- DIAM floor (rgntsw.f:90): `DIAM=[0.3,0.3,0.2, 0.3×11, 0.2×54, 0.3×4]` (72). `if dbh+dg<DIAM → dg=DIAM-dbh`.
- DGSM chain identical to SN `_regent_dg` (DDS=DGSM*(2*bark*d+DGSM)*scale2; DGSM=sqrt((d*bark)^2+DDS)-bark*d).
Then wire small_tree_growth!(::Ontario) (drop the error). Validate: jl per-tree HTG(K)/DG(K) hex vs fort.773
base records (random ±0.1 is RNG-stream — deterministic part bit-exact, random = #206 straddle) + ont_sm .sum cyc0.
Every coefficient has a source line above — NO remaining unknowns.
