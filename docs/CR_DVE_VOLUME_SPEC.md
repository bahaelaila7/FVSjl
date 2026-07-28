# CR volume — DVE method (r3d2hv.f) implementation spec — TURNKEY

Port target: `volume/NVEL/r3d2hv.f` (R3 D2H). Covers 31/38 CR species. Entry point for chunk 8.

## Model (per tree)
`D2H = DBHOB² · HTTOT`. Two species families keyed by `VOLEQU(8:10)` (the 3-digit species code):

### Timber species (093,106,113,122,202,310,314): full cubic + board
Each block (see r3d2hv.f):  ENTIRE = a + b·D2H (total stem cuft);
  UM6 = c + d·((6³·HTTOT)/DBHOB^1.5) + e·DBHOB²;  GCUFT6 = max(ENTIRE−UM6, 0);   (gross cuft to 6″)
  UM4 = c + d·((4³·HTTOT)/DBHOB^1.5) + e·DBHOB²;  GCUFT4 = ENTIRE−UM4;   TWVOL = UM6−UM4;   (topwood 4-6″)
  INTBDFT = GCUFT6·(poly in DBHOB^−1,−2,−3);  SCBDFT = INTBDFT·(poly in DBHOB^−1,−2).
  Coeffs read so far — 093: ENTIRE 0.225466+.002170·D2H, UM* −0.2664752+.006129·(k³H/D^1.5)+.007431·D², INTBDFT
  ·(5.987363−9.847918/D+300.812808/D²−...), SCBDFT·(.878454−15.998458/D²). 113: 0.160889+.002032·D2H; UM
  −0.213005+.004912·(..)+.006061·D². 746(aspen): 0.0327+.002311·D2H; UM −0.236432+.005802·(..)+.006080·D².
  122(special, region-300): SCBDFT/GCUFT6 have their own D2H-breakpoint forms (see r3d2hv.f:67-85).

### Woodland species (015,060,800,999,746-var): DRC-based, GCUFT4 only, NO board
  IF DBHOB>3 OR DRC>3: if DRC>0 D2H=DRC²·HTTOT; D2HA=D2H/1000; branch on FCLASS (1=single else multi) and a
  D2HA breakpoint → GCUFT4 = quadratic (D2HA≤bp) or a+b·D2HA−c/D2HA (>bp). VOLEQU(2:3)=='01' vs '02' sub-variant
  for 800. 999 = default. (060 also has VOLEQ(2:3)=='01'/'02'.) ENTIRE/board = 0 for woodland.

## VOL(15) fill (r3d2hv.f:505-546) — the output mapping
DBHOB<9 ⇒ SCBDFT=INTBDFT=0.  VOL(1)=ENTIRE.
UNT==1 (prod "01", sawtimber): VOL(2)=SCBDFT, VOL(10)=INTBDFT, VOL(4)=GCUFT6, VOL(7)=TWVOL.
UNT==3 (prod "02"): VOL(4)=GCUFT4, VOL(6)=VOL(4)/79.  VOL(15)=UM4.  Clamp VOL(1,2,4,6,7,10,15)≥0.

## jl wiring (volume.jl compute_volumes! ~line 561)
Currently calls `_R8CLARK_VOL(veq[sp],…)` UNCONDITIONALLY. Add dispatch on the method substring:
  m = veq[sp] startswith "NVB" ? :nvb : SubString(veq[sp],4,6)  # "DVE"/"NVB"/"FW2" else Clark
  :DVE → new `_cr_dve_vol(veq[sp], d, h, prod)` returning a 15-vec; read tcf=v[1], scf/mcf=v[4](+v[7]), bf=v[10]
  exactly like the Clark path. prod "01"→UNT1, "02"→UNT3. (HTTOT=height; DRC from t.dbh for woodland or a DRC col.)

## Validation
Per-tree on sp23 (code 800, woodland) in the San Juan cr_calib stand: instrument live r3d2hv.f WRITE(16) VOL(1),
VOL(4), VOL(7), VOL(10) per sp23 tree; compare jl _cr_dve_vol. (Whole-stand .sum needs FW2+NVB too — those species
dominate San Juan; DVE validates per-tree here.) Then a constructed WF/ES timber-DVE stand for the timber blocks.

## Order (validation-driven, from CR_VARIANT_PORT_AUDIT.md)
FW2 (fwinit.f, DF/PP sp3/13) + NVB (nsvb.f, AS sp20) FIRST — they dominate San Juan ⇒ whole-.sum validatable
(live cr_calib TCuFt=4049 MCuFt=3256 BdFt=13487). THEN DVE (this spec, widest coverage).
