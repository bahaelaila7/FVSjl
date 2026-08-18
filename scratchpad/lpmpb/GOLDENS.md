# LPMPB port — oracle + validated goldens (Mountain Pine Beetle, lpmpb/*.f)

## Relink oracle: FVSie_lpmpb  (build_ie_lpmpb.sh)
- IE object set (FVSie_buildDir/*.o), `base/exmpb.f` NO-OP stub DROPPED, real
  `lpmpb/*.o` added (IE block data `mpblkdie.o` only; generic + other-variant
  `mpblkd*.o` excluded to avoid duplicate BLOCK DATA), plus `wsbwe/txnote.o`
  (the ??X.exe warning-flag setter referenced by MPBIN) and `base/exppe.o` stub.
- Build: `bash build_ie_lpmpb.sh [ABS_ovr_dir] [out_binary]`  (gfortran-16
  -std=legacy -w -fno-automatic -O0 + /workspace/.crwork/isoc23_shim.o).
- LPMPB-OFF lp_off.key runs to normal STOP 10 on the pure-lodgepole IE stand.
- LPMPB-ON lp_on.key (`MPB` / `MPBSTART 1` / `DEBUG` / `END`) fires the outbreak in
  cycle 1: lodgepole 113→38 TPA, mortality 152 vs 20 — the .sum-DELTA is large.
- Invoke: `./FVSie_lpmpb --keywordfile=lp_on.key`  (companion lp_on.tre auto-opened).

## Instrumented relink for goldens (g16 single-.o swap)
- `ovr/colmrt.f` = pristine `lpmpb/colmrt.f` + stderr `WRITE(0,…) TRANSFER(x,INDEX)`
  Float32-hex dumps (DBGCM_HDR / DBGCM_CLS START/GREEN/PRKILL / DBGCM_REC
  J/IDX/DBH/PROB/XT/WK2). Build: `bash build_ie_lpmpb.sh ABS/ovr FVSie_lpmpb_dbg`.
- Instrumented .sum data rows BYTE-IDENTICAL to the clean relink (verified
  a.txt==b.txt, timestamp header line excluded) ⇒ dumps are trustworthy.
- Dump captured to `lp_dbg.stderr`.

## RNG identity (mprann.f) — MEASURED, not FFI'd
- MPRANN is the same double-precision MINSTD LCG as DFBRAN, but seeded DIRECTLY to
  55329 by MPBINT (`MPRNSD(LSET=.TRUE.,55329.0)`) — NO DFBSCH `+1128` reset:
  `S1 = DMOD(16807·S0, 2147483647)`, `SEL = REAL(S1/2147483648)`, `S0 = S1`.
  Entries MPRNSD (reseed, forces odd), MPRNGT (get S0), MPRNPT (put S0).
- Only consumer in the whole model: MPBGO's RANSTART branch (`CALL MPRANN(X)`;
  `IF (X.LT.PROTBK)`). MPBMOD/COLMOD draw NO random. ⇒ the entire mortality
  computation is DETERMINISTIC on the MANUAL/MPBSTART path.

## VALIDATED SLICE — deterministic Cole "rate of loss" core (BIT-EXACT, Float32)
Replaying lp_dbg.stderr through src `lpmpb.jl` (check_lpmpb.jl):
  COLDBH START      10/10
  COLMOD GREEN(NUMYRS=10)  10/10   ← includes the full 9-iteration PRNOIN^DEAD pow chain
  COLMRT PRKILL     10/10
  COLIND            18/18
  XT = PRKILL(colind(DBH))·PROB   18/18
  MPBER noer=true (host clears); mpb_idxlp IE=7, CR=11
NOTE: Julia Float32 `^` matches gfortran `REAL**REAL` bit-exact here (no ULP
straddle over the 9-year chain) — measured, re-check on other stands.

## Golden values (class 4..9 nonzero; NUMYRS=10)
START  4=413758B4 5=42031AA4 6=41B276C8 7=41814D1E 8=4143FA1F 9=40226927
GREEN  4=413464DE 5=41C3EA00 6=406A3223 7=3E065EB0 8=36C6E3D7 9=3DCF7CEC
PRKILL 4=3C83EB13 5=3E8172ED 6=3F5601C1 7=3F7DEBEF 8=3F7FFFF8 9=3F75C794
(per-record DBH/PROB/IDX/XT embedded in test_lpmpb.jl)

## MPBINT constants (mpbint.f) — model defaults (overridable by INITMORT/QVALUES)
ZINMOR(XZIN) = 0,0,.0038,.0128,.0206,.0353,.0500,.1429,.1500,.1500
PRNOIN(XPRN) = 1,1,.9935,.982,.965,.909,.743,.309,.285,.285
MPMXYR=10 NCLASS=10 IBOUSE=0 PRBSCL=1 EPIPRB=0.5 FORLAT=44 seed=55329 LPOPDY=F
