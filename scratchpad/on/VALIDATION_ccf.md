# ON CCF — open-grown crown width (cwcalc.f IWHO=1) → ccfcal.f, bit-exact vs FVSon_g16

## Oracle instrument
- `scratchpad/on/cwcalc_ccfdump.f` = canada/on/cwcalc.f patched to dump `ISPC IWHO D CW HILAT HILONG HIELEV HI`
  (hex TRANSFER→Z8.8) when IWHO==1, before RETURN. Build: gfortran-16 -std=legacy -w -fno-automatic
  -finit-local-zero into /workspace/.onwork/instr_cw/cwcalc.o; relink g16obj minus cwcalc.o + this + shim
  → /workspace/.onwork/FVSon_cwdump. Run `echo ont01.key | FVSon_cwdump` in a dir with ont01.key+.tre → fort.772.

## Result — per-tree open-grown CW 8/8 bit-exact (first-cycle dump)
| ISPC | code | US(JSP2) | CWEQ | CW hex (oracle==jl) |
|------|------|----------|------|---------------------|
| 1  | PJ | JP | 10503 | 41C29F7E |
| 5  | PW | WP | 12903 | 41C47486 |
| 6  | SW | WS | 09403 | 4192C3CD |
| 8  | BF | BF | 01203 | 416C9717 |
| 9  | SB | BS | 09503 | 413E9300 |
| 11 | CE | WC | 24101 | 41594FAA |
| 26 | MH | SM | 31803 | 41F54340 |
| 28 | BE | AB | 53101 | 41EB0455 |  ← HI-dependent; nails oracle HILAT/HILONG/HIELEV

## HI inputs (oracle, from the dump)
HILAT=0x423B1EB8=46.78, HILONG=0xC2B83852=-92.11, HIELEV=0x44761020=984.252 ft (=300 m·3.28084),
HI=0x4227AE18=41.920. ⇒ forkod.f US-Superior (KODFOR 915/916) default lat/long (ont01 location lost to
the initre LNOTBK overflow) + metric STDINFO elevation (300 m → 9.8425 hundred-ft). Wired in kw_stdinfo!
(ON metric-elev branch mirroring BC; ON default lat/long).

## .sum
stand_ccf=20844 ⇒ overflows the sumout.f I4 CCF field ⇒ `****` = oracle. Added `_fi` (Fortran Iw edit,
`*×w` on overflow) to write_sum_row — byte-identical to %wd for values that fit (multicycle 339/11 unchanged).
cyc0 .sum row now bit-exact through CCF/AT_CCF; only the trailing FORTYP/size/stock (999/55 vs 125/11) remains.
