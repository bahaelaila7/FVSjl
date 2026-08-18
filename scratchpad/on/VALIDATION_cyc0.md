# ON (Ontario) cyc0 growth .sum — DESIGN plot-count/GROSPC + QMD + merch unblock

Staged: `scratchpad/on/wiring.diff` (git apply-able on kt-variant-port d1fc55a3).
Oracle: FVSon_g16 (== production-object FVSon_clean; ont01.sum byte-identical) on ont01.

## Root cause — DESIGN "plot count", INVYEAR, STDINFO age/aspect all lost
canada/on (and metric/vbase) `initre.f` declares the keyword `ARRAY(7)/LNOTBK(7)/KARD(7)`,
but `base/keyrdr.f` decodes NF=12 fields and its post-decode `DO 50 I=1,NF: LNOTBK(I)=KARD(I).NE.' '`
writes LNOTBK(8:12) PAST the caller's 7-long LNOTBK into the adjacent ARRAY. The gfortran stack
layout puts ARRAY(1:4) exactly there, so those four are overwritten (blank trailing field ⇒
.FALSE. ⇒ 0.0) AFTER the numeric read. LNOTBK(1:4) itself is NOT overflowed, so each field still
reads as PRESENT with value 0.

Measured on FVSon_g16 (instrumented keyrdr.f/initre.f dumps):
- DESIGN  : KARD(4)="      11.0" reads 11.0 (DBGPOSTREAD) then is clobbered to 0.0 by the LNOTBK
            overflow → IPTINV=IFIX(0)=0 → initre.f:318 clamps `IF(IPTINV.LE.0)IPTINV=1`. Field 5
            (NONSTK=1) survives.
- STDINFO : age(f3)=0, aspect(f4)=0 (both clobbered); slope(f5)=30, elev(f6)=300 survive.
- INVYEAR : year(f1)=0.
- NUMCYCLE: f1 clobbered → rejected (FVS04) → default cycle count.
Confirmed byte-identical in the production-object build (ont01_clean.sum == ont01.sum), so this is
genuine production FVS behavior for ON, not a g16 artifact.

## Fixes (all gated to `s.variant isa Ontario`; multicycle 339/11 byte-identical)
1. `src/engine/keyword_dispatch.jl` process_keywords!: for Ontario, after each top-level card read,
   keep `present[1:4]`, zero `values[1:4]` (emulates the ARRAY(1:4) overflow; fields 5+ survive).
   (kw_design! reverted to the generic single-line field-4 read — now fed a 0.)
2. `src/engine/simulate.jl` compute_density!: set `p.qmd = stand_qmd(s)` (RMSQD) for Ontario —
   ON's Penner dgf! reads `p.qmd*ON_INtoCM`; was never populated ⇒ 0. (p.qmd=9.245 in →
   qmdm=23.48 cm = the DGF dump constant 0x41BBDEFD, so this correctly feeds the live DGF.)
3. `src/io/summary.jl`: ON reports the .sum per-HECTARE metric (like BC) — added Ontario to the
   `met` scaling gate + the write_sum_row metric-format gate.
4. `src/engine/volume.jl` init_merch_standards!: Ontario branch (canada/on grinit.f+sitset.f) —
   TOPD=BFTOPD=10cm·CMtoIN, STMP=30cm·CMtoFT, DBHMIN 5/6 & BFMIND 9/11 & BFTOPD 7.6/9.6 by
   softwood(ISPC≤14|>68)/IFOR; scf_* mirror cubic (ON has no Scribner-cubic). Unblocks the .sum.

## cyc0 (year-0) row, FVSjl vs FVSon_g16 (ont01), column-by-column
| col   | JL     | oracle | |
|-------|--------|--------|-|
| year  | 0      | 0      | ✓ |
| age   | 0      | 0      | ✓ |
| TPA   | 84799  | 84799  | ✓ (=34317.87/ac ×2.471; oracle per-tree FVS40 warnings sum to 34317.89) |
| BA    | 3673   | 3673   | ✓ |
| SDI   | 74769  | 74769  | ✓ |
| CCF   | 15     | ****   | ✗ needs ON open-grown crown-width (canada/on/cwcalc.f) + CCFCAL |
| TopHt | 24     | 24     | ✓ |
| QMD   | 23.5   | 23.5   | ✓ (=9.245 in ×2.54 cm) |

Volume columns (cuft/mcuft/bdft) are 0 in JL (ON volume kernel htont/varvol/cubrds/nbolt not yet
ported — a later chunk); accretion/mortality/MAI and the year-N row require height_growth!/morts.

## Gate
- ontario unit tests: 55/55 (test_ontario 11, _dgf_wired 10, _growth_wired 34).
- multicycle regression: 339 pass / 11 broken — byte-identical (no other variant touched).
- wiring.diff applies cleanly on d1fc55a3.

## NEXT (dep-ordered)
1. CCF: port canada/on/cwcalc.f open-grown crown-width (8 ont01 species: PW SW MH BE SB CE PJ BF,
   IWHO=1) + a CCFCAL wrapper; add an Ontario branch to stand_ccf so density is e2e bit-exact.
2. height_growth!(::Ontario) — canada/on/htont.f (run_keyfile currently MethodErrors here at cycle 1).
3. mortality!(::Ontario) — canada/on/morts.f.
4. Full ON volume (htont broken-top + varvol/cubrds/nbolt/TWIGS) → cuft/mcuft/bdft columns +
   revisit the merch scf_* placeholder.
