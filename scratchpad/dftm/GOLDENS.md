# DFTM chunk-0 goldens (from pristine dftm/*.f via gfortran-16 drivers)

## TMRANN (tmrann.f), default seed 55329 — first 6 draws (Float32 hex)
driver_tmrann.f + dftm/tmrann.f :
  3EDDB57A 3F5AB21B 3F630A07 3F2734C9 3EF4FB2C 3F4AF646
Julia dftm_rand! reproduces all 6 BIT-EXACT.

## TMOTPR (tmotpr.f) formula goldens — driver_tmotpr.f (formulas verbatim)
Method 1 (Heller):   elev35 slope.30 aspect1.2 topo2 relden120 reldsp3=40 reldsp4=25 tprob200 -> 3F2DE146
Method 2 (Mika+ash): topo3 ash15.93 pgfba.4 ba180 -> 3DDC07CC
Method 3 (Mika noash):topo3 pgfba.4 ba180 -> 3DF4118C
Julia dftm_otpr matches within <=4 ULP (transcendental exp/cos/sin/sqrt/log straddle).

## Relink oracle: FVSie_dftm  (build_ie_dftm.sh)
- IE object set, exdftm.o dropped, real dftm/*.o (+generic tminit.f) added, +base/exppe.f stub (PPECYC).
- DFTM-OFF ie_off.sum == FVSie_clean clean_off.sum  BYTE-IDENTICAL (stock preserved).
- DFTM-ON  ie_on.key (MANSTART+MANSCHED 2+DEBUG) runs to normal exit(20); .sum DIFFERS from off
  (INSCYC forces the TMBASE=5yr cycle => altered cycle boundaries; full defoliation needs a denser host stand).
- Invoke: `echo ie_off.key | ./FVSie_dftm`  (keyfile NAME on stdin).

## CHUNK 1 DONE (branch dftm-coup, commit ffa66d65) — LIVE dump-replay goldens
Dense DF/GF host stand `dense.key`+`dense.tre` (10 DF + 8 GF) CLEARS DFTMGO (L=T,
INSCYC forces the 5-yr TMBASE outbreak cycle → cycle inserted at 2005). Invoke:
  `./FVSie_dftm --keywordfile=dense.key`  (companion `dense.tre` auto-opened; NOT stdin).
  ⚠ The chunk-0 ie_on.key embedded tree records INLINE → FVS read them as KEYWORDS
    (all-zeros stand). TREEDATA must be EMPTY in the .key with a companion .tre.

Instrumented relink for goldens: `ovr/` = pristine dftm/{tmotpr,dftmgo,tmbmas}.f +
stderr `WRITE(0,...) TRANSFER(x,i)` Float32-hex dumps (DBGOTPR*, DBGGO*, DBGBMAS2).
  `bash build_ie_dftm.sh <ABS ovr dir> FVSie_dftm_dbg` ; instrumented .sum is
  BYTE-IDENTICAL to the clean relink (verified). Binary lands in ieobj/ — cp out.

Ported + VALIDATED (src/engine/dftm.jl; test/unit/test_dftm.jl 59/59):
- dftm_otpr  method 1 PROTBK — LIVE bit-exact (0 ULP). Goldens: IN elev34 slope.3
  aspect5.498 topo1 tmashd15.93 ba72 relden83.08 reldsp3=43.83 reldsp4=39.25 tprob97.49
  -> 3F464F97 ; and ba85.69/relden93.16/... -> 3F4DB2B7.
- dftm_go_gate — CNTDF=424E0DC2 (51.5134, Σ DF PROB serial), CNTGF=421CF227 (39.2365),
  IDF=10 IGF=8 NACLAS=(10,8) L=T. Bit-exact (0 ULP).
- dftm_bmas2_df/_gf — TMBMAS IBMTYP=2 (deterministic; 1/3/4 draw TMBCHL off TMRANN).
  DF FBIOMS + all PCNEWF bit-exact; GF FBIOMS ≤1 ULP (expf straddle). 18-tree goldens
  embedded in the test (DBGBMAS2 rows).

## NEXT chunk oracle work (dependency-ordered)
- TMSCHD (tmschd.f) RANSCHED/MANSCHED auto-scheduler + INSCYC (inscyc.f) cycle-forcing
  — INSCYC MUTATES FVS cycle boundaries (IY/NCYC/IFINT); needs an engine hook, not a
  pure fn. dense.out shows the insertion (IFROM2 IBOUND5 ISPOT3 → IY 1990 2000 2005...).
- garbel.f/grclas.f classification (IPT DBH-sort, CLASS SECTOR POINTERS ISC, JCLASS) +
  the ICON common (Z4=PROB Z5 Z2=PCNEWF Z3=FBIOMS X5=NEWBMAS X6=OLDBMAS X7=EGGS) —
  all in dense.out; instrument garbel.o for hex.
- dftmod.f integrator (g0comp/gfcomp/y0comp/y1comp/z1comp/uv1/uv2/dfole8) → DPCENT1/
  DPCENT2 (dense.out "PERCENT BRANCH DEFOLIATION"). ⚠ EMPTY GF class 11 (PROB 0) makes
  X(1)=NaN at occasion 9 → "RESET TO 0.0"; the NaN reproduces bit-exact only if the
  class build + integrator match — replicate the reset guard.
- tmcoup.f (977 ln) coupling → per-tree WK2 (mortality), DG/HTG loss, top-kill. Final
  "DUMP OF INDIVIDUAL TREE RECORDS: IPT IMC SPECIES WK2 PROB HT HTG DBH DG ..." in
  dense.out is the treelist feedback to validate. Then the .sum-DELTA (on−off) vs
  FVSie_dftm honoring the IE #206 OLDRN growth straddle. NEVER FFI TMRANN.
- Engine seam (simulate.jl): DFTMGO in grincr predict-phase (→ TMBMAS), TMCOUP in
  gradd after MISTOE/before FMMAIN (gradd.f:103). Mirror the DFB seam gating.
