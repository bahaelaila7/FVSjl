# WRD + SVS port + OP loose-end — rollout plan (2026-08-16, USER-directed)

USER directive 2026-08-16: "finish the loose ends, and port WRD and visualization."
Visualization scope (USER-chosen): **Port FVS SVS data path** (per-tree object list FVS emits for
the SVS viewer), validated vs the live FVS `.svs` output. NOT a GUI/viewer.

Branch: kt-variant-port. Doctrine unchanged: bit-exact-or-cornered vs LIVE relinked FVS oracle,
MEASURE don't infer, keep Fortran PRISTINE, commit only validated chunks, never FFI the RNG.

## ★ OP loose-end — ORACLE UNBLOCKED (premise reversed 2026-08-16)
The recorded "FVSop_clean SIGSEGVs at cycle≥1 ⇒ OP multi-cycle unvalidatable" was WRONG. It was an
artifact of (a) the DEBUG/fvsvol dump path and (b) a wrong invocation. FVSop_clean reads FOUR stdin
filenames (keyword/tree/out/treelist); STOP 20 is this build's NORMAL termination. Fed correctly it
produces a clean multi-cycle `.sum`. Reference (stand S248112, NUMCYCLE=2, /workspace/.opwork/op2c.*):
```
YEAR TPA  BA  SDI CCF TopHt QMD TCuFt MCuFt BdFt
1990 536  77  184 114  63   5.1 1472   972   5003
1995 499  98  219 134  71   6.0 2256  1745   9249
2000 470 120  256 154  69   6.8 2988  2419  12427
```
So OP multi-cycle IS validatable. Remaining OP work = wire the growth driver AND port the op-native
crown + mortality (+ small-tree) for the IORG=0 (op-Wykoff) trees:
- OP already has (validated, committed): op native DGF/HTGF (chunk1), full NWO ORGANON engine
  (op_build_organon_buffer!/op_prepare_nwo/op_execute_nwo → DG/HGRO/CR2/DEADEXP/DDS), CCF, cyc0 .sum.
- MISSING (cause of the `small_tree_growth!(::Olympic)` MethodError at cycle 1):
  1. `op_organon_prepare!(s)` LSTART wrapper (mirror oc_organon_prepare!; op_prepare_nwo exists) —
     wire into src/engine/simulate.jl:26 alongside the OregonCoast branch.
  2. single-authority `diameter_growth!(::Olympic)` (mirror oc organon_hook.jl): native dgf! WK2 for
     IORG=0 + op_execute_nwo for IORG=1, apply DBH/HT/CR/MORT copybacks, zero growth fields.
  3. no-op `height_growth!`/`small_tree_growth!`/`mortality!`/`crown_ratio_update!`(::Olympic)
     (ORGANON does it all in one EXECUTE, like OC).
  4. op-native **crown** (op/crown.f) + **mortality** (op/morts.f + op/mortality.f) + small-tree for
     the IORG=0 trees — these ARE the multi-cycle gap (opt01 shows real mortality 536→499). Validate
     bit-exact per-tree via FVSop_g16 (build_g16.sh dbg swap) then the op2c multi-cycle .sum.
Oracle recipe: `printf "op2c.key\nop2c.tre\nop2c.out\nop2c.trl\n" | ./FVSop_clean` in /workspace/.opwork.

## WRD (Western Root Disease) — scope (measured)
- rd/ = 93 .f / 22,632 lines. **Compiled into ONLY 15 variants**: bc bm ci cr ec em ie kt nc pn so
  tt ut wc ws. NOT ca/ak/oc/op/on (only rdpsrt). DEAD (do NOT port): rdinca.f/rdinitca.f/rdppatv.f/
  rdblk1ca.f (~3900 lines), rdpsrt.f (already ported = quickersort.jl).
- Variant-specific = ONLY rdblk1{var}.f block-data (host-species coeff tables). KT uses default
  rdblk1.f. BC/ON metric via rdinit.f:728 LMTRIC (BC = only live metric WRD variant).
- Activation: DORMANT by default; RRINIT→RRMAN, RRTREIN→RRTINV; rdatv.f gate L=RRTINV.OR.RRMAN.
  Keyword reader = rd/rdin.f (42-keyword TABLE, entry RDKEY). RRTYPE selects disease (1 non-host,
  2 P-annosus, 3 S-annosus, 4 armillaria, 5 phellinus; ITOTRR=4 modeled).
- Engine seam (base/): fvs.f:166 RDMN1(1) init, :346 RDMN1(2)/RDPR each cycle; grincr.f:222 RDMN2;
  gradd.f:131 RDTREG (growth-cycle disease); intree.f RDATV/RDTRES; triple/tredel/comprs RD-sync;
  initre.f RDIN (keyword parse). RDPARM.F77 (IRRTRE=1500,ITOTSP=40,ITOTRR=4) + /RDADD/RDCOM/RDCRY/
  RDARRY commons = infra to port first.
- Turnkey oracle EXISTS: tests/FVSie/DBReportTest.key (RDin/RRType 3/RRInit .../SArea 100/RRDOut/End).
- Chunk plan: -1 infra (commons+RDPARM+rdblk1 KT) → 0 mortality path (RDIN{RRTYPE,RRINIT,SAREA,END}+
  rdatv+rdmn1/2+rdtreg+rdcntl reduced to RDSETP→RDINOC→RDAREA→RDMORT+rdpr/rdsum; DF+S-annosus; KT) →
  1 spatial spread (rdsprd/rdcent/rdcloc/rdshrk/rdarea+SPREAD) → 2 inoculum/infprob → 3 RRTREIN
  treelist-init (rdesin/rddam/rdtres+sync) → 4 growth-reduction (rdgrow+TTDMULT/SDIRMULT) → 5
  sporulation/annosus → 6 reporting → 7 bark-beetle (rdbb1-4) → 8 BC metric → 9+ per-variant rdblk1.

## SVS (Stand Visualization System data path) — scope (measured)
- Variant-GENERIC (single copy in base/vbase; NO per-variant SVS code). Coupling only via arrays SVS
  reads: CRWDTH (crown_width.jl, ported), ICR, species codes, #TREEFORM (western cluster → WEST.TRF).
- Keyword: SVS (base/svkey.f). Output: <key>_index.svs (#TREELISTINDEX) + <key>_NNN.svs. Writer =
  vbase/svout.f:243-440 (header + per-tree object record FORMAT at :429-440).
- Needs a SECOND RNG stream: base/svrann.f Park-Miller 16807/2^31-1 but SEL=SVS1/2^31 (divisor 2^31,
  NOT the base rann divisor). Seeded at svstart.f:50 RANNGET(SVS0) from the base stream's current s0
  (tracks RANNSEED). Port as a second FVSRng seeded at the SVSTART seam. NEVER FFI.
- Placement: svestb.f (TPA→objects + fractional lottery, RDPSRT sort, SVRANN draw) → svgtpt.f (x,y:
  rectangle 2-draw for IPLGEM<2, disk for ≥2) → svobol.f overlap reject/retry (≤40). svgtpl.f =
  deterministic subplot layout.
- Engine seam (base/): fvs.f:333 SVSTART (cyc0 picture), grincr.f:277 SVOUT(...,1) each later cycle,
  cuts.f:1868 SVOUT(...,2) post-salvage, fvs.f:453 SVOUT(...,3) final; gradd.f→SVESTB(1) regen add.
- Oracle: FVSkt (relink), stand S248112 (tests/FVSkt/ktt01). No .key uses SVS ⇒ author one.
- Chunk plan: 0 deterministic single-record slice (SVS 0, IPLGEM=0, integer TPA; svkey+svgtpl(0)+
  svgtpt(rect)+svrann+svestb(int)+svout live-tree loop, SVSTART seam, cyc0 diff) → 1 fractional
  lottery+overlap → 2 multi-subplot layout → 3 multi-cycle (grincr/fvs.f seams) → 4 cut/mort removal
  (svrmov+cuts seam) → 5 snags (svsnad/svsnage/svsalv) → 6 (deferred) FFE CWD + ground file.
