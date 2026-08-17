# WWPB (Westwide Pine Beetle) beachhead goldens + architecture + handoff

Branch `wwpb-port` off kt-variant-port 3ec5ba18. Model source: pristine
`/workspace/ForestVegetationSimulator/wwpb/*.f` (53 files, internal prefix `bm`).

## ARCHITECTURE (measured) — WWPB is UNLIKE DFB/DFTM/WPBR

WWPB is a **landscape Parallel-Processing-Extension (PPE)** bark-beetle model
(MPB/WPB/Ips on lodgepole sp 7 + ponderosa sp 10), not a stand/tree-level one.

Two independent keyword entry points:
- **BMIN** — stand-level, keywds.f **option 126** (`base/keywds.f` TABLE(126)='BMIN';
  `initre.f:4839 CALL BMIN(LKECHO)`). TABLE = MAINOUT/TREEOUT/BKPOUT/VOLOUT/END —
  **output scheduling ONLY** (opens .bmm/.bmt/.bmb/.bmv, OPNEW MYACT 2701-2704).
  Does NOT activate the outbreak.
- **BMPPIN** — PPE-landscape, read from PPIN. 40 keywords; **DISPERSE** activates
  the outbreak (GPNEW activity 301 → IBMYR1 start year + duration → IBMYR2). Host
  designations (PBSPEC/HOST), pheromone/salvage mgmt, RANNSEED, etc. all here.

Driver **BMDRV** (bmdrv.f, called from **ALSTD2**) runs only
`IF (LBMSPR .AND. IBMYR1 .GT. 0)`. LBMSPR is set true in bmsetp.f (from ALSTD1)
only when SPLAEX/SPLAAR supply per-stand spatial location+area. FVS mortality
hand-back is **BMKILL** (bmkill.f, from PPMAIN, not BMDRV): converts the TPBK
ledger to per-record `WK2(I)`.

Per-year loop (bmdrv.f DO 30): BMPHER/BMAPH/BMPSTC/BMSALV/BMSANI/BMSMGT →
**BMDRGT** (drought, stochastic) → stand loop {BMCWIN/BMLITE(stoch)/BMFIRE→BMFMRT/
BMOBB/BMDFOL/BMQMRT → **BMMORT(.FALSE.)** fast kills → BMCGRF/BMCBKP/BMCNUM} →
**BMATCT** (landscape attractiveness + BKP redistribution) → stand loop
{BMIPS(stoch)/**BMISTD**(stoch, fills PBKILL) → BMOUT → **BMMORT(.TRUE.)** beetle
kills → BMAGDW}.

**BUILD REALITY**: ALL 24 FVS*_buildDir link the no-op stub `base/exbm.f`
(entries BMIN/BMDRV/BMSETP/BMPPIN/BMKILL/… all RETURN or "**NO BM"). NO
sourceList compiles `wwpb/*.f`; the PPE harness (PPIN/PPMAIN/ALSTD1/ALSTD2/
SPLAEX/GPGET/GPNEW; PPEPRM.F77 lives only under `archive/PPEcommons/`) is ABSENT
— there is no `SUBROUTINE ALSTD2`/`PPMAIN` and no `CALL ALSTD` anywhere in the
tree. ⇒ The beetle OUTBREAK is unreachable via any shipped or straightforwardly-
relinkable single-stand run. A full port must first port/stub the PPE spatial
multi-stand harness. This is materially larger than DFB/DFTM/WPBR.

## RNG — BMRANN (bmrann.f), IDENTICAL LCG to DFTM's TMRANN / DFB's DFBRAN

`BMS1 = DMOD(16807D0*BMS0, 2147483647D0); SEL = BMS1/2147483648D0; BMS0 = BMS1`
(MINSTD a=16807, m=2^31-1, divisor = **2^31 = 2147483648**, NOT m). State in
common /BMRNCM/ (BMS0,BMS1 REAL*8; BMSS REAL*4). Default seed **55329** (odd) set
in every `bmblkd*.f` (`DATA BMS0/55329D0/, BMSS/55329./`). Reseed via RANNSEED
(bmppin.f option 7 → BMRNSD): LSET=true forces odd + stores BMSS/BMS0; LSET=false
restarts stream from BMSS.

### BMRANN golden — seed 55329, first 8 draws (Float32 hex)
`scratchpad/wwpb/driver_bmrann.f` + pristine `wwpb/bmrann.f`, gfortran-16:
```
3EDDB57A 3F5AB21B 3F630A07 3F2734C9 3EF4FB2C 3F4AF646 3F6E7B68 3F67EA59
```
First 6 coincide EXACTLY with DFTM's TMRANN (same LCG + seed). Julia
`wwpb_rand!` reproduces all 8 BIT-EXACT (test_wwpb.jl). Consumers (exact call
order matters for a full port): bmdrgt.f (2 draws→Box-Muller drought normal),
bmlite.f (2), bmips.f (1), bmistd.f (2 sites).

Reproduce: `cd scratchpad/wwpb && gfortran-16 -std=legacy -w -fno-automatic \
  -I/workspace/ForestVegetationSimulator/wwpb driver_bmrann.f \
  /workspace/ForestVegetationSimulator/wwpb/bmrann.f -o drv_bmrann && ./drv_bmrann`

## PORTED + VALIDATED (this beachhead; src/engine/wwpb.jl; test/unit/test_wwpb.jl 33/33)

- **wwpb_rand! / wwpb_seed!** — BMRANN LCG + BMRNSD reseed, bit-exact (8-draw
  golden + odd-force + LSET=false reset).
- **wwpb_defaults!** — BMINIT/BLOCK-DATA defaults: seed 55329, PBSPEC=1 (MPB),
  UPSIZ(10)=/3,6,9,12,15,18,21,25,30,50/, ISCMIN=/3,3,2/.
- **kw_wwpbin!** — BMIN block reader (keywds.f opt 126, dispatched in
  keyword_dispatch.jl on `"BMIN"`): MAINOUT/TREEOUT/BKPOUT/VOLOUT/END → s.wwpb
  (LBMAIN/LBMTRE/LBMBKP/LBMVOL flags + OPNEW output requests with faithful
  IDT(def 1)/nyears(PRMS1 def 100)/incr(PRMS2 def 5) INT-truncated field decode).
- **INERT seam** — a WwpbState projects .sum-byte-identically to no-BMIN
  (validated A/B in test). No simulate.jl per-cycle seam (documented there):
  WWPB applies no mortality without the absent PPE driver.

Wiring: `abstract type AbstractWwpbState` + `wwpb::Union{…,Nothing}` field +
constructor `nothing` (src/core/state.jl); `include("engine/wwpb.jl")` after
wpbr (src/FVSjl.jl); `"BMIN"` dispatch (keyword_dispatch.jl); WWPB no-seam note
(simulate.jl).

## NEXT SUB-CHUNKS (dependency-ordered) — a FULL port is large

The reachable, self-contained deterministic slices below can be dump-replay-
validated WITHOUT the PPE harness by relinking a tiny gfortran-16 driver over the
pristine routine (feed dumped inputs, compare hex) — mirror driver_bmrann.f. The
outbreak itself requires the PPE harness port (hardest, do last).

1. **BMFMRT (bmfmrt.f)** — fire-caused beetle-tree mortality: deterministic
   logistic bark-thickness models per species (codes: 7=lodgepole, 10=ponderosa;
   `bmfmrt.f:154,162`), scorch-height → PMORT. Pure math; driver-testable.
2. **BMQMRT (bmqmrt.f)** — "quick"/extra slow mortality proportion YMORT→OAKILL;
   deterministic given scheduled params (GPGET2 keyword 307).
3. **BMCGRF → BMCBKP/BMCNUM (bmcgrf/bmcbkp/bmcnum.f)** — the deterministic
   GRF/BKP/attractiveness-numerator pipeline. BMCGRF is smallest, no scheduling.
4. **BMMORT (bmmort.f)** — deterministic tree-array decrement given dumped
   PBKILL/OAKILL/ALLKLL (PRDEAD = (PBKILL+OAKILL(,,1)+ALLKLL)/TREE; TREE*=1-PRDEAD)
   + the TPBK ledger.
5. **BMATCT (bmatct.f)** — landscape attractiveness + BKP redistribution
   (calls SPLAAR — needs spatial stubs).
6. **Stochastic dynamics** (exact BMRANN stream + call order): bmdrgt.f, bmlite.f,
   bmips.f, **bmistd.f** (fills PBKILL).
7. **BMKILL (bmkill.f)** — the FVS hand-back: TPBK → WK2(I), bounded
   PROB(I)-WK2(I) ≥ 1e-6; SDWP/DDWP dead-wood pools; SVMORT.
8. **PPE harness** (blocking for any end-to-end run): PPIN/PPMAIN/ALSTD1/ALSTD2 +
   SPLAEX/SPLAAR spatial + GPGET/GPNEW/GPADD activity scheduler + PPEPRM.F77
   (archive/PPEcommons). Without this, no oracle can exercise BMDRV. A single-
   stand degenerate harness (MXSTND=1) + a synthetic DISPERSE (IBMYR1>0) +
   spatial stub (LBMSPR=T) is the minimal path to a relinkable oracle
   (`FVSbm_wwpb` = FVSbm .o − exbm.o + wwpb/*.o + PPE stubs). Verify WWPB-off ≡
   stock first, then BMIN-on ≡ off (confirms our inert model), then feed a
   pine-host stand with DISPERSE.

## Files
- `scratchpad/wwpb/driver_bmrann.f` — BMRANN golden generator (reusable).
- `scratchpad/wwpb/drv_bmrann` — compiled binary (scratch; regenerate as above).
