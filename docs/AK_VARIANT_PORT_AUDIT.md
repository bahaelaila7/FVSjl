# AK (Southeast Alaska) variant port — audit (chunk 0 beachhead)

FVS variant **AK** = "ALASKA", USFS **Region 10**, **MAXSP = 23**, standalone. A standard
**Wykoff DDS** diameter-growth variant (Zeide SDI) with three AK-unique subsystems: a
**PERMAFROST** diameter-growth modifier, the **SEAMRT** mortality-distribution routine, and
**Region-10 volume**. Source: `/workspace/ForestVegetationSimulator/ak/*.f` +
`bin/FVSak_buildDir`. This run is a bounded, validated beachhead — foundation + oracle +
large-tree DGF cyc0 validated bit-exact + a scope note on the three unique subsystems.

Reuse anchor = the ported Utah variant (`src/variants/utah/`), itself a Wykoff-DDS
specialization of the shared engine. AK mirrors that structure.

## Oracle

- Relinked live Fortran: `/workspace/.akwork/FVSak_clean` (596 `bin/FVSak_buildDir/*.o` +
  `/workspace/.crwork/isoc23_shim.o`, gfortran-16). Script: `/workspace/.akwork/relink_ak.sh`
  (`clean` = pristine; `<name> <patched.o>` = single-object swap for instrumentation).
- Reference stand: `tests/FVSak/akt01` (INV 1990, 27 records, species YC/WH/RC/LP/MH/SS).
- DGF per-tree dump obtained with the built-in `DEBUG / DGF / END` keyword block — no custom
  instrumentation needed; the routine's own `IF(DEBUG)WRITE` statements dump every
  intermediate (D, D², ln D, PBAL, PRD, CR, ln CR, TEMEL/TEMSLP/TEMSASP, SSITE, DGCONB1,
  DGCOMP1, DGCOMP2, BASEDG, PFMOD, DGPRED, BRAT, TEMPD1/TEMPD2, COR, DGCON, DDS, WK2).
  Invoke `FVSak_clean --keywordfile=akdbg.key` (auto-derives the `.tre`); dump lands in the
  main `.out` / `fort.16`.

## Chunk verdicts

| Chunk | Scope | Verdict |
|---|---|---|
| 0 Foundation | `SoutheastAlaska <: AbstractVariant`, `variant_from_code("AK")`, 23-species table, grinit constants, `FVSjl.jl` include + export | **DONE** — module loads, resolves, coefficients load |
| 3 Large-tree DGF | `ak/dgf.f` ln(DDS) equation + all coefficient DATA arrays + `ak/bratio.f` bark | **VALIDATED BIT-EXACT** vs live oracle: 54 tree-records, 0 mismatches, WK2 rel-err **0.0** (BASEDG 8e-8 = Float32 ULP on `exp`); AK 3-type bark path also exact (no BRAT drift) |
| 2 Site index / SDImax / forkod | `ak/sitset.f` SITEAR (ISISP=0→11, TEM=70, SLO/SHI interp) + SDIDEF (SDICON) + `ak/forkod.f` (Tongass 1005/IFOR 2 default) | **VALIDATED BIT-EXACT** — akt01 XSITE dump: YC=50, SS=82.5, LP=35, RC=57.5, WH=70, MH=42.5 all match |
| 4 Large-tree HTGF | `ak/htgf.f` single Wykoff HG equation + NOPERM/PERM coeffs + HTLO/HTHI bounding + species mult, wired via `height_growth!(::SoutheastAlaska)` | **VALIDATED BIT-EXACT** vs live FVSak `DEBUG HTGF` on akt01: **27 trees, HTG rel-err 0.0** (POTHTG to F8.4 print precision, ULP-level). Primary gate PASSED. `tools/southeastalaska/validate_htgf.jl` |
| 4b Height-diameter dub | `ak/cratet.f` Curtis-Arney INVENTORY-EQN (LHTDRG=false): H=4.5+HTT11·(1−exp(HTT12·D))^HTT13 ×spmult; HTT11/12/13 MEASURED from `DEBUG CRATET` (all 23 sp) | **DONE** — wired into `dub_missing_heights!`; akt01 dubbed heights match (5-decimal coeff precision) |
| 5 Crown | `ak/crown.f` logistic CR (PRD/HDR/D-QMD) + `ak/dubscr.f` (bachlo RNG) + point-Zeide `ak/sdical.f` SDICAL/SDICLS (XMAXPT/ZRD) | **PORTED** — point-Zeide PRD reproduces the oracle (365.0/592=0.6166 verified); feeds DGF ln(CR) from cyc2+. Not yet independently per-tree-validated |
| 7 Mortality | `ak/morts.f` logistic survival (BM1-5) + SDI/BA iterative-PASS multiplier (NOT SEAMRT — morts.f doesn't call it) | **PORTED** — end-to-end TPA tracks oracle within ~1% mid-run; late-cycle self-thin selection straddle |
| 6 REGENT small-tree | `ak/regent.f` | **STUB (no-op)** — small trees keep large-tree DGF/HTGF; akt01 is mature so bounded. LATER chunk |
| 8 Volume | `ak/sitset.f` VOLEQDEF(VAR='AK',IREGN=10)→NVEL + `ak/logs.f`/`cubrds.f` | **STUB (cuft=0)** — R10 NVEL crosswalk not ported; .sum cuft/bdft columns are 0. LATER chunk |

### End-to-end akt01 `.sum` vs `akt01.sum.save` (unthinned control, NUMCYCLE 10, NOTRIPLE off)

Fix that unblocked it: **AK GRINIT `BAF=62.5`** (ak/grinit.f:173) — jl was using the generic 40 default; akt01's DESIGN omits BAF, so the 62.5 default is load-bearing for the plot expansion (added in `engine/init.jl`). With it, **cycle 0 (1990) is BIT-EXACT**: TPA 669 / BA 118 / SDI 242 / TopHt 64 / QMD 5.7 all match.

| Year | jl TPA/BA/SDI/QMD | oracle TPA/BA/SDI/QMD | ΔBA |
|---|---|---|---|
| 1990 | 669/118/242/5.7 | 669/118/242/5.7 | **0% (exact)** |
| 2000 | 615/144/284/6.6 | 615/147/290/6.6 | −2% |
| 2020 | 545/188/348/8.0 | 553/198/373/8.1 | −5% |
| 2050 | 467/240/416/9.7 | 482/256/453/9.9 | −6% |
| 2070 | 424/266/447/10.7 | 421/270/462/10.8 | −1.5% |
| 2090 | 386/286/467/11.7 | 354/268/446/11.8 | +7% (jl retains more TPA) |

**Residual classification (chunk 1, #209 stream — RESOLVED):** the ~3% YC-COR deficit was a **REAL deterministic bug, now FIXED** (commit below). ROOT CAUSE (MEASURED via `FVSak_g16 DEBUG DGDRIV` + jl per-tree dump): jl's LSTART calibration left the point BA-in-larger-trees (`density.point_bal`, PTBALT) at the **backdated**-dbh value, but FVS's `dense.f` two-pass backdating recomputes WK5 at **current** dbh on the second pass (dense.f:184 sets `D=WK3` only when `LREDO`), and `PTBAL` (dense.f:280 → ptbal.f:148) reads that current WK5 — so the calibration DGF's `PBAL=PTBALT(I)` is **current**, exactly like PTBAA. On akt01 YC (sp 3) tree i=3 jl had PBAL=52.8 vs live 187.5 (=3×BAF 62.5) ⇒ the DGF prediction WK2 over-predicted ⇒ RESLOG under-shot ⇒ **COR[YC]=1.3894 vs live 1.4342** (goal 0.6947 vs 0.7171). FIX: `calibrate_diameter_growth!` now stashes the current-dbh `point_bal` alongside `cur_point_ba` and restores it after backdating (AK-guarded `_ak_pbal_fix`; shared-faithful — same latent gap likely applies to all western variants but is validated here only for AK). After the fix, **COR[YC]=1.4342394 == live 1.4342394 BIT-EXACT** (goal 0.7171197), cyc0 stays bit-exact, and **cyc1 (2000) BA=144 == live 144, TPA=615 == live 615 BIT-EXACT** vs the relinked live `FVSak_clean` (NOTE: the shipped `akt01.sum.save` was produced by an older build — 2000 BA 147 — the LIVE relinked oracle gives 144, which jl now matches; doctrine: validate vs the live oracle, not the shipped save). The remaining multi-cycle BA residual is a **±straddle that OSCILLATES sign** (2010–2030 ~−2…−5%, 2050 exact, 2060–2090 +4…+7%) with NO compounding deterministic bias = the DGSD=2.0 **OLDRN serial-correlation + mortality-selection straddle** = **CORNERED** (the western-cluster bar). Deterministic DDS proven equal: COR and the calibration-pass per-tree WK2 are now bit-exact, and cyc0+cyc1 aggregates match live. Volume + the summary CCF column remain stubs (CCF=0; volume = chunk 8 below).

Reproduce: `julia --project=. tools/southeastalaska/validate_dgf_cyc0.jl <FVSak-DEBUG.out>`
(pulls the *shipped* `src/variants/southeastalaska` coefficients + `ak_bratio`, feeds the
oracle's dumped raw inputs, checks WK2).

### grinit constants (ak/grinit.f, ak/blkdat.f) — all faithful

VARACD `AK`; MAXSP 23; `LZEIDE=.TRUE.` (Zeide SDI, like CR/UT); `DGSD=2.0`; `FINT=10`
(YR=10-yr cycle); `FINTH=5`; `FINTM=5`; `LHTDRG(I)=.FALSE.` (default — **unlike** UT's
`.TRUE.`); RNG seed 55329 (shared blkdat). Species (1..23): SF AF YC TA WS LS BE SS LP RC WH
MH OS AD RA PB AB BA AS CW WI SU OH. FIA: 011 019 042 071 094 094 095 098 108 242 263 264 299
350 351 375 376 741 746 747 920 928 998. (LS "Lutz's/hybrid spruce" has no native FIA; `sitset.f:321`
maps JSP='LS' → 94. `OS` FIA is 299 in blkdat's species-list header but the `dgf.f` header comment
says 298 — a source inconsistency, flagged; irrelevant to DGF, matters only for volume-eq lookup.)

### AK DGF model (single equation, all species — no per-species DDS branches)

```
BASEDG = exp( b1 + b2·D² + b3·ln(D) + b4·PBAL + b5·PRD + b6·ln(CR)
              + b7·ELEV + b8·SLOPE + b9·SLOPE·cos(ASPECT) + b10·ln(SI) )   # annual, outside-bark
DGPRED = YR · BASEDG · PFMOD · ak_dgmult(sp)                              # period, +permafrost, +sp mult
DDS    = ln( (D+DGPRED)²·BRAT² − D²·BRAT² ) + COR(sp) + DGCON(sp)         # inside-bark Δdib², WK2
```
with `ELEV=elevation·100`, `SLOPE=slope·100`, `CR`=crown ratio **percent**, `SI=SITEAR(sp)`,
`PRD`=point ZeideSDI/point maxSDI, `BRAT=ak_bratio(sp,D)` (3 equation types). Species DG
multipliers: AD/WI/OH → 0.45, SU → 0.65, else 1.0. `DGCON` = READCORD term (+ln COR2 if LDCOR2);
`COR` = the shared-engine DGSCOR self-calibration. This is a cleaner form than UT/IE (which
branch DDS per species group); AK converts an annual DOB increment through bark to Δdib².

## Scope notes on the three AK-unique subsystems (LATER runs)

**1. PERMAFROST DG modifier (PFCON) — `ak/dgf.f`.** A multiplicative modifier `PFMOD` on the
annual base DG, applied to permafrost-affected species (SELECT CASE `4:7, 13, 16:23` = TA/WS/LS/BE,
OS, PB/AB/BA/AS/CW/WI/SU/OH). Second regression (own coeffs PFCON/PFDSQ/PFLD/PFDBAL/PFRD/PFLNCR/
PFEL/PFSLOP/PFSASP + presence factor PFPRES): `PFMOD = exp(PFCON + [PFPRES if LPERM] + …)/BASEDG`.
`LPERM` true (PERMAFROST keyword on) ⇒ cap `PFMOD ≤ 1` (permafrost *slows* growth); off ⇒ floor
`PFMOD ≥ 1`. All coefficient arrays are **already transcribed & shipped** in
`dg_coefficients.jl` (AK_PF*) and the modifier logic is **already in `dgf!`**. What remains: (a)
wire the `PERMAFROST` keyword → an `LPERM` control flag, and (b) the `PRD` dependency below.
Not exercised by akt01 (all its species are non-permafrost, `DGRD=0`), so the cyc0 validation
covers the base path exactly; the permafrost path needs a permafrost-species stand to validate.

**2. PRD = point Zeide relative density — the one engine gap.** Both the base `b5·PRD` term and
the permafrost `PFRD·PRD` term use `PRD = point ZeideSDI / point maxSDI`, computed in `ak/dgf.f`
via `SDICAL`(→`XMAXPT`) + `SDICLS`(→`ZRD`) per subplot. The shared FVSjl engine has per-point
BA/BAL/CCF/TPA (`Density`) but **no per-point Zeide-SDI** machinery. `ak_point_zeide_rd` is a
documented stub returning 0 — **exact** for the `DGRD=0` coastal majority (SF/AF/YC/SS/LP/RC/WH/MH,
which dominate SE Alaska), an approximation only for the permafrost/interior species. Porting
`SDICAL`/`SDICLS` point-Zeide is the prerequisite for validating the permafrost + interior path.

**3. SEAMRT mortality — `ak/seamrt.f`.** AK's density-mortality *distribution* routine (analogous
to eastern VARMRT): given a stand kill target `TOKILL`, it distributes deaths across records by a
per-tree efficiency `EFFTR = PEFF·((100−CR)/100)·VARADJ(sp)·0.01`, where `PEFF = 0.84525 −
0.01074·PCT + 2e-7·PCT³` (BA-percentile) and `VARADJ` = a 23-species shade-tolerance array
(0.100 most-tolerant … 0.900 most-intolerant; transcribed below). A geometric-progression
`NPASS` search converges the distributed kill onto `TOKILL`. `TOKILL` itself comes from the AK
`morts.f` SDIMAX self-thin driver. `VARADJ` (SF..OH): 0.300 0.100 0.300 0.500 0.500 0.500 0.500
0.300 0.700 0.100 0.100 0.300 0.700 0.500 0.500 0.500 0.500 0.900 0.900 0.900 0.500 0.500 0.500.

**4. Region-10 volume — `ak/sitset.f` + NVEL.** `sitset.f` assigns per-species NVEL equation
codes via `VOLEQDEF(VAR='AK', IREGN=10, FORST, DIST, IFIASP, PROD, …) → VEQNNC`; volume then runs
through the shared National Volume Estimator library plus AK-local taper (`ak/logs.f` cubic/board
via `cubrds.f` defaults). `IREGN` derives from `KODFOR/100`; BC/Makah/adjacent location codes
(713/720/703, IREGN 74/81) are remapped into Region 10 forests 04/05. This is the standard
`VEQNNC → NVEL` path already used by the other western variants' volume chunks (reuse candidate).

## Files added

- `src/variants/southeastalaska/{southeastalaska,species,dg_coefficients,diameter_growth}.jl`
- `data/southeastalaska/{species_coefficients,species_translation}.csv` (identity + real AK bark;
  other per-species subsystem columns are later chunks)
- wired: `src/variants/variant.jl` (`variant_from_code`), `src/FVSjl.jl` (includes + export)
- `tools/southeastalaska/{validate_dgf_cyc0.jl,akdbg.key}` (reproducible cyc0 validation)
- `/workspace/.akwork/{FVSak_clean,relink_ak.sh}` (oracle, outside the repo)

## What remains for AK (later runs)

Site index / SITEAR (`ak/sitset.f`) · crown ratio (`ak/crown.f`) · height growth
(`ak/htgf.f`, which also carries a permafrost modifier) · small-tree REGENT (`ak/regent.f`) ·
per-point **Zeide-SDI density** (SDICAL/SDICLS) to unlock **PRD** + the **permafrost** DG/HTGF
path · **SEAMRT** mortality + `ak/morts.f` SDIMAX driver · **Region-10 volume** (VEQNNC → NVEL +
ak logs/cubrds) · establishment · full engine-pipeline integration + a multi-cycle `.sum`
differential vs `FVSak_clean` (incl. an interior permafrost-species stand). The DBH-update step
must use `ak_bratio` (3-type), not the engine's `b+a/d` bark surrogate.
