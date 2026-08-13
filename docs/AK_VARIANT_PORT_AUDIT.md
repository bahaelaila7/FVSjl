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
| 8 Volume | `ak/sitset.f` VOLEQDEF(VAR='AK',IREGN=10)→NVEL F32 Flewelling profile (reuses shared `_fw2_*` kernels) + `setcubicdflts.f`/`mrules.f` R10 merch + 32-ft-log board | **VALIDATED BIT-EXACT (per-tree) vs live FVSak_clean TREELIST on akt01 cyc0** — all 29 trees' total cubic + merch cubic + Scribner board match (e.g. LP 21.4/14.5/60, WH 24.9/17.8/60, MH 13.2/10.1/30, YC-dead 240.3/223.9/1010). .sum aggregates: **MCuFt 732 = live 732, BdFt 2417 = live 2417 bit-exact**; TCuFt 2316 vs 2315 (±1 = Float32 summation of the per-acre total, below integer precision). DVE (woodland) + CUR (hardwood) families still 0 (akt01 has none). See §Chunk 8 below |

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
`PFMOD ≥ 1`. **WIRED + VALIDATED BIT-EXACT (#209).** The `PRMFROST` keyword (keywds.f TABLE(146) —
NOT "PERMAFROST"; field 2 = 1 ON / 0 OFF) is parsed (`kw_permafrost!`) into `control.permafrost` and
`dgf!` reads it. **KEY TIMING (MEASURED):** grincr sets LPERM before the GROWTH DGDRIV, but the LSTART
DG/crown calibration runs earlier in **CRATET with LPERM still .FALSE.** — verified from the live dump
(the CRATET DGF passes show PFMOD>1 = the floor≥1 branch even with `PRMFROST` ON). So jl applies LPERM
only on the growth pass, not calibration (gated on `calib_dbh` non-empty). VALIDATED on a white-spruce
(WS, permafrost sp 5) stand vs `FVSak_g16 DEBUG DGF`: (A) same backdated stand, LPERM=false — jl PFCOMP2
(**the PFRD·PRD term**) + PFMOD bit-exact per-tree (max Δ 2.5e-7 / 5.3e-7, 27 trees); (B) LPERM=true —
jl PFMOD = the oracle-derived `min(1, exp(PFCON+PFPRES+PFCOMP1+PFCOMP2)/BASEDG)` from the live-dumped
intermediates, bit-exact (max Δ 2.6e-7, 27 trees). The live LPERM=true GROWTH-pass dump itself is blocked
by the DVE-volume segfault (all permafrost species map to DVE/CUR volume, unported+crashing in the g16
build), so LPERM=true is validated against the live intermediates on the identical stand rather than a
direct growth dump. Inert for akt01 (all `DGRD=0` non-permafrost species: YC COR 1.4342394 + cyc0
unchanged). Composes with the now-live PRD (scope note 2) via the PFRD·PRD term.

**2. PRD = point Zeide relative density — PORTED + WIRED + VALIDATED BIT-EXACT (#209).** Both the base
`b5·PRD` term and the permafrost `PFRD·PRD` term use `PRD = point ZeideSDI / point maxSDI` (`ak/dgf.f`
via `SDICAL`→`XMAXPT` + `SDICLS`→`ZRD` per subplot). `ak_point_zeide!` (crown.jl) ports SDICAL (XMAXPT =
BA-weighted SDIDEF per point) + SDICLS (ZRD = Σ PROB·(PI−NONSTK)·(D/10)^1.605, D≥DBHZEIDE=0), and `dgf!`
now reads `PRD = ZRD(pt)/XMAXPT(pt)` per tree (was a 0 stub). **KEY (MEASURED, same class as the chunk-1
PTBALT fix): during LSTART calibration SDICAL/SDICLS sum the UNCHANGED `DBH(I)` = CURRENT dbh** (FVS
backdates only `DIAM(I)`), so jl stashes the current dbh (`calib.calib_dbh`) for `ak_point_zeide!` —
without it jl computed PRD on the backdated stand (WS D11.5 → PRD 0.2257 vs live 0.2645). VALIDATED vs
`FVSak_g16 DEBUG DGF`: (a) akt01 per-point PRD bit-exact (pt1 0.160517558, pt2 0.441927820, pt4
0.311797887, … all = live to Float32 ULP); (b) a white-spruce (WS, `DGRD=−0.4555`) stand — **all 27
trees' DGF WK2 (deterministic ln-DDS, incl. the b5·PRD term) bit-exact, max rel-err 3.4e-6**. Inert for
akt01 (all-`DGRD=0` species: YC COR 1.4342394 + cyc0 unchanged). The interior/permafrost path now has
its PRD; the PFRD·PRD permafrost term still needs the PERMAFROST keyword wired (scope note 1).

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

## Chunk 8 — Region-10 volume (F32 Flewelling profile) — VALIDATED BIT-EXACT (per-tree, akt01 cyc0)

**Crosswalk (MEASURED from live `FVSak_clean` NVEL equation dump):** `ak/sitset.f` VOLEQDEF (VAR='AK',
IREGN=`KODFOR/100`=10) → `A00F32W###` (R10 Flewelling 2-pt, inside bark; SF/AF/YC/SS/LP/RC/WH/MH),
`A00DVEW###` (DVEST woodland), `A32CURW###` (R10 CUR hardwoods). **akt01 is 100% F32.** DVE/CUR are
not yet ported (return 0 — akt01 has none).

**F32 = the shared engine.** grossvol.f dispatches `MDL='F32'` to the SAME NVEL `PROFILE` (Flewelling)
routine jl already ports (INGY/CR FW2). FWINIT GEOCODE 'A': 042(YC)→JSP31, 242(RC)→32, 098(SS)→33,
260/263/264(WH/MH/LP/SF/AF)→34; `SHP_AK` (f_alaska.f) is byte-for-byte the shared `_fw2_shp_core`
kernel with an AK F-coefficient column (JRSP=JSP−30, hemlock folds onto spruce's F3). So the port
reuses `_fw2_shp_core`/`_fw2_sf_taper`/`_fw2_sf_yhat`/`_fw2_tcubic`/`_fw2_hs`/`_nvb_numlog`/
`_nvb_segmnt`/`_scrib`/`_fw2_dclass` unchanged and adds only the AK F-table
(`data/southeastalaska/volume_coefficients.jl`) + the driver (`src/variants/southeastalaska/volume.jl`).
`COR_AK`/`SF_CORR` is only in the 3-point path — correctly skipped.

**Three MEASURED root-causes closed (instrumented the live NVEL via single-`.o` relinks of `sf_2pt.f`
/ `profile.f`, dumping DBH_IB / F / RFLW/RHFW and LMERCH/MINLEN/NUMSEG/VOL4/VOL2 per tree):**
1. **Bark is the AK VARIANT bark, NOT NVEL FDBT_AK.** FVS's volume driver passes `DBT_USER =
   DBHOB−DBHOB·BRATIO(ISPC)` from `ak/bratio.f`, so SF_SHP's FDBT_AK branch is bypassed. The profile
   is scaled to — and the merch tops converted with — `DBHIB = D·ak_bratio(sp,d)` (bit-exact vs live
   SF2PT `DBH_IB`: SS D5.8→5.16169, LP D11.5→10.93604, …). Using FDBT_AK gave the wrong DBHIB (SS
   5.517) ⇒ the ~1-3% total-cubic error. Fixing it made **total cubic per-tree bit-exact**.
2. **Merch standards** (`setcubicdflts.f`, VARACD 'AK', KODFOR 1005 → AKMERCHCAT 3): DBHMIN=9, TOPD=7,
   SCFMIND=9, SCFTOPD=7, STMP=1. Merch top = TOPD·bark (inside), DBH≥9 gate.
3. **Bucking + board** (`profile.f`): the **Region-10 F3 special** resets `MINLEN=2` when
   `N16SEG=INT(LMERCH/(MAXLEN+TRIM))` is odd (else 8); merch cubic VOL(4) = Σ 0.1-rounded inch-class
   Smalian logs; **board VOL(2) uses the 32-FOOT log rule** (pair 16-ft logs → 32-ft, DIB=INT(raw
   small-end), Σ SCRIB·10) — that 32-ft pairing was the ~1.3× board error. Region-10 NUMLOG/SEGMNT
   params: OPT=23 EVOD=2 MAXLEN=16 MINLEN=8 MERCHL=8 TRIM=0.5.
Also: the total-cubic volume floor is `H>4.5` (not `H≤5`) — a broken-top D8.4/H5.0 tree gets 1.7 cuft
(the H≤5 guard wrongly zeroed it, the −24 cuft aggregate gap).

**Validation vs `FVSak_clean` akt01 cyc0 (TREELIST + .sum):** all **29 trees bit-exact** on total
cubic, merch cubic, and Scribner board (LP 21.4/14.5/60, WH 24.9/17.8/60, LP D9.5 11.6/6.8/20, MH
13.2/10.1/30, YC-dead 240.3/223.9/1010, …). .sum aggregates **MCuFt 732 = live 732** and **BdFt 2417 =
live 2417 bit-exact**; **TCuFt 2316 vs 2315** (±1 = Float32 non-associative summation of the per-acre
0.1-rounded per-tree cuft, below the integer .sum precision — every per-tree value matches). This meets
the CI MCuFt/BdFt exact-at-cyc0 bar. Growth columns unchanged (cyc0 still bit-exact). DVE + CUR
families remain unported (0) — a follow-on chunk; akt01 exercises neither.
