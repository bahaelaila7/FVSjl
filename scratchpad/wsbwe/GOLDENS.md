# WSBWE (Western Spruce Budworm) — scope, goldens, oracle, verdict

6th insect/pathogen model (after DFB/DFTM/WPBR/WWPB/LPMPB). Model source: pristine
`/workspace/ForestVegetationSimulator/wsbwe/*.f` (57 files, internal prefix `bwe`).
FVSjl target branch `kt-variant-port` HEAD 74935cd7.

## ARCHITECTURE (measured) — WSBWE is STAND-LEVEL and end-to-end relinkable

Audit-confirmed **no PPE dependency** in the reachable call graph (grep of the
wsbwe call graph for ALSTD2/PPMAIN/SPLAEX/GPNEW = none; the 3 PPE handoff files
bweppatv/bweppgt/bwepppt are excluded from the variant link and nothing in base
calls them). This is UNLIKE PPE-gated WWPB — a single-stand oracle relinks and runs.

Keyword entry: base keywds.f **option 8 = 'WSBW'** (`base/keywds.f` TABLE(8)='WSBW';
`initre.f:672-680 CALL BWEIN(LKECHO)`). `BWEINT` (init) is called unconditionally
from `initre.f:186`. The extension has its OWN 25-token sub-keyword TABLE (loaded
in each `bwebk<v>.f` BLOCK DATA):

```
END OPEN CLOSE MGMTID RANNSEED STARTYR COMMENT DEFOL SETPRBIO DAMAGE NODAMAGE
PERDAM NOPERDAM GENDEFOL OUTBRLOC RECOVERY WEATHER NEMULT BWOUTPUT PARASITE
FQUALDEV FQUALWT TITLE BWSPRAY OBSCHED
```

Per-cycle gate **BWEGO** (bwego.f) — called from `base/grincr.f:414 CALL BWEGO(LBWEGO)`,
gate consumed at `base/gradd.f:108 IF (IPMODI.EQ.1 .AND. LBWEGO)` — EXACTLY the DFB
LDFBGO pattern. LBWEGO fires when `LDEFOL` (manual DEFOL) OR `LCALBW/LBUDL` (BUDLITE
regional outbreak active this cycle). The annual within-cycle driver is **BWEDR**
(bwedr.f) → BWEDEF (defoliation) → BWEDAM/BWEDIE (growth-loss, top-kill, mortality).

### ACTIVATION vs OUTPUT-ONLY (bwein.f)
- **DEFOL** (opt 8): `LDEFOL=.TRUE.` — user-supplied annual defoliation. **DETERMINISTIC**
  (NO RNG, NO weather). Schedules OPNEW activity 2151. `bwego.f`: LDEFOL ⇒ LCALBW=.FALSE.
- **GENDEFOL** (opt 14): `LBUDL=.TRUE., IBWCHK=1` — the BUDLITE population model.
  **STOCHASTIC** (draws BWERAN) and needs a weather source (WEATHER keyword or model
  block-data weather); without one it dies at `bwein.f:312` reading fort.40 (EOF).
- SETPRBIO(9)=OPNEW 2153, BWSPRAY(24)=OPNEW 2157/2159, OBSCHED(25)=outbreak windows.
- OUTPUT/CONFIG (no growth effect): OPEN CLOSE MGMTID RANNSEED STARTYR COMMENT
  DAMAGE NODAMAGE PERDAM NOPERDAM OUTBRLOC RECOVERY WEATHER NEMULT BWOUTPUT PARASITE
  FQUALDEV FQUALWT TITLE.

### HOST species map IBWSPM (per-variant, bwebk<v>.f) — budworm defoliation classes
1=WF 2=DF 3=GF 4=AF(subalpine fir) 5=ES(Engelmann spruce) 6=WL(western larch) 7=NON-HOST.
Host variants with own coeffs: **BC BM CI EC EM SO TT** (+ generic bwebk.f/bwebms.f).
EM (19 spp): `IBWSPM = [7,6,2,7,7,7,7,5,4,7,7,7,7,7,7,7,7,7,7]` (WL→6, DF→2, ES→5,
AF→4; NHOSTS=6, NCROWN=9). Each other host variant's IBWSPM must be extracted from
its bwebk<v>.f for the effect seam.

## RNG — BWERAN (bweran.f), IDENTICAL LCG to LPMPB/DFTM/WWPB/DFB

`S1 = DMOD(16807D0·S0, 2147483647D0); SEL = S1/2147483648D0; S0 = S1` (MINSTD
a=16807, m=2^31-1, divisor = **2^31 = 2147483648**). Default seed **55329** set
DIRECTLY (`DATA S0/55329D0/` in bweran.f + every bwebk<v>.f) — like LPMPB's MPRANN,
NO +1128 offset. Entries: BWERSD (reseed: LSET=true forces odd + stores SS/S0;
LSET=false restarts from SS), BWERGT/BWERPT (get/put S0). Consumers (BUDLITE only):
BWEBET (beta draws), BWERNP (random defoliation given a mean). The MANUAL DEFOL
path draws NO random numbers.

### BWERAN golden — seed 55329, first 8 draws (Float32 hex) — VALIDATED bit-exact
`scratchpad/wsbwe/driver_bweran.f` + pristine `wsbwe/bweran.f`, gfortran-16:
```
3EDDB57A 3F5AB21B 3F630A07 3F2734C9 3EF4FB2C 3F4AF646 3F6E7B68 3F67EA59
```
IDENTICAL to WWPB(BMRANN)/LPMPB(MPRANN)/DFTM(TMRANN) — same LCG + seed. BWERGT→
BWERPT save/restore reproduces the stream exactly (drv confirms draw2-after-get/put
== draw2). Julia `wsbwe_rand!` reproduces all 8 BIT-EXACT (verified standalone +
test_wsbwe.jl).

Reproduce:
```
cd scratchpad/wsbwe && gfortran-16 -std=legacy -w -fno-automatic \
  -I/workspace/ForestVegetationSimulator/wsbwe driver_bweran.f \
  /workspace/ForestVegetationSimulator/wsbwe/bweran.f -o drv_bweran && ./drv_bweran
```

## ORACLE — FVSem_wsbwe (relinked, RUNS; off ≡ stock byte-identical)

`scratchpad/wsbwe/build_em_wsbwe.sh` (mirror of scratchpad/dftm/build_ie_dftm.sh):
FVSem_buildDir *.o with `base/exbudl.f` no-op stub DROPPED + the real wsbwe/*.o.
File set = authoritative `canada/bin/FVSbcc_notReady.txt` wsbwe block (39 files) with
BC host coeffs swapped to EM (bwebkbc→bwebkem, bwebmsbc→bwebmsem); EXCLUDES generic
bwebk/bwebms, other-variant bwebk*/bwebms*, the 3 PPE files, and txnote. Links clean
(compiled_ok=652 fell_back=2 ; wsbwe compiled=39; LINK_OK, 9.3 MB).

VALIDATION:
- **WSBW-off ≡ stock**: `FVSem_wsbwe em_off.key` vs stock `bin/FVSem em_off.key` →
  `.sum` BYTE-IDENTICAL (`cksum 534720991`, both 1723 bytes). The relink is faithful.
- **Keyword flow live-confirmed**: `em_on.key` (WSBW/GENDEFOL/OBSCHED/END) parses &
  schedules correctly (echoed in `.out`: GENDEFOL "CURRENT OUTBREAK", OBSCHED "START=
  1995, END=2005"). GENDEFOL then dies reading fort.40 (no weather source) — expected;
  BUDLITE needs a WEATHER setup.
- **DEFOL deterministic path runs** (`em_defol.key`, WSBW/DEFOL DF+GF/DAMAGE/END):
  exit 0, `.sum` differs from off ONLY by timestamp (no defoliation effect on this
  LP/PP-dominant stand — little DF host + values/timing not yet tuned). The effect
  path is DEFERRED; a heavy true-fir/spruce host stand is needed to exercise it.

## PORTED + VALIDATED (this beachhead; scratchpad/wsbwe/wsbwe.jl + test_wsbwe.jl)
- **wsbwe_rand! / wsbwe_seed!** — BWERAN LCG + BWERSD reseed, bit-exact (8-draw
  golden + odd-force + LSET=false reset).
- **wsbwe_defaults!** — BWEINT/BLOCK-DATA defaults (seed 55329, IOBLOC=2, IWOPT=1).
- **kw_wsbwe!** — the `WSBW … END` block reader (base keywds.f opt 8), 25 sub-keywords;
  captures activation state (LDEFOL/LBUDL + DEFOL/OBSCHED schedules + seeds).
- **INERT seam** — a WsbweState projects .sum-byte-identically to a no-WSBW run (no
  engine seam wired; the defoliation effect path is deferred). Validated A/B in test.

## VERDICT
RNG bit-exact; oracle relinkable + off≡stock byte-identical; keyword reader faithful
& inert. The deterministic effect path (BWEDR→BWEDAM/BWEDIE) is the next chunk and is
dump-replay-validatable against FVSem_wsbwe on a heavy host stand (see HANDOFF.md).
Expected end-to-end .sum verdict once wired: bit-exact-or-CORNERED by the EM #206
growth straddle (mirrors DFB/LPMPB).

## Files (scratchpad/wsbwe/)
- `wsbwe.jl` — staged src/engine/wsbwe.jl (RNG + reader + defaults; inert).
- `test_wsbwe.jl` — staged test/unit/test_wsbwe.jl.
- `driver_bweran.f` — BWERAN golden generator (reusable; standalone, no commons).
- `build_em_wsbwe.sh` — FVSem_wsbwe oracle relink (regenerates the binary + emobj/).
- `em_off.key / em_on.key / em_defol.key / s248.tre` — oracle test keyfiles.
- `HANDOFF.md` — the 4-seam apply/wire/commit recipe + dependency-ordered next chunks.

## EFFECT KERNEL — DEFOLIATION GROWTH-LOSS + MORTALITY (dump-replay validated)

Next chunk after the beachhead. STAND-LEVEL manual-DEFOL path
(BWEGO→BWECUP→BWESIT→[year loop: BWEDR→BWEDEF→BWEAGE→BWEDAM]→BWEPDM). The two
arithmetic payload kernels are ported into `wsbwe.jl` and validated BIT-EXACT via
Float32-hex dump-replay against the instrumented live oracle FVSem_wsbwe.

### Host stand that FIRES the effect (measurement gate PASSED)
`host.tre` = emt01.tre relabelled DF/S(ES)/AF (the EM budworm hosts; EM has no WF/GF),
read via `OPEN 3 / host.tre / TREEDATA 3`. `host_defol.key` = sustained DEFOL
85/80/70/60% on DF/ES/AF, years 1990-1999. End-to-end `.sum` DELTA (off → defol):

| year | off TPA/BA | defol TPA/BA |
|------|-----------|--------------|
| 1990 | 536 / 77  | 536 / 77 (start) |
| 2000 | 528 / 102 | **52 / 8** |
| 2050 | 459 / 206 | 48 / 31 |
| 2090 | 394 / 264 | 45 / 54 |

(The prior "S248112 shows no effect" was a FIXTURE bug: `TREEDATA\n<file>` without
OPEN yields a 0-TPA degenerate stand — every column zero. `OPEN 3` fixes it.)

### Instrumentation faithfulness
Single-.o swap of `bwedam.f` / `bwepdm.f` (compiled LAST, hex TRANSFER dump to a
private unit 78/79 side file). Instrumented `.sum` == clean-relink `.sum`
BYTE-IDENTICAL (verified `tail -n +2` diff empty) — the dumps are trustworthy.

### Kernel 1 — BWEDAM Nichols growth-loss (bwedam.f:184/188) — `wsbwe_rdds`/`wsbwe_rhtg`
`RDDS = .083861·(AVPRBO_whole·100)^(.4725+.07·RDDSM1)·RDDSM1^.3241`  (cap 1)
`RHTG = .193013·(AVPRBO_top·100)^(.3814−.0212·RHTGM1)·RHTGM1^.5509`  (cap 1)
**90/90 rows bit-exact** (`bwedam_dump.txt`, DF/ES/AF × 3 size classes × 10 years).
Sample goldens (avprbo_top, avprbo_whole, rddsm1, rhtgm1 → rdds, rhtg):
```
3E199998 3E199998 3F800000 3F800000 -> RDDS 3EBA93DE  RHTG 3F030E28   (yr1 carry=1)
3E199998 3E199998 3EBA93DE 3F030E28 -> RDDS 3E6E7E17  RHTG 3EBA67EF   (yr2 carry≠1)
```

### Kernel 2 — BWEPDM Marsden mortality logistic (bwepdm.f:640) — `wsbwe_mort_pr`
`PR = 1/(1+exp( B0+B1·ELEV+B2·ELEV²+B3·PNTBA+B4·PNTHBA+B5·MFT+B6·KTK+B7·MFT·MFM ))`
(host-indexed B0..B7; PR=0 when MFT or MFM ≤ 0, bwepdm.f:645). This is the survival/
mortality probability BEFORE the deterministic period-scaling `PR^IBWYR·(1−BASE)^FA`
and the `max(BASE,PR)` background composition (both in the next chunk with the WK2
apply). **81/81 rows bit-exact** (`bwepdm_dump.txt`). Sample goldens
(ih, elev, pntba, pnthba, mft, mfm, ktk → pr):
```
2 42580000 400A809A 400A809A 41200000 41200000 00000000 -> PR 3E402809  (DF, ELEV=54)
5 42580000 400A809A 400A809A 41200000 41200000 00000000 -> PR 3E402809  (ES)
4 42580000 3FDB330C 3FDB330C 41200000 41200000 00000000 -> PR 3E879871  (AF)
```
(MFT=MFM=10.0=0x41200000 = total missing foliage under this heavy DEFOL; KTK=0 at
cycle 1 since no prior topkill.)

### Seam (wired, INERT)
`wsbwe_go(w)` = BWEGO manual-DEFOL gate (active ∧ ldefol ∧ scheduled). `wsbwe_apply!`
mirrors mpb_apply! and early-returns (feeder chain deferred). INERT proven at the
oracle: off ≡ stock byte-identical; DEFOL-on-non-host (`nonhost_*.key`, all-LP) ≡ off
byte-identical.

### VERDICT (this chunk)
Host stand fires; effect path fully exercised end-to-end by the oracle. The two
deterministic payload kernels (growth-loss + mortality logistic) replay bit-exact.
Seam wired + inert-proven. Remaining to lift inertness = the FOLIAGE/PRBIO feeder
(BWESIT/BWEBMS/BWEAGE) + BWEPDM per-tree apply (draws the damage RNG via
BWERNP/BWEBET) — see HANDOFF.md "NEXT CHUNKS" 1-5.

## FULL DEFOL APPLY — feeder + BWEPDM per-tree + BWERNP/BWEBET RNG (dump-replay bit-exact)

apply! is now LIVE (EM). Every routine of the manual-DEFOL BWECUP chain replays
BIT-EXACT vs instrumented FVSem_wsbwe on host_defol.key (cycle 1, 81 trees, 11 pts;
instrumented .sum ≡ clean relink verified first). Reproduce: build with `ins/`
overrides → `bwe_{sit,age,dam,pdm}.txt` + `bwe_rng.txt`; then `replay_bwepdm.jl`,
`replay_feeder.jl`+`replay_yearloop.jl`, and `test_staged.jl`.

| routine | what | result |
|---|---|---|
| BWEBMS | WK4 per-tree foliage biomass (Moeur, ICVOPT=2) | 81/81 |
| BWESIT | FOLPOT/FOLADJ/BWTPHA/FNEW/FOLD1/FOLD2/FREM/HOSTST | bit-exact |
| BWEAGE | per-year PRBIO + TOTR (10 yrs) | bit-exact |
| BWEDAM | per-year AVPRBO/PEDDS/PEHTG/AVYRMX/CDEF (10 yrs) | bit-exact |
| final | PEDDS/PEHTG/AVYRMX/PRBIO (feed BWEPDM) | bit-exact |
| BWEPDM | topkill+DG+HTG+WK2 + ICR/ITRUNC/NORMHT/IMC/KTK/PR | 81/81 |
| BWERNP/BWEBET | damage RNG stream (seed DSEEDD=55329) | 519/519 draws |

Per-tree draw ORDER (each host tree): BWERNP(AVYRMX,.06)[beta] → BWERAN(topkill) →
BWERNP(PEDDS,.03)[beta]; HTG XH = deterministic branch (PEDDS<0.99, no draw).

Compact unit goldens (test_wsbwe.jl, glibc math):
```
wsbwe_rdds(3E199998,1.0) = 3EBA93DE   wsbwe_rhtg(3E199998,1.0) = 3F030E28
wsbwe_mort_pr(DF,ELEV54,MFT10,MFM10,KTK0) = 3E402809
seed 55329: wsbwe_bernp(0.85,.06)=3F609607 ; next wsbwe_rand!=3F630A07 (topkill roll)
```

⚠ glibc transcendentals: gfortran ALOG/EXP/`**`→glibc logf/expf/powf; Julia native
Float32 math drifts ~1 ULP ~25% of the time. WSBWE routes all transcendentals through
`ccall(:...,"libm.so.6")` (`wsbwe_log/exp/pow`) + `wsbwe_powi` (pow_r4_i4) for
`PR**IBWYR`. Verified 0/60000 mismatch vs gfortran. The RNG LCG stays pure Julia.

Oracle end-to-end target (host_defol.sum vs host_off.sum, EM) — for the main-loop run:
| year | off TPA/BA | defol TPA/BA |
|---|---|---|
| 2000 | 528/102 | 52/8 |
| 2050 | 459/206 | 48/31 |
| 2090 | 394/264 | 45/54 |
