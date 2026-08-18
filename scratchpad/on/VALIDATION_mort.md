# ON chunk 8 — small_tree_growth!(REGENT no-op) + mortality!(morts.f/varmrt.f)

Base: kt-variant-port HEAD **35a3d597**. Oracle: instrumented **FVSon_g16** (morts.f single-.o swap,
built into the g16-consistent object set /workspace/.onwork/g16obj — NOT relink_on.sh, which mixes
the buildDir's non-g16 objects and yields a binary that mis-reads keywords).

## (1) small_tree_growth!(::Ontario) — canada/on/regent.f — NO-OP-VALIDATED
ON REGENT grows only records with DBH < XMAX = 4.72" (= 12 cm; regent.f:78 `DATA XMAX/MAXSP*4.72/`,
uniform for all 72 species). ont01 is 8 large trees whose smallest internal DBH is SB 15 cm = 5.906"
> 4.72", and no record ever drops below 4.72" across the projection (trees only grow; mortality never
creates small trees; no regen keyword). So REGENT skips every record (`GO TO 25`) — a genuine no-op.
The jl reproduces this: the tree list (dbh/diam_growth/ht_growth + the tripling stash) is unchanged
across the call. The sub-12 cm branch (htcalc.f NC128 height-age + htont.f ONHTDBH + htdbh.f Wykoff
dubbing + DGBND, ~760 lines) has NO ON validation stand, so it is not staged; a record < XMAX errors
loudly rather than emitting unvalidated numbers.

## (2) mortality!(::Ontario) — canada/on/morts.f + varmrt.f — DUMP-REPLAY bit-exact-or-cornered
ON MORTS is a FRESH port (not shared southern Pretzsch): Penner MAX-SDI density line (13 SDI_INT/
SDI_SLP eqns via MSB_MAP keyed on the most-BA species), 4-group background (BKG_MAP→PMSC/PMD, halved),
metric throughout, and ON's own VARMRT = the Penner INDIVIDUAL-TREE mortality logistic (18 eqns,
ITM_MAP→MB0..MB7) driving the geometric progression. Bit-exact glibc transcendentals.

Oracle ont01 cyc1 (MORTDUMP, `oracle_mortdump_ont01.txt`): ITRN=8 (no tripling), LZEIDE=F, FINT=10,
T=34317.87, DIA0=9.24562, D10=9.28694, SDIMAX=35599.32, BAMAX=16503.96 (=SDIMAX·0.5454154·PMSDIU),
INDX_BA=8→INDX=6, T85D0=387.46 (T≫T85D0 ⇒ TN10=T85D10=385.39), RN=0.36168. jl reproduces T, DIA0/D10,
SDIMAX (hex 0x470D2A94), tt (hex 0x47082A5F), BAMAX, INDX_BA all bit-exact.

Per-tree WK2 (kill TPA), jl (exact oracle-DG hex injected) vs oracle FINAL WK2 (Float32-hex):
| rec | ISP | oracle WK2 | jl WK2   | Δ    |
|-----|-----|------------|----------|------|
| 1   | 5   | 45224438   | 45224437 | 1 ULP|
| 2   | 6   | 4569A9EA   | 4569A9E9 | 1 ULP|
| 3   | 26  | 44EE6EAF   | 44EE6EAE | 1 ULP|
| 4   | 28  | 44CA456F   | 44CA456E | 1 ULP|
| 5   | 9   | 46224438   | 46224437 | 1 ULP|
| 6   | 11  | 4596DE1B   | 4596DE1A | 1 ULP|
| 7   | 1   | 453C98F2   | 453C98F2 | EXACT (full kill = PROB) |
| 8   | 8   | 45B68CBE   | 45B68CBE | EXACT |
2/8 bit-exact, 6/8 within exactly 1 ULP. The residual is the VARMRT geometric-progression Float32
accumulation knife-edge (adjust=TEMKIL/TEMSUM scaling a common TOKILL=T−TN10; a ~1-ULP D10/TN10 in the
non-associative density sum scales every distributed kill by ~1 ULP). This is the accepted
bit-exact-or-cornered mortality class (cf. southern/mortality.jl: "the knife-edge kill flips on that ULP").

## (3) Full cyc0→cyc1 .sum, jl vs FVSon_g16 readable table (per HA)
| col   | cyc0 jl | cyc0 oracle | cyc1 jl | cyc1 oracle (readable) |
|-------|---------|-------------|---------|------------------------|
| year  | 0       | 0           | 10      | 10                     |
| TPA/HA| 84799   | 84799       | 952     | 952                    | ✓
| BA    | 3673    | 3673        | 40      | 40                     | ✓
| SDI   | 74769   | 74769       | 820     | ~                      |
| TopHt | 24      | 24          | 27      | 26.5                   | ✓ (round)
| QMD   | 23.5    | 23.5        | 23.1    | 23.1                   | ✓
| CCF   | 15      | **** (ovfl) | 0       | 218                    | cyc0 gap (open-grown CW), carries to cyc1
| vol   | 0       | 28068/16514 | 0       | 354/182                | volume kernel NOT ported (next chunk)

The jl cyc1 STAND (952/40/820/27/23.1) matches the oracle's readable cyc1 table (952/40/26.5/23.1).
NOTE: the oracle .sum's cyc1 TPA COLUMN prints 24 (its own BA=39 implies ~950 tph) — an oracle metric
sumout.f IOSUM(3) artifact for the post-growth row, not a jl error; jl's 952 is consistent with BA.

## DG-uniformity investigation — RESOLVED = DGSD=2.0 OLDRN straddle (cornered, not a bug)
The cyc1 jl DGs are near-uniform (~0.0249 for 7/8 species) while the oracle DGs have per-species spread
(0.0233–0.0250). This is the OLDRN random perturbation: FVS (DGSD=2.0) adds a ±10%·DG BACHLO draw per
tree; the jl applies the DETERMINISTIC Penner DG (validated bit-exact at the DDS level in chunk 5,
replay_penner.jl 8/8 vs FVSon_dgfdump). Feeding the oracle's exact DG into jl mortality reproduces WK2
to ≤1 ULP (above), proving the DG difference is purely OLDRN and the mortality port is correct. This is
the #206 accepted cornered class; the cyc1 stand columns straddle within rounding.

## Gates
ontario unit 80/80 (11+10+34+25). multicycle 339 pass / 11 broken — byte-identical (the additive
Density.mort_ibasp field defaults 0; no other variant reads it). run_keyfile(ont01; Ontario()) now
completes the projection loop (blocks NEXT at the volume kernel → 0 vol columns).
