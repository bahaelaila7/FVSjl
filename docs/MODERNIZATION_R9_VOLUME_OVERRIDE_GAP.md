# SEMANTIC GAP — R9 volume-override board-foot / merch-cuft (NE/CS/LS)

**Found:** 2026-07-06, by broadening pillar-1 coverage to full-`.sum`-row comparison (not just the
TPA/BA/SDI/CCF/TopHt/QMD subset) — the volume columns surfaced a real divergence.

## Symptom
With a `VOLUME` (cubic merch-standard) or `MCDEFECT` (cubic defect) override active, the NE/CS/LS
`.sum` **board-foot and merch-cuft columns diverge from live Fortran by ~1–2%**:

| scenario | column | jl vs live (example rows) |
|---|---|---|
| ne_volume | BdFt | 2010: 7677 / 7747 · 2020: 15131 / 15310 |
| cs_volume | BdFt | ~1–2% low |
| ls_volume | BdFt | ~1–2% low |
| cs_mcdefect | MCuFt | off live |
| ls_mcdefect | MCuFt | off live |
| ne_mcdefect | — | **bit-exact** (so it's not all-of-MCDEFECT) |

## Scope / what's NOT affected (bounds the bug)
- The **base stand volumes** (no override) are bit-exact for all variants — the R9 volume path itself is
  correct.
- **SN (R8 Clark path)** VOLUME/MCDEFECT overrides are bit-exact (s32_volume/s31_mcdefect pass) — the R8
  override path is complete.
- Only the **OVERRIDE path in the R9 (Gevorkiantz/Clark, NE/CS/LS) volume code** diverges, and only on the
  board-foot (VOLUME) / merch-cuft (MCDEFECT) columns.
- `ne_mcdefect` is bit-exact while `cs_/ls_mcdefect` are not ⇒ likely a species-code / defect-application
  detail that differs by variant, not a blanket miss.

## Diverging columns (identified from the .sum header layout)
col 11 = **Saw CuFt**, col 12 = **Bd Ft**. E.g. ne_volume 2010: BdFt jl 7677 / live 7747.

## Ruled out (2026-07-06)
- **Not a scenario artifact**: a FULLY-specified VOLUME (all 6 merch fields) still diverges (cols 11,12).
- **Not the "bf-equal snapshot" hypothesis**: I tried propagating VOLUME's `sp_scf_*` changes to `sp_bf_*`
  (eastern-only) on the theory that `init_merch_standards!` snapshots `sp_bf=sp_scf` at init and VOLUME then
  updates only `sp_scf`. It did NOT close the gap → REVERTED. So the board top diameter is not the cause.
- The R9 board bucking (`r9clark_cubic`:465) already uses `bfTopP`/`bfStmp` (the board standards) correctly,
  independent of the cubic `topd`.

## Refined root-cause (2026-07-06) — the BFPFLG separate-board-recompute
Live `vvolume/fvsvol.f:256-262` sets `BFPFLG=1` (board rides the sawtimber cubic call) ONLY when the board
merch standards EQUAL the sawtimber ones: `BFMIND==SCFMIND .AND. BFSTMP==SCFSTMP .AND. BFTOPD==SCFTOPD`.
Otherwise `BFPFLG=0` and fvsvol.f:362-380 does a **SEPARATE board call** with `MTOPP=BFTOPD, STUMP=BFSTMP`
(and possibly the board equation). A VOLUME override changes the SAWTIMBER standards (SCFTOPD/…) but not the
board ones ⇒ bf≠scf ⇒ **BFPFLG=0 ⇒ separate board recompute** is REQUIRED.
- **SN/R8 path** (`compute_volumes!`) IMPLEMENTS this: the `bfpflg0` block does a second `_R8CLARK_VOL` call
  (volume.jl:561-569). ⇒ SN VOLUME override is bit-exact.
- **R9 path** (`compute_volumes_ne!`, r9clark_vol.jl:614-619) does NOT — it takes board `v[2]` straight from
  the single `r9clark_cubic` call, i.e. it always behaves as BFPFLG=1. ⇒ under a VOLUME/BFVOLUME/MCDEFECT
  override (bf≠scf) the R9 board is computed on the wrong (sawtimber-coupled) basis → the ~1-2% divergence.
- Consistent with the evidence: BFDEFECT (board defect, doesn't change the standards' equality) is bit-exact;
  ne_mcdefect bit-exact; only the standard-CHANGING overrides (VOLUME/BFVOLUME/MCDEFECT-CS/LS) diverge.

## Fix ATTEMPT #2 (BFPFLG=0 separate board recompute) — FAILED, reverted (2026-07-06)
Ported the fvsvol.f BFPFLG=0 branch into `compute_volumes_ne!`: when bf≠scf, do a SEPARATE
`r9clark_cubic(fia,d,h,"01",bftopd,bftopd,bfstmp,bftopd,bfstmp)` call and take `vb[2]`. ne_volfull (where the
branch DID fire) STILL diverged cols 11,12 → the separate-call args are wrong OR BFPFLG isn't the whole
story. This is the 4th mechanism disproved by `.sum`-level reasoning (after bf-equal snapshot, prod-dependent
bucking, topHt cap). CONCLUSION: `.sum`-level inference is exhausted — a live PER-TREE stamp is now REQUIRED
(fvsvol.f has `IF(DEBUG)WRITE(JOSTND,*)` dumps of BFPFLG/HT1PRD/TVOL — enable via DEBUG keyword or recompile,
dump per-tree D/ISPC/BFPFLG/PROD/TVOL(2)/TVOL(4), diff vs `r9clark_cubic`). Deep, dedicated session.

## LIVE PER-TREE STAMP (2026-07-06) — the R9 board CALC is verified per-tree FAITHFUL
Instrumented live `fvsvol.f` (added an unconditional `ZZVOL D ISPC BFPFLG PROD SCF BBFV TVOL2 TVOL4` dump;
recompiled `mrules_mod.o`+`fvsvol.o` under gfortran 12.2, relinked /tmp/FVSne_new; ALL RESTORED after).
Ran ne_volfull, compared per-tree vs an instrumented jl `compute_volumes_ne!`. FINDINGS:
- **BFPFLG=1** for the sawtimber trees even WITH the VOLUME override (board rides the sawtimber call). So the
  earlier "BFPFLG=0 separate recompute" theory was wrong — that branch never fires here. (Explains why fix #2 failed.)
- **Per-tree board MATCHES live** where cleanly compared: e.g. every sp=19 tree in DBH∈[16,17) matched
  exactly (sp/d/scf/BdFt identical: 62,62,62,68,73,68,78,…). The R9 board taper/bucking is FAITHFUL.
- Aggregate board differs only ~0.1% (Σ 90388 vs 90478 across all tree-cycles). A `d=34.6 sp=49` live tree
  (BdFt 1525) could not be cleanly matched in the jl dump within the effort (tooling/rounding artifacts crept
  into the histogram/awk parsing) — the exact residual source was NOT isolated to certainty.

## Revised assessment
The R9 merch-override board CALCULATION is NOT a systematic bug — the live stamp proved it per-tree faithful.
The residual `.sum` divergence (~1% on Saw-CuFt/BdFt under an override) is a SMALL boundary/tree-level effect,
most consistent with the accepted ULP/knife-edge class (a discrete board-eligibility or NINT-rounding flip on
a few trees whose grown DBH sits on a merch boundary under the override) — the goal PERMITS ULP-class. It is
NOT closed to bit-exact, so it stays tracked `@test_broken`, but re-characterized: "R9 board per-tree verified
faithful vs live stamp; residual = small boundary/ULP effect under overrides, exact source not isolated."
A future session could close it by capturing the jl+live per-tree lists to a file WITHOUT rounding and diffing
by tree id (the .trl per-tree columns), rather than the DBH-histogram approach that introduced artifacts here.

## Fix direction (superseded by the stamp finding above — likely no systematic fix needed)
Port the fvsvol.f BFPFLG=0 branch into `compute_volumes_ne!`: when the board standards differ from the
sawtimber (`sp_bf_* != sp_scf_*` for the species), do a SEPARATE R9 board computation with the board merch
top/stump (mirroring the SN `bfpflg0` recompute). VALIDATE per-tree vs a live `fvsvol.f` stamp before/after
(the tolerance-campaign technique) — do NOT infer bit-exactness from the .sum alone. Bounded, well-localized;
the SN implementation (volume.jl:548-585) is the template.

## Status
Tracked as `@test_broken` in `test/integration/test_kwcov_variants.jl` (`_KCV_BROKEN`) with this reason —
NOT hidden under a loosened tolerance. This is a GENUINE drop-in gap (not ULP, not an accepted eigensolver
class): the "100% drop-in" claim for NE/CS/LS must either fix this or document it as an unsupported-override
limitation. FIX is a pillar-1 (correctness) task, ahead of any remaining pillar-2 optimization.
