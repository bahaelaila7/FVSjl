# COVER — SHRUB half port (CVBROW→CVSCON→CVSUM-shrub→CVCLAS + report)

Base = `kt-variant-port` HEAD **b04e6801** (canopy half merged).  This agent ran in an
isolated worktree and CANNOT commit.  All work is staged in `scratchpad/cover/`; the
main loop applies the wiring.diff to a real `kt-variant-port` checkout.

## RESULT (validated vs the LIVE oracle by running it)
- **CVSCON scope: MANAGEABLE.** 673-line Fortran but it is *only* DATA coefficient
  tables (PHABOV/PHABUN/PPHYS/PFLOC/PDIST/POTHER site + per-species H/C arrays) plus
  linear index arithmetic → PCON(31)/TCON(2)/HCON(31)/CCON(31).  No branching, no RNG,
  deterministic from SLOPE/ELEV/ASPECT + IOV/IUN/IPHYS/INF/IDIST/ITUN.  ~1500 coeffs
  transcribed verbatim (Fortran column-major DATA → `reshape(vec,N,31)` = TAB[cat,sp]).
- **Kernels dump-replay BIT-EXACT vs FVSem_g16** (DEBUG dump, EM stand S248112, SHRUBS 3.0):
  - CVSCON: **64/64** (TCON(1)=1.71096, TCON(2)=3.65430, HCON[1..31], CCON[1..31]) to the
    F10.5 dump precision.  `replay_shrub.jl`.
  - CVBROW: **164/164** — cyc0 per-species SH/CV/PB/CABHT/PBCV (31×5=155) + PGT0/TCOV/
    TOTLCV, and cyc1/cyc2 PGT0/TCOV/TOTLCV.  `replay_cvbrow.jl`.
  - CVSUM-shrub + CVCLAS: TALLSH/CLOW/CMED/CTALL/TOTLCV/ASHT + the full SHRUB STATISTICS
    text table = **30/30** vs `em_cov.out`.  `replay_report.jl`.
- **End-to-end** `run_keyfile("em_cov.key"; variant=EasternMontana())`:
  - SHRUB STATISTICS table: **byte-exact vs oracle at cyc0 (1990)**; cyc1/cyc2 differ by
    ±0.1 in a handful of F6.1 cells (e.g. VASC prob 8.5 vs 8.4) — the documented EM #206
    OLDRN growth-BA straddle feeding CVBROW's BA (cyc0 BA=85.1313 matches exactly; the
    kernels are proven correct at exact BA, so this is upstream-growth cornered, NOT a
    shrub-kernel bug).
  - CANOPY AND SHRUBS SUMMARY: **understory columns (TIME/PROB/LOW/MED/TALL/TOTAL/AVG
    SHRUB HT/BIOMASS/TWIGS/SUCC STAGE) match the oracle on ALL 3 cycles** (integer
    rounding absorbs the ±0.1).  Overstory STAND AGE/TOP HT/CANOPY COVER/STEMS match at
    cyc0; FOLIAGE BIOMASS + STEM DIAMS carry the **pre-existing canopy-half CVCBMS/SDIAM
    residual** (out of scope for the shrub port — canopy TREES row is bit-exact, proving
    my additive sd2xht/htmax/htmin did not perturb the canopy tallies).
- **INERT:** multicycle regression **339/11** byte-identical; cover unit tests all green.
  A run without the COVER keyword never constructs CoverState → growth untouched.

## THE ORACLE (no relink; DEBUG dump recipe)
`/workspace/.emwork/FVSem_g16` fires COVER on the keyword.  SCREEN routes the *summary*
to stdout; the **detailed report is written to `fort.16`** (JOSTND) in the run dir.
DEBUG dumps also go to fort.16.  Recipe used (built-in DEBUG, no recompile):
```
DEBUG            0.0       1.0     <-- field1 cols 11-20 (=cycle, 0=all), field2 cols 21-30 non-blank
CVSCON CVBROW CVCLAS            <-- one line, space-separated routine names (DBPRSE reads ONE line)
```
then `echo em_dbg.key | FVSem_g16 >/dev/null; cp fort.16 em_dbg_full.out`.
Goldens captured in `em_dbg_full.out`; parsed into the replay harnesses.

## STAGED FILES (scratchpad/cover/)
- `cover.jl.merged`     — the FULL merged src/engine/cover.jl (canopy + shrub).  APPLY = copy this.
- `wiring.diff`         — `git diff b04e6801 -- src/engine/cover.jl` (apply with `git apply`).
- kernel sources (already concatenated into cover.jl.merged, kept for review):
  `cover_shrub.jl` (CVSCON), `cover_shrub_cvbrow.jl` (CVBROW), `cover_shrub_cvsum.jl`
  (CVSUM-shrub + CVCLAS), `cover_shrub_wire.jl` (habitat tables + ShrubRow),
  `cover_shrub_report.jl` (SHRUB STATISTICS + SUMMARY emit).
- replay harnesses: `replay_shrub.jl` (CVSCON), `replay_cvbrow.jl` (CVBROW),
  `replay_report.jl` (SHRUB STATISTICS vs em_cov.out).  Run each with
  `JULIA_DEPOT_PATH=/workspace/.julia_depot julia <file>` (they `include("quickersort.jl")`).
- oracle goldens: `em_cov.golden`/`em_cov.out` (full report), `em_dbg_full.out` (DEBUG),
  `jl_cov.out` (FVSjl end-to-end output), `em_cov.key`/`em_cov.tre`, `em_dbg.key`.

## APPLY (from a real kt-variant-port checkout at b04e6801)
```
cp scratchpad/cover/cover.jl.merged  src/engine/cover.jl
#   OR:  git apply scratchpad/cover/wiring.diff
JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project -e 'using FVSjl'   # loads clean
```
No other file changes are needed — the report seam (keyword_dispatch.jl:2448 kw_coverin!,
io/summary.jl:229 cover_accumulate!, simulate.jl:876 cover_report) was already wired for
the canopy half and the shrub half hooks into the SAME three seams.

## WHAT THE MERGE ADDED TO cover.jl (all inside src/engine/cover.jl)
1. CoverState: 6 new fields (sage0, ihtype, iphys, idist, shrub_const::Any, shrub_rows::Vector)
   with defaults (-1,0,2,1,nothing,[]); constructor + cover_init! updated.
2. kw_coverin! SHRUBS branch: capture fields 1-4 (SAGE/IHTYPE/IPHYS/IDIST).
3. cover_accumulate!: accumulate sd2xht/htmax/htmin in the canopy tree loop (additive,
   INERT for canopy), then `if cv.lbrow` call `_cover_shrub_cycle!`.
4. `_cover_shrub_cycle!`: SAGE = max(sage0<0 ? IAGE : sage0, 3) + (year-year0); CVSCON at
   cycle 0 (habitat→IOV/IUN/ITUN via `_cvb_habitat`, INF via `_cvb_inf(EM=3)`); cvbrow_cycle;
   cvsum_shrub; cvclas; Irwin-Peek SBMASS/TWIGS only if IHTYPE∈{520,530,570}; push ShrubRow.
5. cover_report: gate each table independently (has_canopy/has_shrub/has_sum); add
   `_cover_shrub_stats` (SHRUB STATISTICS) + `_cover_summary` (CANOPY AND SHRUBS SUMMARY).
6. `_cover_stand_habitat(s)`: EM stores habitat_code = IEMTYP (index into EM_JTYPE), so the
   raw ICL5 code (260) is recovered via `EM_JTYPE[iemtyp]`; other variants keep the raw code.
   (Root cause: easternmontana/site_index.jl:166 overwrites habitat_code with IEMTYP.)

## KEY FACTS / GOTCHAS
- Transcendentals: `_f32log`/`_f32exp` = `ccall((:logf/:expf,"libm.so.6"),Float32,...)` —
  the fmath shim is 1-ULP off; libm matches Fortran REAL EXP/ALOG.
- Fortran REAL = Float32 throughout; report fields are IFIX(.5+x) or F6.1/F7.1 so 1-ULP
  kernel diffs never surface.
- RDPSRT (order-sensitive for CABHT & top-3 selection) reuses FVSjl's faithful
  `rdpsrt!(n,a,index,lseq)` (src/engine/quickersort.jl) — NOT sort!/sortperm.
- SHRUBS card fields (cvin.f opt 3): 1 SAGE, 2 IHTYPE (0⇒stand ICL5), 3 IPHYS (default 2),
  4 IDIST (default 1).  LCALIB defaults FALSE (only SHRBLAYR/SHRUBHT/SHRUBPC set it) ⇒
  CVBCAL is skipped and RESIDC/RESIDH=0 (this port's no-calibration path).
- IRWIN&PEEK biomass/twigs = 0 unless IHTYPE∈{520,530,570}; matches oracle (0/0.0).
- Julia `"""..."""` AUTO-DEDENTS — the header text is emitted as per-line string literals
  (`_CV_SHRUB_HDR` etc.), NOT triple-quoted blocks, to preserve exact column alignment.
- The GROHED "FOREST VEGETATION SIMULATOR VERSION ... <timestamp>" page-header line is NOT
  emitted (it carries a wall-clock timestamp); the parity check excludes it (same as canopy).

## NEXT (dep-ordered)
1. **CVBCAL calibration** (covr/cvbcal.f): SHRBLAYR (LCAL1, by layer) + SHRUBHT/SHRUBPC
   (LCAL2, by species) → BHTCF/BPCCF scaling + RESIDC/RESIDH terms in CVBROW.  The
   SH/CV equations here DROP the `+coef*RESIDH(k)` / `+coef*RESIDC` terms (0 without
   calibration); restore them (see cvbrow.f) when wiring CVBCAL.  Add the calibration
   report tables (cvout.f 1000/1050 blocks).  Validate with a SHRBLAYR keyfile.
2. **SHRUB-SMALL CONIFER COMPETITION table** (cvout.f 9090+): SCOV(11)+TRSH(11) by height
   threshold.  SCOV is already computed in ShrubSummary.scov; TRSH (trees/ac by height)
   needs the tree loop's HT-vs-SHTRHT tally (cvsum.f:326-333, 531-537) threaded into the
   canopy loop — cheap add.  Emit deferred (KODE≠3 gate).
3. **Post-thin ITHN=2** path (LTHIN): CVBROW recomputes with post-thin BA + SAGE-reset-to-3
   test (CRAREA/43560 & OCVREM/OCVCUR).  Deferred with the thin seam (same as canopy).
4. **Non-EM habitat raw-code recovery**: `_cover_stand_habitat` handles EM; other western
   variants (IE/CI/KT/…) that overwrite habitat_code need their own reverse-map (or capture
   the raw ICL5 at STDINFO parse into a dedicated plot field).  INF dispatch (`_cvb_inf`)
   covers EM=3, CI (IFOR==1?2:1), IE (INFOR[IFOR]); other variants = Boise/Payette group 1.
5. **Canopy FOLIAGE BIOMASS / SDIAM residual** (pre-existing, canopy half): CVCBMS is
   337/351 @1-ULP but the per-class totals diverge ~7% (jl 4986 vs oracle 5367) — belongs
   to the canopy chunk, surfaces in the SUMMARY overstory column, worth a separate look.
