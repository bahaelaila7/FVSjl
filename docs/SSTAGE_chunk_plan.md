# SSTAGE (stand structural stage) — chunk plan

`SSTAGE` (base/sstage.f, 942 lines) classifies the stand each cycle into a **structural stage
class** (1=SI stand-initiation, 2=SE stem-exclusion, 3=UR understory-reinitiation, 4=young multi-
strata, 5=old single-stratum, 6=old multi-stratum/continuous), plus the per-stratum "Structural
statistics" report. Activated by the **main `STRCLASS` keyword** (ksstag.f) and also called from
FFE each cycle (grincr.f:342, after CUTS). It is a self-contained natural-process classification
(tree list + crown cover + thresholds), so it slots into the C4/output layer.

## ★ It IS validatable against the stripped ground-truth binary
The main `STRCLASS` keyword makes the rebuilt `/tmp/FVSsn_new` print the **Structural statistics**
report to the `.out` — per cycle, per Rm-code (0/1 = before/after thin): for up to 3 strata the
DBH / Nom-Lg-Sm height / Basal / Cover / Sp1 / Sp2 / Crown-diff, then **N-Strata, Tot-Cov, and the
Struct Class** (e.g. `1990 0 … 2 82 3=UR`). So unlike the FFE/Carbon DBS tables (which this stripped
binary can't emit — see [[fvsjl-ground-truth-binary-limits]]), SSTAGE can be diffed column-by-column
against ground truth. **Validate every chunk against this report**, not just the final class.

## Algorithm (sstage.f), most-upstream first
1. **Build the working tree list** (sstage.f:166-274): live trees with HT>0; sort INDEX by HT
   descending; per-tree temp PROB (`TMPPRB`), ICR, DBH, HT, ISP. Thresholds come from the STRCLASS
   keyword / defaults: `TMPGAP` (gap %, the stratum-break height drop), `TMPCCM` (min stratum cover
   %), `TMPTPA` (min TPA to be a stand — default 200), `TMPSSD` (small-tree DBH), `TMPSAW`
   (sawtimber DBH), `TMPPCT` (% of MaxSDI for SE, default 30), `XBAMAX`.
2. **Height-gap stratification** (sstage.f:300-388): walk the HT-sorted list; a GAP is where
   `HT(small) < HT(large) − max(10, HT(large)·TMPGAP·0.01)`, skipping insignificant trees (running
   PROB < 2.0). Track the **two largest gaps** (DIFF1/DIFF2 → boundary index pairs ID1I1/ID1I2,
   ID2I1/ID2I2); swap so the upper gap is on top. ⇒ up to 3 potential strata.
3. **Cover per stratum** (sstage.f:430-465): `COVOLP` (base/covolp.f) on each stratum's trees
   (including the gap trees): `PCCU = CCCOEF·(Σ crown_area/43560)`, `cover = 100·(1−exp(−PCCU))`.
   A stratum is "OK" if cover > `TMPCCM`. `NSTR` = count of OK strata; if 0 and TPA≥TMPTPA, form one
   stratum of all trees. (CCCOEF is the CCADJ coefficient — already in `control.cc_coef`.)
4. **Per-stratum stats** `SSTGHP` (sstage.f:740+): PROB-weighted mean DBH, the nominal/large/small
   heights, mean crown base, and the 1st/2nd dominant species (by stratum basal area). Used for the
   report AND to pick the **dominant stratum's mean DBH `TMPDBH`** (the top OK stratum).
5. **Classify** (sstage.f:539-576) → `TMPSCL`:
   - NSTR=1: DBH<SSD→1; DBH<SAW→2 (but →1 if `XSDI < 0.01·TMPPCT·XBAMAX`); else→5 (→6 if dom min
     DBH<3.0). `XSDI` = SDIBC (before-thin) or SDIAC (after) — the Reineke SDI (see [[fvsjl-keyword-audit]] BSDI).
   - NSTR=2: DBH<SSD→1; DBH<SAW→3; else→6.
   - NSTR=3: DBH<SSD→1; DBH<SAW→4; else→6.
6. **Whole-stand cover** `COVOLP(NTREES)` → `COVER` (the Tot-Cov column).
7. **Report + event-monitor vars** (sstage.f:599-700): print the Structural statistics row; set the
   event-monitor `SSTAGE`/`STRCLASS` variable + `OSTRST` output array. CCADJ (act 444) updates
   CCCOEF via UPDATECCCOEF before this (already verified `.sum`-inert — [[fvsjl-keyword-audit]]).

## Proposed chunks (each validated vs the Structural statistics report)
- **A — helpers:** `covolp` (cover from crown areas + CCCOEF) and `sstghp` (per-stratum PROB-wtd
  DBH/heights/crown/dominant-species). Pure, unit-testable against the report's per-stratum columns.
- **B — stratification:** the HT-sort + two-largest-gap finder → strata boundaries + NSTR. Validate
  the N-Strata column + stratum DBH/cover vs the report.
- **C — classify:** the NSTR×DBH threshold logic → class 1-6. Validate the Struct Class column.
- **D — keyword + wiring + report:** `kw_strclass!` (ksstag.f defaults/overrides), call SSTAGE per
  cycle (after CUTS, like grincr.f:342; emit the before/after Rm rows), the `.out` report writer, the
  event-monitor `SSTAGE` variable, and finally `FVS_StrClass` DBS (only once a fuller binary can
  validate it — the table writer is trivial on top of the class).

## Validation note
Use `fire_early` (FMIN stand, exercises before/after-thin Rm rows) + a managed stand (THIN → the
Rm=1 after-thin row differs). Diff the whole Structural statistics block (strata stats + N-Strata +
Tot-Cov + Struct Class) vs `/tmp/FVSsn_new`. The class is bit-discrete, so it must match exactly;
the cover/DBH are F-format reals (±ULP).
