# COMPRESS (tree-record compression) — chunk plan

`COMPRESS` (initre.f:8000, option 78) schedules **act 250**: reduce the tree list from
`ITRN` records to `NCLAS` classes by clustering similar trees and merging each class into
one representative record. comprs.f is **1010 lines** — an IBM-SSP-eigen PCA clustering —
so it is chunked here rather than ported in one pass. Bit-exactness is *fragile* (it runs
through an iterative eigensolver + several descending sorts with FVS tie-breaks), so each
sub-step must be validated against live Fortran, not just "records reduced."

## Status — ✅ DONE (2026-06-24, with the eigensolver substituted)
- ✅ **Keyword** (`kw_compress!`): recognized + scheduled (icflag 250, params = target
  records `NCLAS`, `PN1`%, date). Defaults: `NCLAS = MAXTRE/2 = 1500`, `PN1 = 50`, date 1.
- ✅ **Algorithm** (`src/engine/compress.jl`, `compress!` + `apply_compress!`): the full
  comprs.f / comcup.f port — standardize → correlation matrix → eigen → PC scores → Method 1
  (gap breaks) → Method 2 (range splits) → PROB-weighted merge — wired at act 250 in
  `grow_cycle!` (suppresses tripling = NOTRIP). **The 1966 IBM-SSP Jacobi eigensolver is replaced
  by `LinearAlgebra.eigen`** (project direction: use Julia's linear algebra rather than re-port
  the routine). Consequences, validated vs live Fortran (`test_compress.jl`):
  - exact: reduces to **NCLAS** records and **conserves total TPA exactly**;
  - at the compression cycle the compressed stand's `.sum` aggregates (TPA/BA/SDI/CCF/TopHt/QMD/
    cubic volume) are **bit-identical** to Fortran — the merge math is correct;
  - the **multi-cycle trajectory diverges** (≈ several % by late cycles) because LAPACK and the
    SSP routine produce slightly different eigenvectors → slightly different PC scores → a
    different class PARTITION among near-identical records (compounded by sort tie-breaks that
    `RDPSRT` resolves differently than a stable sort). This is the accepted cost of the
    eigensolver substitution; bit-exact COMPRESS would require porting `EIGEN` + `RDPSRT` exactly.

### Original chunk plan (kept for reference / if bit-exactness is ever required)

## How it works (comprs.f)
Classification variables (per tree): `HT`, `ICR`, `IMC`, `ln(DBH)`, `DG` (5 vars, NRANK=5).

1. **Standardize** (comprs.f:160-205): `MEANSD` each var (with the floor stddevs 1 / 1e-4 /
   1e-4 / 5e-3 / 0.02); center+scale; accumulate the 5×5 cross-product `XTX` (lower-tri
   vector storage `N(N+1)/2`), subtract `XSUM(i)·XSUM(j)/n`, divide by `n−1` → correlation
   matrix (diagonal forced to 1).
2. **EIGEN** (comprs.f:251): IBM SSP 1966 symmetric eigensolver (Jacobi rotations) →
   eigenvalues on the `XTX` diagonal, eigenvectors in `EIVECT`, sorted. **Port this exactly**
   (the rotation order + convergence threshold drive bit-exactness). Then sign-fix EIVECT(4)
   <0 ⇒ flip PC1, EIVECT(7)>0 ⇒ flip PC2; scale `EIVECT(jk) /= STDDEV(k)`.
3. **Scores** (comprs.f:296-318): `WK3 = 25·(PI·ISP + ITRE)` (species+point base) `+` PC1
   projection `+4`; `WK4 =` PC2 projection. `RDPSRT(WK3, IND, descending)`.
4. **Method 1** (comprs.f:330-365): differences of sorted `WK3`; `NCLS1 = round(NCLAS·PN1)`
   (clamped by the count of non-tied gaps `ISIG`); the `NCLS1−1` largest gaps (via `RDPSRT`
   then `IQRSRT` ascending) are the class breaks; `IND1(NCLS1)=ITRN`.
5. **Method 2** (comprs.f:375-470+): per class, range on PC1 (`WK6`) and PC2 (`CMRANG`→`WK5`);
   `NCLS2 = NCLAS−NCLS1` times, split the class of largest range at its mid-range record
   (re-sorting within-class on PC2 if the PC2 range dominates); stop if largest range
   ≤ `RNGMIN`.
6. **Merge** (comprs.f ~550-800): collapse each class to one record — sum `PROB`+dead, the
   `Nov-1995` averaging uses (PROB + trees-dying) weights for the attribute means; handle
   truncated-tree + subplot classes separately; reset `ITRN=NCLAS`, invalidate the IND/
   species sorts (caller re-sorts).

## Proposed chunks
- **K — keyword:** ✅ DONE (above).
- **A — helpers:** `MEANSD` (mean/var/stddev), `CMRANG` (within-class range), `IQRSRT`
  (integer quicksort ascending). `RDPSRT` already exists (reuse). Pure, unit-testable.
- **B — EIGEN:** port the IBM-SSP symmetric Jacobi eigensolver bit-faithfully; unit-test
  eigenvalues/vectors of a known 5×5 vs the Fortran (a DEBUG dump of `XTX`/`EIVECT`).
- **C — classify:** standardize → corr matrix → scores → Method 1 breaks → Method 2 splits.
  Validate the partition (class membership) vs a Fortran DEBUG run.
- **D — merge + wire:** the class→record averaging; apply at act 250 in `grow_cycle!`
  (COMCUP/GRINCR end); validate the `.sum` bit-exact vs live Fortran on a >NCLAS-record
  stand (e.g. a tripled or large-inventory stand with `COMPRESS … 100`).

## Validation note
Needs a stand with ITRN > NCLAS so compression actually fires (a large inventory, or a
heavily-tripled stand). Expect the merge to conserve total TPA exactly but change BA/volume
(nonlinear in DBH) per the exact partition — so the eigen + sort tie-breaks must match.
