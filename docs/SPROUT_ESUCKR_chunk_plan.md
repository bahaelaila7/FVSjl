# SPROUT / ESUCKR (stump-sprout regeneration) — chunk plan

The last unported *natural process* in SN: stump/root sprouting after a harvest. Verified scoped
from the Fortran (`bin/FVSsn_buildDir/`). It is multi-component and stochastic, so it is chunked
here rather than attempted in one pass. All other SN management/CUTS/volume extensions are ported
and bit-exact; this is the remaining regen piece (alongside C7/C8 fire/insects/econ).

## How it works (Fortran)

1. **Trigger** (`esnutr.f:108-126`): each cycle, if `LSPRUT` (sprouting enabled) **and** `ITRNRM ≥ 1`
   (trees were cut this cycle), call `ESUCKR`. `LSPRUT` defaults `.TRUE.` (`esinit.f:50`), is set by
   the **SPROUT** keyword (`esin.f:598`), and cleared by **NOSPROUT** / when no valid sprouting
   species is supplied (`esin.f:626/655/681/735`).
   - ⚠ FVSjl currently hard-sets `LSPRUT=false` (esinit.jl) — an "untraced COMMON artifact" that
     makes snt01 bit-exact. Resolve this: the real default is `.TRUE.`, but ESUCKR produces no
     sprouts for snt01's cut species (so the net effect matches). Confirm *why* before enabling.

2. **Cut-tracking — `ESTUMP` (`estump.f`, 115 ln)**: called once per cut tree from `cuts.f:1713`
   (`CALL ESTUMP(species, DBH, max(0, PREM−SSNG), plot, ISHAG)`) and from `fmkill.f:80` (fire kill).
   It appends to parallel arrays indexed by `ITRNRM`: `DSTUMP` (stump DBH), `PRBREM` (removed TPA),
   `ISHOOT`, `JSHAGE` (shade/age). When `ITRNRM > MAXTRE` it merges into the best-matching record.
   → **FVSjl gap:** `cuts!` reduces `tpa` but records nothing per-cut. Need a per-cut log
   (species, DBH, removed TPA, plot) populated in every thin method's removal loop.

3. **Sprout generation — `ESUCKR` (`esuckr.f`, 382 ln)**: for each of the `ITRNRM` cut records:
   - `PREM = PRBREM(i)`, `DSTMP = DSTUMP(i)`, `ISSP = ISP(i)`; skip if `PREM < 0.001`.
   - `NSPREC(VARACD, ISSP, NSPRT, DSTMP)` → `NUMSPR` = sprouts per stump (species + stump-DBH).
   - aspen-only `ASSPTN` (INDXAS) — N/A for SN species; skip.
   - `ESSPRT(VARACD, ISSP, PREM, DSTMP)` (`essprt.f`, present in build) → adjusts `PREM` (per-record
     TPA rules); skip if `PREM < 0.001`.
   - create `NUMSPR` sprout records: `IMC=2`, `ISP=ISSP`, `ITRE=plot`, `PROB=PREM·SMULT`, volumes 0.
   - height `SPRTHT(VARACD, ISSP, SITE, ISHAG, HTI)` → `HT = HTI·HMULT`, plus a clamped
     `BACHLO(0, 0.5, ESRANN)` deviation: `HT += BACHLO·HT/5.5`. DBH from the species H-D inverse
     (`HT2/AA/HT1`, `IABFLG`) when `HT > 4.5`, else a small-tree rule.
   - record-list overflow → `ESCPRS` compression.

4. **Sub-routines:** `NSPREC` (sprout count by species/stump DBH — `DMIN/DMAX/SMULT` per species),
   `SPRTHT` (sprout height by species/site/shade), `ESSPRT` (per-record TPA). Each carries
   **per-species coefficient blobs** → CSVs (`data/southern/sprout_*.csv`).

5. **RNG:** `ESRANN` — the establishment stream (FVSjl already has `:estab`; the existing `establish!`
   uses it). The height deviation `BACHLO(0,0.5,ESRANN)` must consume `:estab` in the exact ESUCKR
   order to stay bit-exact.

6. **Keywords:** **SPROUT** (esin.f opt 26): date + species list + `SMULT` (number mult) + `HMULT`
   (height mult), sets `LSPRUT=.TRUE.`. **NOSPROUT** (opt 27): `LSPRUT=.FALSE.`
   - NOSPROUT is a safe no-op in FVSjl *today* (sprouting already off) — could be recognized now.
   - SPROUT is NOT a no-op (it enables a .sum-affecting feature) — needs the full port.

## Proposed chunks

- **A — cut-tracking + keywords:** ✅ **DONE** (Control `cut_log`/`lsprut`/`sprout_smult`/
  `sprout_hmult`; `cuts!` clears the log each cycle and `_thin_sorted!` appends (species, stump DBH,
  removed TPA, plot) per removal, in removal order, gated on `lsprut`; SPROUT/NOSPROUT parsed in
  `kw_estab!`). Inert until C (the log is write-only), so snt01 + all scenarios stay bit-exact
  (2938/2938). NOTE: only `_thin_sorted!` (THINBTA/ATA/BBA/ABA) logs so far — the other thin methods
  (`_thindbh!`/`_thinprsc!`/`_thin_sdi!`/`_thin_rden!`/`_thin_cc!`/`_thin_qfa!`/`_thin_pt!`) each need
  the same one-line append at their removal point when C lands.
- **B — sub-routine coefficients:** ✅ **DONE** (`src/engine/sprout.jl` + `data/southern/sprout_essprt.csv`).
  `nsprec_sn` / `sprtht_sn` / `essprt_sn` ported pure & bit-faithful to essprt.f's SN blocks. The
  ESSPRT per-species blob is the CSV (per-species `essprt_kind`/`essprt_p1`/`essprt_p2`/`essprt_fsp`
  columns, loaded alongside merch_specs); the 5 forest-special species (64/66/70/75/77, forests
  809/810/905/908) keep their distinct formulas in code with the common-forest ELSE form in the CSV.
  NSPREC/SPRTHT are tiny piecewise rules kept inline (NINT via `nint`). 24 unit tests vs hand-computed
  values (`test/unit/test_sprout.jl`); suite 2938→2962. Still inert (no caller) until C.
- **C — ESUCKR sprout-gen + ESRANN:** the generation loop, creating sprout records with the exact
  `:estab` RNG order. The .sum-affecting chunk. Split into:
  - **C1 — cut-log fidelity:** ✅ **DONE.** Revised the Chunk-A cut-log to true ESTUMP semantics
    (estump.f): a `CutRecord` NamedTuple `(species, dstmp, prem, plot, ishag)`; only the 72 SN
    sprouting species (ISPSPE → `is_sprouting` CSV flag) are logged (ESTUMP returns early otherwise);
    `ishag = IFINT` (cycle length, from `plot.cycle_length` set in `cuts!`). Centralized into a single
    `_log_cut!` helper now called at **every** thin method's removal point (`_thindbh!`/`_thinprsc!`/
    `_thin_sorted!` + the `_remove!`/`_rm!`/`_rmc!`/`_rmp!` closures of SDI/RDEN/CC/PT; `_thin_auto!`→
    sorted and `_thin_qfa!`→dbh/sdi are covered by delegation). 4 unit tests (filter, contents,
    lsprut-off); suite 2962→2966. Still inert (write-only) — no live caller until C2.
  - **C2 — generation loop:** port esuckr.f:156-349 (NSPREC count → ESSPRT survival → SPRTHT height +
    clamped `BACHLO(0,.5,ESRANN)` `:estab` deviation → H-D-inverse DBH → CWCALC crown → tree-record
    init), wire into the cycle hook (esnutr.f:112-124, after COMPRS, gated on LSPRUT && ITRNRM≥1),
    handle the SPROUT keyword per-species/DBH multiplier table (OPGET action 450). The .sum chunk.
    Sub-pieces that are already available to reuse: `crown_width(...)` (CWCALC), `bachlo(...;stream=:estab)`
    + `esrann!` (ESRANN), and the establishment tree-record-init / GRADD-order insertion pattern in
    `establish!`.
    ⚠ **PREREQUISITE found (C2a — sprout-DBH H-D coefficients).** ESUCKR's sprout DBH is the **Wykoff**
    inverse `DBH = HT2/(ln(HT−4.5) − AX) − 1`, `AX = (IABFLG==1 ? HT1 : AA)` (esuckr.f:296-307). These
    `HT1/HT2/AA/IABFLG` are the **CRATET** coefficients (cratet.f:303-345: per-species `AA = mean(ln(H−4.5)
    − HT2/(D+1))` over D≥3 trees; `IABFLG=0` when ≥3 obs & `LHTDRG` & `AA≥0`, else default HT1). FVSjl does
    **NOT** compute these — its height dubbing (`dub_missing_heights!`/`_htdbh_height`) uses the **Curtis-Arney**
    form (htdbh.f:300, P2/P3/P4) and is bit-exact to baseline, so SN dubbing is Curtis-Arney while ESUCKR's
    DBH path is a *separate* Wykoff fit. Before C2 can be bit-exact: extract the SN default `HT1`/`HT2` (COEFFS,
    set in coeffs.f/blkdat) into the species CSV, and port the CRATET `AA`/`IABFLG` per-stand fit (or confirm
    that for a fresh sprout stand the K1<3 / `LHTDRG` path keeps `IABFLG=1` ⇒ default `HT1`, sidestepping the
    fit). Resolve empirically vs live Fortran on a SPROUT stand — do not guess the coefficient source.
- **D — validation:** a SPROUT + harvest stand (cut a sprouting species, e.g. an oak), 3-way vs live
  Fortran; resolve the LSPRUT default. Bit-exact bar `:estab` Float32 noise.

## Validation note

Needs a stand that both **cuts a sprouting species** and has **SPROUT** active. snt01 never triggers
ESUCKR (LSPRUT-effectively-false), so a new scenario is required. The sprout heights/DBH are
RNG-driven (`:estab`), so expect the same ±Float32 establishment noise seen in the existing bare-stand
PLANT/NATURAL regen tests.
