# Bugs found in the live Fortran FVS source (the oracle) — jl is CORRECT; do NOT replicate

The live Fortran binaries (`/workspace/ForestVegetationSimulator/bin/FVS{sn,ne,cs}_buildDir`)
are the campaign's oracle. They are normally ground truth, but the FVSjl differential work has
surfaced places where **the real Fortran itself is wrong** — it emits incorrect or physically-
impossible output. Per the campaign doctrine ("the oracle can be WRONG; don't blindly match")
and the user's explicit decision (2026-07-03), **FVSjl stays CORRECT on these paths and does NOT
reproduce the bug.** They are recorded here (separate from the Oracle-A transliteration bugs in
[ORACLE_A_BUGS.md](ORACLE_A_BUGS.md)) so that:
- FVSjl↔live `.sum`/report diffs on these columns are understood as live-side, not jl regressions;
- they can be reported upstream to the FVS maintainers;
- a future "strict byte-fidelity" mode (if ever wanted) knows exactly what to gate.

Ground truth = the FVS source under `bin/FVS{sn,ne,cs}_buildDir/*.f` (cited per entry).

Format per entry:
> ### <title>
> - **Live Fortran:** file.f:line — what it does (wrong)
> - **FVSjl:** correct? (and where)
> - **Evidence it's a bug** (internal inconsistency / source trace)
> - **Impact / how found**
> - **Disposition**

---

### D36 — econ HRVRVN leaves removed-sawlog-cubic (`OSCREM(7)`) stale into a no-cut cycle
- **Live Fortran:** the `.sum` "removed cubic saw volume" column (col 16) = `IOSUM(22,IKNT)` =
  `INT(OSCREM(7)/GROSPC+0.5)` (`disply.f:342`). `ONTREM(7)` (removed trees) is explicitly re-zeroed
  every cycle (`cratet.f:158`, `fvs.f:433`); the full `O*REM(I)` removal arrays are zeroed only inside
  `cuts.f:324-329` **during an actual cut**. There is **no matching per-cycle reset of `OSCREM(7)`**
  (grep-confirmed). With the ECON `HRVRVN` keyword active, the econ path re-populates `OSCREM(7)` after
  the cut and leaves it set, so the NEXT (no-cut) cycle's summary row reports the prior cut's value.
- **FVSjl:** CORRECT — it zeroes all removal columns on a no-cut cycle (`remScuft = 0`).
- **Evidence it's a bug:** the live 2005 row for `econ_strtecon` is `rem: 0 0 0 23 0` — i.e. **23 cuft
  of sawlog "removed" with 0 removed-trees AND 0 removed-total-cubic**, which is physically impossible
  (you cannot remove sawlog cubic volume while removing zero trees and zero total cubic). The thin year
  (2000) is bit-exact both engines (`rem: 20 106 97 23 106`); only the subsequent no-cut cycle diverges.
- **Impact / how found:** found 2026-07-03 by broadening `divergence_sweep.jl` to the REMOVALS columns
  (13–17), which the state-only sweep never checked. Scope: TRIGGERED BY `HRVRVN` — `cut_thinsdi`
  (identical THINSDI, no econ) correctly zeroes `remScuft` at 2005; only econ-harvest stands
  (`econ_strtecon`, `econ_u5`) exhibit it. Affects only the `.sum` removal column; the ECON revenue
  TABLES themselves are bit-exact (validated separately).
- **Disposition:** 📌 jl stays correct; NOT replicated (user decision 2026-07-03). Documented here + as
  D36 in [DIVERGENCE_FIX_CAMPAIGN.md](DIVERGENCE_FIX_CAMPAIGN.md). Upstream-reportable to FVS maintainers.

---

### D37 — htgf.f height-cap on TRIPLED records reads UNINITIALIZED array memory (`HT(ITFN)` before TRIPLE/SVTRIP sets it)
- **Live Fortran:** `htgf.f:297` caps a tripled satellite record's height growth via
  `IF((HT(ITFN)+HTG(ITFN)).GT.SIZCAP(ISPC,4)) HTG(ITFN)=SIZCAP(ISPC,4)-HT(ITFN)`, where
  `ITFN=ITRN+2*I−1` (`htgf.f:292`) is a **not-yet-created** record slot (beyond the live count `ITRN`).
  `HTGF` is `CALL`'d at `grincr.f:443`; `TRIPLE` — whose `SVTRIP` is what actually sets `HT(ITFN)=HT(I)`
  for the new satellite records — is `CALL`'d LATER, at `grincr.f:543`. So at cap time `HT(ITFN)` holds
  **stale/uninitialized array memory** (leftover from a prior compacted record), and the tripled record
  escapes the height cap by a memory-dependent amount.
- **FVSjl:** CORRECT / deterministic — jl caps each satellite record faithfully against the parent height
  it inherits (`copy_tree!`), so its capped tall trees sit ~3 ft lower than the live memory-dependent ones.
- **Evidence it's a bug:** the live TopHt spread among the tripled top trees is **72.0 AND 73.7** (only
  ~1.7 ft apart) — NEITHER a clean `HT(ITFN)=0` full-escape NOR a full cap. A deterministic algorithm on
  the same records would give a single uniform capped height; the non-uniform escape is the fingerprint of
  reading uninitialized memory. `NOTRIPLE` runs are BIT-EXACT vs live (the divergence appears ONLY on the
  tripling path), confirming it is not a growth/mortality logic difference.
- **Impact / how found:** found 2026-07-06 (tolerance-closure campaign, `test_treeszcp.jl:100`). Scope: the
  TREESZCP size-cap × tripling path; TopHt drift ≤3–4 ft. It is **not deterministically reproducible** —
  matching FVS here would mean emulating an uninitialized-array read (undefined behavior). Left as a
  documented `@test_broken` (accepted-irreducible, more fundamentally irreducible than the COMPRESS
  eigensolver, since it is not even a well-defined value).
- **Disposition:** 📌 jl stays correct; NOT replicated. Accepted as a terminal genuinely-irreducible item
  for `TOLERANCE_COMPLETE` (user decision 2026-07-06). **Upstream-reportable to FVS maintainers — confirm
  the `HT(ITFN)` read-before-`SVTRIP`-write with the FVS developer** (likely a real UB bug: the cap should
  run after TRIPLE sets `HT(ITFN)=HT(I)`, or use `HT(I)` directly).

---

### D38 — SN volume init (R9 Clark `r9ht`) SIGFPE-crashes on real FIA stands with short trees (≤17.3 ft: invalid-op; just >17.3 ft: FE_UNDERFLOW — `r9ht` is missing the underflow guard `r9cuft` already has at `r9clark.f:1015`)
- **Live Fortran — EXACT LINE:** `bin/FVSsn_buildDir/r9clark.f:1286`, in subroutine `r9ht`:
  ```fortran
        Y = (1.0 - 17.3/totHt)**p
  ```
  When a tree's total height `totHt < 17.3 ft`, `17.3/totHt > 1`, so the base `(1.0 - 17.3/totHt)` is
  **negative**, and raising a negative real base to the **real (fractional) exponent `p`** yields
  **NaN → SIGFPE (invalid operation)**. The SN binary is built with FP-exception trapping, so it
  **aborts (exit 136)** instead of propagating NaN. Call path (symbolized `-g` backtrace):
  `r9ht (r9clark.f:1286) ← r9clark (r9clark.f:253) ← volinit (volinit.f:414, the `CALL R9CLARK`) ←
  VOLINITNVB (volinit.f:838) ← UPDATE (update.f:108) ← TREGRO (tregro.f:52)`. The whole run aborts with
  **no `.sum` produced**. (Sibling lines 1283/1285 `(1.0-4.5/totHt)**{r,p}` have the same hazard for
  `totHt < 4.5`; 1284 `/(1-G)` and 1287 `/(X-Y)` are related singularities.)
- **Trigger:** any tree with **total height < 17.3 ft** reaching the R9 Clark profile — short
  understory/sapling trees, ubiquitous in dense real FIA stands, and dubbed heights of the many
  `HT`-null records. The R9 Clark profile implicitly assumes `totHt ≥ 17.3` (17.3 ft = the merch/upper
  reference height) and never guards the short-tree case. Confirmed crashing stand `238869289010854`:
  **35 trees, 1 with measured `HT`<17.3 ft and 7 with NULL/zero `HT`** (dense: SDI≈350, RELDEN≈260).
- **FVSjl:** **CORRECT / strictly more robust** — jl's ported volume path (`r9clark_vol.jl`) does NOT
  crash; it produces a valid `.sum`. Measured: jl produced a valid `.sum` on **119/120** treed SN FIA
  stands, **including all 34** that live FVS SIGFPE-killed.
- **Evidence it's a bug:** (1) gfortran backtrace pins the fault at `volinit.f:414` (the R9CLARK call);
  (2) `exit=136` (128+SIGFPE) is uniform — **34/34** no-sum treed SN stands in a stratified sample
  exit 136; (3) reproducible on **both** the extracted subset DB **and** the full 66 GB
  `SQLite_FIADB_ENTIRE.db` (so it is real live-FVS-on-real-data, not an input/harness artifact);
  (4) crashing on valid FIA inventory is by definition wrong — a growth simulator must not abort on a
  stand containing zero-height tree records.
- **Impact / how found:** surfaced by the **stratified FIA-plot mass validation** (task #94,
  `test/harness/fia/sweep/`). Live FVS SN fails to emit a `.sum` on **~28–42%** of *treed* SN FIA COND
  stands (seed-99 sample: 34/120 = 28%). Because these have no live oracle output, they are necessarily
  excluded from the FVSjl↔live `.sum` bit-exact denominator (the 98.4% cycle-0 rate is on both-produced
  stands only); they are NOT jl failures — jl runs them fine.
- **Cross-variant:** the R9 Clark volume library (`r9clark.f`/`r9ht`) is shared, so **all four variants
  crash** on short-tree FIA stands — measured live-FVS SIGFPE rate on 100 treed FIA COND stands/variant:
  **SN 30% · NE 2% · CS 5% · LS 3%** (SN by far the worst; live-OK rates SN 70 / NE 98 / CS 95 / LS 97).
  jl runs all of them.
- **ROOT CAUSE SHARPENED + VERIFIED FIX (2026-07-08):** there are TWO distinct trap conditions in
  `r9ht`, and the second is the more common one:
  1. **Invalid-op** — `totHt < 17.3` ⇒ negative base ⇒ `(neg)**p = NaN` (the line-1286 case above).
  2. **FE_UNDERFLOW** — `totHt` *just above* 17.3 ⇒ base `(1-17.3/totHt)` is tiny and `p` large, so
     `Y = base**p` **underflows to a denormal**, which SIGFPEs under the SN binary's underflow trapping.
     Measured on real crasher `1152014964290487` via a `-g` value stamp: `totHt=17.40, base17=0.00575,
     p=17.81 ⇒ Y≈1e-40` (Float32 denormal) → trap at `r9ht`'s `Y=` line (NOT the `X-Y` divide).
  - **THE FIX IS FVS'S OWN EXISTING GUARD.** FVS *already* guards this identical underflow in the
    **cubic-volume** routine `r9cuft` at **`r9clark.f:1015`**: `IF((1.0-17.3/totht).LT.0.005748.AND.
    p.GT.14)THEN` (and again at `:1082`). The **height** routine `r9ht` is simply **missing the same
    guard**. Applying the verbatim guard to `r9ht`'s `Y=` line — plus a `totHt<=17.3` early return for
    the invalid-op case, and denominator guards on `(1-G)`/`(X-Y)` for very tall stems — is the fix.
  - **Empirically verified** (debug relink `/tmp/FVSsn_fixtest`, since restored): on 300 SN FIA stands
    the fix **cleared 18/18 SIGFPE crashers** (all now emit `.sum`) and left **276/282 non-crashers
    BIT-IDENTICAL**. The **6 changed** non-crashers are cases where `r9cuft` *already* zeroed `Y` (via
    `:1015`) but `r9ht` computed the tiny denormal — i.e. the fix makes `r9ht` **consistent with the
    already-shipping `r9cuft` guard**, arguably a correctness improvement, not a regression.
  - Full fixed source saved at `docs/patches/r9clark_D38_underflow_fix.f` (NOT applied to the oracle —
    the live binaries stay pristine as ground truth per campaign doctrine).
- **Disposition:** live-FVS UB — **jl stays correct (does not crash)**. Upstream fix is the one-guard
  port above (mirror `r9clark.f:1015` into `r9ht`). Recommended alternative for maintainers: since
  gradual underflow to a denormal is well-defined IEEE behavior, simply **not trapping FE_UNDERFLOW**
  clears the underflow case with zero output change (the invalid-op `totHt<17.3` case still needs a
  base guard). Fault lines confirmed via `-g` recompiles: `r9ht` `Y=` (underflow, `totHt=17.4`) and
  `r9clark.f:1286` (invalid, `totHt<17.3`).

## Shared SDI overflow on degenerate ultra-dense micro-stands (SN; FVSjl reproduces)
CN 218434248010854 (SN, 2 tree records, fixed 1/300ac plot ⇒ TPA~7000, AGE 28): at cycle 2026 the reported SDI
jumps to ~4.38 MILLION (physically impossible; SDImax<1500) and TCuFt momentarily reports 0, in BOTH live FVSsn
AND FVSjl. This is a live-FVS numerical pathology on a degenerate micro-stand, not an FVSjl divergence — FVSjl
tracks the (absurd) live SDI within 0.46% (4381035 vs 4361021). Recorded as a FIDELITY success (jl reproduces
even FVS's degenerate behavior). The residual TCuFt/MCuFt wobble (412/350 cuft) is a small-magnitude consequence,
cornered by the vol_max_abs≥300 escalation floor (audit slice 43n).

## Live FVS SIGFPE crash on high-expansion 0.1" seedling records (FVS-UB; FVSjl survives)
Live FVSsn (relinked FS2026.1) CRASHES with a floating-point exception (SIGFPE, exit 136) on stands containing a
tree record representing >1000 TPA at DBH 0.1" (e.g. CN 1224249623290487: sp611 dbh 0.1" TPA 1199). FVS prints its
own guard `FVS40 WARNING: TREE RECORD REPRESENTING GREATER THAN 1000 TPA ENCOUNTERED. MAY CAUSE MATHEMATICAL
ERRORS` and then dies (a div-by-zero in the per-tree small-tree stats). These stands are common in some ecoregions
(the FIA microplot expansion produces them). FVSjl PROJECTS them fine (full .sum). ⇒ FVS-UB, not an FVSjl
divergence — the oracle is un-runnable, so there is nothing to compare against, but FVSjl is strictly MORE robust.
Recorded by the sweep as dig_class `live_crash` (ledger_fia run_live detects termsignal/exit>128), so coverage is
honestly accounted: {comparable = bit_exact+ulp_class} + {live_crash} + {skip}, never silently dropped.
