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

> ★★ **RESOLVED — FIX APPLIED TO LIVE FVS, ALL 4 VARIANTS (2026-07-11).** A crash is a *bug*,
> not a tolerable divergence; the faithful action is to fix FVS at source (for maintainer
> submission) so the stand yields comparable output — never to replicate the crash in FVSjl.
> Sites pinned by `-g` backtraces on real crashers at the **exact original flags** (`-O0
> -ffpe-trap=invalid,zero,underflow,overflow,denormal`; note **`-O2` masks the trap**):
> `r9ht:1286`, `r9dib:~1213`, `r9cuft:1018`, `r9cuft`-V2:1086 — all the Clark taper terms
> `(1-h/totHt)**p`. **Measurement** (e.g. crasher `1224249126290487`: `totHt=17.49, U2=17.30,
> p=19.17, base=0.0108` → `base**p·(totHt-U2)` **subnormal**) showed the dominant failure is a
> **legitimate gradual underflow to a denormal**, not a logic error — and the `YW-2023` fixed
> threshold `0.005748.AND.p.GT.14` can't cover it (it's p-dependent). Since denormal underflow
> is well-defined IEEE and **FVSjl/Julia does not trap it**, the clean FVSjl-exact fix is
> **two-part**:
> 1. **BUILD FLAG** (`bin/makefile_Xbuild`): drop `underflow,denormal` from `-ffpe-trap`
>    (→ `invalid,zero,overflow`). Denormals flow exactly as the port computes them; **fixes
>    every denormal crash site at once with ZERO output change**.
> 2. **SOURCE** (`volume/NVEL/r9clark.f`): guard the residual **invalid-op** (negative base,
>    `totHt<17.3`) in `r9ht`+`r9dib` — `Y=0` is the correct limit (17.3′ term vanishes),
>    mirroring the guard `r9cuft` already has (~:1015). Matches FVSjl (`r9clark_vol.jl` r9ht:236).
> **Validation:** 8/8 first-wave + **162/162 second-wave** SN crashers now emit `.sum`; **cycle-0
> volume bit-exact vs FVSjl**; **40/40 normal stands byte-identical** patched-vs-pristine (data
> rows). **All 4 oracles relinked** (`/tmp/FVS{sn,ne,cs,ls}_new`): SN+CS get the build-flag +
> source guards; **NE/LS get the build-flag only** — their build dirs were compiled with a
> different gfortran (`.mod` ABI mismatch blocks recompiling r9clark.o), but the trap-drop is
> `main.o`-only and clears their (denormal) crashers, verified. Any residual `totHt<17.3`
> invalid-op crash in NE/LS would need the source guard via a full rebuild (tiny: NE 5 / LS 1
> live_crash total). Clean maintainer patch (build flag + source): **`docs/patches/
> r9clark_D38_allsites.patch`**. Previously-crashed SN stands re-run through the patched oracle
> reclassify into the ordinary cornered taxonomy (structure_densephase / print_boundary /
> count_straddle), 0 `live_crash`, 0 UNCLASSIFIED.

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

## `live_crash` sweep category = the D38 `r9ht` SIGFPE bug (already root-caused + fixed + validated)
The full-population sweep's `live_crash` dig_class (ledger_fia run_live detects termsignal/exit>128; e.g. CN
1224249623290487, sp611 dbh 0.1") is the SAME live-FVS SIGFPE documented and RESOLVED above as **D38** — the R9
Clark `r9ht` short-tree underflow/invalid-op crash (r9clark.f:1286 / the missing `r9cuft` guard at :1015). It is
NOT a new bug and NOT "plausible-but-unvalidated": jl carries the D38 fix, and the fix was **empirically validated
against a PATCHED live binary** (`/tmp/FVSsn_fixtest`) — 18/18 crashers cleared (all emit .sum) and 276/282
non-crashers BIT-IDENTICAL. So on a `live_crash` stand jl produces the CORRECT projection (the buggy SHIPPING
oracle just can't run it to confirm; the FIXED oracle does, and jl matches). The sweep records `live_crash` only
for honest COVERAGE ACCOUNTING — {comparable = bit_exact+ulp_class} + {live_crash = D38, jl-correct} + {skip} —
so the crash stands are visible, not silently skipped. The D38 measurement (SN ~30% of treed stands crash live)
explains the region-variable comparable rate. (Meta: this "new" finding was a RE-DISCOVERY of D38 — always grep
FVS_SOURCE_BUGS.md before writing up an FVS crash.)

## D38 addendum — the guard EXISTS in FVS's own tree but is ORPHANED (NVEL vs fvsMod)
Root of the D38 `r9ht` gap: there are TWO copies of the Clark profile in the FVS source tree:
- `volume/NVEL/r9clark.f` — the National Volume Estimator Library (git submodule) copy. **This is what EVERY
  variant compiles** (confirmed in all `bin/FVS*_CmakeDir/FVS*_sourceList.txt`: `../../volume/NVEL/r9clark.f`).
  Its `r9ht` (line 1286) is UNGUARDED — the D38 crash. (It DOES guard `r9cuft`/`r9dib` at line 1015.)
- `volume/r9clark_fvsMod.f` — an older FVS-internal copy (last commit 2025-01-08) that **DOES** carry the
  underflow guard in `r9ht` (at lines 1015/1213/1292). But it is in **NO** variant's source list — orphaned,
  not compiled by anything.
So FVS ships the unguarded NVEL copy while a guarded copy sits unused in the same repo. The D38 fix is therefore
literally FVS's OWN code (the `fvsMod` guard = the `r9cuft` guard). UPSTREAM REPORT: NVEL's `r9ht` should carry
the same guard its sibling `r9cuft` already has (or FVS should compile `r9clark_fvsMod.f`). SN uses this via NVEL
because the Clark profile is the shared national eastern taper model (SN=R8/NE/CS/LS all route to it through
NVEL) — the "R9" filename is the library's origin, not a variant misroute.

## D38 CORRECTION (multi-site) — the crash is NOT single-site, and the fix/patched-oracle are INCOMPLETE
Earlier D38 text (and audit slices 43t/43q) implied the r9ht Y-guard (or my 5-guard patch) fully fixes the crash
and that FVSjl's live_crash stands validate bit-exact against a "D38-patched oracle". THAT WAS OVERCLAIMED.
Measured (backtrace on 40 real live_crash SN stands, -g relink):
- The r9ht Y-guard alone (= FVS's own orphaned fvsMod fix, isolated to a19c41b4, 16 lines) clears 32/40.
- My 5-guard patch (docs/patches/r9clark_D38_underflow_fix.f) clears the SAME 32/40 — NOT more.
- The remaining 8/40 crash at a DIFFERENT site: `r9cuft` (cubic volume), `r9clark.f:1086` — the V2/V3 log-segment
  volume computation `(1-U2/totHt)**p`/`(1-L2/totHt)**p` + `/(totHt-17.3)`,`/(totHt-17.3)**2` divisions. NEITHER
  the fvsMod guards NOR my patch touch this. Trigger: degenerate multi-species 0.1"-seedling stands, NULL heights.
⇒ D38 is a MULTI-SITE hazard (≥ r9ht:1286 AND r9cuft:1086). NO existing guard set (FVS's own or mine) clears all
crashers. Therefore `/tmp/FVSsn_patched` is an INCOMPLETE oracle — it can validate the r9ht-crash subset but NOT
the r9cuft-crash subset. The honest, unaffected ground remains: these are FVS-UB stands the SHIPPING oracle
crashes on (dig_class `live_crash`, honest coverage accounting); FVSjl projects them (its ported r9clark_vol.jl
does not crash), but the "validated bit-exact vs a patched oracle" claim holds ONLY where a COMPLETE guard set
exists — which it does not yet. Root cause of the whole gap unchanged (submodule migration dropped FVS's local
r9clark guards; see D38 addendum). Do NOT claim the crash stands are oracle-validated until r9cuft is also guarded.

---
## CS essprt.f:216-217 — stump-sprout reciprocal SIGFPE under THINBBA (real live-FVS crash, fix proposed, env-blocked)

**Found 2026-07-11 via the Pillar-3 CS management differential** (audit slice 43bw): 3 of 40 sampled CS stands
SIGFPE-crash live FVScs under **THINBBA** (salvage/plant/none do NOT crash). Backtrace (all 3, identical path):
```
essprt.f:217  (SIGFPE, signal 8)
esuckr.f:239 → esnutr.f:119 → gradd.f:229 → tregro.f:52 → fvs.f:376 → main.f
```
**Root cause.** `SUBROUTINE ESSPRT(VAR,ISPC,PREM,DSTMP)` computes the stump-sprout "probability-remaining" multiplier.
The CS branch, CASE(57,58) (essprt.f:215-217):
```fortran
CASE(57,58)
  PREM = PREM * (1. / (1. + EXP(-(-2.8058 + 22.6839 * (1./((DSTMP/0.7788)-0.4403))))))
```
The term `1./((DSTMP/0.7788)-0.4403)` divides by zero when the cut-stump diameter `DSTMP == 0.7788*0.4403 ≈
0.34296"`. Under the shipped build's FPE trap (`-ffpe-trap=invalid,zero,overflow`, confirmed from the .o DWARF
producer string), the `1./0.` divide-by-zero raises SIGFPE. (The DENOM<0 side is also hazardous: `1./(tiny neg)`
→ −∞ → predictor +∞ → `EXP(+∞)` → the `overflow` trap.) THINBBA cuts a species-57/58 hardwood whose stump lands
on/near that diameter ⇒ post-harvest sprouting (esuckr→esnutr→gradd→essprt) hits the singularity ⇒ crash.

**Anomaly.** CASE(57,58) is the ONLY CS branch using the *reciprocal* `1./((DSTMP/0.7788)-0.4403)`. Every other
branch built on the same base term (CASE 47, 54: `((DSTMP/0.7788)-0.4403)*2.54`) uses it as a *linear* predictor
(×2.54, in→cm). The reciprocal form (with the odd 22.6839 coef and no 2.54) looks like a transcription anomaly, but
regardless it is a genuine numerical singularity on valid inputs.

**Both-sides.** FVSjl does NOT crash on these 3 stands (dig_class is `live_crash`, not a jl error ⇒ FVSjl produced
output) — its essprt port evidently uses a safe formulation / guarded limit. So this is a live-FVS-only hazard.

**Proposed maintainer fix (essprt.f, CASE(57,58)) — source guard clamping the LOGISTIC ARGUMENT** (an
epsilon-on-DENOM-only guard is INSUFFICIENT: the DENOM<0 side still overflows EXP):
```fortran
CASE(57,58)
  DENOM = (DSTMP/0.7788) - 0.4403
  IF (ABS(DENOM) .LT. 1.0E-6) THEN
    XARG = SIGN(1.0E30, DENOM)            ! DENOM→0 : predictor → ±∞ limit
  ELSE
    XARG = -2.8058 + 22.6839 * (1./DENOM)
  ENDIF
  XARG = MAX(-80.0, MIN(80.0, XARG))      ! keep EXP finite (no overflow/underflow trap)
  PREM = PREM * (1. / (1. + EXP(-XARG)))  ! → 1 (DENOM>0 limit) or 0 (DENOM<0 limit): the correct limits
```
This preserves the model's value for every non-singular DSTMP and returns the mathematically-correct logistic
limit (1 / 0) at the singularity, with no FP exception. ALTERNATIVE (R9-clark-style build-flag fix): drop `zero`
and `overflow` from main.o's `-ffpe-trap`, letting IEEE ±Inf propagate to the same limits — but that masks
divide-by-zero globally, so the source guard is preferred.

**★ CORRECTION — the BUILD-FLAG fix IS validatable + faithful here (earlier "env-blocked" was WRONG).** The build-dir
`main.o` was compiled with the LOCAL gfortran (Debian **12.2.0**, per its DWARF producer) and does no floating-point
(40-line driver); recompiling it with the FULL exact flags gives BYTE-IDENTICAL `.text/.data/.rodata` (only `-g`
debug sections differ). My earlier "relink non-reproducible" was a two-fold artifact: (1) I dropped the
`-fintrinsic-modules-path`/`-fpre-include` flags (partial-flags → different .text), and (2) I compared the FULL
`.sum` text including the per-run **timestamp header** (a spurious diff). With full flags + timestamp-stripped DATA
comparison: recompiling `main.f` with `-ffpe-trap=invalid` (dropping `zero,overflow`) → relink → **50/50 normal CS
stands DATA-bit-identical** to `/tmp/FVScs_new` AND all 3 crashers FIXED (SIGFPE→valid output, the `1./0.→+Inf`
propagating to the correct logistic limit — same value FVSjl computes). So the crash IS fixable+validatable here.
(NOTE: the *source-recompile* path — recompiling essprt.f itself — remains blocked, since essprt.o was built with
SUSE gfortran **15.2.1** and a 12.2.0 recompile perturbs numerics; only the trap-flag-on-main.o path is faithful.)

**Fix recommendation (unchanged) + deployment decision.** The precise MAINTAINER fix is still the essprt.f source
guard (logistic-arg clamp — fixes ONLY this site). The build-flag fix (`main.o -ffpe-trap=invalid`) is a VALIDATED
alternative here but BROADER: it stops trapping divide-by-zero/overflow GLOBALLY (the R9-clark fix deliberately kept
`zero,overflow`), so it could mask a genuine div/0 elsewhere. Oracle `/tmp/FVScs_new` LEFT PRISTINE and NOT
hot-patched — because the SWEEP runs regime=none, which does NOT hit this crash (it needs THINBBA→stump-sprouting →
essprt CASE 57/58), so the sweep's CS coverage is unaffected; the crash is a Pillar-3 THINBBA-management issue only.
The directive ("fix live FVS for maintainer submission; never tolerate as live_crash") is satisfied: root-caused +
validated working fix + source-guard proposal for submission — not tolerated. Affected sample stands (thinbba):
1910906629290487, 488847180126144, 224864192010661 (+ likely more; the 3-stand sample is not exhaustive). All test
binaries removed; build-dir .o untouched (compiled only to scratch).

## NE summary volume dropped to 0 at extreme-height cycles (FVSjl reproduces the CORRECT nonzero volume)
Found via the FIA-compat sweep (audit 43do), NE stand 207147469020004 (dense regen, 900 TPA seedlings → sawtimber).
The NE .sum reports **0 for ALL volume columns (TCuFt/MCuFt/SCuFt/BdFt) at 2053 and 2063** while the stand carries
135→92 TPA at QMD 16→20 (structure — TPA/BA/SDI/CCF/TopHt/QMD — is bit-exact vs FVSjl every cycle). Volume is
bit-exact 2013-2043, then drops to 0 from 2053 on.
ROOT (proven by instrumenting r9clark.f on a scratch NE build; source restored + oracle pristine): live's
r9clark COMPUTES CORRECT NONZERO per-tree cubic volume for these trees — RCTRACE at the r9cuft gate shows
cfVol = 172-257 cuft/tree with errFlg=0 at all four internal gates (r9Prep/r9dia417/r9totHt/r9cuft). So the
volume is LOST DOWNSTREAM of r9clark (the vollib09 driver or the summary accumulation), NOT in the Clark taper
model. The trigger correlates with extreme tree heights: this stand's NC-128 NE height model yields TopHt 258 ft
(2053) / 295 ft (2063) — 22" trees at ~293 ft, even 0.7" seedlings at ~351 ft (FVSjl reproduces these heights
BIT-EXACTLY, so the height behaviour itself is faithful, not the bug). FVSjl sums the per-tree volumes correctly
(15284 TCuFt @2053) — i.e. FVSjl reports what the .sum SHOULD show; it does NOT replicate this FVS summary bug.
STATUS: cornered as FVS-bug in the FIA-compat campaign (FVSjl is the correct side; no FVSjl change). Exact FVS
loss-location (vollib09.f vs sumout/summary array) not yet pinned — needs a second instrumentation pass on the
driver/summary; low priority since FVSjl is already correct. Repro: run NE stand 207147469020004 to cycle 2053.

## Three SIGFPE crash classes on real FIA inventory — FIXED with minimal source guards (2026-07-18)
Full write-up + patches: `docs/FVS_LIVECRASH_AUDIT.md`, `docs/patches/livecrash_*.patch`. The FIA sweep recorded
60 `live_crash` stands (SIGFPE on the shipped trapping build); the 12 that still reproduced collapse to 3
unguarded-arithmetic root classes (backtraces named the *caller* frame — several `dgdriv:134/:353`/`grincr:437`
are `CALL DGF`, and `grincr:449`→`regent:324` is `CALL HTDBH`). FVSjl runs all clean (correct side). Each trigger
was confirmed by direct instrumentation (not just reproduction, which is FP-precision-sensitive under the
in-container gfortran 12.2.0 ≠ official 15.2.1). Fixes are no-ops off the degenerate path (51/51 normal stands
byte-identical patched-vs-official across all 4 variants) and match FVSjl on the degenerate one.

1. **DGF `ALOG(0)`** — `{ls,cs}/dgf.f` (LS:456 / CS:550): a large tree with tiny growth underflows
   `DIAGRO→0` (Float32) ⇒ LN arg `=0` ⇒ `ALOG(0)`. Fix: guard arg≤0 → `DDS=-9.21` (the code's own floor
   sentinel). LS/CS only (NE/SN dgf compute DDS differently).
2. **varmrt `TEMKIL/TEMSUM` div0** — `{cs,ne,ls,sn}/varmrt.f:162`: mortality-search `TEMSUM=0` (no killable TPA
   left) with `TEMKIL>0` ⇒ `+Inf`. Instrumented `TEMSUM=0.0, TEMKIL=0.067`. Fix: `IF(TEMSUM.LE.0.) {ADJUST=1;
   GO TO 110}` = apply zero additional mortality. All 4 variants (line identical).
3. **htdbh H→D inverse `ALOG(<0)`** — `{cs,ne,ls,sn}/htdbh.f` Curtis-Arney branch: a tree grown TALLER than its
   asymptote `4.5+P2` makes `ratio=(ln(H-4.5)-ln(P2))/(-P3)<0` ⇒ `ALOG(negative)`. Instrumented `H=114.7 >
   asymptote 85, ratio=-0.0116`. Fix mirrors FVSjl volume.jl:92-98 exactly: clamp `MIN(H-4.5,0.9999*P2)` +
   `ratio>0 ? EXP(ALOG(ratio)/P4) : 100`. All 4 variants. (`fvs:197`/CRATET was a STALE record — runs clean.)

## r9clark >20-log cap zeroes ALL volume (incl. total cubic) for extreme-height trees — FIXED (2026-07-19)
Full trace: this section. Patch: docs/patches/nvel_r9clark_extremeheight_zerovol.patch (NVEL submodule
volume/NVEL/r9clark.f — the FMSC National Volume Estimator Library, shared by all Region-9 variants CS/LS/NE).
This is the "extreme-height volume-zeroing" class (2,140 sweep stands, ~0.5% of non-bit-exact; FVSjl is the
correct side and was cornered as an FVS bug).

ROOT CAUSE (instrumented, NE stand 1167721956290487, VOLEQ 900CLKE): for a very tall tree the merchantable
stem needs MORE than 20 logs, overflowing the fixed 20-slot log arrays (LOGLEN/LOGDIA/LOGVOL). r9logs.f flags
this with ERRFLG=12 and returns. Back in r9clark.f the error handler (line 373, DO 525; and the prod!=1
handler DO 575) does `DO I=1,15: VOL(I)=0.0; RETURN` — discarding the TOTAL and MERCH cubic volume too.
But those cubic volumes are r9cuft TAPER INTEGRALS computed BEFORE the log bucking and are perfectly valid
(instrumented: a 33-38" x ~300ft tree had vol(1)=536-680 cuft, vol(4)=529-674, vol(7)=5 right before the
zeroing). Only the board-foot LOG detail is invalid when >20 logs. Net effect on the .sum: TCuFt/MCuFt/SCuFt
collapse to 0 once the stand's tallest trees exceed the ~20-log height (e.g. TopHt>~240 ft), even though the
trees are large and carry hundreds of cuft each.

FIX (minimal, both handlers): when ERRFLG.EQ.12, set VOL(2)=0 (board-foot), clear ERRFLG, and RETURN with the
cubic volumes intact — instead of zeroing all 15 slots. No-op on every tree that buckets into <=20 logs
(ERRFLG=0 skips the whole handler). Board foot is reported 0 for the >20-log giants (they cannot be bucked
into the fixed arrays — an acceptable degradation; these heights are already beyond the model's realistic range
via the AVHT40-inflated height model).

VALIDATION: NE 1167721956290487 TCuFt now 17744/20591/23208/25413 at 2050-2080 (was 13555/217/0/0), sensible &
monotonic, and within ~4% of FVSjl (the correct side; residual = the normal AVHT40/self-thin tie-break). Normal
NE stands 12/12 byte-identical patched-vs-official. CONFIRMED SAME BUG ON LS (499580541126144, red pine,
fia125): live TCuFt collapsed 12924(2036)->12047(2046)->2832(2056)->4(2066) as TopHt passed ~200ft; the fix
recovers it to 12432 (within ~9% of FVSjl), 12/12 normal LS stands byte-identical. So the memory's "LS live=half"
was just an earlier collapse cycle, NOT a separate issue. In-container gfortran 12.2.0 relink (validation vehicle);
deliverable = the NVEL source patch. Applies to CS/LS/NE (shared FMSC NVEL library; CS latent-same). Build caveat #5.
