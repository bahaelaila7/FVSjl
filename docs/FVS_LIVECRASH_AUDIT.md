# FVS live-crash fix — working audit log

Campaign: root-cause + minimally patch every live-FVS SIGFPE on real FIA inventory. Goal/doctrine:
docs/FVS_LIVECRASH_FIX_GOAL.md. Stands: docs/fvs_livecrash_stands.txt (12, recorded 6 sites). FVSjl runs
all of them clean — it is the correct side. Patches: docs/patches/livecrash_*.patch.

## Resolution summary (2026-07-18)
The 12 recorded crash sites collapse to **3 root-cause classes** (the backtrace attributed several to the
*caller* frame because the callee `.o` had no `-g`): a diameter-growth `ALOG(0)`, a mortality-search
divide-by-zero, and a height→DBH inverse `ALOG(negative)`. All three are unguarded arithmetic on
degenerate-but-legal FIA geometry. Each is fixed with a minimal guard that is a **no-op on every
non-degenerate input** and returns the same value FVSjl (the correct side) produces on the degenerate one.

| recorded site(s) | stands | root class | fix file(s) |
|------------------|-------:|-----------|-------------|
| cs/grincr.f:437, ls/dgdriv.f:134, ls/dgdriv.f:353 | 1+3+2 | **DGF `ALOG(0)`** (both `dgdriv:134/:353` are `CALL DGF`) | cs/dgf.f, ls/dgf.f |
| cs/varmrt.f:162 | 1 | **varmrt `TEMKIL/TEMSUM` div0** | cs/ne/ls/sn varmrt.f |
| cs/grincr.f:449, ls/htdbh.f:336 | 3+1 | **htdbh H→D inverse `ALOG(<0)`** (`grincr:449` is `CALL REGENT`→`regent:324` `CALL HTDBH`) | cs/ne/ls/sn htdbh.f |
| cs/fvs.f:197 | 1 | **STALE** — official oracle runs it clean (rc=0); not a live crash | — |

## Build-caveat honesty (doctrine #5)
The shipped oracle build traps FPE (SIGFPE on zero/invalid). In-container relinks use gfortran **12.2.0**
(≠ the official **15.2.1**). Two consequences, both handled:
1. The `ALOG(0)`/underflow triggers are **Float32-precision-sensitive** — whether `DIAGRO` underflows to
   *exactly* 0, or a tree's grown height *exactly* exceeds its asymptote, depends on the compiler's FP
   rounding. Under 12.2.0 some stands that crashed the shipped build now run clean (the invalid value is
   just barely avoided). Reproduction is therefore **not** a reliable validator on its own.
2. So each bug is proven by **direct instrumentation of the degenerate value** (below), independent of
   whether a given relink happens to trap — plus the patched output is shown to match FVSjl bit-exact-or-
   cornered, and normal stands stay byte-identical.

---

## Class 1 — DGF `ALOG(0)`  (cs/dgf.f, ls/dgf.f)
**Sites:** ls/dgdriv.f:134 (`CALL DGF(DBH)`), ls/dgdriv.f:353 (`CALL DGF(WK3)`), cs/grincr.f:437 — all resolve
into `DGF`. LS canonical line 456; CS 550.
**Unguarded op:**
```fortran
      DDS=ALOG(((DBH(I)*BARK+DIAGRI)**2.0)-(DBH(I)*BARK)**2.0)
```
**Trigger (instrumented, all 5 LS DGF-class stands):** `DBH=66.9, DIAGRO=0.0, DIAGRI=0.0, EXP(DDS)=2.4e-4`.
A LARGE tree with TINY predicted growth: `DIAGRO = SQRT(DBH²+EXP(DDS))−DBH` underflows to **exactly 0** in
Float32 because `EXP(DDS)=2.4e-4` is below the ULP of `DBH² = 4476` (ULP ≈ 2⁻¹¹ ≈ 4.9e-4) ⇒ `SQRT(4476)=66.9=
DBH` ⇒ `DIAGRO=0` ⇒ `DIAGRI=0` ⇒ the LN argument `(DBH·BARK+0)² − (DBH·BARK)² = 0` ⇒ `ALOG(0)` ⇒ SIGFPE.
(Directly re-confirmed this session — and notably the 12.2.0 build **does** compute `arg=0` here yet did not
trap on the full run, i.e. the degenerate op genuinely occurs; only whether it SIGFPEs is build/trap-dependent
— see the build caveat above.)
**Fix (minimal, faithful):** floor the degenerate argument to the code's own `-9.21` sentinel (the value the
existing `IF(DDS.LT.-9.21) DDS=-9.21` cap 6 lines below already forces for a vanishingly-small DDS):
```fortran
      DDS=((DBH(I)*BARK+DIAGRI)**2.0)-(DBH(I)*BARK)**2.0
      IF(DDS.LE.0.0)THEN
        DDS=-9.21
      ELSE
        DDS=ALOG(DDS)
      ENDIF
```
No-op when arg>0 (the `ELSE` reproduces the original `ALOG`). Variant scope: only LS/CS `dgf.f` carry this
"outside-bark→inside-bark" DDS formulation; NE/SN `dgf.f` compute DDS differently and have no such line.
**Validation:** patched LS `.sum` bit-exact to FVSjl on 167130809020004 (earlier slice). All 6 LS+CS
DGF-class stands run clean under the patched build.

## Class 2 — varmrt `TEMKIL/TEMSUM` div0  (cs/ne/ls/sn varmrt.f)
**Site:** varmrt.f:162 `ADJUST=TEMKIL/TEMSUM`.
**Trigger (instrumented, unpatched CS on 103630089010661):** `TEMSUM = 0.0, TEMKIL = 0.0672`. In the Newton
mortality-target search, `TEMSUM` accrues only over records with killable TPA (`TPALFT = PROB(I)−WK2(I) > 0`).
When the pass has already consumed all killable TPA, `TEMSUM = 0` while `TEMKIL > 0` ⇒ `TEMKIL/TEMSUM = +Inf`
⇒ SIGFPE.
**Fix:**
```fortran
      IF(TEMSUM .LE. 0.)THEN
        ADJUST=1.0
        GO TO 110
      ENDIF
      ADJUST=TEMKIL/TEMSUM
```
Semantics: no killable TPA ⇒ apply zero additional mortality (all `TEMWK2(I)` are 0 in the apply loop at
label 110, so `ADJUST`'s value is immaterial there; `GO TO 110` just exits the search). No-op when TEMSUM>0.
Variant scope: the line is textually identical in all 4 variants (md5 4129d6af) ⇒ patched in all 4 (CS hit;
NE/LS/SN latent-vulnerable).

## Class 3 — htdbh H→D inverse `ALOG(negative)`  (cs/ne/ls/sn htdbh.f)
**Sites:** ls/htdbh.f:336; cs/grincr.f:449 (`CALL REGENT`) → regent.f:324 (`CALL HTDBH(...,HK,1)`) → the same
`htdbh` inverse. Curtis-Arney branch (`IWYKCA≠0`):
```fortran
      D=EXP(ALOG((ALOG(H-4.5)-ALOG(P2))/(-1.*P3)) * 1./P4)
```
**Trigger (instrumented, LS 1867797854290487):** `H=114.7 ft` but the species asymptotic height is
`4.5+P2 = 85 ft` (P2=80.5, P3=26.98, P4=−2.02). Inverting `H = 4.5 + P2·EXP(−P3·D^P4)` for a tree taller than
its own asymptote gives `ratio = (ALOG(H−4.5)−ALOG(P2))/(−P3) = −0.0116 < 0` ⇒ the outer `ALOG(ratio)` =
`ALOG(negative)` = invalid ⇒ SIGFPE. (Confirmed H>4.5, so it is the *outer* ALOG, not `ALOG(H−4.5)`.)
**Fix (mirrors FVSjl volume.jl:92-98 exactly):**
```fortran
      D=(ALOG(MIN(H-4.5,0.9999*P2))-ALOG(P2))/(-1.*P3)
      IF(D .GT. 0.)THEN
        D=EXP(ALOG(D) * 1./P4)
      ELSE
        D=100.
      ENDIF
```
Clamp `H−4.5` just below the asymptote `P2` (so the inverse is defined for above-asymptote trees, giving a
large-but-finite DBH); fall back to `D=100"` if the ratio is still non-positive. No-op for every below-
asymptote tree (`MIN` returns `H−4.5` unchanged and `ratio>0`, reproducing the original expression). The
Wykoff branch (`:331/:463`) is deliberately **not** guarded — FVSjl leaves it unguarded too (the small-tree
caller guarantees `hk>4.5` and a `db_floor` catches negatives), so guarding it would diverge from the port.
Variant scope: all 4 variants carry the Curtis-Arney line ⇒ patched in all 4 (CS+LS hit; NE/SN latent).

---

## Validation (doctrine #3)
- **Crash-clean:** all 12 recorded stands run to normal completion under the patched builds (no rc=136); the
  `fvs:197` stand was already clean on the shipped oracle (stale record). **Full re-sweep of all 60 recorded
  `live_crash` stands (CS 29 + NE 10 + LS 21) under the patched builds: 60/60 crash-free, 0 SIGFPE.**
- **Matches FVSjl (correct side):** `diff_one` patched-live vs FVSjl on the 4 reproducing CS stands (thinbba):
  103630089010661 **bit-exact all 6 cycles × 6 cols**; the 3 regent/htdbh stands bit-exact early then diverge
  only at late cycles by ULP-scale amounts = the **same dense-stand RDPSRT self-thinning tie-break primitive**
  already cornered in the FIA campaign (docs/fia_divergence_taxonomy.md), NOT a fix artifact.
- **No regression:** patched-vs-official normal-stand summaries **numeric byte-identical** across all 4
  variants — **CS 15/15, NE 12/12, LS 12/12, SN 12/12 (51/51)**. The guards are no-ops off the degenerate
  path, by construction and empirically.

## Oracle state
The shipped `tmp/oracles/FVS{cs,ne,ls,sn}_new` are left **pristine** (official build). The deliverable is the
canonical source patch set (`docs/patches/*.patch`, applied to `ForestVegetationSimulator/{variant}/*.f`,
`git diff`-clean, CRLF-preserving) + this reasoning. Validation oracles were built at `/tmp/FVS*_patched`.
