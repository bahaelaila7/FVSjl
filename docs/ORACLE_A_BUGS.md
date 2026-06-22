# Bugs found in Oracle A (`/workspace/FVSjulia`) that need fixing there

Oracle A is the faithful 1:1 Julia transliteration of Fortran `FVSsn`, used as the
primary reference for validating FVSjl. It is normally trustworthy, but the FVSjl
work has surfaced places where **Oracle A diverges from the real Fortran** — i.e.
transliteration bugs in FVSjulia itself, not in FVSjl. They are recorded here so
they can be fixed in the FVSjulia repo, and so FVSjl validation does not blindly
trust Oracle A on these paths.

Ground truth for each is the original Fortran in
`/workspace/ForestVegetationSimulator_fresh/` (or `.../ForestVegetationSimulator/`).

Format per entry:
> ### <title>
> - **Oracle A (FVSjulia):** file:line — what it does (wrong)
> - **Real Fortran:** file.f:line — what it should do
> - **FVSjl:** correct? (and where)
> - **Impact / how found**
> - **Fix**

---

### HBDECD alpha habitat-code decode drops the match (defaults everything)
- **Oracle A (FVSjulia):** `src/sn/habtyp.jl` `HBDECD` (the `ihb==0` alpha branch,
  ~lines 60–114). Its hand-rolled `temp` construction (double `IOBuffer` with a
  `take!`-based length check) does not reproduce the Fortran fixed-width
  `CHARACTER*8` semantics. For an alpha ecological-unit code such as `232BA`,
  `221HA`, `255AA`, `M222AA`, … it fails to match the SNECU table and falls through
  to the default index 122 (`231DD`).
- **Real Fortran:** `base/hbdecd.f`. `TEMP` and `CNHB(I)` are both `CHARACTER*8`,
  space-padded, and the match is `IF (TEMP(1:8).EQ.CNHB(I)(1:8))`. So `"232BA   "`
  (padded to 8) equals `CNHB` entry `"232BA   "` and **matches** → that ecological
  unit is selected (e.g. `232BA` → p232 physiography).
- **FVSjl:** **correct.** `variants/southern/habitat.jl` `resolve_eco_unit` matches
  the trimmed/uppercased code against SNECU exactly as the real Fortran does
  (verified: `"232BA"` → `232BA` → p232).
- **Impact / how found:** caught while closing the physiography test blindspot —
  alpha-code physiography scenarios resolved differently in Oracle A vs FVSjl. Only
  alpha habitat codes are affected; numeric SNECU-index codes (e.g. `232`→index 232)
  decode correctly in both. snt01 (`231Dd`) is unaffected because its default and
  its intended match are both `231DD`. The physiography test scenarios were
  therefore driven via the numeric-index path (where Oracle A is correct).
- **Fix:** in FVSjulia `HBDECD`, build `temp` as a fixed 8-char space-padded string
  from the first non-blank run of `KARD2` and compare `temp == rpad(code,8)[1:8]`
  (or compare trimmed-vs-trimmed), matching `hbdecd.f`.
- **STATUS: FIXED** in FVSjulia commit `ddfac81` (2026-06-22). Validated vs live
  Fortran: a snt01 stand with habitat `232BA` was 53 `.sum` rows wrong (defaulted to
  231DD), now bit-exact. Found via the FVSjl 3-way harness.

---

### R8 Clark `_r8_remap_spec` operator-precedence bug → zero volume for FIA 123/197
- **Oracle A (FVSjulia):** `src/base/r8clark_vol.jl:583` —
  `spec == 123 || spec == 197 && return 100`. Julia binds `&&` tighter than `||`, so
  this parses as `spec == 123 || (spec == 197 && return 100)`. For `spec==123` the
  expression short-circuits to a discarded `true`; the remap never fires, 123 is
  looked up directly, the Clark coefficient table has no 123 entry, and **cubic
  volume comes out 0**. (197 worked, via the `&&`.)
- **Real Fortran:** `base/r8prep.f` lines 99-116 — `IF (SPEC.EQ.123 .OR. SPEC.EQ.197)
  SPEC=100`; **both** codes remap to species 100's coefficients.
- **FVSjl:** had the **same** bug (inherited from this file) at
  `src/engine/r8clark_vol.jl:583`; fixed in parallel.
- **Impact / how found:** comprehensive 3-way sweep (162 scenarios × with/without
  management) — `all_TM` (Table Mountain pine, FIA 123) reported `cuft=0` every cycle
  in **both** Julia ports while live Fortran gave 1491/2184/2715/…. After the fix
  FVSjulia is bit-exact to Fortran on all_TM.
- **Fix:** parenthesize the `||`: `(spec == 123 || spec == 197) && return 100`.
- **STATUS: FIXED** — FVSjulia commit `06dbc1d`, FVSjl commit `8fa68ea` (2026-06-22).

---

## Open items surfaced by the 3-way sweep (not yet fixed)

- **FFE fire mortality residual (FVSjulia):** `fire_fuel9`/`fire_fuel11` match Fortran
  bit-exact pre-burn, then FVSjulia kills ~10-28 more TPA than Fortran at the burn
  cycle (e.g. fire_fuel9 2010: 133 vs 143 TPA, BA equal). A real but small FFE
  fire-effects divergence in Oracle A's fire extension. Needs tracing in SIMFIRE/FMEFF.
- **Physiography transcendental tail (shared):** `s06_ecounit_232`→p234,
  `s13_phys_p222`→p222, `s16_phys_p255`→p255 each show a small **shared**
  Julia-vs-Fortran cuft drift (±3-10 cuft, ±3-4 TPA) that **both** ports exhibit
  identically (so not a port-specific logic bug). cyc0 is exact; the drift appears in
  growth, is bounded/oscillating (not compounding), and tracks ±0.1 QMD sub-display
  flips → accumulated single-precision transcendental (`exp`/`sqrt`) rounding in the
  DG→volume chain. Pending a definitive live-Fortran per-tree `wk2`/DG trace to confirm.
- **Live Fortran FP-crashes (ground truth):** species `all_AE/EL/RL/SU/WE` abort the
  Fortran `FVSsn` binary itself with SIGFPE — a likely FVS bug or degenerate species
  coefficients; those species cannot be validated against Fortran.
