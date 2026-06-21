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
