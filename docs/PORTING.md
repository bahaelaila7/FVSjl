# Porting guide: Fortran FVS → idiomatic FVSjl

How to translate a subroutine (or a whole new variant) into FVSjl. The source of
truth is `/workspace/ForestVegetationSimulator` (per-variant dirs `base/ common/
sn/ cs/ ne/ ls/`); the already-validated faithful port at `/workspace/FVSjulia`
is the secondary reference (Oracle A) for both behaviour and exact constants.

## Naming conventions

- **Functions**: `snake_case`, mutating ones end in `!`. Keep the Fortran routine
  name in the doc-comment (`# Ported from: sn/dgf.f`).
- **Struct fields**: readable names; **every field carries its original COMMON
  name in a comment** so the two codebases stay cross-referenceable.
- **Pure kernels vs orchestration**: extract the numeric formula into a pure
  function; keep state mutation in the `!` orchestrator. See ARCHITECTURE.md.

## Type rules (bit-exactness)

| Fortran            | Julia        |
|--------------------|--------------|
| `REAL` / `REAL*4`  | `Float32`    |
| `DOUBLE PRECISION` | `Float64`    |
| `INTEGER`          | `Int32`      |
| `LOGICAL`          | `Bool`       |
| `CHARACTER*n`      | `String` (fixed width; pad/trim at edges) |

- Use `Float32` literals (`2.5f0`) in kernels — a stray `Float64` changes results.
- `NINT(x)` → `nint(x)` (`RoundNearestTiesAway`); plain `round` is banker's and is wrong.
- Integer/real division: Fortran `/` on integers truncates → use `div`.

## Control flow

- `GOTO`/labels → restructure into `for`/`while`/`break`/`continue`/early `return`
  or a helper function. **Do not** transliterate with `@goto`; the whole point of
  this rewrite is to remove the 1,131 gotos of the faithful port. If a block is
  genuinely a state machine, write it as one explicitly.
- Computed `GOTO (a,b,c) i` → `if/elseif` or dispatch.
- `ENTRY` points → separate functions sharing a private helper.
- `SAVE` locals → fields on the appropriate state struct (never a Julia global).

## State access

Everything that was a COMMON variable is now a field reached through the
`StandState` passed in. Don't reach for module globals — there are none. If a
routine needs new persistent state, add a field to the matching sub-struct (and
its constructor) rather than introducing a global.

## Constants & coefficient tables

Variant species-coefficient `DATA` tables (e.g. `DGF_INTERC`) become module-level
`const` arrays in the variant directory (`src/variants/southern/species.jl`).
Per-stand *derived* coefficients (calibrated/site-adjusted) live on the state.
Copy the exact float constants from Oracle A to avoid transcription error.

## Adding a new variant (CS/NE/LS)

1. `struct CentralStates <: AbstractVariant end`; `variant_code(::CentralStates)="CS"`.
2. Implement the hooks in `src/variants/variant.jl` for the new singleton
   (`diameter_growth!`, `height_growth!`, `mortality!`, `regenerate!`,
   `load_species_coefficients!`, ...).
3. Add the variant's species/coefficient tables under `src/variants/centralstates/`.
4. The `engine/` is untouched. Add integration baselines under `test/`.

## Verifying a port

Diff against Oracle A at the finest granularity available:
- pure kernel → unit-test the function on representative inputs vs the old routine;
- a phase (growth/mortality) → diff per-tree arrays after that phase;
- a full run → `Oracle.assert_sum_matches("snt01.key")` and DB-table compare.

Bit-exact is the bar in faithful mode. Genuinely unavoidable transcendental ulp
drift is documented in DIVERGENCES.md, not chased.
