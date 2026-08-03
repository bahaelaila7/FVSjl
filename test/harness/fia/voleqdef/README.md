# R4_EQN authoritative dump — validation fixture for `r4_voleq`

`r4_voleq(FORNUM, FIA)` in `src/engine/volume_equations.jl` is a hand-port of the NVEL
`voleqdef.f` `R4_EQN` subroutine (Region-4 Intermountain forest-keyed volume-equation
assignment), shared by the region-4 western variants (CI/UT/TT).

`dump_r4_eqn.f90` is a module-free driver that calls the live `R4_EQN` directly (link
against any `bin/FVS??_buildDir/voleqdef.o`) for all forests × species and prints the
authoritative equation ids. `r4_eqn_authoritative.txt` is its output.

Reproduce / re-validate:
```
gfortran-16 -o dumpr4 dump_r4_eqn.f90 \
    /workspace/ForestVegetationSimulator/bin/FVSut_buildDir/voleqdef.o
./dumpr4 > r4_eqn_authoritative.txt
```
Then compare each line's equation to `FVSjl.r4_voleq(fnum, spec)`.

**Validated 2026-08-03: 555 (FORNUM,FIA) combinations, 0 diffs** — `r4_voleq` reproduces
the authoritative R4_EQN exactly. This is the robust equivalent of CR's
`data/centralrockies/volume_equations_by_forest.csv` (dumped from the same voleqdef.o).
