"""
    FVSjl

An idiomatic, maintainable, thread-safe Julia reimplementation of the USFS Forest
Vegetation Simulator (Southern variant to start). It is a drop-in replacement for
the Fortran `FVSsn`: same `.key`/`.tre` inputs, same SQLite / `.sum` outputs, and
(in the default `faithful=true` mode) bit-exact results.

Design (see docs/ARCHITECTURE.md):
  * all simulation state lives on an explicit `StandState` — no globals;
  * one state per stand/thread → safe parallelism;
  * pure numeric kernels + stateful orchestration → testable;
  * variants are dispatched via `AbstractVariant` singletons.

Ported from the Fortran sources under /workspace/ForestVegetationSimulator and
validated against the faithful port at /workspace/FVSjulia (see test/).
"""
module FVSjl

using Printf
using SQLite
using DBInterface

# --- core (order matters: parameters → rng/units/trees → variant → state) ----
include("core/parameters.jl")
include("core/rng.jl")
include("core/units.jl")
include("core/trees.jl")
include("variants/variant.jl")
include("core/state.jl")

# --- variants ---------------------------------------------------------------
include("variants/southern/southern.jl")

# --- engine, io, extensions, cli are added in later chunks ------------------
# include("io/...")        # C1
# include("engine/...")    # C2–C5
# include("extensions/...")# C6–C8
# include("cli.jl")        # C8

export StandState, Southern, AbstractVariant, variant_code
export FVSRng, rann!, esrann!, bachlo, TreeList, ntrees

end # module FVSjl
