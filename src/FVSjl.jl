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
include("variants/southern/species.jl")

# --- io ---------------------------------------------------------------------
include("io/treedata.jl")
include("io/keyword.jl")
include("io/csv_trees.jl")
include("io/yaml_keywords.jl")
include("io/input.jl")

# --- engine -----------------------------------------------------------------
include("engine/species_translation.jl")

# --- more engine, extensions, cli are added in later chunks -----------------
# include("engine/...")    # C2–C5
# include("extensions/...")# C6–C8
# include("cli.jl")        # C8

export StandState, Southern, AbstractVariant, variant_code
export load_species_coefficients!, init_blockdata!
export resolve_species, translate_species
export FVSRng, rann!, esrann!, bachlo, TreeList, ntrees
export parse_tree_format, parse_tree_record, read_tree_file, TreeRecord, DEFAULT_TREE_FORMAT
export KeywordReader, read_keyword!, KeywordRecord, KeywordStatus, KW_OK, KW_EOF, KW_STOP, KW_PARMS
export read_tree_records, read_trees_csv, write_trees_csv, convert_tre_to_csv, TREE_CSV_HEADER
export read_keyword_records, read_keywords_yaml, write_keywords_yaml, convert_key_to_yaml, read_keyfile_records

end # module FVSjl
