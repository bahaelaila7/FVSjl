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
include("core/coefficients.jl")
include("core/state.jl")

# --- variants ---------------------------------------------------------------
include("variants/southern/southern.jl")
include("variants/southern/coefficients.jl")
include("variants/southern/species.jl")
include("variants/southern/forest_location.jl")
include("variants/southern/habitat.jl")
include("variants/southern/site_index.jl")
include("variants/southern/bark_and_bounds.jl")
include("variants/southern/serial_correlation.jl")
include("variants/southern/diameter_growth.jl")
include("variants/southern/height_growth.jl")
include("variants/southern/small_tree_growth.jl")
include("variants/southern/crown_ratio.jl")
include("variants/southern/mortality.jl")

# --- variants: northeast (NE) — skeleton; equations + data ported chunk by chunk ---
include("variants/northeast/northeast.jl")
include("variants/northeast/species.jl")
include("variants/northeast/site_index.jl")
include("variants/northeast/diameter_growth.jl")    # NE large-tree DG (B1/B2/B3 + BAL competition)

# --- io ---------------------------------------------------------------------
include("io/treedata.jl")
include("io/keyword.jl")
include("io/csv_trees.jl")
include("io/yaml_keywords.jl")
include("io/yaml_stand.jl")
include("io/input.jl")

# --- engine -----------------------------------------------------------------
include("engine/species_translation.jl")
include("engine/treeinput.jl")
include("engine/keyword_dispatch.jl")
include("engine/init.jl")
include("engine/crown_width.jl")
include("engine/forest_type.jl")
include("engine/r8clark_vol.jl")
include("engine/r9clark_vol.jl")        # NE: NVEL Region-9 Clark profile volume
include("engine/volume_equations.jl")
include("engine/volume.jl")
include("engine/standstats.jl")
include("engine/event_monitor.jl")
include("engine/cuts.jl")
include("engine/quickersort.jl")        # RDPSRT/IQRSRT (Scowen 1965) — faithful sorts for COMPRESS
include("engine/compress.jl")            # COMPRESS (act 250) — PC-score tree-record clustering
include("engine/structure_stage.jl")    # SSTAGE — stand structural-stage class (1-6)
include("engine/establishment.jl")
include("engine/sprout.jl")
include("engine/fire/biomass.jl")        # FFE F1 — Jenkins tree biomass
include("engine/fire/crown_biomass.jl")  # FFE F2 — crown biomass by size class (FMCROWE)
include("engine/fire/fuel_loading.jl")   # FFE F3 — initial surface fuel loading (FMCBA)
include("engine/fire/fmcba.jl")          # FFE F3 — per-cycle fuel & cover-type update (FMCBA)
include("engine/fire/fuel_decay.jl")     # FFE F3 — per-cycle surface-fuel decay (FMCWD)
include("engine/fire/fuel_additions.jl") # FFE F3 — annual fuel additions / litterfall (FMCADD)
include("engine/fire/fire_effects.jl")   # FFE F6 — fire-caused tree mortality (FMEFF/FMBRKT)
include("engine/fire/rothermel.jl")      # FFE F5 — Rothermel surface fire behavior (FMFINT)
include("engine/fire/fuel_moisture.jl")  # FFE F5b — fuel-moisture scenario + wind reduction (FMMOIS)
include("engine/fire/fuel_model.jl")     # FFE F5b — dynamic fuel-model construction (FMCFMD3)
include("engine/fire/fmburn.jl")         # FFE F5b — fire event driver (FMBURN/FMEFF) → kill TPA
include("engine/fire/carbon.jl")         # FFE F8 — standing live-tree carbon pools (FMCRBOUT)
include("engine/fire/snag.jl")           # FFE F7 — snag falldown + decay dynamics (FMSFALL)
include("engine/fire/consumption.jl")    # FFE F7/F8 — fire fuel consumption + carbon release (FMCONS)
include("engine/econ.jl")                # C8 — ECON economic-analysis core (eccalc.f)
include("io/summary.jl")
include("io/dbs_output.jl")
include("engine/simulate.jl")

# --- more engine, extensions, cli are added in later chunks -----------------
# include("engine/...")    # C2–C5
# include("extensions/...")# C6–C8
# include("cli.jl")        # C8

export StandState, Southern, Northeast, AbstractVariant, variant_code
export load_species_coefficients!, init_blockdata!
export resolve_species, translate_species
export FVSRng, rann!, esrann!, bachlo, TreeList, ntrees
export parse_tree_format, parse_tree_record, read_tree_file, TreeRecord, DEFAULT_TREE_FORMAT
export KeywordReader, read_keyword!, KeywordRecord, KeywordStatus, KW_OK, KW_EOF, KW_STOP, KW_PARMS
export read_tree_records, read_trees_csv, write_trees_csv, convert_tre_to_csv, TREE_CSV_HEADER
export read_keyword_records, read_keywords_yaml, write_keywords_yaml, convert_key_to_yaml, read_keyfile_records
export write_keyfile, write_tree_file, convert_yaml_to_key, convert_csv_to_tre, translate_io
export initialize, initialize!, each_stand, run_keyfile, process_keywords!, load_trees!, strip_key_ext
export dgcons!, dgf!, notre!, stand_tpa, stand_ba, stand_qmd, stand_sdi, stand_ccf, stand_top_height

end # module FVSjl
