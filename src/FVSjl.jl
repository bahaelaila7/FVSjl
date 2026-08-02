"""
    FVSjl

An idiomatic, maintainable, thread-safe Julia reimplementation of the USFS Forest
Vegetation Simulator. It is a drop-in replacement for the live Fortran FVS — Southern
(`FVSsn`) and Northeast (`FVSne`) are complete and validated, Central States (`FVScs`)
in progress — same `.key`/`.tre` inputs, same SQLite / `.sum` outputs, and (in the
default `faithful=true` mode) bit-exact results.

Design (see docs/ARCHITECTURE.md):
  * all simulation state lives on an explicit `StandState` — no globals;
  * one state per stand/thread → safe parallelism;
  * pure numeric kernels + stateful orchestration → testable;
  * variants are dispatched via `AbstractVariant` singletons.

Ported from the Fortran sources under /workspace/ForestVegetationSimulator and
validated against the live Fortran builds (SN additionally vs the faithful port at
/workspace/FVSjulia; NE/CS have no faithful port — live binary is the sole oracle).
"""
module FVSjl

using Printf
using SQLite
using DBInterface

# --- core (order matters: parameters → rng/units/trees → variant → state) ----
include("core/fmath.jl")                       # gfortran-identical exp/log/pow (doctrine #8)
using .FMath: fexp, flog, fpow, fexp_julia, flog_julia, fpow_julia
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

# CS singleton is registered BEFORE the NE methods because NE's shared TWIGS crown method
# dispatches on `Union{Northeast,CentralStates}` (CS reuses it). The CS *methods* land after.
include("variants/centralstates/centralstates.jl")  # CS singleton + registration (MAXSP 96)
include("variants/lakestates/lakestates.jl")         # LS singleton + registration (MAXSP 68)
include("variants/centralrockies/centralrockies.jl") # CR singleton + registration (MAXSP 38) — first WESTERN variant
include("variants/centralrockies/species.jl")
include("variants/centralrockies/site_index.jl")
include("variants/centralrockies/diameter_growth.jl")
include("variants/centralrockies/height_growth.jl")
include("variants/centralrockies/crown.jl")
include("variants/centralrockies/small_tree_growth.jl")
include("variants/centralrockies/mortality.jl")
include("variants/centralrockies/dwarf_mistletoe_model.jl")  # CR dwarf mistletoe per-cycle spread (mistoe.f)

# --- kootenai (KT) — first of the western Rockies cluster; scaffold, equations chunk by chunk ---
include("variants/kootenai/kootenai.jl")             # KT singleton + registration (MAXSP 11, western Wykoff DDS)
include("variants/kootenai/species.jl")              # KT blkdat init (11 species, seed 55329, Stage SDI) + SPCTRN col 12
include("variants/kootenai/habitat_tables.jl")       # KT habtyp/sitset lookup tables (KOTHAB/JTYPE/KTYPE/MTYPE/BAMAXA)
include("variants/kootenai/site_index.jl")           # KT habtyp 2-level mapping + sitset SDImax (site_setup!)
include("variants/kootenai/dg_coefficients.jl")      # KT large-tree DDS coefficient arrays (kt/dgf.f DATA)
include("variants/kootenai/crown.jl")                # KT per-tree CCF (ccfcal MODE=1) → RELDEN
include("variants/kootenai/diameter_growth.jl")      # KT dgf! (Wykoff DDS) + kt_dgcons!
include("variants/kootenai/height_growth.jl")        # KT htgf! (exp-form height increment)
include("variants/kootenai/regent.jl")               # KT regent (chunk 6): kt_regcons! + small_tree_growth!
include("variants/kootenai/mortality.jl")            # KT mortality (chunk 7): Hamilton MORTS
include("variants/kootenai/volume.jl")               # KT volume (chunk 8): Region-1 Flewelling FW2

# --- inland empire (IE) — 3rd western variant, reuses the KT engine at a discount ---
include("variants/inlandempire/inlandempire.jl")     # IE singleton + registration (MAXSP 23)
include("variants/inlandempire/species.jl")          # IE species block-data init (ie/blkdat.f + grinit.f)
include("variants/inlandempire/site_index.jl")       # IE habtyp/forkod/sitset (chunk 2): ITYPE/IFOR/SITEAR
include("variants/inlandempire/dg_coefficients.jl")  # IE large-tree DDS coefficient arrays (ie/dgf.f)
include("variants/inlandempire/diameter_growth.jl")  # IE large-tree DDS hooks (ie_dgcons! + dgf!)
include("variants/inlandempire/height_growth.jl")    # IE large-tree height growth (ie/htgf.f)
include("variants/inlandempire/crown_coefficients.jl") # IE crown-ratio tables (ie/crown.f, dumped from live)
include("variants/inlandempire/crown.jl")            # IE per-tree CCF (ie/ccfcal.f) + crown ratio (ie/crown.f)
include("variants/inlandempire/mort_coefficients.jl")# IE mortality Hamilton coefficients (ie/morts.f, dumped)
include("variants/inlandempire/mortality.jl")        # IE mortality (ie/morts.f) — Hamilton, reuses KT form
include("variants/inlandempire/mistoe_coefficients.jl") # IE MISTOE (dwarf mistletoe) effect coefficients (mistoe/misintie.f)
include("variants/inlandempire/establishment.jl")      # IE ESSUBH subsequent/planted-tree height model (ie/essubh.f)
include("variants/inlandempire/volume.jl")           # IE volume (ie/sitset.f VOLEQ): FW2 + DVE, reuse CR kernels
include("variants/inlandempire/regent_coefficients.jl") # IE REGENT small-tree coefficients (ie/regent.f, dumped)
include("variants/inlandempire/regent.jl")           # IE REGENT small-tree growth (ie/regent.f): ie_regcons! + hook

# --- easternmontana (EM) — next western Rockies cluster variant (KT-engine discount); chunk 0 scaffold ---
include("variants/easternmontana/easternmontana.jl")  # EM singleton + registration (MAXSP 19, western Wykoff DDS)
include("variants/easternmontana/species.jl")         # EM blkdat init (19 species, seed 55329, Stage SDI) + SPCTRN col 10
include("variants/easternmontana/site_index.jl")      # EM habtyp (JTYPE/NIHMAP→ITYPE) + sitset (site index/SDIDEF) — site_setup!
include("variants/easternmontana/crown.jl")           # EM per-tree CCF (em/ccfcal.f MODE=1 polynomial) — em_tree_ccf
include("variants/easternmontana/dg_coefficients.jl") # EM large-tree DG coefficient arrays (em/dgf.f DATA, 1D+2D)
include("variants/easternmontana/diameter_growth.jl") # EM large-tree DDS (em/dgf.f): em_dgcons! + dgf! (Wykoff main)
include("variants/easternmontana/height_growth.jl") # EM large-tree height growth (em/htgf.f+pothtg.f): POTHTG+RALPH+modifier
include("variants/easternmontana/mortality.jl") # EM mortality (em/morts.f): RI-background + SDI-trend RN
include("variants/easternmontana/regent.jl") # EM small-tree growth (em/regent.f NIVAR): em_regcons! + small_tree_growth!

# --- teton (TT) — next western Rockies cluster variant (EM-DG + CR-Zeide-density discount); chunk 0 scaffold ---
include("variants/teton/teton.jl")            # TT singleton + registration (MAXSP 18, western Wykoff DDS, oracle verified)

# --- variants: northeast (NE) — skeleton; equations + data ported chunk by chunk ---
include("variants/northeast/northeast.jl")
include("variants/northeast/species.jl")
include("variants/northeast/site_index.jl")
include("variants/northeast/diameter_growth.jl")    # NE large-tree DG (B1/B2/B3 + BAL competition)
include("variants/northeast/height_growth.jl")      # NE height growth (NC-128 curve + BAL modifier)
include("variants/northeast/mortality.jl")          # NE VARMRT efficiency (relative-height + VARADJ)
include("variants/northeast/crown_ratio.jl")        # NE crown ratio (TWIGS BA-based model)
include("variants/northeast/small_tree_growth.jl")  # NE small-tree growth (REGENT, dbh<5)

# --- centralstates (CS) methods (singleton already registered above, before NE) -------------
include("variants/centralstates/species.jl")        # CS blkdat init (codes/RNG/YR=10/Zeide SDI)
include("variants/centralstates/site_index.jl")     # CS SITSET (per-species ASITE/BSITE linear)
include("variants/centralstates/diameter_growth.jl")# CS DG: cs_dgf! ln(DDS) regression + cs_dgcons!
include("variants/centralstates/height_growth.jl")  # CS HTGF: NC-128 (MAPCS) + cs_balmod (reuses LTBHEC)
include("variants/centralstates/small_tree_growth.jl") # CS REGENT (d<5): NE shape, XMIN=3, cs_balmod
# CS HT-DBH rides the Southern Curtis-Arney+Wykoff path (cs/htdbh.f ≡ sn/htdbh.f); CS crown
# rides NE's TWIGS method (see northeast/crown_ratio.jl) — both are coefficient-driven, no CS code.

# --- lakestates (LS) methods (singleton already registered above, before NE) -----------------
include("variants/lakestates/species.jl")           # LS blkdat init (codes/RNG/YR=10/Zeide SDI) + spctrn col 5
include("variants/lakestates/site_index.jl")         # LS SITSET (68×68 SICOEF fan-out + SDICON) + FORKOD
include("variants/lakestates/diameter_growth.jl")    # LS DG: ls_dgf! ln(DDS) (== CS model, LS caps) + ls_dgcons!
include("variants/lakestates/height_growth.jl")      # LS HTGF: NC-128 (MAPLS/IVAR=1) + ls_balmod competition
include("variants/lakestates/small_tree_growth.jl")  # LS REGENT (d<5): NE/CS shape, XMIN=3, MAPLS + ls_balmod

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
include("io/fia_database.jl")            # DATABASE/DSNIN input: FIA "FVS-ready" SQLite → stand
include("io/fia_translate.jl")           # raw FIADB (PLOT/COND/TREE/…) → FVS-ready records
include("engine/init.jl")
include("engine/crown_width.jl")
include("engine/forest_type.jl")
include("engine/r8clark_vol.jl")
include("engine/r9clark_vol.jl")        # NE: NVEL Region-9 Clark profile volume
include("engine/cr_dve_vol.jl")         # CR: NVEL R3 D2H (DVE) volume
include("engine/cr_nvb_vol.jl")         # CR: NVEL NSVB (NVB) volume
include("engine/cr_fw2_vol.jl")         # CR: NVEL Flewelling (FW2) stem-profile volume
include("engine/r9vol_gevorkiantz.jl")  # CS: NVEL Region-9 Gevorkiantz '900DVEE' (VOLUME METHC=5)
include("engine/volume_equations.jl")
include("engine/volume.jl")
include("engine/standstats.jl")
include("engine/event_monitor.jl")
include("engine/cuts.jl")
include("engine/quickersort.jl")        # RDPSRT/IQRSRT (Scowen 1965) — faithful sorts for COMPRESS
include("engine/compress.jl")            # COMPRESS (act 250) — PC-score tree-record clustering
include("engine/structure_stage.jl")    # SSTAGE — stand structural-stage class (1-6)
include("../data/centralrockies/establishment/estab_coefs.jl")   # CR establishment coefs (XMIN/HHTMAX/ESSUBH-HHT)
include("../data/centralrockies/dwarf_mistletoe.jl")             # CR dwarf mistletoe coefs + pure kernels (mistoe/misintcr.f)
include("engine/establishment.jl")
include("engine/sprout.jl")
include("engine/fire/biomass.jl")        # FFE F1 — Jenkins tree biomass
include("engine/fire/crown_biomass.jl")  # FFE F2 — crown biomass by size class (FMCROWE)
include("engine/fire/cr_crown_biomass.jl") # FFE F2 — CR western crown biomass (FMCROWW)
include("../data/centralrockies/fire/ffe_live_fuel.jl")  # CR FFE FULIVE/FULIVI live-fuel data
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

export StandState, Southern, Northeast, CentralStates, LakeStates, CentralRockies, Kootenai, InlandEmpire, EasternMontana, Teton, AbstractVariant, variant_code, variant_from_code
export load_species_coefficients!, init_blockdata!
export resolve_species, translate_species
export FVSRng, rann!, esrann!, bachlo, TreeList, ntrees
export parse_tree_format, parse_tree_record, read_tree_file, TreeRecord, DEFAULT_TREE_FORMAT
export KeywordReader, read_keyword!, KeywordRecord, KeywordStatus, KW_OK, KW_EOF, KW_STOP, KW_PARMS
export read_tree_records, read_trees_csv, write_trees_csv, convert_tre_to_csv, TREE_CSV_HEADER
export read_keyword_records, read_keywords_yaml, write_keywords_yaml, convert_key_to_yaml, read_keyfile_records
export write_keyfile, write_tree_file, convert_yaml_to_key, convert_csv_to_tre, translate_io
export write_sum_csv, SUM_CSV_HEADER, yaml_variant_code, yaml_output_format
export initialize, initialize!, each_stand, run_keyfile, process_keywords!, load_trees!, strip_key_ext
export dgcons!, dgf!, notre!, stand_tpa, stand_ba, stand_qmd, stand_sdi, stand_ccf, stand_top_height

end # module FVSjl
