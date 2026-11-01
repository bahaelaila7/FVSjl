#!/usr/bin/env bash
# Launch the Pan-enabled Forest Explorer (branch forest-explorer-pan). Relocatable: every path derives from the workspace
# root, i.e. the folder that holds fe_pan/ (three levels above this script), or from $WORKSPACE if set. Expected layout
# (each entry can be overridden by its own variable, see PAN.md "Run"):
#   $WORKSPACE/fe_pan/apps/forest-explorer        this app
#   $WORKSPACE/FVSclaude/SQLite_FIADB_ENTIRE.db    FIA_DB         FVSjl master FIA DB (read only)
#   $WORKSPACE/FVSclaude/treemap                   TREEMAP_DATA   TreeMap raster/VAT/trees (read only)
#   $WORKSPACE/landis_fe                           PAN_REPO       Pan checkout, branch fe-engine (scenario engine incl. HBAD v2 = plugin hbad2)
#   $WORKSPACE/FIA/truth/fia_truth_r4.duckdb    PAN_FIA_DB     the single curated FIA truth for Pan (read only)
#   $WORKSPACE/FIA/truth/fia_truth_r4.duckdb    FIA_TRAJ_DB     FIA plots / tree-trajectories layer (read only; default PAN_FIA_DB)
#   $WORKSPACE/pan_climate/pan_climate.duckdb      PAN_CLIMATE_DB PRISM baseline + forecasts (PnET climate)
#   $WORKSPACE/landis_region4/data_r4/{plot_soil_r4.csv, pnet_climate_r4_l4.jld2, co2_annual.csv}
#                                                  PAN_SOIL_CSV, PAN_PNET_CLIMATE, PAN_CO2_CSV
#   $WORKSPACE/EPA_ECOREGIONS/us_eco_l4/us_eco_l4_no_st.shp   ECO_L4_SHP  EPA L4 polygons (AOI location → ecoregion for Pan)
#   $WORKSPACE/SOLUS/SOLUS100                      PAN_SOLUS      SOLUS soil rasters (PnET soils at the AOI location)
# Julia packages: private depots inside this app (.julia_depot_fe) and the Pan checkout (.julia_depot_wt), layered over the
# workspace's shared depots when those exist. Missing packages are installed automatically on first start (and again
# whenever a Manifest.toml changes); SETUP=1 forces it, SETUP=0 skips it.
# Usage: ./serve_pan.sh   (HOST/PORT/JULIA_NUM_THREADS/PAN_THREADS/WORKSPACE and the variables above override)
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ws="${WORKSPACE:-$(cd "$here/../../.." && pwd)}"

# --- data paths -------------------------------------------------------------------------------------------------------
export FIA_DB="${FIA_DB:-$ws/FVSclaude/SQLite_FIADB_ENTIRE.db}"
export TREEMAP_DATA="${TREEMAP_DATA:-$ws/FVSclaude/treemap}"
export FE_CACHE_DIR="${FE_CACHE_DIR:-$here/cache/sim_cache}"
export PAN_REPO="${PAN_REPO:-$ws/landis_fe}"
export PAN_FIA_DB="${PAN_FIA_DB:-$ws/FIA/truth/fia_truth_r4.duckdb}"
export FIA_TRAJ_DB="${FIA_TRAJ_DB:-$PAN_FIA_DB}"
export PAN_CLIMATE_DB="${PAN_CLIMATE_DB:-$ws/pan_climate/pan_climate.duckdb}"
export PAN_SOIL_CSV="${PAN_SOIL_CSV:-$ws/landis_region4/data_r4/plot_soil_r4.csv}"
export PAN_PNET_CLIMATE="${PAN_PNET_CLIMATE:-$ws/landis_region4/data_r4/pnet_climate_r4_l4.jld2}"
export PAN_CO2_CSV="${PAN_CO2_CSV:-$ws/landis_region4/data_r4/co2_annual.csv}"
export ECO_L4_SHP="${ECO_L4_SHP:-$ws/EPA_ECOREGIONS/us_eco_l4/us_eco_l4_no_st.shp}"
export PAN_SOLUS="${PAN_SOLUS:-$ws/SOLUS/SOLUS100}"
export PAN_THREADS="${PAN_THREADS:-4}"
export PAN_WORKER_LOG="${PAN_WORKER_LOG:-$here/cache/pan_worker.log}"
mkdir -p "$FE_CACHE_DIR"
for v in FIA_DB TREEMAP_DATA PAN_REPO PAN_FIA_DB FIA_TRAJ_DB PAN_CLIMATE_DB PAN_SOIL_CSV PAN_PNET_CLIMATE PAN_CO2_CSV ECO_L4_SHP PAN_SOLUS; do
    [[ -e "${!v}" ]] || echo "WARNING: $v=${!v} does not exist — set $v (or WORKSPACE) to its location" >&2
done

# --- Julia depots (trailing ':' keeps the bundled stdlibs reachable) ---------------------------------------------------
fe_depot="$here/.julia_depot_fe"
[[ -d "$ws/FVSclaude/.julia_depot" ]] && fe_depot="$fe_depot:$ws/FVSclaude/.julia_depot"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$fe_depot:}"
pan_depot="$PAN_REPO/.julia_depot_wt"
[[ -d "$ws/landis_parametrization_julia/.julia_depot" ]] && pan_depot="$pan_depot:$ws/landis_parametrization_julia/.julia_depot"
# a Pan checkout may ship its own depot layering (wt_env.sh exports JULIA_DEPOT_PATH): use it, so already-installed packages are reused
if [[ -f "$PAN_REPO/wt_env.sh" ]]; then
    wt_depot="$(bash -c "source '$PAN_REPO/wt_env.sh' >/dev/null 2>&1; printf %s \"\${JULIA_DEPOT_PATH:-}\"")"
    [[ -n "$wt_depot" ]] && pan_depot="$wt_depot"
fi
export PAN_DEPOT_PATH="${PAN_DEPOT_PATH:-$pan_depot:}"     # the Pan worker runs with this depot (src/PanEngine.jl)

# --- install packages when needed ----------------------------------------------------------------------------------------
setup() {   # setup <project dir> <depot path> <stamp file>
    local proj="$1" depot="$2" stamp="$3"
    if [[ "${SETUP:-auto}" == 1 || ( "${SETUP:-auto}" == auto && ( ! -f "$stamp" || "$proj/Manifest.toml" -nt "$stamp" ) ) ]]; then
        echo "Installing Julia packages for $proj (depot ${depot%%:*}) ..."
        JULIA_DEPOT_PATH="$depot" julia --project="$proj" -e 'using Pkg; Pkg.instantiate()'
        mkdir -p "$(dirname "$stamp")"; touch "$stamp"
    fi
}
setup "$here" "$JULIA_DEPOT_PATH" "${JULIA_DEPOT_PATH%%:*}/.fe_instantiated"
[[ -f "$PAN_REPO/Project.toml" ]] && setup "$PAN_REPO" "$PAN_DEPOT_PATH" "${PAN_DEPOT_PATH%%:*}/.pan_instantiated"

host="${HOST:-127.0.0.1}"; port="${PORT:-8088}"; threads="${JULIA_NUM_THREADS:-4}"
echo "Forest Explorer (Pan) http://$host:$port  workspace=$ws  depot=$JULIA_DEPOT_PATH  pan=$PAN_REPO (depot $PAN_DEPOT_PATH)"
exec julia --project="$here" --threads="$threads" -e "using ForestExplorer; ForestExplorer.start_server!(; host=\"$host\", port=$port)"
