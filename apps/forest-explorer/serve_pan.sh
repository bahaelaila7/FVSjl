#!/usr/bin/env bash
# Launch the Pan-enabled Forest Explorer (branch forest-explorer-pan) without writing into shared locations:
#   - Julia depot: a private depot layered over the workspace's FVSjl depot (packages are read from there)
#   - FVSjl run caches go to $FE_CACHE_DIR (default: <this app>/cache), not into the TreeMap data dir
#   - Pan worker: PAN_REPO (default /workspace/landis_scen), started in the background on server start
# Usage: ./serve_pan.sh   (HOST/PORT/JULIA_NUM_THREADS/TREEMAP_DATA/PAN_REPO/PAN_THREADS to override)
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$here/.julia_depot_fe:/workspace/FVSclaude/.julia_depot:}"
export FIA_DB="${FIA_DB:-/workspace/FVSclaude/SQLite_FIADB_ENTIRE.db}"   # FVSjl master FIA DB (read-only)
export TREEMAP_DATA="${TREEMAP_DATA:-/workspace/FVSclaude/treemap}"
export FE_CACHE_DIR="${FE_CACHE_DIR:-$here/cache/sim_cache}"
export PAN_REPO="${PAN_REPO:-/workspace/landis_scen}"
export PAN_THREADS="${PAN_THREADS:-4}"
export PAN_WORKER_LOG="${PAN_WORKER_LOG:-$here/cache/pan_worker.log}"
mkdir -p "$FE_CACHE_DIR"
host="${HOST:-127.0.0.1}"; port="${PORT:-8088}"; threads="${JULIA_NUM_THREADS:-4}"
echo "Forest Explorer (Pan) http://$host:$port  depot=$JULIA_DEPOT_PATH  data=$TREEMAP_DATA  pan=$PAN_REPO"
exec julia --project="$here" --threads="$threads" -e "using ForestExplorer; ForestExplorer.start_server!(; host=\"$host\", port=$port)"
