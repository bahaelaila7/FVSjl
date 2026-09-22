#!/usr/bin/env bash
# Launch the Forest Growth Explorer server.
#
# Paths are derived from this script's own location — nothing absolute is hard-coded.
#   - Julia depot: honors an existing $JULIA_DEPOT_PATH, else defaults to the
#     PERSISTENT depot beside the workspace. (~/.julia under $HOME is wiped on
#     restart; only the workspace volume persists, so packages must live there.)
#   - Data dir:   honors $TREEMAP_DATA, else ForestExplorer derives it (see
#     default_datadir); override with e.g. TREEMAP_DATA=/path ./serve.sh
#   - Bind:       HOST (default 127.0.0.1) and PORT (default 8080).
#set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -z "${JULIA_DEPOT_PATH:-}" ]; then
  # here = <ws>/FVSjl/apps/forest-explorer  ->  up 3 = <ws>
  export JULIA_DEPOT_PATH="$(cd "$here/../../.." && pwd)/.julia_depot"
fi

host="${HOST:-127.0.0.1}"
port="${PORT:-8087}"
threads="${JULIA_NUM_THREADS:-auto}"      # parallel FVSjl runs across the AOI's plots

echo "depot: $JULIA_DEPOT_PATH"
echo "serving http://$host:$port  (threads=$threads; HOST/PORT/JULIA_NUM_THREADS/TREEMAP_DATA to override)"
exec julia --project="$here" --threads="$threads" -e "using ForestExplorer; ForestExplorer.start_server!(; host=\"$host\", port=$port)"
