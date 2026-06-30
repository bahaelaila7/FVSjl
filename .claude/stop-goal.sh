#!/usr/bin/env bash
# This session's Stop hook. The CS port is complete (docs/CS_COMPLETE); the active campaign is now the
# non-ULP divergence fix campaign. Delegate to stop-divergence.sh (reads docs/DIVERGENCE_GOAL.md +
# doctrine; off-switch: touch docs/DIVERGENCE_COMPLETE).
exec /workspace/FVSjl/.claude/stop-divergence.sh "$@"
