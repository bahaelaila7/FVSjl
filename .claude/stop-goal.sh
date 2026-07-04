#!/usr/bin/env bash
# This session's Stop hook. SN/NE/CS ports + the non-ULP divergence campaign are complete
# (docs/DIVERGENCE_COMPLETE). The active campaign is now the LS (Lake States) variant port.
# Delegate to stop-ls.sh (reads docs/LS_GOAL.md + doctrine; off-switch: touch docs/LS_COMPLETE).
exec /workspace/FVSjl/.claude/stop-ls.sh "$@"
