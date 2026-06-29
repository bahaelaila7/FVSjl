#!/usr/bin/env bash
# This session loaded THIS path as its Stop hook at startup. The active campaign is the NE port;
# delegate to stop-ne.sh (reads docs/NE_GOAL.md + doctrine; off-switch: touch docs/NE_COMPLETE).
exec /workspace/FVSjl/.claude/stop-ne.sh "$@"
