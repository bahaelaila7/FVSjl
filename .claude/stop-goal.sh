#!/usr/bin/env bash
# This session's Stop hook. The active campaign is the CS (Central States) variant port;
# delegate to stop-cs.sh (reads docs/CS_GOAL.md + doctrine; off-switch: touch docs/CS_COMPLETE).
# (NE port is complete — docs/NE_COMPLETE exists, so stop-ne.sh would no-op anyway.)
exec /workspace/FVSjl/.claude/stop-cs.sh "$@"
