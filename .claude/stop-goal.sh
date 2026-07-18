#!/usr/bin/env bash
# This session's Stop hook dispatcher. The FIA/FVS behaviour-compatibility campaign is COMPLETE
# (docs/FIA_FVS_COMPAT_COMPLETE). The active campaign is now: FIX LIVE-FVS CRASHES on real FIA inventory
# with minimal semantically-plausible upstream source patches (goal docs/FVS_LIVECRASH_FIX_GOAL.md;
# off-switch: touch docs/FVS_LIVECRASH_COMPLETE). (To resume FIA: exec .claude/stop-fia.sh)
exec /workspace/FVSjl/.claude/stop-livecrash.sh "$@"
