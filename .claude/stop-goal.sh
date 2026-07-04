#!/usr/bin/env bash
# This session's Stop hook. SN/NE/CS/LS ports are complete (tags through FVSsn+ne+cs+ls-done).
# The active campaign is now test-tolerance closure: drive every test tolerance to BIT-EXACT or
# PROVEN-IRREDUCIBLE-ULP. Delegate to stop-tolerance.sh (goal docs/TOLERANCE_GOAL.md + doctrine;
# off-switch: touch docs/TOLERANCE_COMPLETE).
exec /workspace/FVSjl/.claude/stop-tolerance.sh "$@"
