#!/usr/bin/env bash
# This session's Stop hook dispatcher. SN/NE/CS/LS ports are complete; the tolerance-closure campaign is
# DONE (docs/TOLERANCE_COMPLETE) and the MODERNIZATION campaign is DONE (docs/MODERNIZATION_COMPLETE,
# user-authorized 2026-07-08). The active campaign is now FIA/FVS BEHAVIOUR COMPATIBILITY on the 4 variants:
# prove FVSjl reproduces live FVS's behaviour on real FIA inventory at scale, over full multi-cycle
# projections and under management — without regressing the closed campaigns' bit-exact-or-cornered floor.
# Delegate to stop-fia.sh (goal docs/FIA_FVS_COMPAT_GOAL.md; off-switch: touch docs/FIA_FVS_COMPAT_COMPLETE).
exec /workspace/FVSjl/.claude/stop-fia.sh "$@"
