#!/usr/bin/env bash
# coverage_check.sh — COVERAGE GATE wrapper (see coverage_gate.jl).
#   coverage_check.sh            check against the committed baseline (fails on lost coverage / unlisted zero-coverage)
#   coverage_check.sh baseline   (re)record the baseline + COVERAGE_GAPS.tsv + docs/TIERED_COVERAGE_<date>.md
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
cd "$ROOT" && exec julia --project=. test/harness/tiered/coverage_gate.jl "${1:-check}"
