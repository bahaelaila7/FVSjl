# FVS live-crash fix — working audit log

Campaign: root-cause + minimally patch every live-FVS SIGFPE on real FIA inventory. Goal/doctrine:
docs/FVS_LIVECRASH_FIX_GOAL.md. Stands: docs/fvs_livecrash_stands.txt (12, 6 sites). FVSjl runs all clean.

## Scope (measured 2026-07-18, current tmp/oracles binaries)
| site | variant | stands | status |
|------|---------|-------:|--------|
| cs/varmrt.f:162 | CS | 1 | root-caused (TEMKIL/TEMSUM div0); fix pending |
| cs/grincr.f:449 | CS | 3 | TBD |
| cs/grincr.f:437 | CS | 1 | TBD |
| cs/fvs.f:197    | CS | 1 | TBD |
| ls/dgdriv.f:134 | LS | 3 | TBD |
| ls/dgdriv.f:353 | LS | 2 | TBD |
| ls/htdbh.f:336  | LS | 1 | TBD |
