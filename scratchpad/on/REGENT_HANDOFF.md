# ON sub-12cm REGENT small-tree growth — NEXT chunk (validation stand READY)

## Status
ON cyc0 .sum FULLY byte-identical vs FVSon_g16 (HEAD 1bdad71b). The last un-exercised growth path
is the REGENT sub-12cm branch (canada/on/regent.f, 632 ln) — currently `small_tree_growth!(::Ontario)`
errors loudly on any DBH < XMAX=4.72" (12 cm). ont01 is all-large so it's a validated no-op there.

## Validation stand (BUILT + oracle output captured)
- `scratchpad/on/ont_sm.{key,tre}` — 8 small trees (DBH 4-11 cm, all < 12 cm ⇒ every record hits REGENT).
- Oracle: `printf "ont_sm.key\nont_sm.tre\n" | /workspace/.onwork/FVSon_g16` (in a run dir) → `ont_sm.sum`:
  cyc0 = TPA 296519 / SDI 44964 / CCF **** / TopHt 7 / QMD 7.8 / FORTYP 122 size 3 stock 1.
  (cyc1 = the degenerate metric year-0/TPA-24 IOSUM artifact, cornered, same as ont01.)
- jl currently ERRORS on record 1 (DBH 1.9685" < 4.72") — confirms the branch is exercised.

## Port scope (regent.f sub-12cm branch, D < XMX)
Dependencies to port/verify: HTCALC (NC128 height-age curve — partially present via findag/volume age dub),
BALMOD (BAL modifier), HTDBH (Wykoff HT→DBH dubbing, htdbh.f), DGBND (DG size-cap), MULTS (XRHMLT/XRDMLT
RNG multipliers). The branch blends the NC128 height-age increment with the large-tree HTG over [XMIN,XMAX]
then dubs DBH via Wykoff HTDBH and caps DG via DGBND. Recipe: instrument FVSon_g16 regent.f per-tree dump
(D/H/HK/DKK/DG/HTG hex, IWHO gate) → single-.o swap into g16obj → relink FVSon_regdump → dump-replay
bit-exact per tree, THEN wire small_tree_growth!(::Ontario) + end-to-end .sum on ont_sm (cyc0 all cols).
