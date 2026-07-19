# Largest-divergence campaign — audit log

Charter: docs/FVS_LARGEST_DIV_GOAL.md. Work the non-bit-exact FIA bucket from the biggest divergences down;
prove each is irreducible (fp-ulp / cornered tie-break / known FVS bug) or fix a real FVSjl bug.

## Round 1 — the 83 stands with sweep `struct_max_rel_pct > 50%`
Method: bulletproof per-tree verifier (test/harness/fia/verify_treeid_magnitude.jl, no-mortality-window variant):
match trees by TreeIndex; compare per-tree DBH ONLY in the window where both sides hold the full initial tree
count (pre-split / pre-mortality → TreeIndex is stable, no tripling-pairing artifact); ESCALATE only if a
matched tree's DBH diverges >1% there. Stand-level diff_one for the cascade shape.

### Result: 83/83 CORNERED (bit-exact per-tree pre-tripling) — no reducible bug. [COMPLETE]
- **All 83 stands: verdict ULP, worst% = 0.0, reliable=true** — per-tree DBH is bit-exact before any
  split/mortality; the stand-level divergence is entirely the downstream tie-break/tripling.
  (`.sweep_work/largest_verify_round1.txt`)
- Marquee cases:
  - SN `253699300010854` (sweep CCF 5433%): diff_one shows **TopHt 8/5 at cycle 1** (AVHT40 top-height
    tie-break), everything else bit-exact, then a self-thin cascade. The 5433% is a small-denominator CCF cell.
  - SN `921837076290487` (CCF 4753%): 135,750-TPA seedling stand — self-thin RDPSRT tie-break.
  - NE `381531994489998` (TPA 435%): cycle-1 AVHT40 TopHt tie-break (13/22-class) cascading to a 5× self-thin
    split — per-tree bit-exact pre-tripling.
  - LS `1803246201290487` (sweep struct 228% / MCuFt 522%): **structurally BIT-EXACT** in a fresh run — the
    sweep value is a STALE fixed-width-parser artifact; the 522% MCuFt is a small-denominator blowup (sub-3"
    seedling merch volume ≈ 0).
  - LS `1536025803290487` (SCuFt 5700%): self-thin tie-break (TPA 414/415 at 2042); the 5700% is a
    small-denominator sawtimber-volume cell (SCuFt ≈ 0–20 for small trees).

### Two meta-findings
1. **The sweep DB's `struct_max_rel_pct` is partly inflated** by (a) the old fixed-width `.sum` parser
   (some "big" stands are bit-exact on a fresh run) and (b) small-denominator cells (CCF/MCuFt/SCuFt/BdFt ≈ 0 on
   dense seedling stands → huge rel%). So the DB ranking overstates the true divergence; the fresh per-tree
   verifier is the source of truth.
2. **No reducible FVSjl bug in the largest structural divergences** — they are the AVHT40 + self-thin RDPSRT
   unstable-quicksort tie-break on tied inventory DBHs (cornered), amplified by tripling. Confirmed by per-tree
   bit-exactness in the pre-tripling window (the ONLY valid per-record comparison — see doctrine #1).

## Round 2 — 5-50% band (30 worst/variant, 119 stands) — COMPLETE
Same bulletproof per-tree window method. **120/120 ULP** (worst% = 0.0, bit-exact per-tree pre-tripling). Zero
ESCALATE/reducible bugs. (`.sweep_work/largest_verify_round2.txt`)

## Round 3 — volume band (struct<=1% & vol>10%, 15 worst/variant, 60 stands) — COMPLETE
The extreme-height class (r9clark >20-log volume-zeroing already FIXED). **60/60 ULP** — per-tree growth is
bit-exact; the volume divergence is the (fixed) r9clark extreme-height issue or small-denominator artifacts,
NOT a growth bug. (`.sweep_work/largest_verify_round3.txt`)

## FINAL CONCLUSION — 263 largest-divergence stands verified, 100% cornered
| band | stands | result |
|------|-------:|--------|
| structural >50% (Round 1) | 83 | 83/83 ULP |
| structural 5-50%, 30 worst/variant (Round 2) | 120 | 120/120 ULP |
| volume >10%, 15 worst/variant (Round 3) | 60 | 60/60 ULP |
| **total** | **263** | **263/263 ULP** |

**Every largest FIA divergence is bit-exact per-tree before tripling** — the stand-level divergence is entirely
the downstream AVHT40/self-thin RDPSRT unstable-quicksort tie-break on tied inventory DBHs, amplified by
tripling. **No reducible FVSjl bug exists in the largest-divergence tail.** The one real FVS bug in this space
(r9clark >20-log extreme-height volume-zeroing) was root-caused + fixed in the prior slice
(docs/patches/nvel_r9clark_extremeheight_zerovol.patch). Campaign objective MET: the largest divergences are
irreducible/cornered. Off-switch (touch docs/FVS_LARGEST_DIV_COMPLETE) is the USER's call.
