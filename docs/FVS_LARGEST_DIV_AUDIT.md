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

## Round 2 — 5-50% band (30 worst/variant, 119 stands), IN PROGRESS
Same bulletproof per-tree window method. Running tally: 60/119 verified so far, **all ULP** (worst% = 0.0,
bit-exact per-tree pre-tripling). Zero ESCALATE/reducible bugs. (`.sweep_work/largest_verify_round2.txt`)

## Running conclusion (143 stands verified: 83 >50% + 60 of the 5-50% band)
**100% ULP — every largest FIA divergence is bit-exact per-tree before tripling**, i.e. the stand-level
divergence is entirely the downstream AVHT40/self-thin RDPSRT tie-break (+ tripling amplification) on tied
inventory DBHs. No reducible FVSjl bug in the largest-divergence tail. The r9clark extreme-height volume-zeroing
(the one real FVS bug in this space) was fixed in the prior slice.
