# ACTIVE GOAL — Verify the LARGEST FIA divergences are irreducible (or fix them)

## Mission
Work down the non-bit-exact FIA bucket from the BIGGEST divergences. For each, prove it is either
(a) irreducible fp-ulp, (b) a cornered primitive (AVHT40/self-thin RDPSRT tie-break, tripling-order), or
(c) a known FVS-side bug where FVSjl is correct — OR find a real FVSjl bug and fix it. Do NOT declare
"irreducible" without a bulletproof measurement.

## DOCTRINE (hard-won this session)
1. **Per-record treelist comparison is INVALID after FVS tripling/splitting.** Tripled sub-records share
   (SpeciesFIA,PtIndex,TreeId); TreeIndex is re-assigned after mortality. Any keyed per-record diff mispairs
   sub-records → FAKE 10-650% "divergences." Trust ONLY: (a) stand-level `.sum` aggregates (tripling-invariant),
   or (b) per-tree DBH strictly in the PRE-split/pre-mortality window (both sides at full initial tree count).
2. **Confirm the trigger by MEASUREMENT, not inference.** Instrument / dump raw records; eyeball before concluding.
3. **Structure bit-exact + volume-only divergence** ⇒ check the volume path (r9clark/vollib), not growth.
4. **Big structural divergence** ⇒ diff_one stand-level: cycle-1 TopHt-only-then-cascade = AVHT40 tie-break
   (cornered); divergence from a self-thin cycle = self-thin RDPSRT (cornered); per-tree DBH divergence in the
   pre-tripling window = a REAL growth bug (dig + fix).
5. **Check input quality.** Extreme heights often trace to garbage SITE_INDEX (e.g. NE SI up to 581). GIGO
   faithfully reproduced by both FVS and FVSjl is not a bug.
6. Document every verdict; fixes validated bit-exact-or-cornered vs live + byte-identical normals.

## Scope (sweep DB data/fia_sweep.db)
Largest first: 101 stands struct>50%, 9,298 struct 5-50%, then the 2,140 vol>10% (r9clark already fixed).
Tooling: test/harness/fia/verify_treeid_magnitude.jl (+ no-mortality-window variant), diff_one.jl.

## Off-switch
`touch docs/FVS_LARGEST_DIV_COMPLETE` (USER's call). Log: docs/FVS_LARGEST_DIV_AUDIT.md.
