# Tolerance & Semantic-Coverage Audit (SN)

## Tolerance taxonomy — every test tolerance is now one of:

1. **Float-ULP / print-resolution** (bit-exact). The `.sum`/`.tre`/carbon-report oracle is a printed
   integer or fixed-decimal column; the tightest possible match is the print granularity:
   - integer columns → `di(x)=trunc(x+0.5)` equality, or `abs(Δ) <= 1` (a sub-ULP float diff can only flip
     the truncation boundary by one — this IS the float-ULP equivalent for a truncated column);
   - F7.1 carbon columns → `atol = 0.05`.
   Examples now tightened to this: snt01 cycle-1 (di==baseline), all carbon_snt LIVE pools + Floor + Shrub
   + Total at inventory (0.05).

2. **Accumulated-ULP "tails"** on multi-cycle VOLUME SUMS (e.g. test_multistand_sum `±8 cuft`): Float32 vs
   Fortran REAL*4 rounding compounding over many trees × cycles. Floating-point in origin (not a semantic
   gap), documented in-line as "ulp tail".

3. **Invariant / range checks** (`0 ≤ p ≤ 1`, `0 ≤ CBD ≤ 0.35`, physical bounds) — not tolerances.

4. **@test_broken** — KNOWN non-bit-exact residuals, tracked honestly instead of hidden behind a passing
   loose tolerance. Current set (10):
   - carbon_snt dead pools BelowD/StandD/DDW (~0.5/0.7/1.2): the crown-lift one-cycle timing lag (applied
     in the next cycle's fuel loop vs FVS same-cycle; inventory + final cycle ARE bit-exact).
   - carbon_jenkins post-mortality DDW (same gap), and the input-snag bole 3.92-vs-3.8 height-dub residual,
     and the age-aware Stand-Dead tracking residual.
   These flip to a (visible) failure the moment the underlying pool reconciles bit-exact.

   NB carbon_jenkins is a NON-bit-exact-GROWTH fixture (synthetic LP diameter-growth tail); its live-pool
   agreement is the growth tail, NOT a carbon property — bit-exact carbon is owned by carbon_snt.

## Coverage audit — SN semantics

A full subsystem sweep (growth/mortality/density/volume/regen/sprout/event-monitor/structure/ECON/IO/FFE)
confirmed bit-exact coverage of the natural-process core and the FFE extension (9/9 carbon+fire DBS tables,
both CARBCALC methods, crown-lift). The audit surfaced ONE active silently-ignored gap and one no-op:

- **YARDLOSS** (cuts.f:1387) — ACTIVE in snt01.key / sn.key but was silently ignored (unrecognized keyword
  fell through to a no-op else). NOW PORTED: reported removed merch/saw/board scaled by (1−PRLOST); total
  cubic + TPA stay at full physical removal (verified cut_specpref 272→136 @ PRLOST=0.5). The loss-routing
  to FFE fuel pools (DSNG/SSNG/CTCRWN) is the remaining C7 coupling. Inert on snt01/sn (their thins remove
  only sub-merch trees), which is why the gap passed unnoticed.
- **SALVAGE** — ABANDONED in the Fortran itself (cuts.f:103); now a RECOGNIZED no-op (not a silent fall-through).

Verified OUT-OF-SCOPE for SN (not gaps): the 8 insect/disease models are Western-only (dfb/lpmpb/wsbwe);
FMSNGHT snag height-loss is a no-op in SN (HTX=0); MORTMSB/MATUREW extra-mortality is keyword-gated and
default-inert (QMDMSB=999); COMPRESS PCA clustering is a recognized keyword (the eigensolver partition is
not bit-identical, documented). The remaining unported items are output-only (C6 report detail) or require
a fuller (non-stripped) DBS validation binary.

KNOWN systemic note: unrecognized keywords currently fall through to a silent no-op (keyword_dispatch.jl
final else). YARDLOSS shows the risk; surfacing unrecognized keywords (collect + warn) would convert future
silent gaps into visible ones — the recommended next safeguard.
