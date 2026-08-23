# MECHPREP / BURNPREP (ESTAB site-preparation) — validation

## TOP-LINE: MECHPREP_CORNERED
The esetpr keyword parse + the WK6 per-plot IPPREP sampler are **transcribed and
proven BIT-EXACT at the kernel level** against the live oracle (FVSie_estabdump).
But the **live end-to-end** wiring cannot be validated bit-exact on the available
IE fixtures, for four independently-measured reasons (below). Staged, not committed
(doctrine: cannot validate live ⇒ do not commit). Gate untouched: src/ and test/
have zero changes ⇒ multicycle 339/11 byte-identical by construction.

Oracle: `/workspace/.iework/FVSie_estabdump` — FVSie_g16 relinked with an
instrumented `estab.f`/`esnutr.f` (dumps ESDRAW, WK6[50], SUMUP, IPPREP,
per-plot IPREP+PROB1+ITPP+NEWTPP, ESETPR parse, and the AUTAL XTES decision).
Fixture: under-stocked stand `753189105290487` (ie_understocked.db).

## What IS bit-exact (measured)

### 1. ESETPR keyword parse (esetpr.f) — BIT-EXACT
Oracle `ZSETPR` on the disturbance tally:
| scenario           | PMECH | PBURN | IALN(2) | IALN(3) |
|--------------------|-------|-------|---------|---------|
| baseline (no kw)   | 0.0   | 0.0   | 0       | 0       |
| MECHPREP 2029 100  | 1.0   | 0.0   | 1       | 0       |
| BURNPREP 2029 100  | 0.0   | 1.0   | 0       | 1       |
| MECHPREP 2029 50   | 0.5   | 0.0   | 1       | 0       |
`ie_esetpr` reproduces all four (test_mechprep.jl, 23/23).

### 2. Per-plot IPPREP assignment (estab.f:373-399 sample-without-replacement) — BIT-EXACT
Replaying the oracle's **exact** normalized SUMUP + WK6[50] through
`ie_esetpr_sample`:
- MECHPREP 100% → all 50 plots IPREP=2 (mech) — matches oracle.
- BURNPREP 100% → all 50 plots IPREP=3 (burn) — matches oracle.
- MECHPREP 50%  → WK6-sampled NONE/MECH mix `1212222211…` — matches oracle bit-for-bit.
- baseline (ESPREP mix 0.799/0.160/0.040) → `1212311111…` — matches oracle bit-for-bit.
The kernel transcription is faithful; given the same (SUMUP, WK6) it is exact.

## The measured obstacle (why LIVE validation is blocked)

### A. Site prep is exercised ONLY on the AUTAL disturbance tally
The estab.f ingrowth path (NTALLY==99 → reset to 1 → `GO TO 44`) **structurally
bypasses the ESETPR call** and forces `IALN(2)=1, PMECH=PBURN=0` ⇒ every plot
IPREP=1. Threshold-free proof: on the plain under-stocked stand (any of MECHPREP
2019, at inventory, or with a non-cutting thin), the oracle dumps `ZSETPR` = absent
and `IPPREP` = all-1 in every scenario, and PROB1 is identical (0.6174). Site prep
needs a **harvest** removing ≥THRES1(=0.10) of the stand → the AUTAL LAUTAL branch
sets NTALLY=1 and reaches ESETPR. Constructed with `THINBTA 2029 10.0` (removes
~99% ⇒ XTES=0.9877), which does drive the oracle onto the ESETPR path.

### B. PROB1 (the designated threshold-free signal) is PREP-INVARIANT here
On the disturbance tally, the oracle's per-plot `PROB1(NCOUNT)` = **0.1870904 for
IPREP = 1, 2, AND 3 alike** (and equal across baseline/MECH100/BURN100/MECH50).
ESTOCK's SPRE site-prep term is 0 for this stand's habitat series, so PROB1 does
not respond to prep. The prep effect lives only in the **species mix** (CPRE/FPRE
via espadv/espxcs) and **heights** (UPRE via essubh) — signals NOT in the
estab.f:1154 per-plot dump. So PROB1 cannot discriminate prep on this fixture; the
task's premise ("per-plot PROB1 should differ by prep") is measurably false here.
(This is bit-exact, not a discrepancy — PROB1 is simply prep-independent.)

### C. jl's disturbance-tally seed stream is DESYNCED from the oracle
On `h_base.key` (THINBTA 2029 + ESTAB):
```
oracle:  cyc1 ingrowth  ESDRAW=43303
         cyc2 DISTURB   ESDRAW=61997   (+ cyc3 continuation 61997)
jl:      cyc2 DISTURB   seed0 =43303   (+ cyc3 continuation 43303)
```
jl does **not** fire the cyc-1 ingrowth tally the oracle fires, so jl's disturbance
tally grabs seed 43303 (the oracle's cyc-1 seed) instead of 61997 — a **one-tally
offset**. The WK6 vector jl would sample IPPREP from therefore does NOT match the
oracle's, so the discriminating cases (baseline mix, MECHPREP 50%) cannot be
bit-exact live even though the kernel is exact given matched inputs. (MECHPREP/
BURNPREP 100% are WK6-independent, but see B/D — their effect is unvalidatable.)

### D. jl's disturbance establishment is itself cornered/garbage on this scenario
jl .sum on `h_base` produces a cyc-3 establishment explosion (2039: TPA 1371,
SDI 310748, QMD 2.4) vs the oracle's ~603 — the known cornered "disturbance-tally
seedling regime" the IE FIA fixtures mask, badly exposed by the near-clearcut.
The species-mix signal that site prep would modulate is buried under this.

### E. jl also lacks ESPREP (esprep.f) entirely
The no-keyword default prep probabilities (XPREP calibration, esblkd.f DATA) are
unported. Needed for the baseline mix; transcribed in `ie_esprep` but XPREP must be
sourced from the FVSie block data. Not on the critical path for the keyword cases.

## Verdict
Faithful, kernel-bit-exact transcription of esprep/esetpr + the WK6 sample-without-
replacement is staged (`esprep_esetpr.jl`) with an oracle-free regression
(`test_mechprep.jl` + `golden_mechprep.jl`, 23/23) and the proposed parse+wiring
patch (`kw_estab_mechprep.patch.jl`). Live wiring is deferred: it cannot be
validated end-to-end (A–E) and would ride the desynced/cornered disturbance stream,
risking regression rather than the required bit-exact-or-cornered proof. The
prerequisite for applying it is resolving the disturbance-tally seed-stream
divergence (C) — a pre-existing establishment-scheduling/timing issue, out of scope
for the site-prep keyword itself.

## Reproduce
```
# oracle (already relinked): dumps → /workspace/.iework/estabdump.txt
cd /workspace/FVSjl/scratchpad/estab/mechprep
echo h_mech100.key | /workspace/.iework/FVSie_estabdump   # ZSETPR/ZIPPREP/ZPROB
# kernel regression (oracle-free):
JULIA_DEPOT_PATH=/workspace/.julia_depot julia -e 'include("staged/test_mechprep.jl")'
# jl live desync evidence:
JULIA_DEPOT_PATH=/workspace/.julia_depot FVSJL_AUTOES_DEBUG=1 \
  julia --project=/workspace/FVSjl -e 'using FVSjl; FVSjl.run_keyfile("h_base.key"; variant=FVSjl.InlandEmpire())'
```
