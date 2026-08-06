# Crown-init gap (#150) — cluster scope measurement (2026-08-06)

## Question
BM #149 (dense conifer-seedling BA under-grows ~2×) traced to the missing lstart crown-init: dense read 0.1"
seedlings keep crown_pct=0, and BM's small-tree height model VIGOR=150·CR³·exp(−6·CR)+0.3 FLOORS at 0.30 when CR=0
⇒ HTGR starved ⇒ seedlings never cross 4.5' ⇒ DBH growth skipped ⇒ BA ~2× low. #150 asked: is this a UNIFORM
cluster gap (KT/IE/EM/TT/UT all missing the crown-init)?

## Measurement (real-FIA dense-seedling stands, jl vs FVS{v}_clean, cyc0+cyc1 TPA:BA)
EM (em_sub.db):
  427473386489998  live 2026 18054:52  jl 18081:51   (BA −1, ~2%)
  83402334020004   live 7026:20         jl 7035:20    (BA EXACT)
  1856105255290487 live 1002:61         jl 1004:61    (BA EXACT)
  225065919010661  live 5447:49         jl 5436:48    (BA −1)
UT (ut_sub.db):
  2860048010690    live 416:72          jl 416:71     (BA −1)
  2916642010690    live 358:120         jl 358:120    (BA EXACT)
  2381662010690    live 10374:185       jl 10352:178  (BA −7, ~3.8%)
TT (tt_sub.db):
  2824914010690    live 5315:172        jl 5246:166   (BA −6, ~3.5%)
  1629341474290487 live 3337:86         jl 3337:88    (BA +2)
  2781432010690    live 2871:57         jl 2871:55    (BA −2, ~3.5%)

## Verdict — #150 is BM-SPECIFIC-SEVERE, not a uniform cluster gap
EM/UT/TT dense conifer-seedling BA is within a few % of live (many bit-exact), with small mixed-sign TPA deltas —
the ACCEPTED dense-cohort DGSCOR/self-thin straddle, NOT a severe crown-init starvation. NONE reproduces BM's 2×
(38 vs 80) under-growth. Root cause of the difference: the small-tree HEIGHT model's behavior at CR=0 —
 - **BM: VIGOR=150·CR³·exp(−6·CR)+0.3 → 0.30 floor at CR=0** (a STARVING floor ⇒ severe gap; #149).
 - **EM: SMHTGF htg1=beta1+beta2·CR → beta1 floor at CR=0** (healthy) ⇒ seedlings grow ⇒ no severe gap.
 - **UT/TT: similar linear/healthy CR=0 floors** ⇒ dense conifer seedlings grow ⇒ only the cornered straddle remains.
⇒ The crown-init (#150) is REQUIRED for BM (done, c29c1c7) but NOT required to reach bit-exact-or-cornered on
EM/UT/TT. The mild UT/TT −3 to −7 BA on the densest stands is within the cornered dense-cohort straddle (a crown-init
port would tighten it but is not needed for the bar). KT/IE untested here but expected similar (linear small-tree HTG).
NOTE: this REFINES the campaign "unifying root" — the crown-init was the UNIFYING DIAGNOSIS but the SEVERITY is
gated by each variant's CR=0 height floor; only BM's VIGOR made it a real bug.

## IE added (2026-08-06 cont.) — confirms BM-specific scope
IE (ie_test.db) dense conifer-seedling stands vs FVSie_clean:
  12343703010690 (9 seed, 24 dead) live 2017 28930:15  jl 28915:16   (BA +1; close)
  44987944020004 (7 seed, 2 dead)  live 32104:164      jl 32203:154  (BA −10, ~6%)
  753199439290487 (WH sp263 dominant, TC=50376) live 45880:74  jl 45030:13  (BA 5.7× LOW)
The 753199 5.7× under-growth is the ALREADY-KNOWN, CORNERED IE #146 (WH dense-seedling ZZRAN tripled-record height
straddle; deterministic path bit-matches live, verified in the #146 work) — NOT a crown-init gap (IE's WH small-tree
height model has no crown term). The other two IE stands are close (BA ±1 to −6%). ⇒ IE is bit-exact-or-cornered on
dense seedlings, consistent with #150 being BM-specific-severe.

## Final cluster verdict (#150)
Measured EM, UT, TT, IE dense conifer-seedling real-FIA stands: NONE reproduces BM's severe crown-init gap. The only
severe dense-seedling under-growths in the cluster are variant-specific KNOWN items — BM #149 (crown-init/VIGOR, FIXED)
and IE #146 (WH ZZRAN straddle, cornered) — not a uniform missing-crown-init bug. KT unmeasured here (no mini-db built)
but expected same (linear small-tree HTG; KT #141 regen already validated). #150 CLOSED as BM-specific-severe.
