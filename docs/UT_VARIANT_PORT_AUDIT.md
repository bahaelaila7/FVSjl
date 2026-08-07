
## 2026-08-06 — #147 VERIFIED FIXED (woodland mortality bit-exact); SEPARATE woodland-DG lead surfaced
Bounded verification of the #147 "UT dense woodland-seedling mortality OVER-KILL ~4×" (tracker still shows pending;
dense-seedling campaign memory says fixed). ut_sub.db has no Gambel-oak(814) but holds the SAME woodland mortality
class — dense pinyon-juniper. Densest stand 434206325489998 (11695 TPA of tiny woodland stems, BA 54), aligned
NUMCYCLE 10, jl vs live FVSut:
  TPA **BIT-EXACT all 9 cycles** (11695→10068 = live exactly; live gently self-thins −14%). ⇒ the ~4× over-kill is
  GONE — #147 mortality is FIXED/verified. Tracker "pending" is STALE.
★ SEPARATE observation (NOT #147, which is mortality): jl BA runs LOW and compounds — 54→**49** vs live 54→**71** by
2095 (jl −31% BA; QMD jl 0.94" vs live 1.14"). TPA identical ⇒ this is a woodland small-tree DIAMETER-GROWTH
under-growth, not mortality. It is a more-severe manifestation of the known UT "multi-cycle BA tail ~−7% by 2090"
(utt01) on a VERY dense PJ stand (11695 TPA). Candidate causes (next session, if pursued): woodland/PJ small-tree DG
(regent/SMDGF) under-prediction on dense sub-1" stems, or a #150/#153-crown-init residual specific to very-dense
woodland. Bounded lead, DISTINCT from #147 (mortality, resolved) — worth a dedicated woodland-DG measurement, low
priority (woodland BA, not mortality; and PJ volume is the DVEW tail already noted). No jl change this session; the UT
port is untouched, so this BA-DG tail is pre-existing.

## SDI-gate tem 35000-cap bug — FIXED (2026-08-07, commit 7ce8f1f, cluster w/ EM 04b15e6)
UT self-thin SDI-in-effect gate `tem` omitted the min(·,35000) cap (ut/morts.f). For UT dq10 IS the
Zeide DR10, so `tem = t55d10` fixes both the missing cap and the correct Zeide d10 basis. On ultra-dense sub-1"
cohorts the uncapped tem ≫ tt → wrong background fallthrough → self-thin under-kill. VALIDATED no-regression:
utt01 TPA bit-exact vs live oracle (fix inert for QMD>~0.9").
