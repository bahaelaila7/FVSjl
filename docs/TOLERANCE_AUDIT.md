# TOLERANCE AUDIT — the closure checklist

Every non-ULP-proven tolerance in the asserted test suite (integration/ + unit/). Each row must end
at **BIT-EXACT** (op fixed, `==`) or **ULP** (op cornered + width proven). Status: ⬜ open · ✅ closed.

Baseline: suite 7658 pass / 2 broken. Doctrine + approach: docs/TOLERANCE_GOAL.md.

## Legend for "kind"
- `%` = relative/percentage bound · `absN` = multi-unit absolute · `±1` = single-unit per-column
- `print` = half-print-width vs a rounded oracle field (may already be provable-exact via `==` on rendered)

---

## C1 — percentage / relative bounds (explicitly empirical)
| ⬜ | file:line | bound | compares | plan |
|----|-----------|-------|----------|------|
| ⬜ | test_allspecies.jl:44-49,72 | max(abs, 1.2-3.0%) | grown .sum cols (4 variants) | re-measure post-QMDGE5; trace residual→op per column |
| ⬜ | test_timeint.jl:47,48 | 0.03·ft+1 | BA, cuft | non-native-cycle DGSCOR? corner/fix |
| ⬜ | test_carbon.jl:50-53 | 0.005·v+0.1 | grown carbon pools | trace |
| ⬜ | test_multicycle.jl:57,59,60,61 | rtol=0.002, absN | TCuFt/BA/SDI/QMD | trace |

## C2 — Scribner board-foot "Float32 noise" (1+pct·x)
| ⬜ | files | bound | plan |
|----|-------|-------|------|
| ⬜ | test_spleave:43, test_voleqnum:50, test_minharv:53, test_sprout_regen:49, test_tfixarea:43, test_cuteff:45, test_bfvolume:47, test_tcondmlt:47, test_fertiliz:47 | 1+0.002·x / 1+0.005·x | BdFt | corner Scribner round op; prefer `==` on rendered BdFt |
| ⬜ | test_dbs_summary:58, test_readcor:57, test_net01:188 | absN | BdFt Scribner | same |

## C2b — multi-unit absolute bounds
| ⬜ | file:line | bound | compares | plan |
|----|-----------|-------|----------|------|
| ⬜ | test_cst01.jl:177-182 | TPA4/SDI4/CCF10/TopHt2/QMD0.2 | grown | re-measure post-QMDGE5 |
| ⬜ | test_cst01.jl:80,81,116,122,123,124 | ≤1, TPA3/CCF3/TopHt2/QMD0.15 | cyc0-vol/2002 | re-measure |
| ⬜ | test_growth.jl:169,170,171 | ≤4 | TCuFt/MCuFt | trace |
| ⬜ | test_cycleat.jl:63,64,65 | TPA8/BA3/SDI6 | grown | trace |
| ⬜ | test_multistand_sum.jl:44, test_multistand.jl:68 | ≤8, ≤3 | cuft | trace |
| ⬜ | test_treeszcp.jl:52,59,85 | ≤5, ≤4 | TPA/TopHt | trace |
| ⬜ | test_dbs_treelist.jl:44 | ≤3 | Σ(TCuFt·TPA) | trace |
| ⬜ | test_mortmsb.jl:64, test_sprout_table.jl:46 | ≤2 | all cols | trace |
| ⬜ | test_resetage.jl:54, test_structure_stage.jl:70,81 | 0.2/0.55/0.1 | MAI/strdbh/DBHNOM | trace |
| ⬜ | test_net01.jl:187,219,311,322,323 | atol 4/2/2/6/2 | TCuFt/BA/TPA | trace |
| ⬜ | test_event_monitor.jl:51,69 | 3, 0.2 | tpa/BSDI | trace |
| ⬜ | test_estab_rng_d10.jl:70,74 | 0.02/0.01 | D10 | trace |
| ⬜ | test_regen_coverage.jl:71 | atol=2 | TPA/space | trace |

## C3 — FFE fire / carbon loosened atols
| ⬜ | file:line | bound | compares | plan |
|----|-----------|-------|----------|------|
| ⬜ | test_carbon.jl:301,302,675-677,702,703,737,777-794 | 0.1-1.0 | FFE AGL/Merch/snag/DDW/Released | trace fire kill-dist + FMCROWE |
| ⬜ | test_lst01_ffe.jl:57,58,72,109 | 0.15/0.5/1/0.6 | flame/scorch/BA/carbon | corner van-Wagner transcendental; fire kill |
| ⬜ | test_fire.jl:38,49,71,92,97,98,119,121,122,143,145,146 | ≤2/3/4 | fire BA/TPA | trace FMEFF kill distribution |

## C4 — ~69 `±1/±2` per-column bounds (labeled ULP, not cornered)
| ⬜ | files | plan |
|----|-------|------|
| ⬜ | test_crnmult:30, test_fixmort:44, test_tripling:34, test_spgroup:34, test_htgstp:41, test_fix_scalers:35, test_hcor_calib:34,45, test_treeszcp:42,53, test_minharv:51, test_cuteff:43, test_bfvolume:45, test_voleqnum:47, test_tcondmlt:45, test_fertiliz:45, test_spleave:41, test_tfixarea:41, test_volume_override:52-54, test_pertree_defect:54-56, test_mcdefect:45, test_setsite:54, test_compute:57, test_estab_pccf:41,56,57, test_multistand:66,67, test_multistand_sum:42,43, test_bamax:66, test_dbs_cutlist:68-70, test_dbs_summary:57, test_growth:169, test_fire:38,71,92,119,143, test_net01:577-579,605, test_structure_stage:69 | for each: prove specific print/sum-order ULP → compare rendered `==`; OR fix op |

## borderline-ULP (verify the traced root holds, else move up)
| ⬜ | test_fire.jl:180,181 | atol 0.005/0.03 van-Wagner ^(7/6)/√ Float32 — confirm irreducible |
| ⬜ | test_net01.jl:41,213,270,363 | print-resolution atols — convert to rendered `==` where possible |

## already ULP-justified (print half-width) — CONFIRM then leave, or convert to rendered `==`
- test_carbon.jl `≤0.05` (carbon 1-dec print), QMD `≤0.05` across ~9 files, test_snt01:29,37, test_init:50.
  Verdict pending: these are print-half-width vs a rounded field. Preferred close = compare the
  RENDERED value with `==` (removes the atol entirely). Mark ULP-with-root only if rendering is infeasible.

## @test_broken (must carry a documented irreducible root)
| ⬜ | test_nohtdreg.jl:87 | WK3/DGSCOR sp33/65 serial-corr tail — re-verify irreducible |
| ⬜ | test_keyword_coverage.jl:160 | COMPRESS s22 eigensolver + s32 R8-VOLUME leak — re-verify |
| ⬜ | test_keyword_coverage.jl:181 | dormant (empty broken-set) — confirm |
