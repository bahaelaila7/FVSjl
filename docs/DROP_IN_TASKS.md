# FVSjl Southern Drop-In — Task Tracker

Persistent progress tracker for the SN drop-in-replacement effort (survives container
restarts). Goal: **FVSjl (Southern) is a bit-exact drop-in for live FVSsn; the only
accepted divergences are ULP floating-point and the COMPRESS eigensolver.**

Canonical test runner: `julia --project=. test/runtests.jl`. Current suite:
**4467 pass + 18 broken** (broken = tracked, documented divergences, NOT failures).
Live-Fortran oracle: `/tmp/FVSsn_new` (rebuilt via `tests/fortran_baseline.sh`).

Keyword-coverage harness: `test/integration/test_keyword_coverage.jl` over
`test/keyword_coverage/scenarios/*` (37 scenarios). See
[memory: fvsjl-keyword-coverage-state] for the gate design.

---

## DONE (committed)
- [x] Wire the 37-scenario keyword-coverage harness into the main suite.
- [x] 10-yr-cycle **mortality** root-caused + fixed (linear FINT-extrapolated SDI trajectory, morts.f:225).
- [x] **s37_thinauto** MAI terminal-row quirk (evtstv.f:414/disply.f:392) — bit-exact.
- [x] Column-aware gate: structural cols strict ULP; merch/board cols 0.3% (Scribner/Behre FP) → **s14/s34** conform.
- [x] Structured-yaml **named schemas** for THINHT/THINRDEN/THINAUTO/VOLUME/BFVOLUME/MCDEFECT/BFDEFECT/SERLCORR/RESETAGE.
- [x] Robust yaml writer: lossless **raw passthrough** on parse failure / unknown keyword → yaml==key 36/37.

## OPEN — engine divergences (must reach ULP/compression-only)
- [ ] **s5/s9** — 10-yr board-foot tail (~4%): calibrated-species HTCALC/AGET + tripling-variance DG/HTG precision under a non-5 period. (TPA already within 3.) [memory: fvsjl-10yr-cycle-mortality]
- [ ] **s26_estab** — CORRECTED (the earlier "245.4 establishment-reduction" was a `.trl` column-parse ERROR — disregard it). Via a targeted GRADD-end LP-sum print in a recompiled FVS, the real per-cycle LP TPA is: FVSjl 2005=270.0/2010=63.31/2015=7.27 ; FVS 2005=270.0/2010=60.31/2015=8.75. **ESTABLISHMENT IS BIT-EXACT** (270 both at 2005). The divergence is the DENSE LP-cohort SELF-THINNING MORTALITY 2005→2010 (FVS 270→60.31, FVSjl 270→63.31, ~5% under-kill). LP are dbh 2.6→3.2 (near dbh_zeide≈3 SDI threshold) ⇒ likely small-tree SDI-inclusion or VARMRT distribution for a dense sub-merch cohort. NOT establishment, NOT the linear-G fix (5-yr=identity). TOOLING WIN: stripped FVS global DEBUG segfaults (fvsvol.f:530/fvs.f:376 uninit WRITEs), but targeted per-cycle prints (recompile one .f at -O0, relink) WORK.
- [ ] **s32_volume** — VOLUME zeroes SCF (cols past 80) → all trees prod=01. Per-tree .trl: FVS scuft=0 below dbh 10, mcuft=0 below ~6; FVSjl leaks small-tree merch (+19 each). mrules.f fixes DON'T apply (sawDib=6 overshoots; merchL=10 BREAKS 3 bit-exact defaults) → SN R8 uses merchL=8/sawDib=7-9; multi-threshold. Default path bit-exact.
- [ ] **s22_compress** — COMPRESS different eigensolver — **ACCEPTED per spec** (keep as documented-broken).

## OPEN — yaml / conversion (user-requested)
- [x] **Conversion tool operational (drop-in paths)** — `test/integration/test_translate.jl`: .key→.yaml→.key engine-equal + .tre→.csv preserves TreeRecords (the engine reads either form). Regression test in suite.
- [ ] **csv→tre legacy re-emission bug** — `write_tree_file` drops an F-field that shares a column with a packed nI1 field in the SN T-specifier layout (e.g. field `T60,F3.1` value 5.0 → 0, overlapping `T54,7I1`). Tracked @test_broken (10). Does NOT affect the drop-in (engine reads the .csv); fix the overlap resolution in write_tree_file for a byte-faithful legacy export.
- [ ] **s20_spgroup** — 2-record SPGROUP round-trip (group name + species-list record) — handle in the hierarchical redesign.
- [ ] **Task 8 — hierarchical, order-aware YAML redesign.** ⚠️ ORDER MATTERS: the current flat form is a `list` so order is preserved. The redesign must NOT become fully order-independent — some stand operations are genuinely order-dependent and must keep their relative order:
      • a SPGROUP / species-group definition must precede a THIN that references the group;
      • multiple same-cycle thinnings/operations apply in sequence;
      • COMPUTE / event-monitor variable definitions vs. their use.
      Make date-scheduled keywords order-free where safe, but preserve declared order within a stand for operations whose effect depends on it (model as ordered sub-lists / explicit sequence, not an unordered map).
- [ ] **Ensure the YML/CSV ↔ .key/.tre conversion tool is operational** (`bin/fvsjl-translate.jl`): verify both directions round-trip on all coverage scenarios + the examples; add a regression test.
- [ ] **Keyword documentation**: update `docs/KEYWORDS.md` to document each keyword in BOTH yml and .key form, how they translate to each other, a good explanation of the keyword + each named parameter, and usage examples.

## OPEN — ill-posed scenarios
- [ ] **s30_thinqfa**, **s36_readcord** — no FVSsn baseline (bad/edge param layout; SCF fields past col 80). Re-author as realistic.

## Deferred
- NE variant port (explicitly deferred while the SN drop-in is in progress).
