# Branch ledger — merged/superseded classification (2026-09-18)

Master tip at classification: `a4a1af39`. Master has since advanced (all proper `git merge`, no
cherry-pick): `bebee0b3` (b94eceda climate AUTOES-ABIRTH) → `7dee9fa1` (this ledger) → `583da79c`
(ie-climate-4759041 = inventory-aspen ABIRTH, the third/final Climate-FVS ABIRTH fix — a real bug
found during the IE close-out, now merged). See `docs/IE_VARIANT_PORT_AUDIT.md` "IE CLOSE-OUT".

**Why this exists:** prior sessions **cherry-picked** dig-branch fixes into master, so the
same fix lives under *different SHAs* on branch vs master — `git branch --merged` can't see it,
and ~58 branches falsely appear "unmerged." Going forward: **proper `git merge` only, never
cherry-pick / rebase-then-ff**, so merged-ness stays verifiable.

**Reliable method used** (not `git cherry` / three-dot diff, both of which give false positives
here — confirmed: they falsely flagged `fix-tt-m331d`, the NC volume/DG branches, and
`fix-ktci-simfire` as KEEP when their code is in fact present in master): for each branch,
count its added `src` **code** lines (non-comment/blank) whose exact text is **absent from
master**. `0 absent` = fully merged. `N absent` on a branch whose feature is present in master
(often confirmed behaviorally at the oracle floor) = **superseded implementation drift** (old /
alternative version of a feature master already has), not an unmerged fix. A coherent block of
absent lines forming a feature master lacks = a genuine KEEP.

## Result: the only genuinely-unmerged code was `b94eceda`. Keep-set = `master`, `forest-explorer`.

- **`b94eceda`** (branch `ie-climate-treemult`): the climate AUTOES-ABIRTH line (1/1 absent) —
  the one real unmerged fix; **merged properly into master this session** (see git log). Branch
  then redundant.
- **`forest-explorer`**: the Forest Growth Explorer webapp (`apps/`) — active, KEEP.
- **All other ~58 branches**: MERGED or SUPERSEDED or STALE → delete-safe (pending owner prune).

## Per-branch (added src code-lines / absent-in-master → verdict)

MERGED (0 absent — code is in master verbatim):
audit-species-major 2/0 · dftm-ode-wip 520/0 · fix-d1b-seedling-selfthin 61/0 · fix-em-181-ccf 3/0 ·
fix-em-estab-volume 2/0 · fix-em-prefix 8/0 · fix-ie-prefix 2/0 · fix-nc-dve-hardwood-vol 32/0 ·
fix-nc-hardwood-dg 54/0 · fix-nc-r6-shorttree-cutoff 2/0 · fix-nc-redwood-brokentop-vol 18/0 ·
fix-op-native-mortality 132/0 · fix-tt-af-selfthin 80/0 · fix-tt-aspen-hcor 52/0 · fix-tt-largetree-dgf 4/0 ·
fix-tt-m331d 4/0 · fix-tt-smdgf-realization 12/0 · fix-tt-sparse-woodland 30/0 · fix-tt-woodland-dg 13/1 ·
fix-ut-woodland 2/0 · fix-ut199-woodland 45/0 · wt/western-sweep 2/0 · worktree-agent-a77fd011 318/0 ·
worktree-agent-a79bc85f 131/0 · worktree-agent-a7e50ffb 20/0

STALE (no src code — data/test/docs/results/probe only):
audit-crosswalks · diag-em-ie-estab · fix-ba-underbuild · fix-harness-seedtest · fix-ie-cohort-rebalance ·
fix-tt-selfthin · fix-ut-sprout · worktree-agent-a9e5ccb4 · worktree-agent-ad16f7908

SUPERSEDED (absent lines are an early/alternative impl of a master-present feature; verified via
symbol-presence + behavioral oracle-floor where a repro exists):
estab-ie-discrete 26/3 · fix-d1-postthin-count 29/5 · fix-d2-regen-height 5/1 · fix-em-aspen-abirth 64/1 ·
fix-em-crownmodel 33/4 · fix-ie-autoes 40/7 · fix-ie-cohort-dist 181/17 · fix-ie-cyc1 19/1 ·
fix-ie-m333-autoes 31/3 · fix-ie-nstore 1/1 · fix-ie-oh-growth 20/9 · fix-ie-singlepoint-autoes 29/6 ·
fix-ktci-simfire 16/2 · fix-nc-largetree-tcubic 78/2 · fix-western-simfire 20/2 · ie-estab-height 53/5 ·
ie-smalltree-dg 10/9 (discarded IEDG_NOSORT probe) · port-op-wscwhr 26/2 · svs-port 228/26 (chunk-0
scaffolding; master has full SVS) · wt/wdig 4/1 · worktree-agent-a238e3041 166/6 (WC foundation) ·
worktree-agent-a7ff266c 40/1 (OC foundation) · worktree-agent-ae6d9269 156/1 (UT htg) ·
worktree-agent-af975d09 127/4 (AK foundation) · worktree-agent-afae4bde 61/3 (OC SWO htg)

Foundations (WC/OC/BM/AK/UT) are early chunks of variants now COMPLETE in master (all variant
modules present in `src/variants/`, full `Pkg.test` 55329/25/75 passes across variants).

## Caveat
"Superseded" rests on: (a) master containing the branch's feature symbols (direct grep on master),
and (b) where a repro exists, master sitting at the oracle floor on it. It does NOT re-litigate
whether master's *variant* is itself fully correct (e.g. NC is a "COMPLETE→retracted→re-sweep"
open campaign item) — that is a master-level sweep question, independent of these stale branches.

---

# 2026-09-19 — BM regime-close campaign branches

Integration branch: **`bm-regime-close`** (worktree `/workspace/.wt-bm`). Every fork branch below was landed
into it by a **proper `git merge --no-ff`** (verified mechanically: `git merge-base --is-ancestor <branch>
bm-regime-close`). `bm-regime-close` itself goes to `master` by proper merge after the final full-suite gate.

**Prune rule:** delete a branch + its worktree only when (1) it is an ancestor of `master` (i.e. after the
`bm-regime-close → master` merge), (2) no ACTIVE fork references its untracked `.dig/` tooling (instrumented
oracles, hex-dump scripts), and (3) the owner approves the prune. Worktree `.dig/` dirs are untracked scratch —
copy anything worth keeping (instrumented-oracle recipes) into the repo or `/workspace` first.

| Branch | Tip | Status | Worktree | After BM→master |
|---|---|---|---|---|
| `bm-regime-close` | (integration) | ACTIVE — integration | `/workspace/.wt-bm` | delete after master merge |
| `forest-explorer` | 9fccd030 | ACTIVE — Forest Growth Explorer webapp; USER 2026-09-20: keep its development on this branch in its OWN worktree so the main checkout stays on master | `/workspace/.wt-explorer` | KEEP (long-lived) |
| `bm-cyc0-vol` | a18802df | merged (R6 MRULES, LHTDRG/HTDBH dub, SF_HS) | `.wt-bm-vol` | delete |
| `bm-thinprsc` | d9e098b2 | no commits (found the g16 `-std=legacy` oracle bug); worktree = **master-equivalent BASELINE** used for base-vs-tip comparisons | `.wt-bm-prsc` | delete (keep worktree until BM closed) |
| `bm-plant-seedling` | 50dded40 | merged (BM strp estab, PLANT cycle-date) | `.wt-bm-plant` | delete |
| `bm-simfire` | 0cfe411b | merged (BM FFE 7 gaps) | `.wt-bm-fire` | delete |
| `bm-mistletoe` | de8d5143 | merged (mistletoe post-TRIPLE order, REGENT per-copy ZZRAN) | `.wt-bm-mist` | delete |
| `bm-rootdis` | e65c2b0f | merged (RDTDEL) | `.wt-bm-rd` | delete |
| `bm-sprout` | b934aa7b | merged (BM ESUCKR/ESSPRT) | `.wt-bm-sprout` | delete |
| `bm-crwdth` | 414ed34f | merged (one variant CRWDTH) | `.wt-bm-crw` | delete |
| `bm-econ` | fb3f1ac5 | merged (ECON engine) | `.wt-bm-econ` | delete |
| `econ-leftovers` | 52e7f988 | merged (PRETEND, THINPT, CUTLIST…) | `.wt-econ2` | delete |
| `ie-addtrees-bridge` | 12075c1d | merged (ADDTREES at ESNUTR head) | `.wt-ie-addt` | delete |
| `estab-stoadj` | 9dbf2623 | merged (STOADJ persistence) | `.wt-stoadj` | delete |
| `ie-estab-draworder` | 793fd960 | merged (IE birth-height draw order, KODFOR) | `.wt-ie-draw` | delete |
| `ieem-postfire-estab` | 635a4bc6 | merged (IE/EM post-fire estab, CUTS-entry TREDEL) | `.wt-postfire` | delete |
| `bm-cover-strclass` | b51975e3 | merged (COVER spmaps, StrClass SSTGHP) | `.wt-bm-cov` | delete |
| `bm-cover2` | f0d6adb5 | merged (cycle-0 COVER, unified CRATET IND) | `.wt-bm-cov2` | delete |
| `bm-base-resid` | b7bf4886 | merged (BM base residuals r1) | `.wt-bm-base` | delete — `.dig/` tooling reused by r3 |
| `bm-base-resid2` | 14f234d4 | merged (BM base residuals r2) | `.wt-bm-base2` | delete — `.dig/` tooling reused by r3 |
| `bm-simfire2` | b67b81ab | merged (post-fire negative ICR) | `.wt-bm-fire2` | delete |
| `fix-wpbr-numtrip` | a9da88f5 | merged (WPBR port, EM plant crown positions) | `.wt-regfix` | delete — `.dig/wpbr` oracle reused by wpbr-tripled |
| `treelist-ulp` | 1fb43056 | merged via `ie-ulp-atrt` (TREELIST gating, birth TPA) | `.wt-trl` | delete — `.dig/ulp` hex oracle reused by ie-ulp3 |
| `ie-ulp-atrt` | 111d7d82 | merged (ATRTList, DSNOUT, IE ULPs cyc≤2) | `.wt-ulp` | delete — `.dig` reused by ie-ulp3 |
| `fix-treeszcp-wpbr` | 0d56d7d5 | merged (fpowi, PCTILE, WPBR activities) | `.wt-tsz` | delete — `.dig_wpbr` reused by wpbr-tripled |
| `bm-base-resid3` | 352a31ef | merged 877064b5 (REAL*4 op forms; BM none 386→398/400) | `.wt-bm-base3` | delete |
| `wpbr-tripled` | 91fd73f6 | merged df8d3c38 — WPBR on tripling cycles, bit-exact on a new tripled fixture; test_wpbr 54/2 → 56/56. The "new" test_root_disease failure was a STALE BASELINE (the tip itself is 1157/2). | `.wt-wpbr3` | delete |
| `list-ht2td` | 50907c14 | merged (Ht2TD columns + NVEL precision) | `.wt-ht2td` | delete |
| `ie-ulp3` | (moving) | **ACTIVE** — IE control-stand divergence (TREEFMT back-tab fixture, print-visible from cycle 1) + per-record growth from cycle 3; also owns the dropped `ie_esnspe` per-plot-BAA finding. | `.wt-ulp3` | merge, then delete |
| `cs-bdft` | 00fcf75f | merged 92d39ea3 — kwcov allowlist generator threw (every stem permanently "broken"); cs_numtrip + cs_treeszcp were stale. | `/workspace/.wt-cs` | delete |
| `tiered-suite` | 98d2b6ef | merged b4ec90f1 (tiered suite + 32MB fixtures, 4 variants) | `.wt-tiered` | delete |
| `ie-climate-4759041` | 29afeea2 | merged (2026-09-18; content == 583da79c) | `/workspace/.wt-regendg` | delete |
| `ie-dfht` | 38091a40 | merged (2026-09-18 bare-plot chain) | `FVSjl/.claude/worktrees/agent-a7c278e41effcaf96` | delete |


Removed already this session: `/workspace/.wt-bisect2` (temp bisect, culprit f26d7cb4), `/workspace/.wt-bm-bisect` (detached, clean), `/workspace/.wt-regfix-bis` (by its fork).

Non-git scratch to clean with the prune: staged oracle binaries `/workspace/.{v}work/FVS{v}_g16.new` (pending the
owner's oracle-swap decision — see `/workspace/ORACLE_SOURCE_AUDIT_2026-09-19.md`), `/workspace/relink_main_nolegacy.sh`,
`/workspace/validate_relink.sh`.

---

## Addendum 2026-09-23 — the westside shared-fix chain

Master advanced to `b56bfeda` by a proper `git merge --no-ff` of the chain below (each branch
built on its predecessor, so ONE merge lands all four). Full `Pkg.test` on the merged tip:
**56533 pass / 23 fail / 45 broken**, by-name identical to the pre-merge tip (removals only).

| branch | tip | worktree | content | verdict |
|---|---|---|---|---|
| `pn-regime-close` | d4ed4e0b | `.wt-pn` | absolute-row WRD tests; PN ATTEN/RELDEN/point-PRD; PN control BA +62 → ±1 | merged in `b56bfeda` |
| `so-dig` | b559337e | `.wt-so` | SO point CCF, SO species-major REGENT, `variant_bratio` unification; SO control BA +14 → ±1 | merged in `b56bfeda` |
| `deadcrown-all` | 1f3fbedb | `.wt-dead` | dead-record crown dub + species-major crown loops for the remaining variants; NC/WS DUBSCR; CI TEMMAI; RMAI grinit defaults; shared CRATET LSTART | merged in `b56bfeda` |
| `oc-dub` | 6a131809 | `.wt-oc` | OC LSTART crown dub (oc/crown.f) + the shared SDIAC the Weibull dub was missing | merged in `8b31b046` |
| `crown-isort` | a6d3716e | `.wt-isort` | one shared `crown_isort` (key = current DBH) + WRD tolerances reconciled | merged in `8b31b046` |
| `em-crown` | 845e593c | `.wt-em` | EM four-class crown restructure, EM REGCAL, LHTDRG defaults, tree-count-stratified sampler | **NOT merged — EM tiered fixture being regenerated** |
| `forest-explorer` | 1b28a1d4 | `.wt-explorer` | Forest Growth Explorer webapp (`apps/forest-explorer`) | KEEP, unmerged by design |

`.wt-embisect` is a throwaway bisect worktree over the `em-crown` chain (detached HEAD); delete
with the prune.

**Six shared-mechanism bugs** landed by this merge — each had been a hand-kept per-variant list
that drifted out of sync with the Fortran, so the fix is one shared implementation, not 23 copies:
dead-record crown dub (`crown.f DO 79`), ATTEN (`dgdriv.f` SIGMA pooling), RELDEN for every
variant, species-major REGENT/crown draw order (`IND1`), `variant_bratio` as the single BRATIO
dispatch, and the per-variant RMAI `grinit` default. See
`docs/PORT_STATUS.md` "Westside shared fixes".

---

## Addendum 2026-09-25 — the crown-initialisation chain

Master advanced `bb30d258` → **`5930ffb5`** by one proper `git merge --no-ff` of the chain below.
Gate: full `Pkg.test` on the tested tip `0b23accd` = **56534 pass / 23 fail / 0 err / 44 broken**,
by-name diff against master **0 new / 0 gone / 0 changed** (the 23 are the pre-existing ON
volume-dump + LPMPB set). The merged tree is byte-identical to that tip.

| branch | tip | worktree | content | verdict |
|---|---|---|---|---|
| `crnmult` | 2e4c47f7 | `.wt-isort` | CRNMULT applies (BM); shared `topkill_icri` (stmt 55); TT SDIAC | merged in `5930ffb5` |
| `kt-lstart-dub` | cac8edba | `.wt-kt` | KT LSTART dub; shared LSTART PCT ordering; IE/AK stmt 55 | merged (via `crnmult`) |
| `crnmult-west` | 7bc24129 | `.wt-crnw` | CRNMULT wired for CR | merged in `5930ffb5` |
| `bc-lstart-dub` | 0b23accd | `.wt-kt` | BC LSTART dub + stmt 55; the KT/simfire golden reconciliations; the empty-stand crash guard | **the merge head** |

Worktrees now merged or throwaway and safe to prune with the owner's approval:
`.wt-isort`, `.wt-kt`, `.wt-crnw`, `.wt-oc`, `.wt-dead`, `.wt-em`, `.wt-embisect`.
`forest-explorer` (`.wt-explorer`) remains the only KEEP.

Three goldens were reconciled **from measurement**, both directions recorded in each case rather
than only the loosening — see the commits: KT's inert-seam self-snapshot, KT's absolute tolerances
(TPA tightened 7→4 and 8→5, BA loosened 6→10, first four cycles now exact), and KT's simfire
pre-fire TREES 282→283, which lost an exact cell while halving the BdFt error (−814 → −479).

## Addendum 2026-09-25 (b) — `crnmult-rest` merge

Branch `crnmult-rest` (worktree `.wt-crn2`), off master `a5825bb8`, landed with one proper `git merge --no-ff`.

It carries:
- CRNMULT at every `crown.f` site for PN, SO, EC, WS, CA, NC, CI, TT, UT, EM, IE, WC and AK.
- The STDINFO habitat decode for WC, PN, SO, CA, NC and UT.
- Faithful SO, CA and NC site setup, including CA Region 5 and the full `nc/forkod`.
- Western volume fixes: inside-bark tops from the start-of-cycle BARK, CFTOPK trims, R5 tables for SO and
  CA, NC merch specs, and the 712 BLMVOL path.
- The CEPMRT/SLPMRT self-thinning latch for EM, NC, UT and TT.
- WS R5CRWD CCF and DVE board-foot.

Details are in `docs/PORT_STATUS.md` "CRNMULT rollout, STDINFO habitat, site setup and western volume".

Gate: the full `Pkg.test` on the final tip is compared by name against master's baseline. The tiered
allowlist was reconciled in both directions (38 tightened, 10 widened, EM latch verified against live), and
UT's WRD row was promoted to bit-exact.

`.wt-crn2` is safe to prune once merged, with the owner's approval.

The oracle-side changes of 2026-09-25 (debug WRITEs removed, CR varmrt guard) are recorded in
`/workspace/ORACLE_SOURCE_AUDIT_2026-09-19.md` §6.

## Addendum 2026-09-25 (c) — `nc-dg` merge

Branch `nc-dg` (worktree `.wt-crn2`, off master `0b09d876`), one proper `git merge --no-ff`. It sets NC's
cycle-1 AUTCOR OLDFNT to grinit's FINT (10) through `dg_measure_period`. This makes NC bit-exact on every
habtest case, and its WRD absolute row is promoted to a passing test. See docs/PORT_STATUS.md item 7.

Gate: the full Pkg.test on the final tip is compared by name against the master `0b09d876` log
(`/workspace/.postswap/pkgtest_master_0b09`).

## Addendum 2026-09-26 (d) — `ws-htg` merge

Branch `ws-htg` (worktree `.wt-ws`, off master `0b09d876`; master `566b618c` merged in as `f22ece52`), landed with one
proper `git merge --no-ff`. It ports the WS `htgf.f` surrogate branches, WS REGENT, WS forkod, WS point PRD,
R5 board MERLEN and the SO Fremont volume table, and gives every variant its own PSIGSQ. See docs/PORT_STATUS.md
("WS height growth, REGENT and site; per-variant PSIGSQ").

Gate: the full Pkg.test on the final tip is compared by name against the master `566b618c` log
(`/workspace/.postswap/pkgtest_master_566b`). The only by-name changes were three PPE MXHRVP EC self-snapshots,
each moving toward the FVSppe oracle; they were re-pinned.

## Addendum 2026-09-26 (e) — `em-vol` merge

Branch `em-vol` (worktree `.wt-em2`, off master `566b618c`; master `f6b7ca75` merged in), landed with one proper
`git merge --no-ff`. It carries:
- EM/KT top-killed trees on NORMHT, and EM Custer PP on 203FW2W122.
- The Fortran-shaped EM REGENT, including ESTAB.
- `em/bratio.f`.
- The cycle-1 WK1 dub.
- LL crown OBA/RDM1/OLDPCT.
- EM HTGMULT and the per-copy LL HTG.

See docs/PORT_STATUS.md, "EM REGENT in the Fortran's shape".

Gate: the full `Pkg.test` on the final tip, compared by name against the master `f6b7ca75` log
(`/workspace/.postswap/pkgtest_wshtg`). The tiered EM allowlist was redrafted from measurement, and the EM WRD row
and the new EM REGENT test pass.

## Addendum 2026-09-26 (g) — `em-regcal` merge

Branch `em-regcal` (worktree `.wt-emr`, off em-vol `27da808a`), landed with one proper `git merge --no-ff` after
master was merged in. It carries the EM and IE LSTART REGCAL fixes:
- Per-species RHCON.
- The cratet backdating-DENSE snapshot.
- IFINTH.
- The cratet.f:153 tie IND.
- CCF at D = 0.
- IE's CR/UT arms.

It also adds the test `test_regcal_em_ie.jl`. See docs/PORT_STATUS.md, "LSTART small-tree height calibration
(REGCAL) for EM and IE".

Gate: the full `Pkg.test` on the final tip, compared by name against the preceding master log. The tiered EM
allowlist was redrafted from measurement, keeping the #247 annotation.

