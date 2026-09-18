# Branch ledger — merged/superseded classification (2026-09-18)

Master tip at classification: `a4a1af39` → (after this) the `b94eceda` climate-completion merge.

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

## Result: the only genuinely-unmerged code was `b94eceda`. Keep-set = `master`, `treemap-explorer`.

- **`b94eceda`** (branch `ie-climate-treemult`): the climate AUTOES-ABIRTH line (1/1 absent) —
  the one real unmerged fix; **merged properly into master this session** (see git log). Branch
  then redundant.
- **`treemap-explorer`**: the Forest Growth Explorer webapp (`apps/`) — active, KEEP.
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
