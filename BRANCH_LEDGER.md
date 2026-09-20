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
| `wpbr-tripled` | 523c2b97 | **WIP, UNVERIFIED** — fork killed by API rate limit mid-verification; it reported test_root_disease 1157/2 vs baseline 1158/1 (new failure unidentified). DO NOT MERGE until traced. | `.wt-wpbr3` | verify → merge or discard |
| `list-ht2td` | 50907c14 | merged (Ht2TD columns + NVEL precision) | `.wt-ht2td` | delete |
| `ie-ulp3` | e30c58ee | **WIP, UNVERIFIED** — fork killed by API rate limit mid-trace (LHTCAL grinit.f:104 / initre.f:2470-2490 lead). DO NOT MERGE until verified. | `.wt-ulp3` | verify → merge or discard |
| `tiered-suite` | 98d2b6ef | merged b4ec90f1 (tiered suite + 32MB fixtures, 4 variants) | `.wt-tiered` | delete |
| `ie-climate-4759041` | 29afeea2 | merged (2026-09-18; content == 583da79c) | `/workspace/.wt-regendg` | delete |
| `ie-dfht` | 38091a40 | merged (2026-09-18 bare-plot chain) | `FVSjl/.claude/worktrees/agent-a7c278e41effcaf96` | delete |


Removed already this session: `/workspace/.wt-bisect2` (temp bisect, culprit f26d7cb4), `/workspace/.wt-bm-bisect` (detached, clean), `/workspace/.wt-regfix-bis` (by its fork).

Non-git scratch to clean with the prune: staged oracle binaries `/workspace/.{v}work/FVS{v}_g16.new` (pending the
owner's oracle-swap decision — see `/workspace/ORACLE_SOURCE_AUDIT_2026-09-19.md`), `/workspace/relink_main_nolegacy.sh`,
`/workspace/validate_relink.sh`.
