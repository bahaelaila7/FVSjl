# Bugs found in the live Fortran FVS source (the oracle) — jl is CORRECT; do NOT replicate

The live Fortran binaries (`/workspace/ForestVegetationSimulator/bin/FVS{sn,ne,cs}_buildDir`)
are the campaign's oracle. They are normally ground truth, but the FVSjl differential work has
surfaced places where **the real Fortran itself is wrong** — it emits incorrect or physically-
impossible output. Per the campaign doctrine ("the oracle can be WRONG; don't blindly match")
and the user's explicit decision (2026-07-03), **FVSjl stays CORRECT on these paths and does NOT
reproduce the bug.** They are recorded here (separate from the Oracle-A transliteration bugs in
[ORACLE_A_BUGS.md](ORACLE_A_BUGS.md)) so that:
- FVSjl↔live `.sum`/report diffs on these columns are understood as live-side, not jl regressions;
- they can be reported upstream to the FVS maintainers;
- a future "strict byte-fidelity" mode (if ever wanted) knows exactly what to gate.

Ground truth = the FVS source under `bin/FVS{sn,ne,cs}_buildDir/*.f` (cited per entry).

Format per entry:
> ### <title>
> - **Live Fortran:** file.f:line — what it does (wrong)
> - **FVSjl:** correct? (and where)
> - **Evidence it's a bug** (internal inconsistency / source trace)
> - **Impact / how found**
> - **Disposition**

---

### D36 — econ HRVRVN leaves removed-sawlog-cubic (`OSCREM(7)`) stale into a no-cut cycle
- **Live Fortran:** the `.sum` "removed cubic saw volume" column (col 16) = `IOSUM(22,IKNT)` =
  `INT(OSCREM(7)/GROSPC+0.5)` (`disply.f:342`). `ONTREM(7)` (removed trees) is explicitly re-zeroed
  every cycle (`cratet.f:158`, `fvs.f:433`); the full `O*REM(I)` removal arrays are zeroed only inside
  `cuts.f:324-329` **during an actual cut**. There is **no matching per-cycle reset of `OSCREM(7)`**
  (grep-confirmed). With the ECON `HRVRVN` keyword active, the econ path re-populates `OSCREM(7)` after
  the cut and leaves it set, so the NEXT (no-cut) cycle's summary row reports the prior cut's value.
- **FVSjl:** CORRECT — it zeroes all removal columns on a no-cut cycle (`remScuft = 0`).
- **Evidence it's a bug:** the live 2005 row for `econ_strtecon` is `rem: 0 0 0 23 0` — i.e. **23 cuft
  of sawlog "removed" with 0 removed-trees AND 0 removed-total-cubic**, which is physically impossible
  (you cannot remove sawlog cubic volume while removing zero trees and zero total cubic). The thin year
  (2000) is bit-exact both engines (`rem: 20 106 97 23 106`); only the subsequent no-cut cycle diverges.
- **Impact / how found:** found 2026-07-03 by broadening `divergence_sweep.jl` to the REMOVALS columns
  (13–17), which the state-only sweep never checked. Scope: TRIGGERED BY `HRVRVN` — `cut_thinsdi`
  (identical THINSDI, no econ) correctly zeroes `remScuft` at 2005; only econ-harvest stands
  (`econ_strtecon`, `econ_u5`) exhibit it. Affects only the `.sum` removal column; the ECON revenue
  TABLES themselves are bit-exact (validated separately).
- **Disposition:** 📌 jl stays correct; NOT replicated (user decision 2026-07-03). Documented here + as
  D36 in [DIVERGENCE_FIX_CAMPAIGN.md](DIVERGENCE_FIX_CAMPAIGN.md). Upstream-reportable to FVS maintainers.
