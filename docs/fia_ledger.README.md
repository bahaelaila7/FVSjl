# FIA per-stand divergence ledger (`docs/fia_ledger.csv`)

Reproducible per-stand record of the FVSjl-vs-live-FVS multi-cycle differential, so a later bug fix can be
checked for STATUS FLIPS by re-running the same stands and diffing the CSV.

## Regenerate / re-run after a fix
```
for v in sn ne cs ls; do V=$(echo $v|tr a-z A-Z)
  LEDGER=/tmp/led_$v.csv julia --project=. test/harness/fia/ledger_fia.jl test/harness/fia/ledger/${v}_1000.stands $V none
done
head -1 /tmp/led_sn.csv > docs/fia_ledger.csv; for v in sn ne cs ls; do tail -n +2 /tmp/led_$v.csv; done >> docs/fia_ledger.csv
git diff docs/fia_ledger.csv   # ← any bit_exact flip or magnitude change shows here
```
Self-contained: the harness builds a temp indexed sub-DB from the read-only master + the committed stand list.
Stand set = 1000 stratified stands/variant (`test/harness/fia/ledger/<v>_1000.stands`); regime = plain (`none`).

## Columns
`variant,regime,cn,n_cycles,bit_exact,div_cols,worst_col,worst_cycle,max_rel_pct,max_abs_diff,struct_max_rel_pct,vol_max_rel_pct,density_bitexact,converges,signature`
The FACTS are primary; `signature` is a deterministic bucket over them (MATERIAL gate = >1 unit AND ≥1% rel,
so a ±1 straddle is not mistaken for a real divergence). `UNCLASSIFIED` / `volume_persistent` = needs a trace.

## Status of this run (1000/variant, plain)
bit-exact: SN 526/824, NE 794/990, CS 870/988, LS 755/999. Signatures: bit_exact 2945, print_boundary 654,
threshold_crossing 75, volume_persistent 49, structure_densephase 43, count_straddle 35.

## ⚠ REAL BUG the ledger surfaced (falsifies the earlier "no masked bug" claim)
`volume_persistent` flagged 2 SN stands (202566908010854, 162992981010854, both Region-8, LOCATION 80216 / 824)
where STRUCTURE is BIT-EXACT (TPA/BA/QMD match live) but jl computes **ZERO volume every cycle** (live TCuFt
1599→3160, jl 0). Includes six 12.2" trees that should carry large cubic. NOT a cornered primitive — a genuine
volume-computation failure on these locations. UNDER INVESTIGATION; the campaign's "complete" status is therefore
PREMATURE until this is root-caused (both-sides-traced) and fixed or genuinely cornered.

## ROOT CAUSE of the zero-volume bug (both-sides-traced)
`setup_volume_equations!` (src/engine/volume_equations.jl:84-92) decodes the national-forest code as
`iregn = kodfor ÷ 10000`, which assumes the LONG LOCATION format `REGION*10000 + FOREST*100 + DIST`
(e.g. 80216 → iregn 8). But FIA `LOCATION` also arrives in the SHORT `REGION*100 + FOREST` form
(e.g. 824 = R8 forest 24, per fia_database.jl:61-62). For a short code, `824 ÷ 10000 = 0 ≠ 8`, so the
`iregn == 8` guard fails and `species.vol_eq[sp]` is left BLANK ⇒ `_R8CLARK_VOL` returns 0 for every tree
(confirmed: a 12.5" loblolly gets cuft=0). Live FVS handles both formats and reports normal volume.
The related stand 80216 decodes to iforst=2, which `_r8_ceqn` doesn't cover — also blank ⇒ 0.
FIX (needs care, both-sides + multi-stand validation, per doctrine): normalize KODFOR to a canonical
(region, forest, district) covering BOTH formats before the `iregn==8` guard — mirroring FVS forkod.f/
sitset.f — and confirm the currently-passing stands (long-format + STDINFO-resolved forest_idx) stay bit-exact.
STATUS: OPEN. Because a real bug is open, `docs/FIA_FVS_COMPAT_COMPLETE` should NOT stand — removing it.

## ✔ RESOLVED (commit 5a4fb9f) — the flip, demonstrated
The zero-volume bug is FIXED (ported forkod.f pseudo-code remap 824→81203). Re-ran the SN ledger; the flagged
stand FLIPPED, exactly the intended workflow:
  OLD: 162992981010854 ... false ... volume_persistent   (all volume 0)
  NEW: 162992981010854 ... true  ... bit_exact            (all 10 cols match live)
SN bit_exact 526→527. This CSV now carries the post-fix SN rows (NE/CS/LS unchanged — the fix is SN-only).
The 2nd "zero-vol" hit (202566908010854, LOCATION 80216) was a false positive of the quick filter = a SCuFt
threshold-crossing (cornered), not a bug. Suite floor 38527/143/0.
