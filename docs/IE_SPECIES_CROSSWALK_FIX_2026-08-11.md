# IE species crosswalk — zero-padded FIA code mis-map (#175) — FIXED 2026-08-11

## Symptom
While triaging #171 (IE larch-heavy stand 373781950489998 BA over-growth) with FVSie_g16 per-tree DG
instrumentation, the subalpine-larch trees (FIA 072) resolved to **live ISPC 14 (LL, Larix lyallii)** but
**jl internal sp2 (WL, Larix occidentalis / western larch)** — a wrong-species assignment (wrong DG /
mortality / volume / crown coefficients).

## Root cause (MEASURED, not inferred)
- The FIA-DB reader `_fia_spcode` (fia_database.jl:49) zero-pads 2-digit FIA codes to 3 chars: `"72"→"072"`.
- Each variant's coef stores its own FIA codes UNPADDED: IE `sp.fia[14]="72"`, `sp.fia[2]="73"`.
- `resolve_species` (intree.f:240) tries a DIRECT match against the variant's own alpha/FIA/PLANTS first,
  then falls back to the SPCTRN "unknown species" crosswalk. The padded `"072"` fails the exact FIA match
  (`"072" != "72"`) → falls through to SPCTRN → data/inlandempire/species_translation.csv row
  `LL,072,LALY,WL,WL,WL,WL` folds subalpine larch → **WL**. That fold is correct for KT (11 species, no LL)
  but wrong for IE (23 species, LL = sp14). The IE CSV was copied wholesale from KT ("target_kt" headers).
- Proven: FVSie_g16 dgdriv.f per-tree dump shows FIA-072 trees at ISPC 14; the jl coef CSV has sp14=LL/FIA72.

## Scope (measured across all 8 western variants)
Padding-normalized "does each variant's own 2-digit FIA code resolve back to itself?" scan:
- **IE only**: 2 mismatches — FIA 072 (LL subalpine larch) → WL(sp2), FIA 066 (RM Rocky-Mtn juniper) → OS(sp23).
- KT/EM/UT/TT/CI/BM/CR: **0 mismatches** (their SPCTRN folds happen to land on the correct own species).

## Fix
`resolve_species` (src/engine/species_translation.jl): leading-zero-normalize the DIRECT FIA match for
all-digit codes, so `"072"` direct-matches the variant's own `"72"` species instead of falling through to
the coarser SPCTRN crosswalk. Non-numeric codes (alpha/PLANTS) are exact-matched exactly as before.

## Validation
- resolve_species("072")→sp14(LL), ("066")→sp16(RM), ("073")→sp2(WL), ("OT")→sp23(OS): all correct.
- Post-fix jl sp14 (LL) per-tree DG on 373781950489998 MATCHES FVSie_g16 (WK2 bit-exact to ~4 decimals,
  DG within the DGSCOR realization straddle) — the larch now grows as LL exactly like live.
- Cyc0 BIT-EXACT jl vs FVSie_clean on 3 larch/RM-juniper FIA stands (31366571010690, 1592831892290487,
  39447478010690); 0 crashes. Stand 373781950489998 cyc0 CCF 86→99 = live (was wrong via WL crowns).
- iet01 BYTE-IDENTICAL (alpha-code treelist → unchanged first branch of resolve_species).
- Regression-free for the other 7 western variants: 0 own-species mismatches ⇒ direct-match lands on the
  same species SPCTRN already gave; alpha/PLANTS branches untouched.

## Relationship to #171
This fix UNMASKS #171: on larch-heavy stands the wrong WL assignment UNDER-grew the larch (~21% of TPA on
373781950489998), partially offsetting the pre-existing small-tree (<3") regent over-growth. With the larch
now correctly growing as the faster LL, the stand's multi-cycle BA over-growth is larger (+7.5%→+14% at
2025) and more clearly attributable to #171 — which FVSie_g16 confirms is NOT the large-tree DG (both AF and
LL large-tree DG are bit-exact) but the small-tree regent cohort (TT #158-class residual). #171 remains open.
