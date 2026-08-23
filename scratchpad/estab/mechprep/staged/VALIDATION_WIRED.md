# MECHPREP / BURNPREP live wiring — RESOLVED BIT-EXACT 2026-08-23 (post AUTOES #143)

## TOP-LINE: MECHPREP_BURNPREP_LIVE_BITEXACT

The 2026-08-23-earlier re-scope (VALIDATION.md) deferred live wiring on FOUR blockers. AUTOES #143
(081e0e6c) resolved the load-bearing one (C, the disturbance-tally seed desync) AND, as a downstream
consequence, blocker (D). With the seed stream now bit-exact, MECHPREP/BURNPREP were WIRED into the live
establishment path and validated **byte-identical vs FVSie_g16** on the AUTAL disturbance tally.

Fixture: under-stocked IE `753189105290487` (ie_understocked.db), `THINBTA 2029 10.0` (near-clearcut ⇒ AUTAL
disturbance tally at cyc2) + `ESTAB` + MECHPREP/BURNPREP. Oracles: FVSie_g16 (.sum) and FVSie_estabdump (IPPREP).

## What the #143 fix unblocked (measured)

| blocker (VALIDATION.md) | prior state | now |
|---|---|---|
| **C** jl omits cyc-1 ingrowth ⇒ disturbance seed 43303 vs oracle 61997 | BLOCKING | **RESOLVED**: jl now fires cyc-1 ingrowth (43303), cyc-2 disturbance grabs **61997** bit-exact, cont 61997, cyc-5 49053 — all seeds aligned |
| **D** jl disturbance establishment cornered (cyc-3 TPA 1371, SDI 310748) | BLOCKING | **RESOLVED** (downstream of C): jl h_base .sum now byte-identical to oracle (2039 TPA 262/BA 99 … 2069 690/190) |
| **A** prep fires only on the AUTAL tally | reachable via THINBTA | reached (cyc-2 NTALLY=1) |
| **B** PROB1 prep-invariant (SPRE=0) | true | true, and irrelevant: the effect is the species-mix (CPRE), which the tally carries |

## Live A/B — jl vs oracle .sum, FULL data rows BYTE-IDENTICAL

| keyfile | prep | jl == oracle |
|---|---|---|
| h_base | none | ✓ byte-identical |
| h_mech100 | MECHPREP 2029 100 (all plots IPREP=2) | ✓ byte-identical |
| h_burn100 | BURNPREP 2029 100 (all plots IPREP=3) | ✓ byte-identical |
| h_mech50 | MECHPREP 2029 50 (WK6-sampled NONE/MECH mix) | ✓ byte-identical |

The prep effect is REAL (not a no-op): jl h_mech100 differs from jl h_base exactly as the oracle does —
2049 TPA 443→441, 2069 BA **190→197**, QMD 68→70 (a species-composition shift toward larger-DBH species;
count near-invariant). h_mech50 being byte-identical proves the WK6 vector jl samples IPPREP from is the SAME
as the oracle's — the AUTOES seed alignment carries all the way through the disturbance-tally site-prep draw.

Direct IPPREP proof (FVSie_estabdump, h_mech50 disturbance tally, SUMUP 0.5/0.5/0):
```
oracle ZIPPREP = 1 2 1 2 2 2 2 2 1 1 2 1 1 1 1 2 1 1 2 1 2 2 2 1 2 1 2 1 2 2 1 2 1 2 2 1 2 1 1 2 2 1 2 2 1 1 1 1 1 2
jl (WK6 off the aligned disturbance seed 61997) == oracle, bit-for-bit.
```

## Wiring (all in the IE establishment path; inert without a MECHPREP/BURNPREP keyword)

- `keyword_dispatch.jl` kw_estab!: MECHPREP/BURNPREP branch schedules activity 493 (MECH) / 491 (BURN) with
  the %-of-plots field (esprin.f), mirroring TALLY.
- `variants/inlandempire/establishment.jl`:
  - ported `ie_esetpr` (keyword→PMECH/PBURN/IALN), `ie_esetpr_normalize` (→SUMUP), `ie_esetpr_sample`
    (estab.f:382-399 sample-without-replacement WK6→IPPREP).
  - `ie_autoes_establish!`: on the DISTURBANCE tally only (`_ntally==1 && !is_ingro`), gathers the scheduled
    493/491 %-plots at the disturbance date → `prep_sumup`, passes it through `ie_autoes_run`.
  - `ie_autoes_tally`: captures the DUPNPT WK6 site-prep draws (previously discarded), samples per-plot IPPREP,
    and uses each plot's IPREP in the species-mix tables (`_prep_tables(ip)` memoizes ie_espadv/ie_espxcs on
    IPREP 1..4). No keyword ⇒ `prep_sumup=nothing` ⇒ every plot IPREP=1 ⇒ byte-identical to the pre-wire path.

Scope notes (faithful, documented): (1) heights stay at the XMIN+0.2 floor (the computed ESSUBH heights fall
below it — bit-exact here, so the UPRE height-prep term is floored out, matching the oracle). (2) TIME is
prep-invariant when the prep date == the harvest date (TIME=FTEMP−ZMECH=FTEMP−ZHARV) — the validation set.
A delayed prep (prep date ≠ harvest) would need the estab.f:626-644 TIME reset (not exercised here). (3) the
estab.f NCOUNT prep-type reordering within a point only matters for idup>1 with a MIXED prep — for the
reachable fixtures (idup=1) each point has one plot so order is unaffected.

## Gate + tests

- `test/integration/test_multicycle.jl` = **339 / 11** byte-identical (change inert without the keyword).
- `test/unit/test_estab_mechprep.jl` (11/11): ie_esetpr parse + normalize + sampler + the LIVE IPPREP off the
  aligned disturbance seed 61997 == oracle h_mech50 vector.
- IE estab suites (test_ie_estock, test_estab_specmult_htadj, test_estab_minplots) still pass — the tally
  refactor is byte-identical on the no-prep path.

## Reproduce
```
d=/workspace/FVSjl/scratchpad/estab/mechprep/staged
for k in h_base h_mech100 h_burn100 h_mech50; do echo $d/$k.key | /workspace/.iework/FVSie_g16; done   # oracle
julia --project -e 'using FVSjl; for k in ["h_base","h_mech100","h_burn100","h_mech50"]; FVSjl.run_keyfile("'$d'/$k.key"; variant=FVSjl.InlandEmpire()); end'
# then diff <(tail -n+2 oracle_h_mech50.sum) <(tail -n+2 jl_h_mech50.sum)  → identical
```
