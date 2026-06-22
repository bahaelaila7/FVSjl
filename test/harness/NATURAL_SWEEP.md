# Natural-process congruence sweep

Verifies FVSjl ≡ Oracle A on **natural dynamics only** by masking every
management/disturbance keyword (cuts, fire/FFE `FMIN…END`, pests) from all
scenarios and diffing the full multi-cycle `.sum`.

    bash mask_natural.sh scenarios /tmp/nat_keys           # strip disturbances
    (julia --project \
        /workspace/FVSjl/test/harness/oracle_sumrow_all.jl /tmp/nat_keys)   # oracle .sum
    julia --project fvsjl_fullsum_all.jl /tmp/nat_keys /tmp/nat_mine_sums   # FVSjl .sum
    julia --project natural_diff.jl /tmp/nat_mine_sums /tmp/nat_keys        # diff
    # DYNONLY=1 restricts to TPA/BA/SDI/TopHt/QMD (drops volume)

## Findings (2026-06-22)
- **147/154 natural scenarios match at 1%.** Growth/mortality/density are
  bit-congruent: the per-cycle RNG state (DGDRIV S0) is IDENTICAL to the oracle on
  every scenario and every cycle, ITRN stays constant (natural mortality reduces PROB
  but never removes records — compaction is thinning-only), and the only residuals are
  Float32 value-ulp (off-by-1 TPA/SDI at late cycles; the REAL*4 baseline has the same).
- **Regen is the one unfinished natural process**: no scenario exercises natural
  establishment (all carry NOAUTOES, and the stocked snt01 stands produce no ingrowth
  even with it off — establishment needs an understocked/BARE stand), and ESTAB is not
  yet ported to FVSjl. Closing regen = add a bare-stand scenario + port the ESTAB model.
- Volume-only residuals (out of the growth/mortality/density/regen scope): Fort Bragg
  cyc0 scuft=0 (sp8 volume-eq selection at forest 701) + all_BK/all_DW board-foot ulp.
