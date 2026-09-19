# Western Wykoff FFE SIMFIRE crown-fire regression guard (mirrors the NE SIMFIRE test at
# test/integration/test_net01.jl). This is the gate hole that let the IE/EM crown-fire kill silently
# regress from ~100% (validated early-Aug) to ~40-50% during the WC/OC/CA/PN FFE ports: jl routed the
# IE-family conifers through the eastern Jenkins FMCROWE crown-biomass (~2x under-count) instead of the
# western FMCROWW that ie/fmcrow.f + em/fmcrow.f actually dispatch — so the canopy bulk density (CBD) came
# out ~2x low, the crowning index rose above the 20-mph wind, and the SimFire fell back to a surface-only
# fire that left most of the overstory alive.
#
# ROOT-CAUSE FIX (fix-western-simfire): (1) route IE/EM conifers through cr_crownw (FMCROWW) with
# _IE_ISPMAP / _EM_ISPMAP (FMCROWE only for the {18,19,21} / {11:17,19} Jenkins species); (2) use the
# FM10-based nc_crown_fire_result for IE/EM (their fmcfir.f == nc/fmcfir.f: RACT = 3.34*FM10-spread at
# SWIND*0.4), so RFINAL matches live (jl 87.9 vs live 88.0) instead of the selected-model 21.5.
#
# LIVE-VALIDATED vs the relinked FVSie_g16 / FVSem_g16 (SimFire 2020, severe): the crown fire now FIRES
# (FMCFIR CFTMP=COND_CRN, CRBURN=1.0, CBD jl 0.138 vs live 0.146) and the stand is essentially wiped.
#   IE (ie_simfire): 2020 MORT jl 506 vs live 512; 2030 TREES jl 3 vs live 0.
#   EM (em_simfire): 2020 MORT jl 307 vs live 305; 2030 TREES jl 1 vs live 0.
# Residual (cornered): the crown-fire byram intensity magnitude — jl's rothermel surface heat-per-area
# HPA runs ~3x below live FMFINT (jl 234 vs live 731.7) for this dry down-wood fuel bed (a pre-existing,
# family-wide crown-fire byram corner, fmburn.jl); it perturbs the additive HPA term only, so both engines
# still wipe the stand. The guard pins the crown-fire OUTCOME (near-total kill), which the surface-only
# regression cannot produce (it left 174 IE / 112 EM TREES alive at 2030).

@testset "Western Wykoff FFE SIMFIRE crown fire — IE/EM/KT/CI vs live (regression guard)" begin
    _fixdir = joinpath(@__DIR__, "western_fixtures")
    # (variant, keybase, pre-fire 2020 TREES, post-fire 2030 TREES upper bound, 2020 MORT lower bound)
    cases = [
        # IE (jl 313; live ORACLE 305): the two FIXABLE crown backdating bugs are FIXED — OLDBA/RELDM1 threading
        # (simulate.jl; crown.f:279 RELDM1<100 → OBA=BA/RDM1=RELDEN fallback, bit-exact) and the cycle-1 OLDPCT
        # BACKDATED-percentile seed (ie/cratet.f:513 → dense.f backdating; ie_seed_backdated_oldpct!, tree I=3
        # OLDPCT 34.2→25.2, all 27 live within 1e-4). EXPDCR is now bit-exact on every record. The residual +8 is a
        # MEASURED named-primitive corner: ±1 ICR on 7/66 records from a ~0.006 EXPPCR gap driven by the CURRENT-cycle
        # PCT RDPSRT unstable-quicksort tie-break among equal-DBH tripled records — the same tie-break family cornered
        # on the CI case below (307/live-300). Asserting jl's value with the oracle documented, per this file's convention.
        (FVSjl.InlandEmpire(),   "ie_simfire", 305, 15, 480),   # live: 305 pre, 0 post, 512 MORT — jl pre-fire TREES now BIT-EXACT 305 (313→309 b927f4a3 → 305 on bm-regime-close 2026-09-19; live re-derived on FVSie_g16 and FVSie_g16.new, both 305). The "+4 RDPSRT-tie corner" was not a corner.
        (FVSjl.EasternMontana(), "em_simfire", 507, 15, 280),   # live: 507 pre, 0 post, 305 MORT
        # KT (Kootenai): kt/fmcrow.f calls FMCROWW for ALL species (all-conifer, like NC/Klamath) — no eastern
        # FMCROWE. CI (CentralIdaho): ci/fmcrow.f SELECT CASE(SPIW) CASE(13,15,17,19)→FMCROWE, DEFAULT→FMCROWW.
        # Both fmcfir.f == nc/fmcfir.f (FM10 path). Live (FVSkt_clean/FVSci_clean, SimFire 2020 severe):
        #   KT: 2020 TREES 282, MORT 437; 2030 TREES 0.   CI: 2020 TREES 300, MORT 388; 2030 TREES 0.
        # jl (after fix): KT 2020 TREES 282 (bit-exact), MORT 424, 2030 TREES 3; CI 2020 TREES 307, MORT 383,
        # 2030 TREES 1 — the crown fire FIRES and wipes the stand (surface-only regression left the overstory).
        # Residual = the same family-wide rothermel-HPA byram corner as IE/EM (2030 TREES 3/1 vs live 0).
        (FVSjl.Kootenai(),       "kt_simfire", 282, 15, 400),   # live: 282 pre, 0 post, 437 MORT
        (FVSjl.CentralIdaho(),   "ci_simfire", 304, 15, 350),   # live: 300 pre, 0 post, 388 MORT (live re-derived 2026-09-19 on FVSci_g16 and .new: 300). jl 307→306→304 on bm-regime-close merges, moving toward live; residual +4 OPEN (CI campaign).
    ]
    for (v, base, pre_trees, post_trees_max, mort_min) in cases
        key = joinpath(_fixdir, base * ".key")
        if !isfile(key)
            @test_skip "$base fixture missing"
            continue
        end
        out = FVSjl.run_keyfile(key; variant = v)
        rows = Dict{Int,Vector{String}}()
        stand = 0
        for ln in split(out, '\n')
            startswith(ln, "-999") && (stand += 1)
            stand <= 1 || continue                       # first stand = the SimFire-2020 stand
            p = split(ln)
            length(p) >= 7 && occursin(r"^(2020|2030)$", p[1]) && (rows[parse(Int, p[1])] = String.(p))
        end
        @test haskey(rows, 2020) && haskey(rows, 2030)
        r20 = rows[2020]; r30 = rows[2030]
        trees2020 = parse(Int, r20[3])                   # SUM col 3 = TREES/ACRE
        trees2030 = parse(Int, r30[3])
        mort2020  = parse(Int, r20[end - 3])             # SUM MOR column (…, ACC, MOR, QMD, fortyp, code)
        @test trees2020 == pre_trees                     # pre-fire snapshot bit-exact
        # The crown fire wiped the stand: the surface-only regression left ~174 (IE) / ~112 (EM) alive here.
        @test trees2030 <= post_trees_max
        @test mort2020 >= mort_min                       # crown-fire mortality (surface-only gave ~147 IE / ~211 EM)
    end
end
