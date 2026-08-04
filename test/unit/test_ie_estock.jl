# IE AUTOES ESTOCK P(stocking) — bit-exact vs live FVSie (task #143 chunk A1).
# Oracle: FVSie_clean on iet01 stand-4 (SHELTERWOOD WITH AUTO REGENERATION) with a DEBUG keyword
# dumps `PN FOR STOCKING= 0.2116` (estab.f:538) alongside its inputs (IHAB=10, IPREP=1, SLO=0.30,
# ASPECT=5.498, ELEV=34, BAA=1, TIME=1, SQREGT=1). IHAB=10 → IEQ=3 (cedar/hemlock series).
# See docs/AUTOES_CHUNK_PLAN.md for the full measured target set.
using Test
using FVSjl

@testset "IE ESTOCK P(stocking) — task #143 chunk A1" begin
    # measured iet01 stand-4 inputs → oracle PN = 0.2116 (cedar/hemlock series, IEQ=3)
    pn = FVSjl.ie_estock(10, 1, 0.30f0, cos(5.498f0), sin(5.498f0), 34.0f0,
                         1.00f0, log(1.00f0), 1.0f0, 1.0f0, 0.0f0, 0.0f0, 4)
    @test isapprox(pn, 0.2116f0; atol = 1f-3)
    # the caller forms the inventory stocking prob FTEMP = 1/(1+exp(-PN)) = 0.5527 (estab.f:540)
    ftemp = 1f0 / (1f0 + exp(-pn))
    @test isapprox(ftemp, 0.5527f0; atol = 1f-3)

    # IEQ dispatch coverage: IHAB=3→DF(1), 6→GF(2), 10→cedar/hemlock(3), 12→subalpine(4).
    # (regression guard: each branch must return a finite logit for nominal inputs)
    for ihab in (3, 6, 10, 12)
        p = FVSjl.ie_estock(ihab, 1, 0.30f0, cos(5.498f0), sin(5.498f0), 34.0f0,
                            50.0f0, log(50.0f0), 1.0f0, 1.0f0, 0.0f0, 0.0f0, 4)
        @test isfinite(p)
    end
    # IPREP>3 "all roads" branch must also compute
    @test isfinite(FVSjl.ie_estock(10, 4, 0.30f0, 0.7f0, -0.7f0, 34.0f0, 1.0f0, 0.0f0, 1.0f0, 1.0f0, 0.0f0, 0.0f0, 4))
end

@testset "IE ESNSPE P(#species) — task #143 chunk A2a" begin
    # iet01 stand-4 plot-1: ISER=4 (WH), ITPP=2, TPP=2, TPPLN=ln2, BAA=1, ELEV=34, REGT=1, BWAF=0.
    # XCOS=cos(asp)·SLO, XSIN=sin(asp)·SLO (SLO-weighted aspect). Oracle PSPE=(0.543,0.393,0,0,0,0).
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    psp = FVSjl.ie_esnspe(4, 2, 2.0f0, log(2.0f0), 1.0f0, 34.0f0, 1.0f0, 0.0f0, xcos, xsin, slo)
    @test isapprox(psp[1], 0.543f0; atol = 1f-3)
    @test isapprox(psp[2], 0.393f0; atol = 1f-3)
    @test psp[3] == 0f0 && psp[4] == 0f0 && psp[5] == 0f0 && psp[6] == 0f0  # ITPP=2 gates ≥3 off
    # ITPP=6 exercises all six count-logits (regression guard: finite & in (0,1))
    psp6 = FVSjl.ie_esnspe(4, 6, 10.0f0, log(10.0f0), 50.0f0, 34.0f0, 1.0f0, 0.0f0, xcos, xsin, slo)
    @test all(0f0 .< collect(psp6) .< 1f0)
end

@testset "IE ESPADV P(advance species) — task #143 chunk A2b" begin
    # iet01 stand-4 plot-1: IHAB=10, IPREP=1, IFO=4, IPHY≠1, TIME=1 (measured), BAA=1, ELEV=34, REGT=1,
    # BWAF=BWB4=0, occupancy=1 (OCURHT/OCURNF/XESMLT all 1 for sp1-9), OVER<9.95 (no bumps). XCOS/XSIN SLO-weighted.
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    # measured OCURHT(grp10,·)=1.0 for sp1-9, 0 for sp10+ (PP) — the occupancy that zeroes PP.
    occ = Float32[1, 1, 1, 1, 1, 1, 1, 1, 1, 0]; over = zeros(Float32, 10)
    padv = FVSjl.ie_espadv(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over)
    # oracle PADV (ie/espadv.f dump, 3 dp): WP WL DF GF WH RC LP ES AF PP
    oracle = (0.062f0, 0.005f0, 0.048f0, 0.485f0, 0.283f0, 0.122f0, 0.001f0, 0.014f0, 0.039f0, 0.0f0)
    for i in 1:10
        @test isapprox(padv[i], oracle[i]; atol = 6f-4)
    end
    # occupancy zeroes the species out
    occ0 = zeros(Float32, 10)
    padv0 = FVSjl.ie_espadv(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ0, over)
    @test all(collect(padv0) .== 0f0)
end

@testset "IE ESPXCS P(excess species) — task #143 chunk A2b" begin
    # iet01 stand-4 plot-1 (same inputs as ESPADV). occ = OCURHT(grp10)=1 sp1-9, 0 PP.
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    occ = Float32[1, 1, 1, 1, 1, 1, 1, 1, 1, 0]; over = zeros(Float32, 10)
    pxcs = FVSjl.ie_espxcs(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over)
    # oracle PXCS (ie/espxcs.f dump, 3 dp): WP WL DF GF WH RC LP ES AF PP
    oracle = (0.045f0, 0.005f0, 0.080f0, 0.327f0, 0.242f0, 0.194f0, 0.043f0, 0.001f0, 0.025f0, 0.0f0)
    for i in 1:10
        @test isapprox(pxcs[i], oracle[i]; atol = 6f-4)
    end
end

@testset "IE ESRANN establishment LCG — task #143 chunk A2c" begin
    # ie/esrann.f Park-Miller LCG from the iet01 stand-4 seed 43303 (live "RANDOM NUMBER SEED= 43303").
    rng = FVSjl.IEEstabRNG(43303.0)
    d = [FVSjl.ie_esrann!(rng) for _ in 1:52]
    @test isapprox(d[1], 0.33890f0; atol = 1f-4)   # first draw
    @test isapprox(d[2], 0.98084f0; atol = 1f-4)
    # draw #52 = the EMSQR magnitude the live run prints (0.219, estab.f:646-650 uses #51 sign + #52 mag)
    @test isapprox(d[52], 0.21862f0; atol = 1f-4)
    # even seed → odd-adjusted (ESRNSD)
    @test FVSjl.IEEstabRNG(43302.0).ess0 == 43303.0
end

@testset "IE OCURHT habitat occupancy — task #143" begin
    # ie/blkdat.f OCURHT(16,MAXSP). Live debug (iet01 stand-4 grp10): sp1-9=1, sp10+ (incl PP)=0.
    @test [FVSjl.ie_ocurht(10, s) for s in 1:23] == Float32[ones(9); zeros(14)]
    @test FVSjl.ie_ocurht(10, 10) == 0f0        # PP not occupant in grp10 (zeroes PADV/PXCS PP)
    @test FVSjl.ie_ocurht(7, 7) == 1f0          # LP occupies all 16 habitats
    @test FVSjl.ie_ocurht(3, 5) == 0f0          # WH absent in DF-series habitats
    @test all(FVSjl.ie_ocurht(h, sp) == 0f0 for h in 1:16, sp in 11:23)  # added species: no natural regen
end

@testset "IE AUTOES species selection END-TO-END — task #143 chunk A2c" begin
    # Reproduces iet01 stand-4 plot-1's species pick straight from the reseeded RNG, exercising the full
    # draw-order model: WK6-fill(IDUP*NPTIDS=50) + EMSQR(2) + ESTPP(1) + NUMSPE-WK6(6) = 59 draws, then
    # draw #60 = the species-selection WK6(1). Oracle: NUMSPE=1, IBEST -> species 5 (WH), from PADV.
    rng = FVSjl.IEEstabRNG(43303.0)
    d = [FVSjl.ie_esrann!(rng) for _ in 1:60]
    @test isapprox(d[52], 0.21862f0; atol = 1f-4)         # EMSQR magnitude checkpoint (draw #52)
    draw60 = d[60]
    @test isapprox(draw60, 0.61708f0; atol = 1f-4)
    # SUMUP = normalized (PADV+PSUB); PSUB=0 here. PADV sp1-9 then 0.
    padv = Float32[0.062, 0.005, 0.048, 0.485, 0.283, 0.122, 0.001, 0.014, 0.039]
    sumup = [padv ./ sum(padv); zeros(Float32, 14)]       # 23-species selection distribution
    @test FVSjl.ie_estab_pick_species(draw60, sumup) == 5  # WH = oracle IBEST
end

@testset "IE AUTOES multi-plot seed chain — task #143 chunk A2c" begin
    # estab.f per-plot reseed: seed_{N+1}=odd(ESAVEGEN_N). Instrument-confirmed on iet01 stand-4.
    seeds = FVSjl.ie_autoes_plot_seeds(43303, 6)
    @test seeds == [43303, 22913, 17231, 97317, 32953, 75193]
    # per-plot EMSQR from each seed (plot-1 offset by wk6=50; plots 2+ EMSQR at body-draw 2)
    emsqr_oracle = (0.219f0, -0.925f0, -0.528f0, 0.863f0, -0.565f0, 0.721f0)
    for (n, sd) in enumerate(seeds)
        rng = FVSjl.IEEstabRNG(sd)
        off = n == 1 ? 50 : 0          # plot-1 one-time WK6 fill precedes its body
        for _ in 1:off; FVSjl.ie_esrann!(rng); end
        sgn = FVSjl.ie_esrann!(rng); mag = FVSjl.ie_esrann!(rng)   # EMSQR sign, magnitude (body draws 1,2)
        emsqr = (sgn < 0.5f0 ? -mag : mag)
        @test isapprox(emsqr, emsqr_oracle[n]; atol = 1f-3)
    end
end

@testset "IE AUTOES per-plot ITPP — task #143 chunk A2c" begin
    # ITPP = INT(ESTPP(body-draw-3)+0.5), capped at MAXING(IHAB) for ingrowth. IHAB=10 → MAXING=7.
    # Oracle iet01 stand-4 plots 1-8: TREES/PLOT = [2,1,7,1,3,1,7,3] (plots 3,7 hit the MAXING=7 cap).
    seeds = FVSjl.ie_autoes_plot_seeds(43303, 8)
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    maxing = FVSjl._IE_MAXING[10]
    @test maxing == 7
    jl = Int[]
    for (n, sd) in enumerate(seeds)
        rng = FVSjl.IEEstabRNG(sd)
        for _ in 1:(n == 1 ? 50 : 0); FVSjl.ie_esrann!(rng); end
        FVSjl.ie_esrann!(rng); FVSjl.ie_esrann!(rng)            # EMSQR sign, mag
        d3 = FVSjl.ie_esrann!(rng)                               # ESTPP draw (body-3)
        push!(jl, clamp(round(Int, FVSjl.ie_estpp(d3, 10, xcos, xsin, slo, 1.0f0, 0.0f0)), 1, maxing))
    end
    @test jl == [2, 1, 7, 1, 3, 1, 7, 3]
end

@testset "IE AUTOES per-plot species selection — task #143 chunk A2c" begin
    # Full per-plot selection from the seed chain: seed → EMSQR → ITPP → NUMSPE → species (IBEST).
    # Oracle iet01 stand-4 plots 1-7 IBEST: [5],[5],[5],[4],[4,6],[4],[8] (plot-5 is NUMSPE=2 GF+RC).
    seeds = FVSjl.ie_autoes_plot_seeds(43303, 7)
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10)
    padv = collect(FVSjl.ie_espadv(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over))
    sumup_base = [Float32.(padv); zeros(Float32, 13)]; sumup_base ./= sum(sumup_base)
    nspnz = count(>(1f-4), sumup_base)
    maxspp = FVSjl._IE_MAXSPP[10]; maxing = FVSjl._IE_MAXING[10]
    oracle = [[5], [5], [5], [4], [4, 6], [4], [8]]
    for (n, sd) in enumerate(seeds)
        rng = FVSjl.IEEstabRNG(sd)
        for _ in 1:(n == 1 ? 50 : 0); FVSjl.ie_esrann!(rng); end
        FVSjl.ie_esrann!(rng); FVSjl.ie_esrann!(rng)
        itpp = clamp(round(Int, FVSjl.ie_estpp(FVSjl.ie_esrann!(rng), 10, xcos, xsin, slo, 1.0f0, 0.0f0)), 1, maxing)
        wk6n = [FVSjl.ie_esrann!(rng) for _ in 1:6]
        wk6s = [FVSjl.ie_esrann!(rng) for _ in 1:6]
        numspe = 1
        if itpp != 1
            pspe = collect(FVSjl.ie_esnspe(4, itpp, Float32(itpp), log(Float32(itpp)), 1.0f0, 34.0f0, 1.0f0, 0.0f0, xcos, xsin, slo))
            cum = cumsum(pspe ./ sum(pspe)); numspe = 6
            for i in 1:5
                if wk6n[i] <= cum[i]; numspe = i; break; end
            end
        end
        numspe = min(numspe, maxspp, nspnz)
        su = copy(sumup_base); picked = Int[]
        for i in 1:numspe
            j = FVSjl.ie_estab_pick_species(wk6s[i], su); push!(picked, j); su[j] = 0f0
            t = sum(su); t > 0 && (su ./= t)
        end
        @test sort(picked) == oracle[n]
    end
end

@testset "IE AUTOES total ingrowth — task #143 chunk A2c" begin
    # Full tally over NCOUNT = NPTIDS·IDUP = 50 plots. Each plot: ITPP trees × PROB1·300/DUPNPT TPA.
    # Measured: PROB1=0.5527 (constant, IHAB=10), DUPNPT=50. Oracle 1999 total ingrowth = 583.7 TPA.
    seeds = FVSjl.ie_autoes_plot_seeds(43303, 50)
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    maxing = FVSjl._IE_MAXING[10]
    sitpp = 0
    for (n, sd) in enumerate(seeds)
        rng = FVSjl.IEEstabRNG(sd)
        for _ in 1:(n == 1 ? 50 : 0); FVSjl.ie_esrann!(rng); end
        FVSjl.ie_esrann!(rng); FVSjl.ie_esrann!(rng)
        d3 = FVSjl.ie_esrann!(rng)
        sitpp += clamp(round(Int, FVSjl.ie_estpp(d3, 10, xcos, xsin, slo, 1.0f0, 0.0f0)), 1, maxing)
    end
    @test sitpp == 176
    total = sitpp * 0.5527f0 * 300f0 / 50f0
    @test isapprox(total, 583.7f0; atol = 0.5f0)   # oracle ingrowth total, bit-exact to print precision
end

@testset "IE AUTOES full per-species tally — task #143 chunk A2c" begin
    # ie_autoes_tally composes the whole tally: seed chain → per-plot selection → excess (PXCS×IBEST) → TPA.
    # Oracle iet01 stand-4 (1999 ingrowth): WP33 WL20 DF7 GF202 WH222 RC50 ES23 AF27, total 583.7.
    occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10); slo = 0.30f0
    t = FVSjl.ie_autoes_tally(seed0 = 43303, nplots = 50, ihab = 10, iser = 4, ifo = 4, iprep = 1, iphy = 3,
        xcos = cos(5.498f0) * slo, xsin = sin(5.498f0) * slo, slo = slo, elev = 34.0f0, baa = 1.0f0,
        regt = 1.0f0, bwaf = 0.0f0, bwb4 = 0.0f0, prob1 = 0.5527f0, dupnpt = 50.0f0, occ = occ, over = over)
    oracle = [33, 20, 7, 202, 222, 50, 0, 23, 27]   # WP WL DF GF WH RC LP ES AF
    for i in 1:9
        @test round(Int, t[i]) == oracle[i]
    end
    @test isapprox(sum(t), 583.7; atol = 0.6)
end

@testset "IE AUTOES seed0 from ESSS + self-contained tally — task #143 chunk A2c" begin
    # The establishment RNG default ESSS=55329 (esblkd.f:30); the first draw derives ESDRAW=seed0.
    @test FVSjl.ie_autoes_seed0(55329) == 43303
    # end-to-end from the DERIVED seed0: species split must still be bit-exact (no hardcoded 43303).
    occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10); slo = 0.30f0
    t = FVSjl.ie_autoes_tally(seed0 = FVSjl.ie_autoes_seed0(55329), nplots = 50, ihab = 10, iser = 4,
        ifo = 4, iprep = 1, iphy = 3, xcos = cos(5.498f0) * slo, xsin = sin(5.498f0) * slo, slo = slo,
        elev = 34.0f0, baa = 1.0f0, regt = 1.0f0, bwaf = 0.0f0, bwb4 = 0.0f0, prob1 = 0.5527f0,
        dupnpt = 50.0f0, occ = occ, over = over)
    @test round(Int, t[4]) == 202 && round(Int, t[5]) == 222   # GF, WH
    @test isapprox(sum(t), 583.7; atol = 0.6)
end

@testset "IE ESADVH advance-regen height — task #143 chunk A2c" begin
    # Bit-exact vs live FVSie (iet01 stand-4 ADVHHT dump). AUTOES advance: AGEL=0 (AGE clamped), BNORM=1.0;
    # DILATE=FIRST(1,sp) order-statistic chain (0.1 → √0.1 → √√0.1 …). Aspect SLO-weighted.
    slo = 0.30f0; xc = cos(5.498f0) * slo; xs = sin(5.498f0) * slo
    kw = (baa = 1.0f0, elev = 34.0f0, xcos = xc, xsin = xs, slo = slo, ihtser = 4, iphy = 3, iprep = 1)
    d1 = 0.1f0; d2 = sqrt(d1); d3 = sqrt(d2)
    # (sp, emsqr, dilate, oracle HHT)
    cases = [(5, 0.21862f0, d1, 0.654103f0),   # WH plot1
             (5, -0.92534f0, d2, 0.550567f0),  # WH plot2
             (5, -0.52763f0, d3, 0.549335f0),  # WH plot3
             (6, -0.56473f0, d1, 0.300050f0),  # RC plot4
             (4, 0.72090f0, d1, 0.176656f0),   # GF plot5
             (4, 0.56212f0, d2, 0.188457f0)]   # GF plot6
    for (sp, em, dil, orc) in cases
        h = FVSjl.ie_esadvh(sp, em, dil, 0.0f0, 1.0f0; kw...)
        @test isapprox(h, orc; atol = 3f-4)
    end
end

@testset "IE ESDLAY germination delay — task #143 chunk A2c" begin
    # Advance (ias=1): DELAY=(Weibull+3)*(-1) is always negative → clamps to 0. Matches live
    # "DELAY TO GERM=0.0000" for iet01 advance regen (BAA=1, BWB4=0). Holds ∀ species, draw.
    @test all(FVSjl.ie_esdlay(sp, 1, dr, 1.0f0, 1.0f0) == 0f0 for sp in 1:10, dr in (0.1f0, 0.5f0, 0.9f0))
    # Subsequent (ias=2): DELAY=Weibull-4, clamped [0,10] — monotone in draw, in range.
    @test 0f0 <= FVSjl.ie_esdlay(4, 2, 0.5f0, 10.0f0, 1.0f0) <= 10f0
end

@testset "IE ESPSUB P(subsequent species) — task #143 chunk A2c" begin
    # Populated at estab.f:774 for the ADV/SUBS ICHOI dispatch. Oracle iet01 stand-4 plot-1 (TIME=1, BAALN=0).
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10)
    ps = FVSjl.ie_espsub(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 0.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over)
    oracle = (0.0628982f0, 0.0310127f0, 0.0701199f0, 0.232294f0, 0.0650079f0,
              0.0100316f0, 0.0276168f0, 0.00603530f0, 0.0104509f0, 0.0f0)
    for i in 1:10
        @test isapprox(ps[i], oracle[i]; atol = 6f-4)
    end
end

@testset "IE AUTOES ICHOI adv/subs dispatch — task #143 chunk A2c" begin
    # Best tree is ADVANCE (→ESADVH) if ADV/SUBS draw ≤ PADV/(PADV+PSUB), else SUBSEQUENT (→ESSUBH).
    # ADV/SUBS draws at body-16..38 (one per species). Oracle ADVHHT dump = advance trees only, in order:
    # WH,WH,WH (plots1-3), RC (plot5), GF (plot6), GF (plot8); plots 4(GF),7(ES) route SUBSEQUENT.
    seeds = FVSjl.ie_autoes_plot_seeds(43303, 8)
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10)
    padv = collect(FVSjl.ie_espadv(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over))
    psub = collect(FVSjl.ie_espsub(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 0.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over))
    sumup_base = [Float32.(padv); zeros(Float32, 13)]; sumup_base ./= sum(sumup_base)
    nspnz = count(>(1f-4), sumup_base); maxspp = FVSjl._IE_MAXSPP[10]; maxing = FVSjl._IE_MAXING[10]
    adv_seq = Int[]                       # advance best species in plot order
    for (n, sd) in enumerate(seeds)
        rng = FVSjl.IEEstabRNG(sd)
        for _ in 1:(n == 1 ? 50 : 0); FVSjl.ie_esrann!(rng); end
        FVSjl.ie_esrann!(rng); FVSjl.ie_esrann!(rng)
        itpp = clamp(round(Int, FVSjl.ie_estpp(FVSjl.ie_esrann!(rng), 10, xcos, xsin, slo, 1.0f0, 0.0f0)), 1, maxing)
        wk6n = [FVSjl.ie_esrann!(rng) for _ in 1:6]; wk6s = [FVSjl.ie_esrann!(rng) for _ in 1:6]
        numspe = 1
        if itpp != 1
            pspe = collect(FVSjl.ie_esnspe(4, itpp, Float32(itpp), log(Float32(itpp)), 1.0f0, 34.0f0, 1.0f0, 0.0f0, xcos, xsin, slo))
            cum = cumsum(pspe ./ sum(pspe)); numspe = 6
            for i in 1:5; wk6n[i] <= cum[i] && (numspe = i; break); end
        end
        numspe = min(numspe, maxspp, nspnz)
        su = copy(sumup_base); best = Int[]
        for i in 1:numspe; j = FVSjl.ie_estab_pick_species(wk6s[i], su); push!(best, j); su[j] = 0f0; t = sum(su); t > 0 && (su ./= t); end
        advsub = [FVSjl.ie_esrann!(rng) for _ in 1:23]
        for j in best
            advsub[j] <= padv[j] / (padv[j] + psub[j]) && push!(adv_seq, j)
        end
    end
    # first 6 advance trees = the captured ADVHHT dump order (plot-8's 2nd advance tree AF follows, not captured).
    @test adv_seq[1:6] == [5, 5, 5, 6, 4, 4]   # WH WH WH RC GF GF
end

@testset "IE AUTOES advance-height integration — task #143 chunk A2c" begin
    # End-to-end: seed chain → selection → ICHOI → ie_esadvh with the FIRST(1,sp) order-statistic chain
    # (0.1 → √0.1 → …) tracked per species. Advance DELAY=0 so heights need only EMSQR + DILATE (no draw idx).
    # Oracle ADVHHT dump: [0.654103, 0.550567, 0.549335, 0.300050, 0.176656, 0.188457].
    seeds = FVSjl.ie_autoes_plot_seeds(43303, 8)
    slo = 0.30f0; xc = cos(5.498f0) * slo; xs = sin(5.498f0) * slo
    occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10)
    padv = collect(FVSjl.ie_espadv(10, 1, 4, 3, xc, xs, slo, 1.0f0, 1.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over))
    psub = collect(FVSjl.ie_espsub(10, 1, 4, 3, xc, xs, slo, 1.0f0, 1.0f0, 0.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over))
    sumup_base = [Float32.(padv); zeros(Float32, 13)]; sumup_base ./= sum(sumup_base)
    nspnz = count(>(1f-4), sumup_base); maxspp = FVSjl._IE_MAXSPP[10]; maxing = FVSjl._IE_MAXING[10]
    kw = (baa = 1.0f0, elev = 34.0f0, xcos = xc, xsin = xs, slo = slo, ihtser = 4, iphy = 3, iprep = 1)
    FIRST = fill(0.1f0, 23); advh = Float32[]
    for (n, sd) in enumerate(seeds)
        rng = FVSjl.IEEstabRNG(sd)
        for _ in 1:(n == 1 ? 50 : 0); FVSjl.ie_esrann!(rng); end
        sgn = FVSjl.ie_esrann!(rng); mag = FVSjl.ie_esrann!(rng); emsqr = sgn < 0.5f0 ? -mag : mag
        itpp = clamp(round(Int, FVSjl.ie_estpp(FVSjl.ie_esrann!(rng), 10, xc, xs, slo, 1.0f0, 0.0f0)), 1, maxing)
        wk6n = [FVSjl.ie_esrann!(rng) for _ in 1:6]; wk6s = [FVSjl.ie_esrann!(rng) for _ in 1:6]
        numspe = 1
        if itpp != 1
            pspe = collect(FVSjl.ie_esnspe(4, itpp, Float32(itpp), log(Float32(itpp)), 1.0f0, 34.0f0, 1.0f0, 0.0f0, xc, xs, slo))
            cum = cumsum(pspe ./ sum(pspe)); numspe = 6
            for i in 1:5; wk6n[i] <= cum[i] && (numspe = i; break); end
        end
        numspe = min(numspe, maxspp, nspnz)
        su = copy(sumup_base); best = Int[]
        for i in 1:numspe; j = FVSjl.ie_estab_pick_species(wk6s[i], su); push!(best, j); su[j] = 0f0; t = sum(su); t > 0 && (su ./= t); end
        advsub = [FVSjl.ie_esrann!(rng) for _ in 1:23]
        for j in best
            if advsub[j] <= padv[j] / (padv[j] + psub[j])
                push!(advh, FVSjl.ie_esadvh(j, emsqr, FIRST[j], 0.0f0, 1.0f0; kw...)); FIRST[j] = sqrt(FIRST[j])
            end
        end
    end
    oracle = [0.654103f0, 0.550567f0, 0.549335f0, 0.300050f0, 0.176656f0, 0.188457f0]
    for i in 1:6
        @test isapprox(advh[i], oracle[i]; atol = 3f-4)
    end
end

@testset "IE AUTOES subsequent + excess heights — task #143 chunk A2c" begin
    # Complete the height model. AUTOES first tally: TIME=1 → essubh AGE=max(1-DELAY-GENTIM,1)=1 (AGEL=0),
    # BNORM=BNORML(1)=1.0; disp = EMSQR·DILATE·BNORM. Excess: ie_esxcsh(sp, TALL, XMIN, TIME=1, draw).
    slo = 0.30f0; xc = cos(5.498f0) * slo; xs = sin(5.498f0) * slo
    # Subsequent best trees (SUBH dump): plot-4 GF EMSQR=0.862929 DILATE=0.1 → 0.143546; plot-7 ES EMSQR=-0.279643 → 0.085346
    @test isapprox(FVSjl.ie_essubh(4, 1.0f0, 1.0f0, 4, 1, 3, xc, xs, slo, 34.0f0, 0.862929f0 * 0.1f0), 0.143546f0; atol = 1f-3)
    @test isapprox(FVSjl.ie_essubh(8, 1.0f0, 1.0f0, 4, 1, 3, xc, xs, slo, 34.0f0, -0.279643f0 * 0.1f0), 0.085346f0; atol = 1f-3)
    # Excess trees (XCSH dump): WH TALL=0.70, XMIN=0.5, DRAW=0.670273 → 0.584198
    @test isapprox(FVSjl.ie_esxcsh(5, 0.70f0, 0.5f0, 1.0f0, 0.670273f0), 0.584198f0; atol = 3f-4)
end

@testset "IE ESTPP trees-per-plot — task #143 chunk A2c" begin
    # draw #53 (after WK6-fill 50 + EMSQR 2) drives ESTPP. Oracle iet01 stand-4 plot-1: TREES/PLOT = ITPP = 2.
    rng = FVSjl.IEEstabRNG(43303.0)
    for _ in 1:52; FVSjl.ie_esrann!(rng); end
    draw53 = FVSjl.ie_esrann!(rng)
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    tpp = FVSjl.ie_estpp(draw53, 10, xcos, xsin, slo, 1.0f0, 0.0f0)
    itpp = clamp(round(Int, tpp), 1, 99)   # INT(TPP+0.5), MAXTPP clamp
    @test itpp == 2   # oracle TREES/PLOT
end
