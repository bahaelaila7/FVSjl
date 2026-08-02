# =============================================================================
# mortality.jl (utah) — UT periodic mortality (ut/morts.f). Chunk 7.
#
# UT mortality = EM's ORIGINAL-species path applied UNIFORMLY to all 24 species (no added-species
# Hamilton branch, no PMDSQ): background RI = 0.5·(1/(1+exp(PMSC[sp]+PMD[sp]·D))), overridden by the
# SDI self-thinning rate RN when the stand is above the mature-stand boundary (T>TEM & RN>0). The
# TN10 self-thinning (SDIMAX/PMSDIU/PMSDIL) is identical to EM (reuse _em_tn10_iter/stand_sdimax);
# UT is a Zeide-SDI variant (LZEIDE=T, set in species.jl), so stand_sdimax takes the Zeide path.
#   WKI = P·(1−(1−RIP)^FINT); SDIMAX<5 ⇒ WKI=P (climate-shift override, morts.f:331).
# =============================================================================

const UT_PMSC = Float32[6.5112, 6.5112, 7.2985, 5.1677, 9.6943, 5.1677, 5.9617, 9.6943, 5.1677, 5.5877,
                        5.1677, 5.1677, 5.1677, 5.1677, 5.1677, 5.1677, 5.1677, 5.9617, 5.9617, 5.9617,
                        5.5877, 5.9617, 5.1677, 5.1677]
const UT_PMD  = Float32[-0.0052485, -0.0052485, -0.0129121, -0.0077681, -0.0127328, -0.0077681, -0.0340128,
                        -0.0127328, -0.0077681, -0.005348, -0.0077681, -0.0077681, -0.0077681, -0.0077681,
                        -0.0077681, -0.0077681, -0.0077681, -0.0052485, -0.0052485, -0.0340128, -0.005348,
                        -0.0052485, -0.0077681, -0.0077681]

function mortality!(s::StandState, ::Utah; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # grown-stand sums (morts.f): T (total tpa), DQ10 (QMD of DBH+DG), DQ0 (QMD of DBH).
    tt = 0f0; sd2sq = 0f0; sd0sq = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bark_ratio(bark_a, bark_b, sp, d)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + 2f0 * d * g + g * g); sd0sq += pr * d * d; tt += pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt); dq0 = sqrt(sd0sq / tt)
    if dq0 < 0.3f0; dq10 = 0.3f0 + dq10 - dq0; dq0 = 0.3f0; end
    # SDI self-thinning boundary (morts.f 320-485) — Zeide SDIMAX for UT.
    sdimax = stand_sdimax(s)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    const_ = sdimax / 0.02483133f0
    tmd10 = const_ * dq10^(-1.605f0); tmd10 > 35000f0 && (tmd10 = 35000f0)
    tmd0  = const_ * dq0^(-1.605f0);  tmd0  > 35000f0 && (tmd0  = 35000f0)
    t85d10 = tmd10 * pmsdiu; t55d10 = tmd10 * pmsdil
    t85d0  = tmd0  * pmsdiu; t55d0  = tmd0  * pmsdil
    local tn10::Float32
    if tt > t85d0
        tn10 = t85d10
    elseif tt > t55d0
        tn10 = abs(t85d0 - tt) <= 5f0 ? t85d10 :
               _em_tn10_iter(tt, dq0, dq10, const_, pmsdil, pmsdiu, t85d10, t55d0, false)
    elseif tt <= t55d10
        tn10 = tt
    else
        tn10 = _em_tn10_iter(tt, dq0, dq10, const_, pmsdil, pmsdiu, t85d10, t55d0, true)
    end
    tn10 > tt && (tn10 = tt); tn10 < 0.1f0 && (tn10 = 0f0)
    rn = 1f0 - (1f0 - (tt - tn10) / tt)^(1f0 / fint)
    tem = const_ * dq10^(-1.605f0) * pmsdil
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]
        ri = 0.5f0 * (1f0 / (1f0 + exp(UT_PMSC[sp] + UT_PMD[sp] * d)))
        rip = rn
        (tt <= tem || rn <= 0f0) && (rip = ri)            # background when SDI not yet limiting
        rip > 1f0 && (rip = 1f0)
        wki = pr * (1f0 - (1f0 - rip)^fint)
        wki > pr && (wki = pr)
        sdimax < 5f0 && (wki = pr)
        killed[i] = wki
    end
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end

# ut/crown.f crown-ratio DCR change model (Weibull, >3" DBH) — PORT PENDING (chunk 5 remainder).
# STUB: no-op for now. Inert for the utt01_s1 single-cycle .sum (only affects cycle 2+ crowns);
# multi-cycle validation needs the real DCR port (mirror EM crown_ratio_update!, UT coefficients).
crown_ratio_update!(s::StandState, ::Utah; fint::Float32 = 10.0f0, lstart::Bool = false, kwargs...) = s
