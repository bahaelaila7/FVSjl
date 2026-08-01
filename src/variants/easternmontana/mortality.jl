# =============================================================================
# mortality.jl (easternmontana) — EM mortality (em/morts.f).
#
# Stand-level (= KT structure): grown-stand DQ10, CONST=SDIMAX/0.02483133, TMD10=CONST·D10^−1.605,
# T85D10=TMD10·PMSDIU, T55D10=TMD10·PMSDIL. SDI self-thinning target TN10 → RN=1−(1−(T−TN10)/T)^(1/FINT).
#   T ≤ T55D10  → TN10=T, RN=0 (below the SDI limit — the emt01 regime; PMSC background dominates).
#   T > T85D10  → TN10=T85D10 (kill to the 85% line).
#   middle (55–85%): iterative linear-fn fit (em/morts.f label 220) — DEFERRED (emt01 is below-SDI).
# Per-tree, ORIGINAL EM species (1-3,7-10,18): RI=0.5/(1+exp(PMSC+PMD·D+PMDSQ·D²)); RIP=RN if SDI limiting
#   (T>TEM and RN>0) else RI; WKI=P·(1−(1−RIP)^FINT). ADDED species (4-6,11-17,19): KT density Hamilton
#   (deferred — not in emt01). Reuses the shared self-thinning RDPSRT + snag booking.
# =============================================================================

const EM_PMSC  = Float32[5.45676, 5.26043, 5.55086, 0.0, 0.2118, 0.0, 3.87794, 6.41265, 5.88697, 5.58766,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.47709, 0.0]
const EM_PMD   = Float32[-0.0118233, -0.0097092, -0.0129121, 0.0, 0.0, 0.0, 0.3078, -0.0127328, -0.0333752,
                         -0.0052485, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0395156, 0.0]
const EM_PMDSQ = Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0174, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

@inline _em_orig_species(sp::Int) = sp <= 3 || (7 <= sp <= 10) || sp == 18

function mortality!(s::StandState, ::EasternMontana; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    ba = p.basal_area
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # grown-stand sums (morts.f): T (total tpa), DQ10 (QMD of DBH+DG), AVED (BA-weighted mean DBH)
    tt = 0f0; sd2sq = 0f0; dsum = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bark_ratio(bark_a, bark_b, sp, d)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + 2f0 * d * g + g * g); tt += pr; dsum += d * pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt)
    aved = dsum / tt
    # SDI self-thinning boundary (morts.f 455-485)
    sdimax = stand_sdimax(s)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    const_ = sdimax / 0.02483133f0
    tmd10 = const_ * dq10^(-1.605f0); tmd10 > 35000f0 && (tmd10 = 35000f0)
    t85d10 = tmd10 * pmsdiu; t55d10 = tmd10 * pmsdil
    # TN10 target (morts.f 502-590). emt01 is below-SDI (T ≤ T55D10) ⇒ TN10=T, RN=0.
    local tn10::Float32
    if tt <= t55d10
        tn10 = tt
    elseif tt > t85d10
        tn10 = t85d10
    else
        error("EM mortality: stand between 55%–85% SDI (iterative TN10 fit) not yet ported (emt01 is below-SDI)")
    end
    tn10 > tt && (tn10 = tt); tn10 < 0.1f0 && (tn10 = 0f0)
    rn = 1f0 - (1f0 - (tt - tn10) / tt)^(1f0 / fint)
    tem = const_ * dq10^(-1.605f0) * pmsdil     # SDI threshold (morts.f 641)
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]
        if _em_orig_species(sp)
            ri = 0.5f0 * (1f0 / (1f0 + exp(EM_PMSC[sp] + EM_PMD[sp] * d + EM_PMDSQ[sp] * d * d)))
            rip = rn
            (tt <= tem || rn <= 0f0) && (rip = ri)              # background when SDI not yet limiting
            rip > 1f0 && (rip = 1f0)
            wki = pr * (1f0 - (1f0 - rip)^fint)
            wki > pr && (wki = pr)
            sdimax < 5f0 && (wki = pr)
            killed[i] = wki
        else
            error("EM mortality: added-species (sp $sp) KT-density Hamilton path not yet ported")
        end
    end
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
