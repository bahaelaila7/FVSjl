# =============================================================================
# mortality.jl (teton) — TT mortality (tt/morts.f).
#
# Stand-level (= KT structure): grown-stand DQ10, CONST=SDIMAX/0.02483133, TMD10=CONST·D10^−1.605,
# T85D10=TMD10·PMSDIU, T55D10=TMD10·PMSDIL. SDI self-thinning target TN10 → RN=1−(1−(T−TN10)/T)^(1/FINT).
#   T ≤ T55D10  → TN10=T, RN=0 (below the SDI limit — the emt01 regime; PMSC background dominates).
#   T > T85D10  → TN10=T85D10 (kill to the 85% line).
#   middle (55–85%): iterative linear-fn fit (tt/morts.f label 220) — DEFERRED (emt01 is below-SDI).
# Per-tree, ORIGINAL EM species (1-3,7-10,18): RI=0.5/(1+exp(PMSC+PMD·D+PMDSQ·D²)); RIP=RN if SDI limiting
#   (T>TEM and RN>0) else RI; WKI=P·(1−(1−RIP)^FINT). ADDED species (4-6,11-17,19): KT density Hamilton
#   (deferred — not in emt01). Reuses the shared self-thinning RDPSRT + snag booking.
# =============================================================================

const TT_PMSC  = Float32[6.5112, 6.5112, 7.2985, 5.1677, 9.6943, 5.1677, 5.9617, 9.6943, 5.1677, 0.2118, 5.1677, 5.1677, 5.5877, 5.1677, 5.9617, 5.9617, 5.1677, 5.9617]
const TT_PMD   = Float32[-0.0052485, -0.0052485, -0.0129121, -0.0077681, -0.0127328, -0.0077681, -0.0340128, -0.0127328, -0.0077681, 0.0, -0.0077681, -0.0077681, -0.005348, -0.0077681, -0.0052485, -0.0340128, -0.0077681, -0.0052485]
const TT_PMDSQ = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]  # TT has no PMDSQ

@inline _tt_mort_default(sp::Int) = sp != 10   # only PP(10) has a special morts form (not in ttt01)

# tt/morts.f label-220 iterative linear-fn fit between the 55%/85% SDI lines → TN10 (target tree count at
# D10). IPATH2=true (came from T≤T55D0, T>T55D10) computes the line once (TEM=T) then goes to 230; IPATH2=false
# (55%<T≤85% at DIA0) Newton-iterates TREEIT (≤100) so exp(CEPT+SLP·ln(DIA0)) ≈ T. TN10 = exp(CEPT+SLP·ln(D10)),
# capped at T85D10.
function _tt_tn10_iter(tt::Float32, dia0::Float32, d10::Float32, const_::Float32, pmsdil::Float32,
                       pmsdiu::Float32, t85d10::Float32, t55d0::Float32, ipath2::Bool)::Float32
    treeit = tt + 0.1f0 * tt
    slp = 0f0; cept = 0f0; knt = 1
    while true
        tem = ipath2 ? tt : treeit
        d55m = (log(tem) - log(pmsdil * const_)) / (-1.605f0)
        t55m = log(tem)
        d85m = d55m * 1.25f0
        while true                                    # tt/morts.f label 221: bump D85M until SLP ≤ −0.5
            d85m > 5f0 && (d85m = 5f0); d85m < 0.125f0 && (d85m = 0.125f0)
            t85m = log(const_ * (exp(d85m)^(-1.605f0)) * pmsdiu)
            slp = (t85m - t55m) / (d85m - d55m)
            (slp > -0.5f0 && d85m < 5f0) ? (d85m += 0.1f0) : break
        end
        cept = t55m - slp * d55m
        (ipath2 || tt <= t55d0) && break              # GOTO 230 (no Newton for the IPATH=2 path)
        tprime = cept + slp * log(dia0)
        diff = tt - exp(tprime)
        (diff <= 5f0 && diff >= -5f0) && break
        treeit += 0.5f0 * diff; knt += 1
        knt > 100 && break
    end
    tn10 = exp(cept + slp * log(d10))
    tn10 >= t85d10 && (tn10 = t85d10)
    return tn10
end

function mortality!(s::StandState, ::Teton; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    ba = p.basal_area
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # grown-stand sums (morts.f): T (total tpa), DQ10 (QMD of DBH+DG), AVED (BA-weighted mean DBH)
    tt = 0f0; sd2sq = 0f0; sd0sq = 0f0; dsum = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bark_ratio(bark_a, bark_b, sp, d)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + 2f0 * d * g + g * g); sd0sq += pr * d * d; tt += pr; dsum += d * pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt)          # DQ10 = QMD of DBH+DG (post-growth)
    dq0  = sqrt(sd0sq / tt)          # DQ0  = QMD of DBH (pre-growth) = DIA0
    aved = dsum / tt
    # DIA0<0.3 reset (morts.f 374-376)
    if dq0 < 0.3f0; dq10 = 0.3f0 + dq10 - dq0; dq0 = 0.3f0; end
    # SDI self-thinning boundary (morts.f 455-485)
    sdimax = stand_sdimax(s)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    const_ = sdimax / 0.02483133f0
    tmd10 = const_ * dq10^(-1.605f0); tmd10 > 35000f0 && (tmd10 = 35000f0)
    tmd0  = const_ * dq0^(-1.605f0);  tmd0  > 35000f0 && (tmd0  = 35000f0)
    t85d10 = tmd10 * pmsdiu; t55d10 = tmd10 * pmsdil
    t85d0  = tmd0  * pmsdiu; t55d0  = tmd0  * pmsdil
    # TN10 target (morts.f 200-271): the SDI mature-stand-boundary self-thinning.
    local tn10::Float32
    if tt > t85d0
        tn10 = t85d10                                              # kill to the 85% line
    elseif tt > t55d0
        tn10 = abs(t85d0 - tt) <= 5f0 ? t85d10 :
               _tt_tn10_iter(tt, dq0, dq10, const_, pmsdil, pmsdiu, t85d10, t55d0, false)
    elseif tt <= t55d10
        tn10 = tt                                                  # below 55% at both — hold (RN=0)
    else
        tn10 = _tt_tn10_iter(tt, dq0, dq10, const_, pmsdil, pmsdiu, t85d10, t55d0, true)   # IPATH=2
    end
    tn10 > tt && (tn10 = tt); tn10 < 0.1f0 && (tn10 = 0f0)
    rn = 1f0 - (1f0 - (tt - tn10) / tt)^(1f0 / fint)
    tem = const_ * dq10^(-1.605f0) * pmsdil     # SDI threshold (morts.f 641)
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]
        if _tt_mort_default(sp)
            ri = 0.5f0 * (1f0 / (1f0 + exp(TT_PMSC[sp] + TT_PMD[sp] * d + TT_PMDSQ[sp] * d * d)))
            rip = rn
            (tt <= tem || rn <= 0f0) && (rip = ri)              # background when SDI not yet limiting
            rip > 1f0 && (rip = 1f0)
            wki = pr * (1f0 - (1f0 - rip)^fint)
            wki > pr && (wki = pr)
            sdimax < 5f0 && (wki = pr)
            killed[i] = wki
        else
            error("TT mortality: PP(10) CI-variant morts form not yet ported — not in ttt01")
        end
    end
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
