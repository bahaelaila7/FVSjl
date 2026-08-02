# =============================================================================
# mortality.jl (bluemountains) — BM mortality (bm/morts.f). Chunk 7.
#
# Uniform Hamilton background (RI = 0.5/(1+exp(PMSC+PMD·D))) + SDI self-thinning (RN), = the UT/EM
# form. RIP = RN when SDI-limiting (T>TEM and RN>0) else RI; WKI = P·(1−(1−RIP)^FINT). The X
# diameter-class modifier (bm/morts.f, FIXMORT) defaults to 1 (no FIXMORT keyword in bmt01).
# =============================================================================

# bm/morts.f DATA PMSC / PMD (species order WP WL DF GF MH WJ LP ES AF PP WB LM PY YC AS CW OS OH).
const BM_PMSC = Float32[6.5112, 6.5112, 7.2985, 5.1677, 9.6943, 5.1677, 5.9617, 9.6943, 5.1677,
                        5.5877, 6.5112, 6.5112, 5.5877, 5.5877, 5.1677, 5.5877, 5.5877, 5.9617]
const BM_PMD  = Float32[-0.0052485, -0.0052485, -0.0129121, -0.0077681, -0.0127328, -0.0077681,
                        -0.0340128, -0.0127328, -0.0077681, -0.005348, -0.0052485, -0.0052485,
                        -0.005348, -0.005348, -0.0077681, -0.005348, -0.005348, -0.0340128]

function mortality!(s::StandState, ::BlueMountains; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    tt = 0f0; sd2sq = 0f0; sd0sq = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bm_bratio(sd, sp, d)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + 2f0 * d * g + g * g); sd0sq += pr * d * d; tt += pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt); dq0 = sqrt(sd0sq / tt)
    if dq0 < 0.3f0; dq10 = 0.3f0 + dq10 - dq0; dq0 = 0.3f0; end
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
        ri = 0.5f0 * (1f0 / (1f0 + exp(BM_PMSC[sp] + BM_PMD[sp] * d)))
        rip = rn
        (tt <= tem || rn <= 0f0) && (rip = ri)
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
