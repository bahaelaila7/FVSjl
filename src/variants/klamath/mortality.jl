# =============================================================================
# mortality.jl (klamath) — NC periodic mortality (nc/morts.f + nc/nwcmrt.f). Chunk 7.
# Background RI = 0.5·(1/(1+exp(PMSC+PMD·D))) (RW/sp12 floored at 0.0001 before the halve, morts.f:480-488),
# overridden by the Zeide-SDI self-thinning rate above the mature-stand boundary. NC is LZEIDE=.TRUE.
# PMSC/PMD sp1-11 are IDENTICAL to EM/UT; sp12 RW (redwood) is special (PMSC=2.5968, PMD=+0.51261).
#
# UNLIKE EM/UT/TT, NC DISTRIBUTES the density (or background) mortality total across records by
# **NWCMRT** (nc/nwcmrt.f) — a percentile/shade-tolerance geometric progression that concentrates kills on
# suppressed (low-BA-percentile) trees — and then runs the morts.f QMD-convergence loop (morts.f:551-598,
# `GO TO 10`): NWCMRT drives QMD up (small trees die) ⇒ TN10 (and the BAMAX target) tighten ⇒ re-iterate.
# The old port applied a UNIFORM per-tree rate, which matches the TOTAL first-pass kill but NOT the
# distribution; combined with the BAMAX residual-BA cap that produced a ~3× TPA UNDER-kill on mature
# over-max-SDI stands (live drives QMD 25→54 killing small trees; uniform kept them). NWCMRT is structurally
# identical to SN's VARMRT (only VARADJ + the 0.1 efficiency factor differ), so it reuses the shared
# `_varmrt!` driver with an NC-specific `_varmrt_efftr!`.
# =============================================================================

const NC_PMSC = Float32[6.5112, 6.5112, 7.2985, 5.1677, 9.6943, 5.1677, 5.9617, 9.6943, 5.1677, 5.5877, 5.1677, 2.59680]
const NC_PMD  = Float32[-0.0052485, -0.0052485, -0.0129121, -0.0077681, -0.0127328, -0.0077681, -0.0340128,
                        -0.0127328, -0.0077681, -0.005348, -0.0077681, 0.51261]

# NWCMRT tree shade tolerance (nc/nwcmrt.f:49-51; species 1=OC 2=SP 3=DF 4=WF 5=M 6=IC 7=BO 8=TO 9=RF
# 10=PP 11=OH 12=RW). 1.0 = most intolerant (highest self-thin efficiency); low = shade-tolerant survivor.
const NC_VARADJ = Float32[0.65, 0.70, 0.65, 0.55, 0.80, 0.60, 1.00, 0.55, 0.50, 0.85, 0.55, 0.80]

# NWCMRT per-record mortality efficiency (nc/nwcmrt.f:97-104): EFFTR = PEFF·VARADJ·0.1, PEFF a cubic of the
# BA percentile PCT (t.crown_ratio, filled by stand_pct!). PCT**3.0 is a gfortran powf ⇒ fpow (doctrine #8).
function _varmrt_efftr!(efftr, s, ::Klamath, t::TreeList, n::Int)
    pct = t.crown_ratio; tpa = t.tpa; sp = t.species
    pass1 = 0f0
    @inbounds for i in 1:n
        pe = clamp(0.84525f0 - 0.01074f0 * pct[i] + 0.0000002f0 * fpow(pct[i], 3f0), 0.01f0, 1f0)
        efftr[i] = pe * NC_VARADJ[sp[i]] * 0.1f0
        pass1 += tpa[i] * efftr[i]
    end
    return pass1
end

function mortality!(s::StandState, ::Klamath; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    yr = htg_period(s.variant)          # NC growth model is 5-yr native (YR=5); the cycle is FINT=10
    # Zeide grown-stand sums (morts.f:186-232, LZEIDE path): T and the Reineke diameters DR0/DR10.
    # G is the OUTSIDE-bark end-of-cycle increment `(DG/BARK)·(FINT/5)` (morts.f:199) — for a 10-yr cycle
    # over a 5-yr-native growth model this is the linear-FINT extrapolation, NOT the raw diam_growth
    # (EM/UT/TT get away with raw because their FINT==YR==10; NC's FINT=10≠YR=5, so the /5 factor matters).
    tt = 0f0; sumdr10 = 0f0; sumdr0 = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bark_ratio(bark_a, bark_b, sp, d)
        g = _mort_traj_g(t.diam_growth[i], d, bark, fint, yr)
        sumdr10 += pr * fpow(d + g, 1.605f0); sumdr0 += pr * fpow(d, 1.605f0); tt += pr
    end
    tt < 1f-6 && return s
    dr10 = fpow(sumdr10 / tt, 1f0 / 1.605f0); dia0 = fpow(sumdr0 / tt, 1f0 / 1.605f0)
    if dia0 < 0.3f0; dr10 = 0.3f0 + dr10 - dia0; dia0 = 0.3f0; end
    sdimax = stand_sdimax(s)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    const_ = sdimax / 0.02483133f0
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    efftr  = @view s.scratch.mort_efftr[1:n]
    temwk2 = @view s.scratch.mort_temwk2[1:n]

    # Background (Hamilton) mortality TOTAL — the target NWCMRT distributes when the SDI self-thin is not yet
    # limiting (morts.f DO 50 loop with RIP=RI, then nwcmrt.f:73-78 TOKILL=ΣWK2). Depends only on start-of-
    # cycle TPA, so compute once. RW (sp12) RI floored at 0.0001 before the halve; MORTMULT (X) hits only the
    # background rate (morts.f:504-508: X=XMORT in [D1,D2], X=1 when the density rate is in effect).
    cur_year = current_cycle_year(s)
    bg_tokill = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; pr <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        ri = 1f0 / (1f0 + fexp(NC_PMSC[sp] + NC_PMD[sp] * d))
        sp == 12 && ri < 0.0001f0 && (ri = 0.0001f0)             # RW floor (morts.f:480-484)
        ri = 0.5f0 * ri                                          # morts.f:488
        x = active_mort_mult(s.control, sp, cur_year, d)         # MORTMULT (background only)
        bg_tokill += min(pr * (1f0 - fpow(1f0 - ri, fint)) * x, pr)
    end

    if sdimax < 5f0
        # Climate shifted enough that the site can no longer support trees — kill everything (morts.f:316-319).
        @inbounds for i in 1:n; killed[i] = t.tpa[i]; end
    else
        # MORTS QMD-convergence iteration (morts.f:242-598, label 10 … GO TO 10): for the assumed end-of-cycle
        # Reineke diameter d10, solve TN10 (self-thin target), distribute the excess T−TN10 (or the background
        # total) by NWCMRT, recompute the post-mortality d10n, and re-iterate with d10=d10n until it converges
        # (|d10−d10n|≤0.1) or the QMD would fall below dia0. Up to 10 passes (morts.f:583 IPASS.EQ.10).
        d10 = dr10
        @inbounds for _ in 1:10
            tmd10 = const_ * fpow(d10, -1.605f0); tmd10 > 35000f0 && (tmd10 = 35000f0)
            tmd0  = const_ * fpow(dia0, -1.605f0); tmd0 > 35000f0 && (tmd0 = 35000f0)
            t85d10 = tmd10 * pmsdiu; t55d10 = tmd10 * pmsdil
            t85d0  = tmd0  * pmsdiu; t55d0  = tmd0  * pmsdil
            local tn10::Float32
            if tt > t85d0
                tn10 = t85d10
            elseif tt > t55d0
                tn10 = abs(t85d0 - tt) <= 5f0 ? t85d10 :
                       _em_tn10_iter(tt, dia0, d10, const_, pmsdil, pmsdiu, t85d10, t55d0, false)
            elseif tt <= t55d10
                tn10 = tt
            else
                tn10 = _em_tn10_iter(tt, dia0, d10, const_, pmsdil, pmsdiu, t85d10, t55d0, true)
            end
            tn10 > tt && (tn10 = tt); tn10 < 0.1f0 && (tn10 = 0f0)
            rn = 1f0 - fpow(1f0 - (tt - tn10) / tt, 1f0 / fint)
            tem = t55d10
            density_on = !(tt <= tem || rn <= 0f0)               # RIP==RN gate (morts.f:497-501)
            if tn10 >= 0.1f0
                # NWCMRT distributes the density excess (or, when the self-thin is off, the background total).
                tokill = density_on ? max(tt - tn10, 0f0) : bg_tokill
                _varmrt!(killed, efftr, temwk2, s, s.variant, t, n, tokill)
            else
                # Whole-stand kill (TN10<0.1): morts.f skips NWCMRT and keeps the uniform per-tree WK2 (RN≈1).
                fill!(killed, 0f0)
                @inbounds for i in 1:n
                    pr = t.tpa[i]; pr <= 0f0 && continue
                    killed[i] = pr * (1f0 - fpow(1f0 - min(rn, 1f0), fint))
                    killed[i] > pr && (killed[i] = pr)
                end
            end
            density_on || break                                 # background ⇒ no d10 dependence, one pass
            # Post-mortality Reineke diameter d10n from the survivors (morts.f:551-582, LZEIDE path).
            ttn = 0f0; sdr = 0f0
            for i in 1:n
                d = t.dbh[i]; pr = t.tpa[i] - killed[i]; pr <= 0f0 && continue
                bark = bark_ratio(bark_a, bark_b, Int(t.species[i]), d)
                g = _mort_traj_g(t.diam_growth[i], d, bark, fint, yr)
                sdr += pr * fpow(d + g, 1.605f0); ttn += pr
            end
            ttn <= 0f0 && break
            d10n = fpow(sdr / ttn, 1f0 / 1.605f0)
            (abs(d10 - d10n) <= 0.1f0 || d10n <= dia0) && break
            d10 = d10n
        end
    end
    # BAMAX residual-BA cap (nc/morts.f:685-754 + sdical.f:203-208). When the user did
    # NOT set BAMAX, live caps residual BA at BAMAX = SDIMAX·0.5454154·PMSDIU (SDI max at
    # 10" DBH) and scales all mortality up proportionally (ADJFAC, iterated ≤100×) until
    # residual BA ≤ BAMAX. This is the QMD≥10"/high-BA regime beyond the Zeide self-thin;
    # with NWCMRT's small-tree-concentrated kill it removes little BA per tree ⇒ drives the
    # residual TPA far below the uniform result (the mature over-max-SDI TPA fix). Inert
    # (immediate break) when residual BA ≤ BAMAX. sdimax<5 (full-kill) handled above.
    if sdimax >= 5f0
        bamax = s.control.ba_max > 0f0 ? s.control.ba_max : sdimax * 0.5454154f0 * pmsdiu
        for _ in 1:100
            banew = 0f0; badead = 0f0
            @inbounds for i in 1:n
                d = t.dbh[i]; sp = Int(t.species[i])
                bark = bark_ratio(bark_a, bark_b, sp, d)
                g = _mort_traj_g(t.diam_growth[i], d, bark, fint, yr)
                ba = 0.0054542f0 * fpow(d + g, 2f0)
                banew  += ba * (t.tpa[i] - killed[i])
                badead += ba * killed[i]
            end
            ((banew - bamax) > 1f0 && badead > 0f0) || break
            adjfac = (banew - bamax) / badead
            @inbounds for i in 1:n
                wki = killed[i] * (1f0 + adjfac)
                wki > t.tpa[i] && (wki = t.tpa[i])
                killed[i] = wki
            end
        end
    end
    apply_fixmort!(s, killed, n, fint)
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
