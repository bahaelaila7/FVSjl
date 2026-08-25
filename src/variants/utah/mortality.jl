# =============================================================================
# mortality.jl (utah) — UT periodic mortality (ut/morts.f + ut/utmrt.f). Chunk 7.
#
# Background RI = 0.5·(1/(1+exp(PMSC[sp]+PMD[sp]·D))), overridden by the Zeide-SDI self-thinning rate above
# the mature-stand boundary. UT is a Zeide-SDI variant (LZEIDE=T, ut/grinit.f:133) ⇒ stand_sdimax takes the
# Zeide path. WKI = P·(1−(1−RIP)^FINT); SDIMAX<5 ⇒ WKI=P (climate-shift override, ut/morts.f:331).
#
# ★ Like NC/CR/BM (and UNLIKE EM/CI/IE which have NO distributor), UT DISTRIBUTES the density (or background)
# mortality TOTAL across records by **UTMRT** (ut/utmrt.f, CALL at ut/morts.f:553) — a BA-percentile (PCT,
# from DENSE) + shade-tolerance (VARADJ) geometric progression that CONCENTRATES kills on suppressed (small)
# trees — then runs the morts.f QMD-convergence loop (ut/morts.f:257-603, `GO TO 10`) and the BAMAX residual-BA
# cap. The old port applied a UNIFORM per-tree rate (matches the TOTAL first-pass kill but NOT the distribution);
# combined with the BAMAX cap that UNDER-killed cycle-1 TPA on mature over-max-SDI stands (live drives QMD up by
# killing small trees; uniform kept them). UTMRT is structurally identical to NC's NWCMRT / CR's VARMRT, so it
# reuses the shared `_varmrt!` driver via a UT-specific `_varmrt_efftr!` (EFFTR = PEFF·VARADJ·0.01 — ★ ·0.01,
# cf. NC's ·0.1; species 17-19,22 use the CR crown-ratio form PEFF·((100−CRI)/100)·VARADJ·0.01, ut/utmrt.f).
# =============================================================================

const UT_PMSC = Float32[6.5112, 6.5112, 7.2985, 5.1677, 9.6943, 5.1677, 5.9617, 9.6943, 5.1677, 5.5877,
                        5.1677, 5.1677, 5.1677, 5.1677, 5.1677, 5.1677, 5.1677, 5.9617, 5.9617, 5.9617,
                        5.5877, 5.9617, 5.1677, 5.1677]
const UT_PMD  = Float32[-0.0052485, -0.0052485, -0.0129121, -0.0077681, -0.0127328, -0.0077681, -0.0340128,
                        -0.0127328, -0.0077681, -0.005348, -0.0077681, -0.0077681, -0.0077681, -0.0077681,
                        -0.0077681, -0.0077681, -0.0077681, -0.0052485, -0.0052485, -0.0340128, -0.005348,
                        -0.0052485, -0.0077681, -0.0077681]

# ut/utmrt.f VARADJ — species shade tolerance (1.0 = most intolerant), the EFFTR scalar.
const UT_VARADJ = Float32[0.80, 0.70, 0.55, 0.50, 0.50, 1.00, 0.90, 0.50, 0.60, 0.85,
                          0.70, 0.70, 0.70, 0.70, 0.70, 0.70, 0.90, 0.90, 0.90, 1.10,
                          0.70, 0.90, 0.75, 0.70]

# ut/utmrt.f:97-131 — per-record mortality efficiency EFFTR = PEFF·VARADJ·0.01, PEFF a cubic of the BA
# percentile PCT (t.crown_ratio, filled by stand_pct!). Species 17-19,22 (GB/NC/FC/BE, CR-variant equations)
# use the CR crown-ratio form PEFF·((100−CRI)/100)·VARADJ·0.01. PCT**3.0 is a gfortran powf ⇒ fpow (doctrine #8).
function _varmrt_efftr!(efftr, s, ::Utah, t::TreeList, n::Int)
    pct = t.crown_ratio; tpa = t.tpa; sp = t.species
    pass1 = 0f0
    @inbounds for i in 1:n
        pe = clamp(0.84525f0 - 0.01074f0 * pct[i] + 0.0000002f0 * fpow(pct[i], 3f0), 0.01f0, 1f0)
        spi = Int(sp[i])
        if (17 <= spi <= 19) || spi == 22
            cri = Float32(t.crown_pct[i])
            efftr[i] = pe * ((100f0 - cri) / 100f0) * UT_VARADJ[spi] * 0.01f0
        else
            efftr[i] = pe * UT_VARADJ[spi] * 0.01f0
        end
        pass1 += tpa[i] * efftr[i]
    end
    return pass1
end

function mortality!(s::StandState, ::Utah; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # grown-stand sums (ut/morts.f): T (total tpa), Reineke DR10/DR0 (LZEIDE path, ut/morts.f:218-219,260-263).
    # Using QMD over-stated D10 on dense sub-1" cohorts ⇒ TMD10 uncapped ⇒ TN10 low ⇒ RN self-thin OVER-KILL
    # (real-FIA sp814 oak seedlings: jl killed 45% vs live 11%, #147). g = DG/BARK (UT FINT==YR==10).
    tt = 0f0; sumdr10 = 0f0; sumdr0 = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bark_ratio(bark_a, bark_b, sp, d)
        g = t.diam_growth[i] / bark
        sumdr10 += pr * (d + g)^1.605f0; sumdr0 += pr * d^1.605f0; tt += pr
    end
    tt < 1f-6 && return s
    dr10 = (sumdr10 / tt)^(1f0 / 1.605f0); dia0 = (sumdr0 / tt)^(1f0 / 1.605f0)
    if dia0 < 0.3f0; dr10 = 0.3f0 + dr10 - dia0; dia0 = 0.3f0; end
    sdimax = stand_sdimax(s)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    const_ = sdimax / 0.02483133f0
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    efftr  = @view s.scratch.mort_efftr[1:n]
    temwk2 = @view s.scratch.mort_temwk2[1:n]

    # Background (Hamilton) mortality TOTAL — the target UTMRT distributes when the SDI self-thin is not yet
    # limiting (ut/morts.f DO 50 loop with RIP=RI, then utmrt.f:73-78 TOKILL=ΣWK2). Depends only on start-of-
    # cycle TPA, so compute once. MORTMULT (X) hits only the background rate (ut/morts.f:510-514: X=XMORT in
    # [D1,D2], X=1 when the density rate is in effect).
    cur_year = current_cycle_year(s)
    bg_tokill = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; pr <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        ri = 0.5f0 * (1f0 / (1f0 + exp(UT_PMSC[sp] + UT_PMD[sp] * d)))
        x = active_mort_mult(s.control, sp, cur_year, d)         # MORTMULT (background only)
        bg_tokill += min(pr * (1f0 - (1f0 - min(ri, 1f0))^fint) * x, pr)
    end

    if sdimax < 5f0
        # Climate shifted enough that the site can no longer support trees — kill everything (ut/morts.f:331).
        @inbounds for i in 1:n; killed[i] = t.tpa[i]; end
    else
        # MORTS QMD-convergence iteration (ut/morts.f:257-603, label 10 … GO TO 10): for the assumed end-of-cycle
        # Reineke diameter d10, solve TN10 (self-thin target), distribute the excess T−TN10 (or the background
        # total) by UTMRT, recompute the post-mortality d10n, and re-iterate with d10=d10n until it converges
        # (|d10−d10n|≤0.1) or the QMD would fall below dia0. Up to 10 passes (ut/morts.f:591 IPASS.EQ.10).
        d10 = dr10
        @inbounds for _ in 1:10
            tmd10 = const_ * d10^(-1.605f0); tmd10 > 35000f0 && (tmd10 = 35000f0)
            tmd0  = const_ * dia0^(-1.605f0); tmd0 > 35000f0 && (tmd0 = 35000f0)
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
            rn = 1f0 - (1f0 - (tt - tn10) / tt)^(1f0 / fint)
            tem = t55d10          # SDI-in-effect gate (ut/morts.f:506-508), the missing CAP fix (#147/#156)
            density_on = !(tt <= tem || rn <= 0f0)               # RIP==RN gate (ut/morts.f:504-508)
            if tn10 >= 0.1f0
                # UTMRT distributes the density excess (or, when the self-thin is off, the background total).
                tokill = density_on ? max(tt - tn10, 0f0) : bg_tokill
                _varmrt!(killed, efftr, temwk2, s, s.variant, t, n, tokill)
            else
                # Whole-stand kill (TN10<0.1): morts.f skips UTMRT and keeps the uniform per-tree WK2 (RN≈1).
                fill!(killed, 0f0)
                @inbounds for i in 1:n
                    pr = t.tpa[i]; pr <= 0f0 && continue
                    killed[i] = pr * (1f0 - (1f0 - min(rn, 1f0))^fint)
                    killed[i] > pr && (killed[i] = pr)
                end
            end
            density_on || break                                 # background ⇒ no d10 dependence, one pass
            # Post-mortality Reineke diameter d10n from the survivors (ut/morts.f:565-582, LZEIDE path).
            ttn = 0f0; sdr = 0f0
            for i in 1:n
                d = t.dbh[i]; pr = t.tpa[i] - killed[i]; pr <= 0f0 && continue
                bark = bark_ratio(bark_a, bark_b, Int(t.species[i]), d)
                g = t.diam_growth[i] / bark
                sdr += pr * (d + g)^1.605f0; ttn += pr
            end
            ttn <= 0f0 && break
            d10n = (sdr / ttn)^(1f0 / 1.605f0)
            (abs(d10 - d10n) <= 0.1f0 || d10n <= dia0) && break
            d10 = d10n
        end
    end
    # BAMAX residual-BA cap (ut/morts.f BA-check + vbase/sdical.f:203-208). When the user did NOT set BAMAX,
    # live caps residual BA at BAMAX = SDIMAX·0.5454154·PMSDIU and scales all mortality up proportionally
    # (ADJFAC, iterated ≤100×) until residual BA ≤ BAMAX — the QMD≥10"/high-BA regime beyond the Zeide self-thin.
    # With UTMRT's small-tree-concentrated kill it removes little BA per tree ⇒ drives residual TPA far below the
    # uniform result (the mature over-max-SDI TPA fix). Inert (immediate break) when residual BA ≤ BAMAX.
    if sdimax >= 5f0
        bamax = s.control.ba_max > 0f0 ? s.control.ba_max : sdimax * 0.5454154f0 * pmsdiu
        for _ in 1:100
            banew = 0f0; badead = 0f0
            @inbounds for i in 1:n
                d = t.dbh[i]; sp = Int(t.species[i])
                bark = bark_ratio(bark_a, bark_b, sp, d)
                g = t.diam_growth[i] / bark
                ba = 0.0054542f0 * (d + g)^2
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
    # FIXMORT (morts.f:781): forced-mortality override applied AFTER the BA-check, before the DM combine.
    # Dwarf-mistletoe mortality (mismrt.f): MAX-combine per-tree DM kill into killed[] (inert without DM ratings).
    apply_fixmort!(s, killed, n, fint)
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
