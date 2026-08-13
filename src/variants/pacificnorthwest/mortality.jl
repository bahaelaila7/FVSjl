# =============================================================================
# mortality.jl (pacificnorthwest) — PN mortality (pn/morts.f). Chunk 7.
#
# pn/morts.f is BYTE-IDENTICAL to WC's ORGANON logistic-RIP model ⇒ reuse `wc_mort_rip` + the shared
# WC_BM0..5 / WC_MCLASS / WC_MORT_MVALUES / WC_MORT_ALPHA / WC_MORT_BETA arrays. The ONLY difference is
# XSITE1 = **raw SITEAR(16)** (VARACD≠'WC' ⇒ NO Curtis→King conversion). Same DQ10/AVED + integer-PASS
# density self-thin (SDI<SDIMAX AND BA<550). Routes the kill through the shared apply_fixmort! /
# book_mortality_snags! / tpa-reduction.
# =============================================================================

function mortality!(s::StandState, ::PacificNorthwest; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    ba = p.basal_area; avh = p.avg_height
    xsite1 = p.sp_site_index[16]                            # PN: RAW DF SITEAR (no King conversion)
    xsite2 = p.sp_site_index[19]                            # WH SI
    dbhstage = s.control.dbh_stage
    sdimax = stand_sdimax(s)
    ciobds = @view s.scratch.mort_efftr[1:n]; g1 = @view s.scratch.mort_temwk2[1:n]
    @inbounds for i in 1:n
        d = t.dbh[i]; sp = Int(t.species[i])
        bark = wc_bratio(sd, sp, d); g = t.diam_growth[i] / bark
        g1[i] = g; ciobds[i] = 2.0f0 * d * g + g * g
    end
    dsum = 0f0; wprob = 0f0
    @inbounds for i in 1:n; wprob += t.tpa[i]; dsum += t.dbh[i] * t.tpa[i]; end
    aved = wprob > 0f0 ? dsum / wprob : 0.0001f0
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        pr = t.tpa[i]; pr <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        d <= 0.5f0 && (d = 0.5f0)
        cr = Float32(t.crown_pct[i]) * 0.01f0
        bal = (1.0f0 - Float32(t.crown_ratio[i]) / 100.0f0) * ba
        ptbal = s.density.point_bal[i]
        rip = wc_mort_rip(sp, d, cr, bal, ptbal, t.height[i], avh, ba, xsite1, xsite2)
        rip < 0.001f0 && (rip = 0.001f0)
        wki = pr * (1.0f0 - (1.0f0 - rip)^fint)
        wki > pr && (wki = pr)
        sdimax < 5.0f0 && (wki = pr)
        killed[i] = wki
    end
    if sdimax >= 5.0f0
        pass = 1
        while pass <= 100
            sd2sqa = 0f0; ta = 0f0
            @inbounds for i in 1:n
                pr = t.tpa[i]; wki = killed[i] * pass; wki > pr && (wki = pr)
                d = t.dbh[i]; d < dbhstage && continue
                sd2sqa += (pr - wki) * (d * d + ciobds[i]); ta += (pr - wki)
            end
            ta <= 0f0 && break
            dq10a = sqrt(sd2sqa / ta); baa = 0.005454154f0 * dq10a * dq10a * ta
            sdia = ta * (dq10a / 10.0f0)^1.605f0
            (sdia < sdimax && baa < 550.0f0) && break
            pass += 1
        end
        pass > 100 && (pass = 100)
        if pass > 1
            @inbounds for i in 1:n
                wki = killed[i] * pass; wki > t.tpa[i] && (wki = t.tpa[i]); killed[i] = wki
            end
        end
    end
    apply_fixmort!(s, killed, n, fint)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
