# =============================================================================
# crown_ratio.jl (northeast) — NE per-cycle crown-ratio update (ne/crown.f)
#
# NE uses the TWIGS crown model (GTR NC-125), NOT SN's Weibull-percentile model:
#
#   CRNEW = 10·( BCR1/(1 + BCR2·BA) + BCR3·(1 − exp(BCR4·D)) )            (crown.f:179)
#
# where BA is the stand basal area per acre and D the tree DBH. BCR1..BCR4 are per
# species (data/northeast/species_coefficients.csv crown_bcr1..4). The change-limit
# (±1%/yr), CRNMLT DBH-window multiplier, crown-length cap, top-kill reduction, and
# the 1/10/95 bounds are the SAME structure as SN's CROWN.
# =============================================================================

# CS shares NE's TWIGS crown model (cs/crown.f ≡ ne/crown.f modulo the BCR4 sign, which CS
# data folds in as a negated `crown_bcr4` so the same `exp(bcr4·d)` kernel reproduces CS's
# `exp(-BCR4·d)`). One method serves both eastern variants.
function crown_ratio_update!(s::StandState, ::Union{Northeast,CentralStates}; fint::Float32 = 10f0,
                             crown_sdi::Float32 = -1f0, relden_override::Float32 = -1f0,
                             ba_override::Float32 = -1f0, lstart::Bool = false)
    t = s.trees; sd = s.coef.species; n = t.n
    n == 0 && return s
    # crown.f COMMON BA = the RAW per-acre total basal area (live-stamped: cst01 ICYC1 BA=109.10 == jl raw
    # basal_area). The FAITHFUL value is raw `basal_area` (no /gross_space), AND it must read the POST-growth
    # BA (the gradd.f DENSE-before-CROWN refresh in simulate.jl fixed the stale pre-growth BA). Together
    # (post-growth + raw) growth_fint10/all_SA go BIT-EXACT. HOWEVER raw currently regresses 7 CS all-species
    # monocultures (TPA off 5-7 — a separate per-species crown issue it unmasks), so raw is DEFERRED pending
    # that investigation; the POST-growth density refresh (the dominant fix) is landed, and `/gross_space` is
    # kept here for now (leaves growth_fint10 at 1.87%, down from 3.72%). ba_override bypasses this (CRATET init).
    ba = ba_override >= 0f0 ? ba_override :
         (s.plot.gross_space > 0f0 ? s.plot.basal_area / s.plot.gross_space : s.plot.basal_area)
    bcr1 = sd[:crown_bcr1]; bcr2 = sd[:crown_bcr2]; bcr3 = sd[:crown_bcr3]; bcr4 = sd[:crown_bcr4]
    cur_year = current_cycle_year(s)
    @inbounds for i in 1:n
        sp = t.species[i]
        icr_old = t.crown_pct[i]
        # crown.f:164 — at init, a tree with an inventory crown is left untouched.
        lstart && icr_old > 0 && continue
        # crown.f:171 — a negative crown means a pest extension already set the change; restore sign.
        if !lstart && icr_old < 0
            t.crown_pct[i] = -icr_old
            continue
        end
        d = t.dbh[i]
        den = 1f0 + bcr2[sp] * ba
        crnew = 10f0 * (bcr1[sp] / den + bcr3[sp] * (1f0 - exp(bcr4[sp] * d)))
        # CRNMLT (crown.f:193-201): the per-species crown-change multiplier, applied only inside its
        # DBH window. active_crn_mult returns 1.0 outside the window / when no CRNMLT keyword (net01).
        cm = active_crn_mult(s.control, sp, cur_year, d)
        change_branch = !(lstart || icr_old == 0)
        if change_branch
            chg = crnew - Float32(icr_old)
            pdifpy = chg / Float32(icr_old) / fint
            pdifpy >  0.01f0 && (chg = Float32(icr_old) *  0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr_old) * -0.01f0 * fint)
            crnew = Float32(icr_old) + chg * cm   # cm==1 outside the CRNMLT window
        end
        icri = trunc(Int32, crnew + 0.5f0)
        if !change_branch                          # crown.f:199-202 — init/dub: scale by CRNMLT in window
            icri = trunc(Int32, Float32(icri) * cm)
        end
        # Crown-length cap (crown.f:206-220) — skipped at init / for dubbed (icr_old==0) trees.
        if change_branch
            crln = t.height[i] * Float32(icr_old) / 100f0
            crmax = (crln + t.ht_growth[i]) / (t.height[i] + t.ht_growth[i]) * 100f0
            Float32(icri) > crmax && (icri = trunc(Int32, crmax + 0.5f0))
            (icri < 10 && cm == 1f0) && (icri = trunc(Int32, crmax + 0.5f0))
        end
        # Top-kill reduction at init (crown.f:224-228).
        if lstart && t.trunc[i] != 0
            hn = Float32(t.norm_ht[i]) / 100f0
            if hn > 0f0
                hd = hn - Float32(t.trunc[i]) / 100f0
                cl = (Float32(icri) / 100f0) * hn - hd
                icri = trunc(Int32, cl * 100f0 / hn + 0.5f0)
            end
        end
        icri > 95 && (icri = Int32(95))
        (icri < 10 && cm == 1f0) && (icri = Int32(10))
        icri < 1 && (icri = Int32(1))
        t.crown_pct[i] = icri
    end
    return s
end
