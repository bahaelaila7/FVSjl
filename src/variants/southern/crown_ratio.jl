# =============================================================================
# crown_ratio.jl — Southern per-cycle crown-ratio update (CROWN)
#
# Ported from: sn/crown.f (CROWN). Runs once per cycle AFTER growth (gradd.f:172,
# i.e. after DG/HTG applied, the tree list re-sorted on dbh, and density recomputed)
# to update each live tree's crown ratio `trees.crown_pct` (ICR). The new crown
# ratio is a Weibull draw at the tree's diameter percentile, recentred on the
# species mean-crown-ratio curve `acrnew(relsdi)`, then limited to a ±1%/yr change
# from the prior crown and capped by the crown-length geometry. Coefficients in
# data/southern/crown_ratio_coeffs.csv (MCREQN mean-CR eqn + WEIBUL params).
#
# Crown ratio is bit-exact at cycle 1 only when this runs AFTER cycle 1's growth
# (so cycle-1 DGF still sees the input crown); without it, the frozen crown made
# every cycle ≥2 diverge (crown → log(ICR) in DGF wk2 → DG → mortality → volume).
# =============================================================================

"""
    eastern_cratet_dead_pct!(s)

The cycle-0 DEAD records' PCT as CRATET's backdating DENSE leaves it (cratet.f:150-195 → dense.f:83-87,244): PCTILE over
the live+dead IND (IND1-seeded, RDPSRT(.FALSE.) on the read DBH) of WK5 = D·(D·P), live D backdated (LBKDEN = IDG<2),
dead D = DBH for HISTORY 6-7 and 0 for 8-9 (IMC 9), dead P expanded ×FINT/FINTM (notre.f). Later DENSEs cover only the
live list, so FVS_TreeList reports this PCT for those rows for good (SN 200267456010854: 0.635/14.38/29.29 for three
IMC-9 snags; jl printed 0). Stored in `calib.cratet_dead_pct`; the live PCT is left as computed.
"""
function eastern_cratet_dead_pct!(s::StandState)
    t = s.trees; nd = Int(t.ndead); nlive = t.n
    nd > 0 || return s
    ntot = nlive + nd
    ind = bm_cratet166_ind(s, view(t.dbh, 1:ntot), nlive, ntot)     # real-DBH IND, dead included
    saved_dbh = t.dbh[1:ntot]; saved_tpa = t.tpa[1:ntot]
    try
        s.control.growth_idg < 2 && _backdate_dbh!(s)                 # live WK3 (t.n = nlive here)
        fintr = s.control.growth_fintm > 0f0 ? s.control.growth_fint / s.control.growth_fintm : 1f0
        @inbounds for j in (nlive + 1):ntot
            t.tpa[j] *= fintr
            (t.history[j] == 8 || t.history[j] == 9) && (t.dbh[j] = 0f0)
        end
        pct = zeros(Float32, ntot)
        t.n = ntot
        _pctile!(pct, t, ind, ntot)
        s.calib.cratet_dead_pct = pct[(nlive + 1):ntot]
        (nlive == 0 && nd == 1) && (s.calib.cratet_dead_pct[1] = 0f0)   # pctile.f N=1 leaves the MAXTRE record at 0
    finally
        t.n = nlive
        @inbounds for j in 1:ntot; t.dbh[j] = saved_dbh[j]; t.tpa[j] = saved_tpa[j]; end
    end
    return s
end

"""
    eastern_cratet_ba(s) -> Float32

COMMON BA as CRATET's `LBKDEN=IDG.LT.2; CALL DENSE` (cratet.f:168-170) leaves it for the LSTART CROWN dub (NE/CS/LS
crown.f DEN=1+BCR2·BA, live AND cycle-0 dead records). At that point IND1/ISCT still hold the input DEAD records
(cratet.f:176-205 drops them only AFTER this DENSE), so dense.f's species loop sums them too: BATREE = 0.005454154·D·(D·P)
with D = the backdated WK3 for live records (dense.f:64-128), the read DBH for HISTORY 6-7 dead (IMC 7) and 0 for 8-9
(IMC 9) (dense.f:83-87); P of a dead record expanded ×FINT/FINTM (notre.f). With LBKDEN the second (current-DBH) pass
is swapped back (dense.f:259-263 BA=OLDBA), so BA = the backdated pass. Without LBKDEN (IDG≥2) one current-DBH pass,
dead included at their DBH. MEASURED live FVScs 1813567613290487 (9 BO snags): the live-only backdated BA dubbed every
missing crown 1-2 points high (PctCr 43 vs live 41).
"""
function eastern_cratet_ba(s::StandState)::Float32
    t = s.trees; n = t.n; nd = Int(t.ndead)
    lbk = s.control.growth_idg < 2
    wk3 = Float32[t.dbh[i] for i in 1:n]
    if lbk
        saved = copy(wk3)
        _backdate_dbh!(s)
        @inbounds for i in 1:n; wk3[i] = t.dbh[i]; t.dbh[i] = saved[i]; end
    end
    fintr = s.control.growth_fintm > 0f0 ? s.control.growth_fint / s.control.growth_fintm : 1f0
    bat = 0f0
    @inbounds for spx in 1:MAXSP                       # dense.f DO 50 ISPC / DO 10 I3 (species-major IND1)
        for i in 1:n
            Int(t.species[i]) == spx || continue
            d = wk3[i]; p = t.tpa[i]
            bat += 0.005454154f0 * (d * (d * p))
        end
        for j in (n + 1):(n + nd)
            Int(t.species[j]) == spx || continue
            d = (lbk && (t.history[j] == 8 || t.history[j] == 9)) ? 0f0 : t.dbh[j]
            p = t.tpa[j] * fintr
            bat += 0.005454154f0 * (d * (d * p))
        end
    end
    return bat
end

"""
    init_crown_ratios!(s)

CRATET: estimate the INITIAL crown ratio for inventory trees that have no input crown
(`crown_pct == 0`), using the CROWN model on the inventory stand (FVS does this in
INITRE/CRATET so the first cycle's DGF and the FFE inventory crown are based on a real
crown, not 0). Trees WITH an input crown are left untouched. No-op if every tree already
has a crown. Idempotent; safe to call once before the first grow.
"""
function init_crown_ratios!(s::StandState)
    t = s.trees
    t.n == 0 && return s
    n = t.n
    # all crowns already set (live AND the cycle-0 dead records, which crown.f DO 79 also dubs)
    (any(@views t.crown_pct[1:n] .== 0) || any(@views t.crown_pct[(n + 1):(n + Int(t.ndead))] .<= 0)) || return s
    compute_density!(s)
    # FVS's CRATET dubs missing crowns with CROWN using DENSE's CCF computed on the BACKDATED dbh
    # (DENSE/LBKDEN: past dbh = sqrt(d²·r), r from the measured DG; unmeasured trees use the stand-average
    # ratio). Replicate that backdated CCF (RELDEN) so the init crown matches — the percentile SCALE depends
    # on it, and using the current (denser) dbh gives a ~5% low crown.
    saved_dbh = copy(@view t.dbh[1:n])
    _backdate_dbh!(s)                                   # dense.f:70-128 backdating (shared w/ DG calibration)
    compute_density!(s)
    bd_relden = stand_ccf(s)                            # CCF on the backdated stand
    bd_ba = s.plot.basal_area                           # RAW backdated BA — the NE/CS crown model's COMMON BA
                                                        # (crown.f uses raw BA, not /gross_space); D39 init dub.
    @inbounds for i in 1:n; t.dbh[i] = saved_dbh[i]; end   # restore current dbh
    # NE: the dead-inclusive CRATET DENSE BA (see eastern_cratet_ba; it backdates the restored current dbh itself)
    s.variant isa Northeast && (bd_ba = eastern_cratet_ba(s))
    compute_density!(s)                                 # rank/SDI use current dbh (as FVS CROWN)
    saved = copy(@view t.crown_pct[1:n])
    crown_ratio_update!(s, s.variant; fint = htg_period(s.variant), relden_override = bd_relden,
                        ba_override = bd_ba, lstart = true)
    @inbounds for i in 1:n
        saved[i] != 0 && (t.crown_pct[i] = saved[i])   # restore input crowns; keep only the estimated 0s
    end
    # crown.f DO 79 — cycle-0 dead records with a missing crown (eastern forms are deterministic, no draw):
    #   SN  (sn/crown.f:386-401): DUBSCR(D,CR): CR = 0.70 − 0.40/24·D (D≤24) else 0.30, bounded [.05,.95]; ICRI=INT(CR·100+.5)
    #   CS/LS/NE (crown.f DO 79): CR = 10·(BCR1/(1+BCR2·BA) + BCR3·(1−EXP(BCR4·D))) (CS's sign folded into crown_bcr4);
    #                             ICRI = INT(CR+.5). BA = the same COMMON BA the live init dub reads.
    v = s.variant
    if v isa Southern
        dub_dead_crowns!(s) do i
            d = t.dbh[i]
            cr = d <= 24f0 ? 0.70f0 - 0.40f0 / 24f0 * d : 0.30f0
            cr < 0.05f0 && (cr = 0.05f0); cr > 0.95f0 && (cr = 0.95f0)
            icri_round(cr)
        end
    elseif v isa Northeast || v isa CentralStates || v isa LakeStates
        sd = s.coef.species
        dub_dead_crowns!(s) do i
            sp = t.species[i]
            den = 1f0 + sd[:crown_bcr2][sp] * bd_ba
            cr = 10f0 * (sd[:crown_bcr1][sp] / den + sd[:crown_bcr3][sp] * (1f0 - fexp(sd[:crown_bcr4][sp] * t.dbh[i])))
            trunc(Int, cr + 0.5f0)
        end
    end
    return s
end

"""
    crown_ratio_update!(state; fint=5f0)

CROWN: update `trees.crown_pct` (ICR, %) for every live record from the post-growth
stand. No-op for an empty stand. Call once per cycle after growth + density.
"""
function crown_ratio_update!(s::StandState, ::Southern; fint::Float32 = 5f0, crown_sdi::Float32 = -1f0,
                             relden_override::Float32 = -1f0, ba_override::Float32 = -1f0, lstart::Bool = false)
    t = s.trees; sd = s.coef.species; n = t.n
    n == 0 && return s
    # CRNMULT: cycle year for the persistent crown-ratio-change multiplier lookup.
    cur_year = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)
    sdiac = crown_sdi >= 0f0 ? crown_sdi : stand_sdi_reineke(s)  # pre-growth Reineke SDIBC (grincr.f:241)
    relden = relden_override >= 0f0 ? relden_override : stand_ccf(s)  # RELDEN — crown competition factor
                                         # (override = the DENSE-backdated CCF, used by CRATET init crown)
    sdidef = s.plot.sp_sdi_def
    # ISORT(IND(JJ)) = ITRN−JJ+1 (sn/crown.f:156-159): the rank in FVS's IND — RDPSRT's DESCENDING diameter order with
    # its own tie order — not a stable ascending sort. LSTART: CRATET's IND (sn/cratet.f:155-157 IND=IND1 + RDPSRT(.FALSE.),
    # :261 RDPSRT(.TRUE.) when dead records exist — the shared crown_isort/bm_cratet_ind!); cycling: gradd.f:186 RDPSRT.
    # MEASURED FVSsn_g16 157577477010854 1972: tied LP pairs (7.1"/49.0', 6.4"/45.4', 5.9"/42.6') had their dubbed CRs
    # 34/35, 31/29, 26/27 swapped against live.
    isort = crown_isort(s; lstart = lstart)
    scale = clamp(1f0 - 0.00167f0 * (relden - 100f0), 0.30f0, 1f0)

    eqn = sd[:mcr_eqn]; ma = sd[:mcr_a]; mc = sd[:mcr_c]; mb = sd[:mcr_b]
    mb2 = sd[:mcr_b2]; mb3 = sd[:mcr_b3]
    wa = sd[:wb_a]; wb0 = sd[:wb_b0]; wb1 = sd[:wb_b1]; wc = sd[:wb_c]
    # Per-species Weibull params (acrnew depends only on species via relsdi). Preallocated Scratch buffers
    # (reused, allocation-free); `seen` is the compute gate, reset false each call. Aw/Bw/Cw need no reset —
    # they are only read where seen[sp] is true, and the lazy per-species compute order is unchanged (bit-exact).
    Aw = s.scratch.crown_aw; Bw = s.scratch.crown_bw; Cw = s.scratch.crown_cw
    seen = s.scratch.crown_seen; fill!(seen, false)
    @inbounds for i in 1:n
        sp = t.species[i]
        if !seen[sp]
            relsdi = sdidef[sp] > 0f0 ? sdiac / sdidef[sp] * 10f0 : 6f0
            relsdi = clamp(relsdi, 1f0, 12f0)
            ie = Int(eqn[sp])
            acrnew = ie == 1 ? exp(ma[sp] + mb[sp] * log(relsdi) + mc[sp] * relsdi) :
                     ie == 2 ? exp(ma[sp] + mb[sp] * log(relsdi)) :
                     ie == 3 ? ma[sp] + mc[sp] * relsdi :
                     ie == 4 ? ma[sp] + mb2[sp] * log10(relsdi) :
                     ie == 5 ? relsdi / (ma[sp] * relsdi + mb3[sp]) : 0f0
            bb = wb0[sp] + wb1[sp] * acrnew; bb < 3f0 && (bb = 3f0)
            cc = wc[sp]; cc < 2f0 && (cc = 2f0)
            Aw[sp] = wa[sp]; Bw[sp] = bb; Cw[sp] = cc; seen[sp] = true
        end
        icr_old = t.crown_pct[i]
        if icr_old < 0   # crown change already computed by the topkill/pest model
            t.crown_pct[i] = -icr_old   # (sn/crown.f:271): restore sign, bypass the recompute
            continue
        end
        d = t.dbh[i]
        # Relative crown position: DBH-rank fraction for live stems, but a RANN draw for dbh ≤ 0 (regen
        # with no DBH) — crown.f:287-292 `IF(DBH>0) X=ISORT/ITRN*SCALE ELSE CALL RANN(RNUMB); X=RNUMB*SCALE`.
        # jl previously used a fixed 0.5, which both mis-set the regen crown AND skipped FVS's RANN draw
        # (desyncing the per-tree DGSCOR RNG stream on the regen path). Drawing in tree-loop order matches FVS.
        x = d > 0f0 ? Float32(isort[i]) / Float32(n) * scale : rann!(s.rng) * scale
        x = clamp(x, 0.05f0, 0.95f0)
        crnew = Aw[sp] + Bw[sp] * ((-log(1f0 - x))^(1f0 / Cw[sp]))
        # Limit change to ±1%/yr of the prior crown (crown.f:442-459).
        if icr_old != 0
            chg = crnew - Float32(icr_old)
            pdifpy = chg / Float32(icr_old) / fint
            pdifpy >  0.01f0 && (chg = Float32(icr_old) *  0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr_old) * -0.01f0 * fint)
            # CRNMULT (crown.f:319): scale the crown-ratio change for trees in the keyword's
            # DBH window (1.0 = no CRNMULT keyword, the common case).
            crnew = Float32(icr_old) + chg * active_crn_mult(s.control, sp, cur_year, d)
        end
        icri = trunc(Int32, crnew + 0.5f0)
        # Crown-length cap: the crown can't exceed (old length + HTG) over new height.
        if icr_old != 0
            crln = t.height[i] * Float32(icr_old) / 100f0
            crmax = (crln + t.ht_growth[i]) / (t.height[i] + t.ht_growth[i]) * 100f0
            (Float32(icri) > crmax || icri < 10) && (icri = trunc(Int32, crmax + 0.5f0))
        end
        # crown.f:55 — at LSTART (init), reduce the dubbed crown of an inventory TOP-KILLED tree:
        # the dead-top portion (NORMHT−ITRUNC) is removed from the crown length, re-expressed over NORMHT.
        # Auto-scoped to dubbed trees: init_crown_ratios! restores input crowns after this call.
        if lstart && t.trunc[i] != 0
            hn = Float32(t.norm_ht[i]) / 100f0
            if hn > 0f0
                hd = hn - Float32(t.trunc[i]) / 100f0
                cl = (Float32(icri) / 100f0) * hn - hd
                icri = trunc(Int32, cl * 100f0 / hn + 0.5f0)
            end
        end
        icri > 95 && (icri = Int32(95))
        icri < 10 && (icri = Int32(10))
        icri < 1  && (icri = Int32(1))
        t.crown_pct[i] = icri
    end
    return s
end
