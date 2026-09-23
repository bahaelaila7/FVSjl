# =============================================================================
# crown_init.jl — cycle-0 (LSTART) crown initialisation shared across variants.
#
# Every variant's crown.f ends with the same `DO 79 I=IREC2,MAXTRE` block: each inventory standing-dead record
# with a missing crown (ICR<=0) is dubbed with the variant's own DUBSCR (a rejection-bounded main-stream BACHLO
# draw), re-expressed on the normal height for top-killed records, and bounded to [10,95]. The draws are consumed
# BEFORE DGDRIV's, so skipping them shifts every later DGSCOR/REGENT draw of the stand (measured on the PN and SO
# WRD fixtures). `dub_dead_crowns!` is that block; the variant supplies DUBSCR as `dub(i) -> CR fraction`.
# =============================================================================

"""
    dub_dead_crowns!(dub, s)

crown.f `DO 79 I=IREC2,MAXTRE` (all variants): dub missing crowns on the cycle-0 dead records. FVS files dead
records from MAXTRE downward in read order, so IREC2→MAXTRE is the REVERSE of jl's dead storage (t.n+1:t.n+ndead).
`dub(i)` returns the variant's ICRI BEFORE the shared top-kill/bounds tail — most variants `INT(CR*100+.5)` after
DUBSCR (use `icri_round`), but the crown-length forms (CR GEMCR, CI/EM/IE CL/H) truncate `INT(CR*100.)`.
"""
@inline icri_round(cr::Real) = trunc(Int, Float32(cr) * 100f0 + 0.5f0)
function dub_dead_crowns!(dub, s::StandState)
    t = s.trees; n = t.n
    t.ndead > 0 || return s
    @inbounds for i in (n + Int(t.ndead)):-1:(n + 1)
        Int(t.crown_pct[i]) > 0 && continue
        icri = Int(dub(i))
        if t.trunc[i] != 0
            hn = Float32(t.norm_ht[i]) / 100f0; hd = hn - Float32(t.trunc[i]) / 100f0
            cl = (Float32(icri) / 100f0) * hn - hd
            icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
        end
        icri > 95 && (icri = 95); icri < 10 && (icri = 10)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# bm/cratet.f LSTART crown-init: DENSE runs over the FULL inventory (live + HISTORY 6-9 standing-dead records),
# and that dead-inclusive density (BA/AVH/point-CCF) is what CROWN→DUBSCR sees when dubbing the D<1 / missing-CR
# LIVE trees. jl partitions the dead into t.n+1:t.n+ndead, so it computes the dead-inclusive scalars by temporarily
# extending the live range, then restores live-only (the grow cycle recomputes density before use). Measured on
# 504443988126144: dead-inclusive BA 55.56 / AVH 85.07 / TPCCF 84.3 (vs live-only 35/45/63) — matches live DUBSCR.
# SHARED (one implementation per Fortran mechanism): bm/, nc/ and pn/cratet.f run the identical LSTART sequence
# (RDPSRT(.FALSE.) IND → `LBKDEN=IDG.LT.2; CALL DENSE` over live+dead → CROWN), so BM, NC and PN all call this.
function crown_init_lstart_dead_inclusive!(s::StandState)
    t = s.trees
    nlive = t.n
    # bm/cratet.f:189-195 `LBKDEN = IDG.LT.2; CALL DENSE` — CROWN (cratet.f:610) dubs against THAT density, whose
    # live WK3 is BACKDATED to the start of the measured-growth period (dense.f:70-128: measured-DG trees by their
    # own increment, the rest by the stand-average BAGR) whenever ≥1 live record carries a measured increment
    # (SN>0; else WK3=DBH). Not backdating put the DUBSCR BA/TPCCF at the current-DBH values (1127576412290487:
    # 3 LP with past DBH ⇒ live BA 16.486 / TPCCF 80.04 vs jl 48.25 / 199.2 ⇒ every seedling crown mis-dubbed).
    # _backdate_dbh! is the shared IDG-faithful dense.f port (live 1:t.n, in place); restored below.
    lbkden = s.control.growth_idg < 2
    # base/fvs.f:193-196 — `SDICLS(0,0.,999.,1,SDIAC,SDIAC2,…)` runs IMMEDIATELY BEFORE `CALL CRATET`,
    # with the source comment "SDICLS IS CALLED HERE SO CROWNS WILL DUB CORRECTLY IN VARIANTS USING THE
    # WEIBULL DISTRIBUTION". So the LSTART dub's RELSDI = SDIAC/SDIDEF is NOT zero: it is the stand
    # Reineke SDI over the LIVE inventory at its READ diameters (before CRATET's backdating DENSE).
    # jl passed no crown_sdi ⇒ RELSDI=0 ⇒ ACRNEW = C0 (its maximum) ⇒ every Weibull-dubbed inventory
    # crown came out too high. MEASURED on the OC control stand (S248112, a missing-crown LP D=11.5):
    # live SDIAC 196.15 ⇒ ICR 76; jl with SDIAC=0 ⇒ 85.
    sdiac = stand_sdi_reineke(s)
    saved_live = lbkden ? t.dbh[1:nlive] : Float32[]
    # #151: dense.f:83-87 — in the CRATET backdating DENSE, standing-dead records get WK3=DBH EXCEPT IMC(I)==9
    # (HISTORY 8,9, older-dead) which LOAD DBH=0 ⇒ they add 0 to BA/CCF/SDI while their HEIGHT still counts
    # toward AVH (measured on 449747082489998: live DUBSCR AVH 67.34 vs jl 1.01 without the dead heights).
    # notre.f:119-124 expands EVERY inventory-dead record (IREC2..MAXTRE) with VP/FP/FP2×(FINT/FINTM), so the dead
    # PROB that CRATET's DENSE sums is TREE_COUNT×FINT/FINTM (BM 10/5 ⇒ ×2; FVS divides it back only where a true
    # density is needed — snag init cratet.f:565, dbstrls MortPA). jl's shared notre! leaves dead TPA unscaled, so
    # scale them here for this pass (24001521010900: 5 HISTORY-8 dead at PROB 12.036/1.998 ⇒ TPCCF 77.788→77.813).
    saved = Tuple{Int,Float32}[]
    saved_tpa = Float32[]
    if t.ndead > 0
        t.n = nlive + Int(t.ndead)
        fintr = s.control.growth_fintm > 0f0 ? s.control.growth_fint / s.control.growth_fintm : 1f0
        @inbounds for i in (nlive + 1):(nlive + Int(t.ndead)); push!(saved_tpa, t.tpa[i]); t.tpa[i] *= fintr; end
    end
    # AVHT40 top height from REAL DBH/HT over live + all dead (IND = real-DBH sort, not WK3), with the dead at their
    # notre-expanded PROB (24001521010900: the ×2 dead fill the top 40 TPA ⇒ AVH 48.466 = live, 41.596 unscaled).
    # dense.f:285-297 AVH walk of the :195 CRATET DENSE, over the :164-166 IND (bm_cratet166_ind: LNKCHN-seeded,
    # dead included, RDPSRT(.FALSE.) on real dbh) — the AVH the :610 CROWN→DUBSCR dub reads. The generic double
    # sort put different equal-DBH records at the 40-TPA boundary (23899355010900 AVH 72.567 vs live 71.606).
    # That DENSE runs BEFORE the missing-height dub (DO 130 :363, DO 145 :464) ⇒ heights are as READ (missing = 0):
    # use the pre-dub snapshot cratet_ht_in (41135212010497: dubbed 1.01 seedlings ⇒ AVH 33.671 vs live 33.160).
    avht_real = let ntot = t.n, ord = bm_cratet166_ind(s, view(t.dbh, 1:t.n), nlive, t.n), hin = s.calib.cratet_ht_in
        use_hin = length(hin) == ntot
        avh = 0f0; ssumn = 0f0
        @inbounds for k in 1:ntot
            ii = Int(ord[k]); p = t.tpa[ii]
            ssumn + p > 40f0 && (p = 40f0 - ssumn)
            ssumn += p; avh += (use_hin ? hin[ii] : t.height[ii]) * p
            ssumn >= 40f0 && break
        end
        ssumn > 0f0 ? avh / ssumn : 0f0
    end
    if lbkden                                 # backdate LIVE WK3 only (after the real-DBH AVH ranking)
        t.n = nlive; _backdate_dbh!(s); t.n = nlive + Int(t.ndead)
    end
    if t.ndead > 0
        @inbounds for i in (nlive + 1):(nlive + Int(t.ndead))
            (t.history[i] == 8 || t.history[i] == 9) || continue
            push!(saved, (i, t.dbh[i])); t.dbh[i] = 0f0
        end
    end
    compute_density!(s)                    # CRATET DENSE: backdated live (+ dead-inclusive) BA / point-CCF
    s.calib.cratet_relden = stand_ccf(s)   # RELDEN after cratet.f:195 DENSE (backdated, dead-inclusive) → REGENT HCOR cal
    @inbounds for (i, d) in saved; t.dbh[i] = d; end
    @inbounds for (k, i) in enumerate((nlive + 1):(nlive + length(saved_tpa))); t.tpa[i] = saved_tpa[k]; end
    t.n = nlive
    lbkden && @inbounds(for i in 1:nlive; t.dbh[i] = saved_live[i]; end)
    s.plot.avg_height = avht_real
    crown_ratio_update!(s, s.variant; lstart = true, crown_sdi = sdiac)   # DUBSCR-dub live D<1 seedlings + Weibull-dub missing-CR overstory
    compute_density!(s; cratet_ind = true)  # restore live-only density; CRATET IND ⇒ cycle-0 PCT/AVH (cratet.f:692 DENSE)
    return s
end


"""
    point_crown_inputs(s) -> (prd, qmdplt, tpccf)

The per-inventory-point inputs crown.f computes for DUBSCR (2021 edits, identical in every variant): PRD = ZRD/XMAXPT
(SDICAL + SDICLS Zeide, live trees, 0 when XMAXPT≤0), QMDPLT = sqrt((PTBAA/PTPA)/0.005454) floored at 1 (1 when the
point has no trees), TPCCF = PCCF(ITRE). Returned as closures over the point index.
"""
function point_crown_inputs(s::StandState)
    xmaxpt, zrd, _ = point_zeide!(s)
    dens = s.density
    prd(pt) = (1 <= pt <= length(xmaxpt) && xmaxpt[pt] > 0f0) ? zrd[pt] / xmaxpt[pt] : 0f0
    function qmdplt(pt)
        baplt = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : 0f0
        tpaplt = (1 <= pt <= length(dens.point_tpa)) ? dens.point_tpa[pt] : 0f0
        q = tpaplt > 0f0 ? sqrt((baplt / tpaplt) / 0.005454f0) : 1f0
        q <= 1f0 ? 1f0 : q
    end
    tpccf(pt) = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
    return prd, qmdplt, tpccf
end
