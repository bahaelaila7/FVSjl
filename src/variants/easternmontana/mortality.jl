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

# em/morts.f MORCON: IPDG/IPDG2(ITYPE,IFOR) → POT index → GMULT/REIN (added-species Hamilton). IFOR 1-6
# collapses to two mappings: forests 1-4 (Beaverhead/Custer/Deerlodge/Gallatin)=Bitterroot, 5-6
# (Helena/Lewis&Clark)=Lolo. POT[k] = the DATA POT literals .25..2.90 (_EM_POT below).
const EM_IPDG_BITT  = Int[7,6,6,6,6,6,5,6,5,6,6,6,6,8,7,7,7,7,7,5,3,4,3,4,4,6,5,1,1,6]
const EM_IPDG_LOLO  = Int[7,7,7,7,6,7,5,6,7,6,6,7,7,9,9,8,9,8,8,5,3,5,4,5,5,7,5,2,4,7]
const EM_IPDG2_BITT = Int[30,29,29,29,28,28,27,31,27,27,28,31,32,32,31,31,32,31,31,25,23,24,23,24,24,27,26,18,16,27]
const EM_IPDG2_LOLO = Int[31,29,30,31,29,30,29,33,31,30,29,34,34,35,34,34,35,33,33,27,23,26,26,27,26,31,26,19,23,31]

# em/morts.f label-220 iterative linear-fn fit between the 55%/85% SDI lines → TN10 (target tree count at
# D10). IPATH2=true (came from T≤T55D0, T>T55D10) computes the line once (TEM=T) then goes to 230; IPATH2=false
# (55%<T≤85% at DIA0) Newton-iterates TREEIT (≤100) so exp(CEPT+SLP·ln(DIA0)) ≈ T. TN10 = exp(CEPT+SLP·ln(D10)),
# capped at T85D10. ALOG/EXP/** are the gfortran libm (flog/fexp/fpow).
function _em_tn10_iter(tt::Float32, dia0::Float32, d10::Float32, const_::Float32, pmsdil::Float32,
                       pmsdiu::Float32, t85d10::Float32, t55d0::Float32, ipath2::Bool;
                       dens::Union{Nothing,Density} = nothing)::Float32
    treeit = tt + 0.1f0 * tt
    slp = 0f0; cept = 0f0; knt = 1
    while true
        tem = ipath2 ? tt : treeit
        d55m = (flog(tem) - flog(pmsdil * const_)) / (-1.605f0)
        t55m = flog(tem)
        d85m = d55m * 1.25f0
        while true                                    # em/morts.f label 221: bump D85M until SLP ≤ −0.5
            d85m > 5f0 && (d85m = 5f0); d85m < 0.125f0 && (d85m = 0.125f0)
            t85m = flog(const_ * fpow(fexp(d85m), -1.605f0) * pmsdiu)
            slp = (t85m - t55m) / (d85m - d55m)
            (slp > -0.5f0 && d85m < 5f0) ? (d85m += 0.1f0) : break
        end
        cept = t55m - slp * d55m
        (ipath2 || tt <= t55d0) && break              # GOTO 230 (no Newton for the IPATH=2 path)
        tprime = cept + slp * flog(dia0)
        diff = tt - fexp(tprime)
        (diff <= 5f0 && diff >= -5f0) && break
        treeit += 0.5f0 * diff; knt += 1
        knt > 100 && break
    end
    # em/morts.f:230 LATCH: the first solved line is kept (SLPMRT/CEPMRT, VARCOM) and reused in later cycles.
    if dens !== nothing
        dens.mort_slope == 0f0 && (dens.mort_slope = slp)
        dens.mort_intercept == 0f0 && (dens.mort_intercept = cept)
        slp = dens.mort_slope; cept = dens.mort_intercept
    end
    tn10 = fexp(cept + slp * flog(d10))
    tn10 >= t85d10 && (tn10 = t85d10)
    return tn10
end

# em/morts.f DATA POT (.25,.30,…,2.90): the REAL*4 decimal literals (not 0.05k+0.20 evaluated in Float32).
const _EM_POT = Float32[Float32(round(0.20 + 0.05 * k; digits = 2)) for k in 1:54]
@inline _em_pot_lit(k::Int) = _EM_POT[clamp(k, 1, 54)]

function mortality!(s::StandState, ::EasternMontana; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t, ctl = s.plot, s.trees, s.control
    n = t.n; n == 0 && return _clim_mort_empty!(s, fint)
    ba = p.basal_area
    icyc = Int(ctl.cycle) + 1                                  # FVS ICYC
    cur_year = current_cycle_year(s)
    zeide = ctl.zeide_sdi
    dthresh = zeide ? ctl.dbh_zeide : ctl.dbh_stage           # em/morts.f:316-317 LZEIDE?DBHZEIDE:DBHSTAGE skip
    fr10 = fint / 10.0f0                                       # G = (DG/BARK)*(FINT/10.0)
    # em/morts.f:280 RMSQD==0 ⇒ reset the latched self-thinning line.
    stand_qmd(s) == 0f0 && (s.density.mort_intercept = 0f0; s.density.mort_slope = 0f0)
    # em/morts.f:305-332 — stand sums in SPECIES-MAJOR IND1 order (DO 20 ISPC / DO 12 I3 / I=IND1(I3)).
    isct = ctl.sp_count_tab; ind1 = s.scratch.idx1
    tt = 0f0; sdq0 = 0f0; sd2sq = 0f0; sumdr0 = 0f0; sumdr10 = 0f0
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        for k in i1:isct[sp, 2]
            i = Int(ind1[k]); (1 <= i <= n) || continue
            d = t.dbh[i]; d < dthresh && continue
            pr = t.tpa[i]
            g = (t.diam_growth[i] / em_bratio(Int(t.species[i]), d)) * fr10
            ciobds = 2f0 * d * g + g * g
            sd2sq += pr * (d * d + ciobds)
            sdq0  += pr * fpow(d, 2f0)                          # P*(D)**2. (REAL exponent ⇒ powf)
            if zeide
                sumdr10 += pr * fpow(d + g, 1.605f0)
                sumdr0  += pr * fpow(d, 1.605f0)
            end
            tt += pr
        end
    end
    # em/morts.f:343-349 trajectory change (ICYC>1 and |T−TPAMRT|>1: thin, ingrowth, fire…) ⇒ reset the line.
    (icyc > 1 && abs(tt - s.density.tpa_mort) > 1f0) && (s.density.mort_intercept = 0f0; s.density.mort_slope = 0f0)
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    tt < 1f0 && @goto morts45                                  # em/morts.f:351 IF(T.LT.1.0) GO TO 45
    dq0  = sqrt(sdq0 / tt)
    dq10 = sqrt(sd2sq / tt)
    dia0 = zeide ? fpow(sumdr0 / tt, 1f0 / 1.605f0) : dq0
    d10  = zeide ? fpow(sumdr10 / tt, 1f0 / 1.605f0) : dq10
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    # SDICAL(0,SDIMAX) (em/morts.f:454): XMAX after CLMAXDEN sets CONST and the SDIMAX<5 kill-all; BAMAX (sdical.f:203-204,
    # BAMAX=XMAX·0.5454154·PMSDIU before CLMAXDEN, unless the user set BAMAX) is what RZ/RIPP and the BA check read.
    sdimax0 = stand_sdimax(s)
    sdimax = clim_sdical_xmax(s, sdimax0, fint)
    const_ = sdimax / 0.02483133f0
    bamax = ctl.ba_max > 0f0 ? ctl.ba_max : sdimax0 * 0.5454154f0 * pmsdiu
    bamax <= 0f0 && (bamax = 1f0)
    # AVED (em/morts.f:394-401): DSUM/WPROB over ALL records in record order (no DBHSTAGE filter).
    dsum = 0f0; wprob = 0f0
    @inbounds for i in 1:n; wprob += t.tpa[i]; dsum += t.dbh[i] * t.tpa[i]; end
    aved = dsum / wprob
    # MORCON (added-species Hamilton): GMULT/REIN from POT(IPDG/IPDG2(ITYPE,IFOR)).
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 6) && (ifor = 1)
    ipdg  = ifor <= 4 ? EM_IPDG_BITT  : EM_IPDG_LOLO
    ipdg2 = ifor <= 4 ? EM_IPDG2_BITT : EM_IPDG2_LOLO
    poten1 = _em_pot_lit(ipdg[itype]);  gmult1 = 0.90f0 / poten1; rein1 = (1f0 - fpow(poten1 / 20f0 + 1f0, -1.605f0)) / 0.06821f0
    poten2 = _em_pot_lit(ipdg2[itype]); gmult2 = 2.50f0 / poten2; rein2 = (1f0 - fpow(poten2 + 1f0, -1.605f0)) / 0.86610f0
    # OLDFNT (grincr.f:60-64): the previous cycle's length; ICYC=1 ⇒ the FINT in force before the projection.
    oldfnt = icyc >= 2 ? Float32(cycle_period_at(ctl, icyc - 2)) : ctl.growth_fint
    sqba = sqrt(ba)
    tn = 0f0; ipass = 0
    # ---- em/morts.f label 10: the QMD-convergence pass (≤10) — each pass recomputes BA10/RZ and TN10 from the
    # CURRENT D10, kills, then re-estimates D10N from the survivors; D10=D10N and repeat while |D10−D10N|>0.1. ----
    while true
        deltba = 0.005454154f0 * d10 * d10 * tt - ba
        if dia0 < 0.3f0                                        # em/morts.f:366-370
            d10 = 0.3f0 + d10 - dia0; dia0 = 0.3f0
        end
        ba10 = ba + ((bamax - ba) / bamax) * deltba
        tb = ba10 / (0.005454154f0 * d10 * d10)
        ttb = (tt - tb) / tt; ttb > 0.9999f0 && (ttb = 0.9999f0)
        rz = 1f0 - fpow(1f0 - ttb, 0.1f0)
        local tn10::Float32
        if sdimax < 5f0                                        # em/morts.f:459-462 ⇒ TN10=0, GO TO 271
            tn10 = 0f0
        else
            tt > 35000f0 && (tt = 35000f0)                     # em/morts.f:465
            tmd0 = const_ * fpow(dia0, -1.605f0); tmd0 > 35000f0 && (tmd0 = 35000f0)
            t85d0 = tmd0 * pmsdiu; t55d0 = pmsdil * tmd0
            tmd10 = const_ * fpow(d10, -1.605f0); tmd10 > 35000f0 && (tmd10 = 35000f0)
            t85d10 = tmd10 * pmsdiu; t55d10 = pmsdil * tmd10
            if tt > t85d0
                tn10 = t85d10                                  # kill to the 85% line
            elseif tt > t55d0
                tn10 = abs(t85d0 - tt) <= 5f0 ? t85d10 :
                       _em_tn10_iter(tt, dia0, d10, const_, pmsdil, pmsdiu, t85d10, t55d0, false; dens = s.density)
            elseif tt <= t55d10
                tn10 = tt                                      # below 55% at both — hold (RN=0)
            else
                tn10 = _em_tn10_iter(tt, dia0, d10, const_, pmsdil, pmsdiu, t85d10, t55d0, true; dens = s.density)
            end
        end
        tn10 > tt && (tn10 = tt); tn10 < 0.1f0 && (tn10 = 0f0)   # label 271
        rn = 1f0 - fpow(1f0 - ((tt - tn10) / tt), 1f0 / fint)
        tem = const_ * fpow(d10, -1.605f0); tem > 35000f0 && (tem = 35000f0); tem = tem * pmsdil   # em/morts.f:639-641
        # per-tree rates, species-major IND1 order (em/morts.f DO 50 ISPC / DO 40 I3)
        @inbounds for sp in 1:MAXSP
            i1 = isct[sp, 1]; i1 == 0 && continue
            for k in i1:isct[sp, 2]
                i = Int(ind1[k]); (1 <= i <= n) || continue
                pr = t.tpa[i]; killed[i] = 0f0
                pr <= 0f0 && continue
                d = t.dbh[i]
                bark = em_bratio(sp, d)
                if _em_orig_species(sp)
                    ri = 1f0 / (1f0 + fexp(EM_PMSC[sp] + EM_PMD[sp] * d + EM_PMDSQ[sp] * d * d))
                    ri = 0.5f0 * ri                            # background rate halved (em/morts.f:622)
                    rip = rn
                    (tt <= tem || rn <= 0f0) && (rip = ri)
                    rip > 1f0 && (rip = 1f0)
                    x = active_mort_mult(ctl, sp, cur_year, d) # XMORT in [XMDIA1,XMDIA2), background rate only
                    rip == rn && (x = 1f0)
                    wki = pr * (1f0 - fpow(1f0 - rip, fint)) * x
                    wki > pr && (wki = pr)
                    killed[i] = wki
                else
                    # ADDED species (4-6,11-17,19) — NI Hamilton potential mortality (em/morts.f:659-725).
                    reldbh = d / aved
                    dd = d <= 0.5f0 ? 0.5f0 : d
                    wk1 = t.dg_prev[i]                          # WK1 = previous cycle's applied DG (dgdriv.f:171)
                    dgt = wk1 / oldfnt
                    dd <= 1f0 && dgt < 0.05f0 && (dgt = 0.05f0)
                    (dd <= 5f0 && dd > 1f0 && dgt < 0.05f0) && (dgt = 0.05f0 * (5f0 - dd) / 4f0)
                    g = wk1 / (bark * oldfnt)
                    wk1 / oldfnt < dgt && (g = dgt / bark)
                    ((icyc == 1 || wk1 == 0f0) && t.diam_growth[i] > 0.5f0) && (g = t.diam_growth[i] / (bark * 10f0))
                    ip = dd <= 5f0 ? 2 : 1
                    g = g * (ip == 1 ? gmult1 : gmult2)
                    rip = 2.76253f0 + 0.222310f0 * sqrt(dd) - 0.0460508f0 * sqba + 11.2007f0 * g -
                          0.554421f0 / dd + EM_PMSC[sp] + 0.246301f0 * reldbh + 6.07129f0 * g / dd
                    rip > 70f0 && (rip = 70f0); rip < -70f0 && (rip = -70f0)
                    rip = 1f0 / (1f0 + fexp(rip))
                    rip = rip * (ip == 1 ? rein1 : rein2)      # ·POTENT
                    ripp = ba * rz
                    ba <= bamax && (ripp = ripp + (bamax - ba) * rip)
                    ripp = ripp / bamax
                    ripp < rip && (ripp = rip); ripp > 1f0 && (ripp = 1f0)
                    x = active_mort_mult(ctl, sp, cur_year, dd)
                    # em/morts.f:705-711 establishment "best"-tree immunity: IY(ICYC)≥IESTAT clears it, else
                    # X·(1−clamp((IESTAT−IY)/FINT,0,1)).
                    if t.iestat[i] > 0
                        iyc = Int32(cur_year)
                        iyc >= t.iestat[i] && (t.iestat[i] = Int32(0))
                        x = x * (1f0 - clamp(Float32(t.iestat[i] - iyc) / fint, 0f0, 1f0))
                    end
                    # species-group share of the NI rate: LL full, RM 20%, the rest 60% (em/morts.f:717-723)
                    wki = pr * (1f0 - fpow(1f0 - ripp, fint)) * x
                    sp == 6 ? (wki = wki * 0.2f0) : (sp != 5 && (wki = wki * 0.6f0))
                    wki > pr && (wki = pr)
                    killed[i] = wki
                end
            end
        end
        # em/morts.f:750-798 — survivor QMD (record order) and the convergence test.
        ipass += 1
        tn = 0f0; sd2sqn = 0f0; sumdr10n = 0f0
        @inbounds for i in 1:n
            d = t.dbh[i]; d < dthresh && continue
            pr = t.tpa[i] - killed[i]
            g = (t.diam_growth[i] / em_bratio(Int(t.species[i]), d)) * fr10
            ciobds = 2f0 * d * g + g * g
            sd2sqn += pr * (d * d + ciobds)
            zeide && (sumdr10n += pr * fpow(d + g, 1.605f0))
            tn += pr
        end
        tn == 0f0 && break
        d10n = zeide ? fpow(sumdr10n / tn, 1f0 / 1.605f0) : sqrt(sd2sqn / tn)
        ipass == 10 && break
        (abs(d10 - d10n) > 0.1f0 && d10n > dia0) || break
        d10 = d10n
    end
    # MORTMSB alternate "mature-stand breakup" mortality (em/morts.f:804-853 + msbmrt.f): inert unless MORTMSB set
    # a non-zero self-thinning slope (default SLPMSB=0 ⇒ QMDMSB=999 never exceeded).
    if ctl.msb_slope != 0f0 && sdimax >= 5f0 && d10 > ctl.msb_qmd && tn > 0f0
        qmd = ctl.msb_qmd; slp = ctl.msb_slope
        cepmsb = flog(const_ * fpow(qmd, -1.605f0)) - slp * flog(qmd)
        tmmsb = fexp(cepmsb + slp * flog(d10)); tmore = tn - tmmsb * pmsdiu; tmore < 0f0 && (tmore = 0f0)
        dlo = ctl.msb_dlo; dhi = ctl.msb_dhi
        tpacls = 0f0                                  # em/morts.f:828-834, TPACLS=TPACLS+PROB-WK2 left to right
        @inbounds for i in 1:n
            dbhend = t.dbh[i] + (t.diam_growth[i] / em_bratio(Int(t.species[i]), t.dbh[i])) * fr10
            (dbhend >= dlo && dbhend < dhi) && (tpacls = (tpacls + t.tpa[i]) - killed[i])
        end
        if tmore > tpacls
            @warn "MORTMSB: additional mortality target exceeds the TPA in the DBH class; alternate mortality cancelled."
        elseif tpacls > 0f0
            mflag = Int(ctl.msb_flag)
            temeff = mflag == 3 ? tmore / tpacls : (tpacls * ctl.msb_eff < tmore ? tmore / tpacls : ctl.msb_eff)
            order = sortperm(view(t.dbh, 1:n); rev = true)
            _msbmrt!(killed, t, order, n, temeff, tmore, dlo, dhi, mflag, s.calib.bark_a, s.calib.bark_b, fint)
        end
    end
    # SIZCAP size-cap mortality floor (em/morts.f:857-869): D+G ≥ SIZCAP(IS,1) and no-mortality flag ≠ 1.
    let sc = ctl.sp_size_cap
        @inbounds for i in 1:n
            sp = Int(t.species[i])
            (sc[sp, 1] >= 999f0 || trunc(Int, sc[sp, 3]) == 1) && continue
            d = t.dbh[i]; pr = t.tpa[i]
            g = (t.diam_growth[i] / em_bratio(sp, d)) * fr10
            if (d + g) >= sc[sp, 1]
                killed[i] = max(killed[i], pr * sc[sp, 2] * fint / 10f0)
                killed[i] > pr && (killed[i] = pr)
            end
        end
    end
    # BA check (em/morts.f:874-943): residual BA ≤ BAMAX, scaling every record's kill by (1+ADJFAC) (≤100 passes).
    begin
        for _ in 1:100
            banew = 0f0; badead = 0f0
            @inbounds for i in 1:n
                d = t.dbh[i]
                g = (t.diam_growth[i] / em_bratio(Int(t.species[i]), d)) * fr10
                ba_ = 0.0054542f0 * fpow(d + g, 2f0)          # (0.0054542*(D+G)**2.)
                banew  += ba_ * (t.tpa[i] - killed[i])
                badead += ba_ * killed[i]
            end
            ((banew - bamax) > 1f0 && badead > 0f0) || break   # (BADEAD=0 would make FVS's ADJFAC Inf/NaN)
            adjfac = (banew - bamax) / badead
            @inbounds for i in 1:n
                wki = killed[i] * (1f0 + adjfac)
                wki > t.tpa[i] && (wki = t.tpa[i])
                killed[i] = wki
            end
        end
    end
    # TPAMRT=TNEW (em/morts.f:944): the post-BA-check survivor TPA over ALL records.
    s.density.tpa_mort = sum(t.tpa[i] - killed[i] for i in 1:n; init = 0f0)
    @label morts45
    # Climate-FVS mortality (em/morts.f:962 CALL CLMORTS — after TPAMRT, before FIXMORT): THISYR = IY(ICYC)+FINT/2.
    (s.climate !== nothing && s.climate.active) &&
        apply_climate_mort!(s, killed, Float32(cur_year) + fint / 2f0, fint)
    # FIXMORT (morts.f): forced-mortality override after the BA check; then the dwarf-mistletoe MAX-combine.
    apply_fixmort!(s, killed, n, fint)
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
