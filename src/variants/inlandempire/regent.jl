# =============================================================================
# regent.jl (inlandempire) — IE small-tree growth (ie/regent.f, chunk 6). Coefficients in
# regent_coefficients.jl (dumped from live). Reuses the KT subcycle scaffold (kper/banext/rdnext/DELMAX);
# the equations are IE-specific (measured — DIFFER from KT).
#
#   ie_regcons!(s)            — REGCON: RHCON(sp) = REGCH + RHSC(sp) + RHHAB(MAPHAB(ITYPE,sp),sp),
#                               REGCH = RHGL(IGL) + (RSAB0+RSAB1·cosASP+RSAB2·sinASP)·SLOPE.
#   small_tree_growth!(...)   — NPER subcycles; NIVAR height HTGRL=CON+RHLH·lnH+RHCCF·RDJ+RHBAL·BAL,
#                               H2=H1+exp(HTGRL)·SCALE·XRHGRO·WK4; DG via AX/BX (=HHT1/HHT2) power H-D model.
# NIVAR (sp1:12,14,23) fully ported+validated; TTVAR(13,17)/CRVAR(19,22)/UTVAR(15,16,18,20,21) faithful
# (stochastic BETA/POTHTG/aspen-FINDAG) — not in iet01, validated later.
# =============================================================================

"""
IE regent small-tree HEIGHT calibration (ie/regent.f:1054-1337, the LHTCAL/mode-40 pass). Computes the RAW
regent HCOR into `c.htg_cor_init[sp]` for every sub-model (NIVAR, TTVAR, CR, UT): for each sub-5" tree with a
measured height increment, accumulate the predicted regent height growth EDH (NIVAR/TTVAR: HK−H grown over the
NPER subcycles; CR/UT: the 10-yr potential·PCTRED·VIGOR or Sheppard increment, halved) and the measured
TERM=HTG·SCALE3; CORNEW = Σ(TERM·P)/Σ(EDH·P); HCOR_raw =
ln(CORNEW), trapped to [0.0821, 12.1825]. NTYR is IFINTH (5 unless the DB HTG_MEASURE set it) — jl had used
FINT (10) ⇒ NPER 2 vs live 1 ⇒ every NIVAR SNX ≈2× live ⇒ CORNEW halved. The stand values are the cratet.f:218
backdating DENSE's (snapshot), and the CR/UT arms (previously unported) are in. Measured: an IE fixture with
NIVAR/TT/CR/UT seedlings carrying HTG — all 14 per-species SUMS equal FVSie_g16's.

ln(CORNEW), trapped to [0.0821, 12.1825]. `calibrate_diameter_growth!`'s shared attenuation (dgdriv.f:188-194,
`htg_cor_small = dg_cor_goal + cormlt_h·(htg_cor_init − dg_cor_goal)`) then produces the applied HCOR.

Without this, IE NIVAR species had htg_cor_init=0 ⇒ the diameter COR (dg_cor_goal) leaked into the regent
height CON, over-growing the small-tree cohort on dense stands where the calibration fires (#171). Inert where
fewer than NCALHT(5) sub-5" trees carry a measured height increment (e.g. iet01) ⇒ htg_cor_init stays 0.
"""
function ie_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    s.control.growth_ifinth == 0 && return s           # regent.f:1064 IF(IFINTH.EQ.0) GOTO 100
    rhcon = ie_regcons!(s)                              # raw RHCON (no HCOR)
    # REGCAL runs inside CRATET (ie/cratet.f:586) on the state its :218 backdating DENSE left — no DENSE or
    # AVHT40 in between — so TEMBA/TEMCCF (=BA/RELDEN), AVH, RELDM1, PCCF and PCT are that DENSE's
    # (crown_init_lstart_dead_inclusive! snapshots them).
    snap = c.cratet_relden > 0f0
    ba = snap ? c.cratet_ba : p.basal_area; relden = snap ? c.cratet_relden : p.relative_density
    avh = snap ? c.cratet_avh : p.avg_height
    reldm1 = snap ? c.cratet_reldm1 : p.relative_density_prev
    pccfv = (snap && !isempty(c.cratet_pccf)) ? c.cratet_pccf : s.density.point_ccf
    pctv = (snap && length(c.cratet_pct) >= t.n) ? c.cratet_pct : t.crown_ratio
    regyr = IE_RG_REGYR
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : Float32(htg_period(s.variant))
    scale3 = regyr / finth                              # regent.f:1065 SCALE3 = REGYR/FINTH
    # regent.f:199 LSTART ⇒ NTYR = IFINTH (5 unless the DB HTG_MEASURE set it; the GROWTH keyword never does).
    ntyr = Int(s.control.growth_ifinth); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1 && length(s.calib.dub_wk3) == t.n
        @inbounds for i in 1:t.n
            # calibration DO 49 I=1,ITRN (regent.f:1081-1100) has NO D>=3 gate (unlike the growth pass DO 6 :246):
            # every live record's increment feeds the subcycle projection. D1=WK3 (backdated), DG = DGDRIV's DO 220
            # (measured, 0 for HT<=4.5, else the dub) — not the raw input DG (-1 when missing).
            d1 = s.calib.dub_wk3[i]
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = ie_bratio(sp, d1)
            d2 = d1 + ie_do220_dg(s, i, saved_dbh[i]) / bark    # REGCAL runs while t.dbh is still backdated
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            c1 = ie_tree_ccf(sp, d1) * pr; c2 = ie_tree_ccf(sp, d2) * pr     # CCFCAL = CCFT·P (regent.f:1089-1090)
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn; banext[j] += Float32(k) * bi * pn
            end
        end
    end
    # regent.f:1104-1111 CR/UT density modifier from X = AVH·(RELDEN/100)
    x = avh * (relden / 100f0); x > 300f0 && (x = 300f0)
    pctred = 1.11436f0 + x*(-0.011493f0 + x*(0.43012f-4 + x*(-0.72221f-7 + x*(0.5607f-10 - x*0.1641f-13))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # regent.f:1116-1132 assigns each species a sub-model: NIVAR (<=12,14,23), TTVAR (13,17), CRVAR (19,22),
    # UTVAR (the rest). The DO 90 loop runs for EVERY species with ISCT≠0 and LHTCAL (default .TRUE.). The TTVAR
    # arm DRAWS ZRAND(I)=BACHLO(0,1,RANN) (regent.f:1181) from the main stream and leaves it on the record.
    ihtg = s.control.growth_ihtg
    @inbounds for sp in 1:23
        nivar = sp <= 12 || sp == 14 || sp == 23
        ttvar = sp == 13 || sp == 17
        crvar = sp == 19 || sp == 22
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        # CR/UT site terms (regent.f:1215-1221): SI clamped to the regent-local SLO/SHI; SJ = SITEAR unclamped.
        sj = p.sp_site_index[sp]
        si = sj; si > IE_RG_SHI[sp] && (si = IE_RG_SHI[sp]); si <= IE_RG_SLO[sp] && (si = IE_RG_SLO[sp] + 0.5f0)
        relsi = IE_RG_SHI[sp] > IE_RG_SLO[sp] ? (si - IE_RG_SLO[sp]) / (IE_RG_SHI[sp] - IE_RG_SLO[sp]) : 0f0
        rsimod = 0.5f0 * (1f0 + relsi)
        snx = 0f0; sny = 0f0; snp = 0f0; nh = 0
        for k in i1:i2
            i = ind1[k]
            hg = t.ht_growth[i]
            hb = t.height[i]; ihtg < 2 && (hb -= hg)      # regent.f:1154 IF(IHTG.LT.2) H=H-HTG(I)
            (saved_dbh[i] >= 5.0f0 || hb < 0.01f0) && continue   # DBH<5, H>=0.01 (regent.f:1155)
            hg < 0.001f0 && continue                      # measured HTG required (regent.f:1156)
            pct = pctv[i]
            cr  = Float32(t.crown_pct[i])                 # CR = FLOAT(ICR(I)) (regent.f:1159)
            ipccf = Int(t.plot_id[i])
            pccf_i = (ipccf >= 1 && ipccf <= length(pccfv)) ? pccfv[ipccf] : 0f0
            hk = hb; edh = 0f0
            for j in 1:nper
                if nivar
                    bal = banext[j] * (100.0f0 - pct) * 0.0001f0
                    htgrl = rhcon[sp] + IE_RG_RHLH[sp]*flog(hk) + IE_RG_RHCCF[sp]*rdnext[j] +
                            IE_RG_RHBAL[sp]*bal
                    hk += fexp(htgrl)                      # regent.f:1167-1168 (NO scale in the calib pass)
                elseif ttvar
                    # TTVAR (regent.f:1170-1206): PPCCF/TPCCF-driven BETA height increment around a persistent
                    # per-tree ZRAND. EDH<=0.1 clamps and resets ZRAND to -999 so the next use redraws.
                    ppccf = (rdnext[j] - reldm1) / reldm1
                    tpccf = pccf_i * ppccf
                    tpccf > 300.0f0 && (tpccf = 300.0f0)
                    tpccf < 25.0f0  && (tpccf = 25.0f0)
                    if t.zrand[i] == -999f0                # regent.f:911 redraw until |ZRAND|<=2
                        while true
                            z = bachlo(s.rng, 0.0f0, 1.0f0)
                            (z >= -2.0f0 && z <= 2.0f0) && (t.zrand[i] = z; break)
                        end
                    end
                    if saved_dbh[i] <= 0f0
                        edh_j = 0f0
                    else
                        beta1 = fexp(1.17527f0 - 0.42124f0*flog(tpccf))
                        beta2 = fexp(-2.56002f0 - 0.58642f0*flog(tpccf))
                        htg1 = beta1 + beta2*cr
                        stddev = htg1*(1.08720f0 - 0.00230f0*cr)
                        edh_j = htg1 + t.zrand[i]*stddev
                        (edh_j <= 0.1f0) && (edh_j = 0.1f0; t.zrand[i] = -999f0)
                    end
                    hk += edh_j
                else
                    # CR/UT (regent.f:1207-1264): EDH is the LAST sub-period's value (not accumulated). CRVAR and
                    # PI/UJ (15,16) take the UT form on SJ (GED 8/16/11); the aspen group (18,20,21) Sheppard's.
                    pothtg = (crvar || sp == 15 || sp == 16) ? ((sj / 5f0) * (sj * 1.5f0 - hb) / (sj * 1.5f0)) * 0.83f0 : 0f0
                    xv = cr / 100f0
                    vigor = 150f0 * xv^3 * fexp(-6f0 * xv) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
                    (sp == 15 || sp == 16) && (vigor = 1f0 - (1f0 - vigor) / 3f0)
                    if sp == 18 || sp == 20 || sp == 21
                        ag1 = (hb * 12f0 * 2.54f0 / 26.9825f0)^0.8509f0
                        h2 = (26.9825f0 * (ag1 + 10f0)^1.1752f0) / (2.54f0 * 12f0)
                        edh = (h2 - hb) * rsimod * rhcon[sp] * 0.75f0
                    else
                        edh = pothtg * pctred * vigor * rhcon[sp]
                    end
                    edh *= 0.5f0
                end
            end
            nivar && (edh = hk - hb)                      # regent.f:1265-1269
            ttvar && (edh = (hk - hb) * rhcon[sp])
            term = hg * scale3
            pr = t.tpa[i]
            snx += edh * pr; sny += term * pr; snp += pr; nh += 1
        end
        nh < 5 && continue                                # NCALHT
        snx /= snp; sny /= snp
        cornew = sny / snx
        cornew <= 0f0 && (cornew = 1f-4)
        (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
        c.htg_cor_init[sp] = flog(cornew)
    end
    return s
end

"""IE REGCON (ie/regent.f:1479): per-species small-tree height constant RHCON."""
function ie_regcons!(s::StandState)
    p = s.plot
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    igl = Int(p.geo_location); (igl < 1 || igl > 3) && (igl = 2)
    asp = p.aspect; slope = p.slope
    regch = IE_RG_RHGL[igl] + (IE_RG_RSAB[1] + IE_RG_RSAB[2]*cos(asp) + IE_RG_RSAB[3]*sin(asp)) * slope
    rhcon = Vector{Float32}(undef, 23)
    @inbounds for sp in 1:23
        # regent.f:1479-1492 — the REGCH+RHSC+RHHAB formula is NIVAR-ONLY; every non-NIVAR (TT/CR/UT)
        # species gets RHCON = 1.0 (their CON = 1.0·EXP(HCOR); the site effect rides in via HCOR calib).
        rc2 = (s.control.regh_cor2_on && sp <= length(s.control.regh_cor2)) ? s.control.regh_cor2[sp] : 0f0
        if sp <= 12 || sp == 14 || sp == 23
            irhhab = clamp(Int(IE_RG_MAPHAB[itype, sp]), 1, 6)
            rhcon[sp] = regch + IE_RG_RHSC[sp] + IE_RG_RHHAB[irhhab, sp]
            rc2 > 0f0 && (rhcon[sp] += log(rc2))          # READCORR: + ln RCOR2
        else
            rhcon[sp] = rc2 > 0f0 ? rc2 : 1.0f0           # READCORR: RHCON = RCOR2
        end
    end
    return rhcon
end


# HTGF per-tripled-copy large-tree height growth (ie/htgf.f:318-347). Each tripled copy's large-tree HTG
# used in the REGENT XWT blend is recomputed from the COPY's OWN diameter growth DG(ITFN)/DG(ITFN+1) —
# NOT the central record's — using the SAME CON (central DBH/HT). NIVAR (recompute) species only; the
# Weibull species (13,17-22) keep the central value (htgf HTG(ITFN)=TEMHTG). `con`,`hdgcof` are the
# central tree's htgf terms; `dg` the copy's large-tree DG; tail = *scale*xht then the HT-based SIZCAP.
@inline function _ie_htgf_copy_large(con::Float32, hdgcof::Float32, dg::Float32,
                                     scale::Float32, xht::Float32, hti::Float32, cap::Float32)
    dg <= 0.0f0 && return 0.1f0 * scale * xht               # ln(DG) undefined ⇒ FVS's pre-tail 0.1 floor
    v = fexp(con + hdgcof * flog(dg)) + IE_HTBIAS
    v < 0.1f0 && (v = 0.1f0)
    v = v * scale * xht
    (hti + v > cap) && (v = max(cap - hti, 0.1f0))          # htgf.f:329-331 SIZCAP on HT(ITFN)=central HT
    return v
end

"""IE `small_tree_growth!` — ie/regent.f's growth section in the Fortran's own shape: one subcycle loop (DO 17 J / DO 16
ISPC / IND1) with the NIVAR/TTVAR/CR/UT arms and the RDNEXT/BANEXT(J+1) density feedback, then one DO-30 assembly with
the 918 tripling loop-back, the DUBSCR crown, and the variables the Fortran carries from tree to tree (BARK, DADJ, the
tripling D). `lestb=true` is REGENT(.TRUE.,ITRNIN), ie_esgent!'s ESTAB entry. It replaced per-species passes that had
TT on TPCCF≈RELDEN (not PCCF·PPCCF) with DLESS3 on it, no CR/PB density feedback, FINDAG's exponent as 0.8509, the
DUBSCR without its TPCCF/AVH/RMAI terms, and no copy crowns. Measured vs FVSie_g16 on the REGCAL fixture
(test/fixtures/inlandempire/regcal): 3 cycles per-tree exact with NOTRIPLE/DGSTDEV 0, tripling, the random component,
and both."""
function small_tree_growth!(s::StandState, stash, ::InlandEmpire; fint::Float32 = 10.0f0,
                            lestb::Bool = false, itrnin::Int = 1,
                            atba::Float32 = -1f0, atccf::Float32 = -1f0, atavh::Float32 = -1f0,
                            ba_now::Float32 = -1f0, relden_now::Float32 = -1f0,
                            pccf_now::Vector{Float32} = Float32[])
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    # PCCF in ESTAB mode = the gradd.f:192 DENSE's (post-growth, PRE-ESNUTR) — the caller's snapshot; jl's density was
    # re-DENSEd with the AUTOES cohort (MEASURED FVSie_g16 3285544010690 thinbba 2012: point-1 PCCF 0.0 live vs 0.18561,
    # point 2 177.69788 vs 178.16805 ⇒ the seedling crown dub CR 0.8950068 → ICR 90 live, 0.8949852 → 89 jl).
    pcv = (lestb && !isempty(pccf_now)) ? pccf_now : dens.point_ccf
    sd = s.coef.species                                 # blkdat HT-DBH :ht1/:ht2 for the aspen log-DK
    t.n == 0 && return s
    n = t.n
    rhcon = ie_regcons!(s)
    # htgf.f HTCONS habitat terms — needed to recompute each tripled COPY's large-tree HTG (its own DG) for
    # the NIVAR XWT blend (ie/htgf.f:318-347). Stand-level; resolved once.
    itype_hg = Int(p.habitat_input); iht_hg = (1 <= itype_hg <= 30) ? IE_HTMAPHAB[itype_hg] : 1
    hghch_hg = IE_HGHC[iht_hg]; h2cof_hg = IE_HGH2[iht_hg]; hdgcof_hg = IE_HGLDD[iht_hg]
    htgf_scale = fint / 10.0f0                              # htgf SCALE=FINT/YR, YR=10 (matches height_growth!)
    # BA/RELDEN = the COMMON values of the last DENSE: the growth pass's; in ESTAB mode the gradd.f:192 DENSE
    # (post-growth, PRE-ESNUTR — the caller snapshots them before the new cohort exists: ba_now/relden_now).
    ba = (lestb && ba_now >= 0f0) ? ba_now : p.basal_area
    relden = (lestb && relden_now >= 0f0) ? relden_now : p.relative_density
    avh = p.avg_height
    dgsd = s.control.dg_sd
    regyr = IE_RG_REGYR
    yr = htg_period(s.variant)                          # /CONTRL/ YR (ie/blkdat.f:57 DATA YR/10.0/), not the TIMEINT cycle length
    ntyr = trunc(Int, fint); iyr = Int(regyr)           # regent.f:198 NTYR=INT(FINT)
    lskiph = false
    if lestb                                            # regent.f:200-201 ESTAB: the rest of the cycle after year 5
        ntyr -= 5; lskiph = ntyr <= 0
    end
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1)
    kper = zeros(Int, 10); itot = ntyr; nn = nper
    @inbounds for k in 1:nper
        nn == 1 && break
        kper[k] = itot ÷ nn; itot -= kper[k]; nn -= 1
    end
    nper > 0 && (kper[nper] = itot)
    banext = zeros(Float32, 11); rdnext = zeros(Float32, 11)    # BANEXT/RDNEXT(10) (+1 for J+1)
    temba = atba > 0f0 ? atba : ba                      # regent.f:218-222 TEMBA/TEMCCF/TEMAHT (after-thin values)
    temccf = atccf > 0f0 ? atccf : relden
    temaht = atavh
    if lestb                                            # regent.f:276-294 (label 8): interpolate from the cycle start
        bayr = 0f0; ccfyr = 0f0
        if !lskiph
            bayr = (ba - temba) / Float32(itot); ccfyr = (relden - temccf) / Float32(itot)
        end
        nyr = 5
        @inbounds for j in 1:nper
            rdnext[j] = temccf + Float32(nyr) * ccfyr; banext[j] = temba + Float32(nyr) * bayr
            nyr += kper[j]
        end
    else
        @inbounds for j in 1:nper; banext[j] = ba; rdnext[j] = relden; end
    end
    if !lestb && nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = ie_bratio(sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            # CCFCAL returns CCFT·P (ccfcal.f:81) ⇒ C1/C2 carry P, and FVS evaluates RDNEXT += ((K*CI)/P)*PN with
            # PN=P*0.985**K (integer power ⇒ __powisf2) — regent.f:262-273. Cancelling the P algebraically (the old
            # `k*ci*pn` on bare CCFT) rounds differently: IE cycle-4 subcycle-2 RDJ was 2 ULP off (bare-plot fixture).
            c1 = ie_tree_ccf(sp, d1) * pr; c2 = ie_tree_ccf(sp, d2) * pr
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]
                pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn
                banext[j] += Float32(k) * bi * pn
            end
        end
    end
    if lestb                                            # regent.f:301-319 DO 13: crown for each new record, STORAGE order
        @inbounds for i in itrnin:n
            ipc = Int(t.plot_id[i])
            pcc = (1 <= ipc <= length(pcv)) ? pcv[ipc] : 0f0
            crn = 0.89722f0 - 0.0000461f0 * pcc
            ran = 0f0
            while true
                ran = bachlo(s.rng, 0f0, 1f0); (ran < -1f0 || ran > 1f0) && continue
                break
            end
            crn = crn + 0.07985f0 * ran
            crn > 0.90f0 && (crn = 0.90f0); crn < 0.20f0 && (crn = 0.20f0)
            t.crown_pct[i] = trunc(Int32, (crn * 100f0) + 0.5f0)
        end
    end
    # DELMAX and the CR/UT PCTRED (regent.f:326-345): AH/R = TEMAHT/TEMCCF from ESTAB, else AVH/RELDEN
    ah = lestb ? temaht : avh; r = lestb ? temccf : relden
    delmax = (ah / 36.0f0) * (0.01232f0 * r - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    xpr = ah * (r / 100.0f0); xpr > 300.0f0 && (xpr = 300.0f0)
    pctred = 1.11436f0 + xpr*(-0.011493f0 + xpr*(0.43012f-4 + xpr*(-0.72221f-7 +
             xpr*(0.5607f-10 - xpr*0.1641f-13))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # per-tree height + diameter accumulators (WK3=H, WK5=D), start at HT/DBH
    wk3 = Float32[t.height[i] for i in 1:n]
    wk5 = Float32[t.dbh[i] for i in 1:n]
    # TTVAR (sp13/17) ZRAND lives in the TREE RECORD (regent.f:513 ZRAND(I)), not in this call: it is seeded at
    # input (intree.f:369/607), re-seeded by the cycle-0 LHTCAL calibration pass (regent.f:1181) and carried
    # across cycles by TRIPLE/tremov, so a tree keeps its deviate until EDH<=0.1 resets it to -999.
    zrand_tt = t.zrand
    cur_year = current_cycle_year(s)
    # WK4(I) = clgmult's per-tree CLIMATE growth multiplier (regent.f:596/598 H2=H1+…·WK4). In the oracle
    # CLGMULT fills WK4 once (dgdriv.f:153) and BOTH the large-tree DDS (dgdriv.f:217, jl apply_climate_dds!)
    # AND the small-tree regent height read it. jl applied it to diameter but OMITTED it in the height path,
    # so under a suppressing climate score the NIVAR/TTVAR seedling over-grew ~1/WK4 (tripping breast-height a
    # cycle early ⇒ escaping self-thinning — the IE climate-regime dense-phase under-kill). WK4≡1 without a
    # CLIMATE keyword (clgmult.f:55,60 return early) ⇒ IEEE-exact no-op on every non-climate stand.
    wk4 = ones(Float32, n)
    if s.climate !== nothing && s.climate.active
        _c = s.climate; _cd = _c.data; _ix = _c.indices
        if !(_ix[:mtcm]==0 || _ix[:mmin]==0 || _ix[:dd0]==0 || _ix[:d100]==0 || _ix[:dd5]==0 || _ix[:gsp]==0)
            _ty = Float32(cur_year) + fint / 2f0                 # THISYR = IY(ICYC)+FINT/2 (matches apply_climate_dds!)
            _A(sym, yr) = algslp(yr, _cd.years, view(_cd.attrs, :, _ix[sym]))
            _smi(yr) = (g = _A(:gsp, yr); g > 0f0 ? _A(:dd5, yr) / g : 0f0)
            _xgsite = _ix[:pSite] > 0 ? clim_xgsite(_A(:pSite, _ty), _A(:pSite, Float32(_c.inv_year))) : 1f0
            _mtcm_now = _A(:mtcm, _ty); _mmin_now = _A(:mmin, _ty); _smi_now = _smi(_ty)
            _ns = length(_c.plant_symbols)
            _vscore = ones(Float32, _ns)
            @inbounds for sp in 1:_ns; _, _vscore[sp] = species_vscore(_cd, _c.plant_symbols[sp], _ty); end
            @inbounds for i in 1:n
                t.dbh[i] <= 0f0 && continue
                sp = Int(t.species[i]); (sp < 1 || sp > _ns) && continue
                _by = _ty - t.birth_age[i]
                _xdf = leites_xdf(_mtcm_now, _A(:mtcm, _by))
                _xwl = leites_xwl(_mmin_now, _A(:mmin, _by), _A(:dd0, _by))
                _xpp = leites_xpp(_smi_now, _smi(_by), _A(:d100, _by))
                _xr = clim_xrelgr(_c.plant_symbols[sp], _xdf, _xpp, _xwl)
                _, _tm = clim_treemult(_xgsite, _xr, _vscore[sp], _c.growmult[sp])
                wk4[i] = _tm
            end
        end
    end
    # ESTAB: WK4 of a new record is the establishment height multiplier HTIMLT (esgent.f; 1 for PLANT, <1 AUTOES).
    lestb && @inbounds(for i in itrnin:n; wk4[i] = t.htimlt[i]; end)
    # FVS walks every REGENT tree loop species-major over IND1 (DO ISPC=1,MAXSP / DO I3=I1,I2 / I=IND1(I3)):
    # the subcycle density feedback RDNEXT/BANEXT(J+1) (regent.f:665-670) is a Float32 running sum, so its
    # accumulation ORDER is part of the result; record order drifted RDJ/BAJ by ULPs from subcycle 2 on.
    _isct = s.control.sp_count_tab; _ind1 = s.scratch.idx1
    _sp_order = Vector{Int}(undef, n); _no = 0
    @inbounds for sp in 1:MAXSP
        i1 = _isct[sp, 1]; i1 == 0 && continue
        i2 = _isct[sp, 2]
        for k in i1:i2
            (1 <= k <= length(_ind1)) || continue
            ii = Int(_ind1[k]); (1 <= ii <= n) || continue
            _no += 1; _sp_order[_no] = ii
        end
    end
    # Fortran variables that CARRY from tree to tree (and across REGENT's two loops): BARK (regent.f:444 subcycle
    # loop; :975/:1004 assembly) feeds the next CR/UT DGK=(DK−DKK)·BARK; DADJ (:463, :866) feeds NIVAR D1; D is
    # the tripling-carried diameter (:746 D=DBH(I), :750 TT D=max(D,DIAM), :1003 D=DBH(K)).
    bark_c = 0f0
    dadj_c = 0f0
    pccf_of(i) = (pt = Int(t.plot_id[i]); (1 <= pt <= length(pcv)) ? pcv[pt] : 0f0)
    # ---- subcycle loop (regent.f:343-679): DO 17 J / DO 16 ISPC / DO 15 I3=I1,I2 (species-major IND1) ----
    ky = 0
    @inbounds for j in 1:nper
        ky += kper[j]
        baj = banext[j]; rdj = rdnext[j]
        ppccf = relden > 0.0f0 ? 1.0f0 + (rdj - relden) / relden : 0.0f0     # regent.f:360-361
        surv = fpowi(0.985f0, ky)                              # 0.985**KY: integer power ⇒ libgcc __powisf2
        for _oi in 1:_no
            i = _sp_order[_oi]
            sp = Int(t.species[i])
            niv = sp <= 12 || sp == 14 || sp == 23; ttv = sp == 13 || sp == 17; crv = sp == 19 || sp == 22
            utv = !(niv || ttv || crv)
            scale = (niv || ttv) ? Float32(kper[j]) / regyr : Float32(ntyr) / yr    # regent.f:384-403
            (ttv || crv) && lskiph && continue                 # regent.f:405
            (crv || utv) && j > 1 && continue                  # regent.f:408 CR/UT: one pass only
            xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
            xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
            con = niv ? rhcon[sp] + c.htg_cor_small[sp] : rhcon[sp] * fexp(c.htg_cor_small[sp])   # regent.f:412-416
            d = t.dbh[i]; h = t.height[i]
            d >= IE_RG_XMAX[sp] && continue                    # regent.f:436
            lestb && i < itrnin && continue                    # regent.f:438
            pr = t.tpa[i]
            h1 = wk3[i]
            # PCT of a new ESTAB record is 0 (estab.f:1252; no DENSE before ESGENT)
            bal = baj * (100.0f0 - (lestb ? 0f0 : t.crown_ratio[i])) * 0.0001f0
            c1 = ie_tree_ccf(sp, d) * pr                       # regent.f:443 CCFCAL(ISPC,D,H1,…) = CCFT·P
            b1 = 0.005454154f0 * d * d
            bark_c = ie_bratio(sp, d)                          # regent.f:445 BARK=BRATIO(ISPC,D,H)
            rsimod = 0f0; pothtg = 0f0
            if crv || utv                                      # regent.f:447-474
                si = p.sp_site_index[sp]
                si > IE_RG_SHI[sp] && (si = IE_RG_SHI[sp])
                si <= IE_RG_SLO[sp] && (si = IE_RG_SLO[sp] + 0.5f0)
                relsi = (si - IE_RG_SLO[sp]) / (IE_RG_SHI[sp] - IE_RG_SLO[sp])
                rsimod = 0.5f0 * (1.0f0 + relsi)
                sj = p.sp_site_index[sp]
                (crv || sp == 15 || sp == 16) && (pothtg = ((sj / 5.0f0) * (sj * 1.5f0 - h) / (sj * 1.5f0)) * 0.83f0)
            end
            pccf = pccf_of(i)
            tpccf = pccf * ppccf; tpccf > 300.0f0 && (tpccf = 300.0f0); tpccf < 25.0f0 && (tpccf = 25.0f0)
            cr = Float32(t.crown_pct[i])
            htgrl = 0f0
            if niv                                             # regent.f:487-501
                relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h1 - 4.5f0) / (ah - 4.5f0)
                relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
                dadj_c = delmax * relh * relh - 2.0f0 * delmax * relh + 0.65f0
                htgrl = con + IE_RG_RHLH[sp] * flog(h1) + IE_RG_RHCCF[sp] * rdj + IE_RG_RHBAL[sp] * bal
            elseif ttv                                         # regent.f:505-538: persistent per-tree ZRAND
                if t.zrand[i] == -999f0
                    while true
                        z = bachlo(s.rng, 0.0f0, 1.0f0); (z < -2.0f0 || z > 2.0f0) && continue
                        t.zrand[i] = z; break
                    end
                end
                if d <= 0.0f0
                    htgrl = 0.0f0
                else
                    beta1 = fexp(1.17527f0 - 0.42124f0 * flog(tpccf))
                    beta2 = fexp(-2.56002f0 - 0.58642f0 * flog(tpccf))
                    htg1 = beta1 + beta2 * cr
                    stddev = htg1 * (1.08720f0 - 0.00230f0 * cr)
                    htgrl = htg1 + t.zrand[i] * stddev
                    htgrl <= 0.1f0 && (htgrl = 0.1f0; t.zrand[i] = -999f0)
                end
            else                                               # CR/UT (regent.f:542-594)
                xv = cr / 100.0f0
                vigor = (150.0f0 * (xv^3.0f0) * fexp(-6.0f0 * xv)) + 0.3f0
                vigor > 1.0f0 && (vigor = 1.0f0)
                (sp == 15 || sp == 16) && (vigor = 1.0f0 - ((1.0f0 - vigor) / 3.0f0))
                if sp == 18 || sp == 20 || sp == 21            # Sheppard aspen; FINDAG (ie/findag.f:17)
                    sitage = lestb ? Float32(t.birth_age[i]) : (h * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0)
                    hite1 = 26.9825f0 * sitage^1.1752f0
                    hite2 = 26.9825f0 * (sitage + 10.0f0)^1.1752f0
                    htgrl = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con
                    htgrl = htgrl * 0.75f0
                else
                    htgrl = pothtg * pctred * vigor * con
                end
            end
            h2 = niv ? h1 + fexp(htgrl) * scale * xrhgro * wk4[i] :              # regent.f:599-605
                 ttv ? h1 + htgrl * scale * xrhgro * con * wk4[i] : h1 + htgrl * scale
            wk3[i] = h2
            # subcycle diameter + density feedback (regent.f:613-676)
            d2 = d
            if niv
                (j >= nper || d >= 3.0f0 || h2 <= 4.5f0) && continue          # GO TO 14 (no WK5, no density)
                ax = IE_RG_HHT1[sp]; bx = IE_RG_HHT2[sp]
                d1v = IE_RG_DIAM[sp] + dadj_c
                h1 > 4.5f0 && (d1v = ax * fpow(h1 - 4.5f0, bx) + dadj_c)
                d2v = ax * fpow(h2 - 4.5f0, bx) + dadj_c
                dgj = (d2v - d1v) * xrdgro; dgj < 0.0f0 && (dgj = 0.0f0)
                d2 = d + dgj
            elseif ttv
                hless4 = h2 - 4.5f0                                                # DLESS3 on the POINT PCCF
                dless3 = 0.000231f0 * hless4 * cr - 0.00005f0 * hless4 * pccf + 0.001711f0 * cr + 0.17023f0 * hless4
                d2 = (dless3 + 0.3f0) * xrdgro
                d2 < IE_RG_DIAM[sp] && (d2 = IE_RG_DIAM[sp])
                wk5[i] = d2
                (j >= nper || d >= 3.0f0 || h2 <= 4.5f0) && continue
            elseif crv
                d >= 1.0f0 && continue                                          # GO TO 15
                h2 <= 4.5f0 && (d2 = d + 0.0001f0 * h2)
            else
                # regent.f:648 `ISPC.EQ.18 .OR. ISPC.EQ.20 .OR. ISPC.EQ.21 .AND. D.GE.3` = 18 ∨ 20 ∨ (21 ∧ D≥3)
                (sp == 18 || sp == 20 || (sp == 21 && d >= 3.0f0)) && continue
                if h2 <= 4.5f0
                    d2 = d + 0.001f0 * h2
                elseif sp == 15 || sp == 16
                    d2 = (h2 - 4.5f0) * 10.0f0 / (p.sp_site_index[sp] - 4.5f0); d2 < 0.1f0 && (d2 = 0.1f0)
                end
            end
            wk5[i] = d2
            c2 = ie_tree_ccf(sp, d2) * pr
            rdnext[j+1] += Float32(ky) * (c2 - c1) / 10.0f0 * surv
            banext[j+1] += (0.005454154f0 * d2 * d2 - b1) * pr * surv
        end
    end
    # ---- assembly (regent.f DO 30 ISPC / DO 25 I3, :683-1047) ----
    rmai_v = _ie_rmai(s)
    scale30 = Float32(ntyr) / yr                       # DO-30 SCALE = NTYR/YR (every group)
    scale2 = Float32(yr) / Float32(ntyr)               # SCALE2 = YR/NTYR (!LSKIPH)
    ltrip = !lestb && stash !== nothing                # regent.f:1026 no tripling from ESTAB
    @inbounds for _oi in 1:_no
        i = _sp_order[_oi]
        sp = Int(t.species[i])
        niv = sp <= 12 || sp == 14 || sp == 23; ttv = sp == 13 || sp == 17; crv = sp == 19 || sp == 22
        utv = !(niv || ttv || crv)
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        xmx = IE_RG_XMAX[sp]; xmn = IE_RG_XMIN[sp]; diam = IE_RG_DIAM[sp]
        ax = IE_RG_HHT1[sp]; bx = IE_RG_HHT2[sp]
        dgmx = ttv ? fint * IE_RG_DGMAX[sp] : crv ? IE_RG_DGMAX[sp] * fint / 10.0f0 : IE_RG_DGMAX[sp]
        d = t.dbh[i]; cr = Float32(t.crown_pct[i])
        d >= xmx && continue                                           # regent.f:748
        lestb && i < itrnin && continue                                # regent.f:749
        h = t.height[i]; hk = wk3[i]
        dk = wk5[i]
        htgr1 = hk - h
        ttv && d < diam && (d = diam)                                  # regent.f:756 IF(TTVAR.AND.D.LT.D1)D=D1
        if ttv && lskiph                                               # regent.f:765-777 (ESTAB, FINT<=5)
            if h >= 4.5f0
                pccf_l = pccf_of(i); hl = h - 4.5f0
                dkl = (0.000231f0*hl*cr - 0.00005f0*hl*pccf_l + 0.001711f0*cr + 0.17023f0*hl + 0.3f0) * xrdgro
                dkl < d && (dkl = d)
                t.dbh[i] = dkl; t.diam_growth[i] = dkl > dgmx ? dgmx : dkl
            end
            continue
        end
        niv && htgr1 < 0.0f0 && (htgr1 = 0.0f0)                        # regent.f:791-794
        icrk = Int(t.crown_pct[i]); icr0 = icrk
        dbh0 = t.dbh[i]                                                # DBH(K) of every slot before REGENT
        pccf = pccf_of(i)
        cap = s.control.sp_size_cap[sp, 4]
        # the large-tree HTG(K) each record blends toward (XWT): central = its htgf value; NIVAR copies = htgf's
        # per-copy value from the copy's own DG (ie/htgf.f:318-347); Weibull species keep the central (TEMHTG).
        large0 = t.ht_growth[i]; largeU = large0; largeL = large0
        if ltrip && niv
            htcon_hg = hghch_hg + IE_HGSC[sp]
            (s.control.htg_cor2_on && s.control.htg_cor2[sp] > 0.0f0) && (htcon_hg += flog(s.control.htg_cor2[sp]))
            con_hg = htcon_hg + h2cof_hg * h * h + IE_HGLD[sp] * flog(dbh0) + IE_HGLH * flog(h)
            xht_hg = active_multiplier(s.control, :htg, sp, cur_year)
            largeU = _ie_htgf_copy_large(con_hg, hdgcof_hg, stash.dgU[i], htgf_scale, xht_hg, h, cap)
            largeL = _ie_htgf_copy_large(con_hg, hdgcof_hg, stash.dgL[i], htgf_scale, xht_hg, h, cap)
        end
        nrec = (ltrip && !ttv) ? 3 : 1
        for l in 0:(nrec - 1)
            htgr = htgr1
            if lskiph                                                  # regent.f:763-780: HTG(K)=0, GO TO 20
                htgr = 0f0
            elseif niv                                                 # regent.f:803-809 (label 18)
                zz = 0.0f0
                while true
                    zz = dgsd >= 1.0f0 ? bachlo(s.rng, 0.0f0, 1.0f0) : 0.0f0
                    (zz > 1.0f0 || zz < -1.5f0) && continue
                    break
                end
                htgr = htgr1 * fexp(zz * IE_RG_HSIGMA)
            elseif !ttv                                                # regent.f:811-825 (label 919)
                zz = 0.0f0
                while true
                    zz = dgsd >= 1.0f0 ? bachlo(s.rng, 0.0f0, 1.0f0) : 0.0f0
                    (zz > 0.5f0 || zz < -2.0f0) && continue
                    break
                end
                utv && (htgr = (htgr1 + zz * 0.1f0) * xrhgro)
                crv && (htgr = (htgr1 + zz * 0.2f0) * xrhgro)
                htgr < 0.1f0 && (htgr = 0.1f0)
            end
            xwt = (d <= xmn || lestb) ? 0.0f0 : (d - xmn) / (xmx - xmn)   # regent.f:833-835 (carried D; 0 from ESTAB)
            largeK = l == 0 ? large0 : (l == 1 ? largeU : largeL)
            htgk = htgr * (1.0f0 - xwt) + xwt * largeK
            if !lskiph
                (crv && htgk < 0.1f0) && (htgk = 0.1f0)
                if h + htgk > cap                                      # regent.f:844-847 SIZCAP
                    htgk = cap - h; htgk < 0.1f0 && (htgk = 0.1f0)
                end
            end
            # label 20 (regent.f:856-859): CR D≥1 ⇒ 23; `UTVAR.AND.ISPC.EQ.15 .OR. ISPC.EQ.16` ⇒ 1020; D≥3 ⇒ 23
            skip23 = (crv && d >= 1.0f0) || (!((utv && sp == 15) || sp == 16) && d >= 3.0f0)
            dbhk = dbh0; dgk_out = 0f0; dbh_set = false; crK = -1
            if !skip23
                d1v = diam
                if niv                                                 # regent.f:863-874
                    relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h - 4.5f0) / (ah - 4.5f0)
                    relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
                    dadj_c = delmax * relh * relh - 2.0f0 * delmax * relh + 0.65f0
                    d1v = diam + dadj_c
                    h > 4.5f0 && (d1v = ax * fpow(h - 4.5f0, bx) + dadj_c)
                end
                hk = h + htgk
                if !ttv && hk < 4.5f0                                   # regent.f:878-882
                    niv && (dbhk = 0.1f0 + diam * 0.01f0 + hk * 0.001f0; dbh_set = true)
                    crv && (dbhk = d + 0.001f0 * hk; dbh_set = true)
                    dgk_out = 0.0f0
                else
                    dkk = 0f0; dgk = 0f0
                    if niv
                        dk = ax * fpow(hk - 4.5f0, bx) + dadj_c
                        dk < diam && (dk = diam)
                        dk = dk + hk * 0.001f0
                    elseif crv || utv
                        if sp == 15 || sp == 16
                            sitear = p.sp_site_index[sp]
                            dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                            dkk = (h - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                            h < 4.5f0 && (dkk = d)
                        else
                            bxw = sd[:ht2][sp]                                            # blkdat HT2
                            axw = c.ht_dbh_iabflg[sp] == 1 ? sd[:ht1][sp] : c.ht_dbh_aa[sp]
                            dk = (bxw / (flog(hk - 4.5f0) - axw)) - 1.0f0; dk < 0.1f0 && (dk = 0.1f0)
                            dkk = h <= 4.5f0 ? d : (bxw / (flog(h - 4.5f0) - axw)) - 1.0f0
                        end
                    else                                               # TTVAR (regent.f:915-929)
                        if hk >= 4.5f0
                            hless4 = h - 4.5f0
                            dless3 = 0.000231f0 * hless4 * cr - 0.00005f0 * hless4 * pccf + 0.001711f0 * cr + 0.17023f0 * hless4
                            dkk = dless3 + 0.3f0; dkk < diam && (dkk = diam)
                            dgk = (dk - dkk) * xrdgro
                        else
                            dgk = 0.0f0
                        end
                    end
                    if lestb                                           # regent.f:931-940 ESTAB diameter
                        dbhk = dk; dbh_set = true
                        dg = ttv ? dgk : 0f0
                        (ttv && dg > dgmx) && (dg = dgmx)
                        (niv || crv) && (dg = dk)
                        if crv || utv
                            dbhk < diam && (dbhk = diam)
                            dbhk = dbhk + 0.001f0 * hk
                            dg = dbhk
                        end
                        (dbhk + dg) < diam && (dg = diam - dbhk)
                        dgk_out = dg
                    else
                    if niv                                             # regent.f:953-966
                        dgk = (dk - d1v) * xrdgro
                    elseif crv || utv
                        if dk < 0.0f0 || dkk < 0.0f0
                            dgk = htgk * 0.2f0 * bark_c * xrdgro        # carried BARK
                            dk = d + dgk
                        else
                            dgk = (dk - dkk) * bark_c * xrdgro          # carried BARK
                        end
                        dgk > dgmx && (dgk = dgmx)
                    end
                    dgk < 0.0f0 && (dgk = 0.0f0)
                    bark_c = ie_bratio(sp, dbhk)                       # regent.f:975 BRATIO(ISPC,DBH(K),HT(K))
                    dg = (niv || ttv) ? dgk * bark_c : dgk
                    dds = dg * (2.0f0 * bark_c * d + dg) * scale2
                    dg = sqrt((d * bark_c)^2.0f0 + dds) - bark_c * d
                    (dbhk + dg) < diam && (dg = diam - dbhk)
                    dgk_out = dg
                    end
                end
                dgk_out = dg_bound(nothing, nothing, sp, dbhk, dgk_out, s.control.sp_size_cap)   # DGBND
                # DUBSCR crown (regent.f:998-1021, not from ESTAB): D=DBH(K) carries into the next tripled record
                if !lestb
                d = dbhk
                bark_c = ie_bratio(sp, d)
                dds2 = (dgk_out * (2.0f0 * bark_c * d + dgk_out)) * scale30
                dg2 = sqrt((d * bark_c)^2 + dds2) - bark_c * d
                dg2 < 0.0f0 && (dg2 = 0.0f0)
                dnew = d + dg2
                if dnew >= 3.0f0
                    if icrk == 0
                        icrk = ie_dubscr(s.rng, sp, dnew, hk, ba, dgsd; tpccf = pccf, avh = avh, tmai = rmai_v)
                    end
                    crK = icrk
                end
                end
            end
            if l == 0
                t.ht_growth[i] = htgk
                if !skip23
                    t.diam_growth[i] = dgk_out
                    dbh_set && (t.dbh[i] = dbhk)
                end
                crK >= 0 && (t.crown_pct[i] = Int32(crK))
            else
                crc = crK >= 0 ? crK : icr0                            # the slot keeps ICR(I) unless DUBSCR set it
                if l == 1
                    stash.htgU[i] = htgk; stash.is_small[i] = true
                    skip23 || (stash.dgU[i] = dgk_out; stash.dbhU[i] = dbhk)
                    stash.crU[i] = Int32(crc)
                else
                    stash.htgL[i] = htgk
                    skip23 || (stash.dgL[i] = dgk_out; stash.dbhL[i] = dbhk)
                    stash.crL[i] = Int32(crc)
                end
            end
        end
        # TTVAR tripling (regent.f:1026-1037): the two copies are exact duplicates of the central record
        if ltrip && ttv
            stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]; stash.is_small[i] = true
            stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        end
    end
    return s
end

# ie/maical.f RMAI (MAI of the site species, capped 128) — the DUBSCR TMAI argument of REGENT. ie/crown.f:180
# has `IF(RMAILM .GT. 128.0)RMAI=128.0` (it caps the COMMON RMAI, not RMAILM), which CROWN runs in CRATET after
# MAICAL, so an LM site giving RMAILM>128 leaves RMAI=128.
function _ie_rmai(s::StandState)::Float32
    isisp = Int(s.plot.site_species); isisp == 0 && (isisp = 3)
    (isisp < 1 || isisp > 23) && (isisp = 3)
    sssi = s.plot.sp_site_index[isisp]; sssi == 0f0 && (sssi = 140f0)
    r = _adjmai(IE_MAI_ISPNUM[isisp], sssi, 10f0); r > 128f0 && (r = 128f0)
    _adjmai(101, s.plot.sp_site_index[13], 10f0) > 128f0 && (r = 128f0)
    return r
end
const IE_MAI_ISPNUM = Int[119, 73, 202, 17, 263, 242, 108, 93, 19, 122, 264,
                          101, 101, 101, 101, 101, 101, 746, 740, 746, 375, 998, 101]

# ie/esgent.f (CALL REGENT(.TRUE.,ITRNIN)) — grow the JUST-ESTABLISHED regen IN its birth cycle (#186).
function ie_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atavh::Float32 = -1.0f0, atba::Float32 = -1.0f0, atrelden::Float32 = -1.0f0,
                    relden_pre::Float32 = -1.0f0, ba_pre::Float32 = -1.0f0, pccf_pre::Vector{Float32} = Float32[])
    t = s.trees
    nstart >= t.n && return s
    # ie/esgent.f (== em/esgent.f): SPESRT, then REGENT(.TRUE.,ITRNIN) — the same Fortran-shaped REGENT as the growth
    # pass in its ESTAB mode (DO-13 crown per new record in storage order; FINT−5 years of subcycling with the density
    # interpolated from the after-thin TEMBA/TEMCCF/TEMAHT; XWT=0; the ESTAB diameter; no tripling or DUBSCR) for
    # EVERY sub-model. jl had grown only the NIVAR conifers here, so planted/natural LM, PI/JU, aspen/PB and CO/OH
    # got no birth-cycle growth and their draws were missing from the stream (IE REGCAL PLANT fixture: planted AS
    # DBH 1.14 vs live 0.78). BA/RELDEN are the gradd.f:192 DENSE's (post-growth, PRE-ESNUTR: ba_pre/relden_pre).
    species_sort!(s)                                   # esgent.f CALL SPESRT
    small_tree_growth!(s, nothing, s.variant; fint = fint, lestb = true, itrnin = nstart + 1,
                       atba = atba, atccf = atrelden, atavh = atavh, ba_now = ba_pre, relden_now = relden_pre,
                       pccf_now = pccf_pre)
    @inbounds for i in (nstart+1):t.n
        sp = Int(t.species[i])
        htemp = t.height[i] + t.ht_growth[i]
        wk4 = t.htimlt[i]
        t.ht_growth[i] = t.ht_growth[i] * wk4
        t.height[i] = t.height[i] + t.ht_growth[i]
        if wk4 < 1f0
            if t.height[i] < 4.5f0
                t.dbh[i] = 0.1f0 + 0.001f0 * t.height[i]; t.diam_growth[i] = 0f0
            else
                t.dbh[i] = t.dbh[i] * (t.height[i] / htemp)
                t.diam_growth[i] = t.diam_growth[i] * (t.height[i] / htemp)
            end
        end
        if t.height[i] > _IE_ES_HHTMAX[sp]
            t.height[i] = _IE_ES_HHTMAX[sp]; t.dbh[i] = 2.95f0
        end
    end
    esgent_add_gentim!(s, nstart, fint)                # estab.f:1504 ABIRTH += GENTIM (after ESGENT)
    return s
end

