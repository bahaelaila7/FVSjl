# =============================================================================
# height_growth.jl (utah) — UT large-tree height growth (ut/htgf.f). Chunk 4.
#
# Conifer/aspen (CASE DEFAULT: sp 1-10,23) = Schreuder-Hafley SBB (same model as TT). utt01 species
# (WB/WF/ES/LP/AS, all sp≤10) use this. UT specifics vs TT:
#   K=(JSPC-1)*5 + (MPCRSI? LSIMAP : KEYCR); JSPC = sp≤10 ? (sp==5?8:sp) : 11;
#   KEYCR: IICR≤2→1, ≤7→2, ≤9→3; LSIMAP=clamp(INT(SITEAR(ISISP)/10−0.5),2,6)−1 (DF/BS/ES/AF, MPCRSI=1);
#   ZBIAS=AZBIAS+BZBIAS·(ELEV−20), active only ELEV∈[80,105], then UNCONDITIONAL ZADJ=.1−.10273Z+.00273Z²;
#   DIA=DBH+DG/BARK; PSI=COF8·((DIA−XI1)/(XI1+COF1−DIA))^COF9·exp(Z·√(1−COF7²)/COF6); H=(PSI/(1+PSI))·COF2+XI2;
#   HTG=max(H−HT,0.1); BS(5)×0.95; ·SCALE·XHMULT(=1)·exp(HTCON)·MISHGF(=1). Small tree (D<1) → REGENT (HTG 0).
# PJ 11:16,24 → REGENT (HTG 0). CR-surrogate 17:19,22 (GB/NC/FC/BE) + MC/BI 20:21 now ported (below).
# Young-LP accelerator (ut/htgf.f:672-707) bypassed (ICYC>1 | IAGE≤0); ISTAGF=0 (grinit).
#
# CR-surrogate CASE(17:19,22) (ut/htgf.f:308-432) — even-aged HHE (Alexander RM-32 spruce-fir site curve
#   IMODTY=4) meshed with uneven-aged HHU (GENGYM diameter curve) via the AGERNG>40 & BA≥70 blend:
#     AP = ABIRTH − 4.5/(−0.22+0.0155·max(TSITE,20));  HHE = 2.75780·TSITE^0.83312·(1−e^(−0.015701·AGETEM))
#       ^(22.71944·TSITE^−0.63557) + 4.5  (AGETEM≥30; if AP<AGETEM linearly interp; ×max(1−BAUTBA,0.728));
#     HHU: sp17 = 42.269377·(1−e^(−0.165687·D))^1.184734 + 4.5 (pinyon/GENGYM);
#          sp18,19,22 = (−2.04+1.4534·TSITE)·(1−e^(−0.058112·D))^(1.894400·BATEM^−0.192979) + 4.5, BATEM=max(BA,10);
#     twin calls at D (current) and DFUT=D+DG/BARK (future age AP+10); HTG=(HHE2−HHE1)·ADJUST;
#     ADJUST: sp17=1.0, sp22=1.3−0.02·(SI−15), else 0.78+0.0023·SI.  + ZZRAN·0.1 (DGSD>0 only), floor 0.1.
#     ABIRTH dubbed at setup via FNDAG (=cr_fndag IMODTY 4, spruce-fir), aged FINT/cycle (gradd.f:205).
# MC/BI CASE(20,21) (ut/htgf.f:439-580) — SO-variant: FINDAG(H)→SITAGE/SITHT/HTMAX/AGMAX, then Chapman-Richards
#   POTHTG = HGUESS(SITAGE+10) − SITHT; HTG = POTHTG·HTGMOD, HTGMOD = .25·HGMDCR + .75·HGMDRH (crown Hoerl +
#   relative-height gen. Chapman-Richards, RH coeffs per sp). H≥HTMAX ⇒ HTG=0.1 (double-scaled). AVH=p.avg_height.
# =============================================================================

@inline _ut_ht_crsurr(sp::Int) = sp == 17 || sp == 18 || sp == 19 || sp == 22

function height_growth!(s::StandState, ::Utah; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    elev = p.elevation
    zon = (elev >= 80f0 && elev <= 105f0)
    isisp = Int(p.site_species); (isisp < 1 || isisp > 24) && (isisp = 7)
    lsimap = clamp(trunc(Int, p.sp_site_index[isisp] / 10f0 - 0.5f0), 2, 6) - 1
    ba = p.basal_area; avh = p.avg_height; dgsd = s.control.dg_sd
    # CR-surrogate (17:19,22) needs the BADIST BA-above-class (BAUTBA) + stand age range (even/uneven blend).
    has_crsurr = any(j -> _ut_ht_crsurr(Int(t.species[j])), 1:t.n)
    ht_bau = has_crsurr ? _cr_badist_bau(t) : nothing
    agerng = has_crsurr ? _cr_agerng(t) : 0f0
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (d <= 0f0 || h <= 0f0) && continue
        if sp in (11,12,13,14,15,16,24)
            continue                                     # PJ/woodland: height from REGENT (HTG 0)
        elseif _ut_ht_crsurr(sp)
            # ut/htgf.f CASE(17:19,22) small gate: DBH<0.5 | HT≤4.5 → REGENT (HTG stays 0).
            (d < 0.5f0 || h <= 4.5f0) && continue
            ssite = p.sp_site_index[sp]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            icls = trunc(Int, d + 1f0); icls > 41 && (icls = 41); icls < 1 && (icls = 1)
            bautba = (ht_bau !== nothing && ba > 0f0) ? ht_bau[icls] / ba : 0f0
            bautba < 0f0 && (bautba = 0f0)
            zzr = 0f0
            if dgsd > 0f0
                while true
                    zzr = bachlo(s.rng, 0f0, 1f0)
                    (zzr <= dgsd && zzr >= -dgsd) && break
                end
            end
            t.ht_growth[i] = _ut_htg_crsurr(sp, ssite, d, h, t.diam_growth[i], bark, bautba,
                t.birth_age[i], agerng, ba, Float32(t.crown_ratio[i]), c.htg_cor[sp], scale,
                s.control.sp_size_cap[sp, 4], zzr)
            continue
        elseif sp == 20 || sp == 21
            # ut/htgf.f CASE(20,21) small gate = DEFAULT: DBH<1.0 → REGENT (HTG stays 0).
            d < 1.0f0 && continue
            ssite = p.sp_site_index[sp]
            t.ht_growth[i] = _ut_htg_mcbi(sp, ssite, d, h, avh, Float32(t.crown_pct[i]),
                c.htg_cor[sp], scale, s.control.sp_size_cap[sp, 4])
            continue
        end
        d < 1.0f0 && continue                            # small tree → REGENT (chunk 6); HTG stays 0
        # SBB coefficient row (ut/htgf.f CASE DEFAULT)
        iicr = trunc(Int, Float32(t.crown_pct[i]) / 10f0 + 0.5f0); iicr > 9 && (iicr = 9)
        keycr = iicr <= 2 ? 1 : (iicr <= 7 ? 2 : 3)
        jspc = sp <= 10 ? (sp == 5 ? 8 : sp) : 11
        col = UT_MPCRSI[sp] == 1 ? lsimap : keycr
        k = (jspc - 1) * 5 + col
        cof1 = UT_HTCOF[k, 1]; cof2 = UT_HTCOF[k, 2]; cof3 = UT_HTCOF[k, 3]
        cof4 = UT_HTCOF[k, 4]; cof5 = UT_HTCOF[k, 5]; cof6 = UT_HTCOF[k, 6]
        cof7 = UT_HTCOF[k, 7]; cof8 = UT_HTCOF[k, 8]; cof9 = UT_HTCOF[k, 9]
        htg = 0.1f0
        if h > 4.5f0 && (_UT_XI1 + cof1) > d && (_UT_XI2 + cof2) > h && d > 0.1f0
            y1 = (d - _UT_XI1) / cof1
            y2 = (h - _UT_XI2) / cof2
            fby1 = log(y1 / (1f0 - y1))
            fby2 = log(y2 / (1f0 - y2))
            z = (cof4 + cof6 * fby2 - cof7 * (cof3 + cof5 * fby1)) * (1f0 - cof7 * cof7)^(-0.5f0)
            zbias = zon ? UT_AZBIAS[sp] + UT_BZBIAS[sp] * (elev - 20f0) : 0f0
            (z - zbias >= 2f0 && zbias < 0f0) && (zbias = 0f0)
            z -= zbias
            zadj = 0.1f0 - 0.10273f0 * z + 0.00273f0 * z * z       # UT: unconditional (ut/htgf.f:661)
            zadj < 0f0 && (zadj = 0f0)
            z += zadj
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            dia = d + t.diam_growth[i] / bark
            if (_UT_XI1 + cof1) > dia
                psi = cof8 * ((dia - _UT_XI1) / (_UT_XI1 + cof1 - dia))^cof9 *
                      exp(z * ((1f0 - cof7 * cof7)^0.5f0) / cof6)
                hh = (psi / (1f0 + psi)) * cof2 + _UT_XI2
                hh < h && (hh = h)
                htg = hh - h
                sp == 5 && (htg *= 0.95f0)                          # BS spruce (ut/htgf.f:727)
                htg < 0.1f0 && (htg = 0.1f0)
            end
        end
        # finalize (ut/htgf.f:730): HTG·SCALE·XHMULT·exp(HTCON)·MISHGF. XHMULT=1, MISHGF=1, HTCON=htg_cor (0 for utt01).
        t.ht_growth[i] = htg * scale * exp(c.htg_cor[sp])
    end
    return s
end

# ut/htgf.f uneven-aged GENGYM diameter curve HHU (sp17 pinyon form; sp18/19/22 aspen/birch form).
@inline function _ut_hhu(sp::Int, d::Float32, tsite::Float32, ba::Float32)::Float32
    if sp == 17
        return 42.269377f0 * (1f0 - exp(-0.165687f0 * d))^1.184734f0 + 4.5f0
    else                                                       # 18,19,22
        batem = ba < 10f0 ? 10f0 : ba
        return (-2.04f0 + 1.4534f0 * tsite) *
               (1f0 - exp(-0.058112f0 * d))^(1.894400f0 * batem^(-0.192979f0)) + 4.5f0
    end
end

# ut/htgf.f even-aged HHE (Alexander RM-32 spruce-fir site curve, IMODTY=4).
@inline function _ut_hhe(age::Float32, tsite::Float32, bautba::Float32)::Float32
    agetem = age < 30f0 ? 30f0 : age
    hhe = (2.75780f0 * tsite^0.83312f0) *
          (1f0 - exp(-0.015701f0 * agetem))^(22.71944f0 * tsite^(-0.63557f0)) + 4.5f0
    age < agetem && (hhe = ((hhe - 4.5f0) / agetem) * age + 4.5f0)
    ratio = 1f0 - bautba; ratio < 0.728f0 && (ratio = 0.728f0)
    return hhe * ratio
end

"""
    _ut_htg_crsurr(sp, ssite, d, hnow, dgi, bark, bautba, ap0, agerng, ba, pct, htcon, scale, cap, zzran) -> HTG

Per-tree ut/htgf.f CASE(17:19,22) CR-surrogate height increment (GB/NC/FC/BE). Even-aged/uneven-aged mesh.
`ap0` = raw ABIRTH (dubbed via FNDAG at setup, aged FINT/cycle). Deterministic unless zzran passed.
"""
function _ut_htg_crsurr(sp::Int, ssite::Float32, d::Float32, hnow::Float32, dgi::Float32, bark::Float32,
                        bautba::Float32, ap0::Float32, agerng::Float32, ba::Float32, pct::Float32,
                        htcon::Float32, scale::Float32, cap::Float32, zzran::Float32)::Float32
    # ADJUST (ut/htgf.f:309-315) uses SITEAR(ISPC) directly (unfloored).
    adjust = sp == 17 ? 1.0f0 :
             sp == 22 ? 1.3f0 - 0.02f0 * (ssite - 15f0) :
                        0.78f0 + 0.0023f0 * ssite             # 18,19
    # AP breast-high-age adjust (ut/htgf.f:333-336): TSITE floored at 20 for the adjust ONLY.
    tsite = ssite < 20f0 ? 20f0 : ssite
    ap = ap0 - (4.5f0 / (-0.22f0 + 0.0155f0 * tsite)); ap < 1f0 && (ap = 1f0)
    tsite = ssite                                             # reset, unfloored (ut/htgf.f:342), used by HHE/HHU
    hhe1 = _ut_hhe(ap, tsite, bautba)
    hhu1 = _ut_hhu(sp, d, tsite, ba)
    agefut = ap + 10f0
    dfut = d + dgi / bark
    hhe2 = _ut_hhe(agefut, tsite, bautba)
    hhu2 = _ut_hhu(sp, dfut, tsite, ba)
    htg = (hhe2 - hhe1) * adjust
    if agerng > 40f0 && ba >= 70f0
        if pct <= 10f0
            htg = hhu2 - hhu1
        elseif pct > 10f0 && pct < 40f0
            hge = (hhe2 - hhe1) * adjust
            hgu = hhu2 - hhu1
            xwt = ((pct - 10f0) * (10f0 / 3f0)) / 100f0
            htg = xwt * hge + (1f0 - xwt) * hgu
        end
    end
    htg += zzran * 0.1f0                                       # ut/htgf.f:423 (ISTAGF=0 ⇒ no DSTAG factor)
    htg < 0.1f0 && (htg = 0.1f0)
    # finalize (ut/htgf.f:730 + 740): HTG·SCALE·XHMULT·exp(HTCON)·MISHGF, then SIZCAP(sp,4). XHMULT=MISHGF=1.
    htg = htg * scale * exp(htcon)
    (hnow + htg > cap) && (htg = max(cap - hnow, 0.1f0))
    return htg
end

# ut/findag.f CASE(20,21) — SO-variant effective-age inversion. Returns (sitage, sitht, htmax, agmax).
function _ut_findag_so(sp::Int, h::Float32, sindx::Float32)
    agmax1 = sp == 20 ? 50f0 : 100f0
    htmax1 = sp == 20 ? 20f0 : 100f0
    if h >= htmax1
        return (agmax1 + (h - htmax1) / 0.10f0, h, htmax1, agmax1)
    end
    ag = 2f0; incrng = 0; hguess = 0f0
    while true
        oldhg = hguess
        s45 = sindx - 4.5f0
        hguess = s45 / (0.6192f0 - 5.3394f0 / s45 + 240.29f0 * ag^(-1.4f0) +
                        (3368.9f0 / s45) * ag^(-1.4f0)) + 4.5f0
        if hguess >= 1f0
            (abs(hguess - h) <= 2f0 || h < hguess) && return (ag, hguess, htmax1, agmax1)
            d2 = hguess - oldhg
            (oldhg != 0f0 && d2 >= 0.05f0) && (incrng = 1)
            (incrng == 1 && d2 < 0.05f0) && return (ag, hguess, htmax1, agmax1)
        end
        ag += 2f0
        ag > agmax1 && return (agmax1, h, htmax1, agmax1)
    end
end

"""
    _ut_htg_mcbi(sp, ssite, d, hnow, avh, icr, htcon, scale, cap) -> HTG

Per-tree ut/htgf.f CASE(20,21) MC/BI (SO-variant) height increment. FINDAG(H)→SITAGE/SITHT/HTMAX/AGMAX,
Chapman-Richards POTHTG, crown + relative-height HTGMOD. `icr` = raw ICR (crown %, 0-100).
"""
function _ut_htg_mcbi(sp::Int, ssite::Float32, d::Float32, hnow::Float32, avh::Float32, icr::Float32,
                      htcon::Float32, scale::Float32, cap::Float32)::Float32
    sc = scale * exp(htcon)                                   # SCALE·XHMULT·exp(HTCON), XHMULT=1
    sitage, sitht, htmax, agmax = _ut_findag_so(sp, hnow, ssite)
    local htg::Float32
    if hnow >= htmax
        htg = 0.1f0 * sc                                      # ut/htgf.f:463 (inner scale) — then 201 scales AGAIN
    else
        local pothtg::Float32
        if sitage > agmax
            pothtg = 0.10f0
        else
            agp10 = sitage + 10f0
            s45 = ssite - 4.5f0
            hguess = s45 / (0.6192f0 - 5.3394f0 / s45 + 240.29f0 * agp10^(-1.4f0) +
                            (3368.9f0 / s45) * agp10^(-1.4f0)) + 4.5f0
            pothtg = hguess - sitht
        end
        # modifiers (ut/htgf.f:489-545): crown Hoerl (CRA=100,CRB=3,CRC=-5) + relative-height Chapman-Richards.
        relht = avh > 0f0 ? hnow / avh : 0f0; relht > 1.5f0 && (relht = 1.5f0)
        crf = icr / 100f0
        hgmdcr = 100f0 * crf^3 * exp(-5f0 * crf); hgmdcr > 1f0 && (hgmdcr = 1f0)
        rhb, rhr, rhm, rhyxs = sp == 20 ? (-1.45f0, 15f0, 1.10f0, 0.10f0) : (-1.10f0, 20f0, 1.10f0, 0.20f0)
        fctrkx = (1f0 / rhyxs)^(rhm - 1f0) - 1f0               # RHK=1
        fctrrb = -1f0 * (rhr / (1f0 - rhb))
        fctrxb = relht^(1f0 - rhb) - 0f0^(1f0 - rhb)          # RHXS=0
        fctrm  = -1f0 / (rhm - 1f0)
        hgmdrh = 1f0 * (1f0 + fctrkx * exp(fctrrb * fctrxb))^fctrm
        htgmod = 0.25f0 * hgmdcr + 0.75f0 * hgmdrh
        htgmod >= 2f0 && (htgmod = 2f0)
        htgmod <= 0f0 && (htgmod = 0.1f0)
        htg = pothtg * htgmod
        (hnow + htg > htmax) && (htg = htmax - hnow)
        htg < 0.1f0 && (htg = 0.1f0)
    end
    # 201/730 finalize (applied to ALL cases — the H≥HTMAX branch is thus double-scaled, faithful to FVS).
    htg = htg * sc
    (hnow + htg > cap) && (htg = max(cap - hnow, 0.1f0))
    return htg
end

"""
    _ut_dub_ages!(s)

ut/cratet.f:662-696 — dub ABIRTH from the current height for inventory trees with no measured age
(birth_age ≤ 0 or > 999), for the aged UT species: aspen/oak closed form (6,13,24); FNDAG spruce-fir
(17:19,22 = cr_fndag IMODTY 4); SO inversion (20,21). Sets AGERNG's input; only CR-surrogate htgf reads
ABIRTH. Runs once at setup (before calibration); gradd.f:205 then ages birth_age by FINT each cycle.
"""
function _ut_dub_ages!(s::StandState)
    p, t = s.plot, s.trees
    any(j -> (sp = Int(t.species[j]);
              sp == 6 || sp == 13 || sp == 24 || _ut_ht_crsurr(sp) || sp == 20 || sp == 21), 1:t.n) || return s
    misscr = false
    @inbounds for i in 1:(t.n + t.ndead)
        t.crown_pct[i] <= 0 && (misscr = true; break)
    end
    bau = misscr ? _cr_badist_bau(t) : nothing
    ba = p.basal_area <= 0f0 ? 25f0 : p.basal_area
    @inbounds for i in 1:t.n
        sp = Int(t.species[i]); h = t.height[i]; h <= 0f0 && continue
        ab = t.birth_age[i]
        (ab > 0f0 && ab <= 999f0) && continue                 # already aged (cratet.f:668)
        ssite = p.sp_site_index[sp]
        sitage = 0f0
        if sp == 6 || sp == 13 || sp == 24
            sitage = (h * 2.54f0 * 12f0 / 26.9825f0)^(1f0 / 1.1752f0)   # aspen/oak (findag.f:97)
        elseif _ut_ht_crsurr(sp)
            bautba = 0f0
            if bau !== nothing
                icls = trunc(Int, t.dbh[i] + 1f0); icls > 41 && (icls = 41); icls < 1 && (icls = 1)
                bautba = ba > 0f0 ? bau[icls] / ba : 0f0; bautba < 0f0 && (bautba = 0f0)
            end
            sitage = cr_fndag(4, ssite, h, bautba, 0f0, sp)   # FNDAG spruce-fir (findag.f:111)
            sitage <= 0f0 && (sitage = 1f0)                   # cratet.f:115
        elseif sp == 20 || sp == 21
            sitage, _, _, _ = _ut_findag_so(sp, h, ssite)
        end
        sitage > 0f0 && (t.birth_age[i] = sitage; t.age_known[i] = true)
    end
    return s
end
