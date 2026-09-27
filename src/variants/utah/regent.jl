# =============================================================================
# regent.jl (utah) — UT small-tree height + diameter growth (ut/regent.f). Chunk 6.
#
# small_tree_growth!(s, stash, ::Utah) overrides HTG/DG for trees below XMAX, blended with the
# large-tree prediction by XWT=(D−XMN)/(XMX−XMN). Height: HTGR = POTHTG·PCTRED·VIGOR·CON (aspen: HTGR=
# POTHTG), + ZZRAN·0.1 (reject-loop until ∈[−2,0.5]), ·XRHGRO(=1)·SCALE(=FINT/REGYR=1). POTHTG per
# SELECT CASE: conifers(1:5,7:10,23)=SJ/5; PJ/GB/MC(11:22,24 non-aspen)=((SJ/5)(SJ·1.5−H)/(SJ·1.5))·0.83;
# aspen(6)=(26.9825·((AGE+10)^1.1752−AGE^1.1752))/(2.54·12)·RSIMOD·CON·0.75 (AGE from FINDAG; here birth_age,
# EM-validated). CON=RHCON(=1)·exp(HCOR). Small-tree DG (D<BKPT): conifer ht_dbh curve (P2/P3/P4 uncal, AX/BX
# calib), PJ linear, MC/BI, PP Wykoff → DK/DKK → DDS → XWT-blended DG. Coefficients ut/regent.f DATA + blkdat.
# utt01 small trees = WF(conifer) + AS(aspen). PJ/GB/MC paths ported faithfully (need pure stands to validate).
# =============================================================================

const UT_RG_DGMAX = Float32[2.8,2.8,2.4,3.6,3.6,2.5,3.5,3.6,3.6,2.8,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.5,2.5,2.0,2.0,2.5,2.8,2.0]
const UT_RG_XMAX  = Float32[4,4,4,4,4,4,5,4,6,6,99,99,99,99,99,99,199,2,2,99,99,2,6,99]
const UT_RG_XMIN  = Float32[2,2,2,2,2,2,1,2,2,2,90,90,90,90,90,90,99,0.5,0.5,90,90,0.5,2,90]
# BKPT = BREAK (ut/blkdat.f DATA BREAK): the DBH at/above which the small-tree HT-DBH diameter increment is
# NOT applied (ut/regent.f:383 IF(D.GE.BKPT) GO TO 23) — the tree keeps its large-tree DGDRIV DG but still
# takes the XWT-blended regent HEIGHT growth. Conifers 3", NC/FC/BE 1", woodland 99 (never large-tree).
const UT_RG_BREAK = Float32[3,3,3,3,3,3,3,3,3,3,99,99,99,99,99,99,99,1,1,99,99,1,3,99]
const UT_RG_DIAM  = Float32[0.4,0.4,0.3,0.3,0.3,0.2,0.4,0.3,0.3,0.5,0.4,0.3,0.3,0.4,0.3,0.3,0.4,0.3,0.3,0.2,0.2,0.3,0.2,0.3]
const UT_RG_AB    = Float32[1.11436, -0.011493, 0.43012f-4, -0.72221f-7, 0.5607f-10, -0.1641f-13]
const _UT_RG_REGYR = 10.0f0

@inline _ut_rg_conifer(sp::Int) = sp <= 5 || (7 <= sp <= 10) || sp == 23
@inline _ut_rg_pjlin(sp::Int)   = (11 <= sp <= 16) || sp == 24 || sp == 17 || (18 <= sp <= 22)  # non-aspen non-conifer

# (The Curtis-Arney P2/P3/P4 ht-dbh inverse formerly here was ut/regent.f's UTBR1 branch, reached only for
#  LHTDRG=.FALSE. species — i.e. MC/BI, which use their own linear DK — so it is NEVER hit for UT conifers.
#  It was wrongly used for the uncalibrated conifer DG; replaced by the Wykoff BX/AX(HT1) curve below.)


function small_tree_growth!(s::StandState, stash, ::Utah; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species; slo = sd[:site_lo]; shi = sd[:site_hi]
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    dgsd = s.control.dg_sd
    scale = fint / _UT_RG_REGYR                       # SCALE = FINT/REGYR (=1 for 10-yr)
    # PCTRED (density modifier) — stand-level, ut/regent.f:171-176
    xd = avh * (relden / 100.0f0)
    ab = UT_RG_AB
    pctred = ab[1] + xd*(ab[2] + xd*(ab[3] + xd*(ab[4] + xd*(ab[5] + xd*ab[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # WK4 = CLGMULT's per-tree climate multiplier (dgdriv.f fills it before REGENT; 1 without CLIMATE).
    wk4 = clim_wk4(s, Float32(current_cycle_year(s)) + fint / 2f0)
    cur_year = current_cycle_year(s)
    # ut/regent.f:183-223 is SPECIES-MAJOR (DO 30 ISPC … I=IND1(I3)); the per-tree ZZRAN draw must follow it.
    @inbounds for i in species_major_order(s)
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= UT_RG_XMAX[sp] || t.tpa[i] <= 0.0f0) && continue
        h = t.height[i]
        sitear = p.sp_site_index[sp]
        si = sitear; si > shi[sp] && (si = shi[sp]); si <= slo[sp] && (si = slo[sp] + 0.5f0)
        relsi = (si - slo[sp]) / (shi[sp] - slo[sp]); rsimod = 0.5f0 * (1.0f0 + relsi)
        sj = sitear
        con = exp(c.htg_cor_small[sp])                # RHCON(=1)·exp(HCOR)
        # POTHTG + HTGR
        if _ut_rg_conifer(sp)
            pothtg = sj / 5.0f0
            xcr = Float32(t.crown_pct[i]) / 100.0f0
            vigor = 150.0f0 * xcr^3 * exp(-6.0f0 * xcr) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
            htgr = pothtg * pctred * vigor * con
        elseif sp == 6                                # aspen — FINDAG site age from HEIGHT (ut/findag.f:97)
            age = (h * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0)   # SITAGE=(H·30.48/26.9825)^(1/1.1752)
            hite1 = 26.9825f0 * age^1.1752f0; hite2 = 26.9825f0 * (age + 10.0f0)^1.1752f0
            htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 0.75f0   # aspen: HTGR=POTHTG (no PCTRED/VIGOR)
        else                                          # PJ/GB/MC (non-aspen non-conifer)
            pothtg = (sj / 5.0f0) * (sj * 1.5f0 - h) / (sj * 1.5f0) * 0.83f0
            xcr = Float32(t.crown_pct[i]) / 100.0f0
            vigor = 150.0f0 * xcr^3 * exp(-6.0f0 * xcr) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
            (11 <= sp <= 17 || sp == 24) && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)
            htgr = pothtg * pctred * vigor * con
        end
        # ut/regent.f:563-567: with tripling the SAME loop body (label 2) runs again for each tripled copy
        # (L=1,2; K=ITRN+2I-2+L): a fresh ZZRAN, the XWT blend with the COPY's large-tree HTG (htgf.f:736-762
        # TEMHTG), its own SIZCAP check, HK=H+HTG(K) and DBH(K)/DG(K) — all from the central's D=DBH(I), H=HT(I).
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xmn = UT_RG_XMIN[sp]; xmx = UT_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        cap = s.control.sp_size_cap[sp, 4]
        large_htg = t.ht_growth[i]
        bark = ut_bratio(s.coef.species, sp, d)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)   # XRDGRO = REGDMULT
        nrec = stash !== nothing && !isempty(stash.htgU) && i <= length(stash.htgU) ? 3 : 1
        for l in 0:(nrec - 1)
            # ZZRAN reject-loop (ut/regent.f:340-342): redraw until ∈[−2,0.5]
            zzran = 0.0f0
            if dgsd >= 1.0f0
                while true
                    zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                    (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                end
            end
            # ut/regent.f:344-350 SELECT CASE(ISPC): the CR-surrogate hardwoods (17:19,22) take ZZRAN·0.2 and the
            # CLGMULT climate multiplier WK4(I); every other species ZZRAN·0.1 and no WK4. XRHGRO = REGHMULT.
            # No floor on HTGR itself — only on the blended HTG(K) (ut/regent.f:360).
            htgk = if 17 <= sp <= 19 || sp == 22
                (htgr + zzran * 0.2f0) * xrhgro * scale * (wk4 === nothing ? 1f0 : wk4[i])
            else
                (htgr + zzran * 0.1f0) * xrhgro * scale
            end
            lh = l == 0 ? large_htg : (stash.htg_copy[i] ? (l == 1 ? stash.htgU[i] : stash.htgL[i]) : large_htg)
            htg = htgk * (1.0f0 - xwt) + xwt * lh; htg < 0.1f0 && (htg = 0.1f0)
            (h + htg > cap) && (htg = max(cap - h, 0.1f0))
            # ut/regent.f:383 IF(D.GE.BKPT) GO TO 23 — for BKPT≤D<XMAX the record gets the regent (XWT-blended)
            # HEIGHT growth but KEEPS its large-tree DGDRIV diameter growth.
            below = d < UT_RG_BREAK[sp]
            dbhk = d; dgk = 0f0; direct = false
            if below
                hk = h + htg
                if hk <= 4.5f0
                    # ut/regent.f:384-386: DG(K)=0.0, DBH(K)=D+0.001·HK — the 0.001·HK nudge every cycle
                    # (see the CCF note in the git history: omitting it left tiny-seedling woodland stands at a ~5×
                    # inflated CCF). diam_growth stays 0 so the apply-loop adds nothing.
                    dgk = 0.0f0; dbhk = d + 0.001f0 * hk; direct = true
                else
                    dgk = _ut_rg_dk_dg(s, c, sd, sp, d, h, hk, htg, bark, xrdgro, sitear, scale, fint)
                end
                dgk = dg_bound(nothing, nothing, sp, dbhk, dgk, s.control.sp_size_cap)   # CALL DGBND(ISPC,DBH(K),DG(K))
            end
            if l == 0
                t.ht_growth[i] = htg
                if below
                    direct && (t.dbh[i] = dbhk)
                    t.diam_growth[i] = dgk
                end
            elseif l == 1
                stash.htgU[i] = htg; stash.is_small[i] = true
                below && (stash.dgU[i] = dgk; stash.dbhU[i] = dbhk)
            else
                stash.htgL[i] = htg
                below && (stash.dgL[i] = dgk; stash.dbhL[i] = dbhk)
            end
        end
    end
    return s
end

# ut/regent.f:388-556 (non-LESTB): DK/DKK from the HT-DBH relation for species ISPC, DG(K) by subtraction with the
# REGDMULT multiplier and bounds, rescaled to the projection length (SCALE2 = YR/FINT). D, H = the central record's.
function _ut_rg_dk_dg(s, c, sd, sp::Int, d::Float32, h::Float32, hk::Float32, htg::Float32, bark::Float32,
                      xrdgro::Float32, sitear::Float32, scale::Float32, fint::Float32)::Float32
    local dk::Float32, dkk::Float32
    if sp == 10                                # PP Wykoff
        dk = (hk - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0
        dkk = h < 4.5f0 ? d : (h - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0
    elseif (11 <= sp <= 17) || sp == 24        # PJ/GB linear-site
        dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
        dkk = h < 4.5f0 ? d : (h - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
    elseif sp == 20 || sp == 21                # MC/BI (SO/WC origin), ut/regent.f:408-465
        dkk = 3.1020f0 + 0.0210f0 * h; dkk < 0.0f0 && (dkk = d)
        dk = 3.1020f0 + 0.0210f0 * hk; dk < dkk && (dk = dkk + 0.01f0)
        # :430-465 inventory Curtis-Arney dub whenever the Wykoff calibration is off or did not happen
        # (`.NOT.LHTDRG .OR. IABFLG==1`) — always for MC (LHTDRG(20)=.FALSE.).
        if !s.control.ht_drag_sp[sp] || c.ht_dbh_iabflg[sp] == 1
            p2, p3, p4 = sp == 20 ? (1709.7229f0, 5.8887f0, -0.2286f0) : (76.5170f0, 2.2107f0, -0.6365f0)
            hat3 = 4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4)
            ca(hh) = hh >= hat3 ? exp(log((log(hh - 4.5f0) - log(p2)) / (-1f0 * p3)) * (1f0 / p4)) :
                                  ((hh - 4.51f0) * 2.7f0) / (4.5f0 + p2 * exp(-1f0 * p3 * (3f0^p4)) - 4.51f0) + 0.3f0
            dk = ca(hk)
            dkk = h <= 4.5f0 ? d : ca(h)
        end
    else                                       # conifers: WYKOFF HT-DBH DK=BX/(ln(HK-4.5)−AX)−1
        # ut/regent.f:398-403 BX=HT2, AX=HT1 (IABFLG=1) or AA (IABFLG=0); the BX/AX branch (466-472) is taken for
        # all LHTDRG=.TRUE. species (every non-MC/BI conifer).
        ax = c.ht_dbh_iabflg[sp] == 0 ? c.ht_dbh_aa[sp] : sd[:ht1][sp]
        bx = sd[:ht2][sp]
        dk = bx / (log(hk - 4.5f0) - ax) - 1.0f0; dk < 0.1f0 && (dk = 0.1f0)
        dkk = h <= 4.5f0 ? d : bx / (log(h - 4.5f0) - ax) - 1.0f0
    end
    dgmx = UT_RG_DGMAX[sp] * scale
    local dgk::Float32
    if sp == 20 || sp == 21                     # ut/regent.f:516-532 CASE(20,21)
        h < 4.5f0 && (dkk = d)
        if dk < 0.0f0 || dkk < 0.0f0
            dgk = htg * 0.2f0 * bark * xrdgro
        else
            dgk = (dk - dkk) * bark * xrdgro
        end
        (s.control.ht_drag_sp[sp] && c.ht_dbh_iabflg[sp] == 0) && (dgk = 0.1f0 * htg * xrdgro)
        dgk < 0.0f0 && (dgk = 0.1f0)
        dgk > dgmx && (dgk = dgmx)
    else
        if dk < 0.0f0 || dkk < 0.0f0           # ut/regent.f:535-540 CASE DEFAULT
            dgk = htg * 0.2f0 * bark * xrdgro
        else
            dgk = (dk - dkk) * bark * xrdgro
        end
    end
    dgk < 0.0f0 && (dgk = 0.0f0)
    dgk > dgmx && (dgk = dgmx)
    scale2 = _UT_RG_REGYR / fint                # YR/FINT (=1 for 10-yr ⇒ transform is identity)
    dds = dgk * (2.0f0 * bark * d + dgk) * scale2
    dgk = sqrt((d * bark)^2 + dds) - bark * d
    (d + dgk) < UT_RG_DIAM[sp] && (dgk = UT_RG_DIAM[sp] - d)
    return dgk
end

# ut_esgent! — ut/esgent.f: SPESRT, REGENT(.TRUE.,ITRNIN) over this cycle's new records (ITRNIN = nstart+1), then the
# WK4=HTIMLT tail. ut/regent.f LESTB (:147-560):
# - FNT=FINT−5 (LSKIPH when FINT≤5), SCALE=FNT/REGYR, SCALE2=YR/FNT, DGMX=DGMAX·SCALE.
# - PCTRED from the mid-period CCF=(5/FINT)·RELDEN+((FINT−5)/FINT)·ATCCF, AVHT likewise from AVH/ATAVH (:162-171), where
#   RELDEN/AVH are GRADD's post-growth PRE-regen DENSE values (establish! has re-DENSEd with the seedlings, so the caller
#   passes pre-regen snapshots).
# - species-major records; each new record first draws its crown (:230-241, PCCF(ITRE(I)), main stream) and then its
#   ZZRAN — the two draws interleave per record. Aspen POTHTG from the Sheppard curve at SITAGE=ABIRTH (:284-287);
#   HTGR=(HTGR+ZZRAN·0.1)·XRHGRO·SCALE (·0.2·WK4 for 17:19,22, WK4=HTIMLT here); XWT=0; HTG floored at 0.1; SIZCAP.
# - HK≤4.5 ⇒ DBH=D+0.001·HK, DG=0; else DBH(K)=DK (PP linear, PJ/WJ/oak (HK−4.5)·10/(SI−4.5), 20/21 SO equations,
#   else Wykoff with AX=HT1|AA), floored at DIAM, +0.001·HK, DG=DBH (:483-497); DGBND.
# esgent.f:51-65 then scales HTG by WK4 (<1 resets sub-4.5' DBH to 0.1+0.001·HT, else rescales DBH/DG by HT/HTEMP) and caps
# HT at HHTMAX; estab.f:707 adds GENTIM to ABIRTH.
function ut_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atavh::Float32 = -1.0f0, atrelden::Float32 = -1.0f0,
                    relden_pre::Float32 = -1.0f0, avh_pre::Float32 = -1.0f0,
                    pccf_pre::Vector{Float32} = Float32[])
    p, t, c = s.plot, s.trees, s.calib
    nstart >= t.n && return s
    sd = s.coef.species; slo = sd[:site_lo]; shi = sd[:site_hi]
    relden = relden_pre >= 0f0 ? relden_pre : p.relative_density
    avh = avh_pre >= 0f0 ? avh_pre : p.avg_height
    pccfv = isempty(pccf_pre) ? s.density.point_ccf : pccf_pre
    dgsd = s.control.dg_sd
    cur_year = current_cycle_year(s)
    n = t.n
    @inbounds for i in (nstart + 1):n; t.diam_growth[i] = 0f0; t.ht_growth[i] = 0f0; end   # estab.f:641-642
    lskiph = false; fnt = fint
    if fint <= 5f0
        lskiph = true
    else
        fnt = fnt - 5f0
    end
    scale = fnt / _UT_RG_REGYR
    ccf = relden; avht = avh
    if fnt > 0f0
        atccf = atrelden >= 0f0 ? atrelden : relden; atah = atavh >= 0f0 ? atavh : avh
        ccf = (5f0 / fint) * relden + ((fint - 5f0) / fint) * atccf
        avht = (5f0 / fint) * avh + ((fint - 5f0) / fint) * atah
    end
    xd = avht * (ccf / 100f0); xd > 300f0 && (xd = 300f0)
    ab = UT_RG_AB
    pctred = ab[1] + xd * (ab[2] + xd * (ab[3] + xd * (ab[4] + xd * (ab[5] + xd * ab[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    order = sort(collect((nstart + 1):n); by = i -> (Int(t.species[i]), i))   # SPESRT order of the new records
    @inbounds for i in order
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(pccfv)) ? pccfv[pt] : 0f0
        cr0 = 0.89722f0 - 0.0000461f0 * pccf                 # ut/regent.f:233-240
        ran = 0f0
        while true; ran = bachlo(s.rng, 0f0, 1f0); (-1f0 <= ran <= 1f0) && break; end
        cr0 = cr0 + 0.07985f0 * ran
        cr0 > 0.90f0 && (cr0 = 0.90f0); cr0 < 0.20f0 && (cr0 = 0.20f0)
        icr = trunc(Int32, cr0 * 100f0 + 0.5f0)
        t.crown_pct[i] = icr; t.crown_ratio[i] = Float32(icr)
        d >= UT_RG_XMAX[sp] && continue
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        con = fexp(c.htg_cor_small[sp])                      # RHCON=1
        sitear = p.sp_site_index[sp]
        si = sitear; si > shi[sp] && (si = shi[sp]); si <= slo[sp] && (si = slo[sp] + 0.5f0)
        rsimod = 0.5f0 * (1f0 + (si - slo[sp]) / (shi[sp] - slo[sp]))
        sj = sitear
        local htg::Float32
        if lskiph
            htg = 0f0
        else
            local htgr::Float32
            if sp == 6                                       # regent.f:273-297 aspen, LESTB ⇒ SITAGE=ABIRTH(I)
                sitage = t.birth_age[i]
                hite1 = 26.9825f0 * fpow(sitage, 1.1752f0)
                hite2 = 26.9825f0 * fpow(sitage + 10f0, 1.1752f0)
                pothtg = (hite2 - hite1) / (2.54f0 * 12f0) * rsimod * con
                htgr = pothtg * 0.75f0
            else
                pothtg = _ut_rg_conifer(sp) ? sj / 5f0 : ((sj / 5f0) * (sj * 1.5f0 - h) / (sj * 1.5f0)) * 0.83f0
                x = Float32(t.crown_pct[i]) / 100f0
                vigor = (150f0 * fpow(x, 3f0) * fexp(-6f0 * x)) + 0.3f0
                vigor > 1f0 && (vigor = 1f0)
                ((11 <= sp <= 17) || sp == 24) && (vigor = 1f0 - ((1f0 - vigor) / 3f0))
                htgr = pothtg * pctred * vigor * con
            end
            zzran = 0f0
            if dgsd >= 1f0
                while true
                    zzran = bachlo(s.rng, 0f0, 1f0)
                    (zzran <= 0.5f0 && zzran >= -2f0) && break
                end
            end
            htgr = (17 <= sp <= 19 || sp == 22) ? (htgr + zzran * 0.2f0) * xrhgro * scale * t.htimlt[i] :
                                                  (htgr + zzran * 0.1f0) * xrhgro * scale
            htg = htgr                                       # XWT=0 under LESTB (:359)
            htg < 0.1f0 && (htg = 0.1f0)
            cap = s.control.sp_size_cap[sp, 4]
            if h + htg > cap; htg = cap - h; htg < 0.1f0 && (htg = 0.1f0); end
        end
        t.ht_growth[i] = htg
        d >= UT_RG_BREAK[sp] && continue                     # regent.f:381 GO TO 23
        hk = h + htg
        diam = UT_RG_DIAM[sp]
        local dbh::Float32, dg::Float32
        if hk <= 4.5f0
            dg = 0f0; dbh = d + 0.001f0 * hk
        else
            local dk::Float32
            dat45 = 0f0
            if sp == 10
                dk = (hk - 8.31485f0 + 0.59200f0 * 7f0) / 3.03659f0
            elseif (11 <= sp <= 17) || sp == 24
                dk = (hk - 4.5f0) * 10f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
            elseif sp == 20 || sp == 21                      # SO (from WC) equations (:402-468)
                dat45 = 3.1020f0 + 0.0210f0 * 4.5f0
                dkk = 3.1020f0 + 0.0210f0 * h; dkk < 0f0 && (dkk = d)
                dk = 3.1020f0 + 0.0210f0 * hk; dk < dkk && (dk = dkk + 0.01f0)
                if !s.control.ht_drag_sp[sp] || c.ht_dbh_iabflg[sp] == 1
                    p2, p3, p4 = sp == 20 ? (1709.7229f0, 5.8887f0, -0.2286f0) : (76.5170f0, 2.2107f0, -0.6365f0)
                    hat3 = 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3.0f0, p4))
                    dk = hk >= hat3 ? fexp(flog((flog(hk - 4.5f0) - flog(p2)) / (-1f0 * p3)) * 1f0 / p4) :
                                      (((hk - 4.51f0) * 2.7f0) / (4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3f0, p4)) - 4.51f0)) + 0.3f0
                end
            else
                bx = sd[:ht2][sp]; ax = c.ht_dbh_iabflg[sp] == 1 ? sd[:ht1][sp] : c.ht_dbh_aa[sp]
                dk = (bx / (flog(hk - 4.5f0) - ax)) - 1f0; dk < 0.1f0 && (dk = 0.1f0)
            end
            if (sp == 20 || sp == 21) && dat45 > 0f0 && hk >= 4.5f0 && s.control.ht_drag_sp[sp] &&
               c.ht_dbh_iabflg[sp] == 0
                dbh = dk - dat45 + diam
            else
                dbh = dk
            end
            dbh < diam && (dbh = diam)
            dbh = dbh + 0.001f0 * hk
            dg = dbh
            (dbh + dg) < diam && (dg = diam - dbh)
        end
        dg = dg_bound(nothing, nothing, sp, dbh, dg, s.control.sp_size_cap)   # DGBND
        t.dbh[i] = dbh; t.diam_growth[i] = dg
    end
    @inbounds for i in (nstart + 1):n                       # esgent.f:51-65
        sp = Int(t.species[i])
        htemp = t.height[i] + t.ht_growth[i]
        t.ht_growth[i] = t.ht_growth[i] * t.htimlt[i]
        t.height[i] = t.height[i] + t.ht_growth[i]
        if t.htimlt[i] < 1f0
            if t.height[i] < 4.5f0
                t.dbh[i] = 0.1f0 + 0.001f0 * t.height[i]; t.diam_growth[i] = 0f0
            else
                t.dbh[i] = t.dbh[i] * (t.height[i] / htemp)
                t.diam_growth[i] = t.dbh[i] * (t.height[i] / htemp)
            end
        end
        t.height[i] > _UT_ES_HHTMAX[sp] && (t.height[i] = _UT_ES_HHTMAX[sp])
    end
    esgent_add_gentim!(s, nstart, fint)                     # estab.f:707 ABIRTH += GENTIM
    return s
end
