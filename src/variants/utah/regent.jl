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
const UT_RG_DIAM  = Float32[0.4,0.4,0.3,0.3,0.3,0.2,0.4,0.3,0.3,0.5,0.4,0.3,0.3,0.4,0.3,0.3,0.4,0.3,0.3,0.2,0.2,0.3,0.2,0.3]
const UT_RG_AB    = Float32[1.11436, -0.011493, 0.43012f-4, -0.72221f-7, 0.5607f-10, -0.1641f-13]
const _UT_RG_REGYR = 10.0f0

@inline _ut_rg_conifer(sp::Int) = sp <= 5 || (7 <= sp <= 10) || sp == 23
@inline _ut_rg_pjlin(sp::Int)   = (11 <= sp <= 16) || sp == 24 || sp == 17 || (18 <= sp <= 22)  # non-aspen non-conifer

# Conifer/MC ht_dbh inverse (ut/regent.f:432-478, uncalibrated IABFLG=1 fixed curve). Returns DBH given HT.
@inline function _ut_htdbh_dbh(sp::Int, ht::Float32)::Float32
    p2, p3, p4 = sp == 20 ? (1709.7229f0, 5.8887f0, -0.2286f0) : (76.5170f0, 2.2107f0, -0.6365f0)
    hat3 = 4.5f0 + p2 * exp(-p3 * 3.0f0^p4)
    if ht >= hat3
        return exp(log((log(ht - 4.5f0) - log(p2)) / (-p3)) * (1.0f0 / p4))
    else
        return ((ht - 4.51f0) * 2.7f0 / (hat3 - 4.51f0)) + 0.3f0
    end
end

# em/regent.f stash push — tripled sub-records use the small-tree DG/HTG, not stale large-tree DG.
@inline function _ut_rg_stash!(stash, t, i::Int)
    if stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
        stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
        !isempty(stash.is_small) && (stash.is_small[i] = true)
    end
end

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
    @inbounds for i in 1:n
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
        elseif sp == 6                                # aspen — FINDAG age (birth_age, EM-validated form)
            age = Float32(t.birth_age[i]); age < 1.0f0 && (age = 1.0f0)
            hite1 = 26.9825f0 * age^1.1752f0; hite2 = 26.9825f0 * (age + 10.0f0)^1.1752f0
            htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 0.75f0   # aspen: HTGR=POTHTG (no PCTRED/VIGOR)
        else                                          # PJ/GB/MC (non-aspen non-conifer)
            pothtg = (sj / 5.0f0) * (sj * 1.5f0 - h) / (sj * 1.5f0) * 0.83f0
            xcr = Float32(t.crown_pct[i]) / 100.0f0
            vigor = 150.0f0 * xcr^3 * exp(-6.0f0 * xcr) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
            (11 <= sp <= 17 || sp == 24) && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)
            htgr = pothtg * pctred * vigor * con
        end
        # ZZRAN reject-loop (ut/regent.f:340-342): redraw until ∈[−2,0.5]
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = (htgr + zzran * 0.1f0) * scale         # XRHGRO=1
        htgr < 0.1f0 && (htgr = 0.1f0)
        # XWT blend with large-tree HTG
        xmn = UT_RG_XMIN[sp]; xmx = UT_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]; htg < 0.1f0 && (htg = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # ---- small-tree DG (ut/regent.f:380-560). DK/DKK are a 10-YR diameter increment; NO XWT blend
        # (the HTG blend enters via HK=H+HTG). Clamp to DGMX, DDS→DG rescale by SCALE2=YR/FINT, DIAM floor.
        hk = h + htg
        bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
        if hk <= 4.5f0
            t.diam_growth[i] = 0.0f0                    # DBH+=0.001·HK is a birth detail; DG=0
        else
            local dk::Float32, dkk::Float32
            if sp == 10                                # PP Wykoff
                dk = (hk - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0
                dkk = h < 4.5f0 ? d : (h - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0
            elseif (11 <= sp <= 17) || sp == 24        # PJ/GB linear-site
                dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                dkk = h < 4.5f0 ? d : (h - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
            elseif sp == 20 || sp == 21                # MC/BI linear
                dk = 3.1020f0 + 0.0210f0 * hk
                dkk = 3.1020f0 + 0.0210f0 * h; dkk < 0.0f0 && (dkk = d)
                dk < dkk && (dk = dkk + 0.01f0)
            elseif c.ht_dbh_iabflg[sp] == 0            # conifers CALIBRATED (ut/regent.f:487): AX=AA, BX=HT2
                ax = c.ht_dbh_aa[sp]; bx = sd[:ht2][sp]
                dk = bx / (log(hk - 4.5f0) - ax) - 1.0f0; dk < 0.1f0 && (dk = 0.1f0)
                dkk = h <= 4.5f0 ? d : bx / (log(h - 4.5f0) - ax) - 1.0f0
            else                                       # conifers UNCALIBRATED — fixed P2/P3/P4 htdbh curve
                dk = _ut_htdbh_dbh(sp, hk)
                dkk = h <= 4.5f0 ? d : _ut_htdbh_dbh(sp, h)
            end
            dgk = (dk - dkk) * bark                     # XRDGRO=1
            dgk < 0.0f0 && (dgk = 0.0f0)
            dgmx = UT_RG_DGMAX[sp] * scale
            dgk > dgmx && (dgk = dgmx)
            scale2 = _UT_RG_REGYR / fint                # YR/FINT (=1 for 10-yr ⇒ transform is identity)
            dds = dgk * (2.0f0 * bark * d + dgk) * scale2
            dgk = sqrt((d * bark)^2 + dds) - bark * d
            (d + dgk) < UT_RG_DIAM[sp] && (dgk = UT_RG_DIAM[sp] - d)
            t.diam_growth[i] = dgk
        end
        _ut_rg_stash!(stash, t, i)
    end
    return s
end
