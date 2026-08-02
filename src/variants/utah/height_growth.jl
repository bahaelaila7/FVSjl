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
# Non-SBB branches DEFERRED (not in utt01): PJ 11:16,24 (→REGENT, HTG 0), CR-surrogate 17:19,22, MC/BI 20:21.
# Young-LP accelerator (ut/htgf.f:672-707) bypassed (ICYC>1 | IAGE≤0); ISTAGF=0 (grinit).
# =============================================================================

function height_growth!(s::StandState, ::Utah; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    elev = p.elevation
    zon = (elev >= 80f0 && elev <= 105f0)
    isisp = Int(p.site_species); (isisp < 1 || isisp > 24) && (isisp = 7)
    lsimap = clamp(trunc(Int, p.sp_site_index[isisp] / 10f0 - 0.5f0), 2, 6) - 1
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (d <= 0f0 || h <= 0f0) && continue
        if sp in (11,12,13,14,15,16,24)
            continue                                     # PJ/woodland: height from REGENT (HTG 0)
        elseif (17 <= sp <= 19) || sp == 22 || sp == 20 || sp == 21
            error("UT height_growth! sp $sp (CR-surrogate/MC-BI) not yet ported — not in utt01")
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
