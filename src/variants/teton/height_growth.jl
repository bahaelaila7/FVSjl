# =============================================================================
# height_growth.jl (teton) — TT large-tree height growth (tt/htgf.f). Chunk 4.
#
# ttt01 species (WB/AS/LP/ES/AF) use the DEFAULT case = the Schreuder-Hafley SBB
# bivariate height-DBH model (tt/htgf.f:590). Height comes AFTER DG (uses DIA=DBH+DG/BARK).
#   IICR=INT(ICR/10+0.5) cap9; KEYCR = IICR≤2?1 : IICR≤7?2 : 3; JSPC = sp≤10?sp : sp14?6 : 11;
#   K=(JSPC−1)·3+KEYCR.  COF(L,K)=TT_HTCOF[K,L].  XI1=0.1, XI2=4.5.
#   OOB (HT≤4.5 | DBH≥XI1+COF1 | HT≥XI2+COF2 | DBH≤0.1) → HTG=0.1.
#   Y1=(DBH−XI1)/COF1; Y2=(HT−XI2)/COF2; FBY=logit; Z=(COF4+COF6·FBY2−COF7·(COF3+COF5·FBY1))·(1−COF7²)^−0.5;
#   ZBIAS=AZBIAS+BZBIAS·ELEV (0 if ELEV<55|>80; sp6/14 use ELEV−20 + a ZADJ quadratic); Z−=ZBIAS.
#   DIA=DBH+DG/BARK; if XI1+COF1≤DIA → HTG=0.1; else PSI=COF8·((DIA−XI1)/(XI1+COF1−DIA))^COF9·exp(Z·√(1−COF7²)/COF6);
#   H=(PSI/(1+PSI))·COF2+XI2; HTG=max(H−HT,0.1).  Finalize: HTG·scale·XHMULT·exp(HTCON) (XHMULT=1; HTCON calib).
# DBH<1.5 (DEFAULT) → 0 here (REGENT small-tree, chunk 6). Non-DEFAULT species error (not in ttt01).
# =============================================================================

const _TT_XI1 = 0.1f0
const _TT_XI2 = 4.5f0
@inline _tt_default_ht(sp::Int) = sp <= 3 || sp == 5 || sp == 6 || (7 <= sp <= 9) || sp == 14 || sp == 17

function height_growth!(s::StandState, ::Teton; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    elev = p.elevation
    zon = (elev >= 55f0 && elev <= 80f0)             # ZBIAS only active for ELEV∈[55,80]
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (d <= 0f0 || h <= 0f0) && continue
        if !_tt_default_ht(sp)
            error("TT height_growth! sp $sp (PP/PM/UJ/RM/NC/OH/BI/MC) not yet ported — not in ttt01")
        end
        d < 1.5f0 && continue                        # small tree → REGENT (chunk 6); HTG stays 0
        # SBB coefficient row
        iicr = trunc(Int, Float32(t.crown_pct[i]) / 10f0 + 0.5f0); iicr > 9 && (iicr = 9)
        keycr = iicr <= 2 ? 1 : (iicr <= 7 ? 2 : 3)
        jspc = sp <= 10 ? sp : (sp == 14 ? 6 : 11)
        k = (jspc - 1) * 3 + keycr
        cof1 = TT_HTCOF[k, 1]; cof2 = TT_HTCOF[k, 2]; cof3 = TT_HTCOF[k, 3]
        cof4 = TT_HTCOF[k, 4]; cof5 = TT_HTCOF[k, 5]; cof6 = TT_HTCOF[k, 6]
        cof7 = TT_HTCOF[k, 7]; cof8 = TT_HTCOF[k, 8]; cof9 = TT_HTCOF[k, 9]
        htg = 0.1f0
        # OOB check (label 180) — else the SBB
        if h > 4.5f0 && (_TT_XI1 + cof1) > d && (_TT_XI2 + cof2) > h && d > 0.1f0
            y1 = (d - _TT_XI1) / cof1
            y2 = (h - _TT_XI2) / cof2
            fby1 = log(y1 / (1f0 - y1))
            fby2 = log(y2 / (1f0 - y2))
            z = (cof4 + cof6 * fby2 - cof7 * (cof3 + cof5 * fby1)) * (1f0 - cof7 * cof7)^(-0.5f0)
            zbias = zon ? (sp == 6 || sp == 14 ? TT_AZBIAS[sp] + TT_BZBIAS[sp] * (elev - 20f0) :
                                                 TT_AZBIAS[sp] + TT_BZBIAS[sp] * elev) : 0f0
            (z - zbias >= 2f0 && zbias < 0f0) && (zbias = 0f0)
            z -= zbias
            if sp == 6 || sp == 14
                zadj = 0.1f0 - 0.10273f0 * z + 0.00273f0 * z * z
                zadj < 0f0 && (zadj = 0f0)
                z += zadj
            end
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            dia = d + t.diam_growth[i] / bark
            if (_TT_XI1 + cof1) > dia
                psi = cof8 * ((dia - _TT_XI1) / (_TT_XI1 + cof1 - dia))^cof9 *
                      exp(z * ((1f0 - cof7 * cof7)^0.5f0) / cof6)
                hh = (psi / (1f0 + psi)) * cof2 + _TT_XI2
                hh < h && (hh = h)
                htg = hh - h
                htg < 0.1f0 && (htg = 0.1f0)
            end
        end
        # finalize (label 201): HTG·SCALE·XHMULT·exp(HTCON). XHMULT=1 (no MULTS); HTCON=0 (uncalibrated ttt01).
        t.ht_growth[i] = htg * scale
    end
    return s
end
