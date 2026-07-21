# =============================================================================
# crown.jl (centralrockies) — CR crown ratio = GENGYM (cr/crown.f + gemcr.f)
#
# gemcr.f predicts crown LENGTH CL (linear in HF/DF/BAT/BAU) and converts to a
# ratio CR=CL/HF, except MH/RC/WL (CASE 6/7/8) which predict CR directly (ICRFLG=1).
# crown.f is the DUB/change wrapper (crown_ratio! hook).
#
# Math: ALOG->flog; RELDEN**2 (int) -> RELDEN*RELDEN.
# =============================================================================

"""
    cr_gemcr(imodty, is, bau, bat, hf, df, h, relden, pcti) -> CR

GENGYM crown ratio (cr/gemcr.f). Returns crown ratio (0-1). `HF`=future height (H+HTG),
`DF`=future dbh, `BAT`=stand BA, `BAU`=BA-above-class, `PCTI`=BA percentile.
"""
function cr_gemcr(imodty::Int, is::Int, bau::Float32, bat::Float32, hf::Float32, df::Float32,
                  h::Float32, relden::Float32, pcti::Float32)::Float32
    icrflg = 0
    cl = 0.0f0; cr = 0.0f0
    if is == 1 || is == 2
        cl = imodty <= 2 ? 0.50706f0 + 0.73070f0 * hf : 0.36135f0 + 0.57085f0 * hf
    elseif is == 3
        cl = 6.47479f0 + 0.54482f0 * df + 0.50703f0 * hf - 0.03326f0 * bat
    elseif is == 4 || is == 5
        cl = 6.22959f0 + 0.67587f0 * hf - 0.03098f0 * bat
    elseif is == 6
        cr = 0.3450f0 - 0.00264f0 * bat + 0.00000512f0 * relden * relden -
             0.25138f0 * flog(h) + 0.05140f0 * flog(pcti)
        icrflg = 1
    elseif is == 7
        cr = -1.6053f0 + 0.17479f0 * flog(bat) - 0.00183f0 * relden -
             0.00560f0 * df + 0.11050f0 * flog(pcti)
        icrflg = 1
    elseif is == 8
        cr = 0.03441f0 - 0.00204f0 * bat + 0.30066f0 * flog(df) - 0.59302f0 * flog(h)
        icrflg = 1
    elseif is == 11
        cl = 5.00215f0 + 0.06334f0 * hf + 0.88236f0 * df - 0.03821f0 * bau
    elseif is == 9 || is == 10 || is == 12 || is == 14 || is == 16 ||
           (23 <= is <= 27) || (29 <= is <= 35) || is == 37
        cl = -0.59373f0 + 0.67703f0 * hf
    elseif is == 13 || is == 36
        if imodty == 1
            cl = 5.63367f0 + 0.56252f0 * hf - 0.06411f0 * bat
        elseif imodty == 2
            cl = 4.35671f0 + 0.84714f0 * df + 0.32549f0 * hf - 0.03802f0 * bat
        else
            cl = 3.49178f0 + 0.80767f0 * df + 0.17421f0 * hf - 0.03272f0 * bat
        end
    elseif is == 15
        cl = 3.03832f0 + 0.65587f0 * hf - 0.01792f0 * bat
    elseif is == 17
        cl = 3.61635f0 + 0.93639f0 * df + 0.61547f0 * hf - 0.02360f0 * bat
    elseif is == 18
        cl = imodty <= 2 ? 1.05857f0 + 0.68442f0 * hf : 3.22244f0 + 0.44315f0 * hf + 0.44755f0 * df
    elseif is == 19
        cl = 0.15768f0 + 0.74697f0 * hf
    elseif is == 20 || is == 21 || is == 22 || is == 28 || is == 38
        cl = 5.17281f0 + 0.32552f0 * hf - 0.01675f0 * bat
    else
        cl = 1.0f0
    end
    if icrflg == 0
        cl < 1.0f0 && (cl = 1.0f0)
        cl > hf && (cl = hf)
        cr = cl / hf
    end
    return cr
end

"""
    _cr_crown_tree(...) -> ICRI (new crown percent, integer)

Per-tree crown-ratio DUB/change from cr/crown.f. `icr_old`=entry ICR (0 ⇒ dub), `htg`/`dg` this cycle's
growth. Change is limited to ±1%/yr; CRNMLT/DLOW/DHI are the CRNMULT-keyword adjustments (defaults
1/0/99 ⇒ inert). Top-kill (itrunc>0, lstart) reduces the crown. Bounds [1 or 10, 95].
"""
function _cr_crown_tree(imodty::Int, sp::Int, d::Float32, h::Float32, pcti::Float32, tbau::Float32,
                        ba::Float32, hf::Float32, df::Float32, relden::Float32, icr_old::Int,
                        htg::Float32, fint::Float32, crnmlt::Float32, dlow::Float32, dhi::Float32,
                        itrunc::Int, normht::Int, lstart::Bool)::Int
    cr = cr_gemcr(imodty, sp, tbau, ba, hf, df, h, relden, pcti)
    crnew = cr * 100.0f0
    if !(lstart || icr_old == 0)
        chg = crnew - Float32(icr_old)
        pdifpy = chg / Float32(icr_old) / fint
        pdifpy > 0.01f0 && (chg = Float32(icr_old) * 0.01f0 * fint)
        pdifpy < -0.01f0 && (chg = Float32(icr_old) * (-0.01f0) * fint)
        crnew = (d >= dlow && d <= dhi) ? Float32(icr_old) + chg * crnmlt : Float32(icr_old) + chg
    end
    icri = trunc(Int, crnew + 0.5f0)
    if lstart || icr_old == 0
        (d >= dlow && d <= dhi) && (icri = trunc(Int, Float32(icri) * crnmlt))
    end
    if !(lstart || icr_old == 0)
        crln = h * Float32(icr_old) / 100.0f0
        crmax = (crln + htg) / (h + htg) * 100.0f0
        icri > crmax && (icri = trunc(Int, crmax + 0.5f0))
        (icri < 10 && crnmlt == 1.0f0) && (icri = trunc(Int, crmax + 0.5f0))
    end
    if lstart && itrunc != 0
        hn = Float32(normht) / 100.0f0
        hd = hn - Float32(itrunc) / 100.0f0
        cl = (Float32(icri) / 100.0f0) * hn - hd
        icri = trunc(Int, (cl * 100.0f0 / hn) + 0.5f0)
    end
    icri > 95 && (icri = 95)
    (icri < 10 && crnmlt == 1.0f0) && (icri = 10)
    icri < 1 && (icri = 1)
    return icri
end

"""
    crown_ratio_update!(s, ::CentralRockies; fint=10, lstart=false, relden_override=-1)

CR crown-ratio update (cr/crown.f): per tree, HF=H+HTG, DF=D+DG/bark, GEMCR → new crown %, with
±1%/yr change limiting and the CRMAX cap. `lstart` fills missing (ICR=0) inventory crowns. CRNMULT
keyword adjustments not yet wired (defaults inert). Writes `t.crown_pct` (ICR).
"""
function crown_ratio_update!(s::StandState, ::CentralRockies; fint::Float32 = 10.0f0,
                             lstart::Bool = false, relden_override::Float32 = -1.0f0, kwargs...)
    p, t, sd = s.plot, s.trees, s.coef.species
    t.n == 0 && return s
    imodty = Int(p.model_type)
    ba = p.basal_area
    relden = relden_override >= 0.0f0 ? relden_override : stand_ccf(s)
    bau = _cr_badist_bau(t)
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0.0f0 && continue
        icr_old = Int(t.crown_pct[i])
        (lstart && icr_old > 0) && continue        # inventory crown present → keep
        icr_old < 0 && (t.crown_pct[i] = Int32(-icr_old); continue)   # pest-computed → restore sign
        sp = Int(t.species[i])
        d = t.dbh[i]; h = t.height[i]; pcti = t.crown_ratio[i]
        icls = trunc(Int, d + 1.0f0); icls > 41 && (icls = 41)
        htg = t.ht_growth[i]; dg = t.diam_growth[i]
        bark = cr_bratio(sd, sp, d, imodty)
        hf = h + htg
        df = d + dg / bark; df < d && (df = d)
        t.crown_pct[i] = Int32(_cr_crown_tree(imodty, sp, d, h, pcti, bau[icls], ba, hf, df, relden,
            icr_old, htg, fint, 1.0f0, 0.0f0, 99.0f0, Int(t.trunc[i]), Int(t.norm_ht[i]), lstart))
    end
    return s
end

# =============================================================================
# cr_crown_width — CR crown competition factor / crown width (cr/ccfcal.f).
# CCF per tree: D≥10 → RD1+D·RD2+D²·RD3; D>0.1 → RDA·D^RDB; else 0.001.
# Crown width = sqrt(CCF/0.001803) (cap 99.9). IMAP = MAP{IMODTY}[sp]. The engine's
# CCF = 0.001803·cw²·tpa then recovers CCFT=CCF·tpa exactly. Ports the Wykoff-Crookston-
# Stage INT-133 Table-8 CCF equations mapped to CR species per model type.
# =============================================================================
const _CR_CCF_MAP = (
    (7,7,2,7,3,1,1,1,1,1,1,1,8,1,1,1,1,6,6,5,4,4,4,4,4,4,4,5,1,1,1,1,1,1,1,8,1,4),  # MAP1 IMODTY 1
    (3,3,2,3,3,1,1,1,1,1,1,7,8,1,1,1,1,1,1,5,4,4,6,6,6,6,6,5,1,1,1,1,7,7,7,8,1,4),  # MAP2 IMODTY 2
    (1,1,2,1,1,1,1,1,1,1,1,1,8,1,1,1,6,6,6,5,1,1,4,4,4,4,4,5,1,1,1,1,1,1,1,8,1,4),  # MAP3 IMODTY 3
    (7,7,2,7,7,1,1,1,1,1,1,1,1,1,1,1,6,6,6,5,4,4,4,4,4,4,4,5,1,1,1,1,1,1,1,1,1,4),  # MAP4 IMODTY 4
    (7,7,2,7,7,1,1,1,1,1,1,1,8,1,1,1,6,6,6,5,4,4,4,4,4,4,4,5,1,1,1,1,1,1,1,8,1,4),  # MAP5 IMODTY ≥5
)
const _CR_CCF_RD1 = (0.01925f0, 0.11f0, 0.04f0, 0.03f0, 0.03f0, 0.03f0, 0.03f0, 0.03f0)
const _CR_CCF_RD2 = (0.01676f0, 0.0333f0, 0.0270f0, 0.0215f0, 0.0238f0, 0.0173f0, 0.0216f0, 0.0180f0)
const _CR_CCF_RD3 = (0.00365f0, 0.00259f0, 0.00405f0, 0.00363f0, 0.00490f0, 0.00259f0, 0.00405f0, 0.00281f0)
const _CR_CCF_RDA = (0.009187f0, 0.017299f0, 0.015248f0, 0.011109f0, 0.008915f0, 0.007875f0, 0.011402f0, 0.007813f0)
const _CR_CCF_RDB = (1.7600f0, 1.5571f0, 1.7333f0, 1.7250f0, 1.7800f0, 1.7360f0, 1.7560f0, 1.7680f0)

@inline function cr_crown_width(sp::Int, d::Float32, imodty::Int)::Float32
    m = 1 <= imodty <= 5 ? imodty : 5
    imap = _CR_CCF_MAP[m][sp]
    ccf = if d >= 10.0f0
        _CR_CCF_RD1[imap] + d * _CR_CCF_RD2[imap] + d * d * _CR_CCF_RD3[imap]
    elseif d > 0.1f0
        _CR_CCF_RDA[imap] * fpow(d, _CR_CCF_RDB[imap])
    else
        0.001f0
    end
    cw = sqrt(ccf / 0.001803f0)
    cw > 99.9f0 && (cw = 99.9f0)
    return cw
end
