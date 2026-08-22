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
    # At LSTART, also dub the DEAD partition's crown (cratet.f processes IREC2..MAXTRE too) so the cycle-0
    # FVS_TreeList dead records carry the crown ratio live shows (PctCr) — which the forest crown width
    # (cr_cwcalc CL term) also needs. Dead crowns feed nothing else (per-tree write; density is live-only).
    nlim = t.n + (lstart ? Int(t.ndead) : 0)
    @inbounds for i in 1:nlim
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

# =============================================================================
# cr_cwcalc — CR forest-grown crown width for the FVS_TreeList CrWidth column
# (base/cwidth.f → cwcalc.f, IWHO=0; western Bechtold/Crookston crown-width library).
# CR is Region 2/3 ⇒ BF=1 and the Region-6 forest-specific BF section is skipped
# (KODFOR<601 GO TO 10). CRMAP maps FVS species index 1-38 → a 5-char CWEQN (FIA+eqn#).
# REPORTING-ONLY: not used in growth/mortality/CCF/volume — the CCF open-grown crown
# width is the separate ccfcal path (cr_crown_width). Ported so the treelist CrWidth
# column matches live (it was the crown_width 0.5 default for all CR rows). Math faithful
# to gfortran: fpow/fexp/flog + left-to-right op order. HI uses the WESTERN Hopkins point.
# =============================================================================
const _CR_CWMAP = ("01905","01801","20205","01703","01505","26403","24205","07303","10201","11301",
                   "10805","10602","12205","10105","11905","06602","09305","09305","09305","74605",
                   "74902","74902","81402","81402","81402","81402","81402","74605","06602","06602",
                   "06602","06602","10602","10602","10602","12205","12205","74902")

@inline function _cr_hopkins(lat::Float32, long::Float32, elev::Float32)::Float32
    hilong = -abs(long)
    hielev = elev * 100f0
    return ((hielev - 5449f0) / 100f0) * 1f0 + (lat - 42.16f0) * 4f0 + (-116.39f0 - hilong) * 1.25f0
end

# Crookston R6 model 2: a·D^b·H^c·CL^dd·(BAREA+1)^e·EXP(EL)^f (BF=1). e=0/f=0 drop the
# BAREA / EL factor (×1.0f0 = bit-exact no-op). OMIND=1 small-tree scaling; EL∈[ello,elhi].
@inline function _cr_r6m2(a::Float32, b::Float32, c::Float32, dd::Float32, e::Float32, f::Float32,
                          d::Float32, h::Float32, cl::Float32, ba1::Float32, el::Float32,
                          ello::Float32, elhi::Float32, cap::Float32)::Float32
    elc = el < ello ? ello : (el > elhi ? elhi : el)
    dm = d >= 1f0 ? d : 1f0
    cw = a * fpow(dm, b) * fpow(h, c) * fpow(cl, dd) * fpow(ba1, e) * fpow(fexp(elc), f)
    d < 1f0 && (cw *= d)                 # ×(D/OMIND), OMIND=1
    cw > cap && (cw = cap)
    return cw
end

# Bechtold 2004 model 2: a + b·D + c·D² + crc·CR + hic·HI. HI∈[hlo,hhi]; MIND=5 small-tree
# scaling; optional D≥25 plateau (dcap25). Coefficients that are 0 drop their term (×0 exact).
@inline function _cr_bech2(a::Float32, b::Float32, c::Float32, crc::Float32, hic::Float32,
                           d::Float32, cr::Float32, hi::Float32, hlo::Float32, hhi::Float32,
                           cap::Float32, dcap25::Bool)::Float32
    hv = hi < hlo ? hlo : (hi > hhi ? hhi : hi)
    dm = d >= 5f0 ? d : 5f0
    cw = a + b * dm + c * dm * dm + crc * cr + hic * hv
    d < 5f0 && (cw *= d / 5f0)
    (dcap25 && d >= 25f0) && (cw = a + b * 25f0 + c * 25f0 * 25f0 + crc * cr + hic * hv)
    cw > cap && (cw = cap)
    return cw
end

# Bechtold 2004 model 1: a + b·D (MIND=5 small-tree scaling).
@inline function _cr_bech1(a::Float32, b::Float32, d::Float32, cap::Float32)::Float32
    dm = d >= 5f0 ? d : 5f0
    cw = a + b * dm
    d < 5f0 && (cw *= d / 5f0)
    cw > cap && (cw = cap)
    return cw
end

function cr_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                   el::Float32, hi::Float32; bf::Float32 = 1f0)::Float32
    (1 <= sp <= 38) || return 0f0
    eqn = _CR_CWMAP[sp]
    cl  = cr * h * 0.01f0
    ba1 = barea + 1f0
    cw = if eqn == "01905"; _cr_r6m2(5.8827f0*bf,0.51479f0,-0.21501f0,0.17916f0,0.03277f0,-0.00828f0, d,h,cl,ba1,el,10f0,85f0,30f0)
    elseif eqn == "20205"; _cr_r6m2(6.0227f0*bf,0.54361f0,-0.20669f0,0.20395f0,-0.00644f0,-0.00378f0, d,h,cl,ba1,el,1f0,75f0,80f0)
    elseif eqn == "01505"; _cr_r6m2(5.0312f0*bf,0.53680f0,-0.18957f0,0.16199f0,0.04385f0,-0.00651f0, d,h,cl,ba1,el,2f0,75f0,35f0)
    elseif eqn == "24205"; _cr_r6m2(6.2382f0*bf,0.29517f0,-0.10673f0,0.23219f0,0.05341f0,-0.00787f0, d,h,cl,ba1,el,1f0,72f0,45f0)
    elseif eqn == "10805"; _cr_r6m2(6.6941f0*bf,0.81980f0,-0.36992f0,0.17722f0,-0.01202f0,-0.00882f0, d,h,cl,ba1,el,1f0,79f0,40f0)
    elseif eqn == "12205"; _cr_r6m2(4.7762f0*bf,0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el,13f0,75f0,50f0)
    elseif eqn == "11905"; _cr_r6m2(5.3822f0*bf,0.57896f0,-0.19579f0,0.14875f0,0f0,-0.00685f0, d,h,cl,ba1,el,10f0,75f0,35f0)
    elseif eqn == "09305"; _cr_r6m2(6.7575f0*bf,0.55048f0,-0.25204f0,0.19002f0,0f0,-0.00313f0, d,h,cl,ba1,el,1f0,85f0,40f0)
    elseif eqn == "10105"; _cr_r6m2(2.2354f0*bf,0.66680f0,-0.11658f0,0.16927f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,40f0)
    elseif eqn == "74605"; _cr_r6m2(4.7961f0*bf,0.64167f0,-0.18695f0,0.18581f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,45f0)
    elseif eqn == "01703"
        dm = d >= 1f0 ? d : 1f0
        v = 1.0303f0 * fexp(1.14079f0 + 0.20904f0*flog(cl) + 0.38787f0*flog(dm))
        d < 1f0 && (v *= d); v > 40f0 && (v = 40f0); v
    elseif eqn == "07303"
        dm = d >= 1f0 ? d : 1f0
        v = 1.02478f0 * fexp(0.99889f0 + 0.19422f0*flog(cl) + 0.59423f0*flog(dm) +
                             (-0.09078f0)*flog(h) + (-0.02341f0)*flog(barea))
        d < 1f0 && (v *= d); v > 40f0 && (v = 40f0); v
    elseif eqn == "26403"
        v = if h < 5f0
            (0.8f0*h*max(0.5f0,cr*0.01f0)) * (1f0-(h-5f0)*0.1f0) * 6.90396f0 *
                fpow(d,0.55645f0) * fpow(h,-0.28509f0) * fpow(cl,0.20430f0) * (h-5f0) * 0.1f0
        elseif h >= 15f0
            6.90396f0 * fpow(d,0.55645f0) * fpow(h,-0.28509f0) * fpow(cl,0.20430f0)
        else
            0.8f0*h*max(0.5f0,cr*0.01f0)
        end
        v > 45f0 && (v = 45f0); v
    elseif eqn == "10602"; _cr_bech2(-5.4647f0,1.9660f0,-0.0395f0,0.0427f0,-0.0259f0, d,cr,hi,-40f0,11f0,25f0,true)
    elseif eqn == "06602"; _cr_bech2(-4.1599f0,1.3528f0,-0.0233f0,0.0633f0,-0.0423f0, d,cr,hi,-37f0,19f0,29f0,true)
    elseif eqn == "74902"; _cr_bech2(4.1687f0,1.5355f0,0f0,0f0,0.1275f0, d,cr,hi,-26f0,-2f0,35f0,false)
    elseif eqn == "81402"; _cr_bech2(0.3309f0,0.8918f0,0f0,0.0510f0,0f0, d,cr,hi,-1f30,1f30,19f0,false)
    elseif eqn == "01801"; _cr_bech1(6.073f0,0.3756f0, d,15f0)
    elseif eqn == "10201"; _cr_bech1(7.4251f0,0.8991f0, d,25f0)
    elseif eqn == "11301"; _cr_bech1(4.0181f0,0.8528f0, d,25f0)
    else 0f0 end
    return cw
end
