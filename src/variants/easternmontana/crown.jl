# =============================================================================
# crown.jl (easternmontana) — per-tree CCF (em/ccfcal.f MODE=1).
#
# EM CCF is the direct per-species polynomial (Paine-Hann / NI INT-133 form), summed
# over the tree list = stand CCF (like KT/IE ccfcal, not the CR crown-width→area path).
#   large:  CCFT = RD1 + D·RD2 + D²·RD3
#   small:  CCFT = RDA · D^RDB
# The large/small DBH threshold is species-group dependent (em/ccfcal.f SELECT CASE):
#   sp 5 (LL):                 D≥10 → poly, else RDA·D^RDB           (no 0.001 floor)
#   sp 1-4,6-10,12,17,18:      D≥1  → poly; D>0.1 → RDA·D^RDB; else 0.001
#   sp 11,13-16,19:            D≥10 → poly; D>0.1 → RDA·D^RDB; else 0.001
# Coefficients transcribed verbatim from em/ccfcal.f DATA (species 1=WB..19=OH).
# =============================================================================

const EM_RD1 = Float32[0.0186, 0.0392, 0.0388, 0.01925, 0.03, 0.01925, 0.01925, 0.03, 0.0172,
                       0.0219, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0204, 0.03]
const EM_RD2 = Float32[0.0146, 0.0180, 0.0269, 0.01676, 0.0216, 0.01676, 0.01676, 0.0173, 0.00876,
                       0.0169, 0.0215, 0.0238, 0.0215, 0.0215, 0.0215, 0.0215, 0.0238, 0.0246, 0.0215]
const EM_RD3 = Float32[0.00288, 0.00207, 0.00466, 0.00365, 0.00405, 0.00365, 0.00365, 0.00259, 0.00112,
                       0.00325, 0.00363, 0.00490, 0.00363, 0.00363, 0.00363, 0.00363, 0.00490, 0.0074, 0.00363]
const EM_RDA = Float32[0.009884, 0.007244, 0.017299, 0.009187, 0.011402, 0.009187, 0.009187, 0.007875,
                       0.011402, 0.007813, 0.011109, 0.008915, 0.011109, 0.011109, 0.011109, 0.011109,
                       0.008915, 0.011109, 0.011109]
const EM_RDB = Float32[1.6667, 1.8182, 1.5571, 1.7600, 1.7560, 1.7600, 1.7600, 1.7360, 1.7560,
                       1.7780, 1.7250, 1.7800, 1.7250, 1.7250, 1.7250, 1.7250, 1.7800, 1.7250, 1.7250]

# em/ccfcal.f MODE=1: per-tree CCF (before the ×P expansion the caller applies).
@inline function em_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    poly()  = EM_RD1[sp] + d * EM_RD2[sp] + d * d * EM_RD3[sp]
    small() = EM_RDA[sp] * d ^ EM_RDB[sp]
    if sp == 5                                   # LL: D≥10 poly, else small (no floor)
        return d >= 10f0 ? poly() : small()
    elseif sp == 11 || (13 <= sp <= 16) || sp == 19   # GA/CW/BA/PW/NC/OH: D≥10 threshold
        return d >= 10f0 ? poly() : (d > 0.1f0 ? small() : 0.001f0)
    else                                          # 1-4,6-10,12,17,18: D≥1 threshold
        return d >= 1f0 ? poly() : (d > 0.1f0 ? small() : 0.001f0)
    end
end

# =============================================================================
# EM crown-ratio update (em/crown.f) — the IE/NI Weibull/DCR change-in-crown model.
# Single PARM(14) vector (IE/NI sp9 form, applied to all species) + habitat intercept
# CRHAB[MAPHAB[ITYPE]] + CRSD=6.35. Mirrors KT's crown_ratio_update! DCR structure
# (crown change = exp(PCR) − exp(DCR-backdated), bounded ±1%/yr; d<3 keeps crown while cycling).
# =============================================================================
# em/crown.f PARM(14): 1-6 density(BA,BA²,lnBA,RELDEN,RELDEN²,lnRELDEN), 7-14(D,D²,lnD,H,H²,lnH,P,lnP)
const EM_CRPARM = Float32[-0.00190, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.23372, 0.0, 0.0, -0.28433, 0.001903, 0.0]
const EM_CRHAB  = Float32[0.09453, -0.07740, 0.07113, 0.2039, 0.06176, 0.1513, 0.09086, 0.1580, 0.09229, 0.01551, 0.0, 0.0, 0.0, 0.0]
const EM_CR_MAPHAB = Int[2,2,2,2,2,2,2,2,2,2,2,2,2,3,4,4,4,4,5,6,6,7,6,1,8,1,9,10,10,6]  # ITYPE→CRHAB idx (13*2,3,4*4,5,6,6,7,6,1,8,1,9,2*10,6)
const _EM_CRSD = 6.35f0

function crown_ratio_update!(s::StandState, ::EasternMontana; fint::Float32 = 10.0f0, lstart::Bool = false, kwargs...)
    p, t = s.plot, s.trees
    t.n == 0 && return s
    itype = Int(p.habitat_input); it = (1 <= itype <= 30) ? itype : 1
    crcon = EM_CRHAB[clamp(EM_CR_MAPHAB[it], 1, 14)]
    ba = p.basal_area; relden = p.relative_density
    lnba = ba > 0f0 ? log(ba) : 0f0; lnrd = relden > 0f0 ? log(relden) : 0f0
    reldm1 = p.relative_density_prev; oba = p.old_ba
    if reldm1 < 100f0; oba = ba; reldm1 = relden; end
    x1 = (!lstart && oba > 0f0) ? log(oba) : 0f0
    x2 = (!lstart && reldm1 > 0f0) ? log(reldm1) : 0f0
    dgsd = s.control.dg_sd
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    P = EM_CRPARM
    nlim = t.n + (lstart ? Int(t.ndead) : 0)
    @inbounds for i in 1:nlim
        t.tpa[i] <= 0f0 && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue
        icr < 0 && (t.crown_pct[i] = Int32(-icr); continue)
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        bark = bark_ratio(ba_a, ba_b, sp, d)
        local icri::Int
        # em/crown.f: the DCR change-in-crown model applies to the lstart DUB (all sizes) and to CYCLING
        # d≥3; cycling d<3 keeps its crown. (Previously d<3 at lstart used a flat-40 placeholder — never
        # exercised because EM wasn't calling the lstart dub at all. Now that the dub is wired (simulate.jl),
        # a 0.1" seedling with ICR=0 dubs via the same model ⇒ CR≈50-79 like live's DUBSCR, activating the
        # _em_smhtgf beta2·cr height term instead of the crown_pct=0 freeze (#137).)
        if d < 3.0f0 && !lstart
            continue                                  # cycling: D<3 keeps its crown
        end
        xcrcon = crcon + P[1]*ba + P[2]*ba*ba + P[3]*lnba + P[4]*relden + P[5]*relden*relden + P[6]*lnrd
        pp = t.crown_ratio[i]; pp < 0.01f0 && (pp = 0.01f0)
        pcr = xcrcon + P[7]*d + P[8]*d*d + P[9]*log(d) + P[10]*h + P[11]*h*h + P[12]*log(h) + P[13]*pp + P[14]*log(pp)
        exppcr = exp(pcr)
        if lstart
            icri = trunc(Int, icr + exppcr*100f0 + 0.50005f0)
            dgsd >= 1.0f0 && (icri = trunc(Int, bachlo(s.rng, Float32(icri), _EM_CRSD)))
        else
            dcrcon = crcon + P[1]*oba + P[2]*oba*oba + P[3]*x1 + P[4]*reldm1 + P[5]*reldm1*reldm1 + P[6]*x2
            db = d - t.diam_growth[i]/bark; db <= 0f0 && (db = d)
            hb = h - t.ht_growth[i]; hb <= 0f0 && (hb = h)
            pb = t.crown_ratio[i]; pb < 0.01f0 && (pb = 0.01f0)
            dcr = dcrcon + P[7]*db + P[8]*db*db + P[9]*log(db) + P[10]*hb + P[11]*hb*hb + P[12]*log(hb) + P[13]*pb + P[14]*log(pb)
            chg = exppcr - exp(dcr)
            if icr > 0
                pdifpy = chg / Float32(icr) / fint * 100f0
                pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint / 100f0)
                pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint / 100f0)
            end
            icri = trunc(Int, Float32(icr) + chg*100f0 + 0.50005f0)
        end
        icri > 95 && (icri = 95); icri < 5 && (icri = 5)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# =============================================================================
# em_cwcalc — EM forest-grown crown width (base/cwidth.f -> cwcalc.f, IWHO=0; the
# western Bechtold/Crookston/Donnelly crown-width library).  EM is Region 1 =>
# KODFOR<601 => BF=1 and the Region-6 forest-BF section is skipped (GO TO 10).
# EMMAP maps FVS EM species 1..19 -> a 5-char CWEQN (FIA code + eqn#).  This is the
# CRWDTH array that COVER's CVCW reads (NOT the em_tree_ccf CCF polynomial above,
# which is the separate CCF path).  Reuses the shared CR kernels (_cr_r6m2 / _cr_bech1
# / _cr_bech2 / _cr_hopkins) plus the R1 (Crookston Region-1) and pure-power helpers
# below.  Math faithful to gfortran: fpow/fexp/flog + left-to-right op order; the
# Hopkins point is the shared WESTERN one.  Dump-replay bit-exact vs FVSem_g16 (all 15
# EM CWEQN forms, incl. small-tree scaling + juniper D>=25 plateau).
# REPORTING-ONLY (COVER); INERT for growth/mortality/CCF/volume.
# =============================================================================
const _EM_CWMAP = ("10105","07303","20203","11301","07204","06602","10803","09303",
                   "01903","12203","74902","74605","74705","74902","74902","74902",
                   "37506","26405","74902")

# Transcendentals routed DIRECTLY through glibc libm (logf/expf/powf), like WSBWE — the
# gfortran-16 oracle's REAL*4 EXP/LOG/** resolve to these, so this is bit-exact vs
# FVSem_g16 (verified 353/353).  The core/fmath shim is built by a *different* gfortran and
# is ~1 ULP off here, so it is intentionally NOT used on this crown-width path.
const _EMCW_LIBM = "libm.so.6"
@inline _emcw_log(x::Float32) = ccall((:logf, _EMCW_LIBM), Float32, (Float32,), x)
@inline _emcw_exp(x::Float32) = ccall((:expf, _EMCW_LIBM), Float32, (Float32,), x)
@inline _emcw_pow(x::Float32, y::Float32) = ccall((:powf, _EMCW_LIBM), Float32, (Float32,Float32), x, y)

# Crookston Region-1 form: mult*EXP(c0 + ccl*ln(CL) + cd*ln(Dm) + ch*ln(H) + cba*ln(BAREA)).
# Only the ln(D) term floors to `dfloor`; CL/H/BAREA use actual values.  Small-tree scale
# ×(D/dfloor).  Terms with a 0 coefficient add +0 (bit-exact no-op).
@inline function _em_r1(mult::Float32, c0::Float32, ccl::Float32, cd::Float32, ch::Float32,
                        cba::Float32, dfloor::Float32, cap::Float32,
                        d::Float32, h::Float32, cl::Float32, barea::Float32)::Float32
    dm = d >= dfloor ? d : dfloor
    v = mult * _emcw_exp(c0 + ccl*_emcw_log(cl) + cd*_emcw_log(dm) + ch*_emcw_log(h) + cba*_emcw_log(barea))
    d < dfloor && (v *= d / dfloor)
    v > cap && (v = cap)
    return v
end

# Pure power form a*D^b (Crookston R6 model-1 07204, Donnelly 37506).  OMIND=1 scaling.
@inline function _em_powf(a::Float32, b::Float32, cap::Float32, d::Float32)::Float32
    dm = d >= 1f0 ? d : 1f0
    v = a * _emcw_pow(dm, b)
    d < 1f0 && (v *= d)
    v > cap && (v = cap)
    return v
end

# Crookston R6 model 2: a·D^b·H^c·CL^dd·(BAREA+1)^e·EXP(EL)^f (BF=1).  OMIND=1 small-tree
# scaling; EL∈[ello,elhi].  e=0/f=0 → the factor is ^0 = ×1.0 (bit-exact no-op).
@inline function _em_r6m2(a::Float32, b::Float32, c::Float32, dd::Float32, e::Float32, f::Float32,
                          d::Float32, h::Float32, cl::Float32, ba1::Float32, el::Float32,
                          ello::Float32, elhi::Float32, cap::Float32)::Float32
    elc = el < ello ? ello : (el > elhi ? elhi : el)
    dm = d >= 1f0 ? d : 1f0
    cw = a * _emcw_pow(dm, b) * _emcw_pow(h, c) * _emcw_pow(cl, dd) *
         _emcw_pow(ba1, e) * _emcw_pow(_emcw_exp(elc), f)
    d < 1f0 && (cw *= d)
    cw > cap && (cw = cap)
    return cw
end

# Bechtold 2004 model 1: a + b·D (MIND=5 small-tree scaling).
@inline function _em_bech1(a::Float32, b::Float32, d::Float32, cap::Float32)::Float32
    dm = d >= 5f0 ? d : 5f0
    cw = a + b * dm
    d < 5f0 && (cw *= d / 5f0)
    cw > cap && (cw = cap)
    return cw
end

# Bechtold 2004 model 2: a + b·D + c·D² + crc·CR + hic·HI.  HI∈[hlo,hhi]; MIND=5 small-tree
# scaling; optional D≥25 plateau (dcap25).
@inline function _em_bech2(a::Float32, b::Float32, c::Float32, crc::Float32, hic::Float32,
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

function em_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                   el::Float32, hi::Float32)::Float32
    (1 <= sp <= 19) || return 0.5f0
    barea <= 1f0 && (barea = 1f0)          # cwcalc.f: IF(BAREA.LE.1.) BAREA=1.
    cl  = cr * h * 0.01f0                   # CL = CR*H*0.01
    ba1 = barea + 1f0
    eqn = _EM_CWMAP[sp]
    cw = if eqn == "10105";     _em_r6m2(2.2354f0,0.66680f0,-0.11658f0,0.16927f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,40f0)
    elseif eqn == "74605";      _em_r6m2(4.7961f0,0.64167f0,-0.18695f0,0.18581f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,45f0)
    elseif eqn == "74705";      _em_r6m2(4.4327f0,0.41505f0,-0.23264f0,0.41477f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,56f0)
    elseif eqn == "26405";      _em_r6m2(3.7854f0,0.54684f0,-0.12954f0,0.16151f0,0.03047f0,-0.00561f0, d,h,cl,ba1,el,10f0,79f0,45f0)
    elseif eqn == "07303";      _em_r1(1.02478f0,0.99889f0,0.19422f0,0.59423f0,-0.09078f0,-0.02341f0, 1f0,40f0, d,h,cl,barea)
    elseif eqn == "20203";      _em_r1(1.01685f0,1.48372f0,0.27378f0,0.49646f0,-0.18669f0,-0.01509f0, 1f0,80f0, d,h,cl,barea)
    elseif eqn == "10803";      _em_r1(1.03992f0,1.58777f0,0.30812f0,0.64934f0,-0.38964f0,0f0, 0.7f0,40f0, d,h,cl,barea)
    elseif eqn == "09303";      _em_r1(1.02687f0,1.28027f0,0.2249f0,0.47075f0,-0.15911f0,0f0, 0.1f0,40f0, d,h,cl,barea)
    elseif eqn == "01903";      _em_r1(1.02886f0,1.01255f0,0.30374f0,0.37093f0,-0.13731f0,0f0, 0.1f0,30f0, d,h,cl,barea)
    elseif eqn == "12203";      _em_r1(1.02687f0,1.49085f0,0.1862f0,0.68272f0,-0.28242f0,0f0, 2f0,46f0, d,h,cl,barea)
    elseif eqn == "07204";      _em_powf(2.2586f0,0.68532f0,33f0, d)
    elseif eqn == "37506";      _em_powf(5.8980f0,0.4841f0,25f0, d)
    elseif eqn == "11301";      _em_bech1(4.0181f0,0.8528f0, d,25f0)
    elseif eqn == "06602";      _em_bech2(-4.1599f0,1.3528f0,-0.0233f0,0.0633f0,-0.0423f0, d,cr,hi,-37f0,19f0,29f0,true)
    elseif eqn == "74902";      _em_bech2(4.1687f0,1.5355f0,0f0,0f0,0.1275f0, d,cr,hi,-26f0,-2f0,35f0,false)
    else 0f0 end
    # cwcalc.f final CRWDTH clamp (after label 9000).
    cw < 0.5f0 && (cw = 0.5f0)
    cw > 99.9f0 && (cw = 99.9f0)
    return cw
end
