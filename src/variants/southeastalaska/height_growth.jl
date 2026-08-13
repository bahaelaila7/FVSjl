# =============================================================================
# height_growth.jl (southeastalaska) — AK large-tree height growth (ak/htgf.f). Chunk 4.
#
#   height_growth!(::SoutheastAlaska; scale) — per-tree HTG (ak/htgf.f main body).
#
# AK is a SINGLE Wykoff-style height equation for ALL species (no per-species branches):
#   BASEHG = exp(B1 + B2·D² + B3·ln(D) + B4·(ELEV·100) + B5·ln(SITEAR) + B6·ln(DG10))   [annual]
#   DG10   = DGLT / BRAT   (outside-bark 10-yr DG; DGLT = the tree's inside-bark DG this cycle);
#            DG10 forced 0.1 when H < 4.5 ft.
# then the same PERMAFROST modifier structure as dgf! (species 4:7,13,16:23):
#   PFHMOD = exp(B1[+B2 if LPERM]+B3·D²+B4·ln(D)+B5·ELEV+B6·ln(DG10))/BASEHG; LPERM⇒cap≤1, else floor≥1.
#   POTHTG = YR·BASEHG·PFHMOD;  then species mult (AD/WI/OH→0.45, SU→0.65, else 1.0).
# then HTLO/HTHI bounding (HGBND∈[0.1,1]), min 0.1, and DG≤0.04 ⇒ HTG=0.1; finally
#   HTG = SCALE·XHT·HTG·EXP(HTCON)·MISHGF, SCALE=FINT/YR (engine `scale`), size-cap.
#
# VALIDATED BIT-EXACT vs the live FVSak `DEBUG HTGF` dump on akt01 (see docs/AK_VARIANT_PORT_AUDIT.md).
# akt01 is all non-permafrost species (PFHMOD=1). The engine supplies DGLT = t.diam_growth[i] (the DDS→DG
# inside-bark increment) exactly as ak/htgf.f reads DG(I) after DGF.
# =============================================================================

# ak/htgf.f DATA — base HG coefficients (NOPERM*), permafrost-mod coefficients (PERM*), HTLO/HTHI bounds.
const AK_HG_B1 = Float32[-2.18236,-2.18236,-2.955092,-0.663746,-0.501134,-0.501134,-0.663746,-2.18236,-3.19512,-2.979661,-2.404892,-2.770793,-0.501134,-0.174905,-0.710501,-0.306839,-0.306839,-0.174905,-0.630114,-0.174905,-0.174905,-0.174905,-0.174905]
const AK_HG_B2 = Float32[-0.000447,-0.000447,-0.000447,-0.003814,-0.003814,-0.003814,-0.003814,-0.000447,-0.000447,-0.000447,-0.000447,-0.000447,-0.003814,-0.003814,-0.003814,-0.003814,-0.003814,-0.003814,-0.003814,-0.003814,-0.003814,-0.003814,-0.003814]
const AK_HG_B3 = Float32[-0.00488,-0.00488,0.125393,0.176358,0.209621,0.209621,0.176358,-0.00488,0.166453,0.130947,0.044809,0.099358,0.209621,-0.051622,0.396143,0.062956,0.062956,-0.051622,0.068496,-0.051622,-0.051622,-0.051622,-0.051622]
const AK_HG_B4 = Float32[0.0,0.0,0.0,-0.000075,-0.000075,-0.000075,-0.000075,0.0,0.0,0.0,0.0,0.0,-0.000075,-0.000075,-0.000075,-0.000075,-0.000075,-0.000075,-0.000075,-0.000075,-0.000075,-0.000075,-0.000075]
const AK_HG_B5 = Float32[0.429243,0.429243,0.429243,0.0,0.0,0.0,0.0,0.429243,0.429243,0.429243,0.429243,0.429243,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const AK_HG_B6 = Float32[0.489281,0.489281,0.376732,0.560124,0.517799,0.517799,0.560124,0.489281,0.342411,0.374658,0.469563,0.408836,0.517799,0.539313,0.195472,0.392915,0.392915,0.539313,0.367413,0.539313,0.539313,0.539313,0.539313]

const AK_HGP_B1 = Float32[0.0,0.0,0.0,-0.661722,-0.509243,-0.509243,-0.661722,0.0,0.0,0.0,0.0,0.0,-0.509243,0.0,0.0,-0.299268,-0.299268,-0.148618,-0.617425,-0.148618,-0.148618,-0.148618,-0.148618]
const AK_HGP_B2 = Float32[0.0,0.0,0.0,-0.130016,-0.130016,-0.130016,-0.130016,0.0,0.0,0.0,0.0,0.0,-0.130016,0.0,0.0,-0.130016,-0.130016,-0.130016,-0.130016,-0.130016,-0.130016,-0.130016,-0.130016]
const AK_HGP_B3 = Float32[0.0,0.0,0.0,-0.003807,-0.003807,-0.003807,-0.003807,0.0,0.0,0.0,0.0,0.0,-0.003807,0.0,0.0,-0.003807,-0.003807,-0.003807,-0.003807,-0.003807,-0.003807,-0.003807,-0.003807]
const AK_HGP_B4 = Float32[0.0,0.0,0.0,0.17915,0.221303,0.221303,0.17915,0.0,0.0,0.0,0.0,0.0,0.221303,0.0,0.0,0.057156,0.057156,-0.075692,0.060586,-0.075692,-0.075692,-0.075692,-0.075692]
const AK_HGP_B5 = Float32[0.0,0.0,0.0,-0.000067,-0.000067,-0.000067,-0.000067,0.0,0.0,0.0,0.0,0.0,-0.000067,0.0,0.0,-0.000067,-0.000067,-0.000067,-0.000067,-0.000067,-0.000067,-0.000067,-0.000067]
const AK_HGP_B6 = Float32[0.0,0.0,0.0,0.537482,0.509459,0.509459,0.537482,0.0,0.0,0.0,0.0,0.0,0.509459,0.0,0.0,0.397696,0.397696,0.556616,0.375558,0.556616,0.556616,0.556616,0.556616]

const AK_HTLO = Float32[131.,83.,92.,47.,78.,78.,47.,131.,70.,102.,130.,83.,78.,38.,86.,77.,77.,84.,76.,84.,38.,55.,38.]
const AK_HTHI = Float32[211.,150.,139.,77.,123.,123.,77.,211.,98.,151.,214.,150.,123.,59.,120.,102.,102.,130.,102.,130.,59.,85.,59.]

# ak/cratet.f INVENTORY-EQN (Curtis-Arney) height-diameter dub for missing heights (LHTDRG=false ⇒ always
# used): H = 4.5 + HTT11·(1−exp(HTT12·D))^HTT13, then ×0.45 (AD/WI/OH) / ×0.65 (SU), floored 4.5; D≤0.1 ⇒ 1.01.
# MEASURED from FVSak DEBUG CRATET (INVENTORY EQN DUBBING dump, all 23 species; 5-decimal print precision).
const AK_HTT11 = Float32[173.57806,173.57806,118.00596,68.53457,115.99145,115.99145,68.53457,173.57806,71.73042,105.74110,141.21040,108.15884,115.99145,115.11522,131.58195,59.01998,59.01998,115.11522,69.44218,115.11522,115.11522,115.11522,115.11522]
const AK_HTT12 = Float32[-0.03470,-0.03470,-0.04487,-0.14965,-0.06808,-0.06808,-0.14965,-0.03470,-0.09667,-0.05192,-0.05410,-0.05487,-0.06808,-0.04711,-0.03474,-0.33004,-0.33004,-0.04711,-0.18146,-0.04711,-0.04711,-0.04711,-0.04711]
const AK_HTT13 = Float32[1.06697,1.06697,1.09078,1.33524,1.05922,1.05922,1.33524,1.06697,1.58941,1.13446,1.24836,1.30982,1.05922,0.79514,0.78304,1.48201,1.48201,0.79514,1.17920,0.79514,0.79514,0.79514,0.79514]

"ak/cratet.f Curtis-Arney missing-height dub (pre-4.5-floor; caller applies the floor for D>0.1)."
@inline function ak_htdbh_dub(sp::Int, d::Float32)::Float32
    d <= 0.1f0 && return 1.01f0
    h = 4.5f0 + AK_HTT11[sp] * (1f0 - exp(AK_HTT12[sp] * d))^AK_HTT13[sp]
    m = (sp == 14 || sp == 21 || sp == 23) ? 0.45f0 : (sp == 22 ? 0.65f0 : 1.0f0)
    return h * m
end

# ak/htgf.f main body — per-tree HTG. `scale` = FINT/YR (engine passes fint/htg_period = 1 for the 10-yr cycle).
function height_growth!(s::StandState, ::SoutheastAlaska; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    yr = 10.0f0                                  # ak/htgf.f YR = 10 (period the annual HG is expanded to)
    temel = p.elevation * 100f0                  # TEMEL = ELEV·100
    lperm = false                                # LPERM (PERMAFROST keyword) — off by default (later chunk)
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); h = t.height[i]; d = t.dbh[i]
        xsite = p.sp_site_index[sp]              # XSITE = SITEAR(ISPC)
        dglt = t.diam_growth[i]                  # DGLT = DG(I) (inside-bark increment this cycle)
        dgchk = dglt <= 0f0 ? 0.0001f0 : dglt    # ak/htgf.f:245 IF(DG.LE.0)DG=0.0001 (used in the DG≤0.04 test)
        brat = ak_bratio(sp, d)
        dg10 = dglt / brat                       # outside-bark DG
        h < 4.5f0 && (dg10 = 0.1f0)
        basehg = exp(AK_HG_B1[sp] + AK_HG_B2[sp]*d*d + AK_HG_B3[sp]*log(d) +
                     AK_HG_B4[sp]*temel + AK_HG_B5[sp]*log(xsite) + AK_HG_B6[sp]*log(dg10))
        pfhmod = 1.0f0
        if sp in AK_PERM_SP
            if lperm
                pfhmod = exp(AK_HGP_B1[sp] + AK_HGP_B2[sp] + AK_HGP_B3[sp]*d*d + AK_HGP_B4[sp]*log(d) +
                             AK_HGP_B5[sp]*temel + AK_HGP_B6[sp]*log(dg10)) / basehg
                pfhmod > 1f0 && (pfhmod = 1f0)
            else
                pfhmod = exp(AK_HGP_B1[sp] + AK_HGP_B3[sp]*d*d + AK_HGP_B4[sp]*log(d) +
                             AK_HGP_B5[sp]*temel + AK_HGP_B6[sp]*log(dg10)) / basehg
                pfhmod < 1f0 && (pfhmod = 1f0)
            end
        end
        pothtg = yr * basehg * pfhmod
        pothtg *= (sp == 14 || sp == 21 || sp == 23) ? 0.45f0 : (sp == 22 ? 0.65f0 : 1.0f0)
        # HTLO/HTHI bounding
        hgbnd = if h >= AK_HTLO[sp] && h < AK_HTHI[sp]
            b = 1f0 - (h - AK_HTLO[sp]) / (AK_HTHI[sp] - AK_HTLO[sp]); b < 0.1f0 ? 0.1f0 : b
        elseif h < AK_HTLO[sp]
            1f0
        else
            0.1f0
        end
        htg = pothtg * hgbnd
        htg <= 0.1f0 && (htg = 0.1f0)
        dgchk <= 0.04f0 && (htg = 0.1f0)
        # SCALE·XHT·HTG·EXP(HTCON)·MISHGF (XHT=1, HTCON=htg_cor, MISHGF=1)
        htg = scale * htg * exp(c.htg_cor[sp])
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
    end
    return s
end
