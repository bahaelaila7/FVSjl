# =============================================================================
# small_tree_growth.jl (oregoncoast) — OC/OP REGENT small-tree height + diameter model.
#
# Ported from oc/regent.f + oc/smhtgf.f + oc/htdbh.f. FVS routes trees with DBH < XMAX(sp) (=4.0",
# sp23/50=10.0) to REGENT, which grows HEIGHT from a species-mapped small-tree height-age model
# (SMHTGF) — modified by a crown/relative-height vigor term that the large-tree HTGF drops when
# PCCF<100 — and then DERIVES the DBH from the resulting height (H-D function). The IORG=0 seedlings
# were previously routed to the large-tree oc_htgf_native (organon_hook.jl), which over-grew their
# height ~7ft/cycle (the RELHT suppression is removed there when PCCF<100) — see docs/OC_VARIANT_PORT_AUDIT.md.
#
# OC/OP are DGSD=0 ⇒ regent.f:234 ZZRAN=0 ⇒ the inventory small-tree path is fully DETERMINISTIC.
# Validated vs the live FVSoc_clean scoped-DEBUG REGENT dump (oct01 stand-1, scratchpad/oc/): DF small
# trees CON=1.0 ⇒ HTGRR=HTGR=1.8110 (D=0.1,H=2.0); the large-tree htgf estimate for those trees = 9.05.
# =============================================================================

# smhtgf.f DATA: SPADJF (species/site adjustment) + MAPSP (species → 1 of 5 height eqs). CA/OC 50-sp order.
const OC_SMHT_SPADJF = Float32[
    1.0,1.0,0.9,1.1,1.1,1.1,1.1,0.8,0.9,0.9, 1.0,1.0,1.0,1.0,1.0,1.1,1.1,1.0,1.1,0.9,
    1.0,0.9,1.0,0.8,1.0,1.1,0.9,1.1,1.1,1.0, 1.1,1.0,1.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    1.1,1.0,1.1,1.2,1.2,1.1,0.8,1.0,1.0,1.0]
const OC_SMHT_MAPSP = Int[
    2,2,2,2,2,2,2,2,2,1, 1,1,1,1,1,1,1,1,1,1,
    1,2,5,2,1,3,3,3,3,3, 3,3,3,4,3,3,4,3,4,3,
    3,4,3,3,3,3,2,4,3,5]

# regent.f DATA: XMAX (large/small routing threshold), XMIN (blend lower), DGMIN (DBH-from-H threshold).
const OC_REG_XMAX  = Float32[i == 23 || i == 50 ? 10.0f0 : 4.0f0 for i in 1:50]
const OC_REG_XMIN  = fill(2.0f0, 50)
const OC_REG_DGMIN = Float32[i == 23 || i == 50 ? 7.0f0 : 3.0f0 for i in 1:50]

"""
    oc_smhtgf(sp, d, h, cr_regent, ba, bal, si, relht) -> HTGR

smhtgf.f — small-tree 5-yr height increment (deterministic). `cr_regent` = ICR/10 (0..10 scale, as
regent.f passes it). One of 5 equations selected by MAPSP(sp): 1=pines 2=firs(incl. DF) 3=blackoak
4=tanoak 5=redwood. FACTOR = SPADJF·(0.80+0.004·(SI−50)); TEMBAL = max(BAL,5). Floored at 0.1.
"""
@inline function oc_smhtgf(sp::Int, d::Float32, h::Float32, cr_regent::Float32,
                           ba::Float32, bal::Float32, si::Float32, relht::Float32)::Float32
    factor = OC_SMHT_SPADJF[sp] * (0.80f0 + 0.004f0*(si - 50f0))
    tembal = bal < 5f0 ? 5f0 : bal
    htgr = 0f0
    branch = OC_SMHT_MAPSP[sp]
    if branch == 1              # PINES
        htgr = fexp(0.7452f0 - 0.003271f0*bal - 0.1632f0*cr_regent + 0.0217f0*cr_regent*cr_regent +
                    0.00536f0*si) * factor * 1.75f0
    elseif branch == 2          # FIRS (incl. Douglas-fir)
        domhtgr = 5f0*(2.2227f0 + 0.4314f0*si)/(29.0f0 - 0.05f0*si)
        crf = cr_regent/10f0
        crmod = 1f0 - fexp(-4.26558f0*crf)
        rhmod = fexp(2.54119f0*(fpow(relht, 0.250537f0) - 1f0))
        htgr = domhtgr*(1.016605f0*crmod*rhmod)
    elseif branch == 3          # BLACK OAK
        htgr = fexp(3.817f0 - 0.7829f0*flog(tembal)) * factor
    elseif branch == 4          # TANOAK
        htgr = fexp(3.385f0 - 0.5898f0*flog(tembal)) * factor
    else                        # REDWOOD (5)
        htmax = 2.242202f0*si
        if htmax - h <= 1f0
            htgr = 0f0
        else
            age1 = 1f0/(-0.010742f0)*flog(1f0 - (h/2.242202f0/si)^(1f0/0.919076f0))
            age2 = age1 + 5f0
            h1 = 2.242202f0*si*fpow(1f0 - fexp(-0.010742f0*age1), 0.919076f0)
            h2 = 2.242202f0*si*fpow(1f0 - fexp(-0.010742f0*age2), 0.919076f0)
            htgr = h2 - h1
        end
    end
    return htgr < 0.1f0 ? 0.1f0 : htgr
end

# blkdat.f HT1/HT2 = the BX/AX Wykoff height-diameter coefficients (regent DBH-from-height). grinit.f sets
# IABFLG(sp)=1 for all ⇒ AX=HT1(sp); LHTDRG(sp)=.FALSE. ⇒ this BX/AX path is used (NOT the Curtis/Arney HTDBH).
const OC_HD_HT1 = Float32[
    4.7874,5.2052,4.7874,5.2180,5.2973,5.2973,5.3076,4.7874,4.7874,4.7874,
    4.6843,4.8358,4.7874,4.7874,5.1419,5.3371,5.2649,5.3820,4.7874,4.6236,
    4.7874,4.7874,5.3401,4.7874,4.7874,4.6618,4.6618,4.6618,4.6618,3.8314,
    4.4907,4.6618,4.6618,4.6618,4.6618,4.6618,4.4809,4.6618,4.6618,4.6618,
    4.6618,4.6618,4.6618,4.6618,4.6618,4.6618,4.6618,4.6618,4.6618,5.3401]
const OC_HD_HT2 = Float32[
    -7.3170,-20.1443,-7.3170,-14.8682,-17.2042,-17.2042,-14.4740,-7.3170,-7.3170,-7.3170,
    -6.5516,-9.2077,-7.3170,-7.3170,-19.8143,-19.3151,-15.5907,-20.4097,-7.3170,-13.0049,
    -7.3170,-7.3170,-15.9354,-7.3170,-7.3170,-8.3312,-8.3312,-8.3312,-8.3312,-4.8221,
    -7.7030,-8.3312,-8.3312,-8.3312,-8.3312,-8.3312,-7.5989,-8.3312,-8.3312,-8.3312,
    -8.3312,-8.3312,-8.3312,-8.3312,-8.3312,-8.3312,-8.3312,-8.3312,-8.3312,-15.9354]
# regent.f DATA DIAM — the minimum DBH floor for small/regenerated trees.
const OC_REG_DIAM = Float32[
    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5, 0.5,0.4,0.5,0.5,0.5,0.5,0.3,0.5,0.5,0.5,
    0.3,0.3,0.3,0.3,0.3,0.2,0.2,0.2,0.2,0.2, 0.2,0.2,0.2,0.3,0.1,0.1,0.2,0.1,0.3,0.4,
    0.2,0.2,0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3]

"regent.f BX/AX height-diameter: DBH at height `h` (h>4.5). AX=HT1(sp) (IABFLG=1), BX=HT2(sp)."
@inline oc_hd_diam(sp::Int, h::Float32)::Float32 = OC_HD_HT2[sp]/(flog(h - 4.5f0) - OC_HD_HT1[sp]) - 1.0f0

"""
    oc_regent_dbh(sp, d, h, htg, bark, dg_large_out, scale2) -> new_dbh

regent.f small-tree DBH-from-height for a tree with DBH < DGMIN(sp). HK = H+HTG. If HK≤4.5 the tree stays
sub-breast-height (DBH = D + 0.001·HK). Else DBH grows by DGSM: DK=diam(HK), DKK=diam(H) (=D if H≤4.5),
DGSM=(DK−DKK)·BARK·XRDGRO → DDS=DGSM·(2·BARK·D+DGSM)·SCALE2 → DGSM=√((D·BARK)²+DDS)−D·BARK, then blended
with the large-tree DG via XDWT (D<1.5 ⇒ pure small-tree). The final DBH increment is added DIRECTLY
(update.f `DBH += DG`, measured: live 9987 DBH-INC = DGSM, NOT DGSM/bark). `dg_large_out` = the large-tree
OUTSIDE-bark DBH increment (jl's dg/bark). `scale2` = YR/FINT (=1 for a 5-yr cycle: the √ is then an
identity). Floored at DIAM(sp). DGSD=0 ⇒ deterministic (regent.f ZZRAN=0).
"""
@inline function oc_regent_dbh(sp::Int, d::Float32, h::Float32, htg::Float32,
                               bark::Float32, dg_large_out::Float32, scale2::Float32)::Float32
    hk = h + htg
    hk <= 4.5f0 && return d + 0.001f0*hk
    dk  = oc_hd_diam(sp, hk)
    dkk = h <= 4.5f0 ? d : oc_hd_diam(sp, h)
    dgsm = (dk - dkk)*bark                        # ·XRDGRO(1)
    dgsm < 0f0 && (dgsm = 0f0)
    dds = dgsm*(2.0f0*bark*d + dgsm)*scale2
    dgsm = sqrt((d*bark)*(d*bark) + dds) - bark*d
    dgsm < 0f0 && (dgsm = 0f0)
    xmn = OC_REG_XMIN[sp]
    xdwt = (sp == 23 || sp == 50) ?
           (d <= xmn ? 0f0 : clamp((d - xmn)/(OC_REG_DGMIN[sp] - xmn), 0f0, 1f0)) :
           (d <= 1.5f0 ? 0f0 : (d >= 3.0f0 ? 1.0f0 : (d - 1.5f0)/1.5f0))
    dg = dgsm*(1f0 - xdwt) + dg_large_out*xdwt     # final DBH increment (added directly)
    newd = d + dg
    newd < OC_REG_DIAM[sp] && (newd = OC_REG_DIAM[sp])
    return newd
end

"""
    oc_regent_htg(sp, d, h, crown_pct, pct, ba, avh, hcor, htg_large) -> HTG

regent.f small-tree HEIGHT increment for an inventory (non-LESTB) tree, DGSD=0 (ZZRAN=0, XRHGRO=SCALE=1).
`crown_pct` = crown ratio % (0..100), `pct` = BA percentile (jl stores it in t.crown_ratio), `hcor` =
small-tree height calibration ln(CORNEW) (0 unless calibrated). `htg_large` = the large-tree oc_htgf_native
estimate, blended in for 2 < D < XMAX via XWT. Returns the final HTG.
"""
@inline function oc_regent_htg(sp::Int, d::Float32, h::Float32, crown_pct::Integer, pct::Float32,
                               ba::Float32, avh::Float32, si::Float32, hcor::Float32,
                               htg_large::Float32)::Float32
    cr_regent = Float32(crown_pct)/10f0                 # regent CR = ICR/10
    bal = (1f0 - pct/100f0)*ba                          # stand BAL (jl pct = t.crown_ratio)
    relht = avh > 0f0 ? h/avh : 1f0
    relht > 1.05f0 && (relht = 1.05f0)
    htgrr = oc_smhtgf(sp, d, h, cr_regent, ba, bal, si, relht)
    con = fexp(hcor)                                    # RHCON=1 (no small-tree-DG calib) · exp(HCOR)
    htgr = htgrr*con                                    # ·XRHGRO(1)·SCALE(1); ZZRAN=0 for DGSD=0
    xmn = OC_REG_XMIN[sp]; xmx = OC_REG_XMAX[sp]
    xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
    htg = if sp == 23 || sp == 50
        lthg = htg_large
        blended = (htgr + lthg)/2f0                     # regent.f:255 GS/RW pre-average (non-LESTB)
        blended*(1f0 - xwt) + xwt*lthg
    else
        htgr*(1f0 - xwt) + xwt*htg_large
    end
    return htg < 0.1f0 ? 0.1f0 : htg
end
