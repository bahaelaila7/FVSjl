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

# regent.f DATA DIAM — the minimum DBH floor for small/regenerated trees.
const OC_REG_DIAM = Float32[
    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5, 0.5,0.4,0.5,0.5,0.5,0.5,0.3,0.5,0.5,0.5,
    0.3,0.3,0.3,0.3,0.3,0.2,0.2,0.2,0.2,0.2, 0.2,0.2,0.2,0.3,0.1,0.1,0.2,0.1,0.3,0.4,
    0.2,0.2,0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3]

# htdbh.f CURARN (Curtis/Arney P2,P3,P4) + SPLINE(Z), the DEFAULT height-diameter model regent uses when
# LHTDRG(sp)=.FALSE. (the default — no HTGROW keyword in oct01). MEASURED: HTDBH(DF,HK=15.63)=2.341 = live
# regent dump DK 2.34; HTDBH(GF,HK=5.14) → DGSM≈0.291 = dump 0.2919. (The BX/AX+AA calibration path applies
# only under HTGROW; not exercised here.)
const OC_CURARN1 = Float32[
    8532.9026,695.4196,487.5415,467.3070,606.3002,606.3002,408.7614,263.1274,233.6987,89.5535,
    101.5170,99.1568,514.1013,514.1013,744.7718,944.9299,422.0948,1267.7589,113.7962,79986.6348,
    60.6009,91.7438,595.1068,127.1698,79986.6348,105.0771,105.0771,59.0941,59.0941,40.3812,
    120.2372,126.7237,55.0,143.9994,55.0,94.5048,117.7410,1176.9704,403.3221,97.7769,
    105.0771,679.1972,55.0,47.3648,179.0706,149.5861,55.0,114.1627,40.3812,595.1068]
const OC_CURARN2 = Float32[
    8.0343,7.5021,5.4444,6.1195,6.2936,6.2936,5.4044,6.9356,6.9059,4.2281,
    4.7066,12.1300,5.5983,5.5983,7.6793,6.2428,6.0404,7.4995,4.7726,9.9284,
    4.1543,17.1081,5.8103,4.8977,9.9284,5.6647,5.6647,6.1195,6.1195,3.7653,
    4.1713,3.1800,5.5,3.5124,5.5,4.0657,4.0764,6.3245,4.3271,8.8202,
    5.6647,5.5698,5.5,15.6276,3.6238,2.4231,5.5,6.0210,3.7653,5.8103]
const OC_CURARN3 = Float32[
    -0.1831,-0.3852,-0.3801,-0.4325,-0.3860,-0.3860,-0.4426,-0.6619,-0.6166,-0.6438,
    -0.9540,-1.3272,-0.2734,-0.2734,-0.3779,-0.3087,-0.4525,-0.3286,-0.7601,-0.1013,
    -0.6277,-1.4429,-0.3821,-0.4668,-0.1013,-0.6822,-0.6822,-1.0552,-1.0552,-1.1224,
    -0.6113,-0.6324,-0.95,-0.5511,-0.95,-0.9592,-0.6151,-0.2739,-0.2422,-1.0534,
    -0.6822,-0.3074,-0.95,-1.9266,-0.5730,-0.1800,-0.95,-0.7838,-1.1224,-0.3821]
const OC_HD_SPLINE = Float32[
    3,6,3,3,3,3,3,3,3,3, 2,5,3,3,5,5,3,2,3,2,
    3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3]

"htdbh.f MODE=1 (Curtis/Arney): DBH at height `h` (h>4.5). The default H-D model (LHTDRG=.FALSE.)."
@inline function oc_hd_diam(sp::Int, h::Float32)::Float32
    p2 = OC_CURARN1[sp]; p3 = OC_CURARN2[sp]; p4 = OC_CURARN3[sp]; z = OC_HD_SPLINE[sp]
    hatz = 4.5f0 + p2*fexp(-p3*fpow(z, p4))
    return h >= hatz ? fexp(flog((flog(h - 4.5f0) - flog(p2))/(-p3)) / p4) :
                       ((h - 4.51f0)*(z - 0.3f0))/(hatz - 4.51f0) + 0.3f0
end
# (oc_htdbh_height MODE=0 — height from DBH — is defined in organon_cratet.jl and reused by the FFE
#  crown-biomass small-unmerch floor; the crown_biomass.jl OC branch calls it.)

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
