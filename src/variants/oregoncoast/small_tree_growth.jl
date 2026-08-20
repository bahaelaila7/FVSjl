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
