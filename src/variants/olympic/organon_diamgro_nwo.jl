# =============================================================================
# organon_diamgro_nwo.jl — OP (Olympic) ORGANON NWO diameter-growth core (validated chunk).
#
# Ported from organon/diagro.f DG_NWO_RUN (VERSION=2). This is the ONE piece of new coefficient data
# OP needs on top of the shared ORGANON engine: the NWO diameter-growth equation. The ln(DDS) form,
# CRADJ crown adjustment, and DIAMGRO_RUN driver are IDENTICAL to SWO (organon/diagro.f DG_SWO);
# only the 11-group × 11-param DGPAR table and the species-group ADJ multiplier differ.
#
#   LNDG = B0 + B1·ln(DBH+K1) + B2·DBH^K2 + B3·ln((CR+0.2)/1.2) + B4·ln(SITE)
#            + B5·(SBAL1^K3 / ln(DBH+K4)) + B6·√SBA1
#   CRADJ = 1 ; if CR≤0.17: CRADJ = 1 − exp(−(25·CR)²)
#   DG = exp(LNDG)·CRADJ·ADJ(group)
#   (DIAMGRO_RUN then multiplies CALIB(3,g)·FERTADJ·THINADJ — all 1.0 on cyc0 no-treatment.)
#
# SITE is ORGANON SI_1 (or SI_2 for WH group 3) minus 4.5. SBAL1 = basal-area-in-larger (GET_BAL),
# SBA1 = stand basal area — both from the shared version-agnostic SSTATS/GET_BAL (organon/statsorg.f,
# organon/diamcal.f), reused from the ORGANON engine in the full-integration chunk.
#
# VALIDATED bit-exact vs the live FVSop_dbg oracle: an instrumented DG_NWO_RUN (WRITE of
# ISPGRP,DBH,CR,SITE,SBAL1,SBA1,DG dumped as Z8 hex for lossless Float32 round-trip; stand S248112,
# 73 tree-calls) reproduced HERE gives 72/73 BIT-EXACT and 1/73 at 1 ULP (the documented irreducible
# libm exp/log intrinsic ULP). Species-group indexing via op_spgroup_nwo (species.jl). buildDir left
# pristine (marker count 0). See docs/OP_VARIANT_PORT_AUDIT.md.
# =============================================================================

# organon/diagro.f DG_NWO_RUN DGPAR(11,11) — groups 1..11 = DF,GF,WH,RC,PY,MD,BL,WO,RA,PD,WI.
const OP_DG_NWO_B0 = Float32[-4.69624,-2.34619,-4.49867,-11.45456097,-9.15835863,-8.84531757,-3.41449922,-7.81267986,-4.39082007,-8.08352683,-8.08352683]
const OP_DG_NWO_B1 = Float32[0.339513,0.594640,0.362369,0.784133664,1.0,1.5,1.0,1.405616529,1.0,1.0,1.0]
const OP_DG_NWO_B2 = Float32[-0.000428261,-0.000976092,-0.00153907,-0.0261377888,-0.00000035,-0.0006,-0.05,-0.0603105850,-0.0945057147,-0.00000035,-0.00000035]
const OP_DG_NWO_B3 = Float32[1.19952,1.12712,1.1557,0.70174783,1.16688474,0.51225596,0.0,0.64286007,1.06867026,0.31176647,0.31176647]
const OP_DG_NWO_B4 = Float32[1.15612,0.555333,1.12154,2.057236260,0.0,0.418129153,0.324349277,1.037687142,0.685908029,0.0,0.0]
const OP_DG_NWO_B5 = Float32[-0.0000446327,-0.0000290672,-0.0000201041,-0.00415440257,0.0,-0.00355254593,0.0,0.0,-0.00586331028,0.0,0.0]
const OP_DG_NWO_B6 = Float32[-0.0237003,-0.0470848,-0.0417388,0.0,-0.02,-0.0321315389,-0.0989519477,-0.0787012218,0.0,-0.0730788052,-0.0730788052]
const OP_DG_NWO_K1 = Float32[1.0,1.0,1.0,5.0,4000.0,110.0,10.0,5.0,5.0,4000.0,4000.0]
const OP_DG_NWO_K2 = Float32[2.0,2.0,2.0,1.0,4.0,2.0,1.0,1.0,1.0,4.0,4.0]
const OP_DG_NWO_K3 = Float32[2.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
const OP_DG_NWO_K4 = Float32[5.0,5.0,5.0,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7]

"organon/diagro.f DG_NWO_RUN ADJ — species-group full adjustment multiplier (VERSION=2)."
@inline function op_dg_nwo_adj(g::Int)
    g == 1 && return 0.7011014f0
    g == 2 && return 0.8722f0
    g == 3 && return 0.7163f0
    g == 6 && return 0.7928f0
    g == 8 && return 1.0f0
    return 0.8f0
end

"""
    op_dg_nwo(g, dbh, cr, site, sbal1, sba1) -> DG

organon/diagro.f DG_NWO_RUN — 5-yr ORGANON diameter growth for NWO species group `g` (1..11).
`site` = ORGANON SI (SITE_1−4.5, or SITE_2−4.5 for WH group 3). REAL*4 throughout (fexp/flog/fpow).
"""
@inline function op_dg_nwo(g::Int, dbh::Float32, cr::Float32, site::Float32,
                           sbal1::Float32, sba1::Float32)
    b0 = OP_DG_NWO_B0[g]; b1 = OP_DG_NWO_B1[g]; b2 = OP_DG_NWO_B2[g]; b3 = OP_DG_NWO_B3[g]
    b4 = OP_DG_NWO_B4[g]; b5 = OP_DG_NWO_B5[g]; b6 = OP_DG_NWO_B6[g]
    k1 = OP_DG_NWO_K1[g]; k2 = OP_DG_NWO_K2[g]; k3 = OP_DG_NWO_K3[g]; k4 = OP_DG_NWO_K4[g]
    lndg = b0 +
           b1*flog(dbh + k1) +
           b2*fpow(dbh, k2) +
           b3*flog((cr + 0.2f0)/1.2f0) +
           b4*flog(site) +
           b5*(fpow(sbal1, k3)/flog(dbh + k4)) +
           b6*sqrt(sba1)
    cradj = 1.0f0
    if cr <= 0.17f0
        cradj = 1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))
    end
    return fexp(lndg)*cradj*op_dg_nwo_adj(g)
end

"""
    op_diamgro_run_nwo(sp, dbh, cr, site, sbal1, sba1; calib3=1f0, fertadj=1f0, thinadj=1f0) -> DGRO

organon/diagro.f DIAMGRO_RUN (VERSION=2) — DGRO = DG_NWO·CALIB(3,g)·FERTADJ·THINADJ. `sp` is the FVS
species sequence number (1..39); the NWO species group is resolved via `op_spgroup_nwo`. On cyc0
inventory with no fertilizer/thinning, CALIB(3,g)=FERTADJ=THINADJ=1.0.
"""
@inline function op_diamgro_run_nwo(sp::Integer, dbh::Float32, cr::Float32, site::Float32,
        sbal1::Float32, sba1::Float32; calib3::Float32=1f0, fertadj::Float32=1f0, thinadj::Float32=1f0)
    g = op_spgroup_nwo(sp)
    g == 0 && error("OP DG_NWO: FVS species $sp maps to no NWO ORGANON group (org FIA $(op_organon_fia(sp)))")
    return op_dg_nwo(g, dbh, cr, site, sbal1, sba1)*calib3*fertadj*thinadj
end
