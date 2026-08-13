# =============================================================================
# organon_cratet.jl — OC (Oregon Coast) setup height/crown dubbing (chunk C8).
#
# Ported from:
#   • bin/FVSoc_buildDir/htdbh.f — HTDBH, the CA-family Curtis-Arney height–diameter INVENTORY
#     equation (CURARN(50,3) + SPLINE(50)); MODE=0 (dbh→height).
#   • oc/cratet.f:650-690 — the missing-height dub. OC has LHTDRG=.FALSE. for every species
#     (oc/grinit.f), so the calibrated-Wykoff branch never fires and the HTDBH INVENTORY equation
#     is the actual dub for a missing-height NON-ORGANON tree (it overwrites the Wykoff H at
#     cratet.f:681). D≤0.1 keeps H=1.01 (handled by the shared `dub_missing_heights!`).
#
# The ORGANON trees' missing HT/CR are instead imputed by ORGANON PREPARE (chunk C2,
# `organon_prepare_swo` PRDHT/PRDCR) — NOT exercised by ocmin (no valid ORGANON tree has a missing
# HT/CR there; `KNTOHT=KNTOCR=0`), so this file ports the load-bearing NON-ORGANON HTDBH dub. Wiring
# PREPARE's dub for a stand that DOES carry an ORGANON missing height is a follow-up.
#
# MEASURED vs the live FVSoc_clean oracle (ocmin CRATET dump): tree-20 LP (sp 12, D=8.5) →
# H=53.3234 via `INVENTORY EQN DUBBING ISPC=12`. See docs/OC_VARIANT_PORT_AUDIT.md (C8).
# =============================================================================

# htdbh.f CURARN(50,3): P2,P3,P4 (Curtis-Arney), species order = OC/CA index 1..50.
const OC_CURARN_P2 = Float32[8532.9026,695.4196,487.5415,467.3070,606.3002,606.3002,408.7614,263.1274,233.6987,89.5535,101.5170,99.1568,514.1013,514.1013,744.7718,944.9299,422.0948,1267.7589,113.7962,79986.6348,60.6009,91.7438,595.1068,127.1698,79986.6348,105.0771,105.0771,59.0941,59.0941,40.3812,120.2372,126.7237,55.0,143.9994,55.0,94.5048,117.7410,1176.9704,403.3221,97.7769,105.0771,679.1972,55.0,47.3648,179.0706,149.5861,55.0,114.1627,40.3812,595.1068]
const OC_CURARN_P3 = Float32[8.0343,7.5021,5.4444,6.1195,6.2936,6.2936,5.4044,6.9356,6.9059,4.2281,4.7066,12.1300,5.5983,5.5983,7.6793,6.2428,6.0404,7.4995,4.7726,9.9284,4.1543,17.1081,5.8103,4.8977,9.9284,5.6647,5.6647,6.1195,6.1195,3.7653,4.1713,3.1800,5.5,3.5124,5.5,4.0657,4.0764,6.3245,4.3271,8.8202,5.6647,5.5698,5.5,15.6276,3.6238,2.4231,5.5,6.0210,3.7653,5.8103]
const OC_CURARN_P4 = Float32[-0.1831,-0.3852,-0.3801,-0.4325,-0.3860,-0.3860,-0.4426,-0.6619,-0.6166,-0.6438,-0.9540,-1.3272,-0.2734,-0.2734,-0.3779,-0.3087,-0.4525,-0.3286,-0.7601,-0.1013,-0.6277,-1.4429,-0.3821,-0.4668,-0.1013,-0.6822,-0.6822,-1.0552,-1.0552,-1.1224,-0.6113,-0.6324,-0.95,-0.5511,-0.95,-0.9592,-0.6151,-0.2739,-0.2422,-1.0534,-0.6822,-0.3074,-0.95,-1.9266,-0.5730,-0.1800,-0.95,-0.7838,-1.1224,-0.3821]
const OC_CURARN_Z  = Float32[3.,6.,3.,3.,3.,3.,3.,3.,3.,3.,2.,5.,3.,3.,5.,5.,3.,2.,3.,2.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.]

"""
    oc_htdbh_height(sp, d) -> H

`bin/FVSoc_buildDir/htdbh.f` HTDBH MODE=0 — the CA-family Curtis-Arney height from DBH `d` for OC
species index `sp`. For `d ≥ Z`: `H = 4.5 + P2·exp(−P3·D^P4)`; for `d < Z`: the linear spline down to
4.51 at D=0.3. (IFOR is unused in this variant.)
"""
@inline function oc_htdbh_height(sp::Int, d::Float32)
    p2 = OC_CURARN_P2[sp]; p3 = OC_CURARN_P3[sp]; p4 = OC_CURARN_P4[sp]; z = OC_CURARN_Z[sp]
    if d >= z
        return 4.5f0 + p2*fexp(-1.0f0*p3*fpow(d, p4))
    else
        return ((4.5f0 + p2*fexp(-1.0f0*p3*fpow(z, p4)) - 4.51f0)*(d - 0.3f0)/(z - 0.3f0)) + 4.51f0
    end
end
