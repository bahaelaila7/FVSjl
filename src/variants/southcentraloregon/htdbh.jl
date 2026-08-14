# =============================================================================
# htdbh.jl (southcentraloregon) — SO height↔DBH (Curtis) missing-height dub. Chunk 4a.
# (so/htdbh.f)
#
# so/cratet.f LHTDRG=.FALSE. for all species ⇒ ALL missing inventory heights are dubbed by the
# FOREST-DEPENDENT Curtis H-D curve (so/htdbh.f MODE=0), like WC/PN/EC. Three physical forest
# coefficient tables (so/htdbh.f DATA): DESCHT (Deschutes, IFOR 1/10 + Warm Springs), FREMNT (Fremont,
# IFOR 2), WINEMA (Winema, default incl. R5 forests). Each is 33 species × 3 params (P2/P3/P4). The
# spline knot is a fixed Z=3.0" (unlike CA's per-species SPLINE). MODE=0 (given D, predict H):
#   D ≥ 3:  H = 4.5 + P2·exp(−P3·D^P4)
#   D < 3:  H = ((4.5 + P2·exp(−P3·3^P4) − 4.51)·(D−0.3)/2.7) + 4.51   (linear taper to the 0.3" seedling)
# MODE=1 (given H, predict D) is the exact inverse (used by regent when LHTDRG — inert for SO).
# =============================================================================

# so/htdbh.f DATA DESCHT/FREMNT/WINEMA (MAXSP,3) → [sp] P2/P3/P4, 33 species in SO order.
const SO_HTDBH_DESCHT_P2 = Float32[582.9947,128.8697,253.2541,235.3340,197.4948,4902.9732,777.9043,290.2790,606.3002,2700.1370,0.0,235.3340,291.4070,171.2219,247.7348,0.0,255.4638,616.3503,317.8257,77.2207,133.7965,484.4591,76.5170,0.0,178.6441,73.3348,40.3812,149.5861,10707.3906,1709.7229,1709.7229,253.2541,1709.7229]
const SO_HTDBH_DESCHT_P3 = Float32[5.4612,6.8868,4.7331,5.7931,6.7218,7.5484,5.2036,6.5834,6.2936,7.1184,0.0,5.7931,5.7543,9.9497,6.1830,0.0,5.5577,5.7620,6.8287,3.5181,6.4050,4.5713,2.2107,0.0,4.5852,2.6548,3.7653,2.4231,8.4670,5.8887,5.8887,4.7331,5.8887]
const SO_HTDBH_DESCHT_P4 = Float32[-0.3435,-0.8701,-0.4843,-0.5972,-0.6528,-0.1783,-0.2843,-0.5753,-0.3860,-0.2312,0.0,-0.5972,-0.5013,-0.9727,-0.6335,0.0,-0.6054,-0.3633,-0.6034,-0.5894,-0.8329,-0.3643,-0.6365,0.0,-0.6746,-1.2460,-1.1224,-0.1800,-0.1863,-0.2286,-0.2286,-0.4843,-0.2286]

const SO_HTDBH_FREMNT_P2 = Float32[391.7346,948.8488,253.2541,705.1903,103.7798,76621.9189,113.7962,290.2790,606.3002,1154.2447,0.0,705.1903,1056.5434,171.2219,247.7348,0.0,255.4638,616.3503,317.8257,77.2207,133.7965,484.4591,76.5170,0.0,178.6441,73.3348,40.3812,149.5861,10707.3906,1709.7229,1709.7229,253.2541,1709.7229]
const SO_HTDBH_FREMNT_P3 = Float32[5.1102,6.1388,4.7331,6.0971,10.6932,10.2682,4.7726,6.5834,6.2936,6.6836,0.0,6.0971,6.6974,9.9497,6.1830,0.0,5.5577,5.7620,6.8287,3.5181,6.4050,4.5713,2.2107,0.0,4.5852,2.6548,3.7653,2.4231,8.4670,5.8887,5.8887,4.7331,5.8887]
const SO_HTDBH_FREMNT_P4 = Float32[-0.3602,-0.2877,-0.4843,-0.3273,-1.0711,-0.1178,-0.7601,-0.5753,-0.3860,-0.2876,0.0,-0.3273,-0.2950,-0.9727,-0.6335,0.0,-0.6054,-0.3633,-0.6034,-0.5894,-0.8329,-0.3643,-0.6365,0.0,-0.6746,-1.2460,-1.1224,-0.1800,-0.1863,-0.2286,-0.2286,-0.4843,-0.2286]

const SO_HTDBH_WINEMA_P2 = Float32[133.7789,222.7080,231.7163,471.6016,130.7104,4518.2601,128.7972,168.9700,606.3002,812.2630,0.0,471.6016,299.3002,171.2219,247.7348,0.0,255.4638,616.3503,317.8257,77.2207,133.7965,484.4591,76.5170,0.0,178.6441,73.3348,40.3812,149.5861,10707.3906,1709.7229,1709.7229,231.7163,1709.7229]
const SO_HTDBH_WINEMA_P3 = Float32[6.9968,6.1735,6.7143,5.7106,7.7823,8.0469,4.9833,13.6848,6.2936,6.4422,0.0,5.7106,6.3401,9.9497,6.1830,0.0,5.5577,5.7620,6.8287,3.5181,6.4050,4.5713,2.2107,0.0,4.5852,2.6548,3.7653,2.4231,8.4670,5.8887,5.8887,6.7143,5.8887]
const SO_HTDBH_WINEMA_P4 = Float32[-0.9072,-0.6122,-0.6647,-0.4035,-0.8830,-0.2090,-0.7463,-1.0635,-0.3860,-0.3348,0.0,-0.4035,-0.5275,-0.9727,-0.6335,0.0,-0.6054,-0.3633,-0.6034,-0.5894,-0.8329,-0.3643,-0.6365,0.0,-0.6746,-1.2460,-1.1224,-0.1800,-0.1863,-0.2286,-0.2286,-0.6647,-0.2286]

# so/htdbh.f SELECT CASE(IFOR): 1/10 → Deschutes, 2 → Fremont, DEFAULT → Winema.
@inline function _so_htdbh_coefs(ifor::Int, sp::Int)
    if ifor == 1 || ifor == 10
        return (SO_HTDBH_DESCHT_P2[sp], SO_HTDBH_DESCHT_P3[sp], SO_HTDBH_DESCHT_P4[sp])
    elseif ifor == 2
        return (SO_HTDBH_FREMNT_P2[sp], SO_HTDBH_FREMNT_P3[sp], SO_HTDBH_FREMNT_P4[sp])
    else
        return (SO_HTDBH_WINEMA_P2[sp], SO_HTDBH_WINEMA_P3[sp], SO_HTDBH_WINEMA_P4[sp])
    end
end

"so/htdbh.f MODE=0: predicted total height (ft) for species `sp` at DBH `d`, forest `ifor`."
@inline function so_htdbh_height(ifor::Int, sp::Int, d::Float32)
    p2, p3, p4 = _so_htdbh_coefs(ifor, sp)
    if d >= 3f0
        return 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(d, p4))
    else
        return ((4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3f0, p4)) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
    end
end

"so/htdbh.f MODE=1: predicted DBH for species `sp` at total height `h` (inverse of so_htdbh_height)."
@inline function so_htdbh_dbh(ifor::Int, sp::Int, h::Float32)
    p2, p3, p4 = _so_htdbh_coefs(ifor, sp)
    hat3 = 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3f0, p4))
    if h >= hat3
        return fexp(log((log(h - 4.5f0) - log(p2)) / (-1f0 * p3)) * (1f0 / p4))
    else
        return (((h - 4.51f0) * 2.7f0) / (4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3f0, p4)) - 4.51f0)) + 0.3f0
    end
end
