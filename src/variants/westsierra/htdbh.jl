# =============================================================================
# htdbh.jl (westsierra) — WS Curtis-Arney height↔diameter (ws/htdbh.f). Chunk 4b.
#
# MODE=1 (H→D) inverse Curtis-Arney, used by regent (ws/regent.f:530-534) for the CA-surrogate regen
# species' DK/DKK breast-height crossing — and it fires even with LHTDRG=false (regent condition
# `.NOT.LHTDRG .OR. (LHTDRG.AND.IABFLG=1)` is TRUE when LHTDRG is false). CURARN(43,3)=P2/P3/P4, SPLINE(43)=Z.
# Nonzero only for the CA-surrogate/GS/RW/MC species (the WS-native main species use the linear DK=AX+BX·HK
# in regent, not this). MODE 0 (D→H, ws_htcalc-style) is a follow-on (not needed by regent's H→D path).
# =============================================================================

# ws/htdbh.f DATA CURARN(43,3) — Curtis-Arney P2 (col1), P3 (col2), P4 (col3); index = ISPC.
const WS_HTDBH_P2 = Float32[
  0,0,0,595.1068,0, 0,0,0,99.1568,99.1568, 0,101.5170,0,101.5170,101.5170,
  101.5170,101.5170,0,101.5170,101.5170, 0,0,595.1068,0,101.5170,
  101.5170,101.5170,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 1709.7229,0,0]
const WS_HTDBH_P3 = Float32[
  0,0,0,5.8103,0, 0,0,0,12.1300,12.1300, 0,4.7066,0,4.7066,4.7066,
  4.7066,4.7066,0,4.7066,4.7066, 0,0,5.8103,0,4.7066,
  4.7066,4.7066,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 5.8887,0,0]
const WS_HTDBH_P4 = Float32[
  0,0,0,-0.3821,0, 0,0,0,-1.3272,-1.3272, 0,-0.9540,0,-0.9540,-0.9540,
  -0.9540,-0.9540,0,-0.9540,-0.9540, 0,0,-0.3821,0,-0.9540,
  -0.9540,-0.9540,0,0,0, 0,0,0,0,0, 0,0,0,0,0, -0.2286,0,0]
const WS_HTDBH_SPLINE = Float32[
  0,0,0,3,0, 0,0,0,5,5, 0,2,0,2,2, 2,2,0,2,2, 0,0,3,0,2,
  2,2,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 3,0,0]

# ws/htdbh.f MODE=0 (D→H). Returns predicted total height (ft) for DBH d. Curtis-Arney with the SPLINE(Z)
# small-tree linear tie to (0.3, 4.51). Used by cratet.f dubbing for the LHTDRG-false surrogate species
# (their CURARN P2/P3/P4 are nonzero); native species dub via the Wykoff HT-DBH (ht1/ht2), not this.
@inline function ws_htdbh_height(ifor::Int, sp::Int, d::Float32)::Float32
    p2 = WS_HTDBH_P2[sp]; p3 = WS_HTDBH_P3[sp]; p4 = WS_HTDBH_P4[sp]; z = WS_HTDBH_SPLINE[sp]
    if d >= z
        return 4.5f0 + p2 * exp(-1f0 * p3 * fpow(d, p4))
    else
        return ((4.5f0 + p2 * exp(-1f0 * p3 * fpow(z, p4)) - 4.51f0) * (d - 0.3f0) / (z - 0.3f0)) + 4.51f0
    end
end

# ws/htdbh.f MODE=1 (H→D). Returns predicted DBH for total height h (ft). ifor unused (matches Fortran).
@inline function ws_htdbh_dbh(ifor::Int, sp::Int, h::Float32)::Float32
    p2 = WS_HTDBH_P2[sp]; p3 = WS_HTDBH_P3[sp]; p4 = WS_HTDBH_P4[sp]; z = WS_HTDBH_SPLINE[sp]
    hatz = 4.5f0 + p2 * exp(-1f0 * p3 * fpow(z, p4))
    if h >= hatz
        return exp(log((log(h - 4.5f0) - log(p2)) / (-1f0 * p3)) * (1f0 / p4))
    else
        return (((h - 4.51f0) * (z - 0.3f0)) / (4.5f0 + p2 * exp(-1f0 * p3 * fpow(z, p4)) - 4.51f0)) + 0.3f0
    end
end
