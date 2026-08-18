# =============================================================================
# ontario/volume.jl — ON (Ontario) per-tree volume (canada/on/vols.f + varvol.f + volont.f).
#
# ON's default volume method is METHC=METHB=8 (the FVS-ON default, vols.f) which routes through
# the Ontario ZAK/HONER special equations in volont.f — NOT R9CLARK/NATCRS (that is method 6).
# For each tree:
#   • VN  = GTV = total merch cubic ft  (pulpwood + sawlog), ground-to-tip.
#   • VM  = GMV = merch sawlog cubic ft  (stump→ZAKHT-merch-height / HONER merch fraction).
#   • BFV = NMV = net merch volume: Mowraski cull applied to GMV (age-based), reported in the
#           board-foot column of the metric .sum.
# ZAK species (varvol.f ONEQ) use the Zakrzewski taper-integral (ZAKVOL) + ZAKHT merch height;
# HON species use the Honer total/merch equations (HONER). The Mowraski cull needs the tree age
# (ABIRTH), which ON dubs from initial height via HTCALC MODE0=0 (Carmean NC-128 SI curves,
# canada/on/htcalc.f + findag.f). We compute-and-freeze that age on the first (cyc0) call, exactly
# as cratet.f seeds ABIRTH once before projection.
#
# Metric constants come from common/METRIC.F77. Transcendentals go through glibc logf/expf/powf
# (on_logf/on_expf/on_powf, defined in diameter_growth.jl) so the Float32 taper integrals and the
# ZAKHT cube-root chain match FVSon_g16 bit-for-bit. SQRT uses Julia's IEEE sqrt (== glibc sqrtf).
#
# VALIDATED: per-tree GTV/GMV/NMV dump-replay BIT-EXACT vs an instrumented canada/on/vols.f
# (single-.o swap into FVSon_g16, /workspace/.onwork/FVSon_voldump) on ont01.
# =============================================================================

const ON_M3toFT3 = 35.314455f0   # common/METRIC.F77  (m³ → ft³)
const ON_CMtoFT  = 0.0328084f0   # common/METRIC.F77  (cm → ft)

# gfortran's `X**4` (__builtin_powi) evaluates as (X²)², which rounds 1 ULP differently from
# Julia's literal `x^4` (power_by_squaring). The Zak taper integrals are catastrophically
# cancelling, so that 1 ULP propagates — use (x²)² everywhere the Fortran writes `**4`.
@inline _on_p4(x::Float32) = (t = x*x; t*t)

# ── species → equation-family maps (volont.f / varvol.f DATA) ────────────────
# ISPMPZ: FVS species → Zak equation (1..12).  ISPMPH: → Honer (1..26).  ISPMPM: → Mowraski (1..14).
const ON_ISPMPZ = Int8[
  3, 3, 2, 2, 1, 5, 4, 6, 4,12, 7,12,12,12,11,11,11, 8, 8,11,
 11,11,11,10,11, 8, 8,11,11,11,11,11,11,11,11,11,11,11,11, 9,
  9, 9,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
 11,11,11,11,11,11,11,11, 3, 1, 5, 4]
const ON_ISPMPH = Int8[
  3, 4, 2, 2, 1, 5, 7, 9, 6,11,10, 8, 1,10,18,18,26,13,13,23,
 20,22,22,14,21,12,12,17,19,16,16,16,16,16,16,16,19,19,19,24,
 26,25,15,12,19,12,22,12,12,12,12,12,22,12,19,22,12,22,12,16,
 23,23,23,12,12,12,12,22, 3, 1, 5, 6]
const ON_ISPMPM = Int8[
  3, 3, 2, 2, 1, 5, 5, 6, 4, 8, 7, 9, 4, 7,14,14,10,13,14,13,
  4,14,14,13,13,12,13,13,13,13,14,13,13,13,14,13,13,13,13,10,
 10,10,11,13,13,13,13,13,14,14,13,13,13,13,13,13,13,13,13,13,
 13,13,13,14,14,14,13,13, 3, 1, 5, 4]
# ONEQ: which family each species uses ('ZAK' vs 'HON').  true = ZAK.
const ON_IS_ZAK = Bool[
 true,false,true,true,true,true,true,true,true,true,
 true,false,true,true,true,true,true,true,true,false,
 false,false,false,true,false,true,true,false,false,false,
 false,false,false,false,false,false,false,false,false,true,
 true,true,true,true,false,true,false,true,true,true,
 true,true,false,true,false,false,true,false,true,false,
 false,false,false,true,true,true,true,false,true,true,
 true,true]

# ── Honer coefficients (volont.f, indexed by ISPMPH 1..26) ───────────────────
const ON_HB1 = Float32[
 0.184,0.151,0.151,0.151,0.176,0.164,0.169,0.155,0.152,0.155,
 0.164,0.145,0.145,0.181,0.176,0.145,0.145,0.145,0.145,0.145,
 0.145,0.145,0.145,0.127,0.127,0.127]
const ON_HC1 = Float32[
 0.691,0.71,0.897,0.897,1.44,1.588,1.226,1.112,2.139,4.167,
 1.588,1.046,1.046,1.449,2.222,1.512,0.959,1.046,1.046,0.634,
 0.948,1.877,0.033,-0.312,0.42,-0.312]
const ON_HC2 = Float32[
 363.676,355.623,348.53,348.53,342.175,333.364,315.832,350.092,301.634,244.906,
 333.364,383.972,383.972,344.754,300.373,336.509,334.829,383.972,383.972,440.496,
 401.456,332.585,393.336,436.683,394.644,436.683]
const ON_HR1 = Float32[
 0.9735,0.9672,0.9635,0.9672,0.9611,0.9526,0.9645,0.9645,0.9352,0.9645,
 0.9526,0.9057,0.9057,0.8778,0.9087,0.9057,0.9057,0.9057,0.9057,0.9057,
 0.9057,0.9057,0.9057,0.9087,0.9354,0.9354]
const ON_HR2 = Float32[
 -0.2348,-0.0393,-0.15,-0.0393,-0.2456,-0.1027,-0.1616,-0.1616,-0.0395,-0.1616,
 -0.1027,-0.0708,-0.0708,-0.2417,-0.3049,-0.0708,-0.0708,-0.0708,-0.0708,-0.0708,
 -0.0708,-0.0708,-0.0708,-0.3049,0.0957,0.0957]
const ON_HR3 = Float32[
 -0.7378,-1.0523,-0.8081,-1.0523,-0.6801,-0.8199,-0.7945,-0.7945,-0.8147,-0.7945,
 -0.8199,-0.8375,-0.8375,-0.5247,-0.5107,-0.8375,-0.8375,-0.8375,-0.8375,-0.8375,
 -0.8375,-0.8375,-0.8375,-0.5107,-1.1613,-1.1613]

# ── Zak coefficients (volont.f, indexed by ISPMPZ 1..12) ─────────────────────
const ON_ZB2 = Float32[
 -1.80769,-1.75095,-1.72170,-1.62979,-1.61060,-1.34130,-1.57624,-1.37784,-1.49699,-1.24598,
 -1.31327,-1.89146]
const ON_ZG2 = Float32[
 0.90926,0.87184,0.84388,0.77509,0.82290,0.50717,0.74247,0.55538,0.63057,0.43204,
 0.45938,1.00000]

# ── Mowraski cull coefficients (volont.f, indexed by ISPMPM 1..14) ───────────
const ON_MAM = Float32[
 0.004286,0.000101,0.012639,0.003614,0.000364,0.003204,0.00894,0.003614,0.000825,0.005207,
 0.00329,0.000324,0.003759,0.015715]
const ON_MBM = Float32[
 2.5398,0.9831,8.3752,3.132,1.1654,1.5584,4.5996,3.132,1.0931,1.4052,
 1.5971,0.4106,1.34313,4.06754]

# ── Age dubbing: LS species map + Carmean NC-128 SI coefficients (htcalc.f) ──
const ON_MAPLS = Int16[
  74,108, 95, 95,104, 68, 64, 55, 70, 59,126,127, 58, 58, 14, 15, 28,  4,  1, 35,
  53, 53, 53,  5, 51,  2,  2, 11, 12, 41, 44, 36, 36, 38, 49, 36, 10, 10, 10, 32,
  32, 25,  8, 16, 16, 16, 58, 50, 58, 58, 58, 58, 58, 53, 19, 58, 58, 58, 27, 25,
  58, 58, 58, 50, 50, 58, 50, 58, 74,104, 68, 70]
# LTBHEC(6,127): B1 B2 B3 B4 B5 BH  (htcalc.f DATA), one row per LS SI curve.
const ON_LTBHEC = Float32[
 2.9435 0.9132 -0.0141 1.6580 -0.1095 0.0
 3.3721 0.8407 -0.0150 2.6208 -0.2661 0.0
 6.1308 0.6904 -0.0195 10.1563 -0.5330 0.0
 1.0645 0.9918 -0.0812 1.5754 -0.0272 0.0
 2.2835 0.9794 -0.0054 0.5819 -0.0281 0.0
 6.0522 0.6768 -0.0217 15.4232 -0.6354 0.0
 1.5980 1.0000 -0.0198 0.9824 0.0000 0.0
 1.7902 0.9522 -0.0173 1.1668 -0.1206 0.0
 2.4321 0.9207 -0.0168 1.5247 -0.1042 0.0
 1.8326 1.0015 -0.0207 1.4080 -0.0005 0.0
 29.7300 0.3631 -0.0127 16.7616 -0.6804 0.0
 1.5768 0.9978 -0.0156 0.6705 0.0182 0.0
 4.1492 0.7531 -0.0269 14.5384 -0.5811 0.0
 4.2286 0.7857 -0.0178 4.6219 -0.3591 0.0
 1.6505 0.9096 -0.0644 125.7045 -0.8908 0.0
 1.2898 0.9982 -0.0289 0.8546 0.0171 0.0
 1.6622 0.8860 -0.1113 66.8363 -1.0584 0.0
 2.2349 0.8420 -0.0808 15.0884 -0.6292 0.0
 1.5932 1.0124 -0.0122 0.6245 0.0130 0.0
 3.5384 0.7932 -0.0244 29.2355 -0.7291 4.5
 1.0902 1.0298 -0.0354 0.7011 0.1178 0.0
 6.0673 0.6774 -0.0391 9.6135 -0.5504 0.0
 1.1798 1.0000 -0.0339 0.8117 -0.0001 0.0
 1.2673 1.0000 -0.0331 1.1149 0.0001 0.0
 1.2941 0.9892 -0.0315 1.0481 -0.0368 0.0
 1.2721 0.9995 -0.0256 0.7447 -0.0019 0.0
 1.3213 0.9995 -0.0254 0.8549 -0.0016 0.0
 1.3615 0.9813 -0.0675 1.5494 -0.0767 0.0
 1.9426 0.8886 -0.0640 21.0753 -0.6786 0.0
 1.2834 0.9571 -0.0680 100.0000 -0.9223 0.0
 1.4475 1.0098 -0.0200 0.8489 0.0146 0.0
 5.2188 0.6855 -0.0301 50.0071 -0.8695 0.0
 2.1567 0.9374 -0.0155 2.4997 -0.2338 4.5
 5.1493 0.6948 -0.0248 20.9210 -0.7143 0.0
 7.1846 0.6781 -0.0222 13.9186 -0.5268 0.0
 2.1037 0.9140 -0.0275 3.7962 -0.2530 0.0
 1.2866 0.9962 -0.0355 1.4485 -0.0316 0.0
 0.4737 1.2905 -0.0236 0.0979 0.6121 0.0
 4.9676 0.7459 -0.0154 7.3640 -0.5144 0.0
 0.5605 1.3105 -0.0145 0.1779 0.4323 0.0
 4.5598 0.8136 -0.0132 2.2410 -0.1880 0.0
 1.6763 0.9837 -0.0220 0.9949 0.0240 0.0
 1.0945 0.9938 -0.0755 2.5601 0.0114 0.0
 1.3466 0.9590 -0.0574 8.9538 -0.3454 0.0
 1.3295 0.9565 -0.0668 16.0085 -0.4157 0.0
 1.9044 0.9752 -0.0162 0.9262 0.0000 0.0
 1.5403 1.0006 -0.0216 1.0616 -0.0044 0.0
 6.1785 0.6619 -0.0241 25.0185 -0.7400 0.0
 2.9989 0.8435 -0.0200 3.4635 -0.3020 0.0
 0.9680 1.0301 -0.0468 0.1639 0.4127 0.0
 4.7633 0.7576 -0.0194 6.5110 -0.4156 0.0
 1.0370 1.1906 -0.0030 0.1391 0.2655 0.0
 6.4362 0.6827 -0.0194 10.9767 -0.5477 0.0
 1.0438 1.0708 -0.0222 1.5915 -0.1068 4.5
 2.0770 0.9303 -0.0285 2.8937 -0.1414 0.0
 8.1200 0.6748 -0.0111 6.4229 -0.4586 4.5
 1.5341 1.0013 -0.0208 0.9986 -0.0012 0.0
 0.9276 1.0591 -0.0424 0.3529 0.3114 0.0
 1.1151 1.0000 -0.0504 1.3076 0.0009 0.0
 1.5470 1.0000 -0.0225 1.1129 0.0000 0.0
 2.5878 0.8755 -0.0220 3.5505 -0.2882 0.0
 1.0716 1.1276 0.0179 1.2495 -0.0029 4.5
 6.7791 0.6876 -0.0280 12.1447 -0.4142 0.0
 8.6744 0.9986 -0.0031 0.8904 -0.0006 0.0
 2.0854 0.9547 -0.0427 3.8109 -0.2288 4.5
 4.4803 0.7382 -0.0289 26.6132 -0.6461 4.5
 11.3079 0.5419 -0.0345 34.1568 -0.6078 0.0
 1.3342 1.0008 -0.0401 1.8068 0.0248 0.0
 20.3317 0.4649 -0.0275 7.2043 -0.5166 4.5
 1.7620 1.0000 -0.0201 1.2307 0.0000 0.0
 16.2120 0.4580 -0.0135 8.0105 -0.5321 0.0
 2.9232 0.8737 -0.0178 4.6529 -0.3319 4.5
 1.3307 1.0442 -0.0496 3.5829 0.0945 0.0
 1.6330 1.0000 -0.0223 1.2419 0.0000 0.0
 2.0141 0.8989 -0.0236 7.9469 -0.5084 4.5
 2.3102 0.8350 -0.0512 152.9926 -1.0625 0.0
 1.2660 1.0034 -0.0365 1.5515 -0.0221 0.0
 1.4232 0.9989 -0.0285 1.2156 0.0088 0.0
 1.4809 1.0145 -0.0142 0.6040 0.0242 0.0
 0.8977 1.0624 -0.0398 0.1419 0.4709 0.0
 2.2746 0.9019 -0.0220 2.5540 -0.2208 0.0
 1.7327 0.9998 -0.0384 1.1430 -0.0004 0.0
 1.4105 0.9977 -0.0747 2.0668 -0.0104 0.0
 1.1557 1.0031 -0.0408 0.9807 0.0314 0.0
 1.6072 0.9192 -0.0390 12.9271 -0.6427 0.0
 1.3379 0.9999 -0.0534 0.9608 -0.0001 0.0
 1.1510 1.0000 -0.0375 0.8394 0.0000 0.0
 2.0185 0.8815 -0.0902 18.8088 -0.5520 0.0
 1.2144 0.9999 -0.0916 1.7968 -0.0002 0.0
 1.8395 1.0001 -0.0335 1.0774 0.0001 0.0
 1.4698 0.9998 -0.0645 1.7575 0.0009 0.0
 1.4210 0.9947 -0.0269 1.1344 -0.0109 0.0
 1.1672 1.0010 -0.0459 1.3460 0.0204 0.0
 1.3196 1.0000 -0.0356 1.4271 0.000017 0.0
 1.8900 1.0000 -0.0198 1.3892 0.0000 0.0
 19.0635 0.5885 -0.00081 3.3922 -0.3418 4.5
 2.0401 1.0003 -0.0361 1.7914 -0.0090 0.0
 2.6359 0.8259 -0.0389 21.5578 -0.6271 0.0
 0.7666 1.0909 -0.0733 3.2335 -0.2947 4.5
 2.0434 0.9978 -0.0147 1.0937 -0.0035 0.0
 13.6713 0.5404 -0.0283 8.7720 -0.5308 4.5
 1.1266 1.0051 -0.0367 0.6780 0.0404 0.0
 1.9660 1.0000 -0.0240 1.8942 0.0000 0.0
 3.2425 0.7980 -0.0435 52.0549 -0.7064 0.0
 0.6689 1.1249 -0.0785 1.9792 -0.1365 4.5
 1.5269 0.9955 -0.0550 1.9118 0.00048 0.0
 1.4638 0.9979 -0.0743 2.2460 -0.0040 0.0
 1.2096 1.0027 -0.0671 1.2282 0.0335 0.0
 1.1727 1.0042 -0.0439 1.3558 0.0240 0.0
 1.1421 1.0042 -0.0374 0.7632 0.0358 0.0
 3.0849 0.8076 -0.0341 26.2342 -0.6702 0.0
 1.1643 0.9999 -0.0413 1.1057 -0.0002 0.0
 1.1579 1.0000 -0.0930 1.4274 0.0001 0.0
 2.9579 0.8274 -0.0581 9.3221 -0.4616 0.0
 1.0060 1.1098 -0.0535 0.5548 0.2433 0.0
 1.0107 1.1384 -0.0393 0.4584 0.2413 0.0
 1.1519 1.0000 -0.1003 1.6640 0.0001 0.0
 1.5177 1.0000 -0.0551 1.4360 -0.0001 0.0
 1.5861 0.9999 -0.0390 0.9753 -0.0017 0.0
 1.1547 0.9973 -0.0915 1.2294 0.0029 0.0
 2.7644 0.9991 -0.0090 0.6293 0.0003 0.0
 1.2697 1.0002 -0.0885 2.0238 0.0006 0.0
 1.2096 1.0140 -0.0380 1.3247 0.0374 0.0
 1.1204 0.9984 -0.0597 2.4448 -0.0284 0.0
 0.7716 1.1087 -0.0348 0.1099 0.5274 0.0
 1.9730 1.0000 -0.0154 1.0895 0.0000 0.0
 2.1493 0.9979 -0.0175 1.4086 -0.0008 0.0]

# ── HONER: total (VN) + merch sawlog (VM) cubic ft.  volont.f ENTRY HONER. ────
@inline function on_honer(sp::Int, d_in::Float32, h_ft::Float32,
                          topd_in::Float32, stmp_ft::Float32)
    j = Int(ON_ISPMPH[sp])
    d_cm = d_in * ON_INtoCM
    h_m  = h_ft * ON_FTtoM
    td_cm = topd_in * ON_INtoCM
    st_m  = stmp_ft * ON_FTtoM
    b1 = ON_HB1[j]
    totcu = 0.0043891f0 * d_cm * d_cm * (1f0 - 0.04365f0*b1)^2 /
            (ON_HC1[j] + 0.3048f0*ON_HC2[j]/h_m)
    x = td_cm^2 * (1f0 + st_m/h_m) / (d_cm^2 * (1f0 - 0.04365f0*b1)^2)
    sawcu = totcu * (ON_HR1[j] + ON_HR2[j]*x + ON_HR3[j]*x*x)
    totcu = totcu * ON_M3toFT3
    sawcu = sawcu * ON_M3toFT3
    return totcu, sawcu   # VN(GTV), VM(GMV)
end

# ── ZAKVOL: Zakrzewski taper-integral cubic volume (ground/stump → U).  ───────
# ltotal=true  → total, zero stump (B=0, U=H);  false → merch (B=stump, U=ZHT).
@inline function on_zakvol(sp::Int, d_in::Float32, h_ft::Float32, u_ft::Float32,
                           bark::Float32, ltotal::Bool, stmp_ft::Float32)
    d = d_in * ON_INtoCM
    h = h_ft * ON_FTtoM
    u = u_ft * ON_FTtoM
    h < 1.40f0 && return 0f0
    b = ltotal ? 0f0 : stmp_ft * ON_FTtoM
    j = Int(ON_ISPMPZ[sp])
    b2 = ON_ZB2[j]; g2 = ON_ZG2[j]
    u > h && (u = h)
    c = d * bark
    a2 = 1f0 + (h / d)
    z0 = 1f0 - 1.3f0 / h
    y11 = z0 - a2
    y22 = z0^2 + b2*z0^3 + g2*_on_p4(z0)
    b12 = y11 / y22
    k = 7.853982f-5 * c^2 * b12
    h2 = h*h; h3 = h2*h; h4 = h3*h
    lu = on_logf(-h + u + a2*h)
    lb = on_logf(-h + b + a2*h)
    tt = (-1f0/12f0 * (12f0 * h4 * a2^3 * lu * b2 +
          3f0 * g2 * _on_p4(u) + 6f0 * h2 * g2 * a2^2 * u^2 -
          12f0 * h * g2 * u^3 - 12f0 * h3 * g2 * a2 * u -
          4f0 * b2 * h * u^3 - 12f0 * h3 * a2 * b2 * u +
          12f0 * b2 * h2 * u^2 - 12f0 * h3 * a2^2 * g2 * u +
          18f0 * g2 * h2 * u^2 - 12f0 * h3 * a2^2 * b2 * u +
          12f0 * h4 * _on_p4(a2) * lu * g2 -
          12f0 * h3 * a2^3 * g2 * u - 4f0 * h * g2 * a2 * u^3 -
          12f0 * h3 * g2 * u + 6f0 * h2 * b2 * a2 * u^2 -
          12f0 * h3 * a2 * u + 12f0 * h4 * a2^2 * lu -
          12f0 * u * h3 - 12f0 * h3 * b2 * u +
          12f0 * h2 * g2 * a2 * u^2 + 6f0 * h2 * u^2) / h3 +
          1f0/12f0 * (-12f0 * h3 * b2 * b -
          4f0 * h * g2 * a2 * b^3 +
          12f0 * h4 * a2^3 * lb * b2 -
          12f0 * h3 * a2 * b2 * b - 12f0 * h3 * a2^2 * g2 * b -
          12f0 * h3 * a2^2 * b2 * b - 12f0 * h3 * a2^3 * g2 * b -
          12f0 * h3 * g2 * a2 * b - 12f0 * h * g2 * b^3 -
          4f0 * b2 * h * b^3 + 12f0 * h2 * g2 * a2 * b^2 +
          18f0 * g2 * h2 * b^2 + 6f0 * h2 * b2 * a2 * b^2 +
          12f0 * h4 * _on_p4(a2) * lb * g2 +
          6f0 * h2 * g2 * a2^2 * b^2 - 12f0 * h3 * a2 * b +
          6f0 * h2 * b^2 + 12f0 * h4 * a2^2 * lb +
          12f0 * b2 * h2 * b^2 + 3f0 * g2 * _on_p4(b) - 12f0 * b * h3 -
          12f0 * h3 * g2 * b) / h3)
    return tt * k * ON_M3toFT3
end

# ── ZAKHT: height (ft) to the merch top diameter for the ZAK merch-vol call. ──
@inline function on_zakht(sp::Int, d_in::Float32, h_ft::Float32, bark::Float32, topd_in::Float32)
    d = d_in * ON_INtoCM
    h = h_ft * ON_FTtoM
    dib = d * bark
    u = max(1.0f0, topd_in * ON_INtoCM)
    j = Int(ON_ISPMPZ[sp])
    b2 = ON_ZB2[j]; g2 = ON_ZG2[j]
    ca_h = 7.853982f-5 * u * u
    z0 = 1f0 - 1.3f0 / h
    s = 1f0 + (h / d)
    k = ((z0 - s) * 7.853982f-5 * dib^2) / (z0^2 + b2*z0^3 + g2*_on_p4(z0))
    t1 = 1f0/g2
    t3 = sqrt(3f0)
    t5 = b2*b2
    t6 = t5*k
    t7 = k*ca_h
    t9 = s*g2
    t11 = ca_h*ca_h
    t13 = s*t5
    t15 = k*k
    t17 = g2*g2
    t21 = k*t11
    t22 = t5*b2
    t26 = ca_h*b2
    t30 = s*s
    t35 = s*t15*k
    t43 = t17*g2
    t48 = t15*ca_h
    t56 = t5*t5
    t60 = 27f0*t11*ca_h*t17 - ca_h*t5*t15 - 4f0*t21*t22 +
          4f0*ca_h*g2*t15 - 80f0*t26*t15*s*g2 +
          128f0*ca_h*t30*t17*t15 - 16f0*t35*g2 + 6f0*t21*t13*g2 -
          192f0*t21*b2*t30*t17 - 256f0*t21*t30*s*t43 +
          18f0*t21*b2*g2 + 18f0*t48*t22*s - 144f0*t21*s*t17 -
          144f0*t48*t30*g2*t5 + 27f0*t48*t30*t56 + 4f0*t35*t5
    t62 = sqrt(ca_h * t60)
    t69 = on_powf((9f0*t7*b2 - 72f0*t7*t9 + 27f0*t11*g2 +
                   27f0*t7*t13 + 2f0*t15 + 3f0*t3*t62)/t15/t43, 1f0/3f0)
    t71 = k*g2
    t73 = on_powf(54f0, 1f0/3f0)
    t74 = t73*t73
    t75 = t69*t69
    t79 = t73*ca_h
    t84 = 1f0/k
    t86 = 1f0/t69
    t88 = sqrt((27f0*t6*t69 - 72f0*t71*t69 + 2f0*t74*t75*k*t17 +
                36f0*t79*b2 + 144f0*t79*t9 + 12f0*t73*k)*t84*t86)
    t90 = sqrt(6f0)
    t92 = t69*t88
    t99 = t88*t73
    t105 = t3*t69
    t118 = sqrt((27f0*t6*t92 - 72f0*t71*t92 - t88*t74*t75*k*t17 -
                 18f0*t99*t26 - 72f0*t99*ca_h*s*g2 - 6f0*t99*k -
                 324f0*t105*b2*k*g2 - 648f0*t105*ca_h*t17 +
                 81f0*t105*t22*k)*t84*t86/t88)
    r = -b2*t1/4f0 - t3*t1*t88/36f0 + t90*t1*t118/36f0
    zht = h*(1f0-r)
    (zht < 0f0 || zht > h) && (zht = 0f0)
    return zht * ON_MtoFT
end

# ── Mowraski cull: NMV = GMV·(1-(1-exp(-AM·age))^BM).  volont.f ENTRY MOWRASKI. ─
@inline function on_mowraski(sp::Int, vm::Float32, age::Float32)
    j = Int(ON_ISPMPM[sp])
    x = 1f0 - on_powf(1f0 - on_expf(-ON_MAM[j]*age), ON_MBM[j])
    return vm * x
end

# ── Age dubbing from initial height (HTCALC MODE0=0 + FINDAG HTMAX-1.1 retry). ─
@inline function _on_htcalc_htmax(sp::Int, si::Float32)
    idx = Int(ON_MAPLS[sp])
    b1 = ON_LTBHEC[idx,1]; b2 = ON_LTBHEC[idx,2]
    return b1 * on_powf(si, b2)
end
@inline function _on_htcalc_age(sp::Int, h::Float32, si::Float32)
    idx = Int(ON_MAPLS[sp])
    b1 = ON_LTBHEC[idx,1]; b2 = ON_LTBHEC[idx,2]; b3 = ON_LTBHEC[idx,3]
    b4 = ON_LTBHEC[idx,4]; b5 = ON_LTBHEC[idx,5]; bh = ON_LTBHEC[idx,6]
    # AGET = 1/B3 * ln(1 - ((H-BH)/B1/SI^B2)^(1/B4/SI^B5))     (htcalc.f, N<=0 branch)
    base = ((h - bh) / b1) / on_powf(si, b2)
    expo = (1f0/b4) / on_powf(si, b5)
    return (1f0/b3) * on_logf(1f0 - on_powf(base, expo))
end
"""
    on_tree_age(sp, h_ft, si) -> age

ON effective tree age dubbed from initial height (findag.f → htcalc.f MODE0=0). `si` is the
species SITEAR. Mirrors FINDAG's HTMAX-1.1 fallback when H is within 1 ft of the species HTMAX
(the age formula is unstable there). Returns 0 if HTMAX-H<=1 even after the fallback (then FVS
leaves ABIRTH unset → Mowraski uses age 0 ⇒ full cull), matching cratet.f `IF(SITAGE>0)ABIRTH=`.
"""
@inline function on_tree_age(sp::Int, h_ft::Float32, si::Float32)
    htmax = _on_htcalc_htmax(sp, si)
    h = h_ft
    if htmax - h <= 1f0
        h = htmax - 1.1f0
        # HTCALC recomputes with the reduced H; if still within 1 ft it returns AGET=0 (GO TO 900).
        (htmax - h <= 1f0) && return 0f0
    end
    aget = _on_htcalc_age(sp, h, si)
    return aget
end

# ── driver: per-tree VN/VM/BFV → tree.cuft_vol/merch_cuft_vol/bdft_vol. ───────
"""
    compute_volumes_on!(s) -> s

ON (Ontario) per-tree volume — canada/on vols.f + varvol.f (METHC=METHB=8 default) + volont.f.
Fills `cuft_vol` (GTV total merch cubic), `merch_cuft_vol`/`saw_cuft_vol` (GMV merch sawlog
cubic), and `bdft_vol` (NMV = Mowraski-culled GMV, the metric .sum board column). Broken-top
top-kill (CFTOPK/BFTOPK) is not exercised by ont01 and is left unported (scoped follow-on).
"""
function compute_volumes_on!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; c = s.control; p = s.plot
    dbhmin = c.sp_dbh_min; bfmind = c.sp_bf_dbhmin; bftopd = c.sp_bf_topd
    topd = c.sp_top_diam; stmp = c.sp_stump_ht
    @inbounds for i in 1:t.n
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
            continue
        end
        bark = on_bratio(sp, d, h)
        # OCFVOL: below the pulpwood-minimum DBH the cubic call returns the 0.0001 sentinel.
        if d < dbhmin[sp]
            t.cuft_vol[i] = 0.0001f0; t.merch_cuft_vol[i] = 0.0001f0
            t.saw_cuft_vol[i] = 0.0001f0; t.bdft_vol[i] = 0f0
            continue
        end
        if ON_IS_ZAK[sp]
            vn = on_zakvol(sp, d, h, h, bark, true, stmp[sp])          # GTV: total, zero stump, U=H
            zht = on_zakht(sp, d, h, bark, topd[sp])
            vm = zht > 0f0 ? on_zakvol(sp, d, h, zht, bark, false, stmp[sp]) : 0f0
        else
            vn, vm = on_honer(sp, d, h, topd[sp], stmp[sp])            # HON: total + merch
        end
        vn < 0f0 && (vn = 0f0)
        vm < 0f0 && (vm = 0f0)
        t.cuft_vol[i] = vn
        t.merch_cuft_vol[i] = vm
        t.saw_cuft_vol[i] = vm
        # Board section (NMV): only when D >= BFMIND and D > BFTOPD (vols.f board gate).
        if d >= bfmind[sp] && d > bftopd[sp]
            # ABIRTH dubbed-and-frozen on first pass (cratet seeds it once before projection).
            if t.birth_age[i] <= 0f0
                si = (sp >= 1 && sp <= length(p.sp_site_index)) ? p.sp_site_index[sp] : 0f0
                age = on_tree_age(sp, h, si)
                age > 0f0 && (t.birth_age[i] = age)
                t.age_known[i] = true
            end
            nmv = on_mowraski(sp, max(vm, 0f0), t.birth_age[i])
            t.bdft_vol[i] = nmv > 0f0 ? nmv : 0f0
        else
            t.bdft_vol[i] = 0f0
        end
    end
    return s
end
