# CI establishment/planted-tree base height (ci/essubh.f — "height of tallest subsequent/planted tree").
# Structurally identical to IE (both Inland-Empire-derived); CI restricts to 19 species. Species 1-10 use
# the Carlson (1988) log-age regressions; 11-19 are fixed (0.5 or 5.0) except OS(18) which reuses the WH eq.
#
# HHT = EXP(PN + disp·SIG), where PN = the deterministic species regression and
#   disp = EMSQR·DILATE·BNORM  (estab.f:646-650 / :1018-1039 / essubh.f:95)
#     EMSQR = ±DRAW2 (sign from DRAW1<0.5), drawn ONCE per stand-establishment (2 ESRANN draws)
#     DILATE = FIRST(2,sp): per-species accumulator, init 0.1, then sqrt-shrinks toward 1 each PLANT of that sp
#     BNORM  = BNORML(IAGE): a DETERMINISTIC age-indexed dispersion table (ci/blkdat.f), NOT a random draw
# So the ONLY stochastic input is the per-stand EMSQR (2 draws); the rest is deterministic given call order.
#
# BNORML(IAGE) — ci/blkdat.f:120-121 (index 1..22). IAGE≥1 clamps into range.
const _CI_BNORML = Float32[1.0, 1.0, 1.0, 1.046, 1.093, 1.139, 1.186, 1.232, 1.278, 1.325,
                           1.371, 1.418, 1.464, 1.510, 1.557, 1.603, 1.649, 1.696, 1.742, 1.789]
@inline _ci_bnorml(iage::Integer) = @inbounds _CI_BNORML[clamp(Int(iage), 1, length(_CI_BNORML))]

# UHAB(H.T.group 1..5, sp) / UPRE(prep 1..4, sp) / UPHY(phys 1..5, sp) — ci/essubh.f DATA (nonzero rows only).
# Coefficients are identical to the shared IE species (CI is IE-derived). Stored [index, sp]; zero elsewhere.
const _CI_ESSUBH_UHAB = let m = zeros(Float32, 5, 19)
    m[:,2]  = Float32[-0.01541, -0.03814,  0.11409,  0.35334, 0.0]
    m[:,3]  = Float32[-0.21858, -0.03354,  0.22756,  0.51988, 0.0]
    m[:,7]  = Float32[-0.29969, -0.15449,  0.04545, -0.00601, 0.0]
    m[:,8]  = Float32[ 0.0,      0.0,      0.18740,  0.26511, 0.0]
    m[:,10] = Float32[-0.02287, -0.14710,  0.19278,  0.13817, 0.0]
    m
end
const _CI_ESSUBH_UPRE = let m = zeros(Float32, 4, 19)
    m[:,2]  = Float32[0.0, -0.11310, -0.06246,  0.009632]
    m[:,3]  = Float32[0.0,  0.06961,  0.19508,  0.17952]
    m[:,4]  = Float32[0.0, -0.08010,  0.01032, -0.05975]
    m[:,6]  = Float32[0.0, -0.41961, -0.22326,  0.15608]
    m[:,7]  = Float32[0.0,  0.11502,  0.02486,  0.13080]
    m[:,8]  = Float32[0.0,  0.10587,  0.27072,  0.16240]
    m[:,10] = Float32[0.0,  0.20729,  0.18491,  0.11864]
    m
end
const _CI_ESSUBH_UPHY = let m = zeros(Float32, 5, 19)
    m[:,1]  = Float32[-0.18731, -0.48682, -0.32160, -0.16113, 0.0]
    m[:,3]  = Float32[-0.27801, -0.20433, -0.12317, -0.26736, 0.0]
    m[:,4]  = Float32[-0.06976, -0.16483, -0.10900, -0.15873, 0.0]
    m[:,7]  = Float32[ 0.32401,  0.14743,  0.22165,  0.24559, 0.0]
    m[:,8]  = Float32[ 0.41120,  0.01164,  0.22217,  0.15834, 0.0]
    m
end

"""
    ci_essubh(sp, age, baa, ihtser, iprep, iphy, xcos, xsin, slo, elev, disp; bwaf=0, bwb4=0) -> Float32

CI subsequent/planted-tree base height (ci/essubh.f). `age` = tree AGE (≥1); `baa` = overstory BA;
`ihtser`/`iprep`/`iphy` = habitat-series/site-prep/physiography index (1-based); `xcos`/`xsin` = aspect
cos/sin; `slo` = slope; `elev` = elevation (100s ft); `disp` = EMSQR·DILATE·BNORM dispersion (drawn/derived
in the estab engine). Returns HHT (ft). Faithful transcription of the per-species CASE(I) branches; special
species → fixed 0.5 / 5.0 ft. Coefficients verbatim from ci/essubh.f.
"""
function ci_essubh(sp::Integer, age::Real, baa::Real, ihtser::Integer, iprep::Integer, iphy::Integer,
                   xcos::Real, xsin::Real, slo::Real, elev::Real, disp::Real; bwaf::Real = 0.0, bwb4::Real = 0.0)::Float32
    a = Float32(age); a < 1f0 && (a = 1f0)
    aln = log(a); baa = Float32(baa); ih = Int(ihtser); ip = Int(iprep); iph = Int(iphy)
    @inbounds uhab(s) = _CI_ESSUBH_UHAB[ih, s]; @inbounds upre(s) = _CI_ESSUBH_UPRE[ip, s]; @inbounds uphy(s) = _CI_ESSUBH_UPHY[iph, s]
    xcos = Float32(xcos); xsin = Float32(xsin); slo = Float32(slo); elev = Float32(elev); disp = Float32(disp)
    pn = 0f0; sig = 0f0; fixed = -1f0
    if sp == 1                                                                # WP
        pn = -1.51302f0 + 1.24537f0*aln - 0.003052f0*baa + uphy(1); sig = 0.46010f0
    elseif sp == 2                                                            # WL
        pn = -1.36257f0 + 1.21548f0*aln - 0.003797f0*baa + uhab(2) + upre(2); sig = 0.52668f0
    elseif sp == 3                                                            # DF
        pn = -2.16416f0 + 1.28151f0*aln - 0.0031363f0*baa + uhab(3) + upre(3) + uphy(3) -
             0.09626f0*xcos - 0.23946f0*xsin - 0.14589f0*slo; sig = 0.55942f0
    elseif sp == 4                                                            # GF
        pn = -2.62001f0 + 1.19408f0*aln - 0.0035489f0*baa + upre(4) + uphy(4) + 0.01871f0*xcos +
             0.09002f0*xsin - 0.37365f0*slo + 0.05070f0*elev - 0.000736f0*elev*elev; sig = 0.52958f0
    elseif sp == 5                                                            # WH
        pn = -2.42379f0 + 1.52366f0*aln - 0.003256f0*baa; sig = 0.54116f0
    elseif sp == 6                                                            # RC
        pn = -0.89895f0 + 1.08584f0*aln - 0.00205f0*baa + upre(6) - 0.01594f0*elev; sig = 0.56107f0
    elseif sp == 7                                                            # LP
        pn = -0.27105f0 + 1.32027f0*aln - 0.008208f0*baa + upre(7) + uphy(7) + uhab(7) -
             0.15385f0*xcos + 0.04156f0*xsin - 0.49186f0*slo - 0.04744f0*elev + 0.0003511f0*elev*elev +
             0.01105f0*Float32(bwaf) + 0.02588f0*Float32(bwb4); sig = 0.47557f0
    elseif sp == 8                                                            # ES
        pn = -2.93213f0 + 1.43503f0*aln - 0.002504f0*baa + upre(8) + uphy(8) + uhab(8); sig = 0.48951f0
    elseif sp == 9                                                            # AF
        pn = -2.06377f0 + 1.18184f0*aln - 0.0044465f0*baa + 0.06615f0*xcos + 0.03085f0*xsin -
             0.37402f0*slo; sig = 0.56740f0
    elseif sp == 10                                                           # PP
        pn = -1.99480f0 + 1.53946f0*aln - 0.00402f0*baa + uhab(10) + upre(10) - 0.01155f0*elev; sig = 0.49076f0
    elseif sp == 18                                                           # OS reuses WH eq
        pn = -2.42379f0 + 1.52366f0*aln - 0.003256f0*baa; sig = 0.54116f0
    elseif sp == 13 || sp == 17 || sp == 19                                   # AS/CW/OH fixed 5.0
        fixed = 5.0f0
    else                                                                      # WB/PY/WJ/MC/LM (11,12,14,15,16) fixed 0.5
        fixed = 0.5f0
    end
    fixed >= 0f0 && return fixed
    return exp(pn + disp*sig)
end
