# =============================================================================
# establishment.jl (inlandempire) — IE ESSUBH subsequent/planted-tree HEIGHT model
# (ie/essubh.f). Assigns heights to trees created by the establishment model
# (PLANT/NATURAL keywords). Per-species: PN = base + b1·ln(AGE) − b2·BAA +
# UHAB[habitat] + UPRE[prep] + UPHY[phys] + aspect/slope/elev terms; then
# HHT = exp(PN + disp·σ_sp), where disp = EMSQR·DILATE·BNORM (the estab-engine
# lognormal dispersion). Special species take fixed heights (0.5 or 5.0 ft).
#
# This is the height model only (coefficients + kernel). Wiring into the shared
# establish! engine (supplying AGE/BAA/IHTSER/IPREP/IPHY + _IE_ES_XMIN/HHTMAX +
# the draw window) + validation vs pure_DF_est is the next step. Not yet
# dispatched — inert until establish! gets an InlandEmpire branch.
# =============================================================================

# ie/blkdat.f establishment constants:
#   XMIN   — per-species minimum (default/natural floor) height (ft).
#   HHTMAX — per-species hard cap on the reported established-tree height (ft).
#   BNORML — age-indexed (IAGE 1..20) lognormal-dispersion multiplier for ESSUBH (BNORM=BNORML(IAGE)).
const _IE_ES_XMIN   = Float32[1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 1.0, 1.0, 0.5, 0.5, 0.5, 1.0, 6.0, 3.0, 6.0, 6.0, 3.0, 0.5]
const _IE_ES_HHTMAX = Float32[23.0, 27.0, 21.0, 21.0, 22.0, 20.0, 24.0, 18.0, 18.0, 17.0, 22.0, 27.0, 27.0, 18.0, 6.0, 6.0, 27.0, 16.0, 16.0, 16.0, 16.0, 16.0, 22.0]
const _IE_ES_BNORML = Float32[1.0, 1.0, 1.0, 1.046, 1.093, 1.139, 1.186, 1.232, 1.278, 1.325, 1.371, 1.418, 1.464, 1.510, 1.557, 1.603, 1.649, 1.696, 1.742, 1.789]

# ie/essubh.f UHAB(5,MAXSP): subsequent-height coef by H.T. group
#   [WET-DF, DRY-DF, GRAND-F, WRC/WH, SAF] × species. Nonzero only sp2,3,7,8,10,12.
const IE_ESSUBH_UHAB = let m = zeros(Float32, 5, 23)
    m[:, 2]  = Float32[-0.01541, -0.03814,  0.11409,  0.35334, 0.0]
    m[:, 3]  = Float32[-0.21858, -0.03354,  0.22756,  0.51988, 0.0]
    m[:, 7]  = Float32[-0.29969, -0.15449,  0.04545, -0.00601, 0.0]
    m[:, 8]  = Float32[ 0.0,      0.0,      0.18740,  0.26511, 0.0]
    m[:, 10] = Float32[-0.02287, -0.14710,  0.19278,  0.13817, 0.0]
    m[:, 12] = Float32[-0.01541, -0.03814,  0.11409,  0.35334, 0.0]   # WB uses WL(2)
    m
end

# ie/essubh.f UPRE(4,MAXSP): subsequent-height coef by site prep [NONE, MECH, BURN, ROAD].
const IE_ESSUBH_UPRE = let m = zeros(Float32, 4, 23)
    m[:, 2]  = Float32[0.0, -0.11310, -0.06246, 0.009632]
    m[:, 3]  = Float32[0.0,  0.06961,  0.19508, 0.17952]
    m[:, 4]  = Float32[0.0, -0.08010,  0.01032, -0.05975]
    m[:, 6]  = Float32[0.0, -0.41961, -0.22326, 0.15608]
    m[:, 7]  = Float32[0.0,  0.11502,  0.02486, 0.13080]
    m[:, 8]  = Float32[0.0,  0.10587,  0.27072, 0.16240]
    m[:, 10] = Float32[0.0,  0.20729,  0.18491, 0.11864]
    m[:, 12] = Float32[0.0, -0.11310, -0.06246, 0.009632]   # WB uses WL(2)
    m
end

# ie/essubh.f UPHY(5,MAXSP): subsequent-height coef by physiographic position
#   [BOTTOM, LOWER, MID, UPPER, RIDGE]. Nonzero only sp1,3,4,7,8.
const IE_ESSUBH_UPHY = let m = zeros(Float32, 5, 23)
    m[:, 1] = Float32[-0.18731, -0.48682, -0.32160, -0.16113, 0.0]
    m[:, 3] = Float32[-0.27801, -0.20433, -0.12317, -0.26736, 0.0]
    m[:, 4] = Float32[-0.06976, -0.16483, -0.10900, -0.15873, 0.0]
    m[:, 7] = Float32[ 0.32401,  0.14743,  0.22165,  0.24559, 0.0]
    m[:, 8] = Float32[ 0.41120,  0.01164,  0.22217,  0.15834, 0.0]
    m
end

"""
    ie_essubh(sp, age, baa, ihtser, iprep, iphy, xcos, xsin, slo, elev, disp; bwaf=0, bwb4=0) -> Float32

IE subsequent/planted-tree height (ie/essubh.f). `age` = tree age (AGE, ≥1); `baa` = stand BA;
`ihtser`/`iprep`/`iphy` = habitat-series/site-prep/physiography index (1-based); `xcos`/`xsin` =
aspect cos/sin; `slo` = slope; `elev` = elevation (100s ft); `disp` = EMSQR·DILATE·BNORM (the
estab-engine lognormal dispersion, drawn there). Returns HHT (ft). Faithful transcription of the
per-species GO TO(I) branches. Special species → fixed 0.5 / 5.0 ft.
"""
function ie_essubh(sp::Integer, age::Real, baa::Real, ihtser::Integer, iprep::Integer, iphy::Integer,
                   xcos::Real, xsin::Real, slo::Real, elev::Real, disp::Real; bwaf::Real = 0.0, bwb4::Real = 0.0)::Float32
    a = Float32(age); a < 1f0 && (a = 1f0)
    aln = log(a); baa = Float32(baa); ih = Int(ihtser); ip = Int(iprep); iph = Int(iphy)
    @inbounds uhab(s) = IE_ESSUBH_UHAB[ih, s]; @inbounds upre(s) = IE_ESSUBH_UPRE[ip, s]; @inbounds uphy(s) = IE_ESSUBH_UPHY[iph, s]
    xcos = Float32(xcos); xsin = Float32(xsin); slo = Float32(slo); elev = Float32(elev); disp = Float32(disp)
    pn = 0f0; sig = 0f0; fixed = -1f0
    if sp == 1                                                                # WP
        pn = -1.51302f0 + 1.24537f0*aln - 0.003052f0*baa + uphy(1); sig = 0.46010f0
    elseif sp == 2 || sp == 12                                               # WL (WB=12 reuses WL)
        pn = -1.36257f0 + 1.21548f0*aln - 0.003797f0*baa + uhab(2) + upre(2); sig = 0.52668f0
    elseif sp == 3                                                           # DF
        pn = -2.16416f0 + 1.28151f0*aln - 0.0031363f0*baa + uhab(3) + upre(3) + uphy(3) -
             0.09626f0*xcos - 0.23946f0*xsin - 0.14589f0*slo; sig = 0.55942f0
    elseif sp == 4                                                           # GF
        pn = -2.62001f0 + 1.19408f0*aln - 0.0035489f0*baa + upre(4) + uphy(4) + 0.01871f0*xcos +
             0.09002f0*xsin - 0.37365f0*slo + 0.05070f0*elev - 0.000736f0*elev*elev; sig = 0.52958f0
    elseif sp == 5 || sp == 11 || sp == 23                                   # WH; MH(11)/OS(23) reuse WH
        pn = -2.42379f0 + 1.52366f0*aln - 0.003256f0*baa; sig = 0.54116f0
    elseif sp == 6                                                           # RC
        pn = -0.89895f0 + 1.08584f0*aln - 0.00205f0*baa + upre(6) - 0.01594f0*elev; sig = 0.56107f0
    elseif sp == 7                                                           # LP
        pn = -0.27105f0 + 1.32027f0*aln - 0.008208f0*baa + upre(7) + uphy(7) + uhab(7) -
             0.15385f0*xcos + 0.04156f0*xsin - 0.49186f0*slo - 0.04744f0*elev + 0.0003511f0*elev*elev +
             0.01105f0*Float32(bwaf) + 0.02588f0*Float32(bwb4); sig = 0.47557f0
    elseif sp == 8                                                           # ES
        pn = -2.93213f0 + 1.43503f0*aln - 0.002504f0*baa + upre(8) + uphy(8) + uhab(8); sig = 0.48951f0
    elseif sp == 9 || sp == 14                                               # AF; LL(14) reuses SAF
        pn = -2.06377f0 + 1.18184f0*aln - 0.0044465f0*baa + 0.06615f0*xcos + 0.03085f0*xsin -
             0.37402f0*slo; sig = 0.56740f0
    elseif sp == 10                                                          # PP
        pn = -1.99480f0 + 1.53946f0*aln - 0.00402f0*baa + uhab(10) + upre(10) - 0.01155f0*elev; sig = 0.49076f0
    elseif sp == 13 || sp == 15 || sp == 16 || sp == 17                      # LM/PI/JU/PY fixed 0.5
        fixed = 0.5f0
    elseif sp == 18 || sp == 19 || sp == 20 || sp == 21 || sp == 22          # AS/CO/MM/PB/OH fixed 5.0
        fixed = 5.0f0
    else
        fixed = 0.5f0
    end
    fixed >= 0f0 && return fixed
    return exp(pn + disp*sig)
end

# ie/esxcsh.f ESXCSH — per-tree height-CLASS distribution: a Weibull inverse-CDF from HTMIN(=XMIN) up to
# HTMAX(=TALL, tallest-subsequent ht), scaled by a random DRAW. Faithful transcription of ie/esxcsh.f.
# ROLE (measured, cont.64): this is the AUTO-ESTABLISHMENT height-class model (estab.f:931). It is NOT
# the NATURAL-keyword height path — live instrumentation showed the single CALL ESXCSH never fires for
# NATURAL regen (XCSHTRC=0 while regen occurs), same as essubh (auto-estab, cont.56). Serves the auto-estab
# mode; the NATURAL-keyword per-tree height source remains unresolved (additive placeholder in establish!).
# Coefficients by time-class (ITIME 1/2/3) × species. sp 13,15-17 → 0.5; sp 18-22 → 5.0 (fixed).
const _IE_ESXCSH_SHIFT = Float32[4.0,4.0,2.0,2.0,2.0,2.0,4.0,2.0,2.0,4.0,2.0, 0,0,0,0,0,0,0,0,0,0,0,0]
# BB/CC (3 time-classes × 23 species). Access [itime, sp].
const _IE_ESXCSH_BB = let m = zeros(Float32, 3, 23)
    m[:,1]=Float32[2.121455,5.060402,5.979549];   m[:,2]=Float32[6.643726,11.422982,19.618871]
    m[:,3]=Float32[3.816083,8.161474,10.987699];  m[:,4]=Float32[3.089571,5.830185,10.105748]
    m[:,5]=Float32[3.347712,6.806825,13.553455];  m[:,6]=Float32[3.169513,4.506403,8.940539]
    m[:,7]=Float32[7.360424,10.928846,25.214411]; m[:,8]=Float32[1.466152,5.159270,9.272780]
    m[:,9]=Float32[2.921356,4.581383,10.333282];  m[:,10]=Float32[2.779221,9.033310,14.131212]
    m[:,11]=Float32[3.347712,6.806825,13.553455]; m[:,12]=Float32[6.643726,11.422982,19.618871]
    m[:,14]=Float32[2.921356,4.581383,10.333282]; m[:,23]=Float32[3.347712,6.806825,13.553455]
    m
end
const _IE_ESXCSH_CC = let m = zeros(Float32, 3, 23)
    m[:,1]=Float32[0.745850,0.782170,0.842171];  m[:,2]=Float32[0.902909,1.166155,1.306380]
    m[:,3]=Float32[0.996732,0.845413,0.948037];  m[:,4]=Float32[0.800681,0.832278,0.954081]
    m[:,5]=Float32[0.567768,0.894628,1.214044];  m[:,6]=Float32[0.640554,0.813543,0.943493]
    m[:,7]=Float32[1.148084,1.232333,1.117025];  m[:,8]=Float32[0.722527,0.739031,1.125510]
    m[:,9]=Float32[0.885137,0.871559,1.043759];  m[:,10]=Float32[0.899325,1.074932,0.930698]
    m[:,11]=Float32[0.567768,0.894628,1.214044]; m[:,12]=Float32[0.902909,1.166155,1.306380]
    m[:,14]=Float32[0.885137,0.871559,1.043759]; m[:,23]=Float32[0.567768,0.894628,1.214044]
    m
end

"""
    ie_esxcsh(sp, htmax, htmin, time, draw) -> Float32

ie/esxcsh.f: per-tree NATURAL height = Weibull inverse-CDF height class from htmin(=XMIN) to
htmax(=TALL), scaled by `draw` (RNG uniform). `time` = plot age (→ITIME 1/2/3). Faithful transcription.
"""
function ie_esxcsh(sp::Integer, htmax::Real, htmin::Real, time::Real, draw::Real)::Float32
    (sp == 13 || (15 <= sp <= 17)) && return 0.5f0
    (18 <= sp <= 22) && return 5.0f0
    itime = time > 12.5 ? 3 : (time > 7.5 ? 2 : 1)
    bb = _IE_ESXCSH_BB[itime, sp]; cc = _IE_ESXCSH_CC[itime, sp]; sh = _IE_ESXCSH_SHIFT[sp]
    (bb <= 0f0) && return Float32(htmin)                     # unaffected species → floor
    class = Float32(htmax) / 0.2f0 - sh
    xuppr = 1f0 - exp(-(((class - Float32(htmin)) / bb)^cc))
    xx = xuppr * Float32(draw)
    hht = ((-log(1f0 - xx))^(1f0 / cc)) * bb + Float32(htmin)
    return 0.2f0 * (hht + sh)
end

# =============================================================================
# ie_estock — AUTOES probability-of-stocking (ie/estock.f, task #143 chunk A1).
# The regeneration establishment model's P(stocking) for the automatic (natural,
# disturbance-triggered) regen tally. Habitat-series dispatch on IHAB → IEQ
# (1=Douglas-fir, 2=grand-fir, 3=cedar/hemlock, 4=subalpine-fir); IPREP>3 uses the
# "all roads" equation. Returns PN (a logit); the caller forms the stocking prob
# FTEMP = 1/(1+exp(-(PN + ESB - ESB1)))·STOADJ (estab.f:579). Coefficients verbatim
# from estock.f DATA statements. VALIDATED bit-exact vs live FVSie instrument on
# iet01 stand-4: IHAB=10→IEQ=3, inputs {SLO=0.30, ASPECT=5.498, ELEV=34, BAA=1,
# IPREP=1, TIME=1, SQREGT=1} → PN=0.2116 (oracle 0.2116). See docs/AUTOES_CHUNK_PLAN.md.
# =============================================================================

const _IE_ESTOCK_SHAB = Float32[1.14975, 0.070196, -0.115832, 0.0, 0.539949, -0.060377, 0.0,
                                -0.178929, 0.0, 0.356001, 0.278596, 0.181755, 0.0, -0.764070,
                                0.403826, -0.837748]
const _IE_ESTOCK_SSER = Float32[0.0, 0.5976365, 1.3547862, 1.5073566, 1.1695859]
# SPRE[iprep, ieq] — site-prep (1=NONE,2=MECH,3=BURN) by habitat series (Fortran SPRE(3,4), col-major).
const _IE_ESTOCK_SPRE = Float32[0.0 0.0 0.0 0.0;
                                -0.180267 -0.161450 -0.185464 0.068146;
                                -0.485674 -0.189325 -0.346349 -0.342800]
const _IE_ESTOCK_FORDF = Float32[0.0, 0.0, 1.077081, 0.0, 0.0, 0.0, 0.0, 0.0, 0.730596, 0.730596,
                                 0.730596, 0.730596, 0.0, 0.0, 0.0, 1.077081, 0.0, 0.0, 0.286133, 0.805539]
const _IE_ESTOCK_FORGF = Float32[0.0, 0.0, 0.0, 0.0, -0.482030, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, -0.415825, 0.0, -1.087558, -1.087558]

"""
    ie_estock(ihab, iprep, slo, xcosas, xsinas, elev, xbaa, xbaaln, time, sqregt, sqbwaf, bwb4, ifo) -> PN

IE regeneration P(stocking) logit (ie/estock.f). `xcosas`/`xsinas` = cos/sin(aspect); `xbaa`/`xbaaln` = BAA
and ln(BAA); `time` = years since disturbance; `sqregt`/`sqbwaf`/`bwb4` = WSBW habitat-code flags; `ifo` =
forest index (1-20). Faithful transcription incl. the IPREP>3 "all roads" branch.
"""
function ie_estock(ihab::Integer, iprep::Integer, slo::Real, xcosas::Real, xsinas::Real, elev::Real,
                   xbaa::Real, xbaaln::Real, time::Real, sqregt::Real, sqbwaf::Real, bwb4::Real,
                   ifo::Integer)::Float32
    xc = Float32(xcosas); xs = Float32(xsinas); el = Float32(elev); ba = Float32(xbaa)
    baln = Float32(xbaaln); srg = Float32(sqregt); sbw = Float32(sqbwaf); b4 = Float32(bwb4)
    elsq = el * el
    shab = (1 <= ihab <= 16) ? _IE_ESTOCK_SHAB[ihab] : 0f0
    fordf = (1 <= ifo <= 20) ? _IE_ESTOCK_FORDF[ifo] : 0f0
    forgf = (1 <= ifo <= 20) ? _IE_ESTOCK_FORGF[ifo] : 0f0
    # IEQ dispatch (estock.f:46-49)
    ieq = 1
    (ihab > 4 && ihab < 9) && (ieq = 2)
    (ihab == 9 || ihab == 10) && (ieq = 3)
    ihab > 10 && (ieq = 4)
    if iprep > 3
        # "all roads" (estock.f:94-100)
        ihab > 9 && (ieq += 1)
        sser = (1 <= ieq <= 5) ? _IE_ESTOCK_SSER[ieq] : 0f0
        return -2.1072788f0 + sser + 0.0111295f0 * ba + 0.4206679f0 * srg -
               0.3558347f0 * b4 + 0.1871430f0 * sbw
    end
    sqslo = sqrt(Float32(slo)); sqsq = sqslo * sqrt(Float32(time))
    spre = _IE_ESTOCK_SPRE[iprep, ieq]
    if ieq == 1                          # Douglas-fir series (estock.f:54-62)
        return -0.829879f0 + shab + 0.074061f0 * xc * sqsq - 0.067207f0 * xs * sqsq -
               0.021187f0 * sqsq + spre - 0.027058f0 * el + 0.213680f0 * srg +
               0.129135f0 * sbw + fordf
    elseif ieq == 2                      # grand-fir series (estock.f:63-73)
        return -2.444807f0 + 0.392834f0 * xc * sqsq + shab + 0.117465f0 * xs * sqsq -
               0.316824f0 * sqsq + spre + forgf + 0.013726f0 * ba - 0.0000822f0 * ba * ba +
               0.076197f0 * el - 0.000839f0 * elsq + 0.555074f0 * sqrt(Float32(time)) -
               0.013542f0 * xc * Float32(slo) * ba - 0.013486f0 * xs * Float32(slo) * ba
    elseif ieq == 3                      # cedar/hemlock series (estock.f:74-84)
        return -6.216694f0 + shab + 0.587967f0 * xc * sqsq + 0.007416f0 * xs * sqsq +
               0.152539f0 * sqsq + spre + 0.007953f0 * ba - 0.0000373f0 * ba * ba +
               0.288594f0 * el - 0.003952f0 * elsq + 0.514807f0 * srg + 0.449858f0 * sbw -
               0.017901f0 * xc * Float32(slo) * ba - 0.006001f0 * xs * Float32(slo) * ba
    else                                 # subalpine-fir series (estock.f:85-93)
        return -0.430349f0 + shab + 0.246248f0 * xc * sqsq - 0.019381f0 * xs * sqsq -
               0.099968f0 * sqsq + spre + 0.136777f0 * baln + 0.239137f0 * srg -
               0.111587f0 * b4 + 0.224696f0 * sbw
    end
end
