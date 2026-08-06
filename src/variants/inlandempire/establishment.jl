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

# =============================================================================
# ie_esnspe — AUTOES P(number of species on a stocked plot) (ie/esnspe.f, task #143 chunk A2a).
# Returns PSPE(1..6) = probability that 1..6 species regenerate on a stocked plot. Six logits in
# BAA/ELEV/REGT/BWAF/TPPLN/aspect(XCOS,XSIN)/SLO/TPP + SPEHAB(ISER,·). PSPE(k) is computed only when
# ITPP≥k (else 0); the caller (estab.f) normalizes to a cumulative and draws the species count via RNG.
# NOTE: XCOS=cos(aspect)·SLO, XSIN=sin(aspect)·SLO (estab.f:480-481) — the SLO-weighted aspect (distinct
# from ESTOCK's plain XCOSAS/XSINAS). VALIDATED bit-exact vs live FVSie (iet01 stand-4, ISER=4/ITPP=2/TPP=2):
# PSPE=(0.543, 0.393, 0, 0, 0, 0) = oracle. SPEHAB verbatim from esblkd.f. See docs/AUTOES_CHUNK_PLAN.md.
# =============================================================================

# SPEHAB[iser, j] (esblkd.f, ESCOM2 SPEHAB(5,4)) — series 1-5 (DF/GF/WRC/WH/SAF) × species-count index 1-4.
const _IE_SPEHAB = Float32[ 0.0       0.0       0.0       0.0;
                           -0.695637  0.436955  0.677341  2.301903;
                           -0.776415  0.426625  0.900422  2.602609;
                           -1.227597  0.363575  1.290210  3.089499;
                           -1.058980  0.721468  0.900652  2.156324]

"""
    ie_esnspe(iser, itpp, tpp, tppln, baa, elev, regt, bwaf, xcos, xsin, slo) -> NTuple{6,Float32}

IE P(k species on a stocked plot), k=1..6 (ie/esnspe.f). `xcos`/`xsin` = SLO-weighted aspect cos/sin
(=cos(asp)·SLO, sin(asp)·SLO); `tpp`/`tppln` = trees-per-plot and ln(tpp); `iser` = habitat series (1-5).
Entries beyond `itpp` are 0 (the count logits are gated on ITPP≥k, matching estab.f's PSPE init to 0).
"""
function ie_esnspe(iser::Integer, itpp::Integer, tpp::Real, tppln::Real, baa::Real, elev::Real,
                   regt::Real, bwaf::Real, xcos::Real, xsin::Real, slo::Real)::NTuple{6,Float32}
    ba = Float32(baa); el = Float32(elev); rg = Float32(regt); bw = Float32(bwaf)
    xc = Float32(xcos); xs = Float32(xsin); sl = Float32(slo); tp = Float32(tpp); tpl = Float32(tppln)
    sh(j) = (1 <= iser <= 5) ? _IE_SPEHAB[iser, j] : 0f0
    p1 = 0f0; p2 = 0f0; p3 = 0f0; p4 = 0f0; p5 = 0f0; p6 = 0f0
    # P(1 species) — always
    pn = 1.399594f0 + 0.002162f0*ba + 0.0213903f0*el - 0.022173f0*rg - 0.039405f0*bw -
         1.017904f0*tpl + sh(1)
    p1 = 1f0 / (1f0 + exp(-pn))
    if itpp >= 2
        pn = -0.441879f0 + 0.577760f0*xc + 0.294070f0*xs - 0.128766f0*sl - 0.010128f0*tp -
             0.011040f0*el + 0.016568f0*rg - 0.000247f0*bw + sh(2)
        p2 = 1f0 / (1f0 + exp(-pn))
    end
    if itpp >= 3
        pn = -2.052740f0 + 0.075063f0*xc - 0.471162f0*xs - 0.840378f0*sl - 0.018585f0*el +
             0.076352f0*bw + 0.004745f0*rg + sh(3) + 0.054661f0*tp - 0.000467f0*tp*tp
        p3 = 1f0 / (1f0 + exp(-pn))
    end
    if itpp >= 4
        pn = -6.24551f0 - 1.277404f0*xc - 0.329897f0*xs + 0.425239f0*sl + 0.071878f0*bw +
             0.001371f0*rg + sh(4) + 0.056784f0*tp - 0.000239f0*tp*tp
        p4 = 1f0 / (1f0 + exp(-pn))
    end
    if itpp >= 5
        pn = -2.043228f0 - 0.052907f0*el + 0.021884f0*tp
        p5 = 1f0 / (1f0 + exp(-pn))
    end
    if itpp >= 6
        pn = -5.938035f0 + 0.020070f0*tp
        p6 = 1f0 / (1f0 + exp(-pn))
    end
    return (p1, p2, p3, p4, p5, p6)
end

# =============================================================================
# ie_espadv — AUTOES P(advance-regen species) (ie/espadv.f, task #143 chunk A2b).
# Per-species probability of ADVANCE regeneration (10 species WP/WL/DF/GF/WH/RC/LP/ES/AF/PP):
# PADV(i) = 1/(1+exp(-PNᵢ)) · occ(i), where occ(i)=OCURHT(IHAB,i)·XESMLT(i)·OCURNF(IFO,i) (the
# habitat/forest occupancy) and PNᵢ is a per-species regression in XCOS/XSIN/SLO/TIME/BAA/BAASQ/BAALN/
# ELEV/ELEVSQ/REGT/BWAF + CHAB(IHAB,i) + CPRE(IPREP,i) + OVER(i)>9.95 & forest & IPHY bumps.
# ★ TIME here = years-since-disturbance passed to the probability call, MEASURED =1.0 for iet01 stand-4
# (NOT the estab.f:606 TIME=10 literal — that is overwritten before ESPADV runs). VALIDATED bit-exact
# (3 dp) vs live FVSie: iet01 stand-4 (IHAB=10,IPREP=1,IFO=4,TIME=1,BAA=1,ELEV=34,occ=1,OVER<9.95,IPHY≠1)
# → PADV=(.062 .005 .048 .485 .283 .122 .001 .014 .039 0) = oracle. Tables verbatim esblkd.f:45-76.
# =============================================================================

# CHAB[ihab, species] (16×10, esblkd.f DATA, col-major parsed). species: WP WL DF GF WH RC LP ES AF PP.
const _IE_CHAB = Float32[
  0.0        0.0        0.0        0.0        0.0  0.0        1.9319803  0.0  0.0       0.0;
  0.0        0.0        0.0        0.0        0.0  0.0        0.0        0.0  0.0       0.0;
  0.0        0.0        0.0        0.0        0.0  0.0        0.0        0.0  0.0       0.0;
  0.0        0.0        0.0        0.0        0.0  0.0        0.0        0.0  0.0       0.0;
  0.0        1.4934765 -0.6200184  0.674118   0.0  0.0        0.0        0.0  0.0      -2.5921190;
 -1.224798   0.0       -0.6200184  0.0        0.0  0.0        1.3270519  0.0  0.0      -2.5921190;
  0.0        0.0       -1.1547343  0.0        0.0  0.0        1.3270519  0.0  0.0      -2.5921190;
 -1.224798   0.0       -1.1547343 -0.3889028  0.0  0.0        0.0        0.0  0.0      -0.7036813;
 -0.4296644  0.0       -2.1232    -0.7004423  0.0  0.0        0.0        0.0  0.0      -2.5921190;
 -1.224798   0.0       -2.3086    -0.7004423  0.0  0.6572366  0.0        0.0  1.404156  0.0;
  0.0        0.0       -1.36119   -1.587698   0.0  0.0        1.9319803  0.0  2.40966   0.0;
 -1.224798   1.4934765 -0.13552   -2.06589    0.0  0.0        1.9319803  0.0  2.40966   0.0;
 -0.4296644  0.0       -1.36119   -0.8935856  0.0  0.0        0.0        0.0  2.40966   0.0;
  0.0        0.0       -1.36119   -0.8935856  0.0  0.0        0.0        0.0  2.40966   0.0;
 -0.4296644  0.0       -1.36119   -2.06589    0.0  0.0        1.3270519  0.0  3.23979   0.0;
  0.0        0.0        0.0        0.0        0.0  0.0        0.0        0.0  3.23979   0.0]

# CPRE[iprep, species] (4×10, esblkd.f DATA). iprep: 1=NONE 2=MECH 3=BURN 4=ROAD.
const _IE_CPRE = Float32[
  0.0        0.0  0.0        0.0        0.0        0.0        0.0  0.0        0.0        0.0;
 -0.8986977  0.0 -0.9135126 -0.9944621 -1.1545026 -0.5937042 0.0 -0.7526562 -0.7226382 -0.5682936;
 -0.8739164  0.0 -1.5960879 -1.9561505 -3.1841443 -1.3850973 0.0 -0.9553895 -1.4396372 -0.7261403;
  0.0        0.0 -1.2667466 -1.3659292 -1.4654944 -1.5701306 0.0 -0.1442125 0.0        0.0]

"""
    ie_espadv(ihab, iprep, ifo, iphy, xcos, xsin, slo, time, baa, elev, regt, bwaf, occ, over) -> NTuple{10,Float32}

IE P(advance-regen species) (ie/espadv.f). `occ[i]` = OCURHT(ihab,i)·XESMLT(i)·OCURNF(ifo,i) occupancy
multiplier; `over[i]` = per-species overstory BA (the >9.95 bumps). `time` = years since disturbance
(MEASURED value at the prob call, not the estab.f:606 literal). Species order WP WL DF GF WH RC LP ES AF PP.
"""
function ie_espadv(ihab::Integer, iprep::Integer, ifo::Integer, iphy::Integer, xcos::Real, xsin::Real,
                   slo::Real, time::Real, baa::Real, elev::Real, regt::Real, bwaf::Real, bwb4::Real,
                   occ::AbstractVector, over::AbstractVector)::NTuple{10,Float32}
    xc = Float32(xcos); xs = Float32(xsin); sl = Float32(slo); tm = Float32(time)
    ba = Float32(baa); el = Float32(elev); rg = Float32(regt); bw = Float32(bwaf); b4 = Float32(bwb4)
    basq = ba*ba; elsq = el*el; baln = ba > 0f0 ? log(ba) : 0f0
    ch(i) = (1 <= ihab <= 16) ? _IE_CHAB[ihab, i] : 0f0
    cp(i) = (1 <= iprep <= 4) ? _IE_CPRE[iprep, i] : 0f0
    ov(i) = over[i] > 9.95f0
    logistic(pn) = 1f0 / (1f0 + exp(-pn))
    p = zeros(Float32, 10)
    # WP(1) — skip if IPREP>3
    if iprep <= 3
        pn = -1.8733029f0 + ch(1) - 1.5893204f0*xc - 3.7865858f0*xs - 3.9780395f0*sl -
             0.228477f0*tm + 0.0251953f0*ba - 0.0003143f0*basq + 0.0385460f0*el + cp(1)
        ov(1) && (pn += 0.6721733f0); (ifo == 7 || ifo == 16) && (pn -= 1.1379313f0)
        p[1] = logistic(pn) * Float32(occ[1])
    end
    # WL(2) — skip if IPREP>2
    if iprep <= 2
        pn = -5.0092864f0 + 0.2560369f0*xc + 1.5965058f0*xs + ch(2) - 0.1719499f0*sl
        ov(2) && (pn += 2.1039474f0)
        p[2] = logistic(pn) * Float32(occ[2])
    end
    # DF(3)
    pn = -0.691118f0 + ch(3) + cp(3) - 0.224932f0*xc - 0.461642f0*xs + 1.390878f0*sl +
         0.0151871f0*ba - 0.0000814f0*basq - 0.0137814f0*el
    ov(3) && (pn += 1.0047915f0)
    p[3] = logistic(pn) * Float32(occ[3])
    # GF(4)
    pn = -2.8426666f0 + ch(4) + cp(4) + 1.3658730f0*sl + 0.0072484f0*ba - 0.0000307f0*basq +
         0.1643756f0*el - 0.0020789f0*elsq - 0.1182024f0*rg - 0.1092534f0*bw
    ov(4) && (pn += 1.0538200f0)
    p[4] = logistic(pn) * Float32(occ[4])
    # WH(5)
    pn = -1.8378316f0 + 3.8978596f0*xc - 0.4192431f0*xs - 0.0317794f0*sl + cp(5)
    ov(5) && (pn += 1.1918564f0)
    p[5] = logistic(pn) * Float32(occ[5])
    # RC(6)
    pn = -0.3039302f0 + 0.6771111f0*xc - 1.0103344f0*xs + 1.7631934f0*sl + 0.0057510f0*ba -
         0.0390010f0*el - 0.0818480f0*tm + cp(6)
    ifo == 4 && (pn -= 1.1554776f0); ov(6) && (pn += 1.4044088f0); iphy == 1 && (pn += 1.1103909f0)
    p[6] = logistic(pn) * Float32(occ[6])
    # LP(7) — skip if IPREP>3
    if iprep <= 3
        pn = -7.8876414f0 + ch(7) - 1.7125410f0*sl + 0.2981326f0*baln + 0.0457937f0*el
        ov(7) && (pn += 2.4799900f0)
        p[7] = logistic(pn) * Float32(occ[7])
    end
    # ES(8)
    pn = -12.2236052f0 + 0.2787733f0*el - 0.0019948f0*elsq - 0.0924103f0*rg + 0.1547716f0*b4 -
         0.0478975f0*bw + cp(8)
    ov(8) && (pn += 1.4522984f0)
    (ifo == 19 || ifo == 20) && (pn += 1.1484108f0); (ifo == 14 || ifo == 16) && (pn += 0.9481026f0)
    ifo == 4 && (pn += 0.8911886f0); ifo == 5 && (pn += 0.9723801f0)
    p[8] = logistic(pn) * Float32(occ[8])
    # AF(9) — skip if IPREP>3
    if iprep <= 3
        pn = -14.9235325f0 + ch(9) + 0.2327975f0*xc - 0.8445729f0*xs + 0.0810013f0*sl + cp(9) +
             0.0053038f0*ba + 0.4097059f0*el - 0.0603046f0*rg - 0.0542131f0*bw - 0.0033011f0*elsq
        ov(9) && (pn += 1.6659029f0)
        p[9] = logistic(pn) * Float32(occ[9])
    end
    # PP(10) — skip if IPREP>3
    if iprep <= 3
        pn = -2.0755525f0 - 2.5323539f0*xc - 0.5925702f0*xs - 0.5511553f0*sl + 0.0227489f0*ba -
             0.0002398f0*basq - 0.1262596f0*rg - 0.0857334f0*bw + ch(10) + cp(10)
        (ifo == 19 || ifo == 20) && (pn += 1.5321411f0); ov(10) && (pn += 1.0523320f0)
        p[10] = logistic(pn) * Float32(occ[10])
    end
    return (p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10])
end

# =============================================================================
# ie_espxcs — AUTOES P(excess-regen species) (ie/espxcs.f, task #143 chunk A2b).
# Same 10-species logistic × occupancy pattern as ie_espadv, with the EXCESS tables FHAB(16,10)/FPRE(4,10)
# (esblkd.f:117-153). Quirks: WH(5) uses FPRE(IPREP,7) & has no FHAB; various per-species IFO/OVER/IPHY bumps.
# VALIDATED BIT-EXACT (3 dp, all 10 sp) vs live FVSie (iet01 stand-4, TIME=1): PXCS=(.045 .005 .080 .327 .242
# .194 .043 .001 .025 0) = oracle. See docs/AUTOES_CHUNK_PLAN.md.
# =============================================================================

# FHAB[ihab, species] (16×10, esblkd.f DATA col-major). species WP WL DF GF WH RC LP ES AF PP.
const _IE_FHAB = Float32[
  0.0        1.0573224  0.7860276  0.0        0.0  0.0         0.0       0.0        0.0        0.0;
  0.0        0.0        0.7860276  0.0        0.0  0.0         0.0       0.0        0.0        0.0;
  0.0        0.0        0.0        0.0        0.0  0.0        -1.675261  0.0        0.0        0.0;
  0.0        0.0        0.0        0.0        0.0  0.0        -0.404111  0.0        0.0        0.0;
 -0.6746051  1.0573224  0.0        0.0        0.0  0.0        -0.404111  0.0        0.0       -0.8002471;
 -0.6746051  1.0573224  0.0       -1.14939    0.0  0.0        -0.404111  0.0        0.0       -0.8002471;
  0.0        0.4331883  0.0       -0.5557902  0.0  0.0        -0.404111  0.9266407  0.0       -1.6319422;
 -0.6746051  0.0       -0.6077863 -1.14939    0.0  0.0        -1.675261  0.0        0.0       -0.8002471;
  0.0        0.433188   0.0        0.0        0.0  0.0        -1.716393  1.859201   0.0       -1.9000533;
  0.0        1.2050298  0.4930977  0.0        0.0 -0.4357035  -0.465368  2.333042   2.101272   0.0;
  0.0        1.556329  -0.6364138 -2.43159    0.0  0.0         0.786260  0.612418   1.8053228  0.0;
 -0.6746051  1.556329   0.5219281 -2.43159    0.0  0.0         0.786260  0.612418   1.8053228  0.0;
 -0.6746051  1.556329   0.0       -0.88282    0.0  0.0        -0.653006  1.852301   3.371353   0.0;
  0.0        0.0       -0.6364138 -2.43159    0.0  0.0         0.786260  0.0        0.0        0.0;
  0.0        1.556329   0.0       -2.43159    0.0  0.0        -0.653006  1.852301   3.371353   0.0;
  0.0        0.0        0.0        0.0        0.0  0.0        -0.653006  1.852301   2.101272   0.0]

# FPRE[iprep, species] (4×10, esblkd.f DATA). iprep 1=NONE 2=MECH 3=BURN 4=ROAD.
const _IE_FPRE = Float32[
  0.0  0.0        0.0        0.0        0.0        0.0        0.0        0.0        0.0        0.0;
  0.0  0.7580304 -0.4115679 -0.4249178 -0.2623001 -0.0647744  0.0321509  0.9904144 -0.5096514  0.0;
  0.0  0.6728193 -0.4939831 -0.8300082 -0.5655652  0.1818555 -0.6495806  0.5408921 -0.8132929  0.0;
  0.0  1.9293740  0.0397986 -0.2121915 -0.1897423  0.5362320 -0.5745220  1.5293173  0.2181278  0.0]

"""
    ie_espxcs(ihab, iprep, ifo, iphy, xcos, xsin, slo, time, baa, elev, regt, bwaf, bwb4, occ, over) -> NTuple{10,Float32}

IE P(excess-regen species) (ie/espxcs.f). Same occupancy/args as [`ie_espadv`](@ref); uses FHAB/FPRE (excess
tables). Note WH(5) reads FPRE(iprep,7) and has no FHAB term (faithful to the Fortran).
"""
function ie_espxcs(ihab::Integer, iprep::Integer, ifo::Integer, iphy::Integer, xcos::Real, xsin::Real,
                   slo::Real, time::Real, baa::Real, elev::Real, regt::Real, bwaf::Real, bwb4::Real,
                   occ::AbstractVector, over::AbstractVector)::NTuple{10,Float32}
    xc = Float32(xcos); xs = Float32(xsin); sl = Float32(slo); tm = Float32(time)
    ba = Float32(baa); el = Float32(elev); rg = Float32(regt); bw = Float32(bwaf); b4 = Float32(bwb4)
    basq = ba*ba; elsq = el*el; baln = ba > 0f0 ? log(ba) : 0f0
    fh(i) = (1 <= ihab <= 16) ? _IE_FHAB[ihab, i] : 0f0
    fp(i) = (1 <= iprep <= 4) ? _IE_FPRE[iprep, i] : 0f0
    ov(i) = over[i] > 9.95f0
    logistic(pn) = 1f0 / (1f0 + exp(-pn))
    p = zeros(Float32, 10)
    # WP(1)
    pn = -2.6601112f0 - 0.2103586f0*xc - 2.1766529f0*xs - 2.6159747f0*sl - 0.0170345f0*ba +
         0.3166977f0*baln + fh(1)
    ov(1) && (pn += 1.0546575f0); (ifo == 7 || ifo == 10 || ifo == 16) && (pn -= 1.8648362f0); ifo == 5 && (pn -= 0.5629999f0)
    p[1] = logistic(pn) * Float32(occ[1])
    # WL(2)
    pn = -15.2532959f0 + 1.5877491f0*xc + 0.4421725f0*xs - 1.5027288f0*sl + 0.4384918f0*el -
         0.0051305f0*elsq + fp(2) + fh(2)
    (ifo == 9 || ifo == 10 || ifo == 14 || ifo == 16) && (pn += 2.3989265f0); ov(2) && (pn += 1.2872436f0)
    p[2] = logistic(pn) * Float32(occ[2])
    # DF(3)
    pn = -2.0080204f0 - 0.1351578f0*xc - 0.3944477f0*xs + 0.8531865f0*sl - 0.0383463f0*el +
         0.0652054f0*tm + fh(3) + fp(3)
    (ifo == 3 || ifo == 9 || ifo == 12 || ifo == 16) && (pn += 1.3155589f0); ov(3) && (pn += 0.7341939f0)
    p[3] = logistic(pn) * Float32(occ[3])
    # GF(4)
    pn = -6.1448393f0 + fh(4) + 1.4774529f0*xc + 0.2616096f0*xs - 0.4769181f0*sl + fp(4) -
         0.0048863f0*ba + 0.2789069f0*el - 0.0036567f0*elsq + 0.0582383f0*rg + 0.0502622f0*bw + 0.1356740f0*baln
    ov(4) && (pn += 0.4348936f0); (ifo == 14 || ifo == 16) && (pn -= 0.5188152f0); (ifo == 19 || ifo == 20) && (pn -= 0.7681130f0)
    p[4] = logistic(pn) * Float32(occ[4])
    # WH(5) — FPRE(iprep,7), no FHAB
    pn = -9.3195372f0 + 4.3167849f0*xc - 0.9011113f0*xs - 0.1339375f0*sl + 0.3964820f0*el -
         0.0055853f0*elsq + 0.0896454f0*tm + fp(7)
    ov(5) && (pn += 1.0252551f0)
    p[5] = logistic(pn) * Float32(occ[5])
    # RC(6)
    pn = -1.2917893f0 + 1.9611115f0*xc - 0.0809641f0*xs - 0.6731737f0*sl + 0.0717763f0*tm + fp(6) + fh(6)
    ov(6) && (pn += 1.3999580f0)
    p[6] = logistic(pn) * Float32(occ[6])
    # LP(7)
    pn = -2.6488557f0 + fh(7) + 0.9309435f0*xc - 0.2925614f0*xs - 2.4103925f0*sl - 0.3734260f0*baln +
         0.0142129f0*el + fp(7)
    ov(7) && (pn += 2.5861046f0); iphy == 1 && (pn += 1.1320554f0)
    p[7] = logistic(pn) * Float32(occ[7])
    # ES(8)
    pn = -26.3272057f0 + 0.7454571f0*el - 0.0064948f0*elsq - 1.9590315f0*sl + 0.0703796f0*tm + fh(8) + fp(8)
    ov(8) && (pn += 1.2224044f0)
    p[8] = logistic(pn) * Float32(occ[8])
    # AF(9)
    pn = -7.4072008f0 + fh(9) + 1.0363630f0*xc + 0.2825538f0*xs - 1.5201763f0*sl + 0.0569783f0*el -
         0.1950940f0*b4 + fp(9)
    (ifo == 3 || ifo == 9 || ifo == 10 || ifo == 11 || ifo == 12 || ifo == 14 || ifo == 16) && (pn += 1.2027771f0)
    ov(9) && (pn += 1.5097677f0)
    p[9] = logistic(pn) * Float32(occ[9])
    # PP(10) — no FPRE
    pn = -18.8911858f0 + 0.6606922f0*el - 0.0068996f0*elsq + fh(10)
    ov(10) && (pn += 0.9117771f0)
    p[10] = logistic(pn) * Float32(occ[10])
    return (p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10])
end

# =============================================================================
# IEEstabRNG / ie_esrann! — AUTOES establishment RNG (ie/esrann.f, task #143 chunk A2c primitive).
# A Park-Miller (Lewis-Goodman-Miller) multiplicative LCG, SEPARATE from the main FVS RNG (seed 55329):
#   ESS1 = mod(16807·ESS0, 2147483647);  SEL = Float32(ESS1 / 2147483648).
# Seeded per stand's establishment (ESRNSD odd-adjusts an even seed). For iet01 stand-4 the live seed is
# 43303. VALIDATED vs live FVSie: from seed 43303, draw #52 = 0.21862 = the oracle EMSQR magnitude (0.219,
# estab.f:646-650 consumes draws #51 sign + #52 magnitude). This exact LCG makes bit-exact end-to-end AUTOES
# feasible — the remaining A2c/A3 work is replicating the driver's ESRANN CALL ORDER, not the RNG itself.
# =============================================================================

mutable struct IEEstabRNG
    ess0::Float64
    function IEEstabRNG(seed::Real = 43303.0)
        s = Float64(seed)
        (mod(s, 2.0) == 0.0) && (s += 1.0)   # ESRNSD: reseed an even seed to odd
        new(s)
    end
end

"""
    ie_esrann!(rng::IEEstabRNG) -> Float32

Next draw from the IE establishment LCG (ie/esrann.f). Advances `rng` state. Returns SEL∈[0,1) as Float32
(matching FVS `REAL(ESS1/2147483648D0)`).
"""
function ie_esrann!(rng::IEEstabRNG)::Float32
    ess1 = mod(16807.0 * rng.ess0, 2147483647.0)
    rng.ess0 = ess1
    return Float32(ess1 / 2147483648.0)
end

# =============================================================================
# ie_ocurht — AUTOES habitat-type-group occupancy (ie/blkdat.f OCURHT(16,MAXSP), task #143).
# 0/1 flag zeroing out species-regen probabilities that cannot occur in a habitat type by definition.
# The `occ` multiplier in ie_espadv/ie_espxcs/ie_espsub = OCURHT(ihab,sp)·XESMLT(sp)·OCURNF(ifo,sp).
# Added species (11-23) have no natural regen ⇒ all 0. VALIDATED vs live FVSie: OCURHT(grp10,·)=[1×sp1-9, 0×sp10+]
# (the debug dump; OCURHT(10,PP=10)=0 is why oracle PADV/PXCS(PP)=0). Cols 1-10 verbatim ie/blkdat.f:80-111.
# =============================================================================

# _IE_OCURHT[ihab, sp] (16×23). Species 1-10 = WP WL DF GF WH RC LP ES AF PP; 11-23 (added) = 0.
const _IE_OCURHT = let m = zeros(Float32, 16, 23)
    # columns 1-10 (habitat rows 1-16), from ie/blkdat.f OCURHT DATA (col-major)
    m[:, 1]  = Float32[0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,0]   # WP
    m[:, 2]  = Float32[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0]   # WL
    m[:, 3]  = Float32[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0]   # DF
    m[:, 4]  = Float32[0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0]   # GF
    m[:, 5]  = Float32[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]   # WH
    m[:, 6]  = Float32[0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0]   # RC
    m[:, 7]  = Float32[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]   # LP
    m[:, 8]  = Float32[0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1]   # ES
    m[:, 9]  = Float32[0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1]   # AF
    m[:, 10] = Float32[1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0]   # PP
    m
end

"""
    ie_ocurht(ihab, sp) -> Float32

IE habitat-type-group occupancy flag (ie/blkdat.f OCURHT). 0 zeroes out a species' regen probability in
a habitat type. Added species (11-23) → 0 (no natural regen).
"""
@inline ie_ocurht(ihab::Integer, sp::Integer)::Float32 =
    (1 <= ihab <= 16 && 1 <= sp <= 23) ? @inbounds(_IE_OCURHT[ihab, sp]) : 0f0

# =============================================================================
# ie_estab_pick_species — AUTOES single species draw from a selection distribution (estab.f:745-753).
# Given a uniform DRAW and per-species selection probs SUMUP (e.g. normalized PADV+PSUB), returns the
# chosen species index: the first J∈1..n-1 with DRAW ≤ cumulative(J), else n. This is the inner pick of
# the NUMSPE-species selection loop (estab.f:742-763). Part of the A2c driver (uses ie_esrann! draws).
# =============================================================================
function ie_estab_pick_species(draw::Real, sumup::AbstractVector)::Int
    d = Float32(draw); n = length(sumup); s = 0f0
    @inbounds for j in 1:n-1
        s += Float32(sumup[j])
        d <= s && return j
    end
    return n
end

# =============================================================================
# ie_estpp — AUTOES trees-per-stocked-plot (ie/estpp.f, task #143 chunk A2c primitive).
# Weibull inverse-CDF: BB = exp(0.697061 + .715825·XCOS − .203452·XSIN − 2.762848·SLO + .041235·REGT +
# .031435·BWAF + habitat-bump); CC=0.6836; TPP = ((−ln(1−VAL))^(1/CC))·BB + 0.9. ITPP=INT(TPP+0.5),
# clamped [1, MAXTPP(IHAB)]. VAL = an ESRANN draw. VALIDATED vs live FVSie: iet01 stand-4 plot-1
# (IHAB=10, draw#53) → ITPP=2 (oracle TREES/PLOT=2). habitat bump: IHAB5-8 +.294473, 9/10 +1.237, >10 +.465788.
# =============================================================================
function ie_estpp(val::Real, ihab::Integer, xcos::Real, xsin::Real, slo::Real, regt::Real, bwaf::Real)::Float32
    bb = 0.697061f0 + 0.715825f0*Float32(xcos) - 0.203452f0*Float32(xsin) - 2.762848f0*Float32(slo) +
         0.041235f0*Float32(regt) + 0.031435f0*Float32(bwaf)
    (ihab > 4 && ihab < 9) && (bb += 0.294473f0)
    (ihab == 9 || ihab == 10) && (bb += 1.237000f0)
    ihab > 10 && (bb += 0.465788f0)
    bb = exp(bb)
    cc = 0.6836f0
    return ((-log(1f0 - Float32(val)))^(1f0/cc)) * bb + 0.9f0
end

# =============================================================================
# ie_autoes_plot_seeds — AUTOES per-plot RNG seed chain (estab.f, task #143 A2c).
# The establishment tally reseeds the ESRANN stream PER PLOT (estab.f:967 ESAVE=INT(DRAW*100000+0.5),
# :1075 CALL ESRNSD). Structure (instrument-confirmed on iet01 stand-4, seed0=43303, wk6=IDUP*NPTIDS=50):
#   • one-time WK6 site-prep fill of `wk6` draws from seed0 (before plot 1);
#   • each plot = a 135-draw body: EMSQR@body-2 (sign,mag), ESTPP@body-3, …, ESAVEGEN@body-135;
#   • seed_{N+1} = odd_adjust(ESAVEGEN_N), where ESAVEGEN_N = INT(draw_last*100000+0.5).
# Plot 1's body shares seed0's stream after the wk6 prefix (so its ESAVEGEN is at draw wk6+135); plots 2+
# start a FRESH reseeded stream (ESAVEGEN at draw 135). VALIDATED BIT-EXACT: seeds [43303,22913,17231,
# 97317,32953,75193] = live FVSie. Returns the per-plot seed vector (each already odd-adjusted).
# =============================================================================
function ie_autoes_plot_seeds(seed0::Integer, nplots::Integer; wk6::Integer = 50, body::Integer = 135)::Vector{Int}
    seeds = Int[]
    s = Int(seed0); iseven(s) && (s += 1)          # ESRNSD odd-adjust of the initial seed
    push!(seeds, s)
    for n in 1:(nplots - 1)
        rng = IEEstabRNG(s)
        ndraw = (n == 1 ? wk6 + body : body)       # plot-1 body is offset by the one-time wk6 fill
        local draw::Float32 = 0f0
        for i in 1:ndraw
            draw = ie_esrann!(rng)
        end
        esave = trunc(Int, draw * 100000f0 + 0.5f0)  # ESAVE = INT(DRAW*100000+0.5)
        iseven(esave) && (esave += 1)                # next-plot ESRNSD odd-adjust
        push!(seeds, esave); s = esave
    end
    return seeds
end

# AUTOES per-habitat caps (estab.f:108-112 DATA), indexed by IHAB (1-16):
#   MAXTPP = max trees per stocked plot; MAXSPP = max species per plot; MAXING = max trees/plot for INGROWTH.
const _IE_MAXTPP = Int[9, 7, 5, 5, 10, 8, 9, 5, 21, 25, 10, 10, 11, 7, 10, 8]
const _IE_MAXSPP = Int[4, 3, 3, 3, 5, 4, 6, 4, 6, 6, 4, 5, 5, 4, 6, 4]
const _IE_MAXING = Int[4, 4, 3, 3, 5, 4, 5, 4, 7, 7, 5, 5, 5, 4, 5, 4]

# =============================================================================
# ie_autoes_tally — AUTOES per-species ingrowth TPA for one tally (ie/estab.f, task #143 chunk A2c).
# Composes the 8 validated primitives + the per-plot seed chain into the full establishment tally, bit-exact
# vs live FVSie. For NCOUNT = nplots plots (=NPTIDS·IDUP), each plot runs a 135-draw body from its chained
# seed: EMSQR@2, ESTPP@3→ITPP (cap MAXING), NUMSPE-WK6@4-9, species-WK6@10-15, ADV/SUBS+heights@16-84,
# excess-WK6@85-134. NUMSPE species (PSPE cumulative) are the "best" (1 tree each, PADV+PSUB SUMUP); the
# ITPP-NUMSPE "excess" trees pick ONLY among the best species weighted by PXCS (estab.f:861-877). Each tree
# = prob1·300/dupnpt TPA. VALIDATED BIT-EXACT on iet01 stand-4 (1999): WP33 WL20 DF7 GF202 WH222 RC50 ES23
# AF27, total 583.7 = oracle. Returns a length-`nsp` per-species TPA vector.
# =============================================================================
function ie_autoes_tally(; seed0::Integer, nplots::Integer, ihab::Integer, iser::Integer, ifo::Integer,
                         iprep::Integer, iphy::Integer, xcos::Real, xsin::Real, slo::Real, elev::Real,
                         baa::Real, regt::Real, bwaf::Real, bwb4::Real, prob1::Real, dupnpt::Real,
                         occ::AbstractVector, over::AbstractVector, time::Real = 1f0,
                         is_ingro::Bool = true, nstore::AbstractVector = Int32[],
                         pnn::AbstractVector = Float32[], nsp::Integer = 23, wk6fill::Integer = 50)
    xc = Float32(xcos); xs = Float32(xsin); sl = Float32(slo); tm = Float32(time)
    padv = collect(ie_espadv(ihab, iprep, ifo, iphy, xc, xs, sl, tm, Float32(baa), Float32(elev),
                             Float32(regt), Float32(bwaf), Float32(bwb4), occ, over))
    pxcs = collect(ie_espxcs(ihab, iprep, ifo, iphy, xc, xs, sl, tm, Float32(baa), Float32(elev),
                             Float32(regt), Float32(bwaf), Float32(bwb4), occ, over))
    sumup_base = zeros(Float32, nsp); sumup_base[1:10] .= padv; sumup_base ./= sum(sumup_base)
    nspnz = count(>(1f-4), sumup_base); maxspp = _IE_MAXSPP[ihab]
    # ITPP cap: FVS uses MAXTPP for disturbances / MAXING for ingrowth (measured: disturbance ITPP=20-25, ingrowth
    # ≤7). BUT jl's per-tree TPA (prob1·300/dupnpt) over-produces at MAXTPP (20-25 trees·prob1 ≫ target) — the real
    # per-tree TPA is ~PLPROB/ITPP so the plot total = PLPROB (prob1 CANCELS via ITPP=PLPROB·DUPNPT/(prob1·300)),
    # INDEPENDENT of ITPP. Until that per-tree-TPA model is ported, MAXING is the empirical best (mean|Δ| 22% vs 36%
    # at MAXTPP). ⇒ TODO: book tpaw = PLPROB/ITPP·scale (total=PLPROB) instead of prob1·300/dupnpt, then use MAXTPP.
    cap = _IE_MAXING[ihab]
    p1 = Float32(prob1); scale = 300f0 / Float32(dupnpt)
    # Per-plot NSTORE/PNN (prior tally's stocked count + PROB1). Empty ⇒ a fresh disturbance (all zeros).
    has_state = length(nstore) == nplots && length(pnn) == nplots
    # Per-plot RNG body = the ACTUAL ESRANN advance per plot (measured live estab.f, EM/IE instrument-replay):
    #   16 [EMSQR(2)+ESTPP(1)+NUMSPE-WK6(6)+species-WK6(6)+1] + 3·nsp [ADV/SUBS(nsp)+heights(2·nsp)]
    #   + 2·MAXTPP[ihab] [excess-WK6]. Validated: IE iet01 (nsp=23,ihab=10,MAXTPP=25)=16+69+50=135 (unchanged);
    #   EM (nsp=19,ihab=3,MAXTPP=5)=16+57+10=83 (measured live per-plot advance = 137−54 = 83, constant). The
    #   old hardcoded 135/69/50 were iet01-specific ⇒ EM (fewer species, smaller MAXTPP) desynced the seed chain.
    #   RESULT: jl ITPP sequence [1,2,1,1,3,3,…] now MATCHES live [1,2,1,1,3,3] bit-exact.
    adv_heights = 3 * nsp
    excess_draws = 2 * _IE_MAXTPP[ihab]
    body_n = 16 + adv_heights + excess_draws
    seeds = ie_autoes_plot_seeds(seed0, nplots; wk6 = wk6fill, body = body_n)
    tally = zeros(Float64, nsp)
    for (n, sd) in enumerate(seeds)
        rng = IEEstabRNG(sd)
        for _ in 1:(n == 1 ? wk6fill : 0); ie_esrann!(rng); end
        ie_esrann!(rng); ie_esrann!(rng)                                     # EMSQR
        itpp = clamp(round(Int, ie_estpp(ie_esrann!(rng), ihab, xc, xs, sl, Float32(regt), Float32(bwaf))), 1, cap)
        ns = has_state ? Int(nstore[n]) : 0
        pn = has_state ? pnn[n] : 0f0
        newtpp = max(0, itpp - ns)
        # ESPROB (estab.f:944-951): a tree at plot-index I gets full PROB1 if new (I>NSTORE); an old tree
        # (I≤NSTORE) gets the increment PROB1-PNN; an ingrowth tally scales ALL trees by NEWTPP/ITPP.
        prob_old = max(p1 - pn, 0.0001f0)
        esprob(i) = is_ingro ? max(p1 * Float32(newtpp) / Float32(itpp), 0.0001f0) : (i <= ns ? prob_old : p1)
        wk6n = ntuple(_ -> ie_esrann!(rng), 6); wk6s = ntuple(_ -> ie_esrann!(rng), 6)
        numspe = 1
        if itpp != 1
            pspe = collect(ie_esnspe(iser, itpp, Float32(itpp), log(Float32(itpp)), Float32(baa),
                                     Float32(elev), Float32(regt), Float32(bwaf), xc, xs, sl))
            cum = cumsum(pspe ./ sum(pspe)); numspe = 6
            for i in 1:5; wk6n[i] <= cum[i] && (numspe = i; break); end
        end
        numspe = min(numspe, maxspp, nspnz)
        su = copy(sumup_base); ibest = zeros(Int, nsp)
        iplot = 0                                                            # tree index within the plot (1..itpp)
        for i in 1:numspe
            iplot += 1; j = ie_estab_pick_species(wk6s[i], su); tally[j] += esprob(iplot) * scale; ibest[j] = 1
            su[j] = 0f0; t = sum(su); t > 0 && (su ./= t)
        end
        we = zeros(Float32, nsp)
        for i in 1:10
            we[i] = pxcs[i] * ibest[i]; (ibest[i] == 1 && we[i] < 0.0001f0) && (we[i] = 0.0001f0)
        end
        twe = sum(we); twe > 0 && (we ./= twe)
        for _ in 1:adv_heights; ie_esrann!(rng); end                         # ADV/SUBS(nsp)+heights(2·nsp)
        wk6e = ntuple(_ -> ie_esrann!(rng), excess_draws)                    # excess-WK6 (2·MAXTPP[ihab])
        nd = 0
        for _ in 1:(itpp - numspe)
            nd += 1; iplot += 1; j = ie_estab_pick_species(wk6e[nd], we); tally[j] += esprob(iplot) * scale; nd += 1
        end
        has_state && (nstore[n] = Int32(itpp); pnn[n] = p1)                  # carry to the next tally
    end
    return tally
end

# =============================================================================
# ie_autoes_seed0 — derive the AUTOES tally seed0 from the establishment RNG default (estab.f:290-295).
# The establishment RNG is a SEPARATE ESRANN stream (ESRNCM), default ESSS=55329 (esblkd.f:30) — numerically
# the same as the main FVS seed but an INDEPENDENT stream. For the FIRST tally the stream is fresh at ESSS, so
# estab.f:291 draws once and ESDRAW=INT(DRAW*100000+0.5) seeds the plot chain (:295 reseed). VALIDATED: from
# ESSS=55329, draw #1 = 0.433025 → ESDRAW = 43303 = the iet01 stand-4 seed0. ⇒ the AUTOES RNG is fully
# self-contained/deterministic (no main-growth-RNG dependency for the first tally).
# =============================================================================
function ie_autoes_seed0(esss::Integer = 55329)::Int
    rng = IEEstabRNG(esss)
    d = ie_esrann!(rng)                       # first draw off the ESSS stream
    return trunc(Int, d * 100000f0 + 0.5f0)   # ESDRAW = INT(DRAW*100000+0.5)
end

# =============================================================================
# ie_esadvh — AUTOES advance-regen tallest-tree height (estb/esadvh.f, task #143 chunk A2c).
# SHARED establishment code (all variants with regen). HHT = EXP(PN + EMSQR·DILATE·BNORM·σ_sp), per-species
# regression PN in AGELN/BAA/ELEV/aspect/slope + THAB(hab)/TPRE(prep)/TPHY(phys) tables. AGE=3-DELAY-GENTIM
# (≥1); for AUTOES advance DELAY=0, GENTIM=5 ⇒ AGE=1, AGELN=0. DILATE=FIRST(1,sp) (order-statistic sqrt chain:
# 0.1→√0.1→…). VALIDATED BIT-EXACT vs live FVSie (iet01 stand-4): plot-1 WH HHT=0.6542 (oracle 0.654103),
# plot-2 WH (DILATE=√.1) 0.5507 (oracle .550567). Coeffs verbatim esadvh.f. (Reported TALLEST = max(HHT,XMIN+0.2).)
# =============================================================================
const _IE_ADVH_THAB = let m = zeros(Float32, 5, 11)   # THAB[ihtser, sp]; only DF(3)/GF(4) rows nonzero
    m[:, 3] = Float32[-0.00683, 0.12521, 0.16327, 0.26886, 0.0]
    m[:, 4] = Float32[0.0, 0.0, 0.20183, 0.31082, 0.0]
    m
end
const _IE_ADVH_TPRE = let m = zeros(Float32, 4, 11)   # TPRE[iprep, sp]; WH(5)/AF(9)/MH(11)
    m[:, 5]  = Float32[0.0, -0.10356, -1.23036, -0.40522]
    m[:, 9]  = Float32[0.0, -0.20770, -0.12903,  0.18322]
    m[:, 11] = Float32[0.0, -0.10356, -1.23036, -0.40522]
    m
end
const _IE_ADVH_TPHY = let m = zeros(Float32, 5, 11)   # TPHY[iphy, sp]; DF(3)/RC(6)/LP(7)/PP(10)
    m[:, 3]  = Float32[0.04770,  0.41224,  0.25028,  0.23537, 0.0]
    m[:, 6]  = Float32[0.32413,  0.39404,  0.25123,  0.23419, 0.0]
    m[:, 7]  = Float32[-0.28223, -0.99702, -0.47684, -0.20872, 0.0]
    m[:, 10] = Float32[-0.18689, 0.27119,  0.70375,  0.65555, 0.0]
    m
end

"""
    ie_esadvh(sp, emsqr, dilate, agel, bnorm; baa, elev, xcos, xsin, slo, ihtser, iphy, iprep, bwaf=0, bwb4=0) -> Float32

IE advance-regen tallest-tree height (estb/esadvh.f). `agel`=ln(AGE); `dilate`=FIRST(1,sp) dispersion;
`bnorm`=BNORML(ITIME). Faithful per-species transcription. Aspect `xcos`/`xsin` = SLO-weighted (=cos·SLO etc.).
"""
function ie_esadvh(sp::Integer, emsqr::Real, dilate::Real, agel::Real, bnorm::Real; baa::Real, elev::Real,
                   xcos::Real, xsin::Real, slo::Real, ihtser::Integer, iphy::Integer, iprep::Integer,
                   bwaf::Real = 0.0, bwb4::Real = 0.0)::Float32
    al = Float32(agel); ba = Float32(baa); el = Float32(elev); xc = Float32(xcos); xs = Float32(xsin)
    sl = Float32(slo); bw4 = Float32(bwb4); bwf = Float32(bwaf); disp = Float32(emsqr)*Float32(dilate)*Float32(bnorm)
    th(s) = (1 <= ihtser <= 5) ? _IE_ADVH_THAB[ihtser, s] : 0f0
    tp(s) = (1 <= iprep <= 4) ? _IE_ADVH_TPRE[iprep, s] : 0f0
    ty(s) = (1 <= iphy <= 5) ? _IE_ADVH_TPHY[iphy, s] : 0f0
    pn = 0f0; sig = 0f0
    if sp == 1
        pn = 0.05585f0 + 0.84765f0*al - 0.003824f0*ba - 0.02835f0*el - 0.79565f0*xc + 0.39278f0*xs - 0.68673f0*sl; sig = 0.51878f0
    elseif sp == 2
        pn = -1.80559f0 + 1.24136f0*al; sig = 0.54325f0
    elseif sp == 3
        pn = -1.15433f0 + 1.09480f0*al + ty(3) + th(3) - 0.04804f0*el + 0.0004225f0*el*el; sig = 0.63678f0
    elseif sp == 4
        pn = -1.96040f0 + 1.02403f0*al - 0.00233f0*ba + th(3) + 0.04315f0*xc + 0.13456f0*xs - 0.21468f0*sl - 0.05224f0*bw4 - 0.01898f0*bwf; sig = 0.61195f0
    elseif sp == 5
        pn = -0.43269f0 + 0.77433f0*al - 0.00378f0*ba + tp(5); sig = 0.54794f0
    elseif sp == 6
        pn = 2.11552f0 + 0.71766f0*al + ty(6) - 0.17259f0*el + 0.12506f0*xc + 0.63747f0*xs - 0.35258f0*sl + 0.0022033f0*el*el; sig = 0.62044f0
    elseif sp == 7
        pn = -0.59267f0 + 0.88997f0*al + ty(7) + 0.79158f0*xc + 0.49060f0*xs + 0.49071f0*sl; sig = 0.68842f0
    elseif sp == 8
        pn = -2.19638f0 + 1.12147f0*al - 0.002270f0*ba; sig = 0.59475f0
    elseif sp == 9
        pn = -1.69509f0 + 0.87242f0*al - 0.001107f0*ba + tp(9) - 0.06402f0*bw4 + 0.02299f0*bwf - 0.01189f0*xc + 0.15379f0*xs + 0.44637f0*sl; sig = 0.59957f0
    elseif sp == 10
        pn = -6.33095f0 + 0.79936f0*al + ty(10) + 0.06347f0*bwf + 0.19305f0*el - 0.0020058f0*el*el; sig = 0.53813f0
    else  # 11 = MH (uses WH eq w/ TPRE(iprep,11))
        pn = -0.43269f0 + 0.77433f0*al - 0.00378f0*ba + tp(11); sig = 0.54794f0
    end
    return exp(pn + disp*sig)
end

# =============================================================================
# ie_esdlay — AUTOES years-to-germination delay (estb/esdlay.f, task #143). SHARED establishment code.
# Weibull: DELAY=((-ln(1-DRAW))^(1/CC))*BB, then ADVANCE(ias=1): DELAY=(DELAY+3)*(-1) [→negative→clamps to 0];
# SUBSEQUENT(ias=2): DELAY=DELAY-4. Clamp [0,10]. BB/CC by species×overstory-BA (advance BADV/CADV) or ×plot-age
# (subsequent BSUB/CSUB); budworm variants (BBW/CBW…) for DF/GF/ES/AF (bwb4/bwaf). VALIDATED: advance DELAY=0
# for iet01 (BAA=1→IBAA=1, BWB4=0; matches live "DELAY TO GERM=0.0000"). Coeffs verbatim esdlay.f.
# =============================================================================
# BADV[ibaa, sp] (2×11) advance, no budworm; CADV likewise. BSUB[it,sp]/CSUB (3×11) subsequent.
const _IE_DLAY_BADV = Float32[6.699826 9.768223 13.121021 11.269182 13.604594 17.779381 7.358880 13.990273 21.962337 8.986115 13.604594;
                             13.431179 27.100242 22.186540 18.245664 19.605485 24.241344 32.955809 21.779362 32.176727 10.312660 19.605485]
const _IE_DLAY_CADV = Float32[1.262533 1.152577 1.043254 1.122139 1.267057 1.337217 0.912295 1.051296 1.101266 1.068472 1.267057;
                             1.279302 1.319692 1.215651 1.082805 1.287445 1.540663 1.230540 1.416222 1.317022 1.470405 1.287445]
const _IE_DLAY_BSUB = Float32[3.52946 5.23792 4.34376 4.17909 4.33094 4.16284 5.33757 5.36466 3.45725 3.81610 4.33094;
                             7.62339 7.38005 6.55916 5.88262 6.30802 7.52536 6.78727 7.44468 6.34975 5.74622 6.30802;
                             12.79801 10.42350 9.16226 8.49857 8.63060 10.20937 9.45827 9.69507 8.65545 9.36345 8.63060]
const _IE_DLAY_CSUB = Float32[1.71621 3.11598 2.33194 2.47058 1.97408 2.03892 4.16994 2.89777 2.06804 3.01975 1.97408;
                             2.72466 3.04038 2.63560 2.25957 2.35053 3.12279 3.59937 2.80504 2.79933 2.09376 2.35053;
                             3.98359 2.87196 2.21663 2.00065 2.00997 2.68340 2.39138 3.27745 2.28687 1.80925 2.00997]

"""
    ie_esdlay(sp, ias, draw, time, baa; bwb4=0, bwaf=0) -> Float32

IE years-to-germination delay (estb/esdlay.f). `ias`=1 advance / 2 subsequent; `draw`=ESRANN uniform;
`time`=plot age; `baa`=overstory BA. Budworm branches (DF/GF/ES/AF) not yet included (bwb4=bwaf=0 here).
"""
function ie_esdlay(sp::Integer, ias::Integer, draw::Real, time::Real, baa::Real; bwb4::Real = 0, bwaf::Real = 0)::Float32
    d = Float32(draw)
    if ias == 1
        ibaa = Float32(baa) > 25.5f0 ? 2 : 1
        bb = _IE_DLAY_BADV[ibaa, sp]; cc = _IE_DLAY_CADV[ibaa, sp]
        delay = ((-log(1f0 - d))^(1f0/cc)) * bb
        delay = (delay + 3f0) * (-1f0)
    else
        it = time > 12.5 ? 3 : (time > 7.5 ? 2 : 1)
        bb = _IE_DLAY_BSUB[it, sp]; cc = _IE_DLAY_CSUB[it, sp]
        delay = ((-log(1f0 - d))^(1f0/cc)) * bb - 4f0
    end
    return clamp(delay, 0f0, 10f0)
end

# =============================================================================
# ie_espsub — AUTOES P(subsequent-regen species) (estb/espsub.f, task #143). SHARED establishment code.
# Same 10-species logistic × occupancy pattern as ie_espadv, with SUBSEQUENT tables DHAB(16,10)/DPRE(4,10)
# (esblkd.f:78-116) + √(TIME)/BAALN/ALOG(ELEV) terms. Called at estab.f:774 to populate PSUB for the ADV/SUBS
# ICHOI dispatch (a best tree is SUBSEQUENT when its ADV/SUBS draw > PADV/(PADV+PSUB)). Coeffs verbatim espsub.f.
# =============================================================================
const _IE_DHAB = Float32[
 0.0        0.6603087  0.0        0.0        0.0  0.0         0.4017175  0.0        0.0       0.0;
 0.0        0.0        0.0        0.0        0.0  0.0         0.4017175  0.0        0.0       0.6557608;
 0.0        0.0        0.0        0.0        0.0  0.0        -0.7714835  0.0        0.0       0.6557608;
 0.0        0.0        0.0        0.0        0.0  0.0        -0.7714835  0.0        0.0       0.6557608;
 0.0        0.6603087  0.0        0.0        0.0  0.0         0.0        0.3424427  0.0      -0.3939905;
 0.0        0.6603087  0.0       -0.2163169  0.0  0.0         0.6710793  0.3424427  0.0       0.0;
 0.0        0.0        0.0        0.0        0.0  0.0         0.0       -0.7935050  0.0      -0.3939905;
 0.0        0.0        0.0        0.0        0.0  0.0        -0.7714835 -0.7935050  0.0       0.6557608;
 0.6572202  0.7437406  0.0        0.2288325  0.0  0.0        -0.543113   0.0        0.0       0.0;
 1.0375612  0.8959712  0.0        0.0        0.0 -0.6403333  -0.0868212  0.2115143  1.285995  0.0;
 0.0        0.7381694 -1.108112  -2.442296   0.0  0.0         1.7321693  0.0954259  1.036427  0.0;
 0.0        0.7381694 -0.1829031 -2.442296   0.0  0.0         0.6710793  0.0954259  1.036427  0.0;
 0.0        0.7381694 -0.1829031 -0.812497   0.0  0.0        -0.6315747  0.0954259  2.449334  0.0;
 0.0        0.0       -1.108112  -2.442296   0.0  0.0         1.7321693  0.0954259  1.036427  0.0;
 0.0        0.7381694 -1.108112  -2.442296   0.0  0.0         0.4315782  0.5586635  2.449334  0.0;
 0.0        0.0        0.0        0.0        0.0  0.0        -0.6315747  0.5586635  2.449334  0.0]
const _IE_DPRE = Float32[
 0.0  0.0        0.0        0.0        0.0        0.0        0.0        0.0        0.0        0.0;
 0.1547952 0.9831529 0.2274977 0.1840748 0.2282684 0.6319225 0.6404021 1.1533183 0.2232152 0.6788768;
 0.0230249 1.1296833 0.4415658 0.0475563 0.6218803 1.0655327 0.3812159 1.0900837 0.1324877 0.7407734;
 0.5517443 1.5132390 0.6000910 0.5660518 0.6568975 1.4428914 0.7928200 1.4772912 0.9357350 0.8536627]

"""
    ie_espsub(ihab, iprep, ifo, iphy, xcos, xsin, slo, time, baa, baaln, elev, sqregt, sqbwaf, bwb4, occ, over) -> NTuple{10,Float32}

IE P(subsequent-regen species) (estb/espsub.f). `time`=years since disturbance (uses √time); `baaln`=ln(BAA);
`occ[i]`=OCURHT·XESMLT·OCURNF. Faithful transcription incl. per-species IFO/OVER/IPHY bumps.
"""
function ie_espsub(ihab::Integer, iprep::Integer, ifo::Integer, iphy::Integer, xcos::Real, xsin::Real,
                   slo::Real, time::Real, baa::Real, baaln::Real, elev::Real, sqregt::Real, sqbwaf::Real,
                   bwb4::Real, occ::AbstractVector, over::AbstractVector)::NTuple{10,Float32}
    xc=Float32(xcos); xs=Float32(xsin); sl=Float32(slo); tm=Float32(time); ba=Float32(baa)
    baln=Float32(baaln); el=Float32(elev); elsq=el*el; srg=Float32(sqregt); sbw=Float32(sqbwaf); b4=Float32(bwb4)
    st=sqrt(tm); dh(i)=(1<=ihab<=16) ? _IE_DHAB[ihab,i] : 0f0; dp(i)=(1<=iprep<=4) ? _IE_DPRE[iprep,i] : 0f0
    ov(i)=over[i]>9.95f0; lg(pn)=1f0/(1f0+exp(-pn)); p=zeros(Float32,10)
    # WP(1)
    pn=-2.4158523f0+dh(1)-0.0190109f0*xc-1.3862606f0*xs-1.8366897f0*sl+dp(1)-0.0096275f0*ba-0.0480603f0*el+0.5816639f0*st
    (ifo==7||ifo==10||ifo==16)&&(pn-=1.5468744f0); ov(1)&&(pn+=1.1525173f0); p[1]=lg(pn)*Float32(occ[1])
    # WL(2)
    pn=-6.6374979f0+dh(2)+1.4214396f0*xc+0.7795442f0*xs-0.4926099f0*sl+dp(2)-0.0097349f0*ba+0.1260484f0*el-0.0016995f0*elsq
    ov(2)&&(pn+=1.1026379f0); (ifo==9||ifo==10||ifo==14||ifo==16)&&(pn+=1.8727359f0); p[2]=lg(pn)*Float32(occ[2])
    # DF(3)
    pn=-2.6235752f0+dh(3)+dp(3)+0.5771993f0*xc-0.0248810f0*xs+0.2794252f0*sl-0.0058675f0*ba-0.0229114f0*el+0.6120327f0*srg+0.5549276f0*sbw
    ov(3)&&(pn+=0.2350508f0); (ifo==19||ifo==20)&&(pn-=0.6635007f0); (ifo in (7,9,10,11,12,3))&&(pn+=0.6195183f0); p[3]=lg(pn)*Float32(occ[3])
    # GF(4)
    pn=-4.2174864f0+0.7765163f0*xc+0.8718039f0*xs-0.4550694f0*sl-0.0045320f0*ba+0.1369318f0*el-0.0017545f0*elsq+dh(4)+0.555869f0*srg-0.249851f0*b4+0.412380f0*sbw+dp(4)
    ov(4)&&(pn+=0.2759684f0); (ifo==3||ifo==10)&&(pn-=1.8014250f0); (ifo in (14,16,19,20))&&(pn-=1.1073681f0); p[4]=lg(pn)*Float32(occ[4])
    # WH(5)
    pn=-12.3053122f0+2.5642173f0*xc-0.4392772f0*xs-1.4051726f0*sl+dp(5)+0.4534289f0*el-0.0058291f0*elsq+0.7455533f0*st
    ov(5)&&(pn+=0.7780439f0); p[5]=lg(pn)*Float32(occ[5])
    # RC(6)
    pn=-3.6279063f0+dh(6)+1.5773680f0*xc+0.7379712f0*xs-1.7963310f0*sl+dp(6)-0.0039867f0*ba-0.0259884f0*el+0.9247357f0*st
    ov(6)&&(pn+=0.9973589f0); p[6]=lg(pn)*Float32(occ[6])
    # LP(7)
    pn=-0.9637001f0+dh(7)-3.1974275f0*sl-0.6840911f0*xc-0.1853668f0*xs-0.0317628f0*ba+dp(7)
    ov(7)&&(pn+=1.5827047f0); (ifo in (4,5,17,19))&&(pn-=1.4140185f0); iphy==1&&(pn+=0.7254971f0); p[7]=lg(pn)*Float32(occ[7])
    # ES(8)
    pn=-17.6667213f0+dh(8)+1.5222952f0*xc+1.0029485f0*xs-2.0576222f0*sl-0.2404502f0*baln+3.173367f0*log(el)+dp(8)+0.5579416f0*srg+0.7190073f0*sbw
    ov(8)&&(pn+=1.3091471f0); (ifo in (3,11,12))&&(pn-=0.5256087f0); ifo==4&&(pn+=1.1098698f0); (ifo==9||ifo==10)&&(pn+=0.3115149f0)
    ifo==14&&(pn+=2.3169837f0); ifo==19&&(pn-=0.3958037f0); ifo==5&&(pn+=0.9730915f0); ifo==17&&(pn+=0.9142514f0); ifo==20&&(pn+=0.7438718f0); p[8]=lg(pn)*Float32(occ[8])
    # AF(9)
    pn=-7.7326031f0+dh(9)+0.8712742f0*xc+0.4329470f0*xs-1.7512965f0*sl-0.0112333f0*ba+0.0526394f0*el+0.5499482f0*srg+0.4826650f0*sbw+dp(9)
    (ifo==14||ifo==16)&&(pn+=0.8137509f0); ov(9)&&(pn+=0.5798956f0); p[9]=lg(pn)*Float32(occ[9])
    # PP(10)
    pn=-3.9285939f0-0.7116187f0*xc-0.2141754f0*xs+dh(10)+dp(10)-0.6540064f0*sl-0.3319028f0*baln-0.0370099f0*el+0.4603494f0*st
    ov(10)&&(pn+=1.0770060f0); (ifo==16||ifo==17)&&(pn+=1.8764150f0); (ifo in (3,19,20,14))&&(pn+=2.4815869f0); (iphy==4||iphy==5)&&(pn+=0.346372f0); p[10]=lg(pn)*Float32(occ[10])
    return (p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10])
end

# =============================================================================
# AUTOES scheduler (esnutr.f) — decides, per cycle, whether the automatic
# establishment tally fires and with what NTALLY. Mirrors the esnutr.f control
# flow: (1) LAUTAL removal path, (2) 20-yr disturbance continuation, (3) LINGRW
# ingrowth. Measured bit-exact on iet01 stand-4 (see docs/AUTOES_CHUNK_PLAN.md A3).
#
# Returns (fire::Bool, ntally::Int) — `ntally` is the value ESTAB receives (99
# for ingrowth, 1 for a removal disturbance, 2 for a continuation). Mutates
# `est.idsdat`/`est.ntally` to the PERSISTENT state carried to the next cycle
# (ingrowth resets to 0 — LONE=T; a heavy removal persists at 1 so the next
# cycle can continue to 2).
#
# `xtes` = max(ONTREM/ONTCUR TPA, OCVREM/OCVCUR cuft) removal fraction from THIS
# cycle's within-cycle thin. `itrn` = the tree count when establish! runs (post-thin,
# post-growth, pre-regen). `next_year` = IY(ICYC+1) (the cycle-end year).
function ie_autoes_schedule!(est::Establishment, icyc::Integer, year::Integer,
                             next_year::Integer, itrn::Integer, xtes::Real,
                             inv_year::Integer)
    kdt = next_year - 1
    # (1) LAUTAL removal path (esnutr.f:264-289). Fires whenever a within-cycle thin removed ≥THRES1 —
    # INCLUDING at the inventory year (measured: iet01 stand-4 cyc1 THINPRSC XTES=0.5526 → NTALLY=1, IDSDAT=1990).
    if est.lautal && Float32(xtes) >= est.thres1
        lone = est.thres1 <= Float32(xtes) < est.thres2   # THRES1..THRES2 = single tally (resets); heavy persists
        est.idsdat = Int32(year)                          # IDSDAT = IY(ICYC)
        est.ntally = lone ? Int32(0) : Int32(1)
        return (true, 1)
    end
    # (2) 20-yr disturbance continuation (esnutr.f:298-306): TALLYONE→TALLYTWO.
    if est.ntally > 0 && (kdt - Int(est.idsdat)) <= 19
        est.ntally += Int32(1)
        return (true, Int(est.ntally))
    end
    # (3) LINGRW ingrowth (esnutr.f:313-343): bare-stand (ICYC=1) or 40-yr gap.
    if est.lingrw && ((itrn == 0 && icyc == 1) || (next_year - Int(est.idsdat) >= 40))
        est.idsdat = Int32(next_year - 20)                # IDSDAT = IY(ICYC+1) - 20
        est.ntally = Int32(0)                             # LONE=T ⇒ reset after firing
        return (true, 99)
    end
    return (false, 0)
end

# =============================================================================
# ESTAB per-plot habitat/forest index derivation (esplt2.f + estab.f). Maps the
# stand's raw habitat code + forest code to the establishment indices the tally
# needs: IHAB (habitat-type-group 1-16), ISER (habitat series), IFO (est. forest
# code), IPHY (physiographic, 3 default), IPREP (site-prep, 1=none default).
# Stand case (no plot-specific data): IPHAB=IHTYPE, IPHYS=3, IPPREP=1 (esplt2.f:262-268).
const _IE_ESTAB_IEND = Int[269,299,319,335,385,394,399,499,509,515,519,522,523,524,529,564,
                           579,584,589,599,634,637,644,649,659,669,689,699,709,719,739,744,799]
const _IE_ESTAB_MYGRUP = Int[3,1,4,2,4,3,4,3,8,6,8,7,5,7,8,9,10,6,8,5,13,16,11,14,16,12,15,12,14,15,11,15,14]
const _IE_ESTAB_MYHABG = Int[1,1,1,1,2,2,2,2,3,4,5,5,5,5,5,5]           # estab.f:111 MYHABG(16)
const _IE_ESTAB_IFORCD = Int[103,104,105,106,621,110,113,114,116,117,118,109,111,112,412,402,108,102,115,0]
const _IE_ESTAB_IFORST = Int[3,4,5,4,7,10,4,14,16,17,4,9,11,12,19,20,11,9,12,4]

"""
    ie_estab_indices(habitat_code, forest_code) -> (ihab, iser, ifo, iphy, iprep)

Derive the AUTOES per-plot indices for a stand (no plot-specific overrides).
`habitat_code` = the raw STDINFO habitat class (KODTYP/ICL5); `forest_code` = the
raw forest-location code (KODFOR). For iet01 stand-4 (570, 118) → (10, 4, 4, 3, 1).
"""
function ie_estab_indices(habitat_code::Integer, forest_code::Integer)
    ihtype = 16                                              # esplt2.f:52 fallback
    @inbounds for i in 1:33                                  # esplt2.f:45-53 IEND/MYGRUP bracket
        if habitat_code <= _IE_ESTAB_IEND[i]
            ihtype = _IE_ESTAB_MYGRUP[i]; break
        end
    end
    iser = _IE_ESTAB_MYHABG[ihtype]                          # estab.f:492 ISER=MYHABG(IHAB)
    ifo = 4                                                  # estab.f:218 default
    @inbounds for i in 1:20                                  # estab.f:212-215 KODFOR match
        if forest_code == _IE_ESTAB_IFORCD[i]
            ifo = _IE_ESTAB_IFORST[i]; break
        end
    end
    return (ihab = ihtype, iser = iser, ifo = ifo, iphy = 3, iprep = 1)
end

# =============================================================================
# ie_autoes_run — compose the full AUTOES tally from STAND-LEVEL inputs. Derives
# the ESTAB indices (ie_estab_indices), the inventory stocking probability PROB1
# (logistic of ie_estock — NOTE the stocking equation uses the UNWEIGHTED aspect
# cos/sin, whereas the species-probability routines use the slope-weighted aspect),
# the occupancy vectors, then runs ie_autoes_tally. Returns the per-species TPA
# (Float64[23]) plus the derived PROB1 and indices. `aspect` is the raw aspect
# angle (radians); `slo` the slope fraction; `baa` the plot basal area (floored to
# ≥1 for a bare stand — TBAAA, estab.f). For a bare stand (BAA=1) this reproduces
# the iet01 stand-4 ingrowth tally (583.65 TPA) purely from (habitat_code=570,
# forest_code=118, ESSS=55329).
function ie_autoes_run(; habitat_code::Integer, forest_code::Integer, seed0::Integer,
                       dupnpt::Real, slo::Real, aspect::Real, elev::Real, baa::Real,
                       time::Real = 1f0, regt::Real = time, bwaf::Real = 0f0, bwb4::Real = 0f0,
                       esb_shift::Real = 0f0, is_ingro::Bool = true,
                       nstore::AbstractVector = Int32[], pnn::AbstractVector = Float32[],
                       nsp::Integer = 23)
    idx = ie_estab_indices(habitat_code, forest_code)
    sl = Float32(slo); asp = Float32(aspect); tm = Float32(time)
    xc_st = cos(asp); xs_st = sin(asp)                       # ESTOCK: unweighted aspect
    xc_sp = xc_st * sl; xs_sp = xs_st * sl                   # species probs: slope-weighted aspect
    ba = max(Float32(baa), 1f0)                              # TBAAA floor (estab.f)
    # ESTOCK/species-prob TIME + REGT = years since disturbance (ESTIME): tally-1 = 10, tally-2 = 20 (a cycle
    # per tally), ingrowth = 1. NOT 1 for the disturbance re-stocking tallies (that hardcoding under-computed PN
    # ~1.6 logit → PROB1 0.55 vs the real 0.88). REGT = TIME, SQREGT = √TIME (measured, stand4_estock_inputs.txt).
    pn = ie_estock(idx.ihab, idx.iprep, sl, xc_st, xs_st, Float32(elev), ba, log(ba), tm,
                   sqrt(Float32(regt)), sqrt(Float32(bwaf)), Float32(bwb4), idx.ifo)
    # PROB1 = logistic(PN + ESB - ESB1)·STOADJ (estab.f:579); esb_shift = ESB-ESB1 (inventory calibration, cyc1/2).
    prob1 = 1f0 / (1f0 + exp(-(pn + Float32(esb_shift))))
    occ = Float32[ie_ocurht(idx.ihab, s) for s in 1:nsp]
    over = zeros(Float32, 10)
    tally = ie_autoes_tally(seed0 = seed0, nplots = Int(dupnpt), ihab = idx.ihab, iser = idx.iser,
                            ifo = idx.ifo, iprep = idx.iprep, iphy = idx.iphy, xcos = xc_sp, xsin = xs_sp,
                            slo = sl, elev = Float32(elev), baa = ba, regt = Float32(regt),
                            bwaf = Float32(bwaf), bwb4 = Float32(bwb4), prob1 = prob1, dupnpt = Float32(dupnpt),
                            occ = occ, over = over, time = tm, is_ingro = is_ingro, nstore = nstore, pnn = pnn,
                            nsp = nsp)
    return (tally = tally, prob1 = prob1, idx = idx)
end

# =============================================================================
# ie_autoes_establish! — the AUTOES cycle hook. Runs the automatic-establishment
# tally for InlandEmpire stands (esnutr.f→estab.f). Unlike establish! (which needs
# a scheduled PLANT/NATURAL activity), AUTOES fires purely off the esnutr scheduler
# (ie_autoes_schedule!). On a firing it runs ie_autoes_run and creates the tallied
# per-species TPA as fresh seedling records (DBH≈0.1, height floored to XMIN — all
# establishment heights are sub-breast-height so esgent.f's nominal DBH applies).
# `xtes` = the removal fraction from THIS cycle's within-cycle thin (see grow_cycle!).
function ie_autoes_establish!(s::StandState; fint::Float32)::Bool
    (s.variant isa InlandEmpire || s.variant isa EasternMontana) || return false
    est = s.estab
    (est.lautal || est.lingrw) || return false
    nsp = nspecies(s.variant)
    # The estb habitat bracket (esplt2.f IEND/MYGRUP, = _IE_ESTAB_IEND/MYGRUP) keys off ICL5 = the FVS/NI habitat
    # CODE. IE's p.habitat_code IS that code; EM's p.habitat_code is IEMTYP (1-118) whose NI code is EM_JTYPE[iemtyp]
    # (em/habtyp.f:123 KODTYP=JTYPE(IEMTYP)). The OCURHT/CHAB/ESTOCK tables + _IE_ES_XMIN are the SHARED estb data
    # (species 1-10 = WP WL DF GF WH RC LP ES AF PP), so no species crosswalk is needed — EM's dry habitats zero the
    # wet-side estb species (OCURHT verified bit-identical to live FVSem for IHAB=3), and EM's establishing species
    # (WL/DF/LP/ES/AF/PP) sit at the matching estb indices 2,3,7,8,9,10.
    ihab_code = s.variant isa EasternMontana ?
        Int(EM_JTYPE[clamp(Int(s.plot.habitat_code), 1, 118)]) : Int(s.plot.habitat_code)
    per = round(Int, fint)
    year = Int(current_cycle_year(s))
    next_year = year + per
    inv_year = Int(s.control.cycle_year[1])
    # esnutr.f:113 — IDSDAT defaults to IY(1)-20 (20 yrs before inventory), NOT -9999. Without this the 40-yr
    # ingrowth-gap rule (next_year - idsdat ≥ 40) would spuriously fire for every stand at cycle 1.
    est.idsdat == Int32(-9999) && (est.idsdat = Int32(inv_year - 20))
    icyc = Int(s.control.cycle) + 1
    itrn = s.trees.n
    fire, _ntally = ie_autoes_schedule!(est, icyc, year, next_year, itrn, est.last_xtes, inv_year)
    est.last_xtes = 0f0                       # consume the removal fraction (one cycle only)
    fire || return false
    year in est.years_done && return false

    p = s.plot
    nptids = max(1, Int(p.points_inv) - Int(p.nonstockable))
    idup   = max(1, cld(_ES_MINREP, nptids))
    dupnpt = Float32(nptids * idup)
    # AUTOES ESRANN seed chain (estab.f:290-295 + the per-plot ESAVE reseed). A NEW disturbance/ingrowth tally
    # (NTALLY==1|99) draws seed0 = ESRANN(es_stream) from the continuing establishment stream, then advances the
    # stream to the post-tally state = the last plot's ESAVE (ie_autoes_plot_seeds[dupnpt+1], validated == the
    # live ESS0 at the next tally: 55329→ESAVE 78807→…). A continuation (NTALLY≥2) reuses this tally's seed0.
    # es_stream starts fresh at ESSS=55329. (Measured seeds 43303/61677/…; the ESAVE chain is bit-exact.)
    if _ntally == 1 || _ntally == 99
        ess0 = est.es_stream == 0f0 ? 55329.0 : Float64(est.es_stream)
        dr = ie_esrann!(IEEstabRNG(ess0))
        seed0 = floor(Int, dr * 100000f0 + 0.5f0)
        est.es_seed = Float32(seed0)                                          # save for the continuation
        est.es_stream = Float32(ie_autoes_plot_seeds(seed0, Int(dupnpt) + 1)[end])  # ESAVE_50 = next stream state
    else
        seed0 = round(Int, est.es_seed)                                       # continuation reuses seed0
    end
    # ESTOCK/species-prob BAA = the per-INVENTORY-POINT basal area BAAA(NNID) (estab.f:482, dense.f:213), NOT the
    # whole-stand BA. After a heavy overstory removal the regen point is bare → BAAA≈0 → TBAAA=max(BAAA,1)=1 →
    # ESTOCK PN high → PROB1 high (the disturbance re-stocking pulse). Using stand_ba (which keeps the residual
    # overstory) collapsed PROB1 to ~0.55 and under-produced ~40%. For the single-inventory-point stand the point
    # is index 1 (validated: jl point_ba[1]=40 vs live BAAA=41.93 at cyc1; 0→1 vs 1 at the bare disturbance tallies).
    baaa = (isempty(s.density.point_ba) ? 0f0 : s.density.point_ba[1])
    # TIME/REGT = years since the disturbance (ESTIME): a disturbance tally is TIME = next_year − IDSDAT (10 for
    # tally-1, 20 for tally-2, …); an ingrowth tally (NTALLY=99) uses TIME=1 (SHORTY, estab.f:252-253).
    time = _ntally == 99 ? 1f0 : Float32(next_year - Int(est.idsdat))
    # Per-plot NSTORE/PNN state (ESPROB weighting). A NEW disturbance/ingrowth tally (NTALLY==1|99) resets the
    # per-plot stocked counts; a continuation (NTALLY≥2) reuses them so it books only the increment ITPP-NSTORE.
    dupnpt_i = Int(dupnpt)
    is_ingro = _ntally == 99
    if _ntally == 1 || _ntally == 99
        est.es_nstore = zeros(Int32, dupnpt_i)
        est.es_pnn = zeros(Float32, dupnpt_i)
    elseif length(est.es_nstore) != dupnpt_i
        est.es_nstore = zeros(Int32, dupnpt_i); est.es_pnn = zeros(Float32, dupnpt_i)
    end
    # ESB inventory calibration (estab.f:319-326,579) applies ONLY at the inventory-year disturbance (INADV=0);
    # the auto/ingrowth tallies (est.idsdat > inv_year) have ESB=ESB1=0. ESB = logit(clamp(logistic(-5.174+0.851·
    # ln(TPACRE)),0.10,0.90)) from the current small-tree (DBH<REGNBK=2.999) TPA; ESB1 = ESTOCK(BAAOLD, TIME=0).
    # Computed once at the first inventory tally, reused by its continuation (est.esb_shift persists ESB-ESB1).
    esb_shift = 0f0
    if Int(est.idsdat) == inv_year
        if isnan(est.esb_shift)
            idx0 = ie_estab_indices(ihab_code, Int(p.user_forest_code))
            tpacre = 0f0
            @inbounds for i in 1:s.trees.n; s.trees.dbh[i] < 2.999f0 && (tpacre += s.trees.tpa[i]); end
            tpacre < 1f0 && (tpacre = 1f0)
            esa = clamp(1f0 / (1f0 + exp(-(-5.17397f0 + 0.85131f0 * log(tpacre)))), 0.10f0, 0.90f0)
            esb = -log(1f0 / esa - 1f0)
            baaold = max(baaa, 1f0)                          # inventory per-point BA (BAAINV); = baaa at cyc1
            asp0 = Float32(p.aspect); sl0 = Float32(p.slope)
            esb1 = ie_estock(idx0.ihab, idx0.iprep, sl0, cos(asp0), sin(asp0), Float32(p.elevation),
                             baaold, log(baaold), 0f0, 0f0, 0f0, 0f0, idx0.ifo)   # ESTOCK(BAAOLD, TIME=0)
            est.esb_shift = esb - esb1
        end
        esb_shift = est.esb_shift
    end
    r = ie_autoes_run(habitat_code = ihab_code, forest_code = Int(p.user_forest_code), nsp = nsp,
                      seed0 = seed0, dupnpt = dupnpt, slo = p.slope, aspect = p.aspect,
                      elev = p.elevation, baa = max(baaa, 1f0), time = time, esb_shift = esb_shift,
                      is_ingro = is_ingro, nstore = est.es_nstore, pnn = est.es_pnn)

    t = s.trees
    xmin = _IE_ES_XMIN
    created = false
    @inbounds for sp in 1:nsp
        tpa_sp = Float32(r.tally[sp])
        tpa_sp > 0f0 || continue
        hht = xmin[sp] + 0.2f0                           # est. height floor TALL=max(HHT,XMIN+0.2) (estab.f:838);
                                                         # the computed ESADVH/ESSUBH heights (0.14-0.65) fall below
                                                         # it, so XMIN+0.2 is the effective seedling height (per-tree
                                                         # ESADVH/ESSUBH refinement pending — this is the floor value)
        dbh = 0.1f0 + 0.001f0 * hht                      # esgent.f:56 sub-breast-height nominal DBH
        n = t.n + 1; n > length(t.dbh) && break
        t.n = n
        t.species[n]     = Int32(sp)
        t.dbh[n]         = dbh
        t.height[n]      = hht
        t.tpa[n]         = tpa_sp
        t.plot_id[n]     = Int32(1)
        # Crown: the REGENT(LESTB) open-grown crown (regent.f:178) CR=0.89722−0.0000461·PCCF, clamped [0.20,0.90].
        # A near-bare regen stand has PCCF≈0 ⇒ CR≈0.90; a nominal deterministic value here (the ±1% RANN draw is
        # a refinement) keeps the downstream crown/DG models well-defined (crown_ratio=0 produced NaN TopHt).
        pccf = s.density.point_ccf[1]
        cr = clamp(0.89722f0 - 0.0000461f0 * pccf, 0.20f0, 0.90f0)
        icr0 = floor(Int32, cr * 100f0 + 0.5f0)
        t.crown_pct[n]   = icr0
        t.crown_ratio[n] = Float32(icr0)
        t.norm_ht[n]     = Int32(0)
        t.sort_key[n]    = Float64(n)
        created = true
    end
    created || return false
    push!(est.years_done, Int32(year))
    compute_density!(s)
    return true
end
