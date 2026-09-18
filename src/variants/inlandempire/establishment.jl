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
# IE AUTOES per-National-Forest occupancy OCURNF(IFO,sp) (ie/blkdat.f) — same gate as EM's, IE's 10 species
# WP/WL/DF/GF/WH/RC/LP/ES/AF/PP. [sp,ifo] 10×20; added species 11-23 → no natural regen (1.0 = inert, gated by
# OCURHT anyway). #143: completes the occ=OCURHT·OCURNF fix for IE (was default 1.0). iet01's forest/species have
# OCURNF=1 (iet01 was bit-exact vs live WITHOUT this ⇒ OCURNF must be 1 there) ⇒ iet01 stays bit-exact.
const IE_AUTOES_OCURNF = let m = ones(Float32, 23, 20)   # default 1.0 for added sp (11-23), OCURHT gates them
    m[1,:]  = Float32[0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0]  # WP
    m[2,:]  = Float32[0,0,1,1,1,0,1,0,1,1,0,0,0,1,0,1,1,0,1,0]  # WL
    m[3,:]  = Float32[0,0,1,1,1,0,1,0,1,1,1,1,0,1,0,1,1,0,1,1]  # DF
    m[4,:]  = Float32[0,0,1,1,1,0,1,0,0,1,0,0,0,1,0,1,1,0,1,1]  # GF
    m[5,:]  = Float32[0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0]  # WH
    m[6,:]  = Float32[0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,1,1,0,0,0]  # RC
    m[7,:]  = Float32[0,0,1,1,1,0,1,0,1,1,1,1,0,1,0,1,1,0,1,1]  # LP
    m[8,:]  = Float32[0,0,1,1,1,0,1,0,1,1,1,1,0,1,0,1,1,0,1,1]  # ES
    m[9,:]  = Float32[0,0,1,1,1,0,1,0,1,1,1,1,0,1,0,1,1,0,1,1]  # AF
    m[10,:] = Float32[0,0,1,1,1,0,1,0,1,0,0,0,0,1,0,1,1,0,1,1]  # PP
    m
end
@inline autoes_ocurnf(::InlandEmpire, ifo::Integer, sp::Integer)::Float32 =
    (1 <= sp <= 23 && 1 <= ifo <= 20) ? @inbounds(IE_AUTOES_OCURNF[sp, ifo]) : 1.0f0

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
# MECHPREP / BURNPREP site preparation (estb/esetpr.f + estab.f:246-399) — ESTAB-packet keyword.
# A MECHPREP/BURNPREP keyword schedules mechanical/broadcast-burn site prep on the AUTAL DISTURBANCE tally
# (NTALLY==1): a fraction of the DUPNPT replicate plots is assigned IPPREP∈{1 NONE,2 MECH,3 BURN} by a
# sample-without-replacement draw off the same WK6 site-prep RNG vector the tally already consumes (estab.f:333
# DO 183). Each prepped plot's regen uses its IPREP index in the advance/excess species-mix (ie_espadv/ie_espxcs
# CPRE prep term) — the prep shifts species COMPOSITION (measured on FVSie: MECHPREP/BURNPREP move the out-year
# BA/QMD, count near-invariant since PROB1's SPRE term is 0 for this habitat series). Heights (UPRE) stay at the
# XMIN+0.2 floor (the computed ESSUBH heights fall below it), so the effect is species-mix only here. TIME is
# prep-invariant when the prep date == the harvest date (TIME=FTEMP−ZMECH=FTEMP−ZHARV).
# =============================================================================

# esetpr.f — resolve the MECHPREP/BURNPREP keyword %-of-plots into (PMECH,PBURN) fractions + the IALN flags.
# pmech_pct / pburn_pct are 0-100 (or `nothing` if that keyword absent). esprin.f always stores ≥1 param so the
# "prep all plots" METH branch is inert (a bare `MECHPREP <date>` gives PMECH=0 — verified on FVSie_estabdump).
function ie_esetpr(pmech_pct, pburn_pct)
    pmech = 0f0; pburn = 0f0; ialn2 = 0; ialn3 = 0
    if pburn_pct !== nothing; pburn = Float32(pburn_pct) / 100f0; ialn3 = 1; end
    if pmech_pct !== nothing; pmech = Float32(pmech_pct) / 100f0; ialn2 = 1; end
    return (pmech = pmech, pburn = pburn, ialn2 = ialn2, ialn3 = ialn3,
            any_kw = (ialn2 == 1 || ialn3 == 1))
end

# estab.f:249-375 — normalize (PNONE,PMECH,PBURN) into cumulative-bucket weights SUMUP(1..3). Keyword path
# (IALN set): mixture is (1−PMECH−PBURN, PMECH, PBURN), renormalized if PMECH+PBURN>1 (estab.f:353-363).
function ie_esetpr_normalize(pnone::Real, pmech::Real, pburn::Real, ialn2::Integer, ialn3::Integer)
    pm = Float32(pmech); pb = Float32(pburn); pn = Float32(pnone)
    if ialn2 == 1 || ialn3 == 1
        s = pm + pb
        if s > 1f0; pm /= s; pb /= s; end
        pn = 1f0 - pm - pb
    end
    tot = pn + pm + pb
    return (pn / tot, pm / tot, pb / tot)
end

# estab.f:382-399 "SAMPLE WITHOUT REPLACEMENT" — assign each of the DUPNPT=nptids·idup replicated plots an
# IPPREP∈{1,2,3} from the normalized bucket weights `sumup` and the WK6 site-prep RNG vector (length DUPNPT,
# already drawn on the :estab stream). BIT-EXACT vs FVSie_estabdump given identical (sumup, WK6).
function ie_esetpr_sample(sumup, wk6::AbstractVector, nptids::Integer, idup::Integer)
    dupnpt = Float32(nptids * idup)
    su = Float32[sumup[1], sumup[2], sumup[3]]
    ipprep = Vector{Int}(undef, nptids * idup)
    n = 0
    @inbounds for _ii in 1:idup, _nn in 1:nptids
        n += 1
        draw = Float32(wk6[n]) * (((dupnpt + 1f0) - Float32(n)) / dupnpt)   # estab.f:388
        s = 0f0; sel = 0
        for i in 1:2                                                        # estab.f:389-393
            s += su[i]
            draw > s && continue
            sel = i; break
        end
        if sel == 0                                                        # estab.f:394-395
            sel = 3; su[3] < 0f0 && (sel = 1)
        end
        ipprep[n] = sel                                                    # estab.f:396
        su[sel] -= 1f0 / dupnpt                                            # estab.f:397 (without replacement)
        su[sel] < 0f0 && (su[sel] = 0f0)
    end
    # estab.f:395-467 DO 203/39/202/201 — the MAIN tally loop does NOT consume IPPREP in the DO-10 build
    # order. For each inventory point NN it reads that point's IDUP replicate preps as a contiguous block
    # (DO 39: ISTART=NN*IDUP-IDUP+1..NN*IDUP) and then PROCESSES them GROUPED BY ASCENDING PREP-TYPE
    # (DO 202 ITYPEP=1,4 → all NONE plots, then MECH, then BURN, then ROAD). NCOUNT — the plot index the
    # per-plot ESTPP/ESNSPE/species draws are booked under — therefore advances in prep-grouped order within
    # each point, NOT in raw replicate order. jl's tally loop consumes ipprep[n] linearly (n's point =
    # div(n-1,idup)+1, so consecutive idup entries = one point), which matches DO 39's contiguous read but
    # NOT DO 202's grouping. Sorting each point's idup-block ascending emits the preps in the exact order the
    # ESTPP draws are paired with them (VALIDATED: bare NOTREES+PLANT-400-DF, 10 pts × 5 reps → jl per-plot
    # IPREP now byte-matches FVSie_g16 for all 50 plots; previously 24/50 diverged, shifting high-ITPP plots
    # onto the GF-favoring NONE prep ⇒ +49 GF / −43 tail-species ⇒ +8 TPA @2002). idup=1 (one replicate per
    # point, the usual FIA case) ⇒ 1-element blocks ⇒ no-op ⇒ byte-identical.
    di = Int(idup)
    if di > 1
        @inbounds for b in 0:(nptids - 1)
            sort!(view(ipprep, b * di + 1 : (b + 1) * di))
        end
    end
    return ipprep
end

# esprep.f — DEFAULT site-prep proportions P(NONE)/P(MECH)/P(BURN) by habitat series, used when the USER
# supplied NO MECHPREP/BURNPREP keyword (estab.f:365-370 → CALL ESPREP(ISER,PNONE,PMECH,PBURN)). Each is a
# per-series logistic in cos/sin(ASPECT)·SLOPE, SLOPE, ln(BA+1) and ELEV. XPREP(prep,iser) coefficients from
# esblkd.f (NONE/MECH/BURN/ROAD × series 1..5). The three probabilities are then NORMALIZED into the SUMUP the
# per-plot IPPREP sample-without-replacement consumes. MEASURED on FVSie_g16 na_def (bare IE dated ESTAB,
# series 4): PNONE/PMECH/PBURN → 0.48/0.26/0.26 (the regen-report SITE PREP SUMMARY), and because the stocking
# logit ie_estock carries a per-IPREP SPRE term (ieq=2 grand-fir: 0/−0.16145/−0.18933), the per-plot PROB1
# splits 0.699/0.663977/0.657729 whose 48/26/26 mean = 0.6792 = the oracle report's stocking probability.
const _IE_XPREP = Float32[  # [iser 1..5, prep(NONE,MECH,BURN,ROAD)] — esblkd.f DATA XPREP (column-major)
    0.0        0.0        0.0        0.0;
    0.085732  -0.226844   0.087087   0.620176;
    0.151164  -0.605840  -0.097998   1.390994;
    0.680760  -0.832692  -0.756665   0.951442;
    0.203387  -0.305596  -0.263820   0.995664]

"""
    ie_esprep(iser, aspect, slope, ba, elev) -> (pnone, pmech, pburn)

IE default site-prep probabilities (estb/esprep.f). `aspect` in radians, `slope` fraction, `ba` stand basal
area, `elev` elevation (hundreds of ft). Faithful transcription; the three are NOT normalized here (the caller
normalizes into the IPPREP SUMUP).
"""
function ie_esprep(iser::Integer, aspect::Real, slope::Real, ba::Real, elev::Real)
    xp(k) = (1 <= iser <= 5) ? _IE_XPREP[iser, k] : 0f0
    asp = Float32(aspect); sl = Float32(slope); el = Float32(elev)
    ca = cos(asp) * sl; sa = sin(asp) * sl; lba = log(Float32(ba) + 1f0)
    pn = 1.043151f0 + xp(1) - 0.220954f0 * ca + 0.369575f0 * sa + 0.769112f0 * sl +
         0.260178f0 * lba - 0.029689f0 * el
    pnone = 1f0 / (1f0 + exp(-pn))
    pn = -1.852031f0 + xp(2) + 0.492668f0 * ca + 0.192020f0 * sa - 0.966674f0 * sl -
         0.085920f0 * lba + 0.024939f0 * el
    pmech = 1f0 / (1f0 + exp(-pn))
    pn = -15.195303f0 + xp(3) + 0.0519477f0 * ca - 0.6135848f0 * sa - 0.0890163f0 * sl -
         0.377915f0 * lba + 0.5303707f0 * el - 0.0049081f0 * el * el
    pburn = 1f0 / (1f0 + exp(-pn))
    return (pnone, pmech, pburn)
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
# MYHTS(IHAB) → IHTSER habitat-SERIES (1-5) for ESADVH/ESSUBH THAB/UHAB tables (estab.f:111 DATA MYHTS/1,3*2,
# 4*3,2*4,6*5/, line 493 IHTSER=MYHTS(IHAB)). Distinct from ISER=MYHABG(IHAB).
const _IE_MYHTS = Int[1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5]

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
                         pnn::AbstractVector = Float32[], nsp::Integer = 23, wk6fill::Integer = 50,
                         idup::Integer = 0, tally_pt::AbstractMatrix = Array{Float64}(undef, 0, 0),
                         point_slope::AbstractVector = Float32[], point_aspect::AbstractVector = Float32[],
                         prep_sumup = nothing, prob1_prep::AbstractVector = Float32[],
                         is_ie::Bool = false, pasmax::Real = Inf32, prob1_pt::AbstractVector = Float32[],
                         emit::Union{Nothing,Vector{NTuple{5,Float64}}} = nothing,
                         ihtser::Integer = 0, gentim::Real = 5f0, call_espadv::Bool = true,
                         point_baa::AbstractVector = Float32[],
                         over_pt::AbstractMatrix = Array{Float32}(undef, 0, 0),
                         plant_sp::AbstractVector = Int[])
    xc = Float32(xcos); xs = Float32(xsin); sl = Float32(slo); tm = Float32(time)
    # Per-INVENTORY-POINT tally accumulation (optional out-param): plot n belongs to point
    # div(n-1,idup)+1 (same NCOUNT order as the NSTORE fill). Filled when caller supplies a sized
    # tally_pt; the summed `tally` return is byte-unchanged. Lets ie_autoes_establish! place each
    # established tree on its true inventory point (plot_id) instead of hardcoding point 1 — so the
    # NEXT ingrowth tally's per-point NSTORE suppresses ALL stocked points, not just point 1 (#172).
    _fillpt = idup > 0 && size(tally_pt, 1) == nsp && size(tally_pt, 2) >= 1
    _npt = _fillpt ? size(tally_pt, 2) : 1
    _ptof(n) = _fillpt ? min(div(n - 1, Int(idup)) + 1, _npt) : 1
    # Per-IPREP species-mix tables (sumup_base = normalized PADV advance-regen weights; pxcs = excess weights;
    # nspnz = # nonzero best species). Site prep (MECHPREP/BURNPREP) modulates these via the CPRE prep-index term
    # in ie_espadv/ie_espxcs, so each plot uses ITS assigned iprep. Memoized on IPREP (only 1..4 possible); the
    # default (scalar `iprep`) is computed eagerly so the no-prep path stays byte-identical.
    # Cache carries: (sumup_base, pxcs, nspnz, padv_raw, psub_raw). padv_raw = ESPADV weights when the internal
    # NTALLY==1 (ESPADV is CALLed; estab.f:773), else zeros; psub_raw = ESPSUB weights (ALWAYS CALLed, estab.f:774).
    # These raw per-species PADV/PSUB drive the ADV/SUBS ICHOI dispatch (DO 63, estab.f:780-788) for the emit path.
    # Per-INVENTORY-POINT species-mix tables (estab.f:482-484 BAA=BAAA(NNID) clamp[1,400] PER POINT; dense.f:214
    # OVER(ISPC,IP) per point). The oracle runs ESPADV/ESPXCS/ESPSUB per point, so a multi-point stand's low-BA
    # points regen the BA²-penalized species (WP/DF/GF) while its dense point regens redcedar; jl formerly used the
    # scalar (point-1, densest) baa/over for ALL plots ⇒ mono-RC over-retention (RC has the lowest mortality). Key
    # the memo on (POINT, IPREP). ★ ESPADV/ESPXCS/ESPSUB draw NO RNG ⇒ per-point weights leave the wk6/pick stream
    # byte-identical. Fall back to the scalar baa / point-1 `over` when the caller supplies no per-point data
    # (non-IE / disturbance / single point): baa_p[1]==scalar baa and over_pt[:,1]==over by construction, so
    # single-point IE and every non-IE caller stay BYTE-IDENTICAL. Slope stays the scalar sl (unchanged).
    _prep_memo = Dict{Tuple{Int,Int},Tuple{Vector{Float32},Vector{Float32},Int,Vector{Float32},Vector{Float32}}}()
    function _prep_tables(ptn::Integer, ip::Integer)
        get!(_prep_memo, (Int(ptn), Int(ip))) do
            baa_p = (!isempty(point_baa) && ptn <= length(point_baa)) ?
                    clamp(Float32(point_baa[ptn]), 1f0, 400f0) : Float32(baa)
            over_p = (size(over_pt, 1) >= 10 && ptn <= size(over_pt, 2)) ? view(over_pt, :, ptn) : over
            # estab.f:604-609 — the ADV/SUBS/EXCESS species probabilities are ALWAYS computed at TIME=10.0
            # (INGRO=0) or SHORTY (INGRO=1), NOT the tally's elapsed time: `TIME=10.0; IF(INGRO.EQ.1) TIME=SHORTY;
            # CALL ESTIME(IDSDAT,IDSDAT+ITIME)` overwrites TIME right before ESPADV/ESPSUB/ESPXCS. A CONTINUATION
            # tally (NTALLY≥2) has elapsed 20/30 yr but STILL calls ESPADV at TIME=10, so REGT=TIME−BWAF=10 too.
            # jl passed the tally `tm`/`regt` (=next_year−IDSDAT: 10 at cycle-1 but 20+ at a cycle-2 continuation),
            # over-aging the logits (advance −0.228·TIME + −0.118·REGT terms) ⇒ GF adv-prob 0.168 vs the oracle's
            # 0.398 ⇒ DF over-book +53 @2012, compounding to +212 TPA @2032. Use the estab.f:606 species-prob
            # TIME/REGT here (10 / SHORTY); ESTOCK/PROB1 (ie_autoes_run) keeps the ACTUAL elapsed time. Cycle-1
            # (tm=regt=10) is unchanged ⇒ byte-identical; only the continuation (tm>10) is corrected.
            tm_sp = is_ingro ? tm : 10f0
            rg_sp = is_ingro ? Float32(regt) : 10f0
            pa = collect(ie_espadv(ihab, ip, ifo, iphy, xc, xs, sl, tm_sp, baa_p, Float32(elev),
                                   rg_sp, Float32(bwaf), Float32(bwb4), occ, over_p))
            px = collect(ie_espxcs(ihab, ip, ifo, iphy, xc, xs, sl, tm_sp, baa_p, Float32(elev),
                                   rg_sp, Float32(bwaf), Float32(bwb4), occ, over_p))
            ps_full = collect(ie_espsub(ihab, ip, ifo, iphy, xc, xs, sl, tm_sp, baa_p, log(baa_p),
                                        Float32(elev), sqrt(rg_sp), sqrt(Float32(bwaf)), Float32(bwb4),
                                        occ, over_p))
            padv_raw = call_espadv ? copy(pa) : zeros(Float32, 10)
            # BEST-species SUMUP = normalized (PADV + PSUB) (estab.f:726-735 FTEMP=PADV(I)+PSUB(I)). PSUB is
            # added ONLY when ITIME=INT(TIME+0.5)>2 (estab.f:606-616): the DATED/DISTURBANCE tally uses TIME=10
            # (ITIME=10>2 ⇒ ESPSUB active), while the AUTOES ingrowth tally uses TIME=SHORTY=1 (ITIME=1≤2 ⇒
            # PSUB stays 0). Omitting PSUB shifted the disturbance-tally best-species mix toward the PADV-dominant
            # species — MEASURED on FVSie_g16 na_def (bare IE dated ESTAB): GF over-booked (jl 436 vs 348), LP
            # under-booked (jl 4 vs 47); per-plot picks diverged exactly where the 2nd best species is a tail
            # species (plots 8→WL, 10→LP). SQREGT/SQBWAF: estime.f SQREGT=√TIME−SQBWAF, REGT=TIME−BWAF; with no
            # WSBW (BWAF=SQBWAF=0) √regt=√TIME=SQREGT — mirrors the ie_estock call. The ingrowth path (tm≤2 ⇒
            # ps=0) stays BYTE-IDENTICAL, so the validated AUTOES-ingrowth tally is unperturbed.
            ps = (is_ie && round(Int, tm_sp + 0.5f0) > 2) ?
                 collect(ie_espsub(ihab, ip, ifo, iphy, xc, xs, sl, tm_sp, baa_p, log(baa_p),
                                   Float32(elev), sqrt(rg_sp), sqrt(Float32(bwaf)), Float32(bwb4),
                                   occ, over_p)) :
                 zeros(Float32, 10)
            sb = zeros(Float32, nsp); sb[1:10] .= pa .+ ps; sb ./= sum(sb)
            (sb, px, count(>(1f-4), sb), padv_raw, ps_full)
        end
    end
    sumup_base, pxcs, nspnz, padv_raw, psub_raw = _prep_tables(1, iprep)   # default (point 1); loop recomputes per point
    maxspp = _IE_MAXSPP[ihab]
    # Per-plot IPPREP (estab.f:382-399): sample-without-replacement from the WK6 site-prep vector when a
    # MECHPREP/BURNPREP keyword supplied prep_sumup (disturbance tally only; ingrowth path forces IPREP=1).
    prep_active = prep_sumup !== nothing && !is_ingro
    ipprep = Int[]
    # ITPP cap (estab.f:681-682): ALWAYS MAXTPP; the MAXING cap applies ONLY when INGRO=1. Instrument-replay
    # (FVSie iet01) CONFIRMED jl's ESTPP draws + TPP are BIT-IDENTICAL to live (0.346302→2.486, 0.835307→14.036,
    # 0.967520→34.5) and jl prob1 (0.60012) matches live PROB1 (0.60122) ⇒ per-tree TPA correct. The sole divergence
    # was jl capping the disturbance path at MAXING(7) too — truncating ITPP 14→7, 25→7 ⇒ UNDER-production. With the
    # split cap jl ITPP = live bit-exact [2,1,14,2,4,2,25,4]. (Old "over-produces at MAXTPP" comment predated the
    # ESTOCK-PROB1 understanding.) NOTE: the ingrowth (is_ingro) dense-stand NSTORE over-book (BA>400 stands) was
    # a MISSING BAAOLD [1,400] clamp on the ESB1 stocking prediction — FIXED in ie_autoes_establish! (see there).
    cap = is_ingro ? _IE_MAXING[ihab] : _IE_MAXTPP[ihab]
    p1s = Float32(prob1); scale = 300f0 / Float32(dupnpt)   # stand-scalar PROB1 (per-point override below)
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
    # PASSALL/PASMAX excess-tree cap (estab.f:1079-1145 NOTE selection + 1288-1370 excess pass). When `pasmax`
    # is finite (only the IE seam passes it; every other caller/variant leaves it Inf ⇒ this block is dead ⇒
    # byte-identical), each plot's per-species EXCESS regen count is capped at PASMAX via es_pasmax_xcsmax. The
    # per-tree booking below is left UNCHANGED (so the uncapped case stays byte-identical); only species whose
    # EXCESS strictly exceeds PASMAX get a deterministic reduction Δ = capped − uncapped applied afterwards.
    docap = isfinite(Float32(pasmax)); pmx = Float32(pasmax)
    exc_cnt = docap ? zeros(Float32, nsp) : Float32[]        # per-species EXCESS count (NOTE=0 trees)
    exc_sum = docap ? zeros(Float32, nsp) : Float32[]        # per-species SUMESP (Σ esprob of those excess trees)
    # ── Establishment-cohort HEIGHT-CLASS / WK4 emission (estab.f DO 99 + DO 33/228). When `emit` is supplied
    # (IE only), the per-plot tally additionally reproduces FVS's per-tree record booking: each best-species
    # first tree (advance ESADVH WK4≈0.60 / subsequent ESSUBH WK4≈0.20/0.00) is one record; the leftover trees
    # accumulate per species into ceil(EXCESS/5) tripled excess records (mean ESXCSH height, WK4=STOMLT). Each
    # record carries (sp, point, height, wk4, tpa). Collapsing to one WK4=0.60 record over-projected the cohort.
    doemit = emit !== nothing
    gtim = Float32(gentim); ihts = Int(ihtser); iphy_i = Int(iphy)
    xmin_e = _IE_ES_XMIN; hhtmax_e = _IE_ES_HHTMAX; bnorml_e = _IE_ES_BNORML
    first1 = doemit ? fill(0.1f0, nsp) : Float32[]           # FIRST(1,sp) advance dilate order-statistic (persists
    first2 = doemit ? fill(0.1f0, nsp) : Float32[]           # across plots within a tally, estab.f:179-183)
    stomlt = doemit ? fill(1f0, nsp) : Float32[]             # per-species WK4 for THIS plot (reset per plot)
    tallh  = doemit ? zeros(Float32, nsp) : Float32[]        # per-species best-tree height (post-floor)
    iasep_e = doemit ? zeros(Int, nsp) : Int[]               # 1=advance 2=subsequent
    # ── Scheduled PLANT/NATURAL trees co-located on EACH plot (estab.f:970-1073 DO 322): each due PLANT/NATURAL
    # activity adds exactly ONE tree per plot at position ITPP+ITODO with species IPNSPE=PRMS(1) (ITODO=#activities,
    # uniform across plots — NOT a TPA distribution). FVS then runs the NBEST best/excess selection (DO 166/168/172)
    # over the COMBINED ITP=ITPP+ITODO set, so a planted NEW-species tree can take a best slot (STEP2, estab.f:1168),
    # pushing a natural tree into EXCESS and tripping the XCSMAX cap (estab.f:1320). jl books PLANT via establish!
    # (separately), so here the planted trees enter ONLY the NBEST ranking as non-booked "phantoms". The natural
    # per-species EXCESS COUNT (what drives XCSMAX) is invariant to the phantom's height — each distinct planted
    # species consumes exactly one best slot whether it is STEP1-tall or STEP2-short (STEP3 then fills the rest with
    # the tallest remaining natural either way), AND estab.f's DO 33 excess accumulation loops only the naturals
    # (N=1,ITPP) — the planted tree at ITPP+ITODO is booked by establish!, not counted as natural excess. So the
    # phantom height is set to 0: it is EXCLUDED from STEP1/STEP3 (both use the FTEMP=0.001 tallest-search floor,
    # and 0<0.001), so it never displaces a natural there; but STEP2 (TALL(J) init 0) still promotes it as the sole
    # representative of a NEW planted species — exactly the slot that pushes a natural species into EXCESS and can
    # trip XCSMAX. When the planted species is ALSO a natural regen species on the plot, the taller natural wins
    # STEP2 and the phantom is inert (matching the oracle, whose planted tree is likewise not among the naturals'
    # best there). Empty plant_sp (no PLANT/NATURAL this tally) ⇒ _nph=0 ⇒ ITP==ITPP ⇒ byte-identical.
    _ph_sp = doemit ? Int[Int(round(x)) for x in plant_sp if 1 <= Int(round(x)) <= nsp] : Int[]
    _nph = length(_ph_sp)
    _ph_ht = zeros(Float32, _nph)
    for (n, sd) in enumerate(seeds)
        rng = IEEstabRNG(sd)
        if n == 1
            if prep_active                                                   # capture the DUPNPT WK6 site-prep draws
                wk6v = Float32[ie_esrann!(rng) for _ in 1:wk6fill]           # (estab.f:333 DO 183) and sample IPPREP
                _id = idup > 0 ? Int(idup) : 1
                ipprep = ie_esetpr_sample(prep_sumup, wk6v, div(nplots, _id), _id)  # nptids = dupnpt/idup
            else
                for _ in 1:wk6fill; ie_esrann!(rng); end
            end
        end
        _emd1 = ie_esrann!(rng); _emd2 = ie_esrann!(rng)                     # EMSQR: sign@1 · magnitude@2
        emsqr = (_emd1 < 0.5f0 ? -1f0 : 1f0) * _emd2                         # estab.f:646-650
        # Per-INVENTORY-POINT slope/aspect for ESTPP (live SLO=PSLO(NNID), XCOS=cos(PASP)·PSLO). Plot n → point
        # div(n-1,idup)+1. When per-point topo is supplied (FIA per-plot SLOPE/ASPECT) each plot's ESTPP uses ITS
        # point's slope; else fall back to the uniform (stand/TREEDATA) xc/xs/sl. FIXES #143 (jl formerly used the
        # stand slope for all plots ⇒ over-suppressed ESTPP on sloped stands). Does NOT draw ⇒ RNG order unchanged.
        local _sln::Float32, _xcn::Float32, _xsn::Float32
        _ptn = idup > 0 ? div(n - 1, Int(idup)) + 1 : 1
        if !isempty(point_slope) && _ptn <= length(point_slope)
            _sln = point_slope[_ptn]                                         # tree-bearing point: its PSLO
            _aspn = _ptn <= length(point_aspect) ? point_aspect[_ptn] : 0f0
            _xcn = cos(_aspn) * _sln; _xsn = sin(_aspn) * _sln
        else
            _sln = sl; _xcn = xc; _xsn = xs                                  # empty point (no tree) → stand slope
        end
        # This plot's species-mix tables use ITS inventory point's BAA + per-species OVERSTORY BA (estab.f runs
        # ESPADV/ESPXCS/ESPSUB per point) and ITS assigned IPREP. Deterministic (no RNG) ⇒ draw stream unchanged;
        # single-point / non-IE fall back to the point-1/scalar tables ⇒ byte-identical (see _prep_tables).
        sumup_base, pxcs, nspnz, padv_raw, psub_raw =
            _prep_tables(_ptn, (prep_active && n <= length(ipprep)) ? ipprep[n] : iprep)
        itpp = clamp(round(Int, ie_estpp(ie_esrann!(rng), ihab, _xcn, _xsn, _sln, Float32(regt), Float32(bwaf))), 1, cap)
        p1 = isempty(prob1_pt) ? p1s : prob1_pt[clamp(_ptn, 1, length(prob1_pt))]   # this plot's inventory-point PROB1
        ns = has_state ? Int(nstore[n]) : 0
        pn = has_state ? pnn[n] : 0f0
        newtpp = max(0, itpp - ns)
        # Per-plot PROB1 (estab.f:572-584): on the DISTURBANCE tally the stocking logit ie_estock carries a
        # per-IPREP SPRE term (estock.f), so each plot's PROB1 depends on its sampled IPPREP. `prob1_prep`
        # (=[PROB1(IPREP=1..3)], supplied only when the default/keyword prep is active) selects it; otherwise
        # the stand-level scalar p1 (IPREP=1). MEASURED FVSie_g16 na_def: PROB1 0.699/0.663977/0.657729 by IPREP.
        # ★ When the per-point prob1_pt is present (the ESB1 correction is active, esb_shift≠0) the SPRE term
        # CANCELS in logistic(PN+ESB−ESB1) — it is time/BA-independent and enters BOTH the end-of-cycle PN and the
        # inventory ESB1 — so PROB1 is per-point AND iprep-independent (the oracle prints one PROB1 per point). Use
        # prob1_pt then and SKIP the per-IPREP override (which would mis-apply point-1's SPRE-uncancelled prob1_prep).
        p1n = (prep_active && !isempty(prob1_prep) && isempty(prob1_pt) && n <= length(ipprep)) ?
              Float32(prob1_prep[ipprep[n]]) : p1
        # ESPROB (estab.f:944-951): a tree at plot-index I gets full PROB1 if new (I>NSTORE); an old tree
        # (I≤NSTORE) gets the increment PROB1-PNN; an ingrowth tally scales ALL trees by NEWTPP/ITPP.
        prob_old = max(p1n - pn, 0.0001f0)
        esprob(i) = is_ingro ? max(p1n * Float32(newtpp) / Float32(itpp), 0.0001f0) : (i <= ns ? prob_old : p1n)
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
            iplot += 1; j = ie_estab_pick_species(wk6s[i], su); _c = esprob(iplot) * scale
            tally[j] += _c; _fillpt && (tally_pt[j, _ptof(n)] += _c); ibest[j] = 1
            su[j] = 0f0; t = sum(su); t > 0 && (su ./= t)
        end
        we = zeros(Float32, nsp)
        for i in 1:10
            we[i] = pxcs[i] * ibest[i]; (ibest[i] == 1 && we[i] < 0.0001f0) && (we[i] = 0.0001f0)
        end
        twe = sum(we); twe > 0 && (we ./= twe)
        # ADV/SUBS + heights (estab.f DO 63 ICHOI[nsp] + DO 122 WK6[2·nsp]). Draw them in-order (identical RNG
        # consumption to the prior bulk discard), then — for the emit path — reproduce ESADVH/ESSUBH per best
        # species: ICHOI (advance if draw ≤ PADV/(PADV+PSUB)) → ESDLAY delay → ESADVH/ESSUBH height + STOMLT.
        ichoi_dr = Vector{Float32}(undef, nsp)
        @inbounds for i in 1:nsp; ichoi_dr[i] = ie_esrann!(rng); end          # DO 63 (one draw per species; used if IBEST)
        wk6d = Vector{Float32}(undef, 2 * nsp)
        @inbounds for i in 1:(2 * nsp); wk6d[i] = ie_esrann!(rng); end        # DO 122 WK6 (ESDLAY draws)
        if doemit
            ip_plot = (prep_active && n <= length(ipprep)) ? Int(ipprep[n]) : Int(iprep)
            @inbounds for i in 1:nsp; stomlt[i] = 1f0; tallh[i] = 0.001f0; iasep_e[i] = 0; end
            ndraw = 0
            @inbounds for sp2 in 1:nsp
                ndraw += 1                                                    # estab.f:803 (every species advances NDRAW)
                ibest[sp2] == 1 || continue                                   # best species are always 1..10
                pav = sp2 <= 10 ? padv_raw[sp2] : 0f0
                psu = sp2 <= 10 ? psub_raw[sp2] : 0f0
                sm = pav + psu
                ftmp = sm > 0f0 ? pav / sm : 0f0                             # PADV/SUM
                isadv = ichoi_dr[sp2] <= ftmp                               # J=1 advance if DRAW ≤ FTEMP (estab.f:785-786)
                local hh::Float32, trage::Float32
                if isadv
                    drw = wk6d[ndraw]
                    delay = ie_esdlay(sp2, 1, drw, tm, baa; bwb4 = bwb4, bwaf = bwaf)
                    nN = trunc(Int, delay + 0.5f0); nN > 2 && (nN = 1)        # esadvh.f: N=INT(DELAY+.5); N>2→1
                    dN = Float32(nN)
                    trage = 3f0 - dN
                    agev = 3f0 - dN - gtim; agev < 1f0 && (agev = 1f0)
                    itime = trunc(Int, tm + 0.5f0); itime < 1 && (itime = 1); itime > 20 && (itime = 20)
                    bnrm = bnorml_e[itime]
                    dil = first1[sp2]; first1[sp2] = sqrt(dil)
                    hh = ie_esadvh(sp2, emsqr, dil, log(agev), bnrm; baa = baa, elev = elev, xcos = xc,
                                   xsin = xs, slo = sl, ihtser = ihts, iphy = iphy_i, iprep = ip_plot,
                                   bwaf = bwaf, bwb4 = bwb4)
                    iasep_e[sp2] = 1
                else
                    ndraw += 1                                                # estab.f:823 (subsequent consumes a 2nd slot)
                    drw = wk6d[ndraw]
                    delay = ie_esdlay(sp2, 2, drw, tm, baa; bwb4 = bwb4, bwaf = bwaf)
                    nN = trunc(Int, delay + 0.5f0); nN < -3 && (nN = -3)      # essubh.f: N=INT(DELAY+.5); N<-3→-3
                    itime = trunc(Int, tm + 0.5f0)
                    dN = (nN > itime) ? tm : Float32(nN)                      # IF(N>ITIME) DELAY=TIME
                    ageA = tm - dN - gtim
                    iage = trunc(Int, ageA + 0.5f0); iage < 1 && (iage = 1); iage > 20 && (iage = 20)
                    agev = ageA < 1f0 ? 1f0 : ageA                           # AGE=AGE+TRAGE(=0); AGE<1→1
                    bnrm = bnorml_e[iage]
                    dil = first2[sp2]; first2[sp2] = sqrt(dil)
                    disp = emsqr * dil * bnrm
                    hh = ie_essubh(sp2, agev, baa, ihts, ip_plot, iphy_i, xc, xs, sl, elev, disp;
                                   bwaf = bwaf, bwb4 = bwb4)
                    trage = tm - dN
                    iasep_e[sp2] = 2
                end
                ft = trage; ft > gtim && (ft = gtim); stomlt[sp2] = ft / (gtim + 0.0001f0)   # STOMLT=min(TRAGE,GENTIM)/…
                tv = hh                                                       # +HTADJ (=0); floor XMIN+0.2; cap HHTMAX
                (tv - xmin_e[sp2] < 0.2f0) && (tv = xmin_e[sp2] + 0.2f0)
                (tv > hhtmax_e[sp2]) && (tv = hhtmax_e[sp2])
                tallh[sp2] = tv
            end
        end
        wk6e = ntuple(_ -> ie_esrann!(rng), excess_draws)                    # excess-WK6 (2·MAXTPP[ihab])
        # Per-plot tree records for the emit path: positions 1..numspe = best-species first trees (numeric order,
        # estab.f DO 146), numspe+1..itpp = excess trees (ESXCSH, filled in the excess loop below).
        e_sp = doemit ? zeros(Int, itpp) : Int[]
        e_ht = doemit ? zeros(Float32, itpp) : Float32[]
        e_wk4 = doemit ? zeros(Float32, itpp) : Float32[]
        e_esp = doemit ? zeros(Float32, itpp) : Float32[]
        if doemit
            _N = 0
            @inbounds for sp2 in 1:nsp
                ibest[sp2] == 1 || continue
                _N += 1; _N > itpp && break
                e_sp[_N] = sp2; e_ht[_N] = tallh[sp2]; e_wk4[_N] = stomlt[sp2]; e_esp[_N] = Float32(esprob(_N))
            end
        end
        # NOTE best/excess split (estab.f:1079-1145). MEASURED bit-exact rule (FVSie estab.f dump, 54/54 plots):
        # NBEST = min(ITP, max(4, #distinct species)) and #distinct == NUMSPE (excess trees pick only among best
        # species). The NBEST trees = the NUMSPE best-species trees (each its species' tallest) + the tallest
        # (NBEST−NUMSPE) excess trees. jl has no per-tree ESXCSH height (the FIRST order-statistic TALL chain is
        # unported), so the promoted excess are taken in generation order — EXACT for single-species-excess plots
        # (the cap is species-independent there) and cornered only for the ~6% multi-species-excess plots (the
        # step-1/3 cross-species height ranking). Inert unless docap.
        npro = docap ? (min(itpp, max(4, numspe)) - numspe) : 0             # excess trees promoted to "best"
        touched = Int[]
        nd = 0; kx = 0
        for _ in 1:(itpp - numspe)
            nd += 1; iplot += 1; kx += 1
            j = ie_estab_pick_species(wk6e[nd], we); _e = esprob(iplot); _c = _e * scale
            tally[j] += _c; _fillpt && (tally_pt[j, _ptof(n)] += _c)
            if doemit                                                      # ESXCSH excess-tree height (estab.f DO 156)
                dhx = wk6e[nd + 1]                                         # 2nd draw = ESXCSH DRAW (NDRAW+1)
                hx = ie_esxcsh(j, tallh[j], xmin_e[j], tm, dhx)           # HTMAX=TALL(II), HTMIN=XMIN(II)
                (hx < xmin_e[j]) && (hx = xmin_e[j]); (hx > hhtmax_e[j]) && (hx = hhtmax_e[j])
                e_sp[iplot] = j; e_ht[iplot] = hx; e_wk4[iplot] = stomlt[j]; e_esp[iplot] = Float32(_e)
            end
            nd += 1
            if docap && kx > npro                                          # a NOTE=0 (true-excess) tree
                (exc_cnt[j] == 0f0) && push!(touched, j)
                exc_cnt[j] += 1f0; exc_sum[j] += _e
            end
        end
        if docap
            pt = _ptof(n)
            for j in touched
                ex = exc_cnt[j]
                if ex > pmx                                                # cap bites (uncapped ⇒ untouched ⇒ byte-identical)
                    ft2 = exc_sum[j] / (ex + 1f-5)                         # FTEMP2 = SUMESP/(EXCESS+1e-5)
                    xcs, ibrk = es_pasmax_xcsmax(ex, pmx)                  # estab.f:1288/1318-1321
                    dlt = Float64(ibrk) * Float64(ft2 * xcs) * Float64(scale) - Float64(exc_sum[j]) * Float64(scale)
                    tally[j] += dlt; _fillpt && (tally_pt[j, pt] += dlt)
                end
                exc_cnt[j] = 0f0; exc_sum[j] = 0f0                         # reset for the next plot
            end
        end
        # ── Best-tree selection (estab.f:1090-1144 DO 166/168/172) + record booking (DO 33 best, DO 228 excess).
        if doemit && itpp >= 1
            pt_e = _ptof(n)
            # Combined pool = natural regen (1..itpp) + co-located PLANT/NATURAL phantoms (itpp+1..itp), matching
            # estab.f's ITP=ITPP+ITODO. NBEST (DO 166/168/172) ranks over the combined set; only the naturals are
            # booked below (phantoms are booked by establish!). _nph=0 ⇒ itp==itpp ⇒ byte-identical to the pre-fix path.
            itp = itpp + _nph
            c_sp = Vector{Int}(undef, itp); c_ht = Vector{Float32}(undef, itp)
            @inbounds for i in 1:itpp; c_sp[i] = e_sp[i]; c_ht[i] = e_ht[i]; end
            @inbounds for k in 1:_nph; c_sp[itpp + k] = _ph_sp[k]; c_ht[itpp + k] = _ph_ht[k]; end
            note_e = zeros(Int, itp)
            if p1n >= 0.00011f0                                              # ISTART=1 unless PROB1<0.00011 (estab.f:1088)
                null_sp = zeros(Int, nsp); nbest = 0
                for _step in 1:2                                             # STEP 1: 2 tallest regardless of species
                    nbest >= itp && break
                    ftv = 0.001f0; itemp = 0
                    @inbounds for i in 1:itp
                        note_e[i] == 1 && continue
                        c_ht[i] < ftv && continue
                        ftv = c_ht[i]; itemp = i
                    end
                    itemp == 0 && break
                    (c_sp[itemp] >= 1) && (null_sp[c_sp[itemp]] = 1); note_e[itemp] = 1; nbest += 1
                end
                talls = zeros(Float32, nsp); ilsp = zeros(Int, nsp)         # STEP 2: tallest of each additional species
                @inbounds for i in 1:itp
                    note_e[i] == 1 && continue
                    s = c_sp[i]; (s < 1) && continue
                    c_ht[i] < talls[s] && continue
                    talls[s] = c_ht[i]; ilsp[s] = i
                end
                @inbounds for s in 1:nsp
                    (null_sp[s] == 1 || ilsp[s] == 0) && continue
                    note_e[ilsp[s]] = 1; nbest += 1
                end
                if nbest < 4 && nbest < itp                                  # STEP 3: tallest remaining until ≥4
                    while true
                        ftv = 0.001f0; itemp = 0
                        @inbounds for i in 1:itp
                            note_e[i] == 1 && continue
                            c_ht[i] < ftv && continue
                            ftv = c_ht[i]; itemp = i
                        end
                        itemp == 0 && break
                        nbest += 1; note_e[itemp] = 1
                        nbest >= itp && break
                        nbest < 4 || break
                    end
                end
            end
            ex_c = zeros(Float32, nsp); ex_h = zeros(Float32, nsp); ex_p = zeros(Float32, nsp)
            @inbounds for N in 1:itpp                                        # DO 33: best trees are individual records
                I = e_sp[N]; (I < 1) && continue; hh = e_ht[N]
                if e_esp[N] >= 0.00011f0 && note_e[N] == 1
                    push!(emit, (Float64(I), Float64(pt_e), Float64(hh), Float64(e_wk4[N]),
                                 Float64(e_esp[N]) * 300.0 / Float64(dupnpt)))
                else                                                        # accumulate excess per species (DO 33:199)
                    ex_c[I] += 1f0; ex_h[I] += hh; ex_p[I] += e_esp[N]
                end
            end
            @inbounds for I in 1:nsp                                         # DO 228: excess → ceil(EXCESS/5) records
                ex_c[I] < 0.5f0 && continue
                if ex_p[I] < 0.00011f0 && I < nsp                            # roll a near-zero-prob species forward
                    ex_p[I+1] += ex_p[I]; ex_p[I] = 0f0; continue
                end
                ft2 = ex_p[I] / (ex_c[I] + 0.00001f0)                        # FTEMP2 = SUMESP/(EXCESS+1e-5)
                xcs, ibrk = es_pasmax_xcsmax(ex_c[I], pmx)                   # XCSMAX + IBRKUP (estab.f:1288/1318-1321)
                hh = ex_h[I] / ex_c[I]                                       # mean excess height (SUMHTS/EXCESS)
                tpa_r = Float64(ft2) * 300.0 * Float64(xcs) / Float64(dupnpt)
                for _ib in 1:ibrk
                    push!(emit, (Float64(I), Float64(pt_e), Float64(hh), Float64(stomlt[I]), tpa_r))
                end
            end
        end
        has_state && (nstore[n] = Int32(itpp); pnn[n] = p1n)                 # carry to the next tally
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
                             inv_year::Integer, npnats::Integer=0)
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
    # (4) PLANT/NATURAL-in-case catch-all (esnutr.f:345-359): a PLANT(430)/NATURAL(431) scheduled THIS cycle forces
    # an ingrowth ESTAB call (NTALLY=99, IDSDAT=IY(ICYC+1)-20) when no earlier rule fired — so a PLANTED stand gets
    # AUTOMATIC natural regen IN the plant cycle (OPFIND(2,[430,431]) = current-cycle PLANT/NATURAL count). Unlike
    # rule 3 this block leaves LONE=.FALSE. (esnutr.f:51 default, never re-set here), so NTALLY stays 99 after ESTAB
    # (the subsequent gap checks still gate the next fire). A NATURAL stand never reaches here — est.lingrw is false
    # (esin.f:1300) so ie_autoes_establish! returns at its 1708 (lautal||lingrw) gate before the scheduler runs.
    if npnats > 0
        est.idsdat = Int32(next_year - 20)                # IDSDAT = IY(ICYC+1) - 20
        est.ntally = Int32(99)                            # LONE=.FALSE. ⇒ NTALLY NOT reset (esnutr.f:416)
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

# ie/habtyp.f plant-association string crosswalk (Colville R6). A 6-char alphanumeric PV_CODE (e.g. "CDS715") is
# matched in PCOML → index i → the NI habitat code = MTYPE[KTYPE[MAPR6[i]]]. (habtyp.f maps NITYPE=MAPR6[i] to
# KODTYP=JTYPE[NITYPE], then the general lookup `first K: KODTYP<JTYPE[K]` gives ITYPE=KTYPE[K-1]=KTYPE[NITYPE]
# since JTYPE (blkdat.f) is strictly ascending; final KODTYP=MTYPE[ITYPE].) ★#143: jl's FIA reader parsed the
# STRING PV_CODE as a number (→0) and fell back to the numeric PV_REF_CODE, selecting the WRONG ESTOCK habitat
# series (ihab=13/subalpine vs live's 3/DF) ⇒ AUTOES over-established up to +108%. This resolves it (CDS715→260).
const _IE_HABTYP_PCOML = ("CCF221","CCF222","CCS311","CDG131","CDG311","CDS632","CDS633","CDS637","CDS715","CDS716",
    "CDS813","CDS814","CEF111","CEF211","CEF421","CEF422","CEF423","CEG311","CEM211","CES210","CES211","CES312",
    "CES313","CES412","CES422","CHF311","CHF312","CHF422","CHF521","CHS711","CLS521","CWF411","CWS214","CWS421",
    "CWS422","CWS821","CCS211","CHS411","CAG112","CDG123")
const _IE_HABTYP_MAPR6 = Int32[45,45,45,19,4,25,18,18,15,17,14,16,68,63,58,60,58,75,27,64,64,72,72,73,61,52,52,
    51,52,54,88,43,43,39,43,55,48,51,84,22]
const _IE_HABTYP_KTYPE = Int32[1,1,1,1,1,2,2,2,2,4,1,1,1,3,4,5,6,7,8,9,8,8,9,7,3,10,10,10,4,11,20,29,11,11,13,14,
    17,12,12,12,12,13,13,13,14,15,14,16,14,16,17,17,17,24,12,24,18,19,21,19,20,20,21,22,19,23,19,24,27,25,25,26,
    27,22,24,27,24,27,28,28,29,28,28,29,29,29,29,27,9,20,24,21,27,24,30]
const _IE_HABTYP_MTYPE = Int32[130,170,250,260,280,290,310,320,330,420,470,510,520,530,540,550,570,610,620,640,
    660,670,680,690,710,720,730,830,850,999]

# ie_pa_habitat_code(code) → NI habitat code (KODTYP), or 0 if `code` isn't a recognized 6-char plant-association
# PV_CODE. Maps e.g. "CDS715"→260 via ie/habtyp.f's PCOML→MAPR6→KTYPE→MTYPE crosswalk. Caller keeps its numeric
# fallback on 0.
function ie_pa_habitat_code(code::AbstractString)::Int
    s = strip(uppercase(code))
    length(s) == 6 || return 0
    idx = findfirst(==(s), _IE_HABTYP_PCOML)
    idx === nothing && return 0
    nitype = Int(_IE_HABTYP_MAPR6[idx])                 # NITYPE = MAPR6[i]
    (1 <= nitype <= length(_IE_HABTYP_KTYPE)) || return 0
    itype = Int(_IE_HABTYP_KTYPE[nitype])               # ITYPE = KTYPE[NITYPE] (JTYPE strictly ascending)
    (1 <= itype <= length(_IE_HABTYP_MTYPE)) || return 0
    return Int(_IE_HABTYP_MTYPE[itype])                 # final KODTYP = MTYPE[ITYPE]
end

# ie_pa_ni_code(code) → the NORTH-IDAHO habitat code ICL5 = JTYPE(NITYPE) for a 6-char plant-association PV_CODE,
# or 0 if unrecognized. ★ This is what ie/habtyp.f STORES AS ICL5 (line 124-125: KODTYP=JTYPE(NITYPE); ICL5=KODTYP)
# BEFORE the final MTYPE(ITYPE) remap — and ICL5, not the MTYPE growth code, is what esplt2.f/estab.f key the
# AUTOES habitat bracket off. (ie_pa_habitat_code above returns the MTYPE growth code, which HBDECD→ITYPE→MTYPE
# yields; for a 6-char code the two DIFFER, e.g. CWS821 → NITYPE 55 → JTYPE(55)=590 (estab ICL5) vs MTYPE 510
# (growth). Feeding 510 to the estab bracket gives group 6 (SHAB −0.06); the oracle's 590 gives group 5
# (SHAB +0.54) — MEASURED live FVSie_g16 645170890126144: ICL5=590, ESTOCK IHAB=5.) Growth is invariant to which
# we store: ie_habtyp(JTYPE(NITYPE)) == ie_habtyp(MTYPE[ITYPE]) == ITYPE (MTYPE codes are the canonical JTYPE
# representative of each group), so the DG/site path re-derives the identical ITYPE either way.
function ie_pa_ni_code(code::AbstractString)::Int
    s = strip(uppercase(code))
    length(s) == 6 || return 0
    idx = findfirst(==(s), _IE_HABTYP_PCOML)
    idx === nothing && return 0
    nitype = Int(_IE_HABTYP_MAPR6[idx])                 # NITYPE = MAPR6[i]
    (1 <= nitype <= length(IE_JTYPE)) || return 0
    return Int(IE_JTYPE[nitype])                        # ICL5 = JTYPE(NITYPE)
end

# ie/pvref1.f: the (PV_CODE, PV_REF_CODE) → HABPVR crosswalk (879 rows, first-match-wins per live's DO-EXIT).
# Live habtyp.f calls PVREF1 whenever a reference code is present; a FULL match (both PVCODE and PVREF hit the
# same row) sets KODTYP=HABPVR; ANY non-full-match (incl. neither found, e.g. PV_CODE "ABR8" + ref 639) leaves
# KODTYP invalid ⇒ habtyp's default ITYPE 4 ⇒ MTYPE(4)=260. #143 ROOT: jl's FIA reader had no PVREF1 and fell
# back to the RAW PV_REF_CODE (639) as the habitat ⇒ ESTOCK ihab 11 (GF-dominant) vs live's ihab 3 (DF/PP);
# 4/5 IE sweep stands are "NOT RECOGNIZED → 260". HABPVR is either a numeric habitat code or a 6-char PCOML
# plant-community code (routed through ie_pa_habitat_code). Built once from the generated data block.
const _IE_PVREF1 = let d = Dict{Tuple{String,String},String}()
    for ln in eachline(IOBuffer(_IE_PVREF1_RAW))
        isempty(strip(ln)) && continue
        f = split(ln)
        length(f) == 3 || continue
        k = (f[1], f[2])
        haskey(d, k) || (d[k] = f[3])        # first-match-wins (live iterates 1..N, EXIT on first)
    end
    d
end

"""
    ie_pvref1(pv_code, pv_ref) -> Int

Resolve the IE FIA habitat KODTYP from a (PV_CODE, PV_REF_CODE) pair via the ported ie/pvref1.f crosswalk.
Returns the HABPVR habitat code on a full match (numeric directly, or a PCOML string mapped through
[`ie_pa_habitat_code`](@ref)), or `0` if the pair is not in the table (caller then applies the live default 260).
"""
function ie_pvref1(pv_code::AbstractString, pv_ref::Integer)::Int
    # CASE-SENSITIVE match: ie/pvref1.f compares ADJUSTL(PVCODE(I)).EQ.ADJUSTL(KARD2T) — it left-justifies but does
    # NOT uppercase, and the table codes are UPPERCASE. So a lowercase DB PV_CODE (e.g. "ces211") does NOT match
    # "CES211" ⇒ no PVREF1 hit ⇒ habtyp.f defaults to 260. (Only 3 IE stands are lowercase: ccf221/cef111/ces211;
    # jl formerly uppercased ⇒ matched ⇒ resolved a habitat the oracle defaults to 260.) Uppercase codes unchanged.
    # ie/habtyp.f:84-96 — BEFORE the PVREF1 lookup the habitat key KARD2 is left-justified and, when it is a
    # ≤2-char PURELY-NUMERIC code, given a leading zero to 3 chars (`IF(LEN_TRIM(KARD2).LE.2 .AND. both chars
    # digits) KARD2='0'//KARD2`). The FIA DB delivers 2-digit PV_CODEs unpadded (e.g. "31"), but the ported
    # ie/pvref1.f table keys are 3-char zero-padded ("031"). Without this pad ie_pvref1("31",110) MISSES the
    # ("031","110")→130 row ⇒ habtyp defaults to 260 ⇒ ITYPE 4 (MTYPE(4)) instead of the correct 130→ITYPE 1,
    # which mis-selects the RHCON habitat term for NIVAR species (e.g. sp10 RHHAB(3)=-0.4345 vs the correct
    # RHHAB(1)=-0.2146, a +0.22 CON deficit ⇒ ~26% under-height ⇒ compounding QMD/BA deficit on bare/pure-
    # establishment IE stands). Mirror habtyp.f: pad a length-2 all-digit code to "0"//code.
    pvc = String(strip(pv_code))
    (length(pvc) == 2 && all(isdigit, pvc)) && (pvc = "0" * pvc)
    hab = get(_IE_PVREF1, (pvc, string(Int(pv_ref))), "")
    isempty(hab) && return 0
    all(isdigit, hab) && return parse(Int, hab)          # numeric HABPVR = habitat code directly
    return ie_pa_habitat_code(hab)                        # PCOML HABPVR (e.g. "CDS715") → NI code, 0 if unmapped
end

"""
    ie_estab_kodtyp(kodtyp) -> Int

Apply ie/habtyp.f's FINAL habitat crosswalk (KODTYP → MTYPE(ITYPE) via the JTYPE/KTYPE lookup,
`ie_habtyp` gives ITYPE) to a raw/pvref1-intermediate habitat code, yielding the FVS "mapped"
habitat code the oracle stores as ICL5 and feeds to esplt2.f/estab.f. ★ ROOT of the IE AUTOES
over-regen population bug: jl's FIA reader stopped at the PVREF1 HABPVR intermediate (e.g.
PV_CODE 591 → pvref1 590) and fed THAT to `ie_estab_indices`, whose IEND/MYGRUP bracket then
gave the wrong ESTAB habitat GROUP (590→group 5, SHAB=+0.54) instead of the oracle's mapped
510→group 6 (SHAB=−0.06) — a ~+0.6 logit that inflated the stocking probability on EVERY plot.
Growth is unaffected (the DG/site path already calls `ie_habtyp`→ITYPE, identical for 590 & 510);
iet01 is inert (570→ITYPE 17→MTYPE 570). Out-of-range codes pass through unchanged.
"""
function ie_estab_kodtyp(kodtyp::Integer)::Int
    it = ie_habtyp(kodtyp)                                    # ie/habtyp.f DO 100/110 → ITYPE (1-30), 0 if out of range
    (1 <= it <= length(_IE_HABTYP_MTYPE)) ? Int(_IE_HABTYP_MTYPE[it]) : Int(kodtyp)
end

"""
    ie_estab_indices(habitat_code, forest_code) -> (ihab, iser, ifo, iphy, iprep)

Derive the AUTOES per-plot indices for a stand (no plot-specific overrides).
`habitat_code` = the MAPPED STDINFO habitat class (KODTYP/ICL5 = `ie_estab_kodtyp` output); `forest_code`
= the raw forest-location code (KODFOR). For iet01 stand-4 (570, 118) → (10, 4, 4, 3, 1).
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
                       esb_shift::Real = 0f0, esb_shift_pt::AbstractVector = Float32[], is_ingro::Bool = true,
                       nstore::AbstractVector = Int32[], pnn::AbstractVector = Float32[],
                       tpacre_ingro::Real = 0f0, point_small_tpa::AbstractVector = Float32[],
                       idup::Integer = 0, nsp::Integer = 23, variant = nothing,
                       point_slope::AbstractVector = Float32[], point_aspect::AbstractVector = Float32[],
                       stoadj::Real = 1f0, spec_mult::AbstractDict = Dict{Int32,Float32}(),
                       prep_sumup = nothing, pasmax::Real = Inf32,
                       point_ba::AbstractVector = Float32[], over_sp::AbstractVector = Float32[],
                       over_pt::AbstractMatrix = Array{Float32}(undef, 0, 0),
                       gentim::Real = 5f0, call_espadv::Bool = true,
                       plant_sp::AbstractVector = Int[], cont_ppba::Bool = false)
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
    # PROB1 = logistic(PN + ESB - ESB1)·STOADJ (estab.f:579-580); esb_shift = ESB-ESB1 (inventory calibration).
    # STOADJ = the STOCKADJ keyword multiplier (default 1.0 ⇒ inert). estab.f:578 clamps STOADJ≥0.001 before the
    # multiply. NOTE: STOADJ<0.0001 (STOCKADJ 0.0, as NATURAL implies) takes a SEPARATE estab.f branch (GO TO
    # 137/229/163 — no stocking model) which jl reaches via its scheduled-NATURAL path, NOT this tally; so the
    # ≈0 case is out of scope for this multiply. The clamp keeps any positive keyword value faithful.
    sa = Float32(stoadj); sa < 0.001f0 && (sa = 0.001f0)
    prob1 = (1f0 / (1f0 + exp(-(pn + Float32(esb_shift))))) * sa
    # estab.f:581-582 clamp FTEMP (=PROB1) to [0.0001, 0.9990] AFTER the STOADJ multiply. Critical when STOADJ>1 (or a
    # high-PN stand) would drive the logistic·STOADJ product past 1 — PROB1 sits in the NSTORE denominator (tpacre/
    # (prob1·300)), so an uncapped prob1>1 shrinks NSTORE and makes ingrowth OVER-book (measured: STOCKADJ 2.0 → jl
    # +588 TPA vs oracle −251 before this clamp). estab.f:583 additionally floors FTEMP at the plot's PNN(NCOUNT)+
    # 0.0001 (existing per-plot stocking prob); jl's prob1 is a stand-level scalar so that per-plot floor isn't
    # represented here — inert wherever prob1 ≥ the plot's PNN (the usual case; verified on the under-stocked fixture).
    prob1 < 0.0001f0 && (prob1 = 0.0001f0)
    prob1 > 0.9990f0 && (prob1 = 0.9990f0)
    # Per-IPREP PROB1 for the DISTURBANCE tally's per-plot IPPREP sampler (estab.f:572 ESTOCK carries a per-IPREP
    # SPRE term). Only needed when site prep is active (default ESPREP or a MECHPREP/BURNPREP keyword); otherwise
    # every plot uses the scalar prob1 (IPREP=1). Mirrors the prob1 pipeline above for IPREP∈{1,2,3}.
    _is_ie = variant !== nothing && variant isa InlandEmpire
    prob1_prep = Float32[]
    if _is_ie && prep_sumup !== nothing && !is_ingro
        prob1_prep = Vector{Float32}(undef, 3)
        for ip in 1:3
            pnp = ie_estock(idx.ihab, ip, sl, xc_st, xs_st, Float32(elev), ba, log(ba), tm,
                            sqrt(Float32(regt)), sqrt(Float32(bwaf)), Float32(bwb4), idx.ifo)
            v = (1f0 / (1f0 + exp(-(pnp + Float32(esb_shift))))) * sa
            v < 0.0001f0 && (v = 0.0001f0); v > 0.9990f0 && (v = 0.9990f0)
            prob1_prep[ip] = v
        end
    end
    # ★ PER-INVENTORY-POINT PROB1 (estab.f:466-582). The oracle recomputes the stocking logit PER PLOT from that
    # plot's inventory-point BAAA(NNID)/PSLO(NNID)/PASP(NNID) via ESTOCK, so on a multi-point FIA stand each point
    # gets its OWN stocking probability (a heavily-stocked point ⇒ near-zero regen, an open point ⇒ high). jl used a
    # single scalar prob1 from point-1's BA/slope for EVERY plot ⇒ on real FIA stands (4 points, heterogeneous BA)
    # the low-BA point-1 value over-produced ingrowth on the high-BA points (MEASURED trace 1627628023290487: pt3
    # BA 266 oracle PROB1 0.097 vs jl 0.60). Compute the per-point vector here (IE; EM stays scalar). nptids==1 ⇒
    # prob1_pt[1]==prob1 (byte-identical). ★ DISTURBANCE tally (NTALLY=1, post-thin) ALSO needs this: a heterogeneous
    # multi-point thin leaves some points bare (BAAA≈0 ⇒ high PROB1) and others stocked, and the oracle recomputes
    # PROB1 per inventory point identically to the ingrowth path (MEASURED FVSie_g16 3291804010690 @2020: oracle
    # per-point PROB1 [0.272,0.317,0.872,0.530], jl scalar point-1 0.272 ⇒ TPA 507 vs 1059, ~half). The per-IPREP
    # SPRE term (estock.f) CANCELS in PROB1=logistic(PN+ESB−ESB1) whenever the ESB1 correction is active (SPRE is
    # time/BA-independent and appears in both the end-of-cycle PN and the inventory ESB1), so a disturbance PROB1
    # with esb_shift≠0 is per-point AND iprep-independent (the oracle prints one PROB1 per point regardless of the
    # plot's IPREP). The caller passes point_ba (the disturbance all-tree per-point BAAA) only when esb_shift≠0.
    _npt_pb = idup > 0 ? max(1, div(Int(dupnpt), Int(idup))) : 1
    prob1_pt = Float32[]
    if _is_ie && !isempty(point_ba)
        prob1_pt = Vector{Float32}(undef, _npt_pb)
        @inbounds for pt in 1:_npt_pb
            ba_pt  = clamp(pt <= length(point_ba)     ? Float32(point_ba[pt])     : ba,  1f0, 400f0)   # BAAA(NNID) clamp [1,400]
            sl_pt  =       pt <= length(point_slope)  ? Float32(point_slope[pt])  : sl                 # PSLO(NNID)
            asp_pt =       pt <= length(point_aspect) ? Float32(point_aspect[pt]) : asp                # PASP(NNID)
            pn_pt = ie_estock(idx.ihab, idx.iprep, sl_pt, cos(asp_pt), sin(asp_pt), Float32(elev), ba_pt,
                              log(ba_pt), tm, sqrt(Float32(regt)), sqrt(Float32(bwaf)), Float32(bwb4), idx.ifo)
            # PER-POINT stocking shift ESB−ESB1(NNID): each point's ingrowth stocking is corrected by ITS OWN
            # inventory prediction (estab.f:557,599). Empty ⇒ the scalar esb_shift (point-1 value) — over-corrects
            # the open points of a heterogeneous multi-point stand (the M333 AUTOES over-production bug).
            shift_pt = pt <= length(esb_shift_pt) ? esb_shift_pt[pt] : Float32(esb_shift)
            v = (1f0 / (1f0 + exp(-(pn_pt + shift_pt)))) * sa                                            # STOADJ already clamped ≥0.001 (sa)
            v < 0.0001f0 && (v = 0.0001f0); v > 0.9990f0 && (v = 0.9990f0)
            prob1_pt[pt] = v
        end
    end
    prob1_of(pt) = isempty(prob1_pt) ? prob1 : prob1_pt[clamp(pt, 1, length(prob1_pt))]
    # #143: INGROWTH NSTORE = the existing small-tree (DBH<REGNBK) stocking (estab.f:589 NSTORE=INT(PLPROB·DUPNPT/
    # (FTEMP·300)+0.5)). PLPROB·DUPNPT = the current DBH<2.999 TPA (measured live: 595·50=29750 = self-thinned
    # cohort), so NSTORE=INT(tpacre/(prob1·300)+0.5) per plot (nptids/idup cancel → uniform; exact single-point,
    # correct uniform-density multi-point). Then NEWTPP=max(0,ITPP−NSTORE)=0 on an already-stocked plot ⇒ ingrowth
    # books ≈0 (live 0.1 TPA), vs the old es_nstore=0 that booked a full fresh cohort (253 TPA, ~2500× over).
    if is_ingro && length(nstore) == Int(dupnpt)
        npt = length(point_small_tpa)
        if npt > 0 && idup > 0 && npt * Int(idup) == Int(dupnpt)
            # PER-POINT NSTORE (estab.f:313/589): NSTORE(pt)=INT(PLPROB·DUPNPT/(PROB1·300)+0.5), PLPROB=point_TPA/DUP
            # ⇒ point_TPA·NPTIDS/(PROB1·300). Fill each point's contiguous IDUP-plot block (FVS NCOUNT order). PROB1 is
            # the POINT's own stocking prob (prob1_of) — a high-BA point has a low PROB1 ⇒ its small-tree stock maps to
            # a LARGER NSTORE (more existing trees per plot) ⇒ less NEWTPP, matching the oracle's per-point suppression.
            fill!(nstore, Int32(0))
            @inbounds for pt in 1:npt
                ns_pt = floor(Int32, Float32(point_small_tpa[pt]) * Float32(npt) / (prob1_of(pt) * 300f0) + 0.5f0)
                ns_pt <= 0 && continue
                base = (pt - 1) * Int(idup)
                for k in 1:Int(idup); nstore[base + k] = ns_pt; end
            end
        elseif Float32(tpacre_ingro) > 0f0
            nsval = floor(Int32, Float32(tpacre_ingro) / (prob1 * 300f0) + 0.5f0)   # fallback: uniform (pre-#143)
            fill!(nstore, nsval)
        end
    end
    # espadv/espxcs occupancy = OCURHT(ihab,sp)·XESMLT(sp)·OCURNF(ifo,sp) (estab.f). jl formerly applied only
    # OCURHT; the OCURNF (per-National-Forest occupancy) gate was omitted (validated inert on iet01 where
    # OCURNF·XESMLT=1). On EM bare-establishment stands this let species EXCLUDED on the stand's NF (e.g. PP on
    # forest 108, OCURNF(ifo,PP)=0) over-establish and over-grow (#143: PP htg1≈4.22 vs DF≈1.50 ⇒ BA 6.5×). Now
    # multiplied in via autoes_ocurnf(variant,ifo,sp) (XESMLT=1 default). EM has its table; other variants→1.0.
    # occ = OCURHT·XESMLT·OCURNF (estab.f); SPECMULT sets the per-species XESMLT (esnutr.f 95), default 1.0.
    occ = Float32[ie_ocurht(idx.ihab, sp) * (variant === nothing ? 1f0 : autoes_ocurnf(variant, Int(idx.ifo), sp)) *
                  (isempty(spec_mult) ? 1f0 : get(spec_mult, Int32(sp), 1f0)) for sp in 1:nsp]
    # OVER(i) = per-species OVERSTORY basal area (D≥REGNBK) at the tally point (dense.f:214 OVER(ISPC,IP)).
    # ESPADV/ESPXCS/ESPSUB add a per-species "seed-source present" bump when OVER(i)>9.95 (a species with
    # overstory ⇒ higher regen probability). jl formerly hardcoded over=0 ⇒ the bump NEVER fired ⇒ species with
    # overstory seed source (esp. LP/GF) were systematically UNDER-picked in the best/excess mix. `over_sp` is the
    # point-1 per-species D≥REGNBK BA supplied by ie_autoes_establish! (IE only; empty ⇒ zeros ⇒ EM/other callers
    # byte-identical). Species order 1-10 = WP WL DF GF WH RC LP ES AF PP (= IE variant indices 1-10).
    over = zeros(Float32, 10)
    @inbounds for i in 1:min(10, length(over_sp)); over[i] = Float32(over_sp[i]); end
    _npt = idup > 0 ? max(1, div(Int(dupnpt), Int(idup))) : 1     # inventory points = dupnpt/idup (=nptids)
    tally_pt = zeros(Float64, nsp, _npt)                          # per-point established TPA (for plot_id placement)
    # IE AND EM emit the faithful per-tree height-class / WK4 records (advance/subsequent/excess); other variants
    # keep the collapsed single-record path (emit=nothing ⇒ byte-identical). ihtser = MYHTS(IHAB).
    # EM shares the estb/estab.f DO-99/33/228 driver and borrows IE/NI coefficients for its 7 establishing
    # species (WB≡WP,WL,DF,LP,ES,AF,PP at positions 1,2,3,7,8,9,10 — verified identical in esadvh/essubh/esxcsh/
    # esdlay + XMIN/HHTMAX/BNORML/MYHTS); OCURNF zeroes the non-shared positions (4,5,6=LM/LL/RM, 11-19), so the
    # emit path is numerically faithful for EM. `is_ie` (below) still gates the IE-ONLY tally internals (PSUB add,
    # prep_active, per-point PROB1) so the validated EM tally TPA + RNG stream stay byte-identical — emit only adds
    # the per-record height/WK4 emission over the SAME draws.
    _is_em = variant !== nothing && variant isa EasternMontana
    emit_recs = (_is_ie || _is_em) ? NTuple{5,Float64}[] : nothing
    ihtser = _IE_MYHTS[clamp(Int(idx.ihab), 1, length(_IE_MYHTS))]
    tally = ie_autoes_tally(seed0 = seed0, nplots = Int(dupnpt), ihab = idx.ihab, iser = idx.iser,
                            ifo = idx.ifo, iprep = idx.iprep, iphy = idx.iphy, xcos = xc_sp, xsin = xs_sp,
                            slo = sl, elev = Float32(elev), baa = ba, regt = Float32(regt),
                            bwaf = Float32(bwaf), bwb4 = Float32(bwb4), prob1 = prob1, dupnpt = Float32(dupnpt),
                            occ = occ, over = over, time = tm, is_ingro = is_ingro, nstore = nstore, pnn = pnn,
                            nsp = nsp, idup = Int(idup), tally_pt = tally_pt, wk6fill = Int(dupnpt),
                            point_slope = point_slope, point_aspect = point_aspect, prep_sumup = prep_sumup,
                            prob1_prep = prob1_prep, is_ie = _is_ie, pasmax = pasmax, prob1_pt = prob1_pt,
                            emit = emit_recs, ihtser = ihtser, gentim = Float32(gentim), call_espadv = call_espadv,
                            # Per-point species tables gated to the INGROWTH tally (the establishment lever); the
                            # DISTURBANCE path keeps the validated scalar/point-1 tables (empty ⇒ fallback).
                            point_baa = (is_ingro || cont_ppba) ? point_ba : Float32[],
                            over_pt = (is_ingro || cont_ppba) ? over_pt : Array{Float32}(undef, 0, 0),
                            plant_sp = plant_sp)
    return (tally = tally, tally_pt = tally_pt, prob1 = prob1, idx = idx, emit = emit_recs)
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
    # INADV=1 bare-stand bypass (estab.f:319,511 `IF(INADV.EQ.1 .OR. NTALLY.NE.1) GO TO …`): on a bare stand FVS
    # runs the establishment scheduler REGARDLESS of the LAUTAL/LINGRW auto-tally flags (so a bare stand with
    # NOAUTALY/NOINGROW still gets the PLANT-forced / bare ingrowth — rule (4)/(3) below). est.inadv persists from
    # the cycle-1 latch (line ~1771); catch cycle-1 inline before it's latched. NON-bare stands (INADV=0) keep the
    # original gate exactly ⇒ byte-identical.
    _bare = est.inadv || (Int(s.control.cycle) + 1 == 1 && s.trees.n == 0)
    (est.lautal || est.lingrw || _bare) || return false
    nsp = nspecies(s.variant)
    # The estb habitat bracket (esplt2.f IEND/MYGRUP, = _IE_ESTAB_IEND/MYGRUP) keys off ICL5 = the FVS/NI habitat
    # CODE. IE's p.habitat_code IS that code; EM's p.habitat_code is IEMTYP (1-118) whose NI code is EM_JTYPE[iemtyp]
    # (em/habtyp.f:123 KODTYP=JTYPE(IEMTYP)). The OCURHT/CHAB/ESTOCK tables + _IE_ES_XMIN are the SHARED estb data
    # (species 1-10 = WP WL DF GF WH RC LP ES AF PP), so no species crosswalk is needed — EM's dry habitats zero the
    # wet-side estb species (OCURHT verified bit-identical to live FVSem for IHAB=3), and EM's establishing species
    # (WL/DF/LP/ES/AF/PP) sit at the matching estb indices 2,3,7,8,9,10.
    # esplt2.f keys the ESTAB habitat bracket off ICL5 = the NORTH-IDAHO habitat code (JTYPE value), which
    # ie/habtyp.f sets to JTYPE(NITYPE) (line 124-125) BEFORE the final MTYPE(ITYPE) remap used only for GROWTH.
    # s.plot.habitat_code already holds that NI code (ie_pvref1's numeric HABPVR / the raw numeric habitat code =
    # ICL5), so feed it DIRECTLY — exactly as the EM branch feeds EM_JTYPE[iemtyp]. ★ MEASURED (FVSie_g16
    # 645170890126144 @2057): oracle ICL5=590 (=JTYPE(55)) → esplt2 group 5 (SHAB=+0.540 ⇒ PROB1 0.726), whereas a
    # PRIOR ie_estab_kodtyp remap sent 590→MTYPE 510 → group 6 (SHAB=−0.060 ⇒ PROB1 0.586), a −0.60 stocking logit
    # that under-produced AUTOES ingrowth by ~100 TPA/cycle. Growth is unaffected: the DG/site path re-derives ITYPE
    # via ie_habtyp(habitat_code), and ie_habtyp(590)==ie_habtyp(510)==ITYPE 12 (identical MTYPE 510).
    ihab_code = s.variant isa EasternMontana ?
        Int(EM_JTYPE[clamp(Int(s.plot.habitat_code), 1, 118)]) :
        Int(s.plot.habitat_code)
    per = round(Int, fint)
    year = Int(current_cycle_year(s))
    next_year = year + per
    inv_year = Int(s.control.cycle_year[1])
    # esnutr.f:113 — IDSDAT defaults to IY(1)-20 (20 yrs before inventory), NOT -9999. Without this the 40-yr
    # ingrowth-gap rule (next_year - idsdat ≥ 40) would spuriously fire for every stand at cycle 1.
    est.idsdat == Int32(-9999) && (est.idsdat = Int32(inv_year - 20))
    icyc = Int(s.control.cycle) + 1
    itrn = s.trees.n
    # EZCRUISE / NOTREES auto-invoke (initre.f:280-289 → esinit.f:79 ESEZCR sets INADV=1): a BARE stand (< 1
    # projectable tree record at inventory) auto-invokes the EZCRUISE regeneration option. INADV=1 then PERSISTS
    # for the whole run (getstd/putstd INTS(33)) and makes estab.f SKIP the ESB inventory-stocking calibration in
    # EVERY cycle (estab.f:319/511). A bare stand has itrn==0 at cycle 1 (pre-regen) — the same condition that
    # triggers the LINGRW ingrowth path — so detect it here and latch the persistent flag.
    icyc == 1 && itrn == 0 && (est.inadv = true)
    # User-scheduled TALLY/TALLYONE/TALLYTWO (esin.f 16/11/12; esnutr.f:163-252): a 427/428/429 whose date falls in
    # THIS cycle's window [year,next_year) and is not stale (KDT+1-IDSDAT ≤ 20) forces a tally, overriding the
    # automatic AUTOES rules. NTALLY = IACTK-427 collapsed to {1 (TALLY/TALLYONE), 2 (TALLYTWO continuation)}. The
    # strict window fires each scheduled tally exactly once; the automatic 20-yr continuation (rule 2) then handles
    # its 2nd pass. The ESTAB-END 427 sits at inv-20 (stale, >20yr) so it is correctly skipped — no regression.
    kdt = next_year - 1
    # STOCKADJ (esin.f opt 13 → activity 440): the stockability multiplier scheduled for THIS cycle's window,
    # default 1.0 (inert). Mirrors the TALLY date→cycle mapping (a cycle number → this cycle's year, else a
    # calendar year). est.stoadj is then consumed by ie_autoes_run (PROB1 = logistic(…)·STOADJ, estab.f:578).
    est.stoadj = 1f0
    for a in s.control.schedule
        a.icflag == Int32(440) || continue
        ay = Int(a.year)
        aidt = (0 < ay < 1000) ? (ay == icyc ? year : -1) : ay
        (year <= aidt < next_year) || continue
        est.stoadj = a.params[1]
        break
    end
    sched_fire = false; sched_ntally = 0
    for a in s.control.schedule
        (a.icflag == 427 || a.icflag == 428 || a.icflag == 429) || continue
        ay = Int(a.year)
        idt = (0 < ay < 1000) ? (ay == icyc ? year : -1) : ay   # cycle-number → this cycle's year; else calendar
        (year <= idt < next_year && (kdt + 1 - idt) <= 20) || continue
        est.idsdat = Int32(idt)
        # NTALLY=IACTK-427 ⇒ 1 for TALLY/TALLYONE; TALLYTWO(429)=2 ONLY if a prior TALLYONE(428) at an earlier
        # date exists, else it collapses to 1 (esnutr.f:204-210). Standalone TALLYTWO ≡ a single tally.
        sched_ntally = (a.icflag == 429 && any(b -> b.icflag == 428 && Int(b.year) < ay, s.control.schedule)) ? 2 : 1
        est.ntally = Int32(sched_ntally)
        sched_fire = true
        break
    end
    # esnutr.f:348 OPFIND(2,[430,431]) — PLANT/NATURAL activities scheduled THIS cycle (same due-window as establish!,
    # engine/establishment.jl:337-339). Drives the rule-6 catch-all so a planted stand gets automatic ingrowth in the
    # plant cycle. (fvscyc = icyc for cycle-number-encoded dates.)
    npnats = count(a -> (a.icflag == Int32(430) || a.icflag == Int32(431)) &&
                        ((year <= Int(a.year) < next_year) || (0 < Int(a.year) < 1000 && Int(a.year) == icyc)),
                   s.control.schedule)
    # Species of the PLANT/NATURAL trees co-located on each plot THIS tally (estab.f:970-985 DO 322 IPNSPE=PRMS(1)),
    # to merge into the AUTOES NBEST best/excess pool (ie_autoes_tally). Same due-window + non-zero-TPA filter as
    # establish! (estab.f:415-421 deletes PRMS(2)/(3)≤0.001). Order = schedule order (DO 322 IDO=1,NTODO). Empty when
    # no PLANT/NATURAL is due ⇒ the tally's NBEST is byte-identical. establish! still books the actual planted trees.
    _plant_sp = Int[]
    for a in s.control.schedule
        ((a.icflag == Int32(430) || a.icflag == Int32(431)) &&
         ((year <= Int(a.year) < next_year) || (0 < Int(a.year) < 1000 && Int(a.year) == icyc))) || continue
        (length(a.params) >= 3 && Float32(a.params[2]) > 0.001f0 && Float32(a.params[3]) > 0.001f0) || continue
        push!(_plant_sp, round(Int, a.params[1]))
    end
    fire, _ntally = sched_fire ? (true, sched_ntally) :
        ie_autoes_schedule!(est, icyc, year, next_year, itrn, est.last_xtes, inv_year, npnats)
    est.last_xtes = 0f0                       # consume the removal fraction (one cycle only)
    fire || return false
    year in est.autoes_years_done && return false   # AUTOES own re-entry guard (NOT establish!'s years_done — a PLANT this cycle must NOT suppress the natural tally, estab.f runs both together)

    p = s.plot
    nptids = max(1, Int(p.points_inv) - Int(p.nonstockable))
    idup   = max(1, cld(Int(est.minrep), nptids))   # MINPLOTS keyword (esin.f MINREP; default 50)
    dupnpt = Float32(nptids * idup)
    # ★ EMPTY-INVENTORY-POINT PSLO/PASP (esplt2.f:168-227). When PROJECTABLE trees exist (IPTKNT>0 ⟺ s.trees.n>0),
    # esplt2.f pads the point list to IPTINV−NONSTK with empty points (no tree record) whose PSLO/PASP are set to −1
    # (:184-193) and then RESOLVED to the STAND slope/aspect XXSLP/XXASP (:216-227, XXSLP=ISLOP·0.01, XXASP=IASPEC·
    # 0.0174533) — NOT 0. The FIA reader (fia_database.jl) fills these empty points with 0, which is correct ONLY for
    # the IPTKNT==0 path (a bare / all-dead stand takes esplt2.f's top block and empty points keep the BLKDAT 0). For
    # IPTKNT>0 the reader's 0 is WRONG: ESTPP's BB carries −2.76·PSLO, so a spurious 0 slope on the empty points of a
    # sloped stand UN-suppresses trees-per-plot ⇒ AUTOES ingrowth OVER-produces (MEASURED FVSie_g16 39598164010690
    # @2039: oracle empty-point PSLO 0.65 ⇒ ITPP 20/21, jl 0 ⇒ ITPP 46/34; +30–120 TPA/cycle across the 4-point
    # sparse `none` stands). The appended empty points are exactly those with p.point_ids[k]==0 (see below). Build a local corrected copy — never mutate
    # p.point_slope (a later cycle's established trees must not turn an inventory-empty point into a tree point; PSLO
    # is frozen at inventory). Bare stands (s.trees.n==0) keep the reader's values (IPTKNT==0 / IPINFO=0 path).
    # An inventory-empty point is exactly one with p.point_ids[k]==0 (the reader's IPVEC, MAXPLT-sized and zero-init;
    # the inventory tree points get their raw plot number ≥1, the appended empty points stay 0). This is frozen at
    # inventory read — established trees never rewrite IPVEC — so it is the correct inventory-emptiness test.
    pslo_es = p.point_slope; pasp_es = p.point_aspect
    if s.trees.n > 0 && !isempty(p.point_slope)
        _need = false
        @inbounds for k in 1:length(p.point_slope)
            (k > length(p.point_ids) || p.point_ids[k] == 0) && (_need = true; break)
        end
        if _need
            pslo_es = copy(p.point_slope); pasp_es = copy(p.point_aspect)
            @inbounds for k in 1:length(pslo_es)
                if k > length(p.point_ids) || p.point_ids[k] == 0
                    pslo_es[k] = Float32(p.slope); pasp_es[k] = Float32(p.aspect)   # esplt2.f XXSLP/XXASP
                end
            end
        end
    end
    # AUTOES ESRANN seed chain (estab.f:290-295 + the per-plot ESAVE reseed). A NEW disturbance/ingrowth tally
    # (NTALLY==1|99) draws seed0 = ESRANN(es_stream) from the continuing establishment stream, then advances the
    # stream to the post-tally state = the last plot's ESAVE (ie_autoes_plot_seeds[dupnpt+1], validated == the
    # live ESS0 at the next tally). A continuation (NTALLY≥2) reuses this tally's seed0.
    # es_stream starts fresh at ESSS=55329.
    # ★#143: the ESAVE chain's per-plot body length = 16 + 3·NOFSPE + 2·MAXTPP(IHAB) (estab.f: EMSQR+ESTPP+
    # NUMSPE-WK6+species-WK6 + ADV/SUBS(NOFSPE) + heights(2·NOFSPE) + excess-WK6(2·MAXTPP) + ESAVE), the SAME
    # count ie_autoes_tally derives per plot. It VARIES with the stand's habitat (MAXTPP) — the old default 135
    # was iet01-specific (ihab 10, MAXTPP 25). Using it for a stand at another habitat (e.g. under-stocked
    # 550→ihab 9, MAXTPP 21 ⇒ body 127) advanced the stream WRONG and desynced every tally after the first
    # (jl 2nd seed 61677 vs live 61997). Deriving body_es from THIS stand's ihab makes the seed chain bit-exact
    # vs live (43303→61997→49053 on the under-stocked fixture).
    if _ntally == 1 || _ntally == 99
        _idx_es = ie_estab_indices(ihab_code, Int(p.user_forest_code))
        body_es = 16 + 3 * nsp + 2 * _IE_MAXTPP[_idx_es.ihab]
        # estab.f:290-295 CALL ESRANN(DRAW) draws seed0 from the SHARED establishment stream ESS0 (ESRNCM),
        # the SAME stream ESUCKR's BACHLO(0,0.5,ESRANN) advances when a disturbance (fire/thin/salvage) removed
        # sprouting-species records earlier this cycle (esuckr.f runs in ESNUTR, BEFORE the AUTOES tally). Model
        # ESS0 as the single canonical `s.rng.es0` (FVSRng :estab; seeded 55329, advanced by esuckr! through
        # bachlo(:estab)), exactly as the SN/base establish! path does (establishment.jl:422-426). Reading the
        # private est.es_stream here instead LOST every ESUCKR draw ⇒ on a disturbance stand AUTOES seeded from a
        # fresh stream (draw #1, ESDRAW 43303) instead of the post-ESUCKR draw (e.g. draw #82, ESDRAW 34087 on
        # the 81-sprout-draw simfire stand 1143093321290487), a systematic post-disturbance seedling miscount.
        # dr is drawn from the raw post-ESUCKR ESS0 (no odd-adjust — estab.f:291 draws directly; ESRNSD odd-forces
        # only on reseed). VALIDATED: 7/8 seed-affected thinbba stands + the simfire stand now match the live
        # oracle's ESDRAW bit-for-bit; on the base `none` regime esuckr! is a no-op so ESS0 stays 55329 and seed0
        # is the previously-validated 43303 chain (no regression). The post-tally ESAVE is written back so the NEXT
        # cycle's esuckr!/tally continue the same stream (both directions of the ESS0 ↔ ESUCKR coupling).
        dr = esrann!(s.rng)                                                   # estab.f:291 — advances shared ESS0
        seed0 = floor(Int, dr * 100000f0 + 0.5f0)
        est.es_seed = Float32(seed0)                                          # save for the continuation
        final_seed = ie_autoes_plot_seeds(seed0, Int(dupnpt) + 1; wk6 = Int(dupnpt), body = body_es)[end]  # ESAVE_last (odd-forced)
        est.es_stream = Float32(final_seed)
        s.rng.es0 = Float64(final_seed)                                       # leave ESS0 at the post-tally ESAVE
    else
        seed0 = round(Int, est.es_seed)                                       # continuation reuses seed0
    end
    # ESTOCK/species-prob BAA = the per-INVENTORY-POINT basal area BAAA(NNID) (estab.f:482, dense.f:213), NOT the
    # whole-stand BA. After a heavy overstory removal the regen point is bare → BAAA≈0 → TBAAA=max(BAAA,1)=1 →
    # ESTOCK PN high → PROB1 high (the disturbance re-stocking pulse). Using stand_ba (which keeps the residual
    # overstory) collapsed PROB1 to ~0.55 and under-produced ~40%. For the single-inventory-point stand the point
    # is index 1 (validated: jl point_ba[1]=40 vs live BAAA=41.93 at cyc1; 0→1 vs 1 at the bare disturbance tallies).
    # ★ #143 (2026-08-06): a D≥REGNBK OVERSTORY-only baaa was TRIED and REFUTED by live measurement — live BAAA is
    # NOT overstory-filtered. On DB stand 753199439290487 (0.2"-QMD, all sub-3" stems) live BAAA=20.66/205.6 counts
    # the tiny trees; the overstory-only version (4.59/12.72) is FARTHER from live than point_ba[1] (13.18/216.12).
    # So BAAA(NNID) is the per-point ALL-tree BA (dense.f REGNBK is ~0 for these stands), and point_ba[1] is the
    # right concept. Residual = a per-point scaling/attribution Δ (live 20.66 vs jl 13.18, ~1.5×) — see AUTOES doc.
    baaa = (isempty(s.density.point_ba) ? 0f0 : s.density.point_ba[1])
    # ★#143 follow-on (MEASURED FVSie_g16): the establishment stocking model (ESTOCK / species-probs) uses the
    # PER-PLOT slope/aspect PSLO(NNID)/PASP(NNID), which FVS reads from the TREE records (esplt1.f:69-70, IPINFO=2),
    # NOT the stand SLOPE/ASPECT. On the under-stocked fixture FVS_TREEINIT_COND carries SLOPE=33/ASPECT=160 while
    # FVS_STANDINIT_COND has SLOPE=0/ASPECT=0 (the stand value drives DGF, but establishment uses the per-plot value).
    # jl fed the stand slope/aspect (0/0) into ESTOCK ⇒ the XCOSAS·SLO·XBAA interaction (+0.30 to ESB1) and the
    # aspect·SQSQ terms vanished ⇒ PROB1 0.463 vs live 0.617. jl already reads these per-plot values into
    # p.point_slope/p.point_aspect (used for ESTPP); use them for the stocking PN + ESB1 too. Empty (TREEDATA / no
    # per-plot topo) ⇒ fall back to the stand slope/aspect (inert where they match, e.g. iet01).
    es_slope  = isempty(pslo_es)  ? Float32(p.slope)  : Float32(pslo_es[1])
    es_aspect = isempty(pasp_es) ? Float32(p.aspect) : Float32(pasp_es[1])
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
    # ESB inventory calibration (estab.f:319-326,579): estab.f applies the ESB-ESB1 actual-vs-predicted stocking
    # correction whenever INADV=0 .AND. NTALLY=1. INADV=0 ⟺ KDT+1-IY(1) ≤ 20 (estab.f:175), i.e. the tally is within
    # 20 yr of inventory (kdt=next_year-1 ⇒ `next_year-inv_year ≤ 20`). ★#143 follow-on (MEASURED FVSie_g16): the
    # FIRST-cycle AUTOES ingrowth IS a NTALLY=1/INADV=0 tally (dump: ntally=1 inadv=0 esb=0.0778 esb1=-0.3583 at ic=1;
    # ntally=1 inadv=1 esb=0 at ic=3/5). jl's old `idsdat==inv_year` gate missed it (ingrowth idsdat=-1) ⇒ esb_shift=0
    # ⇒ PROB1 0.463 vs live 0.617 ⇒ ingrowth (scaled by PROB1·NEWTPP/ITPP) under-booked ~25% ⇒ the .sum TPA residual.
    # ESB = logit(clamp(logistic(-5.174+0.851·ln(TPACRE)),0.10,0.90)) from the current small-tree (DBH<2.999) TPA;
    # ESB1 = ESTOCK(BAAOLD=INVENTORY per-point BA, TIME=0). Computed once (persists in est.esb_shift), applied only
    # while INADV=0. NOTE: baaa here is the PRE-growth (inventory/start-of-cycle) point BA — the right BAAOLD for ESB1;
    # the end-of-cycle stocking PN uses the POST-growth BA refreshed below.
    # INADV=0 AND NTALLY=1-equivalent: the oracle recomputes ESB/ESB1 only on a FRESH tally (NTALLY=1); a
    # continuation (NTALLY≥2) resets them to 0 ⇒ no correction. jl's fresh tallies are ntally∈{1 (disturbance),
    # 99 (ingrowth)}; a continuation is ntally≥2. Guard so a within-20yr continuation doesn't get esb_shift.
    # ★ BARE-STAND FIX (EZCRUISE INADV=1): estab.f gates the ESB calibration on INADV=0 (line 319/511). A bare
    # stand latches est.inadv (above), so its ingrowth uses PROB1 = logistic(PN) with NO shift — MEASURED live
    # FVSem_clean on bare stand 85271557010661: ESB=0 ESB1=0 PROB1=0.3584 (=logistic(PN=-0.5822)), whereas jl
    # previously applied esb-esb1≈-1.39 ⇒ PROB1 0.1222 ⇒ ~2.9× under-production of ingrowth TPA.
    esb_shift = 0f0
    if _ntally == 2 && !isnan(est.esb_shift)
        # DISTURBANCE CONTINUATION (NTALLY=2): reuse the fresh tally's stored ESB−ESB1. FVS resets ESB1(NCOUNT)/ESB
        # ONLY on NTALLY=1 (estab.f:264 GO TO 276 skips the DO-12 reset; :319/:511 skip the recompute for NTALLY≠1),
        # so the persisted ESB1(NCOUNT) is applied at estab.f:579 for the continuation too — it is NOT re-gated by
        # the 20-yr inventory window. jl formerly zeroed it (the `_ntally==1||99` gate excluded 2) ⇒ the cyc-3
        # continuation tally lost its +2.1 stocking shift ⇒ PROB1 0.51 vs oracle 0.90 ⇒ post-thin under-production (D1).
        esb_shift = est.esb_shift
    elseif !est.inadv && (next_year - inv_year) <= 20 && (_ntally == 1 || _ntally == 99)
        # ESB is recomputed from the CURRENT small-tree TPACRE at EVERY fresh tally (estab.f:301-322 reads the live
        # DBH<REGNBK stocking each call; only ESB1's BAAOLD is the frozen inventory value). The `isnan` compute-once
        # guard wrongly FROZE ESB at the first tally: on a stand with a SECOND within-20yr ingrowth (NTALLY=99) — a
        # PLANT-triggered rule-6 tally after the cyc-1 ingrowth — the cyc-1 esb_shift (tpacre≈3, ESB floored, shift
        # ≈−1.23) was reused although cyc-2 tpacre has jumped to ~190 (the cyc-1 cohort), so ESB rises and the shift
        # flips to ≈+0.26. The stale −1.23 collapsed PROB1 0.157 vs oracle 0.63 ⇒ ingrowth NSTORE=INT(tpacre/(PROB1·
        # 300)) inflated 1→4 ⇒ NEWTPP≈0 ⇒ ~7× under-production of plant-cycle natural regen. Recompute on every
        # ingrowth tally (SAFE for none/thinbba/salvage: they have only ONE within-20yr ingrowth ⇒ isnan-equivalent;
        # the NTALLY=2 continuation still reuses the stored value via the branch above). ESB1 stays the frozen
        # inventory prediction (inv_baaold), so this only refreshes the current-stocking ESB half.
        if isnan(est.esb_shift) || _ntally == 99
            idx0 = ie_estab_indices(ihab_code, Int(p.user_forest_code))
            # ESB reads the small-tree (D<REGNBK) TPACRE from the live tree list at the tally (estab.f:301-322 —
            # the current small-tree stocking, NOT frozen at inventory: only ESB1's BAAOLD is the frozen inventory value).
            tpacre = 0f0
            @inbounds for i in 1:s.trees.n; s.trees.dbh[i] < 2.999f0 && (tpacre += s.trees.tpa[i]); end
            tpacre < 1f0 && (tpacre = 1f0)
            esa = clamp(1f0 / (1f0 + exp(-(-5.17397f0 + 0.85131f0 * log(tpacre)))), 0.10f0, 0.90f0)
            esb = -log(1f0 / esa - 1f0)
            # BAAOLD clamped to [1,400] (estab.f:487-489: IF(BAAOLD.LT.1)=1; IF(BAAOLD.GT.400)=400). ★ MEASURED
            # NSTORE over-book fix (dense already-stocked IE stands): on a BA>400 stand the ESB1=ESTOCK(BAAOLD)
            # prediction runs off the un-clamped tail (jl baaold 536 ⇒ esb1 −10.4 vs oracle −6.45 at BAAOLD=400),
            # inflating esb_shift=ESB−ESB1 (jl 8.21 vs oracle 4.26) ⇒ the per-point PROB1 saturates (jl 0.93 vs
            # oracle 0.207) ⇒ the first ingrowth tally over-books ~4.5× (lead 12281578010690: jl 894 vs oracle 199
            # TPA). FVS clamps BOTH BAA (prob1_pt already clamps) AND BAAOLD to 400; jl omitted the BAAOLD clamp.
            # IE-guarded (shared with EM; inert for EM's bare BA<400 stands, but scoped to honor the collision seam).
            # BAAOLD = the INVENTORY per-point OVERSTORY BAAINV (ESFLTR-frozen, snapshot_esb_inputs!), NOT the CURRENT
            # (post-thin) per-point BA `baaa`. ESB1 is defined as the stocking predicted AT INVENTORY (estab.f:506,
            # 536 ESTOCK(BAAOLD=BAAINV(NNID))); using the thinned `baaa` (e.g. 40 vs inventory-overstory 133) collapsed
            # ESB1 (→ +2.4 instead of −0.97) ⇒ esb_shift wrong-signed ⇒ PROB1 0.18 vs oracle 0.87 ⇒ the post-thin
            # re-stocking cohort under-produced ~10× (D1). inv_baaold is overstory-only (D≥REGNBK) per ESFLTR; falls
            # back to the current baaa only if the setup snapshot is absent (never, for IE/EM).
            baaold_raw = isnan(est.inv_baaold) ? baaa : est.inv_baaold
            baaold = (s.variant isa InlandEmpire) ? clamp(baaold_raw, 1f0, 400f0) : max(baaold_raw, 1f0)  # BAAINV
            asp0 = es_aspect; sl0 = es_slope                 # per-plot PSLO/PASP (from tree records), not stand
            esb1 = ie_estock(idx0.ihab, idx0.iprep, sl0, cos(asp0), sin(asp0), Float32(p.elevation),
                             baaold, log(baaold), 0f0, 0f0, 0f0, 0f0, idx0.ifo)   # ESTOCK(BAAOLD, TIME=0)
            est.esb_shift = esb - esb1
            # ★ PER-INVENTORY-POINT stocking shift ESB−ESB1(NNID) (estab.f:507-557 — ESB1 is computed inside the
            # per-plot loop from THAT point's BAAINV(NNID)/PSLO/PASP). The scalar esb_shift above uses point-1's
            # BAAINV, which over-corrects the OPEN points of a heterogeneous multi-point stand: a dense point 1
            # (BAAINV→400 ⇒ ESB1≈−6.2 ⇒ +4 logit) applied to a low-BA point pushes its ingrowth PROB1 to ~0.95
            # (vs the oracle's ~0.25) ⇒ ~4× AUTOES over-production (MEASURED FVSie_g16 1856003217290487 @2043: jl
            # 822 vs oracle 184 TPA). Compute the frozen per-point shift = ESB − ESTOCK(BAAINV(NNID)) so each
            # point's ingrowth stocking is corrected by ITS OWN inventory prediction. IE ingrowth only; empty ⇒
            # scalar fallback (byte-identical). Uses each point's PSLO(NNID)/PASP(NNID) (per-point tree topo).
            if s.variant isa InlandEmpire && !isempty(est.inv_point_baaold)
                npt_e = length(est.inv_point_baaold)
                shpt = Vector{Float32}(undef, npt_e)
                @inbounds for pt in 1:npt_e
                    bo = clamp(est.inv_point_baaold[pt], 1f0, 400f0)                       # BAAINV(NNID) clamp [1,400]
                    slp = pt <= length(p.point_slope)  ? Float32(p.point_slope[pt])  : sl0  # PSLO(NNID)
                    asp = pt <= length(p.point_aspect) ? Float32(p.point_aspect[pt]) : asp0 # PASP(NNID)
                    e1p = ie_estock(idx0.ihab, idx0.iprep, slp, cos(asp), sin(asp), Float32(p.elevation),
                                    bo, log(bo), 0f0, 0f0, 0f0, 0f0, idx0.ifo)             # ESTOCK(BAAINV(NNID), TIME=0)
                    shpt[pt] = esb - e1p
                end
                est.esb_shift_pt = shpt
            end
            # ★ D1b — DISTURBANCE-tally NSTORE/PNN (estab.f:545-559). The existing sub-REGNBK stock — dominated
            # by the stump/root-sprout cohort (esuckr!, AS/PB) just created this cycle — SUPPRESSES new AUTOES regen:
            # NSTORE(pt)=INT((PLPROB·DUPNPT)/(FTEMP·300)+0.5) with FTEMP=logistic(ESB1 PN) and PLPROB·DUPNPT = the
            # point's DBH<2.999 TPA × NPTIDS; PNN(pt)=UNCLAMPED ESA (estab.f:322,554). ie_autoes_run then books only
            # NEWTPP=max(0,ITPP−NSTORE) full-prob trees + the NSTORE "already there" trees at the increment PROB1−PNN.
            # Runs ONLY on the fresh disturbance tally (NTALLY==1, INADV==0, within 20yr) — same gate as the ESB block;
            # the ingrowth (NTALLY==99) NSTORE is handled per-point in ie_autoes_run. Without it the sprout cohort did
            # not suppress AUTOES ⇒ WL/DF/LP over-produced ~9× (flagship 2978686010690: jl 790 vs oracle 84 TPA regen).
            if _ntally == 1
                ftemp1 = 1f0 / (1f0 + exp(-esb1))
                esa_raw = 1f0 / (1f0 + exp(-(-5.17397f0 + 0.85131f0 * log(tpacre))))   # UNCLAMPED ESA → PNN
                psmall = zeros(Float32, nptids)
                @inbounds for i in 1:s.trees.n
                    s.trees.dbh[i] < 2.999f0 || continue
                    pid = Int(s.trees.plot_id[i]); (1 <= pid <= nptids) && (psmall[pid] += s.trees.tpa[i])
                end
                length(est.es_nstore) == dupnpt_i || (est.es_nstore = zeros(Int32, dupnpt_i))
                length(est.es_pnn)    == dupnpt_i || (est.es_pnn    = zeros(Float32, dupnpt_i))
                fill!(est.es_nstore, Int32(0)); fill!(est.es_pnn, esa_raw)
                @inbounds for pt in 1:nptids
                    ns_pt = floor(Int32, psmall[pt] * Float32(nptids) / (ftemp1 * 300f0) + 0.5f0)
                    base = (pt - 1) * idup
                    for k in 1:idup; (base + k) <= dupnpt_i && (est.es_nstore[base + k] = ns_pt); end
                end
            end
        end
        esb_shift = est.esb_shift
    end
    # #143: existing small-tree (DBH<REGNBK=2.999) stocking → the INGROWTH NSTORE (estab.f:305-315,589). Current
    # (self-thinned) DBH<2.999 TPA at the ingrowth cycle; ie_autoes_run turns it into per-plot NSTORE via prob1.
    # #143 ingrowth NSTORE is PER-INVENTORY-POINT, not uniform (measured FVSem_g16 estab.f:313 PLPROB(N=ITRE(I))):
    # the DBH<REGNBK small-tree stock accumulates per point, so points with no small trees leave their IDUP-plot
    # block FREE to establish (NEWTPP=ITPP), while the stocked points suppress. Uniform fill over-suppressed ~5×
    # (em_474…: live 127 vs jl 23). point_small[pt] = the point's DBH<2.999 TPA (t.plot_id = point index 1..nptids).
    tpacre_ingro = 0f0
    point_small = zeros(Float32, nptids)
    if is_ingro
        @inbounds for i in 1:s.trees.n
            s.trees.dbh[i] < 2.999f0 || continue
            tpacre_ingro += s.trees.tpa[i]
            pid = Int(s.trees.plot_id[i])
            (1 <= pid <= nptids) && (point_small[pid] += s.trees.tpa[i])
        end
    end
    # MECHPREP/BURNPREP site prep (estb/esetpr.f): on the DISTURBANCE tally ONLY — estab.f:224-230 cancels prep
    # for NTALLY>1, and the ingrowth path forces IPREP=1. Gather the scheduled prep %-of-plots at the disturbance
    # date → normalized SUMUP for ie_autoes_tally's per-plot IPPREP sampler. No keyword ⇒ prep_sumup=nothing ⇒
    # every plot IPREP=1 (byte-identical to the pre-wire behaviour).
    prep_sumup = nothing
    if _ntally == 1 && !is_ingro
        pmech_pct = nothing; pburn_pct = nothing
        for a in s.control.schedule
            (a.icflag == Int32(493) || a.icflag == Int32(491)) || continue
            ay = Int(a.year)
            idt = (0 < ay < 1000) ? (ay == icyc ? year : -1) : ay
            (year <= idt < next_year) || continue
            a.icflag == Int32(493) ? (pmech_pct = a.params[2]) : (pburn_pct = a.params[2])
        end
        if pmech_pct !== nothing || pburn_pct !== nothing
            es = ie_esetpr(pmech_pct, pburn_pct)
            prep_sumup = ie_esetpr_normalize(0f0, es.pmech, es.pburn, es.ialn2, es.ialn3)
        elseif s.variant isa InlandEmpire
            # estab.f:365-370 — the user supplied NO site-prep keyword ⇒ ESPREP DEFAULT proportions by habitat
            # series drive the per-plot IPPREP (MEASURED FVSie_g16 na_def: NONE/MECH/BURN ≈ 0.48/0.26/0.26). This
            # is NOT inert: ie_estock's per-IPREP SPRE term then splits the per-plot PROB1 (0.699/0.663977/0.657729
            # by IPREP), and the 48/26/26 mean 0.6792 = the oracle stocking prob — closing the bare-IE dated-ESTAB
            # baseline (jl 629→610 internal, ×0.909 = 555). Only the DISTURBANCE tally (NTALLY==1, non-ingrowth);
            # the ingrowth path forces IPREP=1 (estab.f). Species/series from ie_estab_indices; topo/BA/elev = the
            # same stand values fed to ie_estock below.
            _dx = ie_estab_indices(ihab_code, Int(p.user_forest_code))
            pn0, pm0, pb0 = ie_esprep(_dx.iser, es_aspect, es_slope, baaa, p.elevation)
            prep_sumup = ie_esetpr_normalize(pn0, pm0, pb0, 0, 0)
        end
    end
    # ★#143 follow-on: the END-OF-CYCLE stocking PN (estab.f:572 ESTOCK(BAA=BAAA(NNID))) uses the POST-growth
    # per-point BA. jl's s.density was last refreshed PRE-growth (simulate.jl:548, before diameter/height growth),
    # so point_ba lagged one cycle — MEASURED jl 62.51 vs live end-cycle BAAA 77.51 at ic=1 ⇒ PN too low ⇒ PROB1
    # 0.463 vs 0.617. Refresh here so the ingrowth stocking PN sees the grown BA (esb_shift above already captured
    # the PRE-growth baaa as its inventory BAAOLD). Ingrowth-only: the disturbance re-stocking path uses the bare
    # per-point BA deliberately (validated bit-exact) and is inert to a refresh (bare stays bare).
    baaa_pn = baaa
    # ★ OVERSTORY-ONLY per-point BAAA (estb/dense.f:212 `IF(D.LT.REGNBK) GO TO 10`): the end-of-cycle stocking
    # ESTOCK reads BAAA(NNID), which dense.f accumulates from OVERSTORY trees ONLY (D≥REGNBK=2.999, blkdat.f:235).
    # jl fed the ALL-tree point_ba (point_basal_area!, no DBH filter), which on a seedling/regenerating stand counts
    # the sub-REGNBK cohort and inflates BAAA far above the oracle's ~0 (MEASURED FVSie_g16 196386885020004 @2022:
    # oracle BAAA=1.0 all 4 points [all stems <2.999] vs jl point_ba=[13.3,1.1,15.7,5.0] ⇒ PN +0.35 ⇒ PROB1 0.72 vs
    # 0.636 ⇒ NSTORE too low ⇒ AUTOES ingrowth over-books ~+66 TPA). Rebuild the per-point BA counting only D≥REGNBK
    # (same tpa·0.005454·D²·PI/GROSPC scale as point_basal_area!), matching dense.f. On a stocked stand this equals
    # the all-tree value minus the small-tree BA — MEASURED oracle 30194821 @2018 BAAA=59.5 vs jl overstory 60.3
    # (all-tree 68.9) ⇒ closer to the oracle, no regression. Used for the scalar baaa_pn AND the per-point PROB1
    # vector (prob1_pt) below. Ingrowth-only (the disturbance re-stocking path deliberately uses the bare all-tree
    # per-point BA, validated bit-exact, and is not refreshed here).
    # ★ CONTINUATION-tally per-point BAAA (dense.f:212-213). FVS builds BAAA(NNID) = per-point OVERSTORY (D≥REGNBK)
    # POST-growth BA for EVERY establishment tally, not just ingrowth. A CONTINUATION tally (NTALLY≥2 — e.g. the
    # bare-INADV cycle-2 re-tally on a planted cohort) previously fell back to the STALE, ALL-tree scalar
    # s.density.point_ba[1] (measured ≈0.4/pt vs the oracle's overstory 1.6-9.3/pt at cycle 2) ⇒ ESPADV BAA clamped
    # to ~1 vs the oracle's 1.6-3.9 ⇒ GF adv-prob 0.168 vs 0.398 ⇒ DF over-booked +53 @2012, compounding through
    # faithful per-species growth+mortality to +212 TPA @2032. Build the SAME post-growth overstory per-point BAAA
    # the ingrowth path uses (from the current post-growth treelist, D≥REGNBK, same tpa·0.005454·D²·PI/GROSPC scale
    # as dense.f:213) and feed it into the species-prob tables. Gated to NTALLY≥2 so the validated NTALLY=1 initial
    # disturbance (all-tree ≈ overstory on its overstory-dominated residuals) stays byte-identical.
    _ie_cont = (s.variant isa InlandEmpire) && !is_ingro && _ntally >= 2
    overstory_pba = Float32[]
    if is_ingro || _ie_cont
        is_ingro && compute_density!(s)                                   # ingrowth also refreshes PROB1/ESB point stats
        _sc = s.plot.pi / s.plot.gross_space
        overstory_pba = zeros(Float32, nptids)
        @inbounds for i in 1:s.trees.n
            s.trees.dbh[i] < 2.999f0 && continue                          # REGNBK overstory filter (dense.f:212)
            pid = Int(s.trees.plot_id[i]); (1 <= pid <= nptids) || continue
            overstory_pba[pid] += s.trees.tpa[i] * 0.005454154f0 * s.trees.dbh[i]^2 * _sc
        end
        baaa_pn = isempty(overstory_pba) ? baaa : overstory_pba[1]
    end
    # ★ OVER(i) per-species overstory BA (dense.f:214 OVER(ISPC,IP)). The ESPADV/ESPXCS/ESPSUB seed-source bump
    # (OVER(i)>9.95 ⇒ +2.48 LP / +1.05 GF / … logit) is keyed off the TALLY POINT's per-species D≥REGNBK BA.
    # jl formerly passed over=0 ⇒ the bump never fired ⇒ overstory-seed-source species (LP/GF) were massively
    # under-picked (MEASURED FVSie_g16 1856050746290487 @2043: oracle OVER=[GF 11.22, LP 15.60] ⇒ PADV(LP) 0.391,
    # PXCS(LP) 0.301; jl over=0 gave PADV(LP) 0.051, PXCS(LP) 0.031 ⇒ ~1 LP vs oracle 28, +34 AF vs 13). Build the
    # POINT-1 per-species overstory BA (same tpa·0.005454·D²·PI/GROSPC scale as baaa_pn), matching dense.f's OVER,
    # and pass it to ie_autoes_run. Species index = the IE variant index 1-10 (= estab order WP WL DF GF WH RC LP
    # ES AF PP). IE only (EM keeps over=0 pending its own oracle validation). Point 1 is the tally point the port's
    # single ba/occ/species-prob set already uses; single-point IE stands (this bug's regime) are exact.
    # ★ PER-POINT per-species overstory BA `over_pt` (10×nptids): the oracle runs ESPADV/ESPXCS/ESPSUB per
    # inventory point with THAT point's BAAA (clamp[1,400]) and OVER(ISPC,IP), so a multi-point stand's low-BA
    # points regen the BA²-penalized species (WP/DF/GF) while its dense point regens redcedar — the oracle cohort
    # MIXES them. jl formerly used only point-1 (`over_sp`) + the scalar baa for ALL plots ⇒ mono-RC (RC has the
    # lowest mortality ⇒ never culled ⇒ +TPA over-retention). `over_pt[:,pt]` is threaded into ie_autoes_tally's
    # per-point species tables; `over_sp` is kept as the point-1 slice so a single-point stand is byte-identical.
    over_sp = Float32[]; over_pt = Array{Float32}(undef, 0, 0)
    if s.variant isa InlandEmpire
        over_sp = zeros(Float32, 10)
        over_pt = zeros(Float32, 10, nptids)
        _scv = s.plot.pi / s.plot.gross_space
        @inbounds for i in 1:s.trees.n
            s.trees.dbh[i] < 2.999f0 && continue                              # REGNBK overstory filter (dense.f:212)
            pid = Int(s.trees.plot_id[i]); (1 <= pid <= nptids) || continue
            sp = Int(s.trees.species[i]); (1 <= sp <= 10) || continue         # estab species 1-10 carry OVER
            b = s.trees.tpa[i] * 0.005454154f0 * s.trees.dbh[i]^2 * _scv
            over_pt[sp, pid] += b
            pid == 1 && (over_sp[sp] += b)                                    # point-1 slice (byte-identical fallback)
        end
    end
    r = ie_autoes_run(habitat_code = ihab_code, forest_code = Int(p.user_forest_code), nsp = nsp,
                      seed0 = seed0, dupnpt = dupnpt, slo = es_slope, aspect = es_aspect,
                      elev = p.elevation, baa = clamp(baaa_pn, 1f0, 400f0), time = time, esb_shift = esb_shift,
                      # Per-point stocking shift ESB−ESB1(NNID) for the ingrowth per-point PROB1 (fixes the M333
                      # multi-point over-production). Passed only when calibration is active THIS call (scalar
                      # esb_shift≠0 ⇒ the ESB block ran this tally); empty ⇒ ie_autoes_run uses the scalar shift.
                      # ESB−ESB1(NNID) per point: ingrowth AND the IE disturbance tally (both recompute PROB1 per
                      # inventory point). Passed when the ESB block ran this call (esb_shift≠0 ⇒ est.esb_shift_pt filled).
                      esb_shift_pt = (esb_shift != 0f0 && (is_ingro || s.variant isa InlandEmpire)) ? est.esb_shift_pt : Float32[],
                      is_ingro = is_ingro, nstore = est.es_nstore, pnn = est.es_pnn, tpacre_ingro = tpacre_ingro,
                      point_small_tpa = point_small, idup = idup, variant = s.variant,
                      # Per-point slope/aspect (PSLO/PASP) for ESTPP — from the FIA per-plot SLOPE/ASPECT (#143).
                      # Empty (TREEDATA / no per-plot topo) ⇒ ESTPP falls back to the uniform stand slope, inert.
                      point_slope = (isempty(pslo_es) ? Float32[] : @view pslo_es[1:min(nptids, length(pslo_es))]),
                      point_aspect = (isempty(pasp_es) ? Float32[] : @view pasp_es[1:min(nptids, length(pasp_es))]),
                      # Per-inventory-point BAAA(NNID) for the per-point PROB1 stocking logit. INGROWTH: OVERSTORY-ONLY
                      # (D≥REGNBK) post-growth per-point BA (dense.f:212), built above. DISTURBANCE (IE, esb_shift active):
                      # the ALL-tree per-point BA s.density.point_ba (the per-point analog of the validated scalar
                      # baaa=point_ba[1]; the oracle's post-thin BAAA(NNID)=[0,0,162,0] on 3291804010690 matches the
                      # all-tree per-point BA — these points are bare/stocked with negligible sub-REGNBK stock). Empty ⇒
                      # ie_autoes_run keeps the scalar (single-point / non-IE / esb_shift=0 ⇒ byte-identical).
                      point_ba = (is_ingro || _ie_cont) ?
                          (isempty(overstory_pba) ? Float32[] : @view overstory_pba[1:min(nptids, length(overstory_pba))]) :
                          ((s.variant isa InlandEmpire && esb_shift != 0f0 && !isempty(s.density.point_ba)) ?
                              (@view s.density.point_ba[1:min(nptids, length(s.density.point_ba))]) : Float32[]),
                      cont_ppba = _ie_cont,   # IE continuation (NTALLY≥2): use per-point overstory BAAA species tables
                      over_sp = over_sp,        # per-species overstory BA (D≥REGNBK) at point 1 (dense.f OVER); empty ⇒ over=0
                      over_pt = over_pt,        # per-species PER-POINT overstory BA (10×nptids); empty ⇒ scalar/point-1 fallback
                      stoadj = est.stoadj,      # STOCKADJ keyword multiplier (default 1.0 ⇒ inert)
                      spec_mult = est.spec_mult,  # SPECMULT per-species XESMLT (empty ⇒ inert)
                      prep_sumup = prep_sumup,  # MECHPREP/BURNPREP per-plot IPPREP (nothing ⇒ all IPREP=1, inert)
                      # PASSALL/PASMAX excess cap (estab.f:1288). Gated to IE ONLY (EM's cap deferred — task scope);
                      # the AUTOES ingrowth path (NTALLY==99) hardcodes PASMAX=15 (estab.f:249), else the CONFID/
                      # PASSALL value (default 5). Every non-IE caller stays Inf ⇒ dead code ⇒ byte-identical.
                      pasmax = (s.variant isa InlandEmpire) ? (is_ingro ? 15f0 : Float32(est.pasmax)) : Inf32,
                      # Establishment-cohort height-class / WK4 distribution (IE emit path). gentim=FINT−5;
                      # call_espadv = internal-NTALLY==1 (ESPADV CALLed ⇒ advance regen possible): true for a
                      # fresh disturbance (NTALLY=1) AND the ingrowth path (NTALLY 99→INGRO=1,NTALLY=1, estab.f:256);
                      # false only for a continuation (NTALLY≥2) ⇒ PADV=0 ⇒ all-subsequent.
                      gentim = max(fint - 5f0, 0f0), call_espadv = (is_ingro || _ntally == 1),
                      plant_sp = _plant_sp)

    haskey(ENV, "FVSJL_AUTOES_DEBUG") &&
        println(stderr, "AUTOES_IN icyc=$icyc ntally=$(_ntally) seed0=$seed0 es_stream=$(Int(round(est.es_stream))) baaa=$(round(baaa,digits=2)) baa_used=$(round(max(baaa,1f0),digits=2)) time=$time  → total=$(round(sum(r.tally),digits=1))")
    t = s.trees
    xmin = _IE_ES_XMIN
    # Birth-cycle HTG multiplier WK4=HTIMLT (estab.f:801-902/1054-1063). For AUTOES NATURAL regen the per-species
    # WK4=STOMLT(sp)=min(TRAGE,GENTIM)/(GENTIM+1e-4) where GENTIM=FINT-5 and TRAGE is the ESADVH advance-regen age
    # TRAGE=3−DELAY (esadvh.f:86, DELAY=0 at a cycle boundary ⇒ TRAGE=3). So the ADVANCE-regen (dominant) WK4 =
    # min(3,GENTIM)/(GENTIM+1e-4) ⇒ 0.60 at FINT=10 (NOT the PLANT default TRAGE=2 ⇒ 0.40, which under-set the
    # advance cohort). The per-record 0.2/0.0 tail (ESSUBH-subsequent TRAGE=TIME−DELAY + ESXCSH-excess, leftover-
    # STOMLT estab.f:937) is now FULLY PORTED via the ie_autoes_tally emit path (DO 99/33/228), NOT cornered/collapsed.
    # ★ RE-BALANCE MEASUREMENT (2026-09-04, branch fix-ie-cohort-rebalance): the b6edfd9a partition was suspected of
    # over-weighting the slow (WK4≤0.20) classes → under-projection. Instrumenting FVSie_g16 estab.f DO 33/228 (booked
    # WK4+PROB per tree) vs the jl emit records across a POPULATION of bare+sparse establishment stands (ecoregions
    # 331Aa/Ac/Af, 342Ia/b; single- AND multi-point) shows the jl per-class booked TPA MATCHES the oracle within
    # realization noise (advance/subsequent/excess fractions within ~1-3pp; population-mean advance fraction ~equal, jl
    # if anything marginally HIGHER) — the partition is faithful, it does NOT over-weight the slow classes. The residual
    # .sum BA/CCF/SDI cap is a two-sided #206 OLDRN/ZRAND realization straddle (86 stands: BA 38 over / 44 under, mean
    # −0.85; adjacent near-identical stands flip sign; the per-stand advance-fraction discrepancy is UNCORRELATED with
    # the BA sign — jl booking LESS advance yields both +19 and −15). So NO re-balance is warranted (it would break
    # faithfulness to the FVS source without touching the straddle); the residual is legitimately the IE #206 corner.
    # ★ MEASURED (FVSie_g16 esgent trace, 195384161020004): WK4 is applied TWICE — once in REGENT(LESTB) subcycling
    # (regent.f:596 H2=H1+EXP(HTGRL)·SCALE·XRHGRO·WK4) and again in ESGENT (esgent.f:57 HTG=HTG·WK4) ⇒ effective
    # birth multiplier = WK4². ie_esgent! reproduces the square; here we store the single WK4 (per-tree) it reads.
    # TRAGE=3 for IE AND EM: ESADVH advance-regen age TRAGE=3−DELAY (esadvh.f:86, DELAY=0 at a cycle boundary) is
    # variant-general, so the DOMINANT advance cohort's WK4=min(3,GENTIM)/(GENTIM+1e-4)=0.60 (FINT=10), not the
    # PLANT default 0.40. MEASURED directly from the FVSem_g16 esgent trace (stand 11847845010690, esgent.f:57
    # HTG=HTG·WK4 per record): the EM AUTOES cohort books WK4 ∈ {0.60 (advance, 110/138 records), 0.20 (subsequent,
    # 23), 0.00 (excess, 5)} — the advance 0.60 dominates (~80%). jl's collapsed single-record path books that
    # dominant advance class, so it needs 0.60 (EM's em_esgent! applies WK4 ONCE, ie_esgent! squares it — the store
    # is the single per-tree WK4 either way). EM previously kept TRAGE=2 (0.40), UNVALIDATED (left byte-unchanged
    # when the IE emit path landed); 0.60 is the faithful advance value. The 0.20/0.00 subsequent/excess tail
    # (28/138 records) is not separately booked in the collapsed path (that needs the full IE-style emit
    # distribution) — a second-order effect. Measured .sum impact on FIA stands is below integer-BA resolution
    # (the AUTOES ingrowth cohort is a small fraction of tree TPA — ~100 tiny seedlings vs thousands), so this is
    # a faithfulness alignment to the oracle-measured WK4, not a visible .sum mover on the sampled EM stands.
    _autoes_gentim = max(fint - 5f0, 0f0)
    _autoes_trage = (s.variant isa InlandEmpire || s.variant isa EasternMontana) ? 3f0 : 2f0
    _autoes_htimlt = min(_autoes_trage, _autoes_gentim) / (_autoes_gentim + 0.0001f0)
    created = false
    npt_c = size(r.tally_pt, 2)                          # inventory points; established TPA is split per point so
    # Booking table: (sp, pt, height, wk4) → summed TPA. IE uses the faithful per-tree height-class / WK4 emit
    # records (ie_autoes_tally emit path: advance WK4≈0.60 / subsequent 0.20/0.00 first trees + tripled excess),
    # aggregated by identical (sp,pt,height,wk4) — physically equivalent to booking each FVS record separately
    # (same growth/mortality) but keeps the tree list bounded. Other variants (EM) keep the collapsed path.
    book = Vector{NTuple{5,Float32}}()                   # (sp, pt, height, dbh, tpa) with per-record wk4 in a parallel
    bwk4 = Float32[]                                     # array; each is one aggregated seedling record.
    if r.emit !== nothing
        # Book each FVS DO-33 (best) / DO-228 (excess) record 1:1 — NO height/wk4 merge. FVS emits the full
        # per-plot record set (estab.f:1199-1359: best trees individual, excess aggregated per plot at the
        # species-mean height) and ESGENT/REGENT then grows EACH record with its OWN bounded ZZRAN height-
        # growth draw (regent.f:805 CALL BACHLO per record, HTGR=HTGR1·EXP(ZZRAN·HSIGMA)). ★ MEASURED (stand
        # 22963815010497, FVSie_g16 TREELIST): the oracle books 106 records/tally; a prior collapse that merged
        # them by (sp,pt,round(HT·1e3),round(WK4·1e3)) into ~19 gave each merged mass ONE zzran instead of the
        # per-record draws, compressing the birth-cycle HtG spread and — via TPA weighting — systematically
        # under-growing (2057 BA 35 vs oracle 43, TCuFt 944 vs 1199). Booking 1:1 recovers it: BA 42, TCuFt
        # 1178 (residual = the #206 ZRAND realization straddle, ~1.7%). The start heights + per-species TPA +
        # WK4 distribution are already bit-exact to the oracle, so this only restores the record GRANULARITY.
        for rec in r.emit
            sp = Int(rec[1]); pt = Int(rec[2]); hh = Float32(rec[3]); wk4 = Float32(rec[4]); tpa = rec[5]
            (sp < 1 || sp > nsp || tpa <= 0.0) && continue
            dbh = 0.1f0 + 0.001f0 * hh                    # esgent.f:56 sub-breast-height nominal DBH
            push!(book, (Float32(sp), Float32(pt), hh, dbh, Float32(tpa))); push!(bwk4, wk4)
        end
    else
        # Collapsed single-record path (EM / non-emit): one WK4=STOMLT(advance) record per (species, point).
        for sp in 1:nsp
            Float32(r.tally[sp]) > 0f0 || continue
            hh = xmin[sp] + 0.2f0; dbh = 0.1f0 + 0.001f0 * hh
            for pt in 1:npt_c
                tpa_sp = Float32(r.tally_pt[sp, pt]); tpa_sp > 0f0 || continue
                push!(book, (Float32(sp), Float32(pt), hh, dbh, tpa_sp)); push!(bwk4, _autoes_htimlt)
            end
        end
    end
    for (bi, rec) in enumerate(book)
        sp = Int(rec[1]); pt = Int(rec[2]); hht = rec[3]; dbh = rec[4]; tpa_sp = rec[5]
        n = t.n + 1; n + Int(t.ndead) > length(t.dbh) && break   # leave room for the dead block (t.n+1…t.n+ndead)
        t.n = n
        t.species[n]     = Int32(sp)
        t.dbh[n]         = dbh
        t.height[n]      = hht
        t.tpa[n]         = tpa_sp
        t.plot_id[n]     = Int32(pt)
        # ABIRTH = AGEPL+GENTIM (estab.f:628/707); AUTOES natural regen ⇒ AGEPL=0, GENTIM=FINT−5. Read ONLY by
        # Climate-FVS (BIRTHYR=THISYR−ABIRTH → Leites XDF/XPP/XWL transfer distance; apply_climate_dds! +
        # inlandempire/regent.jl clim_treemult), so byte-identical for climate-off IE. Was birth_age=0 ⇒ BIRTHYR=now
        # ⇒ XRELGR≡1 ⇒ under-grown diameter/volume under CLIMATE (matches oracle ABIRTH=GENTIM=5).
        t.birth_age[n]   = _autoes_gentim
        t.htimlt[n]      = bwk4[bi]           # per-tree WK4=HTIMLT (advance 0.60 / subsequent 0.20/0.00 / excess STOMLT)
        # Crown: the REGENT(LESTB) open-grown crown (regent.f:178) CR=0.89722−0.0000461·PCCF, clamped [0.20,0.90].
        pccf = pt <= length(s.density.point_ccf) ? s.density.point_ccf[pt] :
               (isempty(s.density.point_ccf) ? 0f0 : s.density.point_ccf[1])
        cr = clamp(0.89722f0 - 0.0000461f0 * pccf, 0.20f0, 0.90f0)
        icr0 = floor(Int32, cr * 100f0 + 0.5f0)
        t.crown_pct[n]   = icr0
        t.crown_ratio[n] = Float32(icr0)
        t.norm_ht[n]     = Int32(0)
        t.sort_key[n]    = Float64(n)
        # Newly-established AUTOES trees carry NO volume in their birth cycle. Zero the volume fields so the
        # post-ESTAB .sum does not sum a STALE value inherited from this slot's prior occupant (a dead inventory
        # tree or a previously-removed record — compute_volumes! writes cuft_vol over 1:(n+ndead)). Mirrors the
        # esuckr.f sprout path's "the slot may hold a previously-deleted record" zeroing; without it an all-dead-
        # inventory (HISTORY=6) + heavy-AUTOES stand over-reports TCuFt/MCuFt/BdFt (a 0.1" seedling summing tens
        # of cuft). Correct volume is filled by next cycle's compute_volumes!.
        t.cuft_vol[n]       = 0f0; t.merch_cuft_vol[n] = 0f0
        t.saw_cuft_vol[n]   = 0f0; t.bdft_vol[n]       = 0f0
        created = true
    end
    created || return false
    push!(est.autoes_years_done, Int32(year))
    compute_density!(s)
    return true
end

# Climate-FVS: IE PLANTS symbols (PLNJSP, blkdat.f:186) for the viability-column lookup.
climate_plant_symbols(::InlandEmpire) = _IE_PLNJSP
