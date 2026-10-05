# =============================================================================
# organon_diamgro.jl — OC (Oregon Coast) ORGANON SWO diameter growth (chunk C3).
#
# Ported from (ORGANON edition SWO / VERSION=1 path only):
#   • organon/diagro.f   — DIAMGRO_RUN (per-tree DG driver), DG_SWO (the ln(DDS) diameter-growth
#                          equation), DG_THIN / DG_FERT (thinning / fertilizer adjustments),
#                          GET_BAL_RUN (basal-area-in-larger lookup).
#   • organon/statsorg.f — SSTATS (per-cycle stand SBA / BAL / BALL / CCFL / CCFLL / TPA / SCCF).
#   • organon/diamcal.f  — GET_BAL, DIB_SWO / DOB_SWO (bark inside/outside), and the DGCALIB
#                          radial-increment diameter calibration (RAD path; see note below).
#   • organon/submax.f   — SUBMAX (maximum size-density line A1/A2; consumed by mortality C6 and
#                          the stand density index — ported here as the task lists it, inert on DG).
#   • oc/dgdriv.f:444-453 — the FVS-side copy-back seam (DGRO → BARK → DIAGR → DDS → WK2). The
#                          bark ratio there is FVS-native BRATIO (shared engine), NOT ORGANON DIB/
#                          DOB; the ORGANON DIB/DOB bark only feeds the DGCALIB RAD path.
#
# SCOPE: this ports the ORGANON diameter-growth core — the deterministic (DGSD=0) per-tree DGRO
# that overwrites the FVS Wykoff LN(DDS) for the 18 valid ORGANON species. The FVS serial-corr is
# suppressed on the ORGANON path (OLDRN=0, FRM=1; oc/dgdriv.f:549-550), so DGRO is a hard bit-exact
# target with no RNG straddle. VERSION=1 (SWO) ONLY; NWO(2)/SMC(3)/RAP(4) (OP=Olympic) raise a
# clear error — their coefficient branches sit alongside SWO in the same Fortran (thin OP follow-on).
#
# DGCALIB / RAD deferral (from C2), now lifted: on FVS/FIA inventory there are no radial-increment
# cores, so RAD=.FALSE. and TMPCAL(3,*)=1.0 (organon/prepare.f:137-140). `organon_dgcalib_swo`
# returns the all-1.0 CALIB(3,*) directly for RAD=.FALSE. and — now that DG_SWO / SSUM / bark / BAL
# exist — carries the faithful RAD=.TRUE. computation as well (exercised only if radial cores appear).
#
# MEASURED bit-exact vs the live FVSoc_clean oracle (scoped `DEBUG 1 / DGDRIV`, stand S248112 /
# ocmin, 27 records; 17 valid ORGANON trees) — the dgdriv `I,ISPC,DBH,DGRO,BARK,DIAGR,DDS=` growth
# dump. See docs/OC_VARIANT_PORT_AUDIT.md (C3 verdict).
# =============================================================================

# --- DG_SWO DGPAR(18,11): B0,B1,B2,B3,B4,B5,B6,K1,K2,K3,K4 (species groups 1..18) --------------
# organon/diagro.f DATA DGPAR (column-major: each block = 18 group values for one parameter).
const OC_DG_SWO_B0 = Float32[-5.35558894,-5.84904111,-4.51958940,-4.12342552,-2.08551255,-5.70052255,-11.45456097,-9.15835863,-8.84531757,-7.78451344,-3.36821750,-3.59333060,-3.41449922,-7.81267986,-4.43438109,-4.39082007,-8.08352683,-8.08352683]
const OC_DG_SWO_B1 = Float32[0.840528547,1.668196109,0.813998712,0.734988422,0.596043703,0.865087036,0.784133664,1.0,1.5,1.2,1.2,1.2,1.0,1.405616529,0.930930363,1.0,1.0,1.0]
const OC_DG_SWO_B2 = Float32[-0.0427481848,-0.0853271265,-0.0493858858,-0.0425469735,-0.0215223077,-0.0432543518,-0.0261377888,-0.00000035,-0.0006,-0.07,-0.07,-0.07,-0.05,-0.0603105850,-0.0465947242,-0.0945057147,-0.00000035,-0.00000035]
const OC_DG_SWO_B3 = Float32[1.15950313,1.21222176,1.10249641,1.05942163,1.02734556,1.10859727,0.70174783,1.16688474,0.51225596,0.0,0.0,0.51637418,0.0,0.64286007,0.0,1.06867026,0.31176647,0.31176647]
const OC_DG_SWO_B4 = Float32[0.954711126,0.679346647,0.879440023,0.808656390,0.383450822,0.977332597,2.057236260,0.0,0.418129153,1.01436101,0.0,0.0,0.324349277,1.037687142,0.510717175,0.685908029,0.0,0.0]
const OC_DG_SWO_B5 = Float32[-0.00894779670,-0.00809965733,-0.0108521667,-0.0107837565,-0.00489046624,0.0,-0.00415440257,0.0,-0.00355254593,-0.00834323811,0.0,0.0,0.0,0.0,0.0,-0.00586331028,0.0,0.0]
const OC_DG_SWO_B6 = Float32[0.0,0.0,-0.0333706948,0.0,-0.0609024782,-0.0526263229,0.0,-0.02,-0.0321315389,0.0,-0.0339813575,-0.02,-0.0989519477,-0.0787012218,-0.0688832423,0.0,-0.0730788052,-0.0730788052]
const OC_DG_SWO_K1 = Float32[5.0,5.0,5.0,5.0,5.0,5.0,5.0,4000.0,110.0,10.0,10.0,10.0,10.0,5.0,5.0,5.0,4000.0,4000.0]
const OC_DG_SWO_K2 = Float32[1.0,1.0,1.0,1.0,1.0,1.0,1.0,4.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,4.0,4.0]
const OC_DG_SWO_K3 = Float32[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
const OC_DG_SWO_K4 = Float32[2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7,2.7]

# organon/diamcal.f DIB_SWO / DOB_SWO DIBPAR(18,2): B0,B1 (bark inside/outside; DGCALIB RAD path).
const OC_DIB_SWO_B0 = Float32[0.903563,0.904973,0.809427,0.859045,0.837291,0.933707,0.9497,0.97,0.96317,0.94448,0.859151,0.910499,0.97059,0.878457,0.889703,0.947,0.94448,0.94448]
const OC_DIB_SWO_B1 = Float32[0.989388,1.0,1.016866,1.0,1.0,1.0,1.0,1.0,1.0,0.9875170,1.0178109,1.01475,0.993585,1.02393,1.0104062,1.0,0.9875170,0.9875170]

"organon/diagro.f DG_SWO ADJ — species-group full adjustment multiplier (VERSION=1)."
@inline function oc_dg_swo_adj(g::Int)
    g == 1  && return 0.8938f0
    g == 2  && return 0.8722f0
    g == 4  && return 0.7903f0
    g == 9  && return 0.7928f0
    g == 10 && return 0.7259f0
    g == 14 && return 1.0f0
    g == 15 && return 0.7667f0
    return 0.8f0
end

"""
    oc_dg_swo(g, dbh, cr, site, sbal1, sba1) -> DG

organon/diagro.f DG_SWO — 5-yr ln(DDS) diameter growth for species group `g` (SWO). `site` is
SI_1 = SITE_1 − 4.5 (organon/execute2.f:322; DIAMGRO_RUN CASE(1) SITE=SI_1). `sbal1` = basal area
in larger trees (GET_BAL), `sba1` = stand basal area. All REAL*4 via fexp/flog/fpow (doctrine #8).
"""
@inline function oc_dg_swo(g::Int, dbh::Float32, cr::Float32, site::Float32,
                           sbal1::Float32, sba1::Float32)
    b0 = OC_DG_SWO_B0[g]; b1 = OC_DG_SWO_B1[g]; b2 = OC_DG_SWO_B2[g]; b3 = OC_DG_SWO_B3[g]
    b4 = OC_DG_SWO_B4[g]; b5 = OC_DG_SWO_B5[g]; b6 = OC_DG_SWO_B6[g]
    k1 = OC_DG_SWO_K1[g]; k2 = OC_DG_SWO_K2[g]; k3 = OC_DG_SWO_K3[g]; k4 = OC_DG_SWO_K4[g]
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
    return fexp(lndg)*cradj*oc_dg_swo_adj(g)
end

"organon/diamcal.f GET_BAL — basal-area-in-larger lookup at DBH `d` (BALL 1..51, BAL 1..500)."
@inline function oc_get_bal(d::Float32, ball::Vector{Float32}, bal::Vector{Float32})
    if d > 100.0f0
        return 0.0f0
    elseif d > 50.0f0
        k = trunc(Int, d - 49.0f0)
        return ball[k]
    else
        k = trunc(Int, d*10.0f0 + 0.5f0)
        return bal[k]
    end
end

"""
    oc_sstats(spgrp, dbh, ht, expan, ntrees) -> (sba, bal, ball, ccfl, ccfll, tpa, sccf)

organon/statsorg.f SSTATS — stand-level SBA / BAL(500) / BALL(51) / CCFL(500) / CCFLL(51) / TPA /
SCCF from the CURRENT (`TDATAR`) tree list. NPTS=1 in the FVS-coupled path (EXPAN already per-acre),
so the `/FLOAT(NPTS)` divisions are the identity. `MCW**2` and `DBH**2` use integer powers (exact
products), matching the Fortran.
"""
function oc_sstats(spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32},
                   expan::Vector{Float32}, ntrees::Int)
    bal = zeros(Float32, 500); ball = zeros(Float32, 51)
    ccfl = zeros(Float32, 500); ccfll = zeros(Float32, 51)
    sba = 0.0f0; tpa = 0.0f0; sccf = 0.0f0
    @inbounds for i in 1:ntrees
        expan[i] < 0.0001f0 && continue
        g = Int(spgrp[i]); d = dbh[i]; h = ht[i]; ex = expan[i]
        ba = d*d*ex*0.005454154f0
        sba += ba
        tpa += ex
        mcw = oc_mcw_swo(g, d, h)
        ccf = 0.001803f0*(mcw*mcw)*ex
        sccf += ccf
        if d > 50.0f0
            l = trunc(Int, d - 49.0f0); l > 52 && (l = 52)
            for k in 1:500; ccfl[k] += ccf; bal[k] += ba; end
            for k in 1:(l-1); ccfll[k] += ccf; ball[k] += ba; end
        else
            l = trunc(Int, d*10.0f0 + 0.5f0)
            for k in 1:(l-1); ccfl[k] += ccf; bal[k] += ba; end
        end
    end
    return sba, bal, ball, ccfl, ccfll, tpa, sccf
end

"organon/diagro.f DG_THIN — thinning adjustment (VERSION=1). Returns 1.0 with no prior thinning."
function oc_dg_thin(isp::Int32, cyclg::Int, babt::Float32, bart::Vector{Float32}, yt::Vector{Float32})
    if isp == 263
        pt1 = 0.723095045f0; pt2 = 1.0f0; pt3 = -0.2644085320f0
    else # 202 and default (VERSION=1)
        pt1 = 0.6203827985f0; pt2 = 1.0f0; pt3 = -0.2644085320f0
    end
    xtime = Float32(cyclg)*5.0f0
    thinx1 = 0.0f0
    @inbounds for i in 2:5
        thinx1 += bart[i]*fexp((pt3/pt2)*(yt[1]-yt[i]))
    end
    thinx2 = thinx1 + bart[1]
    thinx3 = thinx1 + babt
    prem = thinx3 <= 0.0f0 ? 0.0f0 : thinx2/thinx3
    prem > 0.75f0 && (prem = 0.75f0)
    return 1.0f0 + pt1*fpow(prem, pt2)*fexp(pt3*(xtime - yt[1]))
end

"organon/diagro.f DG_FERT — fertilizer adjustment (VERSION<=3). Returns 1.0 with no fertilizer."
function oc_dg_fert(isp::Int32, cyclg::Int, si_1::Float32, pn::Vector{Float32}, yf::Vector{Float32})
    if isp == 202
        pf1 = 1.368661121f0; pf2 = 0.741476964f0; pf3 = -0.214741684f0; pf4 = -0.851736558f0; pf5 = 2.0f0
    else # 263 and default
        pf1 = 0.0f0; pf2 = 1.0f0; pf3 = 0.0f0; pf4 = 0.0f0; pf5 = 1.0f0
    end
    faldwn = 1.0f0
    xtime = Float32(cyclg)*5.0f0
    fertx1 = 0.0f0
    @inbounds for i in 2:5
        fertx1 += (pn[i]/800.0f0)*fexp((pf3/pf2)*(yf[1]-yf[i]))
    end
    fertx2 = fexp(pf3*(xtime - yf[1]) + pf4*fpow(si_1/100.0f0, pf5))
    return 1.0f0 + (pf1*fpow((pn[1]/800.0f0) + fertx1, pf2)*fertx2)*faldwn
end

"""
    oc_diamgro_run(isp, ispgrp, dbh, cr, si_1, sbal1, sba1, calib3;
                   cyclg, pn, yf, babt, bart, yt) -> DGRO

organon/diagro.f DIAMGRO_RUN (VERSION=1) — DGRO = DG_SWO · CALIB(3,g) · FERTADJ · THINADJ. On cyc0
inventory with no fert/thin, FERTADJ=THINADJ=1.0 and (CALD=.FALSE.) CALIB(3,g)=1.0.
"""
function oc_diamgro_run(isp::Int32, ispgrp::Int, dbh::Float32, cr::Float32, si_1::Float32,
        sbal1::Float32, sba1::Float32, calib3::Float32; cyclg::Int=0,
        pn::Vector{Float32}=zeros(Float32,5), yf::Vector{Float32}=zeros(Float32,5),
        babt::Float32=0.0f0, bart::Vector{Float32}=zeros(Float32,5), yt::Vector{Float32}=zeros(Float32,5))
    dg = oc_dg_swo(ispgrp, dbh, cr, si_1, sbal1, sba1)
    fertadj = oc_dg_fert(isp, cyclg, si_1, pn, yf)
    thinadj = oc_dg_thin(isp, cyclg, babt, bart, yt)
    return dg*calib3*fertadj*thinadj
end

"""
    oc_submax(spgrp, dbh, expan, ntrees, msdi_1, msdi_2, msdi_3) -> (a1, a2)

organon/submax.f SUBMAX (VERSION=1, TRIAL=.FALSE.) — the maximum size-density line. A2=0.62305
(Reineke); A1 = TEMPA1·A1MOD from the DF/TF/PP basal-area proportions of the big-6 groups (1..3).
Consumed by ORGANON mortality (C6) and the RD stand-density index; inert on diameter growth.
"""
function oc_submax(spgrp::Vector{Int32}, dbh::Vector{Float32}, expan::Vector{Float32},
                   ntrees::Int, msdi_1::Float32, msdi_2::Float32, msdi_3::Float32)
    a2 = 0.62305f0
    kb = 0.005454154f0
    tempa1 = msdi_1 > 0.0f0 ? (flog(10.0f0) + a2*flog(msdi_1)) : 6.21113f0
    bagrp = zeros(Float32, 18)
    @inbounds for i in 1:ntrees
        g = Int(spgrp[i])
        bagrp[g] += kb*(dbh[i]*dbh[i])*expan[i]      # EX1 = TDATAR(I,4) (TRIAL=.FALSE.)
    end
    totba = bagrp[1] + bagrp[2] + bagrp[3]
    pdf = ptf = ppp = 0.0f0
    if totba > 0.0f0
        pdf = bagrp[1]/totba
        ptf = bagrp[2]/totba
    end
    tfmod = msdi_2 > 0.0f0 ? (flog(10.0f0) + a2*flog(msdi_2))/tempa1 : 1.03481817f0
    ocmod = msdi_3 > 0.0f0 ? (flog(10.0f0) + a2*flog(msdi_3))/tempa1 : 0.9943501f0
    totba > 0.0f0 && (ppp = bagrp[3]/totba)
    if pdf >= 0.5f0
        a1mod = 1.0f0
    elseif ptf >= 0.6666667f0
        a1mod = tfmod
    elseif ppp >= 0.6666667f0
        a1mod = ocmod
    else
        a1mod = pdf + tfmod*ptf + ocmod*ppp
    end
    a1mod <= 0.0f0 && (a1mod = 1.0f0)
    return tempa1*a1mod, a2
end

"organon/diamcal.f DIB_SWO — diameter inside bark (DGCALIB RAD path)."
@inline oc_dib_swo(g::Int, dob::Float32) = OC_DIB_SWO_B0[g]*fpow(dob, OC_DIB_SWO_B1[g])
"organon/diamcal.f DOB_SWO — diameter outside bark (DGCALIB RAD path)."
@inline oc_dob_swo(g::Int, dib::Float32) = fpow(dib/OC_DIB_SWO_B0[g], 1.0f0/OC_DIB_SWO_B1[g])

# organon/diamcal.f DIB_NWO / DOB_NWO DIBPAR(11,2) (VERSION=2; DIB_SMC carries the same table).
const OP_DIB_NWO_B0 = Float32[0.971330,0.904973,0.933707,0.9497,0.97,0.96317,0.97059,0.878457,0.947,0.94448,0.94448]
const OP_DIB_NWO_B1 = Float32[0.966365,1.0,1.0,1.0,1.0,1.0,0.993585,1.02393,1.0,0.987517,0.987517]

"organon/diamcal.f SSUM — BAL(500)/BALL(51)/SBA (N=1: STDATAR start state; N=2: TDATAR). The CCF term reads
TDATAR(I,4) (the CURRENT expansion) in BOTH modes, `CCF=0.001803*MCW**2*TDATAR(I,4)/NPTS`. Only SBA/BAL/BALL
are consumed by DGCALIB (CCFL feeds CRCALIB, which keeps its own pass)."
function organon_ssum_bal(version::Int, spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32},
                          ex_n::Vector{Float32}, ntrees::Int, npts::Int)
    bal = zeros(Float32, 500); ball = zeros(Float32, 51); sba = 0f0
    fnp = Float32(npts)
    @inbounds for i in 1:ntrees
        ex = ex_n[i] / fnp
        ex <= 0f0 && continue
        d = dbh[i]
        ba = ((d * d) * ex) * 0.005454154f0
        sba += ba
        if d > 50f0
            l = trunc(Int, d - 49f0); l > 52 && (l = 52)
            for k in 1:500; bal[k] += ba; end
            for k in 1:(l - 1); ball[k] += ba; end
        else
            l = trunc(Int, d * 10f0 + 0.5f0)
            for k in 1:(l - 1); bal[k] += ba; end
        end
    end
    return sba, bal, ball
end

"""
    organon_dgcalib(version, nspn, npts, ntrees, species, spgrp, dbh, ht, cr, expan, radgro,
                    sba2, bal2, ball2, og, si_1, si_2, calib1) -> calib3::Vector{Float32}(18)

organon/diamcal.f DGCALIB (identical source in the OC and OP builds) — the ORGANON diameter-growth calibration
CALIB(3,·) from the measured radial increments RADGRO = DG/2 ({oc,op}/cratet.f RADGRO(I)=DG(I)/2.0). Called by
PREPARE (prepare.f:133-141) after CRCALIB/PRDCR, so `dbh/ht/cr/expan` are PREPARE's TDATAR after the HT/CR dub and
`sba2/bal2/ball2/og` are CRCALIB's SSUM(2)/OLDGROWTH on that list. `si_1/si_2` are the RAW site indices (DGCALIB
subtracts 4.5), `calib1` = the HDCALIB height ratios (PHTS). With no RADGRO>0 (RAD=.FALSE.) PREPARE skips DGCALIB
and CALIB(3,·)=1.0 (prepare.f:137-140). VERSION 1 = SWO (OC), 2 = NWO (OP).
Steps (diamcal.f:31-170): predict the last-5-yr DG at the current state (DIAMGRO); per group the measured
DG=DBH−DOB(DIB(DBH)−2·RADGRO) and RATIO=ΣDG/ΣPDG (≤5 measured trees ⇒ the pooled AVRAT); back-date DBH (measured
or RATIO·PDG), height (PHTS, calibrated H-D ratio), crown ratio and expansion (MORTAL: PM=1−(1−XPM)^POW ≤ .5); the
SSUM(1) start-state stand; re-predict DG from the back-dated state and CALTST(ΣDG, ΣPDG, ΣDG²/PDG, RADIN).
"""
function organon_dgcalib(version::Int, nspn::Int, npts::Int, ntrees::Int, species::Vector{Int32},
                         spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32}, cr::Vector{Float32},
                         expan::Vector{Float32}, radgro::Vector{Float32}, sba2::Float32, bal2::Vector{Float32},
                         ball2::Vector{Float32}, og::Float32, si_1::Float32, si_2::Float32, calib1::AbstractVector{Float32})
    calib3 = ones(Float32, 18)
    radin = zeros(Int, 18); rad = false
    @inbounds for i in 1:ntrees                        # prepare.f:484-487 (EDIT)
        if radgro[i] > 0f0
            rad = true; radin[spgrp[i]] += 1
        end
    end
    rad || return calib3
    swo = version == 1
    xsi_1 = si_1 - 4.5f0; xsi_2 = si_2 - 4.5f0
    _dib(g, d) = swo ? oc_dib_swo(g, d) : OP_DIB_NWO_B0[g] * fpow(d, OP_DIB_NWO_B1[g])
    _dob(g, d) = swo ? oc_dob_swo(g, d) : fpow(d / OP_DIB_NWO_B0[g], 1f0 / OP_DIB_NWO_B1[g])
    _hd(g, d)  = swo ? oc_a_hd_swo(g, d) : op_hd_nwo(g, d)
    # DIAMGRO (diamcal.f:172-215): SITE = SI_1 (SWO); NWO: SI_2 for WH (FIA 263), else SI_1
    function _pdg(k, d, c, sba, ball, bal)
        sbal = oc_get_bal(d, ball, bal)
        g = Int(spgrp[k])
        if swo
            return oc_dg_swo(g, d, c, xsi_1, sbal, sba)
        else
            return op_dg_nwo(g, d, c, species[k] == 263 ? xsi_2 : xsi_1, sbal, sba)
        end
    end
    pdg = zeros(Float32, ntrees); growth = zeros(Float32, ntrees)
    sd1 = zeros(Float32, ntrees); sd2 = zeros(Float32, ntrees); sd3 = zeros(Float32, ntrees); sd4 = zeros(Float32, ntrees)
    @inbounds for i in 1:ntrees                        # PREDICT PREVIOUS 5-YR DIAMETER
        pdg[i] = _pdg(i, dbh[i], cr[i], sba2, ball2, bal2)
    end
    ratio = ones(Float32, 18); sum3 = 0f0; sum4 = 0f0
    @inbounds for g in 1:18
        ratio[g] = 1f0
        radin[g] <= 0 && continue
        sum1 = 0f0; sum2 = 0f0
        for j in 1:ntrees
            spgrp[j] == g || continue
            dbh2 = dbh[j]
            dibj = _dib(g, dbh2)
            if radgro[j] > 0f0
                dibck = dibj - 2f0 * radgro[j]
                if dibck > 0f0
                    dbh1 = _dob(g, dibck)
                    dbh1 < 0.1f0 && (dbh1 = 0.1f0)
                else
                    dbh1 = 0.1f0
                end
                dg = dbh2 - dbh1
                sum1 += dg; sum2 += pdg[j]
                sd1[j] = dbh2 - dg
                sd1[j] < 0.1f0 && (sd1[j] = 0.1f0)
                growth[j] = dg
            end
        end
        ratio[g] = sum1 / sum2
        sum3 += sum1; sum4 += sum2
    end
    avrat = sum3 / sum4
    @inbounds for g in 1:18
        radin[g] <= 5 && (ratio[g] = avrat)
    end
    @inbounds for j in 1:ntrees
        g = Int(spgrp[j])
        if radgro[j] <= 0f0
            sd1[j] = dbh[j] - ratio[g] * pdg[j]
            sd1[j] < 0.1f0 && (sd1[j] = 0.1f0)
            growth[j] = dbh[j] - sd1[j]
        end
    end
    @inbounds for i in 1:ntrees                        # PHTS (diamcal.f:217-243)
        g = Int(spgrp[i])
        prsh = _hd(g, sd1[i]); prh = _hd(g, dbh[i])
        cprsh = 4.5f0 + calib1[g] * (prsh - 4.5f0)
        cprh = 4.5f0 + calib1[g] * (prh - 4.5f0)
        sd2[i] = (cprsh / cprh) * ht[i]
        sd2[i] < 4.5f0 && (sd2[i] = 4.5f0)
    end
    @inbounds for i in 1:ntrees                        # PREVIOUS 5-YR CROWN RATIO
        sd3[i] = (cr[i] * ht[i]) / sd2[i]
        sd3[i] > 1f0 && (sd3[i] = 1f0)
    end
    @inbounds for i in 1:ntrees                        # MORTAL (diamcal.f:245-286) → expansion
        pmk = 0f0; pw = 0f0
        if expan[i] > 0f0
            sbal = oc_get_bal(dbh[i], ball2, bal2)
            g = Int(spgrp[i])
            pmk, pw = swo ? oc_pm_swo(g, dbh[i], cr[i], xsi_1, sbal, og) :
                            op_pm_nwo(g, dbh[i], cr[i], xsi_1, xsi_2, sbal)
        end
        xpm = 1f0 / (1f0 + fexp(-pmk))
        ps = fpow(1f0 - xpm, pw)
        pm = 1f0 - ps
        pm > 0.5f0 && (pm = 0.5f0)
        sd4[i] = expan[i] / (1f0 - pm)
    end
    sba1, bal1, ball1 = organon_ssum_bal(version, spgrp, sd1, sd2, sd4, ntrees, npts)
    @inbounds for g in 1:nspn
        sdg = 0f0; spdg = 0f0; yss = 0f0
        for j in 1:ntrees
            (spgrp[j] != g || radgro[j] <= 0f0) && continue
            pdg[j] = _pdg(j, sd1[j], sd3[j], sba1, ball1, bal1)
            sdg += growth[j]; spdg += pdg[j]
            yss += (growth[j] * growth[j]) / pdg[j]
        end
        radin[g] < 2 && continue
        calib3[g] = oc_caltst(sdg, spdg, yss, radin[g])
    end
    return calib3
end

"""
    organon_calib_schedule(acalib, cyclg) -> CALIB(1:3, 1:18)

organon/execute2.f:231-252 — EXECUTE rebuilds CALIB each call from ACALIB with CALH/CALC/CALD (INDS(1..3), all
.TRUE. by default, {oc,op}/grinit.f:356-358): CALIB(k,i) = (1+A)/2 + 0.5**(0.5·CYCLG)·((A−1)/2), so the PREPARE
ratio decays halfway toward 1 with a 2-cycle half-life (= A at CYCLG=0). `0.5**x` is powf.
"""
function organon_calib_schedule(acalib::AbstractMatrix{Float32}, cyclg::Int)
    cal = ones(Float32, 3, 18)
    f = fpow(0.5f0, 0.5f0 * Float32(cyclg))
    @inbounds for i in 1:18, k in 1:3
        a = acalib[k, i]
        cal[k, i] = (1f0 + a) / 2f0 + f * ((a - 1f0) / 2f0)
    end
    return cal
end

"""
    organon_dg_swo(buf; si_1, msdi_1, msdi_2, msdi_3, cyclg=0, radgro=zeros)
        -> (dgro, spgrp, sba1, bal1, ball1, a1, a2)

The ORGANON SWO diameter-growth pass for one cycle — the `GROW` "growth-1" DG sequence
(organon/grow.f:93-109) driven off the `/ORGANON/` input buffer `buf` (chunk C1). Computes species
groups, `SSTATS` stand stats, `SUBMAX` A1/A2, the `DGCALIB` CALIB(3,*) (=1.0 for RAD=.FALSE.), and
returns `dgro[i]` = ORGANON diameter growth (outside bark, 5-yr) for every buffer record — the value
`oc/dgdriv.f:447` bark-converts to `DDS` for the valid ORGANON trees (`buf.iorg[i]==1`).

`si_1` is SITE_1 − 4.5 (the ORGANON SI_1 that DG_SWO's LOG(SITE) consumes). Deterministic (DGSD=0).
The FVS bark conversion (BRATIO → DIAGR → DDS → WK2) and the StandState application are the C7
copy-back seam (FVS-native shared engine), not this chunk.
"""
function organon_dg_swo(buf::OrganonBuffer; si_1::Float32, msdi_1::Float32=0.0f0,
        msdi_2::Float32=0.0f0, msdi_3::Float32=0.0f0, cyclg::Int=0,
        calib3::Union{Nothing,AbstractVector{Float32}}=nothing, stor::AbstractVector{Float32}=zeros(Float32, 6),
        babt::Float32=0f0, bart::Vector{Float32}=zeros(Float32,5), yt::Vector{Float32}=zeros(Float32,5))
    n = buf.ntrees
    spgrp = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        spgrp[i] = oc_spgroup_swo(buf.species[i])
    end
    sba1, bal1, ball1, _, _, _, _ = oc_sstats(spgrp, buf.dbh1, buf.ht1or, buf.expan1, n)
    # execute2.f:327-331: SUBMAX only at CYCLG≤0 (or A1/A2 unset; OSTORY/INGRO/B6THIN = INDS(11..13), 0 here) — later
    # cycles reuse the A1/A2 carried in STOR(3)/STOR(4).
    if cyclg <= 0 || stor[3] == 0f0 || stor[4] == 0f0
        a1, a2 = oc_submax(spgrp, buf.dbh1, buf.expan1, n, msdi_1, msdi_2, msdi_3)
    else
        a1, a2 = stor[3], stor[4]
    end
    stor[3] = a1; stor[4] = a2
    c3 = calib3 === nothing ? ones(Float32, 18) : calib3     # execute2.f CALIB(3,·) (organon_calib_schedule)
    dgro = zeros(Float32, n)
    @inbounds for i in 1:n
        buf.expan1[i] <= 0.0f0 && continue
        g = Int(spgrp[i])
        sbal1 = oc_get_bal(buf.dbh1[i], ball1, bal1)
        dgro[i] = oc_diamgro_run(buf.species[i], g, buf.dbh1[i], buf.cr1[i], si_1,
                                 sbal1, sba1, c3[g]; cyclg=cyclg, babt=babt, bart=bart, yt=yt)
    end
    return dgro, spgrp, sba1, bal1, ball1, a1, a2
end
