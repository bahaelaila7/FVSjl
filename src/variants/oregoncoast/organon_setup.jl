# =============================================================================
# organon_setup.jl — OC (Oregon Coast) ORGANON PREPARE setup calibration (chunk C2).
#
# Ported from (ORGANON edition SWO / VERSION=1 path only):
#   • organon/prepare.f  — PREPARE / EDIT / SPGROUP_EDIT / CKSP_EDIT (edit + SI conversion +
#                          species groups + RAD detection); orchestrates HDCALIB/PRDHT/
#                          CRCALIB/PRDCR/(DGCALIB).
#   • organon/start2.f   — HDCALIB (H-D calibration ratio → TMPCAL(1,*)), PRDHT (missing-height
#                          imputation), CRCALIB (crown-ratio calibration ratio → TMPCAL(2,*)),
#                          PRDCR (missing-crown-ratio imputation), A_HD_SWO, HD40_SWO, HS_H40,
#                          A_HCB_SWO, SPMIX, DFORTY, GET_CCFL_EDIT, CALTST.
#   • organon/diamcal.f  — SSUM (begin/end stand CCFL/BAL/SBA sums), OLDGROWTH (OG indicator),
#                          GET_BAL.
#   • organon/crngrow.f  — MCW_SWO (maximum crown width, needed by SSUM's CCF).
#   • oc/cratet.f:129-401 — the FVS-side driver: TMPCAL init, PREPARE call, and the
#                          TMPCAL→ACALIB load (cratet.f:393-401).
#
# SCOPE: this ports the ORGANON *calibration* seam — the deterministic (DGSD=0) computation of
# the ACALIB/TMPCAL height/crown/diameter calibration multipliers that ORGANON growth (C3-C6)
# consumes — plus the ORGANON HT/CR *dubbing* (PRDHT/PRDCR) for valid ORGANON trees with missing
# height/crown. VERSION=1 (SWO) ONLY. The NWO(2)/SMC(3)/RAP(4) branches (needed for OP=Olympic)
# raise a clear error; their coefficient tables sit alongside SWO in the same Fortran and are a
# thin OP follow-on.
#
# NOT in scope here (measured, see docs/OC_VARIANT_PORT_AUDIT.md C2 verdict):
#   • The DGCALIB (RAD=.TRUE.) diameter-growth calibration branch — it reuses the DG_SWO/bark/BAL
#     machinery that is chunk C3; on inventory data with no radial-increment cores RAD=.FALSE.,
#     so TMPCAL(3,*)=1.0 (organon/prepare.f:137-140) and the branch is not exercised. Deferred to
#     C3 (a loud error fires if RADGRO>0 is ever passed).
#   • The FVS-native Wykoff HT-D dubbing for NON-valid-ORGANON trees (oc/cratet.f:587-804 via
#     HTDBH) — that is shared-engine FVS machinery + OC Wykoff HT1/HT2 coefficients, NOT ORGANON.
#     (This is what actually dubs the ocmin tree-20 LP record to HT=53.32 — measured, correcting
#     the C1 "ORGANON-dubbed" label.)
#
# MEASURED bit-exact vs the live FVSoc_clean oracle (scoped `DEBUG 1 / CRATET`, stand S248112 /
# ocmin, 27 records): the CRATET `ORGANON TMPCAL(k,grp)` dump — see the C2 verdict table.
# =============================================================================

# --- SWO species-group map (organon/prepare.f SPGROUP_EDIT, SCODE1) --------------------------
# FIA code → ORGANON SWO species group (1..18). SCODE1 has 19 entries; the true-fir codes 15 and
# 17 both collapse to group 2 (ISX=J, then ISX=ISX-1 for ISX>2), giving 18 groups:
#   1=DF 2=GF/TF 3=PP 4=SP 5=IC 6=WH 7=RC 8=PY 9=MD/MA 10=GC 11=TA/TO 12=CL 13=BM 14=WO 15=BO
#   16=RA 17=PD/DG 18=WI.
const OC_SCODE1 = Int32[202,15,17,122,117,81,263,242,231,361,431,631,805,312,815,818,351,492,920]

"organon/prepare.f SPGROUP_EDIT (VERSION=1). FIA code → SWO species group 1..18, or -9999."
function oc_spgroup_swo(fia::Integer)
    @inbounds for j in 1:19
        if fia == OC_SCODE1[j]
            isx = j
            isx > 2 && (isx -= 1)
            return isx
        end
    end
    return -9999
end

# --- coefficient tables (all REAL*4 / Float32; species-group order 1..18) --------------------
# organon/start2.f A_HD_SWO HDPAR(18,3): B0,B1,B2 (height-diameter, all species).
const OC_AHD_SWO_B0 = Float32[7.153156143,6.638003799,7.181264435,6.345116767,8.776627288,6.58804,6.14817441,6.402691396,5.42457261,9.21600278,7.398142262,7.762149257,5.02002617,4.69753118,4.907340242,5.59759126,5.252315215,3.862132151]
const OC_AHD_SWO_B1 = Float32[-5.36900835,-5.44399465,-5.90709219,-5.30026188,-7.4383668,-5.35325461,-5.40092761,-4.79802411,-3.56317104,-7.63409138,-5.5099273,-6.04759773,-2.51228202,-3.51586969,-3.18017969,-3.19942952,-3.13509983,-1.5294776]
const OC_AHD_SWO_B2 = Float32[-0.25832512,-0.33929196,-0.27533719,-0.35264183,-0.16906224,-0.31897786,-0.38922036,-0.16317997,-0.36177689,-0.15346440,-0.19080702,-0.16308399,-0.42256497,-0.57665068,-0.46654227,-0.38783403,-0.26979750,-0.62476287]

# organon/start2.f HD40_SWO HD40PAR(3,3): DF,GW,PP (IEQ 1,2,3).
const OC_HD40_SWO_B0 = Float32[-3.485635287,-4.376160718,-4.047994965]
const OC_HD40_SWO_B1 = Float32[-0.255712209,-0.231693907,-0.135864020]
const OC_HD40_SWO_B2 = Float32[-0.001555149,-0.001334070,-0.005647510]

# organon/start2.f A_HCB_SWO HCBPAR(18,7): B0..B6 (height-to-crown-base, all species).
const OC_HCB_SWO_B0 = Float32[1.990155033,4.800089990,2.024723585,3.582314301,3.127730861,0.0,4.49102006,0.0,3.271130882,0.387912505,0.4488479442,1.285465907,1.000364090,1.05786632,2.672850866,0.56713781,0.0,0.0]
const OC_HCB_SWO_B1 = Float32[-0.008180786,0.0,-0.001953589,-0.003256792,-0.004386780,0.0,0.0,0.0,0.0,-0.015000868,-0.009375810,-0.024459278,-0.010636441,0.0,0.0,-0.010377976,0.0,0.0]
const OC_HCB_SWO_B2 = Float32[-0.004696095,-0.003268539,-0.001837480,0.0,-0.003557122,0.0,-0.00132412,0.0,0.0,-0.004098099,-0.001822050,-0.003992574,-0.005950398,-0.00183283,-0.001400851,-0.002066036,-0.004842962,-0.004842962]
const OC_HCB_SWO_B3 = Float32[-0.392033240,-0.858744969,-0.568909853,-0.765250973,-0.637929879,0.0,-1.01460531,0.0,-0.841331291,0.0,0.0,0.0,0.0,-0.28644547,-0.605971926,0.0,-0.567987126,-0.567987126]
const OC_HCB_SWO_B4 = Float32[1.945708371,0.0,4.831886553,3.043845568,0.977816058,3.246352823,0.0,1.225564582,1.791699815,2.104871164,0.0,0.0,0.0,0.0,0.0,1.39796223,0.0,0.0]
const OC_HCB_SWO_B5 = Float32[0.007854260,0.0,0.001653030,0.0,0.005850321,0.0,0.01340624,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0281315332,0.0281315332]
const OC_HCB_SWO_B6 = Float32[0.295593583,0.275679490,0.0,0.0,0.257070387,0.0,0.0,0.0,0.927163029,0.352773356,0.233233237,0.0,0.310672769,0.0,0.430988703,0.0,0.0,0.0]

# organon/crngrow.f MCW_SWO MCWPAR(18,4): B0,B1,B2,PKDBH (maximum crown width).
const OC_MCW_SWO_B0 = Float32[4.6366,6.1880,3.4835,4.6600546,3.2837,4.5652,4.0,4.5652,3.4298629,2.9793895,4.4443,4.4443,4.0953,3.0785639,3.3625,8.0,2.9793895,2.9793895]
const OC_MCW_SWO_B1 = Float32[1.6078,1.0069,1.343,1.0701859,1.2031,1.4147,1.65,1.4147,1.3532302,1.5512443,1.7040,1.7040,2.3849,1.9242211,2.0303,1.53,1.5512443,1.5512443]
const OC_MCW_SWO_B2 = Float32[-0.009625,0.0,-0.0082544,0.0,-0.0071858,0.0,0.0,0.0,0.0,-0.01416129,0.0,0.0,-0.011630,0.0,-0.0073307,0.0,-0.01416129,-0.01416129]
const OC_MCW_SWO_PK = Float32[88.52,999.99,81.35,999.99,83.71,999.99,999.99,999.99,999.99,54.77,999.99,999.99,102.53,999.99,138.93,999.99,54.77,54.77]

# organon/start2.f CALTST t-critical table (df 1..30 direct, then interpolated).
const OC_TVAL = Float32[63.657,9.925,5.841,4.604,4.032,3.707,3.499,3.355,3.250,3.169,3.106,3.055,3.012,2.977,2.947,2.921,2.898,2.878,2.861,2.845,2.831,2.819,2.807,2.797,2.787,2.779,2.771,2.763,2.756,2.750,2.704,2.660,2.617,2.576]

"organon/crngrow.f MCW_SWO — maximum crown width for group `g`, DBH `d`, height `h`."
@inline function oc_mcw_swo(g::Int, d::Float32, h::Float32)
    b0 = OC_MCW_SWO_B0[g]; b1 = OC_MCW_SWO_B1[g]; b2 = OC_MCW_SWO_B2[g]; pk = OC_MCW_SWO_PK[g]
    dbh = d > pk ? pk : d
    if h < 4.501f0
        return h / 4.5f0 * b0
    else
        return b0 + b1*dbh + b2*dbh*dbh
    end
end

"organon/start2.f HS_H40 — Hann-Scrivani dominant top height at age (ISP 1=DF,2=PP)."
@inline function oc_hs_h40(isp::Int, age::Float32, si::Float32)
    s = si - 4.5f0
    if isp == 1
        x1 = 1f0 - fexp(-fexp(-6.21693f0 + 0.281176f0*flog(s) + 1.14354f0*flog(age)))
        x2 = 1f0 - fexp(-fexp(-6.21693f0 + 0.281176f0*flog(s) + 1.14354f0*flog(50.0f0)))
    else
        x1 = 1f0 - fexp(-fexp(-6.54707f0 + 0.288169f0*flog(s) + 1.21297f0*flog(age)))
        x2 = 1f0 - fexp(-fexp(-6.54707f0 + 0.288169f0*flog(s) + 1.21297f0*flog(50.0f0)))
    end
    return 4.5f0 + s*(x1/x2)
end

"organon/start2.f HD40_SWO — predicted height via 40-largest H40/D40 (IEQ 1=DF,2=GW,3=PP)."
@inline function oc_hd40_swo(ieq::Int, ht40::Float32, d40::Float32, dbh::Float32)
    b0 = OC_HD40_SWO_B0[ieq]; b1 = OC_HD40_SWO_B1[ieq]; b2 = OC_HD40_SWO_B2[ieq]
    exd   = fexp(b0*fpow(dbh, b1 + b2*(ht40-4.5f0)))
    exd40 = fexp(b0*fpow(d40, b1 + b2*(ht40-4.5f0)))
    return 4.5f0 + (ht40-4.5f0)*(exd/exd40)
end

"organon/start2.f A_HD_SWO — predicted height from DBH (group `g`)."
@inline function oc_a_hd_swo(g::Int, dbh::Float32)
    b0 = OC_AHD_SWO_B0[g]; b1 = OC_AHD_SWO_B1[g]; b2 = OC_AHD_SWO_B2[g]
    return 4.5f0 + fexp(b0 + b1*fpow(dbh, b2))
end

"organon/start2.f A_HCB_SWO — predicted height to crown base (group `g`)."
@inline function oc_a_hcb_swo(g::Int, ht::Float32, dbh::Float32, ccfl::Float32, ba::Float32,
                              si1::Float32, si2::Float32, og::Float32)
    b0 = OC_HCB_SWO_B0[g]; b1 = OC_HCB_SWO_B1[g]; b2 = OC_HCB_SWO_B2[g]; b3 = OC_HCB_SWO_B3[g]
    b4 = OC_HCB_SWO_B4[g]; b5 = OC_HCB_SWO_B5[g]; b6 = OC_HCB_SWO_B6[g]
    si = g == 3 ? si2 : si1
    return ht/(1f0 + fexp(b0 + b1*ht + b2*ccfl + b3*flog(ba) + b4*(dbh/ht) + b5*si + b6*og*og))
end

"organon/start2.f GET_CCFL_EDIT — crown-competition-factor-in-larger-trees at DBH (1-based)."
@inline function oc_get_ccfl(dbh::Float32, ccfll::Vector{Float32}, ccfl::Vector{Float32})
    if dbh > 100.0f0
        return 0.0f0
    elseif dbh > 50.0f0
        k = trunc(Int, dbh - 49.0f0)
        return ccfll[k]
    else
        k = trunc(Int, dbh*10.0f0 + 0.5f0)
        return ccfl[k]
    end
end

"organon/start2.f CALTST — calibration β with t-test toward 1.0 and [0.5,2.0] clamp."
function oc_caltst(yxs::Float32, xss::Float32, yss::Float32, n::Int)
    beta = yxs/xss
    idegf = n - 1
    degf = Float32(idegf)
    mse = (yss - 2.0f0*beta*yxs + beta*beta*xss)/degf
    varbeta = mse/xss
    if varbeta > 0.0f0
        ttest = abs(beta - 1.0f0)/sqrt(varbeta)
        if idegf <= 30
            critval = OC_TVAL[idegf]
        elseif idegf <= 40
            critval = OC_TVAL[30] + (OC_TVAL[31]-OC_TVAL[30])*((degf-30.0f0)/10.0f0)
        elseif idegf <= 60
            critval = OC_TVAL[31] + (OC_TVAL[32]-OC_TVAL[31])*((degf-40.0f0)/20.0f0)
        elseif idegf <= 120
            critval = OC_TVAL[32] + (OC_TVAL[33]-OC_TVAL[32])*((degf-60.0f0)/60.0f0)
        elseif idegf <= 1000
            critval = OC_TVAL[33] + (OC_TVAL[34]-OC_TVAL[33])*((degf-120.0f0)/880.0f0)
        else
            critval = OC_TVAL[34]
        end
        if ttest < critval
            beta = 1.0f0
        elseif beta > 2.0f0
            beta = 2.0f0
        elseif beta < 0.5f0
            beta = 0.5f0
        end
    else
        if beta > 2.0f0
            beta = 2.0f0
        elseif beta < 0.5f0
            beta = 0.5f0
        end
    end
    return beta
end

"organon/start2.f SPMIX — basal-area species proportions (PDF,PTF,PPP,PWH,PRA)."
function oc_spmix(species::Vector{Int32}, dbh::Vector{Float32}, expan::Vector{Float32},
                  ntrees::Int, npts::Int)
    batot = badf = batf = bapp = bawh = bara = 0.0f0
    @inbounds for i in 1:ntrees
        ba = (dbh[i]*dbh[i]*expan[i])*0.005454154f0/Float32(npts)
        batot += ba
        s = species[i]
        if s == 15 || s == 17
            batf += ba
        elseif s == 122
            bapp += ba
        elseif s == 202
            badf += ba
        elseif s == 263
            bawh += ba
        elseif s == 351
            bara += ba
        end
    end
    if batot > 0.0f0
        return (badf/batot, batf/batot, bapp/batot, bawh/batot, bara/batot)
    end
    return (0.0f0,0.0f0,0.0f0,0.0f0,0.0f0)
end

"organon/start2.f DFORTY — quadratic mean DBH of the 40 largest big-6 trees/acre (VERSION=1)."
function oc_dforty(spgrp::Vector{Int32}, dbh::Vector{Float32}, expan::Vector{Float32},
                   ntrees::Int, npts::Int, ib::Int)
    dcl = zeros(Float32, 1000); trcl = zeros(Float32, 1000)
    @inbounds for i in 1:ntrees
        dbh[i] <= 0.0f0 && continue
        if spgrp[i] <= ib
            id = trunc(Int, dbh[i]*10.0f0)
            id > 1000 && (id = 1000)
            ex = expan[i]/Float32(npts)
            dcl[id]  += dbh[i]*ex
            trcl[id] += ex
        end
    end
    totd = 0.0f0; tottr = 0.0f0
    @inbounds for i in 1000:-1:1
        totd  += dcl[i]
        tottr += trcl[i]
        if tottr > 40.0f0
            trdiff = trcl[i] - (tottr - 40.0f0)
            totd = totd - dcl[i] + ((dcl[i]/trcl[i])*trdiff)
            tottr = 40.0f0
            break
        end
    end
    return tottr > 0.0f0 ? totd/tottr : 0.0f0
end

"organon/diamcal.f OLDGROWTH — old-growth indicator OG from the 5 largest big-6 trees/acre."
function oc_oldgrowth(spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32},
                      expan::Vector{Float32}, ntrees::Int, npts::Int, ib::Int)
    htcl = zeros(Float32, 100); dcl = zeros(Float32, 100); trcl = zeros(Float32, 100)
    @inbounds for i in 1:ntrees
        if spgrp[i] <= ib
            ex = expan[i]/Float32(npts)
            id = trunc(Int, dbh[i]) + 1
            id > 100 && (id = 100)
            htcl[id] += ht[i]*ex
            dcl[id]  += dbh[i]*ex
            trcl[id] += ex
        end
    end
    totht = 0.0f0; totd = 0.0f0; tottr = 0.0f0
    @inbounds for i in 100:-1:1
        totht += htcl[i]; totd += dcl[i]; tottr += trcl[i]
        if tottr > 5.0f0
            trdiff = trcl[i] - (tottr - 5.0f0)
            totht = totht - htcl[i] + ((htcl[i]/trcl[i])*trdiff)
            totd  = totd  - dcl[i]  + ((dcl[i]/trcl[i])*trdiff)
            tottr = 5.0f0
            break
        end
    end
    if tottr > 0.0f0
        ht5 = totht/tottr; dbh5 = totd/tottr
        return dbh5*ht5/10000.0f0
    end
    return 0.0f0
end

"""
organon/diamcal.f SSUM (N=2 / ending stats) — returns (sba, ccfl, ccfll). BAL/BALL are computed
identically but unused by the crown calibration, so only the CCFL sums are returned.
"""
function oc_ssum(spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32},
                 expan::Vector{Float32}, ntrees::Int, npts::Int)
    ccfl = zeros(Float32, 500); ccfll = zeros(Float32, 51)
    sba = 0.0f0
    @inbounds for i in 1:ntrees
        g = Int(spgrp[i])
        ex = expan[i]/Float32(npts)
        ex <= 0.0f0 && continue
        ba = (dbh[i]*dbh[i]*ex)*0.005454154f0
        sba += ba
        mcw = oc_mcw_swo(g, dbh[i], ht[i])
        ccf = 0.001803f0*mcw*mcw*expan[i]/Float32(npts)
        if dbh[i] > 50.0f0
            l = trunc(Int, dbh[i] - 49.0f0)
            l > 52 && (l = 52)
            for k in 1:500
                ccfl[k] += ccf
            end
            for k in 1:(l-1)
                ccfll[k] += ccf
            end
        else
            l = trunc(Int, dbh[i]*10.0f0 + 0.5f0)
            for k in 1:(l-1)
                ccfl[k] += ccf
            end
        end
    end
    return sba, ccfl, ccfll
end

"""
    OrganonCalib

Result of the ORGANON PREPARE setup for one cycle:
  • `acalib` — the 3×18 ACALIB calibration multipliers (row 1=HT, 2=CR, 3=DG) after the
    `oc/cratet.f:393-401` TMPCAL→ACALIB load (keyword-set ACALIB entries are preserved; PREPARE
    values fill only the untouched 1.0 slots that PREPARE moved off 1.0).
  • `tmpcal` — the raw PREPARE calibration ratios (before the ACALIB load), as dumped by CRATET.
  • `dbh`,`ht`,`cr`,`expan` — the (possibly HT/CR-dubbed) tree arrays PREPARE returns.
  • `ierror` — 0 on success (1 if EDIT flagged a stand/tree error and PREPARE bailed).
"""
struct OrganonCalib
    acalib ::Matrix{Float32}
    tmpcal ::Matrix{Float32}
    dbh    ::Vector{Float32}
    ht     ::Vector{Float32}
    cr     ::Vector{Float32}
    expan  ::Vector{Float32}
    ierror ::Int
end

"""
    organon_prepare_swo(species, dbh, ht, cr, expan, radgro, ntrees, npts, stage, bhage,
                        si_1, si_2, msdi_1, msdi_2, msdi_3, pden, ieven; acalib0=nothing)

Port of `organon/prepare.f` PREPARE for ORGANON edition SWO (VERSION=1), plus the `oc/cratet.f`
TMPCAL init + TMPCAL→ACALIB load. Returns an `OrganonCalib`.

Inputs are the `/ORGANON/` buffer arrays (`SPECIES`,`DBH1`,`HT1OR`,`CR1`,`EXPAN1`,`RADGRO`) for
records 1..ntrees, exactly as `oc/cratet.f:228-252` fills them. `acalib0` is the pre-PREPARE
ACALIB (3×18) with any keyword-set multipliers (default all-1.0). Deterministic (DGSD=0).
"""
function organon_prepare_swo(species::Vector{Int32}, dbh0::Vector{Float32}, ht0::Vector{Float32},
        cr0::Vector{Float32}, expan0::Vector{Float32}, radgro::Vector{Float32}, ntrees::Int,
        npts::Int, stage::Int, bhage::Int, si_1::Float32, si_2::Float32, msdi_1::Float32,
        msdi_2::Float32, msdi_3::Float32, pden::Float32, ieven::Int;
        acalib0::Union{Nothing,Matrix{Float32}}=nothing)

    even = ieven == 1
    ib = 5; nspn = 18                              # SWO (organon/prepare.f EDIT)

    # working copies (PREPARE mutates DBH/HT/CR/EXPAN)
    dbh = copy(dbh0); ht = copy(ht0); cr = copy(cr0); expan = copy(expan0)

    # --- EDIT: species groups, missing-value flags, RAD detection, SI conversion --------------
    spgrp = Vector{Int32}(undef, ntrees)
    missht = false; misscr = false; rad = false
    @inbounds for i in 1:ntrees
        spgrp[i] = oc_spgroup_swo(species[i])
        ht[i]  <= 0.0f0 && (missht = true)
        cr[i]  <= 0.0f0 && (misscr = true)
        radgro[i] > 0.0f0 && (rad = true)
    end
    # SI conversion (organon/prepare.f:356-361, VERSION=1)
    if si_1 <= 0.0f0 && si_2 > 0.0f0
        si_1 = 1.062934f0*si_2
    elseif si_2 <= 0.0f0
        si_2 = 0.940792f0*si_1
    end

    tmpcal = ones(Float32, 3, 18)                  # oc/cratet.f:129-134 (ORGANON's ACALIB output)

    if rad
        error("OC ORGANON PREPARE: RAD=.TRUE. (radial-increment) DGCALIB branch is chunk C3 " *
              "(reuses DG_SWO/bark/BAL). No radial cores in FVS/FIA inventory ⇒ RAD=.FALSE. and " *
              "TMPCAL(3,*)=1.0. See docs/OC_ORGANON_PORT_PLAN.md.")
    end

    # --- HDCALIB (organon/start2.f) — height-diameter calibration ratio → tmpcal(1,*) ---------
    entht  = zeros(Int, 18); entdbh = zeros(Int, 18)
    yxs = zeros(Float32, 18); xss = zeros(Float32, 18); yss = zeros(Float32, 18)
    ptrht = zeros(Float32, ntrees)
    d40 = 0.0f0; ht40 = 0.0f0
    pdf = ptf = ppp = pwh = pra = 0.0f0
    if even
        age1 = Float32(bhage)
        (pdf, ptf, ppp, pwh, pra) = oc_spmix(species, dbh, expan, ntrees, npts)
        d40 = oc_dforty(spgrp, dbh, expan, ntrees, npts, ib)
        if pdf >= 0.80f0
            ht40 = oc_hs_h40(1, age1, si_1)
        elseif ptf >= 0.80f0
            ht40 = oc_hs_h40(1, age1, si_1)
        elseif ppp >= 0.80f0
            ht40 = oc_hs_h40(2, age1, si_2)
        end
    end
    @inbounds for i in 1:ntrees
        g = Int(spgrp[i])
        entdbh[g] += 1
        ht[i] <= 0.0f0 && continue
        if pdf >= 0.80f0 && g == 1
            ptrht[i] = oc_hd40_swo(1, ht40, d40, dbh[i])
        elseif ptf >= 0.80f0 && g == 2
            ptrht[i] = oc_hd40_swo(2, ht40, d40, dbh[i])
        elseif ppp >= 0.80f0 && g == 3
            ptrht[i] = oc_hd40_swo(3, ht40, d40, dbh[i])
        else
            ptrht[i] = oc_a_hd_swo(g, dbh[i])
        end
        y = ht[i] - 4.5f0
        x = ptrht[i] - 4.5f0
        wt = dbh[i]
        yxs[g] += y*x/wt
        xss[g] += x*x/wt
        yss[g] += y*y/wt
        entht[g] += 1
    end
    @inbounds for g in 1:nspn
        entdbh[g] == 0 && continue
        entht[g]  <  2 && continue
        tmpcal[1, g] = oc_caltst(yxs[g], xss[g], yss[g], entht[g])
    end

    # --- PRDHT (organon/start2.f) — impute missing heights for valid ORGANON trees -----------
    if missht
        @inbounds for i in 1:ntrees
            g = Int(spgrp[i])
            dbh[i] <= 0.0f0 && continue
            ht[i]  >  0.0f0 && continue
            if pdf >= 0.80f0 && g == 1
                ptrht[i] = oc_hd40_swo(1, ht40, d40, dbh[i])
            elseif ptf >= 0.80f0 && g == 2
                ptrht[i] = oc_hd40_swo(2, ht40, d40, dbh[i])
            elseif ppp >= 0.80f0 && g == 3
                ptrht[i] = oc_hd40_swo(3, ht40, d40, dbh[i])
            else
                ptrht[i] = oc_a_hd_swo(g, dbh[i])
            end
            newht = 4.5f0 + tmpcal[1, g]*(ptrht[i] - 4.5f0)
            newht < 4.6f0 && (newht = 4.6f0)
            ht[i] = newht
        end
    end

    # --- CRCALIB (organon/start2.f) — crown-ratio calibration ratio → tmpcal(2,*) -------------
    sba, ccfl, ccfll = oc_ssum(spgrp, dbh, ht, expan, ntrees, npts)
    og = oc_oldgrowth(spgrp, dbh, ht, expan, ntrees, npts, ib)
    entcr = zeros(Int, 18)
    cyxs = zeros(Float32, 18); cxss = zeros(Float32, 18); cyss = zeros(Float32, 18)
    pcr = zeros(Float32, ntrees)
    xsi_1 = si_1 - 4.5f0; xsi_2 = si_2 - 4.5f0
    @inbounds for i in 1:ntrees
        cr[i] <= 0.0f0 && continue
        g = Int(spgrp[i])
        sccfl = oc_get_ccfl(dbh[i], ccfll, ccfl)
        phcb = oc_a_hcb_swo(g, ht[i], dbh[i], sccfl, sba, xsi_1, xsi_2, og)
        pcr[i] = 1.0f0 - phcb/ht[i]
        x = pcr[i]; y = cr[i]
        cyxs[g] += x*y
        cxss[g] += x*x
        cyss[g] += y*y
        entcr[g] += 1
    end
    @inbounds for g in 1:nspn
        entdbh[g] == 0 && continue
        entcr[g]  <  2 && continue
        tmpcal[2, g] = oc_caltst(cyxs[g], cxss[g], cyss[g], entcr[g])
    end

    # --- PRDCR (organon/start2.f) — impute missing crown ratios for valid ORGANON trees -------
    if misscr
        @inbounds for i in 1:ntrees
            cr[i] > 0.0f0 && continue
            g = Int(spgrp[i])
            dbh[i] <= 0.0f0 && continue
            sccfl = oc_get_ccfl(dbh[i], ccfll, ccfl)
            phcb = oc_a_hcb_swo(g, ht[i], dbh[i], sccfl, sba, xsi_1, xsi_2, og)
            pcr[i] = 1.0f0 - phcb/ht[i]
            newcr = pcr[i]*tmpcal[2, g]
            newcr > 1.0f0  && (newcr = 1.0f0)
            newcr < 0.05f0 && (newcr = 0.05f0)
            cr[i] = newcr
        end
    end

    # RAD=.FALSE. ⇒ tmpcal(3,*)=1.0 (organon/prepare.f:137-140) — already initialized.

    # --- oc/cratet.f:393-401 TMPCAL → ACALIB load --------------------------------------------
    acalib = acalib0 === nothing ? ones(Float32, 3, 18) : copy(acalib0)
    @inbounds for i in 1:18, j in 1:3
        if acalib[j, i] == 1.0f0 && (tmpcal[j, i] > 0.0f0 && tmpcal[j, i] != 1.0f0)
            acalib[j, i] = tmpcal[j, i]
        end
    end

    return OrganonCalib(acalib, tmpcal, dbh, ht, cr, expan, 0)
end

"""
    oc_organon_prepare!(s) -> nothing

Run the OC ORGANON setup (`oc/cratet.f:155-401` — the CRATET ORGANON section) on the initialized
`StandState`, ONCE, before the FVS-native missing-value dubbing. It:

 1. Flags valid ORGANON trees with the CRATET (setup) eligibility gate — `DBH >= 0.1 AND
    (HT == 0 OR HT > 4.5)` AND species ∈ the 18 ORGANON species (`oc/cratet.f:171-201`). NOTE this
    gate INCLUDES blank-height (`HT == 0`) valid records, unlike the grow-time `build_organon_buffer!`
    gate (`HT > 4.5`) — that is the whole point: a blank-height ORGANON tree must be dubbed HERE.
 2. If the stand has ≥1 big-6 tree, marshals ALL live records into the `/ORGANON/` PREPARE buffer
    (`oc/cratet.f:227-252`: `DBH1=max(DBH,0.1)`, `HT1OR=HT` (floored to 4.6 only when HT>0, so a
    missing height stays 0), `CR1=ICR/100`, `EXPAN1=PROB·PI` where `PI=IPTINV`, `SPECIES=ORGSPC`),
    calls `organon_prepare_swo`, and writes the ORGANON-dubbed HT/CR back into the valid ORGANON
    tree records that were MISSING them (`oc/cratet.f:349-365`). The dubbed HT then survives the
    subsequent FVS `dub_missing_heights!` (HT > 0 ⇒ skipped there).
 3. Stores the resulting `ACALIB(3,18)` on `s.calib.organon_acalib` for the growth path (HTGRO2
    consumes row 1; SWO grow-time crown/DG ignore rows 2/3 — only NWO/SMC crown and the RAD DG path
    use them, neither active on OC/FIA inventory).

No-op unless the stand has a big-6 tree (`oc/cratet.f:215-219` `GO TO 261`) — then ACALIB stays
all-1.0 and every tree is dubbed FVS-native. OC is hardcoded even-aged (`oc/grinit.f:359` INDS(4)=1)
and SWO (VERSION=1). Deterministic (DGSD=0).
"""
function oc_organon_prepare!(s::StandState)
    t = s.trees
    n = t.n
    n == 0 && return nothing
    # --- CRATET setup-time eligibility gate (oc/cratet.f:171-201) ---
    iorg   = zeros(Int32, n)
    nbig6  = 0
    @inbounds for i in 1:n
        sp = Int(t.species[i])
        h  = t.height[i]
        ihflag = (h == 0.0f0) || (h > 4.5f0)           # measured-HT lower limit OR missing (cratet.f:178)
        if t.dbh[i] >= 0.1f0 && ihflag
            (sp in OC_ORGANON_BIG6) && (nbig6 += 1)
            iorg[i] = (sp in OC_ORGANON_VALID) ? Int32(1) : Int32(0)
        end
    end
    # --- no big-6 ⇒ ORGANON does not run; FVS-native dubbing/calibration for all (cratet.f:215-219) ---
    nbig6 == 0 && return nothing
    # --- marshal ALL live records into the PREPARE buffer (cratet.f:227-252) ---
    species = Vector{Int32}(undef, n)
    dbh1    = Vector{Float32}(undef, n)
    ht1or   = Vector{Float32}(undef, n)
    cr1     = Vector{Float32}(undef, n)
    expan1  = Vector{Float32}(undef, n)
    pival   = s.plot.pi > 0f0 ? s.plot.pi : 1f0        # PI = FLOAT(IPTINV) (initre.f:336; standstats.jl)
    @inbounds for i in 1:n
        species[i] = orgspc(Int(t.species[i]))
        d = t.dbh[i]; d < 0.1f0 && (d = 0.1f0)
        dbh1[i]  = d
        h = t.height[i]; (h > 0f0 && h < 4.6f0) && (h = 4.6f0)   # floor to 4.6 ONLY when HT>0 (cratet.f:234)
        ht1or[i] = h
        cr1[i]   = Float32(t.crown_pct[i]) / 100f0
        expan1[i] = t.tpa[i] * pival                              # EXPAN1 = PROB·PI (cratet.f:236)
    end
    radgro = zeros(Float32, n)
    stage  = Int(s.plot.stand_age)                     # STAGE = IAGE + IY(ICYC)-IY(1); ICYC=1 ⇒ = IAGE (cratet.f:259)
    bhage  = stage - 6                                  # BREAST HEIGHT AGE (cratet.f:260)
    si_1 = s.plot.sp_site_index[7]                      # DF site (SI_1); PREPARE's own conversion handles ≤0
    si_2 = s.plot.sp_site_index[18]                     # PP site (SI_2)
    npts = max(1, Int(round(pival)))
    res = organon_prepare_swo(species, dbh1, ht1or, cr1, expan1, radgro, n, npts, stage, bhage,
                              si_1, si_2, 0f0, 0f0, 0f0, 0f0, 1)   # IEVEN=1 (OC hardcoded even-aged)
    # --- write ORGANON-dubbed HT/CR back into the valid-ORGANON records that were MISSING them
    #     (oc/cratet.f:349-365: only IORG=1 trees are reloaded; KNTOHT/KNTOCR count the imputed ones) ---
    @inbounds for i in 1:n
        iorg[i] == 1 || continue
        if t.height[i] <= 0.0f0 && res.ht[i] > 0.0f0
            t.height[i] = res.ht[i]
        end
        if t.crown_pct[i] <= 0 && res.cr[i] > 0.0f0
            t.crown_pct[i] = round(Int32, res.cr[i] * 100.0f0, RoundNearestTiesAway)  # NINT (cratet.f:359)
        end
    end
    # --- store ACALIB for the growth path (cratet.f:393-401 already folded into organon_prepare_swo) ---
    copyto!(s.calib.organon_acalib, res.acalib)
    return nothing
end
