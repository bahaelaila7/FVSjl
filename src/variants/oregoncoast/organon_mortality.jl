# =============================================================================
# organon_mortality.jl — OC (Oregon Coast) ORGANON SWO mortality (chunk C6).
#
# Ported from (ORGANON edition SWO / VERSION=1 path only):
#   • organon/mortality.f — MORTAL_RUN (whole-stand mortality driver: individual-tree PM +
#                           SDI-based additional mortality), PM_SWO (the per-tree logistic-mortality
#                           linear predictor), PM_FERT (fertilizer mortality adjustment), QUAD1
#                           (SDI quadratic-mean-diameter target), OLDGRO (old-growth OG; reused C5).
#   • oc/morts.f:498-504 — the FVS-side copy-back: for the ORGANON path `WK2(I)=MORTEXP(I)·(FINT/5)`
#                          (=MORTEXP at FINT=5), MORTEXP(I)=DEADEXP(I)·NPTS (execute2.f:431). The C6
#                          validation target is DEADEXP = EXPAN·PM (trees dying per acre).
#
# SCOPE: the ORGANON SWO mortality core — the deterministic (DGSD=0) per-tree DEADEXP that C4/C5 have
# been consuming as the measured MORTEXP. MORT (INDS(9)=1, oc/grinit.f:364) enables SDI-based
# additional mortality; on a below-carrying-capacity stand (RD ≤ RDCC=0.60) the base individual-tree
# mortality is used and the KR1 density adjustment does not fire, but it is ported faithfully for
# dense stands. VERSION=1 (SWO) ONLY; NWO(2)/SMC(3)/RAP(4) (OP) raise a clear error.
#
# MORTAL_RUN runs inside GROW AFTER growth-1/2 (DGRO/HGRO computed) but BEFORE the DBH/HT update, so
# it sees the ORIGINAL DBH/HT/EXPAN + the C3 dgro (GROWTH(:,2)) for the post-growth BA. It uses the
# start-of-growth BAL (BAL1/BALL1) and the SUBMAX A1/A2. OLDGRO(XIND=0) here is the PRE-growth OG.
#
# MEASURED bit-exact vs the live FVSoc_clean oracle (scoped `DEBUG 1 / DGDRIV HTGF CROWN MORTS`,
# stand S248112 / ocmin, all 27 records) — the dgdriv `MORTEXP` dump. See docs/OC_VARIANT_PORT_AUDIT.md.
# =============================================================================

# --- PM_SWO MPAR(18,9): B0..B7 + POW (species groups 1..18) ------------------------------------
const OC_PM_SWO_B0  = Float32[-4.648483270,-2.215777201,-1.050000682,-1.531051304,-1.922689902,-1.166211991,-0.761609,-4.072781265,-6.089598985,-4.317549852,-2.410756914,-2.990451960,-2.976822456,-6.00031085,-3.108619921,-2.0,-3.020345211,-1.386294361]
const OC_PM_SWO_B1  = Float32[-0.266558690,-0.162895666,-0.194363402,0.0,-0.136081990,0.0,-0.529366,-0.176433475,-0.245615070,-0.057696253,0.0,0.0,0.0,-0.10490823,-0.570366764,-0.5,0.0,0.0]
const OC_PM_SWO_B2  = Float32[0.003699110,0.003317290,0.003803100,0.0,0.002479863,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.018205398,0.015,0.0,0.0]
const OC_PM_SWO_B3  = Float32[-2.118026640,-3.561438261,-3.557300286,0.0,-3.178123293,-4.602668157,-4.74019,-1.729453975,-3.208265570,0.0,-1.049353753,0.0,-6.223250962,-0.99541909,-4.584655216,-3.0,-8.467882343,0.0]
const OC_PM_SWO_B4  = Float32[0.025499430,0.014644689,0.003971638,0.0,0.0,0.0,0.0119587,0.0,0.033348079,0.004861355,0.008845583,0.0,0.0,0.00912739,0.014926170,0.015,0.013966388,0.0]
const OC_PM_SWO_B5  = Float32[0.003361340,0.0,0.005573601,0.0,0.004684133,0.0,0.00756365,0.012525642,0.013571319,0.00998129,0.0,0.002884840,0.0,0.87115652,0.012419026,0.01,0.009461545,0.0]
const OC_PM_SWO_B6  = Float32[0.013553950,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OC_PM_SWO_B7  = Float32[-2.723470950,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OC_PM_SWO_POW = Float32[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]

"""
    oc_pm_swo(g, dbh, cr, si_1, bal, og) -> (pm, pow)

organon/mortality.f PM_SWO — the logistic-mortality linear predictor for species group `g`. `si_1`
is SITE_1−4.5 (so `si_1+4.5` = raw site), `bal` = basal-area-in-larger (start-of-growth), `og` =
old-growth indicator. Returns the linear predictor and POW (=MPAR(g,9)).
"""
@inline function oc_pm_swo(g::Int, dbh::Float32, cr::Float32, si_1::Float32, bal::Float32, og::Float32)
    b0 = OC_PM_SWO_B0[g]; b1 = OC_PM_SWO_B1[g]; b2 = OC_PM_SWO_B2[g]; b3 = OC_PM_SWO_B3[g]
    b4 = OC_PM_SWO_B4[g]; b5 = OC_PM_SWO_B5[g]; b6 = OC_PM_SWO_B6[g]; b7 = OC_PM_SWO_B7[g]
    pow = OC_PM_SWO_POW[g]
    if g == 14   # Oregon white oak
        pm = b0 + b1*dbh + b2*(dbh*dbh) + b3*cr + b4*(si_1 + 4.5f0) + b5*flog(bal + 5.0f0)
    else
        pm = b0 + b1*dbh + b2*(dbh*dbh) + b3*cr + b4*(si_1 + 4.5f0) + b5*bal + b6*bal*fexp(b7*og)
    end
    return pm, pow
end

"organon/mortality.f PM_FERT — fertilizer mortality adjustment (VERSION<=3). 0.0 with no fert."
function oc_pm_fert(g::Int, cyclg::Int, pn::Vector{Float32}, yf::Vector{Float32})
    if g == 1
        pf1 = 0.0000552859f0; pf2 = 1.5f0; pf3 = -0.5f0
    else
        pf1 = 0.0f0; pf2 = 1.0f0; pf3 = 0.0f0
    end
    xtime = Float32(cyclg)*5.0f0
    fertx1 = 0.0f0
    @inbounds for ii in 2:5
        fertx1 += pn[ii]*fexp((pf3/pf2)*(yf[1]-yf[ii]))
    end
    fertx2 = fexp(pf3*(xtime - yf[1]))
    return pf1*fpow(pn[1] + fertx1, pf2)*fertx2
end

"organon/mortality.f QUAD1 (VERSION<=3) — SDI quadratic-mean-diameter target."
@inline function oc_quad1(ni::Float32, no::Float32, rdcc::Float32, a1::Float32)
    a2 = 0.62305f0; a3 = 14.39533971f0
    a4 = -(flog(rdcc)*a2/a1)
    x = a1 - a2*flog(ni) - (a1*a4)*fexp(-a3*(flog(no) - flog(ni)))
    return fexp(x)
end

"""
    organon_mortal_swo(buf, dgro, hgro, spgrp, bal1, ball1, a1, a2; si_1, cyclg=0, mort=true)
        -> deadexp::Vector{Float32}

organon/mortality.f MORTAL_RUN (VERSION=1) — the whole-stand ORGANON mortality pass off the C1
`/ORGANON/` buffer and the C3 `dgro`/SUBMAX `a1`,`a2`/start-of-growth BAL. Returns `deadexp[i]` =
EXPAN·PM (trees dying/acre) — the value `oc/morts.f:499` copies into the FVS mortality (`MORTEXP`).
Runs before the DBH/HT update, so it sees the ORIGINAL DBH/HT/EXPAN + `dgro` for the post-growth BA.
"""
function organon_mortal_swo(buf::OrganonBuffer, dgro::Vector{Float32}, hgro::Vector{Float32},
        spgrp::Vector{Int32}, bal1::Vector{Float32}, ball1::Vector{Float32}, a1::Float32,
        a2::Float32; si_1::Float32, cyclg::Int=0, mort::Bool=true)
    n = buf.ntrees
    ib = 5
    a3 = 14.39533971f0
    rdcc = 0.60f0
    kb = 0.005454154f0
    pn = zeros(Float32, 5); yf = zeros(Float32, 5)
    dbh0 = buf.dbh1; ht0 = buf.ht1or; cr0 = buf.cr1
    expan = copy(buf.expan1)                       # TDATAR(:,4) — mutated by mortality
    deadexp = zeros(Float32, n)
    pmk = zeros(Float32, n)
    pow = ones(Float32, n)

    stba = 0.0f0; stn = 0.0f0
    @inbounds for i in 1:n
        stba += (dbh0[i]*dbh0[i])*kb*expan[i]
        stn  += expan[i]
    end
    sqmda = sqrt(stba/(kb*stn))
    rd = stn/fexp(a1/a2 - flog(sqmda)/a2)

    # OLDGRO(XIND=0) — pre-growth old-growth indicator (growth terms drop out at XIND=0)
    og1 = oc_oldgro(spgrp, dbh0, ht0, expan, dgro, hgro, deadexp, n, ib, 0.0f0)

    # individual-tree mortality linear predictors
    @inbounds for i in 1:n
        expan[i] <= 0.0f0 && continue
        g = Int(spgrp[i])
        fertadj = oc_pm_fert(g, cyclg, pn, yf)
        sbal1 = oc_get_bal(dbh0[i], ball1, bal1)
        cr = cr0[i]                                # SCR=0 ⇒ CR = TDATAR(:,3)
        pm, p = oc_pm_swo(g, dbh0[i], cr, si_1, sbal1, og1)
        pow[i] = p
        pmk[i] = pm + fertadj
    end

    # first-pass NA/BAA (post-growth BA uses DBH+DGRO)
    na = 0.0f0; baa = 0.0f0
    @inbounds for i in 1:n
        cr = cr0[i]
        cradj = cr <= 0.17f0 ? (1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))) : 1.0f0
        xpm = 1.0f0/(1.0f0 + fexp(-pmk[i]))
        ps = fpow(1.0f0 - xpm, pow[i])
        pm = 1.0f0 - ps*cradj
        na  += expan[i]*(1.0f0 - pm)
        baa += kb*((dbh0[i]+dgro[i])*(dbh0[i]+dgro[i]))*expan[i]*(1.0f0 - pm)
    end

    apply_base = true
    if mort
        qmda = sqrt(baa/(kb*na))
        ind = 0; a1max = a1; no = 0.0f0
        # CYCLG==0 initialization (organon/mortality.f:150-173). NOTE: only the CYCLG==0 branch is
        # ported here (ocmin is cyc0); the subsequent-cycle branch (:177-204, needs the carried
        # RD0/PA1MAX/NO state) is a C7-orchestration follow-up.
        rda = na/fexp(a1/a2 - flog(qmda)/a2)
        if cyclg == 0
            if rd >= 1.0f0
                a1max = rda > rd ? (flog(sqmda) + a2*flog(stn)) : (flog(qmda) + a2*flog(na))
                ind = 1
                a1max < a1 && (a1max = a1)
            else
                ind = 0
                if rd > rdcc
                    no = stn*fpow(flog(rd)/flog(rdcc), -1.0f0/a3)
                end
            end
        end
        qmdp = (ind == 0 && no > 0.0f0) ? oc_quad1(na, no, rdcc, a1) : fexp(a1max - a2*flog(na))
        if rd <= rdcc || qmdp > qmda
            apply_base = true
        else
            apply_base = false
            # SDI density adjustment — iterate KR1 up until QMDP≥QMDA (organon/mortality.f:239-298)
            kr1 = 0.0f0
            for kk in 1:7
                nk = 10.0f0/fpow(10.0f0, Float32(kk))
                while true
                    kr1 += nk
                    naa = 0.0f0; baaa = 0.0f0
                    @inbounds for i in 1:n
                        expan[i] < 0.001f0 && continue
                        cr = cr0[i]
                        cradj = cr <= 0.17f0 ? (1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))) : 1.0f0
                        xpm = 1.0f0/(1.0f0 + fexp(-(kr1 + pmk[i])))
                        ps = fpow(1.0f0 - xpm, pow[i])
                        pm = 1.0f0 - ps*cradj
                        naa  += expan[i]*(1.0f0 - pm)
                        baaa += kb*fpow(dbh0[i]+dgro[i], 2.0f0)*expan[i]*(1.0f0 - pm)
                    end
                    qmda = sqrt(baaa/(kb*naa))
                    qmdp = ind == 0 ? oc_quad1(naa, no, rdcc, a1) : fexp(a1max - a2*flog(naa))
                    if qmdp >= qmda
                        kr1 -= nk
                        break
                    end
                end
            end
            @inbounds for i in 1:n
                if expan[i] <= 0.0f0
                    deadexp[i] = 0.0f0; expan[i] = 0.0f0
                else
                    cr = cr0[i]
                    cradj = cr <= 0.17f0 ? (1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))) : 1.0f0
                    xpm = 1.0f0/(1.0f0 + fexp(-(kr1 + pmk[i])))
                    ps = fpow(1.0f0 - xpm, pow[i])
                    pm = 1.0f0 - ps*cradj
                    deadexp[i] = expan[i]*pm
                    expan[i] = expan[i]*(1.0f0 - pm)
                end
            end
        end
    end

    if apply_base
        @inbounds for i in 1:n
            cr = cr0[i]
            cradj = cr <= 0.17f0 ? (1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))) : 1.0f0
            xpm = 1.0f0/(1.0f0 + fexp(-pmk[i]))
            ps = fpow(1.0f0 - xpm, pow[i])
            pm = 1.0f0 - ps*cradj
            deadexp[i] = expan[i]*pm
            expan[i] = expan[i]*(1.0f0 - pm)
        end
    end
    return deadexp
end
