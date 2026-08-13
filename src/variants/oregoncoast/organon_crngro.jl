# =============================================================================
# organon_crngro.jl — OC (Oregon Coast) ORGANON SWO crown growth (chunk C5).
#
# Ported from (ORGANON edition SWO / VERSION=1 path only):
#   • organon/crngrow.f  — CROWGRO (5-yr crown recession → new crown ratio), HCB_SWO (height to
#                          crown base), MAXHCB_SWO (max height to crown base). Reuses the C2/C4
#                          crown-closure primitives (MCW_SWO, GET_CCFL, SSTATS BAL/CCFL).
#   • organon/mortality.f — OLDGRO (old-growth indicator OG, big-6 5-largest; start/end variants).
#   • oc/crown.f:281-288 — the FVS-side copy-back: for IORG=1 trees `CRNEW(I)=ANINT(CR2(I)*100)`.
#                          The C5 validation target is CR2 (the ORGANON crown ratio after growth).
#
# SCOPE: the ORGANON SWO crown-recession core — the deterministic (DGSD=0) per-tree CR2 that FVS
# rounds into CRNEW for the valid ORGANON trees, bypassing native OC crown. CROWGRO runs AFTER DG
# (C3) and HG (C4) advance DBH/HT: it compares the height-to-crown-base at start-of-growth (old
# DBH/HT) and end-of-growth (new DBH/HT) and lets the crown base rise by the difference. VERSION=1
# (SWO) ONLY; NWO(2)/SMC(3)/RAP(4) (OP) raise a clear error.
#
# C6 DEPENDENCY (measured): OLDGRO's start-of-growth OG1 augments each big-6 tree's expansion by its
# DEADEXP (ORGANON mortality, chunk C6) — `EXPAN = TDATAR(4) − XIND·DEADEXP` (mortality.f:763). So
# `organon_cr_swo` takes a `deadexp` vector; until C6 lands it is fed the measured MORTEXP (as C4
# consumed the measured-then-ported DGRO). OG2 (end) uses DEADEXP=0-weighted (XIND=0) so needs none.
#
# MEASURED bit-exact vs the live FVSoc_clean oracle (scoped `DEBUG 1 / DGDRIV HTGF CROWN`, stand
# S248112 / ocmin, 17 valid ORGANON trees) — the `ORG CROWN … CR2 …` dump. See docs/OC_VARIANT_PORT_AUDIT.md.
# =============================================================================

# --- HCB_SWO HCBPAR(18,7): B0..B6 (height to crown base — DISTINCT from C2's A_HCB_SWO) ---------
const OC_HCB2_SWO_B0 = Float32[1.797136911,3.451045887,1.656364063,3.785155749,2.428285297,0.0,4.49102006,0.0,2.955339267,0.544237656,0.833006499,0.5376600543,0.9411395642,1.05786632,2.60140655,0.56713781,0.0,0.0]
const OC_HCB2_SWO_B1 = Float32[-0.010188791,-0.005985239,-0.002755463,-0.009012547,-0.006882851,0.0,0.0,0.0,0.0,-0.020571754,-0.012984204,-0.018632397,-0.00768402,0.0,0.0,-0.010377976,0.0,0.0]
const OC_HCB2_SWO_B2 = Float32[-0.003346230,-0.003211194,0.0,-0.003318574,-0.002612590,0.0,-0.00132412,0.0,0.0,-0.004317523,-0.002704717,0.0,-0.005476131,-0.00183283,-0.002273616,-0.002066036,-0.005666559,-0.005666559]
const OC_HCB2_SWO_B3 = Float32[-0.412217810,-0.671479750,-0.568302547,-0.670270058,-0.572782216,0.0,-1.01460531,0.0,-0.798610738,0.0,0.0,0.0,0.0,-0.28644547,-0.554980629,0.0,-0.745540494,-0.745540494]
const OC_HCB2_SWO_B4 = Float32[3.958656001,3.931095518,6.730693919,2.758645081,2.113378338,4.801329946,0.0,2.030940382,3.095269471,3.132713612,0.0,0.0,0.0,0.0,0.0,1.39796223,0.0,0.0]
const OC_HCB2_SWO_B5 = Float32[0.008526562,0.003115567,0.001852526,0.0,0.008480754,0.0,0.01340624,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.038476613,0.038476613]
const OC_HCB2_SWO_B6 = Float32[0.448909636,0.516180892,0.0,0.841525071,0.506226895,0.0,0.0,0.0,0.700465646,0.483748898,0.2491242765,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

# --- MAXHCB_SWO MAXPAR(18,5): B0,B1,B2,B3,LIMIT -----------------------------------------------
const OC_MAXHCB_SWO_B0  = Float32[0.96,0.96,1.01,1.02,0.97,1.01,0.96,0.85,0.981,1.0,0.98,1.0,1.0,1.0,1.0,0.93,1.0,0.985]
const OC_MAXHCB_SWO_B1  = Float32[0.26,0.31,0.36,0.27,0.22,0.36,0.31,0.35,0.161,0.45,0.33,0.45,0.45,0.3,0.2,0.18,0.45,0.285]
const OC_MAXHCB_SWO_B2  = Float32[-0.987864873,-2.450718394,-1.041915784,-0.922718593,-0.002612590,-0.944528054,-1.059636222,-0.922868139,-1.73666044,-1.219919284,-0.911341687,-0.922025464,-1.020016685,-0.95634399,-1.053892465,-0.928243505,-1.020016685,-0.969750805]
const OC_MAXHCB_SWO_B3  = Float32[1.0,1.0,0.6,0.4,1.0,0.6,1.0,0.8,1.0,1.2,1.0,1.0,1.0,1.1,1.0,1.0,1.0,0.9]
const OC_MAXHCB_SWO_LIM = Float32[0.95,0.95,0.95,0.96,0.95,0.96,0.95,0.80,0.98,0.98,0.97,0.98,0.95,0.98,0.98,0.92,0.95,0.98]

"organon/crngrow.f HCB_SWO — height to crown base (species group `g`)."
@inline function oc_hcb_swo(g::Int, ht::Float32, dbh::Float32, ccfl::Float32, ba::Float32,
                            si_1::Float32, si_2::Float32, og::Float32)
    b0 = OC_HCB2_SWO_B0[g]; b1 = OC_HCB2_SWO_B1[g]; b2 = OC_HCB2_SWO_B2[g]; b3 = OC_HCB2_SWO_B3[g]
    b4 = OC_HCB2_SWO_B4[g]; b5 = OC_HCB2_SWO_B5[g]; b6 = OC_HCB2_SWO_B6[g]
    si = g == 3 ? si_2 : si_1
    return ht/(1.0f0 + fexp(b0 + b1*ht + b2*ccfl + b3*flog(ba) + b4*(dbh/ht) + b5*si + b6*(og*og)))
end

"organon/crngrow.f MAXHCB_SWO — maximum height to crown base (species group `g`)."
@inline function oc_maxhcb_swo(g::Int, ht::Float32, ccfl::Float32)
    b0 = OC_MAXHCB_SWO_B0[g]; b1 = OC_MAXHCB_SWO_B1[g]; b2 = OC_MAXHCB_SWO_B2[g]
    b3 = OC_MAXHCB_SWO_B3[g]; lim = OC_MAXHCB_SWO_LIM[g]
    maxbr = b0 - b1*fexp(b2*fpow(ccfl/100.0f0, b3))
    maxbr > lim && (maxbr = lim)
    return maxbr*ht
end

"""
    oc_oldgro(spgrp, gdbh, ght, expan, dgro, hgro, deadexp, ntrees, ib, xind) -> OG

organon/mortality.f OLDGRO — old-growth indicator from the 5 largest big-6 trees/acre. `gdbh`/`ght`
are the GROWN DBH/HT (as CROWGRO sees `TDATAR`); `xind=-1` recovers the start-of-growth state
(subtract growth, add DEADEXP), `xind=0` the end-of-growth state.
"""
function oc_oldgro(spgrp::Vector{Int32}, gdbh::Vector{Float32}, ght::Vector{Float32},
                   expan::Vector{Float32}, dgro::Vector{Float32}, hgro::Vector{Float32},
                   deadexp::Vector{Float32}, ntrees::Int, ib::Int, xind::Float32)
    htcl = zeros(Float32, 100); dcl = zeros(Float32, 100); trcl = zeros(Float32, 100)
    @inbounds for i in 1:ntrees
        Int(spgrp[i]) <= ib || continue
        ht = ght[i] + xind*hgro[i]
        dbh = gdbh[i] + xind*dgro[i]
        ex = expan[i] - xind*deadexp[i]
        id = trunc(Int, dbh) + 1
        id > 100 && (id = 100)
        htcl[id] += ht*ex
        dcl[id]  += dbh*ex
        trcl[id] += ex
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
        return (totd/tottr)*(totht/tottr)/10000.0f0
    end
    return 0.0f0
end

"""
    organon_cr_swo(buf, dgro, hgro, spgrp, deadexp; si_1, si_2, cyclg=0) -> cr2::Vector{Float32}

The ORGANON SWO crown-recession pass for one cycle — `organon/crngrow.f` CROWGRO driven off the C1
`/ORGANON/` buffer and the C3 `dgro` / C4 `hgro` / C6 `deadexp`. Returns `cr2[i]` — the ORGANON
crown ratio after growth that `oc/crown.f:282` rounds into `CRNEW` for the valid ORGANON trees.

Sequence (VERSION=1, SCR shadow-crown = 0 on inventory): compute start (old DBH/HT) and end (new
DBH/HT, old HT for SSTATS per grow.f) stand stats, OLDGRO OG1/OG2, then per tree HCB_SWO at start &
end + MAXHCB_SWO, let the crown base rise by `HCBG = max(0, HCB2−HCB1)`, and cap the new crown base.
"""
function organon_cr_swo(buf::OrganonBuffer, dgro::Vector{Float32}, hgro::Vector{Float32},
        spgrp::Vector{Int32}, deadexp::Vector{Float32}; si_1::Float32, si_2::Float32, cyclg::Int=0)
    n = buf.ntrees
    ib = 5
    # grown state (TDATAR at CROWGRO time)
    gdbh = Vector{Float32}(undef, n); ght = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        gdbh[i] = buf.dbh1[i] + dgro[i]
        ght[i]  = buf.ht1or[i] + hgro[i]
    end
    # ORGANON mortality (C6) has already reduced TDATAR(4) to survivors before the end-of-growth
    # stats and CROWGRO run: TDATAR(4) = EXPAN·(1−PM) = EXPAN − DEADEXP (mortality.f:233-234).
    surv = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        surv[i] = buf.expan1[i] - deadexp[i]
    end
    # start-of-growth stats (old DBH, old HT) — computed at SOG, BEFORE mortality ⇒ original EXPAN
    sba1, _, _, ccfl1, ccfll1, _, _ = oc_sstats(spgrp, buf.dbh1, buf.ht1or, buf.expan1, n)
    # end-of-growth stats: new DBH, OLD HT (grow.f:169 SSTATS runs before the HT update at :184),
    # SURVIVOR EXPAN (post-mortality)
    sba2, _, _, ccfl2, ccfll2, _, _ = oc_sstats(spgrp, gdbh, buf.ht1or, surv, n)
    # OLDGRO sees the grown TDATAR(4)=survivor; XIND=-1 adds DEADEXP back (→ original), XIND=0 keeps survivor
    og1 = oc_oldgro(spgrp, gdbh, ght, surv, dgro, hgro, deadexp, n, ib, -1.0f0)
    og2 = oc_oldgro(spgrp, gdbh, ght, surv, dgro, hgro, deadexp, n, ib, 0.0f0)

    cr2 = zeros(Float32, n)
    scr = zeros(Float32, n)                        # shadow crown ratio (SCR1B = 0)
    @inbounds for i in 1:n
        g = Int(spgrp[i])
        pht = ght[i] - hgro[i]                     # start height (= old HT)
        pdbh = gdbh[i] - dgro[i]                   # start DBH   (= old DBH)
        sccfl1 = oc_get_ccfl(pdbh, ccfll1, ccfl1)
        hcb1 = oc_hcb_swo(g, pht, pdbh, sccfl1, sba1, si_1, si_2, og1)
        pcr1 = 1.0f0 - hcb1/pht
        phcb1 = (1.0f0 - pcr1)*pht
        ht = ght[i]; dbh = gdbh[i]
        sccfl2 = oc_get_ccfl(dbh, ccfll2, ccfl2)
        hcb2 = oc_hcb_swo(g, ht, dbh, sccfl2, sba2, si_1, si_2, og2)
        maxhcb = oc_maxhcb_swo(g, ht, sccfl2)
        pcr2 = 1.0f0 - hcb2/ht
        phcb2 = (1.0f0 - pcr2)*ht
        hcbg = phcb2 - phcb1
        hcbg < 0.0f0 && (hcbg = 0.0f0)
        ahcb1 = (1.0f0 - buf.cr1[i])*pht
        shcb1 = (1.0f0 - scr[i])*pht
        ahcb2 = ahcb1 + hcbg
        shcb2 = shcb1 + hcbg
        if ahcb1 > shcb1
            if ahcb1 > shcb2
                cr2[i] = 1.0f0 - ahcb1/ht          # (SCR path — inert when SCR=0)
            else
                cr2[i] = 1.0f0 - shcb2/ht
            end
        else
            if ahcb1 >= maxhcb
                cr2[i] = 1.0f0 - ahcb1/ht
            elseif ahcb2 >= maxhcb
                cr2[i] = 1.0f0 - maxhcb/ht
            else
                cr2[i] = 1.0f0 - ahcb2/ht
            end
        end
    end
    return cr2
end
