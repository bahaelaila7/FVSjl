# =============================================================================
# organon_htgf.jl — OC (Oregon Coast) FVS-native height growth for the NON-ORGANON (IORG=0) trees
# (chunk C9 step 4).
#
# Ported from oc/htgf.f (CASE DEFAULT), bin/FVSoc_buildDir/findag.f, bin/FVSoc_buildDir/htcalc.f.
# The IORG=0 trees keep the FVS-native HTGF (oc/htgf.f:105+): FINDAG finds the growth-effective age
# from the tree height + site curve, HTCALC gives the potential height at age+5, and a crown/relative-
# height modifier (SMHMOD) scales POTHTG. IFOR=9 (>5) ⇒ the REGION-6 site curves (HTCALC CASE DEFAULT).
#
# MEASURED bit-exact vs FVSoc_clean DEBUG-HTGF (ocmin): LP tree-1 HTG=6.909, the 0.1"-regen GF/DF
# HTG≈9. See docs/OC_VARIANT_PORT_AUDIT.md (C9 step 4).
# =============================================================================

"""
    oc_htcalc(sp, sindx, ag) -> HGUESS

htcalc.f (REGION 6, JFOR>5) — potential height at age `ag` for FVS species `sp`, site index `sindx`.
Per-species site curves: Hann-Scrivani (DF/GF/…, ×1.05 and JP/PP/… ×1.05), Dahms lodgepole (×1.10),
Dolph red fir, Powers black oak (×0.70), Porter-Wiant tanoak/madrone/red-alder (×0.80).
"""
@inline function oc_htcalc(sp::Int, sindx::Float32, ag::Float32)
    if (1 <= sp <= 4) || sp == 7 || sp == 8 || sp == 16 || sp == 22 || sp == 23 || sp == 25 || sp == 50
        # Hann-Scrivani (PC,IC,RC,WF,DF,WH,SP,BR,GS,OS,RW)
        toptrm = 1f0 - fexp(-fexp(-6.21693f0 + 0.281176f0*flog(sindx - 4.5f0) + 1.14354f0*flog(ag)))
        bottrm = 1f0 - fexp(-fexp(-6.21693f0 + 0.281176f0*flog(sindx - 4.5f0) + 1.14354f0*flog(50f0)))
        return (((sindx - 4.5f0)*toptrm/bottrm) + 4.5f0)*1.05f0
    elseif sp == 15 || (17 <= sp <= 20)
        # Hann-Scrivani (JP,WP,PP,MP,GP)
        toptrm = 1f0 - fexp(-fexp(-6.54707f0 + 0.288169f0*flog(sindx - 4.5f0) + 1.21297f0*flog(ag)))
        bottrm = 1f0 - fexp(-fexp(-6.54707f0 + 0.288169f0*flog(sindx - 4.5f0) + 1.21297f0*flog(50f0)))
        return (((sindx - 4.5f0)*toptrm/bottrm) + 4.5f0)*1.05f0
    elseif sp == 5 || sp == 6 || sp == 9
        # Dolph red fir (RF,SH,MH)
        term  = ag*fexp(ag*(-0.0440853f0))*1.41512f-6
        b     = sindx*term - 3.04951f6*term*term + 5.72474f-4
        term2 = 50f0*fexp(50f0*(-0.0440853f0))*1.41512f-6
        b50   = sindx*term2 - 3.04951f6*term2*term2 + 5.72474f-4
        return ((sindx - 4.5f0)*(1f0 - fexp(-b*fpow(ag, 1.51744f0))))/(1f0 - fexp(-b50*fpow(50f0, 1.51744f0))) + 4.5f0
    elseif (10 <= sp <= 14) || sp == 21
        # Dahms lodgepole (WB,KP,LP,CP,LM,WJ)
        return sindx*(-0.0968f0 + 0.02679f0*ag - 0.00009309f0*ag*ag)*1.10f0
    elseif sp == 24 || (30 <= sp <= 32) || sp == 35 || sp == 39 || sp == 40
        # Powers black oak (PY,WO,BO,VO,BU,DG,FL)
        term = sqrt(ag) - sqrt(50f0)
        return ((sindx*(1f0 + 0.322f0*term)) - 6.413f0*term)*0.70f0
    elseif sp == 42
        return (sindx/(0.204f0 + 39.787f0/ag))*0.80f0          # Porter-Wiant tanoak
    elseif (26 <= sp <= 29) || sp == 33 || sp == 37 || sp == 38
        return (sindx/(0.375f0 + 31.233f0/ag))*0.80f0          # Porter-Wiant madrone
    elseif sp == 34 || sp == 36 || sp == 41 || (43 <= sp <= 49)
        return (sindx/(0.649f0 + 17.556f0/ag))*0.80f0          # Porter-Wiant red alder
    else
        return 0f0
    end
end

"findag.f — growth-effective age/height: iterate AG by 2 until |HGUESS−H|≤2 or H<HGUESS. Returns (sitage, sitht)."
function oc_findag(sp::Int, sindx::Float32, h::Float32)
    agmax = 200f0                                              # IFOR>5 ⇒ AGEMAX=200
    ag = 2f0; incrng = 0; hguess = 0f0
    while true
        oldhg = hguess
        hguess = oc_htcalc(sp, sindx, ag)
        if hguess >= 1f0
            diff = abs(hguess - h)
            if diff <= 2f0 || h < hguess
                return ag, hguess
            end
            d2 = hguess - oldhg
            (oldhg != 0f0 && d2 >= 0.05f0) && (incrng = 1)
            if incrng == 1 && d2 < 0.05f0
                return ag, hguess
            end
        end
        ag += 2f0
        ag > agmax && return agmax, h
    end
end

"""
    oc_htgf_native(sp, h, icr, avh, pccf, sindx) -> HTG

oc/htgf.f CASE DEFAULT — native 5-yr height growth for an IORG=0 tree: FINDAG → SITAGE/SITHT,
POTHTG = HTCALC(SITAGE+5) − SITHT, then `HTG = POTHTG · 1.016605 · CRMOD · RHMOD` with
CRMOD=1−exp(−4.26558·CR), RHMOD=exp(2.54119·(relHt^0.250537−1)), relHt=min(H/AVH,1) forced to 1 when
PCCF<100. Floored at 0.1 (SCALE=XHT=1, HTCON=0 defaults).
"""
function oc_htgf_native(sp::Int, h::Float32, icr::Integer, avh::Float32, pccf::Float32, sindx::Float32)
    sitage, sitht = oc_findag(sp, sindx, h)
    pothtg = sitage > 200f0 ? 0.10f0 : (oc_htcalc(sp, sindx, sitage + 5f0) - sitht)
    cratio = Float32(icr)/100f0
    relht = avh > 0f0 ? h/avh : 0f0
    relht > 1f0 && (relht = 1f0)
    pccf < 100f0 && (relht = 1f0)
    crmod = 1f0 - fexp(-4.26558f0*cratio)
    rhmod = fexp(2.54119f0*(fpow(relht, 0.250537f0) - 1f0))
    htg = pothtg*(1.016605f0*crmod*rhmod)
    htg < 0.1f0 && (htg = 0.1f0)
    return htg
end
