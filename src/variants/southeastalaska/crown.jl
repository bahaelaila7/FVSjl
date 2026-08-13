# =============================================================================
# crown.jl (southeastalaska) — AK crown ratio (ak/crown.f + ak/dubscr.f) + point-Zeide density
# (ak/sdical.f SDICAL/SDICLS). Chunk 5.
#
# AK crown is a LOGISTIC on point relative density (PRD), NOT a Weibull rank model:
#   X = CRINT + CRHDR·ln(H·12/D) + CRRD·PRD + CRDQMD·(D/QMDPLT);  CR = 1/(1+exp(X));  clamp[0.05,0.95].
# with PRD = ZRD(point)/XMAXPT(point) — the per-point Zeide relative density. XMAXPT = BA-weighted
# SDIDEF per point (SDICAL); ZRD = Σ TPACRE·(D/10)^1.605 over the point (SDICLS, Zeide method), with
# TPACRE = PROB·(PI−NONSTK). QMDPLT = point QMD (floored at 1).
#   D<1 (or missing CR): ak/dubscr.f — same logistic + a bounded normal random error (bachlo).
#   Change limit: +3%/yr increase, −1%/yr decrease (ak/crown.f:322-323, faster than the usual 1%/1%).
#   Crown-length cap (CRMAX) so new crown ≤ what all HTG could add. Floor [10,95] (CRNMLT=1).
# =============================================================================

# ak/crown.f CRCONS DATA — logistic crown coefficients (intercept, ln HDR, PRD, D/QMD).
const AK_CRINT  = Float32[-2.818788,-2.818788,0.655733,6.469047,-4.02647,-4.02647,-2.781735,-2.818788,-0.425173,-0.449694,-0.364665,-0.47458,-4.026470,0.977291,-2.799049,-1.268650,-1.268650,0.977291,0.802511,0.977291,-2.406592,-2.406592,0.977291]
const AK_CRHDR  = Float32[0.623633,0.623633,-0.064716,-1.495754,0.655230,0.655230,0.438826,0.623633,0.295262,0.104368,0.116994,0.141604,0.655230,-0.314026,0.660530,0.147556,0.147556,-0.314026,-0.143680,-0.314026,0.347417,0.347417,-0.314026]
const AK_CRRD   = Float32[0.952948,0.952948,0.439411,0.913909,0.914077,0.914077,0.970251,0.952948,0.732675,0.557638,0.573079,0.590263,0.914077,0.848483,0.987503,1.092390,1.092390,0.848483,1.023227,0.848483,0.952753,0.952753,0.848483]
const AK_CRDQMD = Float32[-0.174337,-0.174337,-0.174337,-0.051268,-0.051268,-0.051268,-0.051268,-0.174337,-0.174337,-0.174337,-0.174337,-0.174337,-0.051268,-0.051268,-0.174337,-0.051268,-0.051268,-0.051268,-0.051268,-0.051268,-0.051268,-0.051268,-0.051268]
# ak/dubscr.f DATA CRSD — per-species crown-ratio regression standard error (random-error draw).
const AK_CRSD   = Float32[0.1719,0.1719,0.159,0.1645,0.1542,0.1542,0.1894,0.1719,0.1568,0.1811,0.1582,0.155,0.1542,0.1447,0.1447,0.144,0.144,0.1605,0.1356,0.1605,0.0938,0.0938,0.1605]

"""
    ak_point_zeide!(s) -> (xmaxpt, zrd, xmax)

ak/sdical.f SDICAL(XMAX,XMAXPT) + SDICLS(ZRD). Per point: XMAXPT = BA-weighted SDIDEF; ZRD = Zeide SDI
Σ PROB·(PI−NONSTK)·(D/10)^1.605 (D ≥ DBHZEIDE). XMAX = stand BA-weighted SDIDEF (or 1 if no BA). The
BA weights use raw PROB (0.0054542·D²·PROB) — the PI/GROSPC scaling cancels in the XMAXPT ratio.
"""
function ak_point_zeide!(s::StandState)
    p, t = s.plot, s.trees
    sdidef = p.sp_sdi_def
    # SDICAL/SDICLS sum the UNCHANGED DBH(I) — during DGF calibration the diameters are backdated in
    # DIAM(I) only, so use the stashed CURRENT dbh (calib.calib_dbh) when present. Empty otherwise
    # (growth path + crown init both read the live t.dbh).
    dbharr = isempty(s.calib.calib_dbh) ? t.dbh : s.calib.calib_dbh
    npt = 0
    @inbounds for i in 1:t.n; npt = max(npt, Int(t.plot_id[i])); end
    ptba  = zeros(Float32, npt)                 # PNTBA(pt) = Σ TREEBA
    ptsdi = zeros(Float32, npt)                 # Σ SDIDEF·TREEBA
    zrd   = zeros(Float32, npt)
    totba = 0f0; xsdi = 0f0
    pifac = p.pi - Float32(p.nonstockable)      # (PI − NONSTK)
    dbhzeide = s.control.dbh_zeide
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        d = dbharr[i]; sp = Int(t.species[i]); ip = Int(t.plot_id[i])
        treeba = 0.0054542f0 * d * d * t.tpa[i]
        totba += treeba; xsdi += sdidef[sp] * treeba
        if 1 <= ip <= npt
            ptba[ip]  += treeba
            ptsdi[ip] += sdidef[sp] * treeba
            d >= dbhzeide && (zrd[ip] += t.tpa[i] * pifac * (d / 10f0)^1.605f0)
        end
    end
    xmax = totba <= 0f0 ? 1f0 : xsdi / totba
    xmaxpt = Vector{Float32}(undef, npt)
    @inbounds for ip in 1:npt
        xmaxpt[ip] = ptba[ip] == 0f0 ? xmax : ptsdi[ip] / ptba[ip]
    end
    return xmaxpt, zrd, xmax
end

# ak/dubscr.f — logistic crown for D<1 / missing CR with a bounded normal random error.
@inline function ak_dubscr(rng, sp::Integer, d::Float32, h::Float32, prd::Float32, qmd::Float32)::Float32
    cr = AK_CRINT[sp] + AK_CRHDR[sp] * log(h * 12f0 / d) + AK_CRRD[sp] * prd + AK_CRDQMD[sp] * (d / qmd)
    sd = AK_CRSD[sp]
    fcr = bachlo(rng, 0f0, sd)
    while abs(fcr) > sd                          # dubscr.f label 10: reject |FCR|>SD
        fcr = bachlo(rng, 0f0, sd)
    end
    cr = 1f0 / (1f0 + exp(cr + fcr))
    cr > 0.95f0 && (cr = 0.95f0); cr < 0.05f0 && (cr = 0.05f0)
    return cr
end

function crown_ratio_update!(s::StandState, ::SoutheastAlaska; fint::Float32 = 10.0f0,
                             lstart::Bool = false, crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    xmaxpt, zrd, _ = ak_point_zeide!(s)
    npt = length(xmaxpt)
    pb = s.density.point_ba; ptpa = s.density.point_tpa
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); (sp < 1 || sp > 23) && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue           # inventory crown present ⇒ keep (ak/crown.f:217)
        d = t.dbh[i]; h = t.height[i]; ip = Int(t.plot_id[i])
        baplt  = (1 <= ip <= length(pb))   ? pb[ip]   : 0f0
        tpaplt = (1 <= ip <= length(ptpa)) ? ptpa[ip] : 0f0
        qmdplt = tpaplt > 0f0 ? sqrt((baplt / tpaplt) / 0.005454f0) : 1f0
        qmdplt <= 1f0 && (qmdplt = 1f0)
        xmp = (1 <= ip <= npt) ? xmaxpt[ip] : 0f0
        prd = xmp <= 0f0 ? 0f0 : ((1 <= ip <= npt) ? zrd[ip] : 0f0) / xmp
        local icri::Int
        if d < 1f0 && lstart
            cr = ak_dubscr(s.rng, sp, d, h, prd, qmdplt)
            icri = trunc(Int, cr * 100f0 + 0.5f0)      # CRNMLT=1, DLOW/DHI defaults ⇒ no multiplier
        else
            x = AK_CRINT[sp] + AK_CRHDR[sp] * log(h * 12f0 / d) + AK_CRRD[sp] * prd +
                AK_CRDQMD[sp] * (d / qmdplt)
            x = 1f0 / (1f0 + exp(x))
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = Float32(trunc(Int, x * 100f0 + 0.5f0))
            if !(lstart || icr == 0)               # change limit +3%/yr, −1%/yr (ak/crown.f:319-334)
                chg = crnew - Float32(icr)
                pdifpy = chg / Float32(icr) / fint
                pdifpy > 0.03f0 && (chg = Float32(icr) * 0.03f0 * fint)
                pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
                crnew = Float32(icr) + chg           # CRNMLT=1
            end
            icri = trunc(Int, crnew + 0.5f0)
            if !(lstart || icr == 0)                 # crown-length cap CRMAX (ak/crown.f:346-358)
                crln = h * Float32(icr) / 100f0
                htg = t.ht_growth[i]
                crmax = (crln + htg) / (h + htg) * 100f0
                Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
                (icri < 10) && (icri = trunc(Int, crmax + 0.5f0))   # CRNMLT=1
            end
        end
        icri > 95 && (icri = 95)
        icri < 10 && (icri = 10)                     # CRNMLT=1 (ak/crown.f:393)
        icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
