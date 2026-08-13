# =============================================================================
# mortality.jl (southeastalaska) — AK periodic mortality (ak/morts.f). Chunk 7.
#
# AK mortality is a per-tree LOGISTIC SURVIVAL model with a stand-density iterative multiplier
# (NOT the Zeide self-thin rate used by UT/NC, and NOT SEAMRT — ak/morts.f does not call it):
#   RIP = exp(X)/(1+exp(X)),  X = BM1 + BM2·DT + BM3·DT² + BM4·PTBAL + BM5·PTBAL/DT   (DT = max(D, DSURV))
#   WKI = P·(1 − RIP^FINT)                          [annual survival → period mortality]
# then repeated PASS multiplication: while the post-mortality stand SDI ≥ SDIMAX·PMSDIU or BA ≥ BAMAX,
# scale every record's kill by an increasing integer PASS (≤100) until the stand drops below the maxima.
# SDIMAX = BA-weighted SDIDEF (ak_point_zeide! XMAX); BAMAX = SDIMAX·0.5454154·PMSDIU (no user BAMAX).
# SDIMAX<5 ⇒ kill all. Zeide DQ10/SDI use G = DG/BRAT (outside-bark).
# =============================================================================

# ak/morts.f DATA — logistic survival coefficients + minimum survival diameter.
const AK_BM1 = Float32[5.106086,5.106086,5.334283,5.129246,4.585373,4.585373,5.129246,5.106086,4.935148,6.090056,5.120483,5.440943,4.585373,2.651991,2.987803,3.255821,3.255821,2.651991,2.833541,2.651991,2.651991,2.651991,2.651991]
const AK_BM2 = Float32[0.036136,0.036136,0.036136,0.402707,0.402707,0.402707,0.402707,0.036136,0.036136,0.036136,0.036136,0.036136,0.402707,0.402707,0.402707,0.402707,0.402707,0.402707,0.402707,0.402707,0.402207,0.402207,0.402207]
const AK_BM3 = Float32[-0.000634,-0.000634,-0.001034,-0.051774,-0.020537,-0.020537,-0.051774,-0.000634,-0.000850,-0.000745,-0.00105,-0.000942,-0.020537,-0.011288,-0.014038,-0.018658,-0.018658,-0.011288,-0.015719,-0.011288,-0.011288,-0.011288,-0.011288]
const AK_BM4 = Float32[0.0,0.0,0.0,-0.004365,-0.004414,-0.004414,-0.004365,0.0,0.0,0.0,0.0,0.0,-0.004414,-0.002557,-0.005462,-0.007307,-0.007307,-0.002557,-0.010755,-0.002557,-0.002557,-0.002557,-0.002557]
const AK_BM5 = Float32[-0.010747,-0.010747,-0.004540,0.0,0.0,0.0,0.0,-0.010747,-0.008488,-0.005608,-0.005131,-0.005375,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const AK_DSURV = Float32[0.50,0.50,0.50,0.10,0.10,0.10,0.10,0.50,0.50,0.50,0.50,0.50,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10]

function mortality!(s::StandState, ::SoutheastAlaska; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    dens = s.density
    # SDICAL(0,SDIMAX): stand BA-weighted SDIDEF; BAMAX = SDIMAX·0.5454154·PMSDIU when no user BAMAX.
    _, _, sdimax = ak_point_zeide!(s)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : sdimax * 0.5454154f0 * pmsdiu
    sdiupr = sdimax * pmsdiu
    dbhzeide = s.control.dbh_zeide
    # per-tree outside-bark increment G and the squared-diameter change (for the SDI/BA pass sums)
    g1 = Vector{Float32}(undef, n); ciobds = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        d = t.dbh[i]; sp = Int(t.species[i])
        bark = ak_bratio(sp, d)
        g = t.diam_growth[i] / bark
        g1[i] = g; ciobds[i] = 2f0 * d * g + g * g
    end
    # base per-tree kill wk2 via logistic survival
    wk2 = zeros(Float32, n)
    @inbounds for i in 1:n
        pr = t.tpa[i]; pr <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        ptbal = dens.point_bal[i]
        dt = d < AK_DSURV[sp] ? AK_DSURV[sp] : d
        x = AK_BM1[sp] + AK_BM2[sp]*dt + AK_BM3[sp]*dt*dt + AK_BM4[sp]*ptbal + AK_BM5[sp]*ptbal/dt
        rip = exp(x) / (1f0 + exp(x))
        rip > 0.99999f0 && (rip = 0.99999f0)
        rip < 0.001f0 && (rip = 0f0)
        wki = pr * (1f0 - rip^fint)                       # X (mort mult) = 1
        wki > pr && (wki = pr)
        sdimax < 5f0 && (wki = pr)                        # kill-all under collapsed max SDI
        wk2[i] = wki
    end
    # iterative PASS: scale kills up until stand SDI < SDIUPR and BA < BAMAX (ak/morts.f 55/59 loop)
    pass = 1
    while true
        sd2sqa = 0f0; sumdr10a = 0f0; ta = 0f0
        @inbounds for i in 1:n
            pr = t.tpa[i]
            wki = wk2[i] * pass; wki > pr && (wki = pr)
            surv = pr - wki
            t.dbh[i] >= dbhzeide || continue
            sd2sqa += surv * (t.dbh[i]*t.dbh[i] + ciobds[i])
            sumdr10a += surv * (t.dbh[i] + g1[i])^1.605f0
            ta += surv
        end
        if ta > 0f0
            dq10a = (sumdr10a / ta)^(1f0/1.605f0)
            baa = 0.005454154f0 * dq10a * dq10a * ta
            sdia = ta * (dq10a / 10f0)^1.605f0
            ((sdia < sdiupr && baa < bamax) || pass > 100) && break
        else
            break
        end
        pass += 1
    end
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        pr = t.tpa[i]
        wki = wk2[i] * pass; wki > pr && (wki = pr)
        killed[i] = wki
    end
    apply_fixmort!(s, killed, n, fint)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
