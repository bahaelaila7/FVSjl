# =============================================================================
# mortality.jl (southeastalaska) — AK periodic mortality (ak/morts.f). Chunk 7.
#
# AK mortality is a per-tree LOGISTIC SURVIVAL model with a stand-density iterative multiplier
# (NOT the Zeide self-thin rate used by UT/NC, and NOT SEAMRT — ak/morts.f does not call it):
#   RIP = exp(X)/(1+exp(X)),  X = BM1 + BM2·DT + BM3·DT² + BM4·PTBAL + BM5·PTBAL/DT   (DT = max(D, DSURV))
#   WKI = P·(1 − RIP^FINT)                          [annual survival → period mortality]
# then repeated PASS multiplication: while the post-mortality stand SDI ≥ SDIMAX·PMSDIU or BA ≥ BAMAX,
# scale every record's kill by an increasing integer PASS (≤100) until the stand drops below the maxima.
# SDIMAX = BA-weighted SDIDEF (point_zeide! XMAX); BAMAX = SDIMAX·0.5454154·PMSDIU (no user BAMAX).
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
    iyc = Int(current_cycle_year(s))
    # morts.f:192-193 PMSDIU as a fraction; morts.f:202 SDICAL(0,SDIMAX): stand BA-weighted SDIDEF, and sdical.f:203-209
    # BAMAX = XMAX·0.5454154·PMSDIU unless a user BAMAX is in effect — then XMAX = BAMAX/(0.5454154·PMSDIU).
    xmaxpt_m, _, sdimax = point_zeide!(s)
    ak_es_stash_xmaxpt!(s, xmaxpt_m, sdimax)   # morts.f:202 SDICAL XMAXPT — read by the next DENSE (ESTAB PRDA)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdiu > 1f0 && (pmsdiu = pmsdiu / 100f0)
    local bamax::Float32
    if s.control.ba_max > 0f0
        bamax = s.control.ba_max
        sdimax = bamax / (0.5454154f0 * pmsdiu)
    else
        bamax = sdimax * 0.5454154f0 * pmsdiu
    end
    sdiupr = sdimax * pmsdiu
    dbhzeide = s.control.dbh_zeide
    # morts.f DO 20: G1 = DG/BARK for every record; CIOBDS1 only for the Zeide-eligible (D ≥ DBHZEIDE) ones.
    g1 = Vector{Float32}(undef, n); ciobds = zeros(Float32, n)
    @inbounds for i in 1:n
        d = t.dbh[i]; sp = Int(t.species[i])
        g = t.diam_growth[i] / ak_bratio(sp, d)
        g1[i] = g
        d >= dbhzeide && (ciobds[i] = 2f0 * d * g + g * g)
    end
    # morts.f DO 50/40 (species-major): logistic survival RIP, MORTMULT window, the ESTAB best-tree IESTAT guard,
    # the size-cap mortality, then WKI ≤ P and SDIMAX<5 ⇒ kill all.
    wk2 = zeros(Float32, n)
    @inbounds for i in species_major_order(s)
        pr = t.tpa[i]
        sp = Int(t.species[i]); d = t.dbh[i]
        pr <= 0f0 && continue
        ptbal = dens.point_bal[i]
        dtemp = d < AK_DSURV[sp] ? AK_DSURV[sp] : d
        rip = AK_BM1[sp] + AK_BM2[sp]*dtemp + AK_BM3[sp]*dtemp*dtemp + AK_BM4[sp]*ptbal + AK_BM5[sp]*ptbal/dtemp
        rip = fexp(rip) / (1 + fexp(rip))
        x = active_mort_mult(s.control, sp, iyc, d)
        if t.iestat[i] > 0
            iyc >= t.iestat[i] && (t.iestat[i] = Int32(0))
            xchk = Float32(t.iestat[i] - iyc) / fint
            xchk > 1f0 && (xchk = 1f0); xchk < 0f0 && (xchk = 0f0)
            x = x * (1f0 - xchk)
        end
        rip > 0.99999f0 && (rip = 0.99999f0)
        rip < 0.001f0 && (rip = 0f0)
        wki = pr * (1f0 - fpow(rip, fint)) * x                        # RIP**FINT (REAL**REAL = powf)
        g = (t.diam_growth[i] / ak_bratio(sp, d)) * (fint / 10f0)
        cap = s.control.sp_size_cap
        if (d + g) >= cap[sp, 1] && unsafe_trunc(Int, cap[sp, 3]) != 1
            wki = max(wki, pr * cap[sp, 2] * fint / 10f0)
        end
        wki > pr && (wki = pr)
        sdimax < 5f0 && (wki = pr)
        wk2[i] = wki
    end
    # morts.f 59/60/55 loop: scale every record's kill by PASS until the post-mortality stand is below
    # SDIUPR and BAMAX (BAA from the QUADRATIC DQ10A, before the Zeide override), PASS ≤ 100.
    pass = 1
    while true
        sd2sqa = 0f0; sumdr10a = 0f0; ta = 0f0
        @inbounds for i in 1:n
            wki = wk2[i] * Float32(pass); wki > t.tpa[i] && (wki = t.tpa[i])
            t.dbh[i] >= dbhzeide || continue
            sd2sqa += (t.tpa[i] - wki) * (t.dbh[i] * t.dbh[i] + ciobds[i])
            sumdr10a += (t.tpa[i] - wki) * fpow(t.dbh[i] + g1[i], 1.605f0)
            ta += (t.tpa[i] - wki)
        end
        if !(ta > 0f0)                  # TA=0 ⇒ DQ10A/BAA/SDIA = NaN: only the PASS>100 exit fires (morts.f:371-387)
            pass > 100 && break
            pass += 1; continue
        end
        dq10a = sqrt(sd2sqa / ta)
        baa = 0.005454154f0 * dq10a * dq10a * ta
        dq10a = fpow(sumdr10a / ta, 1f0 / 1.605f0)
        sdia = ta * fpow(dq10a / 10f0, 1.605f0)
        ((sdia < sdiupr && baa < bamax) || pass > 100) && break
        pass += 1
    end
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        pr = t.tpa[i]
        wki = pass > 1 ? wk2[i] * Float32(pass) : wk2[i]; wki > pr && (wki = pr)
        killed[i] = wki
    end
    apply_fixmort!(s, killed, n, fint)
    # MISMRT (mistoe.f:522 → mismrt.f:185-191, misintak.f APMC): WK2=MAX(WK2,PROB·rate); deferred to post-TRIPLE on
    # a tripling cycle (dm_mrt_defer). Inert unless a record carries DMR.
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
