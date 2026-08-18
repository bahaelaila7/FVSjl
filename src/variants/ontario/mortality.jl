# =============================================================================
# ontario/mortality.jl — ON periodic mortality (canada/on/morts.f + canada/on/varmrt.f).
#
# CHUNK 8b. `mortality!(s, ::Ontario; fint)` is a fresh port of ON's MORTS — it is NOT the shared
# southern Pretzsch mortality. ON (Penner 2010, FVS_Mortality_Report) combines:
#   • DENSITY mortality on a Penner MAXIMUM-SDI line, LN(TPHmetric)=A+B·LN(QMDcm), whose (A,B) are
#     chosen by the species with the most basal area (INDX_BA → MSB_MAP → one of 13 SDI_INT/SDI_SLP
#     equations). Density kills back toward the 85%/55%-of-max envelope (PMSDIU/PMSDIL); at cyc1 on a
#     hyper-dense metric-design stand T ≫ T85D0 so TN10 = T85D10 directly.
#   • BACKGROUND mortality RI = 0.5/(1+exp(B0+B1·D)) with B0/B1 from a 4-group map (BKG_MAP → PMSC/PMD).
#   • DISTRIBUTION of the density excess by ON's own VARMRT — the Penner INDIVIDUAL-TREE mortality
#     model (18-equation logistic, ITM_MAP → MB0..MB7, metric DM/DQM/BALM/SIM) driving a geometric
#     progression toward the target TPA (identical progression math to sn/varmrt.f, different EFFTR).
# Everything metric (INtoCM/ACRtoHA/FT2pACRtoM2pHA/FTtoM). G = (DG/BARK)·(FINT/10) — ON is a 10-yr
# variant (YR=10), BARK = on_bratio(sp, DBH, HT) on the CURRENT dims.
#
# Transcendentals via glibc single-precision (on_expf/on_logf/on_powf) to match gfortran.
# VALIDATED bit-exact vs instrumented FVSon_g16 (morts.f g16 single-.o swap, MORTDUMP) on ont01
# cyc1: per-tree WK2 8/8, T/DIA0/D10/SDIMAX/BAMAX/INDX_BA/INDX/TN10/RN all reproduced.
# =============================================================================

const ON_ACRtoHA = 0.4046945f0   # common/METRIC.F77
const ON_HAtoACR = 2.471f0

# --- MAX SDI mortality coefficients (Penner Eqn 9), morts.f DATA; mapped by MSB_MAP -------------
const ON_SDI_INT = Float32[
    7.77029, 6.25277, 5.98521, 5.85622, 4.85176,
    8.99134, 4.49037, 6.90979, 6.97823, 6.32553,
    7.40806, 8.26319, 6.69982]
const ON_SDI_SLP = Float32[
   -0.61566, -0.40357, -0.38783, -0.3828, -0.25916,
   -0.8304, -0.20807, -0.49564, -0.54518, -0.46181,
   -0.58078, -0.74894, -0.51312]
# --- Mature-stand-boundary coefficients (Gary Dixon ~ Penner Eqn 10), mapped by MSB_MAP --------
const ON_MSB_DBH = Float32[
    3.00, 3.30, 2.00, 2.90, 2.75, 3.70, 2.80, 3.10, 3.50, 2.90,
    3.10, 3.30, 2.40]
const ON_MSB_INT = Float32[
   14.64825, 20.51662, 16.27566, 21.06262, 32.03489,
   15.99204, 29.68404, 17.91661, 13.23997, 20.46762,
   15.78771, 18.28047, 13.65976]
const ON_MSB_SLP = Float32[
   -2.30, -4.00, -3.00, -4.60, -8.70, -2.60, -7.70, -3.30, -1.96, -4.50,
   -2.70, -3.60, -2.20]
# 72 FVS species -> one of 13 SDI/MSB equations (morts.f MSB_MAP)
const ON_MSB_MAP = Int[
    3, 2, 2, 2, 1, 4, 4, 6, 5, 3,
    8, 7, 2, 5,11,11,13,10,10,11,
   10,10,10,12,10, 9,10, 9,11,12,
   12,12,12,12,12,12,11,12,12,13,
   13,13,11,12,12,12,10,13,12,10,
   10,10,10,10,12,10,11,11,10,12,
   11,11,11,12,13,13,11,12, 3, 1,
    4, 5]
# --- Background-mortality constants (morts.f), 4 groups, mapped by BKG_MAP ----------------------
const ON_PMSC = Float32[5.1676998, 9.6942997, 5.5876999, 5.9617000]   # B0
const ON_PMD  = Float32[-0.0077681, -0.0127328, -0.0053480, -0.0340128] # B1
const ON_BKG_MAP = Int[
    1, 3, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 3, 3, 4, 1, 4, 1, 1, 4,
    1, 1, 1, 4, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 1, 4, 4, 1, 1, 4, 4, 4, 4, 1, 4, 4, 1, 4, 4, 1, 4,
    4, 4, 4, 1, 1, 4, 1, 4, 1, 3, 1, 1]

# --- Penner individual-tree mortality (varmrt.f), 18 equations, mapped by ITM_MAP --------------
const ON_MB0 = Float32[
    1.56759, -1.75438, -3.67221, -0.3542, -0.83917,
   -0.22749, -2.70186, -1.05177, -3.27602, -2.27416,
   -3.67336, -3.89531, -4.5782, -6.64434, -3.39803,
   -1.84346, -4.63558, -3.60002]
const ON_MB1 = Float32[
   -0.35954, -0.13439, 0.13553, -0.29402, -0.075347,
   -0.17368, -0.15224, 0.18688, -0.06409, -0.75413,
   -0.015301, -0.39851, -0.11127, -0.60664, -0.016078,
   -0.54492, -0.10752, 0.035019]
const ON_MB2 = Float32[
  -18.2223, 0.05, 0.0, -4.41248, -3.61772,
   -2.87478, -1.77768, -7.79339, -2.26039, -5.44546,
   -1.04581, -0.028481, 2.24759, 0.0, -1.73506,
    0.0, -0.17508, -1.79863]
const ON_MB3 = Float32[
    0.004729459, 0.003717782, -0.007633363, 0.007291, 0.00429437,
    0.013202, 0.003566943, -0.000771621, 0.003381, 0.030003,
    0.000577445, 0.00392894, 0.0, 0.009326974, 0.00275289,
    0.005057465, 0.002772803, 0.000206976]
const ON_MB4 = Float32[
   -0.73227, -2.54341, -2.98566, -7.32994, -4.10504,
   -5.39914, -0.53056, -6.43672, -1.18764, -1.28947,
   -0.10899, 2.34899, 1.98514, 4.61152, -1.8512,
    2.90385, 0.25158, -1.3242]
const ON_MB5 = Float32[
   -0.046773, 0.0, 0.029072, 0.04226, 0.018135,
   -0.038558, -0.057377, 0.049045, -0.030958, 0.031923,
   -0.018726, 0.0, 0.0, 0.0, 0.01346,
    0.0, 0.016405, 0.0]
const ON_MB6 = Float32[
    0.0, -1.833, 0.0, 0.5563, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.2,
    0.0, -0.1917, 0.0]
const ON_MB7 = Float32[
    0.0357, 0.034951, 0.026167, 0.085012, 0.064902,
    0.11252, 0.050818, 0.032482, 0.0381, 0.11331,
    0.020787, 0.078816, 0.024981, 0.13838, 0.019736,
    0.0, 0.055896, 0.018238]
const ON_ITM_MAP = Int[
    5, 4, 3, 4, 1, 7, 8,11, 9, 5,
    9,11, 3, 9,18,18,17,13,13,15,
   13,13,13,16,13,12,13,14,18,16,
   16,16,16,16,16,16,15,16,16,17,
   17,17,15,16,16,16,13,17,16,13,
   13,13,13,13,16,13,15,15,13,16,
   15,15,15,16,17,17,15,16, 6, 2,
    8, 10]

"""
    on_varmrt!(killed, efftr, temwk2, s, t, n, tokill, d10) -> sumkil

Ontario VARMRT (canada/on/varmrt.f): distribute `tokill` TPA across the `n` live records by a
geometric progression weighted by the Penner individual-tree mortality efficiency EFFTR (18-eqn
logistic in metric DM/DQM/BALM/SIM, `d10` = the density-loop QMD). Fills `killed[i]` (FVS WK2),
returns the total killed. Bit-exact vs FVSon_g16.
"""
function on_varmrt!(killed::AbstractVector{Float32}, efftr::AbstractVector{Float32},
                    temwk2::AbstractVector{Float32}, s::StandState, t::TreeList, n::Int,
                    tokill::Float32, d10::Float32)
    p = s.plot
    fill!(view(killed, 1:n), 0f0)
    fill!(view(temwk2, 1:n), 0f0)
    tokill <= 0f0 && return 0f0
    ba = p.basal_area
    # --- EFFTR: Penner individual-tree mortality rate (varmrt.f DO I=1,ITRN) ---
    pass1 = 0f0
    dqm = max(1f0, d10 * ON_INtoCM)
    @inbounds for i in 1:n
        jspc = Int(t.species[i])
        sim = (jspc >= 1 && jspc <= length(p.sp_site_index)) ? p.sp_site_index[jspc] * ON_FTtoM : 0f0
        j = ON_ITM_MAP[jspc]
        dm = t.dbh[i] * ON_INtoCM
        balm = (1f0 - (t.crown_ratio[i] / 100f0)) * ba * ON_FT2pACRtoM2pHA
        x = ON_MB0[j] + (ON_MB1[j] * dm) + (ON_MB2[j] / (dm + ON_MB6[j])) +
            (ON_MB3[j] * dm * dm) + (ON_MB4[j] * dm / dqm) + (ON_MB5[j] * balm) + (ON_MB7[j] * sim)
        xc = min(max(x, -88f0), 88f0)
        peff = 1f0 / (1f0 + on_expf(xc))     # SURV/YR
        peff = 1f0 - peff                    # MORT/YR
        peff > 1f0 && (peff = 1f0)
        peff < 0.01f0 && (peff = 0.01f0)
        efftr[i] = peff
        pass1 += t.tpa[i] * peff
    end
    npass = trunc(Int, tokill / pass1) + 1
    sumkil = 0f0; temkil = tokill; short_v = 0f0; jpass = 0
    while true
        jpass += 1
        jpass > 1 && (temkil = short_v)
        iswtch = 0
        temsum = 0f0
        while true                                   # adjust NPASS into [0.8,1.2] (varmrt.f 105)
            temsum = 0f0
            @inbounds for i in 1:n
                tpalft = t.tpa[i] - killed[i]
                if tpalft > 0f0
                    temwk2[i] = -tpalft * ((1f0 - efftr[i])^npass - 1f0)   # NPASS is INTEGER (gfortran integer power)
                    temsum += temwk2[i]
                end
            end
            minstp = npass > 50 ? 5 : (npass > 20 ? 2 : 1)
            adjust = temsum == 0f0 ? 1f0 : temkil / temsum
            if adjust < 0.8f0
                iswtch == 2 && break
                npass -= max(minstp, trunc(Int, (temsum - temkil) / pass1)); iswtch = 1
                npass <= 0 && break
                continue
            elseif adjust > 1.2f0
                iswtch == 1 && break
                npass += max(minstp, trunc(Int, (temkil - temsum) / pass1)); iswtch = 2
                continue
            end
            break
        end
        short_v = 0f0
        adjust = temsum == 0f0 ? 1f0 : temkil / temsum
        @inbounds for i in 1:n
            tpalft = t.tpa[i] - killed[i]
            tpalft < 0.00001f0 && continue
            xkill = temwk2[i] * adjust
            if (t.tpa[i] - killed[i] - xkill) <= 0.00001f0
                xk = t.tpa[i] - killed[i]
                short_v += xkill - xk
                pass1 -= efftr[i]
                killed[i] += xk; sumkil += xk
            else
                killed[i] += xkill; sumkil += xkill
            end
        end
        short_v <= 0f0 && break
        npass = trunc(Int, short_v / pass1) + 1
    end
    return sumkil
end

# NPASS exponent uses a Float32 power to match gfortran `(1.0-EFFTR)**NPASS` (NPASS is INTEGER in
# Fortran, so `**` is a whole-number power — but on ont01 the values agree with on_powf; using the
# integer power `^` is the exact gfortran semantics, so use that instead of on_powf for the base loop).

"""
    mortality!(s, ::Ontario; fint, book_snags) — ON MORTS (canada/on/morts.f).
"""
function mortality!(s::StandState, ::Ontario; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    ctl = s.control
    n = t.n
    n == 0 && return s
    isct = ctl.sp_count_tab; ind1 = s.scratch.idx1
    zeide = ctl.zeide_sdi
    dbhsdi = ctl.dbh_sdi
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    rmsqd = p.qmd
    ba = p.basal_area
    tphmax = 35000f0 * ON_HAtoACR

    dens = s.density
    # morts.f:211 RMSQD==0 -> reset persisted line
    rmsqd == 0f0 && (dens.mort_intercept = 0f0; dens.mort_slope = 0f0)

    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    efftr  = @view s.scratch.mort_efftr[1:n]
    temwk2 = @view s.scratch.mort_temwk2[1:n]

    # --- DO 20/12: per-species BA + stand sums, over DBH>=DBHSDI, IND1 (species-sorted) order ---
    basp = fill(0f0, MAXSP)
    tt = 0f0; sdq0 = 0f0; sd2sq = 0f0; sumdr0 = 0f0; sumdr10 = 0f0
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        for k in i1:i2
            i = ind1[k]
            killed[i] = 0f0
            d = t.dbh[i]
            d < dbhsdi && continue
            pr = t.tpa[i]
            bark = on_bratio(sp, d, t.height[i])
            g = (t.diam_growth[i] / bark) * (fint / 10f0)
            ciobds = 2f0 * d * g + g * g
            tem = pr * (d * d + ciobds)
            sd2sq += tem
            sdq0  += pr * d * d
            if zeide
                sumdr10 += pr * on_powf(d + g, 1.605f0)
                sumdr0  += pr * on_powf(d, 1.605f0)
            end
            tt += pr
            basp[sp] += tem * 0.0054542f0
        end
    end

    # INDX_BA = species with the most BA (morts.f:279-286)
    indx_ba = 1; tem = basp[1]
    @inbounds for sp in 1:MAXSP
        basp[sp] > tem && (indx_ba = sp; tem = basp[sp])
    end
    # dominant-species change -> reset self-thinning line (morts.f:288-292)
    if indx_ba != Int(dens.mort_ibasp)
        dens.mort_ibasp = Int32(indx_ba)
        dens.mort_intercept = 0f0; dens.mort_slope = 0f0
    end
    # trajectory change since last cycle -> reset (morts.f:300-305)
    if Int(ctl.cycle) > 1 && dens.tpa_mort > 0f0 && abs(tt - dens.tpa_mort) > 1f0
        dens.mort_intercept = 0f0; dens.mort_slope = 0f0
    end

    if tt < 1f0
        # morts.f GO TO 45: no density/background mortality; FIXMORT only.
        apply_fixmort!(s, killed, n, fint)
        book_snags && book_mortality_snags!(s, killed, n, fint)
        @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
        return s
    end

    dq0 = sqrt(sdq0 / tt); dq10 = sqrt(sd2sq / tt)
    dr0 = zeide ? on_powf(sumdr0 / tt, 1f0 / 1.605f0) : 0f0
    dr10 = zeide ? on_powf(sumdr10 / tt, 1f0 / 1.605f0) : 0f0
    dia0 = zeide ? dr0 : dq0
    d10 = zeide ? dr10 : dq10

    sdimax = stand_sdimax(s)                          # SDICAL(0,SDIMAX) — Reineke, climate check
    bamax = sdimax * 0.5454154f0 * pmsdiu             # BAMAX (DENSE), morts.f BA-check target

    ipath = 0
    d10cur = d10
    if sdimax < 5f0
        tn10 = 0f0                                    # climate kill: remove all trees
    else
        # --- MORTS QMD-convergence loop (morts.f label 10 .. GO TO 10) ---
        indx = ON_MSB_MAP[indx_ba]
        sdiint = -ON_SDI_INT[indx] / ON_SDI_SLP[indx]
        sdislp = 1f0 / ON_SDI_SLP[indx]
        constv = on_expf(sdiint)
        tn10 = 0f0
        @inbounds for _ in 1:10
            # DIA0<0.3 reset (morts.f:330-334) — inside the loop, uses the current D10
            if dia0 < 0.3f0
                d10cur = 0.3f0 + d10cur - dia0
                dia0 = 0.3f0
            end
            tmd0 = ON_ACRtoHA * on_expf(sdiint + sdislp * on_logf(dia0 * ON_INtoCM))
            tmd0 > tphmax && (tmd0 = tphmax)
            t55d0 = tmd0 * pmsdil; t85d0 = tmd0 * pmsdiu
            tmd10 = ON_ACRtoHA * on_expf(sdiint + sdislp * on_logf(d10cur * ON_INtoCM))
            tmd10 > tphmax && (tmd10 = tphmax)
            t55d10 = tmd10 * pmsdil; t85d10 = tmd10 * pmsdiu
            tn10 = _on_tn10!(dens, tt, dia0, d10cur, constv, sdislp, pmsdil, pmsdiu,
                             t85d0, t55d0, t85d10, t55d10)
            # 271: bound TN10
            tn10 > tt && (tn10 = tt)
            tn10 < 0.1f0 && (tn10 = 0f0)
            rn = 1f0 - on_powf(1f0 - ((tt - tn10) / tt), 1f0 / fint)

            # DO 50: per-tree background/SDI rate into killed (WK2). RIP is uniform per stand.
            tem_g = constv * on_powf(d10cur * ON_INtoCM, sdislp)
            tem_g > tphmax && (tem_g = tphmax)
            tem_g = tem_g * pmsdil * ON_ACRtoHA
            density_on = !(tt <= tem_g || rn <= 0f0)
            @inbounds for sp in 1:MAXSP
                i1 = isct[sp, 1]; i1 == 0 && continue
                i2 = isct[sp, 2]
                b0 = ON_PMSC[ON_BKG_MAP[sp]]; b1 = ON_PMD[ON_BKG_MAP[sp]]
                for k in i1:i2
                    i = ind1[k]
                    killed[i] = 0f0
                    pr = t.tpa[i]; pr <= 0f0 && continue
                    d = t.dbh[i]
                    ri = 1f0 / (1f0 + on_expf(b0 + b1 * d))
                    ri = 0.5f0 * ri
                    rip = density_on ? rn : ri
                    rip > 1f0 && (rip = 1f0)
                    xm = 1f0
                    # MORTMULT (inert without keyword) — density regime forces X=1 (morts.f:609)
                    density_on || (xm = active_mort_mult(ctl, sp, current_cycle_year(s), d))
                    wki = pr * (1f0 - on_powf(1f0 - rip, fint)) * xm
                    wki > pr && (wki = pr)
                    killed[i] = wki
                end
            end

            # distribute (morts.f:637-642): density excess by VARMRT, else the background total
            sumtre = density_on ? (tt - tn10) : 0f0
            sumtre < 0f0 && (sumtre = 0f0)
            tokill = density_on ? sumtre : sum(@view killed[1:n])
            if tn10 >= 0.1f0
                on_varmrt!(killed, efftr, temwk2, s, t, n, tokill, d10cur)
            end
            # background regime has no D10 dependence -> single pass
            density_on || break

            # DO 30: recompute survivor QMD (D10N); iterate if far from D10 (morts.f:647-695)
            ttn = 0f0; sd2sqn = 0f0; sumdr10n = 0f0
            @inbounds for i in 1:n
                d = t.dbh[i]; d < dbhsdi && continue
                pr = t.tpa[i] - killed[i]
                bark = on_bratio(Int(t.species[i]), d, t.height[i])
                g = (t.diam_growth[i] / bark) * (fint / 10f0)
                sd2sqn += pr * (d * d + 2f0 * d * g + g * g)
                zeide && (sumdr10n += pr * on_powf(d + g, 1.605f0))
                ttn += pr
            end
            ttn == 0f0 && break
            d10n = zeide ? on_powf(sumdr10n / ttn, 1f0 / 1.605f0) : sqrt(sd2sqn / ttn)
            (abs(d10cur - d10n) <= 0.2f0) && break
            d10n <= dia0 && (ipath = 0; break)
            d10cur = d10n                             # GO TO 10
        end
    end

    # --- SIZCAP size-cap mortality (morts.f:789-803) ---
    let sc = ctl.sp_size_cap
        @inbounds for i in 1:n
            sp = Int(t.species[i])
            (sc[sp, 1] >= 999f0 || trunc(Int, sc[sp, 3]) == 1) && continue
            d = t.dbh[i]
            bark = on_bratio(sp, d, t.height[i])
            g = (t.diam_growth[i] / bark) * (fint / 10f0)
            if (d + g) >= sc[sp, 1]
                kc = min(t.tpa[i] * sc[sp, 2] * fint / 10f0, t.tpa[i])
                killed[i] < kc && (killed[i] = kc)
            end
        end
    end

    # --- BAMAX enforcement (morts.f:808-873) ---
    if bamax > 0f0
        @inbounds for _ in 1:100
            banew = 0f0; badead = 0f0
            for i in 1:n
                d = t.dbh[i]
                bark = on_bratio(Int(t.species[i]), d, t.height[i])
                g = (t.diam_growth[i] / bark) * (fint / 10f0)
                de2 = 0.0054542f0 * (d + g)^2
                banew  += de2 * (t.tpa[i] - killed[i])
                badead += de2 * killed[i]
            end
            (banew - bamax) > 1f0 || break
            badead <= 0f0 && break
            adjfac = (banew - bamax) / badead
            for i in 1:n
                killed[i] = min(t.tpa[i], killed[i] * (1f0 + adjfac))
            end
        end
    end

    # TPAMRT = surviving TPA (morts.f:874) — locked before FIXMORT
    surv = 0f0
    @inbounds for i in 1:n; surv += t.tpa[i] - killed[i]; end
    dens.tpa_mort = surv

    # FIXMORT (morts.f:880) — inert without the keyword
    apply_fixmort!(s, killed, n, fint)

    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end

"""
    _on_tn10!(dens, t, dia0, d10, constv, sdislp, pmsdil, pmsdiu, t85d0, t55d0, t85d10) -> tn10

ON density-target TN10 (morts.f label 210..270): the number of trees remaining after density
mortality. When the stand is over the 85% envelope at DIA0, kill straight to T85D10; between 55%
and 85%, solve the linear ln-ln self-thinning line (persisted in `dens`) and evaluate it at D10.
"""
function _on_tn10!(dens::Density, t::Float32, dia0::Float32, d10::Float32, constv::Float32,
                   sdislp::Float32, pmsdil::Float32, pmsdiu::Float32,
                   t85d0::Float32, t55d0::Float32, t85d10::Float32, t55d10::Float32)::Float32
    # T > T85D0 -> kill to 85% level at D10 (morts.f:462-464)
    t > t85d0 && return t85d10
    if t > t55d0
        # close to the 85% line at DIA0 (morts.f:475-478)
        abs(t85d0 - t) <= 5f0 && return t85d10
        slp, cept = _on_line_iter(t, dia0, constv, sdislp, pmsdil, pmsdiu, 1)
    else
        # T <= T55D0 (morts.f:240): hold constant if also below 55% at D10; else solve directly.
        t <= t55d10 && return t
        slp, cept = _on_line_iter(t, dia0, constv, sdislp, pmsdil, pmsdiu, 2)
    end
    dens.mort_slope == 0f0 && (dens.mort_slope = slp)
    dens.mort_intercept == 0f0 && (dens.mort_intercept = cept)
    tn10 = on_expf(dens.mort_intercept + dens.mort_slope * on_logf(d10 * ON_INtoCM)) * ON_ACRtoHA
    tn10 >= t85d10 && (tn10 = t85d10)
    return tn10
end

# morts.f label 220-230 linear-fn (ln-ln) solve for (SLP,CEPT). ipath 1 = iterate TREEIT; 2 = direct.
function _on_line_iter(t::Float32, dia0::Float32, constv::Float32, sdislp::Float32,
                       pmsdil::Float32, pmsdiu::Float32, ipath::Int)
    line(tem) = begin
        temp = tem * ON_HAtoACR
        d55m = (on_logf(temp) - on_logf(pmsdil * constv)) / sdislp
        t55m = on_logf(temp)
        d85m = d55m * 1.25f0
        local slp::Float32
        while true
            d85m > 5f0 && (d85m = 5f0)
            d85m < 0.125f0 && (d85m = 0.125f0)
            t85m = on_logf(constv * on_powf(on_expf(d85m), sdislp) * pmsdiu)
            slp = (t85m - t55m) / (d85m - d55m)
            (slp > -0.5f0 && d85m < 5f0) ? (d85m += 0.1f0) : break
        end
        (slp, t55m - slp * d55m)
    end
    if ipath == 1
        treeit = t + 0.1f0 * t
        slp = 0f0; cept = 0f0
        for _ in 1:100
            slp, cept = line(treeit)
            tprime = cept + slp * on_logf(dia0 * ON_INtoCM)
            diff = t - on_expf(tprime) * ON_ACRtoHA
            (-5f0 <= diff <= 5f0) && break
            treeit += 0.5f0 * diff
        end
        return slp, cept
    else
        return line(t)
    end
end
