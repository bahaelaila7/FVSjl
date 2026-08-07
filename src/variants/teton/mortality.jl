# =============================================================================
# mortality.jl (teton) — TT mortality (tt/morts.f).
#
# Stand-level (= KT structure): grown-stand DQ10, CONST=SDIMAX/0.02483133, TMD10=CONST·D10^−1.605,
# T85D10=TMD10·PMSDIU, T55D10=TMD10·PMSDIL. SDI self-thinning target TN10 → RN=1−(1−(T−TN10)/T)^(1/FINT).
#   T ≤ T55D10  → TN10=T, RN=0 (below the SDI limit — the emt01 regime; PMSC background dominates).
#   T > T85D10  → TN10=T85D10 (kill to the 85% line).
#   middle (55–85%): iterative linear-fn fit (tt/morts.f label 220) — DEFERRED (emt01 is below-SDI).
# Per-tree, ORIGINAL EM species (1-3,7-10,18): RI=0.5/(1+exp(PMSC+PMD·D+PMDSQ·D²)); RIP=RN if SDI limiting
#   (T>TEM and RN>0) else RI; WKI=P·(1−(1−RIP)^FINT). ADDED species (4-6,11-17,19): KT density Hamilton
#   (deferred — not in emt01). Reuses the shared self-thinning RDPSRT + snag booking.
# =============================================================================

const TT_PMSC  = Float32[6.5112, 6.5112, 7.2985, 5.1677, 9.6943, 5.1677, 5.9617, 9.6943, 5.1677, 0.2118, 5.1677, 5.1677, 5.5877, 5.1677, 5.9617, 5.9617, 5.1677, 5.9617]
const TT_PMD   = Float32[-0.0052485, -0.0052485, -0.0129121, -0.0077681, -0.0127328, -0.0077681, -0.0340128, -0.0127328, -0.0077681, 0.0, -0.0077681, -0.0077681, -0.005348, -0.0077681, -0.0052485, -0.0340128, -0.0077681, -0.0052485]
const TT_PMDSQ = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]  # TT has no PMDSQ

@inline _tt_mort_default(sp::Int) = sp != 10   # only PP(10) has a special morts form (not in ttt01)

# tt/ttmrt.f VARADJ — species shade tolerance (1.0 = most intolerant), the EFFTR scalar.
const TT_VARADJ = Float32[0.80, 0.70, 0.55, 0.70, 0.50, 1.00, 0.90, 0.50, 0.60, 0.85,
                          0.70, 0.70, 0.70, 1.00, 0.90, 1.10, 0.75, 0.90]

# tt/ttmrt.f — redistribute TOKILL across the treelist by PERCENTILE (PCT = t.crown_ratio) + species tolerance
# (EFFTR), via a geometric progression converging to kill exactly TOKILL. OVERWRITES `killed`. If tokill==0
# (below-SDI/background), TOKILL = Σ current killed (redistribute the background). NC/OH get the extra CRI factor.
function _tt_ttmrt!(killed::AbstractVector{Float32}, tokill::Float32, s::StandState, n::Int)
    t = s.trees
    efftr = Vector{Float32}(undef, n); temwk2 = zeros(Float32, n)
    pass1 = 0f0
    @inbounds for i in 1:n
        pct = Float32(t.crown_ratio[i])                      # PCT (BA percentile)
        peff = 0.84525f0 - 0.01074f0 * pct + 0.0000002f0 * pct * pct * pct
        peff > 1f0 && (peff = 1f0); peff < 0.01f0 && (peff = 0.01f0)
        sp = Int(t.species[i]); cri = Float32(t.crown_pct[i])
        efftr[i] = (sp == 15 || sp == 18) ? peff * ((100f0 - cri) / 100f0) * TT_VARADJ[sp] * 0.01f0 :
                                            peff * TT_VARADJ[sp] * 0.01f0
        pass1 += t.tpa[i] * efftr[i]
    end
    if tokill == 0f0
        @inbounds for i in 1:n; tokill += killed[i]; end
    end
    (pass1 <= 0f0 || tokill <= 0f0) && return                # nothing to distribute
    fill!(killed, 0f0)
    temkil = tokill; short = 0f0
    npass = trunc(Int, tokill / pass1) + 1
    jpass = 0
    while true                                               # ttmrt.f label 100
        jpass += 1
        jpass > 1 && (temkil = short)
        iswtch = 0; adjust = 1f0
        while true                                           # ttmrt.f label 105 (NPASS convergence)
            temsum = 0f0
            @inbounds for i in 1:n
                tpalft = t.tpa[i] - killed[i]
                if tpalft > 0f0
                    temwk2[i] = -tpalft * ((1f0 - efftr[i])^npass - 1f0); temsum += temwk2[i]
                end
            end
            minstp = npass > 50 ? 5 : (npass > 20 ? 2 : 1)
            temsum <= 0f0 && (adjust = 1f0; break)
            adjust = temkil / temsum
            if adjust < 0.8f0
                iswtch == 2 && break
                npass -= max(minstp, trunc(Int, (temsum - temkil) / pass1)); iswtch = 1
                npass <= 0 && break
            elseif adjust > 1.2f0
                iswtch == 1 && break
                npass += max(minstp, trunc(Int, (temkil - temsum) / pass1)); iswtch = 2
            else
                break
            end
        end
        short = 0f0                                          # ttmrt.f label 110 (scale + cap)
        @inbounds for i in 1:n
            tpalft = t.tpa[i] - killed[i]
            tpalft < 0.00001f0 && continue
            xkill = temwk2[i] * adjust
            if (t.tpa[i] - killed[i] - xkill) <= 0.00001f0
                temwk2[i] = t.tpa[i] - killed[i]
                short += xkill - t.tpa[i] + killed[i]; pass1 -= efftr[i]
            else
                temwk2[i] = xkill
            end
            killed[i] += temwk2[i]
        end
        (short > 0f0 && pass1 > 0f0) || break
        npass = trunc(Int, short / pass1) + 1
    end
    return
end

# tt/morts.f MORCON — PP CI-variant REIN (potential mort rate) + GMULT (DG size multiplier) by size class
# IP (1=D>5, 2=D≤5). IPDG=12→POT=0.80, IPDG2=41→POT=2.25 are fixed ⇒ these are constants.
# REIN(1)=(1−(0.80/20+1)^−1.605)/0.06821, REIN(2)=(1−(2.25+1)^−1.605)/0.86610; GMULT(1)=0.90/0.80, GMULT(2)=2.50/2.25.
const TT_PP_REIN  = (Float32(0.89443f0), Float32(0.98042f0))    # (IP=1, IP=2)
const TT_PP_GMULT = (Float32(1.125f0),   Float32(1.111111f0))

# tt/morts.f label-220 iterative linear-fn fit between the 55%/85% SDI lines → TN10 (target tree count at
# D10). IPATH2=true (came from T≤T55D0, T>T55D10) computes the line once (TEM=T) then goes to 230; IPATH2=false
# (55%<T≤85% at DIA0) Newton-iterates TREEIT (≤100) so exp(CEPT+SLP·ln(DIA0)) ≈ T. TN10 = exp(CEPT+SLP·ln(D10)),
# capped at T85D10.
function _tt_tn10_iter(tt::Float32, dia0::Float32, d10::Float32, const_::Float32, pmsdil::Float32,
                       pmsdiu::Float32, t85d10::Float32, t55d0::Float32, ipath2::Bool)::Float32
    treeit = tt + 0.1f0 * tt
    slp = 0f0; cept = 0f0; knt = 1
    while true
        tem = ipath2 ? tt : treeit
        d55m = (log(tem) - log(pmsdil * const_)) / (-1.605f0)
        t55m = log(tem)
        d85m = d55m * 1.25f0
        while true                                    # tt/morts.f label 221: bump D85M until SLP ≤ −0.5
            d85m > 5f0 && (d85m = 5f0); d85m < 0.125f0 && (d85m = 0.125f0)
            t85m = log(const_ * (exp(d85m)^(-1.605f0)) * pmsdiu)
            slp = (t85m - t55m) / (d85m - d55m)
            (slp > -0.5f0 && d85m < 5f0) ? (d85m += 0.1f0) : break
        end
        cept = t55m - slp * d55m
        (ipath2 || tt <= t55d0) && break              # GOTO 230 (no Newton for the IPATH=2 path)
        tprime = cept + slp * log(dia0)
        diff = tt - exp(tprime)
        (diff <= 5f0 && diff >= -5f0) && break
        treeit += 0.5f0 * diff; knt += 1
        knt > 100 && break
    end
    tn10 = exp(cept + slp * log(d10))
    tn10 >= t85d10 && (tn10 = t85d10)
    return tn10
end

function mortality!(s::StandState, ::Teton; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    ba = p.basal_area
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # grown-stand sums (morts.f): T (total tpa), DQ10/DQ0, AVED (BA-weighted mean DBH).
    # ★ #148 (Zeide-QMD fix): TT is a Zeide-SDI variant (tt/grinit.f LZEIDE), so the self-thin diameter is
    # Reineke's DR10=(Σ p·(D+G)^1.605 / T)^(1/1.605), NOT the quadratic mean (tt/morts.f LZEIDE path) — same as UT
    # #147 (6e57347). QMD over-states D10 on dense sub-1" cohorts ⇒ TMD10 uncapped ⇒ TN10 low ⇒ RN self-thin OVER-KILL.
    tt = 0f0; sumdr10 = 0f0; sumdr0 = 0f0; dsum = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = tt_bratio(sp, d)
        g = t.diam_growth[i] / bark
        sumdr10 += pr * fpow(d + g, 1.605f0); sumdr0 += pr * fpow(d, 1.605f0); tt += pr; dsum += d * pr
    end
    tt < 1f-6 && return s
    dq10 = fpow(sumdr10 / tt, 1f0 / 1.605f0)   # Reineke DR10 (Zeide self-thin diameter; tt/morts.f:267 D10=DR10)
    dq0  = fpow(sumdr0 / tt, 1f0 / 1.605f0)     # DR0 = pre-growth Reineke diameter
    aved = dsum / tt
    # DIA0<0.3 reset (morts.f 374-376)
    if dq0 < 0.3f0; dq10 = 0.3f0 + dq10 - dq0; dq0 = 0.3f0; end
    # SDI self-thinning boundary (morts.f 455-485)
    sdimax = stand_sdimax(s)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    const_ = sdimax / 0.02483133f0
    tmd10 = const_ * dq10^(-1.605f0); tmd10 > 35000f0 && (tmd10 = 35000f0)
    tmd0  = const_ * dq0^(-1.605f0);  tmd0  > 35000f0 && (tmd0  = 35000f0)
    t85d10 = tmd10 * pmsdiu; t55d10 = tmd10 * pmsdil
    t85d0  = tmd0  * pmsdiu; t55d0  = tmd0  * pmsdil
    # TN10 target (morts.f 200-271): the SDI mature-stand-boundary self-thinning.
    local tn10::Float32
    if tt > t85d0
        tn10 = t85d10                                              # kill to the 85% line
    elseif tt > t55d0
        tn10 = abs(t85d0 - tt) <= 5f0 ? t85d10 :
               _tt_tn10_iter(tt, dq0, dq10, const_, pmsdil, pmsdiu, t85d10, t55d0, false)
    elseif tt <= t55d10
        tn10 = tt                                                  # below 55% at both — hold (RN=0)
    else
        tn10 = _tt_tn10_iter(tt, dq0, dq10, const_, pmsdil, pmsdiu, t85d10, t55d0, true)   # IPATH=2
    end
    tn10 > tt && (tn10 = tt); tn10 < 0.1f0 && (tn10 = 0f0)
    rn = 1f0 - (1f0 - (tt - tn10) / tt)^(1f0 / fint)
    tem = const_ * dq10^(-1.605f0) * pmsdil     # SDI threshold (morts.f 641)
    # PP CI-variant stand projection (tt/morts.f 273-293): BA forward 10y assuming BA/BAMAX of the BA
    # increment is lost to mortality → an annual TPA-mortality rate RZ. BAMAX defaults from weighted SDImax
    # (sdical.f:204 BAMAX = SDImax·0.5454154·PMSDIU) when not user-set. Only PP (CASE 10) consumes RZ/BAMAX.
    bamax = sdimax * 0.5454154f0 * pmsdiu        # LBAMAX=false default (no user BAMAX keyword in ttpp)
    deltba = 0.005454154f0 * dq10 * dq10 * tt - ba
    ba10 = bamax > 0f0 ? ba + ((bamax - ba) / bamax) * deltba : ba
    tb = ba10 / (0.005454154f0 * dq10 * dq10)
    ttb = (tt - tb) / tt; ttb > 0.9999f0 && (ttb = 0.9999f0)
    rz = 1f0 - (1f0 - ttb)^0.1f0
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]
        if _tt_mort_default(sp)
            ri = 0.5f0 * (1f0 / (1f0 + exp(TT_PMSC[sp] + TT_PMD[sp] * d + TT_PMDSQ[sp] * d * d)))
            rip = rn
            (tt <= tem || rn <= 0f0) && (rip = ri)              # background when SDI not yet limiting
            rip > 1f0 && (rip = 1f0)
            wki = pr * (1f0 - (1f0 - rip)^fint)
            wki > pr && (wki = pr)
            sdimax < 5f0 && (wki = pr)
            killed[i] = wki
        else
            # tt/morts.f CASE(10) — PP from the CI variant. RIP logistic → REIN(IP) potential rate;
            # RIPP blends BA·RZ with a BAMAX-approach term. WK1=prior-cycle DG (dg_prev, 0 at cycle 1 ⇒
            # the ICYC==1 override G=DG/(bark·10) fires for every DG>0.5 tree, matching live).
            dm = d <= 0.5f0 ? 0.5f0 : d
            bark = tt_bratio(10, d)
            reldbh = d / aved
            wk1 = t.dg_prev[i]                      # prior applied DG (inside-bark); 0 at cycle 1
            dgt = wk1 / fint                          # OLDFNT = FINT for the uniform-cycle case
            if dm <= 1f0 && dgt < 0.05f0
                dgt = 0.05f0
            elseif dm > 1f0 && dm <= 5f0 && dgt < 0.05f0
                dgt = 0.05f0 * (5f0 - dm) / 4f0
            end
            g = wk1 / (bark * fint)
            (wk1 / fint < dgt) && (g = dgt / bark)
            dgcur = t.diam_growth[i]                  # current applied DG (inside-bark)
            (wk1 == 0f0 && dgcur > 0.5f0) && (g = dgcur / (bark * 10f0))   # ICYC==1/WK1==0 override
            ip = dm <= 5f0 ? 2 : 1
            g *= TT_PP_GMULT[ip]
            rip = 2.76253f0 + 0.222310f0 * sqrt(dm) - 0.0460508f0 * sqrt(ba) + 11.2007f0 * g -
                  0.554421f0 / dm + TT_PMSC[10] + 0.246301f0 * reldbh + 6.07129f0 * g / dm
            rip = rip > 88.5f0 ? 88.5f0 : (rip < -88.5f0 ? -88.5f0 : rip)
            rip = 1f0 / (1f0 + exp(rip))
            rip *= TT_PP_REIN[ip]                     # POTENT = REIN(IP)
            ripp = ba * rz
            ba <= bamax && (ripp += (bamax - ba) * rip)
            ripp /= bamax
            ripp < rip && (ripp = rip)
            ripp > 1f0 && (ripp = 1f0)
            wki = pr * (1f0 - (1f0 - ripp)^fint)      # X=1 (no MORTMULT); establishment "best-tree" deferred
            wki > pr && (wki = pr)
            sdimax < 5f0 && (wki = pr)
            killed[i] = wki
        end
    end
    # tt/morts.f:682 — REDISTRIBUTE the mortality by percentile+EFFTR (TTMRT). When self-thinning (default trees
    # have rip==rn, i.e. tt>tem & rn>0) TOKILL = T−TN10; else TOKILL=0 ⇒ TTMRT redistributes the background total.
    # Called whenever TN10≥0.1 (morts.f). Overwrites `killed` with the percentile allocation of the same total.
    if tn10 >= 0.1f0
        self_thin = tt > tem && rn > 0f0
        tokill = self_thin ? (tt - tn10) : 0f0
        tokill < 0f0 && (tokill = 0f0)
        _tt_ttmrt!(killed, tokill, s, n)
    end
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
