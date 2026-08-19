# =============================================================================
# lpmpb_lpopdy.jl — LPMPB LPOPDY population-dynamics path (mpbdrv.f/mpbmod.f + helpers).
#
# The alternate MPB mortality path: the POPDYN keyword (mpbin.f opt 26) routes MPBCUP → MPBDRV
# → GARBEL classify → SURFCE → MPBMOD (an 808-line DETERMINISTIC beetle-epidemic simulation) →
# per-class survival → WK2 mortality. The default (Cole rate-of-loss, `mpb_apply!`) is separate.
#
# Every routine below is a faithful port of the Fortran, validated bit-exact (or, for BETIN, to
# 1.56e-8 — sufficient because the epidemic drives class survival to ~1e-11, far below the
# mortality-cap threshold) via g16 single-.o-swap dump-replay against FVSie_lpmpb. The array-based
# signatures let the same code back both the golden-replay regression test (test/fixtures/lpopdy/)
# and the live `mpb_lpopdy!` seam. See scratchpad/lpmpb/LPOPDY_HANDOFF.md for the full derivation.
# =============================================================================

# --- LPMPB constants (mpbint.f defaults) ---
const LPO_EPS    = 1.0e-6
const LPO_BMIN   = 1.0e-10
const LPO_NG     = 2
const LPO_INCRS  = 10
const LPO_IB     = 1
const LPO_MPMXYR = 10
const LPO_CE     = 1.0f0
const LPO_SQFTPA = 43560.0f0
const LPO_CF1    = 0.01f0
const LPO_CF2    = 0.01f0
const LPO_CF3    = 0.5f0
const LPO_CRITAD = 1.5f0
const LPO_AMP1   = 1200.0f0
const LPO_AMP2   = 600.0f0
const LPO_HS     = 1.0f0
const LPO_TAFAC  = 1.7f0
const LPO_TAMIN  = 1.7f0
const LPO_TAMAX  = 3.0f0
const LPO_STRBUG = 500.0f0
const LPO_STRP   = 0.95f0
const LPO_SEXRAT0 = 0.66f0                       # mpb block default (fort.785: 3F28F5C3)
const LPO_EXCON  = 640.0f0
const LPO_DST    = Float32[3000.0, 500.0, 500.0]
const LPO_PGRX   = Float32[0.0, 0.7, 1.0, 1.3, 1.5]   # mpbmod.f resistance table
const LPO_TAY    = Float32[2.0, 2.1, 2.4, 2.7, 3.0]
const LPO_BETTER = Float32[1.0, 4.0]             # GARBEL attr weights (DBH, phloem)

# -----------------------------------------------------------------------------
# BETIN incomplete-beta package (betin/forw/back/pqsml/mpbgam.f), Float64.
# -----------------------------------------------------------------------------
function _lpo_mpbgam(xx::Float64)::Float64
    zz = xx
    local dlgam::Float64
    if xx > 1.0e10
        dlgam = xx < 1.0e70 ? zz*(log(zz)-1.0) : 1.0e75
    elseif xx <= 1.0e-9
        dlgam = -1.0e75
    else
        term = 1.0
        while zz <= 18.0
            term *= zz; zz += 1.0
        end
        rz2 = 1.0/zz^2
        dlgam = (zz-0.5)*log(zz) - zz + Float64(0.9189385332046727f0) - log(term) +
            (1.0/zz)*(0.8333333333333333e-1 - (rz2*(0.2777777777777777e-2 +
            (rz2*(0.7936507936507936e-3 - (rz2*(0.5952380952380952e-3)))))))
    end
    return exp(dlgam)
end

function _lpo_pqsml(x::Float64, p::Float64, q::Float64)::Float64
    u = x^p; s = u/p; xk = 0.0
    while true
        xk += 1.0
        u = (xk-q)*x*u/xk
        v = u/(xk+p)
        s += v
        abs(v) > abs(LPO_EPS*s) || break
    end
    return s*_lpo_mpbgam(p+q)/(_lpo_mpbgam(p)*_lpo_mpbgam(q))
end

function _lpo_forw(x::Float64, p::Float64, q::Float64, x1::Float64, x2::Float64, n::Int)::Float64
    n < 1 && return x1
    n == 1 && return x2
    n1 = n-1; t1 = x1; pq1 = p+q-1.0; xa1 = 1.0-x; forr = x2
    for i in 1:n1
        xi = Float64(i); temp = forr
        forr = forr + (forr-t1)*(pq1+xi)*xa1/(xi+q)
        t1 = temp
    end
    return forr
end

function _lpo_back(x::Float64, p::Float64, q::Float64, z::Float64, n::Int)::Float64
    bk = z
    n == 0 && return bk
    q1 = q-1.0; ap = z; fn = Float64(n); xnu = Float64(n+n+6)
    while true
        xn = xnu; r = 0.0; bk = z
        while true
            xn -= 1.0
            t = p + xn
            u = (t+q1)*x
            r = u/(u+t-t*r)
            xn <= fn && (bk = r*bk)
            bk < LPO_BMIN && return bk
            xn > 1.0 || break
        end
        abs(bk-ap) < abs(LPO_EPS*bk) && return bk
        ap = bk; xnu += 5.0
    end
end

"BETIN (betin.f): regularized incomplete beta I_x(a,b), DOUBLE precision."
function mpb_betin(a::Float64, b::Float64, x::Float64)::Float64
    (a > 0.0 && b > 0.0 && x >= 0.0 && x <= 1.0) || return 0.0
    ia = trunc(Int, a); p = a - ia
    if p == 0.0; p = 1.0; ia -= 1; end
    (x == 0.0 || x == 1.0) && return x
    local bet::Float64
    if x <= 0.5
        q = b; m = trunc(Int, q); q1 = q - m
        if q1 == 0.0; q1 = 1.0; m -= 1; end
        x1 = _lpo_pqsml(x, p, q1)
        x2 = m > 0 ? _lpo_pqsml(x, p, q1+1.0) : 0.0
        bet = _lpo_back(x, p, q, _lpo_forw(x, p, q1, x1, x2, m), ia)
    else
        y = 1.0 - x; q = p; p = b; m = trunc(Int, p); p1 = p - m
        if p1 == 0.0; m -= 1; p1 = 1.0; end
        x1 = _lpo_back(y, p1, q,      _lpo_pqsml(y, p1, q),      m)
        x2 = _lpo_back(y, p1, q+1.0,  _lpo_pqsml(y, p1, q+1.0),  m)
        bet = 1.0 - _lpo_forw(y, p, q, x1, x2, ia)
    end
    bet < LPO_BMIN && (bet = 0.0)
    return bet
end

# -----------------------------------------------------------------------------
# Phloem model (mpbdrv.f:106-116) — XPT per LP tree.
# -----------------------------------------------------------------------------
"XPT phloem thickness proxy: exp(-3.17152 + .12591 ln(BAI5) + .50932 ln(DBH) - .0077 HT), 0 if BAI5≤1e-4."
@inline function mpb_phloem_xpt(dbh::Float32, dg::Float32, ht::Float32)::Float32
    dds5 = (2.0f0*dbh*dg + dg*dg)/2.0f0
    bai5 = dds5*0.7853982f0
    bai5 <= 0.0001f0 && return 0.0f0
    return fexp(-3.17152f0 + 0.12591f0*flog(bai5) + 0.50932f0*flog(dbh) - 0.0077f0*ht)
end

# -----------------------------------------------------------------------------
# GARBEL classifier (garbel.f + grpsum/grclas) — reduce NRECS trees → NACLAS classes.
# Returns (mp1, mp2, ipt, clsprob, avgdbh, avgxpt) over `nclas` classes.
# -----------------------------------------------------------------------------
function mpb_garbel(dbh::Vector{Float32}, xpt::Vector{Float32}, prob::Vector{Float32}, nclass::Int)
    nrecs = length(dbh)
    naclas = min(nrecs, nclass)
    pn1 = 0.5f0
    ipt = collect(1:nrecs)
    Pw = zeros(Float32, nrecs)
    # GRPSUM: standardized weighted score accumulation (DBH wt 1, XPT wt 4).
    function grpsum!(P, atr, wt)
        xn = Float32(nrecs); ave = 0.0f0; stdv = 0.0f0
        for ii in 1:nrecs
            i = ipt[ii]; ave += atr[i]; stdv += atr[i]^2
        end
        stdv = (stdv - ave*ave/xn)/xn
        stdv = stdv > 1.0f-9 ? sqrt(stdv) : 1.0f0
        ave = ave/xn
        for ii in 1:nrecs
            i = ipt[ii]; P[i] += wt*(atr[i]-ave)/stdv
        end
    end
    grpsum!(Pw, dbh, LPO_BETTER[1])
    grpsum!(Pw, xpt, LPO_BETTER[2])
    sort!(ipt, by = i -> -Pw[i])                 # RDPSRT descending
    # method 1: NCL1 max-diff class boundaries (Fortran REAL→INT truncates)
    ncl1 = trunc(Int, naclas*pn1 + 0.5f0); ncl1 < 1 && (ncl1 = 1)
    ncl2 = naclas - ncl1
    if !(nrecs - naclas > 5); ncl1 = naclas; ncl2 = 0; end
    clasd = zeros(Float32, ncl1)
    sz = max(naclas, ncl1) + ncl2 + 4
    isc1 = zeros(Int, sz); isc2 = zeros(Int, sz)
    nc1 = ncl1 - 1
    ipt1 = ipt[1]
    for I in 2:nrecs
        ipt2 = ipt[I]
        diffp = Pw[ipt1] - Pw[ipt2]
        jset = ncl1
        for J in 1:nc1
            if diffp <= clasd[J+1]; jset = J; break; end
            clasd[J] = clasd[J+1]; isc1[J] = isc1[J+1]
        end
        clasd[jset] = diffp; isc1[jset] = I
        ipt1 = ipt2
    end
    isc1[1] = nrecs + 1
    sort!(view(isc1, 1:ncl1))
    I1 = 1
    for J in 1:ncl1
        I2 = isc1[J] - 1
        isc1[J] = I2; isc2[J] = I2 - I1 + 1; I1 = I2 + 1
    end
    # method 2: split the ncl2 largest classes
    for K in 1:ncl2
        mx = 0; ibig = 0
        for J in 1:ncl1
            if isc2[J] > mx; mx = isc2[J]; ibig = J; end
        end
        ncl1 += 1
        mxs = floor(Int, mx/2.0f0 + 0.5f0)
        isc1[ncl1] = isc1[ibig]; isc1[ibig] -= mxs
        isc2[ncl1] = mxs;        isc2[ibig] -= mxs
    end
    sort!(view(isc1, 1:ncl1))
    I1 = 1
    mp1 = zeros(Int, ncl1); mp2 = zeros(Int, ncl1)
    for J in 1:ncl1
        I2 = isc1[J]; mp1[J] = I1; mp2[J] = I2; I1 = I2 + 1
    end
    # STEP7 + GRCLAS: class PROB sum + PROB-weighted avg DBH/XPT
    clsprob = zeros(Float32, ncl1); avgdbh = zeros(Float32, ncl1); avgxpt = zeros(Float32, ncl1)
    for c in 1:ncl1
        sp = 0.0f0
        for jj in mp1[c]:mp2[c]; sp += prob[ipt[jj]]; end
        clsprob[c] = sp
        ad = 0.0f0; ax = 0.0f0
        for jj in mp1[c]:mp2[c]
            i = ipt[jj]; ad += dbh[i]*prob[i]; ax += xpt[i]*prob[i]
        end
        avgdbh[c] = sp > 1.0f-30 ? ad/sp : 0.0f0
        avgxpt[c] = sp > 1.0f-30 ? ax/sp : 0.0f0
    end
    return (mp1, mp2, ipt, clsprob, avgdbh, avgxpt)
end

"SURFLP (surflp.f): LP surface area from class-avg DBH."
@inline mpb_surflp(d::Float32) = d > 5.0f0 ? 8.835f0*d - 40.82f0 : d*0.672f0

# -----------------------------------------------------------------------------
# MPBMOD year-loop (mpbmod.f) — DETERMINISTIC epidemic; LGO path, LAGG/LREP/LPS/LDC=F.
# Returns the per-class survivors TREES(I) (=CLASS(I,IMPROB) after the epidemic).
# -----------------------------------------------------------------------------
function mpb_mpbmod(trees0::Vector{Float32}, surf::Vector{Float32}, diam::Vector{Float32},
                    phloem::Vector{Float32}, ta::Float32, efelev::Float32, eflat::Float32)
    naclas = length(trees0)
    tprob = sum(trees0)
    sexdbh = Float32[0.918f0 - 0.0168f0*diam[i] for i in 1:naclas]
    efphlm = Float32[max(0.0f0, 16.67f0*phloem[i] - 0.667f0) for i in 1:naclas]
    exodus = Float32[1.0f0 - fexp(-LPO_CE*LPO_DST[ig]*LPO_DST[ig]/(LPO_EXCON*LPO_SQFTPA)) for ig in 1:LPO_NG]
    trees = copy(trees0)
    by = LPO_STRBUG; pp = LPO_STRP; sexrat = LPO_SEXRAT0
    mpbyr = 0; sadlpp = 0.0f0
    ninc = LPO_INCRS + 1
    agg = zeros(Float32, naclas, ninc); ad = zeros(Float32, ninc); amp = zeros(Float32, ninc)
    bold = zeros(Float32, LPO_NG)
    ec = Ref(0.0)                                # EMERG persistent C (the -fno-automatic static)
    emerg(byv, t) = (t == 0 ? (ec[] = 1.0/2.0^LPO_INCRS) : (ec[] = ec[]*(LPO_INCRS - t + 1)/t); byv*Float32(ec[]))
    # OS (total live surface) is computed ONCE at MPBMOD entry from the initial CLASS(I,IMPROB)
    # (mpbmod.f:236-243, before the annual loop) — CONSTANT across the epidemic years, NOT recomputed.
    os = 0.0f0
    for I in 1:naclas; os += surf[I]*trees[I]; end
    os += sadlpp
    while true
        mpbyr += 1
        q = 1.0f0 - pp
        geno = LPO_NG == 3 ? Float32[pp*pp, 2*pp*q, q*q] : Float32[pp, q]
        b3sum = zeros(Float32, LPO_NG); fill!(bold, 0.0f0)
        fill!(agg, 0.0f0); fill!(ad, 0.0f0); fill!(amp, 0.0f0)
        for INC in 1:ninc
            inc1 = INC - 1
            e1 = 0.0f0
            for I in 1:naclas; e1 += surf[I]*trees[I]; end
            e3 = 0.0f0; ef3 = 0.0f0; tagg = 0.0f0
            if INC > 1
                for KK in LPO_IB:inc1
                    amp[KK] = LPO_AMP1 - LPO_AMP2*ad[KK]*(1.0f0 - sexrat)
                    amp[KK] < 0.0f0 && (amp[KK] = 0.0f0)
                    for I in 1:naclas
                        tagg += agg[I,KK]; e3 += surf[I]*agg[I,KK]; ef3 += surf[I]*agg[I,KK]*amp[KK]
                    end
                end
            end
            rho1 = 0.0f0
            for I in 1:naclas; rho1 += trees[I]; end
            e2 = os - (e1 + e3)                   # SNOHST=0
            effs = e1 + e2 + ef3
            rho3 = tagg/LPO_SQFTPA
            rho2 = tprob - tagg - rho1; rho2 < 0.0f0 && (rho2 = 0.0f0); rho2 /= LPO_SQFTPA
            rho1 /= LPO_SQFTPA
            bnew = emerg(by, inc1)
            b1inc = 0.0f0; b3inc = 0.0f0
            for ig in 1:LPO_NG
                b0 = geno[ig]*bnew + bold[ig]
                exloss = b0*exodus[ig]
                b0 -= exloss
                b1 = b0*e1/effs; b2 = b0*e2/effs; b3 = b0*ef3/effs
                fm(cf,rho) = (fmx = -cf*2.0f0*sqrt(rho)*LPO_DST[ig]*LPO_DST[ig]/LPO_DST[1]; fmx > -80.0f0 ? fexp(fmx) : 0.0f0)
                b1 -= b1*fm(LPO_CF1,rho1); b2 -= b2*fm(LPO_CF2,rho2); b3 -= b3*fm(LPO_CF3,rho3)
                bold[ig] = b1 + b2
                b1inc += b1; b3inc += b3; b3sum[ig] += b3
            end
            pioden = Float64(b1inc/e1)
            xx = 1.0 - ccall((:exp, "libm.so.6"), Float64, (Float64,), -pioden)  # glibc DEXP
            for I in 1:naclas
                agg[I,INC] = 0.0f0
                dta = Float64(ta); dsmta = Float64(surf[I]) - dta + 1.0
                if xx != 0.0 && dsmta > 0.0
                    agg[I,INC] = Float32(Float64(trees[I])*mpb_betin(dta, dsmta, xx))
                end
                trees[I] -= agg[I,INC]
            end
            if ef3 > 0.0f0
                for KK in LPO_IB:inc1; ad[KK] += amp[KK]*b3inc/ef3; end
            end
        end
        by = 0.0f0; byt = 0.0f0
        for INC in 1:ninc
            xeg = -0.117f0*ad[INC] > -80.0f0 ? fexp(-0.117f0*ad[INC]) : 0.0f0
            eggs = 630.0f0*(1.0f0 - xeg)
            for I in 1:naclas
                psurv = 1.0f0 - fexp(-ad[INC]*0.04328f0)
                ad[INC] >= 2.595f0 && (psurv = psurv*6.812f0*fexp(-1.191f0*sqrt(ad[INC])))
                young = eggs*psurv*efelev*eflat*efphlm[I]
                if ad[INC] >= LPO_CRITAD
                    trkill = agg[I,INC]; skill = surf[I]*trkill
                    by += young*skill*LPO_HS*sexdbh[I]; byt += young*skill*LPO_HS
                else
                    trees[I] += agg[I,INC]; agg[I,INC] = 0.0f0
                end
            end
        end
        byt > 0.0f0 && (sexrat = by/byt)
        b3all = sum(b3sum)
        if b3all != 0.0f0
            pp = LPO_NG == 3 ? (b3sum[1] + 0.5f0*b3sum[2])/b3all : (b3sum[1]/b3all)^2
        end
        (by >= 1.0f0 && mpbyr < LPO_MPMXYR) || break
    end
    return trees
end

"EFELEV/EFLAT (mpbmod.f:150-162) from stand elevation (hundreds-ft) and forest latitude."
@inline function mpb_efelev(elev::Float32)::Float32
    e = elev <= 0.0f0 ? 63.0f0 : elev
    v = 2.62f0 - 2.70f-2*e
    v > 1.0f0 ? 1.0f0 : (v < 0.0f0 ? 0.0f0 : v)
end
@inline function mpb_eflat(forlat::Float32)::Float32
    v = 4.667f0 - 8.333f-2*forlat
    v < 0.0f0 ? 0.0f0 : (v > 1.5f0 ? 1.5f0 : v)
end

"PMSLP (pmslp.f): piecewise-linear interpolation with flat extrapolation."
function mpb_pmslp(xx::Float32, x::Vector{Float32}, y::Vector{Float32})::Float32
    n = length(x)
    for i in 1:(n-1)
        (xx < x[i] || xx > x[i+1]) && continue
        return y[i] + ((y[i+1]-y[i])/(x[i+1]-x[i]))*(xx-x[i])
    end
    return y[n]
end

# -----------------------------------------------------------------------------
# mpb_lpopdy! — the live LPOPDY seam (MPBDRV): phloem → GARBEL → SURFCE → MPBMOD → WK2 mortality.
# Called from mpb_apply! when m.lpopdy and an outbreak is due this cycle. `ta` is the aggregation
# threshold (MPGR resistance); `lpidx` are the LP tree record indices; mutates t.tpa.
# -----------------------------------------------------------------------------
function mpb_lpopdy!(t, old_tpa::AbstractVector{Float32}, lpidx::Vector{Int},
                     ta::Float32, elev::Float32, forlat::Float32)
    nlp = length(lpidx)
    nlp == 0 && return nothing
    dbh = Float32[t.dbh[i] for i in lpidx]
    ht  = Float32[t.height[i] for i in lpidx]
    dg  = Float32[t.diam_growth[i] for i in lpidx]
    prob = Float32[old_tpa[i] for i in lpidx]
    xpt = Float32[mpb_phloem_xpt(dbh[k], dg[k], ht[k]) for k in 1:nlp]
    mp1, mp2, ipt, clsprob, avgdbh, avgxpt = mpb_garbel(dbh, xpt, prob, 10)
    surf = Float32[mpb_surflp(avgdbh[c]) for c in 1:length(mp1)]
    efelev = mpb_efelev(elev); eflat = mpb_eflat(forlat)
    surv = mpb_mpbmod(copy(clsprob), surf, avgdbh, avgxpt, ta, efelev, eflat)
    # per-class SURVIV = survivors/pre; per tree WK2 = max(current, PROB·(1-SURVIV)); cap PROB-1e-6.
    for c in 1:length(mp1)
        surviv = clsprob[c] > 0.0f0 ? surv[c]/clsprob[c] : 0.0f0
        dead = 1.0f0 - surviv
        for jj in mp1[c]:mp2[c]
            j = lpidx[ipt[jj]]                     # tree record index
            prob_j = old_tpa[j]
            wk2 = prob_j * dead
            cur = old_tpa[j] - t.tpa[j]            # current period mortality already applied
            wk2 < cur && (wk2 = cur)               # WK2=max(WK2, current)
            (prob_j - wk2 < 1.0f-6) && (wk2 = prob_j - 1.0f-6)
            t.tpa[j] = prob_j - wk2
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# MPGR (mpgr.f) — periodic growth ratio → aggregation threshold TA (resistance).
# PGR = Σ(FDG/ODG·PROB)/ΣPROB over LP trees with PCT ≥ PCTCO(65); if NPGR≤2 → PGR=0.9.
# FDG = current (projected) DG, ODG = prior/inventory DG rescaled to a yr basis. Then
# TA = clamp(PMSLP(PGR, PGRX, TAY), TAMIN, TAMAX). Bark = BRATIO(7) inside-bark factor.
#
# ⚠ NOT yet wired into mpb_apply! (TA is hardcoded to the captured MPGR value 2.099609 for the
# lp_popdy stand). Wiring it live needs: (1) a MpbState scratch to save each LP tree's MEASURED
# DG (t.diam_growth at input) BEFORE grow_cycle!'s diameter_growth! (simulate.jl:600) overwrites
# it with the projected DG — the MPSVDG-equivalent (mpgr.f:120 ENTRY MPSVDG, saved into XPT);
# (2) the per-tree PCT (BA-percentile-in-larger-trees) for the PCTCO≥65 dominance gate; (3) NPYR
# (last cycle's measurement period) + SCALE=YR/NPYR. The lp_popdy MPGR result is TA=2.0996 (PGR≈0.697).
# -----------------------------------------------------------------------------
const LPO_PCTCO = 65.0f0

"MPGR periodic growth ratio (mpgr.f:60-98): PGR over dominant LP trees; fdg/odg/prob/pct per LP tree."
function mpb_mpgr(fdg::Vector{Float32}, odg::Vector{Float32}, prob::Vector{Float32},
                  pct::Vector{Float32}, bark::Vector{Float32}, dbh::Vector{Float32};
                  scale::Float32 = 1.0f0)::Float32
    sump = 0.0f0; avrat = 0.0f0; npgr = 0
    @inbounds for i in eachindex(fdg)
        pct[i] < LPO_PCTCO && continue
        npgr += 1
        prb = prob[i]; sump += prb
        d2 = dbh[i]*bark[i]                       # inside-bark DBH
        og = odg[i]
        odds = og*(2.0f0*d2 - og)*scale           # old DDS on yr basis
        x = d2*d2 - odds
        ogy = x > 1.0f-6 ? d2 - sqrt(x) : 0.0f0   # ODG on yr basis
        ogy > 0.0f0 && (avrat += fdg[i]/ogy*prb)
    end
    npgr <= 2 && return 0.9f0                      # mpgr.f:53 (≤2 selected LPP)
    return avrat/sump
end

"TA aggregation threshold from PGR (mpbmod.f:228-231): PMSLP(PGR,PGRX,TAY) clamped [TAMIN,TAMAX]."
@inline function mpb_lpopdy_ta(pgr::Float32)::Float32
    ta = mpb_pmslp(pgr, LPO_PGRX, LPO_TAY)
    ta < LPO_TAMIN ? LPO_TAMIN : (ta > LPO_TAMAX ? LPO_TAMAX : ta)
end
