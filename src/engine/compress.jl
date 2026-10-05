# =============================================================================
# compress.jl — COMPRESS (act 250): reduce the tree list to NCLAS representative
# records by principal-component-score clustering, then merge each class into one
# record (PROB-conserving). FAITHFUL transliteration of base/comprs.f + base/comcup.f
# (the act-250 driver): the 1966 IBM-SSP symmetric eigensolver (`_ibm_eigen`, a direct
# port of base/eigen.f) and the Scowen-1965 Quickersort sorts (`rdpsrt!`/`iqrsrt!`,
# quickersort.jl) are ported exactly, so the PC scores AND the discrete class partition
# are bit-reproducible against Fortran — not just the conserved aggregate. (Earlier this
# used LinearAlgebra.eigen + library sorts, whose differing eigenvector basis / tie order
# flipped the partition; that was the accepted COMPRESS divergence, now closed.) comcup.f
# calls this with WK2 (the mortality weight) zeroed, so the per-class weight is just PROB.
#
# Classification variables (comprs.f:161): HT, ICR, IMC, ln(DBH), DG (NRANK=5).
# =============================================================================

const _CMP_RNGMIN = 0.00001
const _CMP_ALGTOL = 0.066
# minimum standard deviations per classification variable (comprs.f:171-188)
const _CMP_STDFLOOR = (1.0, 0.1e-3, 0.1e-3, 5e-3, 0.02)

# =============================================================================
# IBM-SSP EIGEN (1966) — faithful port of base/eigen.f (the routine comprs.f calls).
#
# Real-symmetric eigensolver by cyclic Jacobi diagonalization with threshold sweeps,
# operating on COLUMN-PACKED upper-triangular storage (A[i+(j²−j)/2] = M[i,j], i≤j),
# all in Float64. On return the eigenvalues sit on the diagonal of `a` in DESCENDING
# order and `r` holds the eigenvectors column-wise in the SAME order (column c =
# r[(c-1)n+1 : c·n]). This replaces `LinearAlgebra.eigen`, whose (valid) eigenvector
# basis differs from FVS at the ULP level and so flipped the discrete PC-score class
# partition (the prior accepted COMPRESS divergence). Translated statement-for-statement
# from eigen.f — same Jacobi rotation order, sin/cos formula, and final descending sort —
# so the eigenvectors (and thus the partition) are bit-reproducible against Fortran.
#
# `xtx` is the full 5×5 correlation matrix; only its upper triangle (i≤j) is read.
function _ibm_eigen(xtx::Matrix{Float64}, n::Int)
    a = zeros(Float64, n * (n + 1) ÷ 2)                 # column-packed upper triangle
    @inbounds for j in 1:n, i in 1:j
        a[i + (j * j - j) ÷ 2] = xtx[i, j]
    end
    r = zeros(Float64, n * n)                           # eigenvectors, column-wise
    iq = -n                                             # generate identity (MV=0)
    @inbounds for j in 1:n
        iq += n
        for i in 1:n
            r[iq + i] = (i == j) ? 1.0 : 0.0
        end
    end
    rng = 1.0e-12
    dn = Float64(n)
    anorm = 0.0                                         # off-diagonal norm
    @inbounds for i in 1:n, j in i:n
        i == j && continue
        ia = i + (j * j - j) ÷ 2
        anorm += a[ia] * a[ia]
    end
    if anorm > 0.0
        anorm = sqrt(anorm * 2.0)
        anrmx = anorm * rng / dn
        ind = 0
        thr = anorm
        while true                                      # threshold loop (45)
            thr /= dn
            while true                                  # sweep loop (50)
                @inbounds for l in 1:n-1
                    lq = (l * l - l) ÷ 2
                    for m in l+1:n
                        mq = (m * m - m) ÷ 2
                        lm = l + mq
                        abs(a[lm]) - thr < 0.0 && continue   # below threshold (130)
                        ind = 1
                        ll = l + lq; mm = m + mq
                        x = (a[ll] - a[mm]) / 2.0
                        y = -a[lm] / sqrt(a[lm] * a[lm] + x * x)
                        x <= 0.0 && (y = -y)
                        sinx = y / sqrt(2.0 * (1.0 + sqrt(1.0 - y * y)))
                        sinx2 = sinx * sinx
                        cosx = sqrt(1.0 - sinx2)
                        cosx2 = cosx * cosx
                        sincs = sinx * cosx
                        ilq = n * (l - 1); imq = n * (m - 1)
                        for i in 1:n                    # rotate L and M cols
                            iqi = (i * i - i) ÷ 2
                            if i != l && i != m
                                im = i - m > 0 ? m + iqi : i + mq
                                il = i - l >= 0 ? l + iqi : i + lq
                                xx = a[il] * cosx - a[im] * sinx
                                a[im] = a[il] * sinx + a[im] * cosx
                                a[il] = xx
                            end
                            ilr = ilq + i; imr = imq + i
                            xx = r[ilr] * cosx - r[imr] * sinx
                            r[imr] = r[ilr] * sinx + r[imr] * cosx
                            r[ilr] = xx
                        end
                        xx = 2.0 * a[lm] * sincs
                        yv = a[ll] * cosx2 + a[mm] * sinx2 - xx
                        xv = a[ll] * sinx2 + a[mm] * cosx2 + xx
                        a[lm] = (a[ll] - a[mm]) * sincs + a[lm] * (cosx2 - sinx2)
                        a[ll] = yv
                        a[mm] = xv
                    end
                end
                if ind == 1                             # a rotation occurred → re-sweep (150)
                    ind = 0
                else
                    break
                end
            end
            thr - anrmx > 0.0 || break                  # converged at final norm (160)
        end
    end
    # sort eigenvalues descending + carry eigenvectors (165)
    iq = -n
    @inbounds for i in 1:n
        iq += n
        ll = i + (i * i - i) ÷ 2
        jq = n * (i - 2)
        for j in i:n
            jq += n
            mm = j + (j * j - j) ÷ 2
            a[ll] - a[mm] >= 0.0 && continue
            x = a[ll]; a[ll] = a[mm]; a[mm] = x
            for k in 1:n
                ilr = iq + k; imr = jq + k
                x = r[ilr]; r[ilr] = r[imr]; r[imr] = x
            end
        end
    end
    return r
end

# MEANSD (meansd.f): two-pass mean + sample variance/std of n values from getv(i).
function _cmp_meansd(getv, n::Int)
    apx = 0.0
    for i in 1:n; apx += getv(i); end
    apx /= n
    sm = 0.0; ssq = 0.0
    for i in 1:n
        z = getv(i) - apx; sm += z; ssq += z * z
    end
    abar = sm / n
    var = n > 1 ? (ssq - sm * abar) / (n - 1) : 0.0
    return apx + abar, var, sqrt(max(var, 0.0))
end

# Sort the sub-permutation ind[lo:lo+len-1] by key[ind[*]] (descending if `desc`),
# Range (max − min) of key over the records ind[lo:lo+len-1] (CMRANG, comprs.f).
function _cmp_range(ind::AbstractVector{<:Integer}, lo::Int, len::Int, key::Vector{Float64})
    len < 1 && return 0.0
    mn = key[ind[lo]]; mx = mn
    @inbounds for k in 1:len-1
        v = key[ind[lo+k]]; v < mn && (mn = v); v > mx && (mx = v)
    end
    return mx - mn
end

"""
    compress!(s, nclas, pn1) -> Bool

Reduce the live tree list to `nclas` representative records (comprs.f). `pn1` is the
fraction of classes found by Method 1 (gap clustering); the rest by Method 2 (range
splitting). No-op (returns false) unless `nclas < n`. Conserves total TPA exactly.
"""
function compress!(s::StandState, nclas::Int, pn1::Float64)::Bool
    t = s.trees; n = t.n
    (nclas >= n || n < 2) && return false
    pifac = Float64(s.plot.pi)

    # --- classification variables + standardization (comprs.f:160-236) ----------
    htv(i)  = Float64(t.height[i]);  icrv(i) = Float64(t.crown_pct[i])
    imcv(i) = Float64(t.mort_code[i])
    # comprs.f:30 computes WK4(I)=ALOG(DBH(I)) in SINGLE precision (ALOG of a REAL*4 DBH,
    # stored in the single-precision WK4); match that exactly — log(::Float32)→Float32 then
    # promote — so the ln(DBH) classification variable is bit-identical to Fortran.
    lndv(i) = Float64(log(t.dbh[i]))
    dgv(i)  = Float64(t.diam_growth[i])
    getters = (htv, icrv, imcv, lndv, dgv)
    rmean = zeros(5); stdv = zeros(5)
    for k in 1:5
        m, _, sd = _cmp_meansd(getters[k], n)
        rmean[k] = m; stdv[k] = max(sd, _CMP_STDFLOOR[k])
    end
    # correlation matrix (5×5) via centered/scaled cross-products
    obs = zeros(5); xsum = zeros(5); xtx = zeros(5, 5)
    for i in 1:n
        for k in 1:5; obs[k] = (getters[k](i) - rmean[k]) / stdv[k]; end
        for a in 1:5
            xsum[a] += obs[a]
            for b in 1:a; xtx[a, b] += obs[a] * obs[b]; end
        end
    end
    for a in 1:5, b in 1:a
        xtx[a, b] = (xtx[a, b] - xsum[a] * xsum[b] / n) / (n - 1)
        xtx[b, a] = xtx[a, b]
    end
    for k in 1:5; xtx[k, k] = 1.0; end

    # --- eigen (faithful IBM-SSP EIGEN port) + sign-fix + scale -----------------
    # EIGEN returns eigenvectors column-wise in descending eigenvalue order, so PC1 is
    # column 1 (=R[1:5]) and PC2 is column 2 (=R[6:10]) — matching EIVECT(1:5)/(6:10).
    R = _ibm_eigen(xtx, 5)
    pc1 = R[1:5]; pc2 = R[6:10]
    pc1[4] < 0 && (pc1 .= .-pc1)            # sign-fix on ln(DBH) loading (EIVECT(4))
    pc2[2] > 0 && (pc2 .= .-pc2)            # sign-fix on ICR loading (EIVECT(7))
    e1 = pc1 ./ stdv; e2 = pc2 ./ stdv      # scale by 1/STDDEV(k)

    # --- PC scores: WK3 = species/point base + PC1 + 4 ; WK4 = PC2 (comprs.f:301)
    wk3 = zeros(n); wk4 = zeros(n)
    for i in 1:n
        c1 = htv(i) - rmean[1]; c2 = icrv(i) - rmean[2]; c3 = imcv(i) - rmean[3]
        c4 = lndv(i) - rmean[4]; c5 = dgv(i) - rmean[5]
        base = 25.0 * (pifac * Float64(t.species[i]) + Float64(t.plot_id[i]))
        wk3[i] = base + (c1*e1[1] + c2*e1[2] + c3*e1[3] + c4*e1[4] + c5*e1[5]) + 4.0
        wk4[i] = c1*e2[1] + c2*e2[2] + c3*e2[3] + c4*e2[4] + c5*e2[5]
    end

    # === STEP4: sort IND on WK3 (descending) via RDPSRT (comprs.f:322) ==========
    # Faithful transliteration of comprs.f STEP4–STEP6: identical RDPSRT/IQRSRT sorts
    # (Quickersort tie order) and split mechanics, so the discrete class partition is
    # bit-reproducible against Fortran. WK6/WK5/IND2 carry per-class ranges/lengths;
    # IND1 holds class END positions in the IND ordering. Arrays are sized n (≥ the
    # n−1 / nclas extents Fortran's MAXTRE-sized arrays use here).
    ind = zeros(Int, n)
    rdpsrt!(n, wk3, ind, true)                       # IND ← 1..n, sorted descending on WK3

    # STEP5(A): score differences WK6(J−1)=WK3 gap; count near-ties (≤ ALGTOL).
    wk6 = zeros(n); izers = 0
    x1 = wk3[ind[1]]
    for j in 2:n
        x2 = wk3[ind[j]]
        wk6[j-1] = x1 - x2
        wk6[j-1] <= _CMP_ALGTOL && (izers += 1)
        x1 = x2
    end

    # STEP5(B): number of Method-1 classes.  INT() truncates (comprs.f:343).
    ncls1 = trunc(Int, nclas * pn1 + 0.5)
    isig = n - izers
    isig < ncls1 && (ncls1 = isig)
    ncls1 < 1 && (ncls1 = 1)
    ncls1 > nclas && (ncls1 = nclas)

    # STEP5(C/D): pointer list of the largest gaps (RDPSRT desc), then the first
    # NCLS1−1 of them sorted ascending (IQRSRT); the last class ends at n.
    ind1 = zeros(Int, n)
    rdpsrt!(n - 1, wk6, ind1, true)                  # IND1 ← 1..n−1, sorted descending on WK6
    ncls1 > 1 && iqrsrt!(ind1, ncls1 - 1)
    ind1[ncls1] = n
    ind2 = zeros(Int, n)                             # class lengths
    wk5 = zeros(n)                                    # per-class PC2 (WK4) range
    wk6c = zeros(n)                                   # per-class PC1 (WK3) range (WK6 reused in Fortran)

    # === STEP6: Method 2 — split the widest classes NCLS2 more times ============
    ncls2 = nclas - ncls1
    if ncls2 > 0
        # (B) per-class PC1 range (WK6), length (IND2), PC2 range (WK5/CMRANG).
        let i1 = 1
            for I in 1:ncls1
                i2 = ind1[I]
                jj = ind[i2]; kk = ind[i1]
                wk6c[I] = wk3[kk] - wk3[jj]
                len = i2 - i1 + 1
                ind2[I] = len
                wk5[I] = _cmp_range(ind, i1, len, wk4)
                i1 = i2 + 1
            end
        end
        # (C) split loop: at most NCLS2 new classes.
        for _K in 1:ncls2
            # largest range across classes; LTWO ⇒ PC2 (second score) is the wider axis.
            irec = 0; xrang = 0.0; ltwo = false
            for J in 1:ncls1
                x = wk6c[J]; l2 = wk5[J] > x; l2 && (x = wk5[J])
                x <= xrang && continue
                ltwo = l2; xrang = x; irec = J
            end
            (xrang <= _CMP_RNGMIN || irec == 0) && break
            len = ind2[irec]
            i1 = ind1[irec] - len + 1
            ltwo && rdpsrt!(len, wk4, view(ind, i1:i1+len-1), false)   # DESCENDING on WK4
            i2 = ind1[irec]
            xsmal = ltwo ? wk4[ind[i2]] : wk3[ind[i2]]
            # split point = class median, nudged DOWN off near-equal boundaries.
            if len <= 2
                len = 1
            else
                len = (len + 1) ÷ 2
            end
            jk = ind2[irec] - len - 1
            if jk >= 1
                for _J in 1:jk
                    ii = ind1[irec] - len
                    if ltwo
                        (wk4[ind[ii]] - wk4[ind[ii+1]] > 0.001) && break
                    else
                        (wk3[ind[ii]] - wk3[ind[ii+1]] > 0.001) && break
                    end
                    len += 1
                end
            end
            # perform the split: upper part stays IREC, lower LEN records → class NCLS1+1.
            ncls1 += 1
            ind1[ncls1] = ind1[irec]
            isp = ind1[irec] - len                   # I: last member of the upper class
            ind1[irec] = isp
            ind2[ncls1] = len
            ind2[irec] -= len
            jkk = ind[i1]
            if !ltwo
                j = ind[isp+1]; wk6c[ncls1] = wk3[j] - xsmal
                j = ind[isp];   wk6c[irec]  = wk3[jkk] - wk3[j]
                wk5[irec]  = _cmp_range(ind, i1, ind2[irec], wk4)
                wk5[ncls1] = _cmp_range(ind, isp + 1, ind2[ncls1], wk4)
            else
                j = ind[isp+1]; wk5[ncls1] = wk4[j] - xsmal
                j = ind[isp];   wk5[irec]  = wk4[jkk] - wk4[j]
                len2 = ind2[irec]
                rdpsrt!(len2, wk3, view(ind, i1:i1+len2-1), false)
                len3 = ind2[ncls1]
                rdpsrt!(len3, wk3, view(ind, isp+1:isp+len3), false)
                jkk = ind[i1]; j = ind[isp]
                wk6c[irec] = wk3[jkk] - wk3[j]
                j = ind[i2]; jkk = ind[isp+1]
                wk6c[ncls1] = wk3[jkk] - wk3[j]
            end
        end
        nclas = ncls1
        iqrsrt!(ind1, nclas)                         # (E) resort the class boundaries ascending
    end

    # --- build class member lists from (ind, ind1[1:nclas]) ---------------------
    classes = Vector{Vector{Int}}()
    let i1 = 1
        for c in 1:nclas
            i2 = ind1[c]
            i2 >= i1 && push!(classes, [ind[k] for k in i1:i2])
            i1 = i2 + 1
        end
    end

    # --- merge each class into one PROB-weighted record (comprs.f:648-…) --------
    _merge_classes!(s, classes)
    return true
end

# Merge one class into a single record — comprs.f:600-930, in its REAL*4 arithmetic and member order.
# `mem_in` is the class in IND order; DO 210 (comprs.f:605-612) swaps the smallest record index to the front (IREC1,
# the survivor slot), permuting the rest — that permuted order is the accumulation order of every sum below. The
# weight is WK5 = PROB + WK2 (DO 212; WK2 = 0 on both callers: COMCUP and RDMN1 zero it), WK3 its running sum, TXP
# the class total; one RANN picks the record the nominal attributes come from (DO 216). A one-record class is left
# untouched (GOTO 480). MEASURED FVSbc_instr Fir.20 rootdis RDMN1 COMPRS(400,0.5): the jl Float64 means were 1-2
# ULP off on HT/DBH/PROB of ~130 of 400 records and ICR on 31 (NINT vs a rounded Float64 mean).
function _merge_one!(s::StandState, mem_in::Vector{Int})::Int
    t = s.trees
    mem = copy(mem_in)
    @inbounds for ii in 2:length(mem)                       # DO 210: IND(I1) ↔ IND(II) when IND(I1) > IND(II)
        mem[1] > mem[ii] && ((mem[1], mem[ii]) = (mem[ii], mem[1]))
    end
    dst = mem[1]
    length(mem) == 1 && return dst
    m = length(mem)
    wk5 = Float32[t.tpa[r] + 0f0 for r in mem]            # WK5 = PROB + WK2
    wk3 = similar(wk5); txp = wk5[1]; wk3[1] = txp
    @inbounds for k in 2:m; txp += wk5[k]; wk3[k] = txp; end
    txp == 0f0 && return dst
    x = rann!(s.rng) * txp                                 # CALL RANN(X); X = X*TXP (REAL)
    sel = mem[m]
    @inbounds for k in 1:m
        x <= wk3[k] && (sel = mem[k]; break)
    end
    ltrnk = t.norm_ht[sel] > 0
    # comprs.f:733-745 nominal attributes from the sampled record; IDTREE = IDCMP2 + IY(MAX(1,ICYC)); MISGET/MISPUT.
    t.plot_id[dst] = t.plot_id[sel];   t.species[dst] = t.species[sel]
    t.cut_code[dst] = t.cut_code[sel]; t.special[dst] = t.special[sel]
    t.mort_code[dst] = t.mort_code[sel]; t.iestat[dst] = t.iestat[sel]
    t.tree_id[dst] = Int32(20000000 + cycle_year_at(s.control, Int(s.control.cycle)))   # IY(MAX(1,ICYC)), ICYC = jl cycle+1 (0 ⇒ IY(1) at RDMN1)
    t.dmr[dst] = t.dmr[sel]
    nrm(r) = Float32(t.norm_ht[r]) / 100f0
    if !ltrnk                                              # comprs.f:752-768
        hti = t.norm_ht[dst] > 0 ? nrm(dst) : t.height[dst]
        hti = hti * wk5[1]
        t.norm_ht[dst] = Int32(0); t.trunc[dst] = Int32(0)
        @inbounds for k in 2:m
            r = mem[k]
            xh = t.norm_ht[r] > 0 ? nrm(r) : t.height[r]
            hti = hti + xh * wk5[k]
        end
        t.height[dst] = hti / txp
    else                                                   # comprs.f:770-799
        xx = 0f0; xp = 0f0
        @inbounds for k in 1:m
            r = mem[k]; t.norm_ht[r] > 0 || continue
            xx = xx + Float32(t.trunc[r]) / 100f0 * wk5[k]
            xp = xp + t.height[r] * wk5[k]
        end
        xx = xx / xp
        xit = 0f0; xnr = 0f0; hti = 0f0
        @inbounds for k in 1:m
            r = mem[k]; w = wk5[k]
            hti = hti + t.height[r] * w
            if t.norm_ht[r] > 0
                xnr = xnr + nrm(r) * w
                xit = xit + Float32(t.trunc[r]) / 100f0 * w
            else
                xnr = xnr + t.height[r] * w
                xit = xit + t.height[r] * xx * w
            end
        end
        t.height[dst] = hti / txp
        t.norm_ht[dst] = unsafe_trunc(Int32, xnr / txp * 100f0)
        t.trunc[dst] = unsafe_trunc(Int32, xit / txp * 100f0)
    end
    # comprs.f:800-929 PROB-weighted sums (REAL), first member then the rest in order.
    arrs = (t.bdft_vol, t.cuft_vol, t.merch_cuft_vol, t.saw_cuft_vol, t.cull,
            t.abvgrd_bio, t.merch_bio, t.cubsaw_bio, t.foliage_bio,
            t.abvgrd_carb, t.merch_carb, t.cubsaw_carb, t.foliage_carb, t.carbon_frac,
            t.merch_top_bf, t.merch_top_cf, t.diam_growth, t.ht_growth, t.old_crown_pct, t.crown_ratio,
            t.dg_prev, t.crown_width)
    acc = zeros(Float32, length(arrs))
    dbhi = 0f0; probi = 0f0; xicri = 0f0; decayi = 0; wdldi = 0
    df = zeros(Float32, 4)
    dsplit(dv) = (Float32(dv ÷ 1000000), Float32(dv ÷ 10000 - (dv ÷ 1000000) * 100),
                  Float32(dv ÷ 100 - (dv ÷ 10000) * 100), Float32(dv - (dv ÷ 100) * 100))
    @inbounds for k in 1:m
        r = mem[k]; xp = wk5[k]
        for (a, arr) in enumerate(arrs)
            acc[a] = k == 1 ? arr[r] * xp : acc[a] + arr[r] * xp
        end
        dr = t.dbh[r]
        dbhi = k == 1 ? dr * dr * xp : dbhi + dr * dr * xp
        probi = k == 1 ? t.tpa[r] : probi + t.tpa[r]
        xicri = k == 1 ? Float32(t.crown_pct[r]) * xp : xicri + Float32(t.crown_pct[r]) * xp
        # DECAYI / WDLDSTEMI are INTEGER: each assignment truncates the REAL expression.
        decayi = k == 1 ? unsafe_trunc(Int, Float32(t.decay_code[r]) * xp) :
                          unsafe_trunc(Int, Float32(decayi) + Float32(t.decay_code[r]) * xp)
        wdldi = k == 1 ? unsafe_trunc(Int, Float32(t.woodland_stems[r]) * xp) :
                         unsafe_trunc(Int, Float32(wdldi) + Float32(t.woodland_stems[r]) * xp)
        d1, d2, d3, d4 = dsplit(Int(t.defect[r]))
        df[1] = k == 1 ? d1 * xp : df[1] + d1 * xp
        df[2] = k == 1 ? d2 * xp : df[2] + d2 * xp
        df[3] = k == 1 ? d3 * xp : df[3] + d3 * xp
        df[4] = k == 1 ? d4 * xp : df[4] + d4 * xp
    end
    @inbounds for (a, arr) in enumerate(arrs); arr[dst] = acc[a] / txp; end
    t.decay_code[dst] = Int32(unsafe_trunc(Int, Float32(decayi) / txp))
    t.woodland_stems[dst] = Int32(unsafe_trunc(Int, Float32(wdldi) / txp))
    t.dbh[dst] = sqrt(dbhi / txp)
    t.crown_pct[dst] = Int32(round(Int, xicri / txp, RoundNearestTiesAway))   # NINT
    t.tpa[dst] = probi
    idf = [unsafe_trunc(Int, df[j] / txp + 0.5f0) for j in 1:4]
    t.defect[dst] = Int32(idf[1] * 1000000 + idf[2] * 10000 + idf[3] * 100 + idf[4])
    t.zrand[dst] = -999f0                                  # comprs.f:929 reset the serial correlation
    return dst
end

# Merge all classes and compact the live tree list to one record per class (comprs.f merge + TREDEL).
# FAITHFUL post-COMPRESS ORDER (#29): the merge VALUES incl. the OLDRN serial-corr deviate are correct via
# _merge_one! (debug-FVS dump verified bit-exact). The FINAL RECORD ORDER must match FVS's compaction
# `CALL TREDEL(ITRN-NCLAS, IND)` at comprs.f:1007 — the same cycle's growth draws `rann!` per record IN ORDER, so
# a wrong order desyncs the RNG and drifts the post-compression DG. A debug-FVS dump of IREC1 proved FVS's
# survivor slot = the **minimum** record index of each class (= `dst=minimum(mem)`, where _merge_one! already
# writes the merged record), NOT mem[1]. TREDEL then negates the non-survivor members (comprs.f:863-869) and
# fills the smallest vacancy with the largest survivor (tredel.f:46-120). We simulate that on the dst (=IREC1)
# slots to get FVS's exact order, then permute the merged records into it.
function _merge_classes!(s::StandState, classes::Vector{Vector{Int}})
    t = s.trees
    nc = length(classes)
    nc == 0 && return
    itrn = t.n
    dstc = Int[_merge_one!(s, mem) for mem in classes]   # merged record at dst=min(mem) = FVS IREC1 (survivor)
    # TREDEL order: survivor of class c at slot dstc[c]; every other 1..itrn slot is a vacancy. Fill smallest
    # vacancy with the largest survivor while IVAC ≤ IREC; slot_class[1..nc] after the swaps is the FVS order.
    slot_class = zeros(Int, itrn)
    @inbounds for c in 1:nc; slot_class[dstc[c]] = c; end
    survs = sort(dstc)
    vacs  = Int[i for i in 1:itrn if slot_class[i] == 0]   # ascending (i increasing)
    iv = 1; ir = length(survs)
    @inbounds while iv <= length(vacs) && ir >= 1
        vacs[iv] > survs[ir] && break                      # IVAC > IREC ⇒ done (tredel.f:86)
        slot_class[vacs[iv]]  = slot_class[survs[ir]]      # TREMOV survivor → vacancy
        slot_class[survs[ir]] = 0
        iv += 1; ir -= 1
    end
    fvs_order = Int[slot_class[p] for p in 1:nc]           # class at each final position p (FVS order)
    # gather merged records to 1..nc by ascending dst (safe: dst ≥ target), then permute into FVS order.
    order = sortperm(dstc)
    cur_pos = Vector{Int}(undef, nc)
    @inbounds for j in 1:nc
        c = order[j]; cur_pos[c] = j
        dstc[c] != j && copy_tree!(t, j, dstc[c])
    end
    perm = Int[cur_pos[fvs_order[p]] for p in 1:nc]        # new[p] = old[perm[p]]
    placed = falses(nc); temp = nc + 1                     # temp: a free slot (nc < itrn when compressing)
    @inbounds for p0 in 1:nc
        (placed[p0] || perm[p0] == p0) && (placed[p0] = true; continue)
        copy_tree!(t, temp, p0)
        i = p0
        while true
            j = perm[i]; placed[i] = true
            if j == p0; copy_tree!(t, i, temp); break; end
            copy_tree!(t, i, j); i = j
        end
    end
    if t.ndead > 0                                          # shift the dead block down behind the reps
        @inbounds for k in 1:t.ndead
            copy_tree!(t, nc + k, itrn + k)
        end
    end
    # FVS SPESRT→LNKCHN→SETUP rebuilds the species sort after COMPRESS, re-listing each species in ASCENDING
    # PHYSICAL record index (lnkchn.f is called I=1..IREC1 in order; the tripled lineage order is discarded).
    # Reset sort_key to the compacted position so species_sort! walks each species in physical order — this is
    # the order DGSCOR draws the per-tree RNG, so without it two same-species merged records get SWAPPED bachlo
    # draws (s22: the two sp33 records). Mirrors core/trees.jl compact! (the thin/TREDEL path).
    @inbounds for i in 1:nc; t.sort_key[i] = Float64(i); end
    t.n = nc
    return
end

"""
    apply_compress!(s) -> Bool

COMPRESS (act 250) driver (comcup.f). If a COMPRESS is scheduled for the current cycle and the
record count exceeds the target, compress the live tree list (the LAST request wins if several are
due). Returns whether compression fired — the caller then suppresses record tripling this cycle
(NOTRIP). PN1 is stored as a percent and divided by 100 here (comcup.f:93).
"""
function apply_compress!(s::StandState)::Bool
    isempty(s.control.schedule) && return false
    yr = current_cycle_year(s); fvscyc = Int(s.control.cycle) + 1
    # NOTE: this uses an exact `date == cycle_start` gate (NOT the OPCYCL containing-cycle gate that SETSITE/
    # SIMFIRE/VOLUME use). A calendar-dated mid-cycle COMPRESS (e.g. 2005 in a 10-yr NE cycle) therefore does
    # not fire — but the OPCYCL gate was TRIED (2026-07-07) and REVERTED: it made the mid-cycle compress FIRE
    # density-exact (TPA/BA/SDI bit-exact vs live) yet volume-DIVERGENT (~2-5% by 2040, compounding from the
    # merge), i.e. the same COMPRESS eigensolver near-tie record-distribution sensitivity that corners s22 — so
    # the mid-cycle firing point/merge is NOT bit-exact-validatable here. The common COMPRESS usage is a
    # CYCLE-NUMBER date (<1000, fvscyc path, tested bit-exact), so this exact gate is correct for every exercised
    # case; the mid-cycle-calendar path is a deferred open question (needs a live compress-timing + eigensolver trace).
    due = nothing
    for a in s.control.schedule
        a.icflag == Int32(250) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        due = a                                        # keep the last matching request
    end
    due === nothing && return false
    nclas = round(Int, due.params[1]); pn1 = Float64(due.params[2]) / 100
    # A non-positive target (blank/malformed field, or ≥ the record count) can't compress to fewer
    # classes — no-op, matching live FVS (which produces output rather than crashing). Guards the
    # compress! STEP5 `ind1[ncls1]=n` where ncls1 would clamp to nclas≤0 (comcup.f only fires when
    # ITARG>0 and NCLAS>ITARG).
    (nclas < 1 || nclas >= s.trees.n) && return false
    fired = compress!(s, nclas, pn1)
    if fired
        compute_density!(s)                            # density changed after the merge
        s.control.no_tripling = true                   # COMCUP sets NOTRIP=.TRUE. (comcup.f:126) — but LTRIP for
                                                       # THIS cycle was already latched (grincr.f:74 precedes the
                                                       # :391 COMCUP), so tripling still fires THIS cycle; NOTRIP
                                                       # suppresses only SUBSEQUENT cycles.
    end
    return fired
end
