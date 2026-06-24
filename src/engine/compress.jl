# =============================================================================
# compress.jl — COMPRESS (act 250): reduce the tree list to NCLAS representative
# records by principal-component-score clustering, then merge each class into one
# record (PROB-conserving). Ported from base/comprs.f + base/comcup.f (the act-250
# driver), with the 1966 IBM-SSP symmetric eigensolver replaced by
# `LinearAlgebra.eigen` (project direction: use Julia's linear algebra rather than
# re-port the Jacobi routine). Because the eigenvectors from LAPACK differ from the
# SSP routine at the ULP level, the PC scores — and thus the exact class partition —
# are NOT guaranteed bit-identical to Fortran; the merge still conserves total TPA
# exactly and reduces the record count to NCLAS. comcup.f calls this with WK2 (the
# mortality weight) zeroed, so the per-class weight is just PROB.
#
# Classification variables (comprs.f:161): HT, ICR, IMC, ln(DBH), DG (NRANK=5).
# =============================================================================

using LinearAlgebra

const _CMP_RNGMIN = 0.00001
const _CMP_ALGTOL = 0.066
# minimum standard deviations per classification variable (comprs.f:171-188)
const _CMP_STDFLOOR = (1.0, 0.1e-3, 0.1e-3, 5e-3, 0.02)

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
# stable with the record index as the tie-break (RDPSRT stand-in; see file header on
# why exact Fortran tie-break parity is moot once LAPACK eigenvectors are used).
function _cmp_subsort!(ind::Vector{Int}, lo::Int, len::Int, key::Vector{Float64}, desc::Bool)
    len < 2 && return
    sub = ind[lo:lo+len-1]
    sort!(sub; alg = MergeSort, by = r -> (desc ? -key[r] : key[r], r))
    @inbounds for k in 1:len; ind[lo+k-1] = sub[k]; end
    return
end

# Range (max − min) of key over the records ind[lo:lo+len-1] (CMRANG, comprs.f).
function _cmp_range(ind::Vector{Int}, lo::Int, len::Int, key::Vector{Float64})
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
    imcv(i) = Float64(t.mort_code[i]); lndv(i) = log(Float64(t.dbh[i]))
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

    # --- eigen (Julia, replacing IBM-SSP EIGEN) + sign-fix + scale --------------
    F = eigen(Symmetric(xtx))              # ascending eigenvalues
    pc1 = F.vectors[:, 5]                   # largest eigenvalue ⇒ PC1
    pc2 = F.vectors[:, 4]                   # second ⇒ PC2
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

    # --- Method 1: class breaks at the largest gaps in sorted WK3 (comprs.f:320-362)
    ind = collect(1:n)
    _cmp_subsort!(ind, 1, n, wk3, true)             # descending
    wk6 = zeros(n); izers = 0
    for j in 1:n-1
        wk6[j] = wk3[ind[j]] - wk3[ind[j+1]]
        wk6[j] <= _CMP_ALGTOL && (izers += 1)
    end
    ncls1 = floor(Int, nclas * pn1 + 0.5)
    isig = n - izers
    ncls1 = clamp(min(ncls1, isig), 1, nclas)
    gapord = sortperm(1:n-1; alg = MergeSort, by = j -> (-wk6[j], j))   # gaps, largest first
    ind1 = Int[]                                      # class END positions in `ind`
    if ncls1 > 1
        append!(ind1, sort(gapord[1:ncls1-1]))
    end
    push!(ind1, n)                                    # last class ends at n
    ind2 = Int[]                                      # class lengths
    let prev = 0
        for e in ind1; push!(ind2, e - prev); prev = e; end
    end

    # --- Method 2: split the largest-range classes NCLS2 more times (comprs.f:371-553)
    ncls2 = nclas - ncls1
    if ncls2 > 0
        wk5 = zeros(nclas); wk6c = zeros(nclas)       # per-class PC2 / PC1 ranges
        let i1 = 1
            for c in 1:ncls1
                i2 = ind1[c]
                wk6c[c] = wk3[ind[i1]] - wk3[ind[i2]]
                wk5[c] = _cmp_range(ind, i1, i2 - i1 + 1, wk4)
                i1 = i2 + 1
            end
        end
        for _ in 1:ncls2
            # largest range across classes (PC1 vs PC2 whichever bigger)
            irec = 0; xrang = 0.0; ltwo = false
            for c in 1:ncls1
                x = wk6c[c]; l2 = wk5[c] > x; l2 && (x = wk5[c])
                x <= xrang && continue
                ltwo = l2; xrang = x; irec = c
            end
            (xrang <= _CMP_RNGMIN || irec == 0) && break
            len = ind2[irec]
            i1 = ind1[irec] - len + 1
            ltwo && _cmp_subsort!(ind, i1, len, wk4, false)   # ascending on PC2
            i2 = ind1[irec]
            key = ltwo ? wk4 : wk3
            xsmal = key[ind[i2]]
            # split point: half-way down the (descending) class, nudged off tiny gaps
            local sp::Int
            if len <= 2
                sp = 1
            else
                sp = (len + 1) ÷ 2
                jk = ind2[irec] - sp - 1
                if jk >= 1
                    for _j in 1:jk
                        i = ind1[irec] - sp
                        (key[ind[i]] - key[ind[i+1]] > 0.001) && break
                        sp += 1
                    end
                end
            end
            # perform the split: upper part stays in irec, lower `sp` records → new class
            ncls1 += 1
            push!(ind1, ind1[irec]); push!(ind2, 0); push!(wk5, 0.0); push!(wk6c, 0.0)
            ihi = ind1[irec] - sp
            ind1[irec] = ihi
            ind2[ncls1] = sp
            ind2[irec] -= sp
            jk = ind[i1]
            if !ltwo
                wk6c[ncls1] = wk3[ind[ihi+1]] - xsmal
                wk6c[irec] = wk3[jk] - wk3[ind[ihi]]
                wk5[irec] = _cmp_range(ind, i1, ind2[irec], wk4)
                wk5[ncls1] = _cmp_range(ind, ihi + 1, ind2[ncls1], wk4)
            else
                wk5[ncls1] = wk4[ind[ihi+1]] - xsmal
                wk5[irec] = wk4[jk] - wk4[ind[ihi]]
                _cmp_subsort!(ind, i1, ind2[irec], wk3, false)
                _cmp_subsort!(ind, ihi + 1, ind2[ncls1], wk3, false)
                wk6c[irec] = wk3[ind[i1]] - wk3[ind[ihi]]
                wk6c[ncls1] = wk3[ind[ihi+1]] - wk3[ind[i2]]
            end
        end
        nclas = ncls1
    end
    sort!(ind1)                                      # IQRSRT(IND1,NCLAS): class boundaries ascending

    # --- build class member lists from (ind, ind1) ------------------------------
    classes = Vector{Vector{Int}}()
    let i1 = 1
        for i2 in ind1
            i2 >= i1 && push!(classes, [ind[k] for k in i1:i2])
            i1 = i2 + 1
        end
    end

    # --- merge each class into one PROB-weighted record (comprs.f:648-…) --------
    _merge_classes!(s, classes)
    return true
end

# Merge one class `mem` (record indices) into a single record at its lowest index, PROB-weighted
# (comprs.f:688-966); WK2=0 in the comcup path so the weight is PROB. Conserves the class's TPA.
function _merge_one!(s::StandState, mem::Vector{Int})::Int
    t = s.trees
    dst = minimum(mem)
    length(mem) == 1 && return dst
    txp = 0.0
    for r in mem; txp += t.tpa[r]; end
    txp <= 0.0 && return dst
    # nominal attributes: sample one member proportional to PROB (RANN, comprs.f:725)
    x = Float64(rann!(s.rng)) * txp
    cum = 0.0; sel = mem[end]
    for r in mem
        cum += t.tpa[r]
        if x <= cum; sel = r; break; end
    end
    wmean(f) = (acc = 0.0; for r in mem; acc += f(r) * t.tpa[r]; end; Float32(acc / txp))
    # HT / NORMHT / ITRUNC (truncated-tree handling, comprs.f:756-806)
    ltrnk = t.norm_ht[sel] > 0
    local hnew::Float32; local normnew::Int32; local truncnew::Int32
    if !ltrnk
        acc = 0.0
        for r in mem
            h = t.norm_ht[r] > 0 ? Float64(t.norm_ht[r]) / 100 : Float64(t.height[r])
            acc += h * t.tpa[r]
        end
        hnew = Float32(acc / txp); normnew = Int32(0); truncnew = Int32(0)
    else
        xt = 0.0; xp = 0.0
        for r in mem
            t.norm_ht[r] > 0 || continue
            xt += Float64(t.trunc[r]) / 100 * t.tpa[r]
            xp += Float64(t.height[r]) * t.tpa[r]
        end
        xprop = xp > 0 ? xt / xp : 0.0
        xnr = 0.0; xit = 0.0; hti = 0.0
        for r in mem
            w = t.tpa[r]; hti += Float64(t.height[r]) * w
            if t.norm_ht[r] > 0
                xnr += Float64(t.norm_ht[r]) / 100 * w; xit += Float64(t.trunc[r]) / 100 * w
            else
                xnr += Float64(t.height[r]) * w; xit += Float64(t.height[r]) * xprop * w
            end
        end
        hnew = Float32(hti / txp)
        normnew = round(Int32, xnr / txp * 100); truncnew = round(Int32, xit / txp * 100)
    end
    # quadratic-mean DBH (comprs.f:952)
    dbh2 = 0.0; for r in mem; dbh2 += Float64(t.dbh[r])^2 * t.tpa[r]; end
    # write the merged record into `dst`
    t.height[dst]       = hnew
    t.norm_ht[dst]      = normnew
    t.trunc[dst]        = truncnew
    t.dbh[dst]          = sqrt(Float32(dbh2 / txp))
    t.tpa[dst]          = Float32(txp)                 # PROB summed (conserves TPA)
    t.diam_growth[dst]  = wmean(r -> Float64(t.diam_growth[r]))
    t.ht_growth[dst]    = wmean(r -> Float64(t.ht_growth[r]))
    t.crown_ratio[dst]  = wmean(r -> Float64(t.crown_ratio[r]))
    t.crown_pct[dst]    = round(Int32, wmean(r -> Float64(t.crown_pct[r])))
    t.crown_width[dst]  = wmean(r -> Float64(t.crown_width[r]))
    t.old_crown_pct[dst] = wmean(r -> Float64(t.old_crown_pct[r]))
    t.cuft_vol[dst]       = wmean(r -> Float64(t.cuft_vol[r]))
    t.merch_cuft_vol[dst] = wmean(r -> Float64(t.merch_cuft_vol[r]))
    t.saw_cuft_vol[dst]   = wmean(r -> Float64(t.saw_cuft_vol[r]))
    t.bdft_vol[dst]       = wmean(r -> Float64(t.bdft_vol[r]))
    t.cull[dst]           = wmean(r -> Float64(t.cull[r]))
    t.merch_top_bf[dst]   = wmean(r -> Float64(t.merch_top_bf[r]))
    t.merch_top_cf[dst]   = wmean(r -> Float64(t.merch_top_cf[r]))
    # nominal attributes from the PROB-sampled record (comprs.f:733-741)
    t.species[dst]   = t.species[sel];   t.plot_id[dst]  = t.plot_id[sel]
    t.mort_code[dst] = t.mort_code[sel]; t.cut_code[dst] = t.cut_code[sel]
    t.special[dst]   = t.special[sel];   t.decay_code[dst] = t.decay_code[sel]
    t.defect[dst]    = t.defect[sel];    t.woodland_stems[dst] = t.woodland_stems[sel]
    return dst
end

# Merge all classes and compact the live tree list to one record per class (comcup.f compaction).
function _merge_classes!(s::StandState, classes::Vector{Vector{Int}})
    t = s.trees
    reps = Int[_merge_one!(s, mem) for mem in classes]
    sort!(reps)                                        # final order = ascending record index
    nc = length(reps)
    @inbounds for j in 1:nc
        reps[j] == j && continue
        copy_tree!(t, j, reps[j])                      # gather rep into slot j (reps[j] ≥ j)
    end
    if t.ndead > 0                                     # shift the dead block down behind the reps
        @inbounds for k in 1:t.ndead
            copy_tree!(t, nc + k, t.n + k)
        end
    end
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
    due = nothing
    for a in s.control.schedule
        a.icflag == Int32(250) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        due = a                                        # keep the last matching request
    end
    due === nothing && return false
    nclas = round(Int, due.params[1]); pn1 = Float64(due.params[2]) / 100
    fired = compress!(s, nclas, pn1)
    fired && compute_density!(s)                       # density changed after the merge
    return fired
end
