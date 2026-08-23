# =============================================================================
# PPE (Parallel Processing Extension) — deterministic landscape sort kernels.
# STAGED first chunk (2026-08-23). Bit-exact port of the PPBASE index QuickerSort
# (Scowen 1965, ALGORITHM 271) used by PPMAIN/BYGRPS for master-cycle stand
# ordering. Oracle = per-routine gfortran-16 golden (the full FVSppe landscape exe
# cannot link against modern base-FVS: getstd/putstd marshal a deleted 2014 COMMON
# layout — BB0..BB13/D0/D0MULT — see staged/VALIDATION.md).
#
# C11SRT / C26SRT / CH8SRT are the SAME algorithm at CHARACTER widths 11/26/8; the
# physical array A is never permuted — only INDEX is rearranged so that for each i,
# A[INDEX[i]] <= A[INDEX[i+1]]. Ties are placed by the partition's swap sequence
# (NOT stable), so the exact control flow below is load-bearing for bit-exactness.
#
# ADDITIVE + INERT: staged, wired into nothing. A single-stand projection is
# untouched (the gate multicycle 339/11 is byte-identical).
# =============================================================================

"""
    ppe_index_qsort!(index, keys; lseq) -> index

Faithful bit-exact port of PPBASE C11SRT/C26SRT/CH8SRT (character index
QuickerSort). `keys` is the physical array A (1-based, never permuted). `index`
holds 1-based indices into `keys`; on return its first `n = length(index)`
entries are rearranged so `keys[index[i]] <= keys[index[i+1]]`.

If `lseq` is true, `index` is first loaded 1:n (sort the first n keys). If false,
`index` is taken as the caller-supplied subset/order to sort over (matching the
`.FALSE.` calls in PPMAIN/BYGRPS). Comparison is Fortran CHARACTER lexicographic
order (ASCII); pass keys already space-padded to the fixed field width so the
collation matches the oracle exactly.
"""
function ppe_index_qsort!(index::AbstractVector{<:Integer}, keys::AbstractVector{<:AbstractString}; lseq::Bool)
    n = length(index)
    if lseq
        @inbounds for i in 1:n
            index[i] = i
        end
    end
    n < 2 && return index

    # 1:1 transliteration of C11SRT/C26SRT (Fortran labels preserved as @label Lnn).
    ipush = Vector{Int}(undef, 33)   # IPUSH(33)
    itop = 0
    il = 1
    iu = n
    # locals reused across labels (declare so @goto scoping is happy)
    indil = indiu = indip = indkl = indku = 0
    kl = ku = ip = jl = ju = 0
    t = keys[1]

    @inbounds begin
        @label L30
        if iu <= il ; @goto L40 ; end
        indil = index[il]
        indiu = index[iu]
        if iu > il + 1 ; @goto L50 ; end
        if keys[indil] <= keys[indiu] ; @goto L40 ; end
        index[il] = indiu
        index[iu] = indil

        @label L40
        itop == 0 && return index
        il = ipush[itop-1]
        iu = ipush[itop]
        itop -= 2
        @goto L30

        @label L50
        ip = (il + iu) ÷ 2
        indip = index[ip]
        t = keys[indip]
        index[ip] = indil
        kl = il
        ku = iu

        @label L60
        kl += 1
        if kl > ku ; @goto L90 ; end
        indkl = index[kl]
        if keys[indkl] <= t ; @goto L60 ; end

        @label L70
        indku = index[ku]
        if ku < kl ; @goto L100 ; end
        if keys[indku] < t ; @goto L80 ; end
        ku -= 1
        @goto L70

        @label L80
        index[kl] = indku
        index[ku] = indkl
        ku -= 1
        @goto L60

        @label L90
        indku = index[ku]

        @label L100
        index[il] = indku
        index[ku] = indip
        if ku <= ip ; @goto L110 ; end
        jl = il
        ju = ku - 1
        il = ku + 1
        @goto L120

        @label L110
        jl = ku + 1
        ju = iu
        iu = ku - 1

        @label L120
        itop += 2
        ipush[itop-1] = jl
        ipush[itop] = ju
        @goto L30
    end
end
