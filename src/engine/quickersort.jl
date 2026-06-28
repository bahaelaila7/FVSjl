# =============================================================================
# quickersort.jl — RDPSRT + IQRSRT (Scowen 1965, Algorithm 271 "Quickersort")
#
# Faithful statement-for-statement transliterations of base/rdpsrt.f and base/iqrsrt.f.
# COMPRESS (comprs.f) drives its class partition with these exact routines, and the
# Quickersort permutation of EQUAL keys differs from a stable library sort — so to
# reproduce FVS's discrete tree-record clustering bit-for-bit we must use the same
# algorithm, not `sort!`/`sortperm`. (Same reason rdpsrt.jl was transliterated in the
# FVSjulia port; FVSjl now needs it for COMPRESS faithfulness.)
# =============================================================================

"""
    rdpsrt!(n, a, index, lseq)

Indirect DESCENDING sort (rdpsrt.f): rearrange `index[1:n]` so that
`a[index[1]] ≥ a[index[2]] ≥ … ≥ a[index[n]]`. `a` is indexed by the values held in
`index` (global record numbers), so it is the full key array; `index` may be a view
into a sub-range. The physical `a` is not modified. If `lseq`, `index` is first loaded
with 1..n. Ties resolve by the Quickersort partition order (NOT by index), matching FVS.
"""
function rdpsrt!(n::Int, a::AbstractVector{<:Real}, index::AbstractVector{<:Integer}, lseq::Bool)
    if lseq
        @inbounds for i in 1:n
            index[i] = i
        end
    end
    n < 2 && return nothing

    ipush = zeros(Int, 33)
    itop = 0; il = 1; iu = n
    indil = 0; indiu = 0; indip = 0
    indkl = 0; indku = 0
    ip = 0; kl = 0; ku = 0; jl = 0; ju = 0
    t = zero(eltype(a))

    @label l30
    if iu <= il; @goto l40; end
    indil = Int(index[il]); indiu = Int(index[iu])
    if iu > il + 1; @goto l50; end
    if a[indil] >= a[indiu]; @goto l40; end
    index[il] = indiu; index[iu] = indil

    @label l40
    if itop == 0; return nothing; end
    il = ipush[itop-1]; iu = ipush[itop]; itop -= 2
    @goto l30

    @label l50
    ip    = (il + iu) ÷ 2
    indip = Int(index[ip]); t = a[indip]
    index[ip] = indil
    kl = il; ku = iu

    @label l60
    kl += 1
    if kl > ku; @goto l90; end
    indkl = Int(index[kl])
    if a[indkl] >= t; @goto l60; end

    @label l70
    indku = Int(index[ku])
    if ku < kl; @goto l100; end
    if a[indku] > t; @goto l80; end
    ku -= 1
    @goto l70

    @label l80
    index[kl] = indku; index[ku] = indkl; ku -= 1
    @goto l60

    @label l90
    indku = Int(index[ku])

    @label l100
    index[il] = indku; index[ku] = indip
    if ku <= ip; @goto l110; end
    jl = il; ju = ku - 1; il = ku + 1
    @goto l120

    @label l110
    jl = ku + 1; ju = iu; iu = ku - 1

    @label l120
    itop += 2
    ipush[itop-1] = jl; ipush[itop] = ju
    @goto l30
end

"""
    iqrsrt!(list, n)

In-place ASCENDING integer Quickersort (iqrsrt.f): sort `list[1:n]` so that
`list[1] ≤ list[2] ≤ … ≤ list[n]`. Used by COMPRESS to re-sort the class-boundary
pointer list (IND1).
"""
function iqrsrt!(list::AbstractVector{<:Integer}, n::Int)
    n < 2 && return nothing
    iu = zeros(Int, 33); ilst = zeros(Int, 33)
    m = 1; i = 1; j = n
    k = 0; ij = 0; l = 0
    t = 0; tt = 0

    @label l5
    if i >= j; @goto l70; end

    @label l10
    k = i
    ij = (i + j) ÷ 2
    t = Int(list[ij])
    if Int(list[i]) <= t; @goto l20; end
    list[ij] = list[i]; list[i] = t; t = Int(list[ij])

    @label l20
    l = j
    if Int(list[j]) >= t; @goto l40; end
    list[ij] = list[j]; list[j] = t; t = Int(list[ij])
    if Int(list[i]) <= t; @goto l40; end
    list[ij] = list[i]; list[i] = t; t = Int(list[ij])
    @goto l40

    @label l30
    list[l] = list[k]; list[k] = tt

    @label l40
    l -= 1
    if Int(list[l]) > t; @goto l40; end
    tt = Int(list[l])

    @label l50
    k += 1
    if Int(list[k]) < t; @goto l50; end
    if k <= l; @goto l30; end
    if l - i <= j - k; @goto l60; end
    ilst[m] = i; iu[m] = l; i = k; m += 1
    @goto l80

    @label l60
    ilst[m] = k; iu[m] = j; j = l; m += 1
    @goto l80

    @label l70
    m -= 1
    if m <= 0; return nothing; end
    i = ilst[m]; j = iu[m]

    @label l80
    if j - i >= 11; @goto l10; end
    if i == 1; @goto l5; end
    i -= 1

    @label l90
    i += 1
    if i == j; @goto l70; end
    t = Int(list[i+1])
    if Int(list[i]) <= t; @goto l90; end
    k = i

    @label l100
    list[k+1] = list[k]; k -= 1
    if t < Int(list[k]); @goto l100; end
    list[k+1] = t
    @goto l90
end
