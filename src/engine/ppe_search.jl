# =============================================================================
# PPE (PPBASE) — keyed binary search over a character array via an ascending index.
# Bit-exact ports of C26BSR / CH8BSR (return the position in A) and SPBSRX (returns
# the position in IORD). Oracle = per-routine gfortran-16 golden (see VALIDATION.md;
# the full FVSppe exe won't link — getstd/putstd 2014 COMMON drift). These pair with
# ppe_index_qsort!: IORD is an ascending index over A (as produced by the sort), and
# the search finds F in A without A being physically ordered. ADDITIVE + INERT.
#
# C26BSR/CH8BSR are the SAME algorithm at CHARACTER widths 26/8; SPBSRX shares the
# identical search body and differs ONLY in the final return (index-into-IORD vs
# index-into-A). All three are a 1:1 @goto transliteration (labels L20/L30/L40).
# =============================================================================

# Shared search body. Returns (imid, ip) where imid is the IORD position landed on
# and ip = iord[imid]; the caller maps to the routine's documented return.
@inline function _ppe_bsearch_core(keys::AbstractVector{<:AbstractString},
                                   iord::AbstractVector{<:Integer}, f::AbstractString)
    n = length(iord)
    @inbounds begin
        imid = 1
        i1 = iord[1]
        if f <= keys[i1] ; @goto L40 ; end
        imid = n
        in_ = iord[n]
        if f >= keys[in_] ; @goto L40 ; end
        itop = 1
        ibot = n
        @label L20
        imid = (ibot + itop) ÷ 2
        im = iord[imid]
        if f > keys[im] ; @goto L30 ; end
        ibot = imid - 1
        ib = iord[ibot]
        if f > keys[ib] ; @goto L40 ; end
        @goto L20
        @label L30
        itop = imid + 1
        it = iord[itop]
        if f < keys[it] ; @goto L40 ; end
        @goto L20
        @label L40
        return imid, iord[imid]
    end
end

"""
    ppe_c26bsr(keys, iord, f) -> Int

Bit-exact port of PPBASE C26BSR / CH8BSR (CHARACTER*26 / *8 keyed binary search).
`iord` is an ascending index over `keys` (e.g. from `ppe_index_qsort!`). Returns the
subscript in `keys` where `f` is found, or 0 if `f` is not a member. Pass `f` and the
`keys` space-padded to the same fixed width so the ASCII collation matches the oracle.
"""
function ppe_c26bsr(keys::AbstractVector{<:AbstractString},
                    iord::AbstractVector{<:Integer}, f::AbstractString)
    _, ip = _ppe_bsearch_core(keys, iord, f)
    @inbounds return f != keys[ip] ? 0 : ip
end

"""
    ppe_spbsrx(keys, iord, f) -> Int

Bit-exact port of PPBASE SPBSRX. Same search as C26BSR but returns IMID — the
position in `iord` (so `iord[result]` points into `keys`) — or 0 if `f` is absent.
"""
function ppe_spbsrx(keys::AbstractVector{<:AbstractString},
                    iord::AbstractVector{<:Integer}, f::AbstractString)
    imid, ip = _ppe_bsearch_core(keys, iord, f)
    @inbounds return f != keys[ip] ? 0 : imid
end
