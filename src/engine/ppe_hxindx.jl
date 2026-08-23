# =============================================================================
# PPE (PPBASE) — HXINDX hexagonal-grid neighbor indexing (the spatial-neighbor core
# of PPE's inter-stand disturbance spread; called from HVNEDN). Given a reference
# stand's index in a system of NSTND hexagons (rows/cols as equal as possible, an
# extra row if needed) and a direction NEI (1=N 2=NE 3=SE 4=S 5=SW 6=NW), returns a
# pointer to that neighbor, or 0 if it falls off the grid. Bit-exact port; oracle =
# gfortran-16 golden (see VALIDATION.md). ADDITIVE + INERT.
#
# The grid geometry (odd columns shifted up half a cell) makes the NE/SE/SW/NW row
# depend on the column parity — reproduced exactly below. NR/NC come from
# NR=IFIX(SQRT(FLOAT(NSTND))) computed in REAL*4, so we use Float32 sqrt to match the
# oracle's truncation on/near perfect squares.
# =============================================================================

"""
    ppe_hxindx(nei, istnd, nstnd) -> Int

Bit-exact port of PPBASE HXINDX. `nei` ∈ 1..6 (N/NE/SE/S/SW/NW), `istnd` the 1-based
reference stand, `nstnd` the number of hexagons. Returns the neighbor's 1-based
pointer, or 0 if the neighbor does not exist (off-grid, or any of nei/istnd/nstnd
out of range).
"""
function ppe_hxindx(nei::Integer, istnd::Integer, nstnd::Integer)
    (istnd == 0 || nei < 1 || nei > 6 || nstnd == 0) && return 0
    # NR=IFIX(SQRT(FLOAT(NSTND))) in REAL*4; NC=NR; extra row if NC*NR < NSTND
    nr = trunc(Int, sqrt(Float32(nstnd)))
    nc = nr
    nc * nr < nstnd && (nr += 1)
    # row/col of the reference stand
    ir = istnd ÷ nc
    if mod(istnd, nc) == 0
        ic = nc
    else
        ic = istnd - (ir * nc)
        ir += 1
    end
    local nrow, ncol
    if nei == 1                     # N
        nrow = ir + 1; ncol = ic
    elseif nei == 2                 # NE
        nrow = isodd(ic) ? ir + 1 : ir; ncol = ic + 1
    elseif nei == 3                 # SE
        nrow = isodd(ic) ? ir : ir - 1; ncol = ic + 1
    elseif nei == 4                 # S
        nrow = ir - 1; ncol = ic
    elseif nei == 5                 # SW
        nrow = isodd(ic) ? ir : ir - 1; ncol = ic - 1
    else                            # NW (nei == 6)
        nrow = isodd(ic) ? ir + 1 : ir; ncol = ic - 1
    end
    if ncol <= 0 || nrow <= 0 || ncol > nc || nrow > nr
        return 0
    end
    neiptr = (nrow - 1) * nc + ncol
    neiptr > nstnd && (neiptr = 0)
    return neiptr
end
