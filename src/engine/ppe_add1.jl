# =============================================================================
# PPE (PPBASE) — ADD1: increment one decimal digit of an 11-char internal stand
# number (CISN) at column ICOL (1..7 → character position ICOL+4). Used by PPMAIN's
# internal-stand-number generation. The routine "will not be called if the result
# would exceed 9" (its own comment), so the incremented digit stays a single char.
# Bit-exact port; oracle = gfortran-16 golden. ADDITIVE + INERT.
# =============================================================================

"""
    ppe_add1(cisn::AbstractString, icol::Integer) -> String

Bit-exact port of PPBASE ADD1. Returns `cisn` (11 chars) with the decimal digit at
character position `icol+4` incremented by one. `icol` outside 1..7 returns `cisn`
unchanged (matching the Fortran guard). The digit is assumed 0..8 on entry (result
≤ 9), as the caller guarantees.
"""
function ppe_add1(cisn::AbstractString, icol::Integer)
    s = collect(String(cisn))
    (icol < 1 || icol > 7) && return String(s)
    pos = icol + 4                    # Fortran 1-based char position CISN(ICOL+4:ICOL+4)
    ival = parse(Int, string(s[pos])) # READ (…,'(I1)')
    ival += 1
    s[pos] = Char('0' + ival)         # WRITE (…,'(I1)')  (single digit, guaranteed ≤ 9)
    return String(s)
end
