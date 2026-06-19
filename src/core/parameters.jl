# =============================================================================
# parameters.jl — compile-time sizing constants
#
# Ported from: common/PRGPRM.F77  (Southern variant values)
#
# These are the fixed array dimensions that FVS uses everywhere. In the original
# Fortran they are PARAMETER statements; here they are `const` so the compiler can
# specialize on them (and so a later `juliac --trim` pass can fold them away).
#
# NOTE on variants: a few of these are variant-specific (MAXSP differs: sn=90,
# base=23). For the Southern-only scope they are plain consts. When CS/NE/LS are
# added these will move onto the variant trait (see docs/ARCHITECTURE.md) — until
# then, keeping them as consts is simplest and trim-friendliest.
# =============================================================================

const MAXTRE  = 3000          # max tree records per stand          (PRGPRM)
const MAXTP1  = MAXTRE + 1    # MAXTRE + 1 (slot for the "new" tree)
const MAXPLT  = 500           # max individual plots / point ids
const MAXSP   = 90            # max species (Southern variant)
const MAXCYC  = 40            # max projection cycles
const MAXCY1  = MAXCYC + 1    # cycle-boundary years array length
const MAXSTR  = 20            # max site trees
const MXFRCDS = 20            # max forest codes

# Fortran NINT (round half away from zero) — Julia's `round` is banker's rounding,
# so every NINT in FVS must use this mode to stay bit-exact.
@inline nint(x::Real)::Int32 = Int32(round(x, RoundNearestTiesAway))
