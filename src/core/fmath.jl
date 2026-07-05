# fmath.jl — bit-exact-vs-Fortran elementary transcendentals (TOLERANCE_GOAL.md doctrine #8).
#
# Julia's openlibm `exp`/`log`/`pow(**)` for Float32 differ from gfortran's libm by up to
# 1 ULP (measured: exp ~6.3% of inputs, log ~0.11%, pow ~0.17%). FVS is a gfortran build,
# so in the COMPOUNDING growth/crown/volume paths those sub-ULP differences accumulate
# across cycles/tripled records and surface as the "irreducible" tolerance residuals. To
# DECONFOUND them, FVSjl calls the identical gfortran ops (`deps/fvsmath.f90`, built by the
# same compiler as the oracle) through `ccall`. `sqrt` is IEEE-correctly-rounded in both and
# is NOT wrapped. The pure-Julia forms remain as `fexp_julia`/`flog_julia`/`fpow_julia`
# fallbacks (used verbatim when gfortran/the shim is unavailable — the port still RUNS, just
# with openlibm rounding). CAUTION: only wire these into growth/volume/crown math — never the
# RNG (bachlo/rann), whose bit-exactness is achieved with the current ops.
module FMath

using Libdl

const _SO  = joinpath(@__DIR__, "..", "..", "deps", "libfvsmath.$(Libdl.dlext)")
const _SRC = joinpath(@__DIR__, "..", "..", "deps", "fvsmath.f90")

# cached symbol pointers (C_NULL ⇒ shim unavailable ⇒ Julia fallback)
const _EXP = Ref{Ptr{Cvoid}}(C_NULL)
const _LOG = Ref{Ptr{Cvoid}}(C_NULL)
const _POW = Ref{Ptr{Cvoid}}(C_NULL)
const _ACTIVE = Ref{Bool}(false)

"Compile deps/fvsmath.f90 to a shared lib with gfortran if not already present. Returns true on success."
function _ensure_built()
    isfile(_SO) && return true
    gf = Sys.which("gfortran")
    gf === nothing && return false
    try
        run(pipeline(`$gf -shared -fPIC -O2 -o $_SO $_SRC`; stdout = devnull, stderr = devnull))
    catch
        return false
    end
    return isfile(_SO)
end

function __init__()
    try
        _ensure_built() || return
        lib = Libdl.dlopen(_SO)
        _EXP[] = Libdl.dlsym(lib, :f32_exp)
        _LOG[] = Libdl.dlsym(lib, :f32_log)
        _POW[] = Libdl.dlsym(lib, :f32_pow)
        _ACTIVE[] = true
    catch
        _EXP[] = C_NULL; _LOG[] = C_NULL; _POW[] = C_NULL; _ACTIVE[] = false
    end
    return
end

"true when the gfortran companion is loaded and calls route to it (else Julia fallback)."
@inline is_active() = _ACTIVE[]

# pure-Julia fallbacks (openlibm) — kept verbatim, named `_julia`
@inline fexp_julia(x::Float32) = exp(x)
@inline flog_julia(x::Float32) = log(x)
@inline fpow_julia(x::Float32, p::Float32) = x^p

# the DEFAULT ops: gfortran-identical when the shim is active, else openlibm fallback.
@inline function fexp(x::Float32)
    p = _EXP[]
    p == C_NULL ? fexp_julia(x) : ccall(p, Float32, (Float32,), x)
end
@inline function flog(x::Float32)
    p = _LOG[]
    p == C_NULL ? flog_julia(x) : ccall(p, Float32, (Float32,), x)
end
@inline function fpow(x::Float32, q::Float32)
    p = _POW[]
    p == C_NULL ? fpow_julia(x, q) : ccall(p, Float32, (Float32, Float32), x, q)
end

end # module FMath
