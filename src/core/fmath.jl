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
@inline function fpow_julia(x::Float32, p::Float32)
    # Match C/gfortran powf for a NEGATIVE base: valid only for an integer-valued exponent (sign-preserving,
    # x^p = ±|x|^p by parity), else NaN. Julia's `^` throws DomainError on (negative)^(Float) even when the
    # exponent is integer-valued (e.g. the R4D2H woodland cubic `(a+b·D2H^⅓)**3.` when a+b·c<0 for a tiny tree).
    if x < 0f0
        isinteger(p) || return NaN32
        r = abs(x)^p
        return isodd(round(Int, p)) ? -r : r
    end
    return x^p
end

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

# SIN/COS of REAL*4: gfortran calls glibc sinf/cosf (libm), which can differ from Julia's native Float32
# sin/cos by 1 ULP. Used in the per-species DG/HTG site constants (dgf.f DGCON/SMCON aspect terms).
@inline fsin(x::Float32) = ccall((:sinf, "libm.so.6"), Float32, (Float32,), x)
@inline fcos(x::Float32) = ccall((:cosf, "libm.so.6"), Float32, (Float32,), x)
# TAN/ATAN of REAL*4 likewise (FFE FMCROWE's bole-tip cone angle, fmcrowe.f:381/433)
@inline ftan(x::Float32) = ccall((:tanf, "libm.so.6"), Float32, (Float32,), x)
@inline fatan(x::Float32) = ccall((:atanf, "libm.so.6"), Float32, (Float32,), x)

"""
    fpowi(x::Float32, m::Integer) -> Float32

Fortran `REAL**INTEGER` exactly as gfortran compiles it: a call to libgcc `__powisf2` (libgcc2.c) — binary
square-and-multiply carried out in SINGLE precision (y = x if m odd else 1; repeatedly x=x*x, y=y*x on set bits;
1/y for m<0). Julia's `Float32^Int` instead evaluates in Float64 and rounds once, which differs by 1 ULP on a
large share of inputs; inside a cancellation such as VARMRT's `1-(1-EFFTR)**NPASS` that ULP becomes a ~1e-5
relative error in the kill (measured vs FVSsn, treeszcp_cap cycle 1).
"""
@inline function fpowi(x::Float32, m::Integer)
    n = unsigned(abs(m)); y = isodd(n) ? x : 1f0
    while (n >>= 1) != 0
        x = x * x
        isodd(n) && (y = y * x)
    end
    return m < 0 ? 1f0 / y : y
end

# REAL*8 intrinsics as gfortran emits them: DEXP/DLOG/`**` (REAL*8 operands) are direct calls to glibc libm
# exp/log/pow. Julia's Float64 exp/log/^ are its own implementations and can differ in the last bit, and
# `x^3` lowers to x*x*x (two roundings) where Fortran `X**3.` is one correctly-rounded pow. Used by the
# double-precision NVEL kernels (Flewelling SHP_C2/SHP_OT, SF_TAPER) so their REAL*4 outputs round the same.
const _LIBM = "libm.so.6"
# REAL*4 EXP/ALOG/`**` straight from glibc libm (expf/logf/powf) — what every gfortran build links. THE single binding:
# the named library matters — a bare `ccall(:expf)` resolves to whichever libm Julia loaded first (openlibm under
# Julia 1.13), which differs by 1 ULP (test_lpmpb PROTBK 3C10C78F vs live 3C10C78E; ON Penner DDS). Unlike
# fexp/flog/fpow these never fall back to openlibm when the gfortran companion shim is absent.
@inline expf(x::Float32) = ccall((:expf, _LIBM), Float32, (Float32,), x)
@inline logf(x::Float32) = ccall((:logf, _LIBM), Float32, (Float32,), x)
@inline powf(x::Float32, y::Float32) = ccall((:powf, _LIBM), Float32, (Float32, Float32), x, y)
@inline log10f(x::Float32) = ccall((:log10f, _LIBM), Float32, (Float32,), x)
@inline dexp(x::Float64) = ccall((:exp, _LIBM), Float64, (Float64,), x)
@inline dlog(x::Float64) = ccall((:log, _LIBM), Float64, (Float64,), x)
@inline dpow(x::Float64, y::Float64) = ccall((:pow, _LIBM), Float64, (Float64, Float64), x, y)
@inline dlog10(x::Float64) = ccall((:log10, _LIBM), Float64, (Float64,), x)

end # module FMath
