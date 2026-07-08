# =============================================================================
# rng.jl — faithful Fortran random number generator state + draws
#
# Ported from: base/rann.f, base/bachlo.f, common/RANCOM.F77, common/ESRNCM.F77
#
# FVS uses a Park–Miller multiplicative LCG (multiplier 16807, modulus 2^31-1).
# There are TWO independent streams:
#   * the main stream (RANN)     — used by growth, mortality, tripling, ...
#   * the establishment stream   — seeded 55329, used only by regeneration height
#     (ESRANN), so that turning establishment on/off doesn't perturb the main draws.
#
# To stay bit-exact we keep the exact same Float64 accumulator + Float32 result
# arithmetic as the Fortran. The whole thing is bundled in one mutable struct so
# each `StandState` owns its own RNG — no global state, threads don't contend.
#
# Requirement #8: this is the faithful default. Swapping to a native Julia RNG
# later only means providing another type with `rann!`/`esrann!` methods.
# =============================================================================

const _RANN_MULT = 16807.0
const _RANN_MOD  = 2147483647.0   # 2^31 - 1
const _RANN_DIV  = 2147483648.0   # 2^31

"""
    FVSRng

Holds the two LCG streams. `s*` is the main stream (RANCOM /S0,S1,SS/), `es*` is
the establishment stream (ESRNCM /ESS0,ESS1,ESSS/). Construct fresh per stand.
"""
mutable struct FVSRng
    s0::Float64        # main stream current state          (RANCOM S0)
    s1::Float64        # main stream scratch                 (RANCOM S1)
    ss::Float32        # main stream saved seed              (RANCOM SS)
    es0::Float64       # establishment stream current state  (ESRNCM ESS0)
    es1::Float64       # establishment stream scratch        (ESRNCM ESS1)
    ess::Float32       # establishment stream saved seed     (ESRNCM ESSS)
end

# ESBLKD block-data defaults: establishment stream seeded to 55329.
FVSRng() = FVSRng(0.0, 0.0, 0.0f0, 55329.0, 0.0, 55329.0f0)

"""
    rann!(r::FVSRng) -> Float32

One draw from the main uniform(0,1) stream. Mirrors `rann.f`: advances S0 and
returns S1/2^31. (Caller variants that pass a Ref read the same value back.)
"""
@inline function rann!(r::FVSRng)::Float32
    r.s1 = _RANN_MULT * r.s0 % _RANN_MOD
    sel = Float32(r.s1 / _RANN_DIV)
    r.s0 = r.s1
    return sel
end

"""
    esrann!(r::FVSRng) -> Float32

One draw from the establishment uniform(0,1) stream (separate seed, `estab.f`).
"""
@inline function esrann!(r::FVSRng)::Float32
    r.es1 = _RANN_MULT * r.es0 % _RANN_MOD
    sel = Float32(r.es1 / _RANN_DIV)
    r.es0 = r.es1
    return sel
end

"""
    ranseed!(r, lset, seed) -> Float32

`RANSED` entry of rann.f. `lset=false` restarts the main stream from the saved
seed `ss`; `lset=true` installs `seed` (forced odd) as the new seed.
"""
function ranseed!(r::FVSRng, lset::Bool, seed::Float32)::Float32
    if !lset
        r.s0 = Float64(r.ss)
        return r.ss
    end
    s = seed
    if s % 2.0f0 == 0.0f0
        s += 1.0f0
    end
    r.ss = s
    r.s0 = Float64(s)
    return s
end

rannget(r::FVSRng)::Float64 = r.s0
rannput!(r::FVSRng, s0::Float64) = (r.s0 = s0; nothing)

"""
    bachlo(r, xbar, stdev; stream=:main) -> Float32

`bachlo.f` — a draw from N(xbar, stdev) via Batchelor's composite rejection
(Tocher 1963), using the requested uniform stream. Establishment height uses the
`:estab` stream, everything else `:main`.
"""
function bachlo(r::FVSRng, xbar::Real, stdev::Real; stream::Symbol = :main)::Float32
    xbar_f  = Float32(xbar)
    stdev_f = Float32(stdev)
    stdev_f <= 0.0f0 && return xbar_f
    draw = stream === :estab ? esrann! : rann!
    while true
        u  = draw(r)
        r1 = draw(r)
        r2 = draw(r)
        local x::Float32, z::Float32
        if u > Float32(2.0 / 3.0)
            zval = 3.0f0 * u - 2.0f0
            zval < 0.001f0 && continue
            x = 1.0f0 - 0.5f0 * log(zval)           # native log kept (doctrine: never route the RNG); PROVEN INERT
            z = 0.5f0 * (x - 2.0f0)^2                # on the DGSCOR-divergent record (id 11 hickory, NOTRIPLE) — flog left it unchanged
        else
            x = 1.5f0 * u
            z = 0.5f0 * x * x
        end
        y = -log(r1)
        y <= z && continue
        r2 >= 0.5f0 && (x = -x)
        return x * stdev_f + xbar_f
    end
end
