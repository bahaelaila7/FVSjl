# fire/fuel_decay.jl — FFE down-wood / litter / duff decay (FMCWD, fire/base/fmcwd.f).
#
# Each cycle the surface-fuel pools `fire.cwd[size 1:11, soft(1)/hard(2), decay-class 1:4]` decay by
# size- and decay-class-specific annual rates (DKR), with a fraction (PRDUFF) of the decayed woody
# material routed to duff, and a hard→soft transfer for the woody classes. Soft material decays 10%
# faster than hard (the ·1.1). This is the *decay* half of the FFE fuel dynamics; the *additions*
# half (litterfall + crown breakage + snag falldown, FMCADD) is a companion routine — and the two are
# COUPLED: litter (size 10) decays at 0.65/yr (≈ gone in one 5-yr cycle), so the forest-floor pool only
# holds up because annual litterfall replenishes it. So this routine is faithful on its own but the
# grown-cycle Stand Carbon Report (DDW/Floor) only validates once FMCADD lands too — see
# docs/FFE_FUEL_DYNAMICS_chunk_plan.md.

# DKR — annual decay rate by [size class 1:11, decay class 1:4] (sn/fmvinit.f:70-104). Woody classes
# 1-9 use 0.11 (decay class 1) or 0.11/0.11/0.09/0.07… (classes 2-4, which copy class 2); litter
# (10) = 0.65, duff (11) = 0.002.
const _FM_DKR = Float32[
    0.11 0.11 0.11 0.11     # 1  (<0.25")
    0.11 0.11 0.11 0.11     # 2  (0.25-1")
    0.11 0.09 0.09 0.09     # 3  (1-3")
    0.11 0.07 0.07 0.07     # 4  (3-6")
    0.11 0.07 0.07 0.07     # 5  (6-12")
    0.11 0.07 0.07 0.07     # 6  (12-20")
    0.11 0.07 0.07 0.07     # 7  (20-35")
    0.11 0.07 0.07 0.07     # 8  (35-50")
    0.11 0.07 0.07 0.07     # 9  (>50")
    0.65 0.65 0.65 0.65     # 10 litter
    0.002 0.002 0.002 0.002 # 11 duff
]
const _FM_PRDUFF = 0.02f0   # proportion of decayed woody material that becomes duff (fmvinit.f:112)

"""
    fmcwd!(s, nyrs) -> StandState

Apply `nyrs` years of FFE surface-fuel decay to `fire.cwd` (FMCWD, fmcwd.f:78-134). Duff is decayed
first (so woody-decay duff additions land after), then each woody/litter class 1-10: a PRDUFF fraction
of the decayed amount is moved to duff, the pool is reduced by `(1−DKR·{1.1 soft})^nyrs`, and (woody
classes < 10) a `nyrs·ln(1−DKR)/ln(0.64)` fraction of the hard pool transfers to soft. No-op unless
FFE is active.
"""
function fmcwd!(s::StandState, nyrs::Integer)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    cwd = fs.cwd; n = Float32(nyrs)
    @inbounds for L in 1:4
        # duff (size 11) first, so woody decay can add to it below
        cwd[11, 1, L] *= (1f0 - _FM_DKR[11, L] * 1.1f0)^n
        cwd[11, 2, L] *= (1f0 - _FM_DKR[11, L])^n
        cwd[11, 1, L] < 0f0 && (cwd[11, 1, L] = 0f0)
        cwd[11, 2, L] < 0f0 && (cwd[11, 2, L] = 0f0)
        for J in 1:10
            dk = _FM_DKR[J, L]
            # amount decayed this cycle → a PRDUFF fraction becomes duff (added to the hard duff pool)
            amt = cwd[J, 1, L] - cwd[J, 1, L] * (1f0 - dk * 1.1f0)^n
            amt < 1f-9 && (amt = 0f0); cwd[11, 2, L] += amt * _FM_PRDUFF
            amt = cwd[J, 2, L] - cwd[J, 2, L] * (1f0 - dk)^n
            amt < 1f-9 && (amt = 0f0); cwd[11, 2, L] += amt * _FM_PRDUFF
            # decrease the pools
            cwd[J, 1, L] *= (1f0 - dk * 1.1f0)^n; cwd[J, 1, L] < 1f-9 && (cwd[J, 1, L] = 0f0)
            cwd[J, 2, L] *= (1f0 - dk)^n;        cwd[J, 2, L] < 1f-9 && (cwd[J, 2, L] = 0f0)
            # hard → soft transfer (woody classes only)
            if J < 10
                tosoft = clamp(n * log(1f0 - dk) / log(0.64f0), 0f0, 1f0) * cwd[J, 2, L]
                cwd[J, 1, L] += tosoft; cwd[J, 2, L] -= tosoft
                cwd[J, 1, L] < 1f-9 && (cwd[J, 1, L] = 0f0)
                cwd[J, 2, L] < 1f-9 && (cwd[J, 2, L] = 0f0)
            end
        end
    end
    return s
end
