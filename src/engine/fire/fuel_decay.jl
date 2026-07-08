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
# Lake States annual decay rates (ls/fmvinit.f:72-95). LS is decay-class-INDEPENDENT (DKR(I,J)=DKR(I,1))
# and its litter loss is 0.31/yr — NOT the SN 0.65 — so LS litter equilibrates ~2× higher, which is the
# down-wood loading FMDYN weights the fire fuel model on. Applying the SN table to LS decayed the litter ~2×
# too fast ⇒ SMALL down-wood ~1.47× low ⇒ FMDYN under-weighted the hot model ⇒ under-scorch ⇒ fire under-kill
# (S82-S86; ls_simfire 2020 TPA). Woody classes are also slower than SN (0.06/0.02 vs 0.07). Duff 0.002 == SN.
const _FM_DKR_LS = Float32[
    0.11 0.11 0.11 0.11     # 1
    0.11 0.11 0.11 0.11     # 2
    0.09 0.09 0.09 0.09     # 3
    0.06 0.06 0.06 0.06     # 4
    0.06 0.06 0.06 0.06     # 5
    0.02 0.02 0.02 0.02     # 6
    0.02 0.02 0.02 0.02     # 7
    0.02 0.02 0.02 0.02     # 8
    0.02 0.02 0.02 0.02     # 9
    0.31 0.31 0.31 0.31     # 10 litter (ls/fmvinit.f:94)
    0.002 0.002 0.002 0.002 # 11 duff
]
# Northeast annual decay rates (ne/fmvinit.f:73-96). Decay-class-INDEPENDENT; litter 0.40/yr (NOT SN 0.65);
# woody 0.19/0.19/0.11/0.07/0.03… (Fahey/Arthur/Foster-Lang). Same SN-mismatch class as LS.
const _FM_DKR_NE = Float32[
    0.19 0.19 0.19 0.19     # 1
    0.19 0.19 0.19 0.19     # 2
    0.11 0.11 0.11 0.11     # 3
    0.07 0.07 0.07 0.07     # 4
    0.03 0.03 0.03 0.03     # 5
    0.03 0.03 0.03 0.03     # 6
    0.03 0.03 0.03 0.03     # 7
    0.03 0.03 0.03 0.03     # 8
    0.03 0.03 0.03 0.03     # 9
    0.40 0.40 0.40 0.40     # 10 litter (ne/fmvinit.f:96)
    0.002 0.002 0.002 0.002 # 11 duff
]
# Central States annual decay rates (cs/fmvinit.f:70-93). Litter 0.65 (== SN) but woody is decay-class-
# INDEPENDENT: classes 3-9 use 0.09/0.07 for ALL decay classes (vs SN's 0.11 at decay class 1). Duff 0.002.
const _FM_DKR_CS = Float32[
    0.11 0.11 0.11 0.11     # 1
    0.11 0.11 0.11 0.11     # 2
    0.09 0.09 0.09 0.09     # 3
    0.07 0.07 0.07 0.07     # 4
    0.07 0.07 0.07 0.07     # 5
    0.07 0.07 0.07 0.07     # 6
    0.07 0.07 0.07 0.07     # 7
    0.07 0.07 0.07 0.07     # 8
    0.07 0.07 0.07 0.07     # 9
    0.65 0.65 0.65 0.65     # 10 litter (cs/fmvinit.f:90, == SN)
    0.002 0.002 0.002 0.002 # 11 duff
]
# Variant-default DKR: SN uses `_FM_DKR`; LS/NE/CS use their own faithful tables (ls/ne/cs fmvinit.f).
_fm_dkr_default(::AbstractVariant) = _FM_DKR
_fm_dkr_default(::LakeStates) = _FM_DKR_LS
_fm_dkr_default(::Northeast) = _FM_DKR_NE
_fm_dkr_default(::CentralStates) = _FM_DKR_CS
const _FM_PRDUFF = 0.02f0   # proportion of decayed woody material that becomes duff (fmvinit.f:112)

"""
    apply_fuelmove!(s) -> Bool

FUELMOVE (act 2530, fmtret.f:203-368): transfer surface fuel between size categories at a scheduled cycle.
Each due activity moves XGET = max(amount, proportion·source, source−leave, target−current) tons/ac (capped
at the available source) from size class FROM to TO; size class 0 is the import (FROM=0) / export (TO=0)
sink. The per-size-class totals are then written back by scaling each class's cwd sub-pools (soft/hard ×
decay) by new/old, or dumping into the hard/fast pool if the class was empty. No-op without a due FUELMOVE.
"""
function apply_fuelmove!(s::StandState)::Bool
    fs = s.fire
    (fs === nothing || !fs.active || isempty(s.control.schedule)) && return false
    yr = Int(current_cycle_year(s)); fvscyc = Int(s.control.cycle) + 1
    cwd = fs.cwd
    # FORG/FSRC/FTRG by size class 0:11 (0 = outside sink); +1 offset so idx 1 = class 0, idx j+1 = class j.
    forg = zeros(Float32, 12)
    @inbounds for j in 1:11; forg[j+1] = sum(@view cwd[j, :, :]); end
    fsrc = copy(forg); ftrg = zeros(Float32, 12); altered = false
    for a in s.control.schedule
        a.icflag == Int32(2530) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        ifrm = Int(round(a.params[1])); ito = Int(round(a.params[2]))
        x = a.params[3]; y = a.params[4]; z = a.params[5]; q = a.params[6]
        (0 <= ifrm <= 11 && 0 <= ito <= 11 && ifrm != ito && x >= 0f0 && 0f0 <= y <= 1f0 && z >= 0f0) || continue
        fi = ifrm + 1; ti = ito + 1
        xget = 0f0
        if ifrm > 0
            fsrc[fi] <= 0f0 && continue
            xget = q >= 0f0 ? max(x, y * fsrc[fi], fsrc[fi] - z, q - fsrc[ti]) : max(x, y * fsrc[fi], fsrc[fi] - z)
            xget > fsrc[fi] && (xget = fsrc[fi])
            fsrc[fi] -= xget
        else
            xget = q >= 0f0 ? max(x, q - fsrc[ti]) : x
        end
        ftrg[ti] += xget
        xget > 0f0 && (altered = true)
    end
    altered || return false
    @inbounds for j in 1:11
        ft = fsrc[j+1] + ftrg[j+1]
        abs(forg[j+1] - ft) >= 1f-6 || continue
        if forg[j+1] <= 1f-6
            cwd[j, 2, 3] = ft                              # empty class → all to hard/fast (CWD(1,J1,2,3))
        else
            sc = ft / forg[j+1]
            for k in 1:2, l in 1:4; cwd[j, k, l] *= sc; end
        end
    end
    return true
end

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
    # FUELMULT/FUELDCAY override the DKR matrix; DUFFPROD overrides PRDUFF — both fall back to defaults
    # when unset (size 0×0).
    dkr = size(fs.params.dkr, 1) == 11 ? fs.params.dkr : _fm_dkr_default(s.variant)
    has_pd = size(fs.params.prduff, 1) == 11; pdm = fs.params.prduff
    @inbounds for L in 1:4
        # duff (size 11) first, so woody decay can add to it below
        cwd[11, 1, L] *= (1f0 - dkr[11, L] * 1.1f0)^n
        cwd[11, 2, L] *= (1f0 - dkr[11, L])^n
        cwd[11, 1, L] < 0f0 && (cwd[11, 1, L] = 0f0)
        cwd[11, 2, L] < 0f0 && (cwd[11, 2, L] = 0f0)
        for J in 1:10
            dk = dkr[J, L]
            pd = has_pd ? pdm[J, L] : _FM_PRDUFF       # DUFFPROD-overridable proportion-to-duff
            # amount decayed this cycle → a PRDUFF fraction becomes duff (added to the hard duff pool)
            amt = cwd[J, 1, L] - cwd[J, 1, L] * (1f0 - dk * 1.1f0)^n
            amt < 1f-9 && (amt = 0f0); cwd[11, 2, L] += amt * pd
            amt = cwd[J, 2, L] - cwd[J, 2, L] * (1f0 - dk)^n
            amt < 1f-9 && (amt = 0f0); cwd[11, 2, L] += amt * pd
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
