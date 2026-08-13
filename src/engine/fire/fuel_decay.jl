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
# Central Rockies annual decay rates (cr/fmvinit.f:80-112). Decay class 1: woody 0.12/0.12/0.09 (classes 1-3),
# 0.015 (classes 4-9 = LARGE wood, "agrees with published values for large CWD, Brown et al. 1998"); litter 0.5,
# duff 0.002. Decay classes 2-4 are the class-1 rate × 0.45 (the CR/UT modifier, fmvinit.f:108-112). The SN
# default decayed CR's LARGE down-wood ~7× too fast (0.11 vs 0.015 at decay class 1) ⇒ the down-wood pool
# collapsed over cycles (crt01 large 16.7→1.6) ⇒ FMDYN under-weighted the heavy model on non-MCCT fire stands
# (SFCT/WSCT/ASCT model 8-vs-10) + the crt01 byram residual. This is the CR analogue of the LS DKR fix above.
const _FM_DKR_CR = Float32[
    0.12   0.054   0.054   0.054      # 1  (<0.25")
    0.12   0.054   0.054   0.054      # 2  (0.25-1")
    0.09   0.0405  0.0405  0.0405     # 3  (1-3")
    0.015  0.00675 0.00675 0.00675    # 4  (3-6")   — LARGE CWD 0.015/yr (Brown et al. 1998)
    0.015  0.00675 0.00675 0.00675    # 5  (6-12")
    0.015  0.00675 0.00675 0.00675    # 6  (12-20")
    0.015  0.00675 0.00675 0.00675    # 7  (20-35")
    0.015  0.00675 0.00675 0.00675    # 8  (35-50")
    0.015  0.00675 0.00675 0.00675    # 9  (>50")
    0.5    0.225   0.225   0.225      # 10 litter (cr/fmvinit.f:99)
    0.002  0.0009  0.0009  0.0009     # 11 duff (cr/fmvinit.f:100)
]
# Variant-default DKR: SN uses `_FM_DKR`; LS/NE/CS/CR use their own faithful tables (…/fmvinit.f).
_fm_dkr_default(::AbstractVariant) = _FM_DKR
_fm_dkr_default(::LakeStates) = _FM_DKR_LS
_fm_dkr_default(::Northeast) = _FM_DKR_NE
_fm_dkr_default(::CentralStates) = _FM_DKR_CS
_fm_dkr_default(::CentralRockies) = _FM_DKR_CR
# IE/KT share CR's N-Rockies decay table (ie/fmcwd.f DKR + breakpoints verified identical to cr/fmcwd.f).
# Without this they fell through to the SN `_FM_DKR` default ⇒ wrong 30-yr fuel accumulation ⇒ weak flame.
_fm_dkr_default(::InlandEmpire) = _FM_DKR_CR
_fm_dkr_default(::Kootenai) = _FM_DKR_CR
_fm_dkr_default(::EasternMontana) = _FM_DKR_CR   # em/fmcwd.f DKR verified identical to cr
_fm_dkr_default(::CentralIdaho) = _FM_DKR_CR     # ci/fmcwd.f DKR verified identical to cr
_fm_dkr_default(::Teton) = _FM_DKR_CR            # tt/fmcwd.f DKR == cr
_fm_dkr_default(::Utah) = _FM_DKR_CR             # ut/fmcwd.f DKR == cr
# NC (Klamath) base decay table (nc/fmvinit.f:70-92) — decay-class-INDEPENDENT and MUCH slower than the SN
# default (woody 0.0125-0.025 vs SN 0.07-0.11); subsequently ×DCYMLT (nc/fmcba.f:405, Dunning-code/site index).
# Without this NC fell through to the SN `_FM_DKR` ⇒ LARGE down-wood decayed ~3.7× too fast ⇒ low fuel-model
# weight on the hot model ⇒ weak surface fire. Used as the base that the first-FFE-year DCYMLT then scales.
const _FM_DKR_NC = Float32[
    0.025  0.025  0.025  0.025     # 1  (<0.25")
    0.025  0.025  0.025  0.025     # 2  (0.25-1")
    0.025  0.025  0.025  0.025     # 3  (1-3")
    0.0125 0.0125 0.0125 0.0125    # 4  (3-6")
    0.0125 0.0125 0.0125 0.0125    # 5  (6-12")
    0.0125 0.0125 0.0125 0.0125    # 6  (12-20")
    0.0125 0.0125 0.0125 0.0125    # 7  (20-35")
    0.0125 0.0125 0.0125 0.0125    # 8  (35-50")
    0.0125 0.0125 0.0125 0.0125    # 9  (>50")
    0.5    0.5    0.5    0.5        # 10 litter
    0.002  0.002  0.002  0.002     # 11 duff
]
_fm_dkr_default(::Klamath) = _FM_DKR_NC          # fallback; the first FFE year replaces it with the ×DCYMLT table

# DCYMLT (nc/fmcba.f:395-405 + dunn.f GETDUNN = DUNN50): the site index of the site species → a Dunning
# code (0-7, site↑ ⇒ code↓) → the decay-rate multiplier (Dunning 0-1 → 1.5, 2-4 → 1.0, 5-7 → 0.5).
const _NC_DUNN_XSR = (23f0, 31f0, 39f0, 49f0, 56f0, 75f0, 90f0, 106f0)   # DUNN50, ascending
const _NC_DUNN_YSR = Float32[7, 6, 5, 4, 3, 2, 1, 0]                     # Dunning codes
const _NC_DCY_XD   = (0f0, 1f0, 2f0, 3f0, 4f0, 5f0, 6f0, 7f0)
const _NC_DCY_YD   = Float32[1.5, 1.5, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5]
@inline function nc_dcymlt(si::Float32)::Float32
    dun = _ffe_algslp(si, _NC_DUNN_XSR, _NC_DUNN_YSR)
    return _ffe_algslp(dun, _NC_DCY_XD, _NC_DCY_YD)
end
# NC decay table scaled by DCYMLT (all classes; nct01 has no FUELDCAY ⇒ every SETDECAY<0 ⇒ all scaled).
@inline nc_adjusted_dkr(si::Float32)::Matrix{Float32} = _FM_DKR_NC .* nc_dcymlt(si)

# BM has its OWN base decay table (bm/fmvinit.f:68-113) — NOT the CR table. Faster litter (0.65 vs CR 0.5)
# and different woody rates; used as the base the habitat DKRADJ then scales (see bm_adjusted_dkr).
const _FM_DKR_BM = Float32[
    0.076  0.081  0.090  0.113      # 1  (<0.25")
    0.076  0.081  0.090  0.113      # 2  (0.25-1")
    0.076  0.081  0.090  0.113      # 3  (1-3")
    0.019  0.025  0.033  0.058      # 4  (3-6")
    0.019  0.025  0.033  0.058      # 5  (6-12")
    0.019  0.025  0.033  0.058      # 6  (12-20")
    0.019  0.025  0.033  0.058      # 7  (20-35")
    0.019  0.025  0.033  0.058      # 8  (35-50")
    0.019  0.025  0.033  0.058      # 9  (>50")
    0.65   0.65   0.65   0.65       # 10 litter (bm/fmvinit.f:111)
    0.002  0.002  0.002  0.002      # 11 duff  (bm/fmvinit.f:112)
]
_fm_dkr_default(::BlueMountains) = _FM_DKR_BM     # bm/fmvinit.f — distinct from CR; DKRADJ-scaled at 1st yr
const _FM_PRDUFF = 0.02f0   # proportion of decayed woody material that becomes duff (fmvinit.f:112)

# ── BM decay-rate habitat adjustment (bm/fmcba.f:67-113, 333-368) ──────────────────────────────────
# BM is the lone western variant that DOES NOT use the CR decay table verbatim: at the first FFE year it
# multiplies the base DKR by a habitat-conditioned factor DKRADJ(TEMP,MOIST,K) (from FMR6SDCY). Each BM
# habitat type (ITYPE, 1-92 = the decoded KODTYP / PCOML index) maps to a temperature class BMHMC (1=hot,
# 2=moderate, 3=cold) and a moisture class BMWMD (1=wet, 2=mesic, 3=dry); K is a size-class group
# (1: size 1-3 / <3", 2: size 4-5 / 3-12", 3: size 6-9 / >12"). Without this, BM down-wood decays at the
# generic CR rate and OVER-ACCUMULATES over the projection (measured: bmt01 SMALL/LARGE ~2.4× live), which
# pushes FMDYN onto the hot fuel model → spurious crown fire → mortality over-kill.
const _FM_BMHMC = Int8[   # temperature class by habitat type (bm/fmcba.f:73-84)
    3,3,2,2,2,2,2,2,2,2, 2,2,2,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3,
    2,3,3,2,3,2,2,2,3,3, 3,2,2,2,2,2,3,3,1,1,
    1,1,2,2,2,1,2,2,1,2, 1,1,1,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,3,2,2,2, 2,1]
const _FM_BMWMD = Int8[   # moisture class by habitat type (bm/fmcba.f:88-99)
    3,3,3,3,3,2,3,3,3,3, 3,2,3,2,1,2,3,1,1,2, 1,1,2,2,2,2,2,3,2,3,
    2,3,3,1,1,1,2,1,2,3, 3,3,2,2,2,2,3,3,3,3,
    3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,1,1,2, 2,2,2,1,1,1,3,3,3,2,
    2,2,3,3,2,1,3,2,1,2, 2,2]
# DKRADJ[TEMP, MOIST, K] (bm/fmcba.f:101 — DATA order fills K innermost, then MOIST(J), then TEMP(I)).
const _FM_DKRADJ = let a = Array{Float32}(undef, 3, 3, 3)
    vals = Float32[
        1.7,  2.0,  1.7,   1.49, 1.91, 1.49,  0.75, 0.85, 0.75,    # TEMP=1 (hot):  MOIST 1,2,3
        1.35, 1.85, 1.35,  1.0,  1.7,  1.0,   0.875, 1.2, 0.875,   # TEMP=2 (mod):  MOIST 1,2,3
        1.21, 1.79, 1.21,  1.14, 1.76, 1.14,  0.75, 0.85, 0.75]    # TEMP=3 (cold): MOIST 1,2,3
    n = 0
    for i in 1:3, j in 1:3, k in 1:3
        n += 1; a[i, j, k] = vals[n]
    end
    a
end

"""
    bm_adjusted_dkr(itype) -> Matrix{Float32}

BM habitat-conditioned decay rates (bm/fmcba.f:333-368): the BM base DKR scaled by `DKRADJ(TEMP,MOIST,K)`
for the stand's habitat type `itype`, capped at 1.0, then a second pass (size 9→2) bumps any size class that
would decay slower than the next-larger class up to the larger class's rate. Only the woody classes 1-9 are
adjusted; litter (10) and duff (11) keep the CR default. Applied once, at the first FFE year, when the user
has not overridden decay with FuelDcay.
"""
function bm_adjusted_dkr(itype::Integer)::Matrix{Float32}
    dkr = copy(_FM_DKR_BM)
    (itype < 1 || itype > length(_FM_BMHMC)) && return dkr
    temp = Int(_FM_BMHMC[itype]); moist = Int(_FM_BMWMD[itype])
    @inbounds for i in 1:9
        k = i <= 3 ? 1 : (i <= 5 ? 2 : 3)
        adj = _FM_DKRADJ[temp, moist, k]
        for j in 1:4
            v = dkr[i, j] * adj
            dkr[i, j] = v > 1f0 ? 1f0 : v
        end
    end
    # bump smaller wood up to larger wood's rate where it would otherwise decay more slowly (fmcba.f:359-368)
    @inbounds for i in 9:-1:2, j in 1:4
        (dkr[i, j] - dkr[i-1, j]) > 0f0 && (dkr[i-1, j] = dkr[i, j])
    end
    return dkr
end

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
