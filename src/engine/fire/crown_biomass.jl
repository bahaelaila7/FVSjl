# =============================================================================
# fire/crown_biomass.jl — per-tree crown biomass by size class (FFE chunk F2-fn)
#
# Ported from: bin/FVSsn_buildDir/fmcrowe.f (FMCROWE).
#
# Splits a tree's crown into the FFE size classes XV(0:5):
#   XV[0] = foliage
#   XV[1] = 0–0.25"   XV[2] = 0.25–1"   XV[3] = 1–3"   XV[4] = 3–4"+ (incl. bole tip)
#   XV[5] = (unused, 0)
# returned in tons. These feed canopy bulk density / canopy fuels (FMCBA → crown fire).
#
# Method: Jenkins total aboveground biomass → foliage/bark/wood/branch component
# split → allocation into size classes using one of four published proportion forms
# (red oak, shortleaf pine, maple, aspen) chosen by the tree's Lake-States species
# map (ISPMAP → `ls_spi`), realigned by the unmerchantable bole-tip weight (UMBTW).
#
# Dependencies, all already in FVSjl: BRATIO = `bark_ratio` (identical Clark DIB=a+b·D
# + [0.80,0.99] clamp, bar Fort-Bragg-only overrides), HTDBH = `_htdbh_height`,
# FMSVL2 = the SN cubic-volume model via `_R8CLARK_VOL` (`_fm_cuft`), SG = `v2t`,
# P2T = 0.0005 (lb→ton).
# =============================================================================

const _FM_P2T = 0.0005f0           # FMPARM P2T: pounds → tons

# Cumulative bole-tip cone weight (tons) at bark diameter `dbrk`. Hoisted from the per-call `cone`
# closure inside crown_biomass (pillar-2: no per-tree closure allocation). Bit-identical arithmetic —
# same op order, all Float32 (mypi=3.14159f0, _FM_P2T=0.0005f0).
@inline _cone(dbrk::Float32, ang::Float32, sg::Float32, mp::Float32) =
    sg * (dbrk / 2f0 * tan(ang) * dbrk * dbrk * mp / 1728f0) / _FM_P2T

# Same as _cone but the cone is capped at the tree's own bark diameter (td = min(dbrk, d)) — the
# maple/else-branch `conem` closure, hoisted so `angle` is never captured (it was boxed, de-optimizing
# the whole function). Bit-identical arithmetic, all Float32.
@inline _conem(dbrk::Float32, d::Float32, ang::Float32, sg::Float32, mp::Float32) =
    (td = min(dbrk, d); sg * (td / 2f0 * tan(ang) * td * td * mp / 1728f0) / _FM_P2T)

# Standalone total cubic-foot volume of one (species, dbh, height) tree — FMSVL2's
# `TCF` (fmsvol.f), which is just the SN volume model FVSjl uses in compute_volumes!.
# Mirrors the per-tree cubic block of compute_volumes! (volume.jl:337-344).
@inline function _fm_cuft(s::StandState, sp::Integer, d::Float32, h::Float32;
                          merch::Bool = false)::Float32
    s.control.merch_init || init_merch_standards!(s)
    c = s.control
    if d >= c.sp_scf_dbhmin[sp]
        prod = "01"; stump = c.sp_scf_stump[sp]; mtopp = c.sp_scf_topd[sp]
    else
        prod = "02"; stump = c.sp_stump_ht[sp]; mtopp = c.sp_top_diam[sp]
    end
    v, _, _ = _R8CLARK_VOL(s.species.vol_eq[sp], d, h, mtopp, c.sp_top_diam[sp], stump, prod)
    # `merch=true` returns merch cubic v[4] — FMSVL2's MAX(X,MCF) for SN. v[1] is gross (TCF).
    return merch ? v[4] : v[1]
end

# Jenkins total-aboveground (kg) coefficient pair (b0, b1) selected by the tree's
# Lake-States species index (SPILS) — fmcrowe.f:180-214.
@inline function _fm_totabv_coef(spils::Integer)
    if spils == 17 || (40 <= spils <= 42) || (64 <= spils <= 66)
        (-2.2094f0, 2.3867f0)                         # aspen/alder/cottonwood/willow
    elseif spils == 18 || spils == 19 || spils == 24 || spils == 43 || (49 <= spils <= 52)
        (-1.9123f0, 2.3651f0)                         # soft maple/birch
    elseif spils == 26 || spils == 27 || spils == 28 || (30 <= spils <= 39)
        (-2.0127f0, 2.4342f0)                         # hard maple/oak/hickory/beech
    elseif spils == 10 || spils == 11 || spils == 13 || spils == 14
        (-2.0336f0, 2.2592f0)                         # cedar/larch
    elseif spils == 8 || spils == 12
        (-2.5384f0, 2.4814f0)                         # true fir/hemlock
    elseif 1 <= spils <= 5
        (-2.5356f0, 2.4349f0)                         # pine
    elseif spils == 6 || spils == 7 || spils == 9
        (-2.0773f0, 2.3323f0)                         # spruce
    else                                              # mixed hardwood (15,16,20:23,25,29,44:48,53:63,67,68)
        (-2.4800f0, 2.4835f0)
    end
end

"""
    crown_biomass(s, sp, d, h, ic) -> NTuple{6,Float32}

Per-tree crown biomass by FFE size class `XV(0:5)` (FMCROWE, fmcrowe.f) for species
`sp` (SN code), DBH `d` (in), height `h` (ft), crown ratio `ic` (percent). Returns
`(foliage, 0–0.25", 0.25–1", 1–3", 3–4"+, 0)`.

The arithmetic is a faithful transcription of fmcrowe.f, including its quirks: the
bole-tip cone/frustum weights are scaled by `SG/P2T` (≈×2000) while the sub-breast-
height cylinder is added as a RAW volume (no `SG/P2T`) — a known FMCROWE
inconsistency — so the returned values are in FFE-internal units, not literal tons.

⚠ **Not yet end-to-end validated.** Crown biomass feeds canopy bulk density (FMCBA →
crown fire) and does not appear in the `.sum`, so its magnitude can only be confirmed
against live Fortran once the fire-behavior chunks (F5/F6) are in. The function is
included but **not yet called** in the cycle; the structural tests pin the component
split, size-class ordering, and species-form selection, not absolute bit-exactness.
"""
function crown_biomass(s::StandState, sp::Integer, d::Float32, h::Float32, ic::Integer)::NTuple{6,Float32}
    (d == 0f0 || h == 0f0) && return (0f0, 0f0, 0f0, 0f0, 0f0, 0f0)
    coef = s.coef
    spils = Int(coef_col(coef, :ls_spi)[sp])
    sg    = coef_col(coef, :v2t)[sp] * _FM_P2T   # V2T is rescaled /2000 after init (fmvinit.f:1094);
                                                 # the CSV holds the raw V2T, so apply the /2000 here
    dbhmin = coef_col(coef, :dbh_min)[sp]
    ifor  = Int(s.plot.forest_idx)
    cr    = Float32(ic)
    dx, hx = d, h
    # --- Jenkins total aboveground (lb): trees < 1" use the 1" value, scaled back ---
    dd = d < 1f0 ? 1f0 : d
    b0, b1 = _fm_totabv_coef(spils)
    totabv = exp(b0 + b1 * log(dd * 2.54f0)) * 2.2046f0
    # foliage / bark / wood component fractions (fmcrowe.f:218-232), × TOTABV
    if spils >= 15
        fol  = exp(-4.0813f0 + 5.8816f0 / (dd * 2.54f0)) * totabv
        bark = exp(-2.0129f0 - 1.6805f0 / (dd * 2.54f0)) * totabv
        wood = exp(-0.3065f0 - 5.4240f0 / (dd * 2.54f0)) * totabv
    else
        fol  = exp(-2.9584f0 + 4.4766f0 / (dd * 2.54f0)) * totabv
        bark = exp(-2.0980f0 - 1.1432f0 / (dd * 2.54f0)) * totabv
        wood = exp(-0.3737f0 - 1.8055f0 / (dd * 2.54f0)) * totabv
    end
    branch = totabv - (fol + bark + wood)
    if dx < 1f0                                       # small-tree linear scaling
        fol *= dx; branch *= dx                       # (bark/wood/totabv also scaled in Fortran;
    end                                               #  only fol/branch are used downstream)
    branch < 0f0 && (branch = 0f0)
    ttopw = branch
    # --- small unmerch trees (D < DBHMIN): add the whole-bole weight (FMSVL2) ---
    if dx < dbhmin
        dmin = dbhmin
        hmin = _htdbh_height(coef.species, sp, dmin, ifor; isne = s.variant isa Northeast)
        # FVS uses FMSVL2 = MAX(X, MCF) (merch cubic with the tiny-tree cone floor X=0.005454154·H), NOT
        # the gross cuft — gross over-counted the small-tree bole → crown size-2 over (sp33 d1.5-2.2 1.5-2×).
        vt  = max(0.005454154f0 * hmin, _fm_cuft(s, sp, dmin, hmin; merch = true))
        vt1 = 0.0015f0 * dx * dx * hx                 # cone vol of the actual tree
        vt2 = 0.0015f0 * dmin * dmin * hmin           # cone vol of the DBHMIN tree
        vt  = (vt / vt2) * vt1
        ttopw += sg * vt / _FM_P2T
    end
    mypi = 3.14159f0
    # --- size-class proportions: P1/P2/P3 (or maple F1..F4), fmcrowe.f:320-352 ---
    p1 = p2 = p3 = 0f0; f1 = f2 = f3 = f4 = 0f0
    is_maple = spils == 18 || spils == 19 || spils == 26 || spils == 27 || (49 <= spils <= 52)
    if 30 <= spils <= 39                              # red oak / hickory
        p1 =  6.4735f0 * d^(-1.1313f0) * cr^(-0.5777f0)
        p2 = 36.8351f0 * d^(-0.9345f0) * cr^(-0.7014f0)
        p3 = 28.2916f0 * d^(-0.8658f0) * cr^(-0.4084f0)
    elseif 1 <= spils <= 14                           # conifers (shortleaf pine form)
        p1 = 3.525f0 * d^(-0.778f0) * cr^(-0.412f0)
        p2 = 5.989f0 * d^(-0.565f0) * cr^(-0.346f0)
        p3 = 8.585f0 * d^(-0.517f0) * cr^(-0.223f0)
        d <= 1.5f0 && (p1 = 0.5f0); d <= 1.5f0 && (p2 = 1f0)
        (d <= 10.5f0 || cr <= 35f0) && (p3 = 1f0)
    elseif is_maple
        f1 = 1f0 / (4.6762f0 + 0.1091f0 * d^2.0390f0)
        f2 = 1f0 / (3.3212f0 + 0.0777f0 * d^2.0496f0)
        f3 = 1f0 / (0.9341f0 + 0.0158f0 * d^2.1627f0)
        f4 = 1f0 / (0.8625f0 + 0.0093f0 * d^1.7070f0)
        d < 1.9f0 && (f3 = 1f0); d < 4.8f0 && (f4 = 1f0)
    else                                              # aspen (everything else)
        p1 = 1.856f0 * (d * 2.54f0)^(-0.773f0)
        p2 = 5.317f0 * (d * 2.54f0)^(-0.718f0)
        p3 = 1.793f0 * (d * 2.54f0)^(-0.185f0)
    end
    # --- unmerchantable bole-tip weight by size class (UMBTW) + the missing piece (LILPCE) ---
    # u1=0–.25", u2=0–1", u3=0–3", u4=0–4" (cumulative cone/cylinder weight, tons)
    u1 = u2 = u3 = u4 = 0f0
    lilpce = 0f0
    bark_r = bark_ratio(coef, sp, d)
    dobf = 4f0 / bark_r
    if d > dobf && d > dbhmin
        htf = 4.5f0 + (h - 4.5f0) / d * (d - dobf)
        u4 = (h - htf) > 0f0 ? sg * ((h - htf) * 16f0 * mypi / 1728f0) / _FM_P2T : 0f0
        angle = atan((h - htf) / 2f0)
        u1 = _cone(0.25f0, angle, sg, mypi); u2 = _cone(1f0, angle, sg, mypi); u3 = _cone(3f0, angle, sg, mypi)
        dib = 4f0 * bark_r
        htlp = 4.5f0 + (h - 4.5f0) / d * (d - 4f0)
        lilpce = (htlp - htf) > 0f0 ?
            mypi * (htlp - htf) / 1728f0 * (16f0 + 4f0 * dib + dib * dib) * sg / _FM_P2T : 0f0
        lilpce < 0f0 && (lilpce = 0f0)
        u4 -= lilpce
    else
        if h > 4.5f0
            u4 = sg * ((h - 4.5f0) * d * d * mypi / 1728f0) / _FM_P2T
            angle = atan((h - 4.5f0) / (d / 2f0))
            u1 = _conem(0.25f0, d, angle, sg, mypi)   # j=1 always
            d > 0.25f0 && (u2 = _conem(1f0, d, angle, sg, mypi))   # j>1 only if d > DBRK(j-1)
            d > 1f0    && (u3 = _conem(3f0, d, angle, sg, mypi))
        end
        temp = mypi * d * d / 4f0 / 144f0 * min(4.5f0, h)   # cylinder below 4.5 ft
        k = d <= 0.25f0 ? 1 : d <= 1f0 ? 2 : d <= 3f0 ? 3 : 4
        k <= 1 && (u1 += temp); k <= 2 && (u2 += temp); k <= 3 && (u3 += temp); u4 += temp
    end
    # --- clamps (fmcrowe.f:464-492) ---
    ttopw < 0f0 && (ttopw = 0f0); lilpce < 0f0 && (lilpce = 0f0); fol < 0f0 && (fol = 0f0)
    p1 = clamp(p1, 0f0, 1f0); p2 = clamp(p2, 0f0, 1f0); p3 = clamp(p3, 0f0, 1f0)
    p2 < p1 && (p2 = p1); p3 < p2 && (p3 = p2)
    f1 = clamp(f1, 0f0, 1f0); f2 = clamp(f2, 0f0, 1f0); f3 = clamp(f3, 0f0, 1f0); f4 = clamp(f4, 0f0, 1f0)
    f2 < f1 && (f2 = f1); f3 < f2 && (f3 = f2); f4 < f3 && (f4 = f3)
    u1 < 0f0 && (u1 = 0f0); u2 < 0f0 && (u2 = 0f0); u3 < 0f0 && (u3 = 0f0); u4 < 0f0 && (u4 = 0f0)
    u2 < u1 && (u2 = u1); u3 < u2 && (u3 = u2); u4 < u3 && (u4 = u3)
    # --- assemble XV by species form (fmcrowe.f:499-552) ---
    xv0 = fol; xv5 = 0f0
    if is_maple
        ttopw < u4 && (ttopw = u4); ttopw += fol
        body = ttopw - u4
        xv1 = body * (f2 - f1) + u1
        xv2 = body * (f3 - f2) + (u2 - u1)
        xv3 = body * (f4 - f3) + (u3 - u2)
        xv4 = body * (1f0 - f4) + (u4 - u3)
    elseif 1 <= spils <= 14
        ttopw < u4 && (ttopw = u4)
        body = ttopw - u4
        xv1 = body * p1 + u1
        xv2 = body * (p2 - p1) + (u2 - u1)
        xv3 = body * (p3 - p2) + (u3 - u2)
        xv4 = body * (1f0 - p3) + (u4 - u3)
    elseif 30 <= spils <= 39
        ttopw < (u4 - u2) && (ttopw = u4 - u2)
        body = ttopw - u4 + u2
        xv1 = body * p1
        xv2 = body * (p2 - p1)
        xv3 = body * (p3 - p2) + (u3 - u2)
        xv4 = body * (1f0 - p3) + (u4 - u3)
    else                                              # aspen
        ttopw < (u4 - u1) && (ttopw = u4 - u1)
        body = ttopw - u4 + u1
        xv1 = body * p1
        xv2 = body * (p2 - p1) + (u2 - u1)
        xv3 = body * (p3 - p2) + (u3 - u2)
        xv4 = body * (1f0 - p3) + (u4 - u3)
    end
    (d > dobf && d > dbhmin) && (xv4 += lilpce)       # add the LILPCE back for large trees
    return (xv0, xv1, xv2, xv3, xv4, xv5)
end
