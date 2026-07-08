# =============================================================================
# height_growth.jl (lakestates) — LS periodic height growth (ls/htgf.f + htcalc.f + balmod.f)
#
# Structurally IDENTICAL to NE's HTGF: the NC-128 (Carmean/Hahn/Jacobs) height/age
# curve (the SAME shared 127-row LTBHEC table NE uses, indexed per-species by MAPLS —
# htcalc.f IVAR=1), its 10-yr increment, modulated by a competition modifier + relative
# height + the DG serial-correlation random OLDRN, then HTGMULT × SCALE(=FINT/YR) ×
# exp(HTCON). LS differs from NE in ONLY two places:
#   (a) the species→curve map: MAPLS (IVAR=1) instead of MAPNE (IVAR=3);
#   (b) the competition modifier: ls/balmod.f (an OMEGA/BETA/BAMAX1 QMD-relative model),
#       NOT NE's exp(−dg_b3·BAL).
# Data: data/lakestates/htgf_coeffs.csv (MAPLS + the 8 balmod arrays).
#
#   GMOD   = ls_balmod(sp, DBH, BA, RMSQD)                              (ls/balmod.f)
#   RELHTA = AVH>0 ? min(HT/AVH, 1) : 0
#   GMOD   = (1 − (1−GMOD)·(1−RELHTA))·0.8                              (ls/htgf.f:52)
#   HTG    = HTG1·(1+OLDRN)·GMOD ; floor 0.1                            (:53-55)
#   HTG    = SCALE·XHT·HTG·exp(HTCON) ; size-cap to SIZCAP[sp,4]        (:57-63)
# =============================================================================

# htcalc.f curve coefficients (B1,B2,B3,B4,B5,BH) for LS species `sp` via MAPLS into LTBHEC.
@inline _ls_htcoef(mapls, sp::Integer) = @inbounds LTBHEC[Int(mapls[sp])]

# MAPLS (htcalc IVAR=1 species→LTBHEC row) as a module const, so `ls_htcalc_height` (used by the shared
# establishment path) needs no `sd`. Loaded from data/lakestates/htgf_coeffs.csv column 2 (htcalc_mapls).
const _LS_HTCALC_MAP = let
    path = joinpath(LS_DATADIR, "htgf_coeffs.csv")
    m = Int[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        push!(m, parse(Int, split(strip(line), ',')[2]))
    end
    m
end

"NC-128 height (ft) at `age` for LS species `sp` (htcalc.f MODE 1, IVAR=1 via MAPLS) — the ESSUBH base curve."
@inline function ls_htcalc_height(sp::Integer, si::Real, age::Real)
    b1, b2, b3, b4, b5, bh = @inbounds LTBHEC[_LS_HTCALC_MAP[sp]]; sif = Float32(si)
    return bh + b1 * sif^b2 * (1f0 - exp(b3 * Float32(age)))^(b4 * sif^b5)
end

# sp-based htcalc wrappers (LTBHEC via _LS_HTCALC_MAP) — for the shared establishment REGENT(LESTB) path.
@inline ls_htcalc_htmax(sp::Integer, si::Real) = (b = @inbounds LTBHEC[_LS_HTCALC_MAP[sp]]; b[1] * Float32(si)^b[2])
@inline ls_htcalc_age(sp::Integer, si::Real, h::Real)   = _ls_age(@inbounds(LTBHEC[_LS_HTCALC_MAP[sp]]), si, h)
@inline ls_htcalc_incr(sp::Integer, si::Real, aget::Real) = _ls_incr(@inbounds(LTBHEC[_LS_HTCALC_MAP[sp]]), si, aget)
@inline _ls_htmax(coef, si) = coef[1] * Float32(si)^coef[2]
# Tree age from current height (HTCALC mode 0).
function _ls_age(coef, si, h)
    b1, b2, b3, b4, b5, bh = coef; sif = Float32(si)
    base = (Float32(h) - bh) / (b1 * sif^b2)
    base <= 0f0 && return 0f0
    return (1f0 / b3) * log(1f0 - base^(1f0 / (b4 * sif^b5)))
end
# 10-yr height increment from starting age `aget` (HTCALC mode 9).
function _ls_incr(coef, si, aget)
    b1, b2, b3, b4, b5, bh = coef; sif = Float32(si); a = Float32(aget)
    hmax = b1 * sif^b2; ex = b4 * sif^b5
    h0 = bh + hmax * (1f0 - exp(b3 * a))^ex
    hp = bh + hmax * (1f0 - exp(b3 * (a + 10f0)))^ex
    return hp - h0
end

"""
    ls_balmod(sp, d, ba, rmsqd, check, b1, b2, b3, b4, c1, c2, bamax1) -> Float32

LS height-growth competition modifier (ls/balmod.f). `rmsqd` = stand QMD (RMSQD, dense.f:250),
`ba` = stand basal area, `d` = tree DBH. Returns GM ∈ [0.2, 1].
"""
@inline function ls_balmod(sp::Integer, d::Float32, ba::Float32, rmsqd::Float32,
                           check, b1, b2, b3, b4, c1, c2, bamax1)::Float32
    if rmsqd <= 0f0 || d / rmsqd < check[sp]
        omega = b4[sp]
    else
        expval = b2[sp] * d / rmsqd
        expval > 86f0 && (expval = 86f0)
        omega = b1[sp] * (1f0 - exp(-expval))^b3[sp] + b4[sp]
    end
    beta = c1[sp] * (rmsqd + 1f0)^c2[sp]
    batemp = ba > 1f0 ? ba : 1f0                         # ls/balmod.f:107 IF(BA.LE.1.)BATEMP=1.
    arg = bamax1[sp] / batemp - 1f0
    arg < 0f0 && (arg = 0f0)
    gm = 1f0 - exp(-omega * beta * sqrt(arg))
    gm < 0.2f0 && (gm = 0.2f0)
    return gm
end

"""
    height_growth!(state, ::LakeStates; scale)

LS HTGF: periodic height increment into `trees.ht_growth`. NE's NC-128 curve/increment
(via MAPLS) modulated by the LS balmod competition modifier + relative height + OLDRN,
then HTGMULT × SCALE × exp(HTCON), size-capped. `scale = FINT/YR` (YR=10).
"""
function height_growth!(s::StandState, ::LakeStates; scale::Float32 = 1f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    mapls = sd[:htcalc_mapls]
    check = sd[:balmod_check]; b1 = sd[:balmod_b1]; b2 = sd[:balmod_b2]; b3 = sd[:balmod_b3]
    b4 = sd[:balmod_b4]; c1 = sd[:balmod_c1]; c2 = sd[:balmod_c2]; bamax1 = sd[:balmod_bamax1]
    avh = p.avg_height; ba = p.basal_area; rmsqd = stand_qmd(s)   # RMSQD (dense.f:250) — p.qmd is never stored
    oldrn = t.old_random
    cur_year = current_cycle_year(s)
    @inbounds for i in 1:t.n
        sp = t.species[i]; hti = t.height[i]
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        si = p.sp_site_index[sp]
        coef = _ls_htcoef(mapls, sp)
        htmax = _ls_htmax(coef, si)
        xht = active_multiplier(s.control, :htg, sp, cur_year)
        htcon = c.htg_cor[sp]
        # READCORH/REUSCORH: HTCON += ln(HCOR2) when LHCOR2 on (ls/htgf.f:174-176, HTCONS entry).
        s.control.htg_cor2_on && s.control.htg_cor2[sp] > 0f0 && (htcon += log(s.control.htg_cor2[sp]))
        # htcalc.f:389 — within 1 ft of HTMAX (or below BH): HTG1=0, falls to the 0.1 floor.
        if htmax - hti <= 1f0
            htg1 = 0f0
        else
            aget = _ls_age(coef, si, hti)
            htg1 = _ls_incr(coef, si, aget)
        end
        gmod = ls_balmod(sp, t.dbh[i], ba, rmsqd, check, b1, b2, b3, b4, c1, c2, bamax1)
        relht = avh > 0f0 ? min(hti / avh, 1f0) : 0f0
        gmod = (1f0 - (1f0 - gmod) * (1f0 - relht)) * 0.8f0
        htg = htg1 * (1f0 + oldrn[i]) * gmod
        htg < 0.1f0 && (htg = 0.1f0)
        htg = scale * xht * htg * exp(htcon)
        sc4 = s.control.sp_size_cap[sp, 4]
        (hti + htg) > sc4 && (htg = max(sc4 - hti, 0.1f0))
        t.ht_growth[i] = htg
    end
    return s
end
