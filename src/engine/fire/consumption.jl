# =============================================================================
# fire/consumption.jl — fire fuel consumption + carbon release (FFE F7/F8-rest)
#
# Ported from: bin/FVSsn_buildDir/fmcons.f (FMCONS, the natural-unpiled-fuels path).
#
# When a fire burns, it consumes a moisture-dependent fraction of each surface fuel size
# class; the consumed fuel is removed from the down-wood pools (`fire.cwd`) and, converted
# to carbon at 0.5, is the fire carbon *release*. Only the natural-fuels path is ported
# (the activity-fuels / piled-burn variants need the harvest-year context).
# =============================================================================

# PDIA: midpoints of the 6 large (>3") fuel size classes (fmcons.f:57).
const _FM_PDIA = (4.0f0, 8.0f0, 15.0f0, 15.0f0, 15.0f0, 15.0f0)

"""
    fire_consumption_fractions(mois) -> NTuple{11,Float32}

Consumed fraction of each of the 11 surface fuel classes for natural (unpiled) fuels
(FMCONS, fmcons.f:121-189). The <1" classes burn 90%, 1–3" 65%, litter 100%; the >3"
classes use a moisture-driven diameter-reduction `1 − ((PDIA−DIARED)/PDIA)²` and duff a
moisture-linear `(83.7 − 0.426·m_duff%)/100`. `mois` is the fuel-moisture matrix.
"""
# `so`: so/fmcons.f (SO's own FMCONS) never sets PRBURN(1,3) on the no-activity-fuels path — the 1-3" class is not
# consumed (the shared base fmcons.f sets 0.65). MEASURED private FVSso_g16 FMCONS trace, 15364795010497 2010 FMPOFL:
# PRBURN(1,1:11) = .9 .9 0 … — jl's 0.65 put 3.0 t/ac of 1-3" fuel into Pot_Smoke (0.4470 vs live 0.4354).
# `activity`: (IYR−HARVYR) ≤ 5 — a harvest left slash within 5 years (fmscut.f HARVYR) ⇒ fmcons.f:121-170's ACTIVITY,
# UNPILED branch: PRBURN(1,3) from the 0.25-1" moisture (1.0 below 13.7%, 0 above 34%, else (167−4.89·M%)/100 clamped
# [0,1]) and the >3" diameter reduction DIARED = 4.35−0.096·M% (0 when MOIS(1,4)>0.45 or negative). The returned
# PRBURN(1,1:2) are the BURNZ(1,3)>0 values (0.9, or 1.0 when PRBURN(1,3)>0.9); callers re-apply the BURNZ(1,3)=0 rule.
# SO's own fmcons.f activity branch (HFCORR …) is a different formula and is not ported (`so` keeps its natural path).
function fire_consumption_fractions(mois::AbstractMatrix{Float32}; so::Bool = false, activity::Bool = false)::NTuple{11,Float32}
    m100 = mois[1, 4]                                  # 3+" (100-hr+) moisture drives the large classes
    act = activity && !so
    p3 = so ? 0f0 : 0.65f0
    if act
        m2 = mois[1, 2]
        if m2 < 0.137f0
            p3 = 1f0
        elseif m2 > 0.34f0
            p3 = 0f0                                   # NO BURNING OCCURS (PRBURN stays 0)
        else
            cons = (167f0 - 4.89f0 * m2 * 100f0) / 100f0
            p3 = clamp(cons, 0f0, 1f0)
        end
        diared = 4.35f0 - 0.096f0 * m100 * 100f0
        (m100 > 0.45f0 || diared < 0f0) && (diared = 0f0)
    else
        diared = m100 > 1.25f0 ? 0f0 : max(0f0, 3.38f0 - 0.027f0 * m100 * 100f0)
    end
    big = ntuple(i -> (pd = _FM_PDIA[i]; 1f0 - ((pd - diared) / pd)^2), 6)   # classes 4–9
    prduf = min(1f0, max(0f0, 83.7f0 - 0.426f0 * mois[1, 5] * 100f0) / 100f0)
    psmall = p3 > 0.9f0 ? 1f0 : 0.9f0                  # fmcons.f:139-147 (with BURNZ(1,3) > 0)
    return (psmall, psmall, p3, big[1], big[2], big[3], big[4], big[5], big[6], 1.0f0, prduf)
end

# XDIA / XDIA2 (so/fmcons.f:73-74): the ALGSLP knots of SO's activity-fuels 1000-hr diameter reduction (over PRBURN(1,3)
# on its 0-1 scale) and duff reduction (over the fuel-bed DEPTH).
const _SO_XDIA  = Float32[0, 75, 85, 100]
const _SO_XDIA2 = Float32[0, 1, 2, 999]

"""
    fmcons_prburn(fs, mois; so, activity, fmslop, fwind, depth) -> (prburn::NTuple{11,Float32}, prduf::Float32)

FMCONS's unpiled PRBURN(1,1:11) (before the PSBURN scaling) and PRDUF (the duff percent EXPOSR reads), on the BURNZ(1,·)
of the current `fs.cwd`. The base/IE fmcons.f natural and activity branches go through `fire_consumption_fractions` (the
<1" classes 1.0 when BURNZ(1,3)=0, fmcons.f:139-147). SO in the activity case ((IYR−HARVYR)≤5) runs its own branch
(so/fmcons.f:146-225): 1- and 10-hr fuels all consumed; 100-hr PRBURN = 0.9 − (100·M(1,2) − FMCORR − 12)·0.0535 with
FMCORR = 0.03·LOG(MAX(HFCORR,1E-9))/LOG(2), HFCORR = 0.5·BURNZ(1,3)·(1 + (100·FMSLOP−20)/60 + 0.25·FWIND); the ≥1000-hr
DIARED from M(1,4) (ALGSLP over PRBURN(1,3) below 44%), ×0.78/0.89 when M(1,2)≤0.15; class 4 1−((5.2−DIARED)/5.2)²,
classes 5-9 1−((13.7−DIARED)/5.2)²; duff DFRED by M(1,5), ALGSLP'd over DEPTH, PRBURN(1,11) = DFRED·12.1/BURNZ(1,11);
PRDUF = PRBURN(1,11)·100. `fmslop` = FMSLOP (SLOPE, fraction), `fwind` = the common FWIND (midflame wind of the
fire/potential fire), `depth` = the common DEPTH (the last FMGFMV'd fuel model's depth × DPMOD).
"""
function fmcons_prburn(fs, mois::AbstractMatrix{Float32}; so::Bool = false, activity::Bool = false,
                       fmslop::Float32 = 0f0, fwind::Float32 = 0f0, depth::Float32 = 0f0)
    burnz = ntuple(11) do j
        z = 0f0
        @inbounds for k in 1:2, l in 1:4; z += fs.cwd[j, k, l]; end
        z
    end
    if !(so && activity)
        pr0 = fire_consumption_fractions(mois; so = so, activity = activity)
        small = burnz[3] > 0f0 ? (pr0[3] > 0.9f0 ? 1f0 : 0.9f0) : 1f0
        return (ntuple(i -> i <= 2 ? small : pr0[i], 11), max(0f0, 83.7f0 - 0.426f0 * mois[1, 5] * 100f0))
    end
    hfcorr = 0.5f0 * burnz[3] * (1f0 + (100f0 * fmslop - 20f0) / 60f0 + (0.25f0 * fwind))
    fmcorr = 0.03f0 * flog(max(hfcorr, 1f-9)) / flog(2f0)
    p3 = 0.9f0 - ((100f0 * mois[1, 2]) - fmcorr - 12f0) * 0.0535f0
    p3 = max(0f0, min(p3, 1f0))
    m4 = mois[1, 4]
    if m4 > 0.6f0
        diared = -0.005f0 * (m4 * 100f0) + 0.731f0
    elseif m4 >= 0.44f0 && m4 <= 0.6f0
        diared = -0.0178f0 * (m4 * 100f0) + 1.489f0
    else
        y1 = -0.096f0 * (m4 * 100f0) + 4.6495f0
        y3 = -0.125f0 * (m4 * 100f0) + 6.27f0
        diared = algslp(p3, _SO_XDIA, Float32[y1, y1, y3, y3])
    end
    if mois[1, 2] <= 0.15f0
        if m4 <= 0.4f0
            diared = diared * (1f0 - 0.22f0)
        elseif m4 > 0.4f0 && m4 <= 0.5f0
            diared = diared * (1f0 - 0.11f0)
        end
    end
    p4 = max(0f0, min(1f0 - ((5.2f0 - diared) / 5.2f0)^2, 1f0))
    xdia4 = max(0f0, diared)
    p59 = max(0f0, min(1f0 - ((13.7f0 - diared) / 5.2f0)^2, 1f0))
    big = (p4, p59, p59, p59, p59, p59)
    m5 = mois[1, 5]
    if m5 >= 2f0
        xcons = 0f0
        for il in 4:9; xcons = xcons + (big[il - 3] * burnz[il]); end
        dfred = 0.537f0 + (0.57f0 * xcons)
    elseif m5 >= 1.25f0 && m5 < 2f0
        dfred = 0.323f0 + 1.034f0 + sqrt(xdia4)
    elseif m5 >= 0.5f0 && m5 < 1.25f0
        dfred = 1.323f0 + 1.034f0 + sqrt(xdia4)
    else
        dfred = 2.323f0 + 1.034f0 + sqrt(xdia4)
    end
    dfred = algslp(depth, _SO_XDIA2, Float32[dfred * 0.5f0, dfred * 0.5f0, dfred * 0.75f0, dfred * 0.75f0])
    p11 = burnz[11] <= 1f-20 ? 0f0 : max(0f0, min((dfred * 12.1f0) / burnz[11], 1f0))
    return ((1f0, 1f0, p3, big..., 1f0, p11), p11 * 100f0)
end

"""
    apply_fire_consumption!(fs, mois) -> Float32

Consume the surface fuel pools `fs.cwd` by the natural-fuels consumption fractions
(FMCONS) and return the fire carbon release (tons C/acre). Mutates `fs.cwd` (the unburned
remainder stays for the down-wood pools).

The consumed biomass→carbon conversion splits by pool, mirroring the FVS Carbon Report's
"Carbon Released From Fire" (fmcrbout.f:151 `V(11)=BIOCON(1)·0.37 + BIOCON(2)·0.50`, with
fmdout.f:286-287 `BIOCON(1)=BURNED(litter)+BURNED(duff)`, `BIOCON(2)=rest`): the consumed
FOREST-FLOOR (litter=class 10, duff=class 11) converts at 0.37 (Smith & Heath NE-722) and
all consumed WOODY classes (1–9) at 0.50 — the same split as the standing carbon pools.
"""
function apply_fire_consumption!(fs::FireState, mois::AbstractMatrix{Float32}; psburn::Float32 = 100f0)
    return fire_consumption!(fs, mois; psburn).released
end

# FMCONS PM2.5/PM10 emission factors, lb per ton consumed (fmcons.f DATA EMMFAC(IM,IL,IP,IPM), unpiled IP=1 only):
# [IM moisture type 1..3, IL fuel class 1..11, IPM 1=PM2.5 / 2=PM10]; EMFACL(IL,IPM) live herb/shrub/crown.
const _FM_EMMFAC = let e = zeros(Float32, 3, 11, 2)
    for (ipm, fine, mid, big, lit, duff) in ((1, 7.9f0, 11.9f0, (22.5f0, 18.3f0, 16.2f0), 7.9f0, (23.9f0, 25.8f0, 25.8f0)),
                                             (2, 9.3f0, 14.0f0, (26.6f0, 21.6f0, 19.1f0), 9.3f0, (28.2f0, 30.4f0, 30.4f0)))
        for im in 1:3
            e[im, 1, ipm] = fine; e[im, 2, ipm] = fine; e[im, 3, ipm] = mid
            for il in 4:9; e[im, il, ipm] = big[im]; end
            e[im, 10, ipm] = lit; e[im, 11, ipm] = duff[im]
        end
    end
    e
end
const _FM_EMFACL = (21.3f0, 25.1f0)            # EMFACL(1..4, IPM): the same factor for herb, shrub and crown

"""
    fire_consumption!(fs, mois; psburn=100, burncr=0) -> NamedTuple

FMCONS (fmcons.f, ICALL=0, BTYPE=0 natural unpiled fuels): burn `fs.cwd` and return the per-fire consumption
record — `burned` (tons/ac consumed per fuel class 1..11, BURNED(3,·)), `exposr` (% mineral soil exposed),
`burnlv` (live herb/shrub consumed), `smoke` (PM2.5, PM10 tons/ac; SMOKE·P2T as reported by FMFOUT), and the
carbon `released`. The consumed fractions PRBURN and the live PLVBRN are scaled by PSBURN/100 (fmcons.f:196-200);
the <1" classes burn 100% when the 1-3" class is empty (fmcons.f:125-137). `burncr` = FMEFF's BCROWN (crown
material burned, tons/ac) enters the smoke only. `activity` = a harvest ≤5 yr before the fire (IYR−HARVYR≤5, fmscut.f
HARVYR) selects fmcons.f's activity-fuels fractions.
"""
function fire_consumption!(fs::FireState, mois::AbstractMatrix{Float32}; psburn::Float32 = 100f0, burncr::Float32 = 0f0,
                           so::Bool = false, activity::Bool = false, fmslop::Float32 = 0f0, fwind::Float32 = 0f0,
                           depth::Float32 = 0f0)
    pr0, prduf = fmcons_prburn(fs, mois; so = so, activity = activity, fmslop = fmslop, fwind = fwind, depth = depth)
    pr = ntuple(i -> pr0[i] * psburn / 100f0, 11)   # PRBURN(1,I)*PSBURN/100 (fmcons.f:197)
    burned = zeros(Float32, 11)
    @inbounds for isz in 1:11
        f = pr[isz]
        z = 0f0
        for k in 1:2, l in 1:4; z += fs.cwd[isz, k, l]; end
        z > 0f0 || continue
        burned[isz] = z * f
        for k in 1:2, l in 1:4
            fs.cwd[isz, k, l] *= (1f0 - f)
        end
    end
    # EXPOSR (fmcons.f:189-192): mineral soil exposed from the duff consumption percent PRDUF, ×PSBURN/100.
    exposr = prduf < 10f0 ? 0f0 : (-8.98f0 + 0.899f0 * prduf)
    exposr = exposr * psburn / 100f0
    # Live herb/shrub surface fuels also burn (FMCONS BURNLV, fmcons.f:310-311: herb PLVBRN=1.0, shrub 0.6, ×PSBURN)
    # and release at 0.5 (they're in BIOCON(2), fmdout.f:283/287). FLIVE is recomputed each cycle by fmcba!
    # (the live fuels regrow), so this is RELEASE-ONLY — no pool to mutate.
    burnlv = ((1f0 * psburn / 100f0) * fs.flive[1], (0.6f0 * psburn / 100f0) * fs.flive[2])   # PLVBRN·FLIVE
    # smoke (fmcons.f:325-360): IM from the 3+" moisture; dead classes by EMMFAC, live + crown by EMFACL.
    m4 = mois[1, 4]
    im = m4 <= 0.20f0 ? 3 : m4 <= 0.375f0 ? 2 : 1
    smoke = ntuple(2) do ipm
        ts = 0f0
        @inbounds for il in 1:11; ts += burned[il] * _FM_EMMFAC[im, il, ipm]; end
        ts += burnlv[1] * _FM_EMFACL[ipm]              # fmcons.f:345-347 DO IL=1,2: TSMOKE = TSMOKE + PLVBRN·FLIVE·EMFACL,
        ts += burnlv[2] * _FM_EMFACL[ipm]              # one term at a time (was ts + (herb + shrub))
        ts += burncr * _FM_EMFACL[ipm]
        ts * _FM_P2T
    end
    # Fire carbon release (fmdout.f:266-287 → fmcrbout.f:151): TOTCON = ΣBURNED(3,·) + BURNLV + BURNCR (the crown
    # material FMEFF burned is consumed too), BIOCON(1) = litter+duff at 0.37, BIOCON(2) = the rest at 0.50.
    totcon = 0f0
    @inbounds for ii in 1:11; totcon = totcon + burned[ii]; end
    totcon = totcon + burnlv[1] + burnlv[2] + burncr
    biocon1 = burned[10] + burned[11]
    biocon2 = totcon - biocon1
    released = biocon1 * 0.37f0 + biocon2 * 0.50f0
    return (; burned, exposr, burnlv, smoke, released, totcon)
end
