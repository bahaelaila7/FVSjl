# =============================================================================
# volume.jl — per-tree volume driver (VOLS / CFVOL, the SN R8 Clark path)
#
# Ported from: base/vols.jl (VOLS) + base/fvsvol.jl (CFVOL).
#
# For each live tree, select the merchantability product (sawtimber when
# DBH ≥ SCFMIND, else pulpwood), call the pure `_R8CLARK_VOL` taper model, and
# load the four per-tree volumes the .sum reports:
#   cuft_vol  (CFV)  = vol[1]                     — total cubic
#   merch_cuft_vol (MCFV) = vol[4]+vol[7] if D≥DBHMIN — merch cubic
#   saw_cuft_vol  (SCFV) = vol[4]      if D≥SCFMIND   — sawtimber cubic
#   bdft_vol  (BFV)  = vol[10]                    — board feet
# (snt01 has no defect; defect correction + the two-pass dead handling come later.)
# =============================================================================

# Curtis-Arney H-D parameters for `sp`, with the Fort Bragg (IFOR=20, htdbh.f:145)
# longleaf/loblolly overrides applied. `ifor` defaults to 0 (no override).
@inline function _htdbh_params(sd, sp::Integer, ifor::Integer)
    p2 = sd[:htdbh_p2][sp]; p3 = sd[:htdbh_p3][sp]
    p4 = sd[:htdbh_p4][sp]; db = sd[:htdbh_db][sp]
    if ifor == 20
        sp == 6  && (p2 = 110.3f0; p3 = 7.0670f0; p4 = -1.0420f0)
        sp == 8  && (p2 = 114.6f0; p3 = 4.1840f0; p4 = -0.6940f0)
        sp == 11 && (p2 = 623.9f0; p3 = 4.7396f0; p4 = -0.2763f0)
        sp == 13 && (p2 = 184.3f0; p3 = 4.2660f0; p4 = -0.5496f0)
    end
    return p2, p3, p4, db
end

# NE htdbh has a per-species choice (IWYKCA) of WYKOFF vs CURTIS-ARNEY (htdbh.f:452):
# Wykoff H = exp(HT1 + HT2/(D+1)) + 4.5. SN is Curtis-Arney only (no :htdbh_iwykca
# column → `_uses_wykoff` is false), so this is inert for SN.
@inline _uses_wykoff(sd, sp::Integer) =
    haskey(sd, :htdbh_iwykca) && sd[:htdbh_iwykca][sp] == 0f0

# Wykoff HT-DBH intercept/slope (HT1/HT2) for `sp`, applying the Allegheny NF (IFOR=3)
# Tech-Note-6 (Hough) overrides for 20 hardwood species (sitset.f:428-489). All 20 are
# Wykoff species (IWYKCA=0), so this only ever fires inside the `_uses_wykoff` path —
# i.e. NE-only; SN never reaches it (no :htdbh_iwykca column). `ifor` defaults to 0.
# The IFOR==3 block is the NE Allegheny NF Tech-Note-6 (Hough) override, keyed by NE species indices
# (sp26=RM, 30/31/33=Y/S/P-birch, 41=ash, 55=WO …). It MUST be gated to Northeast — CS (Hoosier) and LS
# also use forest index 3, where the SAME sp index is a DIFFERENT species (CS sp41=YP, LS sp41=QA), so
# applying these NE coefs there was wrong (LS aspen dubbed 64.2 ft vs live 52.2; CS elm — the earlier #99
# birch mis-group). `isne` gates it; non-NE variants fall through to their own base Wykoff coefs.
@inline function _htdbh_wykoff(sd, sp::Integer, ifor::Integer, isne::Bool = false)
    ht1 = sd[:htdbh_ht1][sp]; ht2 = sd[:htdbh_ht2][sp]
    if ifor == 3 && isne
        sp == 26  && (ht1 = 4.6839f0;  ht2 = -4.9622f0)   # red maple
        sp == 27  && (ht1 = 4.6354f0;  ht2 = -4.7168f0)   # sugar maple
        (sp == 30 || sp == 31 || sp == 33) && (ht1 = 4.4635f0; ht2 = -3.6456f0)  # yellow/sweet/paper birch
        sp == 40  && (ht1 = 4.5497f0;  ht2 = -4.6727f0)   # american beech
        (sp == 41 || sp == 42 || sp == 44) && (ht1 = 4.6804f0; ht2 = -4.5561f0)  # ash sp/white/green (=white ash)
        sp == 54  && (ht1 = 4.7614f0;  ht2 = -5.3776f0)   # black cherry
        (sp == 55 || sp == 60 || sp == 64 || sp == 67 || sp == 69) &&
            (ht1 = 4.9100f0; ht2 = -7.2941f0)             # white/scarlet/chestnut/N.red/black oak (=white oak)
        (sp == 71 || sp == 102 || sp == 106) && (ht1 = 4.4393f0; ht2 = -4.0711f0) # other-hw/serviceberry (=hophornbeam)
        sp == 93  && (ht1 = 4.6855f0;  ht2 = -4.8690f0)   # american basswood
        sp == 108 && (ht1 = 4.7614f0;  ht2 = -5.3776f0)   # pin cherry (=black cherry)
    end
    return ht1, ht2
end

"HTDBH mode-0 predicted total height (ft) for `sp` at DBH `d` (htdbh.f)."
@inline function _htdbh_height(sd, sp::Integer, d::Float32, ifor::Integer = 0; isne::Bool = false)
    if _uses_wykoff(sd, sp)
        ht1, ht2 = _htdbh_wykoff(sd, sp, ifor, isne)
        return exp(ht1 + ht2 / (d + 1f0)) + 4.5f0
    end
    p2, p3, p4, db = _htdbh_params(sd, sp, ifor)
    if d >= 3f0
        return 4.5f0 + p2 * exp(-p3 * d ^ p4)
    else
        hat3 = 4.5f0 + p2 * exp(-p3 * 3f0 ^ p4)
        return (hat3 - 4.51f0) * (d - db) / (3f0 - db) + 4.51f0
    end
end

"""
HTDBH mode-1 inverse: dbh (in) from total height `h` (ft) for `sp` (htdbh.f kode 1).

`db_floor=true` applies the per-species DB (budwidth) floor `IF(MODE.NE.0 .AND. D.LT.DB)D=DB`
(NE htdbh.f:475 / CS:436 / LS) — WITHOUT it a Wykoff species just above breast height returns a
NEGATIVE dbh (e.g. LS sugar maple at H≈4.9 → −0.155), which the small-tree DG path then mis-handles
as a fallback. SN's htdbh.f has NO such floor, so SN keeps the default `db_floor=false`.
"""
@inline function _htdbh_dbh(sd, sp::Integer, h::Float32, ifor::Integer = 0; db_floor::Bool = false, isne::Bool = false)
    p2, p3, p4, db = _htdbh_params(sd, sp, ifor)   # db (budwidth) is a per-species array, valid for Wykoff species too
    d = if _uses_wykoff(sd, sp)
        ht1, ht2 = _htdbh_wykoff(sd, sp, ifor, isne)
        ht2 / (log(h - 4.5f0) - ht1) - 1f0          # htdbh.f:463 (:331)
    else
        hat3 = 4.5f0 + p2 * exp(-p3 * 3f0 ^ p4)
        if h >= hat3
            ratio = (log(min(h, 4.5f0 + p2 * 0.9999f0) - 4.5f0) - log(p2)) / (-p3)
            ratio > 0f0 ? exp(log(ratio) * (1f0 / p4)) : 100f0
        else
            ((h - 4.51f0) * (3f0 - db) / (hat3 - 4.51f0)) + db
        end
    end
    return (db_floor && d < db) ? db : d            # htdbh.f:343 IF(MODE.NE.0 .AND. D.LT.DB)D=DB
end

# (The duplicate volume-side bark_ratio + bark_coeffs.csv were removed: bark is now a
# single per-stand source, calib.bark_a/bark_b — see bark_and_bounds.jl / dgcons!.)

# DBH-class breakpoints for the segmented defect curves (vols.f DBHCLS).
const _DBHCLS = (0f0, 5f0, 10f0, 15f0, 20f0, 25f0, 30f0, 35f0, 40f0)

# ALGSLP (algslp.f): segmented-linear interpolation of defect column `sp` of `m` over the
# DBH breakpoints `x` at `xx`, flat-extrapolated beyond the ends. Reads the matrix column
# directly (no SubArray) to keep the volume loop allocation-free.
@inline function _algslp_col(xx::Float32, x::NTuple{9,Float32}, m::Matrix{Float32}, sp::Integer)
    @inbounds begin
        xx < x[1] && return m[1, sp]
        xx >= x[9] && return m[9, sp]
        for i in 1:8
            xx < x[i+1] && return m[i, sp] + (m[i+1, sp] - m[i, sp]) / (x[i+1] - x[i]) * (xx - x[i])
        end
        return m[9, sp]
    end
end

# Behre hyperbola taper (behprm.f / BEHRE) used to redistribute volume after a
# broken/killed top. `behre_params` returns the (AHAT,BHAT) hyperbola constants
# plus a cone flag; `behre` integrates the relative profile between two heights.
@inline function behre_params(vmax::Float32, d::Float32, h::Float32, bark::Float32)
    # Match behprm.f's EXACT left-to-right Float32 multiply order `.00545415*D*D*BARK*BARK*H`
    # (NOT d^2*bark^2 — `(d*d)*(bark*bark)` groups differently and drifts BHAT by 1 ULP, which the
    # 1/BHAT + log + sqrt in AHAT amplify into the broken-top cuft ±1). behprm is broken-top-only,
    # so this is inert for all normal trees.
    bhat = vmax / (0.00545415f0 * d * d * bark * bark * h)
    bhat > 0.95f0 && (bhat = 0.95f0)
    ahat = 0.44277f0 - 0.99167f0 / bhat - 1.43237f0 * flog(bhat) +
           1.68581f0 * sqrt(bhat) - 0.13611f0 * bhat * bhat   # `.13611*BHAT*BHAT` left-to-right; gfortran log (doctrine #8)
    lcone = false
    if abs(ahat) < 0.05f0
        lcone = true
        ahat = ahat < 0f0 ? -0.05f0 : 0.05f0
    end
    bhat = 1f0 - ahat
    bhat < 0.0001f0 && (bhat = 0.0001f0)
    return ahat, bhat, lcone
end

@inline function behre(ahat::Float32, bhat::Float32, l1::Float32, l2::Float32)
    alb1 = ahat * l1 + bhat
    alb2 = ahat * l2 + bhat
    return alb2 - alb1 - 2f0 * bhat * (flog(alb2) - flog(alb1)) -   # gfortran-identical log (doctrine #8)
           bhat * bhat / alb2 + bhat * bhat / alb1
end

"""
    cftopk(merch, sp, d, h, tcf, mcf, scf, vmax, bark, itht) -> (tcf, mcf, scf)

CFTOPK (cftopk.f): reduce total/merch/sawtimber cubic for a broken top at height
`itht/100` ft, using the Behre taper fit to the full-height tree. `merch` carries
the per-stand cubic merch standards (stmp/topd/scfstmp/scftop), so VOLUME keyword
overrides take effect. Pure.
"""
function cftopk(merch, sp::Integer, d::Float32, h::Float32,
                tcf::Float32, mcf::Float32, scf::Float32,
                vmax::Float32, bark::Float32, itht::Integer)
    ahat, bhat, lcone = behre_params(vmax, d, h, bark)
    pht = 0f0; dtrunc = 0f0
    if tcf > 0f0
        volt = behre(ahat, bhat, 0f0, 1f0)
        pht = 1f0 - (Float32(itht) / 100f0) / h
        pht < 0f0 && (pht = 0f0)
        dtrunc = pht / (ahat * pht + bhat)
        if !lcone
            tcf = tcf * behre(ahat, bhat, pht, 1f0) / volt
        else
            tcf = tcf * (1f0 - pht^3)
        end
    end
    if mcf > 0f0
        stump = 1f0 - merch.stmp[sp] / h
        dmrch = merch.topd[sp] / d
        htmrch = (bhat * dmrch) / (1f0 - ahat * dmrch)
        if !lcone
            if dtrunc > dmrch
                mcf = mcf * behre(ahat, bhat, pht, stump) / behre(ahat, bhat, htmrch, stump)
            end
        else
            s3 = stump^3
            dtrunc > dmrch && (mcf = mcf * (s3 - pht^3) / (s3 - htmrch^3))
        end
        mcf > tcf && (mcf = tcf); mcf < 0f0 && (mcf = 0f0)
        if scf > 0f0
            stump = 1f0 - merch.scfstmp[sp] / h
            dmrch = merch.scftop[sp] / d
            htmrch = (bhat * dmrch) / (1f0 - ahat * dmrch)
            if !lcone
                if dtrunc > dmrch
                    scf = scf * behre(ahat, bhat, pht, stump) / behre(ahat, bhat, htmrch, stump)
                end
            else
                s3 = stump^3
                dtrunc > dmrch && (scf = scf * (s3 - pht^3) / (s3 - htmrch^3))
            end
            scf > mcf && (scf = mcf); scf < 0f0 && (scf = 0f0)
        end
    end
    return tcf, mcf, scf
end

"""
    bftopk(merch, sp, d, h, bbfv, vmax, bark, itht) -> bbfv

BFTOPK (bftopk.f): reduce board-foot volume for a broken top at `itht/100` ft,
using the Behre taper and the per-stand board-foot merch limits (BFTOPD/BFSTMP,
overridable by BFVOLUME). Pure.
"""
function bftopk(merch, sp::Integer, d::Float32, h::Float32, bbfv::Float32,
                vmax::Float32, bark::Float32, itht::Integer)
    bbfv <= 0f0 && return bbfv
    ahat, bhat, lcone = behre_params(vmax, d, h, bark)
    pht = 1f0 - (Float32(itht) / 100f0) / h
    dtrunc = pht / (ahat * pht + bhat)
    bftopd = merch.bftopd[sp]
    if dtrunc > bftopd / d
        htmrch = (bhat * bftopd / d) / (1f0 - ahat * bftopd / d)
        stump = 1f0 - merch.bfstmp[sp] / h
        if lcone
            bbfv = bbfv * (stump^3 - pht^3) / (stump^3 - htmrch^3)
        else
            bbfv = bbfv * behre(ahat, bhat, pht, stump) / behre(ahat, bhat, htmrch, stump)
        end
    end
    return bbfv
end

"""
    r4_topkill(t, i, sp, d, h, bark, tcf, mcf, bf, merch) -> (tcf, mcf, bf)

Broken/killed-top volume reduction for the Region-4 western MAT/FW2 paths (fvsvol.f OCFVOL/OBFVOL entries
→ CFTOPK/BFTOPK). A top-killed tree (ITRUNC = `t.trunc[i]` > 0, H ≥ 4.5) has its FULL-height cubic + board
volume trimmed to the standing break height `ITRUNC/100` ft via the Behre taper. `bark` = start-of-cycle
BRATIO (uses the stashed `t.vol_bark[i]` when present). No-op for un-killed trees (the common path). The
region-4 DVE woodland path does NOT call this (fvsvol.f skips CFTOPK for DVE).
"""
# Region-4 cftopk/bftopk use the RAW grinit merch top TOPD=6.0 (fvsvol.f CFTOPK: DMRCH=TOPD(ISPC)/D), NOT the
# species-CSV `top_dib` (which is 4.0 for TT / 0.0 for CI — the secondary-product default, wrong for the broken-
# top reduction). A 6.0-filled top vector keyed by species for the western MAT/FW2 broken-top path.
const _R4_TOPD6 = fill(6.0f0, 64)
const _BM_TOPD45 = fill(4.5f0, 64)                                       # BM grinit TOPD=4.5

@inline function r4_topkill(t, i::Integer, sp::Integer, d::Float32, h::Float32, bark::Float32,
                            tcf::Float32, mcf::Float32, bf::Float32, merch, topdv = _R4_TOPD6)
    (t.trunc[i] > 0 && tcf > 0f0 && h >= 4.5f0) || return (tcf, mcf, bf)
    merch = (stmp = merch.stmp, topd = topdv, scfstmp = merch.scfstmp,
             scftop = topdv, bftopd = topdv, bfstmp = merch.bfstmp)  # grinit TOPD (6.0 R4 / 4.5 BM), not CSV top_dib
    bk = t.vol_bark[i] > 0f0 ? t.vol_bark[i] : bark
    vmx = tcf
    tcf, mcf, _ = cftopk(merch, sp, d, h, tcf, mcf, 0f0, vmx, bk, Int(t.trunc[i]))
    bf = bftopk(merch, sp, d, h, bf, vmx, bk, Int(t.trunc[i]))
    return (max(tcf, 0f0), max(mcf, 0f0), max(bf, 0f0))
end

"""
    dub_missing_heights!(state)

CRATET height resolution (cratet.f:212-265): assign heights to trees missing one
and resolve the full ("normal") height of broken-top trees. Missing-height live
trees get the HTDBH curve height. Topkill trees (norm_ht<0) keep their broken
height but get `norm_ht` = full predicted height ×100 (≥ the standing height), and
a break point `trunc` (80% of standing height when none was supplied).
"""
function dub_missing_heights!(s::StandState)
    t = s.trees; sd = s.coef.species; ifor = Int(s.plot.forest_idx)
    isne = s.variant isa Northeast     # NE-only Allegheny (IFOR=3) HT-DBH overrides (variant-safe gate)
    iscr_dub = s.variant isa CentralRockies   # CR Black Hills (IMODTY 3) no-AA height dub (cratet.f:352)
    # NOHTDREG/LHTDRG (cratet.f:292-335): for each invoked species, fit the Wykoff HT-DBH INTERCEPT from its
    # measured-height trees — `AA = mean(log(H−4.5) − HT2/(D+1))` over trees with H>4.5, NORMHT≥0, D≥3; if ≥3 such
    # trees and AA≥0, set IABFLG=0 so the dub below uses the calibrated Wykoff curve instead of Curtis-Arney.
    # Gated on any LHTDRG species ⇒ fully inert for the default stand (the common path is untouched).
    lhtdrg = s.control.ht_drag_sp; aa = s.calib.ht_dbh_aa; iabflg = s.calib.ht_dbh_iabflg
    # `:wykoff_ht2` is the SN-only NOHTDREG/LHTDRG calibration column (its Wykoff HT2 intercept).
    # Read it ONLY when some species invoked LHTDRG — otherwise it is never used (the calibrated-Wykoff
    # branch below is gated on lhtdrg[sp]), and reading it would KeyError on a variant (NE) that has no
    # such column. The per-tree dub itself uses the variant-generic `_htdbh_height` (htdbh_* coefs).
    # IE's cratet AA-fit uses its blkdat Wykoff HT-DBH HT2 (`:ht2`); `:wykoff_ht2` is IE's separate SPROUT
    # column (≠ blkdat HT2) ⇒ using it gave AA 4.512 vs live 4.2112. Other variants keep `:wykoff_ht2`.
    ht2 = any(lhtdrg) ? coef_col(s.coef, (s.variant isa InlandEmpire || s.variant isa Utah) ? :ht2 : :wykoff_ht2) : nothing
    # TT height-dubbing (tt/cratet.f CASE DEFAULT) uses its OWN Wykoff HT-DBH: H=exp(AX+HT2/(D+1))+4.5,
    # AX=AA(calibrated,IABFLG==0) else HT1(default); PP(sp10,D≤3) linear special. NOT the shared Curtis-Arney
    # `_htdbh_height` (TT defines no htdbh_p2/p3/p4). Load HT1/wykoff_ht2 unconditionally for TT.
    tt_ht1  = s.variant isa Teton ? coef_col(s.coef, :ht1) : nothing
    tt_wht2 = s.variant isa Teton ? coef_col(s.coef, :wykoff_ht2) : nothing
    if any(lhtdrg)
        nmax = length(lhtdrg)
        # FVS accumulates SUMX in REAL (Float32) (cratet.f:292-305); match the dtype.
        sumx = zeros(Float32, nmax); k1 = zeros(Int, nmax)
        @inbounds for i in 1:t.n
            sp = Int(t.species[i]); (1 <= sp <= nmax && lhtdrg[sp]) || continue
            h = t.height[i]; d = t.dbh[i]
            (h > 4.5f0 && t.norm_ht[i] >= 0 && d >= 3f0) || continue       # measured, sound, ≥3" (cratet.f:301)
            sumx[sp] += log(h - 4.5f0) - ht2[sp] / (d + 1f0)               # REAL (Float32), as FVS
            k1[sp]   += 1
        end
        @inbounds for sp in 1:nmax
            (lhtdrg[sp] && k1[sp] >= 3) || continue
            a = sumx[sp] / Float32(k1[sp]); aa[sp] = a
            a >= 0f0 && (iabflg[sp] = Int32(0))                            # IABFLG=0 ⇒ calibrated Wykoff
        end
    end
    # cratet.f dubs missing heights for the LIVE trees (DO loop @337) AND, identically, for the DEAD records
    # (DO 145 @417, II=IREC2..MAXTRE — same AA/HTDBH formula + top-kill handling). jl stores the dead block at
    # t.n+1 : t.n+t.ndead, so dub over BOTH partitions. The AA fit above stays live-only (FVS fits AA from live
    # measured trees, DO 15). The dead-tree heights don't enter the live .sum aggregate but DO feed the DG-
    # calibration backdating (which exposes the dead partition), so dubbing them keeps that calibration faithful.
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; sp = t.species[i]
        tkill = t.norm_ht[i] < 0
        if t.height[i] > 0f0 && !tkill
            continue
        end
        # cratet.f:342-372: calibrated-Wykoff dub when LHTDRG[sp] & IABFLG==0, else the Curtis-Arney HTDBH dub.
        h_v = if d <= 0.1f0
            1.01f0
        elseif lhtdrg[sp] && iabflg[sp] == 0
            exp(aa[sp] + ht2[sp] / (d + 1f0)) + 4.5f0
        elseif iscr_dub && Int(s.plot.model_type) == 3 && lhtdrg[sp] && iabflg[sp] == 1
            # cratet.f:352-360 — Black Hills (IMODTY 3, no AA fit): a distinct SI-driven logistic height
            # curve, NOT the Curtis-Arney/Wykoff dub. Without it, all-missing-height Black Hills ponderosa
            # dub ~8-10 ft low (⇒ TopHt + volume low).
            si = s.plot.sp_site_index[sp]
            if d > 0.5f0
                32.108633f0 * fpow(si, 0.276926f0) *
                    fpow(1f0 - fexp(-0.057766f0 * d), fpow(0.9844340f0, -0.169876f0)) + 4.5f0
            else
                (12.41173f0 + 0.04633f0 * si - 0.000158f0 * si * si) * d
            end
        elseif s.variant isa Teton
            # tt/cratet.f CASE DEFAULT dub: PP(sp10,D≤3) linear; else Wykoff exp(AX+BX/(D+1))+4.5,
            # AX=AA (calibrated) else HT1 (default), BX=wykoff_ht2. (CASE 13,16 P2/P3/P4 forms deferred.)
            if sp == 10 && d <= 3.0f0
                1.74189f0 + 4.17687f0 * d
            else
                ax = (lhtdrg[sp] && iabflg[sp] == 0) ? aa[sp] : tt_ht1[sp]
                exp(ax + tt_wht2[sp] / (d + 1f0)) + 4.5f0
            end
        elseif s.variant isa SoutheastAlaska
            ak_htdbh_dub(Int(sp), d)   # ak/cratet.f Curtis-Arney INVENTORY-EQN dub (LHTDRG=false)
        elseif s.variant isa WestCascades
            # wc/cratet.f:375-377 — WC LHTDRG=.FALSE. for all species ⇒ ALL missing-height dubbing
            # uses the FOREST-DEPENDENT Curtis-Arney HTDBH (MODE=0), not the shared single-table _htdbh_height.
            wc_htdbh_height(_wc_htdbh_ifor(Int(s.plot.forest_idx)), Int(sp), d)
        elseif s.variant isa PacificNorthwest
            pn_htdbh_height(_pn_htdbh_ifor(Int(s.plot.forest_idx)), Int(sp), d)
        elseif s.variant isa EastCascades
            # ec/cratet.f LHTDRG=.FALSE. all species ⇒ forest-dependent Curtis-Arney HTDBH (MODE=0).
            ec_htdbh_height(ec_htdbh_ifor(s.plot), Int(sp), d)
        elseif s.variant isa OregonCoast
            # oc/cratet.f:679-684 — OC LHTDRG=.FALSE. all species ⇒ HTDBH (CA-family Curtis-Arney)
            # is the actual missing-height dub (overwrites the Wykoff H). IFOR unused (C8).
            oc_htdbh_height(Int(sp), d)
        elseif s.variant isa Olympic
            # op/cratet.f:692-694 — OP LHTDRG=.FALSE. all species ⇒ HTDBH MODE=0 (forest-dependent
            # Curtis-Arney, op/htdbh.f == pn/htdbh.f). ORGANON trees are dubbed earlier in PREPARE
            # (op_organon_prepare!); only the FVS-native missing-height records reach here.
            op_htdbh_height(Int(s.plot.forest_idx), Int(sp), d)
        elseif s.variant isa CentralCalifornia
            # ca/cratet.f LHTDRG=.FALSE. all species ⇒ HTDBH (Curtis-Arney) missing-height dub. IFOR unused.
            ca_htdbh_height(Int(sp), d)
        elseif s.variant isa SouthCentralOregon
            # so/cratet.f LHTDRG=.FALSE. all species ⇒ forest-dependent Curtis HTDBH (MODE=0).
            so_htdbh_height(Int(s.plot.forest_idx), Int(sp), d)
        elseif s.variant isa WestSierra
            # ws/cratet.f:447-451 — MIXED LHTDRG (unlike the other westside variants): the NATIVE conifers
            # (LHTDRG=.TRUE.) with a calibrated AA are handled by the top calibrated-Wykoff branch above; every
            # other WS species reaching HERE (LHTDRG=.FALSE. surrogate, OR native with IABFLG==1) has its Wykoff
            # value OVERRIDDEN by HTDBH MODE=0 (`.NOT.LHTDRG .OR. (LHTDRG.AND.IABFLG==1)` ⇒ always TRUE here),
            # so the HTDBH MODE=0 curve IS the dub. (Native-uncalibrated → P2=0 ⇒ H=4.5, faithful to cratet.)
            ws_htdbh_height(0, Int(sp), d)
        elseif s.variant isa Ontario
            # on/cratet.f LHTDRG=.FALSE. all species ⇒ HTDBH MODE=0 (Wykoff/Curtis-Arney) default-coefficient
            # dub. Reaches here for missing / broken-top heights (which drop the measured H).
            _on_htdbh_height(Int(sp), d)
        else
            _htdbh_height(sd, sp, d, ifor; isne = isne)
        end
        # cratet.f:352-354: the D≤0.1 seedling case sets H=1.01 and `GO TO 115`, SKIPPING the 4.5 floor —
        # so a seedling keeps 1.01 ft (sub-breast-height). Only the Curtis-Arney/Wykoff branches (D>0.1) are
        # floored at 4.5. Without the `d>0.1` guard the floor clobbers 1.01→4.5, inflating dubbed seedling
        # heights and the AVH top-height on seedling-heavy stands (live-confirmed: DBH 0.1 → HT 1.01).
        d > 0.1f0 && h_v < 4.5f0 && (h_v = 4.5f0)
        if !tkill
            t.height[i] = h_v
        else
            # cratet.f:381-397: NORMHT/ITRUNC use Fortran INT() = truncate-toward-zero (round-half-UP via +0.5),
            # NOT Julia round() (round-half-to-EVEN) — they diverge by 1 when x is an odd integer.
            t.norm_ht[i] = trunc(Int32, h_v * 100f0 + 0.5f0)
            if t.trunc[i] == 0
                if t.height[i] > 0f0
                    t.trunc[i] = trunc(Int32, 80f0 * t.height[i] + 0.5f0)
                else
                    t.trunc[i] = trunc(Int32, 80f0 * h_v + 0.5f0)
                    t.height[i] = h_v
                end
            else
                if t.height[i] > 0f0
                    t.height[i] < Float32(t.trunc[i]) * 0.01f0 &&
                        (t.height[i] = Float32(t.trunc[i]) * 0.01f0)
                else
                    t.height[i] = Float32(t.trunc[i]) * 0.01f0
                end
            end
            Float32(t.norm_ht[i]) * 0.01f0 < t.height[i] &&
                (t.norm_ht[i] = trunc(Int32, t.height[i] * 100f0))
        end
    end
    return s
end

"""
    init_merch_standards!(state)

Copy the variant's default merch standards (merch_specs.csv) into the per-stand
`Control.sp_*` arrays once, at LSTART. VOLUME / BFVOLUME then overwrite these
per-stand copies (never the shared coefficient tables), and the volume routines
read the copies — that is what makes the override keywords take effect. The copies
are bit-identical to the coef defaults, so an un-overridden stand is unchanged.
"""
function init_merch_standards!(s::StandState)
    s.control.merch_init && return s
    c = s.control
    if s.variant isa Northeast || s.variant isa CentralStates || s.variant isa LakeStates ||
       s.variant isa Southern
        # All four are IFOR-dependent merch code rules (SN setcubicdflts.f, NE/CS/LS sitset.f, via
        # `_sn_merch`/`_ne_merch`/`_cs_merch`/`_ls_merch`) rather than static data. Board-foot mins equal
        # the sawtimber cubic mins (bf-equal) in all four. (SN's merch_specs.csv held exactly the non-NC
        # defaults `_sn_merch` reproduces; routing it here adds the North Carolina IFOR=11 overrides.)
        sn = s.variant isa Southern
        cs = s.variant isa CentralStates
        ls = s.variant isa LakeStates
        # SN uses the resolved IFOR (plot.forest_idx) when set — e.g. STDINFO maps Fort Bragg forest 701 to
        # IFOR=20 even though its KODFOR is remapped to Uwharrie 81110 for VOLEQDEF, and live keys merch on
        # that IFOR, not the remapped code. Only FIA stands leave forest_idx=0, where IFOR+KODIST are decoded
        # from KODFOR the FVS way (IFORST=KODFOR/100−IREGN*100, KODIST=KODFOR mod 100; sitset.f:369/forkod.f:470).
        # KODIST is consulted only when IFOR==11 (North Carolina). NE/CS/LS carry forest_idx directly.
        kodist = 0
        ifor = if sn
            kod = Int(s.plot.user_forest_code)
            kodist = kod % 100
            f = Int(s.plot.forest_idx)
            f != 0 ? f : (kod == 0 ? _SN_DEFAULT_IFOR : (kod ÷ 100 - (kod ÷ 10000) * 100))
        else
            f = Int(s.plot.forest_idx); f == 0 ? (cs ? 1 : (ls ? 5 : _NE_DEFAULT_IFOR)) : f
        end
        @inbounds for j in 1:length(c.sp_dbh_min)
            dbhmin, topd, scfmind, scftopd, stmp, scfstmp =
                sn ? _sn_merch(j, ifor, kodist) : ls ? _ls_merch(j, ifor) :
                cs ? _cs_merch(j, ifor) : _ne_merch(j, ifor)
            c.sp_dbh_min[j] = dbhmin; c.sp_top_diam[j] = topd
            c.sp_scf_dbhmin[j] = scfmind; c.sp_scf_topd[j] = scftopd
            c.sp_stump_ht[j] = stmp; c.sp_scf_stump[j] = scfstmp
            c.sp_bf_dbhmin[j] = scfmind; c.sp_bf_topd[j] = scftopd; c.sp_bf_stump[j] = scfstmp
        end
        c.merch_init = true
        return s
    end
    if s.variant isa WestCascades || s.variant isa PacificNorthwest
        # wc/sitset.f westside merch defaults (IFOR 6 Willamette = CASE DEFAULT): TOPD=BFTOPD=SCFTOPD=4.5,
        # DBHMIN=BFMIND=SCFMIND=7 (LP sp-index 11 = 6), stump=1. WC's species CSV carries no merch columns.
        @inbounds for j in 1:length(c.sp_dbh_min)
            dm = j == 11 ? 6.0f0 : 7.0f0
            c.sp_dbh_min[j] = dm; c.sp_top_diam[j] = 4.5f0; c.sp_stump_ht[j] = 1.0f0
            c.sp_scf_dbhmin[j] = dm; c.sp_scf_topd[j] = 4.5f0; c.sp_scf_stump[j] = 1.0f0
            c.sp_bf_dbhmin[j] = dm; c.sp_bf_topd[j] = 4.5f0; c.sp_bf_stump[j] = 1.0f0
        end
        c.merch_init = true
        return s
    end
    if s.variant isa EastCascades
        # ec/sitset.f westside merch defaults: TOPD=BFTOPD=SCFTOPD=4.5, DBHMIN=7 (LP sp-index 7 = 6), stump=1.
        @inbounds for j in 1:length(c.sp_dbh_min)
            dm = j == 7 ? 6.0f0 : 7.0f0
            c.sp_dbh_min[j] = dm; c.sp_top_diam[j] = 4.5f0; c.sp_stump_ht[j] = 1.0f0
            c.sp_scf_dbhmin[j] = dm; c.sp_scf_topd[j] = 4.5f0; c.sp_scf_stump[j] = 1.0f0
            c.sp_bf_dbhmin[j] = dm; c.sp_bf_topd[j] = 4.5f0; c.sp_bf_stump[j] = 1.0f0
        end
        c.merch_init = true
        return s
    end
    if s.variant isa OregonCoast
        # oc/grinit.f:88-131 DBHMIN=7 (sp-index 11 = 6); oc/sitset.f:241-247 westside (IFOR 6-10)
        # TOPD=BFTOPD=SCFTOPD=4.5; stump=1. OC's species CSV carries no merch columns (like WC).
        @inbounds for j in 1:length(c.sp_dbh_min)
            dm = j == 11 ? 6.0f0 : 7.0f0
            c.sp_dbh_min[j] = dm; c.sp_top_diam[j] = 4.5f0; c.sp_stump_ht[j] = 1.0f0
            c.sp_scf_dbhmin[j] = dm; c.sp_scf_topd[j] = 4.5f0; c.sp_scf_stump[j] = 1.0f0
            c.sp_bf_dbhmin[j] = dm; c.sp_bf_topd[j] = 4.5f0; c.sp_bf_stump[j] = 1.0f0
        end
        c.merch_init = true
        return s
    end
    if s.variant isa Olympic
        # op/sitset.f:228-250 — IFOR 4,5,6 (BLM Salem/Eugene/Coos Bay) ⇒ TOPD=BFTOPD=SCFTOPD=5.0,
        # DBHMIN=BFMIND=SCFMIND=7 (ALL species, no LP special); NF forests (IFOR 1,2,3 = CASE DEFAULT)
        # ⇒ TOPD=4.5, LP(sp-index 11)=6. stump=1. OP's species CSV carries no merch columns.
         blm = 4 <= Int(s.plot.forest_idx) <= 6
        topd = blm ? 5.0f0 : 4.5f0
        @inbounds for j in 1:length(c.sp_dbh_min)
            dm = (!blm && j == 11) ? 6.0f0 : 7.0f0
            c.sp_dbh_min[j] = dm; c.sp_top_diam[j] = topd; c.sp_stump_ht[j] = 1.0f0
            c.sp_scf_dbhmin[j] = dm; c.sp_scf_topd[j] = topd; c.sp_scf_stump[j] = 1.0f0
            c.sp_bf_dbhmin[j] = dm; c.sp_bf_topd[j] = topd; c.sp_bf_stump[j] = 1.0f0
        end
        c.merch_init = true
        return s
    end
    if s.variant isa SouthCentralOregon
        # so/grinit.f DBHMIN=BFMIND=SCFMIND=9.0 ALL species (no sp-11 special, unlike CA);
        # so/sitset.f:167 TOPD=BFTOPD=SCFTOPD = (IFOR∈{1,2,3,10} ? 4.5 : 6.0); stump=1. No merch CSV columns.
        topd = (Int(s.plot.forest_idx) <= 3 || Int(s.plot.forest_idx) == 10) ? 4.5f0 : 6.0f0
        @inbounds for j in 1:length(c.sp_dbh_min)
            c.sp_dbh_min[j] = 9.0f0; c.sp_top_diam[j] = topd; c.sp_stump_ht[j] = 1.0f0
            c.sp_scf_dbhmin[j] = 9.0f0; c.sp_scf_topd[j] = topd; c.sp_scf_stump[j] = 1.0f0
            c.sp_bf_dbhmin[j] = 9.0f0; c.sp_bf_topd[j] = topd; c.sp_bf_stump[j] = 1.0f0
        end
        c.merch_init = true
        return s
    end
    if s.variant isa WestSierra
        # ws/grinit.f:159-170 — UNIFORM all species (no forest/sp-11 special): cubic DBHMIN=7 TOPD=4.5;
        # board BFMIND=10 BFTOPD=6.0; scribner-cubic SCFMIND=10 SCFTOPD=6.0; all stumps=1. No merch CSV columns
        # (WS computes volume via the VEQNNC WO2W/DVEW kernels). Needed by the FFE fire biomass path (_fm_cuft).
        @inbounds for j in 1:length(c.sp_dbh_min)
            c.sp_dbh_min[j]    = 7.0f0; c.sp_top_diam[j]  = 4.5f0; c.sp_stump_ht[j] = 1.0f0
            c.sp_scf_dbhmin[j] = 10.0f0; c.sp_scf_topd[j] = 6.0f0; c.sp_scf_stump[j] = 1.0f0
            c.sp_bf_dbhmin[j]  = 10.0f0; c.sp_bf_topd[j]  = 6.0f0; c.sp_bf_stump[j]  = 1.0f0
        end
        c.merch_init = true
        return s
    end
    if s.variant isa CentralCalifornia
        # ca/grinit.f:85-134 DBHMIN=BFMIND=7 (sp-index 11 = 6); ca/sitset.f:225-238 top diameter is
        # FOREST-dependent — IFOR 6-10 (R6) → 4.5, else (R5) → 6.0; stump=1. No merch CSV columns.
        topd = (6 <= Int(s.plot.forest_idx) <= 10) ? 4.5f0 : 6.0f0
        @inbounds for j in 1:length(c.sp_dbh_min)
            dm = j == 11 ? 6.0f0 : 7.0f0
            c.sp_dbh_min[j] = dm; c.sp_top_diam[j] = topd; c.sp_stump_ht[j] = 1.0f0
            c.sp_scf_dbhmin[j] = dm; c.sp_scf_topd[j] = topd; c.sp_scf_stump[j] = 1.0f0
            c.sp_bf_dbhmin[j] = dm; c.sp_bf_topd[j] = topd; c.sp_bf_stump[j] = 1.0f0
        end
        c.merch_init = true
        return s
    end
    if s.variant isa Ontario
        # canada/on grinit.f + sitset.f merch defaults. grinit seeds TOPD=BFTOPD=10cm·CMtoIN,
        # STMP=BFSTMP=30cm·CMtoFT, DBHMIN=BFMIND=0 (ONMTD is DATA MAXSP*0.0). sitset then fills the
        # zero DBHMIN/BFMIND (TOPD/BFTOPD stay 10cm since they are >0): softwoods (ISPC≤14 or >68)
        # DBHMIN=5, BFMIND=9, BFTOPD=7.6; hardwoods key on IFOR (SELECT CASE, DBHMIN 5/6, BFMIND 9/11,
        # BFTOPD 7.6/9.6) — for ont01 IFOR=9 = CASE DEFAULT → DBHMIN=5, BFMIND=11, BFTOPD=9.6. ON has
        # NO Scribner-cubic (SCF) merch standard, so mirror the cubic values into the scf_* slots
        # (unused until the ON volume kernel — htont/varvol/cubrds/nbolt — lands as a later chunk).
        cmToIn = 0.3937f0; cmToFt = 0.0328084f0
        topd = 10.0f0 * cmToIn                       # 3.937"  (grinit, sitset leaves >0 untouched)
        stmp = 30.0f0 * cmToFt                       # 0.984252 ft (cubic == board stump)
        ifor = Int(s.plot.forest_idx)
        @inbounds for j in 1:length(c.sp_dbh_min)
            sw = (j <= 14 || j > 68)                 # sitset softwood test (ISPC.LE.14 .OR .GT.68)
            hw4042 = (40 <= j <= 42)
            dbhmin = sw ? 5.0f0 :
                     ifor == 2 ? (hw4042 ? 6.0f0 : 5.0f0) :
                     ifor == 6 ? 6.0f0 : 5.0f0
            bfmind = sw ? 9.0f0 :
                     ifor == 2 ? (hw4042 ? 11.0f0 : 9.0f0) :
                     ifor == 5 ? (hw4042 ? 9.0f0 : 11.0f0) : 11.0f0
            bftopd = sw ? 7.6f0 :
                     ifor == 2 ? (hw4042 ? 9.6f0 : 7.6f0) :
                     ifor == 5 ? 7.6f0 : 9.6f0
            c.sp_dbh_min[j]    = dbhmin; c.sp_top_diam[j] = topd; c.sp_stump_ht[j] = stmp
            c.sp_bf_dbhmin[j]  = bfmind; c.sp_bf_topd[j]  = bftopd; c.sp_bf_stump[j] = stmp
            c.sp_scf_dbhmin[j] = dbhmin; c.sp_scf_topd[j] = topd; c.sp_scf_stump[j] = stmp
        end
        c.merch_init = true
        return s
    end
    sd = s.coef.species
    @inbounds for j in 1:length(c.sp_dbh_min)
        c.sp_scf_dbhmin[j] = sd[:scf_min_dbh][j]
        c.sp_scf_topd[j]   = sd[:scf_top_dib][j]
        c.sp_top_diam[j]   = sd[:top_dib][j]
        c.sp_stump_ht[j]   = sd[:stump][j]
        c.sp_scf_stump[j]  = sd[:scf_stump][j]
        c.sp_dbh_min[j]    = sd[:dbh_min][j]
        c.sp_bf_dbhmin[j]  = sd[:bf_min_dbh][j]
        c.sp_bf_topd[j]    = sd[:bf_top_dib][j]
        c.sp_bf_stump[j]   = sd[:bf_stump][j]
    end
    c.merch_init = true
    return s
end

"Assign `val` to the merch array for species `isp`: 0=all, >0=that species, <0=SPGROUP −N."
@inline function _set_merch_sp!(c, arr::Vector{Float32}, isp::Integer, val::Float32)
    if isp == 0
        @inbounds for j in 1:length(arr); arr[j] = val; end
    elseif isp > 0
        isp <= length(arr) && (arr[isp] = val)
    else
        g = -isp
        (1 <= g <= length(c.sp_groups)) || return
        @inbounds for sp in c.sp_groups[g]; arr[sp] = val; end
    end
    return
end

"""
    apply_volume_overrides!(state; fint)

Apply any VOLUME / BFVOLUME merch-standard override whose date has been reached
(volkey.f). Overrides overwrite the per-stand `Control.sp_*` arrays and persist;
re-applying the same event in a later cycle is an idempotent overwrite. Called
only inside `grow_cycle!` — like Fortran VOLKEY's `ICYC.EQ.0` skip, the cycle-0
inventory volume keeps the variant defaults.
"""
function apply_volume_overrides!(s::StandState; fint::Float32 = 5f0)
    isempty(s.control.volume_events) && return s
    s.control.merch_init || init_merch_standards!(s)
    c = s.control
    # OPCYCL containing-cycle bucketing (opcycl.f:58-64, same as the SIMFIRE _fire_due gate): a VOLUME/
    # BFVOLUME/MC-BFDEFECT at date D takes effect in the cycle with IY(i) ≤ D < IY(i+1) and stays applied
    # after — so a mid-cycle date (e.g. 1995 in a 10-yr NE cycle 1990→2000) applies THIS cycle, not one late.
    # The condition is `D < cycle_END`: for D in [cs,ce) it's this (containing) cycle, for D<cs a past cycle
    # already applied. For a BOUNDARY date (SN 5-yr, D=cycle start) this is identical to the old `D ≤ cs` —
    # 1995 ∈ [1995,2000) either way — so SN stays bit-exact; only mid-cycle dates shift one cycle earlier.
    cyc = Int(c.cycle); cs = cycle_year_at(c, cyc)
    ce = cycle_year_at(c, cyc + 1); ce <= cs && (ce = cs + 1)   # degenerate/last cycle ⇒ exact-match
    for ev in c.volume_events
        Int(ev.year) < ce || continue
        isp = round(Int, ev.params[1])
        if ev.icflag == Int32(217)        # VOLUME — cubic merch standards
            _set_merch_sp!(c, c.sp_dbh_min,    isp, ev.params[2])
            _set_merch_sp!(c, c.sp_top_diam,   isp, ev.params[3])
            _set_merch_sp!(c, c.sp_stump_ht,   isp, ev.params[4])
            _set_merch_sp!(c, c.sp_scf_dbhmin, isp, ev.params[5])
            _set_merch_sp!(c, c.sp_scf_topd,   isp, ev.params[6])
            _set_merch_sp!(c, c.sp_scf_stump,  isp, ev.aux)
        elseif ev.icflag == Int32(218)    # BFVOLUME — board-foot merch standards
            _set_merch_sp!(c, c.sp_bf_dbhmin,  isp, ev.params[2])
            _set_merch_sp!(c, c.sp_bf_topd,    isp, ev.params[3])
            _set_merch_sp!(c, c.sp_bf_stump,   isp, ev.params[4])
        elseif ev.icflag == Int32(215)    # MCDEFECT — dated cubic defect curve (sdefet.f)
            _set_defect!(c, c.sp_cf_defect, isp,
                         (ev.params[2], ev.params[3], ev.params[4], ev.params[5], ev.params[6]))
        elseif ev.icflag == Int32(216)    # BFDEFECT — dated board-foot defect curve
            _set_defect!(c, c.sp_bf_defect, isp,
                         (ev.params[2], ev.params[3], ev.params[4], ev.params[5], ev.params[6]))
        end
    end
    return s
end

"""
    compute_volumes!(state)

Fill `trees.{cuft_vol,merch_cuft_vol,saw_cuft_vol,bdft_vol}` for every live tree
from the R8 Clark taper model and the per-stand merch standards (Control.sp_*,
overridable by VOLUME/BFVOLUME). Needs `setup_volume_equations!` to have set
`species.vol_eq`.
"""
# Shared per-tree volume-defect correction (FVS vols.f:285-432, the same ICDF/IBDF block for SN Clark AND
# the NE/CS R9 path — vols.f is the variant-agnostic driver). Cuts the pulpwood part (mcf−scf) by ICDF%
# and board feet + sawtimber cubic by IBDF%, each the MAX of: the per-tree DEFECT input (packed decimal),
# the MCDEFECT/BFDEFECT DBH curve (ALGSLP), and the MCFDLN/BFFDLN log-linear form model. Returns the
# corrected (mcf, scf, bf) with mcf = pulpv + post-board-defect scf. NINT throughout (ties away from zero).
@inline function _apply_tree_defect(mcf::Float32, scf::Float32, bf::Float32, d::Float32, sp::Integer,
                                    dpack::Integer, cfdef, bfdef, cff0, cff1, bff0, bff1,
                                    anydef_cf::Bool, anydef_bf::Bool)
    icdf = dpack ÷ 1000000
    (anydef_cf && mcf > scf) &&     # NINT (vols.f:13,21), not Julia ties-to-even
        (icdf = max(icdf, clamp(round(Int, _algslp_col(d, _DBHCLS, cfdef, sp) * 100f0, RoundNearestTiesAway), 0, 99)))
    temvol = mcf - scf
    if temvol > 0f0 && (cff0[sp] != 0f0 || cff1[sp] != 1f0)
        volcor = exp(cff0[sp] + cff1[sp] * log(temvol))
        icdf = max(icdf, round(Int, (temvol - volcor) / temvol * 100f0, RoundNearestTiesAway))
    end
    icdf = clamp(icdf, 0, 99)
    pulpv = icdf >= 99 ? 0f0 : (mcf - scf) * (1f0 - icdf * 0.01f0)
    # vols.f:352,415-420: the INPUT board-defect is applied to BFV *and* SCFV even when BFV=0 — a
    # too-small-for-boardfeet tree (BFV=0, SCFV>0) still loses sawtimber cubic to its input BF defect.
    # ONLY the curve/form IBDF updates (BFDEFT, log-linear) are gated on BFV>0 (vols.f:393); NINT throughout.
    ibdf = (dpack ÷ 10000) % 100
    if bf > 0f0
        anydef_bf &&
            (ibdf = max(ibdf, clamp(round(Int, _algslp_col(d, _DBHCLS, bfdef, sp) * 100f0, RoundNearestTiesAway), 0, 99)))
        if bff0[sp] != 0f0 || bff1[sp] != 1f0
            volcorb = exp(bff0[sp] + bff1[sp] * log(bf))
            ibdf = max(ibdf, round(Int, (bf - volcorb) / bf * 100f0, RoundNearestTiesAway))
        end
    end
    ibdf = clamp(ibdf, 0, 99)
    if ibdf >= 99
        bf = 0f0; scf = 0f0
    elseif ibdf > 0
        f = 1f0 - ibdf * 0.01f0; bf *= f; scf *= f
    end
    return pulpv + scf, scf, bf
end

function compute_volumes!(s::StandState)
    # Eastern variants (NE + CS + LS) share the NVEL Region-9 Clark cubic + R9LOGS board path,
    # differing only in the IFOR merch standards (_ne_merch / _cs_merch / _ls_merch, dispatched inside)
    # and the per-species METHC=5 DVEE/Gevorkiantz opt-in. LS lst01 defaults METHC=6 ⇒ pure Clark.
    (s.variant isa Northeast || s.variant isa CentralStates || s.variant isa LakeStates) &&
        return compute_volumes_ne!(s)
    s.variant isa CentralRockies && return compute_volumes_cr!(s)
    s.variant isa Kootenai && return compute_volumes_kt!(s)
    s.variant isa EasternMontana && return compute_volumes_em!(s)   # EM = FW2 conifers + DVEW (R1KEMP) non-conifers
    s.variant isa Teton && return compute_volumes_tt!(s)           # TT = Region-4 Matney (r4vol) cubic
    s.variant isa Utah && return compute_volumes_ut!(s)            # UT = MATW r4vol + FW2 + DVEW woodland
    s.variant isa BlueMountains && return compute_volumes_bm!(s)   # BM = FW2W Flewelling conifers (BEHW minor deferred)
    s.variant isa WestCascades && return compute_volumes_wc!(s)   # WC = westside Flewelling (SHP_W3/W4/W5) + INGY FW2 + region-6 Behre
    s.variant isa PacificNorthwest && return compute_volumes_pn!(s)   # PN = westside DF F00 + region-6 Behre (near-clone of WC)
    s.variant isa EastCascades && return compute_volumes_ec!(s)   # EC = forest-8 INGY I11/I12 Flewelling + region-6 Behre
    s.variant isa CentralCalifornia && return compute_volumes_ca!(s)  # CA = R6 Behre + FW2 Flewelling (F06/I00), ca/formcl.f
    s.variant isa CentralIdaho && return compute_volumes_ci!(s)   # CI = MATW r4vol + FW2W Flewelling + DVEW woodland (= UT)
    s.variant isa Klamath && return compute_volumes_nc!(s)         # NC = WO2W R5TAP (Wensel-Krumland) taper + DVEW r5harv CA-hardwood D²H
    s.variant isa WestSierra && return compute_volumes_ws!(s)      # WS = 500WO2W R5TAP + 500DVEW r5harv (reuse NC kernels) — chunk 8
    s.variant isa SoutheastAlaska && return compute_volumes_ak!(s) # AK = R10 VOLEQDEF→NVEL (chunk 8, not yet ported — cuft stubbed 0)
    s.variant isa SouthCentralOregon && return compute_volumes_so!(s) # SO = R6 Behre 616BEHW + INGY FW2 (so/formcl.f) — chunk 8
    s.variant isa OregonCoast && return compute_volumes_oc!(s)     # OC = BLM Behre-taper cubic (blmvol/blmtap) — chunk C10a; board-foot C10b
    s.variant isa Olympic && return compute_volumes_op!(s)         # OP = BLM Behre-taper cubic+board (blmvol/blmtap, reuse OC) — chunk 2
    s.variant isa Ontario && return compute_volumes_on!(s)        # ON = canada ZAK/HONER (volont.f) total+merch cubic + Mowraski NMV
    s.variant isa InlandEmpire && return compute_volumes!(s, InlandEmpire())
    s.variant isa BritishColumbia && return compute_volumes!(s, BritishColumbia())   # BC Kozak taper (total cubic)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; c = s.control
    # R8 board-foot rule: the R8-CLK path reports INTERNATIONAL ¼" board feet (volinit2.f:269-272 VOL(2)=
    # VOL(10)) for GW/JF (IFORST 8), Ouachita (9), Ozark-St Francis (10), and Francis Marion & Sumter (12)
    # except the Andrew Pickens district (IDIST 2); every other R8 forest keeps Scribner. IFORST/IDIST from
    # the KODFOR location code. (D11: s07_forest_808 IFORST=8 & s22_forest_809 IFORST=9 ⇒ Intl (351); snt01/
    # all_GA IFORST=1 ⇒ Scribner (285/174), unchanged.)
    _kodfor = Int(s.plot.user_forest_code)
    _iforst = (_kodfor ÷ 100) - (_kodfor ÷ 10000) * 100
    _idist  = _kodfor - (_kodfor ÷ 100) * 100
    _r8_intl = _iforst in (8, 9, 10) || (_iforst == 12 && _idist != 2)
    # Log-graded HRVRVN (unit 4 = BF_1000_LOG): capture each tree's per-log-DIB gross Scribner BF so
    # the cut path (echarv.f) can bucket it into DIB-class records for the FVS_EconHarvestValue report.
    # Gated on an active ECON with a unit-4 revenue record; otherwise a no-op (the common path).
    log_grade = s.econ !== nothing && s.econ.active && any(r -> r.unit == 4, s.econ.hrv_rev)
    log_grade && empty!(s.econ.tree_log_bf)
    # Cubic log-graded HRVRVN (unit 5 = FT3_100_LOG): the parallel per-log gross-cuft capture (R9LGCFT),
    # bucketed by DIB into FVS_EconHarvestValue cubic columns. Independent of unit 4 — both may be active.
    log_grade_cuft = s.econ !== nothing && s.econ.active && any(r -> r.unit == 5, s.econ.hrv_rev)
    log_grade_cuft && empty!(s.econ.tree_log_ft3)
    scfmin = c.sp_scf_dbhmin; scftop = c.sp_scf_topd; topd = c.sp_top_diam
    stmp = c.sp_stump_ht; scfstmp = c.sp_scf_stump; dbhmin = c.sp_dbh_min
    merch = (stmp = stmp, topd = topd, scfstmp = scfstmp, scftop = scftop,
             bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    # Board-foot equation/standards (BFVOLUME / VOLEQNUM). BFPFLG=1 (fvsvol.f:257) ⇒ board feet rides
    # the cubic call (the default, since SN's board eq+standards equal the sawtimber ones); else a
    # separate board-foot call with the board equation + BFTOPD/BFSTMP is needed. Precompute the
    # per-species flag once so the common all-default path stays a single bool test.
    bfmin = c.sp_bf_dbhmin; bftop = c.sp_bf_topd; bfstm = c.sp_bf_stump; bfeq = c.sp_bf_vol_eq
    bfpflg0 = !isempty(bfeq) && any(k -> bfmin[k] != scfmin[k] || bfstm[k] != scfstmp[k] ||
                                          bftop[k] != scftop[k] || bfeq[k] != veq[k], 1:length(veq))
    cfdef = c.sp_cf_defect; bfdef = c.sp_bf_defect          # MCDEFECT / BFDEFECT defect curves
    cff0 = c.sp_cf_form0; cff1 = c.sp_cf_form1              # MCFDLN cubic log-linear form coefs
    bff0 = c.sp_bf_form0; bff1 = c.sp_bf_form1              # BFFDLN board log-linear form coefs
    anydef_cf = any(!iszero, cfdef); anydef_bf = any(!iszero, bfdef)
    anyform = any(!iszero, cff0) || any(!=(1f0), cff1) || any(!iszero, bff0) || any(!=(1f0), bff1)
    anydef = anydef_cf || anydef_bf || anyform || any(!iszero, t.defect) # gate the no-defect hot path
    # Two reusable 15-vectors for the per-tree R8 volume calls (pillar-2: was a fresh zeros(15) per call).
    # DISTINCT buffers because v (primary cubic) and vb (board recompute) have overlapping lifetimes
    # (v[7] is read after vb is computed). Allocated once here (per compute_volumes! call), not per tree.
    _vbuf = Vector{Float32}(undef, 15); _vbbuf = Vector{Float32}(undef, 15)
    _logbuf = Vector{Float32}(undef, 40)   # sawtimber log-length scratch (both calls run sequentially)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; h = t.height[i]; sp = t.species[i]
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
            continue
        end
        # Broken-top trees: build the volume profile from the full ("normal")
        # height, then truncate it back to the break with CFTOPK (vols.f:60-120).
        tkill = h >= 4.5f0 && t.trunc[i] > 0
        tkill && (h = Float32(t.norm_ht[i]) * 0.01f0)
        if d >= scfmin[sp]
            prod = "01"; stump = scfstmp[sp]; mtopp = scftop[sp]
        else
            prod = "02"; stump = stmp[sp]; mtopp = topd[sp]
        end
        mtops = topd[sp]
        ldref = log_grade ? Base.RefValue{Dict{Int,Float32}}(Dict{Int,Float32}()) : nothing
        lcref = log_grade_cuft ? Base.RefValue{Dict{Int,Float32}}(Dict{Int,Float32}()) : nothing
        v, ht1prd, _ = _R8CLARK_VOL(veq[sp], d, h, mtopp, mtops, stump, prod; log_dib = ldref, log_cuft = lcref, intl_bf = _r8_intl, buf = _vbuf, logbuf = _logbuf)
        tcf = v[1]
        mcf = d >= dbhmin[sp] ? v[4] + v[7] : 0f0
        scf = d >= scfmin[sp] ? v[4] : 0f0
        bf = v[10]
        # Region-8 "≥10 ft of product" rule on the PRIMARY cubic call (fvsvol.f:337-347, the CF
        # section). FVS, inside `IF(IT>0 .AND. D≥DBHMIN)`, runs:
        #   IF(IREGN==8 .AND. PROD=='01' .AND. HT1PRD<10) THEN TVOL(4)=0; TVOL(2)=0; ENDIF
        # i.e. a sawtimber (prod="01") tree whose sawtimber sawlog has < 10 ft of product yields no
        # sawtimber: zero the sawtimber cubic TVOL(4) and the board feet TVOL(2). SCF=TVOL(4)→0,
        # MCF=TVOL(4)+TVOL(7) drops to TVOL(7), BF=TVOL(2)→0. Region 8 is implicit for SN; HT1PRD is
        # this primary call's value. The PROD=='01' gate is FVS's own — default small trees take
        # prod="02" and are untouched. (The BFPFLG=0 board-section copy at fvsvol.f:499 is handled by
        # the board recompute below.)
        if d >= dbhmin[sp] && prod == "01" && ht1prd < 10f0
            scf = 0f0
            mcf = v[7]
            bf  = 0f0
        end
        # Board feet rides the sawtimber call (BFPFLG=1, fvsvol.f:257) — exact by default. When
        # BFVOLUME/VOLEQNUM make the board equation or standards differ from the sawtimber ones
        # (BFPFLG=0), recompute board feet from a separate board call (BFTOPD/BFSTMP + board eq),
        # gated by BFMIND (fvsvol.f:362). That call's Region-8 "≥10 ft of product" rule
        # (fvsvol.f:499) ALSO zeros the sawtimber cubic when the board-top sawlog is < 10 ft, which
        # drops the reported sawtimber + the sawtimber part of merch cubic. (`bf` was set above and
        # may already have been zeroed by the primary-call ≥10-ft rule.)
        # `bfmax` is the BOARD equation's total cubic (VMAX-analog), which the broken-top BOARD top-kill
        # (BFTOPK) uses to fit its Behre taper (FVS vols.f:391 passes BFMAX from BFVOL, NOT the cubic-call
        # VMAX from CFVOL). With the default board (board eq == cubic eq) they're equal, so v[1] serves; when
        # VOLEQNUM/BFVOLUME split the equations (BFPFLG=0), the board top-kill must use the BOARD call's own
        # total (vb[1]), else a broken-top tree's board is scaled by the wrong (cubic) equation's taper.
        bfmax = v[1]
        if bfpflg0 && (bfmin[sp] != scfmin[sp] || bfstm[sp] != scfstmp[sp] ||
                       bftop[sp] != scftop[sp] || bfeq[sp] != veq[sp])
            if d >= bfmin[sp]
                vb, bf_ht1prd, _ = _R8CLARK_VOL(bfeq[sp], d, h, bftop[sp], topd[sp], bfstm[sp], "01"; log_dib = ldref, intl_bf = _r8_intl, buf = _vbbuf, logbuf = _logbuf)
                bf = vb[10]; bfmax = vb[1]                # BFMAX = board-equation total (fvsvol.f BFVOL)
                if bf_ht1prd < 10f0                       # Region-8: a < 10 ft board-top sawlog has
                    bf = 0f0                              # no product — zero board feet (TVOL(2))
                    scf = 0f0                             # and the sawtimber cubic (TVOL(4)), which
                    mcf = d >= dbhmin[sp] ? v[7] : 0f0    # also drops the saw part of merch cubic
                end
            else
                bf = 0f0
            end
        end
        if tkill && tcf > 0f0
            # FVS vols.f:150 computes BARK=BRATIO(D) from the START-of-cycle DBH (before the cycle's
            # DG projection at line 151), and CFTOPK/BFTOPK use THAT bark. Use the value stashed at
            # DG time; fall back to BRATIO(current DBH) at cycle-0 LSTART (no projection yet), which is
            # exactly FVS's LSTART bark. (Was: always recompute from the grown DBH ⇒ broken-top cuft ±1.)
            bark = t.vol_bark[i] > 0f0 ? t.vol_bark[i] :
                   bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, d)
            tcf, mcf, scf = cftopk(merch, sp, d, h, tcf, mcf, scf, v[1], bark, Int(t.trunc[i]))
            bf = bftopk(merch, sp, d, h, bf, bfmax, bark, Int(t.trunc[i]))
        end
        # Volume defect (FVSsn vols.f, SN branch). Two coupled corrections, both keyed off the
        # per-species DBH defect curves (MCDEFECT→CFDEFT, BFDEFECT→BFDEFT) via ALGSLP:
        #   • CUBIC (vols.f:294-325): the pulpwood/topwood part MCFV−SCFV is cut by ICDF% (ICDF≥99
        #     ⇒ all pulpwood gone); sawtimber is left for the board step.
        #   • BOARD (vols.f:419-432): board feet AND sawtimber cubic are cut by IBDF% (≥99 ⇒ both 0),
        #     applied only where board feet exist.
        # Then MCFV = PULPV + (post-board-defect SCFV), so a BFDEFECT also lowers reported merch cubic.
        # ICDF/IBDF are the LARGEST of three sources (vols.f:298): the per-tree DEFECT input, the
        # CFDEFT/BFDEFT DBH curve, and the MCFDLN/BFFDLN log-linear form model VOLCOR=exp(B0+B1·ln(V))
        # (the implied % reduction (V−VOLCOR)/V); the form coefs default to 0/1 ⇒ no-op.
        if anydef
            mcf, scf, bf = _apply_tree_defect(mcf, scf, bf, d, sp, Int(t.defect[i]),
                                              cfdef, bfdef, cff0, cff1, bff0, bff1, anydef_cf, anydef_bf)
        end
        t.cuft_vol[i]       = tcf
        t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i]   = scf
        t.bdft_vol[i]       = bf
        # Stash this tree's per-log-DIB gross BF for the cut path's log-graded revenue accumulation.
        # Only when board feet survived (defect/Region-8 zeroing) so empties don't pollute the lookup.
        log_grade && bf > 0f0 && ldref !== nothing && !isempty(ldref[]) && (s.econ.tree_log_bf[i] = ldref[])
        # Cubic stash: gate on merch cubic surviving (mcf>0), so defect/Region-8-zeroed trees don't pollute.
        log_grade_cuft && mcf > 0f0 && lcref !== nothing && !isempty(lcref[]) && (s.econ.tree_log_ft3[i] = lcref[])
    end
    return s
end
