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

"HTDBH mode-0 predicted total height (ft) for `sp` at DBH `d` (htdbh.f)."
@inline function _htdbh_height(sd, sp::Integer, d::Float32, ifor::Integer = 0)
    p2, p3, p4, db = _htdbh_params(sd, sp, ifor)
    if d >= 3f0
        return 4.5f0 + p2 * exp(-p3 * d ^ p4)
    else
        hat3 = 4.5f0 + p2 * exp(-p3 * 3f0 ^ p4)
        return (hat3 - 4.51f0) * (d - db) / (3f0 - db) + 4.51f0
    end
end

"HTDBH mode-1 inverse: dbh (in) from total height `h` (ft) for `sp` (htdbh.f kode 1)."
@inline function _htdbh_dbh(sd, sp::Integer, h::Float32, ifor::Integer = 0)
    p2, p3, p4, db = _htdbh_params(sd, sp, ifor)
    hat3 = 4.5f0 + p2 * exp(-p3 * 3f0 ^ p4)
    if h >= hat3
        ratio = (log(min(h, 4.5f0 + p2 * 0.9999f0) - 4.5f0) - log(p2)) / (-p3)
        return ratio > 0f0 ? exp(log(ratio) * (1f0 / p4)) : 100f0
    else
        return ((h - 4.51f0) * (3f0 - db) / (hat3 - 4.51f0)) + db
    end
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
    bhat = vmax / (0.00545415f0 * d^2 * bark^2 * h)
    bhat > 0.95f0 && (bhat = 0.95f0)
    ahat = 0.44277f0 - 0.99167f0 / bhat - 1.43237f0 * log(bhat) +
           1.68581f0 * sqrt(bhat) - 0.13611f0 * bhat^2
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
    return alb2 - alb1 - 2f0 * bhat * (log(alb2) - log(alb1)) -
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
    dub_missing_heights!(state)

CRATET height resolution (cratet.f:212-265): assign heights to trees missing one
and resolve the full ("normal") height of broken-top trees. Missing-height live
trees get the HTDBH curve height. Topkill trees (norm_ht<0) keep their broken
height but get `norm_ht` = full predicted height ×100 (≥ the standing height), and
a break point `trunc` (80% of standing height when none was supplied).
"""
function dub_missing_heights!(s::StandState)
    t = s.trees; sd = s.coef.species; ifor = Int(s.plot.forest_idx)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; sp = t.species[i]
        tkill = t.norm_ht[i] < 0
        if t.height[i] > 0f0 && !tkill
            continue
        end
        h_v = d <= 0.1f0 ? 1.01f0 : _htdbh_height(sd, sp, d, ifor)
        h_v < 4.5f0 && (h_v = 4.5f0)
        if !tkill
            t.height[i] = h_v
        else
            t.norm_ht[i] = round(Int32, h_v * 100f0 + 0.5f0)
            if t.trunc[i] == 0
                if t.height[i] > 0f0
                    t.trunc[i] = round(Int32, 80f0 * t.height[i] + 0.5f0)
                else
                    t.trunc[i] = round(Int32, 80f0 * h_v + 0.5f0)
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
                (t.norm_ht[i] = round(Int32, t.height[i] * 100f0))
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
    sd = s.coef.species; c = s.control
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
    yr = Int(c.cycle_year[1]) + Int(c.cycle) * round(Int, fint)
    for ev in c.volume_events
        Int(ev.year) <= yr || continue
        isp = round(Int, ev.params[1])
        if ev.icflag == Int32(217)        # VOLUME — cubic merch standards
            _set_merch_sp!(c, c.sp_dbh_min,    isp, ev.params[2])
            _set_merch_sp!(c, c.sp_top_diam,   isp, ev.params[3])
            _set_merch_sp!(c, c.sp_stump_ht,   isp, ev.params[4])
            _set_merch_sp!(c, c.sp_scf_dbhmin, isp, ev.params[5])
            _set_merch_sp!(c, c.sp_scf_topd,   isp, ev.params[6])
            _set_merch_sp!(c, c.sp_scf_stump,  isp, ev.aux)
        else                               # BFVOLUME — board-foot merch standards
            _set_merch_sp!(c, c.sp_bf_dbhmin,  isp, ev.params[2])
            _set_merch_sp!(c, c.sp_bf_topd,    isp, ev.params[3])
            _set_merch_sp!(c, c.sp_bf_stump,   isp, ev.params[4])
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
function compute_volumes!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; c = s.control
    scfmin = c.sp_scf_dbhmin; scftop = c.sp_scf_topd; topd = c.sp_top_diam
    stmp = c.sp_stump_ht; scfstmp = c.sp_scf_stump; dbhmin = c.sp_dbh_min
    merch = (stmp = stmp, topd = topd, scfstmp = scfstmp, scftop = scftop,
             bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    cfdef = c.sp_cf_defect; anydef = any(!iszero, cfdef)   # MCDEFECT active? (gate the hot path)
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
        v, _, _ = _R8CLARK_VOL(veq[sp], d, h, mtopp, mtops, stump, prod)
        tcf = v[1]
        mcf = d >= dbhmin[sp] ? v[4] + v[7] : 0f0
        scf = d >= scfmin[sp] ? v[4] : 0f0
        # Board feet rides the sawtimber call (BFPFLG=1, fvsvol.f:257) — exact by default since
        # SN's board-foot standards equal the sawtimber ones for every species. When BFVOLUME/
        # VOLUME breaks that coincidence (BFPFLG=0), Fortran recomputes board feet from a separate
        # BF-top call AND that call's Region-8 "10 ft of product" rule (fvsvol.f:499) also zeros
        # the sawtimber cubic — a coupling not yet ported; see DIVERGENCES.md.
        bf = v[10]
        if tkill && tcf > 0f0
            bark = bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, d)  # unified per-stand bark (Fort Bragg)
            tcf, mcf, scf = cftopk(merch, sp, d, h, tcf, mcf, scf, v[1], bark, Int(t.trunc[i]))
            bf = bftopk(merch, sp, d, h, bf, v[1], bark, Int(t.trunc[i]))
        end
        # MCDEFECT cubic defect (vols.f:294-332, SN branch): the pulpwood/topwood part of merch
        # cubic (MCFV−SCFV) is reduced by ICDF%, where ICDF = NINT(ALGSLP(DBH, CFDEFT)·100) clamped
        # to [0,99]; sawtimber (SCFV) is untouched by cubic defect. (Per-tree DEFECT input + the
        # CFLA0/CFLA1 log-linear form model default to no-op and are deferred; see DIVERGENCES.md.)
        if anydef && mcf > scf
            icdf = clamp(round(Int, _algslp_col(d, _DBHCLS, cfdef, sp) * 100f0), 0, 99)
            icdf > 0 && (mcf = scf + (mcf - scf) * (1f0 - icdf * 0.01f0))
        end
        t.cuft_vol[i]       = tcf
        t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i]   = scf
        t.bdft_vol[i]       = bf
    end
    return s
end
