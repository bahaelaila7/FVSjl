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

"HTDBH mode-0 predicted total height (ft) for `sp` at DBH `d` (htdbh.f)."
@inline function _htdbh_height(sd, sp::Integer, d::Float32, ifor::Integer = 0)
    if _uses_wykoff(sd, sp)
        return exp(sd[:htdbh_ht1][sp] + sd[:htdbh_ht2][sp] / (d + 1f0)) + 4.5f0
    end
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
    if _uses_wykoff(sd, sp)
        return sd[:htdbh_ht2][sp] / (log(h - 4.5f0) - sd[:htdbh_ht1][sp]) - 1f0   # htdbh.f:463
    end
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
    yr = cycle_year_at(c, Int(c.cycle))   # IY schedule (TIMEINT/CYCLEAT-aware)
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
function compute_volumes!(s::StandState)
    s.variant isa Northeast && return compute_volumes_ne!(s)   # NVEL R9 Clark cubic (WIP: no board ft)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; c = s.control
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
        v, ht1prd, _ = _R8CLARK_VOL(veq[sp], d, h, mtopp, mtops, stump, prod)
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
        if bfpflg0 && (bfmin[sp] != scfmin[sp] || bfstm[sp] != scfstmp[sp] ||
                       bftop[sp] != scftop[sp] || bfeq[sp] != veq[sp])
            if d >= bfmin[sp]
                vb, bf_ht1prd, _ = _R8CLARK_VOL(bfeq[sp], d, h, bftop[sp], topd[sp], bfstm[sp], "01")
                bf = vb[10]
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
            bark = bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, d)  # unified per-stand bark (Fort Bragg)
            tcf, mcf, scf = cftopk(merch, sp, d, h, tcf, mcf, scf, v[1], bark, Int(t.trunc[i]))
            bf = bftopk(merch, sp, d, h, bf, v[1], bark, Int(t.trunc[i]))
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
            dpack = Int(t.defect[i])
            # ICDF = max(per-tree CF defect, CFDEFT curve, cubic form model on the pulpwood MCFV−SCFV).
            icdf = dpack ÷ 1000000
            (anydef_cf && mcf > scf) &&
                (icdf = max(icdf, clamp(round(Int, _algslp_col(d, _DBHCLS, cfdef, sp) * 100f0), 0, 99)))
            temvol = mcf - scf
            if temvol > 0f0 && (cff0[sp] != 0f0 || cff1[sp] != 1f0)
                volcor = exp(cff0[sp] + cff1[sp] * log(temvol))
                icdf = max(icdf, round(Int, (temvol - volcor) / temvol * 100f0))
            end
            icdf = clamp(icdf, 0, 99)
            pulpv = icdf >= 99 ? 0f0 : (mcf - scf) * (1f0 - icdf * 0.01f0)
            if bf > 0f0
                # IBDF = max(per-tree BF defect, BFDEFT curve, board form model on BFV).
                ibdf = (dpack ÷ 10000) % 100
                anydef_bf &&
                    (ibdf = max(ibdf, clamp(round(Int, _algslp_col(d, _DBHCLS, bfdef, sp) * 100f0), 0, 99)))
                if bff0[sp] != 0f0 || bff1[sp] != 1f0
                    volcorb = exp(bff0[sp] + bff1[sp] * log(bf))
                    ibdf = max(ibdf, round(Int, (bf - volcorb) / bf * 100f0))
                end
                ibdf = clamp(ibdf, 0, 99)
                if ibdf >= 99
                    bf = 0f0; scf = 0f0
                elseif ibdf > 0
                    f = 1f0 - ibdf * 0.01f0; bf *= f; scf *= f
                end
            end
            mcf = pulpv + scf
        end
        t.cuft_vol[i]       = tcf
        t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i]   = scf
        t.bdft_vol[i]       = bf
    end
    return s
end
