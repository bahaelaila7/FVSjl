# =============================================================================
# regent.jl (pacificnorthwest) — PN small-tree growth. Chunk 6.
#
# PN compiles WC's regent.f / esgent.f / vwc/smhgdg.f / dgbnd.f (FVSpn_buildDir == FVSwc_buildDir) ⇒ it runs the
# shared WC REGENT driver. TWO PN data differences: (1) `pn_smhgdg` omits the DF Curtis→King SI line
# (smhgdg.f:260, VARACD≠'WC' ⇒ raw SITEAR); (2) PN's own htdbh.f (4 forest tables on IFOR, DF 5" spline on
# IFOR 2/4/6, SS D≥100 line on IFOR 1/3) for missing-height dubbing + redwood DBH.
# =============================================================================

# pn/htdbh.f — 6-forest × 39-species Curtis-Arney P2/P3/P4 (data/pacificnorthwest/htdbh_coeffs_pn.csv).
let
    P2 = zeros(Float32, 6, 39); P3 = zeros(Float32, 6, 39); P4 = zeros(Float32, 6, 39)
    for l in readlines(joinpath(PN_DATADIR, "htdbh_coeffs_pn.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        fi = parse(Int, f[1]); sp = parse(Int, f[2])
        P2[fi, sp] = parse(Float32, f[3]); P3[fi, sp] = parse(Float32, f[4]); P4[fi, sp] = parse(Float32, f[5])
    end
    global const PN_HTDBH_P2 = P2; global const PN_HTDBH_P3 = P3; global const PN_HTDBH_P4 = P4
end

# PN forkod IFOR (1..6) is used directly as the htdbh forest-table index (pn/htdbh.f takes IFOR, not JFOR):
# the CSV's six blocks are per IFOR — 1/3 OLYMPC, 4 MTHOOD, 5 WILLAM, 2/6 SIUSLW (htdbh.f:205-220).
@inline _pn_htdbh_ifor(ifor::Int)::Int = clamp(ifor, 1, 6)

# pn/htdbh.f:228-229: Douglas-fir on IFOR 2/4/6 splines at 5.0" (label 100) instead of 3.0".
@inline _pn_df_spline5(ifor::Int, sp::Int) = sp == 16 && (ifor == 2 || ifor == 4 || ifor == 6)

"pn/htdbh.f MODE=0 (D→H)."
@inline function pn_htdbh_height(ifor::Int, sp::Int, d::Float32)::Float32
    p2 = PN_HTDBH_P2[ifor, sp]; p3 = PN_HTDBH_P3[ifor, sp]; p4 = PN_HTDBH_P4[ifor, sp]
    if _pn_df_spline5(ifor, sp)                                           # htdbh.f:253-258
        d >= 5.0f0 && return 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(d, p4))
        return ((4.5f0 + p2 * fexp(-1f0 * p3 * fpow(5.0f0, p4)) - 4.51f0) * (d - 0.3f0) / 4.7f0) + 4.51f0
    end
    if d >= 3.0f0                                                         # htdbh.f:230-236
        h = 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(d, p4))
        (d >= 100f0 && sp == 6 && (ifor == 1 || ifor == 3)) && (h = 0.25f0 * d + 248f0)   # SS on the Olympic
        return h
    else
        return ((4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3.0f0, p4)) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
    end
end
"pn/htdbh.f MODE=1 (H→D)."
@inline function pn_htdbh_dbh(ifor::Int, sp::Int, h::Float32)::Float32
    p2 = PN_HTDBH_P2[ifor, sp]; p3 = PN_HTDBH_P3[ifor, sp]; p4 = PN_HTDBH_P4[ifor, sp]
    if _pn_df_spline5(ifor, sp)                                           # htdbh.f:260-266
        hat5 = 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(5.0f0, p4))
        h >= hat5 && return fexp(flog((flog(h - 4.5f0) - flog(p2)) / (-1f0 * p3)) * 1f0 / p4)
        return (((h - 4.51f0) * 4.7f0) / (4.5f0 + p2 * fexp(-1f0 * p3 * fpow(5.0f0, p4)) - 4.51f0)) + 0.3f0
    end
    hat3 = 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3.0f0, p4))                         # htdbh.f:238-246
    if h >= hat3
        d = fexp(flog((flog(h - 4.5f0) - flog(p2)) / (-1f0 * p3)) * 1f0 / p4)
        (h >= 273f0 && sp == 6 && (ifor == 1 || ifor == 3)) && (d = (h - 248f0) / 0.25f0)
        return d
    else
        return (((h - 4.51f0) * 2.7f0) / (4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3.0f0, p4)) - 4.51f0)) + 0.3f0
    end
end

# vwc/smhgdg.f — identical to wc_smhgdg EXCEPT NO DF Curtis→King line (raw SITEAR for PN).
@inline function pn_smhgdg(sp::Int, h::Float32, d::Float32, cr::Float32, ptbal::Float32,
                           ptba::Float32, si::Float32, avht::Float32)
    if sp == 17
        htmax = 2.242202f0 * si
        (htmax - h <= 1.0f0) && return (0.0f0, 0.0f0)
        age1 = (1.0f0 / -0.010742f0) * flog(1.0f0 - fpow((h / 2.242202f0 / si), 1.0f0 / 0.919076f0))
        age2 = age1 + 5.0f0
        h1 = 2.242202f0 * si * fpow((1.0f0 - fexp(-0.010742f0 * age1)), 0.919076f0)
        h2 = 2.242202f0 * si * fpow((1.0f0 - fexp(-0.010742f0 * age2)), 0.919076f0)
        return (h2 - h1, 0.0f0)
    end
    relht = 0.0f0
    avht > 0.0f0 && (relht = h / avht)
    relht > 1.5f0 && (relht = 1.5f0)
    # (PN: NO sp==16 Curtis→King — raw SITEAR)
    ptbal2 = flog(ptbal + 2.71f0); ptba2 = flog(ptba + 2.71f0); relht2 = sqrt(relht)
    boost = ptba < 100.0f0 ? 1.0f0 / (1.0f0 + fexp(-3.1f0 + 0.18f0 * ptba)) : 0.0f0
    hbh = h >= 4.5f0 ? 4.5f0 : h
    a = WC_SMH_ALPHA[sp]; beta = WC_SMH_BETA[sp]
    dgs = WC_SMH_DMAX[sp] / (1.0f0 + fexp(a[1] + a[2]*ptba + a[3]*ptba2 + a[4]*ptbal + a[5]*ptbal2 +
              a[6]*boost + a[7]*cr + a[8]*relht + a[9]*relht2 + a[10]*si))
    hg5 = dgs / beta
    local dg5::Float32
    if h < 4.5f0
        dg5 = (h + hg5) >= 4.5f0 ? beta * (h + hg5 - 4.5f0) : 0.0f0
    else
        dg5 = dgs
    end
    return (hg5, dg5)
end

# pn/regent.f, pn/esgent.f and vwc/smhgdg.f are the WC sources (FVSpn_buildDir == FVSwc_buildDir), so PN runs
# the shared WC driver (_wcpn_small_tree_growth!, wc_esgent!) with these data hooks: smhgdg.f:260 skips the DF
# Curtis→King SI (VARACD≠'WC'), and HTDBH is pn/htdbh.f on IFOR.
@inline _rg_smhgdg(::PacificNorthwest, sp, h, d, cr, ptbal, ptba, si, avht) = pn_smhgdg(sp, h, d, cr, ptbal, ptba, si, avht)
@inline _rg_ifor(::PacificNorthwest, fidx::Int) = _pn_htdbh_ifor(fidx)
@inline _rg_htdbh_dbh(::PacificNorthwest, ifor::Int, sp::Int, h::Float32) = pn_htdbh_dbh(ifor, sp, h)

small_tree_growth!(s::StandState, stash, v::PacificNorthwest; fint::Float32 = 10.0f0) =
    _wcpn_small_tree_growth!(s, stash, v; fint = fint)

regenerate!(s::StandState, ::PacificNorthwest; kwargs...) = s
