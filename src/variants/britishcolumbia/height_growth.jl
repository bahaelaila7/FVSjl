# =============================================================================
# height_growth.jl (britishcolumbia) — BC large-tree height-growth engine hook (canada/bc/htgf.f).
#
# CHUNK 4. `height_growth!(s, ::BritishColumbia; scale)` fills trees.ht_growth per live tree,
# mirroring dgf!/IE height_growth!. V3 regime (LV2ATV=.FALSE., zones ICH/IDF/SBS/SBPS):
#   HTG(i) = V3HTG(sp, HT, DBH, DG) · Y,  Y = SCALE · XHMULT(sp) · MISHGF   (htgf.f:1612-1635)
# V3HTG (bc_v3_htg) already folds MtoFT·MLT (MLT=HCOR2). MISHGF (mistletoe) = 1 for BC.
# HTG uses chunk-3 DG output (t.diam_growth). ⚠ V2 regime (other zones, htgf.f:1618-1622 exp-form) TODO.
# =============================================================================

"""height_growth!(s, ::BritishColumbia; scale) — per-tree periodic height increment → trees.ht_growth (V3)."""
function height_growth!(s::StandState, ::BritishColumbia; scale::Float32 = 1.0f0)
    p, t, ctl = s.plot, s.trees, s.control
    zone, series = bc_stand_zone(s)
    cur_year = current_cycle_year(s)
    nsp = nspecies(BritishColumbia())
    ip = zeros(Int, nsp)
    @inbounds for sp in 1:nsp
        ip[sp] = bc_resolve_lthg(sp, series, zone)
    end
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i])
        (sp < 1 || sp > nsp || ip[sp] < 1) && continue          # LTH%FIT false ⇒ HTG 0
        d = t.dbh[i]; hti = t.height[i]
        (d <= 0f0 || hti <= 0f0) && continue
        dg = t.diam_growth[i]
        mlt = ctl.htg_cor2[sp] > 0f0 ? ctl.htg_cor2[sp] : 1f0    # LTH%MLT = HCOR2 (htgf.f:1810)
        v3  = bc_v3_htg(sp, ip[sp], hti, d, dg; mlt = mlt)
        xht = active_multiplier(ctl, :htg, sp, cur_year)         # XHMULT (HTGMULT keyword; MISHGF=1 for BC)
        t.ht_growth[i] = v3 * scale * xht
    end
    return s
end
