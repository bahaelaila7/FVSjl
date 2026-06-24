# =============================================================================
# small_tree_growth.jl — Southern small-tree (regeneration) growth (REGENT)
#
# Ported from: base/regent.f (growth section, LSTART=false, lestb=false).
#
# Trees below XMAX = 3" DBH grow by a HEIGHT-driven model, not the large-tree
# DGF: a Chapman-Richards height increment (HTCALC) that, over the blend band
# [XMIN=1, XMAX=3], is mixed with the large-tree HTGF increment (weight xwt),
# perturbed by a ±10% random effect (BACHLO, active when DGSD≥1), then converted
# to a diameter increment through the HTDBH height→dbh inverse. This OVERRIDES
# whatever DGF/HTGF produced for those records. Runs after DGDRIV+HTGF, before
# mortality, matching grincr.f's call order (so MORTS sees the corrected small-
# tree DG). Each live small tree spends 3 BACHLO draws (central+upper+lower) when
# tripling is active, 1 otherwise — the upper/lower records are stashed for
# `triple_records!`.
#
# snt01 scope: con = RHCON·exp(HCOR) = 1 (no small-tree height calibration),
# HGADJ = 1, MULTS multipliers = 1, LHTDRG=false ⇒ the HTDBH-inverse branch.
# =============================================================================

const REGENT_XMIN = 1f0
const REGENT_XMAX = 3f0
const REGENT_REGYR = 5f0
const REGENT_DGMAX = 5f0

# Budwidth floor DIAM (in) by species (regent.f DATA / htdbh SNDBAL).
const REGENT_DIAM = Float32[
    0.1,0.3,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
    0.5,0.4,0.5,0.5,0.2,0.2,0.1,0.2,0.2,0.2,
    0.2,0.2,0.3,0.1,0.1,0.2,0.3,0.3,0.1,0.2,
    0.1,0.2,0.1,0.2,0.2,0.2,0.2,0.1,0.2,0.2,
    0.1,0.3,0.4,0.2,0.2,0.2,0.2,0.2,0.2,0.2,
    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.1,0.1,
    0.2,0.1,0.2,0.2,0.1,0.1,0.2,0.1,0.2,0.2,
    0.2,0.1,0.1,0.2,0.2,0.1,0.1,0.2,0.2,0.1,
    0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0.2,0.2]

# Diameter increment for one small-tree record from its height increment `htg`
# (regent.f:195-275, HTDBH-inverse branch). `d`,`h` are the record's pre-growth
# dbh/height; returns the outside-bark DG (before DGBND).
@inline function _regent_dg(sd, bark_a, bark_b, sp::Integer, d::Float32, h::Float32, htg::Float32,
                            scale2::Float32, dgmx::Float32, ifor::Integer = 0, xrdgro::Float32 = 1f0)
    hk = h + htg
    hk <= 4.5f0 && return 0f0                         # (DBH bump path; not hit for snt01)
    dkk = _htdbh_dbh(sd, sp, hk, ifor)
    dk  = h <= 4.5f0 ? d : _htdbh_dbh(sd, sp, h, ifor)
    bark = bark_ratio(bark_a, bark_b, sp, d)          # per-stand bark (Fort Bragg override)
    # REGDMULT: regent.f:347/350 multiply DG by XRDGRO before the DDS-scale conversion.
    dg = ((dk < 0f0 || dkk < 0f0) ? htg * 0.2f0 * bark : (dkk - dk) * bark) * xrdgro
    dg < 0f0 && (dg = 0.1f0)
    dg > dgmx && (dg = dgmx)
    dds = dg * (2f0 * bark * d + dg) * scale2
    dg = sqrt((d * bark)^2 + dds) - bark * d
    return dg
end

"""
    small_tree_growth!(state, stash; fint=5f0)

Apply the REGENT small-tree height/diameter growth, overriding DGF/HTGF for
records with dbh < `REGENT_XMAX`. Writes the central record into
`trees.diam_growth`/`trees.ht_growth`; when `stash !== nothing` it also fills the
tripled upper/lower DG+HTG (`stash.dgU/dgL/htgU/htgL`, `stash.is_small`). Must run
after `height_growth!` (it blends with and reads `trees.ht_growth`).
"""
function small_tree_growth!(s::StandState, stash; fint::Float32 = 5f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    bc = (sd[:ht_curve_b1], sd[:ht_curve_b2], sd[:ht_curve_b3], sd[:ht_curve_b4], sd[:ht_curve_b5])
    montane = !isempty(p.eco_unit) && p.eco_unit[1] == 'M'
    sizcap = s.control.sp_size_cap
    dlo_v = sd[:dg_bound_dbh_lo]; dhi_v = sd[:dg_bound_dbh_hi]
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    species_sort!(s)
    trip = stash !== nothing
    nrec = trip ? 3 : 1
    scale = fint / REGENT_REGYR                       # fnt/REGYR
    scale2 = 1f0                                       # YR/fnt (normal cycle)
    random_on = s.control.dg_stddev_bound >= 1f0   # DGSD (DGSTDEV; default 2, grinit.f) ⇒ random effect on
    # REGHMULT/REGDMULT (MULTS kinds 3/6): per-species regen height/diameter multipliers
    # (regent.f:233 HTGR·XRHGRO, :347/:350 DG·XRDGRO). cur_year = inventory + cycle·period.
    cur_year = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        si = p.sp_site_index[sp]
        htmax = htcalc_htmax(bc, sp, si, montane)
        dgmx = REGENT_DGMAX * scale
        con = exp(c.htg_cor_small[sp])    # RHCON·exp(HCOR); RHCON = 1 (no HCOR2 keyword)
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        for k3 in i1:i2
            i = ind1[k3]
            d = t.dbh[i]
            (d >= REGENT_XMAX || t.tpa[i] <= 0f0) && continue
            h = t.height[i]
            # base small-tree height increment (HGADJ=1; XRHGRO=REGHMULT)
            if htmax - h <= 1f0
                htgr_s = 0.1f0
            else
                aget = htcalc_age(bc, sp, si, h, montane)
                htg1 = htcalc_incr(bc, sp, si, aget, montane)
                htgr_s = max(htg1 * con * scale * xrhgro, 0.1f0)
            end
            xwt = d <= REGENT_XMIN ? 0f0 : (d - REGENT_XMIN) / (REGENT_XMAX - REGENT_XMIN)
            htgr = max(htgr_s * (1f0 - xwt) + xwt * t.ht_growth[i], 0.1f0)
            for l in 0:(nrec - 1)
                ran = 0f0
                if random_on
                    while true
                        ran = bachlo(s.rng, 0f0, 1f0)
                        (-1f0 <= ran <= 1f0) && break
                    end
                end
                htg = max(htgr + ran * 0.1f0 * htgr, 0.1f0)
                (h + htg) > sizcap[sp, 4] && (htg = max(sizcap[sp, 4] - h, 0.1f0))
                dg = _regent_dg(sd, c.bark_a, c.bark_b, sp, d, h, htg, scale2, dgmx, Int(p.forest_idx), xrdgro)
                (d + dg) < REGENT_DIAM[sp] && (dg = REGENT_DIAM[sp] - d)
                dg = dg_bound(dlo_v, dhi_v, sp, d, dg, sizcap)
                if l == 0
                    t.diam_growth[i] = dg; t.ht_growth[i] = htg
                elseif l == 1
                    stash.dgU[i] = dg; stash.htgU[i] = htg; stash.is_small[i] = true
                else
                    stash.dgL[i] = dg; stash.htgL[i] = htg
                end
            end
        end
    end
    return s
end
