# =============================================================================
# establishment.jl — regeneration / establishment (ESNUTR → ESTAB)
#
# Ported from: base/esnutr.f (cycle hook) + base/estab.f (tree creation) +
# base/estab_helpers.f (ESSUBH/ESTIME) + base/esinit.f.
#
# SN's PARTIAL (keyword-driven) establishment model: no auto-ingrowth. When an
# ESTAB packet scheduled PLANT(430)/NATURAL(431) activities are due, ESTAB creates
# the regen trees. A single bare plot (NPTIDS=1) is replicated MINREP=50 times: each
# replicate independently draws an established HEIGHT per species (ESSUBH height-at-
# age + a BACHLO draw on the establishment RNG), and contributes one record per
# species carrying plantedTPA·survival/100 / dupnpt TPA (400/50 = 8) ⇒ 50×2 = 100
# records, 800 TPA. The trees enter AFTER growth+mortality (GRADD order) so they're
# fresh (full TPA) this period; their DBH is derived from the established height.
# =============================================================================

const _ES_MINREP = 50          # MINREP: target plot replication (esinit.f)

# XMIN: per-species establishment min height (blkdat.f; sp 73-90 default 0).
const _ES_XMIN = Float32[
    0.50,2.08,0.50,1.00,1.32,2.51,0.50,2.53,2.75,0.50,
    5.05,0.50,4.70,0.50,1.33,1.33,0.66,2.40,1.35,1.35,
    2.03,0.50,0.50,0.50,0.50,2.08,0.51,0.63,2.08,2.08,
    2.08,2.08,0.50,0.50,0.50,0.92,0.50,5.98,0.94,2.08,
    0.50,3.28,3.28,1.33,0.89,1.53,1.38,3.59,3.59,3.59,
    2.08,2.08,4.15,3.59,3.59,2.08,2.08,2.08,0.89,0.50,
    0.50,0.50,1.38,1.38,1.38,0.50,2.75,2.75,0.50,2.75,
    0.50,1.38, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0]
# HHTMAX: per-species max establishment height (blkdat.f).
const _ES_HHTMAX = Float32[23.0,27.0,21.0,21.0,22.0,20.0,24.0,18.0,18.0,17.0,22.0,
    (20.0 for _ in 12:90)...]

"Establishment-tree disturbance date (TALLY=427 trigger), else inventory − 20 (esnutr.f:63)."
function _es_idsdat(s::StandState)
    for a in s.control.schedule
        a.icflag == Int32(427) && return Int(a.params[1])
    end
    return Int(s.control.cycle_year[1]) - 20
end

"""
    establish!(state; fint=5f0) -> Bool

Create scheduled PLANT/NATURAL regen for the current cycle (ESNUTR/ESTAB). Runs at
the end of `grow_cycle!` (GRADD order). Idempotent per year. Returns whether any
tree was created. No-op unless an ESTAB packet is active.
"""
function establish!(s::StandState; fint::Float32 = 5f0)::Bool
    s.estab.active || return false
    t = s.trees; sd = s.coef.species
    per = round(Int, fint)
    yr = Int32(Int(s.control.cycle_year[1]) + Int(s.control.cycle) * per)
    yr in s.estab.years_done && return false
    due = [a for a in s.control.schedule
           if (a.icflag == Int32(430) || a.icflag == Int32(431)) && yr <= a.year < yr + per]
    isempty(due) && return false

    idsdat = _es_idsdat(s)
    nptids = max(1, Int(s.plot.points_inv))
    idup   = max(1, fld(_ES_MINREP, nptids))
    dupnpt = Float32(nptids * idup)
    bc = (sd[:ht_curve_b1], sd[:ht_curve_b2], sd[:ht_curve_b3], sd[:ht_curve_b4], sd[:ht_curve_b5])
    montane = !isempty(s.plot.eco_unit) && s.plot.eco_unit[1] == 'M'
    ifor = Int(s.plot.forest_idx)
    # gentim/delay/trage timing (esnutr/estab/essubh): age = FINT − delay − gentim + trage.
    gentim = (Int(yr) + per - idsdat) - per; gentim < 0 && (gentim = 0)
    created = false
    @inbounds for rep in 1:idup
        # per-replicate establishment RNG draws (estab.f:216-221): two for emsqr
        # (unused on the no-treeht path), one for esdraw (the re-seed value).
        esrann!(s.rng); esrann!(s.rng)
        esdraw = floor(esrann!(s.rng) * 100000f0 + 0.5f0)
        for a in due
            sp = round(Int, a.params[1]); (1 <= sp <= MAXSP) || continue
            ptree = a.params[2] * (a.params[3] / 100f0) / dupnpt
            ptree <= 0f0 && continue
            delay  = Int(a.year) - Int(yr)
            trage  = a.params[4] < 0.5f0 ? 2f0 : a.params[4]; trage > 10f0 && (trage = 10f0)
            age = Float32(per) - Float32(delay) - Float32(gentim) + trage; age < 1f0 && (age = 1f0)
            si  = s.plot.sp_site_index[sp]
            hht = htcalc_height(bc, sp, si, age, montane)          # ESSUBH base height
            treeht = a.params[5]
            if treeht >= 0.1f0                                      # PLANT specified a height
                hht = treeht; xh = log(hht)
                while true
                    xxh = exp(bachlo(s.rng, xh, 0.5f0; stream = :estab))
                    (0.5f0 * hht <= xxh <= 2f0 * hht) && (hht = xxh; break)
                end
            else                                                   # default ±N(0.5,0.25)
                while true
                    ran = bachlo(s.rng, 0.5f0, 0.25f0; stream = :estab)
                    (0f0 <= ran <= 1.5f0) && (hht += ran; break)
                end
            end
            hht < _ES_XMIN[sp]   && (hht = _ES_XMIN[sp])           # HTADJ=0, floor, cap
            hht > _ES_HHTMAX[sp] && (hht = _ES_HHTMAX[sp])
            ibrkup = floor(Int, ptree / 10f0 + 1f0); brk = Float32(ibrkup)
            dbh = _htdbh_dbh(sd, sp, hht, ifor); dbh < 0.1f0 && (dbh = 0.1f0)
            for _ in 1:ibrkup
                n = t.n + 1; n > length(t.dbh) && break
                t.n = n
                t.species[n]     = Int32(sp)
                t.dbh[n]         = dbh
                t.height[n]      = hht
                t.tpa[n]         = ptree / brk
                t.plot_id[n]     = Int32(rep)
                t.crown_pct[n]   = Int32(40)
                t.crown_ratio[n] = 40f0
                t.norm_ht[n]     = Int32(0)
                t.sort_key[n]    = Float64(n)
                created = true
            end
        end
        es = esdraw; (es % 2f0 == 0f0) && (es += 1f0); s.rng.es0 = Float64(es)  # ESRNSD(true,esdraw)
    end
    push!(s.estab.years_done, yr)
    created && compute_density!(s)
    return created
end
