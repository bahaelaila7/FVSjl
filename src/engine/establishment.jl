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

# XMIN: per-species establishment min height (blkdat.f) lives in
# data/southern/species_coefficients.csv as the `estab_min_ht` column.
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
    es_xmin = sd[:estab_min_ht]   # per-species establishment min height (CSV)
    per = round(Int, fint)
    yr = Int32(current_cycle_year(s))   # IY schedule; yr+per below = next boundary (fint is per-cycle)
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
    pccf = 0f0          # point crown competition factor (≈0 for the sparse established plots)
    created = false
    nstart = t.n        # tree count before establishment (phase-2 crown pass starts here)
    # estab.f outer loop: `for nn in 1:NPTIDS` (each inventory point) × `idup` replicates
    # → NPTIDS·idup records total. For a BARE stand every point is identical (BAAA=0,
    # uniform slope/aspect/habitat), so the per-point variables don't vary; only the
    # record count and the ESRANN draw count scale with NPTIDS. (ptree already divides by
    # dupnpt = NPTIDS·idup, so the planted TPA is conserved across all the records.)
    @inbounds for nn in 1:nptids, rep in 1:idup
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
            hht < es_xmin[sp]   && (hht = es_xmin[sp])           # HTADJ=0, floor, cap
            hht > _ES_HHTMAX[sp] && (hht = _ES_HHTMAX[sp])
            ibrkup = floor(Int, ptree / 10f0 + 1f0); brk = Float32(ibrkup)
            # REGENT establishment dbh (regent.f:331-334, LESTB branch): DBH = HTDBH⁻¹(HK),
            # floored to the species min DIAM, then a small height-proportional add 0.001·HK.
            dbh = _htdbh_dbh(sd, sp, hht, ifor); dbh < 0.1f0 && (dbh = 0.1f0)
            dbh += 0.001f0 * hht
            for _ in 1:ibrkup
                n = t.n + 1; n > length(t.dbh) && break
                t.n = n
                t.species[n]     = Int32(sp)
                t.dbh[n]         = dbh
                t.height[n]      = hht
                t.tpa[n]         = ptree / brk
                # Records go on inventory point `nn` (estab.f:313 ITRE=IPTIDS[nn]).
                # point_ba scales each point's raw BA by PI/GROSPC with PI=NPTIDS, so with
                # the planted TPA spread evenly over NPTIDS points each point_ba comes back
                # to the full stand BA — matching the oracle's pba=ba_v fallback (PTBAA≤0)
                # for fresh establishment, for any NPTIDS (NPTIDS=1 ⇒ this is point 1).
                t.plot_id[n]     = Int32(nn)
                t.crown_pct[n]   = Int32(0)            # crown set in phase 2 (REGENT lestb)
                t.crown_ratio[n] = 0f0
                t.norm_ht[n]     = Int32(0)
                t.sort_key[n]    = Float64(n)
                created = true
            end
        end
        es = esdraw; (es % 2f0 == 0f0) && (es += 1f0); s.rng.es0 = Float64(es)  # ESRNSD(true,esdraw)
    end
    # PHASE 2 — ESGENT → REGENT(lestb): assign each new tree its open-grown crown in
    # SPESRT (species-then-record) order (regent.f:107-116). cr = 0.89722 −
    # 0.0000461·PCCF + 0.07985·N(0,1)[±1], clamp [0.20,0.90]; the crown draw uses the
    # MAIN RANN stream (separate from the ESRANN heights). The per-cycle CROWN
    # (crown_ratio_update!, run after) then applies its ±1%/yr change limit (~85).
    if created
        newidx = sort(collect((nstart + 1):t.n); by = i -> (Int(t.species[i]), i))
        @inbounds for i in newidx
            ran_cr = 0f0
            while true
                ran_cr = bachlo(s.rng, 0f0, 1f0)
                -1f0 <= ran_cr <= 1f0 && break
            end
            cr = clamp(0.89722f0 - 0.0000461f0 * pccf + 0.07985f0 * ran_cr, 0.20f0, 0.90f0)
            icr0 = floor(Int32, cr * 100f0 + 0.5f0)
            t.crown_pct[i]   = icr0
            t.crown_ratio[i] = Float32(icr0)
        end
        compute_density!(s)
    end
    push!(s.estab.years_done, yr)
    return created
end
