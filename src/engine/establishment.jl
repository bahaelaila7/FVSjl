# =============================================================================
# establishment.jl — regeneration / establishment (ESNUTR → ESTAB)
#
# Ported from: base/esnutr.f (cycle hook) + base/estab.f (tree creation).
#
# SN's PARTIAL (keyword-driven) establishment model: no auto-ingrowth. When an
# ESTAB packet scheduled PLANT(430)/NATURAL(431) activities that are due, ESTAB
# creates the regen trees. Per the oracle (bare_plant): a single bare plot
# (NPTIDS=1) is replicated MINREP=50 times, so each PLANT species becomes `idup`
# records each carrying `plantedTPA·survival/100 / dupnpt` TPA (400/50 = 8) — i.e.
# 50 records × 2 species = 100 records, 800 TPA.
#
# STAGE 4a (here): record count + TPA + seedling DBH, so the establishment FIRES and
# the stand regenerates to the right TPA. STAGE 4b (TODO): bit-exact established
# heights/age/CR via ESGENT/ESSUBH driven by the establishment RNG (ESRANN, seed
# 55329), which makes BA/QMD bit-exact.
# =============================================================================

const _ES_MINREP = 50          # MINREP: target plot replication (esinit.f)

"""
    establish!(state; fint=5f0)

Create scheduled PLANT/NATURAL regen for the current cycle (ESNUTR/ESTAB). Runs at
the top of `grow_cycle!` (after CUTS) so the new seedlings grow during the cycle.
Idempotent per year. No-op unless an ESTAB packet is active.
"""
function establish!(s::StandState; fint::Float32 = 5f0)::Bool
    s.estab.active || return false
    t = s.trees
    per = round(Int, fint)
    yr = Int32(Int(s.control.cycle_year[1]) + Int(s.control.cycle) * per)
    yr in s.estab.years_done && return false
    # PLANT(430)/NATURAL(431) activities scheduled within THIS cycle's window [yr, yr+per)
    due = [a for a in s.control.schedule
           if (a.icflag == Int32(430) || a.icflag == Int32(431)) && yr <= a.year < yr + per]
    isempty(due) && return false
    # plot replication (esinit.f MINREP / estab.f dup loop). Single bare plot ⇒
    # idup = MINREP, dupnpt = NPTIDS·idup.
    nptids = max(1, Int(s.plot.points_inv))
    idup   = max(1, fld(_ES_MINREP, nptids))
    dupnpt = Float32(nptids * idup)
    created = false
    @inbounds for a in due
        sp   = round(Int, a.params[1])
        tpa  = a.params[2]; surv = a.params[3]
        (1 <= sp <= MAXSP) || continue
        ptree = tpa * (surv / 100f0) / dupnpt
        ptree <= 0f0 && continue
        for rep in 1:idup
            n = t.n + 1
            n > length(t.dbh) && break          # MAXTRE capacity
            t.n = n
            t.species[n]     = Int32(sp)
            t.dbh[n]         = 0.1f0             # seedling (stage 4b: ESGENT)
            t.height[n]      = 1.0f0
            t.tpa[n]         = ptree
            t.plot_id[n]     = Int32(rep)
            t.crown_pct[n]   = Int32(40)
            t.crown_ratio[n] = 40f0
            t.norm_ht[n]     = Int32(0)
            t.sort_key[n]    = Float64(n)
            created = true
        end
    end
    push!(s.estab.years_done, yr)
    created && compute_density!(s)               # new trees change stand density
    return created
end
