# =============================================================================
# fire/snag.jl — snag falldown + decay dynamics (FFE chunk F7-state core)
#
# Ported from: bin/FVSsn_buildDir/fmsfall.f (FMSFALL snag falldown) + the DECAYX
# hard→soft transition (fmvinit.f / fmsnag.f).
#
# When the fire (or ordinary mortality) kills a tree it becomes a standing snag. Each
# year a fraction of the snags fall (transferring to coarse woody debris) and the
# standing-hard snags decay toward soft. These are the per-snag-record rates that drive
# the stateful snag list (the records + per-cycle loop are F7-state's container, which
# builds on these functions). Rates come from the species' snag class (`snag_fallx`,
# `snag_alldwn`, `snag_decayx` in `fire_species_props.csv`).
# =============================================================================

"""
    snag_fall_density(coef, ksp, d, origden, denttl) -> Float32

Density of snags (stems/acre) that fall in one year for a snag record of species `ksp`
and DBH `d` (FMSFALL, fmsfall.f). `origden` is the record's original density, `denttl`
the density still standing. Small snags (< 12" and not redcedar) fall at a linear
`MODRATE·origden`; large snags use the last-5% logic that ramps the final stems down to
zero by the species' `ALLDWN` year.
"""
function snag_fall_density(coef::SpeciesCoefficients, ksp::Integer, d::Float32,
                           origden::Float32, denttl::Float32)::Float32
    base = max(0.01f0, -0.001679f0 * d + 0.064311f0)
    modrate = min(1f0, base * coef_col(coef, :snag_fallx)[ksp])
    if d < 12f0 && ksp != 2                            # small snag (redcedar=2 keeps last-5% logic)
        return modrate * origden
    end
    alldwn = coef_col(coef, :snag_alldwn)[ksp]
    x = (0.05f0 - 1f0) / (-modrate)                    # year at which 5% remain
    fallm2 = alldwn <= x ? 2f0 : 0.05f0 / (alldwn - x) # final fall rate (last 5%)
    if denttl <= 0.05f0 * origden
        return fallm2 * origden
    end
    dfalln = modrate * origden
    if denttl < dfalln + 0.05f0 * origden              # don't overshoot below 5% in one step
        dfalln = denttl - origden * (0.05f0 - fallm2)
    end
    return dfalln
end

"""
    snag_decay_fraction(coef, ksp) -> Float32

Annual fraction of standing-hard snags of species `ksp` that transition to soft decay
(DECAYX, fmvinit.f — e.g. a 12" tree goes soft in 2/6/10 years for snag class 1/2/3).
"""
@inline snag_decay_fraction(coef::SpeciesCoefficients, ksp::Integer) =
    coef_col(coef, :snag_decayx)[ksp]

"""
    add_snag!(fs, sp, dbh, density, year)

Create a standing-dead snag cohort (FMSADD) for `density` stems/acre of species `sp`,
DBH `dbh`, that died in `year`. New snags start fully hard. No-op for non-positive
density.
"""
function add_snag!(fs::FireState, sp::Integer, dbh::Float32, density::Float32, year::Integer;
                   bolevol::Float32 = 0f0, height::Float32 = 0f0)
    density > 0f0 || return
    sn = fs.snags
    push!(sn.sp, Int32(sp));   push!(sn.dbh, dbh)
    push!(sn.den_hard, density); push!(sn.den_soft, 0f0)
    push!(sn.origden, density);  push!(sn.year, Int32(year)); push!(sn.bolevol, bolevol)
    push!(sn.height, height)
    return
end

"""
    snag_bole_carbon(s) -> Float32

Snag STEM-VOLUME bole carbon in tons C/acre — the faithful FFE Stand-Dead snag basis
(`TOTSNG = (SNVIS+SNVIH)·V2T`, fmdout.f:153): each cohort's death-time stem-volume biomass
(`bolevol`, cuft·V2T) × its still-standing density, × 0.5. This is the bole half of Stand-Dead;
the crown half is CWD2B. (Static here — the snag height-loss that shrinks the bole over time is the
next refinement.) Falls back to Jenkins aboveground for cohorts with `bolevol` unset (e.g. fire snags).
"""
function snag_bole_carbon(s::StandState)::Float32
    fs = s.fire; fs === nothing && return 0f0
    sn = fs.snags; coef = s.coef; c = 0f0
    @inbounds for i in eachindex(sn.sp)
        den = sn.den_hard[i] + sn.den_soft[i]
        den > 0f0 || continue
        b = sn.bolevol[i]
        b <= 0f0 && (b = let (a, _, _) = jenkins_biomass(coef, sn.sp[i], sn.dbh[i]); a end)
        c += b * den
    end
    return c * 0.5f0
end

# CWD down-wood size class (1–9) from a stem diameter, matching the FUINI breakpoints
# (<0.25, .25–1, 1–3, 3–6, 6–12, 12–20, 20–35, 35–50, >50 inches).
@inline _cwd_size_class(d::Float32) =
    d < 0.25f0 ? 1 : d < 1f0 ? 2 : d < 3f0 ? 3 : d < 6f0 ? 4 :
    d < 12f0 ? 5 : d < 20f0 ? 6 : d < 35f0 ? 7 : d < 50f0 ? 8 : 9

# Diameter-class breakpoints BP(0:9), inches (fmcwd.f:56). Index j+1 ↔ BP(j).
const _CWD_BP = (0f0, 0.25f0, 1f0, 3f0, 6f0, 12f0, 20f0, 35f0, 50f0, 9999f0)

"""
    _cwd_cone_fractions(d, ht) -> NTuple{9,Float32}

Fraction of a fallen bole's volume in each CWD size class 1..9, summing to 1 — the cone-taper
distribution of FVS's FMCWD/CWD1 (fmcwd.f label 1000). The standing stem is modeled as a cone
(DBH `d` in at 4.5 ft, total height `ht` ft); each diameter-class breakpoint BP falls at a height
BPH, and the normalized conic volume `P(h)=r(h)²·(ht−h)/(R1²·ht)` between successive breakpoints
(clipped to the [0.1 ft, ht] bole) is that class's share. Replaces the prior single-class dump
(whole bole into the DBH class), which overloaded the 6–12" class. Normalized to sum 1 so the bole
TOTAL is unchanged (carbon DDW preserved); only the per-class split — FMCFMD's (SMALL,LARGE) input —
is corrected. A bole too short to taper (`ht ≤ 4.6`) or with no taperable volume falls back to its
DBH class.
"""
function _cwd_cone_fractions(d::Float32, ht::Float32)
    f = zeros(Float32, 9)
    d <= 0.1f0 && (d = 0.1f0)
    if ht <= 4.6f0
        f[_cwd_size_class(d)] = 1f0
        return f
    end
    htd = ht
    rhrat = ((htd * 12f0) - 54f0) / (0.5f0 * d)
    # BPH(j) = height (ft) where stem diameter = BP(j); index j+1
    bph = ntuple(j -> max(0.10f0, htd - (0.5f0 * _CWD_BP[j] * rhrat) / 12f0), 10)
    loht = 0.10f0; hiht = htd
    r1 = d * 0.0416666667f0                       # radius (ft) at DBH (= d/12 * 0.5)
    r1 = r1 + loht * ((r1 * htd) / (htd - 4.5f0)) # extend cone to the stem base
    r1sq = r1 * r1
    total = 0f0
    @inbounds for j in 1:9
        bphj = bph[j + 1]; bphjm1 = bph[j]        # BPH(j), BPH(j-1)
        (hiht <= bphj || loht > bphjm1) && continue
        hicut = min(hiht, bphjm1); locut = max(loht, bphj)
        locut == hicut && continue
        r2 = r1 * (1f0 - hicut / htd); p1 = (r2 * r2 * (htd - hicut)) / (r1sq * htd)
        r2 = r1 * (1f0 - locut / htd); p2 = (r2 * r2 * (htd - locut)) / (r1sq * htd)
        dif = max(0f0, p2 - p1); f[j] = dif; total += dif
    end
    total <= 0f0 && (f[_cwd_size_class(d)] = 1f0; total = 1f0)
    @inbounds for j in 1:9; f[j] /= total; end
    return f
end

"""
    update_snags!(s, nyears) -> Float32

Advance every snag cohort `nyears` years (FMSNAG): each year the hard snags decay toward
soft (`snag_decay_fraction`) and a `snag_fall_density` share falls — split proportionally
between the hard and soft pools (fmsnag.f:197-221). The fallen snags transfer into the
coarse-woody-debris pools (`fire.cwd`, CWD1): the fallen aboveground biomass (Jenkins ×
fallen density) is added to the down-wood class for the stem DBH and the species' decay
class. Returns the total density (stems/ac) that fell.
"""
function update_snags!(s::StandState, nyears::Integer)::Float32
    fs = s.fire; (fs === nothing) && return 0f0
    sn = fs.snags; coef = s.coef
    cur = Int(current_cycle_year(s))
    fallen = 0f0
    @inbounds for i in eachindex(sn.sp)
        sp = sn.sp[i]
        # Advance the snag only by the years it has actually STOOD since death, capped at the cycle
        # length: a snag that died this cycle (deaths are dated at/near the cycle boundary, not the
        # cycle start) has stood only ~0-1 years, so it must not fall a full cycle's worth — otherwise
        # the cycle's fresh mortality over-falls and Stand-Dead collapses (FMSNAG ages each snag by its
        # own (year − deathyr), not a blanket nyears). Older cohorts fall the full nyears each cycle.
        yrs = clamp(cur - Int(sn.year[i]), 0, Int(nyears))
        yrs > 0 || continue
        # a falling snag transfers its BOLE biomass to down wood; the crown is the separate CWD2B
        # path (so don't double-count it). Fall back to Jenkins for cohorts with bolevol unset.
        a = sn.bolevol[i]
        a <= 0f0 && (a = let (j, _, _) = jenkins_biomass(coef, sp, sn.dbh[i]); j end)
        idc = Int(coef_col(coef, :dkr_cls)[sp])             # decay-rate class
        # Distribute the fallen bole down the cone taper across size classes (FMCWD/CWD1) instead of
        # dumping the whole bole into the DBH class. Fractions depend only on (dbh, height) → compute
        # once per cohort. Height unset (0) ⇒ single-class fallback (no behavior change).
        frac = _cwd_cone_fractions(sn.dbh[i], sn.height[i])
        for _ in 1:yrs
            denttl = sn.den_hard[i] + sn.den_soft[i]
            denttl > 0f0 || break
            shift = min(sn.den_hard[i], sn.den_hard[i] * snag_decay_fraction(coef, sp))
            sn.den_hard[i] -= shift; sn.den_soft[i] += shift
            denttl = sn.den_hard[i] + sn.den_soft[i]
            dfall = min(denttl, snag_fall_density(coef, sp, sn.dbh[i], sn.origden[i], denttl))
            dfis = denttl > 0f0 ? sn.den_soft[i] * dfall / denttl : 0f0
            dfih = denttl > 0f0 ? sn.den_hard[i] * dfall / denttl : 0f0
            sn.den_soft[i] -= dfis; sn.den_hard[i] -= dfih
            add = a * dfall                                 # fallen biomass this step (tons/ac)
            for j in 1:9
                frac[j] > 0f0 && (fs.cwd[j, 2, idc] += add * frac[j])  # spread across CWD size classes
            end
            fallen += dfall
        end
    end
    return fallen
end

"Total standing snag density (stems/ac) currently in the snag list."
snag_standing_density(fs::FireState) = sum(fs.snags.den_hard) + sum(fs.snags.den_soft)

# Snag-summary DBH class lower bounds (SNPRCL, fmvinit.f:46-51): cumulative thresholds — class i counts
# every snag with DBH ≥ SNPRCL[i], so class 1 (≥0) equals the total.
const _FM_SNPRCL = (0f0, 12f0, 18f0, 24f0, 30f0, 36f0)

"""
    snag_summary(s) -> (; hard, soft)

The FFE snag-summary densities (stems/ac) for the FVS_SnagSum table (FMSSUM, fmssum.f:28-53): `hard`
and `soft` are each a 7-tuple — the standing HARD (`den_hard`) and SOFT (`den_soft`) snag density in the
six cumulative DBH classes (≥0/12/18/24/30/36 in) plus the total in slot 7. FVSjl's per-record
hard/soft split is the current-state equivalent of the Fortran's DENIH/DENIS + HARD flag.
"""
function snag_summary(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return (hard = ntuple(_ -> 0f0, 7), soft = ntuple(_ -> 0f0, 7))
    sn = fs.snags; thd = zeros(Float32, 7); tsf = zeros(Float32, 7)
    @inbounds for i in eachindex(sn.sp)
        dh = sn.den_hard[i]; ds = sn.den_soft[i]; d = sn.dbh[i]
        thd[7] += dh; tsf[7] += ds                          # slot 7 = total (all snags)
        for c in 1:6
            d >= _FM_SNPRCL[c] || continue
            thd[c] += dh; tsf[c] += ds
        end
    end
    return (hard = Tuple(thd), soft = Tuple(tsf))
end


"""
    ffe_seed_input_snags!(s) -> StandState

Seed the FFE snag list from the INPUT dead-tree records (the tree-list dead partition `n+1 : n+ndead`,
history codes 6-9) at stand initialization (FMSDIT→FMSADD ITYP=3, fmsdit.f:135). Each becomes a
standing-dead cohort: its STEM-volume bole (`cuft·V2T`) → the Stand-Dead bole, its coarse roots → the
Below-Dead BIOROOT pool. Input snags carry no crown (history ≥7 ⇒ crown already fallen; the records have
`crown_pct = 0`), so there is no CWD2B contribution. Heights/volumes are computed locally here because
`compute_volumes!` covers only the live partition (the dead records arrive with height 0). The snags are
young enough (TSOFT > the inventory age) to be hard, so the static stem bole is exact (no FMSVOL
height-loss yet). Populates the inventory-cycle Stand-Dead / Below-Dead carbon. No-op without FFE / dead records.
"""
function ffe_seed_input_snags!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active || s.trees.ndead <= 0) && return s
    t = s.trees; coef = s.coef; c = s.control; sd = coef.species
    v2t = coef_col(coef, :v2t); ifor = Int(s.plot.forest_idx)
    # Input snags PRE-EXIST the inventory: FVS books the input dead trees at IY(1)−FINTM (a measurement
    # period before the start), so the FFE snag falldown (update_snags!, which ages by current−deathyr)
    # ages them from the inventory rather than holding them frozen-full the first cycle. Use the cycle
    # period as FINTM (the common no-GROWTH-keyword case).
    per = round(Int, c.year); per < 1 && (per = 5)
    yr = Int(current_cycle_year(s)) - per
    s.control.merch_init || init_merch_standards!(s)
    @inbounds for i in (t.n + 1):(t.n + t.ndead)
        den = t.tpa[i]; d = t.dbh[i]
        (den > 0f0 && d >= 1f0) || continue
        sp = Int(t.species[i])
        h = t.height[i] > 0f0 ? t.height[i] : max(4.5f0, _htdbh_height(sd, sp, d, ifor))
        prod, stump, mtopp = d >= c.sp_scf_dbhmin[sp] ?
            ("01", c.sp_scf_stump[sp], c.sp_scf_topd[sp]) : ("02", c.sp_stump_ht[sp], c.sp_top_diam[sp])
        v, _, _ = _R8CLARK_VOL(s.species.vol_eq[sp], d, h, mtopp, c.sp_top_diam[sp], stump, prod)
        # FVS's snag bole (FMDOUT→FMSVOL→CFVOL, fmdout.f:146) is the MERCHANTABLE cubic to the top
        # diameter (v[4]), NOT the gross total-stem cubic (v[1]). For small snags the <top-dia tip is a
        # large fraction (sp27 d7.2: v[1]=5.2 vs v[4]=4.8 = FVS); for large stems v[1]≈v[4]. Using v[4]
        # makes snag_bole_carbon match FVS (3.92→3.77 vs 3.8) and lowers the small-snag falldown into DDW.
        bolevol = v[4] * v2t[sp] / 2000f0
        add_snag!(fs, sp, d, den, yr; bolevol = bolevol, height = h)
        _, _, rbio = jenkins_biomass(coef, sp, d)
        # FVS assumes input snags have been dead 10 years for dead-root decay (fmsadd.f:313-320):
        # XDCAY = (1−CRDCAY)^10. FVSjl was booking the full root biomass (over-counting Below-Dead).
        fs.bioroot += rbio * den * (1f0 - _FM_CRDCAY)^10
    end
    return s
end
