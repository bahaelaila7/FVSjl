# =============================================================================
# keyword_dispatch.jl — process a keyword stream into stand state
#
# Ported from: base/initre.f  (the giant computed-GOTO keyword processor).
#
# INITRE reads keyword records and dispatches each to a handler that mutates the
# stand. We replace the 147-way computed GOTO with named handler functions keyed
# on the keyword. A few keywords pull additional raw lines from the stream
# (TREEFMT reads 2, STDIDENT reads 1); those handlers take the `KeywordReader`.
#
# Scope note: this implements the keywords exercised by the SN test set; others
# are recognized as no-ops (KNOWN_NOOP) so processing doesn't choke. Handlers are
# added here as later chunks need them (thinning, fire, econ, ...).
# =============================================================================

# Keywords recognized but not yet acted on (flags/reports that don't change the
# cycle-0 stand state). Acting on these is added in their respective chunks.
const KNOWN_NOOP = Set([
    "SCREEN", "NOSCREEN", "STATS", "NOAUTOES", "TREELIST", "ECHOSUM", "ECHO",
    "NOECHO", "NOSUM", "NODEBUG", "DEBUG", "CALBSTAT", "COMPRESS", "REWIND",
    "ATRTLIST", "CUTLIST", "MANAGED", "ENDFILE", "FVSSTAND",
    # bare-stand / establishment-adjacent flags (no cycle-0 stand effect yet)
    "NOTREES", "NOHTDREG", "AUTOES",
])

"Read one raw (un-lexed) line from the keyword stream, advancing the record count."
function read_raw_line!(kr::KeywordReader)
    kr.record_count += 1
    return readline(kr.io)
end

# --- individual keyword handlers (each ported from its initre.f label) --------

# OPTION 10 — DESIGN (initre.f:743): plot design.
function kw_design!(s::StandState, rec::KeywordRecord)
    p, v = s.plot, rec.values
    rec.present[1] && (p.baf = v[1])
    rec.present[2] && (p.fixed_plot_inv = v[2])
    rec.present[3] && (p.min_dbh_var_plot = v[3])
    rec.present[4] && (p.points_inv = nint(v[4]))
    rec.present[5] && (p.nonstockable = nint(v[5]))
    rec.present[6] && (p.sample_weight = v[6])
    (p.sample_weight <= 0f0 && p.points_inv > 0) && (p.sample_weight = Float32(p.points_inv))
    g = v[7]
    (g > 1f0 && g <= 100f0) && (g *= 0.01f0)
    (g > 0f0 && g <= 1f0) && (p.gross_space = g)
    return
end

# OPTION 11 — NUMCYCLE (initre.f:771).
function kw_numcycle!(s::StandState, rec::KeywordRecord)
    n = rec.values[1]
    (n >= 1f0 && nint(n) <= MAXCYC) && (s.control.ncycle = nint(n))
    return
end

# OPTION 2 — TIMEINT (vbase/initre.f:533): cycle length (IY). Field 1 = cycle index (0/absent =
# all cycles), field 2 = period length in years (default 10). Sets s.control.year (YR/IFINT), the
# cycle length threaded through the growth models — DDS/HTG scale by FINT/5, autcor/year-age by it.
# UNIFORM path only (all cycles same length); per-cycle lengths (field 1 > 0) deferred.
function kw_timeint!(s::StandState, rec::KeywordRecord)
    v = rec.values; pr = rec.present
    len = pr[2] ? Float32(v[2]) : 10f0
    cyc = pr[1] ? nint(v[1]) : 0
    cyc <= 0 && (s.control.year = len)
    return
end

# OPTION 16 — INVYEAR (initre.f:895).
function kw_invyear!(s::StandState, rec::KeywordRecord)
    rec.present[1] && (s.control.cycle_year[1] = nint(rec.values[1]))
    return
end

# OPTION 94 — SITECODE (initre.f:2546): site index for the site species / all.
# Fields: array[1]=species (or 0=all), array[2]=site index.
function kw_sitecode!(s::StandState, rec::KeywordRecord)
    si = rec.values[2]
    isp = nint(rec.values[1])
    if isp <= 0
        fill!(s.plot.sp_site_index, si)
    elseif 1 <= isp <= MAXSP
        s.plot.sp_site_index[isp] = si
        s.plot.site_species = Int32(isp)        # ISISP — the site species
    end
    rec.present[2] && (s.plot.site_index = si)
    return
end

# OPTION 14 — STDINFO (initre.f:808): stand info. Forest/habitat decoding
# (FORKOD/HABTYP) is deferred; the geometry/site fields are set here.
function kw_stdinfo!(s::StandState, rec::KeywordRecord)
    p, v = s.plot, rec.values
    p.user_forest_code = nint(v[1])
    # SN: STDINFO field 2 is the habitat/ecological-unit field, decoded by HABTYP
    # (numeric → index into SNECU; alpha → matched, uppercased) into the PCOM code.
    rec.present[2] && (p.eco_unit = rpad(resolve_eco_unit(rec.fields[2], rec.values[2]), 10))
    rec.present[3] && (p.stand_age = nint(v[3]))
    rec.present[4] && (p.aspect = v[4] * 0.0174533f0)   # degrees → radians (utils.f)
    rec.present[5] && (p.slope  = v[5] / 100f0)         # percent → fraction (utils.f)
    (rec.present[6] && v[6] > 0f0) && (p.elevation = v[6])
    if rec.present[9]
        org = nint(v[9])
        p.stand_origin = (org < 0 || org > 1) ? Int32(0) : org
    end
    # Fort Bragg (forkod.f:137-140): forest 701 (location 701xx) ⇒ IFOR=20 (special
    # longleaf/loblolly DG + bark equations) AND KODFOR is remapped to NC Uwharrie
    # district 81110 (region 8) for downstream FORTYP + VOLEQDEF. Without the KODFOR
    # remap, VOLEQDEF sees region 7 (70106÷10000) and assigns NO R8 Clark equation ⇒
    # every tree gets zero volume.
    if div(p.user_forest_code, 100) == 701
        p.forest_idx = Int32(20)
        p.user_forest_code = Int32(81110)
    end
    # FORKOD phase 3: default lat/long/elev from the forest code (forkod.f:193).
    lat0, long0, elev0 = forest_location(s.coef, div(p.user_forest_code, 100))
    p.latitude  == 0f0 && (p.latitude  = lat0)
    p.longitude == 0f0 && (p.longitude = long0)
    p.elevation == 0f0 && (p.elevation = elev0)
    return
end

# Thinning/harvest keyword → CUTS method code (icflag). Extended per method as the
# cuts! port lands; THINDBH is the first (milestone 1). (cuts.f label dispatch.)
const _THIN_ICFLAG = Dict("THINBTA" => Int32(3), "THINATA" => Int32(4),
                          "THINBBA" => Int32(5), "THINABA" => Int32(6),
                          "THINPRSC" => Int32(7), "THINDBH" => Int32(8),
                          "THINSDI" => Int32(10), "THINHT" => Int32(12),
                          "THINRDEN" => Int32(14), "THINAUTO" => Int32(1),
                          "THINCC" => Int32(11))

# Parse a THIN* activity: field 1 = calendar year, fields 2-7 = the 6 method params.
# Stores a ScheduledActivity for `cuts!` to apply on the matching cycle.
function kw_thin!(s::StandState, rec::KeywordRecord, icflag::Int32)
    v = rec.values
    yr = nint(v[1])
    params = ntuple(i -> Float32(v[i + 1]), 6)
    push!(s.control.schedule, ScheduledActivity(Int32(yr), icflag, params))
    return
end

# Growth/mortality multipliers (BAIMULT/HTGMULT/MORTMULT/REGHMULT/REGDMULT — MULTS,
# base/mults.f). Field 1 = date, field 2 = species (0 = all), field 3 = multiplier.
function kw_mult!(s::StandState, rec::KeywordRecord, kind::Symbol)
    v = rec.values; pr = rec.present
    # MORTMULT (kind :mort) additionally carries a DBH window (PRM(3)/PRM(4) → XMDIA1/XMDIA2,
    # morts.f:170-171); the others ignore fields 4-5. Defaults: D1=0, D2=99999 (all trees).
    windowed = kind === :mort || kind === :fixdg || kind === :fixhtg || kind === :crn
    d1 = (windowed && pr[4]) ? Float32(v[4]) : 0f0
    d2 = (windowed && pr[5] && Float32(v[5]) > 0f0) ? Float32(v[5]) : 99999f0
    push!(s.control.multipliers,
          GrowthMultiplier(kind, Int32(nint(v[1])), Int32(nint(v[2])), Float32(v[3]), d1, d2))
    return
end

"""
    active_mort_mult(control, sp, year, dbh) -> Float32

The MORTMULT multiplier in effect for species `sp` at cycle `year`, applied only when
the tree's `dbh` falls in the keyword's DBH window [d1, d2) (morts.f:518). Outside the
window — or when no MORTMULT applies — returns 1.0. Same most-recent / species-specific
precedence as `active_multiplier`.
"""
function active_mort_mult(c::Control, sp::Integer, year::Integer, dbh::Real)
    isempty(c.multipliers) && return 1f0
    val = 1f0; d1 = 0f0; d2 = 99999f0; bestyr = typemin(Int32); bestspec = false
    @inbounds for m in c.multipliers
        (m.kind === :mort && sp_field_matches(c, m.species, sp) && m.year <= year) || continue
        spec = m.species != 0
        if m.year > bestyr || (m.year == bestyr && spec && !bestspec)
            val = m.value; d1 = m.d1; d2 = m.d2; bestyr = m.year; bestspec = spec
        end
    end
    return (d1 <= dbh < d2) ? val : 1f0
end

"""
    active_crn_mult(control, sp, year, dbh) -> Float32

The CRNMULT crown-ratio multiplier in effect for species `sp` at cycle `year`, when the
tree's `dbh` is in the keyword's CLOSED DBH window [d1, d2] (sn/crown.f:318). Persists from
the keyword date onward (most recent / species-specific wins). 1.0 when none applies.
"""
function active_crn_mult(c::Control, sp::Integer, year::Integer, dbh::Real)
    isempty(c.multipliers) && return 1f0
    val = 1f0; d1 = 0f0; d2 = 99999f0; bestyr = typemin(Int32); bestspec = false
    @inbounds for m in c.multipliers
        (m.kind === :crn && sp_field_matches(c, m.species, sp) && m.year <= year) || continue
        spec = m.species != 0
        if m.year > bestyr || (m.year == bestyr && spec && !bestspec)
            val = m.value; d1 = m.d1; d2 = m.d2; bestyr = m.year; bestspec = spec
        end
    end
    return (d1 <= dbh <= d2) ? val : 1f0
end

"""
    apply_fix_scalers!(s, stash, kind, fint)

FIXDG/FIXHTG (grincr.f:451-525): a ONE-SHOT per-cycle scaler. In the cycle whose year
range [cyc_start, cyc_start+period) contains the keyword date, multiply DG (kind=:fixdg)
or HTG (:fixhtg) by `value` for every tree of the matching species (0 = all) whose DBH is
in [d1, d2). The tripled upper/lower records (stash dgU/dgL, htgU/htgL) get the same factor
(FVS scales DG(ITFN)/DG(ITFN+1)). Runs after all growth, before mortality (MORTS reads the
scaled DG). Multiple scalers firing the same cycle compound, in keyword order.
"""
function apply_fix_scalers!(s::StandState, stash, kind::Symbol, fint::Float32)
    isempty(s.control.multipliers) && return s
    period = round(Int, s.control.year)
    cyc_start = Int(s.control.cycle_year[1]) + Int(s.control.cycle) * period
    cyc_end = cyc_start + period
    t = s.trees
    nlive = stash === nothing ? t.n : stash.nlive
    isdg = kind === :fixdg
    @inbounds for m in s.control.multipliers
        m.kind === kind || continue
        # one-shot: a fixed-year scaler matches exactly one cycle's [start,end) range. A
        # date before the first cycle fires in cycle 0 (Fortran OPFIND past-date behaviour).
        (cyc_start <= m.year < cyc_end || (m.year < cyc_start && s.control.cycle == 0)) || continue
        for i in 1:nlive
            sp_field_matches(s.control, m.species, t.species[i]) || continue
            d = t.dbh[i]
            (m.d1 <= d < m.d2) || continue
            if isdg
                t.diam_growth[i] *= m.value
                stash !== nothing && (stash.dgU[i] *= m.value; stash.dgL[i] *= m.value)
            else
                t.ht_growth[i] *= m.value
                stash !== nothing && (stash.htgU[i] *= m.value; stash.htgL[i] *= m.value)
            end
        end
    end
    return s
end

# SPGROUP (vbase/initre.f:4726): define a species group. The keyword-line field is the group
# name (optional → auto "GROUPnn"); the NEXT record lists the member species (alpha or numeric
# codes, space-separated, ALL/0/dups skipped, max 90, `&` continues — rare, not handled). Groups
# are numbered in definition order and referenced from a species field by the negative index −N.
# Max 30 groups.
function kw_spgroup!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    line = read_raw_line!(kr)
    length(s.control.sp_groups) >= 30 && return
    members = Int32[]
    for tok in split(line)
        u = uppercase(tok)
        (u == "ALL" || u == "0") && continue
        n = tryparse(Int, tok)
        idx = n !== nothing ? (n > 0 ? Int32(n) : Int32(0)) :
              first(resolve_species(tok, s.variant, s.species, s.coef))
        (idx <= 0 || idx in members) && continue          # skip invalid + duplicates
        push!(members, idx)
        length(members) >= 90 && break
    end
    push!(s.control.sp_groups, members)
    nm = isempty(rec.fields) ? "" : strip(rec.fields[1])
    push!(s.control.sp_group_names,
          isempty(nm) ? "GROUP" * lpad(length(s.control.sp_groups), 2, '0') : uppercase(String(nm)))
    s.control.n_spgroups = Int32(length(s.control.sp_groups))
    return
end

"""
    sp_field_matches(control, isp, sp) -> Bool

Does species `sp` match a keyword species field `isp`? 0 = all species, >0 = that single
species, <0 = every member of SPGROUP group −isp. The one place the ISPCC<0 group branch lives.
"""
@inline function sp_field_matches(c::Control, isp::Integer, sp::Integer)
    isp == 0 && return true
    isp > 0 && return Int(isp) == Int(sp)
    g = -Int(isp)
    return g <= length(c.sp_groups) && (Int32(sp) in c.sp_groups[g])
end

# TREESZCP (base/keywds.f:51 / SIZCAP): a per-species maximum tree size. Field 1 =
# species (0 = all), 2 = cap DBH (SIZCAP[1]), 3 = annual mortality rate applied to
# trees that reach the cap (SIZCAP[2]), 4 = no-mortality flag IDMFLG (SIZCAP[3]),
# 5 = height cap (SIZCAP[4]). Applied immediately to every matching species; the cap
# then governs diameter growth (dgbnd), height growth (htgf) and a size-cap mortality
# floor (morts). No date field — the cap holds for the whole projection.
function kw_treeszcp!(s::StandState, rec::KeywordRecord)
    v = rec.values; pr = rec.present
    isp   = nint(v[1])
    cap   = pr[2] ? Float32(v[2]) : 999f0
    mrate = pr[3] ? Float32(v[3]) : 1f0
    flag  = pr[4] ? Float32(v[4]) : 0f0
    htcap = (pr[5] && Float32(v[5]) > 0f0) ? Float32(v[5]) : 999f0
    sc = s.control.sp_size_cap
    # species selection: 0 = all, >0 = single, <0 = SPGROUP group −isp.
    rng = isp < 0 ? (isp <= -1 && -isp <= length(s.control.sp_groups) ? s.control.sp_groups[-isp] : Int32[]) :
          isp == 0 ? (1:size(sc, 1)) : (isp:isp)
    @inbounds for sp in rng
        sc[sp, 1] = cap; sc[sp, 2] = mrate; sc[sp, 3] = flag; sc[sp, 4] = htcap
    end
    return
end

# HTGSTOP (act 110) / TOPKILL (act 111) top-damage events (htgstp.f). Keyword fields:
# date, species (0/neg = all/group), HT1, HT2, PRB (damage prob), AVEPRB, STDPBR (BACHLO
# mean/sd of the kill proportion). Stored as a ScheduledActivity (icflag = act code).
function kw_htgstp!(s::StandState, rec::KeywordRecord, act::Integer)
    v = rec.values; pr = rec.present
    f(i, dflt) = (pr[i] ? Float32(v[i]) : dflt)
    p = (f(2, 0f0), f(3, 0f0), f(4, 9999f0), f(5, 1f0), f(6, 0f0), f(7, 0f0))
    push!(s.control.htgstp_events, ScheduledActivity(Int32(nint(v[1])), Int32(act), p))
    return
end

# FIXMORT (activity 97, morts.f:781): a one-shot forced-mortality override. Fields: date,
# species (0 = all), rate, DBH window d1/d2, option PRM(5) (0=replace,1=add,2=max,3=mult →
# IP 1/2/3/4), point/size flag PRM(6). Normal (non-concentration) path only — PRM(6)≥10 is
# deferred. The rate is clamped here exactly as Fortran (≤1 when PRM(5) absent; [0,1] when
# PRM(5)<3) so the apply step is a clean per-tree formula.
function kw_fixmort!(s::StandState, rec::KeywordRecord)
    v = rec.values; pr = rec.present
    isp = Float32(nint(v[2]))
    rate = Float32(v[3])
    d1 = (pr[4] ? Float32(v[4]) : 0f0); d1 < 0f0 && (d1 = 0f0)
    d2 = (pr[5] && Float32(v[5]) > 0f0) ? Float32(v[5]) : 999f0
    ip = 1
    if pr[6]                      # PRM(5) given (NP>4)
        p5 = Float32(v[6])
        p5 < 3f0 && (rate = clamp(rate, 0f0, 1f0))
        ip = p5 == 1f0 ? 2 : p5 == 2f0 ? 3 : p5 == 3f0 ? 4 : 1
    else
        rate > 1f0 && (rate = 1f0)
    end
    pflag = pr[7] ? Float32(v[7]) : 0f0
    push!(s.control.fixmort_events,
          ScheduledActivity(Int32(nint(v[1])), Int32(97), (isp, rate, d1, d2, Float32(ip), pflag)))
    return
end

# VOLUME (volkey.f:9915) — override the CUBIC merch standards for a species (0=all, <0=group).
# Keyword fields after the date: 2=species, 3=DBHMIN, 4=TOPD, 5=STMP, 6=FRMCLS, 7=METHC,
# 8=SCFMIND, 9=SCFTOPD, 10=SCFSTMP. FRMCLS/METHC are not used by the SN R8 Clark taper path
# (it has no form-class or method selector), so we carry only the 7 standards the model reads:
# species + DBHMIN/TOPD/STMP/SCFMIND/SCFTOPD in params, SCFSTMP in aux.
function kw_volume!(s::StandState, rec::KeywordRecord)
    v = rec.values
    isp = Float32(nint(v[2]))
    push!(s.control.volume_events,
          ScheduledActivity(Int32(nint(v[1])), Int32(217),
              (isp, Float32(v[3]), Float32(v[4]), Float32(v[5]), Float32(v[8]), Float32(v[9])),
              Float32(v[10])))
    return
end

# BFVOLUME (volkey.f:9905) — override the BOARD-FOOT merch standards for a species.
# Fields after the date: 2=species, 3=BFMIND, 4=BFTOPD, 5=BFSTMP, 6=FRMCLS, 7=METHB
# (FRMCLS/METHB unused by the R8 Clark path).
function kw_bfvolume!(s::StandState, rec::KeywordRecord)
    v = rec.values
    isp = Float32(nint(v[2]))
    push!(s.control.volume_events,
          ScheduledActivity(Int32(nint(v[1])), Int32(218),
              (isp, Float32(v[3]), Float32(v[4]), Float32(v[5]), 0f0, 0f0), 0f0))
    return
end

# Load a per-species defect curve (CFDEFT or BFDEFT rows, sdefet.f:151): the five values at DBH
# 5/10/15/20/25" go to rows 2..6, the 25" value extends flat to rows 7,8,9 (DBH 30/35/40"); row 1
# (DBH 0") stays 0. `isp` selects 0=all / +species / −SPGROUP.
function _set_defect!(c, dmat::Matrix{Float32}, isp::Integer, vals::NTuple{5,Float32})
    setone(sp) = @inbounds begin
        dmat[2, sp] = vals[1]; dmat[3, sp] = vals[2]; dmat[4, sp] = vals[3]
        dmat[5, sp] = vals[4]; dmat[6, sp] = vals[5]
        dmat[7, sp] = vals[5]; dmat[8, sp] = vals[5]; dmat[9, sp] = vals[5]
    end
    if isp == 0
        for sp in 1:size(dmat, 2); setone(sp); end
    elseif isp > 0
        isp <= size(dmat, 2) && setone(isp)
    else
        g = -isp
        (1 <= g <= length(c.sp_groups)) || return
        for sp in c.sp_groups[g]; setone(sp); end
    end
    return
end

_defect_vals(v) = (Float32(v[3]), Float32(v[4]), Float32(v[5]), Float32(v[6]), Float32(v[7]))

# MCDEFECT (sdefet.f, IACTK 215) — per-species CUBIC defect curve (CFDEFT). Fields: 1=date,
# 2=species, 3..7 = defect fractions at DBH 5/10/15/20/25". Applied immediately (undated path,
# the common case — affects cycle 0 like sdefet.f's "date not defined" branch; dated deferred).
kw_mcdefect!(s::StandState, rec::KeywordRecord) =
    _set_defect!(s.control, s.control.sp_cf_defect, nint(rec.values[2]), _defect_vals(rec.values))

# BFDEFECT (sdefet.f, IACTK 216) — per-species BOARD-FOOT defect curve (BFDEFT). Same fields;
# reduces board feet AND sawtimber cubic by IBDF% in the volume path.
kw_bfdefect!(s::StandState, rec::KeywordRecord) =
    _set_defect!(s.control, s.control.sp_bf_defect, nint(rec.values[2]), _defect_vals(rec.values))

"""
    apply_fixmort!(s, killed, n, fint)

FIXMORT (morts.f:781-1042, normal path): a one-shot forced-mortality override applied AFTER
the BA-check, the last word on `killed[]`. In the cycle whose range holds the event date, for
each of the `n` original records of the matching species (0 = all) with d1≤DBH<d2, override the
kill by the option: 1 replace (P·rate), 2 add ((P−killed)·rate), 3 max, 4 multiply (killed·rate,
≤P). Point/size concentration (PRM(6)≥10) is deferred — those events are skipped.
"""
function apply_fixmort!(s::StandState, killed::Vector{Float32}, n::Int, fint::Float32)
    isempty(s.control.fixmort_events) && return
    period = round(Int, s.control.year)
    cyc_start = Int(s.control.cycle_year[1]) + Int(s.control.cycle) * period
    cyc_end = cyc_start + period
    t = s.trees
    for ev in s.control.fixmort_events
        (cyc_start <= ev.year < cyc_end || (ev.year < cyc_start && s.control.cycle == 0)) || continue
        isp = round(Int, ev.params[1]); rate = ev.params[2]
        d1 = ev.params[3]; d2 = ev.params[4]; ip = round(Int, ev.params[5])
        pflag = ev.params[6]
        kbig = (pflag == 10f0 || pflag == 11f0) ? 1 : (pflag == 20f0 || pflag == 21f0) ? 2 : 0
        if kbig == 0
            # NORMAL per-tree override (morts.f:1017): set the kill in place per the IP mode.
            @inbounds for i in 1:n
                sp_field_matches(s.control, isp, t.species[i]) || continue
                d = t.dbh[i]; (d1 <= d < d2) || continue
                p = t.tpa[i]
                if ip == 1
                    killed[i] = p * rate
                elseif ip == 2
                    killed[i] += max(0f0, p - killed[i]) * rate
                elseif ip == 3
                    killed[i] = max(killed[i], p * rate)
                else
                    killed[i] = min(p, killed[i] * rate)
                end
            end
        else
            # SIZE concentration (morts.f:838): pool the window's mortality into XMORE per the IP
            # (zeroing the in-window kills where the IP replaces), then re-impose it WHOLE-RECORD
            # at a time on the trees ranked by grown DBH — KBIG=1 smallest-first (BOTTOM UP),
            # =2 largest-first (TOP DOWN) — until XMORE is spent. KPOINT multi-plot concentration
            # is deferred.
            ba = s.calib.bark_a; bb = s.calib.bark_b
            xmore = 0f0
            @inbounds for i in 1:n
                sp_field_matches(s.control, isp, t.species[i]) || continue
                d = t.dbh[i]; (d1 <= d < d2) || continue
                p = t.tpa[i]
                if ip == 1
                    xmore += p * rate; killed[i] = 0f0
                elseif ip == 2
                    xmore += max(0f0, p - killed[i]) * rate
                elseif ip == 3
                    tmp = max(killed[i], p * rate); tmp > killed[i] && (xmore += tmp - killed[i])
                else
                    xmore += killed[i] * rate; killed[i] = 0f0
                end
            end
            # rank by WORK3 = ∓(DBH+DG/bark), RDPSRT descending (morts.f:889): KBIG=1 negates so the
            # smallest grown trees sort first (bottom up); KBIG=2 keeps the largest first (top down).
            # Reuse the faithful _rdpsrt! so the tie-break among equal grown DBHs matches FVS.
            gdbh(i) = t.dbh[i] + t.diam_growth[i] / bark_ratio(ba, bb, t.species[i], t.dbh[i])
            sgn = kbig == 1 ? -1f0 : 1f0
            key = Float32[sgn * gdbh(i) for i in 1:n]
            ord = Vector{Int32}(undef, n); _rdpsrt!(key, ord)
            credit = 0f0
            @inbounds for ix in ord
                sp_field_matches(s.control, isp, t.species[ix]) || continue
                d = t.dbh[ix]; (d1 <= d < d2) || continue
                avail = t.tpa[ix] - killed[ix]
                if credit + avail <= xmore + 0.0001f0
                    credit += avail; killed[ix] = t.tpa[ix]
                else
                    killed[ix] += xmore - credit; break
                end
            end
        end
    end
    return
end

"""
    htgstp!(s; fint)

HTGSTOP / TOPKILL top-damage events (htgstp.f). Runs once per cycle, after TRIPLE/MORTS and
before UPDATE applies the increments (gradd.f:158). For each event firing this cycle (one-shot,
matched by date), damage the trees of its species whose height is in (HT1, HT2]:
  * act 110 (HTGSTOP): scale the height increment HTG by PKIL∈[0,1].
  * act 111 (TOPKILL): reduce height to H·(1−PKIL) (PKIL≤0.8); for a tall (H≥25, D≥6) tree
    whose Behre top diameter ≥4 set NORMHT/ITRUNC (permanent broken top); cut the crown ratio.
PKIL = BACHLO(AVEPRB, STDPBR) — deterministic (= AVEPRB, no RNG) when STDPBR≤0; a RANN escape
draw skips a tree when PRB<1. Records are visited in species-sorted (IND1) order so the RNG
stream matches FVS when the event is stochastic.
"""
function htgstp!(s::StandState; fint::Float32 = 5f0)
    isempty(s.control.htgstp_events) && return s
    period = round(Int, s.control.year)
    cyc_start = Int(s.control.cycle_year[1]) + Int(s.control.cycle) * period
    cyc_end = cyc_start + period
    t = s.trees; bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    species_sort!(s)
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    for ev in s.control.htgstp_events
        (cyc_start <= ev.year < cyc_end || (ev.year < cyc_start && s.control.cycle == 0)) || continue
        isp = round(Int, ev.params[1]); ht1 = ev.params[2]; ht2 = ev.params[3]
        prb = ev.params[4]; aveprb = ev.params[5]; stdpbr = ev.params[6]
        act = ev.icflag
        sp1 = isp <= 0 ? 1 : isp; sp2 = isp <= 0 ? MAXSP : isp
        @inbounds for sp in sp1:sp2
            # isp<0 (SPGROUP group) loops all species but keeps only members.
            (isp < 0 && !sp_field_matches(s.control, isp, sp)) && continue
            i1 = isct[sp, 1]; i1 == 0 && continue
            i2 = isct[sp, 2]
            for k in i1:i2
                i = ind1[k]
                t.tpa[i] <= 0f0 && continue
                h = t.height[i]
                (h <= ht1 || h > ht2) && continue
                if prb <= 0.99999f0
                    rann!(s.rng) > prb && continue
                end
                pkil = bachlo(s.rng, aveprb, stdpbr)
                pkil <= 0f0 && continue
                if act == 110
                    pkil > 1f0 && (pkil = 1f0)
                    t.ht_growth[i] *= pkil
                else
                    pkil > 0.8f0 && (pkil = 0.8f0)
                    topk = h * pkil; toph = h - topk
                    itrc2 = floor(Int32, toph * 100f0 + 0.5f0)
                    if t.trunc[i] > 0                       # already topkilled: lower the truncation
                        t.trunc[i] > itrc2 && (t.trunc[i] = itrc2)
                        t.height[i] = toph
                        continue
                    end
                    brk = bark_ratio(bark_a, bark_b, sp, t.dbh[i])
                    d = t.dbh[i] * brk
                    if !(h < 25f0 || d < 6f0)              # tall enough to maybe break permanently
                        af = t.cuft_vol[i] / (0.00545415f0 * d * d * h)
                        af = 0.44244f0 - (0.99167f0 / af) - 1.43237f0 * log(af) +
                             1.68581f0 * sqrt(af) - 0.13611f0 * af * af
                        dtk = topk / h
                        dtk = (dtk / (af * dtk + (1f0 - af))) * d
                        if dtk >= 4f0                       # permanent broken top
                            t.trunc[i] = itrc2
                            t.norm_ht[i] = floor(Int32, h * 100f0 + 0.5f0)
                        end
                    end
                    t.height[i] = toph
                    iod = t.crown_pct[i]                    # crown reduction (skip if bug-adjusted <0)
                    if iod >= 0
                        cn = (Float32(iod) / 100f0) * h - h + toph
                        new = floor(Int32, cn / toph * 100f0 + 0.5f0)
                        new < 5 && (new = Int32(5))
                        t.crown_pct[i] = -new
                    end
                end
            end
        end
    end
    return s
end

"""
    active_multiplier(control, kind, sp, year) -> Float32

The growth/mortality multiplier in effect for species `sp` at the cycle `year`
(MULTS): the most recent matching keyword wins, a species-specific one beating an
all-species (`species==0`) one of the same date. Returns 1.0 when none apply (the
common case — short-circuited).
"""
function active_multiplier(c::Control, kind::Symbol, sp::Integer, year::Integer)
    isempty(c.multipliers) && return 1f0
    val = 1f0; bestyr = typemin(Int32); bestspec = false
    @inbounds for m in c.multipliers
        (m.kind === kind && (m.species == 0 || m.species == sp) && m.year <= year) || continue
        spec = m.species != 0
        if m.year > bestyr || (m.year == bestyr && spec && !bestspec)
            val = m.value; bestyr = m.year; bestspec = spec
        end
    end
    return val
end

# THINQFA (option 141, ICFLAG 17): Q-factor diameter-distribution thin — a TWO-record
# keyword (initre.f:5981). Line 1: [year, loDBH, hiDBH, species, Qfactor, classWidth,
# target]; line 2: a single integer for the target units (≤0 BA, ==1 TPA, >1 SDI → 0/1/2).
# Defaults (initre): loDBH=0, hiDBH=24, Q=1.4, classWidth=2. qfatar → activity aux slot.
function kw_thinqfa!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    v = rec.values; pr = rec.present
    yr     = nint(v[1])
    valmin = pr[2] ? Float32(v[2]) : 0f0
    valmax = pr[3] ? Float32(v[3]) : 24f0
    spec   = pr[4] ? Float32(v[4]) : 0f0
    qfac   = pr[5] ? Float32(v[5]) : 1.4f0
    diacw  = pr[6] ? Float32(v[6]) : 2f0
    tarqfa = pr[7] ? Float32(v[7]) : 0f0
    iset   = something(tryparse(Int, strip(read_raw_line!(kr))), 0)   # 2nd record: units switch
    qfatar = iset <= 0 ? 0f0 : (iset <= 1 ? 1f0 : 2f0)
    push!(s.control.schedule,
          ScheduledActivity(Int32(yr), Int32(17), (valmin, valmax, spec, qfac, diacw, tarqfa), qfatar))
    return
end

# OPTION — ESTAB…END establishment packet (esin.f). ESTAB opens the packet (field 1 =
# date of disturbance); subsequent PLANT(430)/NATURAL(431) cards each schedule a regen
# activity (year, species, TPA, %survival, age, height, shade); END closes it. Other
# establishment keywords (TALLY/SPROUT/…) are recognized and skipped for now. At END a
# TALLY(427) trigger is scheduled at the disturbance date and IDSDAT→-9999 (esin.f:100-117).
function kw_estab!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    s.estab.active = true
    idsdat = rec.present[1] ? nint(rec.values[1]) : Int32(-1)   # ESTAB date field
    sched = s.control.schedule
    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            break
        elseif k == "PLANT" || k == "NATURAL"
            ic = k == "PLANT" ? Int32(430) : Int32(431)
            v = r.values
            yr   = r.present[1] ? nint(v[1]) : Int32(1)
            sp   = Float32(v[2])
            tpa  = Float32(v[3])
            tpa <= 0f0 && continue                                # esin.f:143 (TPA must be >0)
            surv = (v[4] < 0.001f0 || v[4] > 100f0) ? 100f0 : Float32(v[4])  # esin.f:149
            push!(sched, ScheduledActivity(yr, ic, (sp, tpa, surv, Float32(v[5]), Float32(v[6]), Float32(v[7]))))
        end
        # other establishment keywords (TALLY/SPROUT/…) not yet ported — skipped
    end
    # END processing (esin.f:100-117): schedule the TALLY(427) establishment trigger at
    # the disturbance date, then mark IDSDAT unset so ESNUTR defaults it.
    if idsdat != Int32(-1)
        push!(sched, ScheduledActivity(max(Int32(1), idsdat), Int32(427),
                                       (Float32(idsdat), 0f0, 0f0, 0f0, 0f0, 0f0)))
    end
    s.estab.idsdat = Int32(-9999)
    return
end

# OPTION 15 — STDIDENT (initre.f:862): stand id from the next raw line.
function kw_stdident!(s::StandState, kr::KeywordReader)
    record = rpad(read_raw_line!(kr), 250)[1:250]
    i1 = findfirst(c -> c != ' ', record[1:26])
    if i1 === nothing
        s.plot.stand_id = rpad(" ", 26)
        s.control.title = rpad("", 72)
        return
    end
    i2 = something(findnext(==(' '), record, i1), 27)
    i2 = min(i2, 27)
    nplt = record[i1:min(i2, 26)]
    s.plot.stand_id = length(nplt) <= 26 ? rpad(nplt, 26) : nplt[1:26]
    rest = i2 + 1 <= length(record) ? record[i2+1:end] : ""
    s.control.title = length(rest) >= 72 ? rest[1:72] : rpad(rest, 72)
    return
end

# OPTION 5 — TREEFMT (initre.f:646): custom tree-record format, 2×80 raw lines.
function kw_treefmt!(s::StandState, kr::KeywordReader)
    line1 = read_raw_line!(kr)
    line2 = read_raw_line!(kr)
    s.control.tree_format = rpad(line1, 80)[1:80] * rpad(line2, 80)[1:80]
    return
end

# OPTION — IF / THEN / ENDIF (event monitor, evmon.f). Reads the condition expression
# (raw lines up to THEN), then the activity keywords up to ENDIF, into a
# ConditionalActivity whose acts fire only in cycles where the condition is true.
# The IF-block activities carry NO date (field 1 blank) — their year is filled in at
# fire time — so they parse with the same (year=v[1]=0, params=v[2:7]) layout.
# OPTION 33 — COMPUTE (vbase/initre.f:1266 → EVUSRV): an event-monitor user-variable block.
# Field 1 = start year (IDT, default 1). The following raw lines are `NAME = expression`
# assignments, re-evaluated each cycle from the start year (a later def may use an earlier one);
# the block ends at END. The values become variables readable by IF/THEN conditions.
function kw_compute!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    date = rec.present[1] ? nint(rec.values[1]) : Int32(1)
    while !eof(kr.io)
        t = strip(read_raw_line!(kr))
        isempty(t) && continue
        uppercase(t) == "END" && break
        eq = findfirst('=', t)
        eq === nothing && continue
        name = uppercase(strip(t[1:eq-1])); expr = strip(t[eq+1:end])
        (isempty(name) || isempty(expr)) && continue
        push!(s.control.compute_defs, (date, String(name), parse_event_condition(String(expr))))
    end
    return
end

function kw_if!(s::StandState, kr::KeywordReader)
    cond = ""
    while true
        line = strip(read_raw_line!(kr))
        (isempty(line) || startswith(uppercase(line), "THEN")) && (isempty(line) ? continue : break)
        cond *= " " * line
        uppercase(line) == "THEN" && break
    end
    acts = ScheduledActivity[]
    while true
        rec = read_keyword!(kr)
        (rec.status == KW_EOF || rec.status == KW_STOP) && break
        kw = strip(rec.name); isempty(kw) && continue
        uppercase(kw) == "ENDIF" && break
        v = rec.values
        if haskey(_THIN_ICFLAG, kw)
            push!(acts, ScheduledActivity(nint(v[1]), _THIN_ICFLAG[kw],
                                          ntuple(i -> Float32(v[i + 1]), 6)))
        elseif kw == "SPECPREF"
            push!(acts, ScheduledActivity(nint(v[1]), Int32(201),
                                          ntuple(i -> Float32(v[i + 1]), 6)))
        end
        # unknown keywords inside IF (e.g. COMPUTE) are skipped for now
    end
    push!(s.control.conditionals,
          ConditionalActivity(parse_event_condition(cond), acts, strip(cond)))
    return
end

"""
    process_keywords!(state, kr, base_path) -> Symbol

Process keyword records from `kr` until PROCESS / STOP / EOF. `base_path` is the
keyword file path with the extension stripped (used to locate the `.tre` file for
TREEDATA). Returns the terminating reason (:process, :stop, :eof).
"""
function process_keywords!(s::StandState, kr::KeywordReader, base_path::AbstractString)
    trees_loaded = false   # an explicit TREEDATA was processed
    notrees      = false   # NOTREES suppresses the default tree read
    nkw          = 0       # real keywords seen (0 ⇒ bare STOP/EOF, not a stand)
    # INITRE end (initre.f:334-336): if no TREEDATA keyword ran, read the tree file once
    # anyway — so a stand without TREEDATA (e.g. snt01's FFE stand, which only REWINDs the
    # shared data unit) still gets its trees. NOTREES (bare stands) suppresses this, and a
    # bare terminator (STOP/EOF with no keywords) is not a stand at all (nkw==0).
    finish(reason) = (nkw > 0 && !trees_loaded && !notrees &&
                      load_trees!(s, base_path * ".tre"); reason)
    while true
        rec = read_keyword!(kr)
        rec.status == KW_EOF && return finish(:eof)
        rec.status == KW_STOP && return finish(:stop)
        kw = strip(rec.name)
        isempty(kw) && continue                      # blank-line record
        nkw += 1
        if     kw == "DESIGN";   kw_design!(s, rec)
        elseif kw == "NUMCYCLE"; kw_numcycle!(s, rec)
        elseif kw == "TIMEINT";  kw_timeint!(s, rec)      # cycle length (period); default 5
        elseif kw == "INVYEAR";  kw_invyear!(s, rec)
        elseif kw == "SITECODE"; kw_sitecode!(s, rec)
        elseif kw == "STDINFO";  kw_stdinfo!(s, rec)
        elseif kw == "STDIDENT"; kw_stdident!(s, kr)
        elseif kw == "TREEFMT";  kw_treefmt!(s, kr)
        elseif kw == "IF";       kw_if!(s, kr)
        elseif kw == "ESTAB";    kw_estab!(s, rec, kr)
        elseif kw == "TREEDATA"; load_trees!(s, base_path * ".tre"); trees_loaded = true
        elseif kw == "NOTREES";  notrees = true       # bare stand — no tree-data read
        elseif kw == "THINQFA"; kw_thinqfa!(s, rec, kr)   # 2-record keyword
        elseif kw == "SPGROUP"; kw_spgroup!(s, rec, kr)   # species group: name + next-record species list
        elseif kw == "COMPUTE"; kw_compute!(s, rec, kr)   # event-monitor variable block (NAME = expr … END)
        elseif kw == "NOTRIPLE"; s.control.icl4 = Int32(0)            # disable record tripling (initre.f:5500)
        elseif kw == "NUMTRIP";  rec.present[1] && (s.control.icl4 = Int32(nint(rec.values[1])))  # set ICL4 (initre.f:2709)
        elseif haskey(_THIN_ICFLAG, kw); kw_thin!(s, rec, _THIN_ICFLAG[kw])
        elseif kw == "SPECPREF"; kw_thin!(s, rec, Int32(201))   # cut modifier: species preference
        elseif kw == "SETPTHIN"; kw_thin!(s, rec, Int32(248))   # point-thin prescription (point, metric)
        elseif kw == "THINPT";   kw_thin!(s, rec, Int32(15))    # point thin (residual + class + dir)
        elseif kw == "BAIMULT";  kw_mult!(s, rec, :bai)   # diameter-growth (BA-increment) multiplier
        elseif kw == "HTGMULT";  kw_mult!(s, rec, :htg)   # height-growth multiplier
        elseif kw == "MORTMULT"; kw_mult!(s, rec, :mort)  # mortality-rate multiplier
        elseif kw == "REGHMULT"; kw_mult!(s, rec, :regh)  # regen height-growth multiplier
        elseif kw == "REGDMULT"; kw_mult!(s, rec, :regd)  # regen diameter-growth multiplier
        elseif kw == "TREESZCP"; kw_treeszcp!(s, rec)     # per-species size cap (SIZCAP)
        elseif kw == "FIXDG";    kw_mult!(s, rec, :fixdg)  # one-shot DG scaler (grincr.f:451)
        elseif kw == "FIXHTG";   kw_mult!(s, rec, :fixhtg) # one-shot HTG scaler (grincr.f:492)
        elseif kw == "HTGSTOP";  kw_htgstp!(s, rec, 110)   # top-damage: scale height growth
        elseif kw == "TOPKILL";  kw_htgstp!(s, rec, 111)   # top-damage: top-kill (htgstp.f)
        elseif kw == "FIXMORT";  kw_fixmort!(s, rec)       # forced-mortality override (morts.f:781)
        elseif kw == "VOLUME";   kw_volume!(s, rec)        # cubic merch-standard override (volkey.f)
        elseif kw == "BFVOLUME"; kw_bfvolume!(s, rec)      # board-foot merch-standard override (volkey.f)
        elseif kw == "MCDEFECT"; kw_mcdefect!(s, rec)      # cubic-volume defect curve (sdefet.f/vols.f)
        elseif kw == "BFDEFECT"; kw_bfdefect!(s, rec)      # board-foot defect curve (sdefet.f/vols.f)
        elseif kw == "CRNMULT";  kw_mult!(s, rec, :crn)    # crown-ratio-change multiplier (crown.f:319)
        elseif kw == "PROCESS";  return finish(:process)
        elseif kw in KNOWN_NOOP
            # recognized, no cycle-0 effect yet
        else
            # unrecognized keyword — ignored for now (handlers added per chunk)
        end
    end
end
