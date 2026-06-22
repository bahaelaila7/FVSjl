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
    "NOTREES", "NOTRIPLE", "NOHTDREG", "AUTOES",
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
                          "THINSDI" => Int32(10), "THINHT" => Int32(12))

# Parse a THIN* activity: field 1 = calendar year, fields 2-7 = the 6 method params.
# Stores a ScheduledActivity for `cuts!` to apply on the matching cycle.
function kw_thin!(s::StandState, rec::KeywordRecord, icflag::Int32)
    v = rec.values
    yr = nint(v[1])
    params = ntuple(i -> Float32(v[i + 1]), 6)
    push!(s.control.schedule, ScheduledActivity(Int32(yr), icflag, params))
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
        elseif kw == "INVYEAR";  kw_invyear!(s, rec)
        elseif kw == "SITECODE"; kw_sitecode!(s, rec)
        elseif kw == "STDINFO";  kw_stdinfo!(s, rec)
        elseif kw == "STDIDENT"; kw_stdident!(s, kr)
        elseif kw == "TREEFMT";  kw_treefmt!(s, kr)
        elseif kw == "IF";       kw_if!(s, kr)
        elseif kw == "ESTAB";    kw_estab!(s, rec, kr)
        elseif kw == "TREEDATA"; load_trees!(s, base_path * ".tre"); trees_loaded = true
        elseif kw == "NOTREES";  notrees = true       # bare stand — no tree-data read
        elseif haskey(_THIN_ICFLAG, kw); kw_thin!(s, rec, _THIN_ICFLAG[kw])
        elseif kw == "SPECPREF"; kw_thin!(s, rec, Int32(201))   # cut modifier: species preference
        elseif kw == "PROCESS";  return finish(:process)
        elseif kw in KNOWN_NOOP
            # recognized, no cycle-0 effect yet
        else
            # unrecognized keyword — ignored for now (handlers added per chunk)
        end
    end
end
