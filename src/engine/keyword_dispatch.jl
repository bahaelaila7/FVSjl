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
    # SN: STDINFO field 2 is the ecological unit code (PCOM), e.g. "231Dd"
    rec.present[2] && (p.eco_unit = rpad(strip(rec.fields[2]), 10))
    rec.present[3] && (p.stand_age = nint(v[3]))
    rec.present[4] && (p.aspect = v[4] * 0.0174533f0)   # degrees → radians (utils.f)
    rec.present[5] && (p.slope  = v[5] / 100f0)         # percent → fraction (utils.f)
    (rec.present[6] && v[6] > 0f0) && (p.elevation = v[6])
    if rec.present[9]
        org = nint(v[9])
        p.stand_origin = (org < 0 || org > 1) ? Int32(0) : org
    end
    # FORKOD phase 3: default lat/long/elev from the forest code (forkod.f:193).
    lat0, long0, elev0 = forest_location(div(p.user_forest_code, 100))
    p.latitude  == 0f0 && (p.latitude  = lat0)
    p.longitude == 0f0 && (p.longitude = long0)
    p.elevation == 0f0 && (p.elevation = elev0)
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

"""
    process_keywords!(state, kr, base_path) -> Symbol

Process keyword records from `kr` until PROCESS / STOP / EOF. `base_path` is the
keyword file path with the extension stripped (used to locate the `.tre` file for
TREEDATA). Returns the terminating reason (:process, :stop, :eof).
"""
function process_keywords!(s::StandState, kr::KeywordReader, base_path::AbstractString)
    while true
        rec = read_keyword!(kr)
        rec.status == KW_EOF && return :eof
        rec.status == KW_STOP && return :stop
        kw = strip(rec.name)
        isempty(kw) && continue                      # blank-line record
        if     kw == "DESIGN";   kw_design!(s, rec)
        elseif kw == "NUMCYCLE"; kw_numcycle!(s, rec)
        elseif kw == "INVYEAR";  kw_invyear!(s, rec)
        elseif kw == "SITECODE"; kw_sitecode!(s, rec)
        elseif kw == "STDINFO";  kw_stdinfo!(s, rec)
        elseif kw == "STDIDENT"; kw_stdident!(s, kr)
        elseif kw == "TREEFMT";  kw_treefmt!(s, kr)
        elseif kw == "TREEDATA"; load_trees!(s, base_path * ".tre")
        elseif kw == "PROCESS";  return :process
        elseif kw in KNOWN_NOOP
            # recognized, no cycle-0 effect yet
        else
            # unrecognized keyword — ignored for now (handlers added per chunk)
        end
    end
end
