# =============================================================================
# fia_database.jl — FVS DBS database INPUT path (DATABASE / DSNIN / StandSQL / TreeSQL)
#
# Faithful port of the FVS "FVS-ready" database reader, so FVSjl consumes the SAME
# keyfile live FVS reads (a DATABASE input block) and simulates an FIA stand
# byte-identically:
#   vdbsqlite/dbsstandin.f  — FVS_STANDINIT_{COND,PLOT} columns → stand/plot state
#   dbsqlite/dbstreesin.f   — FVS_TREEINIT_{COND,PLOT}  columns → tree records
#   base/notre.f            — TPA expansion of the raw TREE_COUNT (already `notre!`)
#
# The stand SQL row plays the role of the STDINFO/SITECODE/DESIGN/GROWTH cards; the
# tree SQL rows play the role of TREEDATA. TREE_COUNT is the raw per-plot count
# (PROB) — `notre!` expands it to trees/acre from BAF/FPA/PI exactly as FVS does.
# =============================================================================

using SQLite
using DBInterface

# Run `sql` (with %StandID% / %Stand_CN% substituted by `sid`) and return every row
# as a Dict{String,Any} keyed by UPPERCASE column name; SQLite NULL stays `missing`.
function _fia_rows(db::SQLite.DB, sql::AbstractString, sid::AbstractString)
    q = sql
    for tok in ("%StandID%", "%Stand_ID%", "%StandCN%", "%Stand_CN%", "%STANDID%")
        q = replace(q, tok => sid)
    end
    out = Dict{String,Any}[]
    for r in DBInterface.execute(db, q)
        d = Dict{String,Any}()
        for nm in propertynames(r)
            d[uppercase(String(nm))] = getproperty(r, nm)
        end
        push!(out, d)
    end
    return out
end

_fia_present(d, k) = haskey(d, k) && d[k] !== missing && d[k] !== nothing
_fia_f32(d, k, dv) = _fia_present(d, k) ? Float32(d[k]) : dv
_fia_int(d, k, dv) = _fia_present(d, k) ? round(Int, Float64(d[k])) : dv
_fia_str(d, k, dv) = _fia_present(d, k) ? strip(string(d[k])) : dv

# FIA numeric species codes arrive un-padded from the DB (e.g. "71"), but the FVS species
# tables key on the 3-digit FIA code ("071"). Zero-pad a purely-numeric 1–2 char code to 3
# digits (dbstreesin.f:353-360 auto-pads a 2-char numeric code with a leading '0'); leave
# alpha codes untouched. Without this, every FIA code < 100 mis-resolves to Other-Hardwood.
_fia_spcode(c::AbstractString) = (s = strip(c); (!isempty(s) && all(isdigit, s) && length(s) <= 2) ? lpad(s, 3, '0') : String(s))

"""
    apply_fia_stand!(s, d)

Map one FVS_STANDINIT_COND/PLOT row `d` (uppercase-keyed Dict) into the stand's
plot/control state — the STDINFO/SITECODE/DESIGN/GROWTH-card equivalent (dbsstandin.f).
"""
function apply_fia_stand!(s::StandState, d::Dict{String,Any})
    p = s.plot; c = s.control
    # INV_YEAR (IY(1)) → first-cycle (start) year
    iy = _fia_int(d, "INV_YEAR", 0); iy > 0 && (c.cycle_year[1] = Int32(iy))
    # LOCATION (KODFOR): direct if present, else composite REGION*100 + FOREST (dbsstandin.f:569)
    loc = _fia_int(d, "LOCATION", 0)
    (loc == 0 && _fia_present(d, "REGION")) && (loc = _fia_int(d, "REGION", 0) * 100 + _fia_int(d, "FOREST", 0))
    loc != 0 && (p.user_forest_code = Int32(loc))
    # Fort Bragg (forkod.f CASE 701): a region-7/forest-1 FIA code (composite LOCATION=701) is remapped to
    # NC Uwharrie 81110 (region 8) — forkod runs for every stand, incl. DB input. Without it VOLEQDEF sees
    # region 7 ⇒ no R8 Clark equation ⇒ zero volume. (Shared with kw_stdinfo!.)
    s.variant isa Southern && sn_fortbragg_remap!(p)
    # AGE (IAGE)
    _fia_present(d, "AGE") && (p.stand_age = Int32(_fia_int(d, "AGE", 0)))
    # ASPECT degrees → radians (TRNASP ×0.0174533); SLOPE percent → fraction (PLOT.F77)
    _fia_present(d, "ASPECT") && (p.aspect = _fia_f32(d, "ASPECT", 0f0) * 0.0174533f0)
    _fia_present(d, "SLOPE")  && (p.slope  = _fia_f32(d, "SLOPE", 0f0) / 100f0)
    # ELEVATION in hundreds of feet; ELEVFT is feet → ×0.01 (dbsstandin.f:710)
    if _fia_present(d, "ELEVFT")
        p.elevation = _fia_f32(d, "ELEVFT", 0f0) * 0.01f0
    elseif _fia_present(d, "ELEVATION")
        p.elevation = _fia_f32(d, "ELEVATION", p.elevation)
    end
    # LATITUDE/LONGITUDE (TLAT/TLONG, dbsstandin.f:254-259) — feed the Hopkins bioclimatic
    # index in the eastern crown-width models.
    _fia_present(d, "LATITUDE")  && (p.latitude  = _fia_f32(d, "LATITUDE", p.latitude))
    _fia_present(d, "LONGITUDE") && (p.longitude = _fia_f32(d, "LONGITUDE", p.longitude))
    # ECOREGION (ecological unit / EUT, e.g. "223Db") → eco_unit. FVS reads it from STANDINIT and adds a
    # per-species ecological-unit DG term (dgf.f EUT categorical coefficients dg_phys_*), plus it drives the
    # montane site/height/estab branches (eco_unit[1]=='M'). Without it the whole EUT DG term is dropped for
    # FIA-DB stands — the FIA/FVS campaign's slice-1 divergence (yellow-poplar large-tree DG ~20% low on a
    # 223Db stand: jl omitted dg_phys_p222 ≈ +0.255 of the −0.344 ln(DDS) deficit). SN-gated: resolve_eco_unit
    # / SNECU is the SOUTHERN ecological-unit table; NE/CS/LS use their own eco-unit handling (left blank as
    # before, a documented follow-up if their FIA differentials show an analogous gap).
    if s.variant isa Southern && _fia_present(d, "ECOREGION")
        p.eco_unit = rpad(resolve_eco_unit(_fia_str(d, "ECOREGION", ""), 0), 10)
    end
    # FORKOD phase-3 default (forkod.f:540-546, mirrored from kw_stdinfo!): fill any geo field the
    # DB left at 0 from the national-forest table. FVS runs forkod BEFORE the DB overrides, and the
    # DB overrides elevation ONLY when >0 (dbsstandin.f:647) — so a null/≤0 ELEVATION keeps the
    # forest default (e.g. LOCATION 80215 → forest 802 → 12.0 hundred-ft, live-confirmed). That
    # elevation drives the Hopkins index for hardwood open-grown crowns, so without it HI (and the
    # reported CCF) drift. Southern-gated: the SN forest_location table is keyed by KODFOR÷100 the
    # same way kw_stdinfo! keys it; NE/CS/LS use a different forkod keying (left as a follow-up).
    if s.variant isa Southern
        lat0, long0, elev0 = forest_location(s.coef, div(Int(p.user_forest_code), 100))
        p.latitude  == 0f0 && (p.latitude  = lat0)
        p.longitude == 0f0 && (p.longitude = long0)
        p.elevation == 0f0 && (p.elevation = elev0)
    end
    # Sampling design (DESIGN card): BAF / FPA / BRK / IPTINV / NONSTK / SAMWT / GROSPC
    _fia_present(d, "BASAL_AREA_FACTOR") && (p.baf = _fia_f32(d, "BASAL_AREA_FACTOR", 0f0))
    _fia_present(d, "INV_PLOT_SIZE")     && (p.fixed_plot_inv = _fia_f32(d, "INV_PLOT_SIZE", 0f0))
    _fia_present(d, "BRK_DBH")           && (p.min_dbh_var_plot = _fia_f32(d, "BRK_DBH", p.min_dbh_var_plot))
    _fia_present(d, "NUM_PLOTS")         && (p.points_inv = Int32(_fia_int(d, "NUM_PLOTS", 1)))
    _fia_present(d, "NONSTK_PLOTS")      && (p.nonstockable = Int32(_fia_int(d, "NONSTK_PLOTS", 0)))
    _fia_present(d, "SAM_WT")            && (p.sample_weight = _fia_f32(d, "SAM_WT", p.sample_weight))
    # GROSPC: STK_PCNT (1..100 → ÷100) if given, else (IPTINV − NONSTK)/IPTINV (dbsstandin.f:740)
    if _fia_present(d, "STK_PCNT")
        g = _fia_f32(d, "STK_PCNT", 1f0)
        (g > 1f0 && g <= 100f0) && (g *= 0.01f0)
        (g > 0f0 && g <= 1f0) && (p.gross_space = g)
    else
        ip = max(1, Int(p.points_inv))
        p.gross_space = Float32(ip - Int(p.nonstockable)) / Float32(ip)
    end
    # Growth calibration transition/measurement (GROWTH card: IDG/FINT/IHTG/FINTH/FINTM).
    # DG_TRANS=1 ⇒ the DG field is a PAST diameter (not an increment) measured DG_MEASURE yrs ago.
    _fia_present(d, "DG_TRANS")     && (c.growth_idg   = Int32(_fia_int(d, "DG_TRANS", 0)))
    _fia_present(d, "DG_MEASURE")   && (c.growth_fint  = _fia_f32(d, "DG_MEASURE", 5f0))
    _fia_present(d, "HTG_TRANS")    && (c.growth_ihtg  = Int32(_fia_int(d, "HTG_TRANS", 0)))
    _fia_present(d, "HTG_MEASURE")  && (c.growth_finth = _fia_f32(d, "HTG_MEASURE", 5f0))
    _fia_present(d, "MORT_MEASURE") && (c.growth_fintm = _fia_f32(d, "MORT_MEASURE", 5f0))
    # SITE_SPECIES (ISISP) + SITE_INDEX (SITEAR): assign to the site species only if given,
    # else to all species (dbsstandin.f:841). ≤7 = Dunning code (not yet handled → direct).
    if _fia_present(d, "SITE_INDEX")
        si = _fia_f32(d, "SITE_INDEX", 0f0)
        isp = 0
        if _fia_present(d, "SITE_SPECIES")
            code = _fia_spcode(_fia_str(d, "SITE_SPECIES", ""))
            if !isempty(code)
                idx, _ = resolve_species(code, s.variant, s.species, s.coef)
                isp = Int(idx)
            end
        end
        if isp >= 1
            p.sp_site_index[isp] = si; p.site_species = Int32(isp); p.site_index = si
        else
            fill!(p.sp_site_index, si); p.site_index = si
        end
    end
    return s
end

"""
    apply_fia_trees!(s, rows) -> Int

Map FVS_TREEINIT_COND/PLOT rows into `TreeRecord`s (dbstreesin.f column reads) and
ingest them via the shared `ingest_tree_records!` (same fixups as the .tre loader).
TREE_COUNT is stored raw as `tpa` (PROB); `notre!` later expands it to trees/acre.
"""
function apply_fia_trees!(s::StandState, rows::Vector{Dict{String,Any}})
    recs = TreeRecord[]
    for d in rows
        dbh = _fia_present(d, "DIAMETER") ? _fia_f32(d, "DIAMETER", 0f0) : _fia_f32(d, "DBH", 0f0)
        dmg = (Int32(_fia_int(d, "DAMAGE1", 0)), Int32(_fia_int(d, "SEVERITY1", 0)),
               Int32(_fia_int(d, "DAMAGE2", 0)), Int32(_fia_int(d, "SEVERITY2", 0)),
               Int32(_fia_int(d, "DAMAGE3", 0)), Int32(_fia_int(d, "SEVERITY3", 0)))
        rec = TreeRecord(
            Int32(_fia_int(d, "PLOT_ID", 1)),           # plot (ITREI) → subplot/IPVEC
            Int32(_fia_int(d, "TREE_ID", 0)),           # id (IDTREE)
            _fia_f32(d, "TREE_COUNT", 1f0),             # tpa (raw PROB; notre! expands)
            Int32(_fia_int(d, "HISTORY", 1)),           # history (ITH; 1 = live default)
            _fia_spcode(_fia_str(d, "SPECIES", "OT")),  # species_code (FIA 3-digit / alpha / PLANTS)
            dbh,                                         # dbh
            _fia_f32(d, "DG", 0f0),                     # diam_growth (PAST dbh when IDG=1)
            _fia_f32(d, "HT", 0f0),                     # height
            _fia_f32(d, "HTTOPK", 0f0),                 # top_height (broken/dead)
            _fia_f32(d, "HTG", 0f0),                    # ht_growth
            Int32(_fia_int(d, "CRRATIO", 0)),           # crown_pct (ICR)
            dmg,
            Int32(_fia_int(d, "TREEVALUE", 0)),         # mort_code (IMC1)
            Int32(_fia_int(d, "PRESCRIPTION", 0)),      # cut_code (KUTKOD)
            (Int32(0), Int32(0), Int32(0), Int32(0), Int32(0)),  # pest_vars
            _fia_f32(d, "AGE", 0f0),                    # birth_age (ABIRTH)
        )
        push!(recs, rec)
    end
    return ingest_tree_records!(s, recs)
end

"""
    load_fia_stand!(s, dbpath, standsql, treesql) -> StandState

Populate stand `s` from an FIA "FVS-ready" SQLite database: run `standsql` (one row →
stand/plot state) and `treesql` (tree records), substituting `%StandID%` with the
stand's id (from STDIDENT). Mirrors the FVS DATABASE/DSNIN input block.
"""
function load_fia_stand!(s::StandState, dbpath::AbstractString,
                         standsql::AbstractString, treesql::AbstractString)
    sid = String(strip(s.plot.stand_id))
    # Open READ-ONLY (URI mode=ro, immutable=1): FVSjl only ever SELECTs from an FIA
    # database — it must never create/modify/journal the source file.
    db = SQLite.DB(startswith(dbpath, "file:") ? dbpath : "file:$(dbpath)?mode=ro&immutable=1")
    try
        srows = _fia_rows(db, standsql, sid)
        isempty(srows) && error("FIA database: no FVS_STANDINIT row for stand '$sid'")
        apply_fia_stand!(s, srows[1])
        isempty(treesql) || apply_fia_trees!(s, _fia_rows(db, treesql, sid))
    finally
        SQLite.close(db)
    end
    return s
end
