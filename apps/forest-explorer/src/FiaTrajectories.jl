"""
    FiaTrajectories

Optional "FIA plots — tree trajectories" layer. Serves
  * `GET  /api/fiaplots/options`        distinct filter values, axes, colour-by choices (builds the plot summary on first use)
  * `POST /api/fiaplots/points`         plots in a bbox that match the filters (dots); also `GET` with query parameters
  * `POST /api/fiaplots/trajectories`   per-tree series of one plot, built as in FIA_DATA_PREP/tree_trajectories.py

Data: table `tree_trajectories` (one row per tree visit) of the curated FIA DuckDB plus PLOT (LAT/LON). The file is opened
READ-ONLY and never written. Path: ENV `FIA_TRAJ_DB`, else `PAN_FIA_DB`, else the workspace default.

Filter semantics follow the notebook: a tree is kept if ANY of its rows satisfies ALL selected filters (so its FULL
trajectory is plotted). The dot filter uses per-plot value sets (a plot is shown when, for every filter, it has at least
one row with a selected value) — a necessary condition; the exact test is applied when a plot's trajectories are loaded.
"""
module FiaTrajectories

using DuckDB, DBInterface, DataFrames, JSON3, HTTP, Oxygen

const DEFAULT_DB = joinpath(get(ENV, "WORKSPACE", normpath(joinpath(@__DIR__, "..", "..", "..", ".."))), "FIA", "truth", "fia_truth_r4.duckdb")   # <workspace root>/FIA/...
db_path() = get(ENV, "FIA_TRAJ_DB", get(ENV, "PAN_FIA_DB", DEFAULT_DB))

const MAX_POINTS = 5000
const OWN_LABELS = Dict(10 => "USFS", 20 => "Other federal", 30 => "State/local", 40 => "Private")
# selectable axes: id => (column, axis label, unit)  (same set as the notebook; X and Y share it)
const AXES = [
    ("drybio_g",   "Biomass (g)",     "AG dry biomass (g / tree)", "g"),
    ("age_calc",   "Age — Selected",  "Age — selected (yr)",       "yr"),
    ("age_totage", "Age — TOTAGE",    "Age — TOTAGE (yr)",         "yr"),
    ("age_stdage", "Age — STDAGE",    "Age — STDAGE (yr)",         "yr"),
    ("age_estage", "Age — estage",    "Age — estage (yr)",         "yr"),
    ("age_model",  "Age — ML model",  "Age — ML model (yr)",       "yr"),
    ("dia",        "Diameter (in)",   "Diameter (in)",             "in"),
    ("ht",         "Height (ft)",     "Height (ft)",               "ft"),
]
const AXIS_COLS = Dict(a[1] => a for a in AXES)
const ECO_FIELDS = [("epa_l3", "EPA L3"), ("epa_l4", "EPA L4"), ("ecosubcd", "Ecosubcd")]
const COLOR_BY = ["species_symbol", "land_use", "terminal_fate", "cut_event_type",
                  "owngrpcd", "sftwd_hrdwd", "cycle", "age_source", "spcd_changed"]
const REPRS = ["all", "first", "last"]

# --- per-plot summary (cached in memory) --------------------------------------------------------------------------------
struct Summary
    key::Vector{String}
    lat::Vector{Float64}
    lon::Vector{Float64}
    n_trees::Vector{Int32}
    n_visits::Vector{Int16}
    y0::Vector{Int16}
    y1::Vector{Int16}
    species::Vector{Vector{String}}    # all species symbols of the plot, most trees first
    sp_n::Vector{Vector{Int32}}
    sp_ids::Vector{Vector{UInt16}}     # sorted ids into `sp_names`
    lu::Vector{UInt8}                   # bit masks (index into the option lists below)
    ce::Vector{UInt8}
    tf::Vector{UInt16}                  # 13 terminal fates (curated censoring / land-use classes)
    og::Vector{UInt8}
    sh::Vector{UInt8}
    cy::Vector{UInt32}
    epa_l3::Vector{String}
    epa_l4::Vector{String}
    ecosubcd::Vector{String}
    sp_names::Vector{String}
    lu_opts::Vector{String}
    ce_opts::Vector{String}
    tf_opts::Vector{String}
    og_opts::Vector{Int}
    sh_opts::Vector{String}
    cy_opts::Vector{Int}
    eco_opts::Dict{String,Vector{String}}
    lon_order::Vector{Int}              # plot indices sorted by lon (bbox pre-selection)
    lons_sorted::Vector{Float64}
    index::Dict{String,Int}
end

const _SUMMARY = Ref{Union{Nothing,Summary}}(nothing)
const _LOCK = ReentrantLock()          # guards the summary build and all DuckDB access
const _DB = Ref{Any}(nothing)
const _CON = Ref{Any}(nothing)

function _con()
    if _CON[] === nothing
        p = db_path()
        isfile(p) || error("FIA trajectories DB not found: $p (set FIA_TRAJ_DB)")
        _DB[] = DuckDB.DB(p; readonly = true)        # read-only: this layer never writes to the FIA database
        _CON[] = DBInterface.connect(_DB[])
    end
    _CON[]
end
_q(sql::AbstractString) = DataFrame(DBInterface.execute(_con(), sql))

_s(x) = (x === missing || x === nothing) ? "" : String(x)
_lst(x) = x === missing ? String[] : String[string(v) for v in x if v !== missing]

function _mask(vals, opts::Vector)
    m = UInt64(0)
    for v in vals
        i = findfirst(==(v), opts)
        i === nothing || (m |= UInt64(1) << (i - 1))
    end
    m
end

function _build_summary()
    t0 = time()
    @info "FIA trajectories: building plot summary" db = db_path()
    sql = """
    WITH sp AS (SELECT plt_cn, species_symbol, count(DISTINCT (subp, tree)) AS n
                FROM tree_trajectories GROUP BY 1, 2),
    spl AS (SELECT plt_cn, list(species_symbol ORDER BY n DESC, species_symbol) AS sps,
                   list(n ORDER BY n DESC, species_symbol) AS spn FROM sp GROUP BY 1),
    pl AS (SELECT plt_cn, any_value(statecd) AS statecd, any_value(unitcd) AS unitcd,
                  any_value(countycd) AS countycd, any_value(plot) AS plot,
                  count(DISTINCT (subp, tree)) AS n_trees, count(DISTINCT measdate) AS n_visits,
                  min(year(measdate)) AS y0, max(year(measdate)) AS y1,
                  list(DISTINCT land_use) AS lu, list(DISTINCT cut_event_type) AS ce,
                  list(DISTINCT terminal_fate) AS tf, list(DISTINCT owngrpcd) AS og,
                  list(DISTINCT sftwd_hrdwd) AS sh, list(DISTINCT cycle) AS cy,
                  any_value(epa_l3) AS epa_l3, any_value(epa_l4) AS epa_l4, any_value(ecosubcd) AS ecosubcd
           FROM tree_trajectories GROUP BY plt_cn),
    xy AS (SELECT statecd, unitcd, countycd, plot, avg(lat) AS lat, avg(lon) AS lon
           FROM PLOT WHERE lat IS NOT NULL AND lon IS NOT NULL GROUP BY ALL)
    SELECT pl.plt_cn, xy.lat, xy.lon, pl.n_trees, pl.n_visits, pl.y0, pl.y1, spl.sps, spl.spn,
           pl.lu, pl.ce, pl.tf, pl.og, pl.sh, pl.cy, pl.epa_l3, pl.epa_l4, pl.ecosubcd
    FROM pl JOIN spl USING (plt_cn) JOIN xy USING (statecd, unitcd, countycd, plot)
    ORDER BY pl.plt_cn"""
    df = _q(sql)
    n = nrow(df)
    uni(col, f = identity) = sort!(unique(vcat([[f(v) for v in _lst(x)] for x in df[!, col]]...)))
    lu_opts = uni(:lu); ce_opts = uni(:ce); tf_opts = uni(:tf); sh_opts = uni(:sh)
    og_opts = sort!(unique(vcat([[round(Int, parse(Float64, v)) for v in _lst(x)] for x in df.og]...)))
    cy_opts = sort!(unique(vcat([[parse(Int, v) for v in _lst(x)] for x in df.cy]...)))
    sp_names = sort!(unique(vcat([_lst(x) for x in df.sps]...)))
    spid = Dict(s => UInt16(i) for (i, s) in enumerate(sp_names))
    species = [_lst(x) for x in df.sps]
    sp_ids = [sort!(UInt16[spid[s] for s in sp]) for sp in species]
    mk(col, opts, conv = identity) = UInt64[_mask([conv(v) for v in _lst(x)], opts) for x in df[!, col]]
    lu = UInt8.(mk(:lu, lu_opts)); ce = UInt8.(mk(:ce, ce_opts)); tf = UInt16.(mk(:tf, tf_opts))
    og = UInt8.(mk(:og, og_opts, v -> round(Int, parse(Float64, v)))); sh = UInt8.(mk(:sh, sh_opts))
    cy = UInt32.(mk(:cy, cy_opts, v -> parse(Int, v)))
    e3 = String[_s(x) for x in df.epa_l3]; e4 = String[_s(x) for x in df.epa_l4]; es = String[_s(x) for x in df.ecosubcd]
    eco_opts = Dict("epa_l3" => sort!(unique(filter(!isempty, e3))), "epa_l4" => sort!(unique(filter(!isempty, e4))),
                    "ecosubcd" => sort!(unique(filter(!isempty, es))))
    keys = String.(df.plt_cn); lat = Float64.(df.lat); lon = Float64.(df.lon)
    ord = sortperm(lon)
    s = Summary(keys, lat, lon, Int32.(df.n_trees), Int16.(df.n_visits), Int16.(df.y0), Int16.(df.y1),
                species, [Int32[Int32(v) for v in x] for x in df.spn], sp_ids, lu, ce, tf, og, sh, cy, e3, e4, es,
                sp_names, lu_opts, ce_opts, tf_opts, og_opts, sh_opts, cy_opts, eco_opts,
                ord, lon[ord], Dict(k => i for (i, k) in enumerate(keys)))
    @info "FIA trajectories: summary ready" plots = n seconds = round(time() - t0; digits = 1)
    s
end

function summary()
    s = _SUMMARY[]
    s === nothing || return s
    lock(_LOCK) do
        _SUMMARY[] === nothing && (_SUMMARY[] = _build_summary())
        _SUMMARY[]
    end
end

# --- filters -------------------------------------------------------------------------------------------------------------
Base.@kwdef struct Filters
    species::Vector{String} = String[]
    land_use::Vector{String} = String[]
    cut_event::Vector{String} = String[]
    fate::Vector{String} = String[]
    own::Vector{Int} = Int[]
    swhw::String = ""
    cycle::Int = 0
    eco_field::String = "epa_l3"
    eco_val::String = ""
end

_strs(x) = x === nothing || x === missing ? String[] :
           x isa AbstractString ? String[strip(s) for s in split(x, ',') if !isempty(strip(s))] :
           String[string(v) for v in x]
_ints(x) = Int[v isa Integer ? Int(v) : parse(Int, string(v)) for v in _strs(x)]
_none_if_all(x) = (x === nothing || x === missing || x in ("", "(all)", "all")) ? "" : String(x)

"Parse filters from a JSON object / Dict-like (`get(o, :key, nothing)`); unknown values are rejected (they are inlined in SQL)."
function parse_filters(o, s::Summary)
    g(k) = o === nothing ? nothing : get(o, k, nothing)
    chk(vals, ok, what) = (bad = [v for v in vals if !(v in ok)]; isempty(bad) || error("unknown $what: $(first(bad))"); vals)
    eco_field = _s(g(:eco_field)); isempty(eco_field) && (eco_field = "epa_l3")
    eco_field in first.(ECO_FIELDS) || error("unknown ecoregion level: $eco_field")
    eco_val = _s(g(:eco_val))
    isempty(eco_val) || eco_val in s.eco_opts[eco_field] || error("unknown ecoregion: $eco_val")
    swhw = _none_if_all(g(:swhw)); isempty(swhw) || swhw in s.sh_opts || error("unknown soft/hardwood: $swhw")
    cyc = _none_if_all(g(:cycle)); cycle = isempty(cyc) ? 0 : parse(Int, string(cyc))
    cycle == 0 || cycle in s.cy_opts || error("unknown cycle: $cycle")
    Filters(; species = chk(_strs(g(:species)), s.sp_names, "species"),
            land_use = chk(_strs(g(:land_use)), s.lu_opts, "land use"),
            cut_event = chk(_strs(g(:cut_event)), s.ce_opts, "cut event"),
            fate = chk(_strs(g(:fate)), s.tf_opts, "terminal fate"),
            own = chk(_ints(g(:own)), s.og_opts, "ownership"),
            swhw, cycle, eco_field, eco_val)
end

_qs(v) = "'" * replace(String(v), "'" => "''") * "'"
_inlist(col, vals) = "$col IN (" * join(_qs.(vals), ", ") * ")"

"SQL WHERE clause (without leading AND) for the row-level filters; empty string if none."
function where_sql(f::Filters)
    w = String[]
    isempty(f.eco_val) || push!(w, "$(f.eco_field) = $(_qs(f.eco_val))")
    isempty(f.species) || push!(w, _inlist("species_symbol", f.species))
    isempty(f.land_use) || push!(w, _inlist("land_use", f.land_use))
    isempty(f.cut_event) || push!(w, _inlist("cut_event_type", f.cut_event))
    isempty(f.fate) || push!(w, _inlist("terminal_fate", f.fate))
    isempty(f.own) || push!(w, "owngrpcd IN (" * join(f.own, ", ") * ")")
    isempty(f.swhw) || push!(w, "sftwd_hrdwd = $(_qs(f.swhw))")
    f.cycle == 0 || push!(w, "cycle = $(f.cycle)")
    join(w, " AND ")
end

function _mask_of(vals, opts)
    m = UInt64(0)
    for v in vals
        m |= UInt64(1) << (findfirst(==(v), opts) - 1)
    end
    m
end

"Indices of the plots inside the bbox that satisfy the per-plot filter sets."
function select_plots(s::Summary, f::Filters, w, so, e, n)
    lo = searchsortedfirst(s.lons_sorted, w); hi = searchsortedlast(s.lons_sorted, e)
    mlu = _mask_of(f.land_use, s.lu_opts); mce = _mask_of(f.cut_event, s.ce_opts); mtf = _mask_of(f.fate, s.tf_opts)
    mog = _mask_of(f.own, s.og_opts); msh = isempty(f.swhw) ? UInt64(0) : _mask_of([f.swhw], s.sh_opts)
    mcy = f.cycle == 0 ? UInt64(0) : _mask_of([f.cycle], s.cy_opts)
    sp = UInt16[UInt16(findfirst(==(x), s.sp_names)) for x in f.species]
    out = Int[]
    for r in lo:hi
        i = s.lon_order[r]
        (so <= s.lat[i] <= n) || continue
        mlu == 0 || (s.lu[i] & mlu) != 0 || continue
        mce == 0 || (s.ce[i] & mce) != 0 || continue
        mtf == 0 || (s.tf[i] & mtf) != 0 || continue
        mog == 0 || (s.og[i] & mog) != 0 || continue
        msh == 0 || (s.sh[i] & msh) != 0 || continue
        mcy == 0 || (s.cy[i] & mcy) != 0 || continue
        if !isempty(f.eco_val)
            v = f.eco_field == "epa_l3" ? s.epa_l3[i] : f.eco_field == "epa_l4" ? s.epa_l4[i] : s.ecosubcd[i]
            v == f.eco_val || continue
        end
        if !isempty(sp)
            ids = s.sp_ids[i]
            any(x -> insorted(x, ids), sp) || continue
        end
        push!(out, i)
    end
    out
end

# --- options ---------------------------------------------------------------------------------------------------------------
function options()
    s = summary()
    cnt = Dict{String,Int}()
    for sp in s.species, x in sp
        cnt[x] = get(cnt, x, 0) + 1
    end
    (; n_plots = length(s.key),
       species = [(; id = x, plots = cnt[x]) for x in s.sp_names],
       land_use = s.lu_opts, cut_event = s.ce_opts, fate = s.tf_opts,
       own = [(; id = c, label = get(OWN_LABELS, c, string(c))) for c in s.og_opts],
       swhw = s.sh_opts, cycle = s.cy_opts,
       eco_fields = [(; id = id, label = lab, values = s.eco_opts[id]) for (id, lab) in ECO_FIELDS],
       axes = [(; id = a[1], label = a[2], title = a[3], unit = a[4]) for a in AXES],
       color_by = COLOR_BY, reprs = REPRS, max_points = MAX_POINTS, db = db_path(),
       defaults = (; x = "age_calc", y = "drybio_g", color_by = "species_symbol", max_traj = 300,
                   repr = "all", imputed = true))
end

# --- points ----------------------------------------------------------------------------------------------------------------
function points(body)
    s = summary()
    f = parse_filters(body, s)
    g(k) = get(body, k, nothing)
    w = Float64(g(:west)); so = Float64(g(:south)); e = Float64(g(:east)); n = Float64(g(:north))
    idx = select_plots(s, f, w, so, e, n)
    total = length(idx)
    truncated = total > MAX_POINTS
    if truncated      # evenly spaced (deterministic) subsample keeps the spatial distribution
        step = total / MAX_POINTS
        idx = idx[unique!(round.(Int, 1:step:total))]
    end
    pts = [(; key = s.key[i], lat = round(s.lat[i]; digits = 5), lon = round(s.lon[i]; digits = 5),
              n_trees = Int(s.n_trees[i]), n_visits = Int(s.n_visits[i]),
              species = s.species[i][1:min(5, end)]) for i in idx]
    (; n_match = total, shown = length(pts), truncated, cap = MAX_POINTS, points = pts)
end

# --- trajectories ----------------------------------------------------------------------------------------------------------
const TAB10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
# Google's polynomial approximation of the "turbo" colormap (matplotlib's turbo sampled at t in [0,1])
function turbo_hex(t::Float64)
    t = clamp(t, 0.0, 1.0)
    r = 0.13572138 + t * (4.61539260 + t * (-42.66032258 + t * (132.13108234 + t * (-152.94239396 + t * 59.28637943))))
    gg = 0.09140261 + t * (2.19418839 + t * (4.84296658 + t * (-14.18503333 + t * (4.27729857 + t * 2.82956604))))
    b = 0.10667330 + t * (12.64194608 + t * (-60.58204836 + t * (110.36276771 + t * (-89.90310912 + t * 27.34824973))))
    c(x) = string(lpad(string(clamp(round(Int, 255x), 0, 255); base = 16), 2, '0'))
    "#" * c(r) * c(gg) * c(b)
end
const GREY = "#808080"

const FATE_MARK = Dict("alive_censored" => ("circle", false), "dead_natural" => ("triangle-down", true),
                       "dead_disturbance" => ("cross", true), "silvicultural_death" => ("square", true),
                       "harvest_inferred" => ("square", false), "disappeared_inferred_removal" => ("diamond", false),
                       "diverted_landuse" => ("diamond", true))   # censored_* -> the default open circle, like alive_censored

_r(x) = round(Float64(x); sigdigits = 7)

function _levelstr(v)
    v === missing && return nothing
    v isa Bool && return v ? "true" : "false"
    v isa AbstractFloat && isinteger(v) && return string(Int(v))
    string(v)
end

function _fit(xs::Vector{Float64}, ys::Vector{Float64}, logx::Bool, logy::Bool)
    length(xs) >= 2 || return nothing
    tx = logx ? log10.(xs) : copy(xs); ty = logy ? log10.(ys) : copy(ys)
    ok = isfinite.(tx) .& isfinite.(ty)
    tx = tx[ok]; ty = ty[ok]
    (length(tx) >= 2 && maximum(tx) > minimum(tx)) || return nothing
    mx = sum(tx) / length(tx); my = sum(ty) / length(ty)
    sxx = sum((tx .- mx) .^ 2); b = sum((tx .- mx) .* (ty .- my)) / sxx; a = my - b * mx
    ssr = sum((ty .- (b .* tx .+ a)) .^ 2); sst = sum((ty .- my) .^ 2)
    r2 = sst > 0 ? 1 - ssr / sst : NaN
    gx = range(minimum(tx), maximum(tx); length = 100); gy = b .* gx .+ a
    (; slope = b, intercept = a, r2 = isfinite(r2) ? r2 : nothing, n = length(tx), logx, logy,
       x = [_r(logx ? 10.0^v : v) for v in gx], y = [_r(logy ? 10.0^v : v) for v in gy])
end

function trajectories(body)
    s = summary()
    key = _s(get(body, :plot, nothing))
    i = get(s.index, key, 0)
    i == 0 && error("unknown plot: $key")
    f = parse_filters(get(body, :filters, nothing), s)
    o = get(body, :options, nothing)
    og(k, d) = (o === nothing || get(o, k, nothing) === nothing) ? d : get(o, k, d)
    xcol = String(og(:x, "age_calc")); ycol = String(og(:y, "drybio_g"))
    haskey(AXIS_COLS, xcol) && haskey(AXIS_COLS, ycol) || error("unknown axis")
    cby = String(og(:color_by, "species_symbol")); cby in COLOR_BY || error("unknown colour-by: $cby")
    repr = String(og(:repr, "all")); repr in REPRS || error("unknown points mode: $repr")
    max_traj = clamp(Int(og(:max_traj, 300)), 1, 5000)
    imputed = Bool(og(:imputed, true)); deriv = Bool(og(:derivative, false))
    logx = Bool(og(:logx, false)); logy = Bool(og(:logy, false)); dofit = Bool(og(:fit, false))

    where_f = where_sql(f)
    where_all = "plt_cn = $(_qs(key))" * (isempty(where_f) ? "" : " AND " * where_f)
    df, n_match = lock(_LOCK) do
        nm = _q("SELECT count(*) AS n FROM (SELECT DISTINCT subp, tree FROM tree_trajectories WHERE $where_all)").n[1]
        d = _q("""
        WITH ids AS (
            SELECT subp, tree FROM (SELECT DISTINCT subp, tree FROM tree_trajectories WHERE $where_all)
            USING SAMPLE reservoir($max_traj ROWS) REPEATABLE(42))
        SELECT t.subp, t.tree, t.measdate, t.age_calc, t.age_totage, t.age_stdage, t.age_estage, t.age_model, t.dia, t.ht,
               t.age_source, t.drybio_g, t.species_symbol, t.spcd_changed, t.intro_type, t.is_terminal, t.terminal_fate,
               t.land_use, t.cut_event_type, t.owngrpcd, t.sftwd_hrdwd, t.cycle
        FROM tree_trajectories t JOIN ids USING (subp, tree)
        WHERE t.plt_cn = $(_qs(key))
        ORDER BY t.subp, t.tree, t.age_calc, t.measdate""")
        (d, Int(nm))
    end

    plot = (; key, lat = round(s.lat[i]; digits = 5), lon = round(s.lon[i]; digits = 5),
              n_trees = Int(s.n_trees[i]), n_visits = Int(s.n_visits[i]), first_year = Int(s.y0[i]), last_year = Int(s.y1[i]),
              epa_l3 = s.epa_l3[i], epa_l4 = s.epa_l4[i], ecosubcd = s.ecosubcd[i],
              land_use = [s.lu_opts[j] for j in eachindex(s.lu_opts) if (s.lu[i] >> (j - 1)) & 1 == 1],
              cut_events = [s.ce_opts[j] for j in eachindex(s.ce_opts) if (s.ce[i] >> (j - 1)) & 1 == 1],
              ownership = [get(OWN_LABELS, s.og_opts[j], string(s.og_opts[j])) for j in eachindex(s.og_opts) if (s.og[i] >> (j - 1)) & 1 == 1],
              cycles = [s.cy_opts[j] for j in eachindex(s.cy_opts) if (s.cy[i] >> (j - 1)) & 1 == 1],
              species = s.species[i][1:min(12, end)])
    _, xlab, xtitle, xunit = AXIS_COLS[xcol]; _, ylab, ytitle, yunit = AXIS_COLS[ycol]
    meta = (; plot, n_match, x = xcol, y = ycol, color_by = cby, repr,
              x_title = xtitle * (logx ? "  (log)" : ""),
              y_title = (deriv && repr == "all" ? "Δ($ytitle) / Δ($xunit)" : ytitle) * (logy ? "  (log)" : ""),
              derivative = deriv && repr == "all", logx, logy, imputed)
    if nrow(df) == 0
        return (; meta..., n_plotted = 0, sampled = false, levels = [], trees = [], fit = nothing)
    end

    # colour levels from ALL fetched rows (sorted by string), as in the notebook
    lvls = sort!(unique(filter(!isnothing, [_levelstr(v) for v in df[!, Symbol(cby)]])))
    nl = max(length(lvls), 1)
    colors = Dict(l => (nl <= 10 ? TAB10[mod1(k, 10)] : turbo_hex(nl == 1 ? 0.0 : (k - 1) / (nl - 1))) for (k, l) in enumerate(lvls))

    xs = df[!, Symbol(xcol)]; ys = df[!, Symbol(ycol)]
    ok = [xs[r] !== missing && ys[r] !== missing for r in 1:nrow(df)]
    # group rows per tree (df is ordered by subp, tree)
    groups = Vector{Vector{Int}}()
    r = 1
    while r <= nrow(df)
        r2 = r
        while r2 < nrow(df) && df.subp[r2+1] == df.subp[r] && df.tree[r2+1] == df.tree[r]
            r2 += 1
        end
        push!(groups, collect(r:r2)); r = r2 + 1
    end
    trees = Any[]; fx = Float64[]; fy = Float64[]; used = Set{String}()
    cl(l) = l === nothing ? GREY : colors[l]
    if repr == "all"
        for gr in groups
            rows = [q for q in gr if ok[q]]
            isempty(rows) && continue
            lv = nothing
            for q in rows
                lv = _levelstr(df[q, Symbol(cby)]); lv === nothing || break
            end
            px = Float64[Float64(xs[q]) for q in rows]; py = Float64[Float64(ys[q]) for q in rows]
            id = "$(df.subp[rows[1]])-$(df.tree[rows[1]])"
            if deriv
                length(rows) < 2 && continue
                dx = diff(px); keep = dx .> 0
                any(keep) || continue
                qx = px[2:end][keep]; qy = (diff(py) ./ dx)[keep]
                lv === nothing || push!(used, lv)
                push!(trees, (; id, lv, color = cl(lv), segs = [(; dash = false, x = _r.(qx), y = _r.(qy))], dot = false,
                                imp = [], term = nothing))
                append!(fx, qx); append!(fy, qy); continue
            end
            lv === nothing || push!(used, lv)
            intro = [String(df.intro_type[q]) for q in rows]
            if length(rows) < 2
                push!(trees, (; id, lv, color = cl(lv), segs = [], dot = true, dotxy = (_r(px[1]), _r(py[1])), imp = [],
                                term = _term(df, rows, px, py, 1)))
                append!(fx, px); append!(fy, py); continue
            end
            k = something(findfirst(!=("phantom_backfill"), intro), length(intro) + 1) - 1   # count of leading phantom rows
            segs = if 0 < k < length(px)
                [(; dash = true, x = _r.(px[1:k+1]), y = _r.(py[1:k+1])), (; dash = false, x = _r.(px[k+1:end]), y = _r.(py[k+1:end]))]
            else
                [(; dash = false, x = _r.(px), y = _r.(py))]
            end
            imp = imputed ? [(_r(px[j]), _r(py[j])) for j in eachindex(intro) if intro[j] in ("phantom_backfill", "gap_interpolated")] : []
            ti = findfirst(q -> df.is_terminal[q] === true, rows)
            push!(trees, (; id, lv, color = cl(lv), segs, dot = false, imp, term = ti === nothing ? nothing : _term(df, rows, px, py, ti)))
            append!(fx, px); append!(fy, py)
        end
    else
        # ONE representative REAL measurement per tree (first/last by age_calc = measdate order)
        for gr in groups
            rows = [q for q in gr if ok[q] && df.intro_type[q] == "real"]
            isempty(rows) && continue
            q = repr == "first" ? rows[1] : rows[end]
            lv = _levelstr(df[q, Symbol(cby)]); lv === nothing || push!(used, lv)
            x = Float64(xs[q]); y = Float64(ys[q])
            push!(trees, (; id = "$(df.subp[q])-$(df.tree[q])", lv, color = cl(lv), segs = [], dot = true,
                            dotxy = (_r(x), _r(y)), imp = [], term = nothing))
            push!(fx, x); push!(fy, y)
        end
    end
    fit = dofit ? _fit(fx, fy, logx, logy) : nothing
    levels = [(; level = l, color = colors[l]) for l in lvls if l in used]
    (; meta..., n_plotted = length(trees), sampled = n_match > max_traj, levels, trees, fit,
       has_phantom = any(x -> x == "phantom_backfill", df.intro_type))
end

function _term(df, rows, px, py, j)
    fate = String(df.terminal_fate[rows[j]])
    shape, filled = get(FATE_MARK, fate, ("circle", false))
    (; x = _r(px[j]), y = _r(py[j]), fate, shape, filled)
end

# --- routes ----------------------------------------------------------------------------------------------------------------
_err(e) = HTTP.Response(400, ["Content-Type" => "application/json"], JSON3.write((; error = sprint(showerror, e))))

function _qparams(req)
    q = HTTP.queryparams(HTTP.URI(req.target))
    d = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in q)
    for k in (:west, :south, :east, :north)
        haskey(d, k) && (d[k] = parse(Float64, d[k]))
    end
    d
end

function register!()
    @get "/api/fiaplots/options" function (req)
        try
            lock(_LOCK) do
                options()
            end
        catch e
            @warn "fiaplots options error" exception = (e, catch_backtrace())
            _err(e)
        end
    end
    @get "/api/fiaplots/points" function (req)
        try
            points(_qparams(req))
        catch e
            @warn "fiaplots points error" exception = (e, catch_backtrace())
            _err(e)
        end
    end
    @post "/api/fiaplots/points" function (req)
        try
            points(JSON3.read(String(req.body)))
        catch e
            @warn "fiaplots points error" exception = (e, catch_backtrace())
            _err(e)
        end
    end
    @post "/api/fiaplots/trajectories" function (req)
        try
            trajectories(JSON3.read(String(req.body)))
        catch e
            @warn "fiaplots trajectories error" exception = (e, catch_backtrace())
            _err(e)
        end
    end
    nothing
end

end # module
