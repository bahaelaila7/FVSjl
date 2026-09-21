# =============================================================================
# init.jl — stand initialization driver (INITRE entry)
#
# Ported from: base/fvs.f (the INITRE call) + base/initre.f setup.
#
# Brings a fresh `StandState` up to the "ready to simulate" point for ONE stand:
# apply variant BLOCK DATA defaults, then process the keyword stream (which loads
# trees, sets plot design, cycles, site, etc.) until PROCESS/STOP/EOF.
#
# A keyword file may hold several stands separated by PROCESS; `initialize!` does
# one stand and returns how it terminated, so a driver can loop for the rest.
# =============================================================================

"Strip a trailing `.key`/`.KEY` (or any extension) to get the run's base path."
function strip_key_ext(keypath::AbstractString)
    # `splitext` removes ONLY the final path component's extension, so a dot in a
    # parent directory (e.g. `/workspace/.ktwork/run/foo.key`) is left intact. The old
    # `findfirst(".k", …)` matched the first `.k` anywhere in the path — truncating the
    # base at a `.ktwork`-style directory and losing the companion `.tre`.
    return first(splitext(keypath))
end

# Open a keyword source as a `KeywordReader` + the run's base path. A `.yaml`/`.yml`
# input (positional OR structured form) is parsed to records and rendered to `.key`
# text so the engine reads it through the same KEYRDR path; `.key`/anything else is
# read verbatim. The base path (used to locate the companion `.tre`/`.csv`) is the
# input minus its extension.
function _keyword_reader(keypath::AbstractString)
    ext = lowercase(splitext(keypath)[2])
    if ext == ".yaml" || ext == ".yml"
        text = keyfile_string(read_keywords_yaml(keypath))
        return KeywordReader(IOBuffer(text)), first(splitext(keypath))
    end
    return open(io -> KeywordReader(io), keypath), strip_key_ext(keypath)
end

# Resolve the variant to run a stand file as: an explicit caller argument always wins;
# otherwise a YAML's top-level `variant:` field selects it; otherwise default to Southern
# (a `.key` has no variant — it's the stock-FVS binary choice). Keeps `.key` runs and
# explicit-variant callers behaving exactly as before.
function _resolve_variant(keypath::AbstractString, variant::Union{AbstractVariant,Nothing})
    variant === nothing || return variant
    ext = lowercase(splitext(keypath)[2])
    if ext == ".yaml" || ext == ".yml"
        vc = yaml_variant_code(keypath)
        vc === nothing || return variant_from_code(vc)
    end
    return Southern()
end

# Resolve the output format (:sum or :csv) for a run: an explicit `output=` argument wins;
# else a YAML's `output_format:`; else :sum (the legacy fixed-column default). A `.key` carries
# no output preference, so for it the format comes only from the argument (else :sum).
function _resolve_output(keypath::AbstractString, output::Union{Symbol,AbstractString,Nothing})
    norm(x) = (s = lowercase(strip(String(x)));
               s in ("sum", "csv") ? Symbol(s) : error("unknown output format '$x' (use :sum or :csv)"))
    output === nothing || return norm(output)
    ext = lowercase(splitext(keypath)[2])
    if ext == ".yaml" || ext == ".yml"
        of = yaml_output_format(keypath)
        of === nothing || return norm(of)
    end
    return :sum
end

"""
    initialize!(state, kr, base_path) -> Symbol

Initialize one stand from an already-open keyword reader. Applies BLOCK DATA
defaults then processes keywords until PROCESS/STOP/EOF (the returned reason).
"""

# Per-variant merchantability-spec defaults, taken verbatim from each variant's `grinit.f`
# (STMP/TOPD/DBHMIN, BFSTMP/BFTOPD/BFMIND, SCFSTMP/SCFTOPD/SCFMIND). FVS sets these in GRINIT
# before INITRE reads keywords, so a VOLUME/MERCH keyword still overrides them. Without them the
# arrays stayed 0 and FVS_InvReference's merch-spec columns were written as 0 for every species.
# Tuple order: (STMP, TOPD, DBHMIN, BFSTMP, BFTOPD, BFMIND, SCFSTMP, SCFTOPD, SCFMIND).
const _MERCH_DEFAULTS = Dict{String,NTuple{9,Float32}}(
    "AK" => (1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "BC" => (30.0f0, 10.0f0, 17.5f0, 30.0f0, 10.0f0, 17.5f0, 0.0f0, 0.0f0, 0.0f0),   # metric
    "BM" => (1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0),
    "CA" => (1.0f0, 0.0f0, 7.0f0, 1.0f0, 0.0f0, 7.0f0, 1.0f0, 0.0f0, 7.0f0),
    "CI" => (1.0f0, 6.0f0, 8.0f0, 1.0f0, 6.0f0, 8.0f0, 1.0f0, 6.0f0, 8.0f0),
    "CR" => (1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "CS" => (0.5f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "EC" => (1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0),
    "EM" => (1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0),
    "IE" => (1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0),
    "KT" => (1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0, 1.0f0, 4.5f0, 7.0f0),
    "LS" => (0.5f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "NC" => (1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "NE" => (0.5f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "OC" => (1.0f0, 0.0f0, 7.0f0, 1.0f0, 0.0f0, 7.0f0, 1.0f0, 0.0f0, 7.0f0),
    "ON" => (30.0f0, 10.0f0, 0.0f0, 30.0f0, 10.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0),    # metric
    "OP" => (1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "PN" => (1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "SN" => (0.5f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "SO" => (1.0f0, 0.0f0, 9.0f0, 1.0f0, 0.0f0, 9.0f0, 1.0f0, 0.0f0, 9.0f0),
    "TT" => (1.0f0, 6.0f0, 8.0f0, 1.0f0, 6.0f0, 8.0f0, 1.0f0, 6.0f0, 8.0f0),
    "UT" => (1.0f0, 6.0f0, 8.0f0, 1.0f0, 6.0f0, 8.0f0, 1.0f0, 6.0f0, 8.0f0),
    "WC" => (1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0),
    "WS" => (1.0f0, 4.5f0, 7.0f0, 1.0f0, 6.0f0, 10.0f0, 1.0f0, 6.0f0, 10.0f0),
)

# Per-SPECIES exceptions that each grinit.f applies after the all-species loop (e.g. BM `DBHMIN(7)=6.0`,
# `BFMIND(7)=6.0`, `SCFMIND(7)=6.0` — lodgepole). Entries: species index => (DBHMIN, BFMIND, SCFMIND).
const _MERCH_SP_OVERRIDES = Dict{String,Dict{Int,NTuple{3,Float32}}}(
    "BC" => Dict(7  => (12.5f0, 12.5f0, -1f0)),      # BC has no SCF arrays in grinit
    "BM" => Dict(7  => (6.0f0, 6.0f0, 6.0f0)),
    "CA" => Dict(11 => (6.0f0, 6.0f0, 6.0f0)),
    "CI" => Dict(7  => (7.0f0, 7.0f0, 7.0f0)),
    "EC" => Dict(7  => (6.0f0, 6.0f0, 6.0f0)),
    "EM" => Dict(7  => (6.0f0, 6.0f0, 6.0f0)),
    "IE" => Dict(7  => (6.0f0, 6.0f0, 6.0f0)),
    "KT" => Dict(7  => (6.0f0, 6.0f0, 6.0f0)),
    "OC" => Dict(11 => (6.0f0, 6.0f0, 6.0f0)),
    "TT" => Dict(7  => (7.0f0, 7.0f0, 7.0f0)),
    "UT" => Dict(7  => (7.0f0, 7.0f0, 7.0f0)),
)

"""
    _init_merch_specs!(s)

GRINIT's per-species merchantability defaults (`grinit.f`): total-cubic stump/top/min-DBH, the
board-foot trio and the sawtimber-cubic trio, applied to every species before keywords are read.
"""
function _init_merch_specs!(s::StandState)
    d = get(_MERCH_DEFAULTS, String(variant_code(s.variant)), nothing)
    d === nothing && return s
    c = s.control
    n = length(c.sp_stump_ht)
    @inbounds for i in 1:n
        c.sp_stump_ht[i]   = d[1]; c.sp_top_diam[i]  = d[2]; c.sp_dbh_min[i]    = d[3]
        c.sp_bf_stump[i]   = d[4]; c.sp_bf_topd[i]   = d[5]; c.sp_bf_dbhmin[i]  = d[6]
        c.sp_scf_stump[i]  = d[7]; c.sp_scf_topd[i]  = d[8]; c.sp_scf_dbhmin[i] = d[9]
    end
    ov = get(_MERCH_SP_OVERRIDES, String(variant_code(s.variant)), nothing)
    if ov !== nothing
        @inbounds for (sp, (dbhmin, bfmind, scfmind)) in ov
            sp <= n || continue
            c.sp_dbh_min[sp] = dbhmin; c.sp_bf_dbhmin[sp] = bfmind
            scfmind >= 0f0 && (c.sp_scf_dbhmin[sp] = scfmind)
        end
    end
    return s
end

function initialize!(s::StandState, kr::KeywordReader, base_path::AbstractString;
                     inherited_format::AbstractString = "")
    load_species_coefficients!(s, s.variant)      # BLOCK DATA (species, TREFMT, RNG seed)
    _init_merch_specs!(s)                         # grinit.f per-variant merch-spec defaults (before keywords)
    # TREFMT persists across stands in FVS (it lives in COMMON; INITRE never resets it —
    # only the TREEFMT keyword changes it). A 2nd+ stand with no TREEFMT keyword inherits
    # the previous stand's format, so re-applying the BLOCK DATA default here would break
    # it. Restore the inherited format; a TREEFMT keyword in this stand still overrides.
    isempty(inherited_format) || (s.control.tree_format = inherited_format)
    # KWDFIL (base/filopn.f:44): the keyword-file path minus its `.k*` extension — the stem the
    # ADDTREES bridge (esaddt.f) builds its `.es1`/`.es2` filenames from. `base_path` is already
    # `strip_key_ext(keypath)`, matching FVS's KWDFIL exactly.
    s.control.keyword_file = String(base_path)
    ranseed!(s.rng, false, s.rng.ss)              # INITRE: RANSED(false,...) → reset to seed
    s.plot.gross_space = -1f0                      # GRINIT reset (sn/grinit.f:156)
    @inbounds for i in 1:MAXSP                      # GRINIT size-cap defaults (sn/grinit.f:62)
        s.control.sp_size_cap[i, 1] = 999f0
        s.control.sp_size_cap[i, 2] = 1f0
        s.control.sp_size_cap[i, 3] = 0f0
        s.control.sp_size_cap[i, 4] = 999f0
    end
    # NE GRINIT default cycle length (ne/grinit.f:172 FINT=10, vs SN=5). Set before keyword
    # processing so TIMEINT/NUMCYCLE-period still overrides. `control.year` is the per-cycle
    # length build_cycle_schedule!/grow drive on; SN keeps its 5 default (bit-exact).
    (s.variant isa Northeast || s.variant isa CentralStates) && (s.control.year = 10f0)
    # NC GRINIT DG-measurement period default (nc/grinit.f:169-171 FINT=10, FINTH=5, FINTM=5). Only FINT differs
    # from the generic 5 default. Set before keyword/DB processing so a GROWTH keyword or DG_MEASURE column still
    # overrides. Drives the FINT/FINTM=2 dead-record PROB inflation for the backdated calibration/crown-init DENSE
    # (notre.f) — the DG SCALE stays the NC 0.5 hardcode (gated on growth_dg_set, unaffected).
    s.variant isa Klamath && (s.control.growth_fint = 10f0)
    # TT grinit.f:194-196 FINT=10 / FINTM=5 ⇒ the FINT/FINTM=2 dead-record PROB inflation (notre.f:122-124) for the
    # backdated calibration + crown-init DENSE. Affects ONLY dead-bearing stands (ndead=0 unchanged); the DG SCALE
    # stays 1 (meas_fint gated on growth_dg_set, else htg_period=10 — unaffected).
    s.variant isa Teton && (s.control.growth_fint = 10f0)
    # AK GRINIT defaults (ak/grinit.f): BAF=62.5 (vs the generic 40), FINT=10-yr cycle. Set before
    # keyword processing so a DESIGN BAF / NUMCYCLE-period still overrides. akt01's DESIGN omits BAF ⇒
    # the 62.5 default is load-bearing for the plot expansion (TPA/BA/QMD).
    s.variant isa SoutheastAlaska && (s.plot.baf = 62.5f0; s.control.year = 10f0)
    reason = process_keywords!(s, kr, base_path)
    finalize_design!(s)                            # INITRE end: PI:=IPTINV, GROSPC
    site_setup!(s, s.variant)                      # SITSET: fan site index to all species (variant-specific)
    return reason
end

"""
    initialize(keypath; variant=Southern(), faithful=true) -> (state, reason)

Convenience: open `keypath`, build a fresh state, and initialize the FIRST stand.
"""
function initialize(keypath::AbstractString; variant::AbstractVariant = Southern(),
                    faithful::Bool = true)
    s = StandState(variant; faithful = faithful)
    kr, base = _keyword_reader(keypath)            # accepts .key OR .yaml/.yml (structured or positional)
    reason = initialize!(s, kr, base)
    return s, reason
end

"""
    each_stand(keypath; variant=Southern(), faithful=true) -> Vector{StandState}

Initialize EVERY stand in a multi-stand keyword file (stands are separated by
`PROCESS` and the run ends at `STOP`/EOF). FVS re-runs INITRE per stand — each stand
gets a fresh `StandState` (so `ITRN` resets) — but the tree-record format (`TREFMT`)
persists in COMMON across stands, so it is carried forward here; a `TREEFMT` keyword
inside a later stand still overrides it. Returns the per-stand initialized states
(cyc0-ready; run `notre!`/`setup_growth!`/`compute_volumes!` to project each).
"""
function each_stand(keypath::AbstractString;
                    variant::Union{AbstractVariant,Nothing} = nothing,
                    faithful::Bool = true)
    variant = _resolve_variant(keypath, variant)   # explicit arg wins; else YAML `variant:`; else SN
    kr, base = _keyword_reader(keypath)
    stands = StandState[]
    fmt = ""                                       # TREFMT carried across stands
    while true
        s = StandState(variant; faithful = faithful)
        reason = initialize!(s, kr, base; inherited_format = fmt)
        fmt = s.control.tree_format                # may have been set by a TREEFMT keyword
        # A bare STOP/EOF after the last stand's PROCESS is the run terminator, not a
        # stand: no STDIDENT/STDINFO ran, no trees, no establishment. Don't emit a phantom
        # stand. Signal = a non-blank stand_id (every real stand carries a STDIDENT); the
        # terminator has none. (NOT user_forest_code — ON's forkod defaults it to 915 in
        # site_setup! even for the empty terminator, which would resurrect the phantom.)
        real = !isempty(strip(s.plot.stand_id)) || s.trees.n > 0 || s.estab.active
        real && push!(stands, s)
        reason in (:stop, :eof) && break
    end
    return stands
end
