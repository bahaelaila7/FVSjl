# =============================================================================
# fia_translate.jl — raw FIADB → "FVS-ready" record translation
#
# Transforms raw FIADB rows (PLOT ⋈ COND ⋈ TREE ⋈ SEEDLING ⋈ SITETREE) into the
# FVS-ready records that `fia_database.jl`'s reader consumes — i.e. it REGENERATES the
# FVS_STANDINIT_COND / FVS_TREEINIT_COND rows that the USFS FIA2FVS utility produces.
#
# LINKAGE (verified against the FIADB): FVS_STANDINIT.STAND_CN = COND.CN,
# PLOT_CN = PLOT.CN, FVS_TREEINIT.TREE_CN = TREE.CN, joined COND.PLT_CN = PLOT.CN.
#
# ── SCOPE / FIDELITY BOUNDARY ────────────────────────────────────────────────
# The tree records and the stand's *direct-measurement* fields are fully derivable
# from raw FIADB and reproduce the FVS-ready rows BIT-EXACT (validated in
# test/harness/fia/validate_translate.jl). But several GROWTH-CRITICAL stand fields
# are produced by the EXTERNAL FIA2FVS tool (R/SQL + reference tables) whose logic is
# NOT in this codebase or the FVS Fortran source, and so CANNOT be reproduced here
# bit-exact:
#   • SITE_INDEX     — a species-specific NON-LINEAR site-curve conversion of
#                      COND.SICOND (e.g. sp95: 38→44, 32→39; sp316: 70→74). Needs the
#                      FIA→FVS site-index equations (external).
#   • LOCATION       — REGION*100 + FOREST, but raw COND.ADFORCD is null for most FIA
#                      plots ⇒ REGION/FOREST come from an external state/county→forest
#                      assignment.
#   • BAF/FPA/BRK/NUM_PLOTS — the sampling design, from a PLOT.DESIGNCD → design-params
#                      lookup table (external).
#   • DG_TRANS/DG_MEASURE/… — the growth-calibration setup (external FIA2FVS policy).
# Therefore the PRODUCTION path is the READER (`fia_database.jl`) over the DB's existing
# FVS-ready tables, which ARE FVS's own authoritative translation. This translator is
# for raw FIA lacking those tables; its external-reference fields are passed through
# when supplied (`extra=`) and otherwise left at the reader's defaults.
# =============================================================================

# Read a value from a raw-row Dict (uppercase keys), missing/nothing → default.
_raw(d, k, dv) = (haskey(d, k) && d[k] !== missing && d[k] !== nothing) ? d[k] : dv

"""
    translate_fia_tree(traw) -> Dict{String,Any}

Raw FIADB `TREE` row (uppercase-keyed Dict) → an FVS_TREEINIT_COND-equivalent record
(the columns `apply_fia_trees!` reads). Direct mappings, verified bit-exact vs FVS-ready:
SPCD→SPECIES, DIA→DIAMETER, HT→HT, CR→CRRATIO, CCLCD→CRCLASS, STATUSCD→HISTORY,
CULL→CULL, SUBP→PLOT_ID, TREE→TREE_ID. TREE_COUNT is the raw per-plot tally (1 for a
variable-radius/prism tree; the plot design does the TPA expansion in `notre!`).
"""
function translate_fia_tree(traw::Dict{String,Any})
    sp = _raw(traw, "SPCD", missing)
    return Dict{String,Any}(
        "SPECIES"    => sp === missing ? missing : string(Int(sp)),
        "DIAMETER"   => _raw(traw, "DIA", missing),
        "HT"         => _raw(traw, "HT", missing),
        "CRRATIO"    => _raw(traw, "CR", missing),
        "CRCLASS"    => _raw(traw, "CCLCD", missing),
        "HISTORY"    => _raw(traw, "STATUSCD", 1),
        "CULL"       => _raw(traw, "CULL", missing),
        "PLOT_ID"    => _raw(traw, "SUBP", 1),
        "TREE_ID"    => _raw(traw, "TREE", missing),
        "TREE_COUNT" => _raw(traw, "TREE_COUNT", 1.0),   # raw tally; design expands (notre!)
    )
end

"""
    translate_fia_stand(craw, praw; extra=Dict()) -> Dict{String,Any}

Raw FIADB `COND` (`craw`) + `PLOT` (`praw`) rows → an FVS_STANDINIT_COND-equivalent
record. DIRECT (bit-exact) fields only; the external-reference fields (SITE_INDEX,
LOCATION, sampling design, DG calibration — see the boundary note above) are taken from
`extra` when the caller supplies them (e.g. from a FIA2FVS reference), else omitted so
the reader's defaults apply.
"""
function translate_fia_stand(craw::Dict{String,Any}, praw::Dict{String,Any};
                             extra::Dict{String,Any} = Dict{String,Any}())
    sisp = _raw(craw, "SISP", missing)
    d = Dict{String,Any}(
        "STAND_CN"           => _raw(craw, "CN", missing),
        "INV_YEAR"           => _raw(praw, "MEASYEAR", _raw(praw, "INVYR", missing)),
        "AGE"                => _raw(craw, "STDAGE", missing),
        "ASPECT"             => _raw(craw, "ASPECT", missing),
        "SLOPE"              => _raw(craw, "SLOPE", missing),
        "ELEVFT"             => _raw(praw, "ELEV", missing),
        "SITE_SPECIES"       => sisp === missing ? missing : string(Int(sisp)),
        "SITE_INDEX_BASE_AG" => _raw(craw, "SIBASE", missing),
        "FOREST_TYPE_FIA"    => _raw(craw, "FORTYPCD", missing),
        "PHYSIO_REGION"      => _raw(craw, "PHYSCLCD", missing),
        "STDORGCD"           => _raw(craw, "STDORGCD", missing),
        "STATE"              => _raw(craw, "STATECD", _raw(praw, "STATECD", missing)),
        "COUNTY"             => _raw(craw, "COUNTYCD", _raw(praw, "COUNTYCD", missing)),
        "LATITUDE"           => _raw(praw, "LAT", missing),
        "LONGITUDE"          => _raw(praw, "LON", missing),
    )
    # External-reference (FIA2FVS) fields: pass through when supplied. LOCATION is
    # REGION*100+FOREST when both are provided in `extra`.
    for k in ("SITE_INDEX", "LOCATION", "REGION", "FOREST", "BASAL_AREA_FACTOR",
              "INV_PLOT_SIZE", "BRK_DBH", "NUM_PLOTS", "DG_TRANS", "DG_MEASURE",
              "HTG_TRANS", "HTG_MEASURE", "MORT_MEASURE")
        haskey(extra, k) && (d[k] = extra[k])
    end
    (!haskey(d, "LOCATION") && haskey(extra, "REGION") && haskey(extra, "FOREST")) &&
        (d["LOCATION"] = Int(extra["REGION"]) * 100 + Int(extra["FOREST"]))
    return d
end
