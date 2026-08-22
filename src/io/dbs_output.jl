# =============================================================================
# dbs_output.jl — DBS (database) output: the SQLite FVS_Summary table.
#
# Ported from: dbssumry.f (the FVS_Summary schema + per-cycle insert) + dbscase.f
# (the FVS_Cases registry). Enabled by the DATABASE block (DSNOUT file + SUMMARY).
#
# The FVS_Summary columns are exactly the `.sum` columns FVSjl already computes
# (`SummaryRow`), so this writes the same per-cycle data into a database instead of
# (in addition to) the text `.sum` — the "modern IO / same SQLite outputs" goal. Only
# the Summary table is emitted so far; the other ~18 DBS tables are the C6 chunk.
# =============================================================================

using SQLite
using DBInterface

# Robust table creation: `CREATE TABLE IF NOT EXISTS` keeps a STALE table from an older jl
# version (different column count) → the current INSERT throws "table X has N columns but M
# values supplied". FVS itself tolerates a pre-existing DB; a faithful drop-in must too. This
# drops a table whose live schema's column count disagrees with the CREATE statement's, then
# (re)creates it. Within one run the schema always matches (no drop); only a cross-version
# stale table (whose data is unusable anyway) is dropped. Extracts the table name + expected
# column count (top-level commas + 1) from `create_sql`.
function _ensure_table!(db, create_sql::AbstractString)
    m = match(r"CREATE TABLE IF NOT EXISTS\s+(\w+)\s*\("i, create_sql)
    name = m === nothing ? nothing : m.captures[1]
    if name !== nothing
        existing = 0
        try
            for _ in DBInterface.execute(db, "PRAGMA table_info($name)"); existing += 1; end
        catch; end
        if existing > 0
            body = create_sql[findfirst('(', create_sql)+1 : findlast(')', create_sql)-1]
            depth = 0; expected = 1
            for ch in body
                ch == '(' && (depth += 1); ch == ')' && (depth -= 1)
                (ch == ',' && depth == 0) && (expected += 1)
            end
            existing != expected && DBInterface.execute(db, "DROP TABLE IF EXISTS $name")
        end
    end
    DBInterface.execute(db, create_sql)
    return
end

# FVS_Summary schema (dbssumry.f:50). Column order matches the INSERT below.
const _FVS_SUMMARY_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Summary(
  CaseID text not null, StandID text not null, Year int, Age int, Tpa int, BA int,
  SDI int, CCF int, TopHt int, QMD real, TCuFt int, MCuFt int, SCuFt int, BdFt int,
  RTpa int, RTCuFt int, RMCuFt int, RSCuFt int, RBdFt int, ATBA int, ATSDI int,
  ATCCF int, ATTopHt int, ATQMD real, PrdLen int, Acc int, Mort int, MAI real,
  ForTyp int, SizeCls int, StkCls int);
"""

# FVS_Cases registry (dbscase.f) — full schema: the per-case run metadata that keys every other
# DBS table. Build/run metadata (Version/RV/RunDateTime/CaseID) is environment-specific.
const _FVS_CASES_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Cases(
  CaseID text primary key, Stand_CN text not null, StandID text not null, MgmtID text,
  RunTitle text, KeywordFile text, SamplingWt real, Variant text, Version text, RV text,
  Groups text, RunDateTime text);
"""

# FVSjl build identifiers written into FVS_Cases (the FVS Version / revision strings).
const FVSJL_VERSION = "FVSjl0.1"
const FVSJL_RV = "20260401"

"""
    write_dbs_cases!(dbpath, caseid, standid; ...)

Register a case in the FVS_Cases table (dbscase.f) — the run metadata (stand/mgmt id, sampling
weight, variant, keyword-file, version, timestamp) that keys every other DBS table. Written once
per stand. The build/run metadata (Version/RV/RunDateTime/CaseID) is FVSjl/environment-specific
and not a Fortran-parity field; the simulation fields (StandID/MgmtID/SamplingWt/Variant) match.
"""
function write_dbs_cases!(dbpath::AbstractString, caseid::AbstractString, standid::AbstractString;
                          mgmt_id::AbstractString = "NONE", variant::AbstractString = "SN",
                          title::AbstractString = "", keyword_file::AbstractString = "",
                          sampling_wt::Real = 1.0, stand_cn::AbstractString = "",
                          groups::AbstractString = "", run_datetime::AbstractString = "")
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_CASES_CREATE)
        DBInterface.execute(db, "INSERT OR REPLACE INTO FVS_Cases VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (caseid, stand_cn, standid, mgmt_id, title, keyword_file, Float64(sampling_wt),
             variant, FVSJL_VERSION, FVSJL_RV, groups, run_datetime))
    finally
        SQLite.close(db)
    end
    return dbpath
end

"The DBS column names, in FVS_Summary order (for the `.sum`-row mapping below)."
const FVS_SUMMARY_COLS = (:Year, :Age, :Tpa, :BA, :SDI, :CCF, :TopHt, :QMD, :TCuFt, :MCuFt,
    :SCuFt, :BdFt, :RTpa, :RTCuFt, :RMCuFt, :RSCuFt, :RBdFt, :ATBA, :ATSDI, :ATCCF, :ATTopHt,
    :ATQMD, :PrdLen, :Acc, :Mort, :MAI, :ForTyp, :SizeCls, :StkCls)

# one FVS_Summary row tuple from a SummaryRow (dbssumry.f arg order)
_dbs_row(r::SummaryRow) = (r.year, r.age, r.tpa, r.ba, r.sdi, r.ccf, r.topht, Float64(r.qmd),
    r.cuft, r.mcuft, r.scuft, r.bdft, r.rem_tpa, r.rem_cuft, r.rem_mcuft, r.rem_scuft,
    r.rem_bdft, r.at_ba, r.at_sdi, r.at_ccf, r.at_topht, Float64(r.at_qmd), r.period,
    r.accretion, r.mortality, Float64(r.mai), r.fortype, r.sizecls, r.stockcls)

"""
    write_dbs_summary!(dbpath, caseid, standid, rows; mgmt_id="NONE", variant="SN", title="")

Append a stand's per-cycle summary `rows` (`SummaryRow`s) to the FVS_Summary table of the
SQLite database `dbpath` (created if absent), registering the case in FVS_Cases. `caseid`
keys the rows (one per run/stand). Mirrors dbssumry.f / dbscase.f.
"""
function write_dbs_summary!(dbpath::AbstractString, caseid::AbstractString,
                            standid::AbstractString, rows::AbstractVector{SummaryRow};
                            mgmt_id::AbstractString = "NONE", variant::AbstractString = "SN",
                            title::AbstractString = "")
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_SUMMARY_CREATE)     # FVS_Cases is registered by write_dbs_cases!
        ins = "INSERT INTO FVS_Summary VALUES (" * join(fill("?", 31), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for r in rows
            DBInterface.execute(stmt, (caseid, standid, _dbs_row(r)...))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_Carbon schema (dbsfmcrpt.f:106-120) — the FFE Stand Carbon Report pools in metric tons/ha.
const _FVS_CARBON_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Carbon(
  CaseID text not null, StandID text not null, Year Int null,
  Aboveground_Total_Live real null, Aboveground_Merch_Live real null,
  Belowground_Live real null, Belowground_Dead real null, Standing_Dead real null,
  Forest_Down_Dead_Wood real null, Forest_Floor real null, Forest_Shrub_Herb real null,
  Total_Stand_Carbon real null, Total_Removed_Carbon real null, Carbon_Released_From_Fire real null)"""

"""
    write_dbs_carbon!(dbpath, caseid, standid, rows) -> dbpath

Write the FFE Stand Carbon Report to the `FVS_Carbon` DBS table (dbsfmcrpt.f, DBSFMCRPT). `rows` is the
`(year, report)` collection from the main simulation (`stand_carbon_report` named tuples) — the same
metric-tons/ha pools as the `.out` carbon report. Total-Removed / Released-from-Fire are 0 (no harvest
/ fire carbon accounting on the carbon-report path yet).
"""
function write_dbs_carbon!(dbpath::AbstractString, caseid::AbstractString,
                           standid::AbstractString, rows::AbstractVector)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_CARBON_CREATE)
        ins = "INSERT INTO FVS_Carbon VALUES (" * join(fill("?", 14), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for row in rows
            yr = row[1]; r = row[2]
            rel = length(row) >= 6 ? Float64(row[6]) : 0.0   # Carbon_Released_From_Fire (fmburn! consumed C)
            DBInterface.execute(stmt, (caseid, standid, Int(yr),
                Float64(r.aboveground), Float64(r.merch), Float64(r.belowground),
                Float64(r.belowground_dead), Float64(r.standing_dead), Float64(r.down_wood),
                Float64(r.forest_floor), Float64(r.shrub_herb), Float64(r.total), 0.0, rel))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_EconHarvestValue schema (dbsecharv.f) — per-species per-DIB-grade log-graded harvest value.
# The BF_1000_LOG (unit 4) and FT3_100_LOG (unit 5) columns are populated by jl per row; the unused
# unit's volume columns + the TPA/Tons/DBH columns are null (FVS's null-when-unset itoc/-1 convention).
const _FVS_ECONHARVEST_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_EconHarvestValue(
  CaseID text, Year int, SpeciesFVS text, SpeciesPLANTS text, SpeciesFIA text,
  Min_DIB real, Max_DIB real, Min_DBH real, Max_DBH real,
  TPA_Removed int, TPA_Value int, Tons_Per_Acre int, Ft3_Removed int, Ft3_Value int,
  Board_Ft_Removed int, Board_Ft_Value int, Total_Value int)"""

"""
    write_dbs_econharvest!(dbpath, caseid, rows) -> dbpath

Write the log-graded harvest-value detail (`FVS_EconHarvestValue`, dbsecharv.f) — one row per
(year, species, unit, DIB grade) from `econ_harvest_value_rows`. Board-foot (unit 4) and cubic-foot
(unit 5) removed/value/total are filled per row; the TPA/Tons/DBH columns are null.
"""
function write_dbs_econharvest!(dbpath::AbstractString, caseid::AbstractString, rows::AbstractVector)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_ECONHARVEST_CREATE)
        ins = "INSERT INTO FVS_EconHarvestValue VALUES (" * join(fill("?", 17), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for r in rows
            DBInterface.execute(stmt, (caseid, r.year, r.sp_fvs, r.sp_plants, r.sp_fia,
                r.min_dib, r.max_dib, missing, missing,
                missing, missing, missing, r.ft3_removed, r.ft3_value,
                r.bf_removed, r.bf_value, r.total_value))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_Fuels schema (dbsfuels.f:64-86) — FFE surface + standing fuel loadings (tons/ac biomass).
const _FVS_FUELS_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Fuels(
  CaseID text not null, StandID text not null, Year Int null,
  Surface_Litter real null, Surface_Duff real null, Surface_lt3 real null, Surface_ge3 real null,
  Surface_3to6 real null, Surface_6to12 real null, Surface_ge12 real null,
  Surface_Herb real null, Surface_Shrub real null, Surface_Total real null,
  Standing_Snag_lt3 real null, Standing_Snag_ge3 real null, Standing_Foliage real null,
  Standing_Live_lt3 real null, Standing_Live_ge3 real null, Standing_Total real null,
  Total_Biomass Int null, Total_Consumed Int null, Biomass_Removed Int null)"""

"""
    write_dbs_fuels!(dbpath, caseid, standid, rows) -> dbpath

Write the FFE fuel loadings to the `FVS_Fuels` DBS table (dbsfuels.f). `rows` is the `(year, …, fuel)`
collection from the main simulation, where `fuel` is an `ffe_fuel_loadings` named tuple (tons/ac biomass).
"""
function write_dbs_fuels!(dbpath::AbstractString, caseid::AbstractString,
                          standid::AbstractString, rows::AbstractVector)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_FUELS_CREATE)
        ins = "INSERT INTO FVS_Fuels VALUES (" * join(fill("?", 22), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for row in rows
            yr = row[1]; f = row[3]
            DBInterface.execute(stmt, (caseid, standid, Int(yr),
                Float64(f.litter), Float64(f.duff), Float64(f.lt3), Float64(f.ge3),
                Float64(f.s3to6), Float64(f.s6to12), Float64(f.ge12),
                Float64(f.herb), Float64(f.shrub), Float64(f.surf_total),
                Float64(f.snag_lt3), Float64(f.snag_ge3), Float64(f.foliage),
                Float64(f.live_lt3), Float64(f.live_ge3), Float64(f.stand_total),
                round(Int, f.total_biomass), round(Int, f.consumed), round(Int, f.removed)))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_SnagSum schema (dbsfmssnag.f:103-121) — standing-snag density (stems/ac) by hard/soft × DBH class.
const _FVS_SNAGSUM_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_SnagSum(
  CaseID text not null, StandID text not null, Year Int null,
  Hard_snags_class1 real null, Hard_snags_class2 real null, Hard_snags_class3 real null,
  Hard_snags_class4 real null, Hard_snags_class5 real null, Hard_snags_class6 real null,
  Hard_snags_total real null,
  Soft_snags_class1 real null, Soft_snags_class2 real null, Soft_snags_class3 real null,
  Soft_snags_class4 real null, Soft_snags_class5 real null, Soft_snags_class6 real null,
  Soft_snags_total real null, Hard_soft_snags_total real null)"""

"""
    write_dbs_snagsum!(dbpath, caseid, standid, rows) -> dbpath

Write the FFE snag-summary densities to the `FVS_SnagSum` DBS table (dbsfmssnag.f). `rows` is the
`(year, …, snags)` collection from the main simulation, where `snags` is a `snag_summary` named tuple.
"""
function write_dbs_snagsum!(dbpath::AbstractString, caseid::AbstractString,
                            standid::AbstractString, rows::AbstractVector)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_SNAGSUM_CREATE)
        ins = "INSERT INTO FVS_SnagSum VALUES (" * join(fill("?", 18), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for row in rows
            yr = row[1]; sg = row[4]
            DBInterface.execute(stmt, (caseid, standid, Int(yr),
                Float64.(sg.hard[1:6])..., Float64(sg.hard[7]),
                Float64.(sg.soft[1:6])..., Float64(sg.soft[7]),
                Float64(sg.hard[7] + sg.soft[7])))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_SnagDet schema (dbsfmdsnag.f) — the DETAILED snag report: one row per (species, death-year, DBH-class)
# cohort. Species in all three code forms (FVS/PLANTS/FIA), like FVS_StrClass.
const _FVS_SNAGDET_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_SnagDet(
  CaseID text not null, StandID text not null, Year int null,
  SpeciesFVS text null, SpeciesPLANTS text null, SpeciesFIA text null,
  DBH_Class int null, Death_DBH real null, Current_Ht_Hard real null, Current_Ht_Soft real null,
  Current_Vol_Hard real null, Current_Vol_Soft real null, Total_Volume real null, Year_Died int null,
  Density_Hard real null, Density_Soft real null, Density_Total real null)"""

"""
    write_dbs_snagdet!(dbpath, caseid, standid, rows, coef) -> dbpath

Write the detailed snag report to the `FVS_SnagDet` DBS table (dbsfmdsnag.f). `rows` is the
`(year, detail)` collection where `detail` is a `snag_detail(s)` vector (per species×death-year×DBH-class
cohort). `coef` supplies the FVS/PLANTS/FIA species-code triplet (as `write_dbs_strclass!` does).
"""
function write_dbs_snagdet!(dbpath::AbstractString, caseid::AbstractString,
                            standid::AbstractString, rows::AbstractVector, coef)
    isempty(rows) && return dbpath
    fia3(x)   = lpad(strip(string(x)), 3, '0')
    fvs(i)    = String(strip(coef.code_alpha[i]))
    plants(i) = String(strip(coef.code_plants[i]))
    fia(i)    = String(fia3(coef.code_fia[i]))
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_SNAGDET_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_SnagDet VALUES (" * join(fill("?", 17), ",") * ")")
        for (yr, detail) in rows
            for r in detail
                DBInterface.execute(stmt, (caseid, standid, Int(yr),
                    fvs(r.sp), plants(r.sp), fia(r.sp), Int(r.jcl),
                    Float64(r.death_dbh), Float64(r.hth), Float64(r.hts),
                    Float64(r.vh), Float64(r.vs), Float64(r.tv), Int(r.yrdied),
                    Float64(r.dh), Float64(r.ds), Float64(r.dt)))
            end
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_StrClass schema (dbsstrclass.f:128-174) — SSTAGE stand-structure classification: up to 3 height strata,
# each with DBHNOM / heights / crown base / cover / the two dominant-crown species (FVS/PLANTS/FIA) / status,
# plus the whole-stand strata count, cover, and structure-class label. Two rows/cycle (Removal_Code 0=before-thin,
# 1=after-thin). All values come from `structure_report` (validated bit-exact vs the sstage.f `.out` report).
const _FVS_STRCLASS_STRATUM_COLS = join(
    ["Stratum_$(k)_DBH real null, Stratum_$(k)_Nom_Ht int null, Stratum_$(k)_Lg_Ht int null, " *
     "Stratum_$(k)_Sm_Ht int null, Stratum_$(k)_Crown_Base int null, Stratum_$(k)_Crown_Cover int null, " *
     "Stratum_$(k)_SpeciesFVS_1 text null, Stratum_$(k)_SpeciesFVS_2 text null, " *
     "Stratum_$(k)_SpeciesPLANTS_1 text null, Stratum_$(k)_SpeciesPLANTS_2 text null, " *
     "Stratum_$(k)_SpeciesFIA_1 text null, Stratum_$(k)_SpeciesFIA_2 text null, Stratum_$(k)_Status_Code int null"
     for k in 1:3], ", ")
const _FVS_STRCLASS_CREATE = "CREATE TABLE IF NOT EXISTS FVS_StrClass(CaseID text not null, StandID text not null, " *
    "Year int null, Removal_Code int null, " * _FVS_STRCLASS_STRATUM_COLS *
    ", Number_of_Strata int null, Total_Cover int null, Structure_Class text null)"

"""
    write_dbs_strclass!(dbpath, caseid, standid, rows, coef) -> dbpath

Write the SSTAGE structure classification to the `FVS_StrClass` DBS table (dbsstrclass.f). `rows` is a
per-report collection of `(year, removal_code, report)` where `report` is a `structure_report` named tuple.
Absent species/strata are stored as "--" / 0 (matching the Fortran emitter).
"""
function write_dbs_strclass!(dbpath::AbstractString, caseid::AbstractString,
                             standid::AbstractString, rows::AbstractVector, coef)
    fia3(x)   = lpad(strip(x), 3, '0')                    # FIAJSP as a 3-char zero-padded code (oracle "017")
    fvs(i)    = i > 0 ? String(strip(coef.code_alpha[i]))  : "--"
    plants(i) = i > 0 ? String(strip(coef.code_plants[i])) : "--"
    fia(i)    = i > 0 ? String(fia3(coef.code_fia[i]))     : "--"
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_STRCLASS_CREATE)
        ins = "INSERT INTO FVS_StrClass VALUES (" * join(fill("?", 46), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for (yr, cd, rep) in rows
            vals = Any[caseid, standid, Int(yr), Int(cd)]
            for k in 1:3
                if k <= length(rep.strata)
                    st = rep.strata[k]
                    push!(vals, Float64(st.dbh), round(Int, st.nomht), round(Int, st.lght), round(Int, st.smht),
                          round(Int, st.crnbase), round(Int, st.cover),
                          fvs(st.sp1), fvs(st.sp2), plants(st.sp1), plants(st.sp2), fia(st.sp1), fia(st.sp2),
                          Int(st.status))
                else
                    push!(vals, 0.0, 0, 0, 0, 0, 0, "--", "--", "--", "--", "--", "--", 0)
                end
            end
            push!(vals, Int(rep.nstr), round(Int, rep.cover), _SS_CLASS_LABEL[rep.class + 1])
            DBInterface.execute(stmt, Tuple(vals))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_CalibStats schema (dbscalib.f:88-100) — the DGSCOR large-tree DG-calibration sample statistics per calibrated
# species (one row per species with NUMCAL>0). Written once per stand after calibration.
const _FVS_CALIBSTATS_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_CalibStats(
  CaseID text not null, StandID text not null, TreeSize text not null,
  SpeciesFVSnum int not null, SpeciesFVS text not null, SpeciesPLANTS text not null, SpeciesFIA text not null,
  NumTrees int null, ScaleFactor real null, StdErrRatio real null, WeightToInput real null, ReadCorMult real null)"""

"""
    write_dbs_calibstats!(dbpath, caseid, standid, calib, coef) -> dbpath

Write the DGSCOR large-tree DG-calibration statistics to the `FVS_CalibStats` DBS table (dbscalib.f). One row per
large-tree-calibrated species (`calib.cal_ntree[sp] > 0`): ScaleFactor = exp(COR), StdErrRatio = STDRAT,
WeightToInput = WC, ReadCorMult = exp(COR/WC).
"""
function write_dbs_calibstats!(dbpath::AbstractString, caseid::AbstractString,
                               standid::AbstractString, calib, coef)
    fia3(x) = lpad(strip(x), 3, '0')
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_CALIBSTATS_CREATE)
        ins = "INSERT INTO FVS_CalibStats VALUES (" * join(fill("?", 12), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for sp in 1:length(calib.cal_ntree)
            n = calib.cal_ntree[sp]
            n > 0 || continue
            wc = calib.cal_wci[sp]
            scale = calib.cal_cortem[sp]                       # CORTEM = exp(COR) at calibration time
            readcor = wc > 0f0 ? exp(log(scale) / wc) : scale  # exp(LOG(CORTEM)/WC)
            DBInterface.execute(stmt, (caseid, standid, "LG", Int(sp),
                String(strip(coef.code_alpha[sp])), String(strip(coef.code_plants[sp])), String(fia3(coef.code_fia[sp])),
                Int(n), Float64(scale), Float64(calib.cal_stdrat[sp]), Float64(wc), Float64(readcor)))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_Climate schema (dbsclsum.f:33-48) — Climate-FVS per-species Viability-and-Effects report.
const _FVS_CLIMATE_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Climate(
  CaseID text not null, StandID text not null, Year Int null,
  SpeciesFVS text null, SpeciesPLANTS text null, SpeciesFIA text null,
  Viability real null, BA real null, TPA real null, ViabMort real null, dClimMort real null,
  GrowthMult real null, SiteMult real null, MxDenMult real null, AutoEstbTPA real null)"""

"""
    write_dbs_climate!(dbpath, caseid, standid, rows, coef) -> dbpath

Write the Climate-FVS per-species Viability-and-Effects report to the `FVS_Climate` DBS table (dbsclsum.f,
CLIMREDB toggle). `rows` is the `(year, report_vector)` collection, where each `report_vector` is a
`climate_report` result (one NamedTuple per reported species); `coef` supplies the FVS/PLANTS/FIA species codes.
"""
function write_dbs_climate!(dbpath::AbstractString, caseid::AbstractString,
                            standid::AbstractString, rows::AbstractVector, coef)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_CLIMATE_CREATE)
        ins = "INSERT INTO FVS_Climate VALUES (" * join(fill("?", 15), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for (yr, rep) in rows
            for r in rep
                sp = r.sp
                DBInterface.execute(stmt, (caseid, standid, Int(yr),
                    strip(coef.code_alpha[sp]), strip(coef.code_plants[sp]), strip(coef.code_fia[sp]),
                    Float64(r.viab), Float64(r.ba), Float64(r.tpa), Float64(r.mort1), Float64(r.mort2),
                    Float64(r.gmult), Float64(r.sitgm), Float64(r.mxden), Float64(r.potestab)))
            end
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_CanProfile schema (dbsfmcanpr.f:97-104) — FFE canopy crown-fuel profile by 1-ft height layer.
const _FVS_CANPROFILE_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_CanProfile(
  CaseID text not null, StandID text not null, Year Int null,
  Height_m real null, Canopy_Fuel_kg_m3 real null, Height_ft real null, Canopy_Fuel_lbs_acre_ft real null)"""

"""
    write_dbs_canprofile!(dbpath, caseid, standid, rows) -> dbpath

Write the FFE canopy crown-fuel profile to the `FVS_CanProfile` DBS table (dbsfmcanpr.f, CANFPROF keyword). `rows`
is the `(year, crfill)` collection where `crfill` is the length-400 `canopy_crfill` vector (lbs/ac-ft by 1-ft
layer). One row per layer `I` with `crfill[I] > 0`: Height_ft=I, Height_m=I·0.3048, fuel in lbs/ac-ft, and the
kg/m³ conversion (·0.45359237/(4046.856422·0.3048)).
"""
function write_dbs_canprofile!(dbpath::AbstractString, caseid::AbstractString,
                               standid::AbstractString, rows::AbstractVector)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_CANPROFILE_CREATE)
        ins = "INSERT INTO FVS_CanProfile VALUES (" * join(fill("?", 7), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        conv = 0.45359237 / (4046.856422 * 0.3048)
        for (yr, crfill) in rows
            for i in eachindex(crfill)
                crfill[i] > 0f0 || continue
                DBInterface.execute(stmt, (caseid, standid, Int(yr),
                    i * 0.3048, Float64(crfill[i]) * conv, Float64(i), Float64(crfill[i])))
            end
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_Down_Wood_Vol schema (dbsfmdwvol.f:61-79) — down-wood volume (cuft/ac) by DBH bin × hard/soft.
const _FVS_DWDVOL_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Down_Wood_Vol(
  CaseID text not null, StandID text not null, Year Int null,
  DWD_Volume_0to3_Hard real null, DWD_Volume_3to6_Hard real null, DWD_Volume_6to12_Hard real null,
  DWD_Volume_12to20_Hard real null, DWD_Volume_20to35_Hard real null, DWD_Volume_35to50_Hard real null,
  DWD_Volume_ge_50_Hard real null, DWD_Volume_Total_Hard real null,
  DWD_Volume_0to3_Soft real null, DWD_Volume_3to6_Soft real null, DWD_Volume_6to12_Soft real null,
  DWD_Volume_12to20_Soft real null, DWD_Volume_20to35_Soft real null, DWD_Volume_35to50_Soft real null,
  DWD_Volume_ge_50_Soft real null, DWD_Volume_Total_Soft real null)"""

# FVS_Down_Wood_Cov schema (dbsfmdwcov.f:59-75) — down-wood percent cover by DBH bin × hard/soft.
const _FVS_DWDCOV_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Down_Wood_Cov(
  CaseID text not null, StandID text not null, Year Int null,
  DWD_Cover_3to6_Hard real null, DWD_Cover_6to12_Hard real null, DWD_Cover_12to20_Hard real null,
  DWD_Cover_20to35_Hard real null, DWD_Cover_35to50_Hard real null, DWD_Cover_ge_50_Hard real null,
  DWD_Cover_Total_Hard real null,
  DWD_Cover_3to6_Soft real null, DWD_Cover_6to12_Soft real null, DWD_Cover_12to20_Soft real null,
  DWD_Cover_20to35_Soft real null, DWD_Cover_35to50_Soft real null, DWD_Cover_ge_50_Soft real null,
  DWD_Cover_Total_Soft real null)"""

"Write the FFE down-wood VOLUME (cuft/ac) to FVS_Down_Wood_Vol (dbsfmdwvol.f). `rows[i][5]` = `ffe_down_wood`."
function write_dbs_dwd_vol!(dbpath, caseid::AbstractString, standid::AbstractString, rows::AbstractVector)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_DWDVOL_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_Down_Wood_Vol VALUES (" * join(fill("?", 19), ",") * ")")
        for row in rows
            dw = row[5]
            DBInterface.execute(stmt, (caseid, standid, Int(row[1]),
                Float64.(dw.vol_hard)..., Float64.(dw.vol_soft)...))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

"Write the FFE down-wood percent COVER to FVS_Down_Wood_Cov (dbsfmdwcov.f). `rows[i][5]` = `ffe_down_wood`."
function write_dbs_dwd_cov!(dbpath, caseid::AbstractString, standid::AbstractString, rows::AbstractVector)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_DWDCOV_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_Down_Wood_Cov VALUES (" * join(fill("?", 17), ",") * ")")
        for row in rows
            dw = row[5]
            DBInterface.execute(stmt, (caseid, standid, Int(row[1]),
                Float64.(dw.cov_hard)..., Float64.(dw.cov_soft)...))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_BurnReport schema (dbsfmburn.f:105-127) — the actual SIMFIRE event's behavior + conditions.
const _FVS_BURNREPORT_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_BurnReport(
  CaseID text not null, StandID text not null, Year int null,
  One_Hr_Moisture real null, Ten_Hr_Moisture real null, Hundred_Hr_Moisture real null,
  Thousand_Hr_Moisture real null, Duff_Moisture real null, Live_Woody_Moisture real null,
  Live_Herb_Moisture real null, Midflame_Wind real null, Slope int null,
  Flame_length real null, Scorch_height real null, Fire_Type text null,
  FuelModl1 int null, Weight1 real null, FuelModl2 int null, Weight2 real null,
  FuelModl3 int null, Weight3 real null, FuelModl4 int null, Weight4 real null)"""

"""
    write_dbs_burnreport!(dbpath, caseid, standid, burns) -> dbpath

Write the SIMFIRE burn events to the `FVS_BurnReport` DBS table (dbsfmburn.f). `burns` is the
`fire.burn_reports` collection (one record per fire) captured by `fmburn!`: dead/live fuel moistures
(×100 = %), midflame wind, flame length, scorch height, and up to four weighted standard fuel models.
Slope is 0 (the SN surface-fire path does not apply a slope term).
"""
function write_dbs_burnreport!(dbpath, caseid::AbstractString, standid::AbstractString, burns::AbstractVector)
    isempty(burns) && return dbpath
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_BURNREPORT_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_BurnReport VALUES (" * join(fill("?", 23), ",") * ")")
        for b in burns
            m = b.mois                                   # 2×5: dead 1/10/100/1000hr+duff, live woody/herb
            fm = b.models                                # vector of (model, weight); pad to 4
            mw(i) = i <= length(fm) ? Int(fm[i][1]) : 0
            ww(i) = i <= length(fm) ? Float64(fm[i][2])*100 : 0.0   # fraction → % (live BurnReport weights are %)
            slp = hasproperty(b, :slope) ? Float64(b.slope)*100 : 0.0  # stand slope 0..1 → % (dbsfmburn.f Slope col)
            ftype = hasproperty(b, :fire_type) ? String(b.fire_type) : "SURFACE"  # Fire_Type (fmcfir.f CFTMP)
            DBInterface.execute(stmt, (caseid, standid, Int(b.year),
                Float64(m[1,1])*100, Float64(m[1,2])*100, Float64(m[1,3])*100, Float64(m[1,4])*100,
                Float64(m[1,5])*100, Float64(m[2,1])*100, Float64(m[2,2])*100,
                Float64(b.wind), slp, Float64(b.flame), Float64(b.scorch), ftype,
                mw(1), ww(1), mw(2), ww(2), mw(3), ww(3), mw(4), ww(4)))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_Mortality schema (dbsfmmort.f:101-122) — fire-killed vs total TPA by DBH class + BA/vol killed, one row
# PER SPECIES (SpeciesFVS/PLANTS/FIA) plus an 'ALL' aggregate row, matching the live FVSne table.
const _FVS_MORTALITY_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Mortality(
  CaseID text not null, StandID text not null, Year Int null,
  SpeciesFVS text null, SpeciesPLANTS text null, SpeciesFIA text null,
  Killed_class1 real null, Total_class1 real null, Killed_class2 real null, Total_class2 real null,
  Killed_class3 real null, Total_class3 real null, Killed_class4 real null, Total_class4 real null,
  Killed_class5 real null, Total_class5 real null, Killed_class6 real null, Total_class6 real null,
  Killed_class7 real null, Total_class7 real null, Bakill real null, Volkill real null)"""

"Write fire mortality (killed vs total TPA by DBH class + BA/vol killed) to FVS_Mortality (dbsfmmort.f): one row
per species (from `b.species_mort`) plus the stand 'ALL' aggregate row."
function write_dbs_mortality!(dbpath, caseid::AbstractString, standid::AbstractString, burns::AbstractVector)
    isempty(burns) && return dbpath
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_MORTALITY_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_Mortality VALUES (" * join(fill("?", 22), ",") * ")")
        clsvals(kil, tot) = (v = Float64[]; for c in 1:7; push!(v, Float64(kil[c]), Float64(tot[c])); end; v)
        for b in burns
            for sm in (hasproperty(b, :species_mort) ? b.species_mort : ())
                DBInterface.execute(stmt, (caseid, standid, Int(b.year), sm.fvs, sm.plants, sm.fia,
                    clsvals(sm.clskil, sm.totcls)..., Float64(sm.bakill), Float64(sm.volkill)))
            end
            DBInterface.execute(stmt, (caseid, standid, Int(b.year), "ALL", "ALL", "ALL",
                clsvals(b.clskil, b.totcls)..., Float64(b.killed_ba), Float64(b.killed_vol)))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_Consumption schema (dbsfuels.f:58, same 22 cols as FVS_Fuels) — fuel CONSUMED by the fire (tons/ac).
const _FVS_CONSUMPTION_CREATE = replace(_FVS_FUELS_CREATE, "FVS_Fuels(" => "FVS_Consumption(")

"Write fuel consumed by the fire (before−after loadings) to FVS_Consumption (dbsfuels.f)."
function write_dbs_consumption!(dbpath, caseid::AbstractString, standid::AbstractString, burns::AbstractVector)
    isempty(burns) && return dbpath
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_CONSUMPTION_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_Consumption VALUES (" * join(fill("?", 22), ",") * ")")
        for b in burns
            f = b.consumed
            DBInterface.execute(stmt, (caseid, standid, Int(b.year),
                Float64(f.litter), Float64(f.duff), Float64(f.lt3), Float64(f.ge3),
                Float64(f.s3to6), Float64(f.s6to12), Float64(f.ge12),
                Float64(f.herb), Float64(f.shrub), Float64(f.surf_total),
                Float64(f.snag_lt3), Float64(f.snag_ge3), Float64(f.foliage),
                Float64(f.live_lt3), Float64(f.live_ge3), Float64(f.stand_total),
                round(Int, f.total_biomass), round(Int, f.consumed), round(Int, f.removed)))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_PotFire schema (dbsfmpf.f:168-195) — potential fire behavior under severe/moderate weather.
const _FVS_POTFIRE_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_PotFire(
  CaseID text not null, StandID text not null, Year int null,
  Surf_Flame_Sev real null, Surf_Flame_Mod real null, Tot_Flame_Sev real null, Tot_Flame_Mod real null,
  PTorch_Sev real null, PTorch_Mod real null, Torch_Index real null, Crown_Index real null,
  Canopy_Ht int null, Canopy_Density real null,
  Mortality_BA_Sev real null, Mortality_BA_Mod real null, Mortality_VOL_Sev real null, Mortality_VOL_Mod real null,
  Pot_Smoke_Sev real null, Pot_Smoke_Mod real null,
  Fuel_Mod1 int null, Fuel_Mod2 int null, Fuel_Mod3 int null, Fuel_Mod4 int null,
  Fuel_Wt1 real null, Fuel_Wt2 real null, Fuel_Wt3 real null, Fuel_Wt4 real null)"""

"""
    write_dbs_potfire!(dbpath, caseid, standid, rows) -> dbpath

Write the Potential Fire report to the `FVS_PotFire` DBS table (dbsfmpf.f). `rows` is the `(year, report)`
collection where `report` is a `potential_fire_report` named tuple (severe/moderate surface fire behavior,
canopy bulk density, torching probabilities, and the severe-case weighted fuel models).
"""
function write_dbs_potfire!(dbpath, caseid::AbstractString, standid::AbstractString, rows::AbstractVector)
    isempty(rows) && return dbpath
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_POTFIRE_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_PotFire VALUES (" * join(fill("?", 27), ",") * ")")
        for (yr, r) in rows
            fm = r.models
            mw(i) = i <= length(fm) ? Int(fm[i][1]) : 0
            ww(i) = i <= length(fm) ? Float64(fm[i][2])*100 : 0.0   # fraction → % (live weights are %)
            DBInterface.execute(stmt, (caseid, standid, Int(yr),
                Float64(r.surf_flame_sev), Float64(r.surf_flame_mod), Float64(r.tot_flame_sev), Float64(r.tot_flame_mod),
                Float64(r.ptorch_sev), Float64(r.ptorch_mod), Float64(r.torch_index), Float64(r.crown_index),
                Int(r.canopy_ht), Float64(r.canopy_density),
                Float64(r.mort_ba_sev), Float64(r.mort_ba_mod), Float64(r.mort_vol_sev), Float64(r.mort_vol_mod),
                # potential smoke is lb/ac in the FFE model; live writes tons/ac (PSMOKE·P2T, fmpofl.f:303)
                Float64(r.smoke_sev) * _FM_P2T, Float64(r.smoke_mod) * _FM_P2T,
                mw(1), mw(2), mw(3), mw(4), ww(1), ww(2), ww(3), ww(4)))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_Hrv_Carbon schema (dbsfmhrpt.f:94-100) — harvested-wood-products carbon fate (metric tons C/ha).
const _FVS_HRVCARBON_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Hrv_Carbon(
  CaseID text not null, StandID text not null, Year int null,
  Products real null, Landfill real null, Energy real null, Emissions real null,
  Merch_Carbon_Stored real null, Merch_Carbon_Removed real null)"""

"""
    write_dbs_hrvcarbon!(dbpath, caseid, standid, rows) -> dbpath

Write the harvested-wood-products carbon report to the `FVS_Hrv_Carbon` DBS table (dbsfmhrpt.f). `rows` is
the `(year, report)` collection where `report` is a `harvested_carbon_report` named tuple (Products /
Landfill / Energy / Emissions / Stored / Removed, metric tons C/ha).
"""
function write_dbs_hrvcarbon!(dbpath, caseid::AbstractString, standid::AbstractString, rows::AbstractVector)
    isempty(rows) && return dbpath
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_HRVCARBON_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_Hrv_Carbon VALUES (" * join(fill("?", 9), ",") * ")")
        for (yr, r) in rows
            DBInterface.execute(stmt, (caseid, standid, Int(yr),
                Float64(r.products), Float64(r.landfill), Float64(r.energy), Float64(r.emissions),
                Float64(r.stored), Float64(r.removed)))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_TreeList schema (dbstrls.f). The columns FVSjl fills directly; the few not yet
# computed (TreeVal/SSCD/PtIndex/MortPA/MistCD/MDefect/BDefect/EstHt/ActPt) are nullable.
const _FVS_TREELIST_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_TreeList(
  CaseID text not null, StandID text not null, Year int, PrdLen int, TreeId text,
  TreeIndex int, SpeciesFVS text, SpeciesPLANTS text, SpeciesFIA text,
  TreeVal int, SSCD int, PtIndex int, TPA real, MortPA real,
  DBH real, DG real, Ht real, HtG real, PctCr int, CrWidth real, MistCD int, BAPctile real,
  PtBAL real, TCuFt real, MCuFt real, SCuFt real, BdFt real, MDefect int, BDefect int, TruncHt int,
  EstHt real, ActPt int, Ht2TDCF real, Ht2TDBF real, TreeAge real);
"""

"""
    treelist_snapshot(s, year, prdlen) -> (year, prdlen, rows)

Capture the start-of-cycle (pre-thin) tree list for the FVS_TreeList table — one tuple per
live record (the columns FVSjl computes directly). Called per cycle by `write_sum_file`'s
`cycle_hook`; the tuples are written later by `write_dbs_treelist!`.
"""
# FVS_TreeList/FVS_CutList CrWidth = CRWDTH(I) (base/cwidth.f → cwcalc.f). The WESTERN variants with a ported,
# per-tree-bit-exact cwcalc kernel compute the forest-grown value here (+ the cwcalc.f [0.5,99.9] final clamp); the
# eastern open-grown crown_width() handles the rest (0.5 default for unknown species). Shared by both the live-tree
# snapshot and the cut-record builder so the two tables stay consistent. (BM/SO/CA/NC excluded — their kernels are not
# per-tree exact; see the recipe. The kernels bake in one forest's Region-6 BF ⇒ bit-exact on the reference forest.)
function _forest_crwdth(s::StandState, sp::Int, d::Float32, h::Float32, crp)::Float32
    p = s.plot
    # WS (WestSierra) is Region-5: cwcalc.f branches to R5CRWD (a function of sp/D/H only — no forest BF,
    # which R5 skips), so ws_r5crwd is per-tree exact for the TreeList (unlike the R6 BF-baked BM/SO kernels).
    s.variant isa WestSierra && return clamp(ws_r5crwd(sp, d, h), 0.5f0, 99.9f0)
    # NC/Klamath (forest 505 = Region-5) uses R5CRWD too — reuse ws_r5crwd via the NC→WS FIA-species map.
    s.variant isa Klamath && return clamp(nc_r5crwd(sp, d, h), 0.5f0, 99.9f0)
    hi = _cr_hopkins(p.latitude, p.longitude, p.elevation)
    # CA/BM: the FVS_TreeList forest-grown CRWDTH applies the R6 forest BF (cwcalc.f IWHO=0), UNLIKE the FFE PERCOV
    # path (fmcba) which is BF-free — so their kernels default to BF-free and the TreeList opts in via forest_bf=true.
    s.variant isa CentralCalifornia &&
        return clamp(ca_cwcalc(sp, d, h, Float32(crp), p.basal_area, p.elevation, hi; forest_bf = true), 0.5f0, 99.9f0)
    s.variant isa BlueMountains &&
        return clamp(bm_cwcalc(sp, d, h, Float32(crp), p.basal_area, p.elevation, hi; forest_bf = true), 0.5f0, 99.9f0)
    wcw = s.variant isa CentralRockies    ? cr_cwcalc :
          s.variant isa OregonCoast       ? oc_cwcalc :
          s.variant isa Olympic           ? op_cwcalc :
          s.variant isa EasternMontana    ? em_cwcalc :
          s.variant isa InlandEmpire      ? ie_cwcalc :
          s.variant isa Kootenai          ? kt_cwcalc :
          s.variant isa CentralIdaho      ? ci_cwcalc :
          s.variant isa Teton             ? tt_cwcalc :
          s.variant isa Utah              ? ut_cwcalc :
          s.variant isa WestCascades      ? wc_cwcalc :
          s.variant isa PacificNorthwest  ? pn_cwcalc :
          s.variant isa EastCascades      ? ec_cwcalc :
          s.variant isa SouthCentralOregon ? so_cwcalc :
          nothing
    wcw === nothing &&
        return crown_width(s.coef, s.species.code2[sp], d, h, 90, 1, p.latitude, p.longitude, p.elevation)
    return clamp(wcw(sp, d, h, Float32(crp), p.basal_area, p.elevation, hi), 0.5f0, 99.9f0)
end

function treelist_snapshot(s::StandState, year::Integer, prdlen::Integer; cycle::Int = -1)
    t = s.trees; c = s.coef; pbal = s.density.point_bal
    g = s.plot.gross_space                      # TPA is per-acre = t.tpa/g (Fortran PROB/GROSPC)
    rows = Vector{Any}[]
    # FVS_TreeList CrWidth = CRWDTH(I). For CR (a western variant) that array is filled by base/cwidth.f →
    # cwcalc.f (IWHO=0, forest-grown Bechtold/Crookston library, actual crown ratio + stand BA/elev/Hopkins),
    # NOT the eastern open-grown crown_width. Precompute the CR stand inputs once; per-tree via cr_cwcalc.
    iscr   = s.variant isa CentralRockies
    # CRWDTH: the shared _forest_crwdth dispatch (western cwcalc kernels + [0.5,99.9] clamp, else eastern crown_width).
    # SpeciesFIA: FVS emits the 3-char zero-padded FIA code (FIAJSP). CR's data has 2-digit western codes
    # unpadded ("15","93") vs live "015"/"093" — pad on output (CR-gated; the DATA stays unpadded so
    # resolve_species still string-matches the unpadded input SPCD). Eastern codes are already 3-char.
    fia3(x) = iscr ? lpad(strip(x), 3, '0') : strip(x)
    @inbounds for i in 1:t.n
        sp = Int(t.species[i])
        # Eastern variants: the OPEN-GROWN crown width (crown_width iwho=1, CR=90). Western: the forest-grown
        # cwcalc.f value (per-variant cwcalc) — matches live's CRWDTH (was the crown_width 0.5 default before).
        cw = _forest_crwdth(s, sp, t.dbh[i], t.height[i], t.crown_pct[i])
        # FVS_TreeList metadata columns (dbstrls.f binds): TreeVal=IMC (mort_code), SSCD=ISPECL (special),
        # PtIndex=ITRE (point), MistCD=IDMR=0 (no dwarf mistletoe in SN), MDefect/BDefect=decoded DEFECT
        # (cubic = (DEF−⌊DEF/1e4⌋·1e4)/100; board = DEF−⌊DEF/100⌋·100), EstHt=normht?(normht+5)/100:HT
        # (dbstrls.f:200-202), ActPt=IPVEC(ITRE) (point id). All sourced from jl state.
        df = Int(t.defect[i]); pid = Int(t.plot_id[i])
        mdef = div(df - div(df, 10000) * 10000, 100); bdef = df - div(df, 100) * 100
        estht = t.norm_ht[i] > 0 ? (Float64(t.norm_ht[i]) + 5) / 100 : Float64(t.height[i])
        actpt = (1 <= pid <= length(s.plot.point_ids)) ? Int(s.plot.point_ids[pid]) : pid
        push!(rows, Any[string(Int(t.tree_id[i])), i, strip(c.code_alpha[sp]),
            strip(c.code_plants[sp]), fia3(c.code_fia[sp]),
            Int(t.mort_code[i]), Int(t.special[i]), pid,           # TreeVal, SSCD, PtIndex
            Float64(t.tpa[i] / g), Float64(t.mort_pa[i] / g),      # TPA, MortPA
            Float64(t.dbh[i]), Float64(t.diam_growth[i]), Float64(t.height[i]),
            Float64(t.ht_growth[i]), Int(t.crown_pct[i]), Float64(cw),
            0,                                                     # MistCD
            Float64(t.crown_ratio[i]), Float64(i <= length(pbal) ? pbal[i] : 0f0),
            Float64(t.cuft_vol[i]), Float64(t.merch_cuft_vol[i]), Float64(t.saw_cuft_vol[i]),
            Float64(t.bdft_vol[i]), mdef, bdef, div(Int(t.trunc[i]) + 5, 100),  # BdFt, MDefect, BDefect, TruncHt
            estht, actpt,                                          # EstHt, ActPt (dbstrls.f: (ITRUNC+5)/100)
            Float64(t.merch_top_cf[i]), Float64(t.merch_top_bf[i]), Float64(t.birth_age[i])])
    end
    # CYCLE-0 DEAD RECORDS (dbstrls.f:308-440): at the inventory year only, FVS appends the input dead
    # trees (HISTORY 6-9) at the bottom of the FVS_TreeList — TPA=0, the mortality expansion in MortPA
    # (P=(PROB/GROSPC)/(FINT/FINTM); FINT/FINTM=1 at cycle 0 ⇒ MortPA = tpa/g), DG=HtG=0, with volume and
    # a point-BAL evaluated against the LIVE stand. jl keeps the dead partition at t.n+1 : t.n+ndead
    # (treeinput.jl). CR-gated (the western validation target); the eastern variants share the same latent
    # gap — enabling them needs their FVS_TreeList treelists re-validated (their compute_volumes cover only
    # live), so this stays scoped to CentralRockies for now.
    if cycle == 0 && s.variant isa CentralRockies && t.ndead > 0
        scale = s.plot.pi / g                            # point_basal_area! per-acre scale (PI/GROSPC)
        @inbounds for i in (t.n + 1):(t.n + t.ndead)
            sp = Int(t.species[i]); dd = t.dbh[i]; pid = Int(t.plot_id[i])
            # PtBAL (dbstrls.f IPTBAL=NINT(PTBALT)): BA of larger records at this dead tree's point. FVS's
            # PTBALT accumulates ALL records (live+dead) in descending DBH, so a dead tree's BAL includes the
            # larger DEAD trees too (only stand BA/SDI excludes dead — that's a separate sum, so the .sum stays
            # bit-exact). Ties: the descending sort is stable on array index, so a record of equal DBH is
            # already accumulated iff its index is lower.
            dbal = 0f0
            for j in 1:(t.n + t.ndead)
                j == i && continue
                (Int(t.plot_id[j]) == pid) || continue
                (t.dbh[j] > dd || (t.dbh[j] == dd && j < i)) || continue
                dbal += t.tpa[j] * BA_PER_TREE * t.dbh[j]^2 * scale
            end
            dbal = Float32(round(Int, dbal))              # NINT(PTBALT(I))
            cw = _forest_crwdth(s, sp, dd, t.height[i], t.crown_pct[i])  # cwcalc.f forest-grown + [0.5,99.9] clamp
            df = Int(t.defect[i])
            mdef = div(df - div(df, 10000) * 10000, 100); bdef = df - div(df, 100) * 100
            estht = t.norm_ht[i] > 0 ? (Float64(t.norm_ht[i]) + 5) / 100 : Float64(t.height[i])
            actpt = (1 <= pid <= length(s.plot.point_ids)) ? Int(s.plot.point_ids[pid]) : pid
            push!(rows, Any[string(Int(t.tree_id[i])), i, strip(c.code_alpha[sp]),
                strip(c.code_plants[sp]), fia3(c.code_fia[sp]),
                Int(t.mort_code[i]), Int(t.special[i]), pid,
                0.0, Float64(t.tpa[i] / g),                # TPA=0, MortPA = mortality expansion
                Float64(dd), 0.0, Float64(t.height[i]),    # DBH, DG=0, Ht
                0.0, Int(t.crown_pct[i]), Float64(cw),     # HtG=0, PctCr, CrWidth
                0,                                         # MistCD
                Float64(t.crown_ratio[i]), Float64(dbal),  # BAPctile, PtBAL
                Float64(t.cuft_vol[i]), Float64(t.merch_cuft_vol[i]), Float64(t.saw_cuft_vol[i]),
                Float64(t.bdft_vol[i]), mdef, bdef, div(Int(t.trunc[i]) + 5, 100),  # TruncHt (ITRUNC+5)/100
                estht, actpt,
                Float64(t.merch_top_cf[i]), Float64(t.merch_top_bf[i]), Float64(t.birth_age[i])])
        end
    end
    return (Int(year), Int(prdlen), rows)
end

# FVS_InvReference schema (dbsinvref.f): a once-per-case dump of the variant's species master
# list — codes, SDI method/max, site index, and the cubic/board volume-equation specs.
const _FVS_INVREF_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_InvReference(
  CaseID text not null, StandID text not null, SpeciesNum int, SpeciesFVS text,
  SpeciesPlants text, SpeciesFIA text, SDIType text, SDIMax int, SiteIndex int,
  CFCruiseType text, CFVolEq text, CFMinDBH real, CFTopDia real, CFStump real,
  CFSawMinDBH real, CFSawTopDia real, CFSawStump real, BFVolEq text, BFMinDBH real,
  BFTopDia real, BFStump real);
"""

"""
    write_dbs_invref!(dbpath, caseid, standid, s)

Write the per-species inventory-reference rows (FVS_InvReference, dbsinvref.f) for stand `s` —
one row per species in the variant master list: FVS/PLANTS/FIA codes, the SDI method + per-species
SDImax and site index, and the cubic/board volume-equation ids and merch specs (min DBH / top
diameter / stump for total, sawtimber, and board). All data the engine already holds after
`compute_volumes!`. A static reference table, so it is written once per stand.
"""
function write_dbs_invref!(dbpath::AbstractString, caseid::AbstractString,
                           standid::AbstractString, s::StandState)
    c = s.control; co = s.coef; p = s.plot; sp_eq = s.species.vol_eq
    nsp = nspecies(s.variant)   # the variant's real species count (code arrays are padded to MAXSP capacity)
    fia3(x) = s.variant isa CentralRockies ? lpad(strip(x), 3, '0') : strip(x)   # 3-char FIAJSP (CR western codes)
    sditype = lpad(c.zeide_sdi ? "ZEIDE" : "REINEKE", 7)   # Fortran right-justifies (e.g. "  ZEIDE")
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_INVREF_CREATE)
        ins = "INSERT INTO FVS_InvReference VALUES (" * join(fill("?", 21), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for sp in 1:nsp
            DBInterface.execute(stmt, (caseid, standid, sp,
                String(strip(co.code_alpha[sp])), String(strip(co.code_plants[sp])),
                String(fia3(co.code_fia[sp])), sditype,
                trunc(Int, p.sp_sdi_def[sp] + 0.5f0), trunc(Int, p.sp_site_index[sp] + 0.5f0),  # FVS NINT (round half up)
                "FVS", String(strip(sp_eq[sp])),
                Float64(c.sp_dbh_min[sp]), Float64(c.sp_top_diam[sp]), Float64(c.sp_stump_ht[sp]),
                Float64(c.sp_scf_dbhmin[sp]), Float64(c.sp_scf_topd[sp]), Float64(c.sp_scf_stump[sp]),
                String(strip(c.sp_bf_vol_eq[sp])),
                Float64(c.sp_bf_dbhmin[sp]), Float64(c.sp_bf_topd[sp]), Float64(c.sp_bf_stump[sp])))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

# FVS_CutList schema (dbscuts.f) — the per-cycle list of REMOVED records (same per-tree columns as
# FVS_TreeList, but TPA = removed trees/acre). The not-yet-computed columns are nullable.
const _FVS_CUTLIST_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_CutList(
  CaseID text not null, StandID text not null, Year int, PrdLen int, TreeId text,
  TreeIndex int, Species text, TreeVal int, SSCD int, PtIndex int, TPA real, MortPA real,
  DBH real, DG real, Ht real, HtG real, PctCr int, CrWidth real, MistCD int, BAPctile real,
  PtBAL real, TCuFt real, MCuFt real, SCuFt real, BdFt real, MDefect int, BDefect int,
  TruncHt int, EstHt real, ActPt int, Ht2TDCF real, Ht2TDBF real, TreeAge real);
"""

# Capture one removed record `i` for FVS_CutList (per-acre removed TPA = prem/GROSPC). The fillable
# per-tree attributes; the rest (TreeVal/SSCD/PtIndex/MortPA/MistCD/MDefect/BDefect/EstHt/ActPt) are
# nullable — exactly as FVS_TreeList. (FVSjl field `crown_ratio` is the BA percentile PCT; `crown_pct`
# is the crown ratio ICR — the confusing names are documented in the TreeList writer.)
function _cut_record(s::StandState, i::Integer, prem::Float32)
    t = s.trees; c = s.coef; g = s.plot.gross_space; pbal = s.density.point_bal
    sp = Int(t.species[i])
    return (treeid = string(Int(t.tree_id[i])), index = Int(i),
            species = String(strip(c.code_alpha[sp])), tpa = Float64(prem / g),
            dbh = Float64(t.dbh[i]), dg = Float64(t.diam_growth[i]), ht = Float64(t.height[i]),
            htg = Float64(t.ht_growth[i]), pctcr = Int(t.crown_pct[i]),
            # CrWidth via the shared forest-grown dispatch (was t.crown_width[i], which is 0 for most variants);
            # TruncHt via dbstrls.f (ITRUNC+5)/100 feet (was the raw hundredths ITRUNC = 100× too large).
            crwidth = Float64(_forest_crwdth(s, sp, t.dbh[i], t.height[i], t.crown_pct[i])),
            bapctile = Float64(t.crown_ratio[i]),
            ptbal = Float64(i <= length(pbal) ? pbal[i] : 0f0), tcuft = Float64(t.cuft_vol[i]),
            mcuft = Float64(t.merch_cuft_vol[i]), scuft = Float64(t.saw_cuft_vol[i]),
            bdft = Float64(t.bdft_vol[i]), truncht = div(Int(t.trunc[i]) + 5, 100),
            ht2tdcf = Float64(t.merch_top_cf[i]), ht2tdbf = Float64(t.merch_top_bf[i]),
            treeage = Float64(t.birth_age[i]))
end

"""
    write_dbs_cutlist!(dbpath, caseid, standid, cycles)

Write the per-cycle removed-record snapshots to FVS_CutList. `cycles` is `[(year, prdlen, recs), …]`
where `recs` are `_cut_record` NamedTuples captured by `_log_cut!` during the cycle's thin.
"""
function write_dbs_cutlist!(dbpath::AbstractString, caseid::AbstractString,
                            standid::AbstractString, cycles)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_CUTLIST_CREATE)
        ins = "INSERT INTO FVS_CutList VALUES (" * join(fill("?", 33), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for (year, prdlen, recs) in cycles, r in recs
            DBInterface.execute(stmt, (caseid, standid, Int(year), Int(prdlen),
                r.treeid, r.index, r.species, missing, missing, missing, r.tpa, missing,
                r.dbh, r.dg, r.ht, r.htg, r.pctcr, r.crwidth, missing, r.bapctile,
                r.ptbal, r.tcuft, r.mcuft, r.scuft, r.bdft, missing, missing,
                r.truncht, missing, missing, r.ht2tdcf, r.ht2tdbf, r.treeage))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

"""
    write_dbs_compute!(dbpath, caseid, standid, var_names, rows)

Write the per-cycle COMPUTE variables to the FVS_Compute table (dbscmpu.f). The schema is
DYNAMIC — one REAL column per COMPUTE variable (`var_names`, in declaration order) — created on
first use. `rows` is `[(year, [(name,value),…]), …]` (the `compute_collect` from `write_sum_file`);
a variable not yet active in a given cycle is written NULL. Only the growing cycles get a row.
"""
function write_dbs_compute!(dbpath::AbstractString, caseid::AbstractString,
                            standid::AbstractString, var_names::Vector{String}, rows)
    isempty(var_names) && return dbpath
    db = SQLite.DB(dbpath)
    try
        cols = join(("\"$v\" real null" for v in var_names), ", ")
        DBInterface.execute(db, "CREATE TABLE IF NOT EXISTS FVS_Compute(" *
            "CaseID text not null, StandID text not null, Year int null, $cols);")
        ins = "INSERT INTO FVS_Compute VALUES (" * join(fill("?", length(var_names) + 3), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for (year, snap) in rows
            vals = Dict{String,Float32}(snap)
            row = Any[caseid, standid, Int(year)]
            for v in var_names
                push!(row, haskey(vals, v) ? Float64(vals[v]) : missing)
            end
            DBInterface.execute(stmt, row)
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end

"""
    write_dbs_treelist!(dbpath, caseid, standid, cycles)

Write the per-cycle tree snapshots (`treelist_snapshot` tuples) to the FVS_TreeList table.
"""
function write_dbs_treelist!(dbpath::AbstractString, caseid::AbstractString,
                             standid::AbstractString, cycles)
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_TREELIST_CREATE)
        ins = "INSERT INTO FVS_TreeList VALUES (" * join(fill("?", 35), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for (year, prdlen, rows) in cycles, r in rows
            # r = [TreeId,TreeIndex,SpFVS,SpPLANTS,SpFIA,TPA,MortPA,DBH,DG,Ht,HtG,PctCr,CrWidth,
            #      BAPctile,PtBAL,TCuFt,MCuFt,SCuFt,BdFt,TruncHt,Ht2TDCF,Ht2TDBF,TreeAge]
            DBInterface.execute(stmt, (caseid, standid, year, prdlen, r...))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end
