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
        DBInterface.execute(db, _FVS_CASES_CREATE)
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
        DBInterface.execute(db, _FVS_SUMMARY_CREATE)     # FVS_Cases is registered by write_dbs_cases!
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
        DBInterface.execute(db, _FVS_CARBON_CREATE)
        ins = "INSERT INTO FVS_Carbon VALUES (" * join(fill("?", 14), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for row in rows
            yr = row[1]; r = row[2]
            DBInterface.execute(stmt, (caseid, standid, Int(yr),
                Float64(r.aboveground), Float64(r.merch), Float64(r.belowground),
                Float64(r.belowground_dead), Float64(r.standing_dead), Float64(r.down_wood),
                Float64(r.forest_floor), Float64(r.shrub_herb), Float64(r.total), 0.0, 0.0))
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
        DBInterface.execute(db, _FVS_FUELS_CREATE)
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
        DBInterface.execute(db, _FVS_SNAGSUM_CREATE)
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
        DBInterface.execute(db, _FVS_DWDVOL_CREATE)
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
        DBInterface.execute(db, _FVS_DWDCOV_CREATE)
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
  Flame_length real null, Scorch_height real null,
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
        DBInterface.execute(db, _FVS_BURNREPORT_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_BurnReport VALUES (" * join(fill("?", 22), ",") * ")")
        for b in burns
            m = b.mois                                   # 2×5: dead 1/10/100/1000hr+duff, live woody/herb
            fm = b.models                                # vector of (model, weight); pad to 4
            mw(i) = i <= length(fm) ? Int(fm[i][1]) : 0
            ww(i) = i <= length(fm) ? Float64(fm[i][2]) : 0.0
            DBInterface.execute(stmt, (caseid, standid, Int(b.year),
                Float64(m[1,1])*100, Float64(m[1,2])*100, Float64(m[1,3])*100, Float64(m[1,4])*100,
                Float64(m[1,5])*100, Float64(m[2,1])*100, Float64(m[2,2])*100,
                Float64(b.wind), 0, Float64(b.flame), Float64(b.scorch),
                mw(1), ww(1), mw(2), ww(2), mw(3), ww(3), mw(4), ww(4)))
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
  TreeIndex int, SpeciesFVS text, SpeciesPLANTS text, SpeciesFIA text, TPA real,
  DBH real, DG real, Ht real, HtG real, PctCr int, CrWidth real, BAPctile real,
  PtBAL real, TCuFt real, MCuFt real, SCuFt real, BdFt real, TruncHt int,
  Ht2TDCF real, Ht2TDBF real, TreeAge real);
"""

"""
    treelist_snapshot(s, year, prdlen) -> (year, prdlen, rows)

Capture the start-of-cycle (pre-thin) tree list for the FVS_TreeList table — one tuple per
live record (the columns FVSjl computes directly). Called per cycle by `write_sum_file`'s
`cycle_hook`; the tuples are written later by `write_dbs_treelist!`.
"""
function treelist_snapshot(s::StandState, year::Integer, prdlen::Integer)
    t = s.trees; c = s.coef; pbal = s.density.point_bal
    g = s.plot.gross_space                      # TPA is per-acre = t.tpa/g (Fortran PROB/GROSPC)
    rows = Vector{Any}[]
    @inbounds for i in 1:t.n
        sp = Int(t.species[i])
        push!(rows, Any[string(Int(t.tree_id[i])), i, strip(c.code_alpha[sp]),
            strip(c.code_plants[sp]), strip(c.code_fia[sp]), Float64(t.tpa[i] / g),
            Float64(t.dbh[i]), Float64(t.diam_growth[i]), Float64(t.height[i]),
            Float64(t.ht_growth[i]), Int(t.crown_pct[i]), Float64(t.crown_width[i]),
            Float64(t.crown_ratio[i]), Float64(i <= length(pbal) ? pbal[i] : 0f0),
            Float64(t.cuft_vol[i]), Float64(t.merch_cuft_vol[i]), Float64(t.saw_cuft_vol[i]),
            Float64(t.bdft_vol[i]), Int(t.trunc[i]), Float64(t.merch_top_cf[i]),
            Float64(t.merch_top_bf[i]), Float64(t.birth_age[i])])
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
    nsp = length(co.code_alpha)
    sditype = lpad(c.zeide_sdi ? "ZEIDE" : "REINEKE", 7)   # Fortran right-justifies (e.g. "  ZEIDE")
    db = SQLite.DB(dbpath)
    try
        DBInterface.execute(db, _FVS_INVREF_CREATE)
        ins = "INSERT INTO FVS_InvReference VALUES (" * join(fill("?", 21), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for sp in 1:nsp
            DBInterface.execute(stmt, (caseid, standid, sp,
                String(strip(co.code_alpha[sp])), String(strip(co.code_plants[sp])),
                String(strip(co.code_fia[sp])), sditype,
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
            crwidth = Float64(t.crown_width[i]), bapctile = Float64(t.crown_ratio[i]),
            ptbal = Float64(i <= length(pbal) ? pbal[i] : 0f0), tcuft = Float64(t.cuft_vol[i]),
            mcuft = Float64(t.merch_cuft_vol[i]), scuft = Float64(t.saw_cuft_vol[i]),
            bdft = Float64(t.bdft_vol[i]), truncht = Int(t.trunc[i]),
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
        DBInterface.execute(db, _FVS_CUTLIST_CREATE)
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
        DBInterface.execute(db, _FVS_TREELIST_CREATE)
        ins = "INSERT INTO FVS_TreeList VALUES (" * join(fill("?", 26), ",") * ")"
        stmt = DBInterface.prepare(db, ins)
        for (year, prdlen, rows) in cycles, r in rows
            # r = [TreeId,TreeIndex,SpFVS,SpPLANTS,SpFIA,TPA,DBH,DG,Ht,HtG,PctCr,CrWidth,
            #      BAPctile,PtBAL,TCuFt,MCuFt,SCuFt,BdFt,TruncHt,Ht2TDCF,Ht2TDBF,TreeAge]
            DBInterface.execute(stmt, (caseid, standid, year, prdlen, r...))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end
