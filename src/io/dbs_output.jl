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

# FVS_Cases registry (dbscase.f) — minimal: identify the case + stand + variant.
const _FVS_CASES_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Cases(
  CaseID text not null, StandID text not null, MgmtID text, RunTitle text, Variant text);
"""

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
        DBInterface.execute(db, _FVS_CASES_CREATE)
        DBInterface.execute(db, _FVS_SUMMARY_CREATE)
        DBInterface.execute(db, "INSERT INTO FVS_Cases VALUES (?,?,?,?,?)",
                            (caseid, standid, mgmt_id, title, variant))
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
