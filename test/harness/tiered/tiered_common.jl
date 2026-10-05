# tiered_common.jl — shared definitions for the tiered integration suite (fixture generator + fast-tier runner).
#
# A tiered fixture is a directory test/fixtures/tiered/<v>/ holding:
#   stands.db                       self-contained FIA sub-DB (FVS_STANDINIT_COND + FVS_TREEINIT_COND) for K stands
#   <cn>_<regime>.key               keyfile (DSNin stands.db / DSNout out.db — RELATIVE, run with cwd = a scratch dir)
#   <cn>_<regime>.live.sum          live-oracle .sum data rows (-999 headers stripped)
#   <cn>_<regime>.<Table>.csv       every live-oracle DBS output table except FVS_Cases (CaseID column dropped)
#   <cn>_<regime>.tables            the sorted list of DBS tables the oracle created (presence/absence golden)
#   PROVENANCE.toml                 oracle binary + sha256 + build flags + FVS source commit + date + generator commit
#
# The keyfile regime blocks are ledger_fia.jl's `regime_block` verbatim (same cards the population sweeps use);
# only the DATABASE plumbing and the per-regime report requests are added here.
import SQLite, DBInterface
using FVSjl

include(joinpath(@__DIR__, "..", "fia", "ledger_fia.jl"))   # BIN, VAR, kwrec, regime_block, build_subdb, parse helpers

const TIERED_ROOT = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "tiered"))

# The regime matrix. plant_cal = PLANT at a calendar year (INV_YEAR+10); plant_cyc = the cycle-number form "2.0"
# (both are real user forms; a cycle-number off-by-one hid behind calendar-only testing until 2026-09-19).
# WPBR is NOT in the matrix: the standard FVS{v}_g16 oracles are not linked with the BRUST extension (its live
# validation uses a privately built FVSie_wpbr_i) — test_wpbr carries its own live goldens.
const REGIMES = ["none", "thinbba", "salvage", "plant_cal", "plant_cyc", "simfire",
                 "mistletoe", "rootdis", "climate", "econ", "cover"]

_ledger_regime(r) = r in ("plant_cal", "plant_cyc") ? "plant" : r

# Per-regime report requests: main-context report keywords (fill the arrays), FFE report keywords (inside the
# regime's FMIn block), and DATABASE block-2 toggles (blank value = default table name; recipe: two DATABASE blocks,
# block 1 = DSNout + DSNin + SQL, block 2 = toggles).
function report_requests(r)
    main = String[]; ffe = String[]; tog = ["SUMMARY"]
    if r == "none"
        push!(main, kwrec("TREELIST", "0")); push!(tog, "TREELIDB")
        push!(main, "STRCLASS"); push!(tog, "STRCLSDB")
    elseif r == "thinbba"
        push!(main, kwrec("CUTLIST", "0")); push!(tog, "CUTLIDB")
        push!(main, "STRCLASS"); push!(tog, "STRCLSDB")
    elseif r in ("simfire", "salvage")
        append!(ffe, ["BURNREPT", "MORTREPT", "CARBREPT", "POTFIRE"])
        append!(tog, ["BURNREDB", "MORTREDB", "CARBREDB", "POTFIRDB"])
    elseif r == "econ"
        push!(tog, "ECONRPTS")
    elseif r == "rootdis"
        push!(tog, "RDSUM")
    elseif r == "mistletoe"
        push!(tog, "MISRPTS")
    elseif r == "climate"
        push!(tog, "CLIMREDB")
    elseif r == "cover"
        push!(main, "STRCLASS"); push!(tog, "STRCLSDB")
    end
    (main, ffe, tog)
end

# Insert FFE report keywords into the regime's FMIn…End block (simfire/salvage), else open a fresh FMIn block.
function _with_ffe(block::String, ffe::Vector{String})
    isempty(ffe) && return block
    if startswith(block, "FMIn\n") && endswith(block, "\nEnd")
        return block[1:end-4] * "\n" * join(ffe, "\n") * "\nEnd"
    end
    return (isempty(block) ? "" : block * "\n") * "FMIn\n" * join(ffe, "\n") * "\nEnd"
end

function tiered_keytext(cn::AbstractString, regime::AbstractString, plantyr::Int)
    main, ffe, tog = report_requests(regime)
    py = regime == "plant_cal" ? plantyr : 0
    blk = _with_ffe(regime_block(_ledger_regime(regime), py, cn), ffe)
    lines = String["STDIDENT", cn, "DATABASE", "DSNout", "out.db", "DSNin", "stands.db",
                   "StandSQL", "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
                   "TreeSQL", "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
                   "END", "NUMCYCLE         5.0"]
    append!(lines, main)
    isempty(blk) || push!(lines, blk)
    push!(lines, "DATABASE"); append!(lines, tog); push!(lines, "END")
    append!(lines, ["ECHOSUM", "PROCESS", "STOP"])
    return join(lines, "\n") * "\n"
end

# ---------------------------------------------------------------------------------------------------------------
# .sum data row, FULL fixed-width layout (sumout.f FORMAT 9014: 2I4,I6,I4,I5,2I4,F5.1,9I6,I4,I5,2I4,F5.1,2X,I6,
# I5,I6,2X,F6.1,1X,I3,1X,2I1). Compared as the PRINTED text of each field — print precision, no tolerance.
const SUM_FIELDS = [
    ("Year",1,4), ("Age",5,8), ("TPA",9,14), ("BA",15,18), ("SDI",19,23), ("CCF",24,27), ("TopHt",28,31),
    ("QMD",32,36), ("TCuFt",37,42), ("MCuFt",43,48), ("SCuFt",49,54), ("BdFt",55,60), ("RTPA",61,66),
    ("RTCuFt",67,72), ("RMCuFt",73,78), ("RSCuFt",79,84), ("RBdFt",85,90), ("ATBA",91,94), ("ATSDI",95,99),
    ("ATCCF",100,103), ("ATTopHt",104,107), ("ATQMD",108,112), ("Period",115,120), ("Accr",121,125),
    ("Mort",126,131), ("MAI",134,139), ("ForTyp",141,143), ("SizeCls",145,145), ("StkCls",146,146)]

# METRIC .sum (BC / ON: metric/vbase/sumout.f FORMAT 9014 = 2I4,I6,I4,I5,2I4,F5.1,7I6,I4,I5,2I4,F5.1,2X,I6,I5,I6,2X,F6.1,
# 1X,I3,1X,2I1 — 134 columns): no SCuFt/RSCuFt; IOSUM(4,5,6)=TCu/MCu/"BdFt" (all m³/ha), IOSUM(7..10)=RTPA/RTCu/RMCu/RBdFt.
const SUM_FIELDS_METRIC = [
    ("Year",1,4), ("Age",5,8), ("TPA",9,14), ("BA",15,18), ("SDI",19,23), ("CCF",24,27), ("TopHt",28,31),
    ("QMD",32,36), ("TCuFt",37,42), ("MCuFt",43,48), ("BdFt",49,54), ("RTPA",55,60), ("RTCuFt",61,66),
    ("RMCuFt",67,72), ("RBdFt",73,78), ("ATBA",79,82), ("ATSDI",83,87), ("ATCCF",88,91), ("ATTopHt",92,95),
    ("ATQMD",96,100), ("Period",103,108), ("Accr",109,113), ("Mort",114,119), ("MAI",122,127), ("ForTyp",129,131),
    ("SizeCls",133,133), ("StkCls",134,134)]
metric_sum(v::AbstractString) = uppercase(v) in ("BC", "ON")
sum_fields(v::AbstractString) = metric_sum(v) ? SUM_FIELDS_METRIC : SUM_FIELDS
sum_width(metric::Bool) = metric ? 134 : 146

"True for a FORMAT-9014 summary data row: year in cols 1-4 and the full 146-column layout (134 metric; the FFE report
tables that also start with a year are ~110 columns and are NOT summary rows)."
function is_sum_row(s::AbstractString; metric::Bool = false)
    length(s) >= sum_width(metric) - 6 || return false
    y = tryparse(Int, strip(s[1:4])); y !== nothing && 1000 <= y <= 3000
end

"Data rows of a .sum, as full lines (rstripped, padded to 146 / 134 metric)."
sum_rows(text::AbstractString; metric::Bool = false) =
    [rpad(rstrip(ln), sum_width(metric)) for ln in split(text, '\n') if is_sum_row(rstrip(ln); metric = metric)]

"Lines of a .sum text that are neither -999 headers, summary rows, nor blank (FVS's .sum has none)."
sum_extra_lines(text::AbstractString; metric::Bool = false) =
    [ln for ln in split(text, '\n') if !isempty(strip(ln)) && !startswith(ln, "-999") && !is_sum_row(rstrip(ln); metric = metric)]
sum_field(row::AbstractString, a::Int, b::Int) = strip(row[a:min(b, length(row))])

# ---------------------------------------------------------------------------------------------------------------
# DBS tables → deterministic text rows. CaseID (random per run) is dropped; FVS_Cases is skipped entirely.
const SKIP_TABLES = Set(["FVS_Cases"])
_cell(v) = v === missing || v === nothing ? "" :
           v isa AbstractFloat ? repr(Float64(v)) : string(v)

function db_tables(path::AbstractString)
    isfile(path) || return String[]
    db = SQLite.DB(path)
    names = sort([String(r.name) for r in DBInterface.execute(db,
                  "SELECT name FROM sqlite_master WHERE type='table'") if !(String(r.name) in SKIP_TABLES)])
    SQLite.close(db); names
end

"(header, rows) of a table, CaseID dropped, rows sorted as full text (content-exact, order-insensitive)."
function db_table_rows(path::AbstractString, tbl::AbstractString)
    db = SQLite.DB(path)
    cols = [String(r.name) for r in DBInterface.execute(db, "PRAGMA table_info(\"$tbl\")")]
    hdr = cols
    rows = Vector{Vector{String}}()
    for r in DBInterface.execute(db, "SELECT * FROM \"$tbl\"")
        push!(rows, [_cell(r[i]) for i in 1:length(cols)])
    end
    SQLite.close(db)
    keep = findall(c -> lowercase(c) != "caseid", hdr)
    hdr = hdr[keep]; rows = [r[keep] for r in rows]
    sort!(rows)
    (hdr, rows)
end

_csvesc(s) = (occursin(',', s) || occursin('"', s) || occursin('\n', s)) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
function write_csv(path, hdr, rows)
    open(path, "w") do io
        println(io, join(_csvesc.(hdr), ","))
        for r in rows; println(io, join(_csvesc.(r), ",")); end
    end
end
function read_csv(path)
    lines = readlines(path)
    parse_line(l) = begin
        f = String[]; buf = IOBuffer(); q = false; i = 1; c = collect(l)
        while i <= length(c)
            ch = c[i]
            if q
                if ch == '"' && i < length(c) && c[i+1] == '"'; write(buf, '"'); i += 1
                elseif ch == '"'; q = false
                else write(buf, ch) end
            elseif ch == '"'; q = true
            elseif ch == ','; push!(f, String(take!(buf)))
            else write(buf, ch) end
            i += 1
        end
        push!(f, String(take!(buf))); f
    end
    hdr = parse_line(lines[1]); rows = [parse_line(l) for l in lines[2:end]]
    (hdr, rows)
end

"Year column value of a DBS row (for allowlist keys), or \"*\" when the table has none."
function row_year(hdr, row)
    i = findfirst(==("Year"), hdr); i === nothing ? "*" : row[i]
end
