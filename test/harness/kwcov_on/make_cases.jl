# make_cases.jl — build the ON (Ontario, metric) keyword-isolation fixtures from the LIVE oracle.
#
# Usage: julia --project=. test/harness/kwcov_on/make_cases.jl [--only=name1,name2] [--stand=<cn>]
#
# Two Ontario FIA stands from test/fixtures/tiered/on/stands.db (DB-read path; ON has no .tre fixture format of its own that
# exercises the DB reader) × one isolated keyword each: base-model keywords (thinnings, multipliers, volume, estab, ...)
# and the FFE keywords the ON build accepts. Every case requests the Summary (+TreeList/StrClass for base cases, every FFE
# report + DBS table for FFE cases — the two-DATABASE-block recipe), so the keyword is compared in the .sum rows AND in
# every DBS table cell. ON is METRIC (BA m2/ha, DBH cm, HT m): keyword magnitudes below are metric.
# Goldens = the stock-flag LIVE oracle /workspace/.onwork/tiered/FVSon_tiered (ORACLE_SOURCE_AUDIT_2026-09-19.md §9), run as
#   cd <dir> && FVSon_tiered --keywordfile=<stem>.key
# Written to test/fixtures/kwcov_on/on/on_<stand>_<case>.{key,live.sum,<Table>.csv,tables} (+ stands.db once).
using SHA
include(joinpath(@__DIR__, "..", "tiered", "tiered_common.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "kwcov_on"))
const SRCDB = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "tiered", "on", "stands.db"))
const ORACLE = "/workspace/.onwork/tiered/FVSon_tiered"
const STANDS = ["21769405010661", "64767919010661"]      # 60 rec/9 sp (inv 1987) and 48 rec/10 sp (inv 1999)
const INVY = Dict("21769405010661" => 1987, "64767919010661" => 1999)

kw(name, vals...) = rpad(name, 10) * join([v === nothing ? " "^10 : lpad(string(v), 10) for v in vals])

const REPORTS = ["FUELOUT", "FUELREPT", "MORTREPT", "BURNREPT", "POTFIRE", "SNAGSUM", "SNAGOUT", "CARBREPT",
                 "DWDVLOUT", "DWDCVOUT"]
const TOGGLES = ["SUMMARY", "FUELSOUT", "FUELREDB", "BURNREDB", "MORTREDB", "SNAGSUDB", "SNAGOUDB", "CARBREDB",
                 "POTFIRDB", "DWDVLDB", "DWDCVDB"]
const BASE_TOGGLES = ["SUMMARY", "TREELIDB", "STRCLSDB"]

# year offsets are applied per stand: y1 = first-cycle boundary + 1 cycle, y2 = +2 cycles (10-yr cycles)
y1(inv) = inv + 10
y2(inv) = inv + 20
sf(inv, yr = y2(inv), w = 10, m = 2, t = 20, mort = 1, pct = 100, seas = 2) = kw("SIMFIRE", yr, w, m, t, mort, pct, seas)

# name => inv -> (main-context lines, FMIN lines (nothing ⇒ base case))
const BASE = Pair{String,Function}[
    "none"      => inv -> ([], nothing),
    "baimult"   => inv -> ([kw("BAIMULT", y1(inv), 0.0, 1.5)], nothing),
    "htgmult"   => inv -> ([kw("HTGMULT", y1(inv), 0.0, 1.5)], nothing),
    "mortmult"  => inv -> ([kw("MORTMULT", y1(inv), 0.0, 1.5)], nothing),
    "crnmult"   => inv -> ([kw("CRNMULT", y1(inv), 0.0, 1.5)], nothing),
    "bamax"     => inv -> ([kw("BAMAX", 25.0)], nothing),
    "sdimax"    => inv -> ([kw("SDIMAX", 0, 600.0)], nothing),
    "sdicalc"   => inv -> ([kw("SDICALC")], nothing),
    "fixdg"     => inv -> ([kw("FIXDG", y1(inv), 0.0, 0.5), kw("FIXHTG", y1(inv), 0.0, 1.5)], nothing),
    "fixmort"   => inv -> ([kw("FIXMORT", y1(inv), 0.0, 0.05)], nothing),
    "dgstdev"   => inv -> ([kw("DGSTDEV", 1.5)], nothing),
    "serlcorr"  => inv -> ([kw("SERLCORR", 0.8, 0.45)], nothing),
    "nocalib"   => inv -> ([kw("NOCALIB")], nothing),
    "notriple"  => inv -> ([kw("NOTRIPLE")], nothing),
    "numtrip"   => inv -> ([kw("NUMTRIP", 4)], nothing),
    "rannseed"  => inv -> ([kw("RANNSEED", 12345.0)], nothing),
    "timeint"   => inv -> ([kw("TIMEINT", 0, 5)], nothing),
    "cycleat"   => inv -> ([kw("CYCLEAT", inv + 13)], nothing),
    "compress"  => inv -> ([kw("COMPRESS", 3.0, 5.0)], nothing),
    "resetage"  => inv -> ([kw("RESETAGE", y1(inv), 50.0)], nothing),
    "managed"   => inv -> ([kw("MANAGED")], nothing),
    "treeszcp"  => inv -> ([kw("TREESZCP", 0.0, 45.0, 1.0)], nothing),
    "tcondmlt"  => inv -> ([kw("TCONDMLT", y1(inv), 0.0, 100.0)], nothing),
    "setsite"   => inv -> ([kw("SETSITE", y1(inv), nothing, 0.0, 0, 20.0)], nothing),
    "htgstop"   => inv -> ([kw("HTGSTOP", y2(inv), 0, 0, 9999, 1.0, 1.0, 0.0)], nothing),
    "topkill"   => inv -> ([kw("TOPKILL", y2(inv), 0, 0, 9999, 1.0, 0.5, 0.0)], nothing),
    "mortmsb"   => inv -> ([kw("MORTMSB", 10.0, -5.0, 0.9, 8.0, 40.0, 3.0)], nothing),
    "strclass"  => inv -> ([kw("STRCLASS", y1(inv))], nothing),
    "volume"    => inv -> ([kw("VOLUME", y1(inv), 0.0, 2.0, 12.0, 20.0)], nothing),
    "bfvolume"  => inv -> ([kw("BFVOLUME", y1(inv), 0.0, 20.0, 12.0, 1.0)], nothing),
    "mcdefect"  => inv -> ([kw("MCDEFECT", y1(inv), 0.0, 5.0, 10.0, 15.0, 20.0)], nothing),
    "bfdefect"  => inv -> ([kw("BFDEFECT", y1(inv), 0.0, 5.0, 10.0, 15.0, 20.0)], nothing),
    "thinbba"   => inv -> ([kw("THINBBA", y1(inv), 20.0)], nothing),
    "thinaba"   => inv -> ([kw("THINABA", y2(inv), 18.0)], nothing),
    "thinbta"   => inv -> ([kw("THINBTA", y2(inv), 20.0)], nothing),
    "thinata"   => inv -> ([kw("THINATA", y1(inv), 20.0)], nothing),
    "thindbh"   => inv -> ([kw("THINDBH", y1(inv), 0.0, 99.0, 1.0, nothing, 10.0)], nothing),
    "thinht"    => inv -> ([kw("THINHT", y2(inv), 0.0, 12.0, 0.5)], nothing),
    "thincc"    => inv -> ([kw("THINCC", y2(inv), 60.0)], nothing),
    "thinrden"  => inv -> ([kw("THINRDEN", y2(inv), 55.0)], nothing),
    "thinsdi"   => inv -> ([kw("THINSDI", y1(inv), 400.0, 1.0)], nothing),
    "thinprsc"  => inv -> ([kw("THINPRSC", y1(inv), 0.999)], nothing),
    "thinpt"    => inv -> ([kw("SETPTHIN", y2(inv), 0.0, 1.0), kw("THINPT", y2(inv), 60.0)], nothing),
    "thinqfa"   => inv -> ([kw("THINQFA", y1(inv), 0.0, 6.0, 0.0, 1.3, 2.0, 20.0), "1"], nothing),
    "thinmult"  => inv -> ([kw("BAIMULT", y1(inv), 0.0, 1.2), kw("THINBTA", y2(inv), 20.0)], nothing),
    "cuteff"    => inv -> ([kw("CUTEFF", 0.5), kw("THINBTA", y2(inv), 20.0)], nothing),
    "minharv"   => inv -> ([kw("MINHARV", 50.0), kw("THINBTA", y2(inv), 20.0)], nothing),
    "yardloss"  => inv -> ([kw("YARDLOSS", y2(inv), 0.3, 0.5, 0.2), kw("THINBTA", y2(inv), 20.0)], nothing),
    "spleave"   => inv -> ([kw("SPLEAVE", y2(inv), 1, 1), kw("THINBBA", y2(inv), 12.0)], nothing),
    "specpref"  => inv -> ([kw("SPECPREF", y2(inv), 1, 9), kw("THINBTA", y2(inv), 20.0)], nothing),
    "estab"     => inv -> ([kw("ESTAB", y1(inv)), kw("PLANT", y1(inv), "PW", 600.0, 90.0), "END"], nothing),
    "natural"   => inv -> ([kw("ESTAB", y1(inv)), kw("NATURAL", y1(inv), "MH", 400.0, 80.0, 2.0), "END"], nothing),
    "thinbba_estab" => inv -> ([kw("THINBBA", y1(inv), 12.0), kw("ESTAB", y1(inv)), kw("PLANT", y1(inv), "SW", 800.0, 95.0), "END"], nothing),
    "eventmon"  => inv -> ([["IF", "(FRAC(CYCLE/2.0) EQ 0.0)", "THEN", kw("THINBBA", 0.0, 18.0), "ENDIF"]...], nothing),
    "compute"   => inv -> ([kw("COMPUTE", 0), "MYBA = BBA", "END"], nothing),
    "fertiliz"  => inv -> ([kw("FERTILIZ", y1(inv), 200.0)], nothing),
    "tfixarea"  => inv -> ([kw("TFIXAREA", 0.04)], nothing),
]

const FFE = Pair{String,Function}[
    "simfire"        => inv -> ([], [sf(inv)]),
    "simfire_s4"     => inv -> ([], [sf(inv, y2(inv), 25, 1, 30, 1, 40, 4)]),
    "simfire_nomort" => inv -> ([], [sf(inv, y2(inv), 10, 3, 15, 0, 100, 1)]),
    "flameadj_mult"  => inv -> ([], [sf(inv), kw("FLAMEADJ", y2(inv), 1.5)]),
    "flameadj_len"   => inv -> ([], [sf(inv), kw("FLAMEADJ", y2(inv), 1, 2.1, 7.5, 5.0)]),
    "moisture"       => inv -> ([], [sf(inv), kw("MOISTURE", y2(inv), 7, 8, 10, 13, 60, 85, 110)]),
    "fuelmove"       => inv -> ([], [kw("FUELMOVE", y1(inv), 6, 11, 3.0), kw("FUELMOVE", y1(inv), 4, 3, 0, 0.5), sf(inv)]),
    "fuelmove_imp"   => inv -> ([], [kw("FUELMOVE", y1(inv), 0, 10, 4.0), kw("FUELMOVE", y1(inv) + 5, 10, 0, 0, 0.3, 0.5), sf(inv)]),
    "pileburn"       => inv -> ([kw("THINBTA", y1(inv), 20.0)], [kw("PILEBURN", y1(inv), 1), sf(inv)]),
    "pileburn_jack"  => inv -> ([kw("THINBTA", y1(inv), 20.0)], [kw("PILEBURN", y1(inv), 2, 100, 30, 60, 0), sf(inv)]),
    "fueltret"       => inv -> ([kw("THINBTA", y1(inv), 20.0)], [kw("FUELTRET", y1(inv), 1, 2), sf(inv, y1(inv) + 3)]),
    "potfmois"       => inv -> ([], [kw("POTFMOIS", 1, 6, 7, 9, 12, 55, 80, 100), kw("POTFMOIS", 2, 12, 13, 15, 18, 80, 120, 140)]),
    "potfwind"       => inv -> ([], [kw("POTFWIND", 40, 20)]),
    "potftemp"       => inv -> ([], [kw("POTFTEMP", 30, 15)]),
    "potfseas"       => inv -> ([], [kw("POTFSEAS", 2, 3)]),
    "potfpab"        => inv -> ([], [kw("POTFPAB", 80, 30)]),
    "snagbrk"        => inv -> ([], [kw("SNAGBRK", "ALL", 12, 10, 30, 25), sf(inv)]),
    "snagfall"       => inv -> ([], [kw("SNAGFALL", "ALL", 2.0, 0.8), sf(inv)]),
    "snagdcay"       => inv -> ([], [kw("SNAGDCAY", "ALL", 1.5), sf(inv)]),
    "snagpsft"       => inv -> ([], [kw("SNAGPSFT", "ALL", 0.9), sf(inv)]),
    "snagpbn"        => inv -> ([], [kw("SNAGPBN", 0.8, 0.5, 3, 12, 1), sf(inv)]),
    "snaginit"       => inv -> ([], [kw("SNAGINIT", "ALL", 30, 40, 40, 5, 20), sf(inv)]),
    "salvage"        => inv -> ([], [sf(inv, y1(inv)), kw("SALVAGE", y1(inv) + 5, 0, 999, 5, 1, 0.9, 0)]),
    "drought"        => inv -> ([], [kw("DROUGHT", y1(inv) + 5, 5), sf(inv)]),
    "fuelinit"       => inv -> ([], [kw("FUELINIT", 0.3, 0.8, 2, 3, 3.5, 4, 1, 0.5, 0.2, 1.5, 10), sf(inv)]),
    "fuelsoft"       => inv -> ([], [kw("FUELSOFT", 0.1, 0.2, 0.5, 1, 1.5, 2, 1, 0.5, 0.2), sf(inv)]),
    "fueldcay"       => inv -> ([], [kw("FUELDCAY", 1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), sf(inv)]),
    "fuelmult"       => inv -> ([], [kw("FUELMULT", 1.2, 1.5), sf(inv)]),
    "duffprod"       => inv -> ([], [kw("DUFFPROD", 1, 0.3, 0.4, 0.5, 0.6, 0.7), sf(inv)]),
    "fuelpool"       => inv -> ([], [kw("FUELPOOL", "ALL", 2), sf(inv)]),
    "carbcalc"       => inv -> ([], [kw("CARBCALC", 1, 1, 0.05, 8, 10), sf(inv)]),
    "carbcalc_t"     => inv -> ([], [kw("CARBCALC", 0, 2)]),
    "carbcut"        => inv -> ([kw("THINBTA", y1(inv), 20.0)], [kw("CARBCUT")]),
    "fmodlist"       => inv -> ([], [kw("FMODLIST", y2(inv), 5, 0), sf(inv)]),
    "fuelmodl"       => inv -> ([], [kw("FUELMODL", y2(inv), 8, 0.6, 10, 0.4), "", sf(inv)]),
    "statfuel"       => inv -> ([], [kw("STATFUEL"), sf(inv)]),
    "soilheat"       => inv -> ([], [kw("SOILHEAT", y2(inv), nothing, 2), sf(inv)]),
    "mortclas"       => inv -> ([], [kw("MORTCLAS", 10, 20, 30, 40, 50, 60, 75), sf(inv)]),
    "fmortmlt"       => inv -> ([], [kw("FMORTMLT", y2(inv), 0.5, "ALL", 0, 999), sf(inv)]),
    "cancalc"        => inv -> ([], [kw("CANCALC", 1, 6, 0, 2, 100), sf(inv)]),
    "canfprof"       => inv -> ([], [kw("CANFPROF"), sf(inv)]),
]

function casekey(cn, name, main, fmin)
    L = String["STDIDENT", cn, "DATABASE", "DSNout", "out.db", "DSNin", "stands.db",
               "StandSQL", "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
               "TreeSQL", "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL", "END",
               "NUMCYCLE         5.0"]
    if fmin === nothing
        push!(L, "TREELIST           0"); push!(L, "STRCLASS")
        append!(L, main)
        push!(L, "DATABASE"); append!(L, BASE_TOGGLES); push!(L, "END")
    else
        append!(L, main)
        push!(L, "FMIN"); append!(L, fmin); append!(L, REPORTS); push!(L, "END")
        push!(L, "DATABASE"); append!(L, TOGGLES); push!(L, "END")
    end
    append!(L, ["ECHOSUM", "PROCESS", "STOP"])
    join(L, "\n") * "\n"
end

function main()
    only = [split(replace(a, "--only=" => ""), ",") for a in ARGS if startswith(a, "--only=")]
    only = isempty(only) ? nothing : only[1]
    stands = [replace(a, "--stand=" => "") for a in ARGS if startswith(a, "--stand=")]
    isempty(stands) && (stands = STANDS)
    dir = joinpath(ROOT, "on"); mkpath(dir)
    cp(SRCDB, joinpath(dir, "stands.db"); force = true)
    # base cases on the first stand only; FFE cases on every stand
    for cn in stands, (cases, isffe) in ((BASE, false), (FFE, true))
        (!isffe && cn != STANDS[1]) && continue
        for (name, f) in cases
            only !== nothing && !(name in only) && continue
            main_l, fmin = f(INVY[cn])
            stem = "on_$(cn)_$(isffe ? "ffe_" : "")$(name)"
            work = mktempdir(); cp(SRCDB, joinpath(work, "stands.db"))
            write(joinpath(work, stem * ".key"), casekey(cn, name, main_l, fmin))
            try
                run(pipeline(Cmd(`$(ORACLE) --keywordfile=$(joinpath(work, stem * ".key"))`; dir = work, ignorestatus = true),
                             stdout = devnull, stderr = devnull))
            catch e
                @warn "oracle failed" stem e
            end
            sumf = joinpath(work, stem * ".sum")
            if !isfile(sumf); @warn "NO live .sum" stem; continue; end
            cp(joinpath(work, stem * ".key"), joinpath(dir, stem * ".key"); force = true)
            cp(sumf, joinpath(dir, stem * ".live.sum"); force = true)
            for f in readdir(dir); startswith(f, stem * ".") && endswith(f, ".csv") && rm(joinpath(dir, f)); end
            tabs = db_tables(joinpath(work, "out.db"))
            for t in tabs
                hdr, rows = db_table_rows(joinpath(work, "out.db"), t)
                write_csv(joinpath(dir, "$(stem).$(t).csv"), hdr, rows)
            end
            write(joinpath(dir, stem * ".tables"), join(tabs, "\n") * "\n")
            outf = joinpath(work, stem * ".out")
            isdir("/workspace/.eastext/outs") && isfile(outf) && cp(outf, "/workspace/.eastext/outs/$(stem).out"; force = true)
            errs = isfile(outf) ? [strip(l) for l in eachline(outf) if (occursin(r"ERROR|WARNING|\*\*\*\*\*\*|FVS\d\d", l) && !occursin("STANDARD ERRORS", l))] : String[]
            println("$stem: $(length(tabs)) tables", isempty(errs) ? "" : "  FVS messages: " * join(unique(first.(errs, 90)), " | "))
            rm(work; recursive = true, force = true)
        end
    end
    open(joinpath(ROOT, "PROVENANCE.toml"), "w") do io
        println(io, "# kwcov_on goldens — LIVE oracle: stock-flag private rebuild (ORACLE_SOURCE_AUDIT §9)")
        println(io, "date = \"", string(Base.Libc.strftime("%Y-%m-%d", time())), "\"")
        println(io, "oracle_on = \"$(ORACLE)\"\nsha256_on = \"", bytes2hex(open(sha256, ORACLE)), "\"")
        println(io, "command = \"cd <dir> && FVSon_tiered --keywordfile=<stem>.key  (generator: test/harness/kwcov_on/make_cases.jl; stands from test/fixtures/tiered/on/stands.db)\"")
    end
end
main()
