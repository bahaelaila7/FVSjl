# make_cases.jl — build the eastern (CS/LS/NE) FFE keyword-isolation fixtures from the LIVE oracle.
#
# Usage: julia --project=. test/harness/kwcov_ffe/make_cases.jl [VARIANT ...] [--only=name1,name2]
#
# One stand per variant (the canonical kwcov stand: cst01/lst01/net01 data, test/fixtures/kwcov/<v>_simfire.tre) × one
# isolated FFE keyword (plus whatever the keyword needs to be observable: a SIMFIRE, a thinning that creates activity
# fuel, ...). Every case also requests the FFE reports and their DBS tables (the two-DATABASE-block recipe), so the
# keyword is compared in the .sum rows AND in every FFE table cell. Goldens = the stock-flag LIVE oracle
#   /workspace/.<v>work/tiered/FVS<v>_tiered   (ORACLE_SOURCE_AUDIT_2026-09-19.md §9; build /workspace/.eastoracle/build_tiered.sh)
# run exactly as test/harness/<v>_oracle.sh runs the shared binary:  cd <dir> && FVS<v>_tiered --keywordfile=<stem>.key
# Written to test/fixtures/kwcov_ffe/<v>/<v>_<case>.{key,tre,live.sum,<Table>.csv,tables}; PROVENANCE.toml records the
# oracle path + sha256 + date.
using SHA
include(joinpath(@__DIR__, "..", "tiered", "tiered_common.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "kwcov_ffe"))
const KWCOV = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "kwcov"))
const ORACLE = Dict(v => "/workspace/.$(v)work/tiered/FVS$(v)_tiered" for v in ("cs", "ls", "ne"))

# keyword card: name in cols 1-10, then 10-wide right-justified fields; `nothing` = blank field
kw(name, vals...) = rpad(name, 10) * join([v === nothing ? " "^10 : lpad(string(v), 10) for v in vals])

const HEAD = Dict(
    "cs" => ["DESIGN                                        11.0       1.0",
             "STDINFO          905                60.0     315.0      30.0      10.0       40.", "SITECODE          19        60"],
    "ls" => ["DESIGN                                        11.0       1.0",
             "STDINFO        903.0       1.0      60.0     315.0      30.0      10.0", "SITECODE           2        60"],
    "ne" => ["DESIGN                                        11.0       1.0",
             "STDINFO        922.0                60.0     315.0      30.0      20.0", "SITECODE          13        56"])
const TREEFMT = ["TREEFMT", "(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,",
                 "T52,I2,T66,5I1,T54,7I1,T75,F3.0)", "TREEDATA"]

# FFE report block appended to every case's FMIN (all reports the eastern builds link) + DBS toggles (block 2)
const REPORTS = ["FUELOUT", "FUELREPT", "MORTREPT", "BURNREPT", "POTFIRE", "SNAGSUM", "SNAGOUT", "CARBREPT",
                 "DWDVLOUT", "DWDCVOUT"]
const TOGGLES = ["SUMMARY", "FUELSOUT", "FUELREDB", "BURNREDB", "MORTREDB", "SNAGSUDB", "SNAGOUDB", "CARBREDB",
                 "POTFIRDB", "DWDVLDB", "DWDCVDB"]

sf(yr = 2010, w = 15, m = 2, t = 70, mort = 1, pct = 100, seas = 2) = kw("SIMFIRE", yr, w, m, t, mort, pct, seas)
thin = kw("THINBTA", 2000, 60)

# name => (main-context lines before FMIN, FMIN lines [keyword cards, supplemental records inline])
const CASES = Pair{String,Tuple{Vector{String},Vector{String}}}[
    "reports"        => ([], [sf()]),
    "simfire_s4"     => ([], [sf(2010, 25, 1, 90, 1, 40, 4)]),
    "simfire_nomort" => ([], [sf(2010, 10, 3, 60, 0, 100, 1)]),
    "flameadj_mult"  => ([], [sf(), kw("FLAMEADJ", 2010, 1.5)]),
    "flameadj_len"   => ([], [sf(), kw("FLAMEADJ", 2010, 1, 7, 25, 15)]),
    "moisture"       => ([], [sf(), kw("MOISTURE", 2010, 7, 8, 10, 13, 60, 85, 110)]),
    "moisture_herb"  => ([], [sf(), kw("MOISTURE", 2010, 7, 8, 10, 13, 60, 85)]),
    "fuelmove"       => ([], [kw("FUELMOVE", 2000, 6, 11, 1.5), kw("FUELMOVE", 2000, 4, 3, 0, 0.5), sf()]),
    "fuelmove_imp"   => ([], [kw("FUELMOVE", 2000, 0, 10, 2.0), kw("FUELMOVE", 2005, 10, 0, 0, 0.3, 0.5), sf()]),
    "pileburn"       => ([thin], [kw("PILEBURN", 2000, 1), sf()]),
    "pileburn_jack"  => ([thin], [kw("PILEBURN", 2000, 2, 100, 30, 60, 0), sf()]),
    "fueltret"       => ([thin], [kw("FUELTRET", 2000, 1, 2), sf(2003)]),
    "fueltret2"      => ([thin], [kw("FUELTRET", 2000, 2, 3), sf(2003)]),
    "potfmois"       => ([], [kw("POTFMOIS", 1, 6, 7, 9, 12, 55, 80, 100), kw("POTFMOIS", 2, 12, 13, 15, 18, 80, 120, 140)]),
    "potfwind"       => ([], [kw("POTFWIND", 25, 12)]),
    "potftemp"       => ([], [kw("POTFTEMP", 85, 60)]),
    "potfseas"       => ([], [kw("POTFSEAS", 2, 3)]),
    "potfpab"        => ([], [kw("POTFPAB", 80, 30)]),
    "snagbrk"        => ([], [kw("SNAGBRK", "ALL", 12, 10, 30, 25), sf()]),
    "snagbrk_sp"     => ([], [kw("SNAGBRK", 19, 8, 6, 20, 15), sf()]),
    "snagfall"       => ([], [kw("SNAGFALL", "ALL", 2.0, 0.8), sf()]),
    "snagfall_sp"    => ([], [kw("SNAGFALL", 19, 0.5, 0.0), sf()]),
    "snagdcay"       => ([], [kw("SNAGDCAY", "ALL", 1.5), sf()]),
    "snagpsft"       => ([], [kw("SNAGPSFT", "ALL", 0.9), sf()]),
    "snagpbn"        => ([], [kw("SNAGPBN", 0.8, 0.5, 3, 12, 1), sf()]),
    "snaginit"       => ([], [kw("SNAGINIT", "ALL", 12, 40, 40, 5, 20), sf()]),
    "salvage"        => ([], [sf(2000), kw("SALVAGE", 2005, 0, 999, 5, 1, 0.9, 0)]),
    "salvage_sp"     => ([], [sf(2000), kw("SALVSP", 2005, 19, 1), kw("SALVAGE", 2005, 0, 999, 5, 2, 0.8, 1)]),
    "drought"        => ([], [kw("DROUGHT", 2005, 5), sf()]),
    "fuelinit"       => ([], [kw("FUELINIT", 0.3, 0.8, 2, 3, 3.5, 4, 1, 0.5, 0.2, 1.5, 10), sf()]),
    "fuelsoft"       => ([], [kw("FUELSOFT", 0.1, 0.2, 0.5, 1, 1.5, 2, 1, 0.5, 0.2), sf()]),
    "fuelfoto"       => ([], [kw("FUELFOTO", 28, 5), sf()]),
    "fuelfoto_1"     => ([], [kw("FUELFOTO", 1, 3), sf()]),                        # reference 1 (total + rotten proportion tables)
    "fuelfoto_12"    => ([], [kw("FUELFOTO", 12, 40), sf()]),
    "fuelfoto_bad"   => ([], [kw("FUELFOTO", 4, 2), sf()]),                        # reference 4 is not a reference: warning, default fuels
    "fuelfoto_big"   => ([], [kw("FUELFOTO", 28, 99), sf()]),                      # code past the reference's last photo
    "fuelfoto_init"  => ([], [kw("FUELFOTO", 2, 10), kw("FUELINIT", nothing, 9.0), kw("FUELSOFT", 0.5), sf()]),   # FUELINIT/FUELSOFT override the photo
    "fueldcay"       => ([], [kw("FUELDCAY", 1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), sf()]),
    "fuelmult"       => ([], [kw("FUELMULT", 1.2, 1.5), sf()]),
    "duffprod"       => ([], [kw("DUFFPROD", 1, 0.3, 0.4, 0.5, 0.6, 0.7), sf()]),
    "fuelpool"       => ([], [kw("FUELPOOL", "ALL", 2), sf()]),
    "carbcalc"       => ([], [kw("CARBCALC", 1, 1, 0.05, 8, 10), sf()]),
    "carbcalc_t"     => ([], [kw("CARBCALC", 0, 2)]),
    "carbcut"        => ([thin], [kw("CARBCUT")]),
    "fmodlist"       => ([], [kw("FMODLIST", 2010, 5, 0), sf()]),
    "fuelmodl"       => ([], [kw("FUELMODL", 2010, 8, 0.6, 10, 0.4), "", sf()]),
    "defulmod"       => ([], [kw("DEFULMOD", 0, 12, 2000, 1800, 1500, 1500), join(lpad.(string.([0.2, 0.6, 0.9, 0.4, 0.25, 1.5, 0.3]), 10)),
                              kw("FUELMODL", 2010, 12, 1.0), "", sf()]),
    "firecalc"       => ([], [kw("FIRECALC", 2000, 1, 0), sf()]),
    "firecalc_ml"    => ([], [kw("FIRECALC", 2000, 2, 1, 2100, 1900, 1600, 0.12, 0.8, 8100), sf()]),
    "firecalc_13"    => ([], [kw("FIRECALC", 2000, 1, 0), sf()]),                  # new logic over the original 13 only (IFMSET 0)
    "firecalc_40"    => ([], [kw("FIRECALC", 2000, 1, 1), sf()]),                  # ... the 40 Scott-Burgan only (IFMSET 1)
    "firecalc_sav"   => ([], [kw("FIRECALC", 2000, 1, 2, 2400, 2000, 1700, 0.08, 0.6), sf()]),
    "firecalc_off"   => ([], [kw("FIRECALC", 2000, 1, 2), kw("FMODLIST", 2000, 8, 1), kw("FMODLIST", 2000, 9, 1), kw("FMODLIST", 2000, 11, 1), sf()]),
    "firecalc_on"    => ([], [kw("FIRECALC", 2000, 1, 0), kw("FMODLIST", 2000, 5, 0), kw("FMODLIST", 2000, 143, 0), sf()]),
    "firecalc_late"  => ([], [kw("FIRECALC", 2010, 2, 1), sf(2010, 15, 2, 70, 1, 100, 2)]),   # modelled loads set the same year the fire burns
    "firecalc_stat"  => ([], [kw("FIRECALC", 2000, 1, 2), kw("STATFUEL"), sf()]),
    "statfuel"       => ([], [kw("STATFUEL"), sf()]),
    "soilheat"       => ([], [kw("SOILHEAT", 2010, nothing, 2), sf()]),
    "mortclas"       => ([], [kw("MORTCLAS", 4, 8, 12, 16, 20, 24, 30), sf()]),
    "fmortmlt"       => ([], [kw("FMORTMLT", 2010, 0.5, "ALL", 0, 999), sf()]),
    "fmortmlt_early" => ([], [kw("FMORTMLT", 2000, 0.5, "ALL", 0, 999), sf()]),
    "fmortmlt_dbh"   => ([], [kw("FMORTMLT", 2010, 0.2, "ALL", 0, 8), kw("FMORTMLT", 2010, 1.5, "ALL", 8, 999), sf()]),
    "cancalc"        => ([], [kw("CANCALC", 1, 20, 0, 6, 100), sf()]),
    "canfprof"       => ([], [kw("CANFPROF"), sf()]),
]

function casekey(v, name, main, fmin)
    L = String["STDIDENT", "$(v) ffe $(name)", "NOAUTOES"]
    append!(L, HEAD[v]); append!(L, ["INVYEAR       1990.0", "NUMCYCLE         5.0", "ECHOSUM"])
    append!(L, ["DATABASE", "DSNout", "out.db", "END"])
    append!(L, main)
    push!(L, "FMIN"); append!(L, fmin); append!(L, REPORTS); push!(L, "END")
    push!(L, "DATABASE"); append!(L, TOGGLES); push!(L, "END")
    append!(L, TREEFMT); append!(L, ["PROCESS", "STOP"])
    join(L, "\n") * "\n"
end

# NB: TREEDATA reads the .tre in the stand's own cwd/stem, so TREEDATA follows TREEFMT then PROCESS, as in kwcov.
function main()
    args = ARGS
    vs = [a for a in args if !startswith(a, "--")]; isempty(vs) && (vs = ["cs", "ls", "ne"])
    only = [split(replace(a, "--only=" => ""), ",") for a in args if startswith(a, "--only=")]
    only = isempty(only) ? nothing : only[1]
    for v in vs
        dir = joinpath(ROOT, v); mkpath(dir)
        for (name, (main, fmin)) in CASES
            only !== nothing && !(name in only) && continue
            stem = "$(v)_$(name)"
            keytxt = casekey(v, name, main, fmin)
            work = mktempdir()
            write(joinpath(work, stem * ".key"), keytxt)
            cp(joinpath(KWCOV, "$(v)_simfire.tre"), joinpath(work, stem * ".tre"))
            try
                run(pipeline(Cmd(`$(ORACLE[v]) --keywordfile=$(joinpath(work, stem * ".key"))`; dir = work, ignorestatus = true),
                             stdout = devnull, stderr = devnull))
            catch e
                @warn "oracle failed" stem e
            end
            sumf = joinpath(work, stem * ".sum")
            if !isfile(sumf)
                @warn "NO live .sum" stem; continue
            end
            cp(joinpath(work, stem * ".key"), joinpath(dir, stem * ".key"); force = true)
            cp(joinpath(work, stem * ".tre"), joinpath(dir, stem * ".tre"); force = true)
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
        end
    end
    open(joinpath(ROOT, "PROVENANCE.toml"), "w") do io
        println(io, "# kwcov_ffe goldens — LIVE oracle: stock-flag private rebuilds (ORACLE_SOURCE_AUDIT §9)")
        println(io, "date = \"", string(Base.Libc.strftime("%Y-%m-%d", time())), "\"")
        for v in ("cs", "ls", "ne")
            isfile(ORACLE[v]) && println(io, "oracle_$(v) = \"$(ORACLE[v])\"\nsha256_$(v) = \"", bytes2hex(open(sha256, ORACLE[v])), "\"")
        end
        println(io, "command = \"cd <dir> && FVS<v>_tiered --keywordfile=<stem>.key  (generator: test/harness/kwcov_ffe/make_cases.jl)\"")
    end
end
main()
