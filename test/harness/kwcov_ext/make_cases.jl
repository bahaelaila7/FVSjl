# make_cases.jl — extension-STUB parity fixtures for CS/LS/NE/ON from the LIVE oracles.
#
# Usage: julia --project=. test/harness/kwcov_ext/make_cases.jl [VARIANT ...] [--only=name1,name2]
#
# FVS{cs,ls,ne,on}_buildDir link the ex*.f stubs for MPB (exmpb), DFTM (exdftm), COVER (excov), BRUST (exbrus), MISTOE (exmist), DFB (exdfb),
# RDIN (exrd), BMIN (exbm), CLIMATE (exclim) and, in ON, FERTILIZ (exfert), BGCIN (exbgc), ECON (execon); WSBW (exbudl) is a silent stub. Each case is
# the canonical stand with ONE extension keyword block (opener [+ unknown sub-keyword] + END), run on the stock-flag live oracle
#   /workspace/.<v>work/tiered/FVS<v>_tiered   (ORACLE_SOURCE_AUDIT_2026-09-19.md §9)
# with a DATABASE block (SUMMARY) so the FVS_Error rows (FVS11 / FVS01 per record) are compared with the .sum. Written to
# test/fixtures/kwcov_ext/<v>/<v>_<case>.{key,tre,live.sum,<Table>.csv,tables}; ON uses stands.db (stand 21769405010661).
using SHA
include(joinpath(@__DIR__, "..", "tiered", "tiered_common.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "kwcov_ext"))
const KWCOV = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "kwcov"))
const SRCDB = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "tiered", "on", "stands.db"))
const ORACLE = Dict(v => "/workspace/.$(v)work/tiered/FVS$(v)_tiered" for v in ("cs", "ls", "ne"))
ORACLE["on"] = "/workspace/.onwork/tiered/FVSon_tiered"
const ONSTAND = "21769405010661"

const HEAD = Dict(
    "cs" => ["DESIGN                                        11.0       1.0",
             "STDINFO          905                60.0     315.0      30.0      10.0       40.", "SITECODE          19        60"],
    "ls" => ["DESIGN                                        11.0       1.0",
             "STDINFO        903.0       1.0      60.0     315.0      30.0      10.0", "SITECODE           2        60"],
    "ne" => ["DESIGN                                        11.0       1.0",
             "STDINFO        922.0                60.0     315.0      30.0      20.0", "SITECODE          13        56"])
const TREEFMT = ["TREEFMT", "(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,",
                 "T52,I2,T66,5I1,T54,7I1,T75,F3.0)", "TREEDATA"]

# name => block lines placed after NUMCYCLE
# BMIN (exbm.f) writes to the unopened JOPPRT unit and the oracle dies with SIGSEGV — no live row to compare, so it is not a case.
const OPENERS = ["MPB", "DFTM", "WSBW", "COVER", "BRUST", "MISTOE", "DFB", "RDIN", "BGCIN", "ECON", "CLIMATE", "ORGANON"]
const CASES = Pair{String,Vector{String}}[]
for o in OPENERS
    lo = lowercase(o)
    push!(CASES, "$(lo)_end" => [o, "END"])
    push!(CASES, "$(lo)_sub" => [o, "FOO       1.0", "END"])
end
push!(CASES, "fertiliz" => ["FERTILIZ      2000.0     200.0"])
push!(CASES, "fertiliz_none" => ["FERTILIZ"])
push!(CASES, "mistoe_real" => ["MISTOE", "MISTPINF", "END"])
push!(CASES, "climate_real" => ["CLIMATE", "CLIMDATA", "END"])
push!(CASES, "rdin_real" => ["RDIN", "RRTYPE           1", "SAREA            20", "RRINIT           1.0", "END"])

function casekey(v, cn, name, body)
    L = String["STDIDENT", v == "on" ? cn : "$(v) ext $(name)"]
    if v == "on"
        append!(L, ["DATABASE", "DSNout", "out.db", "DSNin", "stands.db", "StandSQL",
                    "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL", "TreeSQL",
                    "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL", "END", "NUMCYCLE         5.0"])
        append!(L, body)
        append!(L, ["DATABASE", "SUMMARY", "END", "ECHOSUM", "PROCESS", "STOP"])
    else
        push!(L, "NOAUTOES"); append!(L, HEAD[v])
        append!(L, ["INVYEAR       1990.0", "NUMCYCLE         5.0", "ECHOSUM", "DATABASE", "DSNout", "out.db", "END"])
        append!(L, body)
        push!(L, "DATABASE"); push!(L, "SUMMARY"); push!(L, "END")
        append!(L, TREEFMT); append!(L, ["PROCESS", "STOP"])
    end
    join(L, "\n") * "\n"
end

function main()
    args = ARGS
    vs = [a for a in args if !startswith(a, "--")]; isempty(vs) && (vs = ["cs", "ls", "ne", "on"])
    only = [split(replace(a, "--only=" => ""), ",") for a in args if startswith(a, "--only=")]
    only = isempty(only) ? nothing : only[1]
    for v in vs
        dir = joinpath(ROOT, v); mkpath(dir)
        v == "on" && cp(SRCDB, joinpath(dir, "stands.db"); force = true)
        for (name, body) in CASES
            only !== nothing && !(name in only) && continue
            stem = "$(v)_$(name)"
            work = mktempdir()
            write(joinpath(work, stem * ".key"), casekey(v, ONSTAND, name, body))
            v == "on" ? cp(SRCDB, joinpath(work, "stands.db")) : cp(joinpath(KWCOV, "$(v)_simfire.tre"), joinpath(work, stem * ".tre"))
            try
                run(pipeline(Cmd(`$(ORACLE[v]) --keywordfile=$(joinpath(work, stem * ".key"))`; dir = work, ignorestatus = true),
                             stdout = devnull, stderr = devnull))
            catch e
                @warn "oracle failed" stem e
            end
            sumf = joinpath(work, stem * ".sum")
            if !isfile(sumf); @warn "NO live .sum" stem; continue; end
            cp(joinpath(work, stem * ".key"), joinpath(dir, stem * ".key"); force = true)
            v == "on" || cp(joinpath(work, stem * ".tre"), joinpath(dir, stem * ".tre"); force = true)
            cp(sumf, joinpath(dir, stem * ".live.sum"); force = true)
            for f in readdir(dir); startswith(f, stem * ".") && endswith(f, ".csv") && rm(joinpath(dir, f)); end
            tabs = db_tables(joinpath(work, "out.db"))
            for t in tabs
                hdr, rows = db_table_rows(joinpath(work, "out.db"), t)
                write_csv(joinpath(dir, "$(stem).$(t).csv"), hdr, rows)
            end
            write(joinpath(dir, stem * ".tables"), join(tabs, "\n") * "\n")
            outf = joinpath(work, stem * ".out")
            errs = isfile(outf) ? [strip(l) for l in eachline(outf) if occursin(r"FVS\d\d", l) && !occursin("STANDARD ERRORS", l)] : String[]
            println("$stem: $(length(tabs)) tables  ", join(unique(first.(errs, 60)), " | "))
            rm(work; recursive = true, force = true)
        end
    end
    open(joinpath(ROOT, "PROVENANCE.toml"), "w") do io
        println(io, "# kwcov_ext goldens — LIVE oracles: stock-flag private rebuilds (ORACLE_SOURCE_AUDIT §9)")
        println(io, "date = \"", string(Base.Libc.strftime("%Y-%m-%d", time())), "\"")
        for v in ("cs", "ls", "ne", "on")
            isfile(ORACLE[v]) && println(io, "oracle_$(v) = \"$(ORACLE[v])\"\nsha256_$(v) = \"", bytes2hex(open(sha256, ORACLE[v])), "\"")
        end
        println(io, "command = \"cd <dir> && FVS<v>_tiered --keywordfile=<stem>.key  (generator: test/harness/kwcov_ext/make_cases.jl)\"")
    end
end
main()
