# filter_digworthy.jl — extract the DIG-WORTHY discrepancy rows from a cycle ledger CSV, EXCLUDING stands
# whose (ECOREGION-prefix, signature) has already been cornered as an accepted taxonomy class
# (docs/fia_cornered_clusters.tsv). This lets the full-population coverage sweep advance to NEW strata instead
# of re-pausing every ~200 stands on an already-cornered cluster (e.g. the Appalachian 221H compounded-ULP /
# count-straddle cluster). ESCALATION GUARD: an UNCLASSIFIED signature, or a STRUCTURE/density blow-up
# (worst_col ∈ density/structure cols with max_rel ≥ 15%), is NEVER dropped — it re-surfaces for a manual trace
# even inside a cornered cluster, so a genuinely new bug in that geography still pauses the sweep.
#
# Usage: julia --project=. test/harness/fia/filter_digworthy.jl <cycle.csv> <VARIANT> [cornered.tsv]
#   prints the dig-worthy CSV rows (no header) to stdout; the orchestrator appends them to the dig-queue.
import SQLite, DBInterface
const MASTER = "/workspace/SQLite_FIADB_ENTIRE.db"

# base dig-worthy rule (same as run_expand_cycle.sh's historical awk): a MATERIAL, potentially-real-bug class.
const DIG_SIGS = Set(["UNCLASSIFIED", "volume_persistent", "structure_densephase"])
# columns that force an escalation (never dropped even in a cornered cluster): the true density/structure
# cols. TopHt is DELIBERATELY EXCLUDED — dig-session #2c empirically cornered the AVHT40 top-height tie-break
# as a ULP primitive (no global single/double RDPSRT sort choice is bit-exact; stand-dependent), and TopHt
# divergences with density preserved ARE that primitive's signature. A TopHt-worst row in an ALREADY-CORNERED
# ecoregion is therefore dropped; in a NEW ecoregion it still surfaces (signature not yet in corners). A real
# height bug would also perturb BA/SDI/CCF/QMD or appear in a non-cornered geography, both of which still escalate.
# TCuFt (total cubic, threshold-FREE) is included: a large TCuFt divergence is a real volume-EQUATION bug
# signal (e.g. the slice-41 FORKOD zero-vol), not the threshold-crossing/count-straddle noise that BdFt/SCuFt/
# MCuFt (merch/board step-functions) show. TopHt is EXCLUDED (dig-2c: cornered AVHT40 tie-break ULP primitive).
const ESCALATE_COLS = Set(["TPA", "BA", "SDI", "CCF", "QMD", "TCuFt"])
const ESCALATE_REL = 15.0

is_dig(sig, worst_col, struct_pct, max_rel) =
    sig in DIG_SIGS || (worst_col == "TCuFt" && struct_pct < 1.0 && max_rel >= 5.0)

# escalation: never drop these even in a cornered cluster
is_escalation(sig, worst_col, max_rel) =
    sig == "UNCLASSIFIED" || (worst_col in ESCALATE_COLS && max_rel >= ESCALATE_REL)

function load_cornered(path)
    corners = Set{Tuple{String,String}}()
    isfile(path) || return corners
    for l in eachline(path)
        (isempty(strip(l)) || startswith(strip(l), "#")) && continue
        f = split(l, '\t')
        length(f) >= 2 && push!(corners, (strip(f[1]), strip(f[2])))
    end
    corners
end

function main(csv, v, cornerfile)
    corners = load_cornered(cornerfile)
    # parse candidate rows first (cheap), collect their CNs for a single ecoregion lookup
    rows = String[]; cand = Tuple{Int,String,String,String,Float64}[]  # (rowidx, cn, sig, worst_col, max_rel)
    for (i, l) in enumerate(eachline(csv))
        i == 1 && continue
        f = split(l, ',')
        length(f) < 15 && continue
        cn = f[3]; worst_col = f[7]; sig = strip(f[15])
        max_rel = something(tryparse(Float64, f[9]), 0.0)
        struct_pct = something(tryparse(Float64, f[11]), 0.0)
        is_dig(sig, worst_col, struct_pct, max_rel) || continue
        push!(rows, l); push!(cand, (length(rows), cn, sig, worst_col, max_rel))
    end
    isempty(cand) && return
    # one ecoregion lookup for all candidate CNs
    eco = Dict{String,String}()
    if !isempty(corners)
        inlist = join(["'" * replace(c[2], "'"=>"''") * "'" for c in cand], ",")
        db = SQLite.DB(MASTER)
        for r in DBInterface.execute(db, "SELECT STAND_CN,ECOREGION FROM FVS_STANDINIT_COND WHERE VARIANT='$v' AND STAND_CN IN ($inlist)")
            (r.STAND_CN === missing || r.ECOREGION === missing) && continue
            eco[String(r.STAND_CN)] = String(r.ECOREGION)
        end
        SQLite.close(db)
    end
    for (idx, cn, sig, worst_col, max_rel) in cand
        drop = false
        if !is_escalation(sig, worst_col, max_rel)
            e = get(eco, cn, "")
            for (pfx, csig) in corners
                if csig == sig && startswith(e, pfx); drop = true; break; end
            end
        end
        drop || println(rows[idx])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: filter_digworthy.jl <cycle.csv> <VARIANT> [cornered.tsv]")
    main(ARGS[1], ARGS[2], length(ARGS) >= 3 ? ARGS[3] : "docs/fia_cornered_clusters.tsv")
end
