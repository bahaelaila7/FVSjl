# =============================================================================
# treeinput.jl — load tree records into the stand (INTREE)
#
# Ported from: base/intree.f  (the record→array assignment + per-tree fixups).
#
# Reads the `.tre` data (via the current tree FORMAT), resolves each species code
# to an index, and fills the `TreeList` columns, applying the same per-record
# fixups FVS does: crown-class code → percent, birth-age flag, and the subplot
# index (IPVEC) used later for per-plot expansion.
#
# Scope note: site-tree detection, DBH/species screening and the DBS read path are
# not yet ported (not exercised by the basic SN stands); they slot in here later.
# =============================================================================

"""
    load_trees!(state, trepath) -> Int

Read tree records from `trepath` using `state.control.tree_format`, appending them
to `state.trees`. Returns the number of trees loaded. Mirrors intree.f.
"""
function load_trees!(s::StandState, trepath::AbstractString)
    t = s.trees
    p = s.plot
    plot_ids = Int32[]            # unique record plot numbers (IPVEC)
    dead = Tuple{Any,Int32,Int32}[]  # (record, species idx, subplot) for dead trees
    n0 = t.n

    # Source the records from the modern self-describing `.csv` when one sits beside
    # the `.tre` (named columns → no TREEFMT needed), else parse the fixed-column `.tre`
    # with the current TREEFMT. Both yield the same `Vector{TreeRecord}`.
    csvpath = first(splitext(trepath)) * ".csv"
    local records::Vector{TreeRecord}
    if isfile(csvpath)
        records = read_trees_csv(csvpath)
    else
        path = isfile(trepath) ? trepath :
               (isfile(uppercase(trepath)) ? uppercase(trepath) : trepath)
        isfile(path) || return 0
        fields = parse_tree_format(s.control.tree_format)
        isempty(fields) && return 0
        records = TreeRecord[]
        for line in eachline(path)
            isempty(strip(line)) && continue
            startswith(strip(line), "*") && continue
            occursin("-999", line) && break
            rec = parse_tree_record(fields, line)
            rec === nothing && break
            push!(records, rec)
        end
    end

    for rec in records
        # Subplot index (IPVEC/ITRE) is assigned to EVERY record before the dead /
        # non-stockable exclusions (intree.f: plot counting precedes those checks),
        # so the numbering matches FVS even when some records aren't kept.
        pj = findfirst(==(rec.plot), plot_ids)
        if pj === nothing
            push!(plot_ids, rec.plot); pj = length(plot_ids)
        end

        # IMC1 == 8 marks a non-stockable plot record — not a tree (intree.f:368)
        rec.mort_code == 8 && continue

        idx, fmt = resolve_species(rec.species_code, s.variant, s.species, s.coef)
        p.sp_format[idx] = fmt
        p.sp_format_default <= 0 && (p.sp_format_default = fmt)

        # Dead trees (history/ITH 6-9) are partitioned out of the live stand (intree.f:516):
        # collected here, stored after the live records so live stats use 1:n but the dead
        # remain available (mortality reporting; backdated calibration BA at current dbh).
        if 6 <= rec.history <= 9
            push!(dead, (rec, idx, Int32(pj)))
            continue
        end

        i = t.n + 1
        i > MAXTRE && break
        _store_tree!(t, i, rec, idx, Int32(pj))
        t.n = i
    end

    # append the dead records after the live ones (indices n+1 : n+ndead)
    t.ndead = 0
    for (rec, idx, pj) in dead
        i = t.n + t.ndead + 1
        i > MAXTRE && break
        _store_tree!(t, i, rec, idx, pj)
        t.ndead += 1
    end

    s.control.ntrees_active = Int32(t.n)
    return t.n - n0
end

# Fill all per-tree fields of record `rec` (resolved species `idx`, subplot `pj`)
# into TreeList slot `i`. Shared by the live and dead partitions of the loader.
function _store_tree!(t::TreeList, i::Int, rec, idx::Integer, pj::Int32)
    t.species[i]     = idx
    t.tree_id[i]     = rec.id
    t.tpa[i]         = rec.tpa
    t.history[i]     = rec.history
    t.dbh[i]         = rec.dbh
    t.diam_growth[i] = rec.diam_growth
    t.height[i]      = rec.height
    t.ht_growth[i]   = rec.ht_growth
    t.mort_code[i]   = rec.mort_code
    t.cut_code[i]    = rec.cut_code
    @inbounds for k in 1:6; t.damage[k, i]    = rec.damage[k];    end
    @inbounds for k in 1:5; t.pest_vars[k, i] = rec.pest_vars[k]; end

    # Per-tree volume defect from the damage codes (basdam.f): agent 25 = percent defect for
    # BOTH cubic and board, 26 = cubic only, 27 = board only; the paired value (clamped 0-99)
    # is the percent. Packed as CF·1e6 + BF·1e4 to match DEFECT (digits 11=CF%, 22=BF%); the
    # volume model reads ICDF=DEFECT/1e6 and IBDF=DEFECT/1e4 mod 100.
    defcf = 0; defbf = 0; ispecl = 0
    @inbounds for j in (1, 3, 5)
        ag = rec.damage[j]
        (ag == 25 || ag == 26) && (defcf += clamp(Int(rec.damage[j+1]), 0, 99))
        (ag == 25 || ag == 27) && (defbf += clamp(Int(rec.damage[j+1]), 0, 99))
        ag == 55 && (ispecl = clamp(Int(rec.damage[j+1]), 0, 99))   # special status code (basdam.f)
    end
    t.defect[i] = Int32(defcf * 1000000 + defbf * 10000)
    t.special[i] = Int32(ispecl)

    # Topkill / broken-top detection (intree.f:479-502): damage agent 96 or 97 in
    # any agent slot (odd positions), or an explicit broken-top height THT. Flag
    # with norm_ht=-1 (resolved to the predicted full height later, in CRATET);
    # trunc = break height ×100. A measured height at/below the break is dropped.
    t.trunc[i] = Int32(0); t.norm_ht[i] = Int32(0)
    topkilled = rec.damage[1] == 96 || rec.damage[1] == 97 ||
                rec.damage[3] == 96 || rec.damage[3] == 97 ||
                rec.damage[5] == 96 || rec.damage[5] == 97
    tht = rec.top_height
    if topkilled || tht > 0f0
        t.norm_ht[i] = Int32(-1)
        if tht > 0f0
            t.trunc[i] = round(Int32, tht * 100f0 + 0.5f0)
            t.height[i] <= tht && (t.height[i] = 0f0)
        end
    end

    # crown-class code → crown ratio percent (intree.f:311)
    icr = rec.crown_pct
    if icr > 0
        icr = icr < 10 ? icr * Int32(10) - Int32(5) : min(icr, Int32(99))
    end
    t.crown_pct[i] = icr
    t.crown_ratio[i] = Float32(icr)        # PCT (working crown ratio, percent)

    # birth-age flag (intree.f:190)
    if rec.birth_age <= 0f0
        t.birth_age[i] = 0f0; t.age_known[i] = false
    else
        t.birth_age[i] = rec.birth_age; t.age_known[i] = true
    end
    t.plot_id[i] = pj
    return t
end
