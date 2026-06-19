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
    path = isfile(trepath) ? trepath :
           (isfile(uppercase(trepath)) ? uppercase(trepath) : trepath)
    isfile(path) || return 0

    fields = parse_tree_format(s.control.tree_format)
    isempty(fields) && return 0
    t = s.trees
    p = s.plot
    plot_ids = Int32[]            # unique record plot numbers (IPVEC)
    n0 = t.n

    for line in eachline(path)
        isempty(strip(line)) && continue
        startswith(strip(line), "*") && continue
        occursin("-999", line) && break
        rec = parse_tree_record(fields, line)
        rec === nothing && break

        # Subplot index (IPVEC/ITRE) is assigned to EVERY record before the dead /
        # non-stockable exclusions (intree.f: plot counting precedes those checks),
        # so the numbering matches FVS even when some records aren't kept.
        pj = findfirst(==(rec.plot), plot_ids)
        if pj === nothing
            push!(plot_ids, rec.plot); pj = length(plot_ids)
        end

        # IMC1 == 8 marks a non-stockable plot record — not a tree (intree.f:368)
        rec.mort_code == 8 && continue
        # Dead trees (history/ITH 6-9) are excluded from the live stand; they are
        # partitioned into the dead pool for mortality reporting in C4 (intree.f:516).
        # For now (live cycle-0 statistics) they are skipped.
        (6 <= rec.history <= 9) && continue

        i = t.n + 1
        i > MAXTRE && break

        idx, fmt = resolve_species(rec.species_code, s.variant, s.species)
        p.sp_format[idx] = fmt
        p.sp_format_default <= 0 && (p.sp_format_default = fmt)

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

        t.plot_id[i] = Int32(pj)            # subplot index from the full registration above
        t.n = i
    end

    s.control.ntrees_active = Int32(t.n)
    return t.n - n0
end
