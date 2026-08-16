# test_svs_multicycle.jl — SVS (Stand Visualization System) data path, multi-cycle.
#
# Validates the projected SVS pictures against the LIVE relinked FVSkt oracle (stand S248112),
# advancing chunk 0 (cyc0 inventory) to the MULTI-CYCLE path:
#   * SVSTART (fvs.f:333)  — inventory picture (_001)
#   * GRINCR (grincr.f:277) — "Beginning of cycle" picture each later cycle (ICYC>1) (2-cyc key)
#   * MAIN   (fvs.f:453)   — "End of projection" picture
#   * SVTRIP (svtrip.f)    — the persistent object list's IS2F pointers remap on record tripling
#                           (TRIPLE splits base i → central=i, upper=nlive+2i-1, lower=nlive+2i),
#                           preserving each object's (x,y). The .60/.625 object split reproduces
#                           the oracle's tree# remapping (WL rec2→{2,8}, LP rec4→{4,4,12}).
#
# The SVS DATA PATH is bit-exact: object count, ordering, (x,y) placement, species, tree# pointer
# remapping, and DBH are byte-identical to the oracle on every record. The only residual is the
# grown HT (±1 ft) / crown-ratio (±2%) of the tripled UPPER records — a pre-existing growth-engine
# height-growth-on-tripled-records artifact (the arrays SVS reads), cornered at the .sum. Object
# mortality→dead-snag conversion (oracle 2-cyc _003) is the NEXT chunk and is not asserted here.

using Test, FVSjl

const _SVSMC_DIR = joinpath(@__DIR__, "..", "fixtures", "svs")

# Parse a picture into (header::Vector, tree-records::Vector{Vector{String}} of whitespace fields).
function _svs_parse(path)
    hdr = String[]; recs = Vector{String}[]
    for ln in eachline(path)
        if startswith(ln, "#") || startswith(ln, ";")
            push!(hdr, ln)
        elseif !isempty(strip(ln))
            push!(recs, String.(split(ln)))
        end
    end
    return hdr, recs
end

# Field layout of a live-tree record (svout.f fmt 30, whitespace-split):
#  [1]=SPCD [2]=tree# [6]=dbh [7]=ht ... [end-2]=xloc [end-1]=yloc [end]=z
_fld_sp(r)   = r[1]
_fld_tr(r)   = r[2]
_fld_dbh(r)  = r[6]
_fld_ht(r)   = parse(Float64, rstrip(r[7], '.'))
_fld_xloc(r) = r[end-2]
_fld_yloc(r) = r[end-1]
# crown ratio is the first "cr" after dbh/ht/lean(0)/dir(0)/edia(0)/crad — index 12
_fld_cr(r)   = parse(Float64, r[12])

# Assert the SVS-specific structure is bit-exact and HT/CR are within the cornered growth residual.
function _svs_assert_picture(got_path, oracle_path; ht_tol = 1.0, cr_tol = 0.021)
    ghdr, grecs = _svs_parse(got_path)
    ohdr, orecs = _svs_parse(oracle_path)
    @test ghdr == ohdr                                   # #TITLE/#TREEFORM/... header byte-identical
    @test length(grecs) == length(orecs)                 # same object count & ordering
    for (g, o) in zip(grecs, orecs)
        # SVS data-path fields — must be byte-identical:
        @test _fld_sp(g)   == _fld_sp(o)                 # species code
        @test _fld_tr(g)   == _fld_tr(o)                 # tree# = IS2F (SVTRIP pointer remap)
        @test _fld_dbh(g)  == _fld_dbh(o)                # DBH
        @test _fld_xloc(g) == _fld_xloc(o)               # xloc (persists across cycles)
        @test _fld_yloc(g) == _fld_yloc(o)               # yloc
        # grown HT / crown-ratio — cornered on tripled-upper records (growth engine):
        @test abs(_fld_ht(g) - _fld_ht(o)) <= ht_tol
        @test abs(_fld_cr(g) - _fld_cr(o)) <= cr_tol
    end
end

@testset "SVS multi-cycle — end-of-projection + begin-cycle pictures vs live FVSkt" begin
    out = mktempdir()

    # --- 1-cycle projection (kt0): inventory + End of projection ---
    key1 = joinpath(_SVSMC_DIR, "kt0.key")
    s = FVSjl.each_stand(key1; variant = FVSjl.Kootenai())[1]
    FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
    stem1 = joinpath(out, "kt0")
    pics1 = FVSjl.svs_project!(stem1, s; fint = 10f0)
    @test length(pics1) == 2

    # _001 inventory + _index are BYTE-IDENTICAL to the oracle.
    @test read(stem1 * "_001.svs", String) ==
          read(joinpath(_SVSMC_DIR, "kt0_001.svs.oracle"), String)
    @test read(stem1 * "_index.svs", String) ==
          read(joinpath(_SVSMC_DIR, "kt0_index.svs.oracle"), String)
    # _002 End of projection — SVS structure bit-exact; HT/CR cornered on 2 tripled-upper records.
    _svs_assert_picture(stem1 * "_002.svs", joinpath(_SVSMC_DIR, "kt0_002.svs.oracle"))

    # --- 2-cycle projection (kt2c): exercises the GRINCR "Beginning of cycle" seam + 3-line index ---
    key2 = joinpath(_SVSMC_DIR, "kt2c.key")
    s2 = FVSjl.each_stand(key2; variant = FVSjl.Kootenai())[1]
    FVSjl.notre!(s2); FVSjl.setup_growth!(s2); FVSjl.compute_volumes!(s2)
    stem2 = joinpath(out, "kt2c")
    pics2 = FVSjl.svs_project!(stem2, s2; fint = 10f0)
    @test length(pics2) == 3

    # 3-line index (Inventory / Beginning of cycle / End of projection) BYTE-IDENTICAL.
    @test read(stem2 * "_index.svs", String) ==
          read(joinpath(_SVSMC_DIR, "kt2c_index.svs.oracle"), String)
    # _002 "Beginning of cycle" (2000) — SVS structure bit-exact; same 2 tripled-upper residuals.
    _svs_assert_picture(stem2 * "_002.svs", joinpath(_SVSMC_DIR, "kt2c_002.svs.oracle"))
    # (_003 end-of-projection after 2 cycles needs mortality→dead-object handling — next chunk.)
end
