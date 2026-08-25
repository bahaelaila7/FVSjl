# test_ca_so_forkod_crash.jl — CA/SO forest-index-overflow crash regression guard.
#
# The western full-population FIA sweep SIGSEGV'd (not a MAXTRE overflow) on CA and SO stands whose FVS
# LOCATION code maps, in the JFOR table, to an IFOR index that FVS's forkod.f "FOREST MAPPING CORRECTION"
# remaps DOWN into the 1..10 range that the forest-dimensioned growth arrays (MAPLOC/DGFOR/...) are sized
# for. The jl ports omitted those corrections, so IFOR stayed 11 (CA) / 9,11 (SO) and indexed one column
# past the arrays — an out-of-bounds read that segfaults with @inbounds (a clean BoundsError under
# --check-bounds=yes: "13×10 Matrix at index [1,11]").
#
#   CA ca/forkod.f: TRINITY NF (518, JFOR idx 11) → SHASTA-TRINITY (514, idx 5).   [867 CA stands use 518]
#   SO so/forkod.f: SHASTA NF (514, idx 9) → KLAMATH (505, idx 4);
#                   INDUSTRY LANDS (702, idx 11) → 701 (idx 8).                     [702: 1906 SO stands]
#
# A third, co-located crash on SO whitebark-pine (WB, sp 16) small trees: so/small_tree_growth.jl read
# `t.plot[i]` (no such TreeList field) and `p.point_ccf` (point_ccf lives on s.density) — a FieldError,
# guarded end-to-end by running a WB-bearing stand below.
#
# Guards: (1) the pure forkod remap (no DB), and (2) the three previously-crashing stands run to a real
# projection row through the full engine (self-contained 3-stand FIA sub-DB fixture).
using FVSjl, Test

const _FK_DB = joinpath(@__DIR__, "..", "fixtures", "ca_so_forkod", "ca_so_forkod.db")

# minimal PlotData-like carrier for the pure forkod unit (forkod only touches these fields)
mutable struct _FK
    user_forest_code::Int32
    forest_idx::Int32
    geo_location::Int32
end

function _keytext(cn)
    join(["STDIDENT", cn, "DATABASE", "DSNin", abspath(_FK_DB),
          "StandSQL", "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
          "TreeSQL", "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
          "END", "NUMCYCLE         5.0", "ECHOSUM", "PROCESS", "STOP"], '\n') * '\n'
end

# a projection succeeded iff the .sum has at least one row led by a 4-digit calendar year
_has_data_row(text) = any(let s = strip(ln)
        length(s) >= 4 && (y = tryparse(Int, s[1:4]); y !== nothing && 1000 <= y <= 3000)
    end for ln in split(text, '\n'))

function _run(cn, variant)
    dir = mktempdir(); kf = joinpath(dir, "s.key"); write(kf, _keytext(cn))
    FVSjl.run_keyfile(kf; variant = variant)
end

@testset "CA/SO forkod forest-index-overflow crash" begin
    @testset "pure forkod remap (ca/forkod.f, so/forkod.f)" begin
        # CA: 518 Trinity → IFOR 5 (Shasta-Trinity 514); a normal code is unremapped
        p = _FK(Int32(518), Int32(0), Int32(0)); FVSjl.ca_forkod!(p)
        @test p.forest_idx == 5
        @test p.user_forest_code == 514
        p = _FK(Int32(511), Int32(0), Int32(0)); FVSjl.ca_forkod!(p)   # Plumas — no remap
        @test p.forest_idx == 4
        @test p.user_forest_code == 511

        # SO: 514 → IFOR 4 (Klamath 505); 702 → IFOR 8 (701); a normal code is unremapped
        p = _FK(Int32(514), Int32(0), Int32(0)); FVSjl.so_forkod!(p)
        @test p.forest_idx == 4
        @test p.user_forest_code == 505
        p = _FK(Int32(702), Int32(0), Int32(0)); FVSjl.so_forkod!(p)
        @test p.forest_idx == 8
        @test p.user_forest_code == 701
        p = _FK(Int32(601), Int32(0), Int32(0)); FVSjl.so_forkod!(p)   # Deschutes — no remap
        @test p.forest_idx == 1
        @test p.user_forest_code == 601
    end

    @testset "previously-crashing stands run end-to-end" begin
        @test isfile(_FK_DB)
        # CA 518 (Trinity) dense mega-stand — the original SIGSEGV in ca_dgcons!
        @test _has_data_row(_run("302060823489998", FVSjl.CentralCalifornia()))
        # SO 702 (Industry Lands) — SIGSEGV in so_dgcons! (IFOR 11)
        @test _has_data_row(_run("44907987020004", FVSjl.SouthCentralOregon()))
        # SO stand with whitebark-pine (sp 16) small trees — FieldError in small_tree_growth!
        @test _has_data_row(_run("722687117290487", FVSjl.SouthCentralOregon()))
    end
end
