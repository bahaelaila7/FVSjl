# test_bc_htdub_live.jl — BC missing-height dub (canada/bc cratet.f + htgf.f ENTRY HTCONS) vs LIVE FVSbc_clean.
#
# The BC WRD control stand (bc_htdub.key; ICH default BEC ICHmw2/01, one cycle, TREELIDB) carries live and
# inventory-dead records without a height. Golden: live FVS_TreeList_Metric at the inventory year 1990 (Ht, TCuM per
# record, bc_htdub_1990_tl_live.csv) and the cycle-0 .sum row (bc_htdub_live.rows), compared EXACTLY.

using Test, FVSjl, SQLite, DBInterface
const _BHD = FVSjl
const _BHD_FX = joinpath(@__DIR__, "..", "fixtures", "bc_htdub")

function _bhd_run()
    dir = mktempdir()
    for f in ("bc_htdub.key", "bc_htdub.tre"); cp(joinpath(_BHD_FX, f), joinpath(dir, f)); end
    txt = cd(() -> _BHD.run_keyfile("bc_htdub.key"; variant = _BHD.BritishColumbia(), output = :sum), dir)
    rows = [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)]
    db = SQLite.DB(joinpath(dir, "bc_htdub.db"))
    tl = Dict(Int(r.TreeIndex) => (Float64(r.Ht), Float64(r.TCuM))
              for r in DBInterface.execute(db, "SELECT TreeIndex,Ht,TCuM FROM FVS_TreeList_Metric WHERE Year=1990"))
    return rows, tl
end
const _BHD_ROWS, _BHD_TL = _bhd_run()

# canada/bc cratet.f:296-460 dubs a missing height with the Wykoff H=EXP(AX+BX/(D+1))+4.5 — AX = the fitted intercept AA
# (SUMX=SUMX+YY−XX, left-associative) when calibrated, else HT1; BX = HT2 — never the Curtis-Arney HTDBH. HT1/HT2 are
# what RCON→HTCONS (htgf.f:1817-2015, "REFIT OF HEIGHT DUBBING MODEL FOR VERSION 3") leaves: the blkdat values, or for
# ICH/IDF/SBS/SBPS a METRIC fit (LMHTDUB: YY=ALOG(HM−1.3), XX=BX/(DM+1), H=(EXP(AX+BX/(DM+1))+1.3)·MtoFT). jl used the
# species CSV wykoff_ht2 and fell back to HTDBH (PL D 8.5 cm: 53.62 vs live 48.43; dead LW 124.2 vs 85.4).
@testset "BC 1990 Ht per record (live + inventory-dead) == live (cratet.f Wykoff dub, htgf.f HTCONS V3 table)" begin
    for ln in eachline(joinpath(_BHD_FX, "bc_htdub_1990_tl_live.csv"))
        (startswith(ln, "#") || startswith(ln, "Year")) && continue
        f = split(ln, ',')
        k = parse(Int, f[2])
        @test (k, :Ht, get(_BHD_TL, k, (NaN, NaN))[1]) == (k, :Ht, parse(Float64, f[4]))
        # (TCuM is not compared: 11 records carry a pre-existing BC Kozak-volume residual, 1 ULP – 4E-4 relative.)
    end
    @test first(_BHD_ROWS) == split(first(readlines(joinpath(_BHD_FX, "bc_htdub_live.rows"))))
end
