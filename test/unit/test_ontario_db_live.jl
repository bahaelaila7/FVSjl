# test_ontario_db_live.jl — ON (Ontario) DATABASE-input stand vs LIVE FVSon_g16.
#
# FVSDataHardwood.db stand LD3001 (94 metric trees, 1-ha fixed plot: BASAL_AREA_FACTOR −1, INV_PLOT_SIZE 1, SITE_INDEX
# 12.8 m, LATITUDE but no LONGITUDE, no heights) read through the ON METRIC dbsstandin.f and projected in 5-yr cycles
# (TIMEINT 5 ⇒ FINT=5 ≠ YR=10). Goldens (live FVSon_g16 runs of the fixture keys):
#   hardwood_live.rows        — .sum data rows of hardwood.key (the canonical FVS ON test: 10 cycles, ESTAB/STOCKADJ,
#                               THINBBA 2029)
#   hardwood_2cyc_tl_live.csv — FVS_TreeList_East_Metric of hardwood_2cyc.key (2 cycles), compared EXACTLY per record,
#                               including the three HISTORY-8 inventory-dead records listed in 2004 (TreeIndex 5998-6000).

using Test, FVSjl, SQLite, DBInterface
const _OND = FVSjl
const _OND_FX = joinpath(@__DIR__, "..", "fixtures", "ontario")

function _ond_run(key::AbstractString)
    dir = mktempdir()
    for f in (key, "FVSDataHardwood.db"); cp(joinpath(_OND_FX, f), joinpath(dir, f)); end
    txt = cd(() -> _OND.run_keyfile(key; variant = _OND.Ontario(), output = :sum), dir)
    return [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)], dir
end
function _ond_live_tl()
    tl = Dict{Tuple{Int,Int},NamedTuple}(); hdr = String[]
    for ln in eachline(joinpath(_OND_FX, "hardwood_2cyc_tl_live.csv"))
        startswith(ln, "#") && continue
        f = split(ln, ',')
        if f[1] == "Year"; hdr = String.(f); continue; end
        tl[(parse(Int, f[1]), parse(Int, f[2]))] = (; (Symbol(hdr[k]) => parse(Float64, f[k]) for k in 4:length(f))...)
    end
    return tl
end
function _ond_jl_tl()
    rows, dir = _ond_run("hardwood_2cyc.key")
    db = SQLite.DB(joinpath(dir, "hardwood_2cyc.db"))
    tl = Dict{Tuple{Int,Int},NamedTuple}()
    for r in DBInterface.execute(db, "SELECT Year,TreeIndex,TPH,MortPH,DBH,DG,Ht,HtG,PctCr,CrWidth,BAPctile,PtBAL,TCuM,MCuM,CCum " *
                                     "FROM FVS_TreeList_East_Metric")
        tl[(Int(r.Year), Int(r.TreeIndex))] = (; TPH = Float64(r.TPH), MortPH = Float64(r.MortPH), DBH = Float64(r.DBH),
            DG = Float64(r.DG), Ht = Float64(r.Ht), HtG = Float64(r.HtG), PctCr = Float64(r.PctCr),
            CrWidth = Float64(r.CrWidth), BAPctile = Float64(r.BAPctile), PtBAL = Float64(r.PtBAL),
            TCuM = Float64(r.TCuM), MCuM = Float64(r.MCuM), CCum = Float64(r.CCum))
    end
    return rows, tl
end
const _OND_LIVE = _ond_live_tl()
const _OND_JLROWS, _OND_JL = _ond_jl_tl()
_ond_cmp(col; sel = (k, v) -> true) =
    [((k, col, haskey(_OND_JL, k) ? _OND_JL[k][col] : NaN), (k, col, v[col]))
     for (k, v) in sort(collect(_OND_LIVE)) if sel(k, v)]

@testset "ON Hardwood DB: TreeList live-record set == live" begin
    @test sort([k for k in keys(_OND_JL) if k[2] < 2998]) == sort([k for k in keys(_OND_LIVE) if k[2] < 5998])
end

# (1) canada/on/dbsstandin.f:357-372 converts the METRIC sampling design: BASAL_AREA_FACTOR <0 ⇒ /HAtoACR (≥0 ⇒
# ·M2pHAtoFT2pACR), INV_PLOT_SIZE /HAtoACR, BRK_DBH ·CMtoIN; TREE_COUNT is read raw (dbstreesin.f:100) and notre.f
# expands it (PROB = count·(−BAF)/PI). jl pre-scaled TREE_COUNT by the exact 0.40468564 instead of FVS's 1/2.471 ⇒
# every record 9.999781 TPH vs live 9.999999.
@testset "ON Hardwood DB 2004: TPH per record == live (dbsstandin.f metric BAF / INV_PLOT_SIZE)" begin
    for (j, l) in _ond_cmp(:TPH; sel = (k, v) -> k[1] == 2004 && k[2] < 5998)
        @test j == l
    end
end

# (2) canada/on/dbsstandin.f:419-422: SITE_INDEX (metres) ·MtoFt. Read raw, the 12.8-"ft" site index aged every
# record wrongly in FINDAG (ABIRTH ⇒ Mowraski net-merch cull ⇒ CCum) and slowed the Penner growth. Cycle-0 .sum
# row (NMV 75 vs live 100) and every record's CCum.
@testset "ON Hardwood DB 2004: CCum per record + cycle-0 .sum row == live (dbsstandin.f SITE_INDEX·MtoFt)" begin
    for (j, l) in _ond_cmp(:CCum; sel = (k, v) -> k[1] == 2004 && k[2] < 5998)
        @test j == l
    end
    @test first(_OND_JLROWS) == split(first(readlines(joinpath(_OND_FX, "hardwood_live.rows"))))
end

# (3) canada/on/forkod.f (called by dbsstandin.f BEFORE the DB LATITUDE/LONGITUDE overrides): KODFOR 915/916 sets
# TLAT=46.78, TLONG=92.11, ELEV=16 where still 0. LD3001 carries LATITUDE but no LONGITUDE ⇒ live TLONG=92.11, which
# enters the Hopkins index of the open-grown crown width (CCF). jl left TLONG=0 ⇒ every CrWidth and the CCF column off.
@testset "ON Hardwood DB 2004: CrWidth per record == live (forkod.f KODFOR 915 TLONG default)" begin
    for (j, l) in _ond_cmp(:CrWidth; sel = (k, v) -> k[1] == 2004 && k[2] < 5998)
        @test j == l
    end
end

# (4) canada/on/gradd.f:79-90 rescales DG from the YR=10 basis to FINT only AFTER GRINCR: htgf.f (DBH10=DBH+DG/BARK),
# regent.f (DGGR blend, small-tree DDS·YR/FNT) and morts.f (G=(DG/BARK)·(FINT/10)) all read the 10-yr DG. jl scaled DG
# to FINT inside the DG driver ⇒ with TIMEINT 5 the large-tree HTG and the mortality were ~halved and the small-tree
# DG doubled (2009 record 1: HtG 0.312 vs live 0.608 m, MortPH 0.124 vs 0.277).
@testset "ON Hardwood DB 2009-2014: per record == live (gradd.f DG FINT scaling after GRINCR)" begin
    for col in (:TPH, :MortPH, :DBH, :DG, :Ht, :HtG, :PctCr, :CrWidth, :TCuM, :MCuM, :CCum),
        (j, l) in _ond_cmp(col; sel = (k, v) -> k[1] >= 2009)
        @test j == l
    end
end

# The canonical ON test key: 10 five-year cycles with ESTAB/STOCKADJ and a THINBBA in 2029 — every .sum row == live.
# ESTAB also needs ON's establishment XMIN/HHTMAX (canada/on/esblkd.f): the shared tables had no ON columns ⇒
# KeyError :estab_min_ht, the run crashed at the first ESTAB cycle.
@testset "ON Hardwood DB canonical key: every .sum row == live" begin
    rows, _ = _ond_run("hardwood.key")
    @test length(rows) == 11
    @test rows == [split(l) for l in readlines(joinpath(_OND_FX, "hardwood_live.rows"))]
end

# (6) The cycle-0 inventory-dead records (HISTORY 8 ⇒ IMC 9), listed in 2004 at TreeIndex MAXTRE+1−k (ON MAXTRE=6000):
#  • cratet.f DO 145 dubs a missing dead height with the Penner HTONT(ISPC,D,RMSQD,BA,H) only (the live loop's HTDBH
#    override is absent), RMSQD/BA from the CRATET backdating DENSE (RMSQD over live + FINT/FINTM-inflated dead) and
#    htont.f's SIM = SITEAR(ISISP) (the SITE species');
#  • crown.f DO 79 dubs the dead crown (TWIGS CR at the CRATET BA);
#  • vols.f IPASS=2 volumes them, and OBFVOL's `IMC.GT.3` sends the 0.0001 sentinel to MOWRASKI (CCum);
#  • BAPctile = the CRATET DENSE PCTILE over the IND1-seeded real-DBH IND (dead WK5 = 0), PtBAL = 0.
# jl dubbed them with the live Curtis-Arney HTDBH, left the crown/volume 0 and filed them at TreeIndex 2998-3000.
@testset "ON Hardwood DB 2004 inventory-dead records == live (cratet.f DO 145 HTONT, crown.f DO 79, vols IPASS=2)" begin
    @test sort([k for k in keys(_OND_JL) if k[1] == 2004 && k[2] >= 2998]) == [(2004, 5998), (2004, 5999), (2004, 6000)]
    for col in (:Ht, :PctCr, :CrWidth, :BAPctile, :PtBAL, :TCuM, :MCuM, :CCum, :MortPH),
        (j, l) in _ond_cmp(col; sel = (k, v) -> k[1] == 2004 && k[2] >= 5998)
        @test j == l
    end
end
