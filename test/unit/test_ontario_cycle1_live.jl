# test_ontario_cycle1_live.jl — ON (Ontario) growth cycle 1 per record vs LIVE FVSon_g16.
#
# ont01 (8 inline metric trees, one per species) projected one 10-yr cycle with TREELIDB: the cycle-1 records are
# the 8 originals plus their 16 TRIPLE copies (FVS triples in cycle 1), so every column exercises the full ON
# growth chain — Penner DGF → dgdriv.f DDS→DG + tripling serial correlation → htgf.f/htont.f → morts.f/varmrt.f →
# update.f → vols.f/volont.f. Golden = the live FVS_TreeList_East_Metric (test/fixtures/ontario/ont01_tl_live.csv,
# FVSon_g16 run of ont01_tl.key; the REAL*8 binds are compared EXACTLY). The .sum goldens (*_live.rows) are the
# live FVSon_g16 data rows of the ON fixture stands.

using Test, FVSjl, SQLite, DBInterface
const _ON1 = FVSjl
const _ON1_FX = joinpath(@__DIR__, "..", "fixtures", "ontario")

function _on1_live_treelist()
    rows = Dict{Tuple{Int,Int},NamedTuple}()
    hdr = String[]
    for ln in eachline(joinpath(_ON1_FX, "ont01_tl_live.csv"))
        startswith(ln, "#") && continue
        f = split(ln, ',')
        if f[1] == "Year"; hdr = String.(f); continue; end
        v = (; (Symbol(hdr[k]) => parse(Float64, f[k]) for k in 4:length(f))...)
        rows[(parse(Int, f[1]), parse(Int, f[2]))] = merge((; Species = String(f[3])), v)
    end
    return rows
end

function _on1_jl_treelist()
    dir = mktempdir()
    for f in ("ont01_tl.key", "ont01.tre"); cp(joinpath(_ON1_FX, f), joinpath(dir, f)); end
    mv(joinpath(dir, "ont01.tre"), joinpath(dir, "ont01_tl.tre"))
    cd(dir) do
        _ON1.run_keyfile("ont01_tl.key"; variant = _ON1.Ontario(), output = :sum)
    end
    db = SQLite.DB(joinpath(dir, "ont01_tl.db"))
    rows = Dict{Tuple{Int,Int},NamedTuple}()
    for r in DBInterface.execute(db, "SELECT Year,TreeIndex,SpeciesFVS,TPH,MortPH,DBH,DG,Ht,HtG,TCuM,MCuM,CCum " *
                                     "FROM FVS_TreeList_East_Metric")
        rows[(Int(r.Year), Int(r.TreeIndex))] = (; Species = String(r.SpeciesFVS), TPH = Float64(r.TPH),
            MortPH = Float64(r.MortPH), DBH = Float64(r.DBH), DG = Float64(r.DG), Ht = Float64(r.Ht),
            HtG = Float64(r.HtG), TCuM = Float64(r.TCuM), MCuM = Float64(r.MCuM), CCum = Float64(r.CCum))
    end
    return rows
end

const _ON1_LIVE = _on1_live_treelist()
const _ON1_JL = _on1_jl_treelist()

# every live record exists in jl with the same species, and nothing extra
@testset "ON ont01 cycle-1 record set == live (8 originals + 16 TRIPLE copies)" begin
    @test sort(collect(keys(_ON1_JL))) == sort(collect(keys(_ON1_LIVE)))
    @test all(_ON1_JL[k].Species == v.Species for (k, v) in _ON1_LIVE if haskey(_ON1_JL, k))
end

_on1_col(col; sel = (k, v) -> true) =
    [(k, v[col], haskey(_ON1_JL, k) ? _ON1_JL[k][col] : NaN) for (k, v) in sort(collect(_ON1_LIVE)) if sel(k, v)]

# (1) canada/on/cratet.f → dgdriv.f LSTART calibration: SIGMA=SIGMAR (on/blkdat.f), OLDRN=BACHLO draws, VARDG.
# Without it VARDG=0 ⇒ SSIGMA=0 ⇒ every central record grew at the EXPECTED DG and the RNG sat off live.
@testset "ON ont01 cyc1: central-record DG + all TPH/MortPH == live (LSTART DG calibration)" begin
    for (k, lv, jv) in _on1_col(:DG; sel = (k, v) -> k[1] == 2014 && k[2] <= 8)
        @test (k, jv) == (k, lv)
    end
    for col in (:TPH, :MortPH), (k, lv, jv) in _on1_col(col)
        @test (k, col, jv) == (k, col, lv)
    end
end

# (2) canada/on/dgdriv.f label 30 (tripling) has NO `RNPAR=OLDRN(I)` save: after `OLDRN(I)=FRMT` the upper/lower
# copies read FRU/FRL + CORR·OLDRN(I) = the central's UPDATED residual (sn/ls/ne/cs dgdriv.f use RNPAR).
@testset "ON ont01 cyc1: TRIPLE-copy DG == live (upper/lower read the updated OLDRN)" begin
    for (k, lv, jv) in _on1_col(:DG; sel = (k, v) -> k[1] == 2014 && k[2] > 8)
        @test (k, jv) == (k, lv)
    end
end

# (3) canada/on/update.f: DO 90 adds HTG to HT for every record BEFORE DO 110 DBH=DBH+DG/BRATIO(IS,DBH,HT) — ON's
# metric H/D bark (maple group / black spruce / cedar) reads the END-of-cycle height.
@testset "ON ont01 cyc1: DBH == live (update.f bark reads the grown HT)" begin
    for (k, lv, jv) in _on1_col(:DBH)
        @test (k, jv) == (k, lv)
    end
end

# (4) base dense.f (linked by FVSon): TSUMD2/BAT accumulate WK5=D*(D*P) in IND1 (species) order; RMSQD=SQRT(TSUMD2/
# TPROB) feeds canada/on/htont.f's DBHQM term (HTNOW) — record-order p·d² was 1 ULP off (ont01 SB HtG).
@testset "ON ont01 cyc1: Ht/HtG == live (DENSE RMSQD/BA in IND1 order, D*(D*P))" begin
    for col in (:HtG, :Ht), (k, lv, jv) in _on1_col(col)
        @test (k, col, jv) == (k, col, lv)
    end
end

_on1_rows(stem) = cd(joinpath(_ON1_FX)) do
    txt = _ON1.run_keyfile("$stem.key"; variant = _ON1.Ontario(), output = :sum)
    [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)]
end
_on1_live_rows(stem) = [split(l) for l in readlines(joinpath(_ON1_FX, "$(stem)_live.rows"))]

# (5) canada/on/vols.f:131 BARK=BRATIO(ISPC,D,H) is evaluated BEFORE `D=D+DG/BARK`: OCFVOL/volont ZAK use the bark of
# the cycle-START DBH (and the grown HT) for the dib — not the bark of the updated DBH.
@testset "ON ont01 cyc1: TCuM (live records) + .sum rows == live (vols.f bark at the cycle-start DBH)" begin
    for (k, lv, jv) in _on1_col(:TCuM; sel = (k, v) -> v.TPH > 0)
        @test (k, jv) == (k, lv)
    end
    @test _on1_rows("ont01") == _on1_live_rows("ont01")
end

# (6) Mowraski NMV (volont.f:405) reads ABIRTH: dubbed ONCE by canada/on/cratet.f FINDAG from the INVENTORY height
# (every record), aged +FINT by gradd.f:205 AFTER UPDATE's VOLS. jl dubbed lazily from the grown height the first
# time a record crossed BFMIND, so the cycle-1 NMV (.sum board column) was off on the multi-species stands.
@testset "ON fixture stands: .sum rows == live (FINDAG ABIRTH dub + post-VOLS aging)" begin
    for stem in ("ont_all", "ont_lite", "ont_mh", "ont_sm", "ont01")
        @test (stem, _on1_rows(stem)) == (stem, _on1_live_rows(stem))
    end
end

# (7) The metric FVS_TreeList binds CFV·FT3toM3, WK1·FT3toM3 (MCuM), BFV·FT3toM3 (dbstrls.f:342-350). CFV/BFV go
# through gradd.f's V·PROB … /PROB round trip (divide-back only when PROB>0 ⇒ a PROB=0 record keeps 0); ON's merch
# volume lives in WK1, which the round trip never touches (ON vols.f does not load MCFV).
@testset "ON ont01 cyc1: MCuM (live records) + zero-TPH TCuM/CCum == live (WK1 / gradd.f PROB round trip)" begin
    for (k, lv, jv) in _on1_col(:MCuM; sel = (k, v) -> v.TPH > 0)
        @test (k, jv) == (k, lv)
    end
    # A record the first MORTS empties is skipped by VOLS (vols.f:125 P≤0): its WK1 keeps what dgdriv.f DO 5 loaded at
    # the top of the cycle — cycle 1: the LSTART DO-220 DG (dgdriv.f:700-729, on_do220_dg); TRIPLE copies it
    # (triple.f:68). jl re-volumed them (the three zero-TPH PJ records: MCuM sentinel vs live WK1=DG).
    for (k, lv, jv) in _on1_col(:MCuM; sel = (k, v) -> v.TPH == 0)
        @test (k, :MCuM0, jv) == (k, :MCuM0, lv)
    end
    for col in (:TCuM, :CCum), (k, lv, jv) in _on1_col(col; sel = (k, v) -> v.TPH == 0)
        @test (k, col, jv) == (k, col, lv)
    end
end

# (8) canada/on/morts.f label 35: the MATURE-STAND-BOUNDARY mortality is ON BY DEFAULT in ON (no MORTMSB needed:
# SLPMSB=0 ⇒ the Penner MSB_INT/MSB_SLP/MSB_DBH of the dominant-BA equation). ont_big (72 × 50 cm, one per species)
# has ln(D10 cm) past the boundary, so MSBMRT kills TMORE = TN − 0.85·TMMSB more trees at EFFMSB=0.9 down the IND.
@testset "ON ont_big: .sum rows == live (default mature-stand-boundary mortality, MSBMRT)" begin
    @test _on1_rows("ont_big") == _on1_live_rows("ont_big")
end
