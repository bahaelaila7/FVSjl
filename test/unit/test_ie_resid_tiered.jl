# test_ie_resid_tiered.jl — IE per-record fixes vs the LIVE FVSie_g16 tiered goldens (test/fixtures/tiered/ie: stands.db +
# <stand>_<regime>.key + the live .sum / DBS tables). Each testset names the Fortran it follows and the measured case.
module IEResidTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

function _case(cn, r)
    d = mktempdir()
    txt, db, crashed, _ = run_case("IE", cn, r; dir = d)
    return (txt = txt, db = db, crashed = crashed, ms = compare_case("IE", cn, r, txt, db))
end

# fvsvol.f:90-96 hands VOLINIT IREGN=KODFOR/100; mrules.f then gives region 6 (the Colville, KODFOR 621) COR='N' (raw
# Scribner) and OPT 23, region 1 COR='Y' and OPT 22. jl hard-wired region 1 for IE (MEASURED FVSie_g16 374547584489998
# 2015: RC 11.8" BdFt live 76 / jl 80, LP 12.6" MCuFt 36.9 / 36.6 — 32 inventory-year volume cells).
@testset "IE FW2 volumes follow the KODFOR region's NVEL merch rules (fvsvol.f/mrules.f) vs FVSie_g16" begin
    c = _case("374547584489998", "none")
    @test !c.crashed
    @test count(m -> m.col in ("MCuFt", "BdFt", "TCuFt"), c.ms) == 0
end

# estab.f:545-549 sets PNN(NCOUNT)=ESA on every plot of a calibrated fresh tally — the cycle-1 ingrowth tally too — and
# estab.f:583 floors each plot's PROB1 at PNN+0.0001, which the ingrowth NSTORE (:587-589) then divides by. jl left PNN=0
# on ingrowth (MEASURED FVSie_g16 3285544010690 2011: point 1 logistic 0.3833 < PNN 0.429568 ⇒ live PROB1 0.429668 and
# NSTORE 3, jl 0.3833 and NSTORE 4 ⇒ 2022 TPA live 1028 / jl 949, 20,187 tiered cells).
@testset "IE AUTOES ingrowth PROB1 floored at PNN=ESA (estab.f:545-589) vs FVSie_g16" begin
    for cn in ("3285544010690", "51032748020004")
        c = _case(cn, "none")
        @test !c.crashed
        @test count(m -> m.file in ("sum", "FVS_TreeList") || m.col in ("Tpa", "BA", "SDI", "CCF"), c.ms) == 0
    end
end
# dense.f:179-188 sums TPROB/TSUMD2 over IND1 species-major with WK5=D*(D*P) (dense.f is byte-identical in the IE build);
# RMSQD=SQRT(TSUMD2/TPROB) is FVS_Summary QMD/ATQMD. jl summed IE in record order with P·D² (MEASURED FVSie_g16
# 3285544010690 2012 QMD live 7.00555182, jl 7.00555038 — 1,280 QMD/ATQMD cells across the IE suite).
@testset "IE QMD (RMSQD) over IND1 with dense.f's WK5 (dense.f:179-188) vs FVSie_g16" begin
    for cn in ("3285544010690", "51032748020004", "374547584489998")
        @test count(m -> m.col in ("QMD", "ATQMD"), _case(cn, "none").ms) == 0
    end
end

# The thinbba stand 3027007010690 with a TREELIST added: FVSie_g16's FVS_TreeList (TPA, DBH, Ht, PctCr) 1996-2046.
function _thin_tl()
    fx = joinpath(@__DIR__, "..", "fixtures", "ie_thin")
    d = mktempdir()
    key = replace(read(joinpath(fx, "3027007010690_thinbba_tl.key"), String),
                  "\nout.db\n" => "\n" * joinpath(d, "out.db") * "\n",
                  "\nstands.db\n" => "\n" * joinpath(fixture_dir("IE"), "stands.db") * "\n")
    kp = joinpath(d, "x.key"); write(kp, key)
    FVSjl.run_keyfile(kp; variant = FVSjl.InlandEmpire())
    hdr, rows = db_table_rows(joinpath(d, "out.db"), "FVS_TreeList")
    ix = Dict(c => findfirst(==(c), hdr) for c in ("Year", "TreeId", "TreeIndex", "TPA", "DBH", "Ht", "PctCr"))
    got = Dict((r[ix["Year"]], strip(r[ix["TreeId"]]), r[ix["TreeIndex"]]) => r for r in rows)
    ghdr, grows = read_csv(joinpath(fx, "3027007010690_thinbba_tl.TreeList.csv"))
    return got, ix, grows
end
_nsame(a, b) = a == b || (let x = tryparse(Float64, a), y = tryparse(Float64, b); x !== nothing && x == y end)

# crown.f:479-480 resets OLDPCT to the current PCT when OLDPCT>PCT in a cycle that removed trees (ONTREM(7)>0) — a thin
# lowers every survivor's start-of-cycle percentile. jl omitted the thin branch (MEASURED FVSie_g16 3027007010690
# THINBBA 2006: LP PctCr at 2016 live 47/48/50, jl 44/42/44 on 25 records). Same line in em/kt/bc crown.f.
@testset "IE crown OLDPCT reset after a thin (crown.f:479-480) vs FVSie_g16" begin
    got, ix, grows = _thin_tl()
    @test length(grows) == 1785
    @test count(r -> r[1] == "2016" && !_nsame(r[7], got[(r[1], strip(r[2]), r[3])][ix["PctCr"]]), grows) == 0
end

# cutstk.f CLSSTK sums the class stocking as CSTOCK+TPA*(D*D*0.005454154) — factor first, then ×TPA; jl's (TPA·D²)·factor
# rounded REMOVE differently, so the last record THINBBA partially removes took a different PREM (MEASURED FVSie_g16
# 3027007010690 THINBBA 2006: record 19 cut 5.32047272, jl 5.32045555 ⇒ 2016 TPA, BAAA, the regen PROB1 and on).
# With it and the OLDPCT reset the whole thinned TreeList 1996-2046 (TPA, DBH, Ht, PctCr) equals live.
@testset "THINBBA class stocking at cutstk.f CLSSTK precision vs FVSie_g16 (whole TreeList)" begin
    got, ix, grows = _thin_tl()
    for (j, c) in ((4, "TPA"), (5, "DBH"), (6, "Ht"), (7, "PctCr"))
        @test count(r -> !haskey(got, (r[1], strip(r[2]), r[3])) || !_nsame(r[j], got[(r[1], strip(r[2]), r[3])][ix[c]]), grows) == 0
    end
end
# ---- FFE carbon: IE 11855985010690 SALVAGE (FMIN SALVAGE + CARBREPT) vs FVSie_g16 FVS_Carbon ----------------------------
const _SALV = Ref{Any}(nothing)
_salv() = (_SALV[] === nothing && (_SALV[] = _case("11855985010690", "salvage")); _SALV[].ms)
_crel(m) = (g = parse(Float64, m.gold); abs(parse(Float64, m.got) - g) / max(abs(g), 1e-12))
_ccol(col) = filter(m -> m.file == "FVS_Carbon" && m.col == col, _salv())

# fmcrbout.f:98-146 V(1)=BIOLIVE (fmdout.f TOTFOL+TOTLIV crowns+stems, FMPROB-weighted over FMMAIN's record list),
# V(2)=Σ FMPROB·FMSVL2(merch)·V2T, V(3)=Σ RBIO·FMPROB in REAL*4 at TItoTM/ACRtoHA=0.90718/0.4046945; jl rebuilt the
# report from separate Jenkins/pool helpers (MEASURED IE/EM FVS_Carbon aboveground rows off in the 6th-7th digit).
@testset "FFE carbon V(1)/V(2) at fmdout.f/fmcrbout.f precision vs FVSie_g16" begin
    @test isempty(_ccol("Aboveground_Total_Live"))
    @test isempty(_ccol("Aboveground_Merch_Live"))
end

end # module
