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

# fmcrbout.f:273 BIOROOT·(1−CRDCAY)**NYRS and fmsadd.f:317 (1−CRDCAY)**10 are REAL**INTEGER ⇒ libgcc __powisf2 (fpowi),
# not Julia's Float64-rounded `^` (MEASURED gfortran-16 -O0 `(1.0-0.0425)**10` = Z'3F25D106', Julia 3F25D105).
@testset "Dead-root decay (1-CRDCAY)**N as __powisf2 (fmcrbout.f:273, fmsadd.f:317) vs gfortran" begin
    @test FVSjl.fpowi(1f0 - FVSjl._FM_CRDCAY, 10) === reinterpret(Float32, 0x3F25D106)
    @test all(m -> _crel(m) < 1e-6, _ccol("Belowground_Dead"))
end

# cratet.f:482-488 hands FMSSEE the height HS=ITRUNC·.01 of a top-killed input dead tree, so the class MINHT/MAXHT that
# splits FMSADD's two height classes (fmsadd.f:119-123) sees the broken height, while HTCL still compares HT(I). jl used
# HT for both ⇒ LP 9.5"/75' (top 61') binned low (MEASURED FVSie_g16 2006 snag records: DBHS 9.20/HTDEAD 84.25 ×24/ac,
# jl 9.10/87.33 ×18 ⇒ inventory Standing_Dead 8.9956 vs 8.9772).
@testset "Inventory snag height classes use FMSSEE's broken-top height (cratet.f:482-488) vs FVSie_g16" begin
    @test !any(m -> m.col == "Standing_Dead" && m.year == "2006", _ccol("Standing_Dead"))
end

# fmcwd.f:311-403 (CWD1/CWD2/CWD3): each size class books MAX(0,P2−P1)·TVOLI of the fallen stem, P(h)=(1−h/HTD)³ between
# breakpoints, from LOHT (hard 0.10, soft 1.0) up — NOT renormalized to the whole stem, nothing for HTD≤4.5 ft, and only
# pieces with DIF·density > 1E-6. jl renormalized (+1/P(0.10) ≈ 0.3% per fallen snag) and dumped short snags whole into
# the DBH class (MEASURED FVSie_g16 2006 FMSNAG: hard 3-6" +0.182612 vs jl +0.183183; 2016 hard <0.25" 0.0907 vs 0.0967).
@testset "Snag-fall CWD split at fmcwd.f's un-normalized cone volumes vs FVSie_g16" begin
    @test all(m -> _crel(m) < 1e-6, _ccol("Forest_Down_Dead_Wood"))
    @test all(m -> _crel(m) < 1e-6, _ccol("Forest_Floor"))
end

# ie/fmvinit.f:500-507 TFALL(I,0)=MIN(2,LEAFLF), (I,1)=5, (I,2)=MIN(5,TFALL(I,3)), (I,3..5) per species (10/15/20): FMSCRO
# spreads a dead crown's CWD2B debris over those years. jl used the SN class rows (clamped to 1/1/2/4 yrs) ⇒ LP foliage
# fell in one year (MEASURED 2016 CWD2B(4,0,1..2) 118.08+118.08 live vs 236.16+0 ⇒ 2026 Standing_Dead −0.1%).
@testset "Dead-crown fall years from ie/fmvinit.f TFALL vs FVSie_g16" begin
    @test all(m -> _crel(m) < 1e-6, _ccol("Standing_Dead"))
end

# ie/fmcba.f reads CRWDTH(I) from ie/cwidth.f → cwcalc.f (IEMAP), whose small-tree forms differ from ccfcal MODE=2
# (MEASURED 2026 FMCBA: AF/ES/GF seedlings 0.5/0.55 ft vs MODE=2 1.2/1.09 ⇒ TOTCRA 32619 vs 32643 ⇒ PERCOV ⇒ FLIVE).
@testset "FMCBA canopy cover from cwcalc.f crown widths (ie/fmcba.f:244) vs FVSie_g16" begin
    @test all(m -> _crel(m) < 1e-6, _ccol("Forest_Shrub_Herb"))
end

# fmcrbout.f V(3) is taken in FMMAIN (gradd.f:118), after REGENT set seedling DBH(K)=0.1+DIAM·.01+HK·.001
# (ie/regent.f:882) — jl's pre-growth sample read 0.1 (MEASURED 2006: 0.1060 vs 0.1 ⇒ Belowground_Live 7.19485 vs 7.19215).
@testset "Belowground_Live on FMMAIN's post-REGENT DBH (gradd.f:118) vs FVSie_g16" begin
    @test isempty(_ccol("Belowground_Live"))
end

# ---- FFE fire kill: IE 4769882010690 SIMFIRE (fire 2014) vs FVSie_g16 --------------------------------------------------
# fmeff.f:352-527 pools the burned record's crowns into CWD2B2 in three FMSCRO calls — the crown-fire share CRBURN (no
# foliage, half the 0-0.25" + its OLDCRW), the killed rest ((1−CRBURN)·PMORT, foliage and CRW1BN burned over the scorched
# PROPCR) and the survivors' scorched-dead crown ((1−CRBURN)(1−PMORT), no foliage, PROPCR of each woody size) — and
# FMSADD(IYR,1) (fmeff.f:608) bins the killed trees into class-mean snag records. jl booked every killed tree's crown with
# the scorched-kill form and kept one snag per tree (MEASURED 2014: CWD2B2 sizes 0-3 1899/3399/7488/2582 live vs
# 3248/3748/7427/2553; 67 snag records vs 617 ⇒ Standing_Dead 30.519 vs 31.197).
@testset "Fire-killed crowns (fmeff.f:352-527) + binned fire snags (fmeff.f:608) vs FVSie_g16" begin
    ms = _case("4769882010690", "simfire").ms
    @test !any(m -> m.file == "FVS_Carbon" && m.col == "Standing_Dead" && m.year in ("2004", "2014"), ms)
end

# FMSADD (fmsadd.f, identical in all 24 builds) bins EVERY snag source — inventory ITYP=3, cut ITYP=2, fire/pile ITYP=1,
# mortality ITYP=4, SNAGINIT — into class-mean records in species-major slot order with emptied-record reuse. jl binned
# only the R6 variants' sources (non-R6 inventory snags stayed one per tree, mortality records in first-seen order).
# MEASURED FVSsn_g16 205045340010854 SALVAGE Standing_Dead 2023: live 0.259012, jl 0.259050.
@testset "Snag records via FMSADD binning + slot order in every variant vs FVSsn_g16" begin
    d = mktempdir()
    txt, db, crashed, _ = run_case("SN", "205045340010854", "salvage"; dir = d)
    ms = compare_case("SN", "205045340010854", "salvage", txt, db)
    @test !crashed
    @test !any(m -> m.file == "FVS_Carbon" && m.col == "Standing_Dead" && m.year in ("2018", "2023"), ms)
end

end # module
