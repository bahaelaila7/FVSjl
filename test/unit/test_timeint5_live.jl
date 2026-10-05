# test_timeint5_live.jl — 5-year cycles in YR=10 variants vs LIVE oracles (gradd.f:79-90 DG rescale after GRINCR).
#
# The FVS tests/FVS<v> keys (pnt01, ect01, wct01, bmt01, cit01) re-run with `TIMEINT 5` before every PROCESS, so
# FINT=5 ≠ YR=10. Goldens: the live FVS<v>_g16 .sum data rows of the same keys (<stem>_t5_live.rows).
#
# gradd.f:79-90 rescales DG from the YR-year basis to FINT years only AFTER GRINCR: DGDRIV bounds it, HTGF reads
# DBH+DG/BARK (the 10-year diameter), REGENT blends it, MORTS forms G=(DG/BARK)·(FINT/YR) and TRIPLE copies it, all on
# the 10-year basis. jl scaled DG to FINT inside the DG driver, so under 5-year cycles the height growth and the
# diameter-driven mortality ran on a half-size DG (PN pnt01: 53 of 56 rows off; now 0).

using Test, FVSjl
const _T5 = FVSjl
const _T5_FX = joinpath(@__DIR__, "..", "fixtures", "timeint5")

function _t5_rows(stem, variant)
    dir = mktempdir()
    for ext in ("key", "tre"); cp(joinpath(_T5_FX, "$stem.$ext"), joinpath(dir, "$stem.$ext")); end
    txt = cd(() -> _T5.run_keyfile("$stem.key"; variant = variant, output = :sum), dir)
    return [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)]
end
_t5_live(stem) = [split(l) for l in readlines(joinpath(_T5_FX, "$(stem)_live.rows"))]
# stand index of each row (a new stand starts when the year goes back)
_t5_stand(rows) = (k = 1; [i == 1 ? 1 : (parse(Int, rows[i][1]) < parse(Int, rows[i-1][1]) ? (k += 1) : k) for i in eachindex(rows)])

# Known residuals, named (all outside the DG rescale):
#  • bmt01/cit01 stand 5 (bare ground, PLANT 1992): jl's planted cohort runs one 5-year cycle ahead of live (regent.f LESTB
#    FNT−5 / LSKIPH when FINT≤5 — the establishment-cycle height growth under 5-year cycles is not ported);
#  • one-unit print knife-edges in a single volume column: ect01 row 2018 TCuFt, wct01 rows 2010 MCuFt / 2018 BdFt,
#    cit01 rows 2008/2018 BdFt (its 2060 BdFt edge became exact earlier on this branch).
const _T5_CASES = [("pnt01_t5", _T5.PacificNorthwest(), Int[], Int[]),
                   ("ect01_t5", _T5.EastCascades(),    Int[], Int[]),
                   ("wct01_t5", _T5.WestCascades(),    Int[], Int[]),
                   ("bmt01_t5", _T5.BlueMountains(),   Int[5], Int[]),
                   ("cit01_t5", _T5.CentralIdaho(),    Int[5], Int[])]

@testset "TIMEINT 5, YR=10 variants: .sum rows == live (gradd.f DG rescale after GRINCR)" begin
    for (stem, v, skip_stands, _) in _T5_CASES
        jl = _t5_rows(stem, v); lv = _t5_live(stem)
        @test length(jl) == length(lv)
        st = _t5_stand(lv)
        nfull = 0
        for i in eachindex(lv)
            st[i] in skip_stands && continue
            # Year..QMD (TPA, BA, SDI, CCF, TopHt, QMD) exact in every row
            @test (stem, i, jl[i][1:8]) == (stem, i, lv[i][1:8])
            nfull += (jl[i] == lv[i])
        end
        # every full row exact except the named one-unit volume knife-edges
        nknife = Dict("pnt01_t5" => 0, "ect01_t5" => 1, "wct01_t5" => 2, "bmt01_t5" => 0, "cit01_t5" => 0)[stem]
        @test (stem, nfull) == (stem, count(i -> !(st[i] in skip_stands), eachindex(lv)) - nknife)
    end
end

# KT, CR, WS, CA, SO under TIMEINT 5 (kt/cr/ws/ca/so tests keys; live FVSkt_clean / FVScr_clean / FVS<v>_g16 rows).
# Their DG stays on the YR=10 basis through GRINCR like the variants above, REGENT's SCALE2 is YR/FNT, and KT's MORTS
# reads cycle-1 WK1 from the DO-220 calibration DG with OLDFNT = the previous period. Rows still off are named
# elsewhere (the FFE TEST stands of wst01/cat01, the residuals the 10-year keys also carry, one-unit print edges), so
# this is a ratchet on the exact-row counts measured when the rescale landed (before, full rows / Year..QMD rows:
# KT 1/6, CR 1/28, WS 3/9, CA 32/51, SO 11/27).
const _T5_RATCHET = [("ktt01_t5", _T5.Kootenai(),           29, 46),
                     ("crt01_t5", _T5.CentralRockies(),     29, 35),
                     ("wst01_t5", _T5.WestSierra(),         52, 53),
                     ("cat01_t5", _T5.CentralCalifornia(),  50, 53),
                     ("sot01_t5", _T5.SouthCentralOregon(), 44, 47)]

@testset "TIMEINT 5, KT/CR/WS/CA/SO: exact .sum rows vs live (ratchet)" begin
    for (stem, v, nfull, nhead) in _T5_RATCHET
        jl = _t5_rows(stem, v); lv = _t5_live(stem)
        @test length(jl) == length(lv)
        @test (stem, count(i -> jl[i] == lv[i], eachindex(lv)) >= nfull) == (stem, true)
        @test (stem, count(i -> jl[i][1:8] == lv[i][1:8], eachindex(lv)) >= nhead) == (stem, true)
    end
end

# tt/htgf.f:309-315 + :728: the PP case (CASE 10, the CI equation) falls into the common label-201 tail like every other
# case, HTG·SCALE·XHMULT·EXP(HTCON) with SCALE=FINT/YR. jl returned the 10-year HTG for PP, so under 5-year cycles a PP
# stand doubled its height growth (S248112 PP, TIMEINT 5: 1995 TopHt 67 live; live FVStt_g16 rows).
@testset "TIMEINT 5, TT ponderosa HTG scaled by FINT/YR (tt/htgf.f:728)" begin
    jl = _t5_rows("ttpp_t5", _T5.Teton()); lv = _t5_live("ttpp_t5")
    @test length(jl) == length(lv)
    for i in eachindex(lv); @test (i, jl[i][7]) == (i, lv[i][7]); end     # TopHt
end

# /CONTRL/ YR is the growth models' base period (blkdat.f DATA YR/10.0/), not the TIMEINT cycle length jl kept in
# control.year: em/regent.f:218 SCALE=YR/FINT, ie/regent.f SCALE2=YR/NTYR, ak/dgf.f:521 DGPRED=YR·BASEDG·PFMOD,
# ak/regent.f:200 SCALE2=YR/FNT. emt01 with TIMEINT 5: 35 → 10 rows off live FVSem_g16 (the same 10 the 10-year key has).
@testset "TIMEINT 5, EM REGENT reads YR=10 (em/regent.f:218)" begin
    jl = _t5_rows("emt01_t5", _T5.EasternMontana()); lv = _t5_live("emt01_t5")
    @test length(jl) == length(lv)
    @test count(i -> jl[i] == lv[i], eachindex(lv)) >= 46
end

# AK keeps DG on the YR=10 basis through GRINCR and GRADD rescales it to FINT (gradd.f:79-90), like the other YR=10
# variants, once its DGF/REGENT read YR: akt01 with TIMEINT 5, every row equals live FVSak_g16 (54 of 56 off before).
@testset "TIMEINT 5, AK .sum rows == live (gradd.f:79-90, ak/dgf.f:521)" begin
    jl = _t5_rows("akt01_t5", _T5.SoutheastAlaska()); lv = _t5_live("akt01_t5")
    @test length(jl) == length(lv)
    for i in eachindex(lv); @test (i, jl[i]) == (i, lv[i]); end
end

# tt/morts.f:228,701,812,836 (and ut/morts.f:214,571,682,706) G=(DG(I)/BARK)*(FINT/10.0): MORTS reads the 10-year DG
# (GRADD rescales after GRINCR), so the Zeide D10, the D10N re-pass and the BAMAX check scale it to the cycle. jl used the
# 10-year G under 5-year cycles ⇒ over-thinning (S248112 PP TIMEINT 5: 1995 TPA 516 jl / 521 live).
@testset "TIMEINT 5, TT MORTS G on FINT/10 (tt/morts.f:228) — every row == live" begin
    jl = _t5_rows("ttpp_t5", _T5.Teton()); lv = _t5_live("ttpp_t5")
    @test length(jl) == length(lv)
    for i in eachindex(lv); @test (i, jl[i]) == (i, lv[i]); end
end
