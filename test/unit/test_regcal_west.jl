# test_regcal_west.jl — the LSTART small-tree HEIGHT calibration (REGCAL, {cr,ut,bm,nc}/regent.f label 40) vs
# FVScr_clean / FVSut_g16 / FVSbm_g16 / FVSnc_g16. Sibling of test_regcal_em_ie.jl (same fixture shape).
#
# Fixture (one per variant): the shipped PN inventory (backdated by its measured DG, 2 dead records), base species
# remapped to the variant's codes, plus six seedlings/saplings of 12-15 species carrying a measured HTG, one per DBH
# class 0.6-2.0". Every REGCAL arm is exercised: CR POTHTG (IVFLAG vigor species 9,12,16,23,30) + Sheppard aspen/PB
# (20,28); UT SJ/5 conifers, the woodland/CR/SO-surrogate SJ·(1.5SJ−H) form, aspen; BM SMHTGF, LM SI/5, aspen·2.4,
# WJ vigor; NC HTGR5 on RELHT/BA. The UT seedlings carry PROB 1 so PCTRED (0.2730) keeps most CORNEW inside the trap.
# Golden = the live per-species "SUMS FOR SPECIES n: SNP SNX SNY" (DEBUG REGENT, 2026-09-26) ⇒ CORNEW = SNY/SNX;
# trapped (outside [0.0821,12.1825]) and N<NCALHT species keep HCOR 0.
#
# What had to be right: CR PCTRED from the cratet.f:175 backdating DENSE's AVH and RELDEN (no AVHT40 before REGENT;
# jl recomputed both ad hoc ⇒ every non-aspen CORNEW 4.4% high); UT X = AVHT40·RELDEN(backdated); UT's cratet.f
# 50-yr site conversion run at CRATET time on NOTRE-expanded PROB (TEMCCF 469 ⇒ LP SITEAR 17.09; jl ran it before
# NOTRE ⇒ TEMCCF at the 125 floor ⇒ 30.07); NC XBA = the backdating DENSE's BA and RELHT on the AVHT40 AVH (jl's
# ad-hoc BA/AVH put every SNX ~2.7% off). BM was already exact; it also runs 3 cycles bit-exact on the .sum.
using Test
using FVSjl

const _RCW_FX = joinpath(@__DIR__, "..", "fixtures")

function _rcw_hcor(v, variant)
    key = joinpath(_RCW_FX, v, "regcal", "regcal.key")
    s, _ = FVSjl.initialize(key; variant = variant)
    FVSjl.notre!(s); FVSjl.setup_growth!(s)
    return s.calib.htg_cor_init
end

# species => (SNX, SNY) from the live sums; `nothing` ⇒ HCOR stays 0 (N<5 or CORNEW trapped)
const _RCW_CR = Dict(1 => (1063.52, 1332.00), 3 => (1362.81, 2166.00), 5 => (1063.52, 1332.00), 9 => (589.03, 1332.00),
                     10 => nothing, 11 => (921.72, 1332.00), 12 => (344.69, 1332.00), 13 => (959.90, 1332.00),
                     16 => (268.33, 1332.00), 18 => (1243.50, 1524.00), 20 => (4926.31, 1332.00),
                     21 => (1112.61, 1332.00), 23 => (344.69, 1332.00), 28 => (4926.31, 1332.00), 30 => (268.33, 1332.00))
const _RCW_UT = Dict(1 => (148.12, 444.00), 2 => nothing, 3 => (510.59, 1278.00), 4 => nothing, 6 => (1209.68, 444.00),
                     7 => (253.22, 444.00), 8 => nothing, 10 => nothing, 11 => nothing, 12 => nothing, 13 => nothing,
                     15 => nothing, 17 => (164.40, 444.00), 18 => (300.34, 444.00), 20 => nothing, 21 => nothing,
                     22 => (107.30, 444.00), 23 => (148.12, 444.00), 24 => nothing)
const _RCW_BM = Dict(1 => (384.59, 1332.00), 3 => (391.78, 2166.00), 4 => (256.95, 1332.00), 6 => nothing,
                     7 => (520.69, 1332.00), 8 => (201.94, 1524.00), 9 => (205.04, 1332.00), 10 => (350.97, 1332.00),
                     12 => (165.66, 1332.00), 15 => (8782.90, 1332.00), 16 => (253.34, 1332.00),
                     17 => (351.00, 1332.00), 18 => (253.34, 1332.00))
const _RCW_NC = Dict(1 => (644.03, 666.00), 2 => (685.84, 666.00), 3 => (742.19, 1083.00), 4 => (644.03, 666.00),
                     5 => (1499.04, 666.00), 6 => (644.03, 666.00), 7 => (646.05, 666.00), 8 => (1028.78, 666.00),
                     9 => (705.83, 762.00), 10 => (685.84, 666.00), 11 => (1238.89, 666.00), 12 => (6466.68, 666.00))

@testset "REGCAL small-tree HCOR vs live (CR + UT + BM + NC, all sub-models)" begin
    for (v, variant, gold) in (("centralrockies", FVSjl.CentralRockies(), _RCW_CR),
                               ("utah",           FVSjl.Utah(),           _RCW_UT),
                               ("bluemountains",  FVSjl.BlueMountains(),  _RCW_BM),
                               ("klamath",        FVSjl.Klamath(),        _RCW_NC))
        h = _rcw_hcor(v, variant)
        for (sp, g) in sort(collect(gold); by = first)
            if g === nothing
                @test h[sp] == 0f0
            else
                # live sums print F10.2 ⇒ the smallest SNX (UT 107.30) carries ~5e-5 relative rounding
                @test isapprox(exp(h[sp]), g[2] / g[1]; rtol = 1e-4)
            end
        end
    end
    # BM end-to-end: 3 cycles, every printed .sum field equal to FVSbm_g16
    fx = joinpath(_RCW_FX, "bluemountains", "regcal")
    rows = cd(fx) do
        txt = FVSjl.run_keyfile("regcal.key"; variant = FVSjl.BlueMountains(), output = :sum)
        [split(strip(l)) for l in split(txt, '\n') if occursin(r"^\d{4}\s", strip(l))]
    end
    gold = ["1990  60  6918 134  474 385  63  1.9  1821   978     0  4903     0     0     0     0     0 134  474 385  63  1.9      10   70    34    16.3 267 23",
            "2000  70  5666 140  473 357  71  2.1  2185  1435     0  7337     0     0     0     0     0 140  473 357  71  2.1      10   79    39    20.5 267 12",
            "2010  80  4644 148  475 331  70  2.4  2590  1820     0  9002     0     0     0     0     0 148  475 331  70  2.4      10   89    43    22.8 267 12",
            "2020  90  3816 157  479 311  78  2.7  3053  2269     0 11189     0     0     0     0     0 157  479 311  78  2.7       0    0     0    25.2 267 12"]
    @test length(rows) >= 4
    for (k, g) in enumerate(gold)
        @test rows[k] == split(g)
    end
end
