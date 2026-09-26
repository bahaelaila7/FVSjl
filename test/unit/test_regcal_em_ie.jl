# test_regcal_em_ie.jl — the LSTART small-tree HEIGHT calibration (REGCAL, em/regent.f + ie/regent.f label 40)
# vs FVSem_g16 / FVSie_g16.
#
# Fixture (one per variant): the shipped PN inventory (backdated by its measured DG, with 2 dead records) plus six
# seedlings/saplings of each of 12 species carrying a measured HTG, one per DBH class 0.6-2.0" ⇒ 12-way DBH ties.
# Every REGCAL sub-model is exercised: EM EMVAR/NIVAR(LL)/TTVAR(LM)/CRVAR/UTVAR(aspen, juniper), IE NIVAR/TTVAR/
# CRVAR/UTVAR. Golden = the live per-species "SUMS FOR SPECIES n: SNP SNX SNY" (DEBUG REGENT, 2026-09-26) ⇒
# CORNEW = SNY/SNX; trapped (outside [0.0821,12.1825]) and N<NCALHT species keep HCOR 0.
#
# What had to be right: RHCON = 1 for every non-NI species (jl gave all EM species the NI constant ⇒ CW CORNEW 2.50
# vs 1.004); NTYR = IFINTH (IE used FINT ⇒ NPER 2 vs 1 ⇒ every NIVAR SNX ~2x); the cratet backdating DENSE's
# BA/RELDEN/AVH/RELDM1/PCCF/PCT (the IMC=9 dead's D=0 CCF of 0.001·P; PCT on the cratet.f:153 IND, whose ties
# differ from a stable sort); the IE CR/UT arms (unported before). EM also runs 3 cycles bit-exact on the .sum.
using Test
using FVSjl

const _RC_FX = joinpath(@__DIR__, "..", "fixtures")

function _rc_hcor(v, variant)
    key = joinpath(_RC_FX, v, "regcal", "regcal.key")
    s, _ = FVSjl.initialize(key; variant = variant)
    FVSjl.notre!(s); FVSjl.setup_growth!(s)
    return s.calib.htg_cor_init
end

# species => (SNX, SNY) from the live sums; `nothing` ⇒ HCOR stays 0 (N<5 or CORNEW trapped)
const _RC_EM = Dict(3 => (641.31, 1083.0), 4 => (771.97, 666.0), 5 => (568.28, 666.0), 6 => nothing,
                    7 => (1711.37, 666.0), 8 => nothing, 9 => (581.60, 666.0), 10 => nothing,
                    11 => (167.82, 666.0), 12 => (2401.57, 666.0), 13 => (167.82, 666.0), 14 => (167.82, 666.0),
                    16 => (167.82, 666.0), 17 => (2401.57, 666.0), 18 => nothing)
const _RC_IE = Dict(1 => (668.18, 666.0), 3 => (688.55, 1083.0), 4 => nothing, 5 => (617.31, 666.0),
                    6 => (540.98, 666.0), 7 => (704.42, 666.0), 8 => nothing, 9 => (668.75, 666.0), 10 => nothing,
                    13 => (772.93, 666.0), 14 => (667.73, 666.0), 16 => nothing, 18 => (2401.57, 666.0),
                    19 => (550.78, 666.0), 21 => (2401.57, 666.0), 23 => nothing)

@testset "REGCAL small-tree HCOR vs live (EM + IE, all sub-models)" begin
    for (v, variant, gold) in (("easternmontana", FVSjl.EasternMontana(), _RC_EM),
                               ("inlandempire",   FVSjl.InlandEmpire(),   _RC_IE))
        h = _rc_hcor(v, variant)
        for (sp, g) in sort(collect(gold); by = first)
            if g === nothing
                @test h[sp] == 0f0
            else
                @test isapprox(exp(h[sp]), g[2] / g[1]; rtol = 1e-4)
            end
        end
    end
    # EM end-to-end: 3 cycles, every printed .sum field equal to FVSem_g16
    fx = joinpath(_RC_FX, "easternmontana", "regcal")
    rows = cd(fx) do
        txt = FVSjl.run_keyfile("regcal.key"; variant = FVSjl.EasternMontana(), output = :sum)
        [split(strip(l)) for l in split(txt, '\n') if occursin(r"^\d{4}\s", strip(l))]
    end
    gold = ["1990  60  6427 129  455 286  63  1.9  4930  1006     0  5177     0     0     0     0     0 129  455 286  63  1.9      10   89    54    16.8 201 22",
            "2000  70  5645 164  536 346  67  2.3  5276  1243     0  5987     0     0     0     0     0 164  536 346  67  2.3      10  115    83    17.8 266 22",
            "2010  80  4716 183  565 361  60  2.7  5601  1307     0  6300     0     0     0     0     0 183  565 361  60  2.7      10   74    87    16.3 201 32",
            "2020  90  3932 198  581 368  61  3.0  5475  1252     0  6069     0     0     0     0     0 198  581 368  61  3.0       0    0     0    13.9 901 31"]
    @test length(rows) >= 4
    for (k, g) in enumerate(gold)
        @test rows[k] == split(g)
    end
    # IE end-to-end: 3 cycles, every printed .sum field equal to FVSie_g16 — needs the Fortran-shaped IE REGENT, the
    # DO-220 cycle-1 WK1 dub, the calibration's cratet IND / dead-inclusive RMSQD and the single LM site conversion.
    fx = joinpath(_RC_FX, "inlandempire", "regcal")
    rows = cd(fx) do
        txt = FVSjl.run_keyfile("regcal.key"; variant = FVSjl.InlandEmpire(), output = :sum)
        [split(strip(l)) for l in split(txt, '\n') if occursin(r"^\d{4}\s", strip(l))]
    end
    gold = ["1990  60  6427 129  455 246  63  1.9  4215  1012     0  5233     0     0     0     0     0 129  455 246  63  1.9      10  108    46    16.9 201 23",
            "2000  70  5201 156  506 289  68  2.3  4827  1428     0  7026     0     0     0     0     0 156  506 289  68  2.3      10  109    59    20.4 201 22",
            "2010  80  4107 175  531 306  72  2.8  5320  1848     0  9308     0     0     0     0     0 175  531 306  72  2.8      10   62    62    23.1 267 12",
            "2020  90  3277 192  546 314  75  3.3  5315  2185     0 10546     0     0     0     0     0 192  546 314  75  3.3       0    0     0    24.3 267 12"]
    @test length(rows) >= 4
    for (k, g) in enumerate(gold)
        @test rows[k] == split(g)
    end
end
