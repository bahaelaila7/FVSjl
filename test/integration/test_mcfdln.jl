# test_mcfdln.jl — MCFDLN / BFFDLN log-linear form-model coefficients (sdefln.f → FVSsn vols.f).
#
# MCFDLN sets the cubic form coefs CFLA0/CFLA1, BFFDLN the board coefs BFLA0/BFLA1. The volume
# model then folds an implied defect % into ICDF/IBDF: VOLCOR = exp(B0 + B1·ln(V)), reduction =
# (V−VOLCOR)/V, where V is the pulpwood (MCFV−SCFV) for cubic and BFV for board (vols.f:303-310).
#
# ⚠ NO live-Fortran oracle: the SN Fortran build FPE-crashes whenever the form model is activated —
# `vols.f:306` evaluates `ALOG(TEMVOL)` BEFORE the `TEMVOL==0` guard at `:308`, so any tree with
# MCFV==SCFV (zero pulpwood) crashes it (verified: every non-default B1 cores-dumps). FVSjl
# implements the model robustly (guards `temvol > 0`), so it can't be bit-exact-validated here; see
# docs/DIVERGENCES.md. This test pins the two things that ARE checkable without the oracle:
#   1. the default (B1=1, B0=0) is a no-op — byte-identical to the same stand with no MCFDLN line;
#   2. an activated form model (B1=0.9) deterministically lowers merch cubic in every cycle.

using Test, FVSjl

const _MF_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_mf_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_mfcol(r, c) = parse(Float64, r[c])

@testset "MCFDLN form-model coefficients" begin
    key = joinpath(_MF_DIR, "mcfdln_override.key")
    tre = joinpath(_MF_DIR, "mcfdln_override.tre")
    if !isfile(key)
        @test_skip "mcfdln_override scenario not available"
    else
        act = _mf_rows(FVSjl.run_keyfile(key; faithful = true))
        @test _mf_rows(FVSjl.run_keyfile(key; faithful = true)) == act   # deterministic

        # no-MCFDLN twin (= the default-coefficient path, since B0=0/B1=1 is a no-op)
        nokey = tempname() * ".key"
        write(nokey, join(filter(l -> !startswith(l, "MCFDLN"), readlines(key)), "\n") * "\n")
        cp(tre, replace(nokey, ".key" => ".tre"); force = true)
        none = _mf_rows(FVSjl.run_keyfile(nokey; faithful = true))
        @test length(act) == length(none)

        # 2. activation (B1=0.9 in the scenario) lowers merch cubic (col 10) in every row, by >30 cuft.
        @test all(_mfcol(act[i], 10) <= _mfcol(none[i], 10) + 1 for i in 1:length(act))
        @test all(_mfcol(none[i], 10) - _mfcol(act[i], 10) > 30 for i in 1:length(act))

        # 1. the DEFAULT coefficients (B0=0, B1=1) must be a true no-op vs the no-keyword run.
        defkey = tempname() * ".key"
        write(defkey, join(map(readlines(key)) do l
            startswith(l, "MCFDLN") ? rpad("MCFDLN", 10) * lpad("0.", 10) * lpad("0.0", 10) * lpad("1.0", 10) : l
        end, "\n") * "\n")
        cp(tre, replace(defkey, ".key" => ".tre"); force = true)
        defrows = _mf_rows(FVSjl.run_keyfile(defkey; faithful = true))
        @test defrows == none      # CFLA0=0/CFLA1=1 ⇒ byte-identical to no MCFDLN
    end
end
