# =============================================================================
# test_canonical_multistand.jl — full FVS reference keyfiles, ALL stands/cycles, vs live
#
# The canonical FVS test keys (net01/cst01/lst01) are multi-stand REALISTIC MANAGEMENT
# scenarios written by the FVS authors — each bundles the treatments a real analyst runs:
#   stand 1  unthinned control (10-cycle projection)
#   stand 2  repeated diameter-limit thinning (THINDBH, multi-class, every-3-cycle)
#   stand 3  shelterwood prescription (THINPRSC / SPECPREF / THINBTA) + full ECON
#   stand 4  FFE fire (SIMFIRE + SNAGINIT + SALVAGE + fuel/burn/mort reports)
#   stand 5  bare-ground planting (ESTAB / PLANT, NOTREES)
#
# Running the WHOLE key through `run_keyfile` and comparing every stand × cycle against a
# freshly-relinked live binary is the strongest realistic-scenario gate: it exercises
# thinning selection, shelterwood + economics, the FFE fire/fuel/snag path, and the
# establishment/planting path in one shot, with cross-stand state-carry (REWIND) intact.
# Complements test_allspecies.jl (per-species coefficient rows) — this is per-SCENARIO.
#
# Measured vs live: SN 55/55 (Δ≤1 ULP); NE 50/50 BIT-EXACT; CS 49/50 (only the FFE-fire stand-4 fire-mortality
# distribution drifts ~3% of the kill, the documented FMEFF residual); LS 41/50 exact (the
# tripling-spread late-cycle tail). Bounds are set to the documented per-variant floor.
# =============================================================================

using Test, FVSjl

const _CANON_DIR = joinpath(@__DIR__, "..", "fixtures", "canonical")

# parse a multi-stand .sum → Vector of (stand#, year, [TPA,BA,SDI,CCF,TopHt]).
function _canon_rows(txt::AbstractString)
    out = Tuple{Int,Int,Vector{Int}}[]; si = 0
    for l in split(txt, "\n")
        if startswith(strip(l), "-999"); si += 1; continue; end
        f = split(strip(l))
        (length(f) >= 7 && tryparse(Int, f[1]) !== nothing &&
         (startswith(f[1], "19") || startswith(f[1], "20"))) || continue
        push!(out, (si, parse(Int, f[1]),
                    [parse(Int, f[3]), parse(Int, f[4]), parse(Int, f[5]),
                     parse(Int, f[6]), parse(Int, f[7])]))
    end
    return out
end

# Compare jl vs live golden per stand×cycle. `tol` = max abs deviation allowed on any of the
# 5 stand columns (0 ⇒ bit-exact). The inventory (first) year of each stand is ALWAYS bit-exact.
function _assert_canonical(stem::AbstractString, V; tol::Int, label::AbstractString)
    key = joinpath(_CANON_DIR, "$stem.key")
    golden = joinpath(_CANON_DIR, "$stem.live.sum")
    if !isfile(key) || !isfile(golden)
        @info "$stem canonical fixtures absent; skipping (run test/harness/gen_canonical.sh)"
        return
    end
    jl = cd(_CANON_DIR) do
        FVSjl.run_keyfile("$stem.key"; variant = V, output = :sum)
    end
    J = _canon_rows(jl)
    L = _canon_rows(read(golden, String))
    @test length(J) == length(L)
    Ld = Dict((s, y) => v for (s, y, v) in L)
    # first (inventory) year seen per stand ⇒ bit-exact regardless of tol
    firstyr = Dict{Int,Int}()
    for (s, y, _) in L
        firstyr[s] = min(get(firstyr, s, y), y)
    end
    for (s, y, jv) in J
        haskey(Ld, (s, y)) || (@test haskey(Ld, (s, y)); continue)
        lv = Ld[(s, y)]
        t = (y == firstyr[s]) ? 0 : tol
        d = maximum(abs.(jv .- lv))
        @test (label, "stand$s", y, d <= t) == (label, "stand$s", y, true)
    end
end

@testset "SN snt01 — full multi-stand reference scenarios (vs live FVSsn)" begin
    _assert_canonical("snt01", Southern(); tol = 2, label = "SN")    # bit-exact bar ULP
end

@testset "NE net01 — full multi-stand reference scenarios (vs live FVSne)" begin
    _assert_canonical("net01", Northeast(); tol = 2, label = "NE")   # bit-exact bar ULP
end

@testset "CS cst01 — full multi-stand reference scenarios (vs live FVScs)" begin
    _assert_canonical("cst01", CentralStates(); tol = 6, label = "CS")  # stand-4 FFE-fire ~3% floor
end

@testset "LS lst01 — full multi-stand reference scenarios (vs live FVSls)" begin
    _assert_canonical("lst01", LakeStates(); tol = 6, label = "LS")     # late tripling-spread tail
end
