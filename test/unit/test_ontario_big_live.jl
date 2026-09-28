# test_ontario_big_live.jl — ON (Ontario) stands with ~1000 records vs LIVE FVSon_g16.
#
# ont_c1100 (1100 records, ~49500 TPH) and ont_c950 (950 records, ~64200 TPH): the 72 ON species cycled over a DBH
# ramp (5-45 cm), 2 ten-year cycles. Goldens: live .sum rows (ont_c*_live.rows) and live FVS_TreeList_East_Metric
# TPH/MortPH/BAPctile/MCuM of ont_c950 for 2004-2024 (ont_c950_tl_live.csv), compared EXACTLY.

using Test, FVSjl, SQLite, DBInterface
const _ONB = FVSjl
const _ONB_FX = joinpath(@__DIR__, "..", "fixtures", "ontario")

function _onb_run(stem::AbstractString; tl::Bool = false, extra::Vector{String} = String[], ncyc::Int = 2)
    dir = mktempdir(); out = String[]
    for ln in eachline(joinpath(_ONB_FX, "$stem.key"))
        tl && startswith(ln, "PROCESS") && append!(out, ["DATABASE", "SUMMARY", "TREELIDB", "END"])
        startswith(ln, "TREEDATA") && append!(out, extra)
        startswith(ln, "NUMCYCLE") && (push!(out, "NUMCYCLE" * lpad(string(Float64(ncyc)), 13)); continue)
        push!(out, ln)
        tl && startswith(ln, "TREEDATA") && push!(out, "TREELIST           0")
        tl && length(out) == 2 && append!(out, ["DATABASE", "DSNout", "$stem.db", "END"])
    end
    write(joinpath(dir, "$stem.key"), join(out, "\n") * "\n")
    cp(joinpath(_ONB_FX, "$stem.tre"), joinpath(dir, "$stem.tre"))
    txt = cd(() -> _ONB.run_keyfile("$stem.key"; variant = _ONB.Ontario(), output = :sum), dir)
    rows = [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)]
    tlr = Dict{Tuple{Int,Int},NamedTuple}()
    if tl
        db = SQLite.DB(joinpath(dir, "$stem.db"))
        for r in DBInterface.execute(db, "SELECT Year,TreeIndex,TPH,MortPH,BAPctile,MCuM FROM FVS_TreeList_East_Metric")
            tlr[(Int(r.Year), Int(r.TreeIndex))] = (; TPH = Float64(r.TPH), MortPH = Float64(r.MortPH),
                                                    BAPctile = Float64(r.BAPctile), MCuM = Float64(r.MCuM))
        end
    end
    return rows, tlr
end
_onb_live_rows(stem) = [split(l) for l in readlines(joinpath(_ONB_FX, "$(stem)_live.rows"))]
function _onb_live_tl()
    tl = Dict{Tuple{Int,Int},NamedTuple}()
    for ln in eachline(joinpath(_ONB_FX, "ont_c950_tl_live.csv"))
        (startswith(ln, "#") || startswith(ln, "Year")) && continue
        f = split(ln, ',')
        tl[(parse(Int, f[1]), parse(Int, f[2]))] = (; TPH = parse(Float64, f[4]), MortPH = parse(Float64, f[5]),
                                                   BAPctile = parse(Float64, f[6]), MCuM = parse(Float64, f[7]))
    end
    return tl
end

# (1) canada/on PRGPRM.F77 MAXTRE=6000 (jl had the base 3000 for every variant): grincr.f:74 LTRIP needs
# ITRN ≤ MAXTRE/3 = 2000 in ON, so the 1100-record stand triples in cycle 1 (jl: ≤1000 ⇒ no tripling, 2014 845 vs 830).
@testset "ON 1100-record stand: .sum rows == live (MAXTRE=6000 ⇒ tripling for ITRN ≤ 2000)" begin
    @test first(_onb_run("ont_c1100")) == _onb_live_rows("ont_c1100")
end

const _ONB_ROWS, _ONB_JL = _onb_run("ont_c950"; tl = true)
const _ONB_LIVE = _onb_live_tl()

# (2) canada/on/varmrt.f:246-253: a record whose adjusted kill would exhaust it gets TEMWK2(I)=PROB(I) — the WHOLE PROB,
# added to the WK2 it already carries from earlier JPASSes — so WK2 can exceed PROB; SHORT=SHORT+(XKILL−PROB+WK2).
# The over-kill drives MORTS's survivor sums (TN<0 ⇒ D10N ≤ DIA0 ⇒ stop after one pass), the FVS_TreeList MortPH
# reports WK2 unclamped, and UPDATE clamps the kill to PROB. jl killed only PROB−WK2 (2014: 852 TPH vs live 948).
@testset "ON 950-record hyper-dense stand 2014: TPH/MortPH per record == live (VARMRT exhaustion WK2+PROB)" begin
    for col in (:TPH, :MortPH), (k, v) in sort(collect(_ONB_LIVE))
        k[1] == 2014 || continue
        @test (k, col, haskey(_ONB_JL, k) ? _ONB_JL[k][col] : NaN) == (k, col, v[col])
    end
    @test _ONB_ROWS[1:2] == _onb_live_rows("ont_c950")[1:2]
end

# (3) grincr.f:74 evaluates LTRIP's ITRN.LE.(MAXTRE/3) at GRINCR entry, before CUTS TREDELs the zero-PROB records the
# previous mortality left (cuts.f:255-275): cycle 2 starts with 2850 records (> 2000) ⇒ no tripling, although only
# 1611 survive the TREDEL. jl tested the post-TREDEL count and tripled (2024: 4833 records vs live 1611).
@testset "ON 950-record stand 2024: record set, TPH/MortPH, .sum rows == live (LTRIP uses the pre-CUTS ITRN)" begin
    @test sort([k for k in keys(_ONB_JL) if k[1] == 2024]) == sort([k for k in keys(_ONB_LIVE) if k[1] == 2024])
    for col in (:TPH, :MortPH), (k, v) in sort(collect(_ONB_LIVE))
        k[1] == 2024 || continue
        @test (k, col, haskey(_ONB_JL, k) ? _ONB_JL[k][col] : NaN) == (k, col, v[col])
    end
    @test _ONB_ROWS == _onb_live_rows("ont_c950")
end

# (3b) The same LTRIP latch when a thin is scheduled: the .sum driver runs CUTS (and its zero-PROB TREDEL) before
# grow_cycle!, so the GRINCR-entry ITRN must be latched ahead of that call. THINBTA 2014 (940/ha): live does not
# triple in cycle 2; jl tripled on the post-thin 1544 records (2024: 793 vs live 790 TPH).
@testset "ON 950-record stand + THINBTA 2014: .sum rows == live (LTRIP ITRN latched before the driver's CUTS)" begin
    thin = [rpad("THINBTA", 10) * lpad("2014", 10) * lpad("940.0", 10) * lpad("1.0", 10)]
    @test first(_onb_run("ont_c950"; extra = thin, ncyc = 3)) == _onb_live_rows("ont_c950_thin")
end

# (4) Inventory FVS_TreeList BAPctile = CRATET's PCT: canada/on/cratet.f:128-131 `IND=IND1; RDPSRT(ITRN,DBH,IND,.FALSE.)`
# (species-major read-order seed, no re-seed) feeding the backdating DENSE's PCTILE. ont_c950 repeats every DBH across
# species (400 tie groups); an identity-seeded RDPSRT permuted the tied records' percentiles (553 records off).
@testset "ON 950-record stand 2004: BAPctile per record == live (CRATET IND1-seeded PCT)" begin
    for (k, v) in sort(collect(_ONB_LIVE))
        k[1] == 2004 || continue
        @test (k, haskey(_ONB_JL, k) ? _ONB_JL[k].BAPctile : NaN) == (k, v.BAPctile)
    end
end

# (5) A record MORTS empties is skipped by VOLS (vols.f:125 P≤0), so its FVS_TreeList MCuM (WK1·FT3toM3) is the DG that
# dgdriv.f DO 5 loaded into WK1 at the top of the cycle: cycle 1 the LSTART DO-220 dub (dgdriv.f:700-729), later the
# previous cycle's DG; TRIPLE copies it (triple.f:68). 1239 such records in 2014, 15 in 2024.
@testset "ON 950-record stand 2014-2024: MCuM of the emptied records == live (stale WK1 = DG at DGDRIV entry)" begin
    for (k, v) in sort(collect(_ONB_LIVE))
        (k[1] >= 2014 && v.TPH == 0) || continue
        @test (k, haskey(_ONB_JL, k) ? _ONB_JL[k].MCuM : NaN) == (k, v.MCuM)
    end
end
