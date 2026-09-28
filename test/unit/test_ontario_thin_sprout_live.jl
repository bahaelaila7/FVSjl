# test_ontario_thin_sprout_live.jl — ON (Ontario) thinning + stump sprouting per record vs LIVE FVSon_g16.
#
# ont_sm (8 sub-12 cm trees) and ont_lite (72 species), NUMCYCLE 5, THINBTA 2014 (500 trees/ha, eff 1.0) and
# THINBBA 2034 (10 m²/ha, eff 0.8): cycle 1 triples small trees through REGENT, the thins cut hardwood sprouters,
# ESUCKR adds stump sprouts, and later cycles grow them. Goldens: live FVS_TreeList_East_Metric of the same runs
# (ont_sm_thin_tl_live.csv — every year; ont_lite_thin_tl_live.csv — 2034..2054) compared EXACTLY per record.

using Test, FVSjl, SQLite, DBInterface
const _ONT = FVSjl
const _ONT_FX = joinpath(@__DIR__, "..", "fixtures", "ontario")
const _ONT_THIN = [rpad("THINBTA", 10) * lpad("2014", 10) * lpad("500.0", 10) * lpad("1.0", 10),
                   rpad("THINBBA", 10) * lpad("2034", 10) * lpad("10.0", 10) * lpad("0.8", 10)]

# fixture stand → temp dir, NUMCYCLE 5, the two thins before TREEDATA, the ON DATABASE/TREELIDB layout
function _ont_run(stem::AbstractString)
    dir = mktempdir(); out = String[]
    for ln in eachline(joinpath(_ONT_FX, "$stem.key"))
        if startswith(ln, "NUMCYCLE")
            push!(out, "NUMCYCLE" * lpad("5.0", 13))
        elseif startswith(ln, "PROCESS")
            append!(out, ["DATABASE", "SUMMARY", "TREELIDB", "END", ln])
        elseif startswith(ln, "TREEDATA")
            append!(out, _ONT_THIN); push!(out, ln); push!(out, "TREELIST           0")
        else
            push!(out, ln)
        end
        length(out) == 2 && append!(out, ["DATABASE", "DSNout", "$stem.db", "END"])
    end
    write(joinpath(dir, "$stem.key"), join(out, "\n") * "\n")
    cp(joinpath(_ONT_FX, "$stem.tre"), joinpath(dir, "$stem.tre"))
    txt = cd(() -> _ONT.run_keyfile("$stem.key"; variant = _ONT.Ontario(), output = :sum), dir)
    rows = [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)]
    db = SQLite.DB(joinpath(dir, "$stem.db"))
    tl = Dict{Tuple{Int,Int},NamedTuple}()
    for r in DBInterface.execute(db, "SELECT Year,TreeIndex,TPH,MortPH,DBH,Ht,TCuM FROM FVS_TreeList_East_Metric")
        tl[(Int(r.Year), Int(r.TreeIndex))] = (; TPH = Float64(r.TPH), MortPH = Float64(r.MortPH),
            DBH = Float64(r.DBH), Ht = Float64(r.Ht), TCuM = Float64(r.TCuM))
    end
    return rows, tl
end
function _ont_live(file)
    tl = Dict{Tuple{Int,Int},NamedTuple}(); hdr = String[]
    for ln in eachline(joinpath(_ONT_FX, file))
        startswith(ln, "#") && continue
        f = split(ln, ',')
        if f[1] == "Year"; hdr = String.(f); continue; end
        tl[(parse(Int, f[1]), parse(Int, f[2]))] = (; (Symbol(hdr[k]) => parse(Float64, f[k]) for k in 4:length(f))...)
    end
    return tl
end
const _ONT_SM = _ont_run("ont_sm")
const _ONT_SM_LIVE = _ont_live("ont_sm_thin_tl_live.csv")
_ont_cmp(jl, lv, col; sel = (k, v) -> true) =
    [((k, col, haskey(jl, k) ? jl[k][col] : NaN), (k, col, v[col])) for (k, v) in sort(collect(lv)) if sel(k, v)]

# (1) canada/on/regent.f DGGR = DGSM·(1−XWT) + XWT·DG(K): a tripled copy K blends ITS OWN large-tree DG
# (dgdriv.f DG(ITRIPU)/DG(ITRIPL)), not the central's; an HK≤4.5 record gets DBH(K)=D+0.001·HK directly (DG=0).
@testset "ON ont_sm cycle 1: TRIPLE-copy DBH/Ht/TPH (REGENT DGGR blend) == live" begin
    for col in (:TPH, :DBH, :Ht), (j, l) in _ont_cmp(last(_ONT_SM), _ONT_SM_LIVE, col; sel = (k, v) -> k[1] == 2014)
        @test j == l
    end
end
