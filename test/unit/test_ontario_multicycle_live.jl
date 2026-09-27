# test_ontario_multicycle_live.jl — ON (Ontario) 5-cycle projections vs LIVE FVSon_g16.
#
# The six ON fixture stands (ont01/ont_all/ont_big/ont_lite/ont_mh/ont_sm) re-run with NUMCYCLE 5 (2004-2054): every
# .sum row must equal live's (goldens *_5cyc_live.rows), and the 72-species ont_all stand — whose records thin from
# hyper-dense through the density and background mortality regimes and triple twice — must match live's
# FVS_TreeList per record in every cycle (ont_all_5cyc_tl_live.csv; TPH/MortPH/DBH/Ht compared EXACTLY).

using Test, FVSjl, SQLite, DBInterface
const _ONM = FVSjl
const _ONM_FX = joinpath(@__DIR__, "..", "fixtures", "ontario")

# Copy a fixture stand into a temp dir with NUMCYCLE=ncyc (and, with `tl`, the ON DATABASE/TREELIDB layout).
function _onm_run(stem::AbstractString, ncyc::Int; tl::Bool = false)
    dir = mktempdir()
    out = String[]
    for ln in eachline(joinpath(_ONM_FX, "$stem.key"))
        if startswith(ln, "NUMCYCLE")
            push!(out, "NUMCYCLE" * lpad(string(Float64(ncyc)), 13))
        elseif tl && startswith(ln, "PROCESS")
            append!(out, ["DATABASE", "SUMMARY", "TREELIDB", "END", ln])
        elseif tl && startswith(ln, "TREEDATA")
            append!(out, [ln, "TREELIST           0"])
        else
            push!(out, ln)
        end
        tl && length(out) == 2 && append!(out, ["DATABASE", "DSNout", "$stem.db", "END"])
    end
    write(joinpath(dir, "$stem.key"), join(out, "\n") * "\n")
    cp(joinpath(_ONM_FX, "$stem.tre"), joinpath(dir, "$stem.tre"))
    txt = cd(() -> _ONM.run_keyfile("$stem.key"; variant = _ONM.Ontario(), output = :sum), dir)
    rows = [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)]
    return rows, joinpath(dir, "$stem.db")
end
_onm_live_rows(stem) = [split(l) for l in readlines(joinpath(_ONM_FX, "$(stem)_5cyc_live.rows"))]

function _onm_live_tl()
    rows = Dict{Tuple{Int,Int},NamedTuple}(); hdr = String[]
    for ln in eachline(joinpath(_ONM_FX, "ont_all_5cyc_tl_live.csv"))
        startswith(ln, "#") && continue
        f = split(ln, ',')
        if f[1] == "Year"; hdr = String.(f); continue; end
        rows[(parse(Int, f[1]), parse(Int, f[2]))] =
            (; (Symbol(hdr[k]) => parse(Float64, f[k]) for k in 4:length(f))...)
    end
    return rows
end
function _onm_jl_tl()
    _, dbp = _onm_run("ont_all", 5; tl = true)
    db = SQLite.DB(dbp)
    rows = Dict{Tuple{Int,Int},NamedTuple}()
    for r in DBInterface.execute(db, "SELECT Year,TreeIndex,TPH,MortPH,DBH,Ht FROM FVS_TreeList_East_Metric")
        rows[(Int(r.Year), Int(r.TreeIndex))] = (; TPH = Float64(r.TPH), MortPH = Float64(r.MortPH),
                                                 DBH = Float64(r.DBH), Ht = Float64(r.Ht))
    end
    return rows
end
const _ONM_LIVE_TL = _onm_live_tl()
const _ONM_JL_TL = _onm_jl_tl()
_onm_col(col; sel = (k, v) -> true) =
    [(k, v[col], haskey(_ONM_JL_TL, k) ? _ONM_JL_TL[k][col] : NaN) for (k, v) in sort(collect(_ONM_LIVE_TL)) if sel(k, v)]

@testset "ON ont_all 5-cycle: TreeList record set == live" begin
    @test sort(collect(keys(_ONM_JL_TL))) == sort(collect(keys(_ONM_LIVE_TL)))
end

# (1) canada/on/htgf.f DO 20: X = PROB(I)*DBH10**2 = P*(D*D) feeds BA10/QMD10 (HT10 = HTONT(DBH10,QMD10,BA10));
# jl's (P*D)*D put QMD10 1 ULP off ⇒ HT10 1 ULP ⇒ HTG (a difference of two ~2400-ft Penner heights) off on the
# hyper-dense 72-species stand.
@testset "ON ont_all cycle 1: Ht per record == live (htgf.f QMD10 = P*(D*D))" begin
    for (k, lv, jv) in _onm_col(:Ht; sel = (k, v) -> k[1] == 2014)
        @test (k, jv) == (k, lv)
    end
end

# (2) canada/on/avht40.f: TARG = 100.0/HAtoACR (the 100 largest trees per HECTARE, 40.47/ac) — not the imperial 40/ac.
# AVH is the .sum top height; ont_all 2044 was 25 m vs live 24 with the 40/ac window.
@testset "ON fixture stands, 5 cycles: every .sum row == live (avht40.f metric TARG)" begin
    for stem in ("ont01", "ont_all", "ont_big", "ont_lite", "ont_mh", "ont_sm")
        @test (stem, first(_onm_run(stem, 5))) == (stem, _onm_live_rows(stem))
    end
end

# (3) canada/on/varmrt.f: X = … + MB3(J)*DM**2 (= MB3*(DM*DM)), and a background-regime VARMRT (TOKILL=0) totals the
# MORTS WK2 with a SEQUENTIAL `DO I=1,ITRN: TOKILL=TOKILL+WK2(I)` (Julia's `sum` is pairwise). Either one-ULP drift
# moves every record's kill (EFFTR → TEMWK2 → ADJUST).
@testset "ON ont_all 5-cycle: TPH/MortPH/DBH/Ht per record == live (VARMRT EFFTR + TOKILL)" begin
    for col in (:TPH, :MortPH, :DBH, :Ht), (k, lv, jv) in _onm_col(col)
        @test (k, col, jv) == (k, col, lv)
    end
end
