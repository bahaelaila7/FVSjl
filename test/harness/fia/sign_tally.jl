# sign_tally.jl — signed last-cycle jl−oracle POPULATION sign-tally (any variant, any regime).
# Usage: VARIANT=BM REGIME=none julia --project=. test/harness/fia/sign_tally.jl <standlist>
# Adjudicates straddle-vs-bias: a one-directional material tally = real deterministic bug; mixed/balanced
# with small Σ = realization floor (required for establishment, whose ESRANN stream is seed-fixed).
# Reuses ledger_fia.jl helpers (build_subdb/run_live/keytext/parse_sum10/BIN/VAR/MASTER).
# Records SIGNED residuals (the ledger only keeps abs). Also captures per-stand every-cycle
# signed TPA/BA so we can see WHEN the divergence enters.
include(joinpath(@__DIR__, "ledger_fia.jl"))

const V = get(ENV, "VARIANT", "IE"); const REGIME = get(ENV, "REGIME", "none")
listfile = ARGS[1]
cns = [split(strip(l),'\t')[1] for l in eachline(listfile) if !isempty(strip(l))]
dir = mktempdir(); sub = joinpath(dir,"sub.db")
print(stderr, "building sub-DB ($(length(cns)))..."); flush(stderr)
build_subdb(cns, sub); println(stderr, "ok")
bin = BIN[V]; var = VAR[V]

# column indices into parse_sum10 vals: 1=TPA 2=BA 3=SDI 4=CCF 5=TopHt 6=QMD 7=TCuFt
mat(lv,jv) = (ad=abs(lv-jv); ad>1.0+1e-6 && (lv==0 ? true : ad/abs(lv)>=0.05))  # material: >1 abs & ≥5% rel

# PLANT must use the same CALENDAR-year date as ledger_fia.jl main() (INV_YEAR+10, a cycle boundary). Omitting it
# emits the cycle-number form `PLANT 2.0`, a different scheduler path, so the tally disagreed with the ledger.
invyr = Dict{String,Int}()
let db = SQLite.DB(sub)
    for r in DBInterface.execute(db, "SELECT STAND_CN, INV_YEAR FROM FVS_STANDINIT_COND")
        (r[:STAND_CN] === missing || r[:INV_YEAR] === missing) && continue
        invyr[String(r[:STAND_CN])] = Int(r[:INV_YEAR])
    end
    SQLite.close(db)
end
plantyr_of(cn) = (REGIME == "plant" && haskey(invyr, cn)) ? invyr[cn] + 10 : 0

rows = NamedTuple[]
for cn in cns
    py = plantyr_of(cn)
    live,_ = run_live(bin, cn, sub, REGIME, dir, py)
    keyf = joinpath(dir,"jl.key"); write(keyf, keytext(cn, sub, REGIME, py))
    jl = try FVSjl.run_keyfile(keyf; variant=var) catch e; ""; end
    L = parse_sum10(live); J = parse_sum10(jl)
    (isempty(L) || isempty(J)) && continue
    Jd = Dict(y=>v for (y,v) in J)
    ys = sort([y for (y,_) in L if haskey(Jd,y)]); isempty(ys) && continue
    ly = ys[end]
    lv = Dict(L)[ly]; jv = Jd[ly]
    push!(rows, (cn=cn, year=ly,
        tpa_o=lv[1], tpa_j=jv[1], dtpa=jv[1]-lv[1],
        ba_o=lv[2],  ba_j=jv[2],  dba=jv[2]-lv[2],
        toph_o=lv[5],toph_j=jv[5],dtoph=jv[5]-lv[5],
        qmd_o=lv[6], qmd_j=jv[6], dqmd=jv[6]-lv[6],
        tcuft_o=lv[7],tcuft_j=jv[7],dtcuft=jv[7]-lv[7],
        tpa_mat=mat(lv[1],jv[1]), ba_mat=mat(lv[2],jv[2]), tcuft_mat=mat(lv[7],jv[7])))
end

function tally(rows, key, matkey)
    over = count(r->getfield(r,matkey) && getfield(r,key)>0, rows)
    under= count(r->getfield(r,matkey) && getfield(r,key)<0, rows)
    exact= count(r->!getfield(r,matkey), rows)
    net  = sum(getfield(r,key) for r in rows; init=0.0)
    (over=over, under=under, exact=exact, net=round(net,digits=1))
end
println("\n=== SIGNED SIGN-TALLY  $V $REGIME  (n=$(length(rows)) ran) ===")
for (nm,key,mk) in (("TPA",:dtpa,:tpa_mat),("BA",:dba,:ba_mat),("TCuFt",:dtcuft,:tcuft_mat))
    t = tally(rows,key,mk)
    println(rpad(nm,6), " material OVER(jl>O)=", t.over, "  UNDER(jl<O)=", t.under,
            "  exact/immaterial=", t.exact, "  Σsigned(jl−O)=", t.net)
end
# the material TPA movers, sorted by |dtpa|
println("\n--- material TPA movers (|Δ|≥5 & ≥5%) ---")
movers = sort([r for r in rows if r.tpa_mat], by=r->-abs(r.dtpa))
for r in first(movers, min(30,length(movers)))
    println(rpad(r.cn,20), " ", r.year, "  O_TPA=", round(r.tpa_o), " jl=", round(r.tpa_j),
            " dTPA=", round(r.dtpa), "  dBA=", round(r.dba,digits=1),
            " dTopHt=", round(r.dtoph,digits=1), " dTCuFt=", round(r.dtcuft), " dQMD=", round(r.dqmd,digits=2))
end
println("\nTPA material movers: ", length(movers), " / ", length(rows))
