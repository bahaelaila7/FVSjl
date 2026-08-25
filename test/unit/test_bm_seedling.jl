# test_bm_seedling.jl — BM (Blue Mountains) dense-seedling / establishing-regime regression guard.
#
# Locks in three BM small-tree bugs found by the western full-population FIA sweep (a one-directional
# BA/SDI/CCF/QMD OVER-growth on establishing stands whose live population is dbh=0.1 HT-null regeneration
# records sitting under standing-dead overstory). All three inflate the seedling crown dub / height / DBH:
#
#   1. AVHT40 dead-inclusive top height (bm/dense.f:285-297) — the DUBSCR crown dub at LSTART sees the
#      dead-inclusive AVH, which ranks trees by REAL DBH (IND) and sums real HT; jl zeroed t.dbh for the
#      HISTORY 8/9 dead (BA/CCF WK3 semantics) BEFORE computing AVH, sinking the tall dead below the 0.1"
#      seedlings ⇒ AVH≈seedling height (1.0) not the dead-inclusive top (67) ⇒ seedling CR dubbed ~80 not
#      the capped 95 ⇒ over-vigorous SMHTGF height growth. (bluemountains/crown.jl bm_crown_init_lstart!)
#   2. SMHTGF site index (bm/smhtgf.f:78) — SMHTGF re-reads the RAW SITEAR(I); regent.f's [SLO,SHI] clamp
#      is dead code for SMHTGF species. jl passed the clamped SI to bm_smhtgf, over-growing species whose
#      SITEAR is below SLO (LP 23.9 → clamped 30.5). (bluemountains/regent.jl small_tree_growth!/bm_esgent!)
#   3. LP(7) small-tree DG ht-dbh (bm/regent.f:519) — LHTDRG(7)=.FALSE. ⇒ the fixed CASE(7) -9.8752 curve
#      is unconditionally OVERRIDDEN by HTDBH (Curtis-Arney). jl used the fixed curve (dk≈1.14) ⇒ LP seedling
#      DBH grew ~2.7× vs HTDBH (dk≈0.5). Now LP uses HTDBH like DF/GF/ES. (bluemountains/regent.jl)
#
# Two both-sides-traced establishing stands, diffed vs the live FVSbm_clean .sum golden (10 cycles). Pre-fix
# these ran BA/SDI/CCF up to +40-130% high by the last cycle; post-fix they are bit-exact-or-cornered.
using FVSjl, Test

const _DIR = @__DIR__
const _DB  = joinpath(_DIR, "..", "fixtures", "bm_seedling", "bm_seedling.db")

# .sum is fixed-width (sumout.f 9014): year[1:4] age[5:8] TPA[9:14] BA[15:18] SDI[19:23] CCF[24:27]
# TopHt[28:31] QMD[32:36]. Parse by column (a whitespace split merges SDI+CCF once CCF≥1000).
const _COLS = (:tpa => (9,14), :ba => (15,18), :sdi => (19,23), :ccf => (24,27),
               :topht => (28,31), :qmd => (32,36))
function _rows(text)
    out = Dict{Int,Dict{Symbol,Float64}}()
    for ln in split(text, '\n')
        s = rstrip(ln); length(s) < 36 && continue
        y = tryparse(Int, strip(s[1:4])); (y === nothing || y < 1000 || y > 3000) && continue
        d = Dict{Symbol,Float64}()
        for (name,(a,b)) in _COLS
            v = tryparse(Float64, strip(s[a:b])); v === nothing && (v = 0.0); d[name] = v
        end
        out[y] = d
    end
    out
end

# NOTE: FVS keywords are column-sensitive — every line MUST start at column 1 (no indentation), and NUMCYCLE's
# value must sit in cols 11+. Build the keyfile from left-justified lines, not an indented triple-quoted block.
function _keytext(cn)
    join(["STDIDENT",
          cn,
          "DATABASE",
          "DSNin",
          abspath(_DB),
          "StandSQL",
          "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'",
          "EndSQL",
          "TreeSQL",
          "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'",
          "EndSQL",
          "END",
          "NUMCYCLE         10",
          "ECHOSUM",
          "PROCESS",
          "STOP"], "\n") * "\n"
end

# Per-column relative tolerance. Chosen to ACCEPT the converged cornered state yet FAIL a revert of any of
# the three fixes (which drive BA/SDI/CCF/QMD +40-130% high by the last establishing cycle). Observed cornered
# residuals over the 11-cycle projection: BA ≤0.07, SDI ≤0.08, CCF ≤0.04, QMD ≤0.06; TopHt ≤0.12 (a cycle-boundary
# self-thin/tie straddle, ≤3 ft absolute); TPA ≤0.11 (the documented BM late-cycle self-thin residual, open #140).
const _TOL = Dict(:tpa=>0.15, :ba=>0.12, :sdi=>0.12, :ccf=>0.10, :topht=>0.15, :qmd=>0.10)

@testset "BM dense-seedling establishing-regime (3 small-tree fixes)" begin
    if !isfile(_DB)
        @test_skip "bm_seedling fixture db missing"
    else
        for cn in ("449747082489998", "15144796010497")
            golddf = joinpath(_DIR, "..", "fixtures", "bm_seedling", "$cn.live.sum")
            if !isfile(golddf)
                @test_skip "$cn golden missing"; continue
            end
            G = _rows(read(golddf, String))
            dir = mktempdir(); kf = joinpath(dir, "s.key"); write(kf, _keytext(cn))
            J = _rows(FVSjl.run_keyfile(kf; variant = FVSjl.BlueMountains()))
            @test !isempty(J)
            for y in sort(collect(keys(G)))
                haskey(J, y) || continue
                for (name, _) in _COLS
                    gv = G[y][name]; jv = J[y][name]
                    rel = gv == 0 ? abs(jv) : abs(jv - gv) / abs(gv)
                    @test rel <= _TOL[name]
                end
            end
        end
    end
end
