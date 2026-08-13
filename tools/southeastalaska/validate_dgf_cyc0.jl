# Validate the ported AK DGF equation + coefficients vs the live FVSak oracle DEBUG dump.
# Reproduces WK2 (=DDS) per tree from the ported coefficient arrays + the dumped raw inputs
# (D, PBAL, PRD, CR, TEMEL, TEMSLP, TEMSASP, SSITE, BRAT, COR, DGCON), then compares to the
# oracle's dumped BASEDG / DGPRED / WK2. All Float32 (FVS is single precision).

# --- pull the SHIPPED module coefficients (proves the committed src/variants/southeastalaska code) ---
using FVSjl
const M = FVSjl
DGCONB1 = M.AK_DGCONB1; DGDISQ = M.AK_DGDISQ; DGLD = M.AK_DGLD; DGDBAL = M.AK_DGDBAL
DGRD = M.AK_DGRD; DGLNCR = M.AK_DGLNCR; DGEL = M.AK_DGEL; DGSLOP = M.AK_DGSLOP
DGSASP = M.AK_DGSASP; DGLNSI = M.AK_DGLNSI
YR = 10.0f0
dgmult(sp) = M.ak_dgmult(sp)
akbrat(sp, d) = M.ak_bratio(sp, d)

# --- parse the dump ---
lines = readlines(joinpath(@__DIR__, "akdbg.out"))
num(s) = parse(Float32, s)
# grab last float on a "key= value" segment
function grabf(line, key)
    idx = findfirst(key, line); idx === nothing && return nothing
    rest = line[last(idx)+1:end]
    m = match(r"[-+]?[0-9]*\.?[0-9]+(?:[Ee][-+]?[0-9]+)?", rest)
    m === nothing ? nothing : num(m.match)
end

function main(lines)
cur_ispc = 0; cur_temel=0f0; cur_temslp=0f0; cur_temsasp=0f0; cur_ssite=0f0
maxerr_dds = 0.0; maxerr_basedg = 0.0; n = 0; nbad = 0
tree = Dict{String,Float32}()
function finalize_tree(sp, temel, temslp, temsasp, ssite)
    isempty(tree) && return
    d = tree["D"]; d <= 0 && (empty!(tree); return)
    d2 = d*d
    pbal = tree["PBAL"]; prd = tree["PRD"]; cr = tree["CR"]
    dgcomp1 = DGEL[sp]*temel + DGSLOP[sp]*temslp + DGSASP[sp]*temsasp + DGLNSI[sp]*log(ssite)
    dgcomp2 = DGDISQ[sp]*d2 + DGLD[sp]*log(d) + DGDBAL[sp]*pbal + DGRD[sp]*prd + DGLNCR[sp]*log(cr)
    basedg = exp(DGCONB1[sp] + dgcomp2 + dgcomp1)
    # permafrost: for non-perm species PFMOD=1 (all akt01 species). Perm path deferred.
    pfmod = 1.0f0
    dgpred = YR * basedg * pfmod * dgmult(sp)
    brat = akbrat(sp, d)                       # module ak_bratio (validates AK bark path too)
    abs(brat - tree["BRAT"]) > 1f-6 && println("  BRAT drift sp=$sp D=$d jl=$brat ora=$(tree["BRAT"])")
    tempd1 = d * brat
    dup = d + dgpred
    tempd2 = dup * brat
    cor = tree["COR"]; dgcon = tree["DGCON"]
    dds = log(tempd2^2 - tempd1^2) + cor + dgcon
    dds < -9.21f0 && (dds = -9.21f0)
    # compare
    o_basedg = tree["OBASEDG"]; o_wk2 = tree["OWK2"]
    e_b = abs(basedg - o_basedg) / max(abs(o_basedg), 1f-30)
    e_d = abs(dds - o_wk2) / max(abs(o_wk2), 1f-9)
    maxerr_basedg = max(maxerr_basedg, e_b); maxerr_dds = max(maxerr_dds, e_d)
    n += 1
    if e_d > 1e-5 || e_b > 1e-5
        nbad += 1
        println("  MISMATCH sp=$sp D=$d  DDS jl=$dds ora=$o_wk2 (rel $e_d)  BASEDG jl=$basedg ora=$o_basedg (rel $e_b)")
    end
    empty!(tree)
end

for ln in lines
    if occursin("TOP OF SPECIES LOOP", ln)
        m = match(r"ISPC=\s*(\d+)", ln); m !== nothing && (cur_ispc = parse(Int, m.captures[1]))
    elseif occursin("SPECIES VARIABLES", ln)
        cur_temel = grabf(ln, "TEMEL="); cur_temslp = grabf(ln, "TEMSLP=")
        cur_temsasp = grabf(ln, "TEMSASP="); cur_ssite = grabf(ln, "SSITE=")
    elseif occursin("TOP OF TREE LOOP", ln)
        empty!(tree)
        tree["D"] = grabf(ln, "D=")
    elseif occursin("TREE VARIABLES", ln)
        tree["PBAL"] = grabf(ln, "PBAL="); tree["PRD"] = grabf(ln, "PRD=")
        tree["CR"] = grabf(ln, "CR=")
    elseif occursin("BASEDG EQN", ln)
        tree["OBASEDG"] = grabf(ln, "BASEDG=")
    elseif occursin("END OF TREE LOOP", ln)
        tree["BRAT"] = grabf(ln, "BRAT="); tree["COR"] = grabf(ln, "COR=")
        tree["DGCON"] = grabf(ln, "DGCON(ISPC)")
    elseif occursin("WK2(I)=", ln)
        tree["OWK2"] = grabf(ln, "WK2(I)=")
        finalize_tree(cur_ispc, cur_temel, cur_temslp, cur_temsasp, cur_ssite)
    end
end

println("AK DGF cyc0 validation: $n trees checked, $nbad mismatches")
println("  max rel err BASEDG = $maxerr_basedg")
println("  max rel err WK2/DDS = $maxerr_dds")
println(nbad == 0 ? "  VERDICT: BIT-EXACT-OR-ULP (ported DGF reproduces the live oracle)" : "  VERDICT: MISMATCH")
end

# Usage: julia --project=. tools/southeastalaska/validate_dgf_cyc0.jl <path-to-FVSak-DEBUG-DGF.out>
# Regenerate the dump:  /workspace/.akwork/FVSak_clean --keywordfile=akdbg.key   (akdbg.tre = a copy
# of tests/FVSak/akt01.tre); the `DEBUG / DGF / END` block emits the per-tree DGF dump to fort.16/.out.
dumppath = isempty(ARGS) ? joinpath(@__DIR__, "akdbg.out") : ARGS[1]
main(readlines(dumppath))
