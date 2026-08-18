# Load the DELIVERABLE lpmpb.jl with minimal package stubs and re-validate its
# ACTUAL functions against the FVSie_lpmpb goldens (lp_dbg.stderr).
abstract type AbstractMpbState end
abstract type AbstractVariant end
struct StandState end
struct KeywordReader end
for v in (:CentralRockies,:BlueMountains,:CentralIdaho,:EastCascades,:EasternMontana,
          :InlandEmpire,:SouthCentralOregon,:Teton,:Utah)
    @eval struct $v <: AbstractVariant end
end
const KW_EOF = 0; const KW_STOP = -1
read_keyword!(kr) = nothing
cycle_year_at(c, k) = 0

include("lpmpb.jl")

hex2f32(s) = reinterpret(Float32, parse(UInt32, s; base=16))
f32hex(x::Float32) = uppercase(string(reinterpret(UInt32, x); base=16, pad=8))

recs = Tuple{Int,Int,Float32,Float32,Float32}[]
gs = zeros(Float32,10); gg = zeros(Float32,10); gp = zeros(Float32,10); numyrs = 0
for ln in eachline("lp_dbg.stderr")
    if startswith(ln,"DBGCM_HDR"); global numyrs = parse(Int, match(r"NUMYRS\s+(\d+)",ln)[1])
    elseif startswith(ln,"DBGCM_CLS"); t=split(ln); c=parse(Int,t[2]); gs[c]=hex2f32(t[4]); gg[c]=hex2f32(t[6]); gp[c]=hex2f32(t[8])
    elseif startswith(ln,"DBGCM_REC"); t=split(ln); push!(recs,(parse(Int,t[3]),parse(Int,t[5]),hex2f32(t[7]),hex2f32(t[9]),hex2f32(t[11]))) end
end

lp_dbh = Float32[r[3] for r in recs]; lp_tpa = Float32[r[4] for r in recs]
start = mpb_coldbh_start(lp_dbh, lp_tpa)
green = mpb_colmod(start, numyrs)
prkill = mpb_prkill(start, @view green[numyrs,:])

ok(a,b) = f32hex(a)==f32hex(b)
cs = count(i->ok(start[i],gs[i]),1:10)
cg = count(i->ok(green[numyrs,i],gg[i]),1:10)
cp = count(i->ok(prkill[i],gp[i]),1:10)
ci = count(r->mpb_colind(r[3])==r[2], recs)
cx = count(r->ok(prkill[mpb_colind(r[3])]*r[4], r[5]), recs)
# MPBER minimum condition (BA arbitrary>0 large so pbalpp finite): noer must be true (host clears)
er = mpb_er(lp_dbh, lp_tpa, 300.0f0)
println("lpmpb.jl functions vs FVSie_lpmpb goldens:")
println("  COLDBH START      $cs/10")
println("  COLMOD GREEN(NY)  $cg/10")
println("  COLMRT PRKILL     $cp/10")
println("  COLIND            $ci/$(length(recs))")
println("  XT=PRKILL*PROB    $cx/$(length(recs))")
println("  MPBER noer=$(er.noer) (host clears minimum condition)")
println("  mpb_idxlp(InlandEmpire())=$(mpb_idxlp(InlandEmpire())) (expect 7); CR=$(mpb_idxlp(CentralRockies())) (expect 11)")
@assert cs==10 && cg==10 && cp==10 && ci==length(recs) && cx==length(recs) && er.noer "MISMATCH"
println("ALL BIT-EXACT ✓")
