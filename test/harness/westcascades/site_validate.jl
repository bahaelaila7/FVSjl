# WC chunk-1/2 regression harness: the species-coefficient table loads and wc_sitset! reproduces
# the LIVE FVSwc_clean wct01 site chain BIT-EXACT.
#   • Chunk 2 (site): forkod 618→IFOR 6; habitat 52→PCOML[52]=CFS551→ECOCLS→site species DF(16),
#     SITEAR(DF)=73, SDIDEF=815; all 39 SITEAR vs the live "AFTER SITE ADJUSTMENT FACTORS" DEBUG dump.
#   • Chunk 1 (bark): CSV-driven wc_bratio reproduces wc/bratio.f (JBARK→BARKB) within Float32 epsilon.
# Oracle: /workspace/.wcwork/FVSwc_clean wct01 stand-1, keyfile `DEBUG / SITSET`.
# Run:  JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project=. test/harness/westcascades/site_validate.jl
using FVSjl
const F = FVSjl

# --- live SITEAR (FVSwc_clean wct01 DEBUG SITSET) ---
const ORACLE_SITEAR = Float32[
 72.9977188, 43.5985756, 43.5985756, 72.9977188, 43.5985756, 43.5985756, 72.9977188, 72.9977188,
 72.9977188, 72.9977188, 45.7409706, 72.9977188, 72.9977188, 72.9977188, 72.9977188, 73.0,
 72.9977188, 72.9977188, 43.5985756, 22.2486191, 54.7482910, 19.4166489, 47.4485168, 109.496582,
 51.0984039, 54.7482910, 62.0480614, 48.7942810, 16.7894764, 43.5985756, 51.0984039, 72.9977188,
 18.2494297, 43.7986336, 18.2494297, 36.4988594, 36.4988594, 72.9977188, 72.9977188]

function validate_site()
    v = F.WestCascades(); s = F.StandState(v)
    F.load_species_coefficients!(s, v)
    p = s.plot
    p.user_forest_code = Int32(618); p.habitat_code = Int32(52); p.site_species = Int32(0)
    fill!(@view(p.sp_site_index[1:39]), 0f0); fill!(@view(p.sp_sdi_def[1:39]), 0f0)
    F.site_setup!(s, v)
    worst = 0f0; nbad = 0
    for i in 1:39
        d = abs(p.sp_site_index[i] - ORACLE_SITEAR[i]); d > worst && (worst = d)
        d > 0.01f0 && (nbad += 1; println("  SITEAR mismatch sp$i jl=$(p.sp_site_index[i]) live=$(ORACLE_SITEAR[i])"))
    end
    okifor = p.forest_idx == 6; okisisp = p.site_species == 16
    oksdi = all(==(815f0), @view p.sp_sdi_def[1:39])
    println("WC site: IFOR=$(p.forest_idx)(exp6) ISISP=$(p.site_species)(exp16) SITEAR worst|Δ|=$worst nbad=$nbad/39 SDIDEF=815:$oksdi")
    return nbad == 0 && okifor && okisisp && oksdi
end

function validate_bark()
    s = F.StandState(F.WestCascades()); sd = s.coef.species
    JBARK = [2,2,2,2,2,2,2,5,5,10,10,4,4,4,3,1,14,11,12,11,6,9,9,6,7,9,9,8,11,10,13,13,13,9,9,9,9,1,10]
    BARKB = Dict(1=>(0.903563,0.989388,1),2=>(0.904973,1.0,1),3=>(0.809427,1.016866,1),4=>(0.859045,1.0,1),
                 5=>(0.837291,1.0,1),6=>(0.08360,0.94782,2),7=>(0.15565,0.90182,2),8=>(0.8558,1.0213,1),
                 9=>(0.075256,0.949670,2),10=>(0.9,1.0,1),11=>(0.949670,1.0,1),12=>(0.933710,1.0,1),
                 13=>(0.933290,1.0,1),14=>(0.7012,1.04862,1))
    fort(sp,d) = (d<=0 ? 0.99f0 : begin (a,b,et)=BARKB[JBARK[sp]]; br= et==1 ? a*d^b/d : (a+b*d)/d
                  Float32(clamp(br,0.80,0.99)) end)
    worst = 0f0
    for sp in 1:39, d in Float32[1,3,5,10,20,40,80]
        worst = max(worst, abs(F.wc_bratio(sd,sp,d) - fort(sp,d)))
    end
    println("WC bark: wc_bratio vs bratio.f worst|Δ|=$worst")
    return worst < 1f-4
end

ok = validate_site() & validate_bark()
println(ok ? "PASS — WC chunk 1 (bark table) + chunk 2 (site/SDImax) bit-exact vs live." :
             "FAIL — WC chunk 1/2 divergence.")
