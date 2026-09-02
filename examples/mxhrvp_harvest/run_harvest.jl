# MXHRVP multistand harvest scheduling — worked demo.
#
# Run:  JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project=. examples/mxhrvp_harvest/run_harvest.jl
#
# MXHRVP is a LANDSCAPE feature: instead of treating each stand on its own schedule, you
# define a landscape harvest POLICY (a target resource flow + a per-stand priority + a
# per-stand credit) and the scheduler picks WHICH stands to cut each master cycle to meet
# the target, in priority order. The entry point is `ppe_run_landscape_harvest!` (the
# hvaloc.f coordinator); it returns, per master cycle, the HVSEL selection table.
using FVSjl: ppe_run_landscape_harvest!, PPEStand, EastCascades, hvsel!

const KF = joinpath(@__DIR__, "member_stand.key")
stands3() = [PPEStand(KF; area = 11.0) for _ in 1:3]      # a 3-stand landscape, 11 ac each

showtable(title, res) = begin
    println("\n== $title ==")
    for r in res
        sel = join((s ? "YES" : "no " for s in r.selected), " ")
        println("  $(r.year): target=$(round(r.target,digits=1))  selected=[$sel]  " *
                "resource=$(round(r.selected_resource,digits=1))  " *
                "%oftarget=$(round(r.pct_of_target,digits=1))  hvpart=$(round(r.hvpart,digits=3))")
    end
end

# (a) BASIC landscape target-flow: harvest toward a 1000-unit BA flow, priority & credit = BBA.
res_a = ppe_run_landscape_harvest!(stands3(); variant = EastCascades(),
            labels = ["ALL","ALL","ALL"], mslabel = "ALL",
            target_expr = "1000", priority_expr = "BBA", credit_expr = "BBA",
            master_years = [1990, 2000, 2010])
showtable("(a) basic target-flow (TARGET=1000, PRIORITY=BBA, CREDIT=BBA)", res_a)

# (b) EXACT partial cut: a differential credit (only a SELECTED stand contributes) + a small
#     target forces the marginal stand to be PARTIALLY selected (status 4, HVPART fraction).
res_b = ppe_run_landscape_harvest!(stands3(); variant = EastCascades(),
            labels = ["ALL","ALL","ALL"], mslabel = "ALL",
            target_expr = "150", priority_expr = "BBA", credit_expr = "SELECTED * BBA",
            master_years = [1990, 2000, 2010], lprtct = true)
showtable("(b) EXACT partial cut (CREDIT=SELECTED*BBA, TARGET=150, lprtct=true)", res_b)
println("     → note hvpart > 0 on the cycle where a stand is only partially needed.")

# (c) MXCLRCUT max-contiguous-clearcut constraint (the LHVMXC kernel). Shown at the hvsel!
#     level: three would-be-clearcut stands (status -1), 50 ac each, in a STAR (stand 1
#     borders 2 and 3). Target wants all three, but a 60-acre contiguous-clearcut cap lets
#     only the first (isolated 50 ac) be cut — its two neighbors would create 100 contiguous
#     acres (> 60) and are vetoed.
border = (i, j) -> (Set([(1,2),(1,3),(2,1),(3,1)]) |> s -> ((i,j) in s ? 1.0f0 : -1.0f0))
areas  = Float32[50, 50, 50]
st_unc = Int[-1, -1, -1]
unc, = hvsel!(st_unc, Float32[9,6,3], Float32[50,50,50], Float32[0,0,0], 200f0)
st_mxc = Int[-1, -1, -1]
mxc, = hvsel!(st_mxc, Float32[9,6,3], Float32[50,50,50], Float32[0,0,0], 200f0;
              lhvmxc = true, hvmxcc = 60f0, areas = areas, border = border)
println("\n== (c) MXCLRCUT max-contiguous-clearcut (HVMXCC=60 ac) ==")
println("  unconstrained selection status: $unc   (3 = cut)")
println("  with 60-ac contiguous cap:       $mxc   (2 = vetoed by the contiguous-clearcut cap)")

println("\nDone — MXHRVP harvest scheduling ran on the 3-stand landscape.")
