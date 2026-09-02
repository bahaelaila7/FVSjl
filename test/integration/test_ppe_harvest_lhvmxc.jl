# PPE MXHRVP LHVMXC (max-contiguous-clearcut) — SPCNTG contiguous-acreage bit-exact vs
# gfortran-16 driver-golden (scratchpad/ppe/mxhrvp/driver_spcntg.f, the real spcntg.f with
# stubbed SPLAAR/SPNBBD over a controlled area array + adjacency), plus the hvsel! LHVMXC veto.
using Test
using FVSjl: spcntg!, hvcntg, hvsel!

# border closure from an adjacency Dict{(i,j)=>border}; missing pair ⇒ -1 (not neighbors).
mkborder(adj) = (i, j) -> Float32(get(adj, (min(i,j), max(i,j)), -1.0))

# (label, areas, adjacency, isdwk1, expected contig hex, expected nneig, expected order)
const SPCNTG_GOLDENS = [
    ("A chain 1-2-3",        Float32[10,20,30],    Dict((1,2)=>5.0,(2,3)=>5.0),
        [1,2,3], 0x42700000, 3, [1,2,3]),                       # 60.0
    ("B 1-2 nbr, 3 isolated", Float32[10,20,30],   Dict((1,2)=>5.0),
        [1,2,3], 0x41f00000, 2, [1,2,3]),                       # 30.0
    ("C star 1-{2,3,4}",     Float32[10,20,30,40], Dict((1,2)=>1.0,(1,3)=>1.0,(1,4)=>1.0),
        [1,2,3,4], 0x42c80000, 4, [1,2,3,4]),                   # 100.0
    ("D single stand",       Float32[15],          Dict{Tuple{Int,Int},Float64}(),
        [1], 0x41700000, 1, [1]),                               # 15.0
    ("E two components",     Float32[10,20,30,40], Dict((1,2)=>5.0,(3,4)=>5.0),
        [1,2,3,4], 0x41f00000, 2, [1,2,3,4]),                   # 30.0 (only {1,2})
]

@testset "PPE MXHRVP LHVMXC — SPCNTG contiguous acres bit-exact vs driver-golden" begin
    for (label, areas, adj, isdwk1, exp_contig, exp_nneig, exp_order) in SPCNTG_GOLDENS
        work = copy(isdwk1)
        contig, nneig = spcntg!(work, areas, mkborder(adj))
        @testset "$label" begin
            @test reinterpret(UInt32, contig) == exp_contig
            @test nneig == exp_nneig
            @test work == exp_order
        end
    end
end

@testset "PPE MXHRVP LHVMXC — hvsel! max-contiguous-clearcut veto" begin
    # Three would-be-clearcut candidates (status -1), each yield 50, each 50 acres, in a
    # STAR — stand 1 borders both 2 and 3. Target 200 wants all three, but HVMXCC=60 caps a
    # contiguous clearcut to ~one stand: the 1st selects (contig 50 ≤ 60); the 2nd and 3rd
    # each border the selected stand 1, so selecting either makes 100 contiguous acres
    # (> 60) → both vetoed. Priorities descending so order is 1,2,3.
    areas = Float32[50, 50, 50]
    border = mkborder(Dict((1,2)=>1.0, (1,3)=>1.0))
    prio = Float32[9, 6, 3]
    ysel = Float32[50, 50, 50]
    ynot = Float32[0, 0, 0]

    # WITHOUT the constraint: greedy fills the target → stands 1 & 2 selected (100 ≥ target? no,
    # target 200 needs 1,2,3 fully) — all three selected (status 3).
    st0 = Int[-1, -1, -1]
    s0, = hvsel!(st0, copy(prio), copy(ysel), copy(ynot), 200f0)
    @test s0 == [3, 3, 3]                       # unconstrained: all three cut

    # WITH LHVMXC (hvmxcc=60): stand 1 cut, 2 & 3 vetoed (would exceed 60 contiguous acres).
    st1 = Int[-1, -1, -1]
    s1, = hvsel!(st1, copy(prio), copy(ysel), copy(ynot), 200f0;
                 lhvmxc = true, hvmxcc = 60f0, areas = areas, border = border)
    @test s1[1] == 3                            # first (highest priority) selected
    @test s1[2] == 2 && s1[3] == 2             # neighbors vetoed by the contiguous-acre cap

    # A non-clearcut candidate (+1) is never vetoed by LHVMXC (hvcntg returns 0).
    st2 = Int[1, 1]
    s2, = hvsel!(st2, Float32[9,6], Float32[50,50], Float32[0,0], 40f0;
                 lhvmxc = true, hvmxcc = 1f0, areas = Float32[50,50],
                 border = mkborder(Dict((1,2)=>1.0)))
    @test s2[1] == 3                            # +1 candidate selected despite tiny hvmxcc
end
