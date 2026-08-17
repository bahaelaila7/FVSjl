# test_ontario.jl — ON (Ontario) variant beachhead: Penner large-tree diameter growth
# (canada/on/dgf.f) dump-replay BIT-EXACT vs the relinked/rebuilt oracle FVSon_g16.
#
# ON uses the Penner (2006) annual diameter-increment growth family (NOT Wykoff ln(DDS)).
# Golden = ontario_penner_dgf_dump.txt: per-tree g16 dump of the Penner inputs+outputs,
# 10 Float32-hex cols after I/ISPC/KSP/AGS/ICYC:
#   DIAM(in) DBHM_final(cm) SIM BAM QMDM BALM HTM BARK DIAGR DDS
# on_penner_dds(ksp, ags, diam_in, sim, bam, qmdm, balm, htm, bark) -> (dbhm_final, diagr, dds).
# Transcendentals go through glibc logf/expf (native/fmath drift ~1 ULP over the 10-yr loop).

using Test

const _ON_HERE = @__DIR__
_on_h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))

@testset "ON — Penner large-tree DGF dump-replay (bit-exact vs FVSon_g16)" begin
    seen = Set{Tuple{Int,Int}}()
    nrec = okB = okG = okD = 0
    for ln in eachline(joinpath(_ON_HERE, "ontario_penner_dgf_dump.txt"))
        f = split(strip(ln)); isempty(f) && continue
        I, ISPC, KSP, AGS, _ICYC = parse.(Int, f[1:5])
        (I, ISPC) in seen && continue; push!(seen, (I, ISPC))   # dedup repeated calls
        diam, dbhm_o, sim, bam, qmdm, balm, htm, bark, diagr_o, dds_o = _on_h2f.(f[6:15])
        dbhm_j, diagr_j, dds_j = FVSjl.on_penner_dds(KSP, AGS, diam, sim, bam, qmdm, balm, htm, bark)
        nrec += 1
        okB += reinterpret(UInt32, dbhm_j)  == reinterpret(UInt32, dbhm_o)
        okG += reinterpret(UInt32, diagr_j) == reinterpret(UInt32, diagr_o)
        okD += reinterpret(UInt32, dds_j)   == reinterpret(UInt32, dds_o)
        @test reinterpret(UInt32, dds_j) == reinterpret(UInt32, dds_o)
    end
    @test nrec > 0
    @test (okB, okG, okD) == (nrec, nrec, nrec)   # DBHM_final / DIAGR / DDS all bit-exact
    # registry: the ON code resolves to the Ontario singleton (72 species)
    @test FVSjl.variant_from_code("ON") isa FVSjl.Ontario
end
