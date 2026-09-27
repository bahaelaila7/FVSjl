# test_forkod_default.jl — FORKOD KODFOR → IFOR for the western Prognosis variants, checked against a table
# extracted from each variant's forkod.f (test/fixtures/forkod/expected_ifor.tsv, 346 codes).
#
# Two bugs this guards:
#  * "forest code not found": forkod.f's .NOT.FORFOUND path calls ERRGRO(3) and clears USEIGL but never assigns
#    IFOR, so IFOR keeps grinit.f's default. jl forced IFOR=1 (TT/UT/BM/PN/WC/CA/OP), so an FIA stand whose
#    LOCATION is outside JFOR (TT 11790600010690, LOCATION 404) ran on the wrong forest's DGFOR/MAPLOC.
#  * BIA reservation pseudo-codes (first SELECT CASE (KODFOR)) were missing for CA/OC (18) and SO (3).
using FVSjl, Test

const _FKD_F = Dict("tt"=>FVSjl.tt_forkod!, "ut"=>FVSjl.ut_forkod!, "bm"=>FVSjl.bm_forkod!,
    "pn"=>FVSjl.pn_forkod!, "wc"=>FVSjl.wc_forkod!, "ca"=>FVSjl.ca_forkod!, "op"=>FVSjl.op_forkod!,
    "ec"=>FVSjl.ec_forkod!, "em"=>FVSjl.em_forkod!, "ci"=>FVSjl.ci_forkod!, "ie"=>FVSjl.ie_forkod!,
    "nc"=>FVSjl.nc_forkod!, "so"=>FVSjl.so_forkod!, "ak"=>FVSjl.ak_forkod!, "ws"=>FVSjl.ws_forkod!)

@testset "FORKOD KODFOR → IFOR matches forkod.f" begin
    for l in eachline(joinpath(@__DIR__, "..", "fixtures", "forkod", "expected_ifor.tsv"))
        startswith(l, '#') && continue
        v, k, e = split(l, '\t'); k = parse(Int, k); e = parse(Int, e)
        got = if v == "oc"
            FVSjl.oc_forkod(k)
        else
            p = FVSjl.PlotData(); p.user_forest_code = Int32(k); p.forest_idx = Int32(0)
            _FKD_F[v](p); Int(p.forest_idx)
        end
        @test (v, k, got) == (v, k, e)
    end
end

@testset "FORKOD not-found leaves IGL untouched" begin
    for v in ("tt", "ut", "bm", "pn", "wc", "ca", "op")
        p = FVSjl.PlotData(); p.user_forest_code = Int32(999); p.forest_idx = Int32(0); p.geo_location = Int32(-7)
        _FKD_F[v](p)
        @test (v, p.geo_location) == (v, -7)      # USEIGL = .FALSE.
    end
end
