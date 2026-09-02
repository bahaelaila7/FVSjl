# PPE MXHRVP IHVEXT=1 external-I/O bridge — the PPE_FFERdAccess.txt read-back (SPRDRD/SPRDIS)
# ported as ppe_read_rdaccess, bit-exact vs the gfortran-16 driver-golden
# (scratchpad/ppe/mxhrvp/driver_sprd.f, the real sprdrd.f/sprdis.f). The read-back file is
# the USER-approved staged read: the priorities the external selector would emit, read
# identically by the oracle and FVSjl → the priority-override path is validated end-to-end
# (the selection under those priorities is the already-validated IHVEXT selection core).
using Test
using FVSjl: ppe_read_rdaccess, ppe_run_landscape_harvest!, PPEStand, EastCascades

const RDA = joinpath(@__DIR__, "..", "fixtures", "ppe", "mxhrvp", "rdaccess.txt")

@testset "IHVEXT bridge — PPE_FFERdAccess.txt reader bit-exact vs SPRDRD/SPRDIS golden" begin
    d = ppe_read_rdaccess(RDA)
    @test reinterpret(UInt32, d["STANDA"]) == 0x42480000   # 50.0
    @test reinterpret(UInt32, d["STANDB"]) == 0x41f00000   # 30.0
    @test reinterpret(UInt32, d["STANDC"]) == 0x41200000   # 10.0
    @test length(d) == 3                                    # '-999' terminated the read
end

@testset "IHVEXT bridge — a blank value field maps to the -99999 sentinel" begin
    mktemp() do path, io
        write(io, "STANDA                       50.0\n")   # value present
        write(io, "STANDB                           \n")   # value blank → -99999.0
        write(io, "-999\n")
        close(io)
        d = ppe_read_rdaccess(path)
        @test d["STANDA"] == 50.0f0
        @test d["STANDB"] == -99999.0f0
    end
end

@testset "IHVEXT bridge — coordinator applies the staged priority override" begin
    kf = joinpath(@__DIR__, "..", "fixtures", "ppe", "mxhrvp", "stand.key")
    st = PPEStand(kf; area = 11.0)
    # stage a road-access file keyed to the member stand's id with a positive priority,
    # then run the coordinator in external-selection mode — the priorities are overridden
    # (become the staged value) and HVSEL runs in IHVEXT mode.
    mktemp() do path, io
        write(io, "S248112                      42.0\n-999\n")
        close(io)
        res = ppe_run_landscape_harvest!([st, st, st]; variant = EastCascades(),
                  labels = ["ALL","ALL","ALL"], mslabel = "ALL",
                  target_expr = "1000", priority_expr = "BBA", credit_expr = "BBA",
                  master_years = [1990], rdaccess = path)
        r = res[1]
        @test all(p == 42.0f0 for p in r.priority)    # priority overridden by the staged value
        @test all(r.selected)                          # positive priority → candidates selected
    end
end
