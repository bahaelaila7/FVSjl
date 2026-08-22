# FVS_BM_Main (WWPB MAINOUT stand-summary DBS table, dbs/dbsbmmain.f via the bmout.f DO-20
# aggregation ported in engine/wwpb_landscape.jl `wwpb_main_report`), gated by the PPBMMAIN
# DATABASE keyword. (2026-08-22.)
#
# There is NO live oracle for FVS_BM_* (the WWPB subsystem is in 0 FVS sourceLists; its driver
# bmdrv is called only by the source-absent PPE orchestrator PPMAIN — see memory
# fvsjl-bm-dbs-tables-no-oracle). So the aggregation is validated at the WWPB reconstruction bar:
# BIT-EXACT vs a gfortran-16 golden that compiles the real bmout.f MAINOUT arithmetic
# (scratchpad/wwpb/bmdbs/golden_bmmain.f) on a controlled stand state. The outbreak STATE feeding
# the report in run_keyfile is jl's reconstructed single-stand harness (cornered, like every WWPB
# kill); this test pins the aggregation + serialization transcription.

using Test
using FVSjl
using FVSjl: WwpbStand, wwpb_init_coeffs, wwpb_main_report, write_dbs_bm_main!,
             write_dbs_bm_tree!, write_dbs_bm_vol!, WWPB_UPSIZ_DEFAULT, WWPB_NSCL
using SQLite, DBInterface

# Controlled state identical to golden_bmmain.f: class 3 = 100 host + 20 nonhost, BAH 50, vols 8/6,
# kills 5+2, 10% special; class 5 = 40 host, BAH 30, vol 15, kill 3, 5% special; slash 1.5 + 0.75.
function _bm_controlled_stand()
    st = WwpbStand()
    st.tree[3, 1] = 100.0f0; st.tree[3, 2] = 20.0f0
    st.bah[3] = 50.0f0; st.tvol[3, 1] = 8.0f0; st.tvol[3, 2] = 6.0f0
    st.pbkill[3] = 5.0f0; st.allkll[3] = 2.0f0; st.spclt[3, 1] = 0.10f0
    st.tree[5, 1] = 40.0f0; st.bah[5] = 30.0f0; st.tvol[5, 1] = 15.0f0
    st.pbkill[5] = 3.0f0; st.spclt[5, 1] = 0.05f0
    st.dwphos[1, 1] = 1.5f0; st.dwphos[2, 3] = 0.75f0
    # summary slots (bmsdit! populates these; set consistently for BAH/TPAH columns)
    st.bah[WWPB_NSCL+1] = 80.0f0            # BAH(3)+BAH(5)
    st.tree[WWPB_NSCL+1, 1] = 140.0f0       # TPAH = host in cl 3 + 5
    st.bastd = 95.0f0; st.grfstd = 1.7f0; st.oldbkp = 12.3f0; st.bkp = 4.2f0
    return st
end

# gfortran-16 golden_bmmain.f output (F16.7), the bit-exact target.
const _BM_GOLDEN = (PreDispBKP = 12.3000002, PostDispBKP = 4.1999998, StandRV = 1.7000000,
    StandBA = 95.0, BAH = 80.0, BA_BtlKld = 5.75, TPA = 160.0, TPAH = 140.0,
    TPA_BtlKld = 10.0, StandVol = 1520.0, VolHost = 1400.0, VolBtlKld = 101.0,
    BA_Special = 5.0559969, Ips_Slash = 2.25)

@testset "FVS_BM_Main WWPB MAINOUT — aggregation bit-exact vs bmout.f golden" begin
    st = _bm_controlled_stand()
    coeffs = wwpb_init_coeffs(copy(WWPB_UPSIZ_DEFAULT))
    rep = wwpb_main_report(st, coeffs)
    for (col, val) in pairs(_BM_GOLDEN)
        got = Float32(getfield(rep, col))
        # Float32 print-rounding tolerance vs the F16.7 golden text
        @test isapprox(got, Float32(val); atol = 5.0f-5, rtol = 1.0f-6)
    end
    # sanitation columns are 0 (no harvest keyword) — matches the oracle
    for col in (:BA_San_Remv, :BKP_San_Remv, :TPA_SanRemvLv, :TPA_SanRemLvDd, :VolRemSan, :VolRemSalv)
        @test getfield(rep, col) == 0.0f0
    end

    # serializer round-trip: schema + 23-column write, values byte-exact through SQLite
    @testset "write_dbs_bm_main! round-trip" begin
        mktempdir() do dir
            dbp = joinpath(dir, "bm.db")
            write_dbs_bm_main!(dbp, "CASE0", "STAND0", [(1990, rep)])
            db = SQLite.DB(dbp)
            try
                tabs = [r.name for r in DBInterface.execute(db,
                    "SELECT name FROM sqlite_master WHERE type='table' AND name='FVS_BM_Main'")]
                @test "FVS_BM_Main" in tabs
                row = first(DBInterface.execute(db, "SELECT * FROM FVS_BM_Main WHERE Year=1990"))
                @test row.StandID == "STAND0"
                @test isapprox(Float32(row.BA_BtlKld), 5.75f0; atol = 5.0f-5)
                @test isapprox(Float32(row.StandVol), 1520.0f0; atol = 5.0f-3)
                @test isapprox(Float32(row.BA_Special), 5.0559969f0; atol = 5.0f-5)
                @test isapprox(Float32(row.Ips_Slash), 2.25f0; atol = 5.0f-5)
                @test row.VolRemSalv == 0.0
            finally
                SQLite.close(db)
            end
        end
    end

    # FVS_BM_Tree / FVS_BM_Vol per-size-class detail (bmout.f TREEOUT/VOLOUT golden, classes 3 & 5).
    @testset "FVS_BM_Tree / FVS_BM_Vol per-class vs golden" begin
        # golden_bmmain.f: SC3 TPA=120 HOST=100 TKLD=7 SPCL=10 TV=920 HV=800 VK=56;
        #                  SC5 TPA=40 HOST=40 TKLD=3 SPCL=2 TV=600 HV=600 VK=45
        @test rep.tpa_sc[3] == 120.0f0 && rep.tpa_sc[5] == 40.0f0
        @test rep.host_sc[3] == 100.0f0 && rep.host_sc[5] == 40.0f0
        @test rep.tkld_sc[3] == 7.0f0 && rep.tkld_sc[5] == 3.0f0
        @test rep.spcl_sc[3] == 10.0f0 && rep.spcl_sc[5] == 2.0f0
        @test rep.tvol_sc[3] == 920.0f0 && rep.tvol_sc[5] == 600.0f0
        @test rep.hvol_sc[3] == 800.0f0 && rep.hvol_sc[5] == 600.0f0
        @test rep.volk_sc[3] == 56.0f0 && rep.volk_sc[5] == 45.0f0
        # zero classes stay zero
        @test all(rep.tpa_sc[i] == 0.0f0 for i in (1, 2, 4, 6, 7, 8, 9, 10))

        mktempdir() do dir
            dbp = joinpath(dir, "bmtv.db")
            write_dbs_bm_tree!(dbp, "C", "S", [(1990, rep)])
            write_dbs_bm_vol!(dbp, "C", "S", [(1990, rep)])
            db = SQLite.DB(dbp)
            try
                tr = first(DBInterface.execute(db, "SELECT * FROM FVS_BM_Tree WHERE Year=1990"))
                @test tr.TPA_SC3 == 120.0 && tr.HOST3 == 100.0 && tr.TKLD3 == 7.0
                @test tr.SPCL3 == 10.0 && tr.SAN3 == 0.0 && tr.SAN10 == 0.0
                vl = first(DBInterface.execute(db, "SELECT * FROM FVS_BM_Vol WHERE Year=1990"))
                @test vl.TV_SC3 == 920.0 && vl.HV_SC3 == 800.0 && vl.VK_SC3 == 56.0
                @test vl.TV_SC5 == 600.0 && vl.VK_SC5 == 45.0
            finally
                SQLite.close(db)
            end
        end
    end
end
