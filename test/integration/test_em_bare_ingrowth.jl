# Bare-stand AUTOES ingrowth (EM / EZCRUISE INADV=1) — regression for the dominant western-sweep bug.
#
# The Eastern-Montana FIA population is dominated by BARE plots (0 tree records). FVS auto-invokes the
# NOTREES/EZCRUISE regeneration option (initre.f:280-289 → esinit.f:79 ESEZCR), which sets INADV=1 for the
# WHOLE run. INADV=1 makes estab.f SKIP the ESB inventory-stocking calibration in every cycle (estab.f:319/511
# `IF(INADV.EQ.1.OR.NTALLY.NE.1) GO TO ...`), so a bare-origin ingrowth tally uses PROB1 = logistic(PN) with
# NO ESB-ESB1 shift.
#
# jl formerly applied the esb-esb1 calibration shift unconditionally (its guard checked only the ≤20-yr window,
# not INADV), collapsing PROB1 and UNDER-producing ingrowth ~2.9×. MEASURED live FVSem_clean on the traced bare
# stand 85271557010661 (0 tree records, IHAB=3, IFO=11, SLO=0.05, ELEV=21.9): PN=-0.5822, ESB=0, ESB1=0 ⇒
# PROB1=0.3584 ⇒ 221.5 INGROWTH TREES/ACRE (DF 198 + LP 24) in cycle 1; the ESTPP per-plot ITPP sequence is
# already bit-exact (Σ=103 trees at PROB=2.151). Pre-fix jl gave PROB1=0.1222 ⇒ 75.5 TPA.
#
# Signed 17-stand bare-EM sign-tally (cycle-1 ingrowth TPA, jl−live): BEFORE = 17 UNDER / 0 OVER / 0 TIE
# (one-directional ⇒ real bug); AFTER = 0/0/17 TIE (bit-exact). Gate test_multicycle stays 339/11 byte-identical.
#
# This test is DB-gated (needs the FIA source DB); it is skipped where the DB is absent (CI without /workspace).

using Test
using FVSjl

const _EM_FIA_DB = "/workspace/SQLite_FIADB_ENTIRE.db"

@testset "EM bare-stand AUTOES ingrowth (EZCRUISE INADV=1, no ESB shift)" begin
    if !isfile(_EM_FIA_DB)
        @info "FIA source DB not present ($_EM_FIA_DB) — skipping bare-EM ingrowth regression"
        @test_skip true
    else
        # Traced bare stand: 0 tree records ⇒ EZCRUISE. Live FVSem_clean cyc-1 ingrowth TPA = 222.
        keytext(cn) = """
        STDIDENT
        $cn
        DATABASE
        DSNin
        $_EM_FIA_DB
        StandSQL
        SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
        EndSQL
        TreeSQL
        SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
        EndSQL
        END
        NUMCYCLE         3
        ECHOSUM
        PROCESS
        STOP
        """
        # cyc-1 TPA parsed from the .sum first projection row (cols 9:14).
        function cyc1_tpa(sumtext)
            rows = Int[]
            for ln in split(sumtext, '\n')
                length(ln) < 14 && continue
                y = tryparse(Int, strip(ln[1:4])); (y === nothing || y < 1000 || y > 3000) && continue
                t = tryparse(Int, strip(ln[9:14])); t === nothing && continue
                push!(rows, t)
            end
            length(rows) >= 2 ? rows[2] : -1     # rows[1]=inventory (0), rows[2]=first projection
        end

        dir = mktempdir()
        # (STAND_CN, live cyc-1 ingrowth TPA measured on FVSem_clean)
        cases = (("85271557010661", 222),   # the traced stand (DF 198 + LP 24)
                 ("70100941010661", 238),
                 ("67176330010661", 240),
                 ("2967061010690",  204))
        for (cn, live_tpa) in cases
            k = joinpath(dir, "s.key"); write(k, keytext(cn))
            out = FVSjl.run_keyfile(k; variant = FVSjl.EasternMontana())
            @test cyc1_tpa(out) == live_tpa
        end
    end
end
