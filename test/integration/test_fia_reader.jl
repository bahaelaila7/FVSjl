# test_fia_reader.jl — native FIA "FVS-ready" database INPUT path (src/io/fia_database.jl).
#
# Self-contained: a tiny extracted SQLite fixture (test/fixtures/fia/ls_sample.db, two real
# LS stands) + the golden live-FVS .sum captured from that same fixture. Proves FVSjl's
# DATABASE/DSNIN reader ingests an FIA stand and reproduces live's cycle-0 inventory
# BIT-EXACT on the core columns (TPA/BA/SDI/CCF/TopHt/QMD/CuFt) AND BdFt across all four
# variants — including the R9 per-national-forest board-type gate (Scribner vs International ¼").

using Test
using FVSjl

const FIA_DIR = joinpath(@__DIR__, "..", "fixtures", "fia")
const FIA_DB  = abspath(joinpath(FIA_DIR, "ls_sample.db"))

# Parse the first (cycle-0) data row of a .sum into its whitespace fields.
_sum_cycle0(text) = for ln in split(text, '\n')
    f = split(strip(ln)); length(f) >= 8 || continue
    y = tryparse(Int, f[1]); (y === nothing || y < 1000) && continue
    return f
end

_fia_keyfile(cn) = """
STDIDENT
$cn
DATABASE
DSNin
$FIA_DB
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         1.0
ECHOSUM
PROCESS
STOP
"""

@testset "FIA-ready database reader (DATABASE/DSNIN)" begin
    # Cross-variant: one+ real FIA stand per variant, cycle-0 vs the golden live .sum.
    # Density/size (TPA/BA/SDI/TopHt/QMD) + cubic volume + CCF are BIT-EXACT on every
    # stand/variant — CCF once species resolve (the _fia_spcode 3-digit-FIA fix) AND the
    # stand's lat/long feed the Hopkins index (HI-dependent hardwood crowns). BdFt is
    # bit-exact except the LS hardwood stand, which keeps a small hardwood board-vol residual.
    density = [(3, "TPA"), (4, "BA"), (5, "SDI"), (6, "CCF"), (7, "TopHt"), (8, "QMD"),
               (9, "TotCuFt"), (10, "MerchCuFt")]
    # (stand_cn, variant) — BdFt is BIT-EXACT on every stand once the R9 per-national-forest
    # board-type gate (volinit.f:434-451; _R9_INTL_BDFT_FORESTS) is honored: LS conifer
    # IFORST=10 → Scribner, LS hardwood IFORST=24 → International, NE/CS/SN → their native board.
    cases = [("100180735010661", FVSjl.LakeStates()),    # LS conifer (IFORST=10 → Scribner)
             ("255262523010854", FVSjl.Southern()),      # SN
             ("14173137020004",  FVSjl.CentralStates()), # CS
             ("657546100126144", FVSjl.Northeast()),     # NE (CCF fixed by lat/long)
             ("55482390010661",  FVSjl.LakeStates())]    # LS hardwood (IFORST=24 → R9 International board)
    for (cn, var) in cases
        golden = read(joinpath(FIA_DIR, "ls_$(cn).live.sum"), String)
        lv = _sum_cycle0(golden)
        key = joinpath(mktempdir(), "fia_$(cn).key")
        write(key, _fia_keyfile(cn))
        out = FVSjl.run_keyfile(key; variant = var)
        jl = _sum_cycle0(out)
        @test jl !== nothing && lv !== nothing
        @test jl[1] == lv[1]                       # inventory year
        for (i, name) in density
            @test jl[i] == lv[i]                   # BIT-EXACT rendered inventory stat (incl CCF)
        end
        @test jl[12] == lv[12]                     # BdFt BIT-EXACT (per-forest board type)
    end
end
