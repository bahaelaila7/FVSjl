# =============================================================================
# test_ie_addtrees.jl — ADDTREES establishment keyword (ESTAB opt 28, estb/esaddt.f)
#
# ADDTREES is the external-regeneration-model BRIDGE: the keyword schedules an
# activity 432 (esin.f opt 28) that, when it fires at the ESNUTR establishment seam,
# writes a `.es1` stand summary, runs SYSTEM(CMDLN) (the external regen model), then
# reads back the model's `.es2` activity block (OPRDAT: `IACTK IDT NPRMS PRMS…`) and
# OPADD-schedules those activities — 430/431 route into the already-validated PLANT/
# NATURAL regen path. So the ported unit is the keyword + the activity-injector, and
# the added trees come from the scheduled PLANT activity firing.
#
# Validated with a USER-approved STAGED READ: a trivial CMDLN (`true`) plus a
# pre-staged deterministic `.es2` (one 430/PLANT of 400 TPA sp3), so the live FVSie
# oracle AND FVSjl read the IDENTICAL `.es2` and schedule the identical activity.
#
# LIVE-ORACLE A/B (FVSie_g16, /workspace/.iework; recorded here, reproduced hermetically
# below without the oracle):
#   • ADDTREES(es2:430) is BYTE-IDENTICAL to the direct PLANT keyword — on BOTH the
#     oracle AND FVSjl (the bridge injects exactly what a native PLANT card would).
#   • ADDTREES with NO `.es2` is BYTE-IDENTICAL to the same keyfile with no ADDTREES
#     card — on BOTH sides (a null external model contributes nothing).
#   • Re-measured 2026-09-19 (FVSie_g16, --keywordfile + ECHOSUM): with NOAUTOES at the TOP
#     level (initre.f opt 72 → ESNOAU; it is NOT an ESTAB-packet keyword — inside the packet
#     the IE oracle stops with FVS01 INVALID KEYWORD) the oracle gives bridge == direct PLANT
#     (364 TPA @2002 both). jl also gives bridge == direct, but 911: jl runs the AUTOES tally
#     under NOAUTOES because ESNOAU's STOADJ=0.0 (esin.f:788) is not persisted and the
#     estab.f STOADJ<0.0001 no-stocking branch is not ported — a separate IE establishment
#     gap, NOT a bridge defect. (With AUTOES ON the oracle's bridge ≠ direct — 555 vs 911 —
#     and jl matches both @2002.)
#   • The bridge must run at the TOP of ESNUTR (esnutr.f:59 CALL ESADDT(1)), before the
#     OPFIND 430/431 NPNATS count (:348) and the ESTAB tally (:401): otherwise the AUTOES tally
#     does not see the bridged PLANT (rule-6 catch-all + NBEST plant pool) and diverges.
#
# The two @tests below reproduce the two bit-exact bridge equalities purely in FVSjl
# (no oracle at test time): the bridge == the direct keyword, and the null bridge is
# inert. That is the port's correctness contract.
# =============================================================================

using Test
using FVSjl

# FVS keyword card: name left-justified in cols 1-8, fields right-justified in the
# 10-col fields starting at col 11 (so field n occupies cols 11+10(n-1) .. 20+10(n-1)).
_fld(x) = lpad(x, 10)
_card(name, fields...) = rpad(name, 10) * join(_fld.(fields))

const _IE_STDINFO = _card("STDINFO", "101.0", "520.0", "0.0", "315.0", "30.0", "45.0")

# Common ESTAB packet header for a bare (NOTREES) IE stand, inv year 1992, 4 cycles.
function _base_key(io, extra_estab::Vector{String})
    println(io, "STDIDENT")
    println(io, "IEADDT")
    println(io, _card("DESIGN", "", "", "", "11.0", "1.0"))
    println(io, "NOTREES")
    println(io, _IE_STDINFO)
    println(io, _card("INVYEAR", "1992"))
    println(io, _card("NUMCYCLE", "4"))
    println(io, "NOAUTOES")                  # base keyword (initre.f opt 72) — must precede the ESTAB packet
    println(io, _card("ESTAB", "1992"))
    println(io, "NOINGROW")
    for l in extra_estab; println(io, l); end
    println(io, "END")
    println(io, "PROCESS")
    println(io, "STOP")
end

write_key(path, extra) = open(io -> _base_key(io, extra), path, "w")

# Extract the numeric summary rows (YEAR TREES BA SDI CCF TopHt QMD cuft bdft) from a
# `.sum`-format string for a robust row-by-row equality that ignores the timestamp header.
function _rows(txt)
    out = String[]
    for l in split(txt, "\n")
        occursin(r"^(19|20)\d\d\s", l) && push!(out, rstrip(l))
    end
    return out
end

@testset "IE ADDTREES external-regen bridge (esaddt.f, activity 432)" begin
    # ------------------------------------------------------------------ bridge == PLANT
    mktempdir() do dir
        addt = joinpath(dir, "ie_addtrees.key")
        plant = joinpath(dir, "ie_plant.key")
        # ADDTREES: date 1992, IYR1 offset 0, IMET=1 (run external program); the NEXT card is CMDLN.
        write_key(addt, [_card("ADDTREES", "1992", "0", "1"), "true"])
        # The equivalent native PLANT: date 1992, species 3, 400 TPA, 100% survival, age 2, ht 0.5.
        write_key(plant, [_card("PLANT", "1992", "3", "400", "100.0", "2.0", "0.5")])

        # Pre-stage the `.es2` the bridge will read: <stem>_<NPLT>_<KDT>_<VARACD>.es2, KDT = the
        # end-year of the cycle the 432 fires in (1992 → cycle 1 → IY(2)-1 = 2001). One 430/PLANT
        # activity matching the direct PLANT card above.
        es2 = joinpath(dir, "ie_addtrees_IEADDT_2001_IE.es2")
        write(es2, "1\nIEADDT\n430 1992 6 3 400 100 2 0.5 0\nEnd\n")

        addt_txt  = FVSjl.run_keyfile(addt;  variant = InlandEmpire(), period = 10, output = :sum)
        plant_txt = FVSjl.run_keyfile(plant; variant = InlandEmpire(), period = 10, output = :sum)

        ra, rp = _rows(addt_txt), _rows(plant_txt)
        @test !isempty(ra)
        # The bridge schedules exactly the PLANT the external model returned ⇒ bit-exact vs the
        # native PLANT keyword (mirrors the oracle's ADDTREES==PLANT identity).
        @test ra == rp
        # Sanity: the planted regen actually shows up (non-empty stand by cycle 1).
        @test any(r -> occursin(r"^2002\s+10\s+[1-9]", r), ra)
        # The bridge honored IKEEP=1 (first .es2 line) — the file is kept, not deleted.
        @test isfile(es2)
    end

    # ------------------------------------------------------------------ null bridge is inert
    mktempdir() do dir
        addt   = joinpath(dir, "ie_addtrees.key")
        noaddt = joinpath(dir, "ie_noaddt.key")
        write_key(addt, [_card("ADDTREES", "1992", "0", "1"), "true"])
        write_key(noaddt, String[])   # identical packet, no ADDTREES card
        # No `.es2` staged ⇒ the external model returns nothing ⇒ the 432 schedules no activity.
        addt_txt   = FVSjl.run_keyfile(addt;   variant = InlandEmpire(), period = 10, output = :sum)
        noaddt_txt = FVSjl.run_keyfile(noaddt; variant = InlandEmpire(), period = 10, output = :sum)
        # Byte-identical to the no-ADDTREES packet (mirrors the oracle's null-bridge identity).
        @test _rows(addt_txt) == _rows(noaddt_txt)
    end
end
