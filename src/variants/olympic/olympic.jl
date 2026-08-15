# =============================================================================
# olympic.jl — the OP (Olympic / Pacific-NW ORGANON) variant singleton + registration.
#
# OP is the ORGANON follow-on to OC (Oregon Coast). Like OC its entire per-cycle growth engine is
# the shared ORGANON subsystem (organon/ + vorganon/), but OP selects a DIFFERENT ORGANON edition:
#
#   MEASURED from a live FVSop_clean run (relinked /workspace/.opwork, DEBUG DGF, stand S248112):
#     • .out header reads "VERSION FS2026.1 ORGANON NWO&SMC" (grohed.f).
#     • sitset.f defaults IMODTY=2 and sets VERSION=IMODTY; the run logs "MODEL TYPE NOT RECOGNIZED.
#       BEING SET TO 2= NORTHWEST OREGON" ⇒ OP runs the **NWO (Northwest Oregon, VERSION=2)** ORGANON
#       edition (SMC=3 is keyword-selectable; RAP=4 red-alder-plantation exists but is not the default).
#     • MAXSP = 39 (op/orgspc.f species order; PRGPRM.F77), vs OC's 50.
#
# Reuse verdict (MEASURED against organon/diagro.f): the shared ORGANON *engine structure* — the
# DIAMGRO_RUN driver, GET_BAL/SSTATS stand statistics, SUBMAX, the ln(DDS) diameter-growth form, the
# CRADJ crown adjustment, DGCALIB, PREPARE — is version-agnostic and REUSABLE. What differs for NWO
# and must be ADDED is DATA-only: the DG_NWO coefficient table (11 species groups × 11 params) + the
# NWO full-adjustment ADJ table + the OP species map (39 FVS species → ORGANON FIA code → NWO group).
# The OC Julia port hard-coded VERSION=1 (SWO) and raises on VERSION≠1, so OP carries its own NWO
# data alongside (see species.jl + organon_diamgro_nwo.jl). Growth is deterministic (DGSD=0), so the
# NWO diameter-growth is a genuine BIT-EXACT target, not an RNG straddle.
#
# Oracle: /workspace/.opwork/FVSop_clean (build_g16.sh, gfortran-16 *.o + isoc23 shim). The DGF debug
# dump is fully available; like OC the run later SIGSEGVs in the FVS-native volume path (fvsvol.f) —
# a VOLUME-stage crash that does not block growth validation.
#
# PORT STATUS: FOUNDATION + validated DG_NWO diameter-growth core. The full ORGANON marshalling /
# PREPARE setup / height / crown / mortality / volume integration for NWO is the follow-on (mirrors
# OC chunks C1-C10). Un-ported hooks error loudly (doctrine #5). See docs/OP_VARIANT_PORT_AUDIT.md.
# =============================================================================

"""
    Olympic <: AbstractVariant

The FVS Olympic variant (VARACD "OP", MAXSP = 39) — an ORGANON-coupled variant, edition **NWO**
(Northwest Oregon, VERSION=2). Pass `Olympic()` as the `variant`. PORT IN PROGRESS: foundation +
DG_NWO diameter-growth core validated bit-exact; the rest of the ORGANON engine is unported.
"""
struct Olympic <: AbstractVariant end

variant_code(::Olympic) = "OP"
nspecies(::Olympic) = 39
htg_period(::Olympic) = 5f0     # /CONTRL/ YR = 5 (op FINT=5) — ORGANON's native 5-yr step
organon_version(::Olympic) = 2  # MEASURED: sitset.f IMODTY=2 default ⇒ VERSION=2 (NWO)

const OP_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "olympic"))

# op species table (op/blkdat.f JSP/FIAJSP/PLNJSP, 39 species) + op/bratio.f bark. Chunk 1 carries
# the code columns + bark; the full merch/htdbh/crown/site coefficient columns are follow-ons (the
# FVS-native DGF/HTGF port uses its own hard-coded op/dgf.f + op/htgf.f tables, not this dict).
coefficients(::Olympic) = cached_coefficients(() -> load_species_coefficients(OP_DATADIR), "OP")

# op/sitset.f — ORGANON site-index conversion between DF (species 16) and WH (species 19), Nigh
# (1995) Forest Science 41:84-98. Analogous to OC's DF↔PP; MEASURED from op/sitset.f:144-156. Wired
# here for the foundation; the full R6ADJ site fan + ORGANON SI_1/SI_2 marshalling is a follow-on.
function site_setup!(s::StandState, ::Olympic)
    p = s.plot; si = p.sp_site_index
    if length(si) >= 19 && (si[16] > 0f0 || si[19] > 0f0)
        if si[16] <= 0f0
            si[16] = 0.480f0 + 1.110f0*si[19]      # DF from WH
        elseif si[19] <= 0f0
            si[19] = -0.432f0 + 0.899f0*si[16]     # WH from DF
        end
    end
    return s
end
