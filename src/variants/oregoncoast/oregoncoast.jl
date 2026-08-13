# =============================================================================
# oregoncoast.jl — the OC (Oregon Coast / SW Oregon) variant singleton + registration.
#
# Ported from: oc/*.f (FVS "SORNEC / Oregon-Coast"). OC is NOT a Wykoff-DDS variant like the
# rest of the western cluster — its ENTIRE per-cycle growth engine is ORGANON (the ~13k-line
# organon/ + ~1.3k-line vorganon/ subsystem), edition **SWO (Southwest Oregon)**. MEASURED from
# a live FVSoc_clean run (DEBUG DGF): the .out header reads "VERSION FS2026.1 ORGANON SWO", the
# DGF path loads the FVS tree list into ORGANON via org_intree (`I,PTNO,SPECIES,DBH1,HT1OR,CR1,
# EXPAN1,...` dump), CRATET emits "CRATET ORGANON ERROR CODE, CYCLE= 1 IERROR= 0", and diameter/
# height/crown growth come back through the ORGANON ACALIB/TMPCAL calibration arrays. There is
# effectively NO FVS-native large-tree DGF to validate independently — DGF *is* ORGANON.
#
# Distinguishing infra (MEASURED, oc/grinit.f + oc/blkdat.f + common/PRGPRM.F77):
#   • VARACD='OC', MAXSP=50, seed 55329 (shared).
#   • FINT=FINTH=FINTM=5, IFINT=IFINTH=5  ⇒ 5-YEAR cycle (ORGANON's native step; unlike the
#     western Wykoff variants' 10-yr).
#   • LZEIDE=.FALSE. (Stage SDI).  DGSD=0 (NO DG serial-correlation — ORGANON is deterministic
#     per its own draws).  LHTDRG=.FALSE. all species (height comes from ORGANON, no ht-drag).
# Species JSP (50): PC IC RC GF RF SH DF WH MH WB KP LP CP LM JP SP WP PP MP GP WJ BR GS PY OS
#   LO CY BL EO WO BO VO IO BM BU RA MA GC DG FL WN TO SY AS CW WI CN CL OH RW.  Site species = DF (7).
# The 50 FVS species collapse onto the ~18 SWO ORGANON species via oc/orgspc.f (species chunk).
#
# Oracle: /workspace/.ocwork/FVSoc_clean (relink_oc.sh, gfortran-16 + isoc23 shim). Growth runs
# and the DGF debug dump is fully available; the run SIGSEGVs later in the FVS-native volume path
# (fvsvol.f:530 via natcrs_ ← vols.f:319) — a VOLUME-stage crash that does NOT block growth
# validation (documented in docs/OC_VARIANT_PORT_AUDIT.md).
#
# PORT STATUS: chunk 0 FOUNDATION only. The ORGANON growth subsystem is UNPORTED and is the
# real cost driver — see docs/OC_ORGANON_PORT_PLAN.md. Un-ported hooks error loudly (doctrine #5).
# =============================================================================

"""
    OregonCoast <: AbstractVariant

The FVS Oregon Coast variant (VARACD "OC", MAXSP = 50) — an ORGANON-coupled variant (edition SWO)
whose entire growth engine is the ORGANON subsystem rather than the western Wykoff DDS. Pass
`OregonCoast()` as the `variant`. PORT IN PROGRESS (chunk 0 foundation; ORGANON unported).
"""
struct OregonCoast <: AbstractVariant end

variant_code(::OregonCoast) = "OC"
nspecies(::OregonCoast) = 50
htg_period(::OregonCoast) = 5f0     # /CONTRL/ YR = 5 (oc FINT=5) — ORGANON's native 5-yr step

const OC_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "oregoncoast"))

coefficients(::OregonCoast) = cached_coefficients(() -> load_species_coefficients(OC_DATADIR), "OC")

# SITSET (oc/sitset.f) — fan a per-species site index to species not assigned one by keyword, plus
# the R5/R6-adjusted SDImax (SDIDEF) defaults. In OC the site index and SDImax feed ONLY the ORGANON
# growth/calibration path (the ~9k-line unported engine, chunks C2-C6) and Stage self-thin mortality;
# they are cyc0-INERT for the C1 FVS↔ORGANON boundary marshalling (which reads only DBH/HT/CR/species/
# TPA/ISPECL). Full oc/sitset.f (SITEAR fan, R5SDI/R5ADJ/R6ADJ, SDImax) is a C2 deliverable; this
# minimal stub lets an OC stand initialize without fabricating growth-affecting site coefficients.
site_setup!(s::StandState, ::OregonCoast) = s
