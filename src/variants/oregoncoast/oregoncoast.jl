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
variant_maxtre(::OregonCoast) = 2000   # oc PRGPRM.F77 PARAMETER (MAXTRE=2000)
htg_period(::OregonCoast) = 5f0     # /CONTRL/ YR = 5 (oc FINT=5) — ORGANON's native 5-yr step

const OC_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "oregoncoast"))

coefficients(::OregonCoast) = cached_coefficients("OC") do
    c = load_species_coefficients(OC_DATADIR)
    # FFE crown-biomass columns (the ORGANON growth CSV lacks them): fmvinit.f V2T (raw lb/cuft; the
    # shared crown_biomass applies /2000) + the BLM merch DBH gate. See oregoncoast/ffe_coefficients.jl.
    haskey(c.species, :v2t)          || (c.species[:v2t]          = copy(OC_FFE_V2T))
    haskey(c.species, :dbh_min)      || (c.species[:dbh_min]      = fill(OC_FFE_DBHMIN, 50))
    haskey(c.species, :is_sprouting) || (c.species[:is_sprouting] = copy(OC_FFE_SPROUT))  # ESTUMP sprout filter (blkdat.f ISPSPE)
    c
end

# SITSET (oc/sitset.f) — fan a per-species site index to species not assigned one by keyword, plus
# the R5/R6-adjusted SDImax (SDIDEF) defaults. In OC the site index and SDImax feed ONLY the ORGANON
# growth/calibration path (the ~9k-line engine, C2 calibration + C3-C6 growth) and Stage self-thin
# mortality; they are cyc0-INERT for the C1 FVS↔ORGANON boundary marshalling (which reads only DBH/
# HT/CR/species/TPA/ISPECL).
#
# The one ORGANON-load-bearing SITSET piece needed for the C2 calibration is the DF↔PP site-index
# conversion (oc/sitset.f:181-189): ORGANON's SITE_1=RVARS(1)=SITEAR(7) (DF) and SITE_2=RVARS(2)=
# SITEAR(18) (PP) both feed HDCALIB/CRCALIB, and when only one is set the other is derived here
# (MEASURED: on ocmin SITEAR(7)=92 from the ecoclass, SITEAR(18) unset ⇒ 0.940792·92=86.5528641,
# exactly the RVARS(2) the oracle passes to PREPARE). The R6ADJ site fan and R5SDI/SDImax defaults,
# and the RVARS(3-5)=SDIDEF MSDI marshalling, are wired at the growth entry (C3); the ecoclass-
# derived site species value (SITEAR(7)) itself is supplied by the stand's site-index loader.
# oc/sitset.f:68-73 R6ADJ — per-species site-index adjustment relative to the Hann-Scrivani DF
# site index (SITEAR(I)=HGUESS·R6ADJ(I) for a species with no keyword site index; oc/sitset.f:197).
const OC_R6ADJ = Float32[
    0.90,0.70,0.80,1.00,1.00,1.00,1.00,0.95,0.90,0.90,
    0.90,0.90,0.90,0.90,0.94,1.00,0.94,0.94,0.90,0.90,
    0.76,0.76,1.00,0.40,0.76,0.28,0.42,0.34,0.28,0.40,
    0.56,0.76,0.28,0.76,0.56,0.76,0.76,0.76,0.40,0.70,
    0.40,0.76,0.76,0.40,0.76,0.25,0.25,0.25,0.56,1.00]

# oc/sitset.f = ca/sitset.f + its own R5SDI + the ORGANON DF↔PP site conversion (and no CALCSDI→LZEIDE reset): the R6
# forests (IFOR≥6) ECOCLS the stand's plant association (ca/habtyp.f KODTYP → PCOM; FIA PV_CODE decoded in
# fia_database.jl like CA) for the site species / SI / SDImax, then the Hann-Scrivani R6ADJ fan and the SDIDEF fan
# (BAMAX / the site species' RSDI / R5SDI). The previous OC-only port handled just the no-site default ecoclass
# (CWC221, NSISET==0), so every FIA stand (which carries a SITE_INDEX) kept SDIDEF = 0: FVS_InvReference SDIMax 0
# (live 720/815/923/955 by stand) and ORGANON's MSDI_1/2/3 = SDIDEF(7/18/4) fell back to SUBMAX's built-in A1.
const OC_R5SDI = Float32[
    570, 570, 570, 760, 800, 800, 600, 580, 580, 460,
    430, 580, 430, 460, 430, 430, 460, 430, 430, 430,
    330, 580,1052, 570, 430, 550, 550, 550, 550, 550,
    550, 550, 550, 550, 550, 550, 550, 550, 550, 550,
    550, 550, 550, 550, 550, 550, 550, 550, 550,1052]

function site_setup!(s::StandState, ::OregonCoast)
    p = s.plot
    p.forest_idx = Int32(oc_forkod(Int(p.user_forest_code)))   # oc/forkod.f (= ca/forkod.f) — KODFOR → IFOR (711→9)
    r6_ca_family_sitset!(s; r5sdi = OC_R5SDI, calcsdi_reset = false, organon_si = true)
    return s
end
