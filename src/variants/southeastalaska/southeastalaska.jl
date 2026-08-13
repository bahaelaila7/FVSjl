# =============================================================================
# southeastalaska.jl — the AK (Southeast Alaska) variant singleton + registration
#
# Ported from: ak/*.f (the FVS "ALASKA" variant, USFS Region-10). AK is a STANDARD
# WYKOFF DDS diameter-growth variant (ak/dgf.f reads the classic single-equation
# ln(DDS) form, NO GENGYM) — like UT/IE/EM/KT/TT it is a specialization of the shared
# Wykoff engine, differing in (a) the 23-species list + coefficients, (b) a UNIQUE
# PERMAFROST diameter-growth modifier (ak/dgf.f PFCON..PFSASP; multiplicative PFMOD on
# the annual base DG for permafrost-affected species 4:7,13,16:23), (c) AK's own SEAMRT
# mortality-distribution routine (percentile + shade-tolerance geometric progression),
# and (d) Region-10 volume (sitset VOLEQDEF VAR='AK' IREGN=10 → NVEL + ak logs/cubrds).
#
# AK grinit (ak/grinit.f): VARACD='AK', MAXSP=23, LZEIDE=.TRUE. (Zeide SDI — like CR/UT),
# DGSD=2.0 (DG serial-corr + ZZRAN), FINT=10 (10-yr cycle), FINTH=5, FINTM=5,
# LHTDRG(I)=.FALSE. (default — UNLIKE UT's .TRUE.), RNG seed 55329 (shared blkdat).
#
# Species JSP (23): SF AF YC TA WS LS BE SS LP RC WH MH OS AD RA PB AB BA AS CW WI SU OH
# FIA: 011 019 042 071 094 094 095 098 108 242 263 264 299 350 351 375 376 741 746 747 920 928 998
#   (LS "Lutz's/hybrid spruce" has no native FIA code; sitset.f:321 maps JSP='LS' → IFIASP=94.)
# SPCTRN: AK maps to shared ASPT column 4 (ak spctrn.f `SPCOUT=ASPT(I,4)`).
#
# Validation oracle: live Fortran relinked from bin/FVSak_buildDir/*.o
# (/workspace/.akwork/FVSak_clean via relink_ak.sh, gfortran-16 + isoc23 shim). The
# large-tree DGF ln(DDS) equation + coefficients are VERIFIED BIT-EXACT vs the live
# oracle's `DEBUG DGF` per-tree dump on the akt01 reference stand (54 tree-records,
# 0 mismatches; see docs/AK_VARIANT_PORT_AUDIT.md).
#
# PORT IN PROGRESS (chunk 0 beachhead): singleton + species + grinit constants + the
# validated DGF equation/coefficients. Full engine integration (site_setup!/crown/
# density point-Zeide/mortality/volume) + the PERMAFROST + PRD point-Zeide path are
# later runs. Until a hook is implemented, dispatching it on SoutheastAlaska() errors.
# =============================================================================

"""
    SoutheastAlaska <: AbstractVariant

The FVS Southeast Alaska variant (VARACD "AK", MAXSP = 23, Region-10) — a western
Wykoff-DDS variant (Zeide SDI) with a unique permafrost diameter-growth modifier and
Region-10 volume. Pass `SoutheastAlaska()` as the `variant` to `run_keyfile`/`each_stand`.
PORT IN PROGRESS (chunk 0 beachhead — DGF validated bit-exact).
"""
struct SoutheastAlaska <: AbstractVariant end

variant_code(::SoutheastAlaska) = "AK"
nspecies(::SoutheastAlaska) = 23
htg_period(::SoutheastAlaska) = 10f0     # /CONTRL/ YR = 10 (ak FINT=10), like CR/UT/KT/EM

const AK_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "southeastalaska"))

coefficients(::SoutheastAlaska) = cached_coefficients(() -> load_species_coefficients(AK_DATADIR), "AK")
