# =============================================================================
# klamath.jl — the NC (Klamath Mountains) variant singleton + registration.
#
# Ported from: nc/*.f (FVS "KLAMATH MOUNTAINS"). NC is a western Wykoff-DDS variant — closest
# template is CI (also western-Wykoff + Zeide SDI). nc/dgf.f is the standard western Wykoff DDS
# reading DGLD/DGCR/DGCRSQ/DGDBAL/DGDS/DGFOR/CONSPP (same form as IE/CI), with variant-specific
# DGFASP (aspect DG) + SDICAL/BADIST density handling.
# Distinguishing infra (MEASURED, nc/blkdat.f + nc/grinit.f):
#   • VARACD='NC', MAXSP=12, IFINT=10 (10-yr native period; nct01 runs 5-yr via keyword).
#   • LZEIDE=.TRUE. ⇒ ZEIDE SDI (like CI/UT/TT — NOT Stage). ⇒ reuse the Zeide self-thin path.
#   • DGSD=2.0 ; imperial variant ("FVS VARIANT" banner — standard 9I6 .sum, NOT metric).
# Species JSP (12): OS SP DF WF MA IC BO TO RF PP OH RW
# FIAJSP:          299 117 202 015 361 081 818 631 020 122 998 211
# PLNJSP:          2TN PILA PSME ABCO ARME CADE27 QUKE LIDE3 ABMA PIPO 2TB SESE3
# Oracle: /workspace/.ncwork/FVSnc_clean (straight relink from bin/FVSnc_buildDir/*.o + isoc23 shim, NO stubs).
# Canonical test stand: nct01 (tests/FVSnc/nct01.key + nct01.tre + nct01.sum.save; forest 371).
#
# PORT IN PROGRESS — chunk 0 scaffold. Growth/volume hooks have NO NC method yet, so dispatching
# on Klamath() errors loudly (doctrine #5) until each chunk lands + is validated vs FVSnc_clean.
# Chunk plan (mirrors CI): 1 species block-data · 2 site/SDImax (Zeide) · 3 DG DDS (western Wykoff
# + DGFASP) · 4 height · 5 CCF+crown · 6 regent · 7 mortality (Zeide self-thin) · 8 volume.
# =============================================================================

"""
    Klamath <: AbstractVariant

The FVS Klamath Mountains variant (VARACD "NC", MAXSP = 12) — a western Wykoff-DDS variant with
Zeide SDI (CI template). PORT IN PROGRESS (chunk 0 scaffold).
"""
struct Klamath <: AbstractVariant end

variant_code(::Klamath) = "NC"
nspecies(::Klamath) = 12
htg_period(::Klamath) = 5f0    # NC htgf POTHTG is a 5-yr site-curve rise (HTCALC SITAGE+5) ⇒ YR=5 (scale=fint/5)

const NC_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "klamath"))

coefficients(::Klamath) = cached_coefficients(() -> load_species_coefficients(NC_DATADIR), "NC")
