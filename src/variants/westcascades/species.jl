# =============================================================================
# species.jl (westcascades) — WC species-coefficient table binding + block-data init. Chunk 1.
#
# Chunk 1 = the WC per-species coefficient DATA table (data/westcascades/species_coefficients.csv),
# extracted verbatim from the wc/*.f block data / coefficient sources:
#   • bark1/bark2/bark_imap  — wc/bratio.f JBARK[39]→BARKB[4,14] resolved per species
#     (bark_imap = BARKB eq-type: 1=power a·Dᵇ, 2=linear a+b·D). Consumed by wc_bratio (dgf.f).
#   • dg_resid_sd            — wc/blkdat.f SIGMAR (DG residual SD; drives the DGSD OLDRN draw).
#   • ht1/ht2               — wc/blkdat.f HT1/HT2 (Wykoff HTCALC height-dub coeffs, ≥5" DBH).
#   • sichg_a/b/refage/refloc— wc/sichg.f A/B/REFAGE/REFLOC (site-age reference conversion, chunk 2).
#   • site_redux             — wc/sitset.f misc-hardwood SI reduction factors (chunk 2; sp20 MH & sp28
#                              WO carry special formulas handled in site_index.jl).
#   • crown_imap             — wc/crown.f IMAP[39]→16 crown groups (reference for the crown chunk 5).
# The 39-species code arrays (JSP/FIAJSP/PLNJSP) mirror wc/blkdat.f; slot 6 & 38 are blank ("__").
#
# Distinguishing infra is set by wc_grinit! (westcascades.jl): Reineke SDI, DGSD 1.7, seed 55329.
# =============================================================================

const WC_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "westcascades"))

coefficients(::WestCascades) = cached_coefficients(() -> load_species_coefficients(WC_DATADIR), "WC")

"""
    init_blockdata!(s, ::WestCascades)

WC BLOCK DATA init (wc/blkdat.f + wc/grinit.f): fan the species code arrays into the stand's
`SpeciesData`, set the default tree format, and apply the WC /CONTRL/ + SDI-family flags via
`wc_grinit!` (Reineke SDI, DGSD 1.7, RNG seed 55329, 10-yr native period).
"""
function init_blockdata!(s::StandState, v::WestCascades)
    sd = s.species
    alpha = s.coef.code_alpha; fia = s.coef.code_fia; plants = s.coef.code_plants
    @inbounds for i in 1:nspecies(v)
        sd.alpha[i]  = alpha[i]
        sd.fia[i]    = fia[i]
        sd.plants[i] = plants[i]
        code = rstrip(alpha[i])
        sd.class_codes[i, 1] = code * "1"
        sd.class_codes[i, 2] = code * "2"
        sd.class_codes[i, 3] = code * "3"
        sd.code2[i] = String(rstrip(first(sd.class_codes[i, 1], 2)))
    end
    hab = s.coef.valid_habitat
    copyto!(s.plot.valid_habitat, 1, hab, 1, min(length(hab), length(s.plot.valid_habitat)))

    s.control.tree_format = DEFAULT_TREE_FORMAT
    wc_grinit!(s)                                  # Reineke SDI, DGSD 1.7, seed 55329, 10-yr YR (westcascades.jl)
    return s
end

load_species_coefficients!(s::StandState, v::WestCascades) = init_blockdata!(s, v)

# SPCTRN crosswalk target column (WC CSV mirrors the 7-col layout, target_wc ×4) ⇒ col 4; "other" = OT (39).
spctrn_column(::WestCascades) = 4
other_species(::WestCascades) = Int32(39)
