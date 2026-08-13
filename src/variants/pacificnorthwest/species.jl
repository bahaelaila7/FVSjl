# =============================================================================
# species.jl (pacificnorthwest) — PN species-coefficient table binding + block-data init. Chunk 1.
#
# PN per-species coefficient DATA (data/pacificnorthwest/species_coefficients.csv), extracted verbatim
# from pn/*.f (mirrors the WC layout): bark1/bark2/bark_imap (pn/bratio.f), dg_resid_sd (SIGMAR),
# ht1/ht2 (pn/htcalc.f Wykoff), sichg_a/b/refage/refloc (pn/sichg.f), site_redux (pn/sitset.f),
# crown_imap (pn/crown.f), is_sprouting (pn/blkdat.f ISPSPE). 39 species; slot 6 = SS/Sitka, slot 38 blank.
# =============================================================================

const PN_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "pacificnorthwest"))

coefficients(::PacificNorthwest) = cached_coefficients(() -> load_species_coefficients(PN_DATADIR), "PN")

"""
    init_blockdata!(s, ::PacificNorthwest)

PN BLOCK DATA init (pn/blkdat.f + pn/grinit.f): fan the species code arrays into `SpeciesData`, set the
default tree format, and apply the PN /CONTRL/ + SDI flags via `pn_grinit!` (Reineke SDI, DGSD 1.7).
"""
function init_blockdata!(s::StandState, v::PacificNorthwest)
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
    pn_grinit!(s)
    return s
end

load_species_coefficients!(s::StandState, v::PacificNorthwest) = init_blockdata!(s, v)

spctrn_column(::PacificNorthwest) = 4
other_species(::PacificNorthwest) = Int32(39)
