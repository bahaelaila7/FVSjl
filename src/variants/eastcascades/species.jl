# =============================================================================
# species.jl (eastcascades) — EC species-coefficient table binding + block-data init. Chunk 1.
#
# EC per-species coefficient DATA (data/eastcascades/species_coefficients.csv), extracted verbatim from
# ec/*.f: bark1/bark2/bark_imap (ec/bratio.f), dg_resid_sd (SIGMAR), ht1/ht2 (ec/htcalc.f Wykoff),
# sichg_a/b/refage/refloc (ec/sichg.f), site_redux (ec/sitset.f), crown_imap (ec/crown.f),
# is_sprouting (ec/blkdat.f ISPSPE = {20..30}). 32 species.
# =============================================================================

const EC_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "eastcascades"))

coefficients(::EastCascades) = cached_coefficients(() -> load_species_coefficients(EC_DATADIR), "EC")

"""
    init_blockdata!(s, ::EastCascades)

EC BLOCK DATA init (ec/blkdat.f + ec/grinit.f): fan the species code arrays into `SpeciesData`, set the
default tree format, and apply the EC /CONTRL/ + SDI flags via `ec_grinit!` (Reineke SDI, DGSD 1.7).
"""
function init_blockdata!(s::StandState, v::EastCascades)
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
    ec_grinit!(s)
    return s
end

load_species_coefficients!(s::StandState, v::EastCascades) = init_blockdata!(s, v)

spctrn_column(::EastCascades) = 4
other_species(::EastCascades) = Int32(32)
