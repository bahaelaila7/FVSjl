# =============================================================================
# species.jl (centralcalifornia) — CA species-coefficient table binding + block-data init. Chunk 1.
#
# CA per-species coefficient DATA (data/centralcalifornia/species_coefficients.csv), extracted verbatim
# from ca/*.f: bark1/bark2/bark_imap (ca/bratio.f JBARK→BARKB, 3-method a·DOB^b / a+b·DOB / a·DOB), and
# dg_resid_sd (ca/blkdat.f SIGMAR). 50 species. Later chunks (site/crown/height/mort/vol) add their columns.
# =============================================================================

const CA_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "centralcalifornia"))

coefficients(::CentralCalifornia) = cached_coefficients(() -> load_species_coefficients(CA_DATADIR), "CA")

"""
    init_blockdata!(s, ::CentralCalifornia)

CA BLOCK DATA init (ca/blkdat.f + ca/grinit.f): fan the species code arrays into `SpeciesData`, set the
default tree format, and apply the CA /CONTRL/ + SDI flags via `ca_grinit!` (Zeide SDI, DGSD 1.7).
"""
function init_blockdata!(s::StandState, v::CentralCalifornia)
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
    ca_grinit!(s)
    return s
end

load_species_coefficients!(s::StandState, v::CentralCalifornia) = init_blockdata!(s, v)

spctrn_column(::CentralCalifornia) = 4
other_species(::CentralCalifornia) = Int32(25)   # OS (other softwood); OH(49) is the hardwood catch-all
