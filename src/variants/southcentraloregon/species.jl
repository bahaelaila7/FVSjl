# =============================================================================
# species.jl (southcentraloregon) — SO species-coefficient table binding + block-data init. Chunk 1.
#
# SO per-species coefficient DATA (data/southcentraloregon/species_coefficients.csv), extracted verbatim
# from so/*.f: so_brdat/so_bark1/so_bark2/so_bark_eqtype (so/bratio.f) and dg_resid_sd (so/blkdat.f SIGMAR).
# 33 species. Later chunks (site/dgf/crown/height/mort/vol) add their columns. `so_bratio` = the SO bark
# ratio (so/bratio.f SELECT CASE): constant BRDAT for most species, BARKB eqtype 1/2 for the WC/CA-borrowed
# species, and a special juniper (WJ) form — clamped [0.80,0.99].
# =============================================================================

const SO_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "southcentraloregon"))

coefficients(::SouthCentralOregon) = cached_coefficients(() -> load_species_coefficients(SO_DATADIR), "SO")

# so/bratio.f SELECT CASE — the CASE-1 "from WC/CA variant" species use the BARKB formula (eqtype 1/2);
# WJ(11) uses the UT juniper form; every other species is the constant BRDAT.
const SO_BARK_CASE1 = (9, 15, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 33)

"""
    so_bratio(sd, sp, d) -> bark ratio (DIB/DOB)

so/bratio.f: three-path SELECT CASE, clamped to [0.80, 0.99].
  CASE-1 spp (WC/CA-borrowed): eqtype 1 → DIB=b1·D^b2; eqtype 2 → DIB=b1+b2·D; BRATIO=DIB/D (D>0 else 0.99).
  CASE(11) WJ juniper (from UT): 0.9002 − 0.3089/clamp(D,1,19).
  DEFAULT: the per-species constant BRDAT.
"""
@inline function so_bratio(sd::Dict{Symbol,Vector{Float32}}, sp::Integer, d::Real)
    D = Float32(d)
    local br::Float32
    if sp == 11                                             # WJ juniper (UT form)
        temd = D < 1f0 ? 1f0 : (D > 19f0 ? 19f0 : D)
        br = 0.9002f0 - 0.3089f0 / temd
    elseif sp in SO_BARK_CASE1                              # WC/CA-borrowed: BARKB formula
        if D <= 0f0
            br = 0.99f0
        else
            et = Int(sd[:so_bark_eqtype][sp])
            b1 = sd[:so_bark1][sp]; b2 = sd[:so_bark2][sp]
            br = et == 1 ? (b1 * D^b2) / D :
                 et == 2 ? (b1 + b2 * D) / D : sd[:so_brdat][sp]
        end
    else                                                   # constant BRDAT
        br = sd[:so_brdat][sp]
    end
    br > 0.99f0 && (br = 0.99f0)
    br < 0.80f0 && (br = 0.80f0)
    return br
end

"""
    init_blockdata!(s, ::SouthCentralOregon)

SO BLOCK DATA init (so/blkdat.f + so/grinit.f): fan the species code arrays into `SpeciesData`, set the
default tree format, and apply the SO /CONTRL/ + SDI flags via `so_grinit!` (DGSD 2.0, seed 55329).
"""
function init_blockdata!(s::StandState, v::SouthCentralOregon)
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
    so_grinit!(s)
    return s
end

load_species_coefficients!(s::StandState, v::SouthCentralOregon) = init_blockdata!(s, v)

spctrn_column(::SouthCentralOregon) = 4               # target_so column in species_translation.csv
other_species(::SouthCentralOregon) = Int32(32)       # OS (other softwood); OH(33) is the hardwood catch-all
