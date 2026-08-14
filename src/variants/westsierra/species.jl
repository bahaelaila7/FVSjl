# =============================================================================
# species.jl (westsierra) — WS species-coefficient table binding + block-data init. Chunk 1.
#
# WS per-species coefficient DATA (data/westsierra/species_coefficients.csv), extracted verbatim from
# ws/*.f: bark1/bark2/bark_imap (ws/bratio.f BARK1/BARK2/IMAP) and dg_resid_sd (ws/blkdat.f SIGMAR). 43
# species. Later chunks add their columns. `ws_bratio` = the WS bark ratio (ws/bratio.f IMAP dispatch):
#   IMAP 1 → constant BARK1; IMAP 2 → BARK1 + BARK2/max(D,1) (reciprocal, high-elev pines); IMAP 3 →
#   BARK1·D^BARK2/D (power, GS-sequoia/RW-redwood). Clamped [0.80,0.99] for IMAP 2/3 only.
# ⚠ WS's IMAP encoding differs from the shared wc_bratio's eqtype (wc eqtype 2 is b1+b2·D linear, WS imap 2
#   is b1+b2/D reciprocal) ⇒ WS needs THIS bespoke ws_bratio, NOT wc_bratio.
# =============================================================================

const WS_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "westsierra"))

coefficients(::WestSierra) = cached_coefficients(() -> load_species_coefficients(WS_DATADIR), "WS")

"""
    ws_bratio(sd, sp, d) -> bark ratio (DIB/DOB)

ws/bratio.f: IMAP dispatch, clamped [0.80, 0.99] for IMAP 2/3.
  IMAP 1 → BARK1 (constant, unclamped — values already in range).
  IMAP 2 → BARK1 + BARK2·(1/TEMD), TEMD=max(D,1)  (reciprocal).
  IMAP 3 → DIB=BARK1·D^BARK2; BRATIO=DIB/D          (power; uses D not TEMD).
"""
@inline function ws_bratio(sd::Dict{Symbol,Vector{Float32}}, sp::Integer, d::Real)
    D = Float32(d)
    imap = Int(sd[:bark_imap][sp]); b1 = sd[:bark1][sp]; b2 = sd[:bark2][sp]
    imap == 1 && return b1                                  # constant (unclamped in ws/bratio.f)
    D <= 0f0 && return 0.99f0                               # power form guards against D≤0 (inert for real trees)
    br = if imap == 2
        temd = D < 1f0 ? 1f0 : D
        b1 + b2 * (1f0 / temd)
    else                                                   # IMAP 3 (power)
        (b1 * D^b2) / D
    end
    br > 0.99f0 && (br = 0.99f0)
    br < 0.80f0 && (br = 0.80f0)
    return br
end

"""
    init_blockdata!(s, ::WestSierra)

WS BLOCK DATA init (ws/blkdat.f + ws/grinit.f): fan the species code arrays into `SpeciesData`, set the
default tree format, and apply the WS /CONTRL/ + SDI flags via `ws_grinit!` (DGSD 2.0, Zeide, seed 55329).
"""
function init_blockdata!(s::StandState, v::WestSierra)
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
    ws_grinit!(s)
    return s
end

load_species_coefficients!(s::StandState, v::WestSierra) = init_blockdata!(s, v)

spctrn_column(::WestSierra) = 4               # target_ws column in species_translation.csv
other_species(::WestSierra) = Int32(42)       # OS (other softwood); OH(43) is the hardwood catch-all
