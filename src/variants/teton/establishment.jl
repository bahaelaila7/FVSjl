# =============================================================================
# establishment.jl (teton) — TT establishment constants (tt/blkdat.f ESCOMN).
#
# TT reuses the shared establish! engine (like CR/IE). These are the two per-species
# clamps applied to the established/planted-tree height (tt/estab.f:489-492):
#   HHT < XMIN(sp)   → HHT = XMIN(sp)     (establishment MIN height)
#   HHT > HHTMAX(sp) → HHT = HHTMAX(sp)   (establishment MAX height cap)
# The base height itself uses the IE-style first-cut placeholder (= XMIN) in the
# shared engine's ESSUBH branch until a faithful TT planted-tree height model is ported.
# =============================================================================

# tt/blkdat.f ESCOMN: XMIN (establishment min height, ft) per species 1..18.
const _TT_ES_XMIN = Float32[1.0, 1.0, 1.0, 0.5, 0.5, 6.0, 1.0, 0.5, 0.5, 1.0,
                            0.5, 0.5, 0.5, 6.0, 3.0, 0.5, 0.5, 3.0]

# tt/blkdat.f ESCOMN: HHTMAX (establishment max reported height, ft) per species 1..18.
const _TT_ES_HHTMAX = Float32[23.0, 27.0, 21.0, 6.0, 18.0, 20.0, 24.0, 18.0, 18.0, 17.0,
                              6.0, 6.0, 6.0, 16.0, 16.0, 6.0, 22.0, 16.0]

# tt/essubh.f: subsequent/planted-tree base height per species (fixed constants for all but PP).
# PP (sp10) uses a CI-variant regression (PN=−1.9948+1.53946·lnAGE−0.00402·BAA−0.1471+UPRE(IPREP)−0.01155·ELEV;
# HHT=exp(PN+EMSQR·DILATE·BNORM·0.49076)) — placeholder here (=1.0, PP's XMIN) pending the estab-engine params.
# HHT is then clamped to [XMIN, HHTMAX] in the shared engine (tt/estab.f:489-492).
const _TT_ESSUBH_HHT = Float32[1.0, 0.5, 2.0, 0.5, 1.5, 5.0, 3.0, 1.5, 0.75, 1.0,
                               0.5, 0.5, 1.0, 5.0, 10.0, 1.0, 1.0, 10.0]
