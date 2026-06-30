# =============================================================================
# diameter_growth.jl (centralstates) — CS large-tree DG (cs/dgf.f) [STUB — chunk 3]
#
# The CS diameter-growth model (cs/dgf.f) is the one genuinely-new CS routine: an
# SN-family ln(DDS) regression (DBH / site / crown / BA-percentile / QMD, trees ≥5").
# It lands in chunk 3. For now this file supplies ONLY `cs_dgcons!` — the per-stand
# bark copy that the volume CFTOPK broken-top path and the shared DGDRIV calibration
# read (calib.bark_a/b). The CS bark model is the constant BKRAT (bratio.f), encoded
# as intercept 0 / slope BKRAT — identical structure to NE's ne_dgcons! bark copy.
# =============================================================================

"""
    cs_dgcons!(s)

CS per-stand DG setup — for now just the per-stand bark copy (`calib.bark_a/b` ←
intercept 0 / slope BKRAT). The CS site/slope/forest-type DGCON constant + the
serial-correlation ATTEN clock land with the cs/dgf.f port (chunk 3); until then
dg_const=0 / atten=1000 (inert at cycle-0, where no DG is projected).
"""
function cs_dgcons!(s::StandState)
    c = s.calib; sd = s.coef.species
    ba = sd[:bark_intercept]; bb = sd[:bark_slope]
    @inbounds for sp in 1:MAXSP
        c.bark_a[sp] = ba[sp]; c.bark_b[sp] = bb[sp]
        c.dg_const[sp] = 0f0; c.atten[sp] = 1000f0
    end
    return s
end
