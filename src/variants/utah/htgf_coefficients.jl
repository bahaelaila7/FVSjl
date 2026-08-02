# =============================================================================
# htgf_coefficients.jl (utah) — Schreuder-Hafley SBB height coefficients (ut/htgf.f).
# Loaded at include time from data/utah/htgf_{cof,zbias}_ut.csv (tools/utah/extract_htgf_coeffs.jl,
# verified vs source). UT_HTCOF[K,L] = COF(L,K), K=(JSPC-1)*5+KEYCR|LSIMAP (1..55), L=1..9.
# =============================================================================

let
    cofp = joinpath(UT_DATADIR, "htgf_cof_ut.csv")
    rows = [split(strip(l), ',') for l in readlines(cofp)[2:end]]
    global const UT_HTCOF = Float32[parse(Float32, rows[k][l+1]) for k in 1:55, l in 1:9]   # [K,L]
    zp = joinpath(UT_DATADIR, "htgf_zbias_ut.csv")
    zr = [split(strip(l), ',') for l in readlines(zp)[2:end]]
    global const UT_AZBIAS = Float32[parse(Float32, zr[sp][2]) for sp in 1:24]
    global const UT_BZBIAS = Float32[parse(Float32, zr[sp][3]) for sp in 1:24]
    global const UT_MPCRSI = Int[round(Int, parse(Float32, zr[sp][4])) for sp in 1:24]
end
const _UT_XI1 = 0.1f0
const _UT_XI2 = 4.5f0
