# =============================================================================
# dg_coefficients.jl (utah) — large-tree DG coefficient arrays (ut/dgf.f DATA).
# Loaded at include time from data/utah/dg_coeffs_ut.csv (tools/utah/extract_dg_coeffs.jl,
# verified bit-exact vs the Fortran source). 1D → Float32[24]; 2D → [L, 24] matching Fortran
# ARR(L,K) (the CSV stores species-major, transposed back here).
# =============================================================================

let
    path = joinpath(UT_DATADIR, "dg_coeffs_ut.csv")
    hdr = split(strip(readlines(path)[1]), ',')
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    ci(name) = findfirst(==(name), hdr)
    col1d(name) = Float32[parse(Float32, rows[sp][ci(name)]) for sp in 1:24]
    # 2D: [L, 24], column L of species sp from "name_L"
    col2d(name, L) = Float32[parse(Float32, rows[sp][ci("$(name)_$k")]) for k in 1:L, sp in 1:24]

    global const UT_DGLD   = col1d("DGLD")
    global const UT_DGCR   = col1d("DGCR")
    global const UT_DGCRSQ = col1d("DGCRSQ")
    global const UT_DGPCCF = col1d("DGPCCF")
    global const UT_DGBA   = col1d("DGBA")
    global const UT_DGBAL  = col1d("DGBAL")
    global const UT_DGDBAL = col1d("DGDBAL")
    global const UT_DGCCFA = col1d("DGCCFA")
    global const UT_OBSERV = col1d("OBSERV")
    global const UT_DGCASP = col1d("DGCASP")
    global const UT_DGSASP = col1d("DGSASP")
    global const UT_DGSLOP = col1d("DGSLOP")
    global const UT_DGSLSQ = col1d("DGSLSQ")
    global const UT_DGEL   = col1d("DGEL")
    global const UT_DGEL2  = col1d("DGEL2")
    global const UT_ISMAP  = Int[round(Int, x) for x in col1d("ISMAP")]

    global const UT_IBSERV = [round(Int, x) for x in col2d("IBSERV", 5)]   # [5, 24] (site class × sp)
    global const UT_MAPLOC = [round(Int, x) for x in col2d("MAPLOC", 6)]   # [6, 24] (forest × sp)
    global const UT_MAPDSQ = [round(Int, x) for x in col2d("MAPDSQ", 6)]   # [6, 24]
    global const UT_IDGSIM = [round(Int, x) for x in col2d("IDGSIM", 7)]   # [7, 24] (site-base-sp × sp)
    global const UT_DGFOR  = col2d("DGFOR", 5)                             # [5, 24] (location class × sp)
    global const UT_DGDS   = col2d("DGDS", 4)                              # [4, 24] (dbh² class × sp)
    global const UT_DGSIC  = col2d("DGSIC", 5)                             # [5, 24] (site-index coef × sp)
end
