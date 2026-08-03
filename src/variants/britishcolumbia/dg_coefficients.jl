# BC large-tree DDS coefficient arrays (bc/dgf.f), loaded from data/britishcolumbia/dg_*.csv.
# BC DG = IE western Wykoff DDS (coefficients IDENTICAL to IE for the shared conifers). NI form:
#   DDS = CONSPP + DGLD·lnD + DGBAL·BAL + CR·(DGCR+CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)
#   CONSPP = DGCON + COR + 0.01·DGCCFA·RELDEN;  BAL=(1−PCT/100)·(BA/100).
let
    d = BC_DATADIR
    rows = readlines(joinpath(d, "dg_coeffs_1d.csv"))
    hdr = split(strip(rows[1]), ',')
    ci(n) = findfirst(==(n), hdr)
    col(n) = Float32[parse(Float32, strip(split(strip(rows[1+sp]), ',')[ci(n)])) for sp in 1:15]
    for n in ("DGLD","DGCR","DGCRSQ","DGBAL","DGDBAL","DGSASP","DGCASP","DGSLOP","DGSLSQ","DGEL","DGEL2","OBSERV")
        @eval global const $(Symbol("BC_" * n)) = $(col(n))
    end
    readcols(f) = (ls = readlines(joinpath(d, f));
        hcat([Float32[parse(Float32, strip(x)) for x in split(strip(l), ',')] for l in ls]...))
    global const BC_DGFOR  = readcols("dg_dgfor.csv")         # [6, 15]
    global const BC_DGDS   = readcols("dg_dgds.csv")          # [4, 15]
    global const BC_DGCCFA = readcols("dg_dgccfa.csv")        # [5, 15]
    global const BC_DGHAB  = readcols("dg_dghab.csv")         # [6, 15]
    global const BC_MAPHAB = round.(Int, readcols("dg_maphab.csv"))   # [30, 15]
    global const BC_MAPLOC = round.(Int, readcols("dg_maploc.csv"))   # [11, 15]
    global const BC_MAPDSQ = round.(Int, readcols("dg_mapdsq.csv"))   # [11, 15]
    global const BC_MAPCCF = round.(Int, readcols("dg_mapccf.csv"))   # [11, 15]
end
