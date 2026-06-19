# fvsjl_cycle0.jl — initialize FVSjl on a .key and print the cycle-0 stand summary
# (the quantities that map to the first data row of the Fortran/FVSjulia .sum):
#   TPA  BA  SDI  CCF  TopHt  QMD
#
# Until FVSjl gains the volume + .sum writer (C5), this is the slice of output we
# can three-way against the oracle — it validates init / NOTRE / stand statistics
# across the expanded scenario matrix. Usage: julia fvsjl_cycle0.jl <key>
using FVSjl

key = ARGS[1]
s, _ = initialize(key)
FVSjl.notre!(s)
FVSjl.compute_density!(s)

# per-acre quantities divide by the reciprocal stockable multiplier (gross_space);
# top height and QMD are per-tree and are not scaled.
g    = s.plot.gross_space
tpa  = FVSjl.stand_tpa(s) / g
ba   = FVSjl.stand_ba(s)  / g
sdi  = FVSjl.stand_sdi(s) / g
ccf  = FVSjl.stand_ccf(s) / g
toph = FVSjl.stand_top_height(s)
qmd  = FVSjl.stand_qmd(s)
println(join(round.(Int, (tpa, ba, sdi, ccf, toph)), " "), " ", round(qmd; digits=1))
