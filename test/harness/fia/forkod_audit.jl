import SQLite, DBInterface
using FVSjl
# distinct SN LOCATION codes across the WHOLE FVS-ready population
db = SQLite.DB("/workspace/SQLite_FIADB_ENTIRE.db")
locs = Dict{Int,Int}()
for r in DBInterface.execute(db, "SELECT LOCATION, COUNT(*) c FROM FVS_STANDINIT_COND WHERE VARIANT='SN' AND LOCATION IS NOT NULL GROUP BY LOCATION")
  locs[Int(r.LOCATION)] = Int(r.c)
end
println("distinct SN LOCATION codes: ", length(locs), "  (total stands ", sum(values(locs)), ")")
# apply the forkod remap + the VOLEQDEF iregn decode; flag codes that DON'T resolve to region 8 (=> blank vol_eq => 0 vol)
mutable struct P; user_forest_code::Int32; forest_idx::Int32; end
bad = Tuple{Int,Int,Int}[]   # (location, resolved_code, nstands)
for (loc,c) in locs
  p = P(Int32(loc), Int32(0))
  FVSjl.sn_forkod_remap!(p)
  k = Int(p.user_forest_code); iregn = k ÷ 10000
  iregn != 8 && push!(bad, (loc, k, c))
end
sort!(bad, by=x->-x[3])
println("LOCATION codes still NOT resolving to region 8 (zero-volume risk): ", length(bad), "  (", sum((x->x[3]).(bad); init=0), " stands)")
for (loc,k,c) in bad[1:min(30,end)]; println("  LOCATION=$loc → $k  (iregn=$(k÷10000))  stands=$c"); end
# pre-fix: how many codes/stands would have had iregn!=8 WITHOUT the remap (raw loc/10000)
prebad = [(loc,c) for (loc,c) in locs if loc ÷ 10000 != 8]
println("PRE-FIX (raw loc÷10000): ", length(prebad), " codes / ", sum((x->x[2]).(prebad); init=0), " stands would have blanked vol_eq")
for (loc,c) in sort(prebad, by=x->-x[2])[1:min(12,end)]; println("  was-bad LOCATION=$loc  stands=$c"); end
