using FVSjl
const F = FVSjl
V = F.variant_from_code("IE")
key = "clim_iet.key"
for s in F.each_stand(key; variant=V)
    F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
    c = s.climate
    if c === nothing || !c.active
        println("climate INACTIVE for ", strip(s.plot.stand_id)); break
    end
    t = s.trees
    cd = c.data; ix = c.indices
    ty = Float32(F.current_cycle_year(s))   # cyc0 year (1990)
    A(sym, yr) = F.algslp(yr, cd.years, view(cd.attrs, :, ix[sym]))
    smi(yr) = (g = A(:gsp, yr); g > 0f0 ? A(:dd5, yr) / g : 0f0)
    xgsite = ix[:pSite] > 0 ? F.clim_xgsite(A(:pSite, ty), A(:pSite, Float32(c.inv_year))) : 1f0
    mtcm_now = A(:mtcm, ty); mmin_now = A(:mmin, ty); smi_now = smi(ty)
    ns = length(c.plant_symbols)
    # per-species accumulators
    spba = zeros(Float32, ns); sptpa = zeros(Float32, ns)
    spgm_num = zeros(Float32, ns); spwt = zeros(Float32, ns)
    vscore = ones(Float32, ns); spviab = ones(Float32, ns)
    for sp in 1:ns
        spviab[sp], vscore[sp] = F.species_vscore(cd, c.plant_symbols[sp], ty)
    end
    @inbounds for i in 1:t.n
        t.dbh[i] <= 0f0 && continue
        sp = Int(t.species[i]); (sp < 1 || sp > ns) && continue
        prob = t.tpa[i]; d = t.dbh[i]
        spba[sp]  += d*d*prob*0.005454154f0
        sptpa[sp] += prob
        birthyr = ty - t.birth_age[i]
        xdf = F.leites_xdf(mtcm_now, A(:mtcm, birthyr))
        xwl = F.leites_xwl(mmin_now, A(:mmin, birthyr), A(:dd0, birthyr))
        xpp = F.leites_xpp(smi_now, smi(birthyr), A(:d100, birthyr))
        xr  = F.clim_xrelgr(c.plant_symbols[sp], xdf, xpp, xwl)
        _, tm = F.clim_treemult(xgsite, xr, vscore[sp], c.growmult[sp])
        spgm_num[sp] += tm*prob; spwt[sp] += prob
    end
    println("=== jl cyc0 (", Int(ty), ") per-species climate report, stand ", strip(s.plot.stand_id), " ===")
    for sp in 1:ns
        (spwt[sp] > 0 || spviab[sp] > 0.4f0) || continue
        spgm = spwt[sp] > 0 ? spgm_num[sp]/spwt[sp] : 1f0
        sym = c.plant_symbols[sp]
        println("  sp#", sp, " ", sym, " viab=", round(spviab[sp],digits=4),
                " BA=", round(spba[sp],digits=2), " TPA=", round(sptpa[sp],digits=2),
                " gm=", round(spgm,digits=4), " sm=", round(xgsite^c.growmult[sp],digits=4))
    end
    break
end
