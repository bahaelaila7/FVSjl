# =============================================================================
# diameter_growth.jl (britishcolumbia) — BC large-tree DG engine hooks (canada/bc/dgf.f, V3 metric).
#
# Two hooks (KT/LS/IE pattern):
#   bc_dgcons!(s)             — per-species DGCON setup (dgf.f ENTRY DGCONS): resolve IP/JP by BEC
#                               PrettyName, then DGCON = SSKONST%CON + ZNKONST%CON + elev/aspect/slope.
#   dgf!(s, ::BritishColumbia)— per-tree WK2 = DDS (dgf.f main body, V3 metric path).
#
# Coefficients/matcher/DDS/bark in dg_coefficients.jl (VALIDATED: matcher 15/15, DDS 97.1% bit-exact,
# DGCON = oracle CONSPP on all_BC). Wiring MEASURED: BAL uses FVS PCT = jl `t.crown_ratio` (trees.jl:51
# "(PCT)"); CR = `t.crown_pct`·0.01. ⚠ V2 regime (non ICH/IDF/SBS/SBPS zones) NOT yet ported.
#
# ⚠ PENDING full engine-run validation (needs BEC zone on StandState + RELDEN/BA/PCT from engine).
# TEMP: bc_stand_zone hardcodes all_BC (ICH/ICHmw2/01) until chunk-2 becset habitat→zone is ported.
# =============================================================================

"""Stand BEC zone + site-series. TEMP hardcode (all_BC) until chunk-2 (becset.f) populates StandState."""
function bc_stand_zone(s::StandState)
    # TODO chunk-2: derive from becset habitat→zone map (STDINFO 231Dd → ICH/ICHmw2/01 for all_BC).
    return ("ICH", "ICHmw2/01")
end

"""BC DGCONS — resolve per-species DGCON once/stand (canada/bc/dgf.f:2185-2230)."""
function bc_dgcons!(s::StandState)
    c = s.calib; p = s.plot
    zone, series = bc_stand_zone(s)
    elev = p.elevation; asp = p.aspect; slope = p.slope
    @inbounds for sp in 1:nspecies(BritishColumbia())
        ip, jp = bc_resolve_ipjp(sp, series, zone)
        if ip < 1 || jp < 1
            c.dg_const[sp] = 0f0
            continue
        end
        pelev_ft = sp == 1 ? clamp(Float32(elev) * BC_FTtoM, 5f0, 12f0) / BC_FTtoM : Float32(elev)
        dgcon = bc_dgcon(ip, jp, pelev_ft, asp, slope)
        (s.control.dg_cor2_on && c.dg_cor2[sp] > 0f0) && (dgcon += log(c.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
    end
    return s
end

"""BC `dgf!` hook — fill scratch.wk[2,i] with per-tree DDS (V3 metric, canada/bc/dgf.f main body)."""
function dgf!(s::StandState, ::BritishColumbia)
    p, t, c = s.plot, s.trees, s.calib
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density; ba = p.basal_area
    zone, series = bc_stand_zone(s)
    inICH = occursin("ICH", zone); inIDF = occursin("IDF", zone)
    nsp = nspecies(BritishColumbia())
    ip = zeros(Int, nsp); dgccf1 = zeros(Float32, nsp)
    @inbounds for sp in 1:nsp
        i, _ = bc_resolve_ipjp(sp, series, zone)
        ip[sp] = i
        i > 0 && (dgccf1[sp] = BC_ZNKONST[i].CCFA)
    end
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = Int(t.species[i])
        (sp < 1 || sp > nsp || ip[sp] < 1) && continue
        # RELDN2: per-species CCF floors (dgf.f:1888-1903)
        reldn2 = relden
        (sp == 5 && inICH) && (reldn2 = max(100f0, relden))
        (sp == 7 && inIDF) && (reldn2 = max(100f0, relden))
        (sp == 9 && inICH) && (reldn2 = max(125f0, relden))
        ((sp == 11 || sp == 15) && (inICH || inIDF)) && (reldn2 = max(100f0, relden))
        conspp = c.dg_const[sp] + c.dg_cor[sp] + dgccf1[sp] * reldn2 * 0.01f0
        pct = t.crown_ratio[i]                                 # FVS PCT (trees.jl:51)
        bal = (1f0 - pct / 100f0) * ba * BC_FT2pACRtoM2pHA
        cr  = Float32(t.crown_pct[i]) * 0.01f0
        # crown caps in ICH/IDF for HW/BL/EP-OH (dgf.f:1925-1935)
        ((sp == 5 && inICH) || ((sp == 11 || sp == 15) && (inICH || inIDF))) && (cr = min(0.8f0, cr))
        wk2[i] = bc_v3_dds(sp, ip[sp], zone, d, bal, cr, conspp, bc_bratio(sp))
    end
    return s
end
