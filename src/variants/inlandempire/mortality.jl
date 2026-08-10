# =============================================================================
# mortality.jl (inlandempire) — IE mortality (ie/morts.f). IDENTICAL Hamilton model to KT (same RIP/RIPP/WKI/
# MORCON; IPDG/IPDG2/POT bit-identical to KT, PMSC[1:11]==KT). Differences: IE_MORT_PMSC(23), BAMAX=IE_BAMAXA
# (ie/sitset.f), IFOR = p.forest_idx (ie/forkod.f), bark = ie_bratio (ie/bratio.f). Mirrors mortality!(::Kootenai).
# =============================================================================
function mortality!(s::StandState, ::InlandEmpire; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    ba = p.basal_area
    itype = Int(p.habitat_input)
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : ((1 <= itype <= 30) ? IE_BAMAXA[itype] : 0f0)
    bamax <= 0f0 && (bamax = 1f0)
    sdimax = stand_sdimax(s)
    # stand sums (morts.f:196-210, RAW record order): T, SD2SQ → DQ10; AVED (BA-weighted mean DBH)
    tt = 0f0; sd2sq = 0f0; dsum = 0f0; wprob = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = ie_bratio(sp, d)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + 2f0 * d * g + g * g); tt += pr
        wprob += pr; dsum += d * pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt)
    deltba = 0.005454154f0 * dq10 * dq10 * tt - ba
    ba10 = ba + (bamax - ba) / bamax * deltba
    tb = ba10 / (0.005454154f0 * dq10 * dq10)
    ttb = (tt - tb) / tt; ttb > 0.9999f0 && (ttb = 0.9999f0)
    rz = 1f0 - (1f0 - ttb)^0.1f0
    aved = dsum / wprob
    # MORCON: POTEN → GMULT/REIN per size class (morts.f:646-667)
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 11) && (ifor = 8)
    it = (1 <= itype <= 30) ? itype : 1
    poten1 = IE_MORT_POT[IE_MORT_IPDG[it, ifor]]
    poten2 = IE_MORT_POT[IE_MORT_IPDG2[it, ifor]]
    gmult1 = 0.90f0 / poten1; rein1 = (1f0 - (poten1 / 20f0 + 1f0)^(-1.605f0)) / 0.06821f0
    gmult2 = 2.50f0 / poten2; rein2 = (1f0 - (poten2 + 1f0)^(-1.605f0)) / 0.86610f0
    sqba = sqrt(ba)
    icyc1 = Int(s.control.cycle) == 0
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    sc = s.control.sp_size_cap
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]; bark = ie_bratio(sp, d)
        reldbh = d / aved
        dd = d <= 0.5f0 ? 0.5f0 : d
        dgi = t.diam_growth[i]
        ip = d <= 5f0 ? 2 : 1
        gmult = ip == 1 ? gmult1 : gmult2
        wk1 = t.dg_prev[i]; oldfnt = 10f0                     # WK1 = previous cycle's applied DG
        dgt = wk1 / oldfnt
        d <= 1f0 && dgt < 0.05f0 && (dgt = 0.05f0)
        (1f0 < d <= 5f0) && dgt < 0.05f0 && (dgt = 0.05f0 * (5f0 - d) / 4f0)
        g = wk1 / (bark * oldfnt)
        wk1 / oldfnt < dgt && (g = dgt / bark)
        (icyc1 || wk1 == 0f0) && dgi > 0.5f0 && (g = dgi / (bark * 10f0))
        g = g * gmult
        rip = 2.76253f0 + 0.222310f0 * sqrt(dd) - 0.0460508f0 * sqba + 11.2007f0 * g -
              0.554421f0 / dd + IE_MORT_PMSC[sp] + 0.246301f0 * reldbh + 6.07129f0 * g / dd
        rip > 70f0 && (rip = 70f0); rip < -70f0 && (rip = -70f0)
        rip = 1f0 / (1f0 + exp(rip))
        rip = rip * (ip == 1 ? rein1 : rein2)                 # ·POTENT
        ripp = ba * rz
        ba <= bamax && (ripp += (bamax - ba) * rip)
        ripp /= bamax
        ripp < rip && (ripp = rip); ripp > 1f0 && (ripp = 1f0)
        # ie/morts.f:316-322 species-group rate: NI rate for sp≤12,14,23; 20% for PI/JU
        # (sp15,16); 60% for LM,PY,AS,CO,MM,PB,OH (sp13,17,18,19,20,21,22). X=1 (no MORTMULT).
        smult = (sp <= 12 || sp == 14 || sp == 23) ? 1f0 : (sp == 15 || sp == 16) ? 0.2f0 : 0.6f0
        wki = pr * (1f0 - (1f0 - ripp)^fint) * smult
        gsc = (dgi / bark) * (fint / 10f0)
        if (d + gsc) >= sc[sp, 1] && trunc(Int, sc[sp, 3]) != 1
            wki = max(wki, pr * sc[sp, 2] * fint / 10f0)
        end
        wki > pr && (wki = pr)
        sdimax < 5f0 && (wki = pr)
        killed[i] = wki
    end
    # Climate-FVS mortality (clmorts.f:369 CALL, after base mort, before booking): viability path, THISYR mid-cycle
    # (clmorts.f:78 THISYR=IY(ICYC)+FINT/2). Inert unless a CLIMATE keyword activated s.climate.
    (s.climate !== nothing && s.climate.active) &&
        apply_climate_mort!(s, killed, Float32(current_cycle_year(s)) + fint / 2f0, fint)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
