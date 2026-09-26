# =============================================================================
# mortality.jl (inlandempire) — IE mortality (ie/morts.f). IDENTICAL Hamilton model to KT (same RIP/RIPP/WKI/
# MORCON; IPDG/IPDG2/POT bit-identical to KT, PMSC[1:11]==KT). Differences: IE_MORT_PMSC(23), BAMAX=IE_BAMAXA
# (ie/sitset.f), IFOR = p.forest_idx (ie/forkod.f), bark = ie_bratio (ie/bratio.f). Mirrors mortality!(::Kootenai).
# =============================================================================
function mortality!(s::StandState, ::InlandEmpire; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return _clim_mort_empty!(s, fint)
    ba = p.basal_area
    itype = Int(p.habitat_input)
    # BAMAX: a user BAMAX (keyword/DB ⇒ LBAMAX) is used as given. Otherwise ie/sitset.f:82-84 seeds BAMAXA(ITYPE)
    # WITHOUT setting LBAMAX, so MORTS's CALL SDICAL(0,SDIMAX) (morts.f:193) re-derives BAMAX=XMAX·0.5454154·PMSDIU
    # from the BA-weighted SDIDEF (sdical.f:203-205). The 380→SDIDEF→BAMAX round trip is NOT the identity in single
    # precision (IE habitat 13: 380.00003), and BAMAX divides RIPP (morts.f:291) — using the raw table value moved a
    # kill by ~5 ULP (bare-plot fixture record 241, cycle 4).
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : _ie_sdical_bamax(s, itype)
    bamax <= 0f0 && (bamax = 1f0)
    sdimax = clim_sdical_xmax(s, stand_sdimax(s), fint)   # ie/morts.f:193 SDICAL ⇒ CLMAXDEN-adjusted (SDIMAX<5 kill-all)
    # stand sums (morts.f:196-210, RAW record order): T, SD2SQ → DQ10; AVED (BA-weighted mean DBH)
    tt = 0f0; sd2sq = 0f0; dsum = 0f0; wprob = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = ie_bratio(sp, d)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + (2f0 * d * g + g * g)); tt += pr   # morts.f:198-199 CIOBDS=(2·D·G+G·G); SD2SQ+P·(D·D+CIOBDS)
        wprob += pr; dsum += d * pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt)
    deltba = 0.005454154f0 * dq10 * dq10 * tt - ba
    ba10 = ba + (bamax - ba) / bamax * deltba
    tb = ba10 / (0.005454154f0 * dq10 * dq10)
    ttb = (tt - tb) / tt; ttb > 0.9999f0 && (ttb = 0.9999f0)
    rz = 1f0 - fpow(1f0 - ttb, 0.1f0)                    # morts.f:207 (1−TTB)**0.1 = powf
    aved = dsum / wprob
    # MORCON: POTEN → GMULT/REIN per size class (morts.f:646-667)
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 11) && (ifor = 8)
    it = (1 <= itype <= 30) ? itype : 1
    poten1 = IE_MORT_POT[IE_MORT_IPDG[it, ifor]]
    poten2 = IE_MORT_POT[IE_MORT_IPDG2[it, ifor]]
    gmult1 = 0.90f0 / poten1; rein1 = (1f0 - fpow(poten1 / 20f0 + 1f0, -1.605f0)) / 0.06821f0
    gmult2 = 2.50f0 / poten2; rein2 = (1f0 - fpow(poten2 + 1f0, -1.605f0)) / 0.86610f0
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
        rip = 1f0 / (1f0 + fexp(rip))                         # morts.f:282 EXP = expf
        rip = rip * (ip == 1 ? rein1 : rein2)                 # ·POTENT
        ripp = ba * rz
        ba <= bamax && (ripp += (bamax - ba) * rip)
        ripp /= bamax
        ripp < rip && (ripp = rip); ripp > 1f0 && (ripp = 1f0)
        # ie/morts.f:316-322 species-group rate: NI rate for sp≤12,14,23; 20% for PI/JU
        # (sp15,16); 60% for LM,PY,AS,CO,MM,PB,OH (sp13,17,18,19,20,21,22). X=1 (no MORTMULT).
        smult = (sp <= 12 || sp == 14 || sp == 23) ? 1f0 : (sp == 15 || sp == 16) ? 0.2f0 : 0.6f0
        # ie/morts.f:301-311 establishment "best" trees are immune for 20 yr after the disturbance date: clear IESTAT
        # once IY(ICYC) reaches it, else X·(1−XCHECK) with XCHECK = clamp((IESTAT−IY(ICYC))/FINT, 0, 1).
        xest = 1f0
        if t.iestat[i] > 0
            iyc = Int32(current_cycle_year(s))
            iyc >= t.iestat[i] && (t.iestat[i] = Int32(0))
            xchk = Float32(t.iestat[i] - iyc) / fint
            xchk = clamp(xchk, 0f0, 1f0)
            xest = 1f0 - xchk
        end
        # morts.f:315-322: WKI=P·(1−(1−RIPP)**FINT)·X [·0.2 | ·0.6] — X before the species factor; **FINT = powf
        wki = pr * (1f0 - fpow(1f0 - ripp, fint)) * xest * smult
        gsc = (dgi / bark) * (fint / 10f0)
        if (dd + gsc) >= sc[sp, 1] && trunc(Int, sc[sp, 3]) != 1   # morts.f:326 (D+G) with the D≤0.5→0.5 clamped D
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
    # Dwarf-mistletoe mortality (mismrt.f): MAX-combine per-tree DM kill into killed[] before snags/removal,
    # same as the shared N-Rockies path (southern/mortality.jl). This variant has its own mortality! so it
    # must be wired here; inert on stands with no DM ratings (dmr==0 ⇒ per-tree no-op).
    # FIXMORT (morts.f:781): forced-mortality override applied AFTER the BA-check, before the DM combine —
    # same as southern/mortality.jl:499. This variant has its own mortality! so it must be wired here;
    # inert unless a FIXMORT keyword scheduled events (apply_fixmort! returns early on empty).
    apply_fixmort!(s, killed, n, fint)
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end

"""
    _ie_sdical_bamax(s, itype) -> Float32

SDICAL's non-LBAMAX BAMAX (sdical.f:95-126, 203-205), evaluated exactly as gfortran does: species BA sums over
IND1 (`BAXSP(sp) += ((0.0054542·D)·D)·P`, TOTBA in the same order), `XMAX = (Σ_{sp=1..MAXSP} SDIDEF(sp)·BAXSP(sp))/TOTBA`
(XMAX=1 when TOTBA≤0), then `BAMAX = (XMAX·0.5454154)·PMSDIU` with PMSDIU as a fraction (morts.f:184).
"""
function _ie_sdical_bamax(s::StandState, itype::Int)::Float32
    t = s.trees; p = s.plot
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    baxsp = zeros(Float32, MAXSP); totba = 0f0
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        for k in i1:isct[sp, 2]
            (1 <= k <= length(ind1)) || continue
            i = Int(ind1[k]); (1 <= i <= t.n) || continue
            treeba = 0.0054542f0 * t.dbh[i] * t.dbh[i] * t.tpa[i]
            baxsp[Int(t.species[i])] += treeba
            totba += treeba
        end
    end
    xmax = 0f0
    if totba <= 0f0
        xmax = 1f0
    else
        @inbounds for sp in 1:min(MAXSP, length(p.sp_sdi_def))   # DO 60 I=1,MAXSP (unused species add SDIDEF·0)
            xmax += p.sp_sdi_def[sp] * baxsp[sp]
        end
        xmax = xmax / totba
    end
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    return xmax * 0.5454154f0 * pmsdiu
end
