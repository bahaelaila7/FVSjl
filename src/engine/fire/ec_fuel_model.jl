# =============================================================================
# fire/ec_fuel_model.jl — EC (EastCascades) FMDYN dynamic cover-metagroup fuel-model selection.
#
# Ported from ec/fmcfmd.f (VARACD 'EC' branch). This is NEITHER the California-CWHR classifier (NC/WS/CA/SO)
# NOR the WC/PN FIRE-VPN 6-cover-metagroup — it is EC's own "detailed low fuel model selection": the 32
# species are pooled into 15 cover-type metagroups (DF/LP/SAF/MH/WP/ES/PSF/WL/PP pure + DF-GF/PP-DF/LP-WL
# mixes + SAF-leading + moist/dry mixed conifer), a single dominant metagroup ICT is resolved (>50% BA in a
# 1-sp group, then a 2-sp group, then SAF-leading, then the moist/dry catch-all), and a set of ALGSLP
# weighting rules keyed by ICT, PERCOV, QMD and the single-vs-multi-strata flag LSNGL (from FMSSTAGE) builds
# the candidate fuel-model weights EQWT. Models 10/12/13 are always natural-fuel candidates; model 11 is the
# 5-yr post-activity fuel jump (deferred — no activity ⇒ AFWT=0). Resolved by `_fmdyn` over EC's 13-model
# XPTS (identical to WC's). Values verbatim from ec/fmcfmd.f + ec/fmcba.f (MAPDRY via ec_moist).
#
# NOTE on QMD (ec/fmcfmd.f:171-176): the Fortran loop `DO J=1,ITRN ... FMPROB(I)*DBH(I)*DBH(I)` indexes the
# leftover `I` (=ICLSS+1 after the EQWT-init loop), NOT `J` — a source bug that makes QMD collapse to a single
# record's DBH. It is non-portable (depends on FVS's internal record order) and immaterial to the FMD here:
# EC's QMD cutpoints are 2/4/7/9/19/21, and the proper stand QMD (and the buggy value) both land in the same
# ALGSLP bucket for ect01. We compute the INTENDED QMD (Σprob·D²/Σprob); documented cornered.
# =============================================================================

# ec/fmcfmd.f cover-metagroup enums (PARAMETER block, fmcfmd.f:69-95).
const _EC_DFCT   = 1;  const _EC_LPCT   = 2;  const _EC_SAFCT  = 3;  const _EC_MHCT   = 4
const _EC_WPCT   = 5;  const _EC_ESCT   = 6;  const _EC_PSFCT  = 7;  const _EC_WLCT   = 8
const _EC_PPCT   = 9;  const _EC_DFGFCT = 10; const _EC_PPDFCT = 11; const _EC_LPWLCT = 12
const _EC_SAFMCT = 13; const _EC_MMIXCT = 14; const _EC_DMIXCT = 15

"""
    ec_select_fuel_models(s, mois, sm, lg) -> Vector{(model, weight)}

EC FMDYN dynamic cover-metagroup fuel-model selection (ec/fmcfmd.f). `sm`/`lg` are the SMALL/LARGE down-wood
loads (the FMDYN point); PERCOV comes from the FireState (fmcba). Source-faithful transcription.
"""
function ec_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees
    eqwt = zeros(Float32, 13)                                  # EC ICLSS = 13
    percov = s.fire.percov
    itype = Int(s.plot.habitat_input)

    # per-species FFE basal area (FMTBA, ec/fmcfmd.f:206-208 / fmcba.f:375).
    nsp = length(coef_col(s.coef, :dbh_min))
    fmtba = zeros(Float32, nsp)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i])
        fmtba[sp] += t.tpa[i] * t.dbh[i] * t.dbh[i] * 0.0054542f0
    end
    stndba = 0f0
    @inbounds for sp in 1:nsp; stndba += fmtba[sp]; end

    # QMD (intended; see header note on the ec/fmcfmd.f:171 source bug). Σprob·D² / Σprob over all trees.
    x1 = 0f0; y1 = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        x1 += t.tpa[i] * t.dbh[i] * t.dbh[i]; y1 += t.tpa[i]
    end
    qmd = y1 > 0f0 ? sqrt(x1 / y1) : 0f0

    # FMSSTAGE structure class → LSNGL (ec/fmcfmd.f:188-198): 0,1,2,5 single-stratum; 3,4,6 multi.
    ifmst = structure_class(s).class
    lsngl = (ifmst == 0 || ifmst == 1 || ifmst == 2 || ifmst == 5)

    # cover-metagroup basal areas CTBA[1..15] (ec/fmcfmd.f:210-230).
    ctba = zeros(Float32, 15)
    _g(k) = (1 <= k <= nsp) ? fmtba[k] : 0f0
    ctba[_EC_DFCT]   = _g(3)
    ctba[_EC_LPCT]   = _g(7) + _g(14)
    ctba[_EC_SAFCT]  = _g(9) + _g(15) + _g(16)
    ctba[_EC_MHCT]   = _g(12) + _g(11)
    ctba[_EC_WPCT]   = _g(1)
    ctba[_EC_ESCT]   = _g(8)
    ctba[_EC_PSFCT]  = _g(4)
    ctba[_EC_WLCT]   = _g(2) + _g(17)
    ctba[_EC_PPCT]   = _g(10)
    ctba[_EC_DFGFCT] = _g(3) + _g(6)
    ctba[_EC_PPDFCT] = _g(10) + _g(3)
    ctba[_EC_LPWLCT] = _g(7) + _g(14) + _g(2) + _g(17)
    # ECMOIST: moist-mixed (MMIXCT) iff MAPDRY==1, else dry-mixed (DMIXCT) (ec/fmcfmd.f:225-230).
    ec_moist(itype) == 1 ? (ctba[_EC_MMIXCT] = stndba) : (ctba[_EC_DMIXCT] = stndba)

    # resolve the dominant cover-metagroup ICT (ec/fmcfmd.f:239-278). Falls to OLDICT if no trees/BA.
    ict = 0
    if t.n > 0 && stndba > 0.001f0
        for i in _EC_DFCT:_EC_PPCT                              # (1) one-species dominant >50%
            if ctba[i] / stndba > 0.50f0; ict = i; @goto done; end
        end
        for i in _EC_DFGFCT:_EC_LPWLCT                          # (2) two-species dominant >50%
            if ctba[i] / stndba > 0.50f0; ict = i; @goto done; end
        end
        lsafd = true                                            # (3) SAF leading
        @inbounds for i in 1:nsp
            if i != 9 && fmtba[i] > ctba[_EC_SAFCT]; lsafd = false; break; end
        end
        if lsafd; ict = _EC_SAFMCT; @goto done; end
        for i in _EC_MMIXCT:_EC_DMIXCT                          # (4) moist/dry catch-all >50%
            if ctba[i] / stndba > 0.50f0; ict = i; @goto done; end
        end
    else
        ict = s.fire.covtyp_ict > 0 ? Int(s.fire.covtyp_ict) : _EC_LPCT   # OLDICT (FMVINIT default LPCT... EC=1 DFCT)
    end
    @label done
    s.fire.covtyp_ict = Int32(ict)                              # OLDICT persist

    alg(v, x1p, x2p) = _cr_algslp2(Float32(v), Float32(x1p), Float32(x2p), 0f0, 1f0)

    # ---- SELECT CASE (ICT) → EQWT (ec/fmcfmd.f:294-797) ----
    if ict == _EC_DFCT || ict == _EC_DFGFCT || ict == _EC_LPCT || ict == _EC_LPWLCT ||
       ict == _EC_MHCT || ict == _EC_PSFCT || ict == _EC_SAFCT || ict == _EC_WPCT ||
       ict == _EC_SAFMCT || ict == _EC_MMIXCT
        wt1 = zeros(Float32, 4)
        if percov < 30f0;      wt1[2] = alg(percov, 10, 30); wt1[1] = 1f0 - wt1[2]
        elseif percov < 60f0;  wt1[3] = alg(percov, 40, 60); wt1[2] = 1f0 - wt1[3]
        else;                  wt1[4] = alg(percov, 70, 90); wt1[3] = 1f0 - wt1[4]
        end
        if wt1[1] > 0f0
            if lsngl
                if ict == _EC_MMIXCT
                    w2 = zeros(Float32, 3)
                    if qmd <= 9f0; w2[2] = alg(qmd, 7, 9); w2[1] = 1f0 - w2[2]
                    else;          w2[3] = alg(qmd, 19, 21); w2[2] = 1f0 - w2[3]; end
                    eqwt[1] += wt1[1] * (w2[1] + w2[3]); eqwt[5] += wt1[1] * w2[2]
                else
                    eqwt[1] += wt1[1]
                end
            else
                eqwt[5] += wt1[1]
            end
        end
        if wt1[2] > 0f0
            if ict == _EC_LPCT || ict == _EC_LPWLCT || ict == _EC_WPCT || ict == _EC_DFCT
                w2b = alg(qmd, 2, 4); w2a = 1f0 - w2b
                if ict == _EC_DFCT
                    lsngl ? (eqwt[5] += wt1[2]) : (eqwt[8] += wt1[2])
                else
                    eqwt[1] += wt1[2] * w2a
                    lsngl ? (eqwt[5] += wt1[2] * w2b) : (eqwt[8] += wt1[2] * w2b)
                end
            elseif ict == _EC_DFGFCT || ict == _EC_MHCT || ict == _EC_PSFCT || ict == _EC_MMIXCT
                if lsngl
                    w2 = zeros(Float32, 3)
                    if qmd <= 4f0; w2[2] = alg(qmd, 2, 4); w2[1] = 1f0 - w2[2]
                    else;          w2[3] = alg(qmd, 19, 21); w2[2] = 1f0 - w2[3]; end
                    eqwt[5] += wt1[2] * (w2[1] + w2[3]); eqwt[8] += wt1[2] * w2[2]
                else
                    eqwt[8] += wt1[2]
                end
            elseif ict == _EC_SAFCT || ict == _EC_SAFMCT
                if lsngl
                    w2b = alg(qmd, 2, 4); w2a = 1f0 - w2b
                    eqwt[5] += wt1[2] * w2a; eqwt[8] += wt1[2] * w2b
                else
                    eqwt[8] += wt1[2]
                end
            end
        end
        if wt1[3] > 0f0
            w2b = alg(qmd, 2, 4); w2a = 1f0 - w2b
            eqwt[5] += wt1[3] * w2a; eqwt[8] += wt1[3] * w2b
        end
        if wt1[4] > 0f0
            if lsngl && (ict == _EC_LPCT || ict == _EC_WPCT || ict == _EC_LPWLCT)
                w2b = alg(qmd, 2, 4); w2a = 1f0 - w2b
                eqwt[5] += wt1[4] * w2a; eqwt[8] += wt1[4] * w2b
            else
                eqwt[8] += wt1[4]
            end
        end
    elseif ict == _EC_ESCT || ict == _EC_WLCT
        wt1 = zeros(Float32, 4)
        if qmd <= 4f0;     wt1[2] = alg(qmd, 2, 4);  wt1[1] = 1f0 - wt1[2]
        elseif qmd < 9f0;  wt1[3] = alg(qmd, 7, 9);  wt1[2] = 1f0 - wt1[3]
        else;              wt1[4] = alg(qmd, 19, 21); wt1[3] = 1f0 - wt1[4]
        end
        if wt1[1] > 0f0
            w2 = zeros(Float32, 3)
            if percov <= 30f0; w2[2] = alg(percov, 10, 30); w2[1] = 1f0 - w2[2]
            else;              w2[3] = alg(percov, 70, 90); w2[2] = 1f0 - w2[3]; end
            eqwt[1] += wt1[1] * w2[1]; eqwt[5] += wt1[1] * w2[2]; eqwt[8] += wt1[1] * w2[3]
        end
        if wt1[2] > 0f0
            w2b = alg(percov, 10, 30); w2a = 1f0 - w2b
            lsngl ? (eqwt[1] += wt1[2] * w2a) : (eqwt[5] += wt1[2] * w2a)
            eqwt[8] += wt1[2] * w2b
        end
        if wt1[3] > 0f0
            w2b = alg(percov, 10, 30); w2a = 1f0 - w2b
            if ict == _EC_ESCT || (ict == _EC_WLCT && lsngl)
                eqwt[5] += wt1[3] * w2a; eqwt[8] += wt1[3] * w2b
            else
                eqwt[8] += wt1[3]
            end
        end
        if wt1[4] > 0f0
            w2b = alg(percov, 10, 30); w2a = 1f0 - w2b
            if (ict == _EC_ESCT && lsngl) || (ict == _EC_WLCT && !lsngl)
                eqwt[5] += wt1[4] * w2a
            else
                eqwt[1] += wt1[4] * w2a
            end
            eqwt[8] += wt1[4] * w2b
        end
    elseif ict == _EC_PPCT
        wt1 = zeros(Float32, 4)
        if percov < 30f0;      wt1[2] = alg(percov, 10, 30); wt1[1] = 1f0 - wt1[2]
        elseif percov < 60f0;  wt1[3] = alg(percov, 40, 60); wt1[2] = 1f0 - wt1[3]
        else;                  wt1[4] = alg(percov, 70, 90); wt1[3] = 1f0 - wt1[4]
        end
        wt1[1] > 0f0 && (lsngl ? (eqwt[2] += wt1[1]) : (eqwt[6] += wt1[1]))
        wt1[2] > 0f0 && (lsngl ? (eqwt[6] += wt1[2]) : (eqwt[9] += wt1[2]))
        if wt1[3] > 0f0
            w2b = alg(qmd, 2, 4); w2a = 1f0 - w2b
            eqwt[6] += wt1[3] * w2a; eqwt[9] += wt1[3] * w2b
        end
        wt1[4] > 0f0 && (eqwt[9] += wt1[4])
    elseif ict == _EC_DMIXCT
        wt1 = zeros(Float32, 4)
        if percov < 30f0;      wt1[2] = alg(percov, 10, 30); wt1[1] = 1f0 - wt1[2]
        elseif percov < 60f0;  wt1[3] = alg(percov, 40, 60); wt1[2] = 1f0 - wt1[3]
        else;                  wt1[4] = alg(percov, 70, 90); wt1[3] = 1f0 - wt1[4]
        end
        wt1[1] > 0f0 && (lsngl ? (eqwt[2] += wt1[1]) : (eqwt[6] += wt1[1]))
        wt1[2] > 0f0 && (lsngl ? (eqwt[6] += wt1[2]) : (eqwt[8] += wt1[2]))
        if wt1[3] > 0f0
            w2 = zeros(Float32, 3)
            if qmd <= 4f0; w2[2] = alg(qmd, 2, 4); w2[1] = 1f0 - w2[2]
            else;          w2[3] = alg(qmd, 19, 21); w2[2] = 1f0 - w2[3]; end
            eqwt[6] += wt1[3] * w2[1]
            lsngl ? (eqwt[9] += wt1[3] * w2[2]) : (eqwt[8] += wt1[3] * w2[2])
            eqwt[8] += wt1[3] * w2[3]
        end
        wt1[4] > 0f0 && (eqwt[8] += wt1[4])
    elseif ict == _EC_PPDFCT
        wt1 = zeros(Float32, 4)
        if percov < 30f0;      wt1[2] = alg(percov, 10, 30); wt1[1] = 1f0 - wt1[2]
        elseif percov < 60f0;  wt1[3] = alg(percov, 40, 60); wt1[2] = 1f0 - wt1[3]
        else;                  wt1[4] = alg(percov, 70, 90); wt1[3] = 1f0 - wt1[4]
        end
        wt1[1] > 0f0 && (lsngl ? (eqwt[2] += wt1[1]) : (eqwt[6] += wt1[1]))
        if wt1[2] > 0f0
            if lsngl
                w2 = zeros(Float32, 3)
                if qmd <= 4f0; w2[2] = alg(qmd, 2, 4); w2[1] = 1f0 - w2[2]
                else;          w2[3] = alg(qmd, 19, 21); w2[2] = 1f0 - w2[3]; end
                eqwt[6] += wt1[2] * (w2[1] + w2[3]); eqwt[8] += wt1[2] * w2[2]
            else
                eqwt[8] += wt1[2]
            end
        end
        if wt1[3] > 0f0
            if lsngl
                w2b = alg(qmd, 2, 4); w2a = 1f0 - w2b
                eqwt[6] += wt1[3] * w2a; eqwt[9] += wt1[3] * w2b
            else
                eqwt[8] += wt1[3]
            end
        end
        if wt1[4] > 0f0
            if lsngl
                w2b = alg(qmd, 2, 4); w2a = 1f0 - w2b
                eqwt[8] += wt1[3] * w2a; eqwt[9] += wt1[3] * w2b   # ec/fmcfmd.f:789 uses WT1(3) (source bug — replicated)
            else
                eqwt[8] += wt1[4]
            end
        end
    end

    # post-activity model 11 (AFWT/LATFUEL) — deferred (no activity ⇒ AFWT=0); models 10/12/13 natural fuels.
    eqwt[10] = 1f0
    eqwt[12] = 1f0
    eqwt[13] = 1f0

    return _fmdyn(sm, lg, eqwt, _WC_FMD_XPTS)   # EC XPTS == WC XPTS (5/15…45/100, model 10=(15,30))
end
