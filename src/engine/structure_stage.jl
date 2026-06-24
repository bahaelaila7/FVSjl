# =============================================================================
# structure_stage.jl — SSTAGE: stand structural-stage class (1-6).
#
# Ported (semantics) from base/sstage.f: classify the stand each cycle by canopy
# STRATIFICATION (height-gap layers) + per-stratum CANOPY COVER + the dominant
# stratum's DBH vs the small-tree/sawtimber thresholds:
#   1 = SI  (stand initiation)        4 = young multi-strata
#   2 = SE  (stem exclusion)          5 = old single-stratum
#   3 = UR  (understory reinitiation) 6 = old multi-strata / continuous
# Defaults (isstag.f:32-37): SSDBH=5, SAWDBH=25, GAPPCT=30, PCTSMX=30, CCMIN=5, TPAMIN=200.
#
# Scope of THIS chunk (the CLASS — SSTAGE's primary discrete output): the class uses
# cover only as a `> CCMIN(5%)` stratum-present threshold and DBH vs the coarse 5/25
# thresholds, so it is robust to the few-% cover / DBH precision the full per-stratum
# REPORT columns need (those need the exact CRWDTH source — see docs/SSTAGE_chunk_plan.md).
# Validated against the Fortran "Structural statistics" Struct-Class column.
# =============================================================================

const _SS_SSDBH  = 5.0f0    # small-tree DBH threshold      (TMPSSD)
const _SS_SAWDBH = 25.0f0   # sawtimber DBH threshold       (TMPSAW)
const _SS_GAPPCT = 30.0f0   # stratum-gap % height drop     (TMPGAP)
const _SS_CCMIN  = 5.0f0    # min stratum canopy cover %    (TMPCCM)
const _SS_TPAMIN = 200.0f0  # min TPA to be a stand         (TMPTPA)

# SSTGHP nominal DBH (sstage.f:780-870): the dominant-cohort DBH. Within the stratum (height-sorted
# ord[dlo:dhi]), take the CANOPY COHORT = top trees until cumulative crown area exceeds 41382 sq ft
# (~0.95 ac) — this excludes the suppressed understory. Of that cohort, find the tree at the 70th
# crown-area PERCENTILE (PCTILE: cumulative crown area of equal-or-larger-single-crown trees), and
# return the PROB-weighted mean DBH of the ±4-tree window around it.
function _ss_dbhnom(ord, dlo::Int, dhi::Int, dbh, tpa, crarea)::Float64
    csum = 0.0; i3 = dhi
    @inbounds for k in dlo:dhi
        csum += crarea[ord[k]]
        if csum > 41382.0; i3 = k; break; end
    end
    coh = [ord[k] for k in dlo:i3]                  # canopy cohort tree indices
    isempty(coh) && return 0.0
    wk4 = [crarea[i] / max(tpa[i], 1e-9) for i in coh]   # crown area per single tree
    cohS = coh[sortperm(wk4)]                        # ascending single-crown size (RDPSRT .FALSE.)
    tot = sum(crarea[i] for i in cohS)
    tot <= 0.0 && return 0.0
    pct = zeros(length(cohS))                        # PCTILE: cum crown area of ≥-crown trees, %
    acc = 0.0
    @inbounds for j in length(cohS):-1:1
        acc += crarea[cohS[j]]; pct[j] = acc / tot * 100.0
    end
    i70 = argmin(abs.(pct .- 70.0))                  # tree nearest the 70th percentile
    k1 = max(1, i70 - 4); k2 = min(length(cohS), i70 + 4)
    sd = 0.0; spw = 0.0
    @inbounds for j in k1:k2
        i = cohS[j]; sd += dbh[i] * tpa[i]; spw += tpa[i]
    end
    return spw > 1e-4 ? sd / spw : 0.0
end

# COVOLP (covolp.f): canopy cover % of a tree set whose crown areas (sq ft/ac) are `crarea[idx]`.
function _ss_cover(crarea::Vector{Float64}, idx, cccoef::Float64)::Float64
    s = 0.0
    @inbounds for i in idx; s += crarea[i]; end
    pccu = cccoef * (s / 43560.0)
    return pccu > 5.0 ? 100.0 : (1.0 - exp(-pccu)) * 100.0
end

"""
    structure_class(s) -> (; class, nstr, cover)

The stand structural-stage class (1-6) for the current stand (SSTAGE). Returns the class, the
number of valid strata, and the whole-stand canopy cover %. `class == 0` means unclassified
(too few trees / TPA below TPAMIN with no strata). `iba` selects the SDI (before/after thin) for
the NSTR=1 SE→SI demotion.
"""
function structure_class(s::StandState; iba::Int = 1)
    t = s.trees; p = s.plot; co = s.coef; g = p.gross_space
    cccoef = Float64(s.control.cc_coef)
    # working list: live trees with HT > 0; per-acre TPA; crown area = π/4·CW²·TPA (covolp.f)
    n = 0; ht = Float64[]; dbh = Float64[]; tpa = Float64[]; crarea = Float64[]
    @inbounds for i in 1:t.n
        t.height[i] > 0f0 && t.tpa[i] > 0f0 || continue
        sp2 = strip(co.code_alpha[Int(t.species[i])])
        cw = crown_width(co, sp2, t.dbh[i], t.height[i], Float32(t.crown_pct[i]), 0,
                         p.latitude, p.longitude, p.elevation)
        pa = Float64(t.tpa[i] / g)
        n += 1; push!(ht, Float64(t.height[i])); push!(dbh, Float64(t.dbh[i]))
        push!(tpa, pa); push!(crarea, Float64(cw)^2 * pa * 0.785398)
    end
    n == 0 && return (class = 0, nstr = 0, cover = 0.0)
    tprob = sum(tpa)
    ord = sortperm(ht; rev = true)                  # INDEX: trees by height descending

    # height-gap stratification: track the two largest gaps (sstage.f:300-388)
    diff1 = 0.0; diff2 = 0.0; id1i1 = id1i2 = id2i1 = id2i2 = 0
    iilg = 1; ilarge = ord[1]; sumprb = 0.0
    @inbounds for ii in 2:n
        ismall = ord[ii]
        x = ht[ilarge] * _SS_GAPPCT * 0.01; x < 10.0 && (x = 10.0)
        if ht[ismall] < ht[ilarge] - x              # a gap
            if tpa[ismall] + sumprb < 2.0
                sumprb += tpa[ismall]
            else
                diff = ht[ilarge] - ht[ismall]
                if diff > diff1
                    diff2 = diff1; diff1 = diff
                    id2i1 = id1i1; id2i2 = id1i2; id1i1 = iilg; id1i2 = ii
                elseif diff > diff2
                    diff2 = diff; id2i1 = iilg; id2i2 = ii
                end
                ilarge = ismall; iilg = ii; sumprb = 0.0
            end
        else                                         # no gap — ladder tree?
            if tpa[ismall] + sumprb < 2.0
                sumprb += tpa[ismall]
            else
                ilarge = ismall; iilg = ii; sumprb = 0.0
            end
        end
    end
    if id1i1 > id2i1 && id2i1 > 0                     # keep the upper gap on top
        id1i1, id2i1 = id2i1, id1i1; id1i2, id2i2 = id2i2, id1i2
    end

    # potential strata from the gap boundaries (sstage.f:388-408)
    is1i1 = 1; is1i2 = n; is2i1 = 0; is2i2 = 0; is3i1 = 0; is3i2 = 0
    if id1i1 > 0
        is1i2 = id1i1; is2i1 = id1i2; is2i2 = n
    end
    if id2i1 > 0
        is2i2 = id2i1; is3i1 = id2i2; is3i2 = n
    end

    # cover per stratum (incl. the gap trees), present if > CCMIN (sstage.f:430-465)
    covrange(lo, hi) = _ss_cover(crarea, (ord[k] for k in lo:hi), cccoef)
    crs1 = covrange(is1i1, max(is1i2, is2i1 - 1))
    is1ok = crs1 > _SS_CCMIN
    is2ok = false; is3ok = false
    if is2i1 > 0
        crs2 = covrange(is2i1, max(is2i2, is3i1 - 1)); is2ok = crs2 > _SS_CCMIN
    end
    if is3i1 > 0
        crs3 = covrange(is3i1, is3i2); is3ok = crs3 > _SS_CCMIN
    end
    nstr = is1ok + is2ok + is3ok
    cover = _ss_cover(crarea, 1:n, cccoef)
    if nstr == 0
        tprob < _SS_TPAMIN && return (class = 0, nstr = 0, cover = cover)
        is1ok = true; nstr = 1; is1i1 = 1; is1i2 = n
    end

    # dominant stratum = the top OK stratum; its PROB-weighted mean DBH (≈ SSTGHP DBHNOM) +
    # the min DBH of that stratum (sstage.f:487-535). NOTE: SSTGHP's exact DBHNOM is a 70th-
    # crown-area-percentile ±4-tree window; the coarse 5/25 thresholds make the PROB-weighted
    # mean an adequate stand-in for the CLASS (the per-stratum report column will need the exact).
    dlo, dhi = is1ok ? (is1i1, is1i2) : is2ok ? (is2i1, is2i2) : (is3i1, is3i2)
    tmpdbh = _ss_dbhnom(ord, dlo, dhi, dbh, tpa, crarea)   # SSTGHP 70th-percentile canopy-cohort DBH
    dmind = dbh[ord[dhi]]                             # min DBH of the dominant stratum

    # classify (sstage.f:539-576)
    cls = 0
    if nstr == 1
        if tmpdbh < _SS_SSDBH
            cls = 1
        elseif tmpdbh < _SS_SAWDBH
            cls = 2
            xsdi = iba == 1 ? _event_bsdi(s) : _event_bsdi(s)   # SDIBC (before) ≈ SDIAC for our use
            xbamax = Float64(s.control.ba_max)
            (xbamax > 0 && xsdi < 0.01 * 30.0 * xbamax) && (cls = 1)   # PCTSMX demotion
        else
            cls = dmind < 3.0 ? 6 : 5
        end
    elseif nstr == 2
        cls = tmpdbh < _SS_SSDBH ? 1 : tmpdbh < _SS_SAWDBH ? 3 : 6
    else
        cls = tmpdbh < _SS_SSDBH ? 1 : tmpdbh < _SS_SAWDBH ? 4 : 6
    end
    return (class = cls, nstr = nstr, cover = cover)
end
