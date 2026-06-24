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
# Validated bit-exact vs the Fortran "Structural statistics" report:
#   • the CLASS column (1-6) — every cycle, fire_early + snt01 stand-1;
#   • the Tot-Cov column — every cycle (snt01 10/11 ± a single ULP round, fire stands pre-fire;
#     post-fire diverges only by the known fire kill residual). NB: the crown AREA uses the RAW
#     PROB (`t.tpa`, NOT /GROSPC) — that was the cover fix, not the crown width (which already
#     matched Fortran's CrWidth exactly). CCCOEF=1 default.
#   • the uppermost-stratum DBH (`strdbh`, SSTGHP DBHNOM / the BSTRDBH event var) — 8/11 cycles
#     exact, the rest ≤0.5" (cohort/window-edge boundary). The fix was the WK4 sort DIRECTION:
#     `RDPSRT(.FALSE.)` is DESCENDING (biggest crown first), so the "70th percentile = 30% down from
#     the top" lands on the upper canopy, not mid-cohort.
# =============================================================================

# Default SSTAGE thresholds (isstag.f:32-37), in the `control.strclass_thresh` order:
#   gappct=30 (stratum-gap % height drop), ssdbh=5 (small-tree DBH), sawdbh=25 (sawtimber DBH),
#   ccmin=5 (min stratum cover %), tpamin=200 (min TPA), pctsmx=30 (% MaxSDI for SE).
const SS_THRESH_DEFAULT = (30.0f0, 5.0f0, 25.0f0, 5.0f0, 200.0f0, 30.0f0)

# SSTGHP nominal DBH (sstage.f:780-870): the dominant-cohort DBH. Within the stratum (height-sorted
# ord[dlo:dhi]), take the CANOPY COHORT = top trees until cumulative crown area exceeds 41382 sq ft
# (~0.95 ac) — this excludes the suppressed understory. Of that cohort, find the tree at the 70th
# crown-area PERCENTILE (PCTILE: cumulative crown area of equal-or-larger-single-crown trees), and
# return the PROB-weighted mean DBH of the ±4-tree window around it.
function _ss_dbhnom(ord, dlo::Int, dhi::Int, ht, dbh, tpa, crarea)
    csum = 0.0; i3 = dhi
    @inbounds for k in dlo:dhi
        csum += crarea[ord[k]]
        if csum > 41382.0; i3 = k; break; end
    end
    coh = [ord[k] for k in dlo:i3]                  # canopy cohort tree indices
    isempty(coh) && return (0.0, 0.0)
    wk4 = [crarea[i] / max(tpa[i], 1e-9) for i in coh]   # crown area per single tree
    cohS = coh[sortperm(wk4; rev = true)]            # DESCENDING single-crown size (RDPSRT .FALSE.)
    tot = sum(crarea[i] for i in cohS)
    tot <= 0.0 && return (0.0, 0.0)
    pct = zeros(length(cohS))                        # PCTILE: cum crown area of this-or-larger trees, %
    acc = 0.0
    @inbounds for j in length(cohS):-1:1
        acc += crarea[cohS[j]]; pct[j] = acc / tot * 100.0
    end
    # 70th percentile = "30% down from the top" (the big-crown end) — i70 nearest pct 70 (sstage.f:838)
    i70 = argmin(abs.(pct .- 70.0))
    k1 = max(1, i70 - 4); k2 = min(length(cohS), i70 + 4)
    sd = 0.0; sh = 0.0; spw = 0.0
    @inbounds for j in k1:k2
        i = cohS[j]; sd += dbh[i] * tpa[i]; sh += ht[i] * tpa[i]; spw += tpa[i]
    end
    return spw > 1e-4 ? (sd / spw, sh / spw) : (0.0, 0.0)   # (DBHNOM, nominal height)
end

# COVOLP (covolp.f): canopy cover % of a tree set whose crown areas (sq ft/ac) are `crarea[idx]`.
function _ss_cover(crarea::Vector{Float64}, idx, cccoef::Float64)::Float64
    s = 0.0
    @inbounds for i in idx; s += crarea[i]; end
    pccu = cccoef * (s / 43560.0)
    return pccu > 5.0 ? 100.0 : (1.0 - exp(-pccu)) * 100.0
end

"""
    structure_class(s) -> (; class, nstr, cover, strdbh)

The stand structural-stage class (1-6) for the current stand (SSTAGE). Returns the class, the
number of valid strata, the whole-stand canopy cover %, and the uppermost-stratum dominant DBH
(`strdbh`). `class == 0` means unclassified (too few trees / TPA below TPAMIN with no strata).
`iba` selects the SDI (before/after thin) for the NSTR=1 SE→SI demotion. The STRCLASS keyword can
override the thresholds (`control.strclass_thresh` = gappct/ssdbh/sawdbh/ccmin/tpamin/pctsmx).
"""
# Build the working tree list (live, HT>0, raw PROB, per-tree crown area) and the canopy
# stratification (sstage.f:166-465): up to 3 height-gap strata, each with its cover range
# (incl. gap trees) and OK flag (> CCMIN). Shared by `structure_class` and `structure_report`.
function _ss_strata(s::StandState)
    t = s.trees; p = s.plot; co = s.coef
    cccoef = Float64(s.control.cc_coef)
    gappct = Float64(s.control.strclass_thresh[1]); ccmin = Float64(s.control.strclass_thresh[4])
    n = 0; ht = Float64[]; dbh = Float64[]; tpa = Float64[]; crarea = Float64[]
    species = Int[]; icr = Float64[]
    @inbounds for i in 1:t.n
        t.height[i] > 0f0 && t.tpa[i] > 0f0 || continue
        sp2 = strip(co.code_alpha[Int(t.species[i])])
        cw = crown_width(co, sp2, t.dbh[i], t.height[i], Float32(t.crown_pct[i]), 0,
                         p.latitude, p.longitude, p.elevation)
        pa = Float64(t.tpa[i])                          # PROB (raw, as SSTAGE uses it — NOT /GROSPC)
        n += 1; push!(ht, Float64(t.height[i])); push!(dbh, Float64(t.dbh[i]))
        push!(tpa, pa); push!(crarea, Float64(cw)^2 * pa * 0.785398)
        push!(species, Int(t.species[i])); push!(icr, Float64(t.crown_pct[i]))
    end
    data = (; n, ht, dbh, tpa, crarea, species, icr, cccoef)
    n == 0 && return (data..., ord = Int[], strata = NTuple{4,Int}[], oks = Bool[],
                      covers = Float64[], nstr = 0, tprob = 0.0, cover = 0.0)
    tprob = sum(tpa)
    ord = sortperm(ht; rev = true)                  # INDEX: trees by height descending
    # height-gap stratification: track the two largest gaps (sstage.f:300-388)
    diff1 = 0.0; diff2 = 0.0; id1i1 = id1i2 = id2i1 = id2i2 = 0
    iilg = 1; ilarge = ord[1]; sumprb = 0.0
    @inbounds for ii in 2:n
        ismall = ord[ii]
        x = ht[ilarge] * gappct * 0.01; x < 10.0 && (x = 10.0)
        if ht[ismall] < ht[ilarge] - x
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
        else
            if tpa[ismall] + sumprb < 2.0
                sumprb += tpa[ismall]
            else
                ilarge = ismall; iilg = ii; sumprb = 0.0
            end
        end
    end
    if id1i1 > id2i1 && id2i1 > 0
        id1i1, id2i1 = id2i1, id1i1; id1i2, id2i2 = id2i2, id1i2
    end
    is1i1 = 1; is1i2 = n; is2i1 = 0; is2i2 = 0; is3i1 = 0; is3i2 = 0
    if id1i1 > 0; is1i2 = id1i1; is2i1 = id1i2; is2i2 = n; end
    if id2i1 > 0; is2i2 = id2i1; is3i1 = id2i2; is3i2 = n; end
    # each potential stratum = (lo, hi, cover_lo, cover_hi); cover incl. gap trees (sstage.f:430-465)
    str = NTuple{4,Int}[]
    push!(str, (is1i1, is1i2, is1i1, max(is1i2, is2i1 - 1)))
    is2i1 > 0 && push!(str, (is2i1, is2i2, is2i1, max(is2i2, is3i1 - 1)))
    is3i1 > 0 && push!(str, (is3i1, is3i2, is3i1, is3i2))
    covers = [_ss_cover(crarea, (ord[k] for k in c1:c2), cccoef) for (_, _, c1, c2) in str]
    oks = covers .> ccmin
    nstr = count(oks)
    cover = _ss_cover(crarea, 1:n, cccoef)
    if nstr == 0 && tprob >= Float64(s.control.strclass_thresh[5])   # < TPAMIN ⇒ stays 0
        str = [(1, n, 1, n)]; covers = [cover]; oks = [true]; nstr = 1
    end
    return (data..., ord, strata = str, oks, covers, nstr, tprob, cover)
end

function structure_class(s::StandState; iba::Int = 1)
    th = s.control.strclass_thresh
    ssdbh = Float64(th[2]); sawdbh = Float64(th[3]); pctsmx = Float64(th[6])
    st = _ss_strata(s)
    st.n == 0 && return (class = 0, nstr = 0, cover = 0.0, strdbh = 0.0)
    st.nstr == 0 && return (class = 0, nstr = 0, cover = st.cover, strdbh = 0.0)
    # dominant stratum = the first OK one; its SSTGHP 70th-percentile DBH (sstage.f:487-576)
    di = findfirst(st.oks)
    dlo, dhi = st.strata[di][1], st.strata[di][2]
    tmpdbh, _ = _ss_dbhnom(st.ord, dlo, dhi, st.ht, st.dbh, st.tpa, st.crarea)
    dmind = st.dbh[st.ord[dhi]]
    cls = 0
    if st.nstr == 1
        if tmpdbh < ssdbh
            cls = 1
        elseif tmpdbh < sawdbh
            cls = 2
            xbamax = Float64(s.control.ba_max)
            (xbamax > 0 && _event_bsdi(s) < 0.01 * pctsmx * xbamax) && (cls = 1)   # PCTSMX demotion
        else
            cls = dmind < 3.0 ? 6 : 5
        end
    elseif st.nstr == 2
        cls = tmpdbh < ssdbh ? 1 : tmpdbh < sawdbh ? 3 : 6
    else
        cls = tmpdbh < ssdbh ? 1 : tmpdbh < sawdbh ? 4 : 6
    end
    return (class = cls, nstr = st.nstr, cover = st.cover, strdbh = tmpdbh)
end

"""
    structure_report(s) -> (; class, nstr, cover, strata)

Per-stratum SSTAGE "Structural statistics" data (sstage.f / SSTGHP) for the `.out` report: the
class, number of valid strata, whole-stand canopy cover, and one record per OK stratum (uppermost
first) — `(; dbh, nomht, lght, smht, crnbase, cover, sp1, sp2)`: the 70th-percentile DBHNOM + its
window mean height, the stratum's tallest/shortest height, the mean height-to-crown-base (the
report's "Bas" column = ICRB, not stand basal area), cover, and the two species codes with the most
crown area. All from the same machinery the class uses (validated bit-exact vs the Fortran report).
"""
function structure_report(s::StandState)
    st = _ss_strata(s)
    cls = structure_class(s).class
    strata = NamedTuple[]
    for (k, ok) in enumerate(st.oks)
        ok || continue
        lo, hi, _, _ = st.strata[k]
        dbhnom, nomht = _ss_dbhnom(st.ord, lo, hi, st.ht, st.dbh, st.tpa, st.crarea)
        lght = st.ht[st.ord[lo]]; smht = st.ht[st.ord[hi]]   # tallest / shortest (height-sorted)
        acb = 0.0; sp_ = 0.0; spc = Dict{Int,Float64}()      # ACB = Σ crown-base-ht·PROB (sstage.f:790)
        @inbounds for j in lo:hi
            i = st.ord[j]
            acb += st.ht[i] * (1.0 - st.icr[i] * 0.01) * st.tpa[i]
            sp_ += st.tpa[i]
            spc[st.species[i]] = get(spc, st.species[i], 0.0) + st.crarea[i]   # crown area by species
        end
        crnbase = sp_ > 1e-4 ? acb / sp_ : 0.0               # ICRB = mean height to crown base
        sp = sort(collect(spc); by = x -> -x[2])             # top crown-area species
        sp1 = isempty(sp) ? 0 : sp[1][1]; sp2 = length(sp) >= 2 ? sp[2][1] : 0
        push!(strata, (; dbh = dbhnom, nomht, lght, smht, crnbase, cover = st.covers[k], sp1, sp2))
    end
    return (class = cls, nstr = st.nstr, cover = st.cover, strata = strata)
end
