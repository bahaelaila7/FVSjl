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

using Printf

# Default SSTAGE thresholds (isstag.f:32-37), in the `control.strclass_thresh` order:
#   gappct=30 (stratum-gap % height drop), ssdbh=5 (small-tree DBH), sawdbh=25 (sawtimber DBH),
#   ccmin=5 (min stratum cover %), tpamin=200 (min TPA), pctsmx=30 (% MaxSDI for SE).
const SS_THRESH_DEFAULT = (30.0f0, 5.0f0, 25.0f0, 5.0f0, 200.0f0, 30.0f0)

# SSTGHP (sstage.f:740-873), statement-for-statement in REAL (Float32) with FVS's IFIX(x+.5) rounding.
# For stratum INDEX range ord[i1:i2] (height-descending): IHTL/IHTS = rounded tallest/shortest height; walk the
# stratum summing WK6 cover (SUM), PROB (SP) and crown-base·PROB (ACB), cover by species, and the 95%-cover
# cutoff I3 (41382 = .95·43560); ICRB = IFIX(ACB/SP+.5); top-2 cover species by a strict-> scan over species
# index; then WK4 = single-tree crown area, RDPSRT(.FALSE.) descending on ord[i1:i3], PCTILE on WK6, the record
# nearest the 70th percentile (first minimum), and the PROB-weighted mean DBH / rounded mean height of the ±4
# window. Returns (dbhnom::Float32, iht, ihtl, ihts, icrb, msp1, msp2). `ord` is permuted in place on i1:i3 as
# FVS permutes INDEX (strata ranges don't overlap).
function _ss_sstghp!(ord::Vector{Int}, i1::Int, i2::Int, ht, dbh, tpa, icr, wk6::Vector{Float32}, species)
    (i1 == 0 || i2 == 0) && return (0f0, 0, 0, 0, 0, 0, 0)
    ifix5(x::Float32) = trunc(Int, x + 0.5f0)
    ihtl = ifix5(Float32(ht[ord[i1]])); ihts = ifix5(Float32(ht[ord[i2]]))
    sum_ = 0f0; sp = 0f0; acb = 0f0; i3 = -1
    spcov = Dict{Int,Float32}()
    @inbounds for ii in i1:i2
        i = ord[ii]
        sum_ += wk6[i]; sp += Float32(tpa[i])
        acb += (Float32(ht[i]) * (1f0 - Float32(icr[i]) * 0.01f0)) * Float32(tpa[i])
        spcov[species[i]] = get(spcov, species[i], 0f0) + wk6[i]
        (sum_ > 41382f0 && i3 == -1) && (i3 = ii)
    end
    i3 == -1 && (i3 = i2)
    x1 = 0f0; x2 = 0f0; msp1 = 0; msp2 = 0                  # sstage.f:803-816 (strict >, species-index order)
    for is in sort!(collect(keys(spcov)))
        w = spcov[is]
        if w > x1
            x2 = x1; x1 = w; msp2 = msp1; msp1 = is
        elseif w > x2
            x2 = w; msp2 = is
        end
    end
    icrb = sp > 0.0001f0 ? ifix5(acb / sp) : 0
    m = i3 - i1 + 1
    wk4 = Dict{Int,Float32}()
    @inbounds for ii in i1:i3
        i = ord[ii]; wk4[i] = wk6[i] / Float32(tpa[i])
    end
    key = zeros(Float32, length(ht)); for (i, v) in wk4; key[i] = v; end
    rdpsrt!(m, key, view(ord, i1:i3), false)                 # RDPSRT(I3-I1+1,WK4,INDEX(I1),.FALSE.)
    # PCTILE(N,INDEX(I1),WK6,WK4,T): cumulative WK6 from the bottom, / (TOT/100); top record = 100.
    pct = Dict{Int,Float32}()
    if m > 1
        pct[ord[i3]] = wk6[ord[i3]]
        for j in (i3 - 1):-1:i1
            pct[ord[j]] = pct[ord[j + 1]] + wk6[ord[j]]
        end
        tot = pct[ord[i1]]
        pctin1 = tot / 100f0
        if tot > 0f0
            for j in (i1 + 1):i3; pct[ord[j]] = pct[ord[j]] / pctin1; end
            pct[ord[i1]] = 100f0
        else
            pct[ord[i1]] = pctin1
        end
    end
    i70 = i1                                                 # N==1 ⇒ the lone record (PCTILE returns early)
    if m > 1
        diff = 1f30
        for ii in i1:i3
            d = abs(pct[ord[ii]] - 70f0)
            if d < diff; i70 = ii; diff = d; end
        end
    end
    k1 = max(i70 - 4, i1); k2 = min(i70 + 4, i3)
    sd = 0f0; sh = 0f0; sp = 0f0
    @inbounds for ii in k1:k2
        i = ord[ii]
        sd += Float32(dbh[i]) * Float32(tpa[i]); sh += Float32(ht[i]) * Float32(tpa[i]); sp += Float32(tpa[i])
    end
    dbhnom = 0f0; iht = 0
    if sp > 0.0001f0
        dbhnom = sd / sp; iht = ifix5(sh / sp)
    end
    return (dbhnom, iht, ihtl, ihts, icrb, msp1, msp2)
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
function _ss_strata(s::StandState; thresh = s.control.strclass_thresh)
    t = s.trees; p = s.plot; co = s.coef
    cccoef = Float64(s.control.cc_coef)
    gappct = Float64(thresh[1]); ccmin = Float64(thresh[4])
    n = 0; ht = Float64[]; dbh = Float64[]; tpa = Float64[]; crarea = Float64[]
    species = Int[]; icr = Float64[]; crarea32 = Float32[]
    # CR uses cr_cwcalc (cwcalc.f IWHO=0, the CRWDTH FMSSTAGE reads), not the generic crown_width (0.5 default
    # for CR ⇒ zero cover ⇒ wrong strata/class). Eastern variants keep crown_width. Precompute CR stand inputs.
    _cr_ss = s.variant isa CentralRockies
    _so_ss = s.variant isa SouthCentralOregon    # SO CRWDTH via so_cwcalc (SOMAP Crookston R6, forest-601 BF)
    _oc_ss = s.variant isa OregonCoast           # OC/ORGANON CRWDTH via oc_cwcalc (else generic crown_width→0 cover)
    # Every other western variant: its forest-grown cwcalc CRWDTH (_forest_crwdth — the TreeList CrWidth value);
    # the eastern crown_width returns the 0.5 default for western species ⇒ ~zero cover ⇒ wrong stage/class.
    _west_ss = !(_cr_ss || _so_ss || _oc_ss) && _has_forest_crwdth(s.variant)
    _cr_ba = (_cr_ss || _so_ss || _oc_ss) ? p.basal_area : 0f0
    _cr_el = (_cr_ss || _so_ss || _oc_ss) ? p.elevation : 0f0
    _cr_hi = (_cr_ss || _so_ss || _oc_ss) ? _cr_hopkins(p.latitude, p.longitude, p.elevation) : 0f0
    @inbounds for i in 1:t.n
        t.height[i] > 0f0 && t.tpa[i] > 0f0 || continue
        cw = _cr_ss ? cr_cwcalc(Int(t.species[i]), t.dbh[i], t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi) :
             _so_ss ? so_cwcalc(Int(t.species[i]), t.dbh[i], t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi) :
             _oc_ss ? oc_cwcalc(Int(t.species[i]), t.dbh[i], t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi) :
             _west_ss ? _forest_crwdth(s, Int(t.species[i]), t.dbh[i], t.height[i], t.crown_pct[i]) :  # WK6=CRWDTH (sstage.f:238)
             crown_width(co, strip(co.code_alpha[Int(t.species[i])]), t.dbh[i], t.height[i],
                         Float32(t.crown_pct[i]), 0, p.latitude, p.longitude, p.elevation)
        pa = Float64(t.tpa[i])                          # PROB (raw, as SSTAGE uses it — NOT /GROSPC)
        n += 1; push!(ht, Float64(t.height[i])); push!(dbh, Float64(t.dbh[i]))
        push!(tpa, pa); push!(crarea, Float64(cw)^2 * pa * 0.785398)
        push!(crarea32, ((Float32(cw) * Float32(cw)) * Float32(pa)) * 0.785398f0)   # WK6 in REAL (sstage.f:275-276)
        push!(species, Int(t.species[i])); push!(icr, Float64(t.crown_pct[i]))
    end
    data = (; n, ht, dbh, tpa, crarea, crarea32, species, icr, cccoef)
    n == 0 && return (data..., ord = Int[], strata = NTuple{4,Int}[], oks = Bool[],
                      covers = Float64[], nstr = 0, tprob = 0.0, cover = 0.0)
    tprob = sum(tpa)
    # INDEX: trees by height descending — RDPSRT(NTREES,HT,INDEX,.FALSE.) (sstage.f:269), the Quickersort tie
    # order (not a stable sort): equal-height records' order feeds the SSTGHP 95%-cover cutoff / PCTILE window.
    ord = collect(1:n); rdpsrt!(n, ht, ord, false)
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
    if nstr == 0 && tprob >= Float64(thresh[5])   # < TPAMIN ⇒ stays 0
        str = [(1, n, 1, n)]; covers = [cover]; oks = [true]; nstr = 1
    end
    return (data..., ord, strata = str, oks, covers, nstr, tprob, cover)
end

function structure_class(s::StandState; iba::Int = 1, thresh = s.control.strclass_thresh)
    th = thresh
    ssdbh = Float64(th[2]); sawdbh = Float64(th[3]); pctsmx = Float64(th[6])
    st = _ss_strata(s; thresh = thresh)
    st.n == 0 && return (class = 0, nstr = 0, cover = 0.0, strdbh = 0.0)
    # Single-canopy-tree path (sstage.f:234-268): a stand with only ONE canopy tree is classified by its
    # CROWN-AREA cover (WK6 = CW²·TPA·π/4, already = `st.crarea`), NOT the stratum DBHNOM. Checked BEFORE the
    # stratification (so it fires even when the lone tree forms no OK stratum). Cover < CCMIN% of an acre
    # (435.6 = 0.01·43560 sq ft) ⇒ class 0 (or SI=1 if TPA≥TPAMIN); else SSD→1, SAW→2 (SE→SI demote when the
    # stand SDI < 0.01·PCTSMX·MaxSDI), else 5. (FVS uses SDIAC here; for a lone-tree no-cut stand SDIAC=SDIBC,
    # so the same _event_bsdi as the nstr==1 path.)
    if st.n == 1
        wk6 = st.crarea[1]; tpa1 = st.tpa[1]; dbh1 = st.dbh[1]
        ccmin = Float64(th[4]); tpamin = Float64(th[5])
        cls = if wk6 < 435.6 * ccmin
            tpa1 >= tpamin ? 1 : 0
        elseif dbh1 < ssdbh
            1
        elseif dbh1 < sawdbh
            (_event_bsdi(s) < 0.01 * pctsmx * Float64(stand_sdimax(s))) ? 1 : 2
        else
            5
        end
        # NSTR stays 0 here: the single-canopy branch GOTO-80s before the stratification sets NSTR≥1
        # (sstage.f:235 jumps past :388-478), so the report's N column is 0 for a lone canopy tree.
        return (class = cls, nstr = 0, cover = st.cover, strdbh = dbh1)
    end
    st.nstr == 0 && return (class = 0, nstr = 0, cover = st.cover, strdbh = 0.0)
    # dominant stratum = the first OK one; its SSTGHP 70th-percentile DBH (sstage.f:487-576)
    di = findfirst(st.oks)
    dlo, dhi = st.strata[di][1], st.strata[di][2]
    # SSTGHP is called for every stratum in order (sstage.f:494-513), permuting INDEX within each; DMIND =
    # DBH(INDEX(ISkI2)) is read AFTER those permutations (sstage.f:519-535).
    ordw = copy(st.ord); dbhs = Float32[]
    for (lo, hi, _, _) in st.strata
        push!(dbhs, _ss_sstghp!(ordw, lo, hi, st.ht, st.dbh, st.tpa, st.icr, st.crarea32, st.species)[1])
    end
    tmpdbh = Float64(dbhs[di])
    dmind = st.dbh[ordw[dhi]]
    cls = 0
    if st.nstr == 1
        if tmpdbh < ssdbh
            cls = 1
        elseif tmpdbh < sawdbh
            cls = 2
            # PCTSMX SE→SI demotion (sstage.f:154,544-550): XBAMAX = BTSDIX = the per-cycle SDICAL stand
            # MaxSDI (grincr.f:240), NOT the user BAMAX keyword. SDIBC < 0.01·PCTSMX·MaxSDI → demote to SI.
            xbamax = Float64(stand_sdimax(s))
            (_event_bsdi(s) < 0.01 * pctsmx * xbamax) && (cls = 1)
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

"SSTAGE's NTREES (sstage.f:220-226): the live records carrying more than 0.00001 TPA."
_sstage_ntrees(s::StandState) = (t = s.trees; count(i -> t.tpa[i] > 0.00001f0, 1:t.n))

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
    dom = st.nstr > 0 ? findfirst(st.oks) : 0     # the dominant stratum = first OK
    strata = NamedTuple[]
    # sstage.f:475-476: >1 tree, no valid stratum AND TPROB < TPAMIN ⇒ GOTO 80, skipping every SSTGHP call, so
    # all stratum fields keep their sstage.f:168-200 zero init (DB: 0/"--"/status 0). (_ss_strata already turned
    # the TPROB ≥ TPAMIN case into one all-tree stratum.) NOT the single-record path (sstage.f:238-266,
    # NTREES≤1): that one fills stratum 1 from the lone record before its own GOTO 80 — keep emitting it.
    (st.nstr == 0 && st.n > 1) && return (class = cls, nstr = 0, cover = st.cover, strata = strata)
    ordw = copy(st.ord)
    for k in eachindex(st.strata)
        lo, hi, _, _ = st.strata[k]
        if st.n == 1
            # single-record path (sstage.f:251-265): DBHS1=DBH, IHTS1=IHTSS1=IHTLS1=IFIX(HT+.5),
            # ICRBS1=IFIX(HT·(1-ICR·.01)+.5), SP11 = the record's species, no second species.
            i = st.ord[1]; h32 = Float32(st.ht[i])
            iht = trunc(Int, h32 + 0.5f0)
            icrb = trunc(Int, h32 * (1f0 - Float32(st.icr[i]) * 0.01f0) + 0.5f0)
            push!(strata, (; dbh = Float64(Float32(st.dbh[i])), nomht = iht, lght = iht, smht = iht, crnbase = icrb,
                           cover = st.covers[k], sp1 = st.species[i], sp2 = 0, status = 0))
            continue
        end
        dbhnom, iht, ihtl, ihts, icrb, sp1, sp2 = _ss_sstghp!(ordw, lo, hi, st.ht, st.dbh, st.tpa, st.icr,
                                                             st.crarea32, st.species)
        status = k == dom ? 2 : st.oks[k] ? 1 : 0            # D: 2=dominant, 1=OK, 0=not (IS_OK)
        push!(strata, (; dbh = Float64(dbhnom), nomht = iht, lght = ihtl, smht = ihts, crnbase = icrb,
                       cover = st.covers[k], sp1, sp2, status))
    end
    # single-record path GOTO-80s before NSTR is counted (sstage.f:266) ⇒ Number_of_Strata 0 (as structure_class)
    return (class = cls, nstr = st.n == 1 ? 0 : st.nstr, cover = st.cover, strata = strata)
end

const _SS_CLASS_LABEL = ("0=BG", "1=SI", "2=SE", "3=UR", "4=YM", "5=OS", "6=OM")  # SSCODES (sstage.f:73)

"""
    structure_report_row(s, year, cd) -> String

One "Structural statistics" report row (sstage.f FORMAT 90), byte-for-byte: Year Cd, then for the 3
strata DBH/Nom-Ht/Lg-Ht/Sm-Ht/Bas(crown base)/Cov/Sp1/Sp2/D (zeros + "--" for absent strata), then
N-Strata, Tot-Cov, and the class label. `cd` is the removal code (0 = before-thin, 1 = after).
"""
function structure_report_row(s::StandState, year::Integer, cd::Integer)
    r = structure_report(s); co = s.coef
    spcode(i) = i > 0 ? rpad(strip(co.code_alpha[i]), 3) : "-- "
    blk(st) = @sprintf(" %5.1f %3d %3d %3d %3d %3d", st.dbh, round(Int, st.nomht), round(Int, st.lght),
                       round(Int, st.smht), round(Int, st.crnbase), round(Int, st.cover)) *
              " " * spcode(st.sp1) * " " * spcode(st.sp2) * @sprintf(" %1d", st.status)
    line = @sprintf("%4d %2d", year, cd)
    for k in 1:3
        line *= k <= length(r.strata) ? blk(r.strata[k]) :
                @sprintf(" %5.1f %3d %3d %3d %3d %3d", 0.0, 0, 0, 0, 0, 0) * " -- " * " -- " * " 0"
    end
    return line * @sprintf(" %1d %3d  %s", r.nstr, round(Int, r.cover), _SS_CLASS_LABEL[r.class + 1])
end

# The fixed Structural-statistics column-header lines (sstage.f FORMAT 85, sans the $#*% page marks).
const _SS_REPORT_HEADER = (
    "        ------------ Stratum 1 ------------ ------------ Stratum 2 ------------ ------------ Stratum 3 ------------",
    "     Rm       ---Height-- -Crown- -Major- C       ---Height-- -Crown- -Major- C       ---Height-- -Crown- -Major- C N Tot Struc",
    "Year Cd  DBH  Nom  Lg  Sm Bas Cov Sp1 Sp2 D  DBH  Nom  Lg  Sm Bas Cov Sp1 Sp2 D  DBH  Nom  Lg  Sm Bas Cov Sp1 Sp2 D S Cov Class",
    "---- -- ----- --- --- --- --- --- --- --- - ----- --- --- --- --- --- --- --- - ----- --- --- --- --- --- --- --- - - --- -----")

"""
    write_structure_report(io, stand, ncyc; period=5, stand_id="", mgmt_id="NONE")

Write the SSTAGE "Structural statistics" report (sstage.f) for `stand` over `ncyc+1` cycles to `io`,
byte-for-byte vs the Fortran `.out` block (the page-control marks aside): the header, then per cycle
a before-thin (Rm=0) and after-thin (Rm=1) row. Steps the stand's projection (grow_cycle!), so pass
a stand already through `setup_growth!`/`compute_volumes!` that you don't need afterwards.
"""
function write_structure_report(io::IO, stand::StandState, ncyc::Integer;
                                period::Integer = 5, stand_id::AbstractString = "",
                                mgmt_id::AbstractString = "NONE")
    println(io, "Structural statistics for stand: ", rpad(strip(stand_id), 26), "  MgmtID: ", strip(mgmt_id))
    println(io)
    for h in _SS_REPORT_HEADER; println(io, h); end
    for c in 0:ncyc
        compute_density!(stand)
        yr = Int(current_cycle_year(stand))
        println(io, structure_report_row(stand, yr, 0))     # before-thin (Rm=0): the pre-thin stand
        # Apply this cycle's scheduled thin so the after-thin row reflects the POST-thin stand — its cover,
        # strata, and the after-thin MaxSDI (ATSDIX), per sstage.f:145-155 (IBA≠1 + ONTREM>0). When nothing is
        # cut the stand is unchanged, so the after-thin row equals the before-thin row (sstage.f:146). cuts! is
        # idempotent per year (cuts.jl years_cut guard), so grow_cycle!'s own cut below becomes a no-op.
        cuts!(stand; fint = Float32(period))
        compute_density!(stand)
        println(io, structure_report_row(stand, yr, 1))     # after-thin (Rm=1): the post-thin stand
        c < ncyc && grow_cycle!(stand; fint = Float32(period))
    end
    return io
end
