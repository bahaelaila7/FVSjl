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
    tmpdbh, _ = _ss_dbhnom(st.ord, dlo, dhi, st.ht, st.dbh, st.tpa, st.crarea)
    dmind = st.dbh[st.ord[dhi]]
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
    for k in eachindex(st.strata)
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
        status = k == dom ? 2 : st.oks[k] ? 1 : 0            # D: 2=dominant, 1=OK, 0=not (IS_OK)
        push!(strata, (; dbh = dbhnom, nomht, lght, smht, crnbase, cover = st.covers[k], sp1, sp2, status))
    end
    return (class = cls, nstr = st.nstr, cover = st.cover, strata = strata)
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
