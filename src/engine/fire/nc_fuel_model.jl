# =============================================================================
# fire/nc_fuel_model.jl — NC (Klamath) California-CWHR surface fuel-model selection
#
# Ported from: nc/fmcfmd.f (FMCFMD) + nc/cwhr.f (CWHR + PCNETAVG) — the WS/CA/NC "California
# Wildlife Habitat Relationships" classifier. UNLIKE the interior western variants (IE/CI/EM/BM,
# cover-type × PERCOV) and the SN forest-type path, NC classifies the stand into a CWHR
# SIZE (1..6) × DENSITY (S/P/M/D) structural stage, maps (forest-type IFT × stage JSS) → a base
# Anderson fuel model via a 9×18 table, blends the density sub-models (FMDYN over CWXPTS), and
# adds the natural-fuel candidates 10/12/13 (+ the post-activity 11), resolving the final weighted
# model set with `_fmdyn` over NC's 13-model XPTS. Values transcribed verbatim from nc/*.f.
# =============================================================================

# XPTS (nc/fmcfmd.f:147-160): models 1-9 share (5,15); 10 & 11 share (15,30); 12=(30,60); 13=(45,100).
const _NC_FMD_XPTS = Float32[
    5.0  15.0;  5.0  15.0;  5.0  15.0;  5.0  15.0;  5.0  15.0;  5.0  15.0;  5.0  15.0;
    5.0  15.0;  5.0  15.0; 15.0  30.0; 15.0  30.0; 30.0  60.0; 45.0 100.0]

# CWHR density sub-model lines (nc/fmcfmd.f:124-125) — (S,P,M,D) = tags 10/20/30/40.
const _NC_CWXPTS = Float32[10.3 195.0; 45.0 225.0; 88.3 265.0; 185.0 185.0]
const _NC_CCBP   = (10f0, 40f0, 70f0)        # canopy-cover breakpoints
const _NC_DBHBP  = (1f0, 5f0, 11f0, 21f0)    # DBH class breakpoints
const _NC_CWHR_IPTR = [10, 20, 30, 40]

# FMD_R5 / FMD_R6 CWHR fuel-model classification (nc/fmcfmd.f:89-114): 9 forest-type rows ×
# 18 structural-stage columns. NC's two tables are identical (both blocks match); kept separate
# for faithfulness. Row = IFT (1 pine … 9 other-soft); col = JSS (1, 2S..2D, 3S..3D, …, 6).
const _NC_FMD_R5 = Int[
    5 6 6 6 6 2 2 9 9 2 2 2 9 2 2 9 9 10;   # 1 pine
    5 5 5 8 8 11 11 8 8 8 8 8 8 8 8 8 8 10;  # 2 red fir
    5 5 5 8 8 11 11 11 8 8 8 8 8 8 8 8 8 10; # 3 white fir E
    5 5 5 8 8 11 11 8 8 8 8 8 8 8 8 8 8 10;  # 4 white fir W
    5 5 5 6 6 6 6 8 8 11 11 9 8 11 11 9 8 10; # 5 Douglas-fir
    5 5 5 6 6 11 11 11 9 9 9 9 9 9 9 9 9 10; # 6 hardwoods
    5 5 5 6 6 6 6 6 9 9 9 8 8 8 8 8 8 10;    # 7 mixed pine
    5 5 5 6 6 6 6 6 8 6 6 8 8 6 6 8 8 10;    # 8 mixed conifer
    5 5 5 6 6 6 6 6 8 6 6 8 8 6 6 8 8 10]    # 9 other softwood
const _NC_FMD_R6 = _NC_FMD_R5

# PCNETAVG (nc/cwhr.f:500-511): random-packing-corrected canopy cover.
@inline function _nc_pcnetavg(gross::Float32)::Float32
    pcrnd = (1f0 - 1f0 / exp(gross / 100f0)) * 100f0
    avg1 = pcrnd / 100f0 * pcrnd + (100f0 - pcrnd) / 100f0 * gross
    avg1 > 100f0 && (avg1 = 100f0)
    return (pcrnd + avg1) / 2f0
end

# FINDMOD (nc/fmcfmd.f:384): model NUMBER → index in IPTR=1:13 (identity for 1..13, else 8).
@inline _nc_findmod(m::Integer)::Int = (1 <= m <= 13) ? Int(m) : 8

# FINDJSS (nc/fmcfmd.f:405): CWHR (size, density) → column 1..18 of the classification table.
@inline function _nc_findjss(sz::Char, dn::Char)::Int
    sz == 'X' && return 1
    sz == '1' && return 1
    sz == '6' && return 18
    base = sz == '2' ? 1 : sz == '3' ? 5 : sz == '4' ? 9 : sz == '5' ? 13 : 0
    base == 0 && return 0
    d = dn == 'S' ? 1 : dn == 'P' ? 2 : dn == 'M' ? 3 : dn == 'D' ? 4 : 0
    d == 0 && return 0
    return base + d
end

# CWHR (nc/cwhr.f): classify the stand into a size-class char SZ, density char DN, and the dynamic
# density sub-model tags/weights CWHR_MOD/CWHR_WT (10/20/30/40). `cws` = per-tree crown widths (CRWDTH).
# The base cwhr.f is IDENTICAL across the California variants (ws/cwhr.f == nc/cwhr.f, verified); the only
# per-variant differences are the DBHBP/CCBP/CWXPTS/IPTR constants (passed as DATA args in fmcfmd.f), so the
# core is parameterized here and WS calls it with its own constants (defaults = NC's for bit-identical NC).
nc_cwhr(s::StandState, cws::Vector{Float32}) =
    _ca_cwhr(s, cws, _NC_DBHBP, _NC_CCBP, _NC_CWXPTS, _NC_CWHR_IPTR)

function _ca_cwhr(s::StandState, cws::Vector{Float32}, _NC_DBHBP, _NC_CCBP, _NC_CWXPTS, _NC_CWHR_IPTR)
    t = s.trees
    ht  = zeros(Float32, 6)   # [1]=HT(0) total, [2..6]=HT(1..5) by DBH class
    cc  = zeros(Float32, 6)
    tpa = zeros(Float32, 6)
    basal = 0f0
    @inbounds for i in 1:t.n
        p = t.tpa[i]; p > 0f0 || continue
        d = t.dbh[i]
        basal += p * d * d * 0.0054542f0
        carea = 3.1415927f0 * cws[i] * cws[i] / 4f0
        k = d < _NC_DBHBP[1] ? 1 : d < _NC_DBHBP[2] ? 2 : d < _NC_DBHBP[3] ? 3 : d < _NC_DBHBP[4] ? 4 : 5
        ht[k+1]  += t.height[i] * p
        cc[k+1]  += carea * p
        tpa[k+1] += p
    end
    @inbounds for i in 2:6
        ht[1] += ht[i]; cc[1] += cc[i]; tpa[1] += tpa[i]
    end
    @inbounds for i in 1:6
        cc[i] /= 435.60f0
    end
    ccnetavg = _nc_pcnetavg(cc[1])
    @inbounds for k in 4:6          # HT(3..5) → average heights
        tpa[k] > 0f0 && (ht[k] /= tpa[k])
    end
    ccadj = _nc_pcnetavg(cc[1])

    cwhr_mod = zeros(Int, 4); cwhr_wt = zeros(Float32, 4)

    # (1) non-forest / initiating (cc < 10)
    if ccnetavg < 10f0
        cwhr_mod[1] = 10; cwhr_wt[1] = 1f0
        return (tpa[1] >= 150f0 ? '1' : 'X', '-', cwhr_mod, cwhr_wt)
    end

    # (2) multi-story (dense, all of HT3/4/5 present) → size "6"
    if ccnetavg >= _NC_CCBP[3] && ht[6] > 0f0 && ht[5] > 0f0 && ht[4] > 0f0
        if cc[6] >= 20f0 && cc[6] <= 80f0
            two3 = 2f0 / 3f0
            if (cc[4] >= 20f0 && ht[4] / ht[6] <= two3) ||
               (cc[5] >= 20f0 && ht[5] / ht[6] <= two3) ||
               ((cc[4] + cc[5]) >= 30f0 &&
                ((ht[4] * tpa[4] + ht[5] * tpa[5]) / (tpa[4] + tpa[5])) / ht[6] <= two3)
                cwhr_mod[1] = 40; cwhr_wt[1] = 1f0
                return ('6', '-', cwhr_mod, cwhr_wt)
            end
        end
    end

    local sz::Char
    if ccadj < _NC_CCBP[1]
        # (3) sparse: classify by predominance of cover
        if (cc[4] + cc[5] + cc[6]) < (cc[2] + cc[3])
            sz = cc[3] >= cc[2] ? '2' : '1'
        else
            sz = cc[4] > (cc[5] + cc[6]) ? '3' : (cc[6] >= (cc[5] + cc[4]) ? '5' : '4')
        end
        cwhr_mod[1] = 10; cwhr_wt[1] = 1f0
        return (sz, 'S', cwhr_mod, cwhr_wt)
    else
        # (4) QMD of the largest 75% (or 100% for small-tree stands) of stand BA
        qmdpct = (cc[5] + cc[6] < 10f0) ? 1.00f0 : 0.75f0
        targetba = basal * qmdpct
        sumba = 0f0; ltrees = 0f0; sumd2 = 0f0
        ixs = collect(1:t.n)
        t.n > 0 && rdpsrt!(t.n, t.dbh, ixs, true)     # DBH descending
        @inbounds for j in 1:t.n
            i = ixs[j]; p = t.tpa[i]; d = t.dbh[i]
            treeba = p * d * d * 0.0054542f0
            if sumba + treeba < targetba
                sumba += treeba; sumd2 += p * d * d; ltrees += p
            else
                sng1 = treeba > 0f0 ? (targetba - sumba) / treeba : 0f0
                ltrees += p * sng1; sumd2 += d * d * p * sng1
                break
            end
        end
        qmd = ltrees > 0f0 ? sqrt(sumd2 / ltrees) : 0f0
        if qmd < _NC_DBHBP[1]
            cwhr_mod[1] = 10; cwhr_wt[1] = 1f0
            return ('1', '-', cwhr_mod, cwhr_wt)
        elseif qmd < _NC_DBHBP[2]
            sz = '2'; ccadj = _nc_pcnetavg(cc[1])
        elseif qmd < _NC_DBHBP[3]
            sz = '3'; cc[1] = cc[3] + cc[4] + cc[5] + cc[6]; ccadj = _nc_pcnetavg(cc[1])
        elseif qmd < _NC_DBHBP[4]
            sz = '4'; cc[1] = cc[4] + cc[5] + cc[6]; ccadj = cc[1]
        else
            sz = '5'; cc[1] = cc[4] + cc[5] + cc[6]; ccadj = cc[1]
        end
    end

    # (5) density char (nc/cwhr.f:421-432)
    dn = if ccadj < _NC_CCBP[1]
        'S'
    elseif ccadj < _NC_CCBP[2]
        ccnetavg < _NC_CCBP[1] ? 'S' : 'P'
    elseif ccadj < _NC_CCBP[3]
        ccnetavg < _NC_CCBP[2] ? 'P' : 'M'
    else
        ccnetavg < _NC_CCBP[3] ? 'M' : 'D'
    end

    # (6) density sub-model weights via FMDYN over CWXPTS (the S/P/M/D lines), point=(ccadj, ccnetavg)
    dyn = _fmdyn(ccadj, ccnetavg, ones(Float32, 4), _NC_CWXPTS; iptr = _NC_CWHR_IPTR)
    for (idx, (m, w)) in enumerate(dyn)
        idx > 4 && break
        cwhr_mod[idx] = m; cwhr_wt[idx] = w
    end
    # (7) >2 candidates: drop the low-slope D (tag 40), keep the top-2 by tag, renormalize (cwhr.f:455-484)
    nz = count(>(0), cwhr_mod)
    if nz > 2
        xv = Float32[cwhr_mod[i] == 40 ? 0f0 : Float32(cwhr_mod[i]) for i in 1:4]
        yv = Float32[cwhr_mod[i] == 40 ? 0f0 : cwhr_wt[i] for i in 1:4]
        fill!(cwhr_mod, 0); fill!(cwhr_wt, 0f0)
        ix = collect(1:4); rdpsrt!(4, xv, ix, true)
        for i in 1:2
            j = ix[i]; cwhr_mod[i] = Int(xv[j]); cwhr_wt[i] = yv[j]
        end
        xsum = cwhr_wt[1] + cwhr_wt[2]
        xsum > 0f0 && (cwhr_wt[1] /= xsum; cwhr_wt[2] /= xsum)
    end
    return (sz, dn, cwhr_mod, cwhr_wt)
end

"""
    nc_select_fuel_models(s, mois, sm, lg) -> Vector{(model, weight)}

NC California-CWHR fuel-model selection (nc/fmcfmd.f), the FMBURN candidate-model set + weights.
`sm`/`lg` are the SMALL/LARGE down-wood loads (the FMDYN point); `mois` unused here (kept for the
select_fuel_models signature parity). Region-5 forests (500-599 or ≥705) use FMD_R5, else FMD_R6.
"""
function nc_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees
    eqwt = zeros(Float32, 13)              # NC ICLSS = 13 fuel models (nc/fmcfmd.f)
    # per-tree crown widths (CRWDTH) — same load-time BAREA clamp as fmcba (cwcalc.f:859)
    ba = (s.control.cycle <= Int32(1)) ? 1f0 : s.plot.basal_area
    el = s.plot.elevation
    cws = zeros(Float32, t.n)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        cws[i] = nc_cwcalc(Int(t.species[i]), t.dbh[i], t.height[i],
                           Float32(t.crown_pct[i]), ba, el, 0f0)
    end

    # per-species BA% (BAPCT)
    bapct = zeros(Float32, 12); stndba = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); stndba += t.tpa[i] * t.dbh[i]^2 * 0.0054542f0
    end
    if stndba > 0.001f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            sp = Int(t.species[i])
            bapct[sp] += 100f0 * (t.tpa[i] * t.dbh[i]^2 * 0.0054542f0) / stndba
        end
    end

    # CWHR structural stage
    sz, dn, cwhr_mod, cwhr_wt = nc_cwhr(s, cws)
    jss = _nc_findjss(sz, dn)

    # forest type IFT (nc/fmcfmd.f:241-270)
    ift = 0
    if bapct[10] >= 80f0
        ift = 1                                            # pine (PP)
    elseif (bapct[9] + bapct[4]) >= 80f0                   # RF + WF
        if bapct[9] >= bapct[4]
            ift = 2                                         # red fir
        else
            ss = Int(s.plot.site_species)
            si = (1 <= ss <= length(s.plot.sp_site_index)) ? s.plot.sp_site_index[ss] : 0f0
            ift = si < 55f0 ? 3 : 4                         # white fir E / W
        end
    elseif (bapct[1] + bapct[3]) >= 80f0                   # OS + DF
        ift = 5                                             # Douglas-fir
    elseif (bapct[5] + bapct[7] + bapct[8] + bapct[11]) >= 80f0   # MA+BO+TO+OH
        ift = 6                                             # hardwoods
    elseif (bapct[2] + bapct[6]) >= 80f0                   # SP + IC
        ift = 8                                             # mixed conifer
    else
        ift = bapct[10] > bapct[4] ? 7 : 9                 # mixed pine / other softwood
    end

    # R5 (forests 500-599 or ≥705) vs R6 — identical tables, so the KODFOR split is inert here.
    kodfor = Int(s.plot.user_forest_code)
    cwhrfmd = ((500 <= kodfor < 600) || kodfor >= 705) ? _NC_FMD_R5 : _NC_FMD_R6

    # base model from (IFT, JSS)
    if 1 <= ift <= 9 && 1 <= jss <= 18
        eqwt[_nc_findmod(cwhrfmd[ift, jss])] = 1f0
    end

    # dynamic density blending (LDYNFM; only SS 2-5)
    if sz in ('2', '3', '4', '5')
        imd = zeros(Int, 2)
        for i in 1:2
            if cwhr_wt[i] > 0f0
                dc = cwhr_mod[i] == 10 ? 'S' : cwhr_mod[i] == 20 ? 'P' :
                     cwhr_mod[i] == 30 ? 'M' : cwhr_mod[i] == 40 ? 'D' : ' '
                jc = _nc_findjss(sz, dc)
                imd[i] = (1 <= ift <= 9 && 1 <= jc <= 18) ? _nc_findmod(cwhrfmd[ift, jc]) : 0
            end
        end
        if imd[1] == 0 && imd[2] == 0
            eqwt[8] = 1f0
        elseif imd[1] != 0 && imd[2] != 0 && imd[1] != imd[2]
            eqwt[imd[1]] = cwhr_wt[1]; eqwt[imd[2]] = cwhr_wt[2]
        end
    end

    # FM10 / FM11 sharing (nc/fmcfmd.f:339-357). The 5-yr post-activity FM11 (AFWT/SLCHNG/LATFUEL)
    # needs HARVYR + the small+large fuel-jump tracking, not yet in jl's FFE state — so FM11 stays 0
    # and FM10 = 1 (the AFWT=0 branch), which is exact once a burn is ≥5 yr after the last entry (nct01's
    # 2003 fire is 10 yr after the 1993 THINDBH ⇒ AFWT=0). TODO: wire HARVYR/SLCHNG for the <5-yr window.
    if eqwt[11] == 0f0
        eqwt[10] = 1f0
    else
        eqwt[10] = 0f0
    end
    eqwt[12] = 1f0
    eqwt[13] = 1f0

    return _fmdyn(sm, lg, eqwt, _NC_FMD_XPTS)
end
