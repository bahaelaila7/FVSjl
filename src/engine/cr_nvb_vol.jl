# =============================================================================
# cr_nvb_vol.jl — CR volume, NVB method (NVEL nsvb.f, National-Scale Volume & Biomass).
#
# Total cubic (TCF) = Vtotib (total stem INSIDE-bark cubic) from Table S1 (s1_ib.csv) via
# CalcVOLWT eqn forms 1-5. Merch cubic (MCF) = the NVBC log-bucked sawlog cubic (VOL(4),
# topwood VOL(7)=0 for the CUFT call since MTOPP==MTOPS): buck the stem stump→merch-top into
# NUMLOG/SEGMNT segments, sum per-log Smalian .00272708·(DIBL²+DIBS²)·LEN with NINT-integer
# DIBs, each rounded to 0.1 (nsvb.f NVB_CalcLOGVOL). Merch top MTOPP=MTOPS=TOPD·BARK inside
# bark (fvsvol.f:174), TOPD=4.0 (CR sitset IMODTY≠3), STUMP=1.0. Region-3 merch rule (mrules.f
# PROD 02): MAXLEN=16, MINLEN=MERCHL=10, TRIM=0.5, OPT=22, EVOD=2.
#
# The coefficient row is keyed by (SPCD, DIVISION, STDORG) parsed from VOLEQ the SPEQCOEF way
# (nsvb.f:903): SPCD=VOLEQ[8:10], DIVISION=int(VOLEQ[5:6])·10 (+1000 if VOLEQ[4]=='M'),
# STDORG = (len==11 && VOLEQ[11]=='P') ? 1 : 0 — with fallback to (SPCD, 0, STDORG).
# =============================================================================

const _NVB_S1 = Ref{Union{Nothing,Dict{Tuple{Int,Int,Int},NTuple{10,Float32}}}}(nothing)
const _NVB_S5 = Ref{Union{Nothing,Dict{Tuple{Int,Int,Int},NTuple{2,Float32}}}}(nothing)

function _nvb_load_table(fname::AbstractString)
    d = Dict{Tuple{Int,Int,Int},Vector{Float32}}()   # (spcd,div,stdorg) -> raw fields 4..13
    path = joinpath(@__DIR__, "..", "..", "data", "centralrockies", "nvb", fname)
    for (li, line) in enumerate(eachline(path))
        li == 1 && continue
        f = split(line, ',')
        length(f) < 13 && continue
        spcd = parse(Int, strip(f[1])); div = parse(Int, strip(f[2])); stdorg = parse(Int, strip(f[3]))
        eqn = round(Int, parse(Float32, f[4]))
        d[(spcd, div, stdorg)] = Float32[eqn, (parse(Float32, f[i]) for i in 5:13)...]  # eqn,a,a0,a1,b,b0,b1,b2,c,c1
    end
    return d
end

function _nvb_load_s1()
    _NVB_S1[] !== nothing && return _NVB_S1[]
    raw = _nvb_load_table("s1_ib.csv")
    d = Dict{Tuple{Int,Int,Int},NTuple{10,Float32}}()
    for (k, v) in raw; d[k] = ntuple(i -> v[i], 10); end
    _NVB_S1[] = d
    return d
end

function _nvb_load_s5()
    _NVB_S5[] !== nothing && return _NVB_S5[]
    raw = _nvb_load_table("s5_ratio.csv")
    d = Dict{Tuple{Int,Int,Int},NTuple{2,Float32}}()
    for (k, v) in raw; d[k] = (v[2], v[5]); end   # a (field5), b (field8)
    _NVB_S5[] = d
    return d
end

"SPEQCOEF (nsvb.f:903) — parse (SPCD, DIVISION, STDORG) from a VOLEQ string."
@inline function _nvb_voleq_key(voleq::AbstractString)
    spcd = tryparse(Int, voleq[8:10]); spcd === nothing && return (0, 0, 0)
    dc5 = voleq[5]; dc6 = voleq[6]
    d5 = (dc5 >= '0' && dc5 <= '9') ? (dc5 - '0') : 0
    d6 = (dc6 >= '0' && dc6 <= '9') ? (dc6 - '0') : 0
    div = (d5 * 10 + d6) * 10
    voleq[4] == 'M' && (div += 1000)
    stdorg = (length(voleq) >= 11 && voleq[11] == 'P') ? 1 : 0
    return (spcd, div, stdorg)
end

@inline function _nvb_lookup(tbl, spcd::Int, div::Int, stdorg::Int)
    row = get(tbl, (spcd, div, stdorg), nothing)
    row === nothing && (row = get(tbl, (spcd, 0, stdorg), nothing))
    return row
end

"CalcVOLWT (nsvb.f:862) — volume/weight from the equation form + coefficients."
@inline function _nvb_calcvolwt(spcd::Int, d::Float32, h::Float32, eqn::Int,
                                a, a0, a1, b, b0, b1, b2, c, c1, wdsg::Float32)::Float32
    eqn <= 0 && return 0f0
    if eqn == 1
        return a * fpow(d, b) * fpow(h, c)
    elseif eqn == 2
        k = spcd < 300 ? 9f0 : 11f0
        return d < k ? a0 * fpow(d, b0) * fpow(h, c) : a0 * fpow(k, b0 - b1) * fpow(d, b1) * fpow(h, c)
    elseif eqn == 3
        return a * fpow(d, a1 * fpow(1f0 - fexp(-b1 * d), c1)) * fpow(h, c)
    elseif eqn == 4
        return a * fpow(d, b) * fpow(h, c) * fexp(-(b2 * d))
    elseif eqn == 5
        return a * fpow(d, b) * fpow(h, c) * wdsg / 62.4f0
    end
    return 0f0
end

# --- NVBC taper + log-bucking kernels (nsvb.f) -------------------------------

"CalcRatio (nsvb.f:891): merch-fraction R = (1-(1-h1/H)^a)^b for 0<h1<=H."
@inline function _nvb_calcratio(H::Float32, h1::Float32, a::Float32, b::Float32)::Float32
    (h1 > 0f0 && h1 <= H) ? fpow(1f0 - fpow(1f0 - h1 / H, a), b) : 0f0
end

"NVB_CalcDiaAtHT (nsvb.f:1010): inside-bark diameter at height ht2 given total cubic tcuft."
@inline function _nvb_diaatht(tcuft::Float32, a::Float32, b::Float32, tht::Float32, ht2::Float32)::Float32
    ht2 >= tht && return 0f0
    x = 1f0 - ht2 / tht
    v = (tcuft / 0.005454154f0 / tht) * (a * b * fpow(x, a - 1f0) * fpow(1f0 - fpow(x, a), b - 1f0))
    v <= 0f0 && return 0f0
    return sqrt(v)
end

"NVB_CalcHT2TOPD (nsvb.f:971): bisection for the height where inside-bark dia == topd."
function _nvb_ht2topd(tcuft::Float32, a::Float32, b::Float32, httot::Float32, topd::Float32)::Float32
    low = 0f0; hi = httot; diff = 1f0; mid = 0f0; loop = 0
    while abs(diff) > 0.001f0
        mid = (low + hi) / 2f0
        mid < 0.5f0 && (mid = 0f0; break)
        x = 1f0 - mid / httot
        arg = tcuft / 0.005454154f0 / httot * a * b * fpow(x, a - 1f0) * fpow(1f0 - fpow(x, a), b - 1f0)
        dia = arg <= 0f0 ? 0f0 : fpow(arg, 0.5f0)
        diff = topd - dia
        abs(diff) < 0.001f0 && break
        diff < 0f0 ? (low = mid) : (hi = mid)
        loop += 1; loop > 1000 && break
    end
    return mid
end

@inline _nint(x::Float32) = Float32(floor(Int, x + 0.5f0))    # Fortran NINT (x>=0 here)
@inline _anint10(x::Float32) = floor(x * 10f0 + 0.5f0) / 10f0  # ANINT(x*10)/10

"NUMLOG (numlog.f): number of merch segments in a stem of length lmerch (OPT/EVOD merch rule)."
function _nvb_numlog(opt::Int, evod::Int, lmerch::Float32, maxlen::Float32, minlen::Float32, trim::Float32)::Int
    numseg = trunc(Int, lmerch / (maxlen + trim))
    leftov = lmerch - (maxlen + trim) * numseg
    if numseg > 0 || leftov >= minlen
        if opt < 20
            leftov >= trim + 0.5f0 && (numseg += 1)
        elseif opt == 21 || opt == 22
            evod == 1 && leftov >= trim + 0.5f0 && (numseg += 1)
            evod == 2 && leftov >= trim + 1.0f0 && (numseg += 1)
        elseif opt == 23
            leftov >= trim + minlen && (numseg += 1)
        elseif opt == 24
            leftov >= (maxlen + trim) / 4.0f0 && (numseg += 1)
        end
    else
        numseg = 0
    end
    numseg > 20 && (numseg = 20)
    return numseg
end

"SEGMNT (segmnt.f): per-segment lengths for the nominal-log options (OPT>=20). Returns (loglen, numseg)."
function _nvb_segmnt(opt::Int, evod::Int, lmerch0::Float32, maxlen::Float32, minlen::Float32, trim::Float32, numseg::Int)
    loglen = zeros(Float32, 20)
    numseg == 0 && return (loglen, 0)
    lmerch = lmerch0 - numseg * trim
    lmerch = evod == 1 ? Float32(trunc(Int, lmerch + 0.5f0)) : Float32(trunc(Int, (lmerch + 1.0f0) / 2.0f0) * 2.0f0)
    lmerch > numseg * maxlen && (lmerch = Float32(numseg) * maxlen)
    if numseg == 1
        if opt == 24
            if lmerch < maxlen * 0.25f0
                loglen[1] = 0f0
            elseif lmerch <= maxlen * 0.75f0
                loglen[1] = maxlen / 2f0
            else
                loglen[1] = maxlen
            end
        elseif lmerch >= minlen
            lmerch > maxlen && (lmerch = maxlen)
            loglen[1] = lmerch
        else
            loglen[1] = 0f0
        end
        return (loglen, numseg)
    end
    if opt < 20
        # even/odd distribution (not used for CR OPT=22, ported for completeness)
        avlen = trunc(Int, lmerch / numseg)
        leftov = lmerch - avlen * numseg
        for i in 1:numseg; loglen[i] = avlen; end
        if avlen > trunc(Int, avlen / 2f0) * 2
            for i in 1:numseg
                if numseg - 2i + 1 >= 1
                    loglen[i] += 1f0; loglen[numseg - i + 1] -= 1f0
                end
            end
        end
        # (leftover redistribution omitted: CR never uses OPT<20)
        return (loglen, numseg)
    end
    # nominal-log OPT>=20
    leftov = lmerch - trunc(Int, maxlen) * (numseg - 1)
    for i in 1:numseg; loglen[i] = maxlen; end
    if opt == 21
        if leftov >= maxlen / 2f0
            loglen[numseg] = leftov
        else
            loglen[numseg] = Float32(trunc(Int, (maxlen + leftov) / 2f0))
            loglen[numseg - 1] = maxlen + leftov - loglen[numseg]
            if loglen[numseg] == loglen[numseg - 1] && loglen[numseg] > trunc(Int, loglen[numseg] / 2f0) * 2f0
                loglen[numseg] -= 1f0; loglen[numseg - 1] += 1f0
            end
        end
    elseif opt == 22
        loglen[numseg] = Float32(trunc(Int, (maxlen + leftov) / 2f0))
        loglen[numseg - 1] = maxlen + leftov - loglen[numseg]
        if loglen[numseg] < minlen
            loglen[numseg] = 0f0; loglen[numseg - 1] = maxlen; numseg -= 1
        else
            if loglen[numseg] == loglen[numseg - 1] && loglen[numseg] > trunc(Int, loglen[numseg] / 2f0) * 2f0
                loglen[numseg] -= 1f0; loglen[numseg - 1] += 1f0
            end
        end
    elseif opt == 23
        if leftov >= minlen
            loglen[numseg] = leftov
        else
            loglen[numseg] = 0f0; numseg -= 1
        end
    elseif opt == 24
        if leftov < maxlen * 0.25f0
            loglen[numseg] = 0f0; numseg -= 1
        elseif leftov <= maxlen * 0.75f0
            loglen[numseg] = Float32(trunc(Int, maxlen * 0.5f0 + 0.5f0))
        else
            loglen[numseg] = maxlen
        end
    end
    return (loglen, numseg)
end

"NVB_CalcLOGVOL cubic accumulation (nsvb.f:1044): Σ per-log Smalian cuft (0.1-rounded, NINT DIBs)."
function _nvb_logvol_cuft(numseg::Int, loglen::Vector{Float32}, dibl0::Float32, ht2_0::Float32,
                          tcuft::Float32, trim::Float32, tht::Float32, a::Float32, b::Float32)::Float32
    ht2 = ht2_0; dibl = dibl0; vol4 = 0f0
    @inbounds for i in 1:numseg
        ht2 = ht2 + trim + loglen[i]
        dib = _nvb_diaatht(tcuft, a, b, tht, ht2)
        dibs = _nint(dib)
        logcv = 0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen[i]
        vol4 += _anint10(logcv)
        dibl = dibs
    end
    return vol4
end

# CR region-3 CUFT merch rule (mrules.f, PROD 02): the primary-product cubic bucking constants.
const _NVB_R3_MAXLEN = 16.0f0
const _NVB_R3_MINLEN = 10.0f0
const _NVB_R3_MERCHL = 10.0f0
const _NVB_R3_TRIM   = 0.5f0
const _NVB_R3_OPT    = 22
const _NVB_R3_EVOD   = 2

"NVBC merch cubic (VOL(4)) to the inside-bark top `mtop`, given total inside-bark cubic `vtotib`
and stem taper coefs (a,b). Returns the log-bucked sawlog cubic (VOL(7) topwood is 0 when the
primary and secondary tops coincide, as in the CUFT call)."
function _nvb_merch_cuft(d::Float32, h::Float32, vtotib::Float32, stump::Float32, mtop::Float32,
                         a::Float32, b::Float32)::Float32
    vtotib <= 0f0 && return 0f0
    ht1prd = mtop < d ? _nvb_ht2topd(vtotib, a, b, h, mtop) : 0f0
    ht1prd < stump && (ht1prd = stump)
    lmerch = ht1prd - stump
    lmerch < 0f0 && (lmerch = 0f0)
    lmerch < _NVB_R3_MERCHL && return 0f0
    dibl = _nint(_nvb_diaatht(vtotib, a, b, h, 4.5f0))
    numseg = _nvb_numlog(_NVB_R3_OPT, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, _NVB_R3_MINLEN, _NVB_R3_TRIM)
    numseg == 0 && return 0f0
    loglen, numseg = _nvb_segmnt(_NVB_R3_OPT, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, _NVB_R3_MINLEN, _NVB_R3_TRIM, numseg)
    return _nvb_logvol_cuft(numseg, loglen, dibl, stump, vtotib, _NVB_R3_TRIM, h, a, b)
end

"CR NVB per-tree volume. `voleq`=NVB eq id (division/stdorg parsed from it). `bark`=DIB/DOB ratio,
`topd`=cubic top DOB (4.0), `stump`=stump ht. Returns a 15-vec: VOL[1]=Vtotib=TCF, VOL[4]=merch cubic."
function cr_nvb_vol(voleq::AbstractString, d::Float32, h::Float32; bark::Float32 = 1f0,
                    topd::Float32 = 4f0, stump::Float32 = 1f0, wdsg::Float32 = 0f0)
    vol = zeros(Float32, 15)
    (d < 1f0 || h < 5f0) && return vol
    spcd, div, stdorg = _nvb_voleq_key(voleq)
    spcd == 0 && return vol
    r1 = _nvb_lookup(_nvb_load_s1(), spcd, div, stdorg)
    r1 === nothing && return vol
    eqn = round(Int, r1[1])
    vib = _nvb_calcvolwt(spcd, d, h, eqn, r1[2], r1[3], r1[4], r1[5], r1[6], r1[7], r1[8], r1[9], r1[10], wdsg)
    vol[1] = vib
    # merch cubic (VOL(4)): CUFT call top = TOPD·BARK inside bark (fvsvol.f:174)
    r5 = _nvb_lookup(_nvb_load_s5(), spcd, div, stdorg)
    if r5 !== nothing && vib > 0f0
        mtop = topd * bark
        vol[4] = _nvb_merch_cuft(d, h, vib, stump, mtop, r5[1], r5[2])
    end
    return vol
end
