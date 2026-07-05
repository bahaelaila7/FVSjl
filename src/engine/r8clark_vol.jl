# r8clark_vol.jl — R8 New Clark Volume Equations (R9CLARK path)
# Translated from:
#   r8dib.inc, r8clkcoef.inc, r8cfo.inc  — coefficient tables
#   r8prep.f  — coefficient lookup / DIB preparation
#   r9clark.f — cubic volume and board-feet computation
#   r9logs.f  — log segmentation
#   mrules.f  — R8 merchantability defaults (MAXLEN=8, TRIM=0.5, etc.)
#
# Entry point: R8CLARK_VOL(voleq, dbhOb, htTot, ispc) → vol[15]
# Used by NATCRS in fvsvol.jl when VEQNNC matches "8{d}1CLKE{sss}".

# ---------------------------------------------------------------------------
# Scribner Decimal C board-foot table (scrbnr[120]) from r9clark.f
# ---------------------------------------------------------------------------
# R8 Clark volume coefficient tables — bulky parameter blobs loaded from CSV (data/southern/volume/)
# rather than inlined, to keep them out of the source. Each CSV is row-major; the flat Float32 vector
# preserves the original `(row-1)*ncol + col` access. `value_col` reads only the 2nd column (for the
# 1-D Scribner lookup, whose CSV is index,value).
const _CLARK_DIR = joinpath(@__DIR__, "..", "..", "data", "southern", "volume")
function _load_clark_flat(file; value_col::Bool=false)
    vals = Float32[]
    for (k, line) in enumerate(eachline(joinpath(_CLARK_DIR, file)))
        k == 1 && continue
        toks = split(strip(line), ',')
        isempty(toks[1]) && continue
        if value_col
            push!(vals, parse(Float32, toks[2]))
        else
            for t in toks; push!(vals, parse(Float32, t)); end
        end
    end
    return vals
end

const _SCRBNR = _load_clark_flat("scribner_bdft.csv"; value_col=true)

# R8CF[182,18] — GEOA/SPCD inside-bark taper coefficients from r8dib.inc
const _R8CF = _load_clark_flat("r8clark_cf.csv"; value_col=false)
# Access: _R8CF[(row-1)*18 + col]  (1-based row/col)

# R8CFO[182,9] — outside-bark DIB17 coefficients from r8cfo.inc
const _R8CFO = _load_clark_flat("r8clark_cfo.csv"; value_col=false)

# DIBMEN[49,3] — col1=SPCD, col2=plp fixDI, col3=saw fixDI
const _DIBMEN = _load_clark_flat("r8clark_dibmen.csv"; value_col=false)

# TOTAL[49,7] — inside-bark taper coef: col1=SPCD, col2=R,C,E,P,B,A
const _TOTAL = _load_clark_flat("r8clark_total.csv"; value_col=false)

# OTOTAL[49,7] — outside-bark taper coef
const _OTOTAL = _load_clark_flat("r8clark_ototal.csv"; value_col=false)

# ---------------------------------------------------------------------------
# Accessor helpers (1-based row/col, row-major storage)
# ---------------------------------------------------------------------------
@inline _r8cf(row,col)  = _R8CF[(row-1)*18 + col]
@inline _r8cfo(row,col) = _R8CFO[(row-1)*9  + col]
@inline _dibmen(row,col)= _DIBMEN[(row-1)*3  + col]
@inline _total(row,col) = _TOTAL[(row-1)*7   + col]
@inline _ototal(row,col)= _OTOTAL[(row-1)*7  + col]

# ---------------------------------------------------------------------------
# R8CLARK species remapping (r8prep.f lines 99-116)
# ---------------------------------------------------------------------------
function _r8_remap_spec(spec::Int)
    (spec == 123 || spec == 197)          && return 100   # parens: && binds tighter than ||
    spec == 268                            && return 261
    spec ∈ (313,314,317,650,651,691,711,742,762,920,930,545,546) && return 300
    spec ∈ (521,550,580,601,602,318)      && return 500
    spec ∈ (804,817,820,823,825,826,830,834) && return 800
    return spec
end

# ---------------------------------------------------------------------------
# Coefficient lookup: returns (R,C,E,P,B,A, Ro,Co,Eo,Po,Bo,Ao,
#                               A4,B4,AFI,BFI, A17,B17, A17o,B17o,
#                               SPGRP, sppIdx, ptr, errFlg)
# for the htTot>0 case (always used from NATCRS).
# ---------------------------------------------------------------------------
function _r8clark_lookup(voleq::AbstractString)
    errFlg = 0
    geoa = 0
    spec = 0
    regn = 0
    # Parse voleq: "8{geoa}1CLKE{sss}" → geoa from pos 2, spec from pos 8-10
    length(voleq) >= 10 || return nothing, 1
    geoa_c = voleq[2]
    spec_s = voleq[8:10]
    try
        geoa = parse(Int, string(geoa_c))
        spec = parse(Int, spec_s)
    catch
        return nothing, 6
    end
    (geoa < 1 || geoa > 9 || geoa == 8) && return nothing, 1

    spec = _r8_remap_spec(spec)

    # Sequential search in DIBMEN for sppIdx
    sppIdx = 0
    for k in 1:49
        s = Int(_dibmen(k,1))
        if s == spec
            sppIdx = k; break
        elseif s > spec
            break
        end
    end
    sppIdx == 0 && return nothing, 1

    fixDI4  = _dibmen(sppIdx, 2)   # pulpwood min dib (unused here but kept)
    fixDI79 = _dibmen(sppIdx, 3)   # sawlog min dib

    # Binary search in _R8CF for (geoa*1000 + spec)
    gspec = geoa*1000 + spec
    first_, last_ = 1, 182
    ptr = 0
    done = false
    lastflag = false
    while !done
        if first_ == last_; lastflag = true; end
        half = ((last_ - first_ + 1) ÷ 2) + first_
        check = Int(_r8cf(half,1))*1000 + Int(_r8cf(half,2))
        if gspec == check
            ptr = half; done = true
        elseif gspec > check
            first_ = half
        else
            last_ = max(half - 1, first_)
        end
        if lastflag && !done
            if geoa == 9
                return nothing, 6
            else
                geoa = 9
                gspec = 9000 + spec
                first_ = 135; last_ = 182
                lastflag = false
            end
        end
    end
    ptr == 0 && return nothing, 6

    spgrp = round(Int, _r8cf(ptr, 3))
    (spgrp != 100 && spgrp != 300 && spgrp != 500) && return nothing, 6

    # Inside-bark taper coefficients from TOTAL
    R  = _total(sppIdx, 2)
    C  = _total(sppIdx, 3)
    E  = _total(sppIdx, 4)
    P  = _total(sppIdx, 5)
    B  = _total(sppIdx, 6)
    A  = _total(sppIdx, 7)

    # Outside-bark taper coefficients from OTOTAL
    Ro = _ototal(sppIdx, 2)
    Co = _ototal(sppIdx, 3)
    Eo = _ototal(sppIdx, 4)
    Po = _ototal(sppIdx, 5)
    Bo = _ototal(sppIdx, 6)
    Ao = _ototal(sppIdx, 7)

    # DBH → inside-bark DBH coefficients
    A4  = _r8cf(ptr, 4)
    B4  = _r8cf(ptr, 5)

    # Form class coefficients (inside-bark DIB17 from outside-bark DIB17)
    AFI = _r8cf(ptr, 6)
    BFI = _r8cf(ptr, 7)

    # A17/B17 for htTot > 0 case (inside-bark)
    A17  = _r8cf(ptr, 14)
    B17  = _r8cf(ptr, 15)

    # A17/B17 for htTot > 0 case (outside-bark, from R8CFO col 4-5)
    A17o = _r8cfo(ptr, 4)
    B17o = _r8cfo(ptr, 5)

    return (R=R, C=C, E=E, P=P, B=B, A=A,
            Ro=Ro, Co=Co, Eo=Eo, Po=Po, Bo=Bo, Ao=Ao,
            A4=A4, B4=B4, AFI=AFI, BFI=BFI,
            A17=A17, B17=B17, A17o=A17o, B17o=B17o,
            spgrp=spgrp, sppIdx=sppIdx, ptr=ptr), 0
end

# ---------------------------------------------------------------------------
# r9cuft: 3-segment Clark taper cubic volume integral
# Translated from r9clark.f lines 961-1113.
# Uses inside-bark coefficients: R,C,E,P,B,A,totHt,dbhIb,dib17
# Returns cubic feet.
# ---------------------------------------------------------------------------
function _r9cuft(R::Real, C::Real, E::Real, P::Real, B::Real, A::Real,
                 totHt::Real, dbhIb::Real, dib17::Real,
                 lowrHt::Real, upprHt::Real)
    upprHt <= 0 && return 0.0f0
    dbhIb <= 0  && return 0.0f0
    totHt <= 0  && return 0.0f0

    G = (1 - 4.5/totHt)^R
    W = (C + E/dbhIb^3) / (1 - G)
    X = (1 - 4.5/totHt)^P

    # Avoid underflow: (1-17.3/totHt)^P when near tip
    Y = if (1 - 17.3/totHt) < 0.005748 && P > 14
        0.0f0
    else
        (1 - 17.3/totHt)^P
    end

    Z = (X - Y) > 1e-10 ? (dbhIb^2 - dib17^2) / (X - Y) : 0.0f0
    T = dbhIb^2 - Z * X

    L1 = max(lowrHt, 0.0); U1 = min(upprHt, 4.5)
    L2 = max(lowrHt, 4.5); U2 = min(upprHt, 17.3)
    L3 = max(lowrHt, 17.3); U3 = min(totHt,  upprHt)

    I1 = lowrHt < 4.5
    I2 = lowrHt < 17.3
    I3 = upprHt > 4.5
    I4 = upprHt > 17.3
    I5 = (L3 - 17.3) < A * (totHt - 17.3)
    I6 = (U3 - 17.3) < A * (totHt - 17.3)

    V1 = 0.0f0
    if I1
        t1 = (1 - L1/totHt)^R * (totHt - L1)
        t2 = (1 - U1/totHt)^R * (totHt - U1)
        V1 = dbhIb^2 * ((1 - G*W)*(U1-L1) + W*(t1-t2)/(R+1))
    end

    V2 = 0.0f0
    if I2 && I3
        t1p = (1 - L2/totHt)^P * (totHt - L2)
        t2p_term = (1 - U2/totHt)
        if t2p_term < 0.005748 && P > 14
            V2 = T*(U2-L2) + Z*t1p/(P+1)
        else
            t2p = t2p_term^P * (totHt - U2)
            V2 = T*(U2-L2) + Z*(t1p - t2p)/(P+1)
        end
    end

    V3 = 0.0f0
    if I4
        dth = totHt - 17.3
        V3 = dib17^2 * (
            B*(U3-L3)
            - B*((U3-17.3)^2 - (L3-17.3)^2) / dth
            + (B/3)*((U3-17.3)^3 - (L3-17.3)^3) / dth^2
            + (I5 ? (1/3)*((1-B)/A^2)*(A*dth - (L3-17.3))^3 / dth^2 : 0.0)
            - (I6 ? (1/3)*((1-B)/A^2)*(A*dth - (U3-17.3))^3 / dth^2 : 0.0)
        )
    end

    cfVol = Float32(0.005454154 * (V1 + V2 + V3))
    return max(cfVol, 0.0f0)
end

# ---------------------------------------------------------------------------
# r9ht: height from inside/outside-bark diameter.
# Translated from r9clark.f lines 1244-1358.
# For R8 we call this with outside-bark coefficients (Ro,Co,Eo,Po,Bo,Ao,
#   totHt, dbhOb, dob17) and stmDib = outside-bark top diameter.
# ---------------------------------------------------------------------------
# FVS R9HT (r9clark.f:1267-1360) — height to a given inside-bark diameter, Float32 (`REAL*4`) throughout
# to match Fortran's op sequence (real Clark powers via Float32 `^`, `**(1/r)`/`**(1/p)`/`**0.5`, integer
# exponents as exact multiplications). The `sawHt`/`plpHt` this returns drive the even-foot LOG segmentation
# (LEFTOV → INT rounding in R9LOGLEN); computing in Float64 and rounding once makes sawHt 1 ULP off FVS's
# Float32, which at large trees tips LEFTOV across an even-foot INT knife-edge → different log tops/DIBs →
# a whole Scribner board-foot step. Faithful Float32 keeps the segmentation bit-exact.
function _r9ht(R::Real, C::Real, E::Real, P::Real, B::Real, A::Real,
               totHt::Real, dbhIb::Real, dib17::Real, stmDib::Real)
    R = Float32(R); C = Float32(C); E = Float32(E); P = Float32(P); B = Float32(B); A = Float32(A)
    totHt = Float32(totHt); dbhIb = Float32(dbhIb); dib17 = Float32(dib17); stmDib = Float32(stmDib)
    G = fpow(1f0 - 4.5f0/totHt, R)                       # Clark real powers via gfortran companion (doctrine #8)
    W = (C + E/dbhIb^3) / (1f0 - G)
    X = fpow(1f0 - 4.5f0/totHt, P)
    Y = fpow(1f0 - 17.3f0/totHt, P)
    Z = (X - Y) > 1f-10 ? (dbhIb^2 - dib17^2) / (X - Y) : 0.0f0

    Im = stmDib^2 > B*(A-1f0)^2*dib17^2 ? 1.0f0 : 0.0f0
    Qa =  B + Im*(1f0-B)/A^2
    Qb = -2f0*B - Im*2f0*(1f0-B)/A
    Qc =  B + (1f0-B)*Im - stmDib^2/dib17^2

    Is = stmDib >= dbhIb ? 1.0f0 : 0.0f0
    Ib = (stmDib < dbhIb && stmDib >= dib17) ? 1.0f0 : 0.0f0
    It = stmDib < dib17 ? 1.0f0 : 0.0f0

    stemHt = 0.0f0
    if Is > 0f0
        xxx = (stmDib^2/dbhIb^2 - 1f0)/W + G
        xxx > 0f0 && (stemHt = totHt*(1f0 - fpow(xxx, 1f0/R)))
    elseif Ib > 0f0
        xxx = X - (dbhIb^2 - stmDib^2)/Z
        xxx > 0f0 && (stemHt = totHt*(1f0 - fpow(xxx, 1f0/P)))
    else
        xxx = Qb^2 - 4f0*Qa*Qc
        xxx > 0f0 && (stemHt = 17.3f0 + (totHt - 17.3f0)*((-Qb - fpow(xxx, 0.5f0))/(2f0*Qa)))
    end
    return max(stemHt, 0.0f0)
end

# ---------------------------------------------------------------------------
# R8 form-class minimum DIB17 (r8prep.f lines 346-365)
# ---------------------------------------------------------------------------
function _r8_fcmin_adj(dib17::Float32, dbhOb::Float32, htTot::Float32, spgrp::Int, spec::Int)
    (spec == 221 || spec == 222 || spec == 544) && return dib17
    fcmin = if spgrp == 100
        htTot < 32.5 ? 56 : htTot < 37.5 ? 64 : htTot < 42.5 ? 66 : 67
    elseif spgrp == 300
        htTot < 32.5 ? 57 : htTot < 37.5 ? 60 : htTot < 42.5 ? 64 : 67
    else
        htTot < 32.5 ? 58 : htTot < 37.5 ? 65 : htTot < 42.5 ? 67 : 69
    end
    fcdib = dbhOb * Float32(fcmin) * 0.01f0
    htTot < 47.5f0 && dib17 < fcdib && return fcdib
    return dib17
end

# ---------------------------------------------------------------------------
# Main entry point: compute vol[15] for one tree using R8 Clark equations.
# Called from NATCRS when VEQNNC[ispc] starts with '8' and contains "CLKE".
#
# Arguments:
#   voleq  — e.g. "841CLKE131"
#   dbhOb  — DBH outside bark (inches)
#   htTot  — total height (feet)
#   mTopp  — sawtimber top diam ob (inches); 0 → use default (7 softwood / 9 hardwood)
#   mTopS  — pulpwood top diam ob (inches); 0 → use default (4)
#   stump  — stump height (feet); 0 → use 1.0
#   prod   — '01' sawtimber, '02' pulpwood
# Returns vol[15] array (same indices as NVEL vol array).
# ---------------------------------------------------------------------------
function _R8CLARK_VOL(voleq::AbstractString, dbhOb::Float32, htTot::Float32,
                      mTopp::Float32, mTopS::Float32, stump::Float32,
                      prod::AbstractString; log_dib::Union{Nothing,Base.RefValue{Dict{Int,Float32}}} = nothing,
                      log_cuft::Union{Nothing,Base.RefValue{Dict{Int,Float32}}} = nothing,
                      intl_bf::Bool = false)
    # `intl_bf`: report vol[10] as INTERNATIONAL ¼" board feet instead of Scribner. FVS volinit2.f:291-297
    # replaces vol(2) with vol(10) (International) for R8 National Forests IFORST∈{4,5,8,11,12,14,19,20,21,
    # 22,24,30}; other R8 forests keep Scribner. Same even-foot sawtimber bucking, different per-log rule.
    # `log_dib` (a Ref) is the opt-in per-log-DIB Scribner BF breakdown for the log-graded HRVRVN report
    # (#38): when supplied, the board-feet block below ALSO fills it via `_r8_scribner_bf_by_dib` with the
    # SAME merch params. nothing (the default for all non-econ callers) ⇒ no extra work; return arity unchanged.
    # `log_cuft` is the cubic (HRVRVN unit 5, FT3_100_LOG) analog — per-log gross cuft via `_r8_cuft_by_dib`
    # (R9LGCFT), renormalized to VOL(4)+VOL(7); filled in the same sawtimber block.

    vol = zeros(Float32, 15)
    # Returns (vol, sawHt, plpHt): the merch heights to the sawtimber and pulpwood
    # top diameters (= HT1PRD/HT2PRD in fvsvol.f), used to set HT2TD (merch top ht).
    dbhOb < 1.0f0 && return vol, 0.0f0, 0.0f0   # volinit.f 168: DBH<1 → no volume

    # Short-tree handling (r8prep.f 328-340 / r9totHt): for htTot < 17.4 the
    # taper model uses topHt=17.4 and cubic volumes are scaled by shrtHt/17.3.
    short  = htTot < 17.4f0
    shrtHt = htTot

    coef, err = _r8clark_lookup(voleq)
    err != 0 && return vol, 0.0f0, 0.0f0

    (; R, C, E, P, B, A, Ro, Co, Eo, Po, Bo, Ao,
       A4, B4, AFI, BFI, A17, B17, A17o, B17o, spgrp, sppIdx, ptr) = coef

    spec = Int(_dibmen(sppIdx, 1))

    # Compute inside-bark DBH at 4.5', floored to species min FIXDI4 (r8prep.f 291)
    fixdi4 = _dibmen(sppIdx, 2)
    dbhIb = A4 + B4*dbhOb
    if dbhIb < fixdi4; dbhIb = fixdi4; end

    # totHt = htTot for normal trees; = 17.4 (topHt) for short trees (r9totHt)
    totHt = short ? 17.4f0 : htTot

    # DIB17: inside-bark diameter at 17.3' from dbhOb
    dib17 = dbhOb*(A17 + B17*(17.3f0/totHt)^2)
    dib17 = max(dib17, 0.1f0)

    # Form class minimum adjustment
    dib17 = _r8_fcmin_adj(dib17, dbhOb, totHt, spgrp, spec)

    # Secondary-coefficient DIB17 (COEFFSO%DIB17, r8prep.f:366 + the :507 floor). For 221/222/544
    # (baldcypress/pondcypress/green-ash) live SKIPS the (FCLSS−AFI)/BFI step (r8prep.f:346) so
    # COEFFSO%DIB17 stays 0 and the :507 floor sets it = COEFFS%DIB17 (=dib17); every other species
    # uses (dib17−AFI)/BFI. Applying the :507 floor unconditionally is faithful — it's a no-op for the
    # other species (BFI<1 ⇒ (dib17−AFI)/BFI > dib17) and yields dib17 for the special three.
    dob17 = (spec == 221 || spec == 222 || spec == 544) ? dib17 : (dib17 - AFI)/BFI
    dob17 = max(dob17, dib17)                # r8prep.f:507  COEFFSO%DIB17 = max(.., COEFFS%DIB17)
    dob17 = max(dob17, 0.1f0)

    # Merch defaults from R8 MRULES (mrules.f line 337-369)
    stump  = stump > 0 ? stump : (prod == "01" ? 1.0f0 : 0.5f0)
    maxLen = 8.0f0
    minLen = 2.0f0
    merchL = prod == "08" ? 12.0f0 : 8.0f0   # min merch length (mrules.f 344-345)
    trim   = 0.5f0
    sawDib = mTopp > 0 ? mTopp : (spec < 300 ? 7.0f0 : 9.0f0)
    plpDib = mTopS > 0 ? mTopS : 4.0f0

    isProd1 = prod == "01"

    # 1. Total cubic volume (stump to tip) — r9clark.f 221-227
    vol[1] = _r9cuft(R,C,E,P,B,A, totHt,dbhIb,dib17, stump, totHt)
    if short; vol[1] *= shrtHt/17.3f0; end

    # 2. Pulpwood height (height to plpDib using outside-bark coefs)
    #    r9clark.f 260: zeroed if shorter than merchL+stump+trim
    plpHt = min(_r9ht(Ro,Co,Eo,Po,Bo,Ao, totHt,dbhOb,dob17, plpDib), totHt)
    if plpHt < merchL + stump + trim
        plpHt = 0.0f0
    end

    sawHt = 0.0f0; rawSawHt = 0.0f0
    if isProd1
        # 3. Sawtimber height (height to sawDib using outside-bark coefs). `rawSawHt` is the
        # un-truncated height to the merch top (= VOLINIT's HT1PRD), needed by the caller for the
        # Region-8 "≥10 ft of product" rule (fvsvol.f); `sawHt` is the merch-length-zeroed version
        # used for the cubic/board volumes here.
        rawSawHt = min(_r9ht(Ro,Co,Eo,Po,Bo,Ao, totHt,dbhOb,dob17, sawDib), totHt)
        sawHt = rawSawHt
        if sawHt < merchL + stump + trim   # r9clark.f 320
            sawHt = 0.0f0
        end

        # 4. Sawtimber cubic (vol[4]) — r9clark.f 349-354
        if sawHt > stump
            vol[4] = _r9cuft(R,C,E,P,B,A, totHt,dbhIb,dib17, stump, sawHt)
            if short; vol[4] *= shrtHt/17.3f0; end
        end

        # 5. Topwood cubic (vol[7] = merch to plpHt minus saw part)
        if plpHt > sawHt && plpHt > stump
            plpCf = _r9cuft(R,C,E,P,B,A, totHt,dbhIb,dib17, stump, plpHt)
            if short; plpCf *= shrtHt/17.3f0; end
            vol[7] = max(plpCf - vol[4], 0.0f0)
        end
    else
        # Pulpwood only: vol[4] = total merch cubic — r9clark.f 264-276
        if plpHt - stump >= minLen
            vol[4] = _r9cuft(R,C,E,P,B,A, totHt,dbhIb,dib17, stump, plpHt)
            if short; vol[4] *= shrtHt/17.3f0; end
        end
    end

    # 6. Board feet — Scribner by default, or International ¼" (vol(10)) for the R8 forests that use it.
    if isProd1 && sawHt > stump
        vol[10] = intl_bf ?
            _r8_intlqtr_bf(R,C,E,P,B,A, totHt,dbhIb,dib17, sawHt, stump, minLen, maxLen, trim) :
            _r8_scribner_bf(R,C,E,P,B,A, totHt,dbhIb,dib17,
                                   sawHt, plpHt, stump, minLen, maxLen, trim)
        log_dib !== nothing && (log_dib[] = _r8_scribner_bf_by_dib(R,C,E,P,B,A, totHt,dbhIb,dib17,
                                   sawHt, plpHt, stump, minLen, maxLen, trim))
        log_cuft !== nothing && (log_cuft[] = _r8_cuft_by_dib(R,C,E,P,B,A, totHt,dbhIb,dib17,
                                   sawHt, plpHt, stump, minLen, maxLen, trim, vol[4] + vol[7]))
    end

    # 7. Stump volume (vol[14]) and tip volume (vol[15])
    vol[14] = _r9cuft(R,C,E,P,B,A, totHt,dbhIb,dib17, 0.0f0, stump)
    vol[15] = max(vol[1] - vol[4] - vol[7], 0.0f0)

    # Round to match Fortran NINT() = round-half-AWAY-from-zero (NOT Julia round()'s ties-to-even).
    vol[1]  = Float32(round(vol[1]*10, RoundNearestTiesAway)/10)
    vol[4]  = Float32(round(vol[4]*10, RoundNearestTiesAway)/10)
    vol[7]  = Float32(round(vol[7]*10, RoundNearestTiesAway)/10)
    vol[10] = Float32(round(vol[10], RoundNearestTiesAway))

    return vol, rawSawHt, plpHt
end

# ---------------------------------------------------------------------------
# International ¼" board feet (r9bdft vol(10), r9clark.f:1482) over the SAWTIMBER section [stump, sawHt].
# Same even-foot bucking as `_r8_scribner_bf` (R9LOGLEN) + the R8 Clark log-top DIB, but the per-log board
# is the International rule (`_r9_intl_log`, shared with the R9/NE path) instead of the Scribner table.
# FVS uses this (volinit2.f:296 VOL(2)=VOL(10)) for R8 National Forests {4,5,8,11,12,14,19,20,21,22,24,30}.
# ---------------------------------------------------------------------------
function _r8_intlqtr_bf(R, C, E, P, B, A, totHt, dbhIb, dib17,
                        sawHt, stump, minLen, maxLen, trim)
    lmerch = sawHt - stump
    nlogp  = clamp(floor(Int, lmerch / (maxLen + trim)), 0, 39)
    leftov = lmerch - (maxLen + trim) * nlogp - trim
    logLen = zeros(Float32, 40); tlogs = 0
    if !(lmerch < minLen + trim || (nlogp == 0 && leftov < minLen + trim))
        for i in 1:nlogp; logLen[i] = maxLen; end
        if leftov >= minLen + trim
            nlogp += 1; logLen[nlogp] = leftov
        end
        if nlogp == 1
            logLen[1] = Float32(floor(Int, logLen[1]) ÷ 2 * 2)
        elseif leftov < minLen
            logLen[nlogp] = Float32(floor(Int, logLen[nlogp]) ÷ 2 * 2)
        else
            combined = maxLen + leftov
            logLen[nlogp]   = Float32(floor(Int, combined / 2) ÷ 2 * 2)
            logLen[nlogp-1] = Float32((floor(Int, combined - logLen[nlogp]) ÷ 2) * 2)
        end
        tlogs = nlogp
    end
    tlogs == 0 && return 0f0
    bf = 0f0; ht = stump
    for i in 1:tlogs
        len = logLen[i]
        ht += trim + len                                            # top (small end) of log i
        idib = trunc(Int, _r9dib_clark(R,C,E,P,B,A, totHt,dbhIb,dib17, ht) + 0.499f0)  # r9logdib INT(DIB+0.499)
        bf += _r9_intl_log(len, idib)
    end
    return round(bf, RoundNearestTiesAway)                          # r9bdft:1499 vol(10)=NINT
end

# ---------------------------------------------------------------------------
# Scribner board-foot volume using r9logs segmentation + r9bdft table lookup.
# Translated from r9logs.f and r9clark.f r9bdft subroutine.
# ---------------------------------------------------------------------------
function _r8_scribner_bf(R, C, E, P, B, A, totHt, dbhIb, dib17,
                          sawHt, plpHt, stump, minLen, maxLen, trim)
    # Segment sawtimber section into logs
    lmerch = sawHt - stump
    nlogp  = floor(Int, lmerch / (maxLen + trim))
    # Guard against a degenerate sawHt (no real tree exceeds the 40-slot/340-ft
    # capacity; the Fortran fixed-size LOGLEN array relies on the same). Clamp so
    # the indexing below (incl. the nlogp+1 at line ~960) cannot overflow.
    nlogp  = clamp(nlogp, 0, 39)
    leftov = lmerch - (maxLen + trim)*nlogp - trim

    # Compute log lengths (even-length trimmed); 40-slot array handles trees up to 340 ft sawHt
    logLen = zeros(Float32, 40)
    tlogs = 0

    if !(lmerch < minLen + trim || (nlogp == 0 && leftov < minLen + trim))
        for i in 1:nlogp
            logLen[i] = maxLen
        end
        if leftov >= minLen + trim
            nlogp += 1
            logLen[nlogp] = leftov
        end
        # Round last/second-to-last to even feet
        if nlogp == 1
            logLen[1] = Float32(floor(Int, logLen[1]) ÷ 2 * 2)
        elseif leftov < minLen
            logLen[nlogp] = Float32(floor(Int, logLen[nlogp]) ÷ 2 * 2)
        else
            combined = maxLen + leftov
            logLen[nlogp]   = Float32(floor(Int, combined/2) ÷ 2 * 2)
            logLen[nlogp-1] = Float32((floor(Int, combined - logLen[nlogp]) ÷ 2) * 2)
        end
        tlogs = nlogp
    end

    tlogs == 0 && return 0.0f0

    # vol[2] = Scribner board feet (r9bdft). For SN, METHB=6 → bbfv = vol(2)
    # Scribner (vol(10) International is only used when METHB==9). The log-end
    # DIB is ROUNDED to the nearest inch: INT(DIB+0.499) (r9logdib.f 344), then
    # bf per log = nint(len*scrbnr(iDib)). Height accumulates stump + Σ(trim+len).
    bf_total = 0.0f0
    ht = stump
    for i in 1:tlogs
        len = logLen[i]
        ht_top = ht + trim + len
        dib_top = _r9dib_clark(R,C,E,P,B,A, totHt,dbhIb,dib17, ht_top)
        idib = trunc(Int, dib_top + 0.499f0)   # r9logdib: LOGDIA = INT(DIB+0.499)
        if idib >= 1 && idib <= 120
            bf_total += round(len * Float32(_SCRBNR[idib]), RoundNearestTiesAway)
        end
        ht = ht_top
    end
    return bf_total
end

# R9LOGLEN (r9logs.f:207) — assign even-foot log lengths for one stem segment. Fills `logLen[ilog..]`
# (maxLen logs + a rounded leftover) and returns the last index used. `numseg` is the segment's log
# count (NOLOGP for sawtimber, NOLOGS for topwood); `leftov` the residual merch length past the full logs.
function _r9loglen!(logLen::Vector{Float32}, ilog::Int, jlog::Int, numseg::Int,
                    minLen::Float32, maxLen::Float32, trim::Float32, leftov::Float32)::Int
    if jlog > 0
        for i in ilog:jlog; logLen[i] = maxLen; end
    end
    if leftov >= minLen + trim
        numseg += 1; jlog += 1; logLen[jlog] = leftov
    end
    if numseg == 1                                                   # single log → whole even feet
        logLen[ilog] = Float32(trunc(Int, trunc(Int, logLen[ilog]) / 2) * 2)
    elseif leftov < minLen                                          # short leftover folded into last log
        logLen[jlog] = Float32(trunc(Int, trunc(Int, logLen[jlog]) / 2) * 2)
    else                                                            # resegment top two logs to even feet
        logLen[jlog]   = Float32(trunc(Int, (maxLen + leftov) / 2))
        logLen[jlog]   = Float32(trunc(Int, logLen[jlog] / 2) * 2)
        logLen[jlog-1] = Float32(trunc(Int, maxLen + leftov - logLen[jlog]))
        logLen[jlog-1] = Float32(trunc(Int, logLen[jlog-1] / 2) * 2)
    end
    return jlog
end

# Full-stem per-LOG-DIB GROSS Scribner board feet — ports R9LOGS (sawtimber + topwood segmentation,
# r9logs.f) → R9LOGDIB → R9BDFT (r9clark.f:1433). Returns Dict{idib => Σ gross bf of all logs (saw AND
# pulp, to the pulpwood top) whose rounded end-DIB == idib}. UNLIKE `_r8_scribner_bf` (which sums only
# the sawtimber logs = the net board feet vol[10]), this includes the small topwood logs, so
# `sum(values)` = FVS's ECHARV `treeVol` — the denominator of the log-grade defect proportion
# (defProp = treeVol / netBF). This is what ECHARV (echarv.f) needs to price each log by its DIB grade.
function _r8_scribner_bf_by_dib(R, C, E, P, B, A, totHt, dbhIb, dib17,
                                sawHt, plpHt, stump, minLen, maxLen, trim)::Dict{Int,Float32}
    out = Dict{Int,Float32}()
    sawHt > 0f0 || return out
    logLen = zeros(Float32, 40)
    # --- Sawtimber segmentation (r9logs.f:74-102) ---
    lmerch = sawHt - stump
    nologp = clamp(trunc(Int, lmerch / (maxLen + trim)), 0, 39)
    leftov = lmerch - (maxLen + trim) * nologp - trim
    if !(lmerch < minLen + trim || (nologp == 0 && leftov < minLen + trim))
        nologp = _r9loglen!(logLen, 1, nologp, nologp, minLen, maxLen, trim, leftov)
    else
        nologp = 0
    end
    # --- Topwood (pulpwood) segmentation (r9logs.f:107-145) ---
    nologs = 0
    if plpHt > 0f0 && plpHt > sawHt
        sawTop = stump
        for i in 1:nologp; sawTop += logLen[i] + trim; end
        lmerchp = plpHt - sawTop
        nologs = trunc(Int, lmerchp / (maxLen + trim))
        leftovp = lmerchp - (maxLen + trim) * nologs - trim
        if lmerchp < minLen + trim || (nologs == 0 && leftovp < minLen + trim)
            nologs = 0
        end
        if nologs > 0 || leftovp > minLen + trim
            ilog = nologp + 1
            jlog = ilog + nologs - 1
            if jlog <= 39
                jlog = _r9loglen!(logLen, ilog, jlog, nologs, minLen, maxLen, trim, leftovp)
                nologs = jlog - nologp
            else
                nologs = 0
            end
        end
    end
    tlogs = nologp + nologs
    tlogs == 0 && return out
    # --- R9LOGDIB + R9BDFT: DIB at each log top, round, gross Scribner per log ---
    ht = stump
    for i in 1:tlogs
        len = logLen[i]
        ht_top = ht + trim + len
        ht = ht_top
        len > 0f0 || continue
        dib = Float32(trunc(Int, _r9dib_clark(R,C,E,P,B,A, totHt,dbhIb,dib17, ht_top) + 0.499f0))  # LOGDIA(.,1)
        idib = trunc(Int, dib)                                       # R9BDFT iDib = int(dib)
        if idib >= 1 && idib <= 120
            out[idib] = get(out, idib, 0f0) + round(len * Float32(_SCRBNR[idib]), RoundNearestTiesAway)
        end
    end
    return out
end

# Full-stem per-LOG-DIB GROSS cubic feet — ports R9LOGS (segmentation, shared with the board path) →
# R9LOGDIB (boundary DIBs) → R9LGCFT (r9logs.f:364). Each log's cubic is the Smalian estimate
# 0.00272708·(Dbot²+Dtop²)·len over the predicted (index-2) boundary DIBs, then ALL logs are renormalized
# so Σ = `cfvol` (= VOL(4)+VOL(7), the gross merch cubic; r9clark.f:428). Returns Dict{idib => Σ gross
# cuft of logs whose rounded top scaling-DIB == idib}; `sum(values)` = ECHARV `treeVol` (= cfvol), the
# defProp denominator's numerator. Mirrors `_r8_scribner_bf_by_dib` exactly except for the per-log measure.
function _r8_cuft_by_dib(R, C, E, P, B, A, totHt, dbhIb, dib17,
                         sawHt, plpHt, stump, minLen, maxLen, trim, cfvol)::Dict{Int,Float32}
    out = Dict{Int,Float32}()
    (sawHt > 0f0 && cfvol > 0f0) || return out
    logLen = zeros(Float32, 40)
    # --- segmentation identical to _r8_scribner_bf_by_dib (R9LOGS sawtimber + topwood) ---
    lmerch = sawHt - stump
    nologp = clamp(trunc(Int, lmerch / (maxLen + trim)), 0, 39)
    leftov = lmerch - (maxLen + trim) * nologp - trim
    if !(lmerch < minLen + trim || (nologp == 0 && leftov < minLen + trim))
        nologp = _r9loglen!(logLen, 1, nologp, nologp, minLen, maxLen, trim, leftov)
    else
        nologp = 0
    end
    nologs = 0
    if plpHt > 0f0 && plpHt > sawHt
        sawTop = stump
        for i in 1:nologp; sawTop += logLen[i] + trim; end
        lmerchp = plpHt - sawTop
        nologs = trunc(Int, lmerchp / (maxLen + trim))
        leftovp = lmerchp - (maxLen + trim) * nologs - trim
        if lmerchp < minLen + trim || (nologs == 0 && leftovp < minLen + trim)
            nologs = 0
        end
        if nologs > 0 || leftovp > minLen + trim
            ilog = nologp + 1
            jlog = ilog + nologs - 1
            if jlog <= 39
                jlog = _r9loglen!(logLen, ilog, jlog, nologs, minLen, maxLen, trim, leftovp)
                nologs = jlog - nologp
            else
                nologs = 0
            end
        end
    end
    tlogs = nologp + nologs
    tlogs == 0 && return out
    # --- R9LOGDIB boundary DIBs (index-2 predicted) + R9LGCFT Smalian, then renormalize to cfvol ---
    # boundary 1's DIB is at 4.5 ft (R9LOGDIB:24-26), NOT the stump; boundary I+1 at the cumulative log top.
    dib_bot = _r9dib_clark(R,C,E,P,B,A, totHt,dbhIb,dib17, 4.5f0)            # LOGDIA(1,2)
    sm = zeros(Float32, tlogs); idibs = zeros(Int, tlogs)
    tlogvol = 0f0
    ht = stump
    for i in 1:tlogs
        len = logLen[i]
        ht += trim + len                                                    # R9LOGDIB:33 HT=HT+TRIM+LOGLEN(I)
        dib_top = _r9dib_clark(R,C,E,P,B,A, totHt,dbhIb,dib17, ht)           # LOGDIA(I+1,2)
        if len > 0f0
            sm[i] = 0.00272708f0 * (dib_bot^2 + dib_top^2) * len            # R9LGCFT:392
            tlogvol += sm[i]
            idibs[i] = trunc(Int, dib_top + 0.499f0)                        # LOGDIA(I+1,1) = grading scaling-DIB
        end
        dib_bot = dib_top
    end
    tlogvol > 0f0 || return out
    for i in 1:tlogs
        sm[i] > 0f0 || continue
        idib = idibs[i]
        (idib >= 1 && idib <= 120) || continue
        out[idib] = get(out, idib, 0f0) + sm[i] / tlogvol * cfvol           # R9LGCFT:407 renormalize to cfvol
    end
    return out
end

# ---------------------------------------------------------------------------
# r9dib: inside-bark diameter at height stemHt.
# Translated from r9clark.f lines 1116-1241.
# ---------------------------------------------------------------------------
# FVS R9DIB (r9clark.f:1155-1226) — the Clark stem inside-bark diameter at a height. Computed ENTIRELY
# in Float32 (`REAL*4`), matching Fortran's op sequence: the real-exponent Clark powers (`**R`, `**P`,
# `**0.5`) use Float32 `^`, the integer exponents (`**2`,`**3`) are exact Float32 multiplications, and the
# final `(Ds+Db+Dt)**0.5` is `^0.5f0` (NOT `sqrt`) to match FVS's literal `**0.5`. Computing in Float64
# and rounding once at the end (the earlier form) makes the DIB *more* precise than FVS, which tips the
# `INT(DIB+0.499)` Scribner-bucket at log tops whose raw DIB sits ≈x.50 (the .499 knife-edge) → a whole
# Scribner-row board-foot step (~16 bf on a large tree). Faithful Float32 keeps the bucket bit-exact.
function _r9dib_clark(R::Real, C::Real, E::Real, P::Real, B::Real, A::Real,
                       totHt::Real, dbhIb::Real, dib17::Real, stemHt::Real)
    stemHt = Float32(stemHt)
    stemHt <= 0f0 && return 0.0f0
    R = Float32(R); C = Float32(C); E = Float32(E); P = Float32(P); B = Float32(B); A = Float32(A)
    totHt = Float32(totHt); dbhIb = Float32(dbhIb); dib17 = Float32(dib17)

    # Trap near-tip math errors
    if R < 0f0 && abs(stemHt - totHt) < 0.00001f0
        stemHt = stemHt - 0.1f0
    end

    Is = stemHt < 4.5f0
    Ib = stemHt >= 4.5f0 && stemHt <= 17.3f0
    It = stemHt > 17.3f0
    Im = stemHt < (17.3f0 + A*(totHt - 17.3f0))

    StTot = stemHt/totHt
    if log(max(1f0 - StTot, 1f-20)) < (-20f0/abs(R))
        StTot = 1.0f0
    end

    Ds, Db, Dt = 0.0f0, 0.0f0, 0.0f0

    if Is
        Ds = dbhIb^2 * (1f0 + (C + E/dbhIb^3)*
             (fpow(1f0-StTot, R) - fpow(1f0-4.5f0/totHt, R)) /
             (1f0 - fpow(1f0-4.5f0/totHt, R)))
    end
    if Ib
        Db = dbhIb^2 - (dbhIb^2 - dib17^2)*
             (fpow(1f0-4.5f0/totHt, P) - fpow(1f0-stemHt/totHt, P)) /
             (fpow(1f0-4.5f0/totHt, P) - fpow(1f0-17.3f0/totHt, P))
    end
    if It
        Dt = dib17^2*(B*(((stemHt-17.3f0)/(totHt-17.3f0))-1f0)^2
             + (Im ? ((1f0-B)/A^2)*(A-(stemHt-17.3f0)/(totHt-17.3f0))^2 : 0.0f0))
    end

    val = Ds + Db + Dt
    val > 0f0 || return 0.0f0
    return fpow(val, 0.5f0)
end
