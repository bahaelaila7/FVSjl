# =============================================================================
# r9clark_vol.jl — NVEL Region-9 Clark profile volume model (NE variant)
#
# Translated from volume/r9clark_fvsMod.f (r9Prep / r9dia417 / r9totHt / r9cuft /
# r9ht / r9cor) + ne/mrules.f (REGN.EQ.9) + the fvsvol.f NATCRS driver. NE sets
# METHB=METHC=6 (ne/grinit.f) → every species uses the eastern Clark profile model
# `900CLKE{fia}` (geog 0, topDib digit 0). The taper math here is the Region-9
# sibling of the already-ported R8 Clark (r8clark_vol.jl) — same profile family,
# different coefficient fit (volume/NVEL/r9coeff.inc, 47 species/groups).
#
# Entry: `r9clark_cubic(spp, dbhOb, htTot, prod, mTopP, mTopS, stump)` → vol[15]
# (vol[1]=total cuft, vol[4]=saw cuft, vol[7]=topwood cuft). Board feet (vol[2])
# needs the R9LOGS/r9bdft Scribner path (separate follow-up).
#
# ── VALIDATION STATE (WIP — NOT yet wired into compute_volumes!) ────────────
# Validated PER-TREE against the live FVSne_new cycle-0 COMPLETE TREE LIST (.trl),
# computing on the oracle's OWN (dbh,ht,species) — apples to apples:
#   per-tree CLOSE: SM d=18.1 h=91.9: 65.8 vs 66.0; JP d=11.5: 20.8 vs 20.8;
#   max per-tree total-cuft diff 1.65 (after the d<1 guard).
# Aggregated over the oracle's cyc0 trees (Σ vol·tpa, with the d<1→0 guard that
# compute_volumes! already applies):
#   TOT 1546.5 / 1558.8    MCH 1338.0 / 1346.7    SAW 294.1 / 292.5   — ALL <1%.
# (The earlier "6% / 26%" gaps were ARTIFACTS — a whitespace-split column
# misalignment in the diff harness, and the missing small-tree guard: one YB d=0.1
# seedling gave 2.5 cuft × 27 tpa ≈ the whole apparent total-cuft gap.) The kernels
# (r9dia417/r9totHt/r9cuft/r9ht) + r9Prep group-fallback + MRULES R9 + r9cor are
# faithful. Remaining <1% is per-tree Float32 rounding (Fortran nint vs round) + a
# few medium-tree residuals (e.g. SM d=10.4: 13.8 vs 15.4) + the exact small-tree
# cutoff — refinements, not structural. NEXT: wire into compute_volumes! (variant-
# dispatch + d<1 guard), close the <1% residual, then add board feet (vol[2]).
#
# PER-TREE ORACLE RECIPE (this unblocked the diagnosis): build a clean SINGLE-stand
# keyfile (CR→LF; take lines up to the first PROCESS), insert a `DATABASE/DSNOUT/
# <db>/SUMMARY/TREELIDB/END` block BEFORE `TREEDATA`, and NAME THE TREE FILE TO MATCH
# THE KEYFILE STEM (TREEDATA reads <stem>.tre). Run /tmp/FVSne_new → the `.trl`
# COMPLETE TREE LIST has per-tree TOT CU / MCH CU / SAW CU / SAW BD columns.
#
# Board feet (vol[2]) is also unported. NOT bit-exact → standalone module, not wired.
# =============================================================================

# --- R9 Clark coefficients (volume/NVEL/r9coeff.inc → CSV) -------------------
# Per species/group row: dib4in,a4,b4 (coefA) + a17/b17 at three top-diam classes
# (coef0=topDib 0, coef4=topDib 4, coef79=topDib 7/9) + the total-height profile
# params r,c,e,p,a,b (from coef0 cols 4-9). Species codes ≥1000 are species groups.
struct _R9Coef
    spp::Int
    a4::Float32; b4::Float32
    a17_0::Float32; b17_0::Float32
    a17_4::Float32; b17_4::Float32
    a17_7::Float32; b17_7::Float32
    r::Float32; c::Float32; e::Float32; p::Float32; a::Float32; b::Float32
end

const _R9COEF = let
    path = joinpath(@__DIR__, "..", "..", "data", "northeast", "volume", "r9clark_coef.csv")
    rows = _R9Coef[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        f = split(strip(line), ',')
        isempty(f[1]) && continue
        g(i) = parse(Float32, f[i])
        # cols: spp,dib4in,a4,b4, a17_0,b17_0,r_0,c_0,e_0,p_0,a_0,b_0,
        #       a17_4,b17_4,r_4,c_4,e_4,p_4,q_4, a17_7,b17_7,r_7,c_7,e_7,p_7,q_7
        push!(rows, _R9Coef(parse(Int, f[1]),
            g(3), g(4),            # a4,b4  (coefA cols 3,4; col2 dib4in unused here)
            g(5), g(6),            # a17_0,b17_0
            g(13), g(14),          # a17_4,b17_4
            g(20), g(21),          # a17_7,b17_7
            g(7), g(8), g(9), g(10), g(11), g(12)))  # r_0,c_0,e_0,p_0,a_0,b_0
    end
    sort!(rows; by = r -> r.spp)
    rows
end
const _R9_GRP_START = findfirst(r -> r.spp >= 1000, _R9COEF)  # grpIdx (first group row)

# r9Prep species→group fallback (r9clark_fvsMod.f:592-651): exact species first,
# else a group code (conifer/spruce/pine + the hardwood genus groups).
function _r9_spp_group(spp::Int)
    if spp < 300
        spp in 90:99   && return 1090   # spruces
        spp in 100:199 && return 1100   # pines
        return 1000                     # other conifers
    end
    spp in 310:329 && return 1310       # maples
    spp in 370:379 && return 1370       # birches
    spp in 400:410 && return 1400       # hickories
    spp in 540:549 && return 1540       # ashes
    (spp in (740, 742, 744, 745, 753))      && return 1740   # cottonwoods
    (spp in (741, 743, 746, 752))           && return 1750   # poplars
    spp in 760:769 && return 1760       # cherries/plums
    spp in 800:899 && return 1800       # oaks
    spp in 950:954 && return 1950       # basswoods
    spp in 970:979 && return 1970       # (oaks/elms group)
    return 1300                         # generic hardwood
end

"Look up the R9 Clark coefficient row for FIA species `spp`, applying the group fallback."
function _r9_coef(spp::Int)::_R9Coef
    idx = findfirst(r -> r.spp == spp, _R9COEF)
    idx !== nothing && return _R9COEF[idx]
    grp = _r9_spp_group(spp)
    gidx = findfirst(r -> r.spp == grp, _R9COEF)
    gidx === nothing && error("R9 Clark: no coefficient for species $spp (group $grp)")
    return _R9COEF[gidx]
end

# --- ne/mrules.f REGN.EQ.9 (CLK / NVB) --------------------------------------
# Returns (maxLen, minLen, minLenT, merchL, trim, mTopP, mTopS, stump) after the
# region-9 defaults + the prod-01 stump rule. mTopP/mTopS/stump default only when
# the caller passes ≤0 (the driver passes SCFTOPD/TOPD/SCFSTMP, all set, so the
# fallbacks here only fire for an unset stump).
@inline function _r9_mrules(spp::Int, prod::String, mTopP::Float32, mTopS::Float32, stump::Float32)
    maxLen = 8f0; minLen = 2f0; minLenT = 4f0; merchL = 8f0; trim = 0.3f0
    if mTopP <= 0f0
        mTopP = spp < 300 ? 7.6f0 : 9.6f0
    end
    mTopS <= 0f0 && (mTopS = 4f0)
    if stump <= 0.01f0
        stump = prod == "01" ? 1f0 : 0.5f0
    end
    return maxLen, minLen, minLenT, merchL, trim, mTopP, mTopS, stump
end

# --- mutable per-tree coefficient state (the Fortran CLKCOEF type) -----------
mutable struct _R9State
    r::Float32; c::Float32; e::Float32; p::Float32; a::Float32; b::Float32
    a4::Float32; b4::Float32; a17::Float32; b17::Float32
    dbhIb::Float32; dib17::Float32; totHt::Float32
end

# r9totHt (r9clark_fvsMod.f:920-958): total height from inside-bark dbh, top
# height/diameter and the a/b total-height coefficients.
function _r9_totht(htTot::Float32, dbhIb::Float32, dib17::Float32,
                   topHt::Float32, topDib::Float32, a::Float32, b::Float32)::Float32
    if htTot > 0f0
        return htTot > 17.4f0 ? htTot : topHt
    elseif topHt > 17.3f0
        Im = topDib^2 > b * (a - 1f0)^2 * dib17^2 ? 1f0 : 0f0
        Qa = b + Im * (1f0 - b) / a^2
        Qb = -2f0 * b - Im * 2f0 * (1f0 - b) / a
        Qc = b + (1f0 - b) * Im - topDib^2 / dib17^2
        tot = 17.3f0 + (topHt - 17.3f0) * (2f0 * Qa) /
              (-1f0 * Qb - sqrt(Qb^2 - 4f0 * Qa * Qc))
        tot = max(tot, topHt + topDib * 2f0)
        tot = min(tot, topHt + topDib * 8f0)
        return tot
    else
        return 17.3f0 + dib17 * 3f0
    end
end

# r9dia417 (r9clark_fvsMod.f:817-917): inside-bark diameter at 4.5' and 17.3'.
# Returns (dbhIb, dib17, errFlg). net01 trees carry measured htTot so topDib=0.
function _r9_dia417(st::_R9State, topDib::Float32, dbhOb::Float32, topHt::Float32,
                    ht1Prd::Float32, ht2Prd::Float32, sawDib::Float32, plpDib::Float32)
    dbhIb = st.a4 + st.b4 * dbhOb
    (dbhIb >= dbhOb || dbhIb <= 0f0) && (dbhIb = max(dbhOb - 0.1f0, 0.1f0))
    if (topDib > dbhIb && topHt > 4.5f0) || (topDib < dbhIb && topHt < 4.5f0)
        return dbhIb, 0f0, 11
    end
    local dib17::Float32
    if abs(ht2Prd - 17.3f0) < 1f-5
        dib17 = plpDib
    elseif abs(ht1Prd - 17.3f0) < 1f-5
        dib17 = sawDib
    elseif topHt > 17.3f0
        dib17 = dbhIb * (st.a17 + st.b17 * (17.3f0 / topHt)^2)   # R9 uses dbhIb
        dib17 = max(dib17, topDib + 0.1f0)
    else
        dib17 = topDib - 0.1f0
    end
    if ht1Prd > 17.3f0 && dib17 < sawDib
        dib17 = sawDib + (dbhIb - sawDib) * (ht1Prd - 17.3f0) / (ht1Prd - 4.5f0)
    elseif ht2Prd > 17.3f0 && dib17 < plpDib
        dib17 = plpDib + (dbhIb - plpDib) * (ht2Prd - 17.3f0) / (ht2Prd - 4.5f0)
    end
    dib17 < 0.1f0 && (dib17 = 0.1f0)
    return dbhIb, dib17, 0
end

# r9cuft (r9clark_fvsMod.f:961-1113): cubic foot volume from lowrHt to upprHt
# along the 3-segment Clark profile.
function _r9_cuft(st::_R9State, lowrHt::Float32, upprHt::Float32)::Float32
    upprHt <= 0f0 && return 0f0
    r = st.r; c = st.c; e = st.e; p = st.p; b = st.b; a = st.a
    totht = st.totHt; dbhib = st.dbhIb; dib17 = st.dib17
    G = (1f0 - 4.5f0 / totht)^r
    W = (c + e / dbhib^3) / (1f0 - G)
    X = (1f0 - 4.5f0 / totht)^p
    Y = ((1f0 - 17.3f0 / totht) < 0.005748f0 && p > 14f0) ? 0f0 : (1f0 - 17.3f0 / totht)^p
    Z = (dbhib^2 - dib17^2) / (X - Y)
    T = dbhib^2 - Z * X
    L1 = max(lowrHt, 0f0); U1 = min(upprHt, 4.5f0)
    L2 = max(lowrHt, 4.5f0); U2 = min(upprHt, 17.3f0)
    L3 = max(lowrHt, 17.3f0); U3 = min(totht, upprHt)
    I1 = lowrHt < 4.5f0 ? 1f0 : 0f0
    I2 = lowrHt < 17.3f0 ? 1f0 : 0f0
    I3 = upprHt > 4.5f0 ? 1f0 : 0f0
    I4 = upprHt > 17.3f0 ? 1f0 : 0f0
    I5 = (L3 - 17.3f0) < a * (totht - 17.3f0) ? 1f0 : 0f0
    I6 = (U3 - 17.3f0) < a * (totht - 17.3f0) ? 1f0 : 0f0
    V1 = 0f0; V2 = 0f0; V3 = 0f0
    if I1 > 0f0
        V1 = I1 * dbhib^2 * ((1f0 - G * W) * (U1 - L1) +
             W * ((1f0 - L1 / totht)^r * (totht - L1) -
                  (1f0 - U1 / totht)^r * (totht - U1)) / (r + 1f0))
    end
    if I2 > 0f0 && I3 > 0f0
        if (1f0 - U2 / totht) < 0.005748f0 && p > 14f0
            V2 = T * (U2 - L2) + Z * ((1f0 - L2 / totht)^p * (totht - L2)) / (p + 1f0)
        else
            V2 = T * (U2 - L2) + Z * ((1f0 - L2 / totht)^p * (totht - L2) -
                 (1f0 - U2 / totht)^p * (totht - U2)) / (p + 1f0)
        end
    end
    if I4 > 0f0
        V3 = dib17^2 * (b * (U3 - L3) - b * ((U3 - 17.3f0)^2 - (L3 - 17.3f0)^2) / (totht - 17.3f0) +
             (b / 3f0) * ((U3 - 17.3f0)^3 - (L3 - 17.3f0)^3) / (totht - 17.3f0)^2 +
             I5 * (1f0 / 3f0) * ((1f0 - b) / a^2) * (a * (totht - 17.3f0) - (L3 - 17.3f0))^3 / (totht - 17.3f0)^2 -
             I6 * (1f0 / 3f0) * ((1f0 - b) / a^2) * (a * (totht - 17.3f0) - (U3 - 17.3f0))^3 / (totht - 17.3f0)^2)
    end
    cfVol = 0.005454154f0 * (V1 + V2 + V3)
    return cfVol < 0f0 ? 0f0 : cfVol
end

# r9ht (r9clark_fvsMod.f:1245-1369): height at which inside-bark diameter stmDib occurs.
function _r9_ht(st::_R9State, stmDib::Float32)::Float32
    totHt = st.totHt; dbhIb = st.dbhIb; dib17 = st.dib17
    r = st.r; c = st.c; e = st.e; p = st.p; b = st.b; a = st.a
    G = (1f0 - 4.5f0 / totHt)^r
    W = (c + e / dbhIb^3) / (1f0 - G)
    X = (1f0 - 4.5f0 / totHt)^p
    Y = ((1f0 - 17.3f0 / totHt) < 0.005748f0 && p > 14f0) ? 0f0 : (1f0 - 17.3f0 / totHt)^p
    Z = (dbhIb^2 - dib17^2) / (X - Y)
    Is = stmDib >= dbhIb ? 1f0 : 0f0
    Ib = (stmDib < dbhIb && stmDib >= dib17) ? 1f0 : 0f0
    Im = stmDib^2 > b * (a - 1f0)^2 * dib17^2 ? 1f0 : 0f0
    Qa = b + Im * (1f0 - b) / a^2
    Qb = -2f0 * b - Im * 2f0 * (1f0 - b) / a
    Qc = b + (1f0 - b) * Im - stmDib^2 / dib17^2
    stemHt = 0f0
    if Is == 1f0
        xxx = (stmDib^2 / dbhIb^2 - 1f0) / W + G
        xxx > 0f0 && (stemHt = totHt * (1f0 - xxx^(1f0 / r)))
    elseif Ib == 1f0
        xxx = X - (dbhIb^2 - stmDib^2) / Z
        xxx > 0f0 && (stemHt = totHt * (1f0 - xxx^(1f0 / p)))
    else
        xxx = Qb^2 - 4f0 * Qa * Qc
        xxx > 0f0 && (stemHt = 17.3f0 + (totHt - 17.3f0) * ((-Qb - sqrt(xxx)) / (2f0 * Qa)))
    end
    return stemHt
end

# r9cor (r9clark_fvsMod.f:1507-1574): final per-species/product correction factors.
@inline function _r9_cor!(vol::Vector{Float32}, spp::Int, iProd::Int)
    if spp < 300
        cf1 = cf2 = cf3 = 1.04f0
    elseif (741 <= spp <= 746) || spp == 621
        cf1 = cf2 = cf3 = 1f0
    else
        cf1 = cf2 = cf3 = 1.1f0
    end
    vol[1] *= cf2
    vol[2] *= cf3
    vol[7] *= cf2
    vol[4] *= iProd == 1 ? cf1 : cf2
    return vol
end

# r9dib (r9clark.f:1116): inside-bark diameter at height `h` along the 3-segment Clark profile.
function _r9_dib(st::_R9State, h::Float32)::Float32
    r = st.r; c = st.c; e = st.e; p = st.p; b = st.b; a = st.a
    totHt = st.totHt; dbhIb = st.dbhIb; dib17 = st.dib17
    (r < 0f0 && abs(h - totHt) < 0.00001f0) && (h = h - 0.1f0)
    sttot = h / totHt
    log(1f0 - sttot) < (-20f0 / abs(r)) && (sttot = 1f0)
    ds = 0f0; db = 0f0; dt = 0f0
    if h < 4.5f0
        ds = dbhIb^2 * (1f0 + (c + e / dbhIb^3) *
             ((1f0 - sttot)^r - (1f0 - 4.5f0 / totHt)^r) / (1f0 - (1f0 - 4.5f0 / totHt)^r))
    elseif h <= 17.3f0
        db = dbhIb^2 - (dbhIb^2 - dib17^2) *
             ((1f0 - 4.5f0 / totHt)^p - (1f0 - h / totHt)^p) /
             ((1f0 - 4.5f0 / totHt)^p - (1f0 - 17.3f0 / totHt)^p)
    else
        im = h < (17.3f0 + a * (totHt - 17.3f0))
        dt = dib17^2 * (b * (((h - 17.3f0) / (totHt - 17.3f0)) - 1f0)^2 +
             (im ? ((1f0 - b) / a^2) * (a - (h - 17.3f0) / (totHt - 17.3f0))^2 : 0f0))
    end
    s = ds + db + dt
    return s > 0f0 ? sqrt(s) : 0f0
end

# International ¼-inch board feet per log (r9bdft, r9clark.f:1482) — the NE `.sum` BdFt is FVS's
# vol(10) (International), NOT vol(2) (Scribner). Log-end (small) DIB is the rounded INT(dib+0.499).
@inline function _r9_intl_log(len::Float32, idib::Int)::Float32
    idib < 4 && return 0f0
    d = Float32(idib)
    bd = 0.04976191f0 * len * d^2 + 0.006220239f0 * len^2 * d - 0.1854762f0 * len * d +
         0.0002591767f0 * len^3 + 0.01159226f0 * len^2 + 0.04222222f0 * len
    return round(bd / 5f0, RoundNearestTiesAway) * 5f0
end

# r9bdft (r9clark.f:1380) + R9LOGS bucking (r9logs.f R9LOGLEN/R9LOGDIB): board feet of the sawtimber
# section [stump, sawHt]. Bucks even-foot logs (shared R9LOGLEN rule, identical to `_r8_scribner_bf`),
# takes each log's small-end (top) DIB from the R9 taper `_r9_dib`. Returns vol(10) = International ¼"
# board feet (`nint` of the per-log sum), the value the NE `.sum` BdFt column reports.
function _r9_intlqtr_bf(st::_R9State, sawHt::Float32, stump::Float32,
                        minLen::Float32, maxLen::Float32, trim::Float32)::Float32
    lmerch = sawHt - stump
    nlogp = clamp(floor(Int, lmerch / (maxLen + trim)), 0, 39)
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
        ht += trim + len                        # top (small end) of log i
        idib = trunc(Int, _r9_dib(st, ht) + 0.499f0)   # r9logdib: LOGDIA = INT(DIB+0.499)
        bf += _r9_intl_log(len, idib)
    end
    return round(bf)                            # r9bdft:1499 vol(10)=NINT(vol(10))
end

"""
    r9clark_cubic(spp, dbhOb, htTot, prod, mTopP, mTopS, stump) -> vol::Vector{Float32}

R9 Clark cubic volumes for FIA species `spp`. `prod`="01" (sawtimber) or "02"
(pulpwood). Returns vol[15]: vol[1]=total cuft, vol[4]=sawtimber cuft, vol[7]=
topwood cuft (so merch cuft = vol[4]+vol[7]). Mirrors r9clark's prod-01 cubic
branch + the NATCRS flags (CUTFLG=CUPFLG=SPFLG=1); board feet (vol[2]) is not yet
computed here (needs R9LOGS/r9bdft). NaN/short-tree guards follow the Fortran.
"""
function r9clark_cubic(spp::Int, dbhOb::Float32, htTot::Float32, prod::String,
                       mTopP::Float32, mTopS::Float32, stump::Float32)::Vector{Float32}
    vol = zeros(Float32, 15)
    iProd = prod == "01" ? 1 : 2
    co = _r9_coef(spp)
    maxLen, minLen, minLenT, merchL, trim, mTopP, mTopS, stump =
        _r9_mrules(spp, prod, mTopP, mTopS, stump)
    sawDib = mTopP != 0f0 ? mTopP : (spp < 300 ? 7.6f0 : 9.6f0)
    plpDib = mTopS != 0f0 ? mTopS : 4f0

    # topHt / topDib from measured total height (r9Prep:700-753). net01 trees all
    # carry htTot, so this is the htTot>0 branch.
    short = false; shrtHt = 0f0; topDib = 0f0; topHt = 0f0
    if htTot > 0f0
        topDib = 0f0
        if htTot >= 17.4f0
            topHt = htTot
        else
            short = true; topHt = 17.4f0; shrtHt = htTot
        end
    else
        return vol      # height required for net01 path
    end
    dbhOb <= topDib && return vol

    st = _R9State(co.r, co.c, co.e, co.p, co.a, co.b, co.a4, co.b4,
                  co.a17_0, co.b17_0, 0f0, 0f0, 0f0)
    # A17/B17 select on topDib (=0 here → coef0 class).
    st.a17 = co.a17_0; st.b17 = co.b17_0

    dbhIb, dib17, err = _r9_dia417(st, topDib, dbhOb, topHt, 0f0, 0f0, sawDib, plpDib)
    err != 0 && return vol
    st.dbhIb = dbhIb; st.dib17 = dib17
    st.totHt = _r9_totht(htTot, dbhIb, dib17, topHt, topDib, st.a, st.b)
    st.totHt <= 17.3f0 && return vol

    # total cuft to the tip (CUTFLG=1)
    cfVol = _r9_cuft(st, stump, st.totHt)
    short && (cfVol *= shrtHt / 17.3f0)
    vol[1] = cfVol

    # secondary-product flag: SPFLG=1 when D≥SCFMIND, i.e. the sawtimber (prod-01) path.
    spFlg = iProd == 1

    # Pulpwood / merch-cubic section (r9clark:237-280): runs for prod≠1 OR spFlg.
    # tcfVol = cubic stump→pulp-top. For prod-02 trees this IS the merch cubic (vol[4]);
    # for prod-01 trees it is carried so topwood = tcfVol − saw cubic.
    tcfVol = 0f0
    if iProd != 1 || spFlg
        plpHt = _r9_ht(st, plpDib)
        (topDib <= plpDib && topHt < plpHt) && (plpHt = topHt)
        plpHt < merchL + stump + trim && (plpHt = 0f0)
        if plpHt - stump >= minLen
            tcfVol = _r9_cuft(st, stump, plpHt)
            short && (tcfVol *= shrtHt / 17.3f0)
            iProd != 1 && (vol[4] = tcfVol)            # prod-02 merch cubic (CUPFLG=1)
        end
    end

    # Sawtimber section (r9clark:286-367): iProd==1 → saw cubic + topwood.
    if iProd == 1
        sawHt = _r9_ht(st, sawDib)
        (topDib <= sawDib && topHt < sawHt) && (sawHt = topHt)
        sawHt < merchL + trim + stump && (sawHt = 0f0)
        if sawHt > 0f0
            scfVol = _r9_cuft(st, stump, sawHt)
            short && (scfVol *= shrtHt / 17.3f0)
            vol[4] = scfVol
            spFlg && (vol[7] = max(tcfVol - scfVol, 0f0))     # topwood
            # Board feet of the saw section (r9bdft International ¼"); cf correction by _r9_cor!
            # (vol(10)*=cf4, and cf4≡cf3 — so storing it in vol[2] gets the right factor).
            vol[2] = _r9_intlqtr_bf(st, sawHt, stump, minLen, maxLen, trim)
        end
    end

    # rounding (r9clark:456-461) then r9cor corrections
    vol[1] > 0f0 && (vol[1] = round(vol[1] * 10f0) / 10f0)
    vol[4] > 0f0 && (vol[4] = round(vol[4] * 10f0) / 10f0)
    vol[7] > 0f0 && (vol[7] = round(vol[7] * 10f0) / 10f0)
    _r9_cor!(vol, spp, iProd)
    return vol
end

# --- NE merch standards (ne/sitset.f:505-560) + the per-tree volume driver -----
# IFOR-dependent merch rules. Softwood = species index ≤ 25. Stumps from grinit.f
# (STMP 0.5 / SCFSTMP 1.0). Returns the per-tree (dbhmin, topd, scfmind, scftopd,
# stmp, scfstmp) used by the cubic columns. Board-foot mins are bf-equal here.
@inline function _ne_merch(spi::Integer, ifor::Integer)
    if spi <= 25                          # softwoods
        return (5f0, 4f0, 9f0, 7.6f0, 0.5f0, 1f0)
    else                                  # hardwoods
        dbhmin = (ifor == 1 || ifor == 3) ? 6f0 : (ifor == 4 ? 8f0 : 5f0)
        topd   = ifor == 3 ? 5f0 : 4f0
        return (dbhmin, topd, 11f0, 9.6f0, 0.5f0, 1f0)
    end
end

"""
    compute_volumes_ne!(s)

NE per-tree volume driver (the Region-9 analogue of the SN `compute_volumes!`
body): NVEL R9 Clark cubic via `r9clark_cubic`, with the NE IFOR-dependent merch
standards. Loads cuft/merch-cuft/saw-cuft/board-feet per tree (board feet = Scribner
via R9LOGS bucking + r9bdft). Broken-top (CFTOPK) reuse is TODO.
"""
function compute_volumes_ne!(s::StandState)
    t = s.trees; co = s.coef
    ifor = Int(s.plot.forest_idx); ifor == 0 && (ifor = _NE_DEFAULT_IFOR)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
            continue
        end
        fias = strip(string(co.code_fia[sp]))
        fia = isempty(fias) ? 0 : parse(Int, fias)
        dbhmin, topd, scfmind, scftopd, _stmp, _scfstmp = _ne_merch(sp, ifor)
        prod = d >= scfmind ? "01" : "02"
        mtopp = d >= scfmind ? scftopd : topd
        v = r9clark_cubic(fia, d, h, prod, mtopp, topd, 0f0)
        t.cuft_vol[i]      = v[1]
        t.merch_cuft_vol[i] = d >= dbhmin  ? v[4] + v[7] : 0f0
        t.saw_cuft_vol[i]   = d >= scfmind ? v[4] : 0f0
        t.bdft_vol[i]       = d >= scfmind ? v[2] : 0f0   # Scribner board feet (R9LOGS/r9bdft)
    end
    return s
end
