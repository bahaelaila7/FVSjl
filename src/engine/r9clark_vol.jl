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

# --- per-tree coefficient state (the Fortran CLKCOEF type) -------------------
# IMMUTABLE + isbits (all Float32) ⇒ stack-allocated, no per-tree heap. r9clark_cubic rebuilds it once
# with the resolved dbhIb/dib17/totHt (the only fields set after construction); the taper functions
# (_r9_dia417/_r9_cuft/_r9_ht/_r9_dib/_r9_bucked_bf) only READ it.
struct _R9State
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
    G = fpow(1f0 - 4.5f0 / totht, r)                       # Clark real powers via gfortran companion (doctrine #8)
    W = (c + e / dbhib^3) / (1f0 - G)
    X = fpow(1f0 - 4.5f0 / totht, p)
    Y = ((1f0 - 17.3f0 / totht) < 0.005748f0 && p > 14f0) ? 0f0 : fpow(1f0 - 17.3f0 / totht, p)
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
             W * (fpow(1f0 - L1 / totht, r) * (totht - L1) -
                  fpow(1f0 - U1 / totht, r) * (totht - U1)) / (r + 1f0))
    end
    if I2 > 0f0 && I3 > 0f0
        if (1f0 - U2 / totht) < 0.005748f0 && p > 14f0
            V2 = T * (U2 - L2) + Z * (fpow(1f0 - L2 / totht, p) * (totht - L2)) / (p + 1f0)
        else
            V2 = T * (U2 - L2) + Z * (fpow(1f0 - L2 / totht, p) * (totht - L2) -
                 fpow(1f0 - U2 / totht, p) * (totht - U2)) / (p + 1f0)
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
    G = fpow(1f0 - 4.5f0 / totHt, r)                       # Clark real powers **R/**P → FFI companion (doctrine #8)
    W = (c + e / dbhIb^3) / (1f0 - G)
    X = fpow(1f0 - 4.5f0 / totHt, p)
    Y = ((1f0 - 17.3f0 / totHt) < 0.005748f0 && p > 14f0) ? 0f0 : fpow(1f0 - 17.3f0 / totHt, p)
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
        xxx > 0f0 && (stemHt = totHt * (1f0 - fpow(xxx, 1f0 / r)))
    elseif Ib == 1f0
        xxx = X - (dbhIb^2 - stmDib^2) / Z
        xxx > 0f0 && (stemHt = totHt * (1f0 - fpow(xxx, 1f0 / p)))
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
    flog(1f0 - sttot) < (-20f0 / abs(r)) && (sttot = 1f0)   # ALOG + **R/**P → FFI companion (doctrine #8)
    ds = 0f0; db = 0f0; dt = 0f0
    if h < 4.5f0
        ds = dbhIb^2 * (1f0 + (c + e / dbhIb^3) *
             (fpow(1f0 - sttot, r) - fpow(1f0 - 4.5f0 / totHt, r)) / (1f0 - fpow(1f0 - 4.5f0 / totHt, r)))
    elseif h <= 17.3f0
        # r9clark.f r9dib: guard the near-zero denominator term for a very short bole with a
        # steep taper — Y=0 when (1−17.3/totHt) < 0.005748 AND p > 14 (else the usual (…)^p).
        y = ((1f0 - 17.3f0 / totHt) < 0.005748f0 && p > 14f0) ? 0f0 : fpow(1f0 - 17.3f0 / totHt, p)
        db = dbhIb^2 - (dbhIb^2 - dib17^2) *
             (fpow(1f0 - 4.5f0 / totHt, p) - fpow(1f0 - h / totHt, p)) /
             (fpow(1f0 - 4.5f0 / totHt, p) - y)
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

# Scribner Decimal-C board feet per log (r9clark.f R9BDFT `nint(len·scrbnr(int(dib)))`, scrbnr = the
# 120-entry FACTOR table in data/lakestates/scribner_factor.csv). The LS `.sum` BdFt is FVS's vol(2)
# (Scribner), NOT vol(10) (International) — cf. ls/vols.f:348-387 "COMPUTE SCRIBNER BOARD FOOT VOLUME".
# The r9cor correction (cf3 == cf4) and the DIB rounding (INT(dib+0.499)) are identical to International,
# so Scribner differs ONLY in this per-log kernel (and each log is nint-rounded before the sum).
const _R9_SCRBNR = let
    path = joinpath(@__DIR__, "..", "..", "data", "lakestates", "scribner_factor.csv")
    v = zeros(Float32, 120)
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        f = split(strip(line), ',')
        v[parse(Int, f[1])] = parse(Float32, f[2])
    end
    v
end
@inline function _r9_scrib_log(len::Float32, idib::Int)::Float32
    idib < 1 && return 0f0
    round(len * _R9_SCRBNR[min(idib, 120)], RoundNearestTiesAway)   # nint(len·scrbnr(idib))
end

# r9bdft (r9clark.f:1380) + R9LOGS bucking (r9logs.f R9LOGLEN/R9LOGDIB): board feet of the sawtimber
# section [stump, sawHt]. Bucks even-foot logs (shared R9LOGLEN rule, identical to `_r8_scribner_bf`),
# takes each log's small-end (top) DIB from the R9 taper `_r9_dib`. `logfn` maps (len, idib)→board feet:
# `_r9_intl_log` = vol(10) International ¼" (NE/CS `.sum`), `_r9_scrib_log` = vol(2) Scribner (LS `.sum`).
function _r9_bucked_bf(st::_R9State, sawHt::Float32, stump::Float32,
                       minLen::Float32, maxLen::Float32, trim::Float32, logfn::F,
                       logbuf::Union{Vector{Float32},Nothing} = nothing)::Float32 where {F}
    lmerch = sawHt - stump
    nlogp = clamp(floor(Int, lmerch / (maxLen + trim)), 0, 39)
    leftov = lmerch - (maxLen + trim) * nlogp - trim
    # logLen entries are written before read (1:tlogs, tlogs≤nlogp), so a reused dirty buffer is value-safe.
    logLen = logbuf === nothing ? zeros(Float32, 40) : logbuf; tlogs = 0
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
        bf += logfn(len, idib)
    end
    return round(bf)                            # r9bdft:1499 vol(10)=NINT(vol(10)) / vol(2)=NINT(vol(2))
end

# International ¼" (vol10 — NE/CS) and Scribner Decimal-C (vol2 — LS) board-foot wrappers over the
# shared R9LOGS bucking. Same log lengths + DIBs; only the per-log kernel differs.
@inline _r9_intlqtr_bf(st::_R9State, sawHt, stump, minLen, maxLen, trim, logbuf = nothing) =
    _r9_bucked_bf(st, sawHt, stump, minLen, maxLen, trim, _r9_intl_log, logbuf)
@inline _r9_scribner_bf(st::_R9State, sawHt, stump, minLen, maxLen, trim, logbuf = nothing) =
    _r9_bucked_bf(st, sawHt, stump, minLen, maxLen, trim, _r9_scrib_log, logbuf)

"""
    r9clark_cubic(spp, dbhOb, htTot, prod, mTopP, mTopS, stump) -> vol::Vector{Float32}

R9 Clark cubic volumes for FIA species `spp`. `prod`="01" (sawtimber) or "02"
(pulpwood). Returns vol[15]: vol[1]=total cuft, vol[4]=sawtimber cuft, vol[7]=
topwood cuft (so merch cuft = vol[4]+vol[7]). Mirrors r9clark's prod-01 cubic
branch + the NATCRS flags (CUTFLG=CUPFLG=SPFLG=1); board feet (vol[2]) is not yet
computed here (needs R9LOGS/r9bdft). NaN/short-tree guards follow the Fortran.
"""
function r9clark_cubic(spp::Int, dbhOb::Float32, htTot::Float32, prod::String,
                       mTopP::Float32, mTopS::Float32, stump::Float32,
                       bfTopP::Float32 = -1f0, bfStmp::Float32 = -1f0;
                       board_scribner::Bool = false,
                       vbuf::Union{Vector{Float32},Nothing} = nothing,
                       logbuf::Union{Vector{Float32},Nothing} = nothing)::Vector{Float32}
    vol = vbuf === nothing ? zeros(Float32, 15) : fill!(vbuf, 0f0)
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

    # A17/B17 select on topDib (=0 here → coef0 class) — i.e. the construction values, so no post-set.
    st = _R9State(co.r, co.c, co.e, co.p, co.a, co.b, co.a4, co.b4,
                  co.a17_0, co.b17_0, 0f0, 0f0, 0f0)
    dbhIb, dib17, err = _r9_dia417(st, topDib, dbhOb, topHt, 0f0, 0f0, sawDib, plpDib)
    err != 0 && return vol
    totHt = _r9_totht(htTot, dbhIb, dib17, topHt, topDib, st.a, st.b)
    totHt <= 17.3f0 && return vol
    # immutable _R9State ⇒ rebuild with the resolved dbhIb/dib17/totHt (stack-allocated, no heap)
    st = _R9State(co.r, co.c, co.e, co.p, co.a, co.b, co.a4, co.b4,
                  co.a17_0, co.b17_0, dbhIb, dib17, totHt)

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
    # FVS books vol(4)=cfVol and vol(7)=max(tcfVol−cfVol,0) inside `if(cupFlg.eq.1 .or. spFlg.eq.1)` —
    # NOT gated on sawHt>0. When the saw bole is too short (sawHt→0) FVS's r9cuft(stump,0)=0, so vol(4)
    # (saw cubic)=0 BUT vol(7) (topwood) = tcfVol−0 = the FULL merch cubic. Gating the whole block on
    # sawHt>0 (the old bug) dropped the merch bole for a sawtimber-SIZED tree with no valid sawlog — e.g.
    # a hardwood right at SCFMIND whose saw top (9.6") sits just below DBH so the sawlog can't make MERCHL.
    if iProd == 1
        sawHt = _r9_ht(st, sawDib)
        (topDib <= sawDib && topHt < sawHt) && (sawHt = topHt)
        sawHt < merchL + trim + stump && (sawHt = 0f0)
        scfVol = 0f0
        if sawHt > 0f0
            scfVol = _r9_cuft(st, stump, sawHt)
            short && (scfVol *= shrtHt / 17.3f0)
        end
        vol[4] = scfVol                                     # saw cubic (0 when no valid sawlog)
        spFlg && (vol[7] = max(tcfVol - scfVol, 0f0))       # topwood = merch − saw (= full merch when sawHt=0)
    end

    # Board feet (r9bdft International ¼") — computed for ANY board-eligible tree, INDEPENDENT of the sawtimber
    # prod class. FVS vols.f:354 books BFV whenever `D ≥ BFMIND .AND. D > BFTOPD` (the caller applies that gate),
    # NOT gated on SCFMIND/prod — so a tree with a raised SCFMIND (VOLUME) that is prod-02 for the cubic still
    # gets board feet if it clears BFMIND (the volume_override case). The board bole runs to the BOARD merch top
    # + stump (BFVOLUME sp_bf_topd/sp_bf_stump); when those equal the sawtimber values `bfHt≡sawHt`, so a prod-01
    # bf==saw tree is BIT-EXACT (same _r9_ht + topHt-cap + short-log checks). cf correction (_r9_cor! vol(10),
    # cf4≡cf3) applies via vol[2].
    bfDib = bfTopP > 0f0 ? bfTopP : sawDib
    bfSt  = bfStmp > 0.01f0 ? bfStmp : stump
    bfHt = _r9_ht(st, bfDib)
    (topDib <= bfDib && topHt < bfHt) && (bfHt = topHt)
    bfHt < merchL + trim + bfSt && (bfHt = 0f0)
    bfHt > 0f0 && (vol[2] = board_scribner ? _r9_scribner_bf(st, bfHt, bfSt, minLen, maxLen, trim, logbuf) :
                                             _r9_intlqtr_bf(st, bfHt, bfSt, minLen, maxLen, trim, logbuf))

    # FVS r9clark.f:454-462 applies the correction (r9cor) FIRST, THEN nint-rounds. The port had these
    # REVERSED (round-half-even then cor), so each printed value was cor×(0.1-multiple) instead of the clean
    # nint(cor×raw) FVS emits — a consistent ~0.3% high bias on TCuFt/BdFt (per-tree within .trl print
    # precision, but it accumulated ×TPA). nint = round-half-AWAY-from-zero (Fortran NINT), not Julia's
    # default round-half-to-even. Board feet (vol[2]) rounds to a whole number (r9clark.f:457).
    _r9_cor!(vol, spp, iProd)
    rnd(x) = round(x, RoundNearestTiesAway)
    vol[1] > 0f0 && (vol[1] = rnd(vol[1] * 10f0) / 10f0)
    vol[2] > 0f0 && (vol[2] = rnd(vol[2]))
    vol[4] > 0f0 && (vol[4] = rnd(vol[4] * 10f0) / 10f0)
    vol[7] > 0f0 && (vol[7] = rnd(vol[7] * 10f0) / 10f0)
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

# CS merch standards (cs/sitset.f:130-227). Softwoods = species 1-7 (RC/JU/SP/VP/LP/OS/WP);
# eastern redcedar (sp 1) gets a lower sawtimber min, but only on forest IFOR=1 (Mark Twain).
# Board-foot mins equal the sawtimber-cubic mins (BFMIND==SCFMIND, BFTOPD==SCFTOPD), so the
# returned (scfmind, scftopd) cover both — bf-equal, like _ne_merch. Returns the 6-tuple
# (dbhmin, topd, scfmind, scftopd, stmp, scfstmp) the eastern volume driver consumes.
@inline function _cs_merch(spi::Integer, ifor::Integer)
    if spi <= 7                           # softwoods
        scfmind = (ifor == 1 && spi == 1) ? 6f0   : 9f0
        scftopd = (ifor == 1 && spi == 1) ? 5f0   : 7.6f0
        return (5f0, 4f0, scfmind, scftopd, 0.5f0, 1f0)
    else                                  # hardwoods
        dbhmin  = ifor == 1 ? 5f0   : 6f0
        topd    = ifor == 2 ? 5f0   : 4f0
        scfmind = ifor == 1 ? 9f0   : 11f0
        scftopd = ifor == 1 ? 7.6f0 : 9.6f0
        return (dbhmin, topd, scfmind, scftopd, 0.5f0, 1f0)
    end
end

# LS merch standards (ls/sitset.f:305-385). Softwoods = species index ≤ 14. Aspen/poplar
# (species 40-42 = BT/QA/BP) get raised sawtimber mins on some forests. Board-foot mins equal
# the sawtimber-cubic mins for LS in every IFOR case (BFMIND==SCFMIND, BFTOPD==SCFTOPD), so the
# returned (scfmind, scftopd) cover both — bf-equal. Stumps from grinit.f (STMP 0.5 / SCFSTMP 1.0).
# TOPD (merch-cubic top) is a flat 4 for LS (no IFOR dependence, unlike CS). Returns the 6-tuple
# (dbhmin, topd, scfmind, scftopd, stmp, scfstmp) the eastern volume driver consumes.
@inline function _ls_merch(spi::Integer, ifor::Integer)
    if spi <= 14                          # softwoods
        return (5f0, 4f0, 9f0, 7.6f0, 0.5f0, 1f0)
    else                                  # hardwoods
        aspen = (40 <= spi <= 42)
        dbhmin  = ifor == 2 ? (aspen ?  6f0 :  5f0) : (ifor == 6 ? 6f0 : 5f0)
        scfmind = ifor == 2 ? (aspen ? 11f0 :  9f0) : (ifor == 5 ? (aspen ? 9f0 : 11f0) : 11f0)
        scftopd = ifor == 2 ? (aspen ? 9.6f0 : 7.6f0) : (ifor == 5 ? 7.6f0 : 9.6f0)
        return (dbhmin, 4f0, scfmind, scftopd, 0.5f0, 1f0)
    end
end

# Fallback SN forest when a stand carries no forest code — a neutral non-special forest so `_sn_merch`
# yields the region-8 defaults (the North Carolina IFOR=11 and IFOR=10 branches never fire on it).
const _SN_DEFAULT_IFOR = 1

# SN merch standards (setcubicdflts.f:350-414, the region-8 block). Softwoods = species index ≤ 17
# (plus sp 88); hardwoods otherwise. FVS overwrites every species from this block keyed on
# softwood/hardwood + IFOR — the only non-default forests are North Carolina (IFOR=11) and a single
# IFOR=10 softwood case. The two NC coastal districts (KODIST 3/10) get yet-lower tops, but KODIST is
# not plumbed through the FIA reader, so the common non-coastal NC branch is used. Board-foot mins are
# bf-equal (like the eastern variants); stumps STMP 0.5 / SCFSTMP 1.0 (grinit.f). Returns the 6-tuple
# (dbhmin, topd, scfmind, scftopd, stmp, scfstmp). For every non-NC forest this reproduces the SN
# merch_specs.csv defaults exactly (softwood 7/10, hardwood 9/12, top 4), so non-NC stays bit-exact.
@inline function _sn_merch(spi::Integer, ifor::Integer, kodist::Integer)
    nc = ifor == 11                          # North Carolina overrides (setcubicdflts.f:363-413)
    coastal = nc && (kodist == 3 || kodist == 10)   # NC coastal districts get lower tops (:368,388,411)
    softwood = spi <= 17 || spi == 88
    topd = nc ? 3.5f0 : 4f0                  # :364-365 (NC) / :360
    dbhmin = if nc                           # :368-373
        coastal ? (softwood ? 5.6f0 : 6f0) : 8f0
    else
        spi in (7,13,39,43,44,52,53,55,63) ? 6f0 : 4f0
    end
    if softwood                              # :377-402
        if nc
            scfmind, scftopd = coastal ? (11f0, 6.3f0) :
                               (spi in (2,12,15,16,17) ? (12f0, 9f0) : (10f0, 6.3f0))
        elseif ifor == 10 && spi == 2
            scfmind, scftopd = 9f0, 7f0
        else
            scfmind, scftopd = 10f0, 7f0
        end
    else                                     # hardwoods :404-413
        scfmind, scftopd = nc ? (coastal ? (13f0, 8f0) : (15f0, 11f0)) : (12f0, 9f0)
    end
    return (dbhmin, topd, scfmind, scftopd, 0.5f0, 1f0)
end

# Region-9 Clark `.sum` board type is per-national-forest (volinit.f:434-451, the R9 branch):
# after R9CLARK fills both vol2 (Scribner) and vol10 (International ¼"), FVS overwrites
# vol2←vol10 for the listed R9 forests, so those report International and every other R9
# forest reports Scribner. IFORST is the 2-digit forest within KODFOR/LOCATION (LOCATION mod
# 100, or the middle two digits when KODFOR>10000, per fvsvol.f:89-96). Confirmed against live
# FVSls (debug-stamped volinit): FIA stand LOCATION=924→IFORST=24∈list→International; lst01
# LOCATION=903→IFORST=3∉list→Scribner. LS defaults to Scribner ⇒ only these forests flip.
const _R9_INTL_BDFT_FORESTS = (4, 5, 8, 11, 12, 14, 19, 20, 21, 22, 24, 30)
_r9_iforst(kodfor::Integer) = (k = Int(kodfor); k > 10000 ? (k ÷ 100) % 100 : k % 100)

"""
    compute_volumes_ne!(s)

NE per-tree volume driver (the Region-9 analogue of the SN `compute_volumes!`
body): NVEL R9 Clark cubic via `r9clark_cubic`, with the NE IFOR-dependent merch
standards. Loads cuft/merch-cuft/saw-cuft/board-feet per tree (board feet = Scribner
via R9LOGS bucking + r9bdft). Broken-top (CFTOPK) reuse is TODO.
"""
function compute_volumes_ne!(s::StandState)
    # LS `.sum` BdFt is Scribner (vol2) EXCEPT on the R9 forests that FVS maps to International
    # ¼" (vol10) per _R9_INTL_BDFT_FORESTS; NE/CS report International ¼" (vol10). ls/vols.f:348-387.
    board_scribner = (s.variant isa LakeStates) &&
                     !(_r9_iforst(s.plot.user_forest_code) in _R9_INTL_BDFT_FORESTS)
    t = s.trees; co = s.coef
    cs = s.variant isa CentralStates
    ifor = Int(s.plot.forest_idx); ifor == 0 && (ifor = cs ? 1 : _NE_DEFAULT_IFOR)
    # Merch standards live in Control.sp_* (init_merch_standards! seeds them from _ne_merch/_cs_merch, the
    # IFOR code rules — so values are IDENTICAL for a stand with no override), and VOLUME/BFVOLUME override
    # those same fields (apply_volume_overrides!). Read them here (not a fresh _ne_merch call) so the R9 path
    # honors VOLUME/BFVOLUME just like the SN Clark path — FVS's merch standards are one overridable common.
    init_merch_standards!(s)
    md = s.control
    # Volume defect (FVS vols.f:285-432) — the SAME ICDF/IBDF block the SN path uses; vols.f is the shared
    # driver, so the R9 (NE/CS) volumes get defect-corrected too. Gated on any defect source being present.
    ctl = s.control
    cfdef = ctl.sp_cf_defect; bfdef = ctl.sp_bf_defect        # MCDEFECT / BFDEFECT DBH curves
    cff0 = ctl.sp_cf_form0; cff1 = ctl.sp_cf_form1            # MCFDLN cubic log-linear form coefs
    bff0 = ctl.sp_bf_form0; bff1 = ctl.sp_bf_form1            # BFFDLN board log-linear form coefs
    anydef_cf = any(!iszero, cfdef); anydef_bf = any(!iszero, bfdef)
    anyform = any(!iszero, cff0) || any(!=(1f0), cff1) || any(!iszero, bff0) || any(!=(1f0), bff1)
    anydef = anydef_cf || anydef_bf || anyform || any(!iszero, t.defect)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
            continue
        end
        # Broken-top trees (FVS vols.f:146 `IF(TKILL) H=NORMHT/100`): build the volume profile from the
        # NORMAL (predicted full) height — `norm_ht`, resolved by dub_missing_heights! — then truncate it
        # back to the break with CFTOPK/BFTOPK (vols.f:193). Without this a top-killed tree's cubic was
        # built on the SHORT broken height (net01 SM d10.4 HTTOPK49: jl TOT 13.8 vs live 15.4).
        tkill = h >= 4.5f0 && t.trunc[i] > 0
        tkill && (h = Float32(t.norm_ht[i]) * 0.01f0)
        fias = strip(string(co.code_fia[sp]))
        fia = isempty(fias) ? 0 : parse(Int, fias)
        # D35: VOLUME field-7 METHC==5 selects the CS DVEE/Gevorkiantz model ('900DVEE', r9vol.f R9VOL)
        # instead of the Clark taper. iforst = KODFOR−900 (the R9 forest number → LS/CS/NE region in R9_MHTS).
        if md.sp_methc[sp] == 5
            tcf, mcf, scf, _bf = r9vol_gevorkiantz(fia, d, h, Int(s.plot.user_forest_code) - 900)
            # DVEE merch DBHMIN gate (fvsvol.f:512). CS default (sitset.f:130-141) = 5 for softwoods AND, for
            # the region-9 DVEE forests (905/908 ⇒ IFOR=1 case), hardwoods too — a stamp of live's DBHMIN(ISPC)
            # gave 5 for SP, and DBHMIN=6 for BH REGRESSED Mcuft (3093→2881 vs live 3090), confirming 5. Fall
            # back to 5 when unset (the planted DVEE species carry sp_dbh_min=0).
            _dv_dbhmin = md.sp_dbh_min[sp] > 0f0 ? md.sp_dbh_min[sp] : 5f0
            d < _dv_dbhmin && (mcf = 0f0)
            # BOARD feet: `VOLUME …5` sets only METHC=5 (cubic); METHB stays 6 ⇒ VEQNNB='900CLKE' (sitset.f:260-
            # 270). So the board is the CLARK model, NOT DVEE — compute it via r9clark like the CLKE path below.
            _scfm = md.sp_scf_dbhmin[sp]; _topd = md.sp_top_diam[sp]
            _prod = d >= _scfm ? "01" : "02"
            _mtopp = d >= _scfm ? md.sp_scf_topd[sp] : _topd
            _vc = r9clark_cubic(fia, d, h, _prod, _mtopp, _topd, 0f0, md.sp_bf_topd[sp], md.sp_bf_stump[sp];
                                board_scribner = board_scribner,
                                vbuf = s.scratch.r9_vol, logbuf = s.scratch.r9_logbuf)
            bf = (d >= md.sp_bf_dbhmin[sp] && d > md.sp_bf_topd[sp]) ? _vc[2] : 0f0
            if anydef
                mcf, scf, bf = _apply_tree_defect(mcf, scf, bf, d, sp, Int(t.defect[i]),
                                                  cfdef, bfdef, cff0, cff1, bff0, bff1, anydef_cf, anydef_bf)
            end
            t.cuft_vol[i] = tcf; t.merch_cuft_vol[i] = mcf
            t.saw_cuft_vol[i] = scf; t.bdft_vol[i] = bf
            continue
        end
        dbhmin = md.sp_dbh_min[sp]; topd = md.sp_top_diam[sp]
        scfmind = md.sp_scf_dbhmin[sp]; scftopd = md.sp_scf_topd[sp]
        stmp = md.sp_stump_ht[sp]; scfstmp = md.sp_scf_stump[sp]
        bfmind = md.sp_bf_dbhmin[sp]                           # BFVOLUME board-foot min DBH (bf-equal by default)
        prod = d >= scfmind ? "01" : "02"
        mtopp = d >= scfmind ? scftopd : topd
        v = r9clark_cubic(fia, d, h, prod, mtopp, topd, 0f0, md.sp_bf_topd[sp], md.sp_bf_stump[sp];
                          board_scribner = board_scribner,
                          vbuf = s.scratch.r9_vol, logbuf = s.scratch.r9_logbuf)
        tcf = v[1]
        mcf = d >= dbhmin  ? v[4] + v[7] : 0f0
        scf = d >= scfmind ? v[4] : 0f0
        bf  = (d >= bfmind && d > md.sp_bf_topd[sp]) ? v[2] : 0f0   # board feet: FVS vols.f:354 D≥BFMIND & D>BFTOPD
        if tkill && tcf > 0f0
            bark = bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, d)
            # _ne_merch returns per-species SCALARS; wrap as 1-tuples so cftopk/bftopk's `merch.x[sp]`
            # works with sp=1. The BOARD top-kill (bftopk) uses the BOARD's own top/stump (sp_bf_topd/sp_bf_stump),
            # NOT the sawtimber scftopd/scfstmp — bf-equal by default so this is inert, but a VOLUME card with a
            # blank SCFTOPD zeroes scftopd (keeping sp_bf_topd) ⇒ a broken-top tree's board must not follow the
            # sawtimber-cubic top. (The non-tkill path at r9clark_cubic already passes sp_bf_topd for the board.)
            mk = (stmp = (stmp,), topd = (topd,), scfstmp = (scfstmp,), scftop = (scftopd,),
                  bftopd = (md.sp_bf_topd[sp],), bfstmp = (md.sp_bf_stump[sp],))
            tcf, mcf, scf = cftopk(mk, 1, d, h, tcf, mcf, scf, v[1], bark, Int(t.trunc[i]))
            bf = bftopk(mk, 1, d, h, bf, v[1], bark, Int(t.trunc[i]))
        end
        if anydef
            mcf, scf, bf = _apply_tree_defect(mcf, scf, bf, d, sp, Int(t.defect[i]),
                                              cfdef, bfdef, cff0, cff1, bff0, bff1, anydef_cf, anydef_bf)
        end
        t.cuft_vol[i]      = tcf
        t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i]   = scf
        t.bdft_vol[i]       = bf
    end
    return s
end
