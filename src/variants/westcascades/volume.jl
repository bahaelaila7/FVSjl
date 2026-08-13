# =============================================================================
# volume.jl (westcascades) — WC volume (R6 NVEL / VOLEQDEF). Chunk 8.
#
# WC forest-618 (Willamette) VOLEQ (dumped from FVSwc_clean; wc/sitset.f:232 → NVEL voleqdef.f R6_EQN
# westside branch, VAR='WC', FORNUM=KODFOR%100):
#   • WESTSIDE Flewelling  F05FW2W202 (DF), F03FW2W263 (WH), …W242 (RC) — f_west.f SHP_W3/W4/W5,
#     INSIDE-bark profile calibrated to DBHIB=D−FDBT_C1 (mirror of the INGY dibat path). NEW port.
#   • INGY Flewelling  I00FW2W017 (GF), I00FW2W108 (NF), I00FW2W073 (IC) — reuse cr_fw2_vol (SHP_C2).
#   • Behre  616BEHW<fia>  — everything else — reuse bm_r6vol3/r6dibs/r6vol1 + WC form-class (wc/formcl.f).
# Merch (wc/sitset.f westside CASE DEFAULT / IFOR 6): TOPD=BFTOPD=4.5, STUMP=1, DBHMIN=BFMIND=7 (LP sp11=6).
#   .sum: TCuFt=VOL(1), MerchCuFt=VOL(4) (D≥DBHMIN), MerchBdFt=VOL(2) Scribner (D≥BFMIND).
#
# wct01 species present: DF(F05FW2W202, westside), LP/PP/SP/WF/ES(616BEHW, Behre). GF/NF/IC route to INGY.
# =============================================================================

# --- WC form class (wc/formcl.f): 6 forest tables (39×5) + 4 BLM tables (39×1). ---
# forest key = forkod JFOR index (plot.forest_idx): 1..6 = GIFPFC/MBSNFC/MTHDFC/ROGRFC/UMPQFC/WILLFC;
# 7..10 = BLM708/709/710/711. IFCDBH=INT((D−1)/10+1), clamp≥1, D>40.9→5.
const WC_FORMCL = let
    T = Dict{Tuple{Int,Int,Int},Int}()   # (ifor, sp, dbhclass) → fc
    for l in readlines(joinpath(WC_DATADIR, "formcl_wc.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        T[(parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]))] = parse(Int, f[4])
    end
    # densify into a (10, 39, 5) array (BLM tables 7..10 replicate their single class into all 5 slots)
    A = zeros(Int, 10, 39, 5)
    for ((fi, sp, dc), v) in T
        if fi <= 6
            A[fi, sp, dc] = v
        else
            for c in 1:5; A[fi, sp, c] = v; end
        end
    end
    A
end

function wc_formcl(sp::Integer, ifor::Int, d::Real)::Int
    (sp < 1 || sp > 39) && return 80
    (ifor < 1 || ifor > 10) && (ifor = 6)       # default Willamette
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1)
    Float32(d) > 40.9f0 && (ifcdbh = 5)
    v = WC_FORMCL[ifor, sp, ifcdbh]
    v == 0 ? 80 : v
end

# --- WC R6_EQN westside VOLEQ (voleqdef.f VAR='WC'). Implemented for FORNUM 18 (Willamette, wct01);
#     other westside forests fall back to region-6 Behre 616BEHW<fia> (their FW2 overrides deferred). ---
function _wc_r6_eqn(fornum::Int, fia::Int)::String
    beh = "616BEHW" * lpad(string(fia), 3, '0')
    if fornum == 18                                  # Willamette (validated vs FVSwc_clean)
        fia == 17  && return "I00FW2W017"            # GF → INGY
        fia == 22  && return "I00FW2W108"            # NF → INGY (lodgepole eqn)
        fia == 81  && return "I00FW2W073"            # IC → INGY (incense-cedar eqn 073)
        fia == 202 && return "F05FW2W202"            # DF → westside Flewelling
        fia == 263 && return "F03FW2W263"            # WH → westside Flewelling
    end
    return beh
end

# ---------------------------------------------------------------------------
# Westside Flewelling shape (f_west.f SHP_W3/W4/W5) + breast-height bark (FDBT_C1).
# Fills RFLW(6)/RHFW(4) consumed UNCHANGED by _fw2_sf_taper / _fw2_sf_yhat / _fw2_tcubic. Double-precision
# kernel (REAL*4 out), matching the proven cr_fw2 pattern. GEOSUB '01'..'08' → regional coefficient index.
# ---------------------------------------------------------------------------
@inline function _wc_geoidx(geosub::AbstractString)::Int
    (length(geosub) == 2 && geosub[1] == '0' && '1' <= geosub[2] <= '8') ? Int(geosub[2] - '0') : 0
end
@inline _logit(u::Float64) = exp(u) / (1.0 + exp(u))

const _W3_R25 = (-2.0262, -1.7945, -2.0366, -2.0811, -1.9868, -2.0151, -1.9475, -2.0151)
const _W3_R34 = ( 5.2132,  5.1417,  5.1768,  5.1156,  5.1654,  5.2413,  5.1427,  5.2413)

function _fw2_shp_w3(d::Float32, h::Float32, geosub::AbstractString)
    D = Float64(d); H = Float64(h)
    F13=-0.07799267; F14=0.8096211; F15=-0.9247244
    F17=0.166749;    F18=-9.20884;  F19=0.1212094
    F21=5.0719922;   F22=-0.12555733; F23=-1.6408313; F24=-0.3005388
    F25=-2.000088;   F26=-0.410677
    F29=-1.428671;   F30=-0.8660924; F31=0.0543758
    F33=0.80981083;  F34=5.1785820;  F35=1.8935219;  F36=-2.420031; F37=-0.009232401
    F38=-11.49305;   F39=2.037438;   F40=-0.0114178
    F42=2.246288;    F43=-1.1410822; F44=1.7293166
    F45=-3.18221;    F46=-2.119838;  F47=1.3756086
    fr25 = F25; fr34 = F34
    if geosub != "00"
        id = _wc_geoidx(geosub)
        id != 0 && (fr25 = _W3_R25[id]; fr34 = _W3_R34[id])
    end
    dmedian = 0.566 * (H - 4.5)^(0.634 + 0.00074 * H); dform = D / dmedian - 1.0
    u7 = F13 + F14*log(D + 1.0) + F15*log(H)
    u9a = clamp(F18 + F19*H, -7.0, 7.0); u9 = F17 * _logit(u9a)
    u8 = F21 + F22*log(H) + F23*dform + F24*(D/10.0)^1.5
    u1 = fr25 + F26*log(H)
    u2 = F29 + F30*log(D) + F31*D
    u3 = fr34 + F35*log(D + 1.0) + F36*log(H) + F37*dform*H + F33*dform
    u4 = F38 + F39*D + F40*H*D
    u5 = F42 + F43*H + F44*dform
    u1=clamp(u1,-7.,7.); u2=clamp(u2,-7.,7.); u3=clamp(u3,-7.,7.); u4=clamp(u4,-7.,7.)
    u5=clamp(u5,-7.,7.); u7=clamp(u7,-7.,7.); u8=clamp(u8,-7.,7.)
    u6 = clamp(F45 + F46*log(D + 1.0) + F47*log(H), -6.0, 6.0); u6 = 1.0 + exp(u6)
    u6 = u6 < 1.005 ? 1.005 : (u6 > 100.0 ? 100.0 : u6)
    _wc_assemble_shp(u1,u2,u3,u4,u5,u6,u7,u8,u9)
end

function _fw2_shp_w4(d::Float32, h::Float32, geosub::AbstractString)
    D = Float64(d); H = Float64(h)
    F13=-3.1137977; F14=1.1996084; F15=-0.01195901
    F17=-0.066829984; F18=0.026990398; F19=0.24661021
    F21=-8.3415555; F22=2.4274384; F23=-5.6918023; F24=0.56487213
    F25=-7.6464554; F26=5.1709049; F27=-2.7133381; F28=0.38918349
    F29=-7.0
    F34=-1.2882571; F35=35.688410; F36=0.17995769; F37=1.5565605
    F38=6.4397446; F39=-1.3439736; F40=6.3442558
    F42=11.898092; F43=-3.6789851; F44=0.15168209
    F45=-1.3248733; F46=-0.11788962; F47=-0.015909154
    R25 = (-7.511,-7.687,-7.224,-7.355,-7.632,-7.646,-7.687,-7.911)
    R34 = (-1.215,-1.355,-1.177,-1.373,-1.398,-1.188,-1.355,-1.449)
    fr25 = F25; fr34 = F34
    if geosub != "00"
        id = _wc_geoidx(geosub); id != 0 && (fr25 = R25[id]; fr34 = R34[id])
    end
    dmedian = 0.2855*(H-4.5)^(0.307 - 0.00505*H + 0.00001745*H*H + 0.19*log(H)); dform = D/dmedian - 1.0
    u7 = clamp(F13 + F14*(1.0 - exp(F15*H)), -7.0, 7.0); rhi1 = _logit(u7); rhi1 > 0.5 && (rhi1 = 0.5)
    u9a = F17 + F18*log(H) + F19*dform
    (rhi1 + u9a > 0.75) && (u9a = 0.75 - rhi1); u9 = max(u9a, 0.0)
    u8 = F21 + F22*log(H) + F23*dform + F24*log(H)*dform
    u1 = fr25 + F26/(H/100) + F27/(H/100)^2 + F28/(H/100)^3
    u2 = F29
    u3 = fr34 + F35/H + F36*H*D/1000 + F37*dform
    u4 = F38 + F39*log(H) + F40*dform
    u5 = F42 + F43*log(H) + F44*D
    u1=clamp(u1,-7.,7.); u2=clamp(u2,-7.,7.); u3=clamp(u3,-7.,7.); u4=clamp(u4,-7.,7.); u5=clamp(u5,-7.,7.); u8=clamp(u8,-7.,7.)
    u6 = clamp(F45 + F46*log(D) + F47*H, -6.0, 6.0); u6 = 1.0 + exp(u6)
    u6 = u6 < 1.005 ? 1.005 : (u6 > 100.0 ? 100.0 : u6)
    _wc_assemble_shp_rhi(u1,u2,u3,u4,u5,u6,rhi1,u8,u9)
end

function _fw2_shp_w5(d::Float32, h::Float32, geosub::AbstractString)
    D = Float64(d); H = Float64(h)
    F13=-1.1917515; F15=0.15632952
    F17=0.44805624; F18=-0.10269221
    F21=0.828142; F22=17.750205
    F25=-17.959260; F26=3.8308966; F27=45.163788
    F29=-1.7546090; F30=-7.9230813; F31=-112.53761
    F34=-0.17943249; F35=0.15534485; F36=-0.032777026
    F38=8.5305448; F39=-0.83350599; F40=12.100013
    F42=8.3021283; F43=-150.0
    F45=-7.1594841; F46=0.11263465; F47=22.396603
    dmedian = 0.11*(H-4.5)^(1.08 + 0.0006*H); dform = D/dmedian - 1.0
    u7 = F13 + F15*dform; u7 = clamp(u7, -7.0, 1.0); rhi1 = _logit(u7)   # W5: upper clamp 1.0, NO 0.5 cap
    u9a = F17 + F18*log(H); (rhi1 + u9a > 0.75) && (u9a = 0.75 - rhi1); u9 = max(u9a, 0.0)
    u8 = F21 + F22/D
    u1 = F25 + F26*log(D) + F27/D
    u2 = F29 + F30/D + F31/(D*D)
    u3 = F34 + F35*D + F36*H
    u4 = F38 + F39*D + F40*dform
    u5 = F42 + F43/D
    u1=clamp(u1,-7.,7.); u2=clamp(u2,-7.,7.); u3=clamp(u3,-7.,7.); u4=clamp(u4,-7.,7.); u5=clamp(u5,-7.,7.); u8=clamp(u8,-7.,7.)
    u6 = clamp(F45 + F46*D + F47/D, -6.0, 6.0); u6 = 1.0 + exp(u6)
    u6 = u6 < 1.005 ? 1.005 : (u6 > 100.0 ? 100.0 : u6)
    _wc_assemble_shp_rhi(u1,u2,u3,u4,u5,u6,rhi1,u8,u9)
end

# Shared tail: build (RFLW, RHFW) from the U-transforms. Two variants differ only in whether RHI1's
# 0.5 cap is applied to the raw logit (W3) or already applied (W4/W5 pass rhi1 in).
@inline function _wc_assemble_shp(u1,u2,u3,u4,u5,u6,u7,u8,u9)
    rhi1 = _logit(u7); rhi1 > 0.5 && (rhi1 = 0.5)
    _wc_assemble_shp_rhi(u1,u2,u3,u4,u5,u6,rhi1,u8,u9)
end
@inline function _wc_assemble_shp_rhi(u1,u2,u3,u4,u5,u6,rhi1,u8,u9)
    r1=_logit(u1); r2=_logit(u2); r3=_logit(u3); r4=_logit(u4)
    r5 = 0.5 + 0.5*_logit(u5); a3 = u6
    rhlongi = u9; rhi2 = rhi1 + rhlongi
    rhc = rhi2 + (1.0 - rhi2)*_logit(u8)
    ((Float32(r1),Float32(r2),Float32(r3),Float32(r4),Float32(r5),Float32(a3)),
     (Float32(rhi1),Float32(rhi2),Float32(rhc),Float32(rhlongi)))
end

# FDBT_C1 (f_west.f:895) — double bark thickness at BH; DBHIB = D − FDBT_C1. JSP 3=DF,4=WH,5=RC.
function _wc_fdbt_c1(jsp::Int, geosub::AbstractString, d::Float32, h::Float32)::Float32
    D = Float64(d); H = Float64(h); id = _wc_geoidx(geosub)
    if jsp == 3
        dmedian = 0.566*(H-4.5)^(0.634 + 0.00074*H); dform = D/dmedian - 1.0
        ratio = if geosub == "00" || id == 0
            exp(-2.4641 + 0.04393*log(D) - 0.2922*dform + 0.05964*dform*log(D))
        else
            roff = (0.0,0.117,0.121,0.133,0.025,0.028,0.088,0.028)[id]
            exp(-2.5087 + 0.03600*log(D) - 0.4086*dform + 0.10120*dform*log(D) + roff)
        end
        return Float32(ratio * D)
    elseif jsp == 4
        ratio = if geosub == "00" || id == 0
            0.04504*(1.0 + 0.8307*exp(-0.2048*D))
        else
            roff1 = (0.0,0.118,0.204,0.105,0.058,0.061,0.118,0.071)[id]
            0.04221*(1.0 + 0.8836*exp(-0.2145*D))*(1.0 + roff1)
        end
        return Float32(ratio * D)
    else                                          # jsp==5 RC
        du = max(D, 3.8); ratio = 0.01949*(1.0 + 15.599/du - 29.212/du^2)
        return Float32(ratio * D)
    end
end

# R6 westside log-bucking constants (NVEL mrules.f REGN 6): OPT=23 (full-log-first SEGMNT, ≠ R3's 22),
# EVOD=2, MAXLEN=16, MINLEN=2, TRIM=0.5, MERCHL=8, board COR='N' (full board feet AINT(len·volfac+.5), not
# the R3 decimal-C ×10). Inside-bark log small-end diameters (LOGDIA(:,2)), merch top = TOPD·BARK (IB).
const _WC_R6_OPT = 23; const _WC_R6_EVOD = 2
const _WC_R6_MAXLEN = 16.0f0; const _WC_R6_MINLEN = 2.0f0; const _WC_R6_TRIM = 0.5f0; const _WC_R6_MERCHL = 8.0f0

"R6 westside merch cubic VOL(4): buck stump→mtop (IB) with the R6 (OPT=23) SEGMNT, per-log 0.1-rounded Smalian."
function _wc_fw2_merch_cuft(dibat, h::Float32, mtop::Float32, stump::Float32)::Float32
    hs = _fw2_hs(dibat, mtop, h); lmerch = hs - stump
    lmerch < _WC_R6_MERCHL && return 0f0
    ns = _nvb_numlog(_WC_R6_OPT, _WC_R6_EVOD, lmerch, _WC_R6_MAXLEN, _WC_R6_MINLEN, _WC_R6_TRIM); ns == 0 && return 0f0
    loglen, ns = _nvb_segmnt(_WC_R6_OPT, _WC_R6_EVOD, lmerch, _WC_R6_MAXLEN, _WC_R6_MINLEN, _WC_R6_TRIM, ns)
    dibl = _fw2_dclass(dibat(4.5f0)); ht2 = stump; vol4 = 0f0
    @inbounds for i in 1:ns
        ht2 += _WC_R6_TRIM + loglen[i]; dib = dibat(ht2)
        (i == ns && dib < mtop) && (dib = mtop)
        dibs = _fw2_dclass(dib)
        vol4 += floor(0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen[i] * 10f0 + 0.5f0) / 10f0
        dibl = dibs
    end
    return vol4
end

"R6 westside Scribner board VOL(2): R6 (OPT=23) SEGMNT + SCRIB COR='N' (full board feet per log)."
function _wc_fw2_board(dibat, h::Float32, bftop::Float32, stump::Float32)::Float32
    hs = _fw2_hs(dibat, bftop, h); lmerch = hs - stump
    lmerch < _WC_R6_MERCHL && return 0f0
    ns = _nvb_numlog(_WC_R6_OPT, _WC_R6_EVOD, lmerch, _WC_R6_MAXLEN, _WC_R6_MINLEN, _WC_R6_TRIM); ns == 0 && return 0f0
    loglen, ns = _nvb_segmnt(_WC_R6_OPT, _WC_R6_EVOD, lmerch, _WC_R6_MAXLEN, _WC_R6_MINLEN, _WC_R6_TRIM, ns)
    ht2 = stump; vol2 = 0f0
    @inbounds for i in 1:ns
        ht2 += _WC_R6_TRIM + loglen[i]; dib = dibat(ht2)
        (i == ns && dib < bftop) && (dib = bftop)
        vol2 += _scrib(_fw2_dclass(dib), loglen[i], 'N')   # COR='N': AINT(len·volfac+.5), no ×10, no exception
    end
    return vol2
end

# VOLEQ(1)='F' westside per-tree volume. Returns (VOL1 tcuft, VOL4 merch-cuft, VOL2 Scribner-bf).
# INSIDE-bark profile calibrated to DBHIB = D·BARK (the fvsvol variant bark, wc_bratio — passed as DBTBH=
# D·(1−BARK) into sf_shp's DBT_USER, so FDBT_C1 is bypassed): the SAME convention as the INGY dibat path.
# Merch/board buck with the R6 rules (mrules.f REGN 6: OPT=23 SEGMNT + board COR='N'); tops = TOPD·BARK (IB).
function wc_fw2_westside_vol(voleq::AbstractString, d::Float32, h::Float32, bark::Float32;
                            topd::Float32 = 4.5f0, bftopd::Float32 = 4.5f0, stump::Float32 = 1f0)
    (d < 1f0 || h <= 5f0) && return (0f0, 0f0, 0f0)
    spec = voleq[8:10]
    jsp = (spec == "202" || spec == "205" || spec == "204") ? 3 : spec == "263" ? 4 : spec == "242" ? 5 : 0
    jsp == 0 && return (0f0, 0f0, 0f0)
    geosub = length(voleq) >= 3 ? voleq[2:3] : "00"
    rflw, rhfw = jsp == 3 ? _fw2_shp_w3(d, h, geosub) :
                 jsp == 4 ? _fw2_shp_w4(d, h, geosub) : _fw2_shp_w5(d, h, geosub)
    tapcoe = _fw2_sf_taper(rhfw, rflw)
    dbhib = d * bark                                  # DBHIB from the WC variant bark (== INGY convention)
    yhat_bh = _fw2_sf_yhat(4.5f0 / h, tapcoe, rhfw, rflw, 1.0f0)
    yhat_bh == 0f0 && return (0f0, 0f0, 0f0)
    f = dbhib / yhat_bh
    dibat = ht -> _fw2_sf_yhat(ht / h, tapcoe, rhfw, rflw, f)
    stump_dib = h <= 15f0 ? _fw2_fwsmall(jsp, h, dibat(1.0f0), dbhib) : -1f0
    v1 = Float32(round(_fw2_tcubic(dibat, h; stump_dib = stump_dib) * 10.0f0)) / 10.0f0
    v4 = _wc_fw2_merch_cuft(dibat, h, topd * bark, stump)
    v2 = _wc_fw2_board(dibat, h, bftopd * bark, stump)
    return (max(v1, 0f0), max(v4, 0f0), max(v2, 0f0))
end

# --- Behre (616BEHW) per-tree volume — reuse the BM R6 machinery + WC form class. ---
function wc_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32)
    fclass = wc_formcl(sp, ifor, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    vol2 = 0f0; vol4 = 0f0
    v1 = if h <= 17.3f0
        0.00272708f0 * dbhib * dbhib * h            # R6VOL short-tree cylinder (R6DIBS/R6VOL1 skipped)
    else
        v = bm_r6vol3(d, dbtbh, fclass, h, 1)
        mtopp = 4.5f0 * bark                         # TOPDIAM = TOPD·BARK
        xlogs, ld1 = bm_r6dibs(d, fclass, mtopp, h)
        lv1, lv4 = bm_r6vol1(d, fclass, xlogs, ld1)
        nlog = Int(floor(xlogs)); nacc = (xlogs - nlog) > 0f0 ? nlog + 1 : nlog
        for k in 1:nacc
            vol2 += bm_anint(lv1[k])
            vol4 += bm_anint(lv4[k] * 10f0) / 10f0
        end
        v
    end
    return (max(v1, 0f0), max(vol4, 0f0), max(vol2, 0f0))
end

function compute_volumes_wc!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species
    ifor = Int(s.plot.forest_idx)                   # forkod JFOR index (1..10) — form-class table + INGY iregn
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 39
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = wc_bratio(sd, sp, d)
        dbhmin = sp == 11 ? 6.0f0 : 7.0f0            # wc/sitset.f westside: LP(sp11)=6, else 7
        # Top-killed trees: full cubic uses NORMAL height (norm_ht); the profile truncates at the break.
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, bf = wc_fw2_westside_vol(eq, d, hv, bark)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = 6)
            tcf = max(v[1], 0f0); mcf = max(v[4] + v[7], 0f0); bf = max(v[2], 0f0)
        else                                          # 616BEHW
            tcf, mcf, bf = wc_behre_vol(sp, ifor, d, hv, bark)
        end
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = d >= dbhmin ? max(bf, 0f0) : 0f0
    end
    return s
end
