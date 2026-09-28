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

# --- WC VOLEQDEF (wc/sitset.f:218-243 → voleqdef.f R6_EQN westside branch for REGN 6 / R7_EQN for the BLM
#     REGN 7 forests; VAR='WC', DIST='  ', KODFOR = wc/forkod.f's JFOR(IFOR)). The complete WC input domain is
#     the 10 JFOR forests × 39 species, so this is the full decision table — read off FVSwc_g16's "NATIONAL
#     VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS" block for each forest (cubic == board for every entry). jl had
#     only Willamette (618) and sent every other forest's DF/WH/GF/… to region-6 Behre: +25-45% cubic on
#     603/605/615 FIA stands (DF F03/F08/F00FW2W202, INGY I11-I13 subregions). ---
const _WC_VOLEQ_BY_FOREST = Dict{Int,NTuple{39,String}}(
    603 => ("I12FW2W017", "616BEHW015", "616BEHW017", "I00FW2W108", "616BEHW020", "616BEHW000", "616BEHW022", "616BEHW042", "616BEHW081", "616BEHW093", "616BEHW108", "616BEHW116", "616BEHW117", "616BEHW119", "616BEHW122", "F03FW2W202", "616BEHW211", "616BEHW242", "F00FW2W263", "616BEHW264", "616BEHW312", "A16CURW351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "616BEHW103", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    605 => ("616BEHW011", "616BEHW015", "616BEHW017", "616BEHW019", "616BEHW020", "616BEHW000", "616BEHW022", "616BEHW042", "616BEHW081", "616BEHW093", "616BEHW108", "616BEHW116", "616BEHW117", "616BEHW119", "616BEHW122", "F08FW2W202", "616BEHW211", "616BEHW242", "F03FW2W263", "616BEHW264", "616BEHW312", "616BEHW351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "616BEHW103", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    606 => ("I12FW2W017", "616BEHW015", "I13FW2W017", "616BEHW019", "616BEHW020", "616BEHW000", "I13FW2W017", "616BEHW042", "616BEHW081", "I11FW2W093", "I11FW2W108", "616BEHW116", "616BEHW117", "616BEHW119", "I12FW2W122", "F03FW2W202", "616BEHW211", "616BEHW242", "I11FW2W260", "616BEHW264", "616BEHW312", "616BEHW351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "616BEHW103", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    610 => ("616BEHW011", "I00FW2W093", "616BEHW017", "616BEHW019", "616BEHW020", "616BEHW000", "616BEHW022", "616BEHW042", "616BEHW081", "616BEHW093", "616BEHW108", "616BEHW116", "616BEHW117", "616BEHW119", "I00FW2W073", "F06FW2W202", "616BEHW211", "616BEHW242", "F06FW2W263", "616BEHW264", "616BEHW312", "616BEHW351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "616BEHW103", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    615 => ("I13FW2W017", "I00FW2W017", "616BEHW017", "616BEHW019", "I00FW2W012", "616BEHW000", "616BEHW022", "616BEHW042", "I00FW2W073", "I00FW2W093", "I00FW2W108", "616BEHW116", "616BEHW117", "616BEHW119", "I00FW2W073", "F00FW2W202", "616BEHW211", "I00FW2W012", "I11FW2W260", "I00FW2W242", "616BEHW312", "616BEHW351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "I00FW2W260", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    618 => ("616BEHW011", "616BEHW015", "I00FW2W017", "616BEHW019", "616BEHW020", "616BEHW000", "I00FW2W108", "616BEHW042", "I00FW2W073", "616BEHW093", "616BEHW108", "616BEHW116", "616BEHW117", "616BEHW119", "616BEHW122", "F05FW2W202", "616BEHW211", "616BEHW242", "F03FW2W263", "616BEHW264", "616BEHW312", "616BEHW351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "616BEHW103", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    708 => ("B00BEHW011", "B00BEHW015", "B00BEHW017", "B00BEHW015", "B00BEHW021", "B00BEHW999", "B00BEHW022", "B00BEHW042", "B00BEHW081", "B00BEHW093", "B00BEHW108", "B00BEHW116", "B00BEHW117", "B00BEHW119", "B00BEHW122", "B01BEHW202", "B00BEHW211", "B00BEHW242", "B00BEHW260", "B00BEHW260", "B00BEHW312", "B00BEHW351", "B00BEHW361", "B00BEHW999", "B00BEHW431", "B00BEHW999", "B00BEHW747", "B00BEHW800", "B00BEHW242", "B00BEHW073", "B00BEHW119", "B00BEHW108", "B00BEHW231", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999"),
    709 => ("B00BEHW011", "B00BEHW015", "B00BEHW017", "B00BEHW015", "B00BEHW021", "B00BEHW999", "B00BEHW022", "B00BEHW042", "B00BEHW081", "B00BEHW093", "B00BEHW108", "B00BEHW116", "B00BEHW117", "B00BEHW119", "B00BEHW122", "B01BEHW202", "B00BEHW211", "B00BEHW242", "B00BEHW260", "B00BEHW260", "B00BEHW312", "B00BEHW351", "B00BEHW361", "B00BEHW999", "B00BEHW431", "B00BEHW999", "B00BEHW747", "B00BEHW800", "B00BEHW242", "B00BEHW073", "B00BEHW119", "B00BEHW108", "B00BEHW231", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999"),
    710 => ("B00BEHW011", "B00BEHW015", "B00BEHW017", "B00BEHW015", "B00BEHW021", "B00BEHW999", "B00BEHW022", "B00BEHW042", "B00BEHW081", "B00BEHW093", "B00BEHW108", "B00BEHW116", "B00BEHW117", "B00BEHW119", "B00BEHW122", "B01BEHW202", "B00BEHW211", "B00BEHW242", "B00BEHW260", "B00BEHW260", "B00BEHW312", "B00BEHW351", "B00BEHW361", "B00BEHW999", "B00BEHW431", "B00BEHW999", "B00BEHW747", "B00BEHW800", "B00BEHW242", "B00BEHW073", "B00BEHW119", "B00BEHW108", "B00BEHW231", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999"),
    711 => ("B00BEHW011", "B00BEHW015", "B00BEHW017", "B00BEHW015", "B00BEHW021", "B00BEHW999", "B00BEHW022", "B00BEHW042", "B00BEHW081", "B00BEHW093", "B00BEHW108", "B00BEHW116", "B00BEHW117", "B00BEHW119", "B00BEHW122", "B01BEHW202", "B00BEHW211", "B00BEHW242", "B00BEHW260", "B00BEHW260", "B00BEHW312", "B00BEHW351", "B00BEHW361", "B00BEHW999", "B00BEHW431", "B00BEHW999", "B00BEHW747", "B00BEHW800", "B00BEHW242", "B00BEHW073", "B00BEHW119", "B00BEHW108", "B00BEHW231", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999"),
)

"VOLEQ for WC species `sp` on forest `kodfor` (post-FORKOD JFOR code); unknown codes use the grinit 618 row."
_wc_voleq(kodfor::Int, sp::Int)::String =
    (1 <= sp <= 39) ? get(_WC_VOLEQ_BY_FOREST, kodfor, _WC_VOLEQ_BY_FOREST[618])[sp] : "           "

# ---------------------------------------------------------------------------
# Westside Flewelling shape (f_west.f SHP_W3/W4/W5) + breast-height bark (FDBT_C1).
# Fills RFLW(6)/RHFW(4) consumed UNCHANGED by _fw2_sf_taper / _fw2_sf_yhat / _fw2_tcubic. Double-precision
# kernel (REAL*4 out), matching the proven cr_fw2 pattern. GEOSUB '01'..'08' → regional coefficient index.
# ---------------------------------------------------------------------------
@inline function _wc_geoidx(geosub::AbstractString)::Int
    (length(geosub) == 2 && geosub[1] == '0' && '1' <= geosub[2] <= '8') ? Int(geosub[2] - '0') : 0
end
@inline _logit(u::Float64) = exp(u) / (1.0 + exp(u))

# f_west.f precision (SHP_W3/W4/W5 are IMPLICIT DOUBLE PRECISION with REAL*4 DBHOB/HTTOT/outputs): the F(·) DATA are
# REAL*8 D-literals, but the M1-M3 / R25 / R34 DATA are REAL*4, the DMEDIAN power is evaluated in SINGLE precision
# (REAL operands ⇒ powf), LOG of a REAL argument is logf, LOG/EXP of a REAL*8 argument are glibc log/exp, and each
# logistic is written `DEXP(U)/(1.0D0+DEXP(U))` (a coefficient in front multiplies BEFORE the division). RHI1/RHI2/
# RHLONGI and the RFLW/RHFW outputs are REAL*4. Every piece below follows that statement order.
# (_r4 = REAL*4 DATA widened to REAL*8 — shared helper, data/centralrockies/fw2/ingy_coefs.jl)
@inline _lf(x::Float32) = Float64(flog(x))                     # LOG(REAL) → logf, then widened
@inline _lgt(u::Float64) = dexp(u) / (1.0 + dexp(u))
@inline _clp(u::Float64, lo::Float64, hi::Float64) = u < lo ? lo : (u > hi ? hi : u)

const _W3_R25 = (-2.0262, -1.7945, -2.0366, -2.0811, -1.9868, -2.0151, -1.9475, -2.0151)   # REAL*4 DATA
const _W3_R34 = ( 5.2132,  5.1417,  5.1768,  5.1156,  5.1654,  5.2413,  5.1427,  5.2413)   # REAL*4 DATA

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
        id != 0 && (fr25 = _r4(_W3_R25[id]); fr34 = _r4(_W3_R34[id]))
    end
    # DMEDIAN = m1*(HTTOT-4.5)**(m2+m3*HTTOT), m1..m3 REAL*4 ⇒ all single
    dmedian = Float64(0.566f0 * fpow(h - 4.5f0, 0.634f0 + 0.00074f0 * h))
    dform = D / dmedian - 1.0
    u7 = F13 + F14 * _lf(d + 1f0) + F15 * _lf(h)
    u9a = _clp(F18 + F19 * H, -7.0, 7.0)
    u9 = F17 * dexp(u9a) / (1.0 + dexp(u9a))
    u8 = F21 + F22 * _lf(h) + F23 * dform + F24 * dpow(D / 10.0, 1.5)
    u1 = fr25 + F26 * _lf(h)
    u2 = F29 + F30 * _lf(d) + F31 * D
    u3 = fr34 + F35 * dlog(D + 1.0) + F36 * _lf(h) + F37 * dform * H + F33 * dform
    u4 = F38 + F39 * D + F40 * H * D
    u5 = F42 + F43 * H + F44 * dform
    u1 = _clp(u1, -7.0, 7.0); u2 = _clp(u2, -7.0, 7.0); u3 = _clp(u3, -7.0, 7.0); u4 = _clp(u4, -7.0, 7.0)
    u5 = _clp(u5, -7.0, 7.0); u7 = _clp(u7, -7.0, 7.0); u8 = _clp(u8, -7.0, 7.0)
    u6 = _clp(F45 + F46 * dlog(D + 1.0) + F47 * _lf(h), -6.0, 6.0)
    u6 = 1.0 + dexp(u6)
    u6 < Float64(1.005f0) && (u6 = 1.005); u6 > 100.0 && (u6 = 100.0)       # `U6 .LT. 1.005` (REAL literal)
    rhi1 = Float32(_lgt(u7)); rhi1 > 0.5f0 && (rhi1 = 0.5f0)
    _wc_assemble_shp_rhi(u1, u2, u3, u4, u5, u6, rhi1, u8, u9)
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
    F45=-1.3248733; F46=-0.11788962; F47=-0.015909154                    # F(48) has no DATA ⇒ 0 (static storage)
    R25 = (-7.511,-7.687,-7.224,-7.355,-7.632,-7.646,-7.687,-7.911)        # REAL*4 DATA
    R34 = (-1.215,-1.355,-1.177,-1.373,-1.398,-1.188,-1.355,-1.449)
    fr25 = F25; fr34 = F34
    if geosub != "00"
        id = _wc_geoidx(geosub); id != 0 && (fr25 = _r4(R25[id]); fr34 = _r4(R34[id]))
    end
    lh = _lf(h)
    dmedian = Float64(0.2855f0 * fpow(h - 4.5f0, 0.307f0 - 0.00505f0 * h + 0.00001745f0 * h * h + 0.19f0 * flog(h)))
    dform = D / dmedian - 1.0
    u7 = _clp(F13 + F14 * (1.0 - dexp(F15 * H)), -7.0, 7.0)
    rhi1 = Float32(_lgt(u7)); rhi1 > 0.5f0 && (rhi1 = 0.5f0)
    u9a = F17 + F18 * lh + F19 * dform
    (Float64(rhi1) + u9a > 0.75) && (u9a = Float64(0.75f0 - rhi1))
    u9 = u9a > 0.0 ? u9a : 0.0
    u8 = F21 + F22 * lh + F23 * dform + F24 * lh * dform
    hh = H / 100.0
    u1 = fr25 + F26 / hh + F27 / (hh * hh) + F28 / (hh * (hh * hh))
    u2 = F29
    u3 = fr34 + F35 / H + F36 * H * D / 1000.0 + F37 * dform
    u4 = F38 + F39 * lh + F40 * dform
    u5 = F42 + F43 * lh + F44 * D
    u1 = _clp(u1, -7.0, 7.0); u2 = _clp(u2, -7.0, 7.0); u3 = _clp(u3, -7.0, 7.0); u4 = _clp(u4, -7.0, 7.0)
    u5 = _clp(u5, -7.0, 7.0); u8 = _clp(u8, -7.0, 7.0)
    u6 = _clp(F45 + F46 * _lf(d) + F47 * H + 0.0 * lh, -6.0, 6.0)
    u6 = 1.0 + dexp(u6)
    u6 < 1.005 && (u6 = 1.005); u6 > 100.0 && (u6 = 100.0)                 # W4 compares 1.005D0
    _wc_assemble_shp_rhi(u1, u2, u3, u4, u5, u6, rhi1, u8, u9)
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
    dmedian = Float64(0.11f0 * fpow(h - 4.5f0, 1.08f0 + 0.0006f0 * h))
    dform = D / dmedian - 1.0
    u7 = _clp(F13 + F15 * dform, -7.0, 1.0)                                  # W5: upper clamp 1.0, NO 0.5 cap
    rhi1 = Float32(_lgt(u7))
    u9a = F17 + F18 * _lf(h)
    (Float64(rhi1) + u9a > 0.75) && (u9a = Float64(0.75f0 - rhi1))
    u9 = u9a > 0.0 ? u9a : 0.0
    u8 = F21 + F22 / D
    u1 = F25 + F26 * _lf(d) + F27 / D
    u2 = F29 + F30 / D + F31 / Float64(d * d)                                # DBHOB*DBHOB is a REAL product
    u3 = F34 + F35 * D + F36 * H
    u4 = F38 + F39 * D + F40 * dform
    u5 = F42 + F43 / D
    u1 = _clp(u1, -7.0, 7.0); u2 = _clp(u2, -7.0, 7.0); u3 = _clp(u3, -7.0, 7.0); u4 = _clp(u4, -7.0, 7.0)
    u5 = _clp(u5, -7.0, 7.0); u8 = _clp(u8, -7.0, 7.0)
    u6 = _clp(F45 + F46 * D + F47 / D, -6.0, 6.0)
    u6 = 1.0 + dexp(u6)
    u6 < 1.005 && (u6 = 1.005); u6 > 100.0 && (u6 = 100.0)
    _wc_assemble_shp_rhi(u1, u2, u3, u4, u5, u6, rhi1, u8, u9)
end

# Shared tail (identical in W3/W4/W5): R1-R4 = DEXP(U)/(1+DEXP(U)), R5 = 0.5D0+0.5D0*DEXP(U5)/(1+DEXP(U5)),
# A3=U6, RHLONGI=U9 (REAL*4), RHI2=RHI1+RHLONGI (single), RHC = RHI2 + (1.0-RHI2)*DEXP(U8)/(1+DEXP(U8)).
@inline function _wc_assemble_shp_rhi(u1, u2, u3, u4, u5, u6, rhi1::Float32, u8, u9)
    r1 = Float32(_lgt(u1)); r2 = Float32(_lgt(u2)); r3 = Float32(_lgt(u3)); r4 = Float32(_lgt(u4))
    r5 = Float32(0.5 + 0.5 * dexp(u5) / (1.0 + dexp(u5))); a3 = Float32(u6)
    rhlongi = Float32(u9); rhi2 = rhi1 + rhlongi
    rhc = Float32(Float64(rhi2) + Float64(1f0 - rhi2) * dexp(u8) / (1.0 + dexp(u8)))
    ((r1, r2, r3, r4, r5, a3), (rhi1, rhi2, rhc, rhlongi))
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
function _wc_fw2_merch_cuft(dibat, h::Float32, mtop::Float32, stump::Float32; hs = nothing)::Float32
    hs === nothing && (hs = _fw2_hs(dibat, mtop, h)); lmerch = hs - stump
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
function _wc_fw2_board(dibat, h::Float32, bftop::Float32, stump::Float32; hs = nothing)::Float32
    hs === nothing && (hs = _fw2_hs(dibat, bftop, h)); lmerch = hs - stump
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
                            topd::Float32 = 4.5f0, bftopd::Float32 = 4.5f0, stump::Float32 = 1f0, ht2td = nothing)
    # `ht2td` (optional 2-slot buffer) ← [HT1PRD of the cubic call, HT1PRD of the board call]
    ht2td === nothing || (ht2td[1] = 0f0; ht2td[2] = 0f0)
    (d < 1f0 || h <= 5f0) && return (0f0, 0f0, 0f0)
    spec = voleq[8:10]
    jsp = (spec == "202" || spec == "205" || spec == "204") ? 3 : spec == "263" ? 4 : spec == "242" ? 5 : 0
    jsp == 0 && return (0f0, 0f0, 0f0)
    geosub = length(voleq) >= 3 ? voleq[2:3] : "00"
    rflw, rhfw = jsp == 3 ? _fw2_shp_w3(d, h, geosub) :
                 jsp == 4 ? _fw2_shp_w4(d, h, geosub) : _fw2_shp_w5(d, h, geosub)
    tapcoe = _fw2_sf_taper(rhfw, rflw)
    # fvsvol.f:153 DBTBH = D*(1-BARK); sf_shp.f DBHIB = DBHOB-DBTBH (NOT D*BARK — differs in the last bit);
    # sf_2pt.f F = DBH_IB/SF_YHAT(4.5/TOTALH) with SF_YHAT at sf_yhat.f's own REAL/REAL*8 mix, and TCUBIC/MERLEN
    # read the profile through SF_DS → the same SF_YHAT (BRK_UP leaves DIB alone for JSP 3-5).
    dbhib = d - d * (1f0 - bark)
    yhat_bh = _fw2_sf_yhat_f(4.5f0 / h, tapcoe, rhfw, rflw, 1.0f0, h, false)[1]
    yhat_bh == 0f0 && return (0f0, 0f0, 0f0)
    f = dbhib / yhat_bh
    dibat = ht -> (Float32(ht) > h ? 0f0 : _fw2_sf_yhat_f(Float32(ht) / h, tapcoe, rhfw, rflw, f, h, false)[1])
    stump_dib = h <= 15f0 ? _fw2_fwsmall(jsp, h, dibat(1.0f0), dbhib) : -1f0
    v1 = _nint(_fw2_tcubic(dibat, h; stump_dib = stump_dib) * 10.0f0) * 1f-1   # profile.f:293 VOL(1)=NINT(TCVOL*10.0)*1E-1
    # profile.f MERLEN (VOLEQ(4:4)='F'): the merch-top height is SF_HS itself (the Newton/bisection solver, not
    # a diameter-tolerance bisection — they straddle SEGMNT even-foot boundaries); HT1PRD = LMERCH + STUMP with
    # LMERCH = MAX(HS − STUMP, 0) (profile.f:342, 1052) → fvsvol.f HT2TD.
    hsc = _fw2_sf_hs(tapcoe, rhfw, rflw, f, h, topd * bark)
    hsb = bftopd == topd ? hsc : _fw2_sf_hs(tapcoe, rhfw, rflw, f, h, bftopd * bark)
    v4 = _wc_fw2_merch_cuft(dibat, h, topd * bark, stump; hs = hsc)
    v2 = _wc_fw2_board(dibat, h, bftopd * bark, stump; hs = hsb)
    if ht2td !== nothing
        lc = hsc - stump; lc < 0f0 && (lc = 0f0); ht2td[1] = lc + stump
        lb = hsb - stump; lb < 0f0 && (lb = 0f0); ht2td[2] = lb + stump
    end
    return (max(v1, 0f0), max(v4, 0f0), max(v2, 0f0))
end

# r6vol.f:118-124: after R6DIBS, HT1PRD = 1.0 + Σ XLEN(1:20) — only reached when TTH > FC_HT (17.3) and
# DBHIB ≥ MTOPP (the short-tree cylinder and the total-cubic-only exits leave the caller's HT1PRD = 0).
function r6vol_ht1prd(d::Float32, fclass, mtopp::Float32, h::Float32, dbhib::Float32)::Float32
    (h <= 17.3f0 || dbhib < mtopp) && return 0f0
    _, _, xl = bm_r6dibs(d, fclass, mtopp, h)
    ht1 = 1.0f0
    for k in 1:20; ht1 += xl[k]; end
    return ht1
end

# --- Behre (616BEHW) per-tree volume — reuse the BM R6 machinery + WC form class. ---
function wc_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32; topd::Float32 = 4.5f0,
                      ht2td = nothing)
    fclass = wc_formcl(sp, ifor, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    if ht2td !== nothing                               # cubic and board calls share MTOPP = TOPD·BARK (BFTOPD = TOPD)
        ht2td[1] = r6vol_ht1prd(d, fclass, topd * bark, h, dbhib); ht2td[2] = ht2td[1]
    end
    vol2 = 0f0; vol4 = 0f0
    v1 = if h <= 17.3f0
        0.00272708f0 * dbhib * dbhib * h            # R6VOL short-tree cylinder (R6DIBS/R6VOL1 skipped)
    else
        v = bm_r6vol3(d, dbtbh, fclass, h, 1)
        mtopp = topd * bark                          # TOPDIAM = TOPD·BARK
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

# --- A16CURW351: red alder on Gifford Pinchot (603). volinit.f:454-472 sends MDL='CUR' to R10VOL only for REGN 10
#     (or A01/A02), so on REGN 6 it is PROFILE with the Region-10 red-alder taper R10TAP (DVREDA, shared with AK:
#     `_ak_cur_dib`), the merch length from R10HTS (profile.f:225, `_ak_cur_lmerch`; MERLEN only if that is 0,
#     profile.f:322-339), and MRULES region 6 (mrules.f:302: OPT 23, EVOD 2, MAXLEN 16, MINLEN 2, TRIM 0.5,
#     MERCHL 8, COR 'N'). TCUBIC 4-ft Smalian, VOL(1)=NINT(TCVOL·10)·0.1 (profile.f:293); GETDIB inch classes
#     (INT, +1 above .501), top log floored at MTOPP, DEM/CUR stump D/36 for D>36 (profile.f:1143).

"MERLEN (profile.f:978) non-Flewelling branch: 0.1-ft bisection to the height where INT((DIB+.005)·10) ≥ TOP1."
function _cur_merlen(dibat, h::Float32, top::Float32, stump::Float32)::Float32
    top1 = Float32(round((top + 0.005f0) * 10.0f0, RoundNearestTiesAway))      # ANINT((DS+.005)*10)
    first = 1; last = trunc(Int, h + 0.5f0) * 10
    for _ in 1:last
        first == last && break
        half = (first + last + 1) ÷ 2
        dib = Float32(trunc(Int, (dibat(Float32(half) / 10.0f0) + 0.005f0) * 10.0f0))
        top1 <= dib ? (first = half) : (last = half - 1)
    end
    lm = Float32(first) / 10.0f0 - stump
    return lm < 0f0 ? 0f0 : lm
end

"PROFILE log products for A16CURW351: merch cubic VOL(4) (`board=false`) or Scribner VOL(2) (`board=true`)."
function _cur_profile_logs(d::Float32, h::Float32, dibat, mtop::Float32, stump::Float32, board::Bool)::Float32
    lmerch = _ak_cur_lmerch(d, h, mtop, stump)
    lmerch <= 0f0 && (lmerch = _cur_merlen(dibat, h, mtop, stump))
    lmerch < 8f0 && return 0f0                                        # MERCHL
    lmerch / (16f0 + 0.5f0) > 20f0 && return 0f0                      # ERRFLAG 12
    numseg = _nvb_numlog(23, 2, lmerch, 16f0, 2f0, 0.5f0)
    numseg == 0 && return 0f0
    loglen, numseg = _nvb_segmnt(23, 2, lmerch, 16f0, 2f0, 0.5f0, numseg)
    raw = zeros(Float32, numseg + 1)
    raw[1] = dibat(4.5f0)                                             # STUMP ≤ 4.5: butt large end at DBH
    ht2 = d > 36f0 ? d / 36f0 : stump
    @inbounds for i in 1:numseg
        ht2 += 0.5f0 + loglen[i]
        raw[i + 1] = dibat(ht2)
    end
    raw[numseg + 1] < mtop && (raw[numseg + 1] = mtop)
    v = 0f0
    if board
        @inbounds for i in 1:numseg
            v += Float32(round(_scrib(_fw2_dclass(raw[i + 1]), loglen[i], 'N'), RoundNearestTiesAway))
        end
    else
        dibl = _fw2_dclass(raw[1])
        @inbounds for i in 1:numseg
            dibs = _fw2_dclass(raw[i + 1])
            logv = 0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen[i]
            v += Float32(round(logv * 10.0f0, RoundNearestTiesAway)) / 10.0f0
            dibl = dibs
        end
    end
    return v
end

"NATCRS for A16CURW351: cubic call (MTOPP=TOPD·BARK → VOL(1), VOL(4)); board call for D ≥ BFMIND (BFTOPD·BARK)."
function wc_cur_vol(d::Float32, h::Float32, bark::Float32; topd::Float32, bftopd::Float32, stump::Float32,
                    bfmind::Float32, ht2td = nothing)
    ht2td === nothing || (ht2td[1] = 0f0; ht2td[2] = 0f0)
    (d < 1f0 || h < 5f0) && return (0f0, 0f0, 0f0)                  # profile.f:117
    if ht2td !== nothing
        # profile.f HT1PRD of the R10HTS/MERLEN merch length (see _cur_ht1prd, shared with AK).
        ht2td[1] = _cur_ht1prd(d, h, topd * bark, stump)
        d >= bfmind && (ht2td[2] = _cur_ht1prd(d, h, bftopd * bark, stump))
    end
    dibat = ht -> _ak_cur_dib(d, h, Float32(ht))
    tcf = _nint(_fw2_tcubic(dibat, h) * 10.0f0) * 1f-1
    mcf = _cur_profile_logs(d, h, dibat, topd * bark, stump, false)
    bf = d >= bfmind ? _cur_profile_logs(d, h, dibat, bftopd * bark, stump, true) : 0f0
    return (tcf, mcf, bf)
end

function compute_volumes_wc!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species; c = s.control
    ifor = Int(s.plot.forest_idx)                   # forkod JFOR index (1..10) — form-class table + INGY iregn
    # vols.f:86-90 zeroes HT2TD for every record; NATCRS (fvsvol.f:337-339 / 484-487) stores the NVEL HT1PRD.
    fill!(t.merch_top_cf, 0f0); fill!(t.merch_top_bf, 0f0)
    htb = zeros(Float32, 2)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 39
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        htb[1] = 0f0; htb[2] = 0f0
        # wc/vols.f:150-151: BARK=BRATIO(ISPC,D_start,H) before `D=D+DG(I)/BARK` ⇒ projected cycles use the stashed
        # start-of-cycle bark (t.vol_bark) for the merch tops / DBHIB / CFTOPK; grown-DBH bark at cycle 0 / dead records.
        bark = (i <= t.n && t.vol_bark[i] > 0f0) ? t.vol_bark[i] : wc_bratio(sd, sp, d)
        # wc/sitset.f:192-213 merch specs (init_merch_standards!): BLM IFOR 7-10 TOPD=BFTOPD=5, DBHMIN=7 for all;
        # otherwise TOPD=BFTOPD=4.5, DBHMIN=BFMIND=7 (LP sp11 = 6).
        dbhmin = c.sp_dbh_min[sp]; bfmind = c.sp_bf_dbhmin[sp]
        topd = c.sp_top_diam[sp]; stmp = c.sp_stump_ht[sp]
        bftopd = c.sp_bf_topd[sp]; bfstmp = c.sp_bf_stump[sp]
        # wc/vols.f:145-146: a top-killed tree (H≥4.5, ITRUNC>0) is volumed at its NORMAL height NORMHT …
        tkill = h >= 4.5f0 && t.trunc[i] > 0
        hv = tkill ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, bf = wc_fw2_westside_vol(eq, d, hv, bark; topd = topd, bftopd = bftopd, stump = stmp, ht2td = htb)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = topd, bftopd = bftopd, stump = stmp, iregn = 6, board_cor = 'N',
                           merch_opt = 23, sf_hs = true, ht2td = htb)
            tcf = max(v[1], 0f0); mcf = max(v[4] + v[7], 0f0); bf = max(v[2], 0f0)
        elseif se[1] == 'B' || se[1] == 'b'           # BLM forests (R7_EQN B00BEHW/B01BEHW202) → NVEL BLMVOL
            tcf, mcf, bf = _blm_natcrs(eq, wc_formcl(sp, ifor, d), d, hv, bark, topd, bftopd, bfmind; ht2td = htb)
        elseif mdl == "CUR"                           # A16CURW351 (603 RA) → PROFILE + R10TAP red alder
            tcf, mcf, bf = wc_cur_vol(d, hv, bark; topd = topd, bftopd = bftopd, stump = stmp, bfmind = bfmind, ht2td = htb)
        else                                          # 616BEHW
            tcf, mcf, bf = wc_behre_vol(sp, ifor, d, hv, bark; topd = topd, ht2td = htb)
        end
        d >= dbhmin && (t.merch_top_cf[i] = htb[1])  # fvsvol.f:337-339 HT2TD(IT,2) = HT1PRD (HT2PRD = 0 in R6)
        d >= bfmind && (t.merch_top_bf[i] = htb[2])  # fvsvol.f:484-487 board call HT2TD(IT,1) = HT1PRD
        tcf = max(tcf, 0f0)
        mcf = d >= dbhmin ? max(mcf, 0f0) : 0f0      # fvsvol.f:513 MCF only for D≥DBHMIN
        bf  = d >= bfmind ? max(bf, 0f0) : 0f0       # fvsvol.f:517 BBFV=0 for D<BFMIND
        # … then trimmed to the standing broken stem by CFTOPK/BFTOPK (vols.f:191-193, 390-391; NATCRS sets
        # CTKFLG=BTKFLG=.TRUE., VMAX=TCF, BFMAX=VMAX). jl volumed the full NORMHT stem and never truncated.
        if tkill && tcf > 0f0
            vmax = tcf
            tcf, mcf = cr_cftopk(tcf, mcf, d, hv, vmax, bark, Int(t.trunc[i]), stmp, topd)
            bf = cr_bftopk(bf, d, hv, vmax, bark, Int(t.trunc[i]), bfstmp, bftopd)
        end
        t.cuft_vol[i] = tcf
        t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = bf
    end
    return s
end
