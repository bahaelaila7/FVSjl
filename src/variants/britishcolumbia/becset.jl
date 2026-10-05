# =============================================================================
# becset.jl (britishcolumbia) — BEC-zone parsing (habtyp.f) + site-series BAMAX table (sitset.f). CHUNK 2.
#
# Removes the ICHmw2/01 hardcode: the STDINFO habitat field (KARD(1)//KARD(2)) is parsed by HABTYP into
# Region/Zone/SubZone/Series; a Region (CAR/KAM/NEL) is REQUIRED, else the grinit default ICHmw2/01 is kept
# (so all_BC — '231Dd', no Region — stays ICHmw2/01). BAMAX comes from the sitset.f Zone→SubZone→Series table.
# The DG/HTG/mortality coefficient matchers (bc_resolve_*) are already zone/series-general. Validated on the
# all_BC_idf.key stand (STDINFO 'KAMIDFDK1/01' ⇒ IDFdk1/01, BAMAX 60 m²/ha vs ICH's 89).
# =============================================================================

const BC_HAB_RGN = ("CAR", "KAM", "NEL")
const BC_HAB_ZN  = ("ESSF", "ICH", "IDF", "MS", "PP", "SBS", "SBPS")   # habtyp INDEX order (first match wins)
const BC_HAB_SZ  = ("DC","DK","DM","DW","MC","MH","MK","MM","MW","VK","WC","WK","WM","WW","XC","XH","XK","XM")

const BC_M2pHAtoFT2pACR = 4.3560773f0   # METRIC.F77 PARAMETER M2pHAtoFT2pACR (sitset.f:260); NOT 1/FT2pACRtoM2pHA (4.356078, 2 ULP high)

"""Parse a raw BEC/STDINFO string → (Zone, SubZone, Series, PrettyName) or `nothing` (habtyp.f). Region required."""
function bc_habtyp(raw::AbstractString)
    full = uppercase(filter(!isspace, String(raw)))
    any(r -> occursin(r, full), BC_HAB_RGN) || return nothing      # Region required (GOTO 3 → default)
    zone = ""
    for z in BC_HAB_ZN; if occursin(z, full); zone = z; break; end; end
    isempty(zone) && return nothing
    subzone = ""
    for sz in BC_HAB_SZ
        j1 = findfirst(sz, full)
        if j1 !== nothing
            subzone = lowercase(sz)
            js = findfirst('/', full)
            if js !== nothing && js == last(j1) + 2                 # digit between SubZone and '/' (habtyp.f:116)
                subzone *= lowercase(string(full[last(j1) + 1]))
            end
            break
        end
    end
    isempty(subzone) && return nothing
    js = findfirst('/', full); js === nothing && return nothing     # Series after '/'
    series = lowercase(strip(full[(js + 1):min(lastindex(full), js + 5)]))
    return (zone, subzone, series, zone * subzone * "/" * series)
end

"""iSeries (sitset.f:51-59): numeric series, or 1 for special '-'/'a'/'b' suffixes."""
function bc_iseries(series::AbstractString)
    (occursin("-", series) || occursin("a", series) || occursin("b", series)) && return 1
    something(tryparse(Int, strip(series)), 0)
end

"""BC site-series maximum basal area (sitset.f:61-251), Zone→SubZone→Series → m²/ha → ft²/ac. Fallback 10 m²/ha."""
function bc_sitset_bamax(zone::AbstractString, subzone::AbstractString, iser::Integer)
    m = 0f0
    if zone == "IDF"
        if     subzone == "dk1"; m = iser == 1 ? 60f0 : iser == 3 ? 55f0 : (iser in (2,4,5,6,7)) ? 58f0 : 0f0
        elseif subzone == "dk2"; m = (iser == 1 || iser in 3:7) ? 89f0 : 0f0
        elseif subzone == "dk3"; m = 53f0
        elseif subzone == "dk4"; m = 40f0
        elseif subzone == "dm1"; m = (iser == 1 || iser in 3:7) ? 43f0 : 0f0
        elseif subzone == "dm2"; m = iser == 4 ? 52f0 : (iser in (1,3,5,6,7)) ? 53f0 : 0f0
        elseif subzone in ("mw1","mw2"); m = 47f0
        elseif subzone == "xh1"; m = iser == 4 ? 32f0 : 48f0
        elseif subzone == "xh2"; m = iser == 1 ? 42f0 : 48f0
        elseif subzone == "xm";  m = iser in (2,3,4) ? 25f0 : iser == 8 ? 40f0 : (iser in (1,5,6,7)) ? 30f0 : 0f0
        end
    elseif zone == "ICH"
        if     subzone == "mk1"; m = iser == 1 ? 59f0 : iser == 2 ? 50f0 : iser == 3 ? 46f0 : iser == 4 ? 42f0 : (iser in 5:7) ? 50f0 : 0f0
        elseif subzone == "mw2"; m = (iser in 1:7) ? 89f0 : 0f0
        elseif subzone == "mw3"; m = iser == 1 ? 81f0 : (iser in (2,3,6,7)) ? 76f0 : (iser in (4,5)) ? 69f0 : 0f0
        elseif subzone == "wk1"; m = (iser == 1 || iser in 3:7) ? 67f0 : 0f0
        elseif subzone == "dw";  m = 50f0
        end
    elseif zone == "ESSF"
        if     subzone == "dk"; m = iser == 1 ? 64f0 : iser == 2 ? 62f0 : iser == 3 ? 61f0 : (iser in 4:7) ? 62f0 : 0f0
        elseif subzone in ("wc4","wm"); m = 62f0
        end
    elseif zone == "MS"
        if     subzone == "dk";  m = iser == 1 ? 48f0 : (iser in (2,3)) ? 55f0 : iser == 4 ? 38f0 : (iser in 5:7) ? 55f0 : 0f0
        elseif subzone == "dm1"; m = (iser in (1,4)) ? 63f0 : (iser in (2,3,7)) ? 62f0 : iser == 5 ? 57f0 : iser == 6 ? 64f0 : 0f0
        end
    elseif zone == "PP"
        if     subzone == "dh2"; m = (iser == 1 || iser in 3:7) ? 32f0 : 0f0
        elseif subzone == "xh2"; m = (iser in 1:4) ? 32f0 : (iser in (6,7)) ? 21f0 : 0f0
        end
    elseif zone == "SBS"
        m = subzone == "mh" ? 60f0 : 55f0
    elseif zone == "SBPS"
        m = subzone == "mk" ? 55f0 : 50f0
    end
    m <= 0f0 && (m = 10f0)                       # sitset.f:244 undefined fallback
    return m * BC_M2pHAtoFT2pACR                 # m²/ha → ft²/ac
end

"""Resolve the stand's BEC (parse eco_unit; default ICHmw2/01 like grinit). Returns (Zone,SubZone,Series,PrettyName)."""
function bc_becset(s::StandState)
    parsed = bc_habtyp(s.plot.eco_unit)
    parsed === nothing ? ("ICH", "mw2", "01", "ICHmw2/01") : parsed
end

# canada/bc blkdat.f:120-127 DATA HT1 / HT2 — the Wykoff HT-DBH intercept/slope CRATET dubs missing heights with (the
# species CSV ht1/ht2 carry other values for the hardwoods EP/AT/AC/OH, and wykoff_ht2 is a different table).
const BC_BLK_HT1 = Float32[5.19988, 4.97407, 4.81519, 5.00233, 4.97331, 4.89564, 4.62171, 4.92190, 4.76537, 4.92880,
                           5.152, 5.152, 5.152, 4.81519, 5.152]
const BC_BLK_HT2 = Float32[-9.26718, -6.78347, -7.29306, -8.19365, -8.19730, -8.39057, -5.32481, -8.30289, -7.61062,
                           -9.32795, -13.576, -13.576, -13.576, -7.29306, -13.576]

"""
    bc_ht_coefs(s) -> (ht1, ht2, lmhtdub)

canada/bc htgf.f ENTRY HTCONS step 2 (:1817-2015, "REFIT OF HEIGHT DUBBING MODEL FOR VERSION 3"): RCON calls it at
LSTART (cratet.f), AFTER sitset's BECSET — it restarts from the compile-time HT1/HT2 (XHT1B) and, for the ICH/IDF/SBS/SBPS
zones and the SBSdw1/SBSmh/SBSdw2 units, replaces a species' HT1/HT2 with a METRIC Wykoff fit (C0,C1 on ln(m−1.3) vs
cm), flagged LMHTDUB. (becset.f's V2 ESSF/MS/PP table is overwritten by this reset before CRATET dubs.) Other species
and zones keep blkdat HT1/HT2 (imperial form).
"""
function bc_ht_coefs(s::StandState)
    zone, sz, _, pretty = bc_becset(s)
    ht1 = copy(BC_BLK_HT1); ht2 = copy(BC_BLK_HT2)
    z(x...) = any(y -> occursin(y, zone), x); hs(x...) = any(y -> occursin(y, sz), x); pn(x...) = any(y -> occursin(y, pretty), x)
    set!(i, a, b) = (ht1[i] = a; ht2[i] = b)
    for i in 1:15
        if i == 1
            if z("ICH")
                if hs("mw");      set!(i, 3.9738f0, -20.8459f0)
                elseif hs("wk");  set!(i, 3.7771f0, -15.8162f0)
                end
            end
        elseif i == 2
            if z("ICH")
                if hs("dw", "mk");  set!(i, 3.9045f0, -22.0149f0)
                elseif hs("mw");    set!(i, 3.7825f0, -14.5352f0)
                end
            elseif z("IDF")
                hs("dm", "xh") && set!(i, 3.9012f0, -19.1924f0)
            end
        elseif i == 3 || i == 14
            if z("ICH") || pn("SBSdw1", "SBSmh")
                if hs("dk", "mk") || pn("SBSdw1", "SBSmh"); set!(i, 3.9432f0, -25.9091f0)
                elseif hs("mw");  set!(i, 3.9274f0, -20.1322f0)
                elseif hs("wk");  set!(i, 4.1218f0, -29.6135f0)
                elseif hs("dw");  set!(i, 3.7859f0, -22.0166f0)
                end
            elseif z("IDF", "SBPS") || pn("SBSdw2")
                if hs("dk", "mw", "xm") || z("SBPS") || pn("SBSdw2"); set!(i, 3.7049f0, -22.9781f0)
                elseif hs("dm", "ww");  set!(i, 3.7877f0, -23.1796f0)
                elseif hs("xh");        set!(i, 3.6337f0, -18.2998f0)
                end
            end
        elseif i == 5
            if z("ICH")
                if hs("vc", "vk");      set!(i, 3.9685f0, -30.4439f0)
                elseif hs("mw", "wk");  set!(i, 3.9491f0, -26.1417f0)
                elseif hs("mc");        set!(i, 3.8170f0, -21.3197f0)
                end
            elseif z("IDF")
                if hs("mm");                  set!(i, 3.2884f0, -12.9475f0)
                elseif hs("wv", "ww", "wc");  set!(i, 4.0577f0, -28.1854f0)
                end
            end
        elseif i == 6
            if z("ICH")
                if hs("mw");                  set!(i, 3.8565f0, -25.7777f0)
                elseif hs("mm", "vk");        set!(i, 4.1251f0, -41.8833f0)
                elseif hs("wk", "mc", "mk");  set!(i, 3.9432f0, -30.5647f0)
                end
            elseif z("IDF")
                if hs("dm");      set!(i, 3.2692f0, -10.1899f0)
                elseif hs("ww");  set!(i, 4.0537f0, -36.7598f0)
                end
            end
        elseif i == 7
            if z("ICH") || pn("SBSdw1", "SBSmh")
                (hs("mk", "mw", "mc") || pn("SBSdw1", "SBSmh")) && set!(i, 3.8470f0, -18.5836f0)
            elseif z("IDF", "SBPS") || pn("SBSdw2")
                if hs("xh", "dm");                            set!(i, 3.5683f0, -13.0736f0)
                elseif hs("dk") || z("SBPS") || pn("SBSdw2"); set!(i, 3.1192f0, -8.2208f0)
                end
            end
        elseif i == 8
            if z("ICH", "IDF", "SBPS", "SBS")
                if hs("dk") || z("SBPS") || pn("SBSdw2");          set!(i, 3.6519f0, -17.8872f0)
                elseif hs("dm", "mk", "mc") || z("SBS");           set!(i, 4.0115f0, -25.0730f0)
                elseif hs("mw");                                   set!(i, 3.9729f0, -22.6984f0)
                elseif hs("wk");                                   set!(i, 4.1277f0, -29.2780f0)
                end
            end
        elseif i == 9 || i == 4
            if z("ICH", "IDF", "SBS")
                if hs("mk", "vk", "dk") || z("SBS");  set!(i, 3.7656f0, -20.1044f0)
                elseif hs("mw", "wk", "mc");          set!(i, 3.9598f0, -24.0969f0)
                end
            end
        elseif i == 10
            if z("ICH")
                if hs("dw");      set!(i, 4.0178f0, -27.1633f0)
                elseif hs("mk");  set!(i, 3.7798f0, -29.9875f0)
                end
            elseif z("IDF")
                if hs("xh", "dm");  set!(i, 3.7948f0, -26.8113f0)
                elseif hs("dk");    set!(i, 3.8113f0, -29.7283f0)
                end
            end
        end
    end
    lmhtdub = Bool[(ht1[i] != BC_BLK_HT1[i]) || (ht2[i] != BC_BLK_HT2[i]) for i in 1:15]
    return ht1, ht2, lmhtdub
end
