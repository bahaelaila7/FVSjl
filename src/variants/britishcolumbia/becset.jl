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

const BC_M2pHAtoFT2pACR = 1f0 / 0.2295643f0   # inverse of FT2pACRtoM2pHA (sitset.f)

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
