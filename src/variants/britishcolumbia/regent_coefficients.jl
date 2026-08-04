# =============================================================================
# regent_coefficients.jl (britishcolumbia) — BC small-tree growth coefficients + V3 HTG regression.
#
# CHUNK 6. `canada/bc/regent.f` V3 small-tree height increment. Zone-matched ST_COEF blocks (27),
# same matcher pattern as DG/HTG. VALIDATED: V3 HTG regression 91.7% bit-exact / 100% ≤1 ULP vs
# oracle (test/harness/britishcolumbia/regent_htg_validate.jl). ⚠ CON is SITE-ADJUSTED (like DGCON).
#
# Data: data/britishcolumbia/regent_stcoef.csv (27 ST_COEF: SPP, ZONE, CON/CASP/SASP/SLP/HT/LNHT/
# CCF/BAL/EXPMLT/HCCFL/CCFSL/SD). ⚠ V2 regime + SBS logistic branch + DBH-from-ht/blend/SIZCAP TODO.
# =============================================================================

struct BCStCoef
    spp::Vector{Int}
    zones::Vector{String}
    CON::Float32; CASP::Float32; SASP::Float32; SLP::Float32
    HT::Float32; LNHT::Float32; CCF::Float32; BAL::Float32
    EXPMLT::Float32; HCCFL::Float32; CCFSL::Float32; SD::Float32
end

function _bc_load_stcoef()
    rows = _bc_readcsv("regent_stcoef.csv"); h = rows[1]
    ci(n) = findfirst(==(n), h); f(r,n) = parse(Float32, strip(r[ci(n)]))
    out = BCStCoef[]
    for r in rows[2:end]
        isempty(strip(r[1])) && continue
        push!(out, BCStCoef(_bc_splitlist(r[ci("spp")], Int), _bc_splitstr(r[ci("zones")]),
            f(r,"CON"), f(r,"CASP"), f(r,"SASP"), f(r,"SLP"), f(r,"HT"), f(r,"LNHT"),
            f(r,"CCF"), f(r,"BAL"), f(r,"EXPMLT"), f(r,"HCCFL"), f(r,"CCFSL"), f(r,"SD")))
    end
    out
end
const BC_STCOEF = _bc_load_stcoef()

# --- small-tree DBH-from-height (regent.f:186-195): D = HHT1·(H-4.5)^HHT2 + DADJ (H>4.5) else DIAM+DADJ ---
const BC_RG_HHT1 = Float32[0.0781,0.0751,0.0828,0.1155,0.0729,0.0730,0.0988,0.0658,0.0658,0.2160,
                           0.0729,0.0729,0.0729,0.0828,0.0729]
const BC_RG_HHT2 = Float32[1.1645,1.1176,1.1713,1.0688,1.1988,1.2343,1.0807,1.3817,1.3817,1.0049,
                           1.1988,1.1988,1.1988,1.1713,1.1988]
const BC_RG_DIAM = Float32[0.4,0.3,0.3,0.3,0.2,0.2,0.4,0.3,0.3,0.5,0.2,0.2,0.2,0.3,0.2]
const BC_RG_HSIGMA = 0.59f0                            # V2 ZZRAN mult (V3 uses ST_COEF%SD)
# small/large blend bounds XMINV3/XMAXV3 (per sp; sp14 measured=2/4; ⚠ others TODO via one instrument run)
const BC_RG_XMIN = fill(2.0f0, 15)                     # TODO: instrument all sp (sp14=2.0 confirmed)
const BC_RG_XMAX = fill(4.0f0, 15)                     # TODO: instrument all sp (sp14=4.0 confirmed)

"""small-tree DBH from height (regent.f:1361-1363): power form + DELMAX/DADJ density adj."""
@inline function bc_st_dbh(sp::Integer, h::Real, dadj::Real)
    h = Float32(h)
    h > 4.5f0 ? BC_RG_HHT1[sp]*(h - 4.5f0)^BC_RG_HHT2[sp] + Float32(dadj) : BC_RG_DIAM[sp] + Float32(dadj)
end

"""Resolve the ST_COEF block for `sp` (regent.f 2-pass PrettyName match, fallback `zone*"/all "`)."""
function bc_resolve_stcoef(sp::Integer, series::AbstractString, zone::AbstractString)
    for pat in (series, zone * "/all ")
        for (idx, b) in enumerate(BC_STCOEF)
            (sp in b.spp) || continue
            for z in b.zones
                occursin(z, pat) && return idx
            end
        end
    end
    0
end

"""Site-adjusted small-tree constant (regent HTCONS): ST_COEF%CON + aspect/slope terms (like DGCON)."""
function bc_stcon(ip::Integer, aspect::Real, slope::Real)
    ip < 1 && return 0f0
    b = BC_STCOEF[ip]; asp = Float32(aspect); slp = Float32(slope)
    return b.CON + b.CASP*cos(asp)*slp + b.SASP*sin(asp)*slp + b.SLP*slp
end

"""
    bc_v3_sthg(sp, ip, h1_ft, bal, rdj, rhcon, aspect, slope) -> small-tree HTG (ft)

Port of the V3 non-SBS small-tree height regression (canada/bc/regent.f:1396-1402). VALIDATED
91.7% bit-exact. `h1_ft`=working height, `bal`=tree BAL (ft²/ac), `rdj`=RDNEXT (CCF term),
`rhcon`=RHCON small-tree calib. ⚠ SBS-zone logistic branch NOT ported (all_BC is ICH).
"""
function bc_v3_sthg(sp::Integer, ip::Integer, h1_ft::Real, bal::Real, rdj::Real, rhcon::Real,
                    aspect::Real, slope::Real)
    ip < 1 && return 0f0                                    # STG%FIT false
    b = BC_STCOEF[ip]
    htm = max(2f0, Float32(h1_ft) * BC_FTtoM)
    con = bc_stcon(ip, aspect, slope)
    pre = con + Float32(rhcon) + b.LNHT*log(htm) + b.HT*htm +
          b.BAL*(Float32(bal) * BC_FT2pACRtoM2pHA * 0.01f0) + b.CCF*Float32(rdj)
    return exp(max(-88f0, pre)) * 3.28084f0                 # MtoFT
end
