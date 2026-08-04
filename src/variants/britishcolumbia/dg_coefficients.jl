# =============================================================================
# dg_coefficients.jl (britishcolumbia) — BC large-tree diameter growth coefficients + model.
#
# ⚠ SUPERSEDES the earlier IE-premise scaffold. MEASURED from canada/bc/dgf.f (doctrine #2):
# BC DG is NOT the IE Wykoff DDS. It is a dual-regime, BEC-biogeoclimatic-zone-driven model
# (LV2ATV toggle, becset.f:855): zones ICH/IDF/SBS/SBPS ⇒ "V3" metric species/site-specific model
# (ported here); other zones ⇒ "V2" IE-like form (TODO). Coefficients are per-(species,site-series)
# in ZNKONST (zone) + SSKONST (site) DATA blocks, resolved by matching the stand's BEC PrettyName.
#
# Data: data/britishcolumbia/dg_znkonst.csv (28), dg_sskonst.csv (25) — extracted+verified from
# dgf.f (provenance/extract_dg_coeffs.py). Validated: PrettyName→IP/JP matcher 15/15 vs oracle
# (test/harness/britishcolumbia/dg_match_validate.jl); V3 default DDS 97.1% bit-exact / 100% ≤1 ULP
# (dg_dds_validate.jl) on all_BC (ICH/ICHmw2/01, sp14→ZNK(5)).
# =============================================================================

# --- metric constants (common/METRIC.F77 — EXACT truncated literals, critical for bit-exactness) ---
const BC_CMtoIN         = 0.3937f0        # NOT 1/2.54 — using 1/2.54 breaks bit-exactness
const BC_INtoCM         = 2.54f0
# Per-species population residual variance for the DGSCOR empirical-Bayes COR weight (bc/dgdriv.f:104
# DATA PSIGSQ; EP/AT/AC taken from WC). Indexed by species 1..15. Used by calibrate_diameter_growth!.
const BC_PSIGSQ = Float32[0.0408, 0.0586, 0.1556, 0.0970, 0.0858, 0.1433, 0.0636, 0.0970, 0.0970,
                          0.0636, 0.0898, 0.0898, 0.0898, 0.1556, 0.0898]
const BC_FTtoM          = 0.3048f0
const BC_FT2pACRtoM2pHA = 0.2295643f0

# --- V3 zone-constant record (MD_STR, canada/bc/dgf.f:85) ---
struct BCZoneCoef
    spp::Vector{Int}          # species this block serves (primary + surrogate), e.g. [3,14]
    zones::Vector{String}     # ZONE(P_ST) BEC/site-series PrettyName list
    OBSERV::Float32
    CON::Float32; CASP::Float32; SASP::Float32; EL::Float32; EL2::Float32; CCFA::Float32
    LD::Float32; DSQ::Float32; DBAL1::Float32; DBAL2::Float32; CR::Float32; BAL::Float32; SIGMAR::Float32
end

# --- V3 site record (SS_STR, canada/bc/dgf.f:108) ---
struct BCSiteCoef
    spp::Vector{Int}
    prettynames::Vector{String}
    CON::Float32
end

_bc_readcsv(f) = (rows = split(strip(replace(read(joinpath(BC_DATADIR, f), String), "\r" => "")), '\n');
                  [split(l, ',') for l in rows])
_bc_splitlist(s, T) = T[parse(T, x) for x in split(s, ';') if !isempty(strip(x))]
_bc_splitstr(s)     = String[strip(x) for x in split(s, ';') if !isempty(strip(x))]

function _bc_load_znkonst()
    rows = _bc_readcsv("dg_znkonst.csv"); h = rows[1]
    ci(n) = findfirst(==(n), h)
    f32(r, n) = parse(Float32, strip(r[ci(n)]))
    out = BCZoneCoef[]
    for r in rows[2:end]
        isempty(strip(r[1])) && continue
        push!(out, BCZoneCoef(
            _bc_splitlist(r[ci("spp")], Int), _bc_splitstr(r[ci("zones")]),
            f32(r,"OBSERV"), f32(r,"CON"), f32(r,"CASP"), f32(r,"SASP"), f32(r,"EL"), f32(r,"EL2"),
            f32(r,"CCFA"), f32(r,"LD"), f32(r,"DSQ"), f32(r,"DBAL1"), f32(r,"DBAL2"),
            f32(r,"CR"), f32(r,"BAL"), f32(r,"SIGMAR")))
    end
    out
end

function _bc_load_sskonst()
    rows = _bc_readcsv("dg_sskonst.csv"); h = rows[1]
    ci(n) = findfirst(==(n), h)
    out = BCSiteCoef[]
    for r in rows[2:end]
        isempty(strip(r[1])) && continue
        push!(out, BCSiteCoef(_bc_splitlist(r[ci("spp")], Int),
                              _bc_splitstr(r[ci("prettynames")]), parse(Float32, strip(r[ci("CON")]))))
    end
    out
end

const BC_ZNKONST = _bc_load_znkonst()
const BC_SSKONST = _bc_load_sskonst()

# --- bark ratio (canada/bc/bratio.f): IMAP=2 & BARK2=0 for all sp ⇒ BRATIO = BARK1(sp), clamp [0.80,0.99] ---
const BC_BARK1 = Float32[0.964, 0.851, 0.867, 0.915, 0.934, 0.950, 0.969,
                         0.956, 0.937, 0.890, 0.939, 0.916, 0.897, 0.867, 0.939]
@inline bc_bratio(sp::Integer) = clamp(BC_BARK1[sp], 0.80f0, 0.99f0)

"""
    bc_resolve_ipjp(sp, series, zone) -> (ip, jp)

Port of canada/bc/dgf.f:2140-2181 PrettyName resolution. Two passes (site-series `PAT2`, then
zone fallback `zone*"/all "`); first block (in DATA order) whose `spp` contains `sp` AND some
PrettyName is a substring of the pattern. Returns (0,0) if unresolved. VALIDATED 15/15 on all_BC.
"""
function bc_resolve_ipjp(sp::Integer, series::AbstractString, zone::AbstractString)
    _find(blocks, names) = begin
        for pat in (series, zone * "/all ")
            for (idx, b) in enumerate(blocks)
                (sp in b.spp) || continue
                for z in names(b)
                    occursin(z, pat) && return idx
                end
            end
        end
        0
    end
    ip = _find(BC_ZNKONST, b -> b.zones)
    jp = _find(BC_SSKONST, b -> b.prettynames)
    return (ip, jp)
end

"""
    bc_dgcon(ip, jp, elev_ft, aspect, slope) -> DGCON

Port of DGCONS (canada/bc/dgf.f:2200): site-dependent constant, resolved once/stand/species.
PELEV = elev_ft*FTtoM (m); white-pine (sp1) caller clamps PELEV∈[5,12] — pass clamped elev for sp1.
"""
function bc_dgcon(ip::Integer, jp::Integer, elev_ft::Real, aspect::Real, slope::Real)
    (ip < 1 || jp < 1) && return 0f0
    z = BC_ZNKONST[ip]; ss = BC_SSKONST[jp]
    pelev = Float32(elev_ft) * BC_FTtoM
    return ss.CON + z.CON + z.EL*pelev + z.EL2*(pelev^2) +
           z.SASP*sin(Float32(aspect))*Float32(slope) + z.CASP*cos(Float32(aspect))*Float32(slope)
end

"""
    bc_v3_dds(sp, ip, zone, d_in, bal_m2ha, cr_frac, conspp, brat) -> WK2 (ln DDS, inches)

Port of the V3 metric DDS (canada/bc/dgf.f:1943-1996). `d_in`=DBH inches, `bal_m2ha`=BAL (m²/ha),
`cr_frac`=crown ratio, `conspp`=DGCON+COR+DGCCF1·RELDN2·0.01, `brat`=BRATIO(sp, d_cm, 0).
Default branch VALIDATED 97.1% bit-exact / 100% ≤1 ULP. Birch(11,15)/aspen(12,13) forms per source.
"""
function bc_v3_dds(sp::Integer, ip::Integer, zone::AbstractString,
                   d_in::Real, bal_m2ha::Real, cr_frac::Real, conspp::Real, brat::Real)
    z = BC_ZNKONST[ip]
    D   = Float32(d_in) * BC_INtoCM       # cm
    BAL = Float32(bal_m2ha)
    CR  = Float32(cr_frac)
    D2  = max(D, 1f0)
    bald1 = BAL / D2
    _power() = z.LD * (D2 ^ (z.DBAL1 + z.DBAL2 * bald1)) * exp(z.DSQ * D2 * D2)
    _linexp() = exp(max(-9.21f0,
        Float32(conspp) + z.LD*D + z.DSQ*D*D + z.BAL*BAL +
        z.DBAL1*(BAL/D) + z.DBAL2*(BAL/log(D+1f0)) + z.CR*CR))
    dds = if (sp == 11 || sp == 15)
        occursin("IDF", zone) ? _power() : _linexp()
    elseif (sp == 12 || sp == 13)
        _power()
    else
        _linexp()
    end
    # CM DG → ln(DDS) inches via DIB (dgf.f:1987)
    dds = (dds*dds + 2f0*dds*D*Float32(brat)) * BC_CMtoIN * BC_CMtoIN
    return max(-9.21f0, log(max(0.001f0, dds)))
end
