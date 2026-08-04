# =============================================================================
# height_growth_coefficients.jl (britishcolumbia) — BC large-tree HEIGHT growth (canada/bc/htgf.f, V3).
#
# CHUNK 4. Mirrors chunk-3 DG structure: dual-regime (LV2ATV V2-old / V3-new), HTCONS site setup,
# BEC-zone PrettyName resolution (same 2-pass matcher). HTG predicted from species/HT/DBH/DG
# (⇒ needs chunk-3 DG output). V3 core = V3HTG (htgf.f:2094-2113); gated by FIT (zone matched).
#
# Coefficients: `dg... `no — `htg_lthg.csv` (57 LTHG blocks: SPP, PRETTYNAME[zones], SI/B3/C [old],
# CN/CNSI/DG/DG2/DBH1/DBH2/HT2 [new]) + LTSP[15,3] = A1/B1/B2 per species (htgf.f:112-126, hardcoded
# below). MLT = HCOR2 (height calib, engine state; default 1). ⚠ V2 regime (htgf.f:1620/1796) TODO.
# =============================================================================

# LTSP[sp] = (A1, B1, B2) old-v3 per-species constants (htgf.f:112-126). New-model sp (7,8,11,12,13,15) = 0.
const BC_LTSP = Float32[
    4.1060   0.0      0.0;      # 1 PW
    2.2341  -0.0532  -0.00222;  # 2 LW
    2.6937  -0.2034  -0.00387;  # 3 FD
    2.0204   0.0      0.0;      # 4 BG [use BL]
    1.4113   0.4019  -0.00760;  # 5 HW [use CW]
    1.4113   0.4019  -0.00760;  # 6 CW
    0.0      0.0      0.0;      # 7 PL  (new model)
    0.0      0.0      0.0;      # 8 SE  (new model)
    2.0204   0.0      0.0;      # 9 BL
    5.1763  -0.3013   0.0;      # 10 PY
    0.0      0.0      0.0;      # 11 EP (new model)
    0.0      0.0      0.0;      # 12 AT (new model)
    0.0      0.0      0.0;      # 13 CT (new model)
    2.6937  -0.2034  -0.00387;  # 14 OC [FD]
    0.0      0.0      0.0 ]     # 15 OH [EP] (new model)

# LTHG zone-coefficient record (htgf.f:90 TYPE LTHG_STR) — the site-series block.
struct BCLthgCoef
    spp::Vector{Int}
    zones::Vector{String}       # PRETTYNAME(P_SS)
    SI::Float32; B3::Float32; C::Float32                    # old-v3 site terms
    CN::Float32; CNSI::Float32; DG::Float32; DG2::Float32   # new-v3
    DBH1::Float32; DBH2::Float32; HT2::Float32
end

function _bc_load_lthg()
    rows = _bc_readcsv("htg_lthg.csv"); h = rows[1]
    ci(n) = findfirst(==(n), h); f(r,n) = parse(Float32, strip(r[ci(n)]))
    out = BCLthgCoef[]
    for r in rows[2:end]
        isempty(strip(r[1])) && continue
        push!(out, BCLthgCoef(_bc_splitlist(r[ci("spp")], Int), _bc_splitstr(r[ci("zones")]),
            f(r,"SI"), f(r,"B3"), f(r,"C"), f(r,"CN"), f(r,"CNSI"), f(r,"DG"), f(r,"DG2"),
            f(r,"DBH1"), f(r,"DBH2"), f(r,"HT2")))
    end
    out
end
const BC_LTHG = _bc_load_lthg()

"""Resolve the LTHG block index for species `sp` (htgf.f:1740-1762, 2-pass PrettyName match).
Fallback pattern is `zone*"/all"` (NO trailing space — cf. DG's `"/all "`). 0 ⇒ FIT=false."""
function bc_resolve_lthg(sp::Integer, series::AbstractString, zone::AbstractString)
    for pat in (series, zone * "/all")
        for (idx, b) in enumerate(BC_LTHG)
            (sp in b.spp) || continue
            for z in b.zones
                occursin(z, pat) && return idx
            end
        end
    end
    0
end

"""
    bc_v3_htg(sp, ip, ht_ft, dbh_in, dg_in; mlt=1) -> HTG (ft)

Port of V3HTG (canada/bc/htgf.f:2094-2113). `ip`=LTHG index (0 ⇒ FIT=false ⇒ 0). NEW branch for
sp∈{7,8,11,12,13,15}, OLD branch (A1/B1/B2 from LTSP + SI/B3/C from LTHG) otherwise. `mlt`=HCOR2.
"""
function bc_v3_htg(sp::Integer, ip::Integer, ht_ft::Real, dbh_in::Real, dg_in::Real; mlt::Real = 1f0)
    ip < 1 && return 0f0                       # LTH%FIT false
    L = BC_LTHG[ip]
    HTM = Float32(ht_ft)  * BC_FTtoM
    D   = Float32(dbh_in) * BC_INtoCM
    DGM = Float32(dg_in)  * BC_INtoCM
    htgm = if sp in (7, 8, 11, 12, 13, 15)     # NEW v3
        (L.CN + L.CNSI + L.DG*DGM + L.DG2*DGM*DGM) *
            (HTM ^ (L.DBH1 + L.DBH2*D)) * exp(L.HT2*HTM*HTM)
    else                                        # OLD v3
        A1, B1, B2 = BC_LTSP[sp,1], BC_LTSP[sp,2], BC_LTSP[sp,3]
        (A1 + L.SI) * (HTM ^ (B1 + B2*D + L.B3*DGM)) * exp(L.C*HTM*HTM)
    end
    return htgm * 3.28084f0 * Float32(mlt)      # MtoFT
end
