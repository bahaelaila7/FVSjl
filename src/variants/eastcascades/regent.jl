# =============================================================================
# regent.jl (eastcascades) — EC small-tree growth (ec/regent.f + ec/smhtgf.f) + htdbh. Chunk 6.
#
# EC's small-tree model is BESPOKE (not the WC smhgdg reuse): ec/smhtgf.f is HEIGHT-growth-only with
# inline per-species Chapman-Richards/curve literals (WP/RC Chapman-Richards; MH/OS ×3.281 metric), and
# small-tree DIAMETER growth lives in ec/regent.f (per-species H→DBH inversion bounded by DGMAX·SCALE,
# SLO/SHI SI clamp, no DGMIN). Missing-height dubbing + small-tree DBH inversion use the forest-dependent
# Curtis-Arney htdbh (ec/htdbh.f: MTHOOD/OKANOG/GIFFPC/WENATC by IFOR/IGL).
# =============================================================================

# ec/htdbh.f — 4 physical forest tables (MTHOOD/OKANOG/GIFFPC/WENATC), stored in htdbh_coeffs_ec.csv by
# resolved CSV forest_idx {1=MTHOOD, 2/4/7=OKANOG, 3=WENATC(IGL1), 5/6=GIFFPC(IGL2,3)} × 32 species.
let
    P2 = zeros(Float32, 7, 32); P3 = zeros(Float32, 7, 32); P4 = zeros(Float32, 7, 32)
    for l in readlines(joinpath(EC_DATADIR, "htdbh_coeffs_ec.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        fi = parse(Int, f[1]); sp = parse(Int, f[2])
        P2[fi, sp] = parse(Float32, f[3]); P3[fi, sp] = parse(Float32, f[4]); P4[fi, sp] = parse(Float32, f[5])
    end
    global const EC_HTDBH_P2 = P2; global const EC_HTDBH_P3 = P3; global const EC_HTDBH_P4 = P4
end

# ec/htdbh.f SELECT CASE(IFOR)/CASE(IGL) → CSV forest_idx. forkod passes the CORRECTED IFOR (1..4) + IGL.
@inline function _ec_htdbh_fidx(ifor::Int, igl::Int)::Int
    ifor == 1 && return 1                 # MTHOOD
    (ifor == 2 || ifor == 4) && return 2  # OKANOG
    (ifor == 3 && (igl == 2 || igl == 3)) && return 5  # GIFFPC
    return 3                              # WENATC (IFOR 3, IGL 1) / default
end

@inline function ec_htdbh_height(fidx::Int, sp::Int, d::Float32)::Float32
    p2 = EC_HTDBH_P2[fidx, sp]; p3 = EC_HTDBH_P3[fidx, sp]; p4 = EC_HTDBH_P4[fidx, sp]
    if d >= 3.0f0
        return 4.5f0 + p2 * exp(-1f0 * p3 * d^p4)
    else
        return ((4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
    end
end
@inline function ec_htdbh_dbh(fidx::Int, sp::Int, h::Float32)::Float32
    p2 = EC_HTDBH_P2[fidx, sp]; p3 = EC_HTDBH_P3[fidx, sp]; p4 = EC_HTDBH_P4[fidx, sp]
    hat3 = 4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4)
    if h >= hat3
        return exp(log((log(h - 4.5f0) - log(p2)) / (-1f0 * p3)) * (1f0 / p4))
    else
        return (((h - 4.51f0) * 2.7f0) / (hat3 - 4.51f0)) + 0.3f0
    end
end

@inline ec_htdbh_ifor(p)::Int = _ec_htdbh_fidx(Int(p.forest_idx), Int(p.geo_location))

regenerate!(s::StandState, ::EastCascades; kwargs...) = s
