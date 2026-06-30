# =============================================================================
# diameter_growth.jl (northeast) — NE large-tree diameter increment (ne/dgf.f)
#
# Structurally different from SN: a Monserud-style potential basal-area growth
# modulated by BAL (basal-area-in-larger-trees) competition — the one genuinely
# new NE growth mechanism.
#
#   POTBAG = B1·SITEAR·(1 − exp(−B2·D))·0.7        (ne/dgf.f:131-132)
#   GMOD   = exp(−B3·BAL),  clamped ≥ 0.5           (ne/balmod.f)
#   one cycle = 10 ANNUAL steps; each step recomputes POTBAG & GMOD at the
#   current (growing) D, grows basal area: QTRBA = DELD + D²·0.0054542,
#   D ← √(QTRBA/0.0054542)                          (ne/dgf.f:127-151)
#   DIAGR = (Dfinal − D0)·bark ; DDS = DIAGR·(2·D0·bark + DIAGR)
#   WK2   = ln(DDS) + COR[sp]                        (COR = DG calibration)
#
# SITEAR(ISPC) is the per-species site index — already produced by the SICOEF
# fan-out (`ne_site_index_setup!` → plot.sp_site_index). BAL comes from BADIST.
# B1/B2/B3 are loaded from data/northeast/dg_coeffs.csv.
# =============================================================================

"""
    ne_badist!(ebau, s)

BADIST (ne/badist.f): basal-area-in-larger-trees by 1-inch DBH class. Fills the
50-element `ebau`: bin each tree's per-acre basal area into class `min(⌊D+1⌋,50)`,
then cumulate from the top so `ebau[c]` = Σ BA in classes ≥ c.
"""
function ne_badist!(ebau::Vector{Float32}, s::StandState)
    fill!(ebau, 0f0)
    tr = s.trees
    # During DG calibration the stand dbh is backdated for the per-tree PREDICTION, but FVS NE computes the
    # BAL competition on the CURRENT stand (verified live: badist.f EBAU=52 on every ICYC=1 call). So when the
    # calibration stash is set, build the BAL array from the current dbh; otherwise use the live tree dbh.
    dbh = isempty(s.calib.calib_dbh) ? tr.dbh : s.calib.calib_dbh
    @inbounds for i in 1:tr.n
        d = dbh[i]; d <= 0f0 && continue
        icls = min(floor(Int, d + 1f0), 50)
        # BADIST (ne/badist.f:45-47) floors DBH at 1.0 for the BA contribution (TDBH), but bins
        # by the ACTUAL DBH; PROB is already per-acre, so NO gross_space division (matching the
        # bit-exact `stand_ba`). The earlier `/gross_space` scaled every class by 10/11 — masked
        # for high-BAL trees by the GMOD≥0.5 floor, but it under-grew the large/low-competition
        # trees whose GMOD>0.5 reads the BAL directly.
        tdbh = d < 1f0 ? 1f0 : d
        ebau[icls] += 0.0054542f0 * tdbh * tdbh * tr.tpa[i]   # per-acre basal area (PROB·BA/tree)
    end
    @inbounds for i in 49:-1:1
        ebau[i] += ebau[i + 1]
    end
    return ebau
end

"BALMOD (ne/balmod.f): BAL competition modifier for species `sp`, DBH `d`."
@inline function ne_balmod(b3::Float32, ebau::Vector{Float32}, d::Float32)::Float32
    icls = floor(Int, d + 1f0) - 2          # competition from same-or-larger neighbours
    icls < 1 && (icls = 1); icls > 50 && (icls = 50)
    bal = @inbounds ebau[icls]
    bal <= 0f0 && return 1f0
    g = exp(-b3 * bal)
    return g < 0.5f0 ? 0.5f0 : g
end

"""
    ne_diameter_increment(s, i, ebau) -> Float32

Outside-bark DBH increment (in) for tree `i` over one cycle — the 10 annual
DGF steps. Pure (no calibration COR / bark here; the caller applies those).
"""
function ne_diameter_increment(s::StandState, i::Integer, ebau::Vector{Float32})::Float32
    tr = s.trees; sp = Int(tr.species[i])
    sd = s.coef.species
    b1 = sd[:dg_b1][sp]; b2 = sd[:dg_b2][sp]; b3 = sd[:dg_b3][sp]
    sitear = s.plot.sp_site_index[sp]
    d0 = tr.dbh[i]; d0 <= 0f0 && return 0f0
    d = d0
    @inbounds for _ in 1:10
        potbag = b1 * sitear * (1f0 - exp(-(b2 * d))) * 0.7f0
        gmod = ne_balmod(b3, ebau, d)
        deld = potbag * gmod
        qtrba = deld + d * d * 0.0054542f0
        d = sqrt(qtrba / 0.0054542f0)
    end
    return d - d0
end

"""
    ne_dgf!(s)

NE variant fill of `scratch.wk[2, i] = ln(DDS) + COR` — the `ne/dgf.f` analog of SN's `dgf!`,
i.e. the only variant-specific part of the (otherwise shared) DGDRIV pipeline. DDS is the
inside-bark Δ(d²) implied by the BAL model: `diagr = ne_diameter_increment · bark` (inside-bark
increment), `dib = d·bark`, `DDS = diagr·(2·dib + diagr)`. COR = `calib.dg_cor[sp]` (the DG
calibration adjustment). Clamped `≥ −9.21` exactly like `dgf!`. BADIST (the BAL array) is the
cycle-start basis — computed once. Bark = the constant `BKRAT` via the shared `bark_ratio`
(NE coefs: intercept 0, slope BKRAT).
"""
function ne_dgf!(s::StandState)
    t = s.trees; c = s.calib; sd = s.coef.species; ctl = s.control
    wk2 = view(s.scratch.wk, 2, :)
    ba = sd[:bark_intercept]; bb = sd[:bark_slope]
    ebau = zeros(Float32, 50); ne_badist!(ebau, s)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        bark  = bark_ratio(ba, bb, sp, d)
        dib   = d * bark
        diagr = ne_diameter_increment(s, i, ebau) * bark
        # COR2 bark-growth compensation (dgf.f:158, LDCOR2-gated) THEN the DIAGR≥.0001 floor (dgf.f:159) —
        # the floor is NOT calibration-gated and guarantees DDS>0, so NE dgf has NO −9.21 clamp (dgf.f:169),
        # unlike SN's dgf!. (Was wrongly carrying the SN −9.21 isms; the floor bottoms WK2 near −8.6.)
        ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0 && (diagr *= ctl.dg_cor2[sp])
        diagr < 0.0001f0 && (diagr = 0.0001f0)
        dds   = diagr * (2f0 * dib + diagr)
        wk2[i] = log(dds) + c.dg_cor[sp]
    end
    return s
end

# Variant `dgf!` hook (the only variant-specific step of the shared DGDRIV pipeline): NE fills
# wk2 from the BAL model. SN's `dgf!(s, ::Southern)` is the linear ln(DDS) model. Both are called
# as `dgf!(s, s.variant)` by the (to-be-generic) calibrate / diameter_growth! driver.
dgf!(s::StandState, ::Northeast) = ne_dgf!(s)

"""
    ne_dgcons!(s)

NE per-stand DG setup (the `dgcons!` analog). ne/dgf.f:188 zeros DGCON/ATTEN/SMCON — NE's DDS comes
straight from B1/B2/B3 + SITEAR + BAL, with no SN-style site/slope/forest-type linear constant. So
the only state this needs to populate is the per-stand bark copy (`calib.bark_a/b` ← the constant
`BKRAT`, encoded as intercept 0 / slope BKRAT) that the shared calibrate / DGDRIV downstream reads.
"""
function ne_dgcons!(s::StandState)
    c = s.calib; sd = s.coef.species
    ba = sd[:bark_intercept]; bb = sd[:bark_slope]
    @inbounds for sp in 1:MAXSP
        c.bark_a[sp] = ba[sp]; c.bark_b[sp] = bb[sp]
        c.dg_const[sp] = 0f0; c.atten[sp] = 1000f0   # DGCONS ATTEN=1000 (dgf.f:195) — prior-obs weight
    end
    return s
end
