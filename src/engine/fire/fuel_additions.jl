# fire/fuel_additions.jl — FFE annual fuel additions (FMCADD, fire/base/fmcadd.f).
#
# Each YEAR (the FFE fuel update runs annually, fmmain.f:226-259 with NYRS=1) live trees add debris to
# the surface-fuel pools `fire.cwd`:
#   - litterfall: foliage biomass ÷ leaf lifespan → the size-10 litter pool (this routine);
#   - woody crown breakage (LIMBRK) + crown-lift dead material → woody pools (chunk 2b);
#   - snag-crown falldown from the CWD2B debris-in-waiting pool (chunk 2b, couples to snag.jl).
# Material enters the tree species' decay class (`dkr_cls`). Foliage comes from `crown_biomass`
# (FMCROWE) and is in FFE-internal pounds; P2T = 1/2000 converts to tons/acre.

const _FM_LIMBRK = 0.01f0   # woody-crown limb-breakage fraction per year (sn/fmvinit.f:126)
const _FM_P2T = 1f0 / 2000f0

"""
    crown_lift_rate(oldht, oldcrl, ht, crown_pct, cyclen) -> Float32

The annual crown-base-rise fraction **X** of the FFE crown-LIFT term (FMSDIT, fmsdit.f:103-117). As a
tree's crown base rises over `cyclen` years from the previous-cycle base `OLDBOT = oldht − oldcrl` to the
current base `NEWBOT = ht − ht·crown_pct/100`, the fraction `X = (NEWBOT − OLDBOT)/oldcrl/cyclen` of the
OLD crown is shed (dies) into down wood EACH year — the dominant post-mortality down-wood addition
(`X·CROWNW·TPA·P2T`, validated at ~0.39 t/ac/yr vs an instrumented FMCADD dump; see
FFE_FUEL_DYNAMICS_chunk_plan.md). Returns 0 if the old crown length is ~0 or the base did not rise.

This is the standalone, faithful semantic of the crown-lift rate; wiring it into `fmcadd_woody!` needs
the PREVIOUS-cycle per-tree `oldht`/`oldcrl`, which require tree-record tracking across the
regen/mortality-changing tree list (the remaining plumbing — kept separate so the formula is locked and
tested independently).
"""
@inline function crown_lift_rate(oldht::Real, oldcrl::Real, ht::Real, crown_pct::Real, cyclen::Real)::Float32
    (oldcrl > 0.001f0 && cyclen > 0) || return 0f0
    newbot = Float32(ht) - Float32(ht) * Float32(crown_pct) / 100f0
    oldbot = Float32(oldht) - Float32(oldcrl)
    rise = newbot - oldbot
    return rise > 0f0 ? rise / Float32(oldcrl) / Float32(cyclen) : 0f0
end

"""
    ffe_fuel_update!(s, nyrs) -> StandState

The per-cycle FFE surface-fuel update (the deterministic core of fmmain.f:226-259): refresh the
cover type + live fuels and load the initial dead fuels once (`fmcba!`), then run the ANNUAL fuel
loop `nyrs` times — decay (`fmcwd!`) + litterfall and woody breakage (`fmcadd_*`). The tree crowns
are held at the cycle's start (their end-of-previous-cycle state, as FVS records them). No-op unless
FFE is active. This is what makes the Stand Carbon Report's dead/down-wood/floor pools evolve across
grown cycles (validated bit-exact on carbon_jenkins). Snag-debris falldown (CWD2B) is still pending.
"""
const _FM_CRDCAY = 0.0425f0   # dead coarse-root decay rate per year (fminit.f:918)

# TFALL — years for a crown component to fall, by `tfall_cls` (1..6) and crown size 0..5
# (sn/fmvinit.f:1018-1058). Foliage(0)=1; branch(1,2)=row1; size 3=row3; size 4,5=row4.
const _FM_TFALL1 = (5f0, 3f0, 2f0, 1f0, 1f0, 1f0)
const _FM_TFALL3 = (10f0, 6f0, 5f0, 4f0, 3f0, 2f0)
const _FM_TFALL4 = (25f0, 12f0, 10f0, 8f0, 6f0, 4f0)
@inline function _fm_tfall(cls::Int, sz::Int)::Float32
    sz == 0 && return 1f0
    (sz == 1 || sz == 2) && return _FM_TFALL1[cls]
    sz == 3 && return _FM_TFALL3[cls]
    return _FM_TFALL4[cls]
end

"""
    fmscro!(s, sp, dbh, xv, density, dkcl)

Schedule a dying tree's crown debris into the CWD2B debris-in-waiting pool (FMSCRO, fmscro.f): each
crown component `xv[size]·density` is spread EQUALLY over years 1..min(TSOFT, TFALL(sp,size)), where
TSOFT = `(1.24·dbh + 13.82)·DECAYX` (FMSNGDK). The un-fallen CWD2B is the Stand-Dead crown; it flows
to the down-wood pool as it falls. `xv` is the `crown_biomass` tuple (foliage, woody1-5), in lb.
"""
function fmscro!(s::StandState, sp::Integer, dbh::Float32, xv, density::Float32, dkcl::Integer)
    fs = s.fire; coef = s.coef
    cls = clamp(Int(coef_col(coef, :tfall_cls)[sp]), 1, 6)
    tsoft = (1.24f0 * dbh + 13.82f0) * coef_col(coef, :snag_decayx)[sp]
    @inbounds for sz in 0:5
        amt = xv[sz + 1] * density
        amt > 0f0 || continue
        ilife = clamp(round(Int, min(tsoft, _fm_tfall(cls, sz))), 1, 60)
        annual = amt / ilife
        for yr in 1:ilife
            fs.cwd2b[dkcl, sz + 1, yr] += annual
        end
    end
    return s
end

"""
    snag_crown_carbon(s) -> Float32

The Stand-Dead CROWN carbon (tons C/acre): the crown debris still in the CWD2B waiting pool
(not yet fallen to down wood), summed × P2T × 0.5 (fmdout.f:173). The other half of Stand-Dead
is `snag_bole_carbon`.
"""
snag_crown_carbon(s::StandState)::Float32 =
    (s.fire === nothing ? 0f0 : sum(s.fire.cwd2b) * _FM_P2T) * 0.5f0

# One year of CWD2B falldown → the down-wood pools (FMCADD, fmcadd.f:122-135): the year-1 pool of each
# crown size flows to cwd (foliage size-0 → litter cwd[10]; woody 1-5 → cwd[1-5]) at P2T, then shift.
function _cwd2b_fall!(fs::FireState)
    c2 = fs.cwd2b
    @inbounds for dkcl in 1:4, sz in 0:5
        down = c2[dkcl, sz + 1, 1]
        down > 0f0 || continue
        fs.cwd[sz == 0 ? 10 : sz, 2, dkcl] += down * _FM_P2T
    end
    @inbounds for dkcl in 1:4, sz in 1:6
        for yr in 1:59
            c2[dkcl, sz, yr] = c2[dkcl, sz, yr + 1]
        end
        c2[dkcl, sz, 60] = 0f0
    end
    return fs
end

function ffe_fuel_update!(s::StandState, nyrs::Integer)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    fmcba!(s)
    @inbounds for _ in 1:nyrs
        _cwd2b_fall!(fs)                       # FMCADD: CWD2B crown debris → down wood
        fmcwd!(s, 1); fmcadd_litterfall!(s); fmcadd_woody!(s)
    end
    fs.bioroot *= (1f0 - _FM_CRDCAY)^nyrs    # dead-root decay (fmcrbout.f:273)
    return s
end

"""
    fmcadd_litterfall!(s) -> StandState

Add ONE year of foliage litterfall to the FFE litter pool (FMCADD, fmcadd.f:72-76): for each live
tree, `foliage · TPA / LEAFLF(sp) · P2T` tons/ac into `cwd[10, hard, dkr_cls(sp)]`. This is the term
that offsets the 0.65/yr litter decay so the forest-floor pool holds up. No-op unless FFE is active.
(Woody breakage + snag falldown are the companion FMCADD terms — chunk 2b.)
"""
function fmcadd_litterfall!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    t = s.trees; coef = s.coef
    leaflf = coef_col(coef, :leaf_life); dkrcls = coef_col(coef, :dkr_cls)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i])
        ll = leaflf[sp]; ll <= 0f0 && continue
        xv = crown_biomass(s, sp, t.dbh[i], t.height[i], Int(round(t.crown_pct[i])))
        dkcl = clamp(Int(dkrcls[sp]), 1, 4)
        fs.cwd[10, 2, dkcl] += xv[1] * t.tpa[i] / ll * _FM_P2T
    end
    return s
end

"""
    fmcadd_woody!(s) -> StandState

Add ONE year of woody crown-breakage debris to the FFE down-wood pools (FMCADD, fmcadd.f:78-84): for
each live tree and crown size class SIZE 1..5, `LIMBRK · CROWNW(SIZE) · TPA · P2T` tons/ac into
`cwd[SIZE, hard, dkr_cls(sp)]`. `CROWNW(SIZE)` is the woody crown component from `crown_biomass`
(`xv[2..6]`), now correctly scaled after the V2T /2000 fix. No-op unless FFE is active.

REMAINING — the crown-LIFT term (fmcadd.f:95-102): `X · CROWNW(SIZE) · TPA · P2T`, the lower crown that
dies as the crown base rises, where `X = (NEWBOT−OLDBOT)/OLDCRL/CYCLEN` is the annual crown-base-rise
fraction (FMSDIT, fmsdit.f:103-117). An instrumented FMCADD dump showed this is NOT small — it is the
DOMINANT post-mortality down-wood addition (~0.39 vs 0.15 t/ac/yr for breakage on carbon_jenkins), the
source of the last DDW residual (Fortran 3.8 vs FVSjl 2.1 @2000). It needs the PREVIOUS-cycle per-tree
height + crown length, tracked across the regen/mortality-changing tree list (FVS does this with the
OLDCRW record-maintenance in FMTDEL/FMTRIP/FMCMPR) — a stable tree-record id, not a naive index
snapshot (a first index-snapshot attempt failed because regen grows the list 6→18 each cycle). That is
the focused next step for the DDW column; see FFE_FUEL_DYNAMICS_chunk_plan.md.
"""
function fmcadd_woody!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    t = s.trees; coef = s.coef; dkrcls = coef_col(coef, :dkr_cls)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i])
        xv = crown_biomass(s, sp, t.dbh[i], t.height[i], Int(round(t.crown_pct[i])))
        dkcl = clamp(Int(dkrcls[sp]), 1, 4)
        for sz in 1:5
            fs.cwd[sz, 2, dkcl] += _FM_LIMBRK * xv[sz + 1] * t.tpa[i] * _FM_P2T
        end
    end
    return s
end


