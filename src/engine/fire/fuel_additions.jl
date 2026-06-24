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
(`xv[2..6]`), now correctly scaled after the V2T /2000 fix. The crown-LIFT term (fmcadd.f:95-102,
dead material shed as the crown base rises) needs previous-cycle crown tracking and is deferred —
small for a closing-canopy stand. No-op unless FFE is active.
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
