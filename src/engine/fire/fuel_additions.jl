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

# ⛔ Woody crown-breakage additions (FMCADD, fmcadd.f:78-84: `LIMBRK·CROWNW(SIZE)·TPA·P2T` into the
# woody cwd pools) are NOT ported yet — they are BLOCKED on the `crown_biomass` WOODY components
# (`xv[2..6]`), which are in FMCROWE's "FFE-internal units, not literal tons" (a known fmcrowe.f
# scaling inconsistency: bole-tip frusta ×SG/P2T, sub-breast cylinder raw). For an 8" loblolly they
# come out as (8.7, 545, 14394, 20484) — clearly not tons — and applying P2T gives DDW ≈ 120 t/ac vs
# the report's 2.5. The FOLIAGE component (`xv[1]`, Jenkins-derived pounds) IS right (it lands the
# forest floor bit-exact above), so only the woody side is blocked. ⇒ chunk 2b needs FMCROWE woody-
# component validation/fix first (the F5/F6 crown-biomass chunk). See FFE_FUEL_DYNAMICS_chunk_plan.md.
