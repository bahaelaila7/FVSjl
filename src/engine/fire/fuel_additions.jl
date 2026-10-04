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
# (sn/fmvinit.f:1017-1058). Foliage(0)=1 EXCEPT redcedar=3; branch(1,2)=row1; size 3=row3; size 4,5=row4.
const _FM_TFALL1 = (5f0, 3f0, 2f0, 1f0, 1f0, 1f0)
const _FM_TFALL3 = (10f0, 6f0, 5f0, 4f0, 3f0, 2f0)
const _FM_TFALL4 = (25f0, 12f0, 10f0, 8f0, 6f0, 4f0)
# `sp` is the SN species index; fmvinit.f:1017 `IF(I.EQ.2)` gives eastern redcedar a 3-yr foliage fall vs 1 for
# all others. (SN-scoped, like the _FM_TFALL tables themselves — a variant porting NE would re-source these.)
@inline function _fm_tfall(cls::Int, sz::Int, sp::Integer)::Float32
    sz == 0 && return sp == 2 ? 3f0 : 1f0          # foliage (redcedar = 3, fmvinit.f:1018)
    (sz == 1 || sz == 2) && return _FM_TFALL1[cls]
    sz == 3 && return _FM_TFALL3[cls]
    return _FM_TFALL4[cls]
end

# {ie,em,kt}/fmvinit.f:500-513 (byte-identical block): TFALL(I,1)=5, TFALL(I,2)=MIN(5,TFALL(I,3)), TFALL(I,3) per species
# (the CSV tfall_cls column, verified against each fmvinit CASE), TFALL(I,4)=TFALL(I,5)=TFALL(I,3), and the dead-leaf
# fall TFALL(I,0)=MIN(2,LEAFLF(I)). The SN class rows above clamp these variants to row 6 (foliage 1 yr, branches 1 yr,
# size 3-5 2/4 yr) — MEASURED FVSie_g16 11855985010690: LP foliage CWD2B(4,0,·) 118.08 in slots 1 AND 2 (TFALL(7,0)=
# MIN(2,3)=2), jl 236.16 all in slot 1.
# so/fmcba.f:1033-1099 (SO-FFE initialization at FMCBA): TFALL(I,0:2) = 3/10/15 on the California forests (KODFOR 500-599
# or 701), 2/5/5 on the Oregon ones; TFALL(I,0)=1 for IC (6), WL (17), the hardwoods 21-31 and 33; TFALL(I,3)=20 for
# IC/RC/PY (6,18,20), 10 for ES/PP (8,10), else 15; TFALL(4:5)=TFALL(3), TFALL(0)<=LEAFLF, TFALL(2)<=TFALL(3).
@inline function _fm_tfall_so(coef, sp::Int, sz::Int, kodfor::Int)::Float32
    ca = (500 <= kodfor < 600) || kodfor == 701
    t0, t1, t2 = ca ? (3f0, 10f0, 15f0) : (2f0, 5f0, 5f0)
    (sp == 6 || sp == 17 || 21 <= sp <= 31 || sp == 33) && (t0 = 1f0)
    t3 = (sp == 6 || sp == 18 || sp == 20) ? 20f0 : (sp == 8 || sp == 10) ? 10f0 : 15f0
    t0 = min(t0, Float32(coef_col(coef, :leaf_life)[sp]))
    t2 > t3 && (t2 = t3)
    return sz == 0 ? t0 : sz == 1 ? t1 : sz == 2 ? t2 : t3
end
_fm_tfall_iestyle(v) = v isa InlandEmpire || v isa EasternMontana || v isa Kootenai
@inline function _fm_tfall_ie(coef, sp::Integer, sz::Int)::Float32
    t3 = coef_col(coef, :tfall_cls)[sp]
    sz == 0 && return min(2f0, coef_col(coef, :leaf_life)[sp])
    sz == 1 && return 5f0
    sz == 2 && return min(5f0, t3)
    return t3
end

"""
    fmscro!(s, sp, dbh, xv, density, dkcl)

Schedule a dying tree's crown debris into the CWD2B debris-in-waiting pool (FMSCRO, fmscro.f): each
crown component `xv[size]·density` is spread EQUALLY over years 1..min(TSOFT, TFALL(sp,size)), where
TSOFT = `(1.24·dbh + 13.82)·DECAYX` (FMSNGDK). The un-fallen CWD2B is the Stand-Dead crown; it flows
to the down-wood pool as it falls. `xv` is the `crown_biomass` tuple (foliage, woody1-5), in lb.
"""
function fmscro!(s::StandState, sp::Integer, dbh::Float32, xv, density::Float32, dkcl::Integer; icall::Int = 1)
    fs = s.fire; coef = s.coef
    cls = clamp(Int(coef_col(coef, :tfall_cls)[sp]), 1, 6)
    dcx = get(fs.params.snag_decayx_ovr, Int32(sp), coef_col(coef, :snag_decayx)[sp])  # SNAGDCAY override
    # FMSNGDK (fmscro.f:99): the R6 variants' JYRSOFT·DECAYX (_snag_dktime); others the (1.24·D+13.82)·DECAYX form.
    tsoft = r6_ffe_code(s.variant) === :none ? (1.24f0 * dbh + 13.82f0) * dcx : _snag_dktime(s, Int(sp), dbh, dcx)
    @inbounds for sz in 0:5
        amt = xv[sz + 1] * density
        amt > 0f0 || continue
        # ILIFE = ceil(RLIFE), floor 1 (fmscro.f:126-131: INT(RLIFE) then +1 if truncated or ≤0) — NOT round.
        # TFALL(SP,SIZE) (fmscro.f:124): EC/WC/PN/OP read their own fmvinit table; the western variants' CSV tfall_cls
        # holds TFALL(I,3) (10/15/20), which the SN class lookup clamps to row 6 (open for IE/EM/CR/BM/… too).
        tft = _fm_tfall_table(s.variant)
        tf = tft !== nothing ? tft[sp, sz + 1] : _fm_tfall_iestyle(s.variant) ? _fm_tfall_ie(coef, sp, sz) :
             s.variant isa SouthCentralOregon ? _fm_tfall_so(coef, Int(sp), sz, Int(s.plot.user_forest_code)) :
             _fm_tfall(cls, sz, sp)
        ilife = clamp(ceil(Int, min(tsoft, tf)), 1, 60)
        annual = amt / ilife
        # fmscro.f:160-170: mortality reconciliation (ICALL=4, FMKILL after the annual loop) books straight into CWD2B;
        # every other caller (fire, pile burn, cut/new snags) into CWD2B2, which FMMAIN merges after that year's FMCADD
        # — so a fire's crowns start falling the year AFTER the burn.
        c2 = icall == 4 ? fs.cwd2b : fs.cwd2b2
        for yr in 1:ilife
            c2[dkcl, sz + 1, yr] += annual
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
    (s.fire === nothing ? 0f0 : (sum(s.fire.cwd2b) + sum(s.fire.cwd2b2)) * _FM_P2T) * 0.5f0   # fmdout.f:173 CWD2B+CWD2B2

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

"""
    compute_crown_lift!(s, old_ht, old_dbh, old_cr, cyclen) -> StandState

Compute this cycle's crown-LIFT down-wood input (FMSDIT/FMCADD, fmsdit.f:103-119 + fmcadd.f:86-102) into
`fire.crown_lift_annual` (the per-YEAR addition, added each year by `ffe_fuel_update!`). As a live tree
grows its crown base rises; the lower woody crown left below the new base dies and falls to down wood at
`FMPROB · X · OLDCRW(size) · P2T` per year, where `X = crown_lift_rate(oldht, oldcrl, ht, ICR, cyclen)`
and `OLDCRW(size)` is the PREVIOUS cycle's woody crown weight (sizes 1-5; foliage size-0 is excluded —
leaf-lifespan already books it). Both FVS thresholds (raw OLDCRW and FMPROB·OLDCRW < 0.0000625) are
applied. `old_*` are the per-tree previous-cycle height / DBH / crown-% snapshots, aligned by record to
the current tree list; if the list changed (regen/compaction) the term is skipped (left zero) — faithful
only while records are stable (the general case needs OLDHT to travel with the record, FMOLDC). This is
the DOMINANT post-mortality down-wood source (~2.5 t/ac/cycle on carbon_snt, instrumented vs Fortran).
"""
function compute_crown_lift!(s::StandState, cyclen::Real)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    cl = fs.crown_lift_annual; fill!(cl, 0f0)
    t = s.trees; coef = s.coef; dkrcls = coef_col(coef, :dkr_cls)
    ocrw = t.ffe_oldcrw
    @inbounds for i in 1:t.n
        for k in 1:5; ocrw[k, i] = 0f0; end                # reset this record's OLDCRW (recomputed below)
        t.tpa[i] > 0f0 || continue
        oldht = t.ffe_oldht[i]
        oldht > 0f0 || continue                            # prev-cycle state not set (1st cycle / regen)
        sp = Int(t.species[i])
        oldcr = t.ffe_oldcr[i]
        oldcrl = oldht * oldcr / 100f0
        x = crown_lift_rate(oldht, oldcrl, t.height[i], Float32(t.crown_pct[i]), cyclen)
        x > 0f0 || continue
        # OLDCRW = the PREVIOUS-cycle woody crown weights (recomputed from the old tree state, = FMOLDC)
        # FMOLDC (fmoldc.f, one file in every build) saved the record's CROWNW itself (fmsdit.f:106 OLDCRW = X·OLDCRW) —
        # the array ffe_fmcrow! filled at the last FMSDIT with its PCTILE height percentile (and any fire reduction), not a
        # recompute from the old dims (only a stand FMCROW never ran on falls back to that).
        xvold = fs.fmcrow_on ? ntuple(k -> t.ffe_crownw[k, i], 6) :
                crown_biomass(s, sp, t.ffe_olddbh[i], oldht, Int(round(oldcr)))
        dkcl = clamp(Int(dkrcls[sp]), 1, 4)
        for sz in 1:5
            ocw = xvold[sz + 1]
            ocw < 0.0000625f0 && continue                  # FMSDIT raw-OLDCRW threshold
            lift = x * ocw                                 # OLDCRW after the X scaling
            ocrw[sz, i] = lift                             # store per-tree OLDCRW (for the at-death FMSCRO term)
            t.tpa[i] * lift < 0.0000625f0 && continue      # FMCADD FMPROB·OLDCRW threshold (down-wood only)
            cl[sz, dkcl] += t.tpa[i] * lift * _FM_P2T
        end
    end
    return s
end

"""
    crown_lift_at_death(t, i, cyclen) -> NTuple{6,Float32}

The crown-lift the dying record `i` would have shed this cycle, that instead joins its snag crown at death
(FVS FMSCRO, fmscro.f:147: `ANNUAL += YRSCYC·OLDCRW(SIZE)·X`, X=1 for a kill) — `YRSCYC · OLDCRW(size)`,
read from the per-tree `ffe_oldcrw` (the X-scaled crown-lift stored by `compute_crown_lift!` last cycle).
Sizes 1–5 only (foliage size-0 excluded). Returns the per-size addition to the dying record's `crown_biomass`.
"""
@inline function crown_lift_at_death(t::TreeList, i::Integer, cyclen::Real)::NTuple{6,Float32}
    yrs = Float32(cyclen)
    @inbounds ntuple(sz -> sz == 1 ? 0f0 : yrs * t.ffe_oldcrw[sz - 1, i], 6)
end

"""
    snapshot_ffe_oldcrown!(s) -> StandState

Store each live record's current height / DBH / crown-% into its `ffe_old*` fields (FMOLDC): the
previous-cycle crown state that next cycle's `compute_crown_lift!` reads to size the crown-base rise.
Called at the END of a cycle's fuel processing; the values then travel through the next grow's record
tripling/compaction (`copy_tree!`), keeping the crown-lift aligned. No-op without FFE.
"""
function snapshot_ffe_oldcrown!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    t = s.trees
    @inbounds for i in 1:t.n
        t.ffe_oldht[i]  = t.height[i]
        t.ffe_olddbh[i] = t.dbh[i]
        t.ffe_oldcr[i]  = Float32(t.crown_pct[i])
    end
    return s
end

function ffe_fuel_update!(s::StandState, nyrs::Integer; vtrip::Bool = false)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    fmcba!(s; vtrip = vtrip)
    # FVS FMMAIN year-loop ORDER: FMSNAG (snag fall → bole into down wood) → FMCWD (decay) → FMCADD
    # (cwd2b crown fall + litterfall + woody breakage + crown-lift). The snag falldown MUST precede the
    # decay so the freshly-fallen bole is decayed in the same year it falls (else it over-accumulates by
    # ~a cycle's worth of decay — the DDW size-4/5 overshoot). The snag DENSITY falldown is identical
    # whether stepped 1yr×nyrs here or nyrs at once (it compounds the same), so Stand-Dead is unchanged.
    cur0 = Int(current_cycle_year(s))
    @inbounds for k in 1:nyrs
        # FMSNAG: snag bole → down wood (this year). Pass the ACTUAL annual year so a fire snag created
        # this cycle (before this loop) ages across the loop and falls in the years after the burn —
        # ordinary-mortality snags are created after the loop, so they're absent this cycle regardless.
        isempty(fs.snags.sp) || update_snags!(s, 1; at_year = cur0 + (k - 1))
        ffe_snag_height_loss!(s, 1; at_year = cur0 + (k - 1))   # SNAGBRK bole breakage (no-op unless HTX set)
        fmcwd!(s, 1)                                   # FMCWD: decay (now also decays this year's bole)
        fmcadd!(s; vtrip = vtrip)                      # FMCADD: litterfall, breakage, crown lift, then CWD2B year-1 fall
        fs.cwd2b .+= fs.cwd2b2; fill!(fs.cwd2b2, 0f0)  # fmmain.f:243-257 CWD2B += CWD2B2; CWD2B2 = 0
    end
    fs.bioroot *= fpowi(1f0 - _FM_CRDCAY, nyrs)    # dead-root decay (fmcrbout.f:273, REAL**INTEGER ⇒ __powisf2)
    return s
end

"""
    fmcadd!(s) -> StandState

One year of FMCADD (fmcadd.f, the same file in every variant build), in its own order and REAL*4 association: per
live record I (FMPROB>0) the foliage litterfall (CROWNW(I,0)·FMPROB/LEAFLF)·P2T, then for SIZE 1..5 the limb breakage
(LIMBRK·FMPROB·CROWNW(I,SIZE))·P2T and the crown lift (FMPROB·OLDCRW(I,SIZE))·P2T (0 when FMPROB·OLDCRW < 6.25E-5),
all into the hard pool of DKRCLS(SP); then the year-1 CWD2B crown-debris slot falls (DOWN/2000) and the slots shift.
The crown-lift term reads the CURRENT FMPROB and OLDCRW each year, so a fire's kill and its OLDCRW(I,1) halving
(fmeff.f:498-503) shrink it at once — jl added a cycle-start total (MEASURED FVSie_g16 4769882010690 SIMFIRE 2014
FMCADD: hard <0.25" +0.0183 live vs +0.0398 jl).
"""
function fmcadd!(s::StandState; vtrip::Bool = false)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    t = s.trees; coef = s.coef
    leaflf = coef_col(coef, :leaf_life); ocrw = t.ffe_oldcrw; cwd = fs.cwd
    # DO I=1,ITRN over FMMAIN's list: in a tripling cycle that is the tripled list (originals ×.60, then each record's
    # ×.25/×.15 copies with the same CROWNW/OLDCRW), and jl's is still untripled here — walk it the FVS way.
    _fm_record_walk(t, vtrip) do i, p
        p > 0f0 || return
        sp = Int(t.species[i])
        dkcl = clamp(ffe_dkr_cls(s, sp), 1, 4)
        xv = _ffe_crownw(s, i, sp, t.dbh[i], t.height[i], Int(round(t.crown_pct[i])))
        ll = leaflf[sp]
        @inbounds begin
            ll > 0f0 && (cwd[10, 2, dkcl] += (xv[1] * p / ll) * _FM_P2T)
            for sz in 1:5
                cwd[sz, 2, dkcl] += (_FM_LIMBRK * p * xv[sz + 1]) * _FM_P2T
                po = p * ocrw[sz, i]
                po < 0.0000625f0 || (cwd[sz, 2, dkcl] += po * _FM_P2T)
            end
        end
    end
    c2 = fs.cwd2b
    @inbounds for dkcl in 1:4
        for sz in 0:5
            down = c2[dkcl, sz + 1, 1]                 # PDOWN·CWD2B(DKCL,SIZE,1)/NYRS, NYRS=1 in the annual loop
            down > 0f0 && (cwd[sz == 0 ? 10 : sz, 2, dkcl] += down / 2000f0)
        end
    end
    @inbounds for dkcl in 1:4, sz in 1:6
        for yr in 1:59
            c2[dkcl, sz, yr] = c2[dkcl, sz, yr + 1]
        end
        c2[dkcl, sz, 60] = 0f0
    end
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
        xv = _ffe_crownw(s, i, sp, t.dbh[i], t.height[i], Int(round(t.crown_pct[i])))
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
        xv = _ffe_crownw(s, i, sp, t.dbh[i], t.height[i], Int(round(t.crown_pct[i])))
        dkcl = clamp(Int(dkrcls[sp]), 1, 4)
        for sz in 1:5
            fs.cwd[sz, 2, dkcl] += _FM_LIMBRK * xv[sz + 1] * t.tpa[i] * _FM_P2T
        end
    end
    return s
end



"""
    apply_pileburn!(s) -> Bool

PILEBURN (act 2523, fmtret.f): a scheduled jackpot/pile burn. The NET effect on the down-wood pools (jl
does not model FVS's transient piled/unpiled CWD dimension) is to CONSUME the staged fraction of each fuel
size class — size 1-9: AFFECT·FULCON, litter/duff (10-11): AFFECT·ATREAT — times the FMCONS consumption
fraction at medium moisture (FMOIS=3). Optionally kills TRMORT of every tree's TPA (→ snags + crown debris,
the same booking as a fire kill). Params: type / AFFECT% / ATREAT% / FULCON% / TRMORT%. No-op without a due
PILEBURN.
"""
function apply_pileburn!(s::StandState)::Bool
    fs = s.fire
    (fs === nothing || !fs.active || isempty(s.control.schedule)) && return false
    yr = Int(current_cycle_year(s)); fvscyc = Int(s.control.cycle) + 1
    fired = false
    for a in s.control.schedule
        a.icflag == Int32(2523) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        affect = a.params[2] / 100f0; atreat = a.params[3] / 100f0
        fulcon = a.params[4] / 100f0; trmort = clamp(a.params[5] / 100f0, 0f0, 1f0)
        # consume the staged (piled) fraction of each fuel size class (FMCONS net, FMOIS=3 medium)
        mois = fuel_moisture(3, s.variant); fr = fire_consumption_fractions(mois); cwd = fs.cwd
        @inbounds for sz in 1:11
            stage = sz <= 9 ? affect * fulcon : affect * atreat
            cf = stage * fr[sz]
            cf > 0f0 || continue
            for k in 1:2, l in 1:4; cwd[sz, k, l] *= (1f0 - cf); end
        end
        # optional uniform tree mortality → snags + crown debris (mirrors the fmburn! fire kill→snag path)
        if trmort > 0f0
            t = s.trees; coef = s.coef
            pend = Tuple{Int,Float32,Float32,Float32,Float32,Float32,Float32}[]
            @inbounds for i in 1:t.n
                t.tpa[i] > 0f0 || continue
                trkil = t.tpa[i] * trmort; t.tpa[i] -= trkil
                sp = Int(t.species[i]); d = t.dbh[i]
                push!(pend, (sp, d, t.height[i], t.height[i], t.height[i], trkil, -1f0))   # FMSSEE (fmtret.f:154)
                xvc = _ffe_crownw(s, i, sp, d, t.height[i], Int(t.crown_pct[i]))
                fmscro!(s, sp, d, xvc, trkil, clamp(ffe_dkr_cls(s, sp), 1, 4))           # fmtret.f:160
                _, _, rbio = jenkins_biomass(coef, sp, d)
                fs.bioroot += rbio * trkil                                               # fmsadd.f:320 (XDCAY 1)
            end
            # FMSADD(IYR,1) (fmtret.f:166): the killed trees are binned into snag records like any other source
            isempty(pend) || fmsadd_bin!(s, pend, yr; ityp = 1, bolefn = _r6_snag_bolefn(s))
            compute_density!(s)
        end
        fired = true
    end
    return fired
end
