# =============================================================================
# fire/carbon.jl — standing live-tree carbon pools (FFE chunk F8 — carbon)
#
# Ported from: bin/FVSsn_buildDir/fmcrbout.f (the live-tree pools of the FFE Stand
# Carbon Report), which sums the Jenkins biomass (FMCBIO, ported as `jenkins_biomass`)
# over the tree list and converts biomass → carbon at 0.5 (fmcrbout.f:89/158).
#
# These are the live aboveground / merchantable / belowground carbon pools (tons C/ac).
# The dead pools (snags, down wood, forest floor) build on the F7 snag/CWD model; this
# is the live-tree foundation. Carbon appears in the DBS Carbon report, not the `.sum`.
# =============================================================================

"""
    stand_live_carbon(s) -> (; aboveground, merch, belowground)

Standing live-tree carbon pools in tons C/acre (FFE Stand Carbon Report, fmcrbout.f):
the per-tree Jenkins aboveground / merchantable / belowground (root) biomass summed over
the tree list (weighted by TPA) and converted to carbon at the 0.5 biomass→carbon ratio.
"""
# FMCBIO's merch gate DBHMIN(KSP): the variant's cubic merch standard once set, else the species CSV value.
_jenkins_dbhmin(s::StandState, sp::Int) = _fm_dbhmin(s, sp)

function stand_live_carbon(s::StandState)
    t = s.trees; coef = s.coef
    above = 0f0; merch = 0f0; root = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        a, m, r = jenkins_biomass(coef, t.species[i], t.dbh[i]; dbhmin = _jenkins_dbhmin(s, Int(t.species[i])))
        above += a * t.tpa[i]
        merch += m * t.tpa[i]
        root  += r * t.tpa[i]
    end
    return (aboveground = above * 0.5f0, merch = merch * 0.5f0, belowground = root * 0.5f0)
end

"""
    ffe_down_wood(s) -> (; vol_hard, vol_soft, cov_hard, cov_soft)

Down-wood VOLUME (cuft/ac) and percent COVER for the FVS_Down_Wood_Vol / FVS_Down_Wood_Cov DBS tables
(fmdout.f:312-373). Volume = `cwd`(biomass tons/ac)·2000/CWDDEN, with CWDDEN = 18.72 (soft) / 24.96
(hard) lbs/cuft (SG 0.3/0.4). `vol_*` are 8-tuples: the 7 DBH bins (0-3=sizes 1-3, then 3-6/6-12/12-20/
20-35/35-50/≥50 = size classes 4-9) + total. Cover = `a·vol^b` per size (sizes 1-3 have 0 cover); `cov_*`
are 7-tuples (6 bins 3-6…≥50 + total). All derived from the same validated `cwd` pools.
"""
function ffe_down_wood(s::StandState)
    z8 = ntuple(_ -> 0f0, 8); z7 = ntuple(_ -> 0f0, 7)
    fs = s.fire
    (fs === nothing || !fs.active) && return (vol_hard = z8, vol_soft = z8, cov_hard = z7, cov_soft = z7)
    cw = fs.cwd; den = (18.72f0, 24.96f0)                    # CWDDEN by cwd hardness index (1=soft, 2=hard)
    # fmdout.f:316-343: CWDVOL(I,J,K,L) = CWD(I,J,K,L)·2000/CWDDEN per decay class, then summed over I (piles) and L —
    # the volume of each pool, not the volume of the pooled biomass.
    vol(sz, K) = (let v = 0f0; for L in 1:4; v += (cw[sz, K, L] * 2000f0) / den[K]; end; v end)
    function vbins(K)
        b = (vol(1, K) + vol(2, K) + vol(3, K), vol(4, K), vol(5, K), vol(6, K), vol(7, K), vol(8, K), vol(9, K))
        (b..., sum(b))
    end
    # cover power-law (a, b) for size classes 4-9 (3-6 … ≥50 in); sizes 1-3 contribute 0 (fmdout.f:362-376)
    ccf = ((4, 0.0166f0, 0.8715f0), (5, 0.0092f0, 0.8795f0), (6, 0.0063f0, 0.8728f0),
           (7, 0.0069f0, 0.8134f0), (8, 0.0033f0, 0.8617f0), (9, 0.0949f0, 0.5f0))
    function cbins(K)
        c = ntuple(i -> (let (sz, a, b) = ccf[i]; a * vol(sz, K)^b end), 6)
        (c..., sum(c))
    end
    return (vol_hard = vbins(2), vol_soft = vbins(1), cov_hard = cbins(2), cov_soft = cbins(1))
end

"""
    ffe_fuel_loadings(s) -> NamedTuple

The FFE fuel loadings (tons/ac **biomass**, NOT carbon — no ×0.5) that feed the FVS_Fuels DBS table
(DBSFUELS, fmdout.f:399). Surface pools come from `fire.cwd` (litter=10, duff=11, woody by size class
1-9) and `fire.flive` (herb/shrub) — the same pools that give the validated DDW / Forest-Floor. Standing
pools: snag biomass split by DBH ≤3/>3 (bole `bolevol·density` + the CWD2B crown, as in `standing_dead`)
and live biomass (foliage + woody crown + stem `_fm_cuft·v2t`, split by tree DBH ≤3/>3, = the CARBCALC=0
`BIOLIVE` components, fmdout.f:218-258). Consumed / removed are 0 without a fire / harvest this cycle.
"""
function ffe_fuel_loadings(s::StandState; vtrip::Bool = false)
    fs = s.fire
    z = (litter=0f0, duff=0f0, lt3=0f0, ge3=0f0, s3to6=0f0, s6to12=0f0, ge12=0f0, herb=0f0, shrub=0f0,
         surf_total=0f0, snag_lt3=0f0, snag_ge3=0f0, foliage=0f0, live_lt3=0f0, live_ge3=0f0,
         stand_total=0f0, total_biomass=0f0, consumed=0f0, removed=0f0)
    (fs === nothing || !fs.active) && return z
    cw = fs.cwd
    sumc(rng) = sum(@view cw[rng, :, :])
    litter = sumc(10:10); duff = sumc(11:11)
    lt3 = sumc(1:3); s3to6 = sumc(4:4); s6to12 = sumc(5:5); ge12 = sumc(6:9); ge3 = s3to6 + s6to12 + ge12
    herb = fs.flive[1]; shrub = fs.flive[2]
    surf_total = litter + duff + lt3 + ge3 + herb + shrub
    # standing snags (TOTSNG(1|2): bole by DBHS + the CWD2B/CWD2B2 crowns) and live trees (TOTFOL, TOTLIV(1|2)) are
    # FMDOUT's accumulators (fmdout.f:132-258), shared with the carbon report — `fmdout_bio`, on FMMAIN's record list
    fb = fmdout_bio(s; vtrip = vtrip)
    snag_lt3 = fb.totsng1; snag_ge3 = fb.totsng2; foliage = fb.totfol; live_lt3 = fb.totliv1; live_ge3 = fb.totliv2
    stand_total = live_lt3 + live_ge3 + foliage + snag_lt3 + snag_ge3        # TOTSTD (fmdout.f:261)
    return (; litter, duff, lt3, ge3, s3to6, s6to12, ge12, herb, shrub, surf_total,
            snag_lt3, snag_ge3, foliage, live_lt3, live_ge3, stand_total,
            total_biomass = surf_total + stand_total, consumed = 0f0, removed = 0f0)
end

# FFE live-tree stem merch cubic (NATCRS MCF = v[4]+v[7]) as FMSVL2/FMDOUT compute it (fmsvol.f:130-150).
# The common path returns the tree's cached `merch_cuft_vol`. The ONLY exception is a broken-top tree on
# the SN R8-Clark volume path: there the .sum's `merch_cuft_vol` was built from the NORMAL height (norm_ht)
# + CFTOPK truncation, but FMSVL2 for a LIVE tree calls NATCRS with the ACTUAL (broken) height as the total
# height and LTKIL=.FALSE. (no top-kill), giving a different merch cubic. Recompute it there to match FVS.
# Gated to the Southern R8-Clark path (NE/CS/LS use a separate NVEL routine and are out of scope / already
# matched); a NON-broken tree is never recomputed, preserving the bit-exact `merch_cuft_vol` for 299/300.
@inline function _ffe_stem_mcf(s::StandState, i::Int, sp::Int, d::Float32, h::Float32)::Float32
    (s.variant isa Southern && h >= 4.5f0 && s.trees.trunc[i] > 0) || return s.trees.merch_cuft_vol[i]
    c = s.control
    if d >= c.sp_scf_dbhmin[sp]
        prod = "01"; stump = c.sp_scf_stump[sp]; mtopp = c.sp_scf_topd[sp]
    else
        prod = "02"; stump = c.sp_stump_ht[sp];  mtopp = c.sp_top_diam[sp]
    end
    v, _, _ = _R8CLARK_VOL(s.species.vol_eq[sp], d, h, mtopp, c.sp_top_diam[sp], stump, prod)
    return d >= c.sp_dbh_min[sp] ? v[4] + v[7] : 0f0     # NATCRS MCF; no CFTOPK (LTKIL=.FALSE. for live)
end

"""
    ffe_live_carbon(s) -> (; aboveground, merch)

Live aboveground / merchantable carbon by the **FFE-fuel** method (CARBCALC=0; fmcrbout.f:120-141 +
fmdout.f:225-258), in tons C/acre. Aboveground `BIOLIVE` = the FFE crown biomass (foliage + woody sizes
1-5, `crown_biomass`, lb→tons via P2T) **plus** the stem biomass (`_fm_cuft` cubic volume × `v2t`);
merch = the stem biomass alone. Both × 0.5 for carbon. Belowground (roots) stays the Jenkins value in
both methods (fmcrbout.f:144-146), so it is not recomputed here.

NB the OLDCRW crown-lift term that fmdout.f adds to BIOLIVE is `X·CROWNW` with `X` the per-year crown-
base-rise fraction (~7e-4/yr) — i.e. <0.1% of the crown — so it is omitted here (negligible, and the
live FFE oracle is unavailable in the stripped validation binary; see FFE_FUEL_DYNAMICS_chunk_plan.md).
"""
function ffe_live_carbon(s::StandState)
    t = s.trees; coef = s.coef; v2t = coef_col(coef, :v2t)
    above = 0f0; merch = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        xv = _ffe_crownw(s, i, sp, d, h, Int(round(t.crown_pct[i])))   # (foliage, woody 1..5), lb
        crown = xv[1]; for sz in 1:5; crown += xv[sz + 1]; end         # foliage + all woody (lb)
        # Stem volume = FMSVL2 with LMERCH=.FALSE., which for SN (VARACD∈{CS,LS,NE,SN}) returns MAX(X,MCF)
        # (fmsvol.f:149-151), where X = 0.005454154·H is the tiny-tree cone floor — NOT gross/TCF. SN's MCF is the
        # NATCRS merch cubic (v[4]+v[7]). For a NON-broken tree this equals the tree's `merch_cuft_vol`, so that
        # is used directly (the common path). For a BROKEN-TOP tree the two DIVERGE: FMSVL2/FMDOUT calls NATCRS
        # with the tree's ACTUAL (broken) height as total height and LTKIL=.FALSE. (no CFTOPK), whereas the .sum's
        # `merch_cuft_vol` builds the profile from the NORMAL height (norm_ht) and then truncates via CFTOPK
        # (vols.f:191-193). These give different merch cubics — e.g. carbon_ffe sp22 D10.4 H55 (norm 64.77):
        # FMSVL2 MCF = 11.2 (verified live via DEBUG FMDOUT), .sum merch_cuft_vol = 13.2. So for broken-top trees
        # (SN R8-Clark path) the stem merch is RECOMPUTED at the actual height to match FMSVL2; non-broken trees
        # (299/300) keep merch_cuft_vol bit-exact. See docs/TOLERANCE_AUDIT.md 2026-07-05y.
        if _ffe_west_vol(s.variant)
            # Western FMSVOL ({v}/fmsvol.f:146-151): BIOLIVE's stem is FMSVL2 LMERCH=.FALSE. ⇒ MAX(X,TCF) (fmdout.f:247)
            # while the merch pool uses LMERCH=LVWEST=.TRUE. ⇒ MCF (fmcrbout.f:127) — both NATCRS on the actual
            # height, no CFTOPK. jl used the merch cubic for both (EM Aboveground_Total_Live 11.67 vs live 14.70;
            # IE 12.04 vs 15.57 on S248112 1990).
            tcf, mcf, _, _ = ffe_west_nocut(s, sp, d, h)
            above += t.tpa[i] * (crown + max(0.005454154f0 * h, tcf) * v2t[sp]) * _FM_P2T
            merch += t.tpa[i] * (mcf * v2t[sp]) * _FM_P2T
            continue
        end
        mcf = _ffe_stem_mcf(s, i, sp, d, h)
        stem = max(0.005454154f0 * h, mcf) * v2t[sp]                  # MAX(X,MCF) × V2T = stem biomass (lb)
        above += t.tpa[i] * (crown + stem) * _FM_P2T                   # BIOLIVE = crown + stem, lb→tons
        merch += t.tpa[i] * stem * _FM_P2T                             # merch = stem only (FMSVL2·V2T/2000)
    end
    return (aboveground = above * 0.5f0, merch = merch * 0.5f0)        # biomass → carbon (×0.5)
end

"""
    standing_dead_carbon(s) -> Float32

Standing-dead (snag) carbon pool in tons C/acre: the Jenkins aboveground biomass of each
snag cohort (`fire.snags`, F7) weighted by its still-standing density (hard + soft),
converted to carbon at 0.5 (fmcrbout.f). Zero when FFE is off / no snags.
"""
function standing_dead_carbon(s::StandState)::Float32
    # Stand-Dead = snag stem-volume BOLE (SNVIS·V2T) + the crown debris still in CWD2B (fmdout.f:153/173).
    return snag_bole_carbon(s) + snag_crown_carbon(s)
end

"""
    down_wood_carbon(s) -> Float32

Down dead wood carbon pool in tons C/acre: the 9 woody surface-fuel size classes of `fire.cwd`
(the down-wood loadings in tons/ac, F3) converted to carbon at 0.5 (fmcrbout.f BIODDW). Litter
and duff are the SEPARATE forest-floor pool (`forest_floor_carbon`). Zero when FFE is off.
"""
function down_wood_carbon(s::StandState)::Float32
    fs = s.fire; fs === nothing && return 0f0
    return sum(@view fs.cwd[1:9, :, :]) * 0.5f0
end

"""
    forest_floor_carbon(s) -> Float32

Forest-floor carbon pool in tons C/acre: litter + duff (the last two `fire.cwd` size classes)
converted to carbon at **0.37** (Smith & Heath, NE-722; fmcrbout.f:90/160 — the litter/duff
carbon fraction differs from the 0.5 used for all woody/live pools). Zero when FFE is off.
"""
function forest_floor_carbon(s::StandState)::Float32
    fs = s.fire; fs === nothing && return 0f0
    return sum(@view fs.cwd[10:11, :, :]) * 0.37f0
end

"""
    belowground_dead_carbon(s) -> Float32

Dead coarse-root carbon pool in tons C/acre: the FFE `BIOROOT` accumulator (Jenkins root biomass of
trees as they die, decayed at CRDCAY each cycle; fmsadd.f:320 / fmcrbout.f:273) × 0.5. Zero when FFE
is off / nothing has died.
"""
belowground_dead_carbon(s::StandState)::Float32 =
    (s.fire === nothing ? 0f0 : s.fire.bioroot) * 0.5f0

"""
    shrub_herb_carbon(s) -> Float32

Live shrub + herb carbon pool in tons C/acre: `BIOSHRB = FLIVE(1) + FLIVE(2)` (the FFE live
surface-fuel loadings, fmdout.f:283) converted to carbon at 0.5. Zero when FFE is off.
"""
function shrub_herb_carbon(s::StandState)::Float32
    fs = s.fire; fs === nothing && return 0f0
    return (fs.flive[1] + fs.flive[2]) * 0.5f0
end

"""
    stand_carbon(s) -> (; live_above, live_below, standing_dead, down_wood, forest_floor, shrub_herb, total)

All the main stand carbon pools (tons C/acre): live aboveground + belowground (trees), standing
dead (snags), down dead wood, forest floor (litter+duff, at the 0.37 fraction), and live shrub+herb
(FFE Stand Carbon Report, fmcrbout.f). The FFE pools are zero unless `fmcba!` has run this cycle.
"""
function stand_carbon(s::StandState)
    lc = stand_live_carbon(s)
    sd = standing_dead_carbon(s)
    dw = down_wood_carbon(s)
    ff = forest_floor_carbon(s)
    sh = shrub_herb_carbon(s)
    return (; live_above = lc.aboveground, live_below = lc.belowground,
            standing_dead = sd, down_wood = dw, forest_floor = ff, shrub_herb = sh,
            total = lc.aboveground + lc.belowground + sd + dw + ff + sh)
end

# short tons/acre → metric tons/hectare (the Stand Carbon Report's units, fmcrbout.f METRIC.F77).
const _TONAC_TO_MTHA = 0.90718474f0 / 0.40468564f0

# FAPROP — harvested-wood-products fate fractions [habitat 1:2, year-since-harvest 1:101, fate 1:3
# (in-use/landfill/energy), product 1:2 (pulp/saw), group 1:2 (sw/hw)] (fmcblk.f, Smith et al.). Loaded
# once from data/southern/fire_hwp_fate.csv. Const lookup table (not mutable state).
const _FAPROP = let
    path = joinpath(@__DIR__, "..", "..", "..", "data", "southern", "fire_hwp_fate.csv")
    a = zeros(Float32, 2, 101, 3, 2, 2)
    fate_i = Dict("in_use" => 1, "landfill" => 2, "energy" => 3)
    prod_i = Dict("pulpwood" => 1, "sawtimber" => 2)
    grp_i  = Dict("softwood" => 1, "hardwood" => 2)
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue                          # header
        f = split(strip(line), ',')
        isempty(f[1]) && continue
        a[parse(Int, f[1]), parse(Int, f[5]), fate_i[f[2]], prod_i[f[3]], grp_i[f[4]]] = parse(Float32, f[6])
    end
    a
end

"""
    accrue_harvest_carbon!(s, sp, dbh, removed_tpa, merch_biomass, year)

Bucket one cut tree's harvested merch biomass into the FFE harvested-wood-products FATE accumulator
(FMSCUT, fmscut.f:147-151): product = pulpwood if DBH ≤ `CDBRK` else sawtimber (CDBRK = 9″ softwood /
11″ hardwood), group = softwood if `biogrp ≤ 5` else hardwood. `merch_biomass` is the per-acre harvested
merch biomass (tons/ac, merch-cuft × V2T × removed TPA, or Jenkins merch × TPA). No-op without FFE.
"""
function accrue_harvest_carbon!(s::StandState, sp::Integer, dbh::Float32,
                                merch_biomass::Float32, year::Integer)
    fs = s.fire
    (fs === nothing || !fs.active || merch_biomass <= 0f0) && return
    grp = Int(coef_col(s.coef, :biogrp)[sp]) > 5 ? 2 : 1     # 1 = softwood, 2 = hardwood
    cdbrk = grp == 1 ? 9f0 : 11f0                            # pulp/saw DBH break (fminit.f:919-920)
    prod = dbh > cdbrk ? 2 : 1                               # 1 = pulpwood, 2 = sawtimber
    key = (Int(year), prod, grp)
    fs.hwp_fate[key] = get(fs.hwp_fate, key, 0f0) + merch_biomass
    return
end

"""
    harvested_carbon_report(s, year, habitat) -> (; products, landfill, energy, emissions, stored, removed)

The harvested-wood-products carbon pools (metric tons C/ha) at report `year` (FMCHRVOUT, fmchrvout.f:83-103).
For each past cut, the harvested merch biomass is distributed by the FAPROP year-since-harvest fate curves
into In-use (`products`), Landfill, Energy, and Emissions (= the residual `1 − Σ fate`). `stored` = products
+ landfill; `removed` = energy + emissions + stored. All × 0.5 carbon, × the tons/ac→t/ha factor.
"""
function harvested_carbon_report(s::StandState, year::Integer, habitat::Integer)
    fs = s.fire
    (fs === nothing || !fs.active) && return (; products = 0f0, landfill = 0f0, energy = 0f0,
                                              emissions = 0f0, stored = 0f0, removed = 0f0)
    h = clamp(Int(habitat), 1, 2)
    v = zeros(Float64, 4)                                    # in-use, landfill, energy, emissions
    for ((cyr, prod, grp), bio) in fs.hwp_fate
        kyr = year - cyr + 1; kyr < 1 && continue
        kyr > 101 && (kyr = 101)
        xtmp = 0.0
        for fate in 1:3
            f = _FAPROP[h, kyr, fate, prod, grp]
            xtmp += f
            v[fate] += bio * f
        end
        v[4] += bio * (1.0 - xtmp)                           # emissions = the un-accounted residual
    end
    k = 0.5f0 * _TONAC_TO_MTHA                               # biomass → carbon → metric t/ha
    products = Float32(v[1]) * k; landfill = Float32(v[2]) * k
    energy = Float32(v[3]) * k; emissions = Float32(v[4]) * k
    stored = products + landfill                             # V(5) = in-use + landfill
    removed = energy + emissions + stored                   # V(6) = energy + emissions + stored
    return (; products, landfill, energy, emissions, stored, removed)
end

# FMMAIN (gradd.f:118) runs after grincr.f:543's TRIPLE, so in a tripling cycle FMDOUT/FMCRBOUT loop over the tripled list:
# records 1..ITRN at PROB·0.60, then per record the copies ITRN+2I-1 (PROB·0.25) and ITRN+2I (PROB·0.15) (triple.f) — same
# D/H/ICR, different FMPROB. jl samples the cycle-top report before its own triple_records!, so `vtrip` replays that
# record walk (the fire-cycle sample runs after jl's triple and passes vtrip=false).
@inline function _fm_record_walk(f, t::TreeList, vtrip::Bool)
    n = t.n
    if vtrip
        @inbounds for i in 1:n; f(i, t.tpa[i] * 0.60f0); end
        @inbounds for i in 1:n; f(i, t.tpa[i] * 0.25f0); f(i, t.tpa[i] * 0.15f0); end
    else
        @inbounds for i in 1:n; f(i, t.tpa[i]); end
    end
    return nothing
end

# Is the cycle about to start a TRIPLE (grincr.f:74 LTRIP=(ICYC.LE.ICL4 .AND. ITRN.LE.MAXTRE/3 .AND. .NOT.NOTRIP))?
_fm_will_triple(s::StandState) = !s.control.no_tripling && Int(s.control.cycle) < Int(s.control.icl4) &&
                                 s.trees.n > 0 && s.trees.n <= (variant_maxtre(s.variant) - Int(s.trees.ndead)) ÷ 3

"""
    fmdout_bio(s; vtrip=false) -> NamedTuple

fmdout.f:98-286, the accumulators FMCRBOUT reads, in FVS's own loop order and REAL*4 arithmetic:
SMALL2/LARGE2 over CWD(3,J,K,5) (= Σ_I Σ_L CWD(I,J,K,L), no piles here), TOTLIT/TOTDUF, the snag stems then the CWD2B/
CWD2B2 crown debris (DO ISZ / IDC / ITM) into TOTSNG(1|2), and per live record TOTFOL, TOTLIV(1)/(2) crowns
((CROWNW+OLDCRW)·P2T·FMPROB, P2T·FMPROB·(CROWNW+OLDCRW)) plus the FMSVL2 stem FMPROB·VT·V2T (V2T already /2000,
fmvinit.f:481). BIOLIVE=TOTFOL+TOTLIV(1)+TOTLIV(2), BIOSNAG=TOTSNG(1)+TOTSNG(2), BIODDW=SMALL2+LARGE2,
BIOFLR=TOTLIT+TOTDUF, BIOSHRB=FLIVE(1)+FLIVE(2).
"""
function fmdout_bio(s::StandState; vtrip::Bool = false)
    fs = s.fire; t = s.trees; coef = s.coef
    cw = fs.cwd
    c35(j, k) = (a = 0f0; @inbounds(for l in 1:4; a += cw[j, k, l]; end); a)
    small2 = 0f0; large2 = 0f0
    for isz in 1:3
        jsz = isz + 3; ksz = isz + 6
        small2 = small2 + c35(isz, 1) + c35(isz, 2)
        large2 = large2 + c35(jsz, 1) + c35(jsz, 2) + c35(ksz, 1) + c35(ksz, 2)
    end
    totduf = c35(11, 1) + c35(11, 2); totlit = c35(10, 1) + c35(10, 2)
    sn = fs.snags; tsng1 = 0f0; tsng2 = 0f0
    west = _ffe_west_vol(s.variant)
    @inbounds for i in eachindex(sn.sp)
        (sn.den_hard[i] + sn.den_soft[i]) > 0f0 || continue
        if west && sn.height[i] > 0f0
            # fmdout.f:139-155: SNVIH = FMSVOL(I,HTIH)·DENIH, SNVIS = FMSVOL(I,HTIS)·DENIS — the TOTAL cubic of (DBHS,HTDEAD)
            # with CFTOPK at the current height on every snag (XHT>−1 ⇒ LTKIL), fresh each report — then (SNVIS+SNVIH)·V2T
            sp = Int(sn.sp[i])
            vv = ffe_west_snag_vol_at(s, sp, sn.dbh[i], sn.height[i], sn.htcur[i]; always = true)
            snvih = sn.den_hard[i] > 0f0 ? vv * sn.den_hard[i] : 0f0
            snvis = sn.den_soft[i] > 0f0 ? vv * sn.den_soft[i] : 0f0
            v = (snvis + snvih) * (coef_col(coef, :v2t)[sp] / 2000f0)
        elseif _snag_east_vol(s.variant) && sn.height[i] > 0f0
            # same FMSVOL(I,HTIx)·DENIx, (SNVIS+SNVIH)·V2T with V2T pre-divided by 2000 (fmvinit.f:1094) — CS/LS/NE/SN
            sp = Int(sn.sp[i])
            vv = ffe_east_snag_vol_at(s, sp, sn.dbh[i], sn.height[i], sn.htcur[i])
            snvih = sn.den_hard[i] > 0f0 ? vv * sn.den_hard[i] : 0f0
            snvis = sn.den_soft[i] > 0f0 ? vv * sn.den_soft[i] : 0f0
            v = (snvis + snvih) * (coef_col(coef, :v2t)[sp] / 2000f0)
        else
            b = _snag_bole_tons(s, i)
            v = b * sn.den_soft[i] + b * sn.den_hard[i]      # (SNVIS+SNVIH)·V2T with the bole already in tons/stem
        end
        sn.dbh[i] <= 3f0 ? (tsng1 += v) : (tsng2 += v)
    end
    c2 = fs.cwd2b; c22 = fs.cwd2b2; tfm = size(c2, 3)
    @inbounds for isz in 0:3, idc in 1:4, itm in 1:tfm
        tsng1 += _FM_P2T * (c2[idc, isz + 1, itm] + c22[idc, isz + 1, itm])
        (isz > 0 && isz < 3) && (tsng2 += _FM_P2T * (c2[idc, isz + 4, itm] + c22[idc, isz + 4, itm]))
    end
    v2t = coef_col(coef, :v2t); ocw = t.ffe_oldcrw
    snfam = variant_code(s.variant) in ("CS", "LS", "NE", "SN")
    totfol = 0f0; tl1 = 0f0; tl2 = 0f0
    _fm_record_walk(t, vtrip) do i, pr
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        xv = _ffe_crownw(s, i, sp, d, h, Int(round(t.crown_pct[i])))
        totfol += xv[1] * pr * _FM_P2T
        for j in 1:3
            tl1 += (xv[j + 1] + ocw[j, i]) * _FM_P2T * pr
            j < 3 && (tl2 += _FM_P2T * pr * (xv[j + 4] + ocw[j + 3, i]))
        end
        # FMSVL2('L', LMERCH=.FALSE., no top-kill): MAX(X,MCF) for CS/LS/NE/SN, MAX(X,TCF) for the western variants
        vt = snfam ? max(0.005454154f0 * h, _ffe_stem_mcf(s, i, sp, d, h)) :
             _ffe_west_vol(s.variant) ? max(0.005454154f0 * h, ffe_west_nocut(s, sp, d, h)[1]) :
                     max(0.005454154f0 * h, t.cuft_vol[i])
        v2tp = v2t[sp] / 2000f0
        d <= 3f0 ? (tl1 += pr * vt * v2tp) : (tl2 += pr * vt * v2tp)
    end
    return (totfol = totfol, totliv1 = tl1, totliv2 = tl2, totsng1 = tsng1, totsng2 = tsng2,
            small2 = small2, large2 = large2, totlit = totlit, totduf = totduf,
            biolive = totfol + tl1 + tl2, biosnag = tsng1 + tsng2, bioddw = small2 + large2,
            bioflr = totlit + totduf, bioshrb = fs.flive[1] + fs.flive[2])
end

"""
    stand_carbon_report(s; vtrip=false) -> (; aboveground, merch, belowground, belowground_dead, standing_dead,
                                             down_wood, forest_floor, shrub_herb, total)

The Stand Carbon Report pools V(1..9) exactly as fmcrbout.f:93-183 forms them: per record (the FMMAIN list, see
`_fm_record_walk`) the Jenkins FMCBIO ABIO/MBIO/RBIO × FMPROB (V(3) always; V(1)/V(2) for CARBCALC method 1), or for the
FFE method V(2) += FMPROB·VT·V2T with FMSVL2 LMERCH=LVWEST (western MCF, eastern MAX(X,MCF)) and V(1)=BIOLIVE;
V(4..8)=BIOROOT/BIOSNAG/BIODDW/BIOFLR/BIOSHRB (`fmdout_bio`); ×0.5 (×0.37 for the floor); ICMETRC 1 ⇒ ×TItoTM/ACRtoHA,
2 ⇒ ×TItoTM; V(9)=V(1)+V(3)+V(5)+V(6)+V(7)+V(8) (+V(4) when LDCAY).
"""
function stand_carbon_report(s::StandState; vtrip::Bool = false)
    t = s.trees; coef = s.coef; v2t = coef_col(coef, :v2t)
    ffe = s.control.carbon_method == 0
    west = _ffe_west_vol(s.variant)
    v1 = 0f0; v2 = 0f0; v3 = 0f0
    _fm_record_walk(t, vtrip) do i, pr
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        abio, mbio, rbio = jenkins_biomass(coef, sp, d; dbhmin = _jenkins_dbhmin(s, sp))
        abio = abio * pr; mbio = mbio * pr; rbio = rbio * pr
        v3 += rbio
        if !ffe
            v1 += abio; v2 += mbio
        else
            vt = west ? ffe_west_nocut(s, sp, d, h)[2] : max(0.005454154f0 * h, _ffe_stem_mcf(s, i, sp, d, h))
            v2 += pr * vt * (v2t[sp] / 2000f0)
        end
    end
    fb = (s.fire !== nothing && s.fire.active) ? fmdout_bio(s; vtrip = vtrip) : nothing
    ffe && (v1 = fb === nothing ? 0f0 : fb.biolive)
    v4 = s.fire === nothing ? 0f0 : s.fire.bioroot
    v5 = fb === nothing ? 0f0 : fb.biosnag; v6 = fb === nothing ? 0f0 : fb.bioddw
    v7 = fb === nothing ? 0f0 : fb.bioflr;  v8 = fb === nothing ? 0f0 : fb.bioshrb
    v1 *= 0.50f0; v2 *= 0.50f0; v3 *= 0.50f0; v4 *= 0.50f0; v5 *= 0.50f0; v6 *= 0.50f0; v8 *= 0.50f0
    v7 *= 0.37f0
    if s.control.carbon_units == 1                         # V(I)·TItoTM/ACRtoHA (METRIC.F77 0.90718, 0.4046945)
        cv(x) = x * _TITOTM / _ACRTOHA
        v1 = cv(v1); v2 = cv(v2); v3 = cv(v3); v4 = cv(v4); v5 = cv(v5); v6 = cv(v6); v7 = cv(v7); v8 = cv(v8)
    elseif s.control.carbon_units == 2
        v1 *= _TITOTM; v2 *= _TITOTM; v3 *= _TITOTM; v4 *= _TITOTM; v5 *= _TITOTM; v6 *= _TITOTM; v7 *= _TITOTM; v8 *= _TITOTM
    end
    total = v1 + v3 + v5 + v6 + v7 + v8
    _FM_CRDCAY > 0f0 && (total = total + v4)               # LDCAY (dead-root decay on)
    return (; aboveground = v1, merch = v2, belowground = v3, belowground_dead = v4,
            standing_dead = v5, down_wood = v6, forest_floor = v7, shrub_herb = v8, total = total)
end

"""
    carbon_report_fmmain_v3(rep, s, stash, vtrip) -> NamedTuple

fmcrbout.f:98-112 V(3) = Σ RBIO(DBH(I))·FMPROB(I) is taken inside FMMAIN (gradd.f:118), after GRINCR's REGENT has
reset the DBH of the small trees it sizes directly — e.g. ie/regent.f:882 `DBH(K)=0.1+DIAM(ISPC)*.01+HK*0.001` for a
seedling still under 4.5 ft, per tripled record K. jl samples the non-fire report before growth, so re-derive V(3)
(and the stand total) from the post-REGENT DBH once `small_tree_growth!` has run: the central record's `t.dbh`, a
tripled copy's `stash.dbhU/dbhL` (−1 ⇒ unchanged). Every other pool is REGENT-independent and kept.
(MEASURED FVSie_g16 11855985010690 2006: seedlings DBH 0.1060/0.1054 vs 0.1 ⇒ Belowground_Live 7.19485 vs 7.19215.)
"""
function carbon_report_fmmain_v3(rep, s::StandState, stash, vtrip::Bool)
    t = s.trees; coef = s.coef; n = t.n
    trip = vtrip && stash !== nothing && hasproperty(stash, :dbhU) && length(stash.dbhU) >= n
    v3 = 0f0
    rb(i, d, pr) = (sp = Int(t.species[i]); (_, _, r) = jenkins_biomass(coef, sp, d; dbhmin = _jenkins_dbhmin(s, sp)); r * pr)
    @inbounds if vtrip
        for i in 1:n; v3 += rb(i, t.dbh[i], t.tpa[i] * 0.60f0); end
        for i in 1:n
            du = trip && stash.dbhU[i] >= 0f0 ? stash.dbhU[i] : t.dbh[i]
            dl = trip && stash.dbhL[i] >= 0f0 ? stash.dbhL[i] : t.dbh[i]
            v3 += rb(i, du, t.tpa[i] * 0.25f0); v3 += rb(i, dl, t.tpa[i] * 0.15f0)
        end
    else
        for i in 1:n; v3 += rb(i, t.dbh[i], t.tpa[i]); end
    end
    v3 *= 0.50f0
    if s.control.carbon_units == 1
        v3 = v3 * _TITOTM / _ACRTOHA
    elseif s.control.carbon_units == 2
        v3 *= _TITOTM
    end
    total = rep.aboveground + v3 + rep.standing_dead + rep.down_wood + rep.forest_floor + rep.shrub_herb
    _FM_CRDCAY > 0f0 && (total = total + rep.belowground_dead)
    return merge(rep, (belowground = v3, total = total))
end

# METRIC.F77 parameters the carbon reports convert with (TItoTM, ACRtoHA).
const _TITOTM  = 0.90718f0
const _ACRTOHA = 0.4046945f0

# The fixed header block of the Stand Carbon Report exactly as the Fortran prints it to the `.out`
# (fmcrbout.f FORMATs 700-709, with the FVS `1X,I5,1X` line-prefix stripped as it is in the file).
const _CARBON_SEP = "-"^110
const _CARBON_HEADER = (
    "                              ******  CARBON REPORT VERSION 1.0 ******",
    "                                         STAND CARBON REPORT (BASED ON STOCKABLE AREA)",
    "                         ALL VARIABLES ARE REPORTED IN METRIC TONS/HECTARE",
    "",
)
const _CARBON_COLHDR = (
    "      Aboveground Live    Belowground                        Forest             Total    Total     Carbon",
    "     ----------------- -----------------    Stand  -------------------------    Stand  Removed   Released",
    "YEAR    Total    Merch     Live     Dead     Dead      DDW    Floor  Shb/Hrb   Carbon   Carbon  from Fire",
)

"""
    carbon_report_row(s, year; removed=0f0, released=0f0) -> String

One data row of the Stand Carbon Report, byte-for-byte as the Fortran `.out` (fmcrbout.f FORMAT 800:
`I4` year, then ten `2X,F7.1` pool columns and a final `4X,F7.1`). Columns are the metric-tons/ha pools
from `stand_carbon_report`: Aboveground Total / Merch, Belowground Live / Dead, Stand Dead, DDW, Forest
Floor, Shrub-Herb, Total Stand Carbon, then Total Removed and Carbon Released-from-Fire (`removed` /
`released`, 0 without a harvest/fire this cycle).
"""
function _format_carbon_row(year::Integer, r; removed::Real = 0f0, released::Real = 0f0)
    vals = (r.aboveground, r.merch, r.belowground, r.belowground_dead, r.standing_dead,
            r.down_wood, r.forest_floor, r.shrub_herb, r.total, Float32(removed))
    io = IOBuffer()
    @printf(io, "%4d", year)
    for v in vals; @printf(io, "  %7.1f", v); end       # 10 × (2X, F7.1)
    @printf(io, "    %7.1f", Float32(released))          # final 4X, F7.1
    return String(take!(io))
end

carbon_report_row(s::StandState, year::Integer; removed::Real = 0f0, released::Real = 0f0) =
    _format_carbon_row(year, stand_carbon_report(s); removed = removed, released = released)

"""
    write_carbon_report(io, stand, ncyc; period=5, stand_id="", mgmt_id="NONE") -> IO

Write the FFE Stand Carbon Report to `io` exactly as FVS prints it to the `.out` (CARBREPT; fmcrbout.f),
for the inventory cycle plus `ncyc` grown cycles. Drives the per-cycle FFE fuel update + growth (the
same loop as the multi-cycle carbon test). The header block and the per-row format are byte-for-byte
vs the Fortran; the pool values are the validated metric-tons/ha pools (8/9 columns bit-exact, the
post-mortality DDW tracking within the LP growth tail — see FFE_FUEL_DYNAMICS_chunk_plan.md).
"""
function write_carbon_report(io::IO, stand::StandState, ncyc::Integer;
                             period::Integer = 5, stand_id::AbstractString = "",
                             mgmt_id::AbstractString = "NONE")
    println(io, _CARBON_SEP)
    for h in _CARBON_HEADER; println(io, h); end
    println(io, "STAND ID: ", rpad(strip(stand_id), 26), "    MGMT ID: ", strip(mgmt_id))
    println(io, _CARBON_SEP)
    for h in _CARBON_COLHDR; println(io, h); end
    println(io, _CARBON_SEP)
    fs = stand.fire
    # seed the inventory snags from the input dead-tree records (no-op when there are none); the
    # per-cycle snag falldown then runs inside grow_cycle! (update_snags!, simulate.jl:211).
    fs !== nothing && fs.active && (ffe_seed_input_snags!(stand); ffe_add_snaginit!(stand))
    # Snapshot the INVENTORY crown as the OLD state (FVS calls FMOLDC in the inventory FMMAIN, before the
    # first grow), so the FIRST cycle's crown-lift has a valid OLDCRW. Without this the 1st cycle's
    # crown-lift is skipped (ffe_oldht=0), losing ~1.9 t/ac of fine down-wood — the DDW sizes-1-3 gap.
    fs !== nothing && fs.active && snapshot_ffe_oldcrown!(stand)
    for c in 0:ncyc
        compute_density!(stand)
        # Refresh cover type + live herb/shrub fuels (FLIVE) from the CURRENT (post-growth) stand at
        # each report point — FVS reports the cycle's own live fuels, so computing them only pre-growth
        # lags the Shrub/Herb column one cycle. fmcba! loads the initial dead fuels just once (fuels_init).
        if fs !== nothing && fs.active
            compute_forest_type!(stand); fmcba!(stand)
        end
        println(io, carbon_report_row(stand, Int(current_cycle_year(stand))))
        if c < ncyc
            fs !== nothing && fs.active && ffe_fuel_update!(stand, period)
            grow_cycle!(stand; fint = Float32(period))
            # crown-lift: the lower crown shed as the base rises during THIS cycle's growth (FMSDIT);
            # the per-year input is applied in the NEXT cycle's fuel loop (FMCADD), then snapshot the
            # post-growth crown as next cycle's OLD state (FMOLDC). Both no-op on the first grow
            # (ffe_old* unset ⇒ zero), matching FVS's ICYC>1 gate.
            if fs !== nothing && fs.active
                compute_crown_lift!(stand, period)
                snapshot_ffe_oldcrown!(stand)
            end
        end
    end
    return io
end

"""
    write_carbon_report_block(io, rows; stand_id="", mgmt_id="NONE") -> IO

Write the Stand Carbon Report header block + the per-cycle `rows` to `io`, byte-for-byte as the Fortran
`.out`. Each element of `rows` is a `(year, report)` tuple (`report` = a `stand_carbon_report` named
tuple) collected during the main simulation loop (`write_sum_file`'s `carbon_collect`).
"""
function write_carbon_report_block(io::IO, rows::AbstractVector;
                                   stand_id::AbstractString = "", mgmt_id::AbstractString = "NONE")
    println(io, _CARBON_SEP)
    for h in _CARBON_HEADER; println(io, h); end
    println(io, "STAND ID: ", rpad(strip(stand_id), 26), "    MGMT ID: ", strip(mgmt_id))
    println(io, _CARBON_SEP)
    for h in _CARBON_COLHDR; println(io, h); end
    println(io, _CARBON_SEP)
    for row in rows
        println(io, _format_carbon_row(row[1], row[2]; released = length(row) >= 6 ? row[6] : 0f0))
    end
    return io
end
