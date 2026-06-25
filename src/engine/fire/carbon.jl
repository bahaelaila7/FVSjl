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
function stand_live_carbon(s::StandState)
    t = s.trees; coef = s.coef
    above = 0f0; merch = 0f0; root = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        a, m, r = jenkins_biomass(coef, t.species[i], t.dbh[i])
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
    vol(sz, K) = (let v = 0f0; for L in 1:4; v += cw[sz, K, L]; end; v end) * 2000f0 / den[K]
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
function ffe_fuel_loadings(s::StandState)
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
    # standing snags: bole biomass by DBH + the CWD2B crown (sizes 0-3 → ≤3, 4-5 → >3), tons/ac
    coef = s.coef; sn = fs.snags; snag_lt3 = 0f0; snag_ge3 = 0f0
    @inbounds for i in eachindex(sn.sp)
        den = sn.den_hard[i] + sn.den_soft[i]; den > 0f0 || continue
        b = sn.bolevol[i]; b <= 0f0 && (b = let (a,_,_) = jenkins_biomass(coef, sn.sp[i], sn.dbh[i]); a end)
        (sn.dbh[i] <= 3f0 ? (snag_lt3 += b*den) : (snag_ge3 += b*den))
    end
    snag_lt3 += sum(@view fs.cwd2b[:, 1:4, :]) * _FM_P2T     # CWD2B crown sizes 0-3 (idx 1-4)
    snag_ge3 += sum(@view fs.cwd2b[:, 5:6, :]) * _FM_P2T     # CWD2B crown sizes 4-5 (idx 5-6)
    # standing live: foliage + woody crown + stem, split by tree DBH (fmdout.f:218-258)
    t = s.trees; v2t = coef_col(coef, :v2t); foliage = 0f0; live_lt3 = 0f0; live_ge3 = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); d = t.dbh[i]
        xv = crown_biomass(s, sp, d, t.height[i], Int(round(t.crown_pct[i])))
        foliage += xv[1] * t.tpa[i] * _FM_P2T
        woody = 0f0; for sz in 1:5; woody += xv[sz+1]; end
        stem = _fm_cuft(s, sp, d, t.height[i]) * v2t[sp]
        (d <= 3f0 ? (live_lt3 += (woody + stem) * t.tpa[i] * _FM_P2T) :
                    (live_ge3 += (woody + stem) * t.tpa[i] * _FM_P2T))
    end
    stand_total = snag_lt3 + snag_ge3 + foliage + live_lt3 + live_ge3
    return (; litter, duff, lt3, ge3, s3to6, s6to12, ge12, herb, shrub, surf_total,
            snag_lt3, snag_ge3, foliage, live_lt3, live_ge3, stand_total,
            total_biomass = surf_total + stand_total, consumed = 0f0, removed = 0f0)
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
        xv = crown_biomass(s, sp, d, h, Int(round(t.crown_pct[i])))   # (foliage, woody 1..5), lb
        crown = xv[1]; for sz in 1:5; crown += xv[sz + 1]; end         # foliage + all woody (lb)
        stem = _fm_cuft(s, sp, d, h) * v2t[sp]                         # stem biomass (lb; v2t is raw lb/ft³)
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

"""
    stand_carbon_report(s) -> (; aboveground, merch, belowground, standing_dead,
                                 down_wood, forest_floor, shrub_herb, total)

The Stand Carbon Report pools in **metric tons C / hectare** (CARBREPT, CARBCALC method 1 = Jenkins;
fmcrbout.f), matching the report columns. The live aboveground / merchantable / belowground (root)
pools come from `stand_live_carbon` (Jenkins biomass × 0.5 × TPA) — bit-exact vs the Fortran report.
The standing-dead / down-wood / forest-floor / shrub-herb pools come from the FFE surface-fuel model
(`fmcba!`) and are zero unless it has populated `fire.cwd`/`fire.flive`; down-wood and forest floor
(at the 0.5 / 0.37 carbon fractions) reconcile bit-exact, while shrub-herb tracks `FLIVE`, which
carries the FFE live-fuel-loading residual. `total` = above + below + snag + ddw + floor + shrub
(fmcrbout.f:178). NB the FFE pools require `fmcba!` to have run this cycle (the per-cycle FFE fuel
update — only triggered on a fire event in the current main path; that lifecycle wiring is the
remaining increment).
"""
function stand_carbon_report(s::StandState)
    lc = stand_live_carbon(s)
    # CARBCALC method: 1 = JENKINS national biomass (default), 0 = FFE crown+stem biomass. Belowground
    # (roots) is the Jenkins value in BOTH methods (fmcrbout.f:144-146); only Above/Merch differ.
    if s.control.carbon_method == 0
        fc = ffe_live_carbon(s)
        above = fc.aboveground * _TONAC_TO_MTHA
        merch = fc.merch       * _TONAC_TO_MTHA
    else
        above = lc.aboveground * _TONAC_TO_MTHA
        merch = lc.merch       * _TONAC_TO_MTHA
    end
    below = lc.belowground * _TONAC_TO_MTHA
    sd  = standing_dead_carbon(s) * _TONAC_TO_MTHA
    bd  = belowground_dead_carbon(s) * _TONAC_TO_MTHA
    dw  = down_wood_carbon(s)     * _TONAC_TO_MTHA
    ff  = forest_floor_carbon(s)  * _TONAC_TO_MTHA
    sh  = shrub_herb_carbon(s)    * _TONAC_TO_MTHA
    return (; aboveground = above, merch = merch, belowground = below, belowground_dead = bd,
            standing_dead = sd, down_wood = dw, forest_floor = ff, shrub_herb = sh,
            total = above + below + sd + dw + ff + sh)
end

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
    fs !== nothing && fs.active && ffe_seed_input_snags!(stand)
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
    for row in rows; println(io, _format_carbon_row(row[1], row[2])); end
    return io
end
