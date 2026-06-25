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
    above = lc.aboveground * _TONAC_TO_MTHA
    merch = lc.merch       * _TONAC_TO_MTHA
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
function carbon_report_row(s::StandState, year::Integer; removed::Real = 0f0, released::Real = 0f0)
    r = stand_carbon_report(s)
    vals = (r.aboveground, r.merch, r.belowground, r.belowground_dead, r.standing_dead,
            r.down_wood, r.forest_floor, r.shrub_herb, r.total, Float32(removed))
    io = IOBuffer()
    @printf(io, "%4d", year)
    for v in vals; @printf(io, "  %7.1f", v); end       # 10 × (2X, F7.1)
    @printf(io, "    %7.1f", Float32(released))          # final 4X, F7.1
    return String(take!(io))
end

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
