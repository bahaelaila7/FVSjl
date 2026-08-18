# =============================================================================
# ontario/diameter_growth.jl — Penner large-tree diameter growth (canada/on/dgf.f).
#
#   on_penner_dds(...)  — the core annual-iteration Penner DG.  *** VALIDATED BIT-EXACT ***
#       replaying the g16 per-tree dump reproduces DBHM_final, DIAGR AND DDS 8/8 bit-exact
#       (Float32-hex) across 8 distinct equations incl. both AGS/UGS quality species
#       (scratchpad/on/replay_penner.jl vs /workspace/.onwork/FVSon_dgfdump).
#   dgf!(s, ::Ontario) — per-tree wrapper writing WK2 = DDS.  Stand-context inputs
#       (SIM/BAM/QMDM metric, from dgcons/sitset) are wired to StandState here; that outer
#       wiring is PENDING full-engine integration on kt-variant-port (the core math above is
#       what the go/no-go proves).  Penner shade-tolerant maple/beech HYBRID adjustment
#       (dgf.f after the species loop) is NOT YET ported — scoped as a follow-on slice.
#
# Transcendentals MUST be glibc single-precision logf/expf: the 10-iteration compounding makes
# native Julia log/exp and the core/fmath shim drift ~1 ULP; measured glibc = 0-mismatch.
# =============================================================================

@inline on_logf(x::Float32) = ccall(:logf, Float32, (Float32,), x)
@inline on_expf(x::Float32) = ccall(:expf, Float32, (Float32,), x)

"""
    on_penner_dds(ksp, ags, diam_in, sim, bam, qmdm, balm, htm, bark) -> (dbhm_final, diagr, dds)

Core Ontario Penner large-tree DG (canada/on/dgf.f main body). `ksp` is the 1-based Penner
equation (OSPMAP value, +1 for a quality tree); `ags` is 1 for an AGS/quality tree else 0.
All stand inputs are METRIC (cm, m, m²/ha); `diam_in` is DBH in inches. Bit-exact vs FVSon_g16.
"""
@inline function on_penner_dds(ksp::Int, ags::Int, diam_in::Float32, sim::Float32, bam::Float32,
                               qmdm::Float32, balm::Float32, htm::Float32, bark::Float32)
    dbhm = diam_in * ON_INtoCM
    x1 = -(ON_B0[ksp]) - (ON_BBAL[ksp]*balm) - (ON_BHT[ksp]*htm) - (ON_BSI[ksp]*sim) -
         (ON_BBA[ksp]*bam) - (ON_BDBHQ[ksp]*qmdm) - (ON_BAGS[ksp]*Float32(ags))
    @inbounds for _ in 1:10
        dgln = x1 + (ON_B1[ksp]*on_logf(dbhm)) - (ON_B2[ksp]*dbhm)
        deld = min(max(dgln, -5f0), 5f0)
        deld = on_expf(deld)
        deld = min(max(deld, 0.0001f0), ON_B95[ksp])
        dbhm = dbhm + deld
    end
    d = dbhm * ON_CMtoIN
    diagr = (d - diam_in) * bark
    dds = diagr * (2f0*diam_in*bark + diagr)
    return dbhm, diagr, dds
end

@inline on_powf(x::Float32, y::Float32) = ccall(:powf, Float32, (Float32, Float32), x, y)

"""
    on_bratio(is, d, h) -> Float32

Ontario bark ratio (canada/on/bratio.f). `is` = FVS species (1-based), `d` = diameter
(inches, the GROWN diameter D = DBHM·CMtoIN as passed at dgf.f:362), `h` = height (ft).
Two sequential override blocks over the `ON_BKRAT` fallback: species-constant bark, then
(for a non-zero diameter) the Zakrzewski metric H/D equations. Bit-exact vs FVSon_g16.
"""
@inline function on_bratio(is::Int, d::Float32, h::Float32)::Float32
    br = is == 0 ? 0.93f0 : ON_BKRAT[is]
    # --- constant-value overrides (bratio.f) ---
    if is == 3 || is == 4
        br = 0.92608f0                                   # red pine
    elseif is == 1 || is == 69 || is == 10 || is == 13
        br = 0.9609f0                                    # jack pine (nat/plant), tamarack, other softwoods
    elseif is == 6 || is == 71 || is == 7
        br = 0.955f0                                     # white / Norway spruce
    elseif is == 8
        br = 0.9456f0                                    # balsam fir
    elseif is == 24 || is == 43
        br = 0.9451f0                                    # birch
    elseif is == 15 || is == 16 || is == 17 || is == 40 || is == 41 || is == 42
        br = 0.9365f0                                    # aspen / ash / cottonwood
    end
    d == 0f0 && return br
    # --- Ontario (Zakrzewski) metric H/D equations, applied after the constants ---
    d_cm = d * ON_INtoCM
    h_m  = h * ON_FTtoM
    hdrat = h_m / d_cm
    if is == 5 || is == 70                               # white pine
        br = 1f0 - 0.04513f0 * on_powf(d_cm, 1.168567f0 - 1.0f0)
    elseif is == 9 || is == 72                           # black spruce (clamp HDRAT<14)
        br = hdrat < 14f0 ? 1f0 - 0.067259f0 * hdrat : 1f0 - 0.067259f0 * 14f0
    elseif is == 11 || is == 14                          # white cedar / N. red cedar (clamp HDRAT<7)
        br = hdrat < 7f0 ? 1f0 - 0.13665f0 * hdrat : 1f0 - 0.13665f0 * 7f0
    elseif is == 18 || is == 19 || is == 26 || is == 27 || is == 44 || is == 46 ||
           is == 48 || is == 49 || is == 50 || is == 51 || is == 52 || is == 54 ||
           is == 57 || is == 59 || is == 64 || is == 65 || is == 66 || is == 67  # sugar/red maple group (clamp HDRAT<8)
        br = hdrat < 8f0 ? 1f0 - 0.111681f0 * hdrat : 1f0 - 0.111681f0 * 8f0
    end
    return br
end

"""
    dgf!(s::StandState, ::Ontario)

Ontario large-tree DG hook (canada/on/dgf.f). Faithful two-pass transcription:

 1. species-sorted loop (ISCT/IND1, ascending species — the order set by `species_sort!`
    which the driver runs first): resolve the Penner equation via OSPMAP (+1 / AGS=1 for the
    8 LQUAL species when the tree quality code IMC==1), iterate the annual Penner increment to
    the grown metric diameter, compute per-tree BARK = on_bratio(sp, D, HT), and store
    DIAGR = (D−DBH)·BARK in `temd`.  The stand context is assembled exactly as dgf.f does:
    BAM = BA·FT2pACRtoM2pHA, QMDM = RMSQD·INtoCM, SIM = SITEAR(sp)·FTtoM,
    BALM = (1−PCT/100)·BAM (PCT lives in `crown_ratio`, FVS convention), HTM = HT·FTtoM.

 2. final loop (dgf.f DO 303): DDS = DIAGR·(2·DBH·BARK + DIAGR) then WK2 = ln(DDS)+COR.
    FAITHFUL QUIRK: the FVS DDS term reuses the SCALAR `BARK` left over from the LAST tree of
    the last species processed in pass 1 (dgf.f never re-indexes BARK in DO 303), NOT each
    tree's own bark. `barklast` reproduces this bit-for-bit (measured: fort.772 dumps BARK =
    beech 0.94 for every tree in ont01). The LHYBRID maple/beech stand adjustment is hard-
    disabled in dgf.f (LHYBRID=.FALSE. at dgf.f:402, "Temporarily disabled for BlueSource"),
    so it is intentionally not ported.

`species_sort!(s)` MUST have run so `sp_count_tab`/`idx1` (ISCT/IND1) are current; the shared
`diameter_growth!(::AbstractVariant)` driver does this before calling `dgf!`.
"""
function dgf!(s::StandState, ::Ontario)
    p, t, c = s.plot, s.trees, s.calib
    ctl = s.control
    wk2  = view(s.scratch.wk, 2, :)
    isct = ctl.sp_count_tab
    ind1 = s.scratch.idx1
    bam  = p.basal_area * ON_FT2pACRtoM2pHA   # BA (ft²/ac) -> m²/ha
    qmdm = p.qmd * ON_INtoCM                  # RMSQD (in) -> cm

    # TEMD holds DIAGR per tree (dgf.f DO 4 seeds it with DIAM; overwritten for every grown tree).
    temd = Vector{Float32}(undef, t.n)
    @inbounds for i in 1:t.n
        temd[i] = t.dbh[i]
        wk2[i]  = 0f0
    end

    # --- Pass 1 (dgf.f DO 20 / DO 10): per-tree DIAGR, tracking the leftover BARK scalar ---
    barklast = 0.93f0
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]
        i1 == 0 && continue
        (sp < 1 || sp > length(ON_OSPMAP)) && continue
        base = ON_OSPMAP[sp]
        base == 0 && continue                 # dgf.f:299 OSPMAP==0 -> GOTO 20 (no ON large-tree DG)
        i2 = isct[sp, 2]
        sim = p.sp_site_index[sp] * ON_FTtoM  # SIM = SITEAR(sp)·FTtoM
        islqual = sp in ON_LQUAL_SPP
        for i3 in i1:i2
            i = ind1[i3]
            d = t.dbh[i]
            d <= 0f0 && continue              # temd stays = DIAM (DO 4 seed); wk2 stays 0
            if islqual && t.mort_code[i] == 1  # IMC==1 -> AGS equation (+1)
                ksp = base + 1; ags = 1
            else
                ksp = base; ags = 0
            end
            pct  = t.crown_ratio[i]           # FVS PCT (percentile of BA in larger trees)
            balm = (1f0 - pct / 100f0) * bam
            htm  = t.height[i] * ON_FTtoM
            dbhm_final, _, _ = on_penner_dds(ksp, ags, d, sim, bam, qmdm, balm, htm, 1f0)
            dgrown = dbhm_final * ON_CMtoIN
            bark = on_bratio(sp, dgrown, t.height[i])   # BRATIO uses the GROWN diameter (dgf.f:361-362)
            temd[i] = (dgrown - d) * bark
            barklast = bark
        end
    end

    # --- Pass 2 (dgf.f DO 303): DDS with the leftover BARK, then WK2 = ln(DDS)+COR ---
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]
        i1 == 0 && continue
        i2 = isct[sp, 2]
        cor = c.dg_cor[sp]
        cor2on = ctl.dg_cor2_on
        cor2 = cor2on ? ctl.dg_cor2[sp] : 0f0
        for i3 in i1:i2
            i = ind1[i3]
            diagr = temd[i]
            (cor2on && cor2 > 0f0) && (diagr *= cor2)
            diagr <= 0.0001f0 && (diagr = 0.0001f0)
            d = t.dbh[i]
            dds = diagr * (2f0 * d * barklast + diagr)
            dds > 0f0 && (wk2[i] = on_logf(dds) + cor)
        end
    end
    return s
end
