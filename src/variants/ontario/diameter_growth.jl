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

"""
    dgf!(s::StandState, ::Ontario)

Ontario large-tree DG hook: per tree, resolve the Penner equation via OSPMAP (+1 AGS for the
8 LQUAL species when IMC==1) and write WK2 = DDS. Mirrors dgf.f's per-species tree loop.
NOTE: SIM/BAM/QMDM sourcing + the maple/beech HYBRID adjustment are pending engine integration.
"""
function dgf!(s::StandState, ::Ontario)
    p, t = s.plot, s.trees
    wk2 = view(s.scratch.wk, 2, :)
    bam  = p.basal_area            # metric BA (m²/ha) — confirm units vs dgf.f DGCONS on wire-up
    qmdm = p.qmd * ON_INtoCM       # QMD (in) -> cm ; PLACEHOLDER: match dgf.f DGCONS QMDM source
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = Int(t.species[i])
        (sp < 1 || sp > length(ON_OSPMAP)) && continue
        base = ON_OSPMAP[sp]
        base == 0 && continue
        imc = t.imc[i]             # tree quality code (1 = AGS); trees.jl field TBD on wire-up
        if (sp in ON_LQUAL_SPP) && imc == 1
            ksp = base + 1; ags = 1
        else
            ksp = base; ags = 0
        end
        pct  = t.crown_ratio[i]                       # FVS PCT
        balm = (1f0 - pct/100f0) * bam
        htm  = t.ht[i] * ON_FTtoM                     # ON_FTtoM from metric consts
        sim  = s.calib.site_index_sp[sp]              # per-species SI (sitset); field TBD on wire-up
        bark = bark_ratio(s, Ontario())               # bratio.f metric bark; TBD
        _, _, dds = on_penner_dds(ksp, ags, d, sim, bam, qmdm, balm, htm, bark)
        wk2[i] = dds > 0f0 ? on_logf(dds) : 0f0       # WK2 = ln(DDS)+COR ; +COR pending (dgcons)
    end
    return s
end

diameter_growth!(s::StandState, v::Ontario) = dgf!(s, v)
