# Dwarf-mistletoe infection & mortality summary report (mistoe/misprt.f MISPRT).
#
# Produces the two per-cycle DM summary aggregations that FVS prints to the .out (auto, whenever
# mistletoe damage codes are present) AND writes to the FVS_DM_Spp_Sum / FVS_DM_Stnd_Sum DBS tables
# (via DBSMIS1/DBSMIS2, gated by the MISRPTS database keyword). Also produces the by-DBH-class
# aggregation for FVS_DM_Sz_Sum (DBSMIS3, gated additionally by the MISTPRT report keyword).
#
# The aggregation is fully deterministic in the (validated) DM model quantities: per-tree TPA (PROB),
# DMR (IMIST, seeded from damage codes 30-34), DBH, total cubic volume (CFV), and the per-tree DM
# mortality rate (mismrt.f — ie_dm_mortality_rate / cr_dm_mortality_rate). Validated bit-exact vs the
# relinked FVSie_clean .out mistletoe tables on scratchpad/dm/dm6 (cyc0: species DF DMR 2.4/DMI 3.1/
# INF 1110/MORT 56/%INF 75/%MORT 4/COMP 85; stand T1729/BA660/VOL14666/INF 1110·353·7854/MORT 56·17·382/
# %64·54·3·3/DMR 2.0/DMI 3.1 — every column matches).

"minimum tree DBH for DMR/DMI statistics (misin0.f:79 DMRMIN default; MISTMIN keyword not yet wired)"
const _DM_DMRMIN = 1.0f0

# Variants that run the base mistoe.f DM model + its MISPRT report: the N-Rockies cluster
# (_ie_mis_variant), CentralRockies (its own cr_mistoe!/cr_dm_mortality — always-on, per-tree DMR),
# and BritishColumbia (NEWSPRED, keyword-activated). Distinct from `_dm_effects_variant` (the
# effect-application gate, which CR bypasses via its own cr_ gates).
@inline _dm_report_variant(v)::Bool = _ie_mis_variant(v) || v isa CentralRockies || v isa BritishColumbia

"true when any tree carries a dwarf-mistletoe rating (misprt.f DMFLAG) on a DM-report variant"
function _dm_report_active(s::StandState)::Bool
    _dm_report_variant(s.variant) || return false
    t = s.trees
    @inbounds for i in 1:t.n
        t.dmr[i] > 0 && return true
    end
    return false
end

"per-tree DM mortality proportion for the cycle (mismrt.f), dispatched to the active DM model"
@inline function _dm_mortality_rate(s::StandState, sp::Int, dmr::Int, dbh::Float32, fint::Float32)::Float32
    dmr == 0 && return 0f0
    if s.variant isa CentralRockies
        return cr_dm_mortality_rate(sp, dmr, dbh, fint)
    end
    _, _, pmc, maxsp = _mis_tables(s.variant)
    return ie_dm_mortality_rate(pmc, maxsp, sp, dmr, dbh, fint)
end

"""
    mistletoe_report(s; fint, top4, nage) -> NamedTuple

Aggregate the dwarf-mistletoe summary (misprt.f MISPRT) at the current (start-of-cycle) stand state.
`fint` is the cycle length used for the per-tree DM mortality projection (MISMRT). `top4` is a
persistent cache of the top-4 most-infected species indices (FVS LSORT4: sorted ONCE, on the first
cycle with infection, then frozen) — pass the same `Vector{Int}` across cycles; it is filled in place.
`nage` is the current stand age (misprt.f NAGE = IAGE + year − inventory year).

Returns `(; nage, stand, species, dbhclass)` where `stand` carries the per-acre composite columns,
`species` is a length-≤4 vector of per-species top-infected rows, and `dbhclass` the 10 2-inch DBH
classes for the size table.
"""
function mistletoe_report(s::StandState; fint::Float32, top4::Vector{Int}, nage::Integer)
    t = s.trees
    n = t.n
    maxsp = nspecies(s.variant)
    SPTPAI = zeros(Float64, maxsp); SPTPAX = zeros(Float64, maxsp)
    SPTPAT = zeros(Float64, maxsp); SPDMRS = zeros(Float64, maxsp)
    SPTPAM = zeros(Float64, maxsp)
    STVOL = 0.0; STVOLI = 0.0; STVOLM = 0.0; STTPAM = 0.0
    # 20 2-inch DBH classes (0-2.9, 3-4.9, ...); the report lumps 11-20 into class 10.
    DCTPA = zeros(Float64, 20); DCTPAX = zeros(Float64, 20); DCINF = zeros(Float64, 20)
    DCMRT = zeros(Float64, 20); DCSUM = zeros(Float64, 20)
    dmrmin = _DM_DMRMIN
    @inbounds for i in 1:n
        sp = Int(t.species[i]); (sp < 1 || sp > maxsp) && continue
        p = Float64(t.tpa[i]); idmr = Int(t.dmr[i]); dbh = t.dbh[i]
        cfv = Float64(t.cuft_vol[i])
        # species infection totals (DBH >= DMRMIN gate on the DMR statistics)
        if dbh >= dmrmin
            idmr > 0 && (SPTPAI[sp] += p)
            SPDMRS[sp] += idmr * p
            SPTPAX[sp] += p
        end
        SPTPAT[sp] += p
        # stand volumes
        STVOL += cfv * p
        (idmr > 0 && dbh >= dmrmin) && (STVOLI += cfv * p)
        # per-tree DM mortality
        dmm = p * Float64(_dm_mortality_rate(s, sp, idmr, dbh, fint))
        SPTPAM[sp] += dmm
        if dmm > 0.0
            STTPAM += dmm
            STVOLM += cfv * dmm
        end
        # DBH-class buckets (0-2.9,3-4.9,...): IDBH = IVAL/2 + MOD(IVAL,2), clamped 1..20
        ival = trunc(Int, dbh); idbh = clamp(ival ÷ 2 + ival % 2, 1, 20)
        DCTPA[idbh] += p
        dbh >= dmrmin && (DCTPAX[idbh] += p)
        DCMRT[idbh] += dmm
        if idmr > 0 && dbh >= dmrmin
            DCINF[idbh] += p
            DCSUM[idbh] += idmr * p
        end
    end
    STTPAI = sum(SPTPAI); STTPAT = sum(SPTPAT); STTPAX = sum(SPTPAX); STDMRS = sum(SPDMRS)
    ba = Float64(s.plot.basal_area)
    STBAI = STVOL != 0 ? ba / STVOL * STVOLI : 0.0
    STBAM = STVOL != 0 ? ba / STVOL * STVOLM : 0.0
    STDMR = STTPAX != 0 ? STDMRS / STTPAX : 0.0
    STDMI = STTPAI != 0 ? STDMRS / STTPAI : 0.0
    STPIT = STTPAT != 0 ? STTPAI / STTPAT * 100.0 : 0.0
    STPMT = STTPAT != 0 ? STTPAM / STTPAT * 100.0 : 0.0
    STPIV = STVOL != 0 ? STVOLI / STVOL * 100.0 : 0.0
    STPMV = STVOL != 0 ? STVOLM / STVOL * 100.0 : 0.0

    # Top-4 most-infected species (misprt.f LSORT4: sort ONCE by SPTPAI, descending, then freeze).
    if isempty(top4)
        order = sortperm(SPTPAI; rev = true)
        for k in 1:min(4, maxsp)
            SPTPAI[order[k]] > 0.0 ? push!(top4, order[k]) : push!(top4, 0)
        end
    end
    species = NamedTuple[]
    for infno in top4
        infno == 0 && continue
        SPTPAI[infno] <= 0.0 && continue
        spdmr = SPTPAX[infno] != 0 ? SPDMRS[infno] / SPTPAX[infno] : 0.0
        spdmi = SPTPAI[infno] != 0 ? SPDMRS[infno] / SPTPAI[infno] : 0.0
        sppin = SPTPAT[infno] != 0 ? SPTPAI[infno] / SPTPAT[infno] * 100.0 : 0.0
        sppmr = SPTPAT[infno] != 0 ? SPTPAM[infno] / SPTPAT[infno] * 100.0 : 0.0
        sppoc = STTPAT != 0 ? SPTPAT[infno] / STTPAT * 100.0 : 0.0
        push!(species, (sp = infno, mean_dmr = spdmr, mean_dmi = spdmi,
                        inf_tpa = SPTPAI[infno], mort_tpa = SPTPAM[infno],
                        inf_pct = sppin, mort_pct = sppmr, comp_pct = sppoc))
    end

    # DBH-class table (10 printed classes; 11-20 lumped into 10), all-trees + infected-only DMRs.
    @inbounds for c in 11:20
        DCTPA[10] += DCTPA[c]; DCTPAX[10] += DCTPAX[c]; DCINF[10] += DCINF[c]
        DCMRT[10] += DCMRT[c]; DCSUM[10] += DCSUM[c]
    end
    dcdmr = zeros(Float64, 10); dcdmi = zeros(Float64, 10)
    @inbounds for c in 1:10
        DCTPAX[c] != 0 && (dcdmr[c] = DCSUM[c] / DCTPAX[c])
        DCINF[c]  != 0 && (dcdmi[c] = DCSUM[c] / DCINF[c])
    end

    stand = (sttpat = STTPAT, ba = ba, stvol = STVOL, sttpai = STTPAI, stbai = STBAI,
             stvoli = STVOLI, sttpam = STTPAM, stbam = STBAM, stvolm = STVOLM,
             stpit = STPIT, stpiv = STPIV, stpmt = STPMT, stpmv = STPMV,
             stdmr = STDMR, stdmi = STDMI)
    dbhclass = (dctpa = DCTPA[1:10], dcinf = DCINF[1:10], dcmrt = DCMRT[1:10],
                dcdmr = dcdmr, dcdmi = dcdmi)
    return (nage = Int(nage), stand = stand, species = species, dbhclass = dbhclass)
end
