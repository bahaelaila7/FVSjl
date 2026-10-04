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

"""
true when MISPRT's DMFLAG is set on a DM-report variant: any tree carries a dwarf-mistletoe rating (misprt.f:375), or —
after a projection cycle (DMFLAG is reset only at ICYC=0, misprt.f:265) — the cycle's MISTOE/MISINF left it set
(mistoe.f:267, misinf.f:182: a MISTPINF card on a host species with no trees still prints a zero-infection row;
MEASURED FVSem_g16 231908428020004 MISTPINF: one 2023 FVS_DM_Stnd_Sum row, Inf_TPA 0, and an empty FVS_DM_Spp_Sum).
"""
function _dm_report_active(s::StandState)::Bool
    _dm_report_variant(s.variant) || return false
    (_ie_mis_variant(s.variant) && s.control.cycle > 0 && s.control.dm_flag) && return true
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
    # misprt.f MISPRT, REAL*4 throughout (P, IDMR·P, the SP*/ST*/DC* accumulators and the ratios are REAL), and in its
    # own loop orders: the per-species sums walk each species' records in IND1 order (ISCT/IND1, the first DO 90 loop),
    # the stand totals add the species subtotals in species order, SPTPAM walks IND1 again (DO 265), and the DBH-class /
    # volume sums walk the records 1..ITRN. jl summed in Float64 in record order (FVS_DM_* Mean_DMR/Mean_DMI a few ULP
    # to ~1e-7 relative off).
    t = s.trees
    maxsp = nspecies(s.variant)
    SPTPAI = zeros(Float32, maxsp); SPTPAX = zeros(Float32, maxsp)
    SPTPAT = zeros(Float32, maxsp); SPDMRS = zeros(Float32, maxsp)
    SPTPAM = zeros(Float32, maxsp)
    STVOL = 0f0; STVOLI = 0f0; STVOLM = 0f0; STTPAM = 0f0
    DCTPA = zeros(Float32, 20); DCTPAX = zeros(Float32, 20); DCINF = zeros(Float32, 20)
    DCMRT = zeros(Float32, 20); DCSUM = zeros(Float32, 20)
    dmrmin = _DM_DMRMIN
    ind1 = _ind1_order(s)
    # DMMTPA(I) (MISMRT(.FALSE.), mismrt.f:141-195): WKI=PTPA·DMMORT, 0 when PROB≤0 or IDMR=0
    dmm = zeros(Float32, t.n)
    @inbounds for i in 1:t.n
        sp = Int(t.species[i]); (sp < 1 || sp > maxsp) && continue
        p = t.tpa[i]; idmr = Int(t.dmr[i])
        (p <= 0f0 || idmr == 0) && continue
        dmm[i] = p * _dm_mortality_rate(s, sp, idmr, t.dbh[i], fint)
    end
    # DO 90: per-species sums in IND1 order; the stand totals add each species' subtotal (species order)
    STTPAI = 0f0; STTPAT = 0f0; STTPAX = 0f0; STDMRS = 0f0
    @inbounds for i in ind1
        sp = Int(t.species[i]); (sp < 1 || sp > maxsp) && continue
        p = t.tpa[i]; idmr = Int(t.dmr[i])
        if t.dbh[i] >= dmrmin
            idmr > 0 && (SPTPAI[sp] += p)
            SPDMRS[sp] += Float32(idmr) * p
            SPTPAX[sp] += p
        end
        SPTPAT[sp] += p
    end
    @inbounds for sp in 1:maxsp
        STTPAI += SPTPAI[sp]; STTPAT += SPTPAT[sp]; STTPAX += SPTPAX[sp]; STDMRS += SPDMRS[sp]
    end
    # tree loop 1..ITRN: DBH classes and the stand volume / mortality sums
    @inbounds for i in 1:t.n
        sp = Int(t.species[i]); (sp < 1 || sp > maxsp) && continue
        p = t.tpa[i]; idmr = Int(t.dmr[i]); dbh = t.dbh[i]; cfv = t.cuft_vol[i]
        ival = trunc(Int, dbh); idbh = clamp(ival ÷ 2 + ival % 2, 1, 20)
        DCTPA[idbh] += p
        dbh >= dmrmin && (DCTPAX[idbh] += p)
        DCMRT[idbh] += dmm[i]
        if idmr > 0 && dbh >= dmrmin
            DCINF[idbh] += p
            DCSUM[idbh] += Float32(idmr) * p
        end
        STVOL += cfv * p
        (idmr > 0 && dbh >= dmrmin) && (STVOLI += cfv * p)
        if dmm[i] > 0f0
            STTPAM += dmm[i]
            STVOLM += cfv * dmm[i]
        end
    end
    ba = Float32(s.plot.basal_area)
    STBAI = STVOL != 0f0 ? ba / STVOL * STVOLI : 0f0
    STBAM = STVOL != 0f0 ? ba / STVOL * STVOLM : 0f0
    STDMR = 0f0; STPIT = 0f0; STPMT = 0f0
    if STTPAX != 0f0                                   # misprt.f: STPIT/STPMT sit inside the STTPAX test
        STDMR = STDMRS / STTPAX
        STPIT = STTPAI / STTPAT * 100f0
        STPMT = STTPAM / STTPAT * 100f0
    end
    STDMI = STTPAI != 0f0 ? STDMRS / STTPAI : 0f0
    STPIV = STVOL != 0f0 ? STVOLI / STVOL * 100f0 : 0f0
    STPMV = STVOL != 0f0 ? STVOLM / STVOL * 100f0 : 0f0
    # DO 265: SPTPAM walks IND1 again
    @inbounds for i in ind1
        sp = Int(t.species[i]); (sp < 1 || sp > maxsp) && continue
        SPTPAM[sp] += dmm[i]
    end

    if isempty(top4)
        # misprt.f LSORT4 block (DO 180): species in index order enter slot 4 when their infected TPA beats it,
        # then a pairwise exchange pass re-sorts the four slots (strict >, so ties keep the lower species index).
        sortsp = zeros(Float32, 4); isv = zeros(Int, 4)
        @inbounds for ispc in 1:maxsp
            if SPTPAI[ispc] > sortsp[4]
                sortsp[4] = SPTPAI[ispc]; isv[4] = ispc
                for i in 1:3, j in (i + 1):4
                    if sortsp[j] > sortsp[i]
                        sortsp[i], sortsp[j] = sortsp[j], sortsp[i]
                        isv[i], isv[j] = isv[j], isv[i]
                    end
                end
            end
        end
        append!(top4, isv)
    end
    species = NamedTuple[]
    # misprt.f:616-632 fills every frozen top-4 slot with INFNO≠0 and DBSMIS1 skips only the '**' (INFNO=0) slots, so a
    # top-4 species whose infection has since dropped to 0 still gets its (zero-infection) row.
    for infno in top4
        infno == 0 && continue
        spdmr = SPTPAX[infno] != 0f0 ? SPDMRS[infno] / SPTPAX[infno] : 0f0
        spdmi = SPTPAI[infno] != 0f0 ? SPDMRS[infno] / SPTPAI[infno] : 0f0
        sppin = SPTPAT[infno] != 0f0 ? SPTPAI[infno] / SPTPAT[infno] * 100f0 : 0f0
        sppmr = SPTPAT[infno] != 0f0 ? SPTPAM[infno] / SPTPAT[infno] * 100f0 : 0f0
        sppoc = STTPAT != 0f0 ? SPTPAT[infno] / STTPAT * 100f0 : 0f0
        push!(species, (sp = infno, mean_dmr = spdmr, mean_dmi = spdmi,
                        inf_tpa = SPTPAI[infno], mort_tpa = SPTPAM[infno],
                        inf_pct = sppin, mort_pct = sppmr, comp_pct = sppoc))
    end

    @inbounds for c in 11:20
        DCTPA[10] += DCTPA[c]; DCTPAX[10] += DCTPAX[c]; DCINF[10] += DCINF[c]
        DCMRT[10] += DCMRT[c]; DCSUM[10] += DCSUM[c]
    end
    dcdmr = zeros(Float32, 10); dcdmi = zeros(Float32, 10)
    @inbounds for c in 1:10
        DCTPAX[c] != 0f0 && (dcdmr[c] = DCSUM[c] / DCTPAX[c])
        DCINF[c]  != 0f0 && (dcdmi[c] = DCSUM[c] / DCINF[c])
    end

    stand = (sttpat = STTPAT, ba = ba, stvol = STVOL, sttpai = STTPAI, stbai = STBAI,
             stvoli = STVOLI, sttpam = STTPAM, stbam = STBAM, stvolm = STVOLM,
             stpit = STPIT, stpiv = STPIV, stpmt = STPMT, stpmv = STPMV,
             stdmr = STDMR, stdmi = STDMI)
    dbhclass = (dctpa = DCTPA[1:10], dcinf = DCINF[1:10], dcmrt = DCMRT[1:10],
                dcdmr = dcdmr, dcdmi = dcdmi)
    return (nage = Int(nage), stand = stand, species = species, dbhclass = dbhclass)
end
