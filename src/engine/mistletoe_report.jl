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
    # misprt.f is REAL*4 throughout; the species sums walk IND1 by species (DO 90/80, ISCT), the stand volumes and
    # DBH classes walk the records (DO ITREE=1,ITRN), and DMMTPA=PROB·DMMORT (mismrt.f:193). jl accumulated in
    # Float64 record order (MEASURED FVSem_g16 196378260020004 mistletoe 2032: Mean_DMR live 2.8820932 = REAL*4,
    # jl 2.8820885801961484).
    SPTPAI = zeros(Float32, maxsp); SPTPAX = zeros(Float32, maxsp)
    SPTPAT = zeros(Float32, maxsp); SPDMRS = zeros(Float32, maxsp)
    SPTPAM = zeros(Float32, maxsp)
    STTPAI = 0f0; STTPAT = 0f0; STTPAX = 0f0; STDMRS = 0f0
    STVOL = 0f0; STVOLI = 0f0; STVOLM = 0f0; STTPAM = 0f0
    # 20 2-inch DBH classes (0-2.9, 3-4.9, ...); the report lumps 11-20 into class 10.
    DCTPA = zeros(Float32, 20); DCTPAX = zeros(Float32, 20); DCINF = zeros(Float32, 20)
    DCMRT = zeros(Float32, 20); DCSUM = zeros(Float32, 20)
    dmrmin = _DM_DMRMIN
    DMMTPA = zeros(Float32, n)
    @inbounds for i in 1:n                                  # mismrt.f:141-195 WKI=PTPA*DMMORT (0 if PROB≤0/IDMR=0)
        sp = Int(t.species[i]); (sp < 1 || sp > maxsp) && continue
        p = t.tpa[i]; idmr = Int(t.dmr[i])
        (p <= 0f0 || idmr == 0) && continue
        DMMTPA[i] = p * _dm_mortality_rate(s, sp, idmr, t.dbh[i], fint)
    end
    ord = _ind1_order(s)
    @inbounds for sp in 1:maxsp                             # misprt.f DO 90 ISPC / DO 80 I3=ISCT(ISPC,1..2)
        any_sp = false
        for i in ord
            Int(t.species[i]) == sp || continue
            any_sp = true
            p = t.tpa[i]; idmr = Int(t.dmr[i])
            if t.dbh[i] >= dmrmin
                idmr > 0 && (SPTPAI[sp] += p)
                SPDMRS[sp] += Float32(idmr) * p
                SPTPAX[sp] += p
            end
            SPTPAT[sp] += p
        end
        any_sp || continue
        STTPAI += SPTPAI[sp]; STTPAT += SPTPAT[sp]; STTPAX += SPTPAX[sp]; STDMRS += SPDMRS[sp]
    end
    # Top-4 most-infected species (misprt.f LSORT4, :448-464): sorted ONCE, on the first call, by the insertion scheme
    # below (strict > — a tie keeps the earlier species), then frozen. Slot 0 ⇒ fewer than 4 infected species.
    if isempty(top4)
        sortsp = zeros(Float32, 4); isv = zeros(Int, 4)
        for sp in 1:maxsp
            if SPTPAI[sp] > sortsp[4]
                sortsp[4] = SPTPAI[sp]; isv[4] = sp
                for a in 1:3, b in (a + 1):4
                    if sortsp[b] > sortsp[a]
                        sortsp[a], sortsp[b] = sortsp[b], sortsp[a]; isv[a], isv[b] = isv[b], isv[a]
                    end
                end
            end
        end
        append!(top4, isv)
    end
    @inbounds for i in 1:n                                  # misprt.f DO ITREE=1,ITRN (stand volumes, DBH classes)
        sp = Int(t.species[i]); (sp < 1 || sp > maxsp) && continue
        p = t.tpa[i]; idmr = Int(t.dmr[i]); dbh = t.dbh[i]; cfv = t.cuft_vol[i]
        ival = unsafe_trunc(Int, dbh); idbh = clamp(ival ÷ 2 + ival % 2, 1, 20)
        DCTPA[idbh] += p
        dbh >= dmrmin && (DCTPAX[idbh] += p)
        DCMRT[idbh] += DMMTPA[i]
        if idmr > 0 && dbh >= dmrmin
            DCINF[idbh] += p
            DCSUM[idbh] += Float32(idmr) * p
        end
        STVOL += cfv * p
        (idmr > 0 && dbh >= dmrmin) && (STVOLI += cfv * p)
        if DMMTPA[i] > 0f0
            STTPAM += DMMTPA[i]
            STVOLM += cfv * DMMTPA[i]
        end
    end
    ba = s.plot.basal_area
    STBAI = 0f0; STBAM = 0f0
    if STVOL != 0f0
        STBAI = ba / STVOL * STVOLI
        STBAM = ba / STVOL * STVOLM
    end
    STDMR = 0f0; STDMI = 0f0; STPIT = 0f0; STPMT = 0f0; STPIV = 0f0; STPMV = 0f0
    if STTPAX != 0f0                                        # misprt.f:578-588
        STDMR = STDMRS / STTPAX
        STPIT = STTPAI / STTPAT * 100f0
        STPMT = STTPAM / STTPAT * 100f0
    end
    STTPAI != 0f0 && (STDMI = STDMRS / STTPAI)
    if STVOL != 0f0
        STPIV = STVOLI / STVOL * 100f0
        STPMV = STVOLM / STVOL * 100f0
    end
    @inbounds for sp in 1:maxsp, i in ord                   # misprt.f DO 265: SPTPAM over IND1
        Int(t.species[i]) == sp && (SPTPAM[sp] += DMMTPA[i])
    end
    species = NamedTuple[]
    for infno in top4
        infno == 0 && continue
        spdmr = SPTPAX[infno] != 0f0 ? SPDMRS[infno] / SPTPAX[infno] : 0f0
        spdmi = SPTPAI[infno] != 0f0 ? SPDMRS[infno] / SPTPAI[infno] : 0f0
        sppin = 0f0; sppmr = 0f0
        if SPTPAT[infno] != 0f0
            sppin = SPTPAI[infno] / SPTPAT[infno] * 100f0
            sppmr = SPTPAM[infno] / SPTPAT[infno] * 100f0
        end
        sppoc = STTPAT != 0f0 ? SPTPAT[infno] / STTPAT * 100f0 : 0f0
        push!(species, (sp = infno, mean_dmr = spdmr, mean_dmi = spdmi,
                        inf_tpa = SPTPAI[infno], mort_tpa = SPTPAM[infno],
                        inf_pct = sppin, mort_pct = sppmr, comp_pct = sppoc))
    end

    # DBH-class table (10 printed classes; 11-20 lumped into 10), all-trees + infected-only DMRs.
    @inbounds for c in 11:20
        DCTPA[10] += DCTPA[c]; DCTPAX[10] += DCTPAX[c]; DCINF[10] += DCINF[c]
        DCMRT[10] += DCMRT[c]; DCSUM[10] += DCSUM[c]
    end
    dcdmr = zeros(Float32, 10); dcdmi = zeros(Float32, 10)
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
