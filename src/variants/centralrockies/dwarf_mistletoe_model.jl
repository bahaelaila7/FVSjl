# =============================================================================
# dwarf_mistletoe_model.jl (centralrockies) — the per-cycle DM spread/intensification
#
# Ported from mistoe/mistoe.f (the MISTOE driver). Runs each cycle for CentralRockies
# AFTER small_tree_growth! and BEFORE mortality (the FVS gradd.f:96 position — after
# GRINCR's growth increments/MORTS, before UPDATE; uses the cycle's HTG). Updates each
# tree's DMR (Hawksworth 0-6) via the Hawksworth spread model, drawing from the shared
# main RNG stream (rann!) in the SAME species-sorted order (ISCT/IND1) as FVS, so the
# stream stays bit-exact-aligned (verified jl[1..105]==FVS[1..105] up to this point).
#
# Gated on infection: a species with stand-mean DMR (SMR) == 0 draws NO random numbers
# (mistoe.f:263) — so mistletoe-free stands are completely untouched (RNG-safe).
#
# Coefficients: spread probs (CR_DM_B*/CR_DM_D*) are variant-uniform (mistoe.f DATA);
# YPLMLT/YNGMLT/DMMMLT default 1.0 (misin0.f, no MISTMULT keyword). DM mortality (MISMRT
# max-combine) and diameter growth-loss (DGPDMR) are wired separately (mortality/DG).
# =============================================================================

"""
    cr_mistoe!(s; fint)

CR dwarf-mistletoe spread/intensification for the cycle (mistoe.f MISTOE). Mutates
`s.trees.dmr`. No-op for non-CentralRockies. Draws `rann!` per host-species tree only
when that species carries infection (SMR>0), matching FVS's draw count/order.
"""
function cr_mistoe!(s::StandState; fint::Float32)
    s.variant isa CentralRockies || return s
    t = s.trees
    t.n == 0 && return s
    species_sort!(s)                              # ensure ISCT/IND1 (species-sorted, FVS mistoe.f ISCT order)
    isct = s.control.sp_count_tab
    ind1 = s.scratch.idx1
    rng  = s.rng
    fscale = fint / 10f0
    @inbounds for ispc in 1:38
        CR_DM_MISFIT[ispc] == 0 && continue          # MISFIT: non-host species skip (no spread)
        i1 = isct[ispc, 1]; i1 == 0 && continue
        i2 = isct[ispc, 2]
        # stand total TPA and mean DMR for this species (mistoe.f:250-262)
        tottpa = 0f0; smr = 0f0
        for i3 in i1:i2
            i = Int(ind1[i3]); p = t.tpa[i]
            tottpa += p; smr += Float32(t.dmr[i]) * p
        end
        tottpa <= 0f0 && continue
        smr /= tottpa
        smr == 0f0 && continue                        # mistletoe-free species ⇒ NO draws (mistoe.f:263)
        # tallest infected tree height per point (DMTALL, mistoe.f:280-290)
        dmtall = Dict{Int32,Float32}()   # only allocated for an actually-infected host species
        for i3 in i1:i2
            i = Int(ind1[i3])
            if t.dmr[i] > 0
                pl = t.plot_id[i]; h = t.height[i]
                (get(dmtall, pl, 0f0) < h) && (dmtall[pl] = h)
            end
        end
        # per-tree spread / intensification (mistoe.f:305-498)
        for i3 in i1:i2
            i = Int(ind1[i3])
            idmr = Int(t.dmr[i])
            htgr10 = t.ht_growth[i]
            pplus = 0f0
            if idmr < 6
                pplus = CR_DM_BCONST + CR_DM_BDMR[idmr + 1] +
                        CR_DM_BHTG * (htgr10 * 10f0 / fint) + CR_DM_BTPA * tottpa
                pplus != 0f0 && (pplus = 1f0 / (1f0 + fexp(-pplus)))     # ×YPLMLT(=1)
                pplus = pplus >= 1f0 ? 1f0 : 1f0 - fpow(1f0 - pplus, fscale)
            end
            if idmr != 0
                pminus = CR_DM_DCONST + CR_DM_DDMR * idmr +
                         CR_DM_DHTG * (htgr10 * 10f0 / fint) + CR_DM_DTPA * tottpa
                pminus != 0f0 && (pminus = 1f0 / (1f0 + fexp(-pminus)))  # ×YNGMLT(=1)
                pminus = pminus >= 1f0 ? 1f0 : 1f0 - fpow(1f0 - pminus, fscale)
                xnum = rann!(rng)
                dtall = get(dmtall, t.plot_id[i], 0f0)
                if idmr != 6 && pplus > xnum
                    if dtall * 0.7f0 > t.height[i]
                        x2 = rann!(rng)                                  # overstory intensification magnitude
                        m = Int(t.dmr[i])
                        inc = if m == 1
                            x2 < 0.61f0 ? 1 : (x2 < 0.83f0 ? 2 : 3)
                        elseif m == 2 || m == 3
                            x2 < 0.34f0 ? 1 : (x2 < 0.67f0 ? 2 : 3)
                        elseif m == 4
                            x2 < 0.55f0 ? 1 : 2
                        else
                            1
                        end
                        t.dmr[i] += Int32(inc)
                    else
                        t.dmr[i] += Int32(1)
                    end
                    t.dmr[i] > 6 && (t.dmr[i] = Int32(6))
                end
                (pminus > xnum) && (t.dmr[i] -= Int32(1))               # decrease (mistoe.f:462)
            else
                xnum = rann!(rng)                                        # uninfected: introduce infection
                dtall = get(dmtall, t.plot_id[i], 0f0)
                if dtall * 0.7f0 > t.height[i]
                    if xnum < 0.55f0
                        x2 = rann!(rng)
                        t.dmr[i] = x2 < 0.69f0 ? Int32(1) : (x2 < 0.87f0 ? Int32(2) : Int32(3))
                    end
                else
                    (pplus > xnum) && (t.dmr[i] = Int32(1))
                end
            end
        end
        # highly-infected trees flagged (mistoe.f:503-507 IMC=3)
        for i3 in i1:i2
            i = Int(ind1[i3])
            t.dmr[i] >= 4 && (t.mort_code[i] = Int32(3))
        end
    end
    return s
end

"""
    cr_dm_mortality_combine!(killed, s, fint, n)

MAX-combine the dwarf-mistletoe-induced per-tree mortality into `killed[]` (mismrt.f:185-191:
`WKI = PTPA·DMMORT; IF WK2<WKI, WK2=WKI`). Uses the post-spread DMR and start-of-cycle DBH.
No-op for non-CentralRockies. Order-independent (per-tree max), so iterate physical order.
"""
function cr_dm_mortality_combine!(killed::AbstractVector{Float32}, s::StandState, fint::Float32, n::Int)
    s.variant isa CentralRockies || return
    t = s.trees
    @inbounds for i in 1:n
        dmr = Int(t.dmr[i]); dmr == 0 && continue
        pr = t.tpa[i]; pr <= 0f0 && continue
        rate = cr_dm_mortality_rate(Int(t.species[i]), dmr, t.dbh[i], fint)
        wki = pr * rate
        killed[i] < wki && (killed[i] = wki)
    end
    return
end

"""
    cr_dm_growth_loss!(s, stash)

Multiply each infected tree's diameter growth by its DM growth-potential proportion
(misdgf.f: `DG = DG·DGPDMR(sp,DMR+1)`), for the central record AND the tripled sub-records
(dgdriv.f:230/276/290). Uses the START-of-cycle DMR (growth runs before this cycle's spread,
FVS GRINCR before GRADD's MISTOE). No-op for non-CentralRockies. Deterministic (no RNG).
"""
function cr_dm_growth_loss!(s::StandState, stash)
    s.variant isa CentralRockies || return
    t = s.trees
    n = stash === nothing ? t.n : stash.nlive
    @inbounds for i in 1:n
        dmr = Int(t.dmr[i]); dmr == 0 && continue
        m = Float32(cr_dm_dg_mult(Int(t.species[i]), dmr))
        m == 1f0 && continue
        t.diam_growth[i] *= m
        stash !== nothing && (stash.dgU[i] *= m; stash.dgL[i] *= m)
    end
    return
end
