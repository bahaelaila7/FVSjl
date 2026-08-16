# =============================================================================
# organon_hook_op.jl — OP (Olympic) COOPERATING ORGANON/Wykoff growth driver.
#
# OP couples ORGANON-NWO (grows the IORG=1 valid-ORGANON trees) with op-native Wykoff growth (grows
# the IORG=0 trees). Faithful to op/dgdriv.f: DGF loads WK2=ln(DDS) for EVERY tree on the ORIGINAL
# DBH, then the ORGANON EXECUTE overwrites WK2 for the IORG=1 trees (op/dgdriv.f:432-441), and the
# shared FVS diameter path grows everyone from WK2. This is COOPERATING HOOKS (each of DG/HT/crown
# runs natively for IORG=0 and consumes stashed ORGANON outputs for IORG=1); MORTALITY is the one
# exception — op/morts.f uses ORGANON MORTEXP for ALL trees whenever ORGANON ran (SMORMT>0).
#
#   op_organon_prepare!(s)          — LSTART PREPARE (op/cratet.f) : dub missing HT/CR for the
#                                     valid-ORGANON trees, compute + store ACALIB. Mirrors OC.
#   diameter_growth!(s,::Olympic)   — the per-cycle DGDRIV driver: op_dgcons!+dgf! (IORG=0 WK2),
#                                     op_execute_nwo (IORG=1), fold DGRO→WK2, DDS→DG for all + DGBND;
#                                     stash HGRO/CR2/MORTEXP on s.calib for the other hooks.
#
# OP is DGSD=0 AND ICL4=0 (op/grinit.f) ⇒ FULLY DETERMINISTIC: no OLDRN serial-correlation, no record
# tripling, and DGSCOR returns FRM=1 (base/dgscor.f, DGSD<1). So the DDS→DG conversion is the plain
# DG=sqrt(d_ib²+exp(WK2+ln XDMULT))−d_ib with d_ib=DBH·op_bratio. Validated vs live FVSop_clean/g16.
# =============================================================================

"op/dgdriv.f RVARS: SITE_1=SITEAR(16=DF), SITE_2=SITEAR(19=WH). Returns the RAW site indices."
@inline function _op_organon_site(s::StandState)
    si = s.plot.sp_site_index
    si1 = length(si) >= 16 ? si[16] : 0f0
    si2 = length(si) >= 19 ? si[19] : 0f0
    return si1, si2
end

"op/sitset.f:325-327 RVARS(3/4/5) = SDIDEF(16=DF), SDIDEF(3=GF), SDIDEF(19=WH) (MSDI_1/2/3)."
@inline function _op_organon_msdi(s::StandState)
    sd = s.plot.sp_sdi_def
    m1 = length(sd) >= 16 ? sd[16] : 0f0
    m2 = length(sd) >= 3  ? sd[3]  : 0f0
    m3 = length(sd) >= 19 ? sd[19] : 0f0
    return m1, m2, m3
end

"op DGBND (op/dgbnd.f): DGMAX=7.92·exp(−0.03·min(DBH,150)) cap (RW sp17 exempt), then the SIZCAP cap."
@inline function op_dgbnd(sp::Int, dbh::Float32, ddg::Float32, sizcap::AbstractMatrix)::Float32
    if sp != 17
        temdbh = dbh > 150f0 ? 150f0 : dbh
        dgmax = 7.92f0 * exp(-0.03f0 * temdbh)
        ddg > dgmax && (ddg = dgmax)
        ddg < 0f0 && (ddg = 0f0)
    end
    if (dbh + ddg) > sizcap[sp, 1] && sizcap[sp, 3] < 1.5f0
        ddg = sizcap[sp, 1] - dbh
        ddg < 0.01f0 && (ddg = 0.01f0)
    end
    return ddg
end

"""
    op_organon_prepare!(s)

op/cratet.f LSTART ORGANON PREPARE (VERSION=2) — the setup-time pass that (1) gates on a big-6 tree,
(2) marshals ALL live records into the /ORGANON/ PREPARE buffer, (3) calls `op_prepare_nwo` to dub
missing HT/CR for the valid-ORGANON (IORG=1) records and compute the ACALIB calibration, and (4)
stores ACALIB on `s.calib.organon_acalib` for the growth path. No-op unless the stand has a big-6
(GF/DF) tree. Mirrors `oc_organon_prepare!`. Deterministic (DGSD=0). Even-aged (op/grinit.f INDS(4)=1).
"""
function op_organon_prepare!(s::StandState)
    t = s.trees; n = t.n
    n == 0 && return nothing
    # --- CRATET setup-time eligibility gate (op/cratet.f:160-209) ---
    iorg = zeros(Int32, n); nbig6 = 0
    @inbounds for i in 1:n
        sp = Int(t.species[i]); h = t.height[i]
        ihflag = (h == 0f0) || (h > 4.5f0)                 # measured-HT lower limit OR missing (cratet.f:167)
        if t.dbh[i] >= 0.1f0 && ihflag
            (sp in OP_ORGANON_BIG6) && (nbig6 += 1)
            iorg[i] = (sp in OP_ORGANON_VALID) ? Int32(1) : Int32(0)
        end
    end
    nbig6 == 0 && return nothing                            # no big-6 ⇒ PREPARE skipped (cratet.f:204-208)
    # --- marshal ALL live records into the PREPARE buffer (op/cratet.f:216-241) ---
    species = Vector{Int32}(undef, n); dbh1 = Vector{Float32}(undef, n)
    ht1or   = Vector{Float32}(undef, n); cr1 = Vector{Float32}(undef, n)
    expan1  = Vector{Float32}(undef, n); radgro = zeros(Float32, n)
    pival = s.plot.pi > 0f0 ? s.plot.pi : 1f0               # PI = FLOAT(IPTINV)
    @inbounds for i in 1:n
        species[i] = op_organon_fia(Int(t.species[i]))      # ORGSPC
        d = t.dbh[i]; d < 0.1f0 && (d = 0.1f0); dbh1[i] = d
        h = t.height[i]; (h > 0f0 && h < 4.6f0) && (h = 4.6f0)   # floor to 4.6 ONLY when HT>0 (cratet.f:223)
        ht1or[i] = h
        cr1[i] = Float32(t.crown_pct[i]) / 100f0
        expan1[i] = t.tpa[i] * pival                        # EXPAN1 = PROB·PI (cratet.f:225)
        radgro[i] = t.diam_growth[i] / 2f0                  # RADGRO = DG/2 (cratet.f:226)
    end
    stage = Int(s.plot.stand_age)                           # STAGE = IAGE (ICYC=1)
    bhage = stage - 6                                       # BREAST HEIGHT AGE (cratet.f:249)
    si_1, si_2 = _op_organon_site(s)                        # RAW SITEAR(16)/(19); PREPARE's own DF↔WH conv handles ≤0
    npts = max(1, Int(round(pival)))
    res = op_prepare_nwo(species, dbh1, ht1or, cr1, expan1, n, npts, stage, bhage, si_1, si_2)
    # --- write ORGANON-dubbed HT/CR back into the valid-ORGANON records that were MISSING them
    #     (op/cratet.f:338-354: only IORG=1 trees reloaded) ---
    @inbounds for i in 1:n
        iorg[i] == 1 || continue
        (t.height[i] <= 0f0 && res.ht[i] > 0f0)   && (t.height[i]   = res.ht[i])
        (t.crown_pct[i] <= 0 && res.cr[i] > 0f0)  && (t.crown_pct[i] = round(Int32, res.cr[i] * 100f0, RoundNearestTiesAway))
    end
    copyto!(s.calib.organon_acalib, res.acalib)             # ACALIB → growth path (cratet.f:382-390)
    return nothing
end

"""
    diameter_growth!(s, ::Olympic; sfint=5, tripling=false, kwargs...) -> nothing

op/dgdriv.f per-cycle DGDRIV — the COOPERATING diameter driver. Fills `t.diam_growth` for EVERY tree
from WK2 (op-native ln(DDS) for IORG=0; ORGANON DGRO-folded ln(DDS) for IORG=1) and stashes the
ORGANON HGRO/CR2/MORTEXP on `s.calib` for height_growth!/crown_ratio_update!/mortality!. Returns
`nothing` (no tripling: ICL4=0 on OP) so the shared `triple_records!` no-ops.
"""
function diameter_growth!(s::StandState, ::Olympic; sfint::Float32 = 5f0,
                          tripling::Bool = false, kwargs...)
    t, c, p = s.trees, s.calib, s.plot
    n = t.n
    # COR attenuation toward the calibration goal (dgdriv.f:508-513); dgf! reads c.dg_cor. Elapsed-time
    # clock from the IY schedule (cormlt=1 at the first projection cycle). Inert (COR=0) when no calibration.
    elapsed = Float32(current_cycle_year(s) - Int(s.control.cycle_year[1]))
    cormlt = exp(-0.02773f0 * elapsed)
    @inbounds for sp in 1:MAXSP
        c.dg_cor[sp] = c.dg_cor_goal[sp] + cormlt * c.dg_cor_goal[sp]
    end
    op_dgcons!(s)                                           # ENTRY DGCONS — per-species site DGCON
    dgf!(s, s.variant)                                      # WK2 = op-native ln(DDS) for IORG=0 trees
    wk2 = view(s.scratch.wk, 2, :)
    buf = op_build_organon_buffer!(s)                       # IORG gate + big-6 gate + /ORGANON/ buffer
    si1r, si2r = _op_organon_site(s)
    m1, m2, m3 = _op_organon_msdi(s)
    cyclg = Int(s.control.cycle)                            # CYCLG = ICYC-1 = jl cycle (0-based)
    isp_fvs = Int[Int(t.species[i]) for i in 1:n]
    calib1 = Float32[s.calib.organon_acalib[1, g] for g in 1:11]
    calib2 = Float32[s.calib.organon_acalib[2, g] for g in 1:11]
    g = buf.runs ? op_execute_nwo(buf, isp_fvs; si_1 = si1r - 4.5f0, si_2 = si2r - 4.5f0,
                        msdi_1 = m1, msdi_2 = m2, msdi_3 = m3, calib1 = calib1, calib2 = calib2,
                        cyclg = cyclg) : nothing
    # Fold the ORGANON DGRO→DDS into WK2 for the IORG=1 trees (op/dgdriv.f:432-441).
    if g !== nothing
        @inbounds for i in 1:n
            buf.iorg[i] == 1 && (wk2[i] = g.dds[i])
        end
    end
    # DDS→DG for EVERY tree (op/dgdriv.f:521-554): DG=sqrt(d_ib²+EXP(WK2+ln XDMULT))−d_ib, then DGBND.
    sizcap = s.control.sp_size_cap
    cur_year = current_cycle_year(s)
    @inbounds for i in 1:n
        d = t.dbh[i]
        if d <= 0f0
            t.diam_growth[i] = 0f0; continue
        end
        sp = isp_fvs[i]
        xbai = active_multiplier(s.control, :bai, sp, cur_year)
        xdgrow = flog(xbai)
        bark = op_bratio(sp, d)
        dib = d * bark
        dds = fexp(wk2[i] + xdgrow)
        dg = sqrt(dib * dib + dds) - dib
        t.diam_growth[i] = op_dgbnd(sp, d, dg, sizcap)
    end
    # Stash the ORGANON per-tree outputs for the cooperating hooks (op has no tripling ⇒ index-stable).
    if length(c.op_iorg) != n
        resize!(c.op_iorg, n); resize!(c.op_hgro, n); resize!(c.op_cr2, n); resize!(c.op_mortexp, n)
    end
    @inbounds for i in 1:n
        c.op_iorg[i] = (g !== nothing) ? buf.iorg[i] : Int32(0)
        c.op_hgro[i] = (g !== nothing) ? g.hgro[i] : 0f0
        c.op_cr2[i]  = (g !== nothing) ? g.cr2[i]  : 0f0
        c.op_mortexp[i] = (g !== nothing) ? g.deadexp[i] : 0f0
    end
    c.op_org_ran = (g !== nothing)
    return nothing
end
