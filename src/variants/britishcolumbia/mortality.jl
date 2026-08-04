# =============================================================================
# mortality.jl (britishcolumbia) — BC mortality driver (canada/bc/morts.f). CHUNK 7.
#
# DUAL-REGIME (LV2ATV, becset.f:855):
#   • V2 (non-ICH/IDF/SBS/SBPS zones): shared Hamilton RIP regression (like KT/IE). NOT PORTED here —
#     errors loudly; needs a V2 fixture to validate (all_BC is V3/ICH). TODO.
#   • V3 (ICH/IDF/SBS/SBPS incl. all_BC): per-species branch —
#       – TABULAR species (PW LW FD BG HW CW PL SE PY; sp3/sp14→FD, sp4/sp9→BG): RIP = TMRT(mrtcls,ba,dbh)
#         = FMRT(beccls,mrtcls,ba,dbh)/100  (BECCLS: ICH/SBS→1, IDF/SBPS/SBSdw2→2; MORCON morts.f:1109-1117);
#       – LMRT hardwoods (EP AT AC OH = sp11/12/13/15): logistic survival→mortality (LMRT(1) coefs);
#       – else (MRTCLS=0, not in LMRT): RIP = 99 (complete mortality) — no BC species hits this.
# The BAMAX-approach tail (RZ / RIPP / SIZCAP / climate-death SDIMAX<5) is IDENTICAL to the shared
# western driver (mirrors KT src/variants/kootenai/mortality.jl). Reuses the shared apply-tail
# (book_mortality_snags! + TPA removal). ⚠ Establishment/FIXMORT/MORTMULT windows omitted (X=1).
# =============================================================================

# --- V3 tabular mortality table (morts.f:207-354). Annual %/yr; /100 → annual proportion in TMRT. ---
# BC_FMRT[bec, mrtcls, ba(1-5), dbh(1-4)];  bec: 1=ICH, 2=IDF.  mrtcls: 1PW 2LW 3FD 4BG 5HW 6CW 7PL 8SE 9PY.
const BC_DBHBRK = Float32[10f0, 25f0, 40f0]           # cm DBH-class breaks (DX ≤ brk)
const BC_BABRK  = Float32[20f0, 30f0, 40f0, 50f0]     # m²/ha BA-class breaks (BAX ≤ brk)

# One 5×4 (BA×DBH) block per (bec,spp). Rows = BA class, cols = DBH class.
const _BC_FMRT_BLOCKS = Dict{Tuple{Int,Int},Matrix{Float32}}(
    # ICH,PW (IDF PW = ICH)
    (1,1) => Float32[0.153 0.141 0.088 0.062; 0.305 0.281 0.176 0.124; 0.305 0.281 0.176 0.124; 0.364 0.266 0.188 0.124; 0.462 0.266 0.188 0.124],
    # ICH,LW (IDF LW = ICH)
    (1,2) => Float32[0.196 0.185 0.131 0.058; 0.393 0.370 0.263 0.115; 0.393 0.350 0.263 0.115; 0.456 0.334 0.263 0.115; 0.456 0.334 0.263 0.115],
    # ICH,FD
    (1,3) => Float32[0.281 0.072 0.136 0.123; 0.563 0.145 0.272 0.245; 0.592 0.297 0.272 0.245; 0.400 0.326 0.181 0.326; 0.159 0.385 0.233 0.252],
    # IDF,FD (distinct)
    (2,3) => Float32[0.265 0.259 0.305 0.213; 0.529 0.518 0.609 0.427; 0.515 0.443 0.332 0.358; 0.409 0.647 0.689 0.163; 0.409 0.510 0.247 0.163],
    # ICH,BG(BL) (IDF = ICH)
    (1,4) => Float32[0.103 0.108 0.064 0.020; 0.205 0.216 0.127 0.039; 0.480 0.392 0.216 0.039; 0.510 0.223 0.187 0.079; 0.202 0.092 0.070 0.070],
    # ICH,HW (IDF = ICH)
    (1,5) => Float32[0.125 0.101 0.076 0.045; 0.251 0.202 0.153 0.089; 0.198 0.330 0.210 0.089; 0.288 0.330 0.210 0.089; 0.394 0.200 0.151 0.103],
    # ICH,CW (IDF = ICH)
    (1,6) => Float32[0.078 0.068 0.065 0.059; 0.157 0.136 0.130 0.119; 0.225 0.142 0.136 0.130; 0.295 0.223 0.151 0.151; 0.202 0.092 0.070 0.049],
    # ICH,PL
    (1,7) => Float32[0.218 0.317 0.196 0.116; 0.435 0.633 0.391 0.232; 0.457 0.607 0.276 0.232; 0.391 0.468 0.405 0.232; 0.256 0.388 0.128 0.232],
    # IDF,PL (distinct)
    (2,7) => Float32[0.282 0.226 0.130 0.130; 0.563 0.452 0.260 0.260; 0.493 0.445 0.195 0.195; 0.413 0.592 0.245 0.245; 0.413 0.287 0.433 0.245],
    # ICH,SE (combined ICH+IDF)
    (1,8) => Float32[0.126 0.101 0.122 0.058; 0.253 0.202 0.244 0.116; 0.253 0.202 0.244 0.116; 0.253 0.242 0.242 0.116; 0.219 0.182 0.152 0.116],
    # ICH,PY (same as IDF)
    (1,9) => Float32[0.168 0.116 0.089 0.063; 0.335 0.231 0.179 0.127; 0.320 0.222 0.092 0.046; 0.662 0.434 0.206 0.220; 0.662 0.434 0.206 0.220],
)

"""Assemble BC_FMRT[bec,spp,ba,dbh] from the per-block table, filling IDF=ICH where IDF is not distinct."""
function _bc_build_fmrt()
    fmrt = zeros(Float32, 2, 9, 5, 4)
    for spp in 1:9, bec in 1:2
        blk = get(_BC_FMRT_BLOCKS, (bec, spp), _BC_FMRT_BLOCKS[(1, spp)])  # IDF → ICH fallback
        @inbounds for iba in 1:5, idbh in 1:4
            fmrt[bec, spp, iba, idbh] = blk[iba, idbh]
        end
    end
    fmrt
end
const BC_FMRT = _bc_build_fmrt()

# species → tabular MRTCLS (morts.f:512-530); 0 = not tabular (LMRT hardwood or complete-mort).
const BC_MORT_MRTCLS = Int[1, 2, 3, 4, 5, 6, 7, 8, 4, 9, 0, 0, 0, 3, 0]  # sp3/14→FD(3), sp4/9→BG(4)

# LMRT(1) logistic (morts.f:366-386, applies to sp 11,12,13,15 in all V3 zones).
# COEF order (morcon.f copy): CON, INVDBH, DBH, DBHSQ, RDBH, BAL, SQRTBA, SPH.
const BC_LMRT_SPP  = (11, 12, 13, 15)
const BC_LMRT_COEF = Float32[-2.774444f0, 0.684626f0, -0.180211f0, 0.004550f0,
                             -0.634960f0, 0.063153f0, -0.276905f0, 0.000100f0]
@inline _bc_lmrt_fit(sp::Integer) = sp in BC_LMRT_SPP

"""BEC-class for the tabular table (MORCON morts.f:1100-1107): ICH/SBS→1, IDF/SBPS→2, SBSdw2→2."""
function bc_beccls(zone::AbstractString, series::AbstractString)
    bec = 1
    occursin("ICH", zone)  && (bec = 1)
    occursin("IDF", zone)  && (bec = 2)
    occursin("SBPS", zone) && (bec = 2)
    occursin("SBS", zone)  && (bec = 1)
    occursin("SBSdw2", series) && (bec = 2)
    bec
end

"""LV2ATV (becset.f:855): V3 (false) for ICH/IDF/SBS/SBPS zones, else V2 (true)."""
function bc_lv2atv(zone::AbstractString)
    v3 = occursin("ICH", zone) || occursin("IDF", zone) ||
         occursin("SBS", zone) || occursin("SBPS", zone)
    !v3
end

"""
    mortality!(s, ::BritishColumbia; fint, book_snags)

BC mortality (canada/bc/morts.f). V3 per-species RIP (tabular FMRT / LMRT logistic) → RIPP (BAMAX/RZ
weighting) → WKI = P·(1−(1−RIPP)^FINT); climate-death when weighted SDImax < 5. V2 Hamilton path errors
(not yet ported — needs a V2 fixture). Reuses the shared apply-tail.
"""
function mortality!(s::StandState, ::BritishColumbia; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    zone, series = bc_stand_zone(s)
    bc_lv2atv(zone) && error("BC mortality: V2 (Hamilton) regime not yet ported (zone=$zone); " *
                             "only V3 (ICH/IDF/SBS/SBPS) validated. Needs a V2 fixture.")
    beccls = bc_beccls(zone, series)
    ba = p.basal_area
    sdimax = stand_sdimax(s)
    # BAMAX: user (control.ba_max) or SDICAL Stage default = SDImax·0.5454154·PMSDIU (PMSDIU=0.85). sdical.f:204
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : sdimax * 0.5454154f0 * 0.85f0
    bamax <= 0f0 && (bamax = 1f0)
    # stand sums (morts.f:436-472): T, SD2SQ → DQ10; AVED (BA-weighted mean DBH, inches)
    tt = 0f0; sd2sq = 0f0; dsum = 0f0; wprob = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bc_bratio(sp)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + 2f0 * d * g + g * g); tt += pr
        wprob += pr; dsum += d * pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt)
    deltba = 0.005454154f0 * dq10 * dq10 * tt - ba
    ba10 = ba + (bamax - ba) / bamax * deltba
    tb = ba10 / (0.005454154f0 * dq10 * dq10)
    ttb = (tt - tb) / tt; ttb > 0.9999f0 && (ttb = 0.9999f0)
    rz = 1f0 - (1f0 - ttb)^0.1f0
    aved = dsum / wprob
    # metric stand values (morts.f:476-492): BAX (m²/ha), SPH (stems/ha), BA-class
    bax = ba * BC_FT2pACRtoM2pHA
    sph = tt * 2.471f0                                    # HAtoACR
    bacls = 5
    @inbounds for i in 1:4
        if bax <= BC_BABRK[i]; bacls = i; break; end
    end
    sqrtbax = sqrt(bax)
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    sc = s.control.sp_size_cap
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]; bark = bc_bratio(sp)
        dx = max(2.5f0, d * BC_INtoCM)                    # cm, ≥2.5
        reldbh = dx / max(2.5f0, aved * BC_INtoCM)
        mrtcls = BC_MORT_MRTCLS[sp]
        if mrtcls > 0                                     # TABULAR (morts.f:594-601)
            dbhcls = 4
            for j in 1:3
                if dx <= BC_DBHBRK[j]; dbhcls = j; break; end
            end
            rip = BC_FMRT[beccls, mrtcls, bacls, dbhcls] / 100f0
        elseif _bc_lmrt_fit(sp)                           # LMRT logistic (morts.f:610-625)
            bal = (1f0 - t.crown_ratio[i] / 100f0) * bax
            c = BC_LMRT_COEF
            rip = c[1] + c[2] / dx + c[3] * dx + c[4] * dx * dx +
                  c[5] * reldbh + c[6] * bal + c[7] * sqrtbax + c[8] * sph
            rip = max(-70f0, min(70f0, rip))
            rip = 1f0 - (1f0 / (1f0 + exp(rip)))
        else                                              # no estimate → complete mortality
            rip = 1f0 - (1f0 / (1f0 + exp(min(70f0, 99f0))))
        end
        # BAMAX-approach RIPP (morts.f:632-637) — identical to shared western driver
        ripp = ba * rz
        ba <= bamax && (ripp += (bamax - ba) * rip)
        ripp /= bamax
        ripp < rip && (ripp = rip); ripp > 1f0 && (ripp = 1f0)
        wki = pr * (1f0 - (1f0 - ripp)^fint)              # X=1 (no MORTMULT/estab window)
        gsc = (t.diam_growth[i] / bark) * (fint / 10f0)
        if (d + gsc) >= sc[sp, 1] && trunc(Int, sc[sp, 3]) != 1
            wki = max(wki, pr * sc[sp, 2] * fint / 10f0)
        end
        wki > pr && (wki = pr)
        sdimax < 5f0 && (wki = pr)                        # climate death (morts.f:676-678)
        killed[i] = wki
    end
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
