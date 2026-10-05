# =============================================================================
# regent.jl (eastcascades) — EC small-tree growth (ec/regent.f + ec/smhtgf.f) + htdbh. Chunk 6.
#
# EC's small-tree model is BESPOKE (not the WC smhgdg reuse): ec/smhtgf.f is HEIGHT-growth-only with
# inline per-species Chapman-Richards/curve literals (WP/RC Chapman-Richards; MH/OS ×3.281 metric), and
# small-tree DIAMETER growth lives in ec/regent.f (per-species H→DBH inversion bounded by DGMAX·SCALE,
# SLO/SHI SI clamp, no DGMIN). Missing-height dubbing + small-tree DBH inversion use the forest-dependent
# Curtis-Arney htdbh (ec/htdbh.f: MTHOOD/OKANOG/GIFFPC/WENATC by IFOR/IGL).
# =============================================================================

# ec/htdbh.f — 4 physical forest tables (MTHOOD/OKANOG/GIFFPC/WENATC), stored in htdbh_coeffs_ec.csv by
# resolved CSV forest_idx {1=MTHOOD, 2/4/7=OKANOG, 3=WENATC(IGL1), 5/6=GIFFPC(IGL2,3)} × 32 species.
let
    P2 = zeros(Float32, 7, 32); P3 = zeros(Float32, 7, 32); P4 = zeros(Float32, 7, 32)
    for l in readlines(joinpath(EC_DATADIR, "htdbh_coeffs_ec.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        fi = parse(Int, f[1]); sp = parse(Int, f[2])
        P2[fi, sp] = parse(Float32, f[3]); P3[fi, sp] = parse(Float32, f[4]); P4[fi, sp] = parse(Float32, f[5])
    end
    global const EC_HTDBH_P2 = P2; global const EC_HTDBH_P3 = P3; global const EC_HTDBH_P4 = P4
end

# ec/htdbh.f SELECT CASE(IFOR)/CASE(IGL) → CSV forest_idx. forkod passes the CORRECTED IFOR (1..4) + IGL.
@inline function _ec_htdbh_fidx(ifor::Int, igl::Int)::Int
    ifor == 1 && return 1                 # MTHOOD
    (ifor == 2 || ifor == 4) && return 2  # OKANOG
    (ifor == 3 && (igl == 2 || igl == 3)) && return 5  # GIFFPC
    return 3                              # WENATC (IFOR 3, IGL 1) / default
end

@inline function ec_htdbh_height(fidx::Int, sp::Int, d::Float32)::Float32
    p2 = EC_HTDBH_P2[fidx, sp]; p3 = EC_HTDBH_P3[fidx, sp]; p4 = EC_HTDBH_P4[fidx, sp]
    if d >= 3.0f0
        return 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(d, p4))
    else
        return ((4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3.0f0, p4)) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
    end
end
@inline function ec_htdbh_dbh(fidx::Int, sp::Int, h::Float32)::Float32
    p2 = EC_HTDBH_P2[fidx, sp]; p3 = EC_HTDBH_P3[fidx, sp]; p4 = EC_HTDBH_P4[fidx, sp]
    hat3 = 4.5f0 + p2 * fexp(-1f0 * p3 * fpow(3.0f0, p4))
    if h >= hat3
        return fexp(flog((flog(h - 4.5f0) - flog(p2)) / (-1f0 * p3)) * 1f0 / p4)   # htdbh.f:354 ALOG(..)*1./P4 = (x*1.)/P4
    else
        return (((h - 4.51f0) * 2.7f0) / (hat3 - 4.51f0)) + 0.3f0
    end
end

@inline ec_htdbh_ifor(p)::Int = _ec_htdbh_fidx(Int(p.forest_idx), Int(p.geo_location))

regenerate!(s::StandState, ::EastCascades; kwargs...) = s
# mortality!(::EastCascades) = the shared Stage/Reineke driver mortality!(::AbstractVariant)
# (southern/mortality.jl); EC's only variant hook is _varmrt_efftr!(::EastCascades) in mortality.jl.


# ── ec/regent.f + ec/smhtgf.f small-tree growth (chunk 6).
const EC_RG_DGMAX = Float32[2.8,2.8,2.4,3.6,2.5,2.5,3.5,3.6,3.6,2.8,5.0,2.8,5.0,5.0,5.0,2.5,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,2.8,5.0]
const EC_RG_XMAX  = Float32[4.0,4.0,4.0,4.0,10.0,4.0,5.0,4.0,6.0,6.0,4.0,6.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,6.0,4.0]
const EC_RG_XMIN  = Float32[2.0,2.0,2.0,2.0,2.0,2.0,1.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0]
const EC_RG_SLO   = Float32[20,50,50,50,15,50,30,40,50,70,0,15,0,0,0,50,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,0]
const EC_RG_SHI   = Float32[80,110,110,110,30,110,70,120,150,140,999,30,999,999,999,110,999,999,999,999,999,999,999,999,999,999,999,999,999,999,30,999]
const EC_RG_DIAM  = Float32[0.4,0.3,0.3,0.3,0.2,0.3,0.4,0.3,0.3,0.5,0.2,0.2,0.2,0.4,0.3,0.3,0.3,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
const EC_RG_AB    = Float32[1.11436,-0.011493,0.43012e-4,-0.72221e-7,0.5607e-10,-0.1641e-13]   # ec/regent.f:89 REAL AB(9) — single-precision Horner
const EC_RG_REGYR = 10.0f0

# ec/smhtgf.f — small-tree potential HEIGHT increment over DTIME years (MODE 1 = established, effective-age).
@inline function ec_smhtgf(ispc::Int, h::Float32, dtime::Float32, si::Float32)::Float32
    if ispc == 1                                       # WP — Chapman-Richards
        c1=0.375045f0; c2=0.92503f0; c3=-0.020796f0; c4=2.48811f0
        effage = flog((1f0 - fpow(c1/si*h, 1f0/c4))/c2)/c3; agepdt = effage + dtime
        return (si/c1)*fpow(1f0-c2*fexp(c3*agepdt), c4) - (si/c1)*fpow(1f0-c2*fexp(c3*effage), c4)
    elseif ispc == 5                                   # RC — Chapman-Richards
        c1=0.752842f0; c2=1.0f0; c3=-0.0174f0; c4=1.4711f0
        e = (1f0 - fpow(c1/si*h, 1f0/c4))/c2
        effage = e > 0f0 ? flog(e)/c3 : 100f0; agepdt = effage + dtime
        return (si/c1)*fpow(1f0-c2*fexp(c3*agepdt), c4) - (si/c1)*fpow(1f0-c2*fexp(c3*effage), c4)
    elseif ispc == 2 || ispc == 17
        return ((-3.97245f0 + 0.50995f0*si)/(28.11668f0-0.05661f0*si))*dtime
    elseif ispc == 3
        return ((2.0f0 + 0.420f0*si)/(28.5f0 - 0.05f0*si))*dtime
    elseif ispc == 4
        return ((-0.6667f0 + 0.4333f0*si)/(28.5f0 - 0.05f0*si))*dtime
    elseif ispc == 6 || ispc == 16
        return ((-1.0470f0 + 0.4220f0*si)/(28.7739f0 - 0.0597f0*si))*dtime
    elseif ispc == 7
        return (0.3277f0 + 0.01296f0*si)*dtime
    elseif ispc == 8
        return ((-8.0f0 + 0.35f0*si)/(53.72545f0 - 0.274509f0*si))*dtime
    elseif ispc == 9
        return ((6.0f0 + 0.14f0*si)/(33.882f0 - 0.06588f0*si))*dtime
    elseif ispc == 10
        return ((-1.0f0 + 0.32857f0*si)/(28.0f0 - 0.042857f0*si))*dtime
    elseif ispc == 11
        return ((-5.74874f0 + 0.54576f0*si)/(26.15767f0 - 0.03596f0*si))*dtime
    elseif ispc == 12 || ispc == 31
        return ((0.965758f0 + 0.082969f0*si)/(55.249612f0 - 1.288852f0*si))*dtime*3.280833f0
    elseif ispc == 15
        return ((11.26677f0 + 0.12027f0*si)/(27.93806f0 - 0.02873f0*si))*dtime
    elseif ispc == 22
        return (-0.007025f0 + 0.056794f0*si)*dtime
    elseif ispc == 28
        return (((-37.60812f0*flog(1f0-fpow(si/114.24569f0, 0.44444f0)))*0.01f0)-0.1f0)*dtime
    else                                               # 13,14,18:21,23:27,29,30,32
        return ((1.47043f0 + 0.23317f0*si)/(31.56252f0 - 0.05586f0*si))*dtime
    end
end

# ec/essubh.f — planted/subsequent-tree base height: SMHTGF called with MODE=0, DTIME=AGE. For the two
# Chapman-Richards species (WP=1, RC=5) MODE=0 forces EFFAGE=0 (HHT = curve(AGE) − curve(0)); every other
# species is linear (coef·AGE, H unused) so ec_smhtgf already gives the MODE=0 value. SI = the species' SITEAR.
@inline function ec_essubh_hht(sp::Int, si::Float32, age::Float32)::Float32
    age <= 0f0 && return 0f0                            # SMHTGF: DTIME≤0 → no growth (ec/smhtgf.f:97)
    if sp == 1
        c1=0.375045f0; c2=0.92503f0; c3=-0.020796f0; c4=2.48811f0
        return (si/c1)*fpow(1f0-c2*fexp(c3*age), c4) - (si/c1)*fpow(1f0-c2, c4)     # EFFAGE=0
    elseif sp == 5
        c1=0.752842f0; c2=1.0f0; c3=-0.0174f0; c4=1.4711f0
        return (si/c1)*fpow(1f0-c2*fexp(c3*age), c4) - (si/c1)*fpow(1f0-c2, c4)     # EFFAGE=0
    else
        return ec_smhtgf(sp, 0f0, age, si)             # linear species: H unused ⇒ = coef·AGE
    end
end

# ec/regent.f:380-458 — DK (diameter at the grown height) / DKK (at the start height) by species form, then the
# inventory HTDBH override when height-diameter calibration is off for the species or did not happen
# (`.NOT.LHTDRG .OR. IABFLG==1`). jl used HTDBH for every species, so a DF whose Wykoff intercept AA was
# calibrated (IABFLG=0) got the wrong seedling DBH.
function _ec_regent_dk_dkk(s::StandState, ifor::Int, sp::Int, d::Float32, h::Float32, hk::Float32)
    sd = s.coef.species; c = s.calib
    local dk::Float32, dkk::Float32
    if sp == 11 || (13 <= sp <= 15) || (17 <= sp <= 30) || sp == 32      # WC logic
        if sp == 11
            dkk = -0.674f0 + 1.522f0*flog(h); dk = -0.674f0 + 1.522f0*flog(hk)
        elseif (13 <= sp <= 15) || sp == 17
            dkk = -2.089f0 + 1.980f0*flog(h); dk = -2.089f0 + 1.980f0*flog(hk)
        elseif sp == 18 || sp == 19
            dkk = -0.532f0 + 1.531f0*flog(h); dk = -0.532f0 + 1.531f0*flog(hk)
        else
            dkk = 3.102f0 + 0.021f0*h; dk = 3.102f0 + 0.021f0*hk
        end
        dkk < 0f0 && (dk = d)                                              # regent.f: IF(DKK.LT.0.0) DK=D (sic)
    elseif sp == 10                                                        # PP (EC logic)
        dk = (hk - 8.31485f0 + 0.59200f0*7f0) / 3.03659f0
        dkk = (h - 8.31485f0 + 0.59200f0*7f0) / 3.03659f0
        h < 4.5f0 && (dkk = d)
    else                                                                   # 1-9,12,16,31: inverse Wykoff
        bx = sd[:ht2][sp]
        ax = c.ht_dbh_iabflg[sp] == 1 ? sd[:ht1][sp] : c.ht_dbh_aa[sp]
        dk = (bx / (flog(hk - 4.5f0) - ax)) - 1f0
        dkk = h <= 4.5f0 ? d : (bx / (flog(h - 4.5f0) - ax)) - 1f0
    end
    if !s.control.ht_drag_sp[sp] || c.ht_dbh_iabflg[sp] == 1
        dk = ec_htdbh_dbh(ifor, sp, hk)
        dkk = h <= 4.5f0 ? d : ec_htdbh_dbh(ifor, sp, h)
    end
    return (dk, dkk)
end

function small_tree_growth!(s::StandState, stash, ::EastCascades; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    cw = clim_wk4(s, Float32(current_cycle_year(s)) + fint / 2f0)   # CLGMULT WK4 (ec/regent.f:339); nothing ⇒ 1
    sd = s.coef.species
    avh = p.avg_height; dgsd = s.control.dg_sd
    relden = p.relative_density
    ifor = ec_htdbh_ifor(p)
    isisp = Int(p.site_species); (isisp < 1 || isisp > 32) && (isisp = 10)
    scale = fint / EC_RG_REGYR; scale2 = EC_RG_REGYR / fint
    cur_year = current_cycle_year(s)
    # density modifier PCTRED (ec/regent.f): X = AVH·(RELDEN/100), capped 300; 5th-order polynomial.
    xden = avh * (relden / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = Float32(EC_RG_AB[1] + xden*(EC_RG_AB[2] + xden*(EC_RG_AB[3] + xden*(EC_RG_AB[4] + xden*(EC_RG_AB[5] + xden*EC_RG_AB[6])))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # ec/regent.f:251-272 is SPECIES-MAJOR (DO 30 ISPC … I=IND1(I3)); the per-tree ZZRAN draw must follow it.
    @inbounds for i in species_major_order(s)
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= EC_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        h = t.height[i]
        # ec/smhtgf.f:98 reads SI=SITEAR(ISPC) — the tree's OWN site index, UNCLAMPED. ec/regent.f's clamped
        # SI (=SITEAR(ISISP) bounded to the site species' [SLO+0.5, SHI]) is a dead local. jl clamped the tree's SI
        # to the site species' range, so on Mt Hood (site species GF, SHI 110) ES's SI 148 became 110 and its
        # seedlings grew 7.9 ft instead of live's 17.4 ft in one cycle.
        si = p.sp_site_index[sp]
        con = fexp(c.htg_cor_small[sp])                  # RHCON=1 default ⇒ CON = exp(HCOR)
        x = Float32(t.crown_pct[i]) / 100f0
        vigor = 150f0 * fpow(x, 3f0) * fexp(-6f0*x) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
        pothtg = ec_smhtgf(sp, h, 10f0, si)
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        wcform = sp in EC_HT_WCFORM
        xmn = EC_RG_XMIN[sp]; xmx = EC_RG_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        cap = s.control.sp_size_cap[sp, 4]
        bark = wc_bratio(sd, sp, d)
        large0 = t.ht_growth[i]
        ltrip = stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
        # ec/regent.f label 2 … `L=L+1; K=ITRN+2*I-2+L; GO TO 2`: the central record and each tripled copy get
        # their OWN ZZRAN, HTG blend (toward that record's large-tree HTG) and, below 3", their own DBH increment;
        # at or above 3" a copy keeps its own DGDRIV DG. jl ran one pass and copied the central onto both copies.
        for l in 0:(ltrip ? 2 : 0)
            htgr = pothtg * pctred * vigor * con
            zzran = 0f0
            if dgsd >= 1f0
                while true
                    zzran = bachlo(s.rng, 0f0, 1f0)
                    (zzran <= 0.5f0 && zzran >= -2f0) && break
                end
            end
            if wcform                                    # regent.f:336-340 ·WK4(I) = the CLGMULT climate multiplier
                htgr = (htgr + zzran*0.1f0) * xrhgro * scale * (cw === nothing ? 1f0 : cw[i]); htgr < 0.1f0 && (htgr = 0.1f0)
            else
                htgr = (htgr + zzran*0.1f0) * xrhgro * scale; htgr < 0f0 && (htgr = 0f0)
            end
            large = l == 0 ? large0 : (stash.htg_copy[i] ? (l == 1 ? stash.htgU[i] : stash.htgL[i]) : large0)
            htg = htgr*(1f0 - xwt) + xwt*large
            htg < 0.1f0 && (htg = 0.1f0)
            (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))
            dg = 0f0; dbh_set = -1f0
            if d < 3f0
                hk = h + htg
                if hk <= 4.5f0
                    dbh_set = d + 0.001f0*hk
                else
                    dk, dkk = _ec_regent_dk_dkk(s, ifor, sp, d, h, hk)
                    if dk < 0f0 || dkk < 0f0
                        dg = htg * 0.2f0 * bark * xrdgro
                        dk = d + dg
                    else
                        dg = (dk - dkk) * bark * xrdgro
                    end
                    # regent.f:489-493 hardwoods (20:30,32): the linear height-DBH form degenerates to 0.021·HG, so
                    # with a calibrated Wykoff fit (LHTDRG and IABFLG==0) FVS substitutes the rule of thumb 0.1·HTG.
                    ((20 <= sp <= 30 || sp == 32) && s.control.ht_drag_sp[sp] && c.ht_dbh_iabflg[sp] == 0) &&
                        (dg = 0.1f0 * htg * xrdgro)
                    dg < 0f0 && (dg = 0f0)
                    dgmx = EC_RG_DGMAX[sp] * scale
                    dg > dgmx && (dg = dgmx)
                    dds = dg*(2f0*bark*d + dg)*scale2
                    dg = sqrt(fpow(d*bark, 2f0) + dds) - bark*d
                    (d + dg) < EC_RG_DIAM[sp] && (dg = EC_RG_DIAM[sp] - d)
                end
                dg = wc_dgbnd(sp, dbh_set >= 0f0 ? dbh_set : d, dg, s.control.sp_size_cap[sp, 1], s.control.sp_size_cap[sp, 3])
            end
            if l == 0
                t.ht_growth[i] = htg
                if d < 3f0
                    t.diam_growth[i] = dg
                    dbh_set >= 0f0 && (t.dbh[i] = dbh_set)
                end
            else
                stash.is_small[i] = true
                if l == 1
                    stash.htgU[i] = htg
                    d < 3f0 && (stash.dgU[i] = dg; stash.dbhU[i] = dbh_set >= 0f0 ? dbh_set : d)
                else
                    stash.htgL[i] = htg
                    d < 3f0 && (stash.dgL[i] = dg; stash.dbhL[i] = dbh_set >= 0f0 ? dbh_set : d)
                end
            end
        end
    end
    return s
end


# ec/blkdat.f DATA HHTMAX — ESGENT's cap on a newly established tree's height.
const EC_HHTMAX = Float32[23,27,21,21,22,20,24,18,18,17, 20,22,20,20,20,20,20,20,20,20,
                          20,50,20,20,20,20,20,20,20,20, 22,20]
# ec/regent.f:385-397 DAT45 (the WC-logic diameter at 4.5 ft; 0 for the EC-logic species).
@inline function _ec_dat45(sp::Int)::Float32
    sp == 11 && return -0.674f0 + 1.522f0*flog(4.5f0)
    ((13 <= sp <= 15) || sp == 17) && return -2.089f0 + 1.980f0*flog(4.5f0)
    (sp == 18 || sp == 19) && return -0.532f0 + 1.531f0*flog(4.5f0)
    ((20 <= sp <= 30) || sp == 32) && return 3.102f0 + 0.021f0*4.5f0
    return 0f0
end

"""
    ec_esgent!(s, nstart; fint, atavh, atrelden, relden_pre, avh_pre)

ec/esgent.f → REGENT(.TRUE.,ITRNIN): grow the records established this cycle (nstart+1:n) over the rest of the
cycle, FNT = FINT−5 (LSKIPH when FINT≤5). Per new record in SPESRT order: the crown draw
(CR = 0.89722 − 0.0000461·PCCF + 0.07985·RAN), then SMHTGF·PCTRED·VIGOR·CON with its ZZRAN draw (XWT=0), then
the ESTAB diameter DBH=DK (or DK−DAT45+DIAM for a calibrated WC-logic fit), floored at DIAM, +0.001·HK, DG=DBH.
PCTRED reads the mid-period CCF/top height (5/FINT·current + (FINT−5)/FINT·start-of-cycle). ESGENT then adds
HTG to HT and caps it at HHTMAX. EC had no birth-cycle growth, so planted stands lagged a cycle (ect01 PLANT
stand: 2002 BA 0 vs live 22).
"""
function ec_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atavh::Float32 = -1.0f0, atrelden::Float32 = -1.0f0,
                    relden_pre::Float32 = -1.0f0, avh_pre::Float32 = -1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    nstart >= t.n && return s
    sd = s.coef.species
    relden = relden_pre >= 0f0 ? relden_pre : p.relative_density
    avh = avh_pre >= 0f0 ? avh_pre : p.avg_height
    dgsd = s.control.dg_sd
    lskiph = fint <= 5f0
    fnt = lskiph ? fint : fint - 5f0
    scale = fnt / EC_RG_REGYR
    cur_year = current_cycle_year(s)
    ifor = ec_htdbh_ifor(p)
    ccf = relden; avht = avh
    if fnt > 0f0 && atrelden >= 0f0 && atavh >= 0f0      # regent.f:232-235
        ccf = (5f0/fint)*relden + ((fint - 5f0)/fint)*atrelden
        avht = (5f0/fint)*avh + ((fint - 5f0)/fint)*atavh
    end
    xden = avht * (ccf / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = Float32(EC_RG_AB[1] + xden*(EC_RG_AB[2] + xden*(EC_RG_AB[3] + xden*(EC_RG_AB[4] + xden*(EC_RG_AB[5] + xden*EC_RG_AB[6])))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    newidx = sort(collect((nstart+1):t.n); by = i -> (Int(t.species[i]), i))   # SPESRT species-then-record
    @inbounds for i in newidx
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= EC_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        ran = 0f0
        while true
            ran = bachlo(s.rng, 0f0, 1f0); (-1f0 <= ran <= 1f0) && break
        end
        pccf = s.density.point_ccf[Int(t.plot_id[i])]
        cr0 = 0.89722f0 - 0.0000461f0 * pccf + 0.07985f0 * ran
        cr0 > 0.90f0 && (cr0 = 0.90f0); cr0 < 0.20f0 && (cr0 = 0.20f0)
        icr0 = trunc(Int32, cr0 * 100f0 + 0.5f0)
        t.crown_pct[i] = icr0; t.crown_ratio[i] = Float32(icr0)
        h = t.height[i]
        local htg::Float32
        if lskiph
            htg = 0f0
        else
            con = fexp(c.htg_cor_small[sp])
            x = Float32(icr0) / 100f0
            vigor = 150f0 * fpow(x, 3f0) * fexp(-6f0*x) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
            htgr = ec_smhtgf(sp, h, 10f0, p.sp_site_index[sp]) * pctred * vigor * con
            zzran = 0f0
            if dgsd >= 1f0
                while true
                    zzran = bachlo(s.rng, 0f0, 1f0); (zzran <= 0.5f0 && zzran >= -2f0) && break
                end
            end
            xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
            htgr = (htgr + zzran*0.1f0) * xrhgro * scale
            if sp in EC_HT_WCFORM
                htgr < 0.1f0 && (htgr = 0.1f0)
            else
                htgr < 0f0 && (htgr = 0f0)
            end
            htg = htgr                                   # XWT = 0 under LESTB
            htg < 0.1f0 && (htg = 0.1f0)
            cap = s.control.sp_size_cap[sp, 4]
            (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))
        end
        if d < 3f0
            hk = h + htg
            if hk <= 4.5f0
                t.dbh[i] = d + 0.001f0*hk; t.diam_growth[i] = 0f0
            else
                dk, _ = _ec_regent_dk_dkk(s, ifor, sp, d, h, hk)
                dat45 = (!s.control.ht_drag_sp[sp] || c.ht_dbh_iabflg[sp] == 1) ? 0f0 : _ec_dat45(sp)
                dbh = (dat45 > 0f0 && hk >= 4.5f0 && s.control.ht_drag_sp[sp] && c.ht_dbh_iabflg[sp] == 0) ?
                      dk - dat45 + EC_RG_DIAM[sp] : dk
                (dbh < EC_RG_DIAM[sp] || hk < 4.5f0) && (dbh = EC_RG_DIAM[sp])
                dbh = dbh + 0.001f0*hk
                t.dbh[i] = dbh; t.diam_growth[i] = dbh
            end
        end
        # ec/esgent.f:51-65: HTG·WK4 (WK4 = HTIMLT, ec/estab.f:510-516,652 — 5/5.0001 for a start-of-cycle PLANT after
        # ec/essubh.f:52-58 rewrote DELAY/TRAGE), the WK4<1 DBH/DG rescale (DG from the NEW DBH, as written), HHTMAX cap.
        # MEASURED: every EC plant_cal/plant_cyc case had a QMD 1e-6-rel high (the WK4=1 increment).
        esgent_finish!(t, i, htg, EC_HHTMAX[sp])
    end
    return s
end

"""ec/regent.f LSTART small-tree HEIGHT calibration (label 40; called from ec/cratet.f:726). PCTRED from
X=AVH·RELDEN/100 (AVH = the cratet.f:684 AVHT40, RELDEN = the cratet.f:231 backdating DENSE's). Per species (IND1
order, LHTCAL default .TRUE.): current DBH<5, backdated H=HT−HTG (IHTG<2) ≥0.01, measured HTG≥0.001 ⇒
EDH = SMHTGF(ISPC, H, REGYR)·PCTRED·VIGOR·RHCON floored at 0.1 (SMHTGF reads SITEAR(ISPC); the SLO/SHI-clamped
SITEAR(ISISP) in regent.f is a dead local, as in the growth path), TERM = HTG·REGYR/FINTH; CORNEW = mean TERM /
mean EDH when N≥NCALHT(5), ≤0 ⇒ 1E−4, trapped to [0.0821, 12.1825] else 1. HCOR = ln(CORNEW) → htg_cor_init."""
function ec_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector, avh::Float32)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    s.control.growth_ifinth == 0 && return s                  # regent.f IF(IFINTH.EQ.0) GOTO 95
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : 5f0
    scale3 = EC_RG_REGYR / finth
    relden = c.cratet_relden
    x = avh * (relden / 100f0); x > 300f0 && (x = 300f0)
    pctred = Float32(EC_RG_AB[1] + x*(EC_RG_AB[2] + x*(EC_RG_AB[3] + x*(EC_RG_AB[4] + x*(EC_RG_AB[5] + x*EC_RG_AB[6])))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for sp in 1:size(isct, 1)
        i1 = Int(isct[sp, 1]); i1 == 0 && continue
        i2 = Int(isct[sp, 2])
        rhcon = (s.control.regh_cor2_on && s.control.regh_cor2[sp] > 0f0) ? s.control.regh_cor2[sp] : 1f0
        snp = 0f0; snx = 0f0; sny = 0f0; nh = 0
        for k in i1:i2
            i = Int(ind1[k])
            h = t.height[i]
            s.control.growth_ihtg < 2 && (h = h - t.ht_growth[i])
            (saved_dbh[i] >= 5f0 || h < 0.01f0) && continue
            xv = Float32(t.crown_pct[i]) / 100f0
            vigor = 150f0 * fpow(xv, 3f0) * fexp(-6f0*xv) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
            edh = ec_smhtgf(sp, h, EC_RG_REGYR, p.sp_site_index[sp]) * pctred * vigor * rhcon
            edh < 0.1f0 && (edh = 0.1f0)
            hg = t.ht_growth[i]; hg < 0.001f0 && continue
            pr = t.tpa[i]
            snp += pr; snx += edh * pr; sny += hg * scale3 * pr; nh += 1
        end
        nh < 5 && continue
        snx /= snp; sny /= snp
        cornew = sny / snx
        cornew <= 0f0 && (cornew = 1f-4)
        (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
        c.htg_cor_init[sp] = flog(cornew)
    end
    return s
end
