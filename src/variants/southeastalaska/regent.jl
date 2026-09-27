# =============================================================================
# regent.jl (southeastalaska) — AK small-tree height/diameter growth (ak/regent.f + ak/htcalc.f),
# the AK ESTAB birth-cycle grower (ak/esgent.f → REGENT(.TRUE.,ITRNIN)) and the LSTART small-tree
# height calibration (REGENT(.FALSE.,1) from cratet.f:629).
#
# Height: HTCALC site curves — Payandeh (SS/WH/MH, sp 8/11/12) or Hegyi (all others) — MODE 0 inverts the
# tree height to an age, MODE 2 gives the 10-yr increment H(AGET+10)−H(AGET) (0.10 when H is within 1 ft
# of the asymptote HTMAX). HTGR = (HTGR + ZZRAN·0.1)·XRHGRO·SCALE1·CON, ZZRAN~N(0,1) redrawn outside
# [−2,0.5] (DGSD≥1). Growth (non-ESTAB) pass averages it with the large-tree HTG, HTGR=(HTGR+LTHG)/2, then
# blends HTG = HTGR·(1−XWT) + XWT·LTHG, XWT=(D−XMIN)/(XMAX−XMIN) (0 when D≤XMIN).
# Diameter (D<DGMAX): the Chapman-Richards HT-DBH inverse (HTT11/12/13) at H and H+HTG gives DK/DKK;
# DGSM=(DK−DKK)·BARK, rescaled in DDS space by SCALE2=YR/FNT, blended with the large-tree DG by XDWT.
# Tripling: each record's two copies (K=ITRN+2I−2+L) re-run labels 2..23 with their own ZZRAN draw and
# their own large-tree HTG/DG; H and D stay the record's own.
# =============================================================================

const AK_RG_REGYR = 10.0f0
const AK_RG_XMIN  = fill(1.0f0, 23)
const AK_RG_XMAX  = Float32[5,5,8,3,3, 3,3,5,8,8, 7,8,3,3,3, 3,3,3,3,3, 3,3,3]
const AK_RG_DGMAX = Float32[5,5,5,3,3, 3,3,5,5,5, 5,5,3,3,3, 3,3,3,3,3, 3,3,3]
const AK_RG_DIAM  = Float32[0.3,0.3,0.3,0.2,0.4, 0.4,0.2,0.4,0.4,0.3, 0.3,0.3,0.4,0.2,0.2,
                            0.1,0.1,0.2,0.3,0.2, 0.1,0.1,0.2]
# ak/blkdat.f HT1/HT2 (Wykoff HT-DBH) — only read by the LHTDRG (HTDREG keyword) branch of REGENT.
const AK_HT1 = Float32[5.047089,5.047089,4.683932,4.350320,4.633182,4.633182,4.350320,5.047089,4.513771,4.716522,
                       5.007972,4.701229,4.633182,4.622461,4.542755,4.393426,4.393426,4.622461,4.479633,4.622461,
                       4.622461,4.622461,4.622461]
const AK_HT2 = Float32[-12.629014,-12.629014,-10.690737,-5.776563,-6.81926,-6.81926,-5.776563,-12.629014,-10.853785,
                       -11.426736,-12.085418,-12.133655,-6.81926,-6.696442,-6.654068,-3.968868,-3.968868,-6.696442,
                       -5.03023,-6.696442,-6.696442,-6.696442,-6.696442]
# ak/blkdat.f ESCOMN XMIN / HHTMAX — establishment minimum and maximum (ESGENT cap) seedling heights.
const AK_ES_XMIN   = Float32[1.0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.0,0.5, 0.5,0.5,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                             1.0,1.0,1.0]
const AK_ES_HHTMAX = Float32[20.0,20.0,12.0,13.0,13.0,13.0,13.0,20.0,9.0,15.0, 19.0,10.0,13.0,7.6,30.0,20.0,20.0,
                             16.9,18.0,16.9,7.6,11.0,7.6]

# ak/htcalc.f DATA HEGYI1..3 (Hegyi site curve) + the Payandeh constants (SS/WH/MH).
const AK_HEGYI1 = Float32[1.1945,1.3832,1.1243,1.1637,1.2883,1.2883,1.1637,1.0458,1.0236,1.1243,1.1514,1.1514,
                          1.2883,1.1318,1.0142,1.1580,1.1580,1.1318,1.2025,1.1318,1.1318,1.1318,1.1318]
const AK_HEGYI2 = Float32[-0.0236,-0.0155,-0.0263,-0.0215,-0.0181,-0.0181,-0.0215,-0.0380,-0.0465,-0.0263,-0.0237,
                          -0.0237,-0.0181,-0.0226,-0.0421,-0.0175,-0.0175,-0.0226,-0.0158,-0.0226,-0.0226,-0.0226,-0.0226]
const AK_HEGYI3 = Float32[1.7918,1.3597,1.5662,1.2243,1.4177,1.4177,1.2243,1.9804,2.4269,1.5662,1.4365,1.4365,
                          1.4177,1.1233,0.9422,0.7687,0.7687,1.1233,0.7994,1.1233,1.1233,1.1233,1.1233]
const _AK_PY_B1 = 1.5469f0; const _AK_PY_B2 = 1.0018f0; const _AK_PY_B3 = -0.0114f0
const _AK_PY_B4 = 1.0883f0; const _AK_PY_B5 = 0.0072f0
@inline _ak_payandeh(sp::Integer) = sp == 8 || sp == 11 || sp == 12

"ak/htcalc.f HTMAX (computed on every call): Payandeh B1·SI^B2, else Hegyi C1·SI."
@inline ak_htcalc_htmax(sp::Integer, xsite::Float32)::Float32 =
    _ak_payandeh(sp) ? _AK_PY_B1 * fpow(xsite, _AK_PY_B2) : AK_HEGYI1[sp] * xsite

"ak/htcalc.f MODE 0 — tree age from height (caller has checked HTMAX−H>1)."
@inline function ak_htcalc_age(sp::Integer, xsite::Float32, h::Float32)::Float32
    if _ak_payandeh(sp)
        return 1.0f0 / _AK_PY_B3 * flog(1f0 - fpow(h / _AK_PY_B1 / fpow(xsite, _AK_PY_B2),
                                               1.0f0 / _AK_PY_B4 / fpow(xsite, _AK_PY_B5)))
    else
        c1 = AK_HEGYI1[sp]; c2 = AK_HEGYI2[sp]; c3 = AK_HEGYI3[sp]
        return 1f0 / c2 * flog(1f0 - fpow(h / c1 / xsite, 1f0 / c3))
    end
end

"ak/htcalc.f MODE 1 — height at age AGET."
@inline function ak_htcalc_height(sp::Integer, xsite::Float32, aget::Float32)::Float32
    if _ak_payandeh(sp)
        return _AK_PY_B1 * fpow(xsite, _AK_PY_B2) * fpow(1f0 - fexp(_AK_PY_B3 * aget), _AK_PY_B4 * fpow(xsite, _AK_PY_B5))
    else
        return AK_HEGYI1[sp] * xsite * fpow(1f0 - fexp(AK_HEGYI2[sp] * aget), AK_HEGYI3[sp])
    end
end

"ak/htcalc.f MODE 2 — 10-yr height increment H(AGET+10) − H(AGET)."
@inline ak_htcalc_htg(sp::Integer, xsite::Float32, aget::Float32)::Float32 =
    ak_htcalc_height(sp, xsite, aget + 10.0f0) - ak_htcalc_height(sp, xsite, aget)

"REGENT labels 2..4 height potential: HTCALC MODE 0 then MODE 2 (0.10 within 1 ft of HTMAX)."
@inline function ak_regent_htgr(sp::Integer, xsite::Float32, h::Float32)::Float32
    htmax = ak_htcalc_htmax(sp, xsite)
    htmax - h <= 1f0 && return 0.10f0
    return ak_htcalc_htg(sp, xsite, ak_htcalc_age(sp, xsite, h))
end

"Chapman-Richards HT-DBH inverse (regent.f DK): 1/HTT12·LOG(1−((H−4.5)/HTT11)^(1/HTT13))."
@inline ak_regent_dinv(sp::Integer, h::Float32)::Float32 =
    1f0 / AK_HTT12[sp] * flog(1f0 - fpow((h - 4.5f0) / AK_HTT11[sp], 1f0 / AK_HTT13[sp]))

"ak/dgbnd.f — SIZCAP diameter cap only."
@inline function ak_dgbnd(dbh::Float32, ddg::Float32, sizcap1::Float32, sizcap3::Float32)::Float32
    if (dbh + ddg) > sizcap1 && sizcap3 < 1.5f0
        ddg = sizcap1 - dbh
        ddg < 0.01f0 && (ddg = 0.01f0)
    end
    return ddg
end

# regent.f keeps DGSM in a SAVEd local (gfortran -fno-automatic); the DK<0 branch reads the stale value.
const _AK_RG_DGSM = Ref(0f0)

"""
One REGENT pass for slot K (labels 2..23 of ak/regent.f): returns (htg, dbh_direct (NaN ⇒ unchanged), dg, dg_set).
`lthg`/`dglt` are the slot's large-tree HTG/DG, `dbhk` the slot's DBH(K) (=D for the record and its copies).
"""
function _ak_regent_slot(s::StandState, sp::Int, h::Float32, d::Float32, bark::Float32, htgr0::Float32,
                         lthg::Float32, dglt::Float32, dbhk::Float32, con::Float32, scale1::Float32,
                         scale2::Float32, lestb::Bool, lskiph::Bool, xrhgro::Float32, xrdgro::Float32)
    local htg::Float32
    if lskiph
        htg = 0f0
    else
        zzran = 0f0
        if s.control.dg_sd >= 1f0
            while true
                zzran = bachlo(s.rng, 0f0, 1f0)
                (zzran > 0.5f0 || zzran < -2.0f0) || break
            end
        end
        htgr = (htgr0 + zzran * 0.1f0) * xrhgro * scale1 * con
        xmn = AK_RG_XMIN[sp]
        xwt = (d <= xmn || lestb) ? 0f0 : (d - xmn) / (AK_RG_XMAX[sp] - xmn)
        lestb || (htgr = (htgr + lthg) / 2.0f0)
        htg = htgr * (1.0f0 - xwt) + xwt * lthg
        htg < 0.1f0 && (htg = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        if (h + htg) > cap
            htg = cap - h; htg < 0.1f0 && (htg = 0.1f0)
        end
    end
    d >= AK_RG_DGMAX[sp] && return (htg, NaN32, dglt, false)            # GO TO 23
    xmn = AK_RG_XMIN[sp]
    xdwt = (d <= xmn || lestb) ? 0f0 : (d - xmn) / (AK_RG_DGMAX[sp] - xmn)
    hk = h + htg
    local dg::Float32
    dbh_set = NaN32
    if hk <= 4.5f0
        dg = 0f0
        dbh_set = d + 0.001f0 * hk
        dbhk = dbh_set
    else
        local dk::Float32, dkk::Float32
        if !s.control.ht_drag_sp[sp] || s.calib.ht_dbh_iabflg[sp] == 1
            dk = ak_regent_dinv(sp, hk)
            dkk = h <= 4.5f0 ? d : ak_regent_dinv(sp, h)
        else
            bx = AK_HT2[sp]
            ax = s.calib.ht_dbh_iabflg[sp] == 1 ? AK_HT1[sp] : s.calib.ht_dbh_aa[sp]
            dk = bx / (flog(hk - 4.5f0) - ax) - 1.0f0
            dkk = h <= 4.5f0 ? d : bx / (flog(h - 4.5f0) - ax) - 1.0f0
        end
        if lestb
            dbh_set = dk < 0.1f0 ? 0.1f0 : dk
            dbh_set = dbh_set + 0.001f0 * hk
            dbhk = dbh_set
            dg = dbh_set
        else
            dgk = dglt
            if dk < 0f0 || dkk < 0f0
                dgk = htg * 0.2f0 * bark * xrdgro
            else
                _AK_RG_DGSM[] = (dk - dkk) * bark * xrdgro
            end
            dgsm = _AK_RG_DGSM[]
            dgsm < 0f0 && (dgsm = 0f0)
            dds = (dgsm * (2.0f0 * bark * d + dgsm)) * scale2
            dgsm = sqrt((d * bark)^2 + dds) - bark * d
            _AK_RG_DGSM[] = dgsm
            dg = dgsm * (1.0f0 - xdwt) + xdwt * dgk
        end
        (dbhk + dg) < AK_RG_DIAM[sp] && (dg = AK_RG_DIAM[sp] - dbhk)
    end
    dg = ak_dgbnd(dbhk, dg, s.control.sp_size_cap[sp, 1], s.control.sp_size_cap[sp, 3])
    return (htg, dbh_set, dg, true)
end

@inline function _ak_rg_con(s::StandState, sp::Int)::Float32
    rhcon = (s.control.regh_cor2_on && s.control.regh_cor2[sp] > 0f0) ? s.control.regh_cor2[sp] : 1f0
    return rhcon * fexp(s.calib.htg_cor_small[sp])           # regent.f CON = RHCON·EXP(HCOR)
end

"ak/regent.f REGENT(.FALSE.,1) — the cycling small-tree growth (grincr.f:449, after DGDRIV/HTGF)."
function small_tree_growth!(s::StandState, stash, ::SoutheastAlaska; fint::Float32 = 10.0f0)
    p, t = s.plot, s.trees
    t.n == 0 && return s
    scale1 = fint / AK_RG_REGYR                        # SCALE1 = FNT/REGYR
    scale2 = s.control.year / fint                     # SCALE2 = YR/FNT
    trip = stash !== nothing && !isempty(stash.dgU)
    yr_now = current_cycle_year(s)
    @inbounds for i in species_major_order(s)          # DO 30 ISPC / DO 25 I3 → I=IND1(I3)
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= AK_RG_XMAX[sp] && continue
        h = t.height[i]
        bark = ak_bratio(sp, d)
        xsite = p.sp_site_index[sp]
        xrhgro = active_multiplier(s.control, :regh, sp, yr_now)
        xrdgro = active_multiplier(s.control, :regd, sp, yr_now)
        con = _ak_rg_con(s, sp)
        htgr0 = ak_regent_htgr(sp, xsite, h)
        nrec = trip ? 3 : 1
        # HTGF gave each tripled copy HTG(ITFN)=TEMHTG, the record's own large-tree HTG (htgf.f:419-434) — read it
        # before the L=0 pass overwrites HTG(I) with the regent value.
        lthg0 = t.ht_growth[i]
        for l in 0:(nrec - 1)
            lthg = l == 0 ? lthg0 :
                   (stash.htg_copy[i] ? (l == 1 ? stash.htgU[i] : stash.htgL[i]) : lthg0)
            dglt = l == 0 ? t.diam_growth[i] : (l == 1 ? stash.dgU[i] : stash.dgL[i])
            htg, dbhs, dg, dgset = _ak_regent_slot(s, sp, h, d, bark, htgr0, lthg, dglt, d, con, scale1, scale2,
                                                   false, false, xrhgro, xrdgro)
            if l == 0
                t.ht_growth[i] = htg
                if dgset
                    isnan(dbhs) || (t.dbh[i] = dbhs)
                    t.diam_growth[i] = dg
                end
            else
                if dgset
                    if l == 1
                        stash.dbhU[i] = isnan(dbhs) ? d : dbhs; stash.dgU[i] = dg
                    else
                        stash.dbhL[i] = isnan(dbhs) ? d : dbhs; stash.dgL[i] = dg
                    end
                end
                l == 1 ? (stash.htgU[i] = htg) : (stash.htgL[i] = htg)
                stash.is_small[i] = true
            end
        end
    end
    return s
end

"""
    ak_esgent!(s, nstart; fint, pccf_pre)

ak/esgent.f for the records ak_establish! created this cycle (nstart+1:n): SPESRT, then REGENT(.TRUE.,ITRNIN) in
IND1 species order — the open-grown crown draw CR=0.89722−0.0000461·PCCF(IP)+0.07985·RANC (RANC~N(0,1) in [−1,1]),
FNT=FINT−5 (LSKIPH when FINT≤5), no large-tree blend, DBH from the HT-DBH inverse at H+HTG — then HT += HTG·WK4,
the WK4<1 DBH rescale and the HHTMAX cap. Returns ΣPROB·HT / ΣPROB by species over IMC=1 records (SUMPX/SUMPI).
"""
function ak_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0, pccf_pre::Vector{Float32} = Float32[],
                    sumpx::Vector{Float32} = zeros(Float32, 23), sumpi::Vector{Float32} = zeros(Float32, 23))
    p, t = s.plot, s.trees
    nstart >= t.n && return s
    lskiph = fint <= 5f0
    fnt = lskiph ? fint : fint - 5f0
    scale1 = fnt / AK_RG_REGYR
    scale2 = s.control.year / fnt
    yr_now = current_cycle_year(s)
    spesrt_reorder!(t)                                   # esgent.f:52 SPESRT
    newidx = filter(>(nstart), species_major_order(s))
    @inbounds for i in newidx
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= AK_RG_XMAX[sp] && continue
        h = t.height[i]
        bark = ak_bratio(sp, d)
        ip = Int(t.plot_id[i])
        pccf = (1 <= ip <= length(pccf_pre)) ? pccf_pre[ip] :
               (1 <= ip <= length(s.density.point_ccf) ? s.density.point_ccf[ip] : 0f0)
        cr = 0.89722f0 - 0.0000461f0 * pccf
        ranc = 0f0
        while true
            ranc = bachlo(s.rng, 0f0, 1f0)
            (ranc < -1.0f0 || ranc > 1.0f0) || break
        end
        cr = cr + 0.07985f0 * ranc
        cr > 0.90f0 && (cr = 0.90f0); cr < 0.20f0 && (cr = 0.20f0)
        icr = unsafe_trunc(Int32, cr * 100.0f0 + 0.5f0)
        t.crown_pct[i] = icr; t.crown_ratio[i] = Float32(icr)
        xsite = p.sp_site_index[sp]
        xrhgro = active_multiplier(s.control, :regh, sp, yr_now)
        xrdgro = active_multiplier(s.control, :regd, sp, yr_now)
        con = _ak_rg_con(s, sp)
        htgr0 = lskiph ? 0f0 : ak_regent_htgr(sp, xsite, h)
        htg, dbhs, dg, dgset = _ak_regent_slot(s, sp, h, d, bark, htgr0, 0f0, t.diam_growth[i], d, con, scale1,
                                               scale2, true, lskiph, xrhgro, xrdgro)
        if dgset
            isnan(dbhs) || (t.dbh[i] = dbhs)
            t.diam_growth[i] = dg
        end
        # esgent.f:59-74
        wk4 = t.htimlt[i]
        htemp = h + htg
        htg = htg * wk4
        hn = h + htg
        if wk4 < 1.0f0
            if hn < 4.5f0
                t.dbh[i] = 0.1f0 + 0.001f0 * hn; t.diam_growth[i] = 0f0
            else
                t.dbh[i] = t.dbh[i] * (hn / htemp); t.diam_growth[i] = t.diam_growth[i] * (hn / htemp)
            end
        end
        hn > AK_ES_HHTMAX[sp] && (hn = AK_ES_HHTMAX[sp])
        t.ht_growth[i] = htg
        t.height[i] = hn
        if t.mort_code[i] == 1
            sumpx[sp] += t.tpa[i] * hn; sumpi[sp] += t.tpa[i]
        end
    end
    return s
end

"""
ak/regent.f LSTART (cratet.f:629 REGENT(.FALSE.,1)): the small-tree height calibration. Per LHTCAL species, over
records with DBH<3 and a measured HTG (≥0.001): EDH = HTCALC-predicted 10-yr increment on the start-of-period
height (H−HTG when IHTG<2), TERM = HTG·SCALE3 (SCALE3=REGYR/FINTH); CORNEW = ΣTERM·P/ΣEDH·P with ≥NCALHT obs,
trapped to [0.0821,12.1825] (else 1). HCOR_init = ln(CORNEW); dgdriv attenuates it per cycle (DIFH).
"""
function ak_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector, saved_dbh::AbstractVector)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    s.control.growth_ifinth == 0 && return s                       # IF(IFINTH.EQ.0) GOTO 95
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : 10f0
    scale3 = AK_RG_REGYR / finth
    ncalht = 5                                                      # NCALHT (grinit.f; CALBSTAT keyword unported)
    @inbounds for sp in 1:23
        i1 = isct[sp, 1]; i1 == 0 && continue            # LHTCAL all .TRUE. (grinit.f:110)
        i2 = isct[sp, 2]
        xsite = p.sp_site_index[sp]
        snp = 0f0; snx = 0f0; sny = 0f0; n = 0
        for k in i1:i2
            i = ind1[k]
            h = t.height[i]
            s.control.growth_ihtg < 2 && (h = h - t.ht_growth[i])
            (saved_dbh[i] >= 3.0f0 || h < 0.01f0) && continue
            edh = ak_regent_htgr(sp, xsite, h)
            pr = t.tpa[i]
            t.ht_growth[i] < 0.001f0 && continue
            term = t.ht_growth[i] * scale3
            snp += pr; snx += edh * pr; sny += term * pr; n += 1
        end
        n < ncalht && continue
        snx /= snp; sny /= snp
        cornew = sny / snx
        cornew <= 0f0 && (cornew = 1f-4)
        (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
        c.htg_cor_init[sp] = flog(cornew)
    end
    return s
end
