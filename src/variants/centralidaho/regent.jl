# =============================================================================
# regent.jl (centralidaho) — CI small-tree growth (ci/regent.f). Chunk 6. Overrides DG/HTG for
# small trees (D < XMAX). NIVAR path (CI conifers): per-subcycle height via the HTGRL regression,
# accumulated, then DBH from HTDBH; DG = new−old DBH. Coeffs from data/centralidaho/regent_coeffs_1d.csv.
#   HTGRL = CNST + HTBA·RELHT·PTBAA + PBAL·PTBALI + STBA·BA + RLHT·RELHT + CRSQ·RCR² + CR·RCR
#         + BAL·TBAL + HDM1·RHDM1 + HDM2·RHDM2 + PTBA·PTBAA   (ci/regent.f:641)
#   CON = RHCON·exp(HCOR)  (RHCON≈1, HCOR = c.htg_cor_small); H2 = H1 + HTGRL·SCALE·XRHGRO·CON
#   DK = exp(DHCN + DHHT·ln(HK) + DHCR·ln(RCR))  (HTDBH); DG = DK − D_old
# SMHTGF (stochastic BETA/FINDAG) is the establishment/aspen path — deferred; cit01 small trees are conifers.
# =============================================================================

let
    path = joinpath(CI_DATADIR, "regent_coeffs_1d.csv")
    hdr = split(strip(readlines(path)[1]), ',')
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    ci(n) = findfirst(==(n), hdr)
    col(n) = Float32[parse(Float32, strip(rows[sp][ci(n)])) for sp in 1:19]
    for n in ("CNST","CR","CRSQ","DGMAX","DHCN","DHCR","DHHT","HDM1","HDM2","HTBA","RLHT","STBA",
              "XMAX","XMIN","BAL","PBAL","PTBA")
        @eval global const $(Symbol("CI_RG_" * n)) = $(col(n))
    end
end
const CI_RG_REGYR = 5.0f0
const CI_RG_DIAM = Float32[0.4,0.3,0.3,0.3,0.2,0.2,0.4,0.3,0.3,0.5,0.3,0.3,0.2,0.3,0.2,0.3,0.2,0.2,0.2]  # ci/regent.f DATA DIAM (min DBH)
const CI_RG_AB = Float32[1.11436, -0.011493, 0.43012f-4, -0.72221f-7, 0.5607f-10, -0.1641f-13]  # PCTRED poly (regent.f:448)
@inline _ci_ut_species(sp::Int) = sp == 13 || sp == 14 || sp == 15 || sp == 17 || sp == 19  # ci/regent.f UTVAR branch

"""ci/regent.f sp17/19 (CW/OH) Curtis-Arney height→DBH inverse (P2=1709.7229,P3=5.8887,P4=−0.2286)."""
@inline function _ci_ut_htdbh(ht::Float32)::Float32
    p2 = 1709.7229f0; p3 = 5.8887f0; p4 = -0.2286f0
    hat3 = 4.5f0 + p2 * exp(-p3 * 3.0f0^p4)
    ht >= hat3 ? exp(log((log(ht - 4.5f0) - log(p2)) / (-p3)) * (1.0f0 / p4)) :
                 ((ht - 4.51f0) * 2.7f0 / (hat3 - 4.51f0)) + 0.3f0
end

# Push the regent growth into the tripling stash so the upper/lower sub-records use the REGENT growth, not the
# stale/explosive large-tree dgf DG/HTG (diameter_growth! balloons tiny DBH). FVS computes HTG(K) for EVERY tripled
# record (ci/regent.f DO 25 L-loop) ⇒ copies always get the regent HTG (is_small); the D<BKPT diameter dub (DG(K))
# is applied per-record too, so dgU/dgL carry it only when the central got dubbed (small=true, D<BKPT=3 for CIVAR).
# ZZRAN stays deferred (intentional, regent.f:934 = grow-phase RNG desync class) ⇒ copies == central deterministically.
@inline function _ci_rg_stash!(stash, t, i::Int, small::Bool)
    if stash !== nothing && !isempty(stash.htgU) && i <= length(stash.htgU)
        stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
        !isempty(stash.is_small) && (stash.is_small[i] = true)
        if small && !isempty(stash.dgU) && i <= length(stash.dgU)
            stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        end
    end
end

function small_tree_growth!(s::StandState, stash, ::CentralIdaho; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species
    t.n == 0 && return s
    n = t.n
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    kodtyp = Int(p.habitat_code)
    rhdm1 = (500 <= kodtyp < 600) ? 1.0f0 : 0.0f0
    rhdm2 = (600 <= kodtyp < 700) ? 1.0f0 : 0.0f0
    regyr = CI_RG_REGYR
    ntyr = Int(round(fint))                               # NTYR (grow cycle; LSTART/ESTAB partial-cycle deferred)
    scale_h = Float32(ntyr) / regyr                       # ci/regent.f:518 CIVAR SCALE = NTYR/REGYR (single pass)
    scale2 = 1.0f0                                        # SCALE2 = YR/NTYR = 1 for a full cycle ⇒ DDS round-trip identity
    scale_ut = 1.0f0                                      # UTVAR SCALE = NTYR/YR = 1 for a full cycle (regent.f:514)
    cur_year = current_cycle_year(s)
    slo_a = sd[:site_lo]; shi_a = sd[:site_hi]
    xd = avh * (relden / 100.0f0); xd > 300.0f0 && (xd = 300.0f0)   # PCTRED input X = AH·(R/100) (regent.f:446)
    pctred = CI_RG_AB[1] + xd*(CI_RG_AB[2] + xd*(CI_RG_AB[3] + xd*(CI_RG_AB[4] + xd*(CI_RG_AB[5] + xd*CI_RG_AB[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # CIVAR (sp 1-10,18) small-tree height+DBH — single pass (regent.f:528 J>1 skips). Verified vs the PRISTINE
    # FVSci_clean via the DEBUG keyword: seedling I=24 H=1.01→HK=5.267, DK=0.5179, DKK=0.1 ⇒ DG=0.387 (matches).
    # ZZRAN random height (regent.f:934) DEFERRED = deterministic mean HTGR=HTGR1 (grow-phase RNG desync class).
    # TT(11,12,16)=SMHTGF and UT(13,14,15,17,19)=aspen/juniper get the CIVAR form as a fallback (rare/deferred).
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d0 = t.dbh[i]
        d0 >= CI_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        (sp < 1 || sp > 19) && continue
        h0 = t.height[i]
        # ---- UTVAR path: aspen(13)/juniper(14)/MC(15)/cottonwood(17)/hardwood(19) borrowed from UT variant
        # (ci/regent.f UTVAR branch). POTHTG·PCTRED·VIGOR·CON height; per-species DBH. ZZRAN deferred (det. mean).
        if _ci_ut_species(sp)
            sitear = p.sp_site_index[sp]; sj = sitear
            con = exp(c.htg_cor_small[sp])                # RHCON=1
            if sp == 13                                   # aspen: FINDAG site-age from height (Sheppard)
                slo = slo_a[sp]; shi = shi_a[sp]
                si = sitear; si > shi && (si = shi); si <= slo && (si = slo + 0.5f0)
                relsi = (si - slo) / (shi - slo); rsimod = 0.5f0 * (1.0f0 + relsi)
                age = (h0 * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0)
                hite1 = 26.9825f0 * age^1.1752f0; hite2 = 26.9825f0 * (age + 10.0f0)^1.1752f0
                htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 0.75f0
            else                                          # 14,15,17,19: POTHTG·PCTRED·VIGOR·CON
                pothtg = (sj / 5.0f0) * (sj * 1.5f0 - h0) / (sj * 1.5f0) * 0.83f0
                xcr = Float32(t.crown_pct[i]) / 100.0f0
                vigor = 150.0f0 * xcr^3 * exp(-6.0f0 * xcr) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
                sp == 14 && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)   # juniper
                htgr = pothtg * pctred * vigor * con
            end
            htgr = htgr * scale_ut                        # ZZRAN deferred; XRHGRO=1
            htgr < 0.1f0 && (htgr = 0.1f0)
            xmn = CI_RG_XMIN[sp]; xmx = CI_RG_XMAX[sp]
            xwt = d0 <= xmn ? 0.0f0 : (d0 - xmn) / (xmx - xmn)
            htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]; htg < 0.1f0 && (htg = 0.1f0)
            hcap = s.control.sp_size_cap[sp, 4]           # SIZCAP(sp,4) max height (regent.f:955)
            (hcap > 0f0 && h0 + htg > hcap) && (htg = max(hcap - h0, 0.1f0))
            t.ht_growth[i] = htg
            hk = h0 + htg
            bark = ci_bratio(sd, sp, d0)
            if hk <= 4.5f0
                t.diam_growth[i] = 0.0f0
            else
                local dk::Float32, dkk::Float32
                if sp == 14                               # juniper linear-site
                    dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                    dkk = h0 < 4.5f0 ? d0 : (h0 - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                elseif sp == 15                           # mountain-mahogany linear
                    dk = 3.1020f0 + 0.0210f0 * hk
                    dkk = 3.1020f0 + 0.0210f0 * h0; dkk < 0.0f0 && (dkk = d0); dk < dkk && (dk = dkk + 0.01f0)
                else                                      # 17,19 cottonwood/hardwood Curtis-Arney
                    dk = _ci_ut_htdbh(hk)
                    dkk = h0 <= 4.5f0 ? d0 : _ci_ut_htdbh(h0)
                end
                dg = (dk - dkk) * bark; dg < 0.0f0 && (dg = 0.0f0)
                dgmx = CI_RG_DGMAX[sp] * scale_ut; dg > dgmx && (dg = dgmx)
                dds = dg * (2.0f0 * bark * d0 + dg) * scale2
                dg = sqrt((d0 * bark)^2 + dds) - bark * d0
                (d0 + dg) < CI_RG_DIAM[sp] && (dg = CI_RG_DIAM[sp] - d0)
                t.diam_growth[i] = dg
            end
            _ci_rg_stash!(stash, t, i, d0 < 3.0f0)   # copies use the regent growth (UTVAR); dgU/dgL only if dubbed
            continue
        end
        pt = Int(t.plot_id[i])
        ptba = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : ba
        pct = t.crown_ratio[i]
        relht = avh > 0f0 ? h0 / avh : 1.5f0; relht > 1.5f0 && (relht = 1.5f0)
        ptbali = ptba * (1.0f0 - pct / 100.0f0)
        tbal = (1.0f0 - pct / 100.0f0) * ba
        iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9); iicr < 1 && (iicr = 1)
        rcr = Float32(iicr)
        htgrl = CI_RG_CNST[sp] + CI_RG_HTBA[sp]*relht*ptba + CI_RG_PBAL[sp]*ptbali +
                CI_RG_STBA[sp]*ba + CI_RG_RLHT[sp]*relht + CI_RG_CRSQ[sp]*rcr*rcr +
                CI_RG_CR[sp]*rcr + CI_RG_BAL[sp]*tbal + CI_RG_HDM1[sp]*rhdm1 +
                CI_RG_HDM2[sp]*rhdm2 + CI_RG_PTBA[sp]*ptba
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        con = exp(c.htg_cor_small[sp])                    # CON = RHCON·exp(HCOR); RHCON=1
        h2 = h0 + htgrl * scale_h * xrhgro * con
        htgr1 = h2 - h0; htgr1 < 0.0f0 && (htgr1 = 0.0f0) # HTGR1 = HK−H (regent.f:907-908)
        htgr = htgr1                                      # + ZZRAN*0.1*SCALE deferred
        htgr < 0.1f0 && (htgr = 0.1f0)                    # Dixon 3/4/09 floor (regent.f:961)
        # Blend the small-tree (regent) HTGR with the large-tree HTG(K) via XWT (regent.f:973-979): for a tree
        # with XMN<D<XMAX the final height increment is a diameter-weighted mix of the regent prediction and the
        # large-tree htgf increment already in t.ht_growth[i]. jl previously used pure HTGR (XWT=0) — correct only
        # for D≤XMN; trees between XMN and XMAX (e.g. 2-5" CIVAR conifers) then under-grew (the large-tree htgf
        # increment is much bigger for young vigorous stems). Both HK and the DBH below flow from the blended HTG
        # (regent.f:998 HK=H+HTG(K)), matching live.
        xmn = CI_RG_XMIN[sp]; xmx = CI_RG_XMAX[sp]
        xwt = d0 <= xmn ? 0.0f0 : (d0 - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]
        hcap = s.control.sp_size_cap[sp, 4]               # SIZCAP(sp,4) max height (regent.f:983-986)
        (hcap > 0f0 && h0 + htg > hcap) && (htg = max(hcap - h0, 0.1f0))
        t.ht_growth[i] = htg
        hk = h0 + htg                                     # HK = H + HTG(K) (blended, regent.f:998)
        if hk < 4.5f0
            t.diam_growth[i] = 0.0f0
            _ci_rg_stash!(stash, t, i, d0 < 3.0f0)
            continue
        end
        dhcn = CI_RG_DHCN[sp]; dhht = CI_RG_DHHT[sp]; dhcr = CI_RG_DHCR[sp]
        dk = exp(dhcn + dhht * log(hk) + dhcr * log(rcr))
        dkk = h0 < 4.5f0 ? d0 : exp(dhcn + dhht * log(h0) + dhcr * log(rcr))   # regent.f:1009-1012
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        bark = ci_bratio(sd, sp, d0)                      # BRATIO(ISPC,DBH,HT)
        dg = (dk - dkk) * bark * xrdgro; dg < 0.0f0 && (dg = 0.0f0)            # regent.f:1245
        dds = dg * (2.0f0 * bark * d0 + dg) * scale2                          # regent.f:1251
        dg = sqrt((d0 * bark)^2 + dds) - bark * d0                            # regent.f:1253
        dg < 0.0f0 && (dg = 0.0f0)
        (d0 + dg) < CI_RG_DIAM[sp] && (dg = CI_RG_DIAM[sp] - d0)              # regent.f:1256 min-DBH floor
        t.diam_growth[i] = dg
        _ci_rg_stash!(stash, t, i, d0 < 3.0f0)   # copies use the regent growth (CIVAR/NIVAR); dgU/dgL only if dubbed
    end
    return s
end

# ci/esgent.f (CALL REGENT(.TRUE.,ITRNIN)) — grow the JUST-ESTABLISHED regen IN its birth cycle. CI was OMITTED
# from the esgent dispatch (simulate.jl had CR/TT/EM/UT), so planted/established CI seedlings never got their
# first-cycle height growth (CI BARE-PLANT: TopHt 2 vs live 12, BA ~half thru age 40). Same class as EM #137 /
# UT #184. Mirrors small_tree_growth!'s UTVAR + CIVAR height/DBH over the birth subperiod (subyr=FINT−GENTIM=5),
# applying HT/DBH directly (esgent.f HT(I)=HT(I)+HTG(I)·WK4). Gated to the new records nstart+1:n.
function ci_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0, avh_pre::Float32 = -1.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species
    nstart >= t.n && return s
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    # ci/regent.f:624 CIVAR RELHT uses ATAVH (grincr.f:318 ATAVH=AVH, grinit.f:242 ATAVH=0) = the PRE-regen stand
    # avg height; on a bare-plant stand ATAVH=AVH=0 ⇒ RELHT=1.5 (not HT/avh over the just-established seedlings).
    # #194: establish! recomputes avg_height INCLUDING the new regen before ci_esgent! runs, so p.avg_height is the
    # WRONG (post-regen) denominator — jl relht 0.634 vs live 1.5 ⇒ htgrl 4.08 vs live 5.26 (~22% under). Use the
    # captured pre-establishment avg height (avh_pre); its 0-on-bare ⇒ RELHT=1.5 = live.
    atavh = avh_pre >= 0.0f0 ? avh_pre : avh
    kodtyp = Int(p.habitat_code)
    rhdm1 = (500 <= kodtyp < 600) ? 1.0f0 : 0.0f0
    rhdm2 = (600 <= kodtyp < 700) ? 1.0f0 : 0.0f0
    regyr = CI_RG_REGYR
    gentim = max(fint - 5.0f0, 0.0f0)
    bscale = (fint - gentim) / regyr                     # birth-cycle fraction (WK4; =0.5 for fint=10)
    scale2 = 1.0f0
    cur_year = current_cycle_year(s)
    slo_a = sd[:site_lo]; shi_a = sd[:site_hi]
    xd = avh * (relden / 100.0f0); xd > 300.0f0 && (xd = 300.0f0)
    pctred = CI_RG_AB[1] + xd*(CI_RG_AB[2] + xd*(CI_RG_AB[3] + xd*(CI_RG_AB[4] + xd*(CI_RG_AB[5] + xd*CI_RG_AB[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for i in (nstart+1):t.n
        sp = Int(t.species[i]); d0 = t.dbh[i]
        d0 >= CI_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        (sp < 1 || sp > 19) && continue
        h0 = t.height[i]
        if _ci_ut_species(sp)                            # UTVAR aspen/juniper/MC/CW/hardwood
            sitear = p.sp_site_index[sp]; sj = sitear
            con = exp(c.htg_cor_small[sp])
            if sp == 13
                slo = slo_a[sp]; shi = shi_a[sp]
                si = sitear; si > shi && (si = shi); si <= slo && (si = slo + 0.5f0)
                relsi = (si - slo) / (shi - slo); rsimod = 0.5f0 * (1.0f0 + relsi)
                age = (h0 * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0)
                hite1 = 26.9825f0 * age^1.1752f0; hite2 = 26.9825f0 * (age + 10.0f0)^1.1752f0
                htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 0.75f0
            else
                pothtg = (sj / 5.0f0) * (sj * 1.5f0 - h0) / (sj * 1.5f0) * 0.83f0
                xcr = Float32(t.crown_pct[i]) / 100.0f0
                vigor = 150.0f0 * xcr^3 * exp(-6.0f0 * xcr) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
                sp == 14 && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)
                htgr = pothtg * pctred * vigor * con
            end
            htgr = htgr * bscale                         # birth-cycle subperiod (was scale_ut)
            htgr < 0.1f0 && (htgr = 0.1f0)
            xmn = CI_RG_XMIN[sp]; xmx = CI_RG_XMAX[sp]
            xwt = d0 <= xmn ? 0.0f0 : (d0 - xmn) / (xmx - xmn)
            htg = htgr * (1.0f0 - xwt); htg < 0.1f0 && (htg = 0.1f0)   # new tree: large-tree HTG(K)=0
            hcap = s.control.sp_size_cap[sp, 4]
            (hcap > 0f0 && h0 + htg > hcap) && (htg = max(hcap - h0, 0.1f0))
            hk = h0 + htg
            t.height[i] = hk; t.ht_growth[i] = htg
            bark = ci_bratio(sd, sp, d0)
            if hk > 4.5f0
                local dk::Float32, dkk::Float32
                if sp == 14
                    dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                    dkk = h0 < 4.5f0 ? d0 : (h0 - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                elseif sp == 15
                    dk = 3.1020f0 + 0.0210f0 * hk
                    dkk = 3.1020f0 + 0.0210f0 * h0; dkk < 0.0f0 && (dkk = d0); dk < dkk && (dk = dkk + 0.01f0)
                else
                    dk = _ci_ut_htdbh(hk)
                    dkk = h0 <= 4.5f0 ? d0 : _ci_ut_htdbh(h0)
                end
                dg = (dk - dkk) * bark; dg < 0.0f0 && (dg = 0.0f0)
                dgmx = CI_RG_DGMAX[sp] * bscale; dg > dgmx && (dg = dgmx)
                dds = dg * (2.0f0 * bark * d0 + dg) * scale2
                dg = sqrt((d0 * bark)^2 + dds) - bark * d0
                (d0 + dg) < CI_RG_DIAM[sp] && (dg = CI_RG_DIAM[sp] - d0)
                dg > 0.0f0 && (t.dbh[i] = d0 + dg; t.diam_growth[i] = dg)
            end
            continue
        end
        # CIVAR conifers (HTGRL regression)
        pt = Int(t.plot_id[i])
        ptba = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : ba
        pct = t.crown_ratio[i]
        relht = atavh > 0f0 ? h0 / atavh : 1.5f0; relht > 1.5f0 && (relht = 1.5f0)   # ci/regent.f:624 ATAVH (pre-regen)
        ptbali = ptba * (1.0f0 - pct / 100.0f0)
        tbal = (1.0f0 - pct / 100.0f0) * ba
        iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9); iicr < 1 && (iicr = 1)
        rcr = Float32(iicr)
        htgrl = CI_RG_CNST[sp] + CI_RG_HTBA[sp]*relht*ptba + CI_RG_PBAL[sp]*ptbali +
                CI_RG_STBA[sp]*ba + CI_RG_RLHT[sp]*relht + CI_RG_CRSQ[sp]*rcr*rcr +
                CI_RG_CR[sp]*rcr + CI_RG_BAL[sp]*tbal + CI_RG_HDM1[sp]*rhdm1 +
                CI_RG_HDM2[sp]*rhdm2 + CI_RG_PTBA[sp]*ptba
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        con = exp(c.htg_cor_small[sp])
        h2 = h0 + htgrl * bscale * xrhgro * con          # birth-cycle fraction (was scale_h)
        htgr1 = h2 - h0; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        htgr = htgr1 < 0.1f0 ? 0.1f0 : htgr1
        xmn = CI_RG_XMIN[sp]; xmx = CI_RG_XMAX[sp]
        xwt = d0 <= xmn ? 0.0f0 : (d0 - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt)                       # new tree: large-tree HTG(K)=0
        hcap = s.control.sp_size_cap[sp, 4]
        (hcap > 0f0 && h0 + htg > hcap) && (htg = max(hcap - h0, 0.1f0))
        hk = h0 + htg
        t.height[i] = hk; t.ht_growth[i] = htg
        if hk >= 4.5f0
            dhcn = CI_RG_DHCN[sp]; dhht = CI_RG_DHHT[sp]; dhcr = CI_RG_DHCR[sp]
            dk = exp(dhcn + dhht * log(hk) + dhcr * log(rcr))
            dkk = h0 < 4.5f0 ? d0 : exp(dhcn + dhht * log(h0) + dhcr * log(rcr))
            xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
            bark = ci_bratio(sd, sp, d0)
            dg = (dk - dkk) * bark * xrdgro; dg < 0.0f0 && (dg = 0.0f0)
            dds = dg * (2.0f0 * bark * d0 + dg) * scale2
            dg = sqrt((d0 * bark)^2 + dds) - bark * d0; dg < 0.0f0 && (dg = 0.0f0)
            (d0 + dg) < CI_RG_DIAM[sp] && (dg = CI_RG_DIAM[sp] - d0)
            dg > 0.0f0 && (t.dbh[i] = d0 + dg; t.diam_growth[i] = dg)
        end
    end
    return s
end
