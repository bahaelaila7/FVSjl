# =============================================================================
# height_growth.jl (centralrockies) — CR height growth = GENGYM (cr/htgf.f + gemht.f)
#
# gemht.f computes an EVEN-AGED estimate HHE (site-index curves, by IMODTY 1-5)
# and an UNEVEN-AGED estimate HHU (GENGYM diameter-based, SELECT CASE(species));
# surrogate species (GF/MH/RC/WL/LM) instead return a DIRECT increment HTGI with
# IHTG=1. htgf.f meshes HHE/HHU (twin calls at current + future age) into HTG.
#
# Math: ALOG->flog, EXP->fexp, **realexp->fpow, **2 (int)->x*x, AMAX1->max.
# =============================================================================

"""
    cr_gemht(imodty, si, dp, h, is, bat, ap, bautba, ccf, dgi, bark) -> (hhe, hhu, ihtg, htgi)

GENGYM height (cr/gemht.f). `HHE`=even-aged height estimate, `HHU`=uneven-aged; `IHTG=1` + `HTGI`
for surrogate species that return a direct increment. `BATEM` is set in the IMODTY block and reused
(then possibly reset) per species case — faithfully carried.
"""
function cr_gemht(imodty::Int, si::Float32, dp::Float32, h::Float32, is::Int, bat::Float32,
                  ap::Float32, bautba::Float32, ccf::Float32, dgi::Float32, bark::Float32)
    ihtg = 0; htgi = 0.0f0; hhe = 0.0f0; hhu = 0.0f0
    batem = 0.0f0

    # ---- EVEN-AGED HHE by model type ----
    if imodty == 1                                            # SW mixed conifer
        batem = bat * 0.01f0; batem < 0.8f0 && (batem = 0.8f0)
        agetem = ap; agetem > 210.0f0 && (agetem = 210.0f0)
        si < 80.0f0 && (agetem = max(agetem, 110.0f0 - si))
        hg = 109.559129f0 * fpow(1.0f0 - 0.975884f0 * fexp(-0.014377f0 * agetem), 1.289266f0) + 4.5f0
        hl =  72.512644f0 * fpow(1.0f0 - 0.876961f0 * fexp(-0.020066f0 * agetem), 2.016632f0)
        hhe = -((hg - hl) * ((82.488f0 - si) / 26.279f0)) + hg
        if si < 80.0f0
            ap < agetem && (hhe = ((hhe - 4.5f0) / agetem) * ap + 4.5f0)
        else
            ap < 20.0f0 && (hhe = (0.02348f0 * si - 0.93429f0) * ap + 4.5f0)
        end
        ratio = 1.0f0 - bautba; ratio < 0.768f0 && (ratio = 0.768f0)
        hhe = hhe * ratio
    elseif imodty == 2                                        # SW ponderosa pine
        batem = bat * 0.01f0; batem < 0.8f0 && (batem = 0.8f0)
        agetem = ap; agetem > 210.0f0 && (agetem = 210.0f0)
        si < 80.0f0 && (agetem = max(agetem, 120.0f0 - si))
        hg = 106.493954f0 * fpow(1.0f0 - 0.938775f0 * fexp(-0.016066f0 * agetem), 1.550720f0) + 4.5f0
        hl =  78.078735f0 * fpow(1.0f0 - 0.843715f0 * fexp(-0.020412f0 * agetem), 2.280435f0)
        hhe = -((hg - hl) * ((81.5585f0 - si) / 21.7149f0)) + hg
        if si < 80.0f0
            ap < agetem && (hhe = ((hhe - 4.5f0) / agetem) * ap + 4.5f0)
        else
            ap < 20.0f0 && (hhe = (0.02463f0 * si - 1.1025f0) * ap + 4.5f0)
        end
        ratio = 1.0f0 - bautba; ratio < 0.768f0 && (ratio = 0.768f0)
        hhe = hhe * ratio
    elseif imodty == 3                                        # Black Hills ponderosa pine
        batem = bat; batem < 20.0f0 && (batem = 20.0f0)
        agetem = ap
        htmax = (si + 0.3846f0) * 1.2999886f0
        if h >= htmax
            hhe = h + 0.1f0
        elseif agetem > 165.0f0
            hhe = h + 0.1f0
        else
            temht = (si + 0.3846f0) * (-0.5234f0 + 1.8234f0 * fexp(-fpow(1.0989f0 - 0.006105f0 * agetem, 2.35f0)))
            age10 = agetem + 10.0f0
            hhe = (si + 0.3846f0) * (-0.5234f0 + 1.8234f0 * fexp(-fpow(1.0989f0 - 0.006105f0 * age10, 2.35f0)))
            ratio = 1.0f0 - bautba; ratio < 0.793f0 && (ratio = 0.793f0)
            phg = hhe - temht; phg < 0.0f0 && (phg = 0.0f0)
            hhe = h + phg * ratio
        end
    elseif imodty == 4                                        # spruce-fir
        batem = bat; batem < 20.0f0 && (batem = 20.0f0)
        agetem = ap; agetem < 30.0f0 && (agetem = 30.0f0)
        hhe = (2.75780f0 * fpow(si, 0.83312f0)) *
              fpow(1.0f0 - fexp(-0.015701f0 * agetem), 22.71944f0 * fpow(si, -0.63557f0)) + 4.5f0
        ap < agetem && (hhe = ((hhe - 4.5f0) / agetem) * ap + 4.5f0)
        ratio = 1.0f0 - bautba; ratio < 0.728f0 && (ratio = 0.728f0)
        hhe = hhe * ratio
    elseif imodty == 5                                        # lodgepole pine
        batem = bat; batem < 20.0f0 && (batem = 20.0f0)
        agetem = ap; aptem = agetem
        agetem < 30.0f0 && (agetem = 30.0f0)
        agetem > 200.0f0 && (agetem = 200.0f0)
        ccftem = ccf - 125.0f0; ccftem < 0.0f0 && (ccftem = 0.0f0)
        hhe = 9.89331f0 - 0.19177f0 * agetem + 0.00124f0 * agetem * agetem -
              0.00082f0 * ccftem * si + 0.01387f0 * agetem * si -
              0.0000455f0 * agetem * agetem * si
        aptem <= 30.0f0 && (hhe = (hhe / agetem) * aptem)
        ratio = 1.0f0 - bautba; ratio <= 0.742f0 && (ratio = 0.742f0)
        hhe = hhe * ratio
    else
        hhe = 0.0f0
    end

    # ---- UNEVEN-AGED HHU / direct HTGI by species ----
    if is == 1 || is == 2                                     # sub-alpine / corkbark fir
        if imodty <= 2
            hhu = 4.514294f0 * fpow(si, 0.755380f0) *
                  fpow(1.0f0 - fexp(-0.080869f0 * dp), 1.409884f0 * fpow(batem, 0.003919f0)) + 4.5f0
        else
            hhu = (15.5f0 + 1.1f0 * si) *
                  fpow(1.0f0 - fexp(-0.097152f0 * dp), 4.698567f0 * fpow(batem, -0.252630f0)) + 4.5f0
        end
    elseif is == 3                                            # Douglas-fir
        batem = bat; batem < 0.8f0 && (batem = 0.8f0)
        hhu = 13.096420f0 * fpow(si, 0.480509f0) *
              fpow(1.0f0 - fexp(-0.077408f0 * dp), 1.237589f0 * fpow(batem, -0.063297f0)) + 4.5f0
    elseif is == 4                                            # grand fir (direct)
        con = 2.03035f0 - 0.6458f0 - 0.00013358f0 * h * h - 0.09775f0 * flog(dp) + 0.23315f0 * flog(h)
        htgi = fexp(con + 0.62144f0 * flog(dgi)) + 0.4809f0
        temsi = si; temsi < 40.0f0 && (temsi = 40.0f0); temsi > 120.0f0 && (temsi = 120.0f0)
        adjfac = temsi <= 80.0f0 ? 0.1f0 + 0.0125f0 * temsi : -0.7f0 + 0.0225f0 * temsi
        htgi *= adjfac; ihtg = 1
    elseif is == 5                                            # white fir
        batem = bat; batem < 0.8f0 && (batem = 0.8f0)
        hhu = 13.822088f0 * fpow(si, 0.462393f0) *
              fpow(1.0f0 - fexp(-0.075766f0 * dp), 1.312638f0 * fpow(batem, -0.040708f0)) + 4.5f0
    elseif is == 6                                            # mountain hemlock (direct)
        con = 1.74090f0 - 0.6458f0 - 0.0000446f0 * h * h - 0.09775f0 * flog(dp) + 0.23315f0 * flog(h)
        htgi = fexp(con + 0.34003f0 * flog(dgi)) + 0.4809f0
        temsi = si; temsi < 40.0f0 && (temsi = 40.0f0); temsi > 70.0f0 && (temsi = 70.0f0)
        adjfac = 0.36f0 + 0.012f0 * temsi
        htgi *= adjfac; ihtg = 1
    elseif is == 7                                            # western red cedar (direct)
        con = 2.21104f0 - 0.9941f0 - 0.00003631f0 * h * h - 0.1219f0 * flog(dp) + 0.23315f0 * flog(h)
        htgi = fexp(con + 0.37042f0 * flog(dgi)) + 0.4809f0
        temsi = si; temsi < 40.0f0 && (temsi = 40.0f0); temsi > 110.0f0 && (temsi = 110.0f0)
        adjfac = 0.0875f0 + 0.01375f0 * temsi; adjfac < 1.0f0 && (adjfac = 1.0f0)
        htgi *= adjfac; ihtg = 1
    elseif is == 8                                            # western larch (direct)
        con = 1.81759f0 + 0.1433f0 - 0.00002607f0 * h * h - 0.3899f0 * flog(dp) + 0.23315f0 * flog(h)
        htgi = fexp(con + 0.75756f0 * flog(dgi)) + 0.4809f0
        temsi = si; temsi < 40.0f0 && (temsi = 40.0f0); temsi > 120.0f0 && (temsi = 120.0f0)
        adjfac = 0.23337f0 + 0.008333f0 * temsi
        htgi *= adjfac; ihtg = 1
    elseif is == 10                                           # limber pine (direct, logistic)
        ihtg = 1
        d10 = dp + dgi / bark
        if dp > 45.1f0 || h > 94.5f0 || d10 > 45.1f0
            htgi = 0.1f0
        else
            y1 = (dp - 0.1f0) / 45.0f0
            y2 = (h - 4.5f0) / 90.0f0
            fby1 = flog(y1 / (1.0f0 - y1))
            fby2 = flog(y2 / (1.0f0 - y2))
            z = (0.30546f0 + 0.94823f0 * fby2 - 0.70453f0 * (1.64770f0 + 1.35015f0 * fby1)) *
                fpow(1.0f0 - 0.70453f0 * 0.70453f0, -0.5f0)
            psi = 2.46480f0 * fpow((d10 - 0.1f0) / (45.1f0 - d10), 1.00316f0) *
                  fexp(z * fpow(1.0f0 - 0.70453f0 * 0.70453f0, 0.5f0) / 0.94823f0)
            h2 = ((psi / (1.0f0 + psi)) * 90.0f0) + 4.5f0
            h2 < h && (h2 = h + 0.1f0)
            htgi = h2 - h
            temsi = si; temsi < 20.0f0 && (temsi = 20.0f0); temsi > 60.0f0 && (temsi = 60.0f0)
            adjfac = temsi < 40.0f0 ? 0.2f0 + 0.015f0 * temsi : -0.1f0 + 0.0225f0 * temsi
            htgi = htgi * adjfac
        end
    elseif is == 11 || is == 14                               # lodgepole / whitebark pine
        batem = bat; batem < 20.0f0 && (batem = 20.0f0)
        hhu = (8.5f0 + 1.1f0 * si) *
              fpow(1.0f0 - fexp(-0.085004f0 * dp), 1.709643f0 * fpow(batem, -0.163186f0)) + 4.5f0
    elseif is == 9 || is == 12 || is == 16 || (23 <= is <= 27) || (29 <= is <= 35) || is == 37
        hhu = 42.269377f0 * fpow(1.0f0 - fexp(-0.165687f0 * dp), 1.184734f0) + 4.5f0
    elseif is == 13 || is == 36                               # ponderosa / Chihuahua pine
        if imodty == 1
            hhu = 24.244690f0 * fpow(si, 0.343864f0) *
                  fpow(1.0f0 - fexp(-0.069180f0 * dp), 1.251384f0 * fpow(batem, -0.272018f0)) + 4.5f0
        elseif imodty == 2
            hhu = 40.78321f0 * fpow(si, 0.332614f0) *
                  fpow(1.0f0 - fexp(-0.021471f0 * dp), 0.922811f0 * fpow(batem, -0.133923f0)) + 4.5f0
        else
            batem = bat * 0.01f0; batem < 1.0f0 && (batem = 1.0f0)
            hhu = 32.108633f0 * fpow(si, 0.276926f0) *
                  fpow(1.0f0 - fexp(-0.057766f0 * dp), 0.984340f0 * fpow(batem, -0.169876f0)) + 4.5f0
        end
    elseif is == 15                                           # western white pine
        batem = bat; batem < 0.8f0 && (batem = 0.8f0)
        hhu = 18.967185f0 * fpow(si, 0.379790f0) *
              fpow(1.0f0 - fexp(-0.071482f0 * dp), 1.159608f0 * fpow(batem, -0.099449f0)) + 4.5f0
    elseif is == 17 || is == 19                               # blue / white spruce
        batem = bat; batem < 0.8f0 && (batem = 0.8f0)
        hhu = 54.180173f0 * fpow(si, 0.177962f0) *
              fpow(1.0f0 - fexp(-0.089253f0 * dp), 1.533535f0 * fpow(batem, -0.028852f0)) + 4.5f0
    elseif is == 18                                           # Engelmann spruce
        if imodty <= 2
            hhu = 10.616238f0 * fpow(si, 0.549461f0) *
                  fpow(1.0f0 - fexp(-0.087283f0 * dp), 1.488355f0 * fpow(batem, -0.027226f0)) + 4.5f0
        else
            hhu = (15.5f0 + 1.1f0 * si) *
                  fpow(1.0f0 - fexp(-0.110383f0 * dp), 6.262899f0 * fpow(batem, -0.286055f0)) + 4.5f0
        end
    elseif is == 20 || is == 21 || is == 22 || is == 28 || is == 38   # aspen/birch/cottonwood/OH
        if imodty <= 2
            hhu = 14.187987f0 * fpow(si, 0.416525f0) *
                  fpow(1.0f0 - fexp(-0.126806f0 * dp), 1.310744f0 * fpow(batem, -0.245126f0)) + 4.5f0
        else
            batem = bat; batem < 10.0f0 && (batem = 10.0f0)
            hhu = (-2.04f0 + 1.4534f0 * si) *
                  fpow(1.0f0 - fexp(-0.058112f0 * dp), 1.894400f0 * fpow(batem, -0.192979f0)) + 4.5f0
        end
    else
        hhu = 0.0f0
    end

    return hhe, hhu, ihtg, htgi
end

# HTOSI(6, MAXSP) — height-to-site adjust factor by (IMODTY 1-6, species) (cr/htgf.f DATA).
const _CR_HTOSI = Float32[
  1.25 1.25 1.18 1.00 1.20 1.00 1.00 1.00 1.00 1.00 1.20 1.00 1.20 1.20 1.20 1.00 1.20 1.20 1.20 1.00 1.15 1.15 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.20 1.15 1.15
  1.25 1.25 1.18 1.00 1.20 1.00 1.00 1.00 1.00 1.00 1.15 1.00 1.20 1.15 1.20 1.00 1.20 1.20 1.20 1.00 1.15 1.15 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.20 1.15 1.15
  1.13 1.13 1.10 1.10 1.13 1.10 1.10 1.10 1.00 1.07 1.13 1.00 1.10 1.10 1.10 1.00 1.15 1.15 1.15 1.00 1.05 1.05 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.10 1.05 1.15
  1.15 1.15 1.18 1.00 1.15 1.00 1.00 1.00 1.00 1.00 1.10 1.00 1.15 1.05 1.15 1.00 1.10 1.10 1.10 1.00 1.05 1.05 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.15 1.05 1.05
  1.25 1.25 1.25 1.10 1.25 1.10 1.10 1.10 1.00 1.07 1.30 1.00 1.25 1.20 1.20 1.00 1.20 1.20 1.20 1.00 1.18 1.18 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.25 1.15 1.20
  1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ]

# Breast-high-age adjustment for AP (cr/htgf.f:172-190), IMODTY 1/2/4 only.
@inline function _cr_bhage_adjust(ap::Float32, imodty::Int, is::Int, ssite::Float32)
    (imodty == 1 || imodty == 2 || imodty == 4) || return ap
    if is == 20 || is == 28 || is == 38 || is == 14
        ap -= 4.5f0 / (0.1f0 + ssite / 50.0f0)
    elseif imodty == 1
        tsite = ssite < 30.0f0 ? 30.0f0 : ssite
        ap -= 4.5f0 / (-0.642f0 + 0.02285f0 * tsite)
    elseif imodty == 2
        ap -= 4.5f0 / (0.25f0 + 0.00467f0 * ssite)
    elseif imodty == 4
        tsite = ssite < 20.0f0 ? 20.0f0 : ssite
        ap -= 4.5f0 / (-0.22f0 + 0.0155f0 * tsite)
    end
    ap < 1.0f0 && (ap = 1.0f0)
    return ap
end

"""
    height_growth!(s, ::CentralRockies; scale=1)

CR periodic height growth (cr/htgf.f). Per tree ≥0.5" DBH and >4.5 ft: twin GEMHT calls (current + future
age) meshed into HTG (even-aged, or even/uneven blend when AGERNG>40 & BA≥70), aspen/birch ASPFAC, HTCON
calibration × SCALE × XHMULT, mistletoe (=1), SIZCAP(·,4) cap. Small trees defer to regent (HTG=0 for now).
The stochastic ZZRAN term fires only when DGSD>0 (crt01 default 0 ⇒ deterministic). Writes `t.ht_growth`.
"""
function height_growth!(s::StandState, ::CentralRockies; scale::Float32 = 1.0f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    imodty = Int(p.model_type)
    ba = p.basal_area
    dens = s.density
    agerng = _cr_agerng(t)
    cur_year = current_cycle_year(s)
    dgsd = s.control.dg_sd

    # BADIST BAU (BA-above-class) — same pre-pass as dgf! (cr/htgf.f reads BAU from DGF's BADIST).
    bau = _cr_badist_bau(t)

    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0.0f0
        t.tpa[i] <= 0.0f0 && continue
        d = t.dbh[i]; hnow = t.height[i]
        (d < 0.5f0 || hnow <= 4.5f0) && continue          # small trees → regent (not yet ported)
        sp = Int(t.species[i])
        ssite = p.sp_site_index[sp]
        icls = trunc(Int, d + 1.0f0); icls > 41 && (icls = 41)
        bark = cr_bratio(sd, sp, d, imodty)
        xhmult = active_multiplier(s.control, :htg, sp, cur_year)
        # ZZRAN stochastic increment (htgf.f:288-294): rejection-sampled BACHLO(0,1,RANN) with |z|≤DGSD.
        # NOTE: FVS draws in species-sorted IND1 order; this loop is tree-index order — the per-tree draw
        # VALUES therefore only bit-match live once the cycle RNG sequence is reconciled (chunk-9 concern).
        zzran = 0.0f0
        if dgsd > 0.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= dgsd && zzran >= -dgsd) && break
            end
        end
        t.ht_growth[i] = _cr_htg_tree(imodty, sp, ssite, _CR_HTOSI[imodty, sp], d, hnow, bark,
            bau[icls] / ba, dens.point_ccf[Int(t.plot_id[i])], t.diam_growth[i], t.birth_age[i],
            agerng, ba, t.crown_ratio[i], c.htg_cor[sp], scale, xhmult, s.control.sp_size_cap[sp, 4], zzran)
    end
    return s
end

"""
    _cr_htg_tree(...) -> HTG

The per-tree cr/htgf.f height-increment computation (single source of truth for `height_growth!`
and the validation replay). `birth_age` = raw ABIRTH (the BH-age adjust + LPP-H30 use it directly).
Deterministic: the ZZRAN stochastic term (DGSD>0) is omitted — the caller adds it when enabled.
"""
function _cr_htg_tree(imodty::Int, sp::Int, ssite::Float32, adjust::Float32, d::Float32, hnow::Float32,
                      bark::Float32, bautba::Float32, pccfi::Float32, dgi::Float32, birth_age::Float32,
                      agerng::Float32, ba::Float32, pct::Float32, htcon::Float32, scale::Float32,
                      xhmult::Float32, cap::Float32, zzran::Float32 = 0.0f0)::Float32
    ap = _cr_bhage_adjust(birth_age, imodty, sp, ssite)
    hhe1, hhu1, ihtg, htgi = cr_gemht(imodty, ssite, d, hnow, sp, ba, ap, bautba, pccfi, dgi, bark)
    local htg::Float32
    if ihtg == 1
        htg = htgi * adjust
    else
        agefut = ap + 10.0f0
        dfut = d + dgi / bark
        hhe2, hhu2, _, _ = cr_gemht(imodty, ssite, dfut, hnow, sp, ba, agefut, bautba, pccfi, dgi, bark)
        if imodty == 3
            hhe2 = hhe1; hhe1 = hnow
        end
        htg = (hhe2 - hhe1) * adjust
        if imodty == 5 && birth_age <= 31.0f0
            tccf = pccfi; tccf <= 125.0f0 && (tccf = 0.0f0)
            h30 = 5.25621f0 + 0.37515f0 * ssite - 0.00082f0 * tccf * ssite
            h30 = h30 * birth_age / 31.0f0
            thtg = h30 - hnow
            thtg > htg && (htg = thtg)
        end
        if agerng > 40.0f0 && ba >= 70.0f0
            if pct <= 10.0f0
                htg = hhu2 - hhu1
            elseif pct > 10.0f0 && pct < 40.0f0
                hge = (hhe2 - hhe1) * adjust
                hgu = hhu2 - hhu1
                xwt = ((pct - 10.0f0) * (10.0f0 / 3.0f0)) / 100.0f0
                htg = xwt * hge + (1.0f0 - xwt) * hgu
            end
        end
        if sp == 20 || sp == 28
            temsi = ssite; temsi < 30.0f0 && (temsi = 30.0f0); temsi > 90.0f0 && (temsi = 90.0f0)
            htg = htg * (0.6253f0 + 0.00583f0 * temsi)
        end
    end
    htg = htg + zzran * 0.1f0                          # ZZRAN stochastic increment (htgf.f:295)
    htg < 0.1f0 && (htg = 0.1f0)                        # floor AFTER the random add (htgf.f:301)
    htg = htg * fexp(htcon) * scale * xhmult            # HTCON calib × SCALE × XHMULT; MISHGF=1
    if hnow + htg > cap
        htg = cap - hnow; htg < 0.1f0 && (htg = 0.1f0)
    end
    return htg
end

# BADIST BA-above-dbh-class (cr/badist.f) — same as dgf!'s inline pre-pass; summed in tree-index
# order (i=1:n) for Float32 bit-exactness. htgf reads BAU from DGF's earlier BADIST.
function _cr_badist_bau(t)
    bau = zeros(Float32, 41); totba = 0.0f0
    @inbounds for i in 1:t.n
        t.height[i] < 4.5f0 && continue
        tdbh = t.dbh[i]; icls = trunc(Int, tdbh + 1.0f0); icls > 41 && (icls = 41)
        tdbh < 1.0f0 && (tdbh = 1.0f0)
        treeba = 0.0054542f0 * tdbh * tdbh * t.tpa[i]
        totba += treeba; bau[icls] += treeba
    end
    bau[1] = totba - bau[1]; bau[1] < 0.0f0 && (bau[1] = 0.0f0)
    @inbounds for j in 2:41
        bau[j] = bau[j-1] - bau[j]; bau[j] < 0.0f0 && (bau[j] = 0.0f0)
    end
    return bau
end

# AGERNG (cr/dgf.f + htgf.f): birth-age range over established trees (ABIRTH>1, HT>4.5); else 1000.
function _cr_agerng(t)
    young = 1000.0f0; old = 0.0f0; anyage = false
    @inbounds for i in 1:t.n
        ab = t.birth_age[i]
        (ab <= 1.0f0 || t.height[i] <= 4.5f0) && continue
        anyage = true; ab < young && (young = ab); ab > old && (old = ab)
    end
    return anyage ? abs(old - young) : 1000.0f0
end

# cr_fndag (cr/fndag.f) — invert the even-aged site height curve to find the total age (SITAGE) at which a
# tree of species `sp`/site `site` reaches height `rht`. Linear search AP=10 step 5 up to AGEMAX (fndag.f:100),
# stop when HH≥RHT or |HH−RHT|<2, + breast-high adjust (fndag.f:9990). Aspen/birch/pinyon use a closed form.
# Used by _cr_dub_ages! to dub ABIRTH for un-aged inventory trees (cratet.f:552) — else htgf's AP floors to 1
# and tall trees over-grow height 2-3× (the TopHt divergence).
function cr_fndag(imodty::Int, site::Float32, rht::Float32, bautba::Float32, relden::Float32, sp::Int)::Float32
    if sp == 20 || sp == 28 || sp == 38 || sp == 14 ||
       (imodty == 3 && (sp == 21 || sp == 22 || sp == 16 || (29 <= sp <= 32) || sp == 37))
        agetem = fpow(rht * 12.0f0 * 2.54f0 / 26.9825f0, 0.8509f0)
        return agetem + 4.5f0 / (0.1f0 + site / 50.0f0)
    end
    agemax = imodty == 3 ? 165.0f0 : imodty == 5 ? 200.0f0 : 210.0f0
    m = 1 <= imodty <= 5 ? imodty : 1                          # IMODTY 6 → 1 (fndag.f:6000 GO TO 1000)
    ap = 10.0f0; agetem = ap
    while true
        agetem = ap; agetem > agemax && (agetem = agemax)
        local hh::Float32
        if m == 3
            htmax = (site + 0.3846f0) * 1.2999886f0
            rht >= htmax && return 165.0f0
            hh = (site + 0.3846f0) * (-0.5234f0 + 1.8234f0 * fexp(-fpow(1.0989f0 - 0.006105f0 * agetem, 2.35f0)))
            hh *= max(1.0f0 - bautba, 0.793f0)
        elseif m == 4
            agetem < 30.0f0 && (agetem = 30.0f0)
            hh = (2.75780f0 * fpow(site, 0.83312f0)) *
                 fpow(1.0f0 - fexp(-0.015701f0 * agetem), 22.71944f0 * fpow(site, -0.63557f0)) + 4.5f0
            hh *= max(1.0f0 - bautba, 0.728f0)
        elseif m == 5
            agetem < 30.0f0 && (agetem = 30.0f0)
            ccftem = relden - 125.0f0; ccftem < 0.0f0 && (ccftem = 0.0f0)
            hh = 9.89331f0 - 0.19177f0 * agetem + 0.00124f0 * agetem * agetem -
                 0.00082f0 * ccftem * site + 0.01387f0 * agetem * site - 0.0000455f0 * agetem * agetem * site
            hh *= max(1.0f0 - bautba, 0.742f0)
        elseif m == 1
            hg = 109.559129f0 * fpow(1.0f0 - 0.975884f0 * fexp(-0.014377f0 * agetem), 1.289266f0) + 4.5f0
            hl = 72.512644f0 * fpow(1.0f0 - 0.876961f0 * fexp(-0.020066f0 * agetem), 2.016632f0)
            hh = -((hg - hl) * ((82.488f0 - site) / 26.279f0)) + hg
            hh *= max(1.0f0 - bautba, 0.768f0)
        else                                                   # m == 2
            hg = 106.493954f0 * fpow(1.0f0 - 0.938775f0 * fexp(-0.016066f0 * agetem), 1.550720f0) + 4.5f0
            hl = 78.078735f0 * fpow(1.0f0 - 0.843715f0 * fexp(-0.020412f0 * agetem), 2.280435f0)
            hh = -((hg - hl) * ((81.5585f0 - site) / 21.7149f0)) + hg
            hh *= max(1.0f0 - bautba, 0.768f0)
        end
        (abs(hh - rht) < 2.0f0 || hh > rht) && break
        ap += 5.0f0
        ap > agemax && return agemax
    end
    tage = agetem
    if imodty == 1
        tsite = site < 30.0f0 ? 30.0f0 : site; tage += 4.5f0 / (-0.642f0 + 0.02285f0 * tsite)
    elseif imodty == 2
        tage += 4.5f0 / (0.25f0 + 0.00467f0 * site)
    elseif imodty == 4
        tsite = site < 20.0f0 ? 20.0f0 : site; tage += 4.5f0 / (-0.22f0 + 0.0155f0 * tsite)
    end
    return tage
end

# _cr_dub_ages! (cr/cratet.f:540-552 + FINDAG) — dub ABIRTH from the current height for inventory trees with no
# measured age (birth_age ≤ 0 or > 999). Runs once at setup (before the first height growth); gradd.f:205 then
# increments birth_age by FINT each cycle. Site = per-species SITEAR; BAUTBA = BAU(dbh-class)/BA.
function _cr_dub_ages!(s::StandState)
    p, t = s.plot, s.trees
    imodty = Int(p.model_type)
    bau = _cr_badist_bau(t); ba = p.basal_area <= 0.0f0 ? 25.0f0 : p.basal_area
    relden = stand_ccf(s)
    @inbounds for i in 1:t.n
        ab = t.birth_age[i]
        (ab > 0.0f0 && ab <= 999.0f0) && continue              # already aged (cratet.f:541)
        sp = Int(t.species[i]); h = t.height[i]
        h <= 0.0f0 && continue
        icls = trunc(Int, t.dbh[i] + 1.0f0); icls > 41 && (icls = 41)
        bautba = ba > 0.0f0 ? bau[icls] / ba : 0.0f0; bautba < 0.0f0 && (bautba = 0.0f0)
        sitage = cr_fndag(imodty, p.sp_site_index[sp], h, bautba, relden, sp)
        sitage > 0.0f0 && (t.birth_age[i] = sitage; t.age_known[i] = true)
    end
    return s
end
