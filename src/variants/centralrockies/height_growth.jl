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
