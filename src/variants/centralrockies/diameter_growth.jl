# =============================================================================
# diameter_growth.jl (centralrockies) — CR large-tree DDS = GENGYM (cr/gemdg.f)
#
# Faithful transcription of GEMDG: dispatch on IMODTY (model type 1-5) then
# SELECT CASE(species) with per-species regression equations. Two output paths:
#   IDDS=0  -> compute DF (forecast dbh);  DIAGR = (DF - DPP)*BARK
#   IDDS=1  -> compute DDS (ln change);    DIAGR = sqrt(exp(DDS)+(DPP*BARK)^2) - DPP*BARK
# common tail: DSTAG (ISTAGF, inactive by default), BH-ponderosa x0.80, then
# convert DIAGR back to ln(DDS) floored at -9.21. Returns DDS (ln change in
# squared inside-bark diameter). The dgf.f wrapper adds COR + DGCON (next step).
#
# Math: ALOG->flog, EXP->fexp, **2.0->fpow(.,2f0), SQRT/SIN/COS->Base (bit-exact).
# =============================================================================

const _LN2F = 2.0f0

"GENGYM DDS (cr/gemdg.f). `dbhmax` = per-species DBHMAX (cr/sitset.f); `agerng` (dgf.f=1000, cratet age range);
`istagf` per-species stagnation flag (0 by default => DSTAG inert). Returns ln(DDS)."
function cr_gemdg(imodty::Int, is::Int, bautba::Float32, spba::Float32, si::Float32, dp::Float32,
                  bat::Float32, bark::Float32, cr::Float32, slope::Float32, aspect::Float32,
                  pbal::Float32, pccfi::Float32, relden::Float32, bal::Float32;
                  dbhmax::Float32, agerng::Float32 = 1000.0f0, istagf::Int = 0, dstag::Float32 = 1.0f0)::Float32
    l2(x) = fpow(flog(x), _LN2F)                    # (ALOG(x))**2.0
    dpp = dp < 1.0f0 ? 1.0f0 : dp
    batem = bat
    bgttba = bautba
    idds = 0
    bal100 = bal / 100.0f0
    df = 0.0f0; dds = 0.0f0

    # --- IMODTY bounds (cr/gemdg.f:59-83) ---
    if imodty == 1
        batem < 65.0f0 && (bgttba = 0.0f0)
        batem < 1.0f0  && (batem = 1.0f0)
    elseif imodty == 2
        batem < 21.0f0 && (batem = 21.0f0)
    elseif imodty == 3
        batem < 5.0f0  && (batem = 5.0f0)
    elseif imodty == 4
        batem < 5.0f0  && (batem = 5.0f0)
    elseif imodty == 5
        batem < 14.0f0 && (batem = 14.0f0)
    end

    # --- SELECT CASE(species) (cr/gemdg.f:86-345) ---
    if is == 1 || is == 2                       # subalpine/corkbark fir
        if imodty <= 2
            df = 0.50283f0 + 0.93999f0*dpp - 0.0016822f0*bat + 0.00791f0*si + 0.16235f0*l2(dpp)
        else
            df = 2.19459f0 + 0.94659f0*dpp - 0.20410f0*bgttba - 0.35970f0*flog(batem) +
                 0.11326f0*l2(dpp) + 0.00615f0*si
        end
    elseif is == 3                              # Douglas-fir
        df = 0.58234f0 + 0.94542f0*dpp - 0.0025092f0*bat + 0.0070474f0*si -
             0.14044f0*bgttba + 0.15738f0*l2(dpp)
    elseif is == 4                              # grand fir (CI variant eqn)
        dds = -0.057982f0 + 0.977362f0 + 0.005045f0*60.0f0 - 0.000093f0*60.0f0*60.0f0 +
              0.009335f0*sin(aspect)*slope - 0.004469f0*cos(aspect)*slope - 0.033374f0*slope -
              0.418343f0*slope*slope + 1.286963f0*flog(dpp) - 0.217923f0*flog(bat) +
              1.175105f0*cr + 0.219013f0*cr*cr - 0.0004408f0*dpp*dpp -
              0.000578f0*pbal/flog(dpp+1.0f0) - 0.000512f0*pccfi
        idds = 1
    elseif is == 5                              # white fir
        df = 0.46172f0 + 1.04089f0*dpp - 0.0019984f0*bat + 0.0065947f0*si -
             0.09197f0*bgttba - 0.0011355f0*dpp*dpp
    elseif is == 6                              # mountain hemlock (NI variant eqn)
        dds = -1.52111f0 + 0.08518f0*60.0f0 - 0.000943f0*60.0f0*60.0f0 +
              0.13363f0*sin(aspect)*slope + 0.17935f0*cos(aspect)*slope + 0.07628f0*slope +
              0.89778f0*flog(dpp) - 0.10744f0*0.01f0*relden + 1.28403f0*cr - 0.000484f0*dpp*dpp -
              0.66110f0*bal100/flog(dpp+1.0f0)
        idds = 1
    elseif is == 7                              # western redcedar (NI variant eqn)
        dds = 1.61452f0 - 0.00175f0*60.0f0 - 0.000067f0*60.0f0*60.0f0 +
              0.05534f0*sin(aspect)*slope - 0.06625f0*cos(aspect)*slope + 0.11931f0*slope +
              0.58705f0*flog(dpp) - 0.15356f0*0.01f0*relden + 0.74596f0*bal100 + 1.29360f0*cr -
              2.28375f0*bal100/flog(dpp+1.0f0)
        idds = 1
    elseif is == 8                              # western larch (NI variant eqn)
        dds = 0.51291f0 + 0.20004f0 + 0.03730f0*60.0f0 - 0.000433f0*60.0f0*60.0f0 +
              0.03430f0*sin(aspect)*slope - 0.21337f0*cos(aspect)*slope + 0.33523f0*slope -
              0.70216f0*slope*slope - 0.05438f0*0.01f0*relden + 0.54140f0*flog(dpp) +
              0.43637f0*bal100 + 1.03478f0*cr + 0.07509f0*cr*cr - 0.000310f0*dpp*dpp -
              2.03256f0*bal100/flog(dpp+1.0f0)
        idds = 1
    elseif is == 10                             # limber pine (UT variant eqn)
        con1 = sin(aspect-0.7854f0)*slope*(-0.01752f0)
        con2 = cos(aspect-0.7854f0)*slope*(-0.609774f0)
        dds = 1.911884f0 + 0.001766f0*si + con1 + con2 - 2.05706f0*slope +
              2.113263f0*slope*slope - 0.199592f0*0.01f0*relden + 0.213947f0*flog(dpp) -
              0.358634f0*bal100 + 1.523464f0*cr - 0.0006538f0*dpp*dpp
        idds = 1
    elseif is == 11                             # lodgepole pine
        df = 1.32652f0 + 0.96279f0*dpp - 0.17138f0*bgttba - 0.21827f0*flog(batem) +
             0.06298f0*l2(dpp) + 0.00590f0*si
    elseif is == 9 || is == 12 || is == 16 || (23 <= is <= 27) || (29 <= is <= 35) || is == 37
        # bristlecone/pinyons/junipers/oaks/other-softwood
        df = 0.25897f0 + 1.03129f0*dpp - 0.0002025464f0*batem + 0.00177f0*si
        (df - dpp) > 1.0f0 && (df = dpp + 1.0f0)
    elseif is == 13 || is == 36                 # ponderosa pine / Chihuahua pine
        if imodty == 1                          # SWMC: blend MCPP/SWPP
            dfmcpp = 0.26265f0 + 0.94001f0*dpp - 0.0016405f0*bat + 0.0060074f0*si + 0.16328f0*l2(dpp)
            bgttba = bautba
            batem = bat
            batem < 21.0f0 && (batem = 21.0f0)
            dfswpp = 4.10552f0 + 0.88872f0*dpp - 0.83531f0*flog(batem) + 0.013781f0*si -
                     0.52784f0*bgttba + 0.23834f0*fpow(flog(dp), _LN2F)
            adjb = 0.0f0
            if bat > 60.0f0 && bat < 100.0f0
                adjb = -0.025f0*bat + 2.5f0
            elseif bat < 60.1f0
                adjb = 1.0f0
            end
            bapp = spba < 1.0f0 ? 1.0f0 : spba
            barat = bapp / bat
            adjr = 0.0f0
            if barat > 0.3f0 && barat < 0.5f0
                adjr = 5.0f0*barat - 1.5f0
            elseif barat > 0.499f0
                adjr = 1.0f0
            end
            adj = max(adjb, adjr)
            df = (1.0f0 - adj)*dfmcpp + adj*dfswpp
        elseif imodty == 2                      # SWPP
            df = 4.10552f0 + 0.88872f0*dpp - 0.83531f0*flog(batem) + 0.013781f0*si -
                 0.52784f0*bgttba + 0.23834f0*l2(dpp)
        elseif imodty == 3                      # Black Hills PP
            df = 2.50428f0 + 0.94573f0*dpp - 0.45765f0*flog(batem) + 0.01354f0*si -
                 0.46350f0*bgttba + 0.05940f0*l2(dpp)
            if agerng > 40.0f0 && bgttba > 0.3f0
                siminr = 0.57593f0 + 0.92561f0*si
                batem < 21.0f0 && (batem = 21.0f0)
                dfswpp = 4.10552f0 + 0.88872f0*dpp - 0.83531f0*flog(batem) + 0.013781f0*siminr -
                         0.52784f0*bgttba + 0.23834f0*l2(dpp)
                if df > dfswpp
                    if bgttba >= 0.5f0
                        df = dfswpp
                        df < (dpp + 0.1f0) && (df = dpp + 0.1f0)
                    else
                        df = (2.5f0 - 5.0f0*bgttba)*df + (-1.5f0 + 5.0f0*bgttba)*dfswpp
                    end
                end
            end
        else                                    # lodgepole/spruce-fir types
            df = 2.50428f0 + 0.94573f0*dpp - 0.45765f0*flog(batem) + 0.01354f0*si -
                 0.46350f0*bgttba + 0.05940f0*l2(dpp)
        end
    elseif is == 14                             # whitebark pine (EM variant eqn)
        xslope = slope / 10.0f0
        dds = 0.01545f0 + 1.5675f0 - 0.00565f0*80.0f0 - 0.01606f0*sin(aspect)*xslope +
              0.00270f0*cos(aspect)*xslope - 0.20011f0*xslope + 0.80110f0*flog(dpp) +
              0.00064f0*bal + 1.02878f0*cr - 0.45448f0*cr*cr - 0.00328f0*bal/flog(dpp+1.0f0) -
              0.25717f0*flog(relden)
        idds = 1
    elseif is == 15                             # western white pine
        df = 0.89451f0 + 0.95435f0*dpp - 0.0012183f0*bat + 0.00159f0*si -
             0.20815f0*bgttba + 0.09987f0*l2(dpp)
    elseif is == 17                             # blue spruce
        df = 0.50941f0 + 1.02697f0*dpp - 0.0035783f0*bat + 0.01022f0*si - 0.00098825f0*dpp*dpp
    elseif is == 18                             # Engelmann spruce
        if imodty <= 2
            df = 0.67074f0 + 1.03068f0*dpp - 0.0020055f0*bat + 0.00633f0*si - 0.0010875f0*dpp*dpp
        else
            df = 2.28652f0 + 0.94475f0*dpp - 0.48056f0*bgttba - 0.41632f0*flog(batem) +
                 0.10778f0*l2(dpp) + 0.01295f0*si
        end
    elseif is == 19                             # white spruce
        df = 2.94734f0 + 0.88423f0*dpp - 0.64449f0*flog(batem) + 0.01207f0*si + 0.25363f0*l2(dpp)
    elseif is == 20 || is == 21 || is == 22 || is == 28 || is == 38
        # aspen/cottonwoods/paper-birch/other-hardwood
        if imodty <= 2
            df = 0.24506f0 + 1.01291f0*dpp - 0.00084659f0*bat + 0.00631f0*si
        else
            df = 1.55986f0 + 1.01825f0*dpp - 0.29342f0*flog(batem) + 0.00672f0*si - 0.00073f0*bgttba
        end
        df = df * 1.05f0
    else                                        # CASE DEFAULT
        dds = -9.21f0
        idds = 1
    end

    # --- common tail (cr/gemdg.f:355-377) ---
    if idds > 0
        diagr = sqrt(fexp(dds) + fpow(dpp*bark, _LN2F)) - dpp*bark
    else
        df > dbhmax && (df = dbhmax)
        df < dpp && (df = dpp)
        diagr = (df - dpp) * bark
    end
    istagf != 0 && (diagr = diagr * dstag)
    (imodty == 3 && (is == 13 || is == 36)) && (diagr = diagr * 0.80f0)

    if diagr <= 0.0f0
        dds = -9.21f0
    else
        dds = flog(diagr * (2.0f0*dpp*bark + diagr))
        dds < -9.21f0 && (dds = -9.21f0)
    end
    return dds
end
