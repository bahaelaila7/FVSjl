# =============================================================================
# crown_width.jl — open-grown / forest-grown crown width (CWCALC)
#
# Ported from: base/cwcalc.jl (the eastern-US crown-width equation library, 144
# species/equation variants). Used by the crown competition factor (CCF).
#
# Pure: maps the species 2-char code to an equation, applies the crown-width
# formula (Bechtold / Bragg / Ek / Smith / Krajicek forms), clamps to [0.5, 99.9].
# `iwho==1` selects the open-grown form; `hi` is the Hopkins bioclimatic index
# (needs stand lat/long/elevation). Globals TLAT/TLONG/ELEV/NSP are now arguments.
# =============================================================================

"""
    crown_width(sp2, d, h, cr, iwho, lat, long, elev) -> Float32

Crown width (ft) for a tree of species 2-char code `sp2`, DBH `d`, height `h`,
crown ratio `cr`. `iwho==1` → open-grown crown.
"""
function crown_width(sp2::AbstractString, d::Real, h::Real, cr::Real, iwho::Integer,
                     lat::Real, long::Real, elev::Real)::Float32
    df    = Float32(d)
    crf   = Float32(cr)
    mind  = 5.0f0
    omind = 3.0f0
    # Hopkins bioclimatic index (Bechtold 2003)
    hilat  = Float32(lat)
    hilong = -abs(Float32(long))
    hielev = Float32(elev) * Float32(100)
    hi     = ((hielev - Float32(887)) / Float32(100)) * Float32(1) +
             (hilat  - Float32(39.54)) * Float32(4) +
             (Float32(-82.52) - hilong) * Float32(1.25)


    # ---------------------------------------------------------------------------
    # 1. Map species 2-char code → equation number
    # ---------------------------------------------------------------------------
    cweq = ""
    if     sp2 == "AB"; cweq = "53101"
    elseif sp2 == "AC"; cweq = iwho == 1 ? "40703" : "40701"
    elseif sp2 == "AE"; cweq = iwho == 1 ? "97203" : "97201"
    elseif sp2 == "AH"; cweq = "39101"
    elseif sp2 == "AI"; cweq = "49101"
    elseif sp2 == "AP"; cweq = "76102"
    elseif sp2 == "AS"; cweq = iwho == 1 ? "54403" : "54401"
    elseif sp2 == "AW"; cweq = "24101"
    elseif sp2 == "BA"; cweq = "54301"
    elseif sp2 == "BB"; cweq = "37301"
    elseif sp2 == "BC"; cweq = iwho == 1 ? "76203" : "76201"
    elseif sp2 == "BE"; cweq = "31301"
    elseif sp2 == "BF"; cweq = iwho == 1 ? "01203" : "01201"
    elseif sp2 == "BG"; cweq = "69301"
    elseif sp2 == "BH"; cweq = "40201"
    elseif sp2 == "BI"; cweq = "40801"
    elseif sp2 == "BJ"; cweq = "82401"
    elseif sp2 == "BK"; cweq = "90101"
    elseif sp2 == "BL"; cweq = iwho == 1 ? "97203" : "97201"
    elseif sp2 == "BM"; cweq = iwho == 1 ? "31803" : "31801"
    elseif sp2 == "BN"; cweq = "60201"
    elseif sp2 == "BO"; cweq = iwho == 1 ? "83704" : "83701"
    elseif sp2 == "BP"; cweq = "74101"
    elseif sp2 == "BR"; cweq = iwho == 1 ? "82303" : "82301"
    elseif sp2 == "BS"; cweq = iwho == 1 ? "09503" : "09501"
    elseif sp2 == "BT"; cweq = "74301"
    elseif sp2 == "BU"; cweq = iwho == 1 ? "40703" : "40701"
    elseif sp2 == "BW" || sp2 == "BD"; cweq = "95101"
    elseif sp2 == "BY"; cweq = "22101"
    elseif sp2 == "CA"; cweq = "93101"
    elseif sp2 == "CB"; cweq = "81201"
    elseif sp2 == "CC"; cweq = "76102"
    elseif sp2 == "CK"; cweq = "82601"
    elseif sp2 == "CM"; cweq = iwho == 1 ? "31803" : "31801"
    elseif sp2 == "CO"; cweq = "83201"
    elseif sp2 == "CT"; cweq = "65101"
    elseif sp2 == "CW"; cweq = iwho == 1 ? "74203" : "74201"
    elseif sp2 == "DM"; cweq = iwho == 1 ? "97203" : "97201"
    elseif sp2 == "DO"; cweq = "83501"
    elseif sp2 == "DP"; cweq = "83501"
    elseif sp2 == "DW"; cweq = "49101"
    elseif sp2 == "EC"; cweq = iwho == 1 ? "74203" : "74201"
    elseif sp2 == "EH"; cweq = "26101"
    elseif sp2 == "EL"; cweq = iwho == 1 ? "97203" : "97201"
    elseif sp2 == "FM"; cweq = iwho == 1 ? "31803" : "31801"
    elseif sp2 == "FR"; cweq = iwho == 1 ? "01203" : "01201"
    elseif sp2 == "GA"; cweq = iwho == 1 ? "54403" : "54401"
    elseif sp2 == "GB"; cweq = iwho == 1 ? "37503" : "37501"
    elseif sp2 == "HA"; cweq = "49101"
    elseif sp2 == "HB"; cweq = "46201"
    elseif sp2 == "HH"; cweq = "70101"
    elseif sp2 == "HI"; cweq = iwho == 1 ? "40703" : "40701"
    elseif sp2 == "HK"; cweq = "46201"
    elseif sp2 == "HL"; cweq = "55201"
    elseif sp2 == "HM"; cweq = "26101"
    elseif sp2 == "HT"; cweq = "49101"
    elseif sp2 == "HY"; cweq = "59101"
    elseif sp2 == "JP"; cweq = iwho == 1 ? "10503" : "10501"
    elseif sp2 == "JU"; cweq = "06801"
    elseif sp2 == "KC"; cweq = "90101"
    elseif sp2 == "LB"; cweq = "65301"
    elseif sp2 == "LK"; cweq = "82001"
    elseif sp2 == "LL"; cweq = iwho == 1 ? "12105" : "12101"
    elseif sp2 == "LO"; cweq = "83801"
    elseif sp2 == "LP"; cweq = iwho == 1 ? "13105" : "13101"
    elseif sp2 == "MA"; cweq = "55201"
    elseif sp2 == "MB"; cweq = "68201"
    elseif sp2 == "MG"; cweq = "65301"
    elseif sp2 == "MH"; cweq = "40901"
    elseif sp2 == "ML"; cweq = "65301"
    elseif sp2 == "MM"; cweq = "31301"
    elseif sp2 == "MP"; cweq = iwho == 1 ? "31803" : "31801"
    elseif sp2 == "MS"; cweq = "65301"
    elseif sp2 == "MV"; cweq = "65301"
    elseif sp2 == "NK"; cweq = "81201"
    elseif sp2 == "NP"; cweq = "80901"
    elseif sp2 == "NS"; cweq = iwho == 1 ? "09104" : "09101"
    elseif sp2 == "OB"; cweq = iwho == 1 ? "40703" : "40701"
    elseif sp2 == "OC"; cweq = "06801"
    elseif sp2 == "OG"; cweq = "69101"
    elseif sp2 == "OH"; cweq = "93101"
    elseif sp2 == "OK"; cweq = iwho == 1 ? "80204" : "80201"
    elseif sp2 == "OO"; cweq = "93101"
    elseif sp2 == "OP"; cweq = iwho == 1 ? "12903" : "12901"
    elseif sp2 == "OS"; cweq = "06801"
    elseif sp2 == "OT"; cweq = iwho == 1 ? "31603" : "31601"
    elseif sp2 == "OV"; cweq = iwho == 1 ? "82303" : "82301"
    elseif sp2 == "PA"; cweq = "54101"
    elseif sp2 == "PB"; cweq = iwho == 1 ? "37503" : "37501"
    elseif sp2 == "PC"; cweq = "22101"
    elseif sp2 == "PD"; cweq = "12801"
    elseif sp2 == "PE"; cweq = "40201"
    elseif sp2 == "PH"; cweq = "40301"
    elseif sp2 == "PI"; cweq = iwho == 1 ? "09403" : "09401"
    elseif sp2 == "PL"; cweq = "76102"
    elseif sp2 == "PN"; cweq = "83001"
    elseif sp2 == "PO"; cweq = "83501"
    elseif sp2 == "PP"; cweq = "12601"
    elseif sp2 == "PR"; cweq = "76102"
    elseif sp2 == "PS"; cweq = "52101"
    elseif sp2 == "PU"; cweq = "13201"
    elseif sp2 == "PW"; cweq = "93101"
    elseif sp2 == "PY"; cweq = iwho == 1 ? "74203" : "74201"
    elseif sp2 == "PZ"; cweq = iwho == 1 ? "13105" : "13101"
    elseif sp2 == "QA"; cweq = iwho == 1 ? "74603" : "74601"
    elseif sp2 == "QI"; cweq = "81701"
    elseif sp2 == "QN"; cweq = "81901"
    elseif sp2 == "QS"; cweq = "81201"
    elseif sp2 == "RA"; cweq = "72101"
    elseif sp2 == "RB"; cweq = "37301"
    elseif sp2 == "RC"; cweq = "06801"
    elseif sp2 == "RD"; cweq = "49101"
    elseif sp2 == "RE"; cweq = iwho == 1 ? "97203" : "97201"
    elseif sp2 == "RL"; cweq = "97501"
    elseif sp2 == "RM"; cweq = iwho == 1 ? "31603" : "31601"
    elseif sp2 == "RN"; cweq = iwho == 1 ? "12503" : "12501"
    elseif sp2 == "RO"; cweq = iwho == 1 ? "83303" : "83301"
    elseif sp2 == "RP"; cweq = iwho == 1 ? "12503" : "12501"
    elseif sp2 == "RS"; cweq = "09701"
    elseif sp2 == "RY"; cweq = "68201"
    elseif sp2 == "SA"; cweq = "11101"
    elseif sp2 == "SB"; cweq = "37201"
    elseif sp2 == "SC"; cweq = "13001"
    elseif sp2 == "SD"; cweq = "71101"
    elseif sp2 == "SE"; cweq = "35601"
    elseif sp2 == "SG"; cweq = "46201"
    elseif sp2 == "SH"; cweq = iwho == 1 ? "40703" : "40701"
    elseif sp2 == "SI"; cweq = iwho == 1 ? "97203" : "97201"
    elseif sp2 == "SK"; cweq = "81201"
    elseif sp2 == "SL"; cweq = iwho == 1 ? "40703" : "40701"
    elseif sp2 == "SM"; cweq = iwho == 1 ? "31803" : "31801"
    elseif sp2 == "SN"; cweq = "83201"
    elseif sp2 == "SO"; cweq = "80601"
    elseif sp2 == "SP"; cweq = iwho == 1 ? "11005" : "11001"
    elseif sp2 == "SR"; cweq = iwho == 1 ? "11005" : "11001"
    elseif sp2 == "SS"; cweq = "93101"
    elseif sp2 == "ST"; cweq = "31301"
    elseif sp2 == "SU"; cweq = "61101"
    elseif sp2 == "SV"; cweq = "31701"
    elseif sp2 == "SW"; cweq = iwho == 1 ? "80204" : "80201"
    elseif sp2 == "SY"; cweq = "73101"
    elseif sp2 == "TA"; cweq = iwho == 1 ? "07103" : "07101"
    elseif sp2 == "TL"; cweq = "69301"
    elseif sp2 == "TM"; cweq = "12601"
    elseif sp2 == "TO"; cweq = "81901"
    elseif sp2 == "TS"; cweq = "69401"
    elseif sp2 == "UA"; cweq = "54101"
    elseif sp2 == "VP"; cweq = "13201"
    elseif sp2 == "WA"; cweq = "54101"
    elseif sp2 == "WB"; cweq = "95101"
    elseif sp2 == "WC"; cweq = "24101"
    elseif sp2 == "WE"; cweq = "97101"
    elseif sp2 == "WH"; cweq = "40201"
    elseif sp2 == "WI"; cweq = iwho == 1 ? "97203" : "97201"
    elseif sp2 == "WK"; cweq = "82701"
    elseif sp2 == "WL"; cweq = "83101"
    elseif sp2 == "WM"; cweq = "68201"
    elseif sp2 == "WN"; cweq = "60201"
    elseif sp2 == "WO"; cweq = iwho == 1 ? "80204" : "80201"
    elseif sp2 == "WP"; cweq = iwho == 1 ? "12903" : "12901"
    elseif sp2 == "WR"; cweq = "37301"
    elseif sp2 == "WS"; cweq = iwho == 1 ? "09403" : "09401"
    elseif sp2 == "WT"; cweq = "69101"
    elseif sp2 == "YB"; cweq = "37101"
    elseif sp2 == "YP"; cweq = "62101"
    elseif sp2 == "YY"; cweq = iwho == 1 ? "40703" : "40701"
    end

    # ---------------------------------------------------------------------------
    # 2. Apply crown width equation
    # ---------------------------------------------------------------------------
    cw_val = Float32(0)

    # --- SOFTWOODS ---

    if cweq == "01201"  # Balsam fir — Bechtold Model 2
        cw_val = df >= mind ? Float32(0.6564) + Float32(0.8403)*df + Float32(0.0792)*crf :
                              (Float32(0.6564) + Float32(0.8403)*mind + Float32(0.0792)*crf)*(df/mind)
        if cw_val > Float32(34); cw_val = Float32(34); end
    elseif cweq == "01202"  # Balsam fir — Bragg
        cw_val = Float32(-3.746931) + Float32(7.122778)*df^Float32(0.396998)
    elseif cweq == "01203"  # Balsam fir — Ek
        cw_val = df >= omind ? Float32(0.3270) + Float32(5.1160)*df^Float32(0.5035) :
                               (Float32(0.3270) + Float32(5.1160)*omind^Float32(0.5035))*(df/omind)
        if cw_val > Float32(34); cw_val = Float32(34); end

    elseif cweq == "06801"  # Eastern redcedar — Bechtold Model 2
        cw_val = df >= mind ? Float32(1.2359) + Float32(1.2962)*df + Float32(0.0545)*crf :
                              (Float32(1.2359) + Float32(1.2962)*mind + Float32(0.0545)*crf)*(df/mind)
        if cw_val > Float32(33); cw_val = Float32(33); end

    elseif cweq == "07101"  # Tamarack — Bechtold Model 2
        cw_val = df >= mind ? Float32(-0.3276) + Float32(1.3865)*df + Float32(0.0517)*crf :
                              (Float32(-0.3276) + Float32(1.3865)*mind + Float32(0.0517)*crf)*(df/mind)
        if cw_val > Float32(29); cw_val = Float32(29); end
    elseif cweq == "07102"  # Tamarack — Bragg
        cw_val = Float32(2.503585) + Float32(1.100883)*df^Float32(1.056165)
    elseif cweq == "07103"  # Tamarack — Ek
        cw_val = df >= omind ? Float32(2.205) + Float32(3.475)*df^Float32(0.7506) :
                               (Float32(2.205) + Float32(3.475)*omind^Float32(0.7506))*(df/omind)
        if cw_val > Float32(29); cw_val = Float32(29); end

    elseif cweq == "09101"  # Norway spruce — Bechtold Model 3
        cw_val = df >= mind ? Float32(1.8336) + Float32(0.9932)*df + Float32(0.0431)*crf + Float32(0.1012)*hi :
                              (Float32(1.8336) + Float32(0.9932)*mind + Float32(0.0431)*crf + Float32(0.1012)*hi)*(df/mind)
        if cw_val > Float32(27); cw_val = Float32(27); end
    elseif cweq == "09104"  # Norway spruce — Krajicek
        cw_val = df >= omind ? Float32(5.0570) + Float32(1.1313)*df :
                               (Float32(5.0570) + Float32(1.1313)*omind)*(df/omind)
        if cw_val > Float32(47); cw_val = Float32(47); end

    elseif cweq == "09401"  # White spruce — Bechtold Model 2
        cw_val = df >= mind ? Float32(0.3789) + Float32(0.8658)*df + Float32(0.0878)*crf :
                              (Float32(0.3789) + Float32(0.8658)*mind + Float32(0.0878)*crf)*(df/mind)
        if cw_val > Float32(30); cw_val = Float32(30); end
    elseif cweq == "09402"  # White spruce — Bragg
        cw_val = Float32(3.067563) + Float32(1.944947)*df^Float32(0.718583)
    elseif cweq == "09403"  # White spruce — Ek
        cw_val = df >= omind ? Float32(3.5940) + Float32(1.9630)*df^Float32(0.8820) :
                               (Float32(3.5940) + Float32(1.9630)*omind^Float32(0.8820))*(df/omind)
        if cw_val > Float32(37); cw_val = Float32(37); end

    elseif cweq == "09501"  # Black spruce — Bechtold Model 2
        cw_val = df >= mind ? Float32(-0.8566) + Float32(0.9693)*df + Float32(0.0573)*crf :
                              (Float32(-0.8566) + Float32(0.9693)*mind + Float32(0.0573)*crf)*(df/mind)
        if cw_val > Float32(27); cw_val = Float32(27); end
    elseif cweq == "09502"  # Black spruce — Bragg
        cw_val = Float32(4.281343) + Float32(0.153325)*df^Float32(1.787982)
    elseif cweq == "09503"  # Black spruce — Ek
        cw_val = df >= omind ? Float32(3.6550) + Float32(1.398)*df^Float32(1.000) :
                               (Float32(3.6550) + Float32(1.398)*omind^Float32(1.000))*(df/omind)
        if cw_val > Float32(27); cw_val = Float32(27); end

    elseif cweq == "09701"  # Red spruce — Bechtold Model 3
        if df >= mind
            dlim = min(df, Float32(30))
            cw_val = Float32(-1.2151) + Float32(1.6098)*dlim - Float32(0.0277)*dlim*dlim +
                     Float32(0.0674)*crf - Float32(0.0474)*hi
        else
            cw_val = (Float32(-1.2151) + Float32(1.6098)*mind - Float32(0.0277)*mind*mind +
                      Float32(0.0674)*crf - Float32(0.0474)*hi) * (df/mind)
        end

    elseif cweq == "10501"  # Jack pine — Bechtold Model 2
        cw_val = df >= mind ? Float32(0.7478) + Float32(0.8712)*df + Float32(0.0913)*crf :
                              (Float32(0.7478) + Float32(0.8712)*mind + Float32(0.0913)*crf)*(df/mind)
        if cw_val > Float32(25); cw_val = Float32(25); end
    elseif cweq == "10502"  # Jack pine — Bragg
        cw_val = Float32(3.382473) + Float32(0.529126)*df^Float32(1.269473)
    elseif cweq == "10503"  # Jack pine — Ek
        cw_val = df >= omind ? Float32(0.2990) + Float32(5.6440)*df^Float32(0.6036) :
                               (Float32(0.2990) + Float32(5.6440)*omind^Float32(0.6036))*(df/omind)
        if cw_val > Float32(30); cw_val = Float32(30); end

    elseif cweq == "11001"  # Shortleaf pine — Bechtold Model 3
        cw_val = df >= mind ? Float32(-2.2564) + Float32(1.3004)*df + Float32(0.1031)*crf - Float32(0.0562)*hi :
                              (Float32(-2.2564) + Float32(1.3004)*mind + Float32(0.1031)*crf - Float32(0.0562)*hi)*(df/mind)
        if cw_val > Float32(34); cw_val = Float32(34); end
    elseif cweq == "11005"  # Shortleaf pine — Smith 1992
        if df >= omind
            dcm = df * Float32(2.54)
            cw_val = (Float32(0.5830) + Float32(0.2450)*dcm + Float32(0.0009)*dcm*dcm) * Float32(3.28084)
        else
            dcm = omind * Float32(2.54)
            cw_val = (Float32(0.5830) + Float32(0.2450)*dcm + Float32(0.0009)*dcm*dcm) * Float32(3.28084) * (df/omind)
        end
        if cw_val > Float32(45); cw_val = Float32(45); end

    elseif cweq == "11101"  # Slash pine — Bechtold Model 3
        if df >= mind
            dlim = min(df, Float32(30))
            cw_val = Float32(-6.9659) + Float32(2.1192)*dlim - Float32(0.0333)*dlim*dlim +
                     Float32(0.0587)*crf - Float32(0.0959)*hi
        else
            cw_val = (Float32(-6.9659) + Float32(2.1192)*mind - Float32(0.0333)*mind*mind +
                      Float32(0.0587)*crf - Float32(0.0959)*hi) * (df/mind)
        end

    elseif cweq == "12101"  # Longleaf pine — Bechtold Model 3
        cw_val = df >= mind ? Float32(-12.2105) + Float32(1.3376)*df + Float32(0.1237)*crf - Float32(0.2759)*hi :
                              (Float32(-12.2105) + Float32(1.3376)*mind + Float32(0.1237)*crf - Float32(0.2759)*hi)*(df/mind)
        if cw_val > Float32(50); cw_val = Float32(50); end
    elseif cweq == "12105"  # Longleaf pine — Smith 1992
        if df >= omind
            cw_val = (Float32(0.113) + Float32(0.259)*(df*Float32(2.54))) * Float32(3.28084)
        else
            cw_val = (Float32(0.113) + Float32(0.259)*(omind*Float32(2.54))) * Float32(3.28084) * (df/omind)
        end
        if cw_val > Float32(50); cw_val = Float32(50); end

    elseif cweq == "12501"  # Red pine — Bechtold Model 2
        if df >= mind
            dlim = min(df, Float32(24))
            cw_val = Float32(-3.6548) + Float32(1.9565)*dlim - Float32(0.0409)*dlim*dlim + Float32(0.0577)*crf
        else
            cw_val = (Float32(-3.6548) + Float32(1.9565)*mind - Float32(0.0409)*mind*mind + Float32(0.0577)*crf)*(df/mind)
        end
    elseif cweq == "12502"  # Red pine — Bragg
        cw_val = Float32(3.499341) + Float32(0.806186)*df^Float32(1.090937)
    elseif cweq == "12503"  # Red pine — Ek
        cw_val = df >= omind ? Float32(4.2330) + Float32(1.4620)*df^Float32(1.0000) :
                               (Float32(4.2330) + Float32(1.4620)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(39); cw_val = Float32(39); end

    elseif cweq == "12601"  # Pitch pine — Bechtold Model 3
        cw_val = df >= mind ? Float32(-0.9442) + Float32(1.4531)*df + Float32(0.0543)*crf - Float32(0.1144)*hi :
                              (Float32(-0.9442) + Float32(1.4531)*mind + Float32(0.0543)*crf - Float32(0.1144)*hi)*(df/mind)
        if cw_val > Float32(34); cw_val = Float32(34); end

    elseif cweq == "12801"  # Pond pine — Bechtold Model 1
        if df >= mind
            dlim = min(df, Float32(18))
            cw_val = Float32(-8.7711) + Float32(3.7252)*dlim - Float32(0.1063)*dlim*dlim
        else
            cw_val = (Float32(-8.7711) + Float32(3.7252)*mind - Float32(0.1063)*mind*mind)*(df/mind)
        end

    elseif cweq == "12901"  # Eastern white pine — Bechtold Model 2
        cw_val = df >= mind ? Float32(0.3914) + Float32(0.9923)*df + Float32(0.1080)*crf :
                              (Float32(0.3914) + Float32(0.9923)*mind + Float32(0.1080)*crf)*(df/mind)
        if cw_val > Float32(45); cw_val = Float32(45); end
    elseif cweq == "12902"  # Eastern white pine — Bragg
        cw_val = Float32(3.874199) + Float32(1.062309)*df^Float32(0.969580)
    elseif cweq == "12903"  # Eastern white pine — Ek
        cw_val = df >= omind ? Float32(1.6200) + Float32(3.1970)*df^Float32(0.7981) :
                               (Float32(1.6200) + Float32(3.1970)*omind^Float32(0.7981))*(df/omind)
        if cw_val > Float32(58); cw_val = Float32(58); end

    elseif cweq == "13001"  # Scotch pine — Bechtold Model 2
        cw_val = df >= mind ? Float32(3.5522) + Float32(0.6742)*df + Float32(0.0985)*crf :
                              (Float32(3.5522) + Float32(0.6742)*mind + Float32(0.0985)*crf)*(df/mind)
        if cw_val > Float32(27); cw_val = Float32(27); end

    elseif cweq == "13101"  # Loblolly pine — Bechtold Model 2
        cw_val = df >= mind ? Float32(-0.8277) + Float32(1.3946)*df + Float32(0.0768)*crf :
                              (Float32(-0.8277) + Float32(1.3946)*mind + Float32(0.0768)*crf)*(df/mind)
        if cw_val > Float32(55); cw_val = Float32(55); end
    elseif cweq == "13105"  # Loblolly pine — Smith 1992
        if df >= omind
            dcm = df * Float32(2.54)
            cw_val = (Float32(0.7380) + Float32(0.2450)*dcm + Float32(0.000809)*dcm*dcm) * Float32(3.28084)
        else
            dcm = omind * Float32(2.54)
            cw_val = (Float32(0.7380) + Float32(0.2450)*dcm + Float32(0.000809)*dcm*dcm) * Float32(3.28084) * (df/omind)
        end
        if cw_val > Float32(66); cw_val = Float32(66); end

    elseif cweq == "13201"  # Virginia pine — Bechtold Model 2
        cw_val = df >= mind ? Float32(-0.1211) + Float32(1.2319)*df + Float32(0.1212)*crf :
                              (Float32(-0.1211) + Float32(1.2319)*mind + Float32(0.1212)*crf)*(df/mind)
        if cw_val > Float32(34); cw_val = Float32(34); end

    elseif cweq == "22101"  # Baldcypress — Bechtold Model 2
        cw_val = df >= mind ? Float32(-1.0183) + Float32(0.8856)*df + Float32(0.1162)*crf :
                              (Float32(-1.0183) + Float32(0.8856)*mind + Float32(0.1162)*crf)*(df/mind)
        if cw_val > Float32(37); cw_val = Float32(37); end

    elseif cweq == "24101"  # Northern white-cedar — Bechtold Model 2
        cw_val = df >= mind ? Float32(-0.0634) + Float32(0.7057)*df + Float32(0.0837)*crf :
                              (Float32(-0.0634) + Float32(0.7057)*mind + Float32(0.0837)*crf)*(df/mind)
        if cw_val > Float32(27); cw_val = Float32(27); end
    elseif cweq == "24102"  # Northern white-cedar — Bragg
        cw_val = Float32(2.123722) + Float32(1.898797)*df^Float32(0.764193)
    elseif cweq == "24103"  # Northern white-cedar — Ek
        cw_val = df >= omind ? Float32(3.4960) + Float32(1.0930)*df^Float32(1.0000) :
                               (Float32(3.4960) + Float32(1.0930)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(27); cw_val = Float32(27); end

    elseif cweq == "26101"  # Eastern hemlock — Bechtold Model 3
        if df >= mind
            dlim = min(df, Float32(40))
            cw_val = Float32(6.1924) + Float32(1.4491)*dlim - Float32(0.0178)*dlim*dlim - Float32(0.0341)*hi
        else
            cw_val = (Float32(6.1924) + Float32(1.4491)*mind - Float32(0.0178)*mind*mind - Float32(0.0341)*hi)*(df/mind)
        end
    elseif cweq == "26102"  # Eastern hemlock — Bragg
        cw_val = Float32(0.868672) + Float32(4.526525)*df^Float32(0.589487)
    elseif cweq == "26103"  # Eastern hemlock — Ek
        cw_val = df >= omind ? Float32(0.5230) + Float32(1.6320)*df^Float32(1.0000) :
                               (Float32(0.5230) + Float32(1.6320)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(39); cw_val = Float32(39); end

    # --- HARDWOODS ---

    elseif cweq == "31301"  # Boxelder — Bechtold Model 3
        cw_val = df >= mind ? Float32(6.4741) + Float32(1.0778)*df + Float32(0.0719)*crf - Float32(0.0637)*hi :
                              (Float32(6.4741) + Float32(1.0778)*mind + Float32(0.0719)*crf - Float32(0.0637)*hi)*(df/mind)
        if cw_val > Float32(57); cw_val = Float32(57); end

    elseif cweq == "31601"  # Red maple — Bechtold Model 3
        if df >= mind
            dlim = min(df, Float32(50))
            cw_val = Float32(2.7563) + Float32(1.4212)*dlim - Float32(0.0143)*dlim*dlim +
                     Float32(0.0993)*crf - Float32(0.0276)*hi
        else
            cw_val = (Float32(2.7563) + Float32(1.4212)*mind - Float32(0.0143)*mind*mind +
                      Float32(0.0993)*crf - Float32(0.0276)*hi)*(df/mind)
        end
    elseif cweq == "31602"  # Red maple — Bragg
        cw_val = Float32(5.394872) + Float32(1.844592)*df^Float32(0.875755)
    elseif cweq == "31603"  # Red maple — Ek
        cw_val = df >= omind ? Float32(0.0) + Float32(4.7760)*df^Float32(0.7656) :
                               (Float32(0.0) + Float32(4.7760)*omind^Float32(0.7656))*(df/omind)
        if cw_val > Float32(55); cw_val = Float32(55); end

    elseif cweq == "31701"  # Silver maple — Bechtold Model 3
        cw_val = df >= mind ? Float32(3.3576) + Float32(1.1312)*df + Float32(0.1011)*crf - Float32(0.1730)*hi :
                              (Float32(3.3576) + Float32(1.1312)*mind + Float32(0.1011)*crf - Float32(0.1730)*hi)*(df/mind)
        if cw_val > Float32(45); cw_val = Float32(45); end

    elseif cweq == "31801"  # Sugar maple — Bechtold Model 3
        cw_val = df >= mind ? Float32(4.9399) + Float32(1.0727)*df + Float32(0.1096)*crf - Float32(0.0493)*hi :
                              (Float32(4.9399) + Float32(1.0727)*mind + Float32(0.1096)*crf - Float32(0.0493)*hi)*(df/mind)
        if cw_val > Float32(54); cw_val = Float32(54); end
    elseif cweq == "31802"  # Sugar maple — Bragg
        cw_val = Float32(5.588697) + Float32(2.302860)*df^Float32(0.828795)
    elseif cweq == "31803"  # Sugar maple — Ek
        cw_val = df >= omind ? Float32(0.8680) + Float32(4.1500)*df^Float32(0.7514) :
                               (Float32(0.8680) + Float32(4.1500)*omind^Float32(0.7514))*(df/omind)
        if cw_val > Float32(54); cw_val = Float32(54); end

    elseif cweq == "35601"  # Serviceberry — Bechtold Model 1
        cw_val = df >= mind ? Float32(6.9814) + Float32(1.6032)*df :
                              (Float32(6.9814) + Float32(1.6032)*mind)*(df/mind)
        if cw_val > Float32(27); cw_val = Float32(27); end

    elseif cweq == "37101"  # Yellow birch — Bechtold Model 3
        if df >= mind
            dlim = min(df, Float32(24))
            cw_val = Float32(-1.1151) + Float32(2.2888)*dlim - Float32(0.0493)*dlim*dlim +
                     Float32(0.0985)*crf - Float32(0.0396)*hi
        else
            cw_val = (Float32(-1.1151) + Float32(2.2888)*mind - Float32(0.0493)*mind*mind +
                      Float32(0.0985)*crf - Float32(0.0396)*hi)*(df/mind)
        end
    elseif cweq == "37102"  # Yellow birch — Bragg
        cw_val = Float32(2.374661) + Float32(4.110366)*df^Float32(0.677280)

    elseif cweq == "37201"  # Sweet birch — Bechtold Model 2
        cw_val = df >= mind ? Float32(4.6725) + Float32(1.2968)*df + Float32(0.0787)*crf :
                              (Float32(4.6725) + Float32(1.2968)*mind + Float32(0.0787)*crf)*(df/mind)
        if cw_val > Float32(54); cw_val = Float32(54); end

    elseif cweq == "37301"  # River birch — Bechtold Model 1
        cw_val = df >= mind ? Float32(11.6634) + Float32(1.0028)*df :
                              (Float32(11.6634) + Float32(1.0028)*mind)*(df/mind)
        if cw_val > Float32(68); cw_val = Float32(68); end

    elseif cweq == "37501"  # Paper birch — Bechtold Model 3
        cw_val = df >= mind ? Float32(2.8399) + Float32(1.2398)*df + Float32(0.0855)*crf - Float32(0.0282)*hi :
                              (Float32(2.8399) + Float32(1.2398)*mind + Float32(0.0855)*crf - Float32(0.0282)*hi)*(df/mind)
        if cw_val > Float32(42); cw_val = Float32(42); end
    elseif cweq == "37502"  # Paper birch — Bragg
        cw_val = Float32(6.342475) + Float32(0.552092)*df^Float32(1.325344)
    elseif cweq == "37503"  # Paper birch — Ek
        cw_val = df >= omind ? Float32(3.6390) + Float32(1.9530)*df^Float32(1.0000) :
                               (Float32(3.6390) + Float32(1.9530)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(42); cw_val = Float32(42); end

    elseif cweq == "39101"  # American hornbeam — Bechtold Model 3
        cw_val = df >= mind ? Float32(0.9219) + Float32(1.6303)*df + Float32(0.1150)*crf - Float32(0.1113)*hi :
                              (Float32(0.9219) + Float32(1.6303)*mind + Float32(0.1150)*crf - Float32(0.1113)*hi)*(df/mind)
        if cw_val > Float32(42); cw_val = Float32(42); end

    elseif cweq == "40201"  # Bitternut hickory — Bechtold Model 1
        cw_val = df >= mind ? Float32(8.0118) + Float32(1.4212)*df :
                              (Float32(8.0118) + Float32(1.4212)*mind)*(df/mind)
        if cw_val > Float32(41); cw_val = Float32(41); end

    elseif cweq == "40301"  # Pignut hickory — Bechtold Model 2
        cw_val = df >= mind ? Float32(3.9234) + Float32(1.5220)*df + Float32(0.0405)*crf :
                              (Float32(3.9234) + Float32(1.5220)*mind + Float32(0.0405)*crf)*(df/mind)
        if cw_val > Float32(53); cw_val = Float32(53); end

    elseif cweq == "40701"  # Shagbark hickory — Bechtold Model 2
        cw_val = df >= mind ? Float32(4.5453) + Float32(1.3721)*df + Float32(0.0430)*crf :
                              (Float32(4.5453) + Float32(1.3721)*mind + Float32(0.0430)*crf)*(df/mind)
        if cw_val > Float32(54); cw_val = Float32(54); end
    elseif cweq == "40703"  # Shagbark hickory — Ek
        cw_val = df >= omind ? Float32(2.3600) + Float32(3.5480)*df^Float32(0.7986) :
                               (Float32(2.3600) + Float32(3.5480)*omind^Float32(0.7986))*(df/omind)
        if cw_val > Float32(54); cw_val = Float32(54); end
    elseif cweq == "40704"  # Shagbark hickory — Krajicek
        cw_val = df >= omind ? Float32(1.9310) + Float32(1.9990)*df :
                               (Float32(1.9310) + Float32(1.9990)*omind)*(df/omind)
        if cw_val > Float32(54); cw_val = Float32(54); end

    elseif cweq == "40801"  # Black hickory — Bechtold Model 1
        if df >= mind
            dlim = min(df, Float32(15))
            cw_val = Float32(-5.8749) + Float32(4.1555)*dlim - Float32(0.1343)*dlim*dlim
        else
            cw_val = (Float32(-5.8749) + Float32(4.1555)*mind - Float32(0.1343)*mind*mind)*(df/mind)
        end

    elseif cweq == "40901"  # Mockernut hickory — Bechtold Model 2
        cw_val = df >= mind ? Float32(1.5838) + Float32(1.6318)*df + Float32(0.0721)*crf :
                              (Float32(1.5838) + Float32(1.6318)*mind + Float32(0.0721)*crf)*(df/mind)
        if cw_val > Float32(55); cw_val = Float32(55); end

    elseif cweq == "46201"  # Hackberry — Bechtold Model 2
        cw_val = df >= mind ? Float32(7.1043) + Float32(1.3041)*df + Float32(0.0456)*crf :
                              (Float32(7.1043) + Float32(1.3041)*mind + Float32(0.0456)*crf)*(df/mind)
        if cw_val > Float32(51); cw_val = Float32(51); end

    elseif cweq == "49101"  # Flowering dogwood — Bechtold Model 2
        cw_val = df >= mind ? Float32(2.9646) + Float32(1.9917)*df + Float32(0.0707)*crf :
                              (Float32(2.9646) + Float32(1.9917)*mind + Float32(0.0707)*crf)*(df/mind)
        if cw_val > Float32(36); cw_val = Float32(36); end

    elseif cweq == "52101"  # Common persimmon — Bechtold Model 2
        cw_val = df >= mind ? Float32(3.5393) + Float32(1.3939)*df + Float32(0.0625)*crf :
                              (Float32(3.5393) + Float32(1.3939)*mind + Float32(0.0625)*crf)*(df/mind)
        if cw_val > Float32(36); cw_val = Float32(36); end

    elseif cweq == "53101"  # American beech — Bechtold Model 3
        cw_val = df >= mind ? Float32(3.9361) + Float32(1.1500)*df + Float32(0.1237)*crf - Float32(0.0691)*hi :
                              (Float32(3.9361) + Float32(1.1500)*mind + Float32(0.1237)*crf - Float32(0.0691)*hi)*(df/mind)
        if cw_val > Float32(80); cw_val = Float32(80); end

    elseif cweq == "54101"  # White ash — Bechtold Model 2
        cw_val = df >= mind ? Float32(1.7625) + Float32(1.3413)*df + Float32(0.0957)*crf :
                              (Float32(1.7625) + Float32(1.3413)*mind + Float32(0.0957)*crf)*(df/mind)
        if cw_val > Float32(62); cw_val = Float32(62); end
    elseif cweq == "54102"  # White ash — Bragg
        cw_val = Float32(5.715288) + Float32(0.914942)*df^Float32(1.113606)
    elseif cweq == "54103"  # White ash — Ek
        cw_val = df >= omind ? Float32(2.3260) + Float32(2.8390)*df^Float32(1.0000) :
                               (Float32(2.3260) + Float32(2.8390)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(62); cw_val = Float32(62); end

    elseif cweq == "54301"  # Black ash — Bechtold Model 1
        cw_val = df >= mind ? Float32(5.2824) + Float32(1.1184)*df :
                              (Float32(5.2824) + Float32(1.1184)*mind)*(df/mind)
        if cw_val > Float32(34); cw_val = Float32(34); end
    elseif cweq == "54302"  # Black ash — Bragg
        cw_val = Float32(2.761995) + Float32(2.560977)*df^Float32(0.742525)

    elseif cweq == "54401"  # Green ash — Bechtold Model 2
        cw_val = df >= mind ? Float32(2.9672) + Float32(1.3066)*df + Float32(0.0585)*crf :
                              (Float32(2.9672) + Float32(1.3066)*mind + Float32(0.0585)*crf)*(df/mind)
        if cw_val > Float32(61); cw_val = Float32(61); end
    elseif cweq == "54403"  # Green ash — Ek
        cw_val = df >= omind ? Float32(0.0000) + Float32(4.7550)*df^Float32(0.7381) :
                               (Float32(0.0000) + Float32(4.7550)*omind^Float32(0.7381))*(df/omind)
        if cw_val > Float32(61); cw_val = Float32(61); end

    elseif cweq == "55201"  # Honeylocust — Bechtold Model 2
        cw_val = df >= mind ? Float32(4.1971) + Float32(1.5567)*df + Float32(0.0880)*crf :
                              (Float32(4.1971) + Float32(1.5567)*mind + Float32(0.0880)*crf)*(df/mind)
        if cw_val > Float32(46); cw_val = Float32(46); end

    elseif cweq == "59101"  # American holly — Bechtold Model 2
        cw_val = df >= mind ? Float32(4.5803) + Float32(1.0747)*df + Float32(0.0661)*crf :
                              (Float32(4.5803) + Float32(1.0747)*mind + Float32(0.0661)*crf)*(df/mind)
        if cw_val > Float32(31); cw_val = Float32(31); end

    elseif cweq == "60201"  # Black walnut — Bechtold Model 2
        cw_val = df >= mind ? Float32(3.6031) + Float32(1.1472)*df + Float32(0.1224)*crf :
                              (Float32(3.6031) + Float32(1.1472)*mind + Float32(0.1224)*crf)*(df/mind)
        if cw_val > Float32(37); cw_val = Float32(37); end
    elseif cweq == "60203"  # Black walnut — Ek
        cw_val = df >= omind ? Float32(4.901) + Float32(2.480)*df^Float32(1.0000) :
                               (Float32(4.901) + Float32(2.480)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(37); cw_val = Float32(37); end

    elseif cweq == "61101"  # Sweetgum — Bechtold Model 3
        cw_val = df >= mind ? Float32(1.8853) + Float32(1.1625)*df + Float32(0.0656)*crf - Float32(0.0300)*hi :
                              (Float32(1.8853) + Float32(1.1625)*mind + Float32(0.0656)*crf - Float32(0.0300)*hi)*(df/mind)
        if cw_val > Float32(50); cw_val = Float32(50); end

    elseif cweq == "62101"  # Yellow-poplar — Bechtold Model 2
        cw_val = df >= mind ? Float32(3.3543) + Float32(1.1627)*df + Float32(0.0857)*crf :
                              (Float32(3.3543) + Float32(1.1627)*mind + Float32(0.0857)*crf)*(df/mind)
        if cw_val > Float32(61); cw_val = Float32(61); end

    elseif cweq == "65101"  # Cucumbertree — Bechtold Model 1
        cw_val = df >= mind ? Float32(4.1711) + Float32(1.6275)*df :
                              (Float32(4.1711) + Float32(1.6275)*mind)*(df/mind)
        if cw_val > Float32(39); cw_val = Float32(39); end

    elseif cweq == "65301"  # Sweetbay — Bechtold Model 1
        cw_val = df >= mind ? Float32(8.2119) + Float32(0.9708)*df :
                              (Float32(8.2119) + Float32(0.9708)*mind)*(df/mind)
        if cw_val > Float32(41); cw_val = Float32(41); end

    elseif cweq == "68201"  # Red mulberry — Bechtold Model 1
        cw_val = df >= mind ? Float32(13.3255) + Float32(1.0735)*df :
                              (Float32(13.3255) + Float32(1.0735)*mind)*(df/mind)
        if cw_val > Float32(46); cw_val = Float32(46); end

    elseif cweq == "69101"  # Water tupelo — Bechtold Model 2
        cw_val = df >= mind ? Float32(5.3409) + Float32(0.7499)*df + Float32(0.1047)*crf :
                              (Float32(5.3409) + Float32(0.7499)*mind + Float32(0.1047)*crf)*(df/mind)
        if cw_val > Float32(37); cw_val = Float32(37); end

    elseif cweq == "69301"  # Blackgum — Bechtold Model 3
        cw_val = df >= mind ? Float32(5.5037) + Float32(1.0567)*df + Float32(0.0880)*crf + Float32(0.0610)*hi :
                              (Float32(5.5037) + Float32(1.0567)*mind + Float32(0.0880)*crf + Float32(0.0610)*hi)*(df/mind)
        if cw_val > Float32(50); cw_val = Float32(50); end

    elseif cweq == "69401"  # Swamp tupelo — Bechtold Model 2
        cw_val = df >= mind ? Float32(1.3564) + Float32(1.0991)*df + Float32(0.1243)*crf :
                              (Float32(1.3564) + Float32(1.0991)*mind + Float32(0.1243)*crf)*(df/mind)
        if cw_val > Float32(41); cw_val = Float32(41); end

    elseif cweq == "70101"  # Eastern hophornbeam — Bechtold Model 3
        cw_val = df >= mind ? Float32(7.8084) + Float32(0.8129)*df + Float32(0.0941)*crf - Float32(0.0817)*hi :
                              (Float32(7.8084) + Float32(0.8129)*mind + Float32(0.0941)*crf - Float32(0.0817)*hi)*(df/mind)
        if cw_val > Float32(39); cw_val = Float32(39); end
    elseif cweq == "70102"  # Eastern hophornbeam — Bragg
        cw_val = Float32(-33.898790) + Float32(38.731332)*df^Float32(0.152718)

    elseif cweq == "71101"  # Sourwood — Bechtold Model 3
        cw_val = df >= mind ? Float32(7.9750) + Float32(0.8303)*df + Float32(0.0423)*crf - Float32(0.0706)*hi :
                              (Float32(7.9750) + Float32(0.8303)*mind + Float32(0.0423)*crf - Float32(0.0706)*hi)*(df/mind)
        if cw_val > Float32(36); cw_val = Float32(36); end

    elseif cweq == "72101"  # Redbay — Bechtold Model 3
        cw_val = df >= mind ? Float32(4.2756) + Float32(1.0773)*df + Float32(0.1526)*crf + Float32(0.1650)*hi :
                              (Float32(4.2756) + Float32(1.0773)*mind + Float32(0.1526)*crf + Float32(0.1650)*hi)*(df/mind)
        if cw_val > Float32(25); cw_val = Float32(25); end

    elseif cweq == "73101"  # American sycamore — Bechtold Model 2
        cw_val = df >= mind ? Float32(-1.3973) + Float32(1.3756)*df + Float32(0.1835)*crf :
                              (Float32(-1.3973) + Float32(1.3756)*mind + Float32(0.1835)*crf)*(df/mind)
        if cw_val > Float32(66); cw_val = Float32(66); end

    elseif cweq == "74101"  # Balsam poplar — Bechtold Model 1
        cw_val = df >= mind ? Float32(6.2498) + Float32(0.8655)*df :
                              (Float32(6.2498) + Float32(0.8655)*mind)*(df/mind)
        if cw_val > Float32(25); cw_val = Float32(25); end
    elseif cweq == "74102"  # Balsam poplar — Bragg
        cw_val = Float32(7.522796) + Float32(0.125282)*df^Float32(1.855258)

    elseif cweq == "74201"  # Eastern cottonwood — Bechtold Model 1
        cw_val = df >= mind ? Float32(3.4375) + Float32(1.4092)*df :
                              (Float32(3.4375) + Float32(1.4092)*mind)*(df/mind)
        if cw_val > Float32(80); cw_val = Float32(80); end
    elseif cweq == "74203"  # Eastern cottonwood — Ek
        cw_val = df >= omind ? Float32(2.934) + Float32(2.538)*df^Float32(0.8617) :
                               (Float32(2.934) + Float32(2.538)*omind^Float32(0.8617))*(df/omind)
        if cw_val > Float32(80); cw_val = Float32(80); end

    elseif cweq == "74301"  # Bigtooth aspen — Bechtold Model 3
        cw_val = df >= mind ? Float32(0.6847) + Float32(1.1050)*df + Float32(0.1420)*crf - Float32(0.0265)*hi :
                              (Float32(0.6847) + Float32(1.1050)*mind + Float32(0.1420)*crf - Float32(0.0265)*hi)*(df/mind)
        if cw_val > Float32(43); cw_val = Float32(43); end
    elseif cweq == "74302"  # Bigtooth aspen — Bragg
        cw_val = Float32(4.031684) + Float32(1.132992)*df^Float32(1.024800)
    elseif cweq == "74303"  # Bigtooth aspen — Ek
        cw_val = df >= omind ? Float32(0.0750) + Float32(5.5770)*df^Float32(0.5996) :
                               (Float32(0.0750) + Float32(5.5770)*omind^Float32(0.5996))*(df/omind)
        if cw_val > Float32(43); cw_val = Float32(43); end

    elseif cweq == "74601"  # Quaking aspen — Bechtold Model 2
        cw_val = df >= mind ? Float32(0.7315) + Float32(1.3180)*df + Float32(0.0966)*crf :
                              (Float32(0.7315) + Float32(1.3180)*mind + Float32(0.0966)*crf)*(df/mind)
        if cw_val > Float32(39); cw_val = Float32(39); end
    elseif cweq == "74602"  # Quaking aspen — Bragg
        cw_val = Float32(2.303376) + Float32(2.371714)*df^Float32(0.807622)
    elseif cweq == "74603"  # Quaking aspen — Ek
        cw_val = df >= omind ? Float32(4.203) + Float32(2.129)*df^Float32(1.0000) :
                               (Float32(4.203) + Float32(2.129)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(43); cw_val = Float32(43); end

    elseif cweq == "76102"  # Pin cherry — Bragg (with MIND floor)
        cw_val = df >= mind ? Float32(4.102718) + Float32(1.396006)*df^Float32(1.077474) :
                              (Float32(4.102718) + Float32(1.396006)*mind^Float32(1.077474))*(df/mind)
        if cw_val > Float32(52); cw_val = Float32(52); end

    elseif cweq == "76201"  # Black cherry — Bechtold Model 3
        cw_val = df >= mind ? Float32(3.0237) + Float32(1.1119)*df + Float32(0.1112)*crf - Float32(0.0493)*hi :
                              (Float32(3.0237) + Float32(1.1119)*mind + Float32(0.1112)*crf - Float32(0.0493)*hi)*(df/mind)
        if cw_val > Float32(52); cw_val = Float32(52); end
    elseif cweq == "76202"  # Black cherry — Bragg
        cw_val = Float32(1.304425) + Float32(4.592688)*df^Float32(0.526895)
    elseif cweq == "76203"  # Black cherry — Ek
        cw_val = df >= omind ? Float32(0.621) + Float32(7.059)*df^Float32(0.5441) :
                               (Float32(0.621) + Float32(7.059)*omind^Float32(0.5441))*(df/omind)
        if cw_val > Float32(52); cw_val = Float32(52); end

    elseif cweq == "80201"  # White oak — Bechtold Model 3
        cw_val = df >= mind ? Float32(3.2375) + Float32(1.5234)*df + Float32(0.0455)*crf - Float32(0.0324)*hi :
                              (Float32(3.2375) + Float32(1.5234)*mind + Float32(0.0455)*crf - Float32(0.0324)*hi)*(df/mind)
        if cw_val > Float32(69); cw_val = Float32(69); end
    elseif cweq == "80203"  # White oak — Ek
        cw_val = df >= omind ? Float32(3.689) + Float32(1.838)*df^Float32(1.0000) :
                               (Float32(3.689) + Float32(1.838)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(69); cw_val = Float32(69); end
    elseif cweq == "80204"  # White oak — Krajicek
        cw_val = df >= omind ? Float32(1.8000) + Float32(1.8830)*df :
                               (Float32(1.8000) + Float32(1.8830)*omind)*(df/omind)
        if cw_val > Float32(69); cw_val = Float32(69); end

    elseif cweq == "80601"  # Scarlet oak — Bechtold Model 2
        cw_val = df >= mind ? Float32(0.5656) + Float32(1.6766)*df + Float32(0.0739)*crf :
                              (Float32(0.5656) + Float32(1.6766)*mind + Float32(0.0739)*crf)*(df/mind)
        if cw_val > Float32(66); cw_val = Float32(66); end

    elseif cweq == "80901"  # Northern pin oak — Bechtold Model 1
        cw_val = df >= mind ? Float32(4.8935) + Float32(1.6069)*df :
                              (Float32(4.8935) + Float32(1.6069)*mind)*(df/mind)
        if cw_val > Float32(44); cw_val = Float32(44); end

    elseif cweq == "81201"  # Southern red oak — Bechtold Model 2
        cw_val = df >= mind ? Float32(2.1517) + Float32(1.6064)*df + Float32(0.0609)*crf :
                              (Float32(2.1517) + Float32(1.6064)*mind + Float32(0.0609)*crf)*(df/mind)
        if cw_val > Float32(56); cw_val = Float32(56); end

    elseif cweq == "81701"  # Shingle oak — Bechtold Model 1
        cw_val = df >= mind ? Float32(9.8187) + Float32(1.1343)*df :
                              (Float32(9.8187) + Float32(1.1343)*mind)*(df/mind)
        if cw_val > Float32(54); cw_val = Float32(54); end

    elseif cweq == "81901"  # Turkey oak — Bechtold Model 1
        cw_val = df >= mind ? Float32(5.8858) + Float32(1.4935)*df :
                              (Float32(5.8858) + Float32(1.4935)*mind)*(df/mind)
        if cw_val > Float32(29); cw_val = Float32(29); end

    elseif cweq == "82001"  # Laurel oak — Bechtold Model 1
        cw_val = df >= mind ? Float32(6.3149) + Float32(1.6455)*df :
                              (Float32(6.3149) + Float32(1.6455)*mind)*(df/mind)
        if cw_val > Float32(54); cw_val = Float32(54); end

    elseif cweq == "82301"  # Bur oak — Bechtold Model 2
        cw_val = df >= mind ? Float32(1.7827) + Float32(1.6549)*df + Float32(0.0343)*crf :
                              (Float32(1.7827) + Float32(1.6549)*mind + Float32(0.0343)*crf)*(df/mind)
        if cw_val > Float32(61); cw_val = Float32(61); end
    elseif cweq == "82303"  # Bur oak — Ek
        cw_val = df >= omind ? Float32(0.942) + Float32(3.539)*df^Float32(0.7952) :
                               (Float32(0.942) + Float32(3.539)*omind^Float32(0.7952))*(df/omind)
        if cw_val > Float32(78); cw_val = Float32(78); end

    elseif cweq == "82401"  # Blackjack oak — Bechtold Model 2
        cw_val = df >= mind ? Float32(0.5443) + Float32(1.4882)*df + Float32(0.0565)*crf :
                              (Float32(0.5443) + Float32(1.4882)*mind + Float32(0.0565)*crf)*(df/mind)
        if cw_val > Float32(37); cw_val = Float32(37); end

    elseif cweq == "82601"  # Chinkapin oak — Bechtold Model 3
        cw_val = df >= mind ? Float32(0.5189) + Float32(1.4134)*df + Float32(0.1365)*crf - Float32(0.0806)*hi :
                              (Float32(0.5189) + Float32(1.4134)*mind + Float32(0.1365)*crf - Float32(0.0806)*hi)*(df/mind)
        if cw_val > Float32(45); cw_val = Float32(45); end

    elseif cweq == "82701"  # Water oak — Bechtold Model 3
        cw_val = df >= mind ? Float32(1.6349) + Float32(1.5443)*df + Float32(0.0637)*crf - Float32(0.0764)*hi :
                              (Float32(1.6349) + Float32(1.5443)*mind + Float32(0.0637)*crf - Float32(0.0764)*hi)*(df/mind)
        if cw_val > Float32(57); cw_val = Float32(57); end

    elseif cweq == "83001"  # Pin oak — Bechtold Model 3
        cw_val = df >= mind ? Float32(-5.6268) + Float32(1.7808)*df + Float32(0.1231)*crf + Float32(0.1578)*hi :
                              (Float32(-5.6268) + Float32(1.7808)*mind + Float32(0.1231)*crf + Float32(0.1578)*hi)*(df/mind)
        if cw_val > Float32(63); cw_val = Float32(63); end

    elseif cweq == "83101"  # Willow oak — Bechtold Model 2
        cw_val = df >= mind ? Float32(1.6477) + Float32(1.3672)*df + Float32(0.0846)*crf :
                              (Float32(1.6477) + Float32(1.3672)*mind + Float32(0.0846)*crf)*(df/mind)
        if cw_val > Float32(74); cw_val = Float32(74); end

    elseif cweq == "83201"  # Chestnut oak — Bechtold Model 2
        if df >= mind
            dlim = min(df, Float32(50))
            cw_val = Float32(2.1480) + Float32(1.6928)*dlim - Float32(0.0176)*dlim*dlim + Float32(0.0569)*crf
        else
            cw_val = (Float32(2.1480) + Float32(1.6928)*mind - Float32(0.0176)*mind*mind + Float32(0.0569)*crf)*(df/mind)
        end

    elseif cweq == "83301"  # Northern red oak — Bechtold Model 2
        cw_val = df >= mind ? Float32(2.8908) + Float32(1.4077)*df + Float32(0.0643)*crf :
                              (Float32(2.8908) + Float32(1.4077)*mind + Float32(0.0643)*crf)*(df/mind)
        if cw_val > Float32(82); cw_val = Float32(82); end
    elseif cweq == "83302"  # Northern red oak — Bragg
        cw_val = Float32(2.280575) + Float32(2.679718)*df^Float32(0.830741)
    elseif cweq == "83303"  # Northern red oak — Ek
        cw_val = df >= omind ? Float32(2.8500) + Float32(3.7820)*df^Float32(0.7968) :
                               (Float32(2.850) + Float32(3.782)*omind^Float32(0.7968))*(df/omind)
        if cw_val > Float32(82); cw_val = Float32(82); end
    elseif cweq == "83304"  # Northern red oak — Krajicek
        cw_val = df >= omind ? Float32(4.5100) + Float32(1.6700)*df :
                               (Float32(4.5100) + Float32(1.6700)*omind)*(df/omind)
        if cw_val > Float32(82); cw_val = Float32(82); end

    elseif cweq == "83501"  # Post oak — Bechtold Model 2
        cw_val = df >= mind ? Float32(1.6125) + Float32(1.6669)*df + Float32(0.0536)*crf :
                              (Float32(1.6125) + Float32(1.6669)*mind + Float32(0.0536)*crf)*(df/mind)
        if cw_val > Float32(45); cw_val = Float32(45); end

    elseif cweq == "83701"  # Black oak — Bechtold Model 2
        cw_val = df >= mind ? Float32(2.8974) + Float32(1.3697)*df + Float32(0.0671)*crf :
                              (Float32(2.8974) + Float32(1.3697)*mind + Float32(0.0671)*crf)*(df/mind)
        if cw_val > Float32(52); cw_val = Float32(52); end
    elseif cweq == "83703"  # Black oak — Ek
        cw_val = df >= omind ? Float32(4.5040) + Float32(2.4170)*df^Float32(1.0000) :
                               (Float32(4.5040) + Float32(2.4170)*omind^Float32(1.0000))*(df/omind)
        if cw_val > Float32(52); cw_val = Float32(52); end
    elseif cweq == "83704"  # Black oak — Krajicek
        cw_val = df >= omind ? Float32(4.5100) + Float32(1.6700)*df :
                               (Float32(4.5100) + Float32(1.6700)*omind)*(df/omind)
        if cw_val > Float32(52); cw_val = Float32(52); end

    elseif cweq == "83801"  # Live oak — Bechtold Model 1
        cw_val = df >= mind ? Float32(5.6694) + Float32(1.6402)*df :
                              (Float32(5.6694) + Float32(1.6402)*mind)*(df/mind)
        if cw_val > Float32(66); cw_val = Float32(66); end

    elseif cweq == "90101"  # Black locust — Bechtold Model 2
        cw_val = df >= mind ? Float32(3.0012) + Float32(0.8165)*df + Float32(0.1395)*crf :
                              (Float32(3.0012) + Float32(0.8165)*mind + Float32(0.1395)*crf)*(df/mind)
        if cw_val > Float32(48); cw_val = Float32(48); end

    elseif cweq == "93101"  # Sassafras — Bechtold Model 2
        cw_val = df >= mind ? Float32(4.6311) + Float32(1.0108)*df + Float32(0.0564)*crf :
                              (Float32(4.6311) + Float32(1.0108)*mind + Float32(0.0564)*crf)*(df/mind)
        if cw_val > Float32(29); cw_val = Float32(29); end

    elseif cweq == "95101"  # American basswood — Bechtold Model 3
        cw_val = df >= mind ? Float32(1.6871) + Float32(1.2110)*df + Float32(0.1194)*crf - Float32(0.0264)*hi :
                              (Float32(1.6871) + Float32(1.2110)*mind + Float32(0.1194)*crf - Float32(0.0264)*hi)*(df/mind)
        if cw_val > Float32(61); cw_val = Float32(61); end
    elseif cweq == "95102"  # American basswood — Bragg
        cw_val = Float32(7.172413) + Float32(0.662556)*df^Float32(1.127814)
    elseif cweq == "95103"  # American basswood — Ek
        cw_val = df >= omind ? Float32(0.1350) + Float32(3.7030)*df^Float32(0.7307) :
                               (Float32(0.1350) + Float32(3.7030)*omind^Float32(0.7307))*(df/omind)
        if cw_val > Float32(61); cw_val = Float32(61); end

    elseif cweq == "97101"  # Winged elm — Bechtold Model 2
        cw_val = df >= mind ? Float32(4.3649) + Float32(1.6612)*df + Float32(0.0643)*crf :
                              (Float32(4.3649) + Float32(1.6612)*mind + Float32(0.0643)*crf)*(df/mind)
        if cw_val > Float32(40); cw_val = Float32(40); end

    elseif cweq == "97201"  # American elm — Bechtold Model 3
        cw_val = df >= mind ? Float32(1.7296) + Float32(2.0732)*df + Float32(0.0590)*crf - Float32(0.0869)*hi :
                              (Float32(1.7296) + Float32(2.0732)*mind + Float32(0.0590)*crf - Float32(0.0869)*hi)*(df/mind)
        if cw_val > Float32(50); cw_val = Float32(50); end
    elseif cweq == "97202"  # American elm — Bragg
        cw_val = Float32(-53.239079) + Float32(61.327257)*df^Float32(0.060166)
    elseif cweq == "97203"  # American elm — Ek
        cw_val = df >= omind ? Float32(2.8290) + Float32(3.4560)*df^Float32(0.8575) :
                               (Float32(2.8290) + Float32(3.4560)*omind^Float32(0.8575))*(df/omind)
        if cw_val > Float32(72); cw_val = Float32(72); end

    elseif cweq == "97501"  # Slippery elm — Bechtold Model 3
        cw_val = df >= mind ? Float32(9.0023) + Float32(1.3933)*df - Float32(0.0785)*hi :
                              (Float32(9.0023) + Float32(1.3933)*mind - Float32(0.0785)*hi)*(df/mind)
        if cw_val > Float32(49); cw_val = Float32(49); end
    end

    # ---------------------------------------------------------------------------
    # 3. Final clamp [0.5, 99.9]
    # ---------------------------------------------------------------------------
    if cw_val < Float32(0.5); cw_val = Float32(0.5); end
    if cw_val > Float32(99.9); cw_val = Float32(99.9); end
    return cw_val
end
