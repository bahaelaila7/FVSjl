# =============================================================================
# volume.jl (bluemountains) — BM volume (bm VEQNNC). Chunk 8.
#   FW2W  → Flewelling FW2 (cr_fw2_vol, same as KT/EM/UT) — the main conifers WL/DF/GF/LP/ES/AF/PP.
#   616BEHW → region-6 Behre (minor species WP/MH/WJ/WB/LM/PY/YC/AS/CW/OS/OH — not in bmt01, DEFERRED).
# Merch (bm/grinit.f): TOPD=4.5, DBHMIN=7 (LP sp7=6), BFTOPD=4.5, BFMIND=7 (sp7=6).
# =============================================================================

function compute_volumes_bm!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq
    sd = s.coef.species
    iforst = bm_kodfor_remap(Int(s.plot.user_forest_code)) % 100   # forkod-remapped R6 forest (619→616, 8117→614)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = bm_bratio(sd, sp, d)
        dbhmin = sp == 7 ? 6.0f0 : 7.0f0
        if mdl == "FW2"
            v = cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = 6)
            t.cuft_vol[i] = max(v[1], 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin ? max(v[2], 0f0) : 0f0
        else                                                 # 616BEHW (region-6 Behre)
            fclass = bm_formcl(sp, iforst, d)                # form class keyed by BM species index (formcl.f)
            dbtbh = d * (1f0 - bark)                          # double bark thickness (fvsvol.f:153)
            dbhib = d - dbtbh
            vol2 = 0f0; vol4 = 0f0
            v1 = if h <= 17.3f0                               # R6VOL short-tree guard (TTH≤FC_HT):
                0.00272708f0 * dbhib * dbhib * h             # cylinder VOL(1); R6DIBS/R6VOL1 SKIPPED
            else
                v = bm_r6vol3(d, dbtbh, fclass, h, 1)        # ZONE 1 total cubic → VOL(1)
                mtopp = 4.5f0 * bark                         # TOPDIAM = TOPD·BARK (fvsvol.f)
                xlogs, ld1 = bm_r6dibs(d, fclass, mtopp, h)  # log bucking → small-end diams
                lv1, lv4 = bm_r6vol1(d, fclass, xlogs, ld1)  # per-log Scribner (VOL2) + merch cubic (VOL4)
                nlog = Int(floor(xlogs)); nacc = (xlogs - nlog) > 0f0 ? nlog + 1 : nlog
                for k in 1:nacc
                    vol2 += bm_anint(lv1[k])                  # r6vol.f:176 VOL(2)=Σ ANINT(LOGVOL(1))
                    vol4 += bm_anint(lv4[k] * 10f0) / 10f0    # r6vol.f:174 VOL(4)=Σ round(LOGVOL(4)·10)/10
                end
                v
            end
            t.cuft_vol[i] = max(v1, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(vol4, 0f0) : 0f0   # MCF=VOL(4), D≥DBHMIN
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin ? max(vol2, 0f0) : 0f0         # BdFt=VOL(2) Scribner, D≥BFMIND
        end
    end
    return s
end

# bm/NVEL formclas.f FORMCL_BM — R6 Blue Mountains form-class lookup. Binary search of the FIA
# code SPEC in FIAJSP (11 sorted codes), IFCDBH = (D−1)/10+1 clamped [1,5] (D>40.9→5), then
# FC = forest_table[ISPC, IFCDBH]. IFORST 4/7/14/16 → Malheur/Ochoco/Umatilla/Wallowa-Whitman;
# FULL 18-species form-class tables (bm/formcl.f DATA MALHFC/OCHOFC/UMATFC/WLWHFC, dimensioned MAXSP=18 × 5).
# Indexed by the BM SPECIES INDEX (1=WP 2=WL 3=DF 4=GF 5=MH 6=WJ 7=LP 8=ES 9=AF 10=PP 11=WB 12=LM 13=PY 14=YC
# 15=AS 16=CW 17=OS 18=OH), NOT a FIA subset — the earlier 11-col FIAJSP table left woodland minors (WJ/PY/YC/
# WB/LM/AS/CW/OH) defaulting to 80 (bmt01-DEFERRED). Stored [ifcdbh(1..5), sp(1..18)]; column-major from Fortran.
# Forest→table (formcl.f): 604→MALH 607→OCHO 614→UMAT 616→WLWH (iforst = kodfor%100 = 4/7/14/16).
const BM_FCL_MALH = Int[78 78 78 76 75 60 80 77 78 78 80 80 56 56 77 76 60 77;
                        78 79 77 78 79 60 83 80 80 78 81 81 60 66 77 78 60 77;
                        79 80 77 77 79 60 83 82 80 80 81 81 60 68 77 78 60 77;
                        81 82 80 76 79 60 80 84 82 82 82 82 60 68 77 78 60 77;
                        78 77 77 76 78 60 80 84 82 83 82 82 60 68 76 78 60 76]
const BM_FCL_OCHO = Int[78 78 79 76 75 60 70 82 78 76 80 80 56 56 77 76 60 77;
                        80 78 79 78 78 60 75 82 76 78 81 81 60 66 77 78 60 77;
                        80 80 76 77 79 60 75 82 74 78 81 81 60 68 77 78 60 77;
                        82 80 76 74 79 60 75 82 74 80 82 82 60 68 77 78 60 77;
                        80 80 76 74 78 60 75 82 74 80 82 82 60 68 76 78 60 76]
const BM_FCL_UMAT = Int[78 78 77 76 75 60 86 77 74 78 80 80 56 56 77 76 60 77;
                        78 78 77 78 75 60 86 77 74 78 81 81 60 66 77 78 60 77;
                        80 78 77 77 75 60 86 75 74 80 81 81 60 68 77 78 60 77;
                        81 78 77 76 79 60 86 75 75 81 82 82 60 68 77 78 60 77;
                        81 78 77 76 78 60 86 75 75 81 82 82 60 68 76 78 60 76]
const BM_FCL_WLWH = Int[78 78 78 76 75 60 85 84 78 78 80 80 56 56 77 76 60 77;
                        78 82 77 78 79 60 86 84 79 78 81 81 60 66 77 78 60 77;
                        78 77 77 77 79 60 85 84 79 80 81 81 60 68 77 78 60 77;
                        78 75 77 76 79 60 85 84 79 82 82 82 60 68 77 78 60 77;
                        78 75 77 76 78 60 85 84 79 83 82 82 60 68 76 78 60 76]

function bm_formcl(sp::Integer, iforst::Int, d::Real)::Int
    (sp < 1 || sp > 18) && return 80
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1)
    Float32(d) > 40.9f0 && (ifcdbh = 5)
    tbl = iforst == 4 ? BM_FCL_MALH : iforst == 7 ? BM_FCL_OCHO :
          iforst == 14 ? BM_FCL_UMAT : BM_FCL_WLWH
    return tbl[ifcdbh, sp]
end

# bm/NVEL r6vol3.f — Behre total-cubic profile (VOLEQ 616BEH*** → ZONE 1). Smalian-integrated taper
# DR = HRATIO/(0.62·HRATIO+0.38). Inputs: DBHOB, DBTBH (=D·(1-bark)), FCLASS (form class), HTTOT.
# Used for the BM minor species (WP/MH/WJ/WB/LM/PY/YC/AS/CW/OS/OH) whose VEQNNC = 616BEHW.
function bm_r6vol3(dbhob::Float32, dbtbh::Float32, fclass::Int, httot::Float32, zone::Int)::Float32
    topd = 4.0f0
    d17 = Float32(fclass) / 100.0f0 * dbhob
    h17 = zone == 1 ? 17.3f0 : 33.6f0
    dbhib = dbhob - dbtbh
    (dbhib <= 0f0 || dbhib > dbhob) && (dbhib = dbhob)
    if dbhib < topd                                            # whole stem below top diam
        return 0.00272708f0 * dbhib * dbhib * httot
    elseif d17 < topd                                          # two-log: butt (H17) + upper cylinder
        v = 0.00272708f0 * (dbhib * dbhib + d17 * d17) * h17
        v += 0.00272708f0 * (d17 * d17) * (httot - h17)
        return v
    end
    a = 0.62f0; b = 1.0f0 - a
    D = zeros(Float32, 21)
    vol = 0.00272708f0 * (dbhib * dbhib + d17 * d17) * h17     # butt log stump→17.3ft
    htup = httot - h17
    D[1] = d17
    s = 0.0f0; goto100 = false; goto130 = false; iexit = 21
    i = 2
    while i <= 20
        hratio = (htup - ((i - 1) * 16.3f0)) / htup
        if hratio <= 0.0f0; iexit = i; goto100 = true; break; end
        dr = hratio / (a * hratio + b)
        D[i] = dr * D[1]
        if D[i] < topd; iexit = i; goto100 = true; break; end
        vol += 0.00272708f0 * (D[i-1] * D[i-1] + D[i] * D[i]) * 16.3f0
        if D[i] == topd; s = 16.3f0; iexit = i; goto130 = true; break; end
        i += 1
    end
    if !goto130                                                # label 100 (incl. natural loop end, iexit=21)
        goto100 || (iexit = 21)
        dr = topd / d17
        hx = (dr * b * htup) / (1.0f0 - (a * dr))
        hh = (iexit - 2) * 16.3f0
        s = htup - hx - hh
        vol += 0.00272708f0 * (D[iexit-1] * D[iexit-1] + topd * topd) * s
    end
    htup2 = httot - (16.3f0 * (iexit - 2) + h17) - s           # label 130: top-of-tree cone
    vol += 0.00272708f0 * (topd * topd) * htup2
    return vol
end

# Fortran ANINT (round half away from zero); all volumes here are ≥0.
bm_anint(x::Real)::Float32 = Float32(floor(Float32(x) + 0.5f0))

# bm/NVEL r6vol1.f IFTR — Scribner board-foot table (132 log small-end-diameter entries), loaded once.
const BM_IFTR = let
    rows = readlines(joinpath(BM_DATADIR, "r6vol1_iftr.csv"))
    Int[parse(Int, strip(l)) for l in rows[2:end]]
end

# bm/NVEL r6dibs.f — ZONE-1 (16.3-ft log), A=0.62, total-height (TLH=0) log-bucking path
# (labels 70→140): Behre taper DR=HR/(A·HR+B), bucked to merch top MTOPP. Returns (XLOGS number
# of logs incl. fractional, LOGDIA(:,1) integer small-end diameters, length 21). DBHOB=outside-bark.
function bm_r6dibs(dbhob::Float32, fclass::Int, mtopp::Float32, th::Float32)
    a = 0.62f0; b = 1.0f0 - a; fc16 = Float32(fclass) / 100.0f0
    ld2 = zeros(Float32, 21); sl = zeros(Float32, 20); xl = zeros(Float32, 20)
    ld2[1] = dbhob * fc16
    sl[1] = 16.3f0; xl[1] = 16.3f0
    xlogs = 0.0f0; goto130 = false; goto100 = false; iat = 0
    if ld2[1] <= mtopp                                         # single log to top
        ld2[1] = mtopp; xlogs = 1.0f0; goto130 = true
    else
        h1 = th - 16.3f0
        i = 2
        while i <= 19
            hx = h1 - ((i - 1) * 16.3f0); hr = hx / h1
            if hr <= 0.0f0
                ld2[i] = 1.0f0; iat = i; goto100 = true; break
            end
            dr = hr / (a * hr + b); ld2[i] = dr * ld2[1]
            sl[i] = 16.3f0; xl[i] = 16.3f0
            if ld2[i] < mtopp
                iat = i; goto100 = true; break
            elseif ld2[i] > mtopp
                i += 1; continue
            else                                              # exactly at top → whole logs
                xlogs = Float32(i - 1); sl[i] = 16.3f0; xl[i] = 16.3f0; goto130 = true; break
            end
        end
        if !goto130 && !goto100                               # loop ran to i=19 → fall through, i=20
            iat = 20; goto100 = true
        end
        if goto100                                            # label 100: fractional top log
            i = iat
            dr = mtopp / ld2[1]
            hx = (dr * b * h1) / (1.0f0 - (a * dr))
            hh = Float32(i - 2) * 16.3f0
            s = h1 - hx - hh
            if s < 4.0f0                                      # <4ft stub → drop
                xlogs = Float32(i - 1); ld2[i] = 0.0f0; sl[i] = 0.0f0; xl[i] = s
            elseif s <= 12.0f0                                # 4–12ft → half log
                xlogs = Float32(i - 1) + 0.5f0; ld2[i] = mtopp; sl[i] = 8.0f0; xl[i] = s
            else                                              # >12ft → full log
                xlogs = Float32(i); ld2[i] = mtopp; sl[i] = 16.3f0; xl[i] = s
            end
        end
    end
    ld1 = zeros(Int, 21)
    for i in 1:20; ld1[i] = Int(floor(ld2[i] + 0.5f0)); end    # LOGDIA(:,1)=INT(LOGDIA(:,2)+0.5)
    return xlogs, ld1
end

# bm/NVEL r6vol1.f — ZONE-1 (IAPZ=1) per-log volumes from log small-end diameters:
#   LOGVOL(1,i) Scribner board (IFTR integer table, (IFTR·16+500)÷1000), LOGVOL(4,i) merch cubic
#   (butt log = 0.06239·D²·(FC/100)²+0.025624·D²; rest Smalian F=0.005454154). Skips INTL14
#   International (feeds VOL(10), not the .sum). Returns (logvol1, logvol4) length 20.
function bm_r6vol1(dbhob::Float32, fclass::Int, xlogs::Float32, ld1::Vector{Int})
    logs = Int(floor(xlogs)); lv1 = zeros(Float32, 20); lv4 = zeros(Float32, 20)
    for i in 1:logs
        kd = ld1[i]; k = kd <= 11 ? kd : kd + 6
        (k >= 1 && k <= 132) && (lv1[i] = Float32((BM_IFTR[k] * 16 + 500) ÷ 1000))
    end
    x = xlogs - Float32(logs)
    if x != 0.0f0                                             # fractional (8-ft) top log
        kd = ld1[logs+1]
        k = kd < 6 ? kd : (kd >= 12 ? kd + 6 : kd + 121)      # sequential-IF mapping (r6vol1.f:42-44)
        (k >= 1 && k <= 132) && (lv1[logs+1] = Float32((BM_IFTR[k] * 8 + 500) ÷ 1000))
    end
    f = 0.005454154f0
    lv4[1] = 0.06239f0 * dbhob^2 * (Float32(fclass) / 100f0)^2 + 0.025624f0 * dbhob^2
    if x != 0.0f0
        lv4[logs+1] = (Float32(ld1[logs+1])^2 * f + Float32(ld1[logs])^2 * f) / 2f0 * 8f0
    end
    for i in 2:logs
        lv4[i] = (Float32(ld1[i])^2 * f + Float32(ld1[i-1])^2 * f) / 2f0 * 16f0
    end
    return lv1, lv4
end
