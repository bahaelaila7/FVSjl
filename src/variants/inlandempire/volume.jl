# =============================================================================
# volume.jl (inlandempire) — IE volume (chunk 8). Region-1 NVEL. Merch standards IDENTICAL to KT
# (ie/grinit.f: TOPD=4.5/DBHMIN=7, sp7=6; BFTOPD=4.5/BFMIND=7, sp7=6; stump=1). sp1-14,23 use Flewelling
# FW2 (I00FW2W<fia>, reuse cr_fw2_vol); sp15-22 (PM/RM/PY/AS/CO/MM/PB/OH) use DVE/Behre (reuse cr_dve_vol;
# Behre sp17 deferred). VOLEQ strings DUMPED from the live binary (VOLEQDEF VAR='IE', forest 118).
# NOTE: VOLEQDEF is forest-keyed; this is the forest-118 (iet01) default — a full voleqdef port covers all forests.
# =============================================================================

# IE per-species VOLEQ (dumped from live FVSie sitset VEQNNC, forest 118).
const IE_VOL_EQ = String[
    "I00FW2W119",  # sp1
    "I00FW2W073",  # sp2
    "I00FW2W202",  # sp3
    "I00FW2W017",  # sp4
    "I00FW2W260",  # sp5
    "I00FW2W242",  # sp6
    "I00FW2W108",  # sp7
    "I00FW2W093",  # sp8
    "I00FW2W019",  # sp9
    "I00FW2W122",  # sp10
    "I00FW2W260",  # sp11
    "I00FW2W012",  # sp12
    "I00FW2W073",  # sp13
    "I00FW2W019",  # sp14
    "102DVEW106",  # sp15
    "102DVEW106",  # sp16
    "616BEHW231",  # sp17
    "102DVEW746",  # sp18
    "102DVEW740",  # sp19
    "200DVEW746",  # sp20
    "101DVEW375",  # sp21
    "200DVEW746",  # sp22
    "I00FW2W260",  # sp23
]

# IE Colville NF (region-6, forest 621, IFOR==5) VOLEQ — VOLEQDEF is forest-keyed and Colville is IE's ONLY
# region-6 forest (ie/forkod.f JFOR: only 621 is R6; the other 10 are region-1). Region-6 assigns the INGY
# Flewelling subregion (I11/I12/I13 → GFSUB/PPSUB/WLSUB coefficients in cr_fw2_vol, applied via _fw2_ingy_frow)
# for the conifers and the region-6 GIRARD-form-class BEHRE equation (616BEHW###, → ie_behre_vol) for the
# minor species — NOT the region-1 FW2/DVE strings above. DUMPED verbatim from the live FVSie_g16 "NATIONAL
# VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS" table for a forest-621 stand. (DF/LP/ES/AF resolve to the FW2 base
# row either way, so those were already bit-exact; WL/GF/WH/PP/MH pick up INGY subregion coefs, and every minor
# species — WB/LM/LL/PM/RM/PY/AS/CO/MM/PB/OH/OS — switches from FW2/DVE to Behre.)
const IE_VOL_EQ_COLVILLE = String[
    "I00FW2W119",  # sp1  WP
    "I11FW2W073",  # sp2  WL (WLSUB 11)
    "I11FW2W202",  # sp3  DF (base)
    "I11FW2W017",  # sp4  GF (GFSUB 11)
    "I11FW2W017",  # sp5  WH → GF equation (GFSUB 11)
    "I11FW2W242",  # sp6  RC
    "I11FW2W108",  # sp7  LP (base)
    "I13FW2W093",  # sp8  ES (base)
    "I11FW2W202",  # sp9  AF → DF equation (base)
    "I12FW2W122",  # sp10 PP (PPSUB 12)
    "I11FW2W017",  # sp11 MH → GF equation (GFSUB 11)
    "616BEHW101",  # sp12 WB Behre
    "616BEHW113",  # sp13 LM Behre
    "616BEHW072",  # sp14 LL Behre
    "616BEHW000",  # sp15 PM Behre
    "616BEHW066",  # sp16 RM Behre
    "616BEHW231",  # sp17 PY Behre
    "616BEHW746",  # sp18 AS Behre
    "616BEHW740",  # sp19 CO Behre
    "616BEHW321",  # sp20 MM Behre
    "616BEHW375",  # sp21 PB Behre
    "616BEHW998",  # sp22 OH Behre
    "616BEHW299",  # sp23 OS Behre
]

# ie/formcl.f COLVFC — Colville NF (IFOR=5, forest 621, a Region-6 forest) Girard form-class table,
# [IFCDBH 1..5, ISPC 1..23]. Column-major Fortran DATA transcribed to [ifcdbh, sp]. Region-1 IE forests
# ignore form class (their FW2/DVE routines don't use it) and get FC=80; Colville feeds it to the R6 Behre
# routine. FRMCLS keyword override (ie/formcl.f) not modeled — FIA stands carry none.
const IE_COLVFC = permutedims(Float32[
    78 78 78 76 76 64 80 77 78 78 75 81 75 78 56 56 56 77 76 70 70 70 75
    80 78 76 78 78 65 82 79 76 80 78 82 75 76 56 56 60 77 78 70 70 70 78
    80 80 75 77 80 66 82 80 74 80 79 82 75 74 56 56 60 77 78 70 70 70 79
    82 80 74 76 80 66 80 80 74 82 79 80 74 74 56 56 60 77 78 70 70 70 79
    80 80 74 76 82 66 80 81 74 80 78 80 74 74 56 56 60 77 78 70 70 70 78])

# ie/formcl.f FORMCL: Region-6 Colville (IFOR=5) uses COLVFC; all other IE forests default FC=80.
@inline function ie_formcl(sp::Integer, ifor::Int, d::Real)::Int
    (ifor == 5 && sp >= 1 && sp <= 23) || return 80
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1); Float32(d) > 40.9f0 && (ifcdbh = 5)
    return Int(IE_COLVFC[sp, ifcdbh])
end

# IE sp17 (PY, VEQNNC 616BEHW231) region-6 Behre volume — reuses the validated BM/SO R6 machinery
# (bm_r6vol3/bm_r6dibs/bm_r6vol1) with the IE form class. IE TOPD=BFTOPD=4.5 ⇒ MTOPP=4.5·BARK (fvsvol.f).
function ie_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32)
    fclass = ie_formcl(sp, ifor, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    vol2 = 0f0; vol4 = 0f0; ht1 = 0f0
    v1 = if h <= 17.3f0                                    # R6VOL short-tree guard (TTH≤FC_HT): cylinder VOL(1)
        0.00272708f0 * dbhib * dbhib * h
    else
        v = bm_r6vol3(d, dbtbh, fclass, h, 1)             # ZONE 1 total cubic → VOL(1)
        mtopp = 4.5f0 * bark                              # TOPDIAM = TOPD·BARK
        xlogs, ld1, xl = bm_r6dibs(d, fclass, mtopp, h)   # log bucking → small-end diameters
        # r6vol.f:121-125 HT1PRD = 1.0 + Σ XLEN(1:20) (not reached when DBHIB<MTOPP ⇒ GO TO 1000)
        if dbhib >= mtopp
            ht1 = 1.0f0
            for k in 1:20; ht1 += xl[k]; end
        end
        lv1, lv4 = bm_r6vol1(d, fclass, xlogs, ld1)       # per-log Scribner (VOL2) + merch cubic (VOL4)
        nlog = Int(floor(xlogs)); nacc = (xlogs - nlog) > 0f0 ? nlog + 1 : nlog
        for k in 1:nacc
            vol2 += bm_anint(lv1[k]); vol4 += bm_anint(lv4[k] * 10f0) / 10f0
        end
        v
    end
    return (max(v1, 0f0), max(vol4, 0f0), max(vol2, 0f0), ht1)   # 4th = HT1PRD (fvsvol.f → HT2TD)
end

# ---------------------------------------------------------------------------
# ie_snag_bole_cuft — the FFE snag-bole TOTAL cubic (FVS FMSVL2 'D' ⇒ TCF) for IE/KT. IE/KT's `vol_eq`
# is a Region-1 NVEL string (Flewelling FW2 / Gevorkiantz DVE / R6 Behre), NOT an R8-Clark string, so the
# shared _R8CLARK_VOL snag path returns 0 ⇒ the input/SNAGINIT snag bole collapses to 0 and update_snags!
# falls back to the Jenkins ABOVEGROUND stem biomass — ~2.3× the true total-stem cubic (the 18" input snag:
# Jenkins 1.146 t vs FVS TVOLI·V2T 0.505 t) — over-loading the coarse (3"+) down-wood ⇒ over-fire. FVS
# FMSVOL/FMSVL2 (fmsvol.f:152-155) returns MAX(X,TCF) for the western NVEL variants, so the western snag
# bole (report AND CWD1 falldown) == the total cubic VOL(1). Mirror `cr_snag_bole_cuft`/`nc_snag_bole_cuft`.
function ie_snag_bole_cuft(s::StandState, sp::Int, d::Float32, h::Float32)::Float32
    (d < 1f0 || h <= 0f0 || sp < 1) && return 0f0
    veq = s.species.vol_eq
    sp > length(veq) && return 0f0
    eq = veq[sp]
    if s.variant isa Kootenai                        # KT: FW2 only (kt/bratio.f bark), no DVE/Behre
        bark = bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, d)
        startswith(eq, "I") || return 0f0
        return max(cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = 1)[1], 0f0)
    end
    bark = ie_bratio(sp, d)
    occursin("BEH", eq) && return max(ie_behre_vol(sp, Int(s.plot.forest_idx), d, h, bark)[1], 0f0)
    v = startswith(eq, "I") ?
        cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = 1) :
        cr_dve_vol(eq, d, h)
    return max(v[1], 0f0)                             # VOL(1) total cubic (= FMSVL2 'D' TCF)
end

# ---------------------------------------------------------------------------
# IE minor-species (aspen/cottonwood/mtn-mahogany/paper-birch/pinyon) volume.
# IE assigns these species NVEL "DVE"-method VOLEQ strings that are REGION-1
# ("102…"/"101…") and REGION-2 ("200…") — NOT the region-3 r3d2hv that
# `cr_dve_vol` implements. FVS routes them (vols.f METHC=6 → NATCRS → VOLINIT →
# dvest.f) by VOLEQ prefix: region-1 '02'-'06' → R1KEMP (r1kemp.f), region-1
# '01' → R1ALLENC (r1allen.f), region-2 → R2OLDV (r2oldv.f). The total-cubic
# call carries PROD='02' (fvsvol.f:184), which fixes R1KEMP's KLASS (live→3,
# dead→1). Sending these to `cr_dve_vol` made every aspen/cottonwood/birch tree
# return 0 cuft (r3d2hv has no branch for spc 746/740/375) — the IE TCuFt=0 bug.
# All three routines return the FVS VOL(1..15) convention: [1]=total cubic,
# [4]=merch cubic, [2]=board, [7]=topwood. (Measured vs FVSie_g16 dvest/vols dumps.)

# r1kemp.f species index (VOLEQ(8:10) → ISPEC) for the region-1 '02'-'06' D2H tables.
function _ie_r1kemp_ispec(spc::AbstractString)
    spc == "746" ? 1  : spc == "740" ? 2  : spc == "017" ? 3  : spc == "019" ? 4  :
    (spc == "070" || spc == "073") ? 5 : (spc == "090" || spc == "093") ? 6 :
    spc == "101" ? 7  : spc == "108" ? 8  : spc == "119" ? 9  : spc == "122" ? 10 :
    spc == "202" ? 11 : (spc == "240" || spc == "242") ? 12 :
    (spc == "260" || spc == "263") ? 13 : spc == "060" ? 14 : spc == "106" ? 15 : 0
end

# r1kemp.f CBVOLE(ISPEC,1..11) cubic-foot coefficient table (region-1 Kemp D2H).
const _IE_CBVOLE = (
    (0.3482f0,-0.0384f0,0.001427f0,-0.842503f0,0.224f0,-0.343f0,0.217f0,1.071f0,0f0,0f0,0f0),
    (0.1064f0,-0.00778f0,0.000176f0,-0.265342f0,0.204f0,-0.749f0,0.194f0,4.285f0,0f0,0f0,0f0),
    (0.3386f0,-0.03359f0,0.001109f0,-0.918645f0,0.219f0,-0.563f0,0.197f0,9.969f0,0.2153f0,-0.00167f0,0.50f0),
    (0.4529f0,-0.052f0,0.002003f0,-1.113416f0,0.183f0,1.449f0,0.117f0,26.222f0,0.2153f0,0.00167f0,0.67f0),
    (0.4172f0,-0.04693f0,0.001782f0,-1.086592f0,0.17f0,-0.056f0,0.132f0,19.409f0,0.1922f0,0.09023f0,0.35f0),
    (0.2619f0,-0.02345f0,0.000671f0,-0.716502f0,0.214f0,0.48f0,0.174f0,19.041f0,0.2306f0,0.14528f0,0.35f0),
    (0.6808f0,-0.07974f0,0.003113f0,-1.692512f0,0.221f0,1.052f0,0.197f0,5.369f0,0.2306f0,0.14528f0,0.35f0),
    (0.6808f0,-0.07974f0,0.003113f0,-1.692512f0,0.221f0,1.052f0,0.197f0,5.369f0,0.2306f0,0.14528f0,0.35f0),
    (0.4544f0,-0.05119f0,0.001945f0,-1.14765f0,0.206f0,0.166f0,0.194f0,4.508f0,0.2306f0,0.14528f0,0.35f0),
    (0.4041f0,-0.04535f0,0.001726f0,-1.054732f0,0.203f0,-1.656f0,0.218f0,-9.637f0,0.2306f0,0.14528f0,0.25f0),
    (0.5125f0,-0.05817f0,0.002208f0,-1.320519f0,0.178f0,0.437f0,0.165f0,7.702f0,0.1795f0,0.16949f0,0.47f0),
    (0.3349f0,-0.03565f0,0.001273f0,-0.851441f0,0.174f0,1.141f0,0.146f0,8.931f0,0.1922f0,0.09023f0,0.67f0),
    (0.2213f0,-0.01913f0,0.000533f0,-0.635045f0,0.209f0,-0.991f0,0.210f0,2.544f0,0.2153f0,-0.00167f0,0.43f0),
    (0f0,0f0,0f0,0f0,0.211f0,-0.597f0,0.211f0,-0.597f0,0f0,0f0,0f0),
    (0f0,0f0,0f0,0f0,0.211f0,-0.597f0,0.211f0,-0.597f0,0f0,0f0,0f0))

# r1kemp.f BFVOL(ISPEC,1,1..4) — board-foot "01" table (JTAB=1; region-1 '02'⇒JTAB=1).
const _IE_BFVOL01 = (
    (1.197f0,-18.544f0,1.216f0,-21.309f0),(1.046f0,-15.966f0,1.140f0,-46.735f0),
    (1.293f0,-34.127f0,1.218f0,10.603f0),(1.011f0,-11.403f0,0.694f0,124.425f0),
    (0.997f0,-29.790f0,0.841f0,85.150f0),(1.149f0,-11.851f0,1.158f0,1.620f0),
    (1.208f0,-8.085f0,1.103f0,14.111f0),(1.208f0,-8.085f0,1.103f0,14.111f0),
    (1.189f0,-26.729f0,1.181f0,-32.516f0),(1.201f0,-50.340f0,1.595f0,-298.784f0),
    (1.003f0,-25.332f0,1.011f0,-9.522f0),(0.878f0,-10.742f0,0.799f0,-4.064f0),
    (1.203f0,-37.314f0,1.306f0,-50.680f0),(1.208f0,-8.085f0,1.103f0,14.111f0),
    (0f0,0f0,0f0,0f0))

# R1KEMP (r1kemp.f) — region-1 Kemp D2H cubic + board. Returns (CBGRS, BFGRS) = (VOL(1)=VOL(4), VOL(2)).
# PROD='02' always for the total-cubic call (fvsvol.f:184) ⇒ KLASS=3 (live) / KLASS=1 (dead, ispec≠8).
function ie_r1kemp_vol(ispec::Int, d::Float32, h::Float32, islive::Bool)
    (ispec < 1 || ispec > 15) && return (0f0, 0f0)
    cb = _IE_CBVOLE[ispec]
    d2h100 = d * d * h / 100f0
    # KLASS (r1kemp.f:253): live PROD='02'⇒3; dead (ispec≠8)⇒1.
    klass = islive ? 3 : 1
    # CBGRS — cubic foot (r1kemp.f:345-366).
    cbgrs = if ispec == 14 || ispec == 15
        d < 5f0 ? 0f0 : (d <= 20.5f0 ? cb[5]*d2h100 + cb[6] : cb[7]*d2h100 + cb[8])
    else
        if d < 5f0
            (cb[9]*d2h100 + cb[10]) * cb[11]
        elseif d <= 9.5f0
            d2h100 * (cb[1]*d + cb[2]*d*d + cb[3]*d*d*d + cb[4])
        elseif d <= 20.5f0
            cb[5]*d2h100 + cb[6]
        else
            cb[7]*d2h100 + cb[8]
        end
    end
    # Cubic minimums (r1kemp.f:373-384): KLASS=3⇒2.4, KLASS≤2 (not dead-LP/WP)⇒1.6.
    if klass == 3
        cbgrs < 2.4f0 && (cbgrs = 2.4f0)
    else
        cbgrs < 1.6f0 && (cbgrs = 1.6f0)
    end
    # BFGRS — board foot, JTAB=1 "01" table (r1kemp.f:272-276, 336).
    bf = _IE_BFVOL01[ispec]
    bfgrs = d < 21f0 ? bf[1]*d2h100 + bf[2] : bf[3]*d2h100 + bf[4]
    bfgrs < 10f0 && (bfgrs = 10f0)
    return (cbgrs, bfgrs)
end

# R2OLDV (r2oldv.f) — region-2 "200DVEW746" aspen (RM-232). Returns (TCUFT, GCUFT, GRSBDT) = VOL(1),VOL(4),VOL(2).
function ie_r2oldv_746(d::Float32, h::Float32)
    d2h = d * d * h
    tcuft = d2h <= 12470f0 ? 0.002219f0*d2h : 0.001896f0*d2h + 4.0267f0
    gcuft = d2h <= 11800f0 ? 0.002195f0*d2h - 0.9076f0 : 0.001837f0*d2h + 3.3075f0
    grsbdt = 0f0
    if d > 7f0
        grsbdt = d2h <= 2500f0 ? 8f0 : (d2h <= 8850f0 ? 0.011389f0*d2h - 20.5112f0 : 0.010344f0*d2h - 11.2615f0)
    end
    return (max(tcuft, 0f0), max(gcuft, 0f0), max(grsbdt, 0f0))
end

# R1ALLENC (r1allen.f) — region-1 '01' cubic. IE uses only "101DVEW375" paper birch (ISPC=12); returns VOL(1)=VOL(4)=CUVOL.
function ie_r1allenc_375(d::Float32, h::Float32)
    d2h = d * d * h
    cuvol = d < 5f0 ? 0f0 : (d < 11f0 ? 0.988264f0 + 0.002732f0*d2h : 2.512836f0 + 0.002446f0*d2h)
    return max(cuvol, 0f0)
end

# IE minor-species volume dispatch (dvest.f): returns FVS VOL(1..15) with [1]=total,[4]=merch,[2]=board,[7]=topwood.
function ie_dve_vol(eq::AbstractString, d::Float32, h::Float32, islive::Bool)
    vol = zeros(Float32, 15)
    length(eq) >= 10 || return vol
    spc = eq[8:10]
    if eq[1] == '1'
        r23 = eq[2:3]
        if r23 == "01"                                   # region-1 '01' → R1ALLENC (birch 375)
            spc == "375" && (v = ie_r1allenc_375(d, h); vol[1] = v; vol[4] = v)
        else                                             # region-1 '02'-'06' → R1KEMP
            cb, bf = ie_r1kemp_vol(_ie_r1kemp_ispec(spc), d, h, islive)
            vol[1] = cb; vol[4] = cb; vol[2] = bf
        end
    elseif eq[1] == '2'                                  # region-2 → R2OLDV
        if spc == "746"
            tc, gc, bd = ie_r2oldv_746(d, h)
            vol[1] = tc; vol[4] = gc; vol[2] = bd
        end
    end
    return vol
end

function compute_volumes!(s::StandState, ::InlandEmpire)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq
    ifor = Int(s.plot.forest_idx)
    topd = 4.5f0; bftopd = 4.5f0; stump = 1.0f0; iregn = 1
    # vols.f:86-90 zeroes HT2TD for every record; FVSVOL/NATCRS then fills the FW2 merch-top heights.
    fill!(t.merch_top_cf, 0f0); fill!(t.merch_top_bf, 0f0)
    htb = zeros(Float32, 2)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        dbhmin = sp == 7 ? 6f0 : 7f0
        bfmind = sp == 7 ? 6f0 : 7f0
        bark = ie_bratio(sp, d)
        eq = veq[sp]
        if occursin("BEH", eq)                            # region-6 Behre (sp17 PY, 616BEHW231)
            # Broken/dead-top trees (trunc = break-ht·100, norm_ht = predicted full ht·100): FVS builds the
            # full-bole volume at NORMHT then trims to the break via the Behre form-class taper (cftopk/bftopk),
            # exactly as the FW2 variants do. (Validated vs FVSie_g16 treelist: D=11.7 broken tree 9.16/7.44/29.5
            # vs live 9.2/7.4/29.5.)
            broken = t.trunc[i] > 0 && t.norm_ht[i] > 0
            hbase = broken ? Float32(t.norm_ht[i]) / 100f0 : h
            tcf, mcf, bf, ht1 = ie_behre_vol(sp, ifor, d, hbase, bark)
            d >= dbhmin && (t.merch_top_cf[i] = ht1)           # fvsvol.f:337-339 (CF call, same MTOPP)
            d >= bfmind && (t.merch_top_bf[i] = ht1)           # fvsvol.f:484-487 (BF call)
            if broken && tcf > 0f0 && hbase >= 4.5f0
                vmax = tcf
                tcf, mcf = cr_cftopk(tcf, mcf, d, hbase, vmax, bark, Int(t.trunc[i]), 1f0, 4.5f0)
                bf = cr_bftopk(bf, d, hbase, vmax, bark, Int(t.trunc[i]), 1f0, 4.5f0)
            end
            t.cuft_vol[i] = max(tcf, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= bfmind ? max(bf, 0f0) : 0f0
            continue
        end
        if startswith(eq, "I")                            # Flewelling FW2 (sp1-14,23)
            # Broken/dead-top trees (vols.f:145-146,193): the full cubic/board is built at the
            # PREDICTED FULL HEIGHT (NORMHT), then trimmed to the standing broken stem via the
            # Behre taper (CFTOPK/BFTOPK, TOPD=4.5). Mirrors the Behre/DVE branches below and KT's
            # FW2 path. (Validated vs FVSie_g16 vols dump: D=15.04 NORMHT=57.81 ITRUNC=16 broken
            # redcedar VMAX 27.3 → TCF 14.84 / MCF 10.61.)
            broken = t.trunc[i] > 0 && t.norm_ht[i] > 0
            hbase = broken ? Float32(t.norm_ht[i]) / 100f0 : h
            # sf_hs: MERLEN's merch-top height via the faithful SF_HS Newton (profile.f MERLEN → sf_hs.f), which also
            # supplies HT1PRD → HT2TD (fvsvol.f:337-339 cubic, :484-487 board).
            v = cr_fw2_vol(eq, d, hbase; bark = bark, topd = topd, bftopd = bftopd, stump = stump, iregn = iregn,
                           sf_hs = true, ht2td = htb)
            d >= dbhmin && (t.merch_top_cf[i] = htb[1])
            d >= bfmind && (t.merch_top_bf[i] = htb[2])
            tcf = max(v[1], 0f0)
            mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
            bf  = d >= bfmind ? v[2] : 0f0
            if broken && tcf > 0f0 && hbase >= 4.5f0
                vmax = tcf
                tcf, mcf = cr_cftopk(tcf, mcf, d, hbase, vmax, bark, Int(t.trunc[i]), stump, topd)
                bf = cr_bftopk(bf, d, hbase, vmax, bark, Int(t.trunc[i]), stump, bftopd)
            end
            t.cuft_vol[i] = tcf; t.merch_cuft_vol[i] = max(mcf, 0f0)
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        else                                              # region-1/2 NVEL DVE (aspen/cottonwood/mm/birch/pinyon)
            # FVS routes these by VOLEQ prefix (R1KEMP / R1ALLENC / R2OLDV), PROD='02' total-cubic call.
            broken = t.trunc[i] > 0 && t.norm_ht[i] > 0
            hbase = broken ? Float32(t.norm_ht[i]) / 100f0 : h
            v = ie_dve_vol(eq, d, hbase, i <= t.n)
            tcf = max(v[1], 0f0)
            mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
            bf  = d >= bfmind ? max(v[2], 0f0) : 0f0
            # CFTOPK/BFTOPK (vols.f:191-193): broken-top trees have their full-height cubic + board
            # reduced to the standing broken stem via the Behre taper (TOPD=4.5).
            if broken && tcf > 0f0 && hbase >= 4.5f0
                vmax = tcf
                tcf, mcf = cr_cftopk(tcf, mcf, d, hbase, vmax, bark, Int(t.trunc[i]), 1f0, 4.5f0)
                bf = cr_bftopk(bf, d, hbase, vmax, bark, Int(t.trunc[i]), 1f0, 4.5f0)
            end
            t.cuft_vol[i] = max(tcf, 0f0); t.merch_cuft_vol[i] = max(mcf, 0f0)
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        end
    end
    return s
end
