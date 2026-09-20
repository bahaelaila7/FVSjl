# =============================================================================
# crown.jl (bluemountains) — CCF / crown width (bm/ccfcal.f). Chunk 5 (partial: the CCF
# piece is a prerequisite for chunk-3 DG, since RELDEN = stand CCF feeds CONSPP).
#
# bm/ccfcal.f MODE=1 per-tree CCFT (Paine & Hann, Oregon RP46):
#   CASE(1:12,15,17):  D>=1  → RD1 + D·RD2 + D²·RD3 ; 0.1<D<1 → RDA·D^RDB ; D<=0.1 → 0.001
#   CASE(13,14,16,18): D<1   → D·(RD1+RD2+RD3)      ; D>=1    → RD1 + RD2·D + RD3·D²
# Stand CCF = Σ CCFT·TPA = RELDEN.  Coeffs from data/bluemountains/ccf_coeffs_bm.csv (verified vs source).
# =============================================================================

let
    path = joinpath(BM_DATADIR, "ccf_coeffs_bm.csv")
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    col(j) = Float32[parse(Float32, rows[sp][j]) for sp in 1:18]
    global const BM_RD1 = col(2); global const BM_RD2 = col(3); global const BM_RD3 = col(4)
    global const BM_RDA = col(5); global const BM_RDB = col(6)
end

@inline function bm_tree_ccf(sp::Integer, d::Real)::Float32
    # No D≤0 guard (bm/ccfcal.f:38-56): a zero diameter falls to the D≤0.1 ⇒ 0.001 branch for sp 1-12,15,17.
    # The CRATET backdating DENSE zeroes IMC-9 (HISTORY 8/9) dead WK3 and still counts them at 0.001·P in the
    # point CCF the cycle-0 crown dub reads (1127576412290487: TPCCF 79.977 → 80.037 = live, 10 dead × 6.018 TPA).
    dd = Float32(d)
    if sp == 13 || sp == 14 || sp == 16 || sp == 18
        return dd < 1f0 ? dd * (BM_RD1[sp] + BM_RD2[sp] + BM_RD3[sp]) :
                          BM_RD1[sp] + BM_RD2[sp] * dd + BM_RD3[sp] * fpow(dd, 2f0)   # ccfcal.f:53 RD3*D**2.0 (powf)
    else
        return dd >= 1f0 ? BM_RD1[sp] + dd * BM_RD2[sp] + dd * dd * BM_RD3[sp] :
               dd > 0.1f0 ? BM_RDA[sp] * fpow(dd, BM_RDB[sp]) : 0.001f0     # ccfcal.f:44 RDA*(D**RDB) (powf)
    end
end

# --- bm/crown.f WEIBULL crown-ratio change model (= TT/UT/CR form). Coeffs from crown_coeffs_bm.csv. ---
let
    path = joinpath(BM_DATADIR, "crown_coeffs_bm.csv")
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    col(name) = Float32[parse(Float32, rows[sp][findfirst(==(name), split(strip(readlines(path)[1]), ','))]) for sp in 1:18]
    global const BM_WEIBA  = col("WEIBA");  global const BM_WEIBB0 = col("WEIBB0"); global const BM_WEIBB1 = col("WEIBB1")
    global const BM_WEIBC0 = col("WEIBC0"); global const BM_WEIBC1 = col("WEIBC1")
    global const BM_CRC0   = col("C0");     global const BM_CRC1   = col("C1")
    global const BM_CRNMLT = col("CRNMLT"); global const BM_CR_DLOW = col("DLOW"); global const BM_CR_DHI = col("DHI")
end

# --- bm/dubscr.f DUBSCR — crown-ratio dub for read-inventory D<1 (missing CR) + regen inserts. ---
# Logistic (sp 1-12,15,17) / linear-rescale (sp 13,14,16,18) crown model. Coeffs verified vs bm/dubscr.f DATA.
const BM_BCR0  = Float32[-1.669490,-1.669490,-.426688,-.426688,-.426688,-2.19723,-1.669490,-.426688,-.426688,-1.669490,-1.66949,-1.66949,6.489813,7.558538,-.426688,5.0,-1.669490,5.0]
const BM_BCR1  = Float32[-.209765,-.209765,-.093105,-.093105,-.093105,0.0,-.209765,-.093105,-.093105,-.209765,-.209765,-.209765,0.0,0.0,-.093105,0.0,-.209765,0.0]
const BM_BCR2  = Float32[0.0,0.0,.022409,.022409,.022409,0.0,0.0,.022409,.022409,0.0,0.0,0.0,-0.029815,-0.015637,.022409,0.0,0.0,0.0]
const BM_BCR3  = Float32[.003359,.003359,.002633,.002633,.002633,0.0,.003359,.002633,.002633,.003359,.003359,.003359,-0.009276,-0.009064,.002633,0.0,.003359,0.0]
const BM_BCR5  = Float32[.011032,.011032,0.0,0.0,0.0,0.0,.011032,0.0,0.0,.011032,.011032,.011032,0.0,0.0,0.0,0.0,.011032,0.0]
const BM_BCR6  = Float32[0.0,0.0,-.045532,-.045532,-.045532,0.0,0.0,-.045532,-.045532,0.0,0.0,0.0,0.0,0.0,-.045532,0.0,0.0,0.0]
const BM_BCR8  = Float32[.017727,.017727,0.0,0.0,0.0,0.0,.017727,0.0,0.0,.017727,.017727,.017727,0.0,0.0,0.0,0.0,.017727,0.0]
const BM_BCR9  = Float32[-.000053,-.000053,.000022,.000022,.000022,0.0,-.000053,.000022,.000022,-.000053,-.000053,-.000053,0.0,0.0,.000022,0.0,-.000053,0.0]
const BM_BCR10 = Float32[.014098,.014098,-.013115,-.013115,-.013115,0.0,.014098,-.013115,-.013115,.014098,.014098,.014098,0.0,0.0,-.013115,0.0,.014098,0.0]
const BM_DUBSD = Float32[.5000,.5000,.6957,.6957,.6957,0.200,.6124,.6957,.6957,.4942,.5000,.5000,2.0426,1.9658,.9310,0.500,.4942,0.500]

# so/adjmai.f ADJMAI — adjusted MAI site-index polynomial (10 eq-groups keyed by FIA site-species code).
# maical.f: ISPNUM maps BM species idx 1-18 → FIA code for the site-index lookup; RMAI capped at 128.
const _BM_ISPNUM = Int32[119,117,202,15,264,101,108,93,21,122,101,101,101,42,746,746,122,746]
const _BM_ADJ_ISP = Int32[19,41,42,71,81,94,95,242,263,264,98,298, 101,103,108,109,120,124, 119, 116,117,122, 201,202, 73, 92,93, 11,15,17,19,20,21,22, 211,212]
const _BM_ADJ_IMAP = Int32[fill(Int32(1),12); fill(Int32(2),6); Int32(3); fill(Int32(4),3); fill(Int32(5),2); Int32(6); fill(Int32(7),2); fill(Int32(9),7); fill(Int32(10),2)]

@inline function _bm_adjmai_group(grp::Integer, si::Float32)::Float32
    a = 0f0
    if grp == 1;      si < 33f0 && return 0f0; a = -63.689706f0 + 1.9402941f0*si
    elseif grp == 2;  si < 11f0 && return 0f0; a = -12.0388f0 + 1.18672f0*si
    elseif grp == 3;  a = 5.972615f0 + 1.857675f0*si
    elseif grp == 4;  a = 2.305357f0 + 0.033890056f0*si + 0.0090108543f0*si*si
    elseif grp == 5;  si < 29f0 && return 0f0; a = -10.303313f0 + .032929911f0*si + .012207163f0*si*si - .00003543129f0*si*si*si
    elseif grp == 6;  si < 11f0 && return 0f0; a = -6.0892857f0 + .45178571f0*si + .014464286f0*si*si
    elseif grp == 7;  si < 10f0 && return 0f0; a = -18.4f0 + 1.92f0*si
    elseif grp == 8;  si < 32f0 && return 0f0; a = -53.892857f0 + 1.7178571f0*si
    elseif grp == 9;  a = -4.89001f0 + 311.29546f0*((exp((si/170f0-1f0)^3/0.343f0)-0.055f0)/0.95f0)
    elseif grp == 10; si < 62f0 && return 0f0; a = 157.94643f0 - 1.78125f0*si + .014330357f0*si*si
    end
    a < 0f0 && (a = 0f0)                                  # ADJMAI = ADJMAI*POINTS/10 (POINTS=10 ⇒ ×1), floored 0
    return a
end

# adjmai.f DO 2 I=1,32 — only the first 32 ISP entries are searched (a live quirk: codes 21,22,211,212 as SITE
# species fall through to illegal ⇒ ADJMAI=0). Replicated faithfully.
@inline function bm_adjmai(inspec::Integer, sindex::Real)::Float32
    si = Float32(sindex)
    inspec >= 300 && return _bm_adjmai_group(8, si)
    grp = 0
    @inbounds for i in 1:32
        if inspec == _BM_ADJ_ISP[i]; grp = Int(_BM_ADJ_IMAP[i]); break; end
    end
    grp == 0 && return 0f0
    return _bm_adjmai_group(grp, si)
end

# maical.f RMAI — stand-scalar adjusted MAI for the site species (default DF/idx3; SITEAR default 140), capped 128.
@inline function bm_rmai(p)::Float32
    isisp = Int(p.site_species); isisp == 0 && (isisp = 3)
    (isisp < 1 || isisp > 18) && (isisp = 3)
    sssi = p.sp_site_index[isisp]; sssi == 0f0 && (sssi = 140f0)
    rmai = bm_adjmai(_BM_ISPNUM[isisp], sssi)
    rmai > 128f0 && (rmai = 128f0)
    return rmai
end

@inline function bm_dubscr(rng, sp::Integer, d, h, ba, tpccf, avh, rmai)::Float32
    hf = Float32(h); hf <= 0f0 && (hf = 0.1f0)
    cr = BM_BCR2[sp]*hf + BM_BCR1[sp]*Float32(d) + BM_BCR5[sp]*Float32(tpccf) +
         BM_BCR6[sp]*(Float32(avh)/hf) + BM_BCR8[sp]*Float32(avh) + BM_BCR3[sp]*Float32(ba) +
         BM_BCR9[sp]*(Float32(ba)*Float32(tpccf)) + BM_BCR10[sp]*Float32(rmai) + BM_BCR0[sp]
    sd = BM_DUBSD[sp]
    fcr = 0f0
    while true                                             # dubscr.f label 10: reject |FCR|>SD
        fcr = bachlo(rng, 0f0, sd)
        abs(fcr) > sd && continue
        break
    end
    if sp == 13 || sp == 14 || sp == 16 || sp == 18        # CASE(13,14,16,18): linear rescale
        cr = cr + fcr
        cr = ((cr - 1f0)*10f0 + 1f0)/100f0
    else                                                   # CASE(1:12,15,17): logistic
        abs(cr + fcr) >= 86f0 && (cr = 86f0)               # faithful: sets +86 regardless of sign
        cr = 1f0/(1f0 + exp(cr + fcr))
    end
    cr < 0.05f0 && (cr = 0.05f0); cr > 0.95f0 && (cr = 0.95f0)
    return cr
end

# CRATET's :164-166 IND over the FULL inventory record set (live 1:n AND the not-yet-deleted dead n+1:n+ndead):
# IND=IND1 — species-major, each species in LNKCHN (read) order with live and dead interleaved (intree.f:629,
# lnkchn.f, setup.f) — then RDPSRT(ITRN,DBH,IND,.FALSE.) on the REAL dbh. This is the order the CRATET :195 DENSE
# (LBKDEN backdating pass) walks for its PCTILE percentile (dense.f:241-244) and its AVH (dense.f:285-297) — i.e.
# the calibration PCT the :630 DGDRIV reads and the AVH the :610 CROWN→DUBSCR crown dub reads. `dbhv` = real dbh of
# records 1:ntot. Without the ingest read order (s.calib.input_seq) falls back to IND1(live) + dead by index.
function bm_cratet166_ind(s::StandState, dbhv::AbstractVector{Float32}, nlive::Int, ntot::Int)
    t = s.trees; iseq = s.calib.input_seq
    ord = Vector{Int32}(undef, ntot)
    k = 0
    if length(iseq) == ntot
        @inbounds for sp_o in 1:MAXSP
            mem = Int32[j for j in 1:ntot if Int(t.species[j]) == sp_o]
            sort!(mem; by = j -> iseq[j])
            for j in mem; k += 1; ord[k] = j; end
        end
    else
        nsave = t.n; t.n = nlive
        species_sort!(s)
        isct_c = s.control.sp_count_tab; ind1_c = s.scratch.idx1
        @inbounds for sp_o in 1:MAXSP
            if isct_c[sp_o, 1] != 0
                for i3 in Int(isct_c[sp_o, 1]):Int(isct_c[sp_o, 2]); k += 1; ord[k] = ind1_c[i3]; end
            end
            for j in (nlive + 1):ntot
                Int(t.species[j]) == sp_o && (k += 1; ord[k] = Int32(j))
            end
        end
        t.n = nsave
    end
    _rdpsrt!(dbhv, ord; lseq = false)
    return ord
end

# CRATET's cycle-0 IND (bm/cratet.f:164-166, :199, :270): IND is SEEDED from IND1 (species-major SPESRT order) and
# re-sorted by DBH with RDPSRT(.FALSE.) — Scowen's UNSTABLE quicksort, so equal-DBH ties fall by the seed order. Only
# when the inventory has standing-dead records (IREC2<MAXTP1) does CRATET later re-sort RDPSRT(.TRUE.) from identity.
# That IND feeds CROWN's ISORT (cratet.f:610) and the initial DENSE→PCTILE PCT (cratet.f:692) that the first cycle's
# DGDRIV reads. A fresh identity sort inverted a tied pair on 171243999020004 (two PP 11.7"/56' records: PCT
# 22.24/17.94 swapped ⇒ BAL swapped ⇒ DG 0.9602/0.8192 vs jl 0.9412/0.8359).
function bm_cratet_ind!(s::StandState, idx::AbstractVector{Int32})
    t = s.trees; n = t.n
    dbhv = view(t.dbh, 1:n)
    if t.ndead > 0
        _rdpsrt!(dbhv, idx)                                   # cratet.f:270 RDPSRT(ITRN,DBH,IND,.TRUE.)
    else
        idx .= bm_cratet166_ind(s, dbhv, n, n)               # cratet.f:164-166 IND=IND1; RDPSRT(.FALSE.)
    end
    return idx
end

# bm/crown.f — Weibull crown-ratio (all species; small trees D<1 → REGENT). RELSDI=SDIAC/SDIDEF,
# ACRNEW=C0+C1·RELSDI·100, Weibull A/B/C (B<1→1, C<2→2), SCALE=1−0.00167·(RELDEN−100), rank-based X.
function crown_ratio_update!(s::StandState, ::BlueMountains; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    relden = p.relative_density; sdiac = crown_sdi
    p_pccf = s.density.point_ccf
    rmai = lstart ? bm_rmai(p) : 0f0                        # maical RMAI (stand scalar) for the DUBSCR crown dub
    # ISORT(IND(JJ)) = ITRN−JJ+1 over FVS's IND (bm/crown.f:172-175). Cycling: IND is the gradd.f:186 / esnutr.f:325
    # RDPSRT(DBH,.TRUE.) of the ALREADY-GROWN DBH (UPDATE precedes CROWN; jl applies DBH before crown too), so the
    # key is t.dbh — the old dbh+DG/BARK re-added this cycle's growth. LSTART: CRATET's IND (bm_cratet_ind!).
    idx = Vector{Int32}(undef, n)
    if lstart
        bm_cratet_ind!(s, idx)
    else
        _rdpsrt!(view(t.dbh, 1:n), idx)
    end
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    # Visit order = bm/crown.f:185-213 `DO 70 ISPC=1,MAXSP; DO 60 I3=ISCT(ISPC,1),ISCT(ISPC,2); I=IND1(I3)`
    # (species-major, IND1). Only the LSTART DUBSCR dub draws RNG (one rejection-bounded BACHLO per D<1
    # missing-crown tree), so the order fixes which seedling gets which draw; record order handed a
    # mixed-species seedling cohort the wrong draws (171243999020004: PP 974-TPA seedling took DF's −0.68
    # ⇒ CR 32 vs live 17). Cycling (no draws, order-independent) keeps record order.
    order = if lstart
        species_sort!(s)
        isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
        Int[Int(ind1[i3]) for sp_o in 1:MAXSP
            for i3 in (isct[sp_o, 1] == 0 ? (1:0) : (Int(isct[sp_o, 1]):Int(isct[sp_o, 2])))]
    else
        collect(1:n)
    end
    @inbounds for i in order
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (lstart && t.crown_pct[i] > 0) && continue
        if d < 1f0 && lstart                               # bm/crown.f:336 label 58 — D<1 missing-CR at LSTART → DUBSCR
            pt = Int(t.plot_id[i])
            tpccf = (1 <= pt <= length(p_pccf)) ? p_pccf[pt] : 0f0
            cr = bm_dubscr(s.rng, sp, d, h, p.basal_area, tpccf, p.avg_height, rmai)
            icri = trunc(Int, cr*100f0 + 0.5f0)
            (d >= BM_CR_DLOW[sp] && d <= BM_CR_DHI[sp]) && (icri = trunc(Int, Float32(icri) * BM_CRNMLT[sp]))
            icri > 95 && (icri = 95)
            (icri < 10 && BM_CRNMLT[sp] == 1f0) && (icri = 10)
            icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri)
            continue
        end
        icr = Int(t.crown_pct[i])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = BM_CRC0[sp] + BM_CRC1[sp] * relsdi * 100f0
        A = BM_WEIBA[sp]
        # bm/crown.f:202-206 — B floor is 1.0 only for WJ(6)/LM(12)/AS(15); every other species floors at 3.0.
        B = BM_WEIBB0[sp] + BM_WEIBB1[sp] * acrnew
        bfloor = (sp == 6 || sp == 12 || sp == 15) ? 1f0 : 3f0
        B < bfloor && (B = bfloor)
        C = BM_WEIBC0[sp] + BM_WEIBC1[sp] * acrnew; C < 2f0 && (C = 2f0)
        scale = 1f0 - 0.00167f0 * (relden - 100f0)
        scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
        x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
        x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
        crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = (d >= BM_CR_DLOW[sp] && d <= BM_CR_DHI[sp]) ? Float32(icr) + chg * BM_CRNMLT[sp] :
                    Float32(icr) + chg
        end
        icri = trunc(Int, crnew + 0.5f0)
        if lstart || icr == 0
            (d >= BM_CR_DLOW[sp] && d <= BM_CR_DHI[sp]) && (icri = trunc(Int, Float32(icri) * BM_CRNMLT[sp]))
        else
            # CRMAX cap (bm/crown.f:301-314): crown length can't exceed the height-growth-adjusted max.
            crln = h * Float32(icr) / 100f0; htg = t.ht_growth[i]
            crmax = (crln + htg) / (h + htg) * 100f0
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
            (icri < 10 && BM_CRNMLT[sp] == 1f0) && (icri = trunc(Int, crmax + 0.5f0))
        end
        # final clamps (bm/crown.f:347-349)
        icri > 95 && (icri = 95)
        (icri < 10 && BM_CRNMLT[sp] == 1f0) && (icri = 10)
        icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    # ---- CYCLE-0 DEAD-TREE CROWN DUB (bm/crown.f:360-380 `DO 79 I=IREC2,MAXTRE`) ----
    # After the live loop FVS dubs MISSING crowns on the standing-dead records with DUBSCR (one main-stream
    # rejection-bounded BACHLO each, regardless of DBH), iterating IREC2→MAXTRE = the REVERSE of jl's dead
    # storage (FVS files dead from MAXTRE downward in read order; same layout as ie/crown.f DO 79). TPCCF =
    # PCCF(ITRE(I)); BA/AVH/RMAI are the same stand scalars the live dub saw. Broken-top (ITRUNC>0) crowns are
    # re-expressed on the normal height (crown.f:373-376); bounds [10,95] for all species (crown.f:378-379).
    if lstart && t.ndead > 0
        @inbounds for i in (n + Int(t.ndead)):-1:(n + 1)
            Int(t.crown_pct[i]) > 0 && continue
            sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
            pt = Int(t.plot_id[i])
            tpccf = (1 <= pt <= length(p_pccf)) ? p_pccf[pt] : 0f0
            cr = bm_dubscr(s.rng, sp, d, h, p.basal_area, tpccf, p.avg_height, rmai)
            icri = trunc(Int, cr*100f0 + 0.5f0)
            if t.trunc[i] != 0
                hn = Float32(t.norm_ht[i]) / 100f0
                hd = hn - Float32(t.trunc[i]) / 100f0
                cl = (Float32(icri) / 100f0) * hn - hd
                icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
            end
            icri > 95 && (icri = 95); icri < 10 && (icri = 10)
            t.crown_pct[i] = Int32(icri)
        end
    end
    return s
end

# bm/cratet.f LSTART crown-init: DENSE runs over the FULL inventory (live + HISTORY 6-9 standing-dead records),
# and that dead-inclusive density (BA/AVH/point-CCF) is what CROWN→DUBSCR sees when dubbing the D<1 / missing-CR
# LIVE trees. jl partitions the dead into t.n+1:t.n+ndead, so it computes the dead-inclusive scalars by temporarily
# extending the live range, then restores live-only (the grow cycle recomputes density before use). Measured on
# 504443988126144: dead-inclusive BA 55.56 / AVH 85.07 / TPCCF 84.3 (vs live-only 35/45/63) — matches live DUBSCR.
function bm_crown_init_lstart!(s::StandState)
    t = s.trees
    nlive = t.n
    # bm/cratet.f:189-195 `LBKDEN = IDG.LT.2; CALL DENSE` — CROWN (cratet.f:610) dubs against THAT density, whose
    # live WK3 is BACKDATED to the start of the measured-growth period (dense.f:70-128: measured-DG trees by their
    # own increment, the rest by the stand-average BAGR) whenever ≥1 live record carries a measured increment
    # (SN>0; else WK3=DBH). Not backdating put the DUBSCR BA/TPCCF at the current-DBH values (1127576412290487:
    # 3 LP with past DBH ⇒ live BA 16.486 / TPCCF 80.04 vs jl 48.25 / 199.2 ⇒ every seedling crown mis-dubbed).
    # _backdate_dbh! is the shared IDG-faithful dense.f port (live 1:t.n, in place); restored below.
    lbkden = s.control.growth_idg < 2
    saved_live = lbkden ? t.dbh[1:nlive] : Float32[]
    # #151: dense.f:83-87 — in the CRATET backdating DENSE, standing-dead records get WK3=DBH EXCEPT IMC(I)==9
    # (HISTORY 8,9, older-dead) which LOAD DBH=0 ⇒ they add 0 to BA/CCF/SDI while their HEIGHT still counts
    # toward AVH (measured on 449747082489998: live DUBSCR AVH 67.34 vs jl 1.01 without the dead heights).
    # notre.f:119-124 expands EVERY inventory-dead record (IREC2..MAXTRE) with VP/FP/FP2×(FINT/FINTM), so the dead
    # PROB that CRATET's DENSE sums is TREE_COUNT×FINT/FINTM (BM 10/5 ⇒ ×2; FVS divides it back only where a true
    # density is needed — snag init cratet.f:565, dbstrls MortPA). jl's shared notre! leaves dead TPA unscaled, so
    # scale them here for this pass (24001521010900: 5 HISTORY-8 dead at PROB 12.036/1.998 ⇒ TPCCF 77.788→77.813).
    saved = Tuple{Int,Float32}[]
    saved_tpa = Float32[]
    if t.ndead > 0
        t.n = nlive + Int(t.ndead)
        fintr = s.control.growth_fintm > 0f0 ? s.control.growth_fint / s.control.growth_fintm : 1f0
        @inbounds for i in (nlive + 1):(nlive + Int(t.ndead)); push!(saved_tpa, t.tpa[i]); t.tpa[i] *= fintr; end
    end
    # AVHT40 top height from REAL DBH/HT over live + all dead (IND = real-DBH sort, not WK3), with the dead at their
    # notre-expanded PROB (24001521010900: the ×2 dead fill the top 40 TPA ⇒ AVH 48.466 = live, 41.596 unscaled).
    # dense.f:285-297 AVH walk of the :195 CRATET DENSE, over the :164-166 IND (bm_cratet166_ind: LNKCHN-seeded,
    # dead included, RDPSRT(.FALSE.) on real dbh) — the AVH the :610 CROWN→DUBSCR dub reads. The generic double
    # sort put different equal-DBH records at the 40-TPA boundary (23899355010900 AVH 72.567 vs live 71.606).
    # That DENSE runs BEFORE the missing-height dub (DO 130 :363, DO 145 :464) ⇒ heights are as READ (missing = 0):
    # use the pre-dub snapshot cratet_ht_in (41135212010497: dubbed 1.01 seedlings ⇒ AVH 33.671 vs live 33.160).
    avht_real = let ntot = t.n, ord = bm_cratet166_ind(s, view(t.dbh, 1:t.n), nlive, t.n), hin = s.calib.cratet_ht_in
        use_hin = length(hin) == ntot
        avh = 0f0; ssumn = 0f0
        @inbounds for k in 1:ntot
            ii = Int(ord[k]); p = t.tpa[ii]
            ssumn + p > 40f0 && (p = 40f0 - ssumn)
            ssumn += p; avh += (use_hin ? hin[ii] : t.height[ii]) * p
            ssumn >= 40f0 && break
        end
        ssumn > 0f0 ? avh / ssumn : 0f0
    end
    if lbkden                                 # backdate LIVE WK3 only (after the real-DBH AVH ranking)
        t.n = nlive; _backdate_dbh!(s); t.n = nlive + Int(t.ndead)
    end
    if t.ndead > 0
        @inbounds for i in (nlive + 1):(nlive + Int(t.ndead))
            (t.history[i] == 8 || t.history[i] == 9) || continue
            push!(saved, (i, t.dbh[i])); t.dbh[i] = 0f0
        end
    end
    compute_density!(s)                    # CRATET DENSE: backdated live (+ dead-inclusive) BA / point-CCF
    s.calib.cratet_relden = stand_ccf(s)   # RELDEN after cratet.f:195 DENSE (backdated, dead-inclusive) → REGENT HCOR cal
    @inbounds for (i, d) in saved; t.dbh[i] = d; end
    @inbounds for (k, i) in enumerate((nlive + 1):(nlive + length(saved_tpa))); t.tpa[i] = saved_tpa[k]; end
    t.n = nlive
    lbkden && @inbounds(for i in 1:nlive; t.dbh[i] = saved_live[i]; end)
    s.plot.avg_height = avht_real
    crown_ratio_update!(s, s.variant; lstart = true)   # DUBSCR-dub live D<1 seedlings + Weibull-dub missing-CR overstory
    compute_density!(s; cratet_ind = true)  # restore live-only density; CRATET IND ⇒ cycle-0 PCT/AVH (cratet.f:692 DENSE)
    return s
end
