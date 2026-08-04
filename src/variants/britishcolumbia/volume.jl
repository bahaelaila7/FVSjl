# =============================================================================
# volume.jl (britishcolumbia) — BC total cubic volume (canada/bc: VOLS→CFVOL→MIN→LOG). CHUNK 8.
#
# BC volume is a BC-specific Kozak-2002 TAPER model, NOT the shared NVEL driver. LOG integrates the
# taper dib(h)=FF·X^EXPO in METRIC (m³); MIN wraps it (imperial↔metric); CFVOL gates merch. all_BC
# exercises TOTAL cubic only (merch=0). Coefficients: data/britishcolumbia/vol_taper_iz2.csv (the IZ=2
# / FIZ=2 column, the only one BC uses) + vol_sptr.csv (BC species → Kozak species). Thresholds from
# grinit.f: STMP=30cm, TOPD=10cm, DBHMIN=17.5cm (12.5 for PL sp7). ⚠ merch/board (BFVOL) TODO.
# =============================================================================

const BC_M3toFT3 = 35.314455f0
const BC_VOL_CONS  = 0.00007854f0   # log.f:188 (stump element π/4·1e-4)
const BC_VLM_CONS  = 0.00005236f0   # log.f:332 (Simpson composite in VLM)
const BC_VOL_GOL   = 3.0f0          # log length (m)

struct BCTaperCoef
    A::NTuple{8,Float32}   # A1..A8
    PP::Float32
end

function _bc_load_taper()
    rows = _bc_readcsv("vol_taper_iz2.csv"); h = rows[1]
    ci(n) = findfirst(==(n), h)
    out = BCTaperCoef[]
    for r in rows[2:end]
        isempty(strip(r[1])) && continue
        a = ntuple(k -> parse(Float32, strip(r[ci("A$k")])), 8)
        push!(out, BCTaperCoef(a, parse(Float32, strip(r[ci("PP")]))))
    end
    out
end
const BC_TAPER = _bc_load_taper()               # indexed by Kozak species (1..16)

function _bc_load_sptr()
    rows = _bc_readcsv("vol_sptr.csv"); h = rows[1]
    ci(n) = findfirst(==(n), h)
    sptr = zeros(Int, 15)
    for r in rows[2:end]
        isempty(strip(r[1])) && continue
        sptr[parse(Int, strip(r[ci("bc_sp")]))] = parse(Int, strip(r[ci("kozak_sp")]))
    end
    sptr
end
const BC_SPTR = _bc_load_sptr()                 # BC species (1..15) → Kozak taper species

# per-species merch thresholds (grinit.f:99-142), stored METRIC (cm/cm) — converted in cfvol like the oracle.
bc_vol_stmp(sp) = 30.0f0                          # stump height (cm)
bc_vol_topd(sp) = 10.0f0                          # top diameter (cm)
bc_vol_dbhmin(sp) = sp == 7 ? 12.5f0 : 17.5f0    # min merch DBH (cm); PL=12.5

"""Kozak taper exponent (log.f EX): D4·D²+D5·ln(D+0.001)+D6·√D+D7·E+D8·exp(D)."""
@inline bc_vol_ex(b4,b5,b6,b7,b8,d,e) = b4*d*d + b5*log(d+0.001f0) + b6*sqrt(d) + b7*e + b8*exp(d)

"""Newton/Simpson-composite volume of a stem segment [hl,hu] (log.f VLM), dib endpoints d1,d2 (cm)."""
function bc_vol_vlm(hl, hu, d1, d2, ff, per, dh, B, htm)
    dis = hu - hl
    dis < 0.01f0 && return 0f0
    n = dis > 5.01f0 ? 6 : 4
    d = zeros(Float32, 7); d[1] = d1; d[n+1] = d2
    ds = dis / n; x1 = hl
    @inbounds for i in 2:n
        x1 += ds; y = x1 / htm
        expo = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], y, dh)
        d[i] = ff * ((1f0 - sqrt(y))/per) ^ expo
    end
    @inbounds for i in 1:(n+1); d[i] *= d[i]; end
    n == 4 ? BC_VLM_CONS*ds*(d[1]/2f0 + 2f0*d[2] + d[3] + 2f0*d[4] + d[5]/2f0) :
             BC_VLM_CONS*ds*(d[1]/2f0 + 2f0*d[2] + d[3] + 2f0*d[4] + d[5] + 2f0*d[6] + d[7]/2f0)
end

"""
    bc_vol_log(is, dbh_cm, htm, shm, td_cm, htrunc_m) -> total volume GROS (m³)

Port of canada/bc/log.f (Kozak taper integration). `is` = Kozak species; all lengths metric.
Returns GROS = merch-log volume (TVOL) + stump (STMV) + top (TOPV). Total-cubic path (merch is
computed but for all_BC the tree falls below TD ⇒ TVOL=0, GROS=STMV+TOPV = whole stem).
"""
function bc_vol_log(is::Int, dbh::Float32, htm::Float32, shm::Float32, td::Float32, htrunc::Float32)
    (htm <= 0f0 || dbh <= 0f0) && return 0f0
    c = BC_TAPER[is]; B = c.A
    zzz = c.PP; per = 1f0 - sqrt(zzz)
    dh = dbh / htm
    ff = B[1] * dbh^B[2] * B[3]^dbh
    y = shm / htm; x = (1f0 - sqrt(y))/per
    expo = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], y, dh)
    dbt = ff * x^expo                              # dib at stump
    tvol = 0f0; nl = 0; di3 = dbt; hm = shm
    diam = 0f0
    if htrunc < htm                                # top-kill diameter
        y = htrunc/htm; x = (1f0 - sqrt(y))/per
        expo = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], y, dh)
        diam = ff * x^expo
    end
    if dbt > td                                    # MERCH present: find merch height + accumulate logs
        x = 0.9f0; nn = 0
        while true
            expo = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], x, dh)
            yv = (td/ff)^(1f0/expo); yv = (1f0 - yv*per)^2
            abs(x - yv) < 0.0001f0 && break
            x = x + (yv - x)/2f0
            x > 1f0 && (x = (shm/htm)*(di3/td))
            nn += 1; nn > 9 && break
        end
        hmm = x*htm; hmm < shm && (hmm = shm + 0.01f0)
        htrunc < hmm && (hmm = htrunc)
        hm = hmm
        nl = trunc(Int, (hm - shm)/BC_VOL_GOL + 1f0)
        # log accumulation (log.f:263-299). K1=1 unless a partial first log below 0.3m is peeled off.
        x1 = shm; k1 = 1; tdlprev = 0f0; skip_main = false
        if x1 < 0.3f0
            x1 += BC_VOL_GOL; x2 = 0.3f0; x1 > hm && (x1 = hm); k1 = 2
            yv = x2/htm; e1 = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], yv, dh); di3b = ff*((1f0 - sqrt(yv))/per)^e1
            yv = x1/htm; e2 = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], yv, dh); tdl1 = ff*((1f0 - sqrt(yv))/per)^e2
            tvol += BC_VOL_CONS*di3b*di3b*(0.3f0 - shm) + bc_vol_vlm(x2, x1, di3b, tdl1, ff, per, dh, B, htm)
            tdlprev = tdl1
            nl < 2 && (skip_main = true)                 # log.f: IF (NL.LT.2) GO TO 15
        end
        if !skip_main
            for i in k1:nl
                x1 += BC_VOL_GOL; x2 = x1 - BC_VOL_GOL; x1 > hm && (x1 = hm)
                yv = x1/htm; e = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], yv, dh); tdli = ff*((1f0 - sqrt(yv))/per)^e
                d1v = i < 2 ? dbt : tdlprev               # log.f:292-296 (DI3=DBT here); i≥2 → TDL(i-1)
                tvol += bc_vol_vlm(x2, x1, d1v, tdli, ff, per, dh, B, htm)
                tdlprev = tdli
            end
        end
    end
    # stump (log.f:303-312)
    if shm <= 0.3f0
        stmv = BC_VOL_CONS*di3*di3*shm
    else
        yv = 0.3f0/htm; xv = (1f0 - sqrt(yv))/per
        e = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], yv, dh); di3s = ff*xv^e
        stmv = BC_VOL_CONS*di3s*di3s*0.3f0 + bc_vol_vlm(0.3f0, shm, di3s, dbt, ff, per, dh, B, htm)
    end
    dd = di3 < td ? di3 : td
    topv = hm != htrunc ? bc_vol_vlm(hm, htrunc, dd, diam, ff, per, dh, B, htm) : 0f0
    return tvol + stmv + topv
end

"""BC total cubic-foot volume for one tree (VOLS→CFVOL→MIN→LOG). `d` in, `h` ft, `itht` topkill (cm·100)."""
function bc_tree_cuft(sp::Integer, d::Real, h::Real, itht::Integer)
    h < 4.5 && return 0f0
    is = BC_SPTR[Int(sp)]; is < 1 && return 0f0
    htrunc_ft = itht > 0 ? Float32(itht)/100f0 : Float32(h)
    # MIN: convert imperial → metric for LOG (min.f:114-118)
    dbh_cm = Float32(d) * BC_INtoCM
    ht_m   = Float32(h) * BC_FTtoM
    sh_m   = (bc_vol_stmp(sp) * 0.01f0)                     # STMP cm → m
    td_cm  = bc_vol_topd(sp)
    htr_m  = htrunc_ft * BC_FTtoM
    gros = bc_vol_log(is, dbh_cm, ht_m, sh_m, td_cm, htr_m)  # m³
    return gros * BC_M3toFT3                                 # → ft³
end

"""Fill `t.cuft_vol` (total cubic ft) for all trees; BC Kozak taper. Merch/board TODO."""
function compute_volumes!(s::StandState, ::BritishColumbia)
    t = s.trees
    @inbounds for i in 1:t.n
        d = t.dbh[i]; h = t.height[i]
        t.cuft_vol[i] = (d <= 0f0 || h <= 0f0) ? 0f0 :
            bc_tree_cuft(Int(t.species[i]), d, h, Int(t.trunc[i]))
    end
    return s
end
