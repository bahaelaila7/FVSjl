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

"""Kozak taper exponent (log.f EX statement function): D4·D**2+D5·ALOG(D+0.001)+D6·SQRT(D)+D7·E+D8·EXP(D) —
ALOG/EXP are gfortran's logf/expf (glibc, `flog`/`fexp`), D**2 = D·D."""
@inline bc_vol_ex(b4,b5,b6,b7,b8,d,e) = b4*(d*d) + b5*flog(d+0.001f0) + b6*sqrt(d) + b7*e + b8*fexp(d)

"""dib at relative height y (log.f): FF·X**EXPO, X = (1−SQRT(Y))/PER, REAL**REAL = powf."""
@inline function _bc_dib(y::Float32, ff::Float32, per::Float32, dh::Float32, B)
    x = (1f0 - sqrt(y)) / per
    expo = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], y, dh)
    return ff * fpow(x, expo)
end

"""Newton/Simpson-composite volume of a stem segment [hl,hu] (log.f FUNCTION VLM), dib endpoints d1,d2 (cm)."""
function bc_vol_vlm(hl, hu, d1, d2, ff, per, dh, B, htm)
    dis = hu - hl
    dis < 0.01f0 && return 0f0
    n = dis > 5.01f0 ? 6 : 4
    d = zeros(Float32, 7); d[1] = d1; d[n+1] = d2
    ds = dis / n; x1 = hl
    @inbounds for i in 2:n
        x1 += ds; y = x1 / htm
        d[i] = _bc_dib(y, ff, per, dh, B)
    end
    @inbounds for i in 1:(n+1); d[i] *= d[i]; end
    n == 4 ? BC_VLM_CONS*ds*(d[1]/2f0 + 2f0*d[2] + d[3] + 2f0*d[4] + d[5]/2f0) :
             BC_VLM_CONS*ds*(d[1]/2f0 + 2f0*d[2] + d[3] + 2f0*d[4] + d[5] + 2f0*d[6] + d[7]/2f0)
end

"""
    bc_vol_log(is, dbh_cm, htm, shm, td_cm, htrunc_m) -> (GROS, TVOL) (m³)

Port of canada/bc/log.f SUBROUTINE LOG (Kozak taper integration). `is` = Kozak species; all lengths metric.
GROS = merch logs (TVOL) + stump (STMV) + top (TOPV). DI3 is the log.f variable: the stump dib DBT, re-assigned to the
0.3-m dib when the stump height is under 0.3 m (log.f:260-271) — which the stump volume (log.f:304) then uses.
"""
function bc_vol_log(is::Int, dbh::Float32, htm::Float32, shm::Float32, td::Float32, htrunc::Float32)
    (htm <= 0f0 || dbh <= 0f0) && return (0f0, 0f0)
    c = BC_TAPER[is]; B = c.A
    zzz = c.PP; per = 1f0 - sqrt(zzz)
    dh = dbh / htm
    ff = B[1] * fpow(dbh, B[2]) * fpow(B[3], dbh)            # FF=B(1)*DBH**B(2)*B(3)**DBH
    dbt = _bc_dib(shm / htm, ff, per, dh, B)                  # dib at stump
    tvol = 0f0; nl = 0; di3 = dbt; hm = shm
    diam = 0f0
    htrunc < htm && (diam = _bc_dib(htrunc / htm, ff, per, dh, B))   # top-kill diameter
    if dbt > td                                               # log.f:239 IF (DBT.LE.TD) GO TO 10
        x = 0.9f0; nn = 0
        while true                                            # log.f:3 merch-height iteration
            expo = bc_vol_ex(B[4],B[5],B[6],B[7],B[8], x, dh)
            yv = fpow(td/ff, 1f0/expo); yv = (1f0 - yv*per)^2
            abs(x - yv) < 0.0001f0 && break
            x = x + (yv - x)/2f0
            x > 1f0 && (x = (shm/htm)*(di3/td))
            nn += 1; nn > 9 && break
        end
        hm = x*htm; hm < shm && (hm = shm + 0.01f0)
        htrunc < hm && (hm = htrunc)
        nl = Int(unsafe_trunc(Int16, (hm - shm)/BC_VOL_GOL + 1f0))   # INT(…,2)
        x1 = shm; k1 = 1; tdlprev = 0f0; skip_main = false
        if x1 < 0.3f0                                         # log.f:258 IF (X1.GE.0.3) GO TO 9
            x1 += BC_VOL_GOL; x2 = 0.3f0; x1 > hm && (x1 = hm); k1 = 2
            di3 = _bc_dib(x2/htm, ff, per, dh, B)
            tdl1 = _bc_dib(x1/htm, ff, per, dh, B)
            tvol += BC_VOL_CONS*(di3*di3)*(0.3f0 - shm) + bc_vol_vlm(x2, x1, di3, tdl1, ff, per, dh, B, htm)
            tdlprev = tdl1
            nl < 2 && (skip_main = true)                     # IF (NL.LT.2) GO TO 15
        end
        if !skip_main
            for i in k1:nl
                x1 += BC_VOL_GOL; x2 = x1 - BC_VOL_GOL; x1 > hm && (x1 = hm)
                tdli = _bc_dib(x1/htm, ff, per, dh, B)
                d1v = i < 2 ? (shm >= 0.3f0 ? dbt : di3) : tdlprev   # log.f:290-294
                tvol += bc_vol_vlm(x2, x1, d1v, tdli, ff, per, dh, B, htm)
                tdlprev = tdli
            end
        end
    end
    # stump (log.f:298-305)
    if shm <= 0.3f0
        stmv = BC_VOL_CONS*(di3*di3)*shm
    else
        di3 = _bc_dib(0.3f0/htm, ff, per, dh, B)
        stmv = BC_VOL_CONS*(di3*di3)*0.3f0 + bc_vol_vlm(0.3f0, shm, di3, dbt, ff, per, dh, B, htm)
    end
    dd = di3 < td ? di3 : td
    topv = hm != htrunc ? bc_vol_vlm(hm, htrunc, dd, diam, ff, per, dh, B, htm) : 0f0
    return tvol + stmv + topv, tvol            # (GROS total, TVOL merch-logs)
end

"""
    bc_tree_vol(sp, d, h, itht) -> (total_cuft, merch_cuft)

BC cubic volume for one tree (VOLS→CFVOL→MIN→LOG). `d` in, `h` ft, `itht` topkill (cm·100). Merch (VM)
is the merch-log volume, gated by D≥DBHMIN and D≥TOPD (cfvol.f:55); else 0. Total (VN) always computed.
"""
function bc_tree_vol(sp::Integer, d::Real, h::Real, itht::Integer)
    h < 4.5 && return (0f0, 0f0)
    is = BC_SPTR[Int(sp)]; is < 1 && return (0f0, 0f0)
    htrunc_ft = itht > 0 ? Float32(itht)/100f0 : Float32(h)   # cfvol.f HTRUNC=ITHT/100, 0 ⇒ H
    # CFVOL passes the grinit.f imperial thresholds STMP = 30·CMtoFT (ft), TOPD = 10·CMtoIN (in), DBHMIN = 17.5|12.5·CMtoIN;
    # MIN (min.f) converts them back with FTtoM / INtoCM — SH = 0.29999998 m, TD = 9.99998 cm, not 0.3 / 10.
    stmp_ft = bc_vol_stmp(sp) * BC_CMtoFT; topd_in = bc_vol_topd(sp) * BC_CMtoIN
    dbh_cm = Float32(d) * BC_INtoCM                        # MIN metric conversion (min.f:114-118)
    ht_m   = Float32(h) * BC_FTtoM
    sh_m   = stmp_ft * BC_FTtoM
    td_cm  = topd_in * BC_INtoCM
    htr_m  = htrunc_ft * BC_FTtoM
    gros, tvol = bc_vol_log(is, dbh_cm, ht_m, sh_m, td_cm, htr_m)   # m³
    vn = gros * BC_M3toFT3
    # merch gate (cfvol.f:55): D ≥ DBHMIN(ISPC) .AND. D ≥ TOPD(ISPC), compared in INCHES (a 12.5-cm PL is D = 12.5·CMtoIN
    # = DBHMIN exactly, merchantable; the metric compare 4.92125·2.54 = 12.49997 < 12.5 dropped it).
    vm = (Float32(d) >= bc_vol_dbhmin(sp) * BC_CMtoIN && Float32(d) >= topd_in) ? tvol * BC_M3toFT3 : 0f0
    return (vn, vm)
end

bc_tree_cuft(sp, d, h, itht) = bc_tree_vol(sp, d, h, itht)[1]   # total-only (validation helper)

"""Fill `t.cuft_vol` (total) + `t.merch_cuft_vol` (merch) for all trees; BC Kozak taper. Board deferred (0)."""
function compute_volumes!(s::StandState, ::BritishColumbia)
    t = s.trees
    @inbounds for i in 1:t.n
        d = t.dbh[i]; h = t.height[i]
        if d <= 0f0 || h <= 0f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
        else
            # bc/vols.f:137-138: a top-killed tree (H≥4.5, ITRUNC>0) is volumed at its NORMAL height NORMHT
            # (CRATET's dubbed full height), truncated at ITRUNC — not at the recorded HT.
            (h >= 4.5f0 && t.trunc[i] > 0) && (h = Float32(t.norm_ht[i]) / 100f0)
            vn, vm = bc_tree_vol(Int(t.species[i]), d, h, Int(t.trunc[i]))
            t.cuft_vol[i] = vn; t.merch_cuft_vol[i] = vm
        end
    end
    return s
end
