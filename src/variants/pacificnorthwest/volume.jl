# =============================================================================
# volume.jl (pacificnorthwest) — PN volume (R6 NVEL). Chunk 8.
#
# PN forest-612 (Siuslaw) VOLEQ (dumped from FVSpn_clean): DF → F00FW2W202 (westside Flewelling, geosub
# '00' = base, no regional override); WH → F03FW2W263; everything else → Behre 616BEHW<fia> (PN does NOT
# use the INGY overrides WC had for GF/NF/IC; RA → NVBM240351, not reachable on pnt01). Reuses the WC
# westside SHP_W3/W4/W5 (wc_fw2_westside_vol handles geosub '00') + BM's R6 Behre machinery + PN form class
# (formcl_pn.csv). Merch (pn/sitset.f non-BLM) = the WC westside defaults (TOPD=BFTOPD=4.5, DBHMIN=7, LP=6).
# =============================================================================

const PN_FORMCL = let
    T = Dict{Tuple{Int,Int,Int},Int}()
    for l in readlines(joinpath(PN_DATADIR, "formcl_pn.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        T[(parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]))] = parse(Int, f[4])
    end
    A = zeros(Int, 10, 39, 5)
    for ((fi, sp, dc), v) in T
        if fi <= 6; A[fi, sp, dc] = v else; for c in 1:5; A[fi, sp, c] = v; end; end
    end
    A
end
function pn_formcl(sp::Integer, ifor::Int, d::Real)::Int
    (sp < 1 || sp > 39) && return 80
    (ifor < 1 || ifor > 10) && (ifor = 2)
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1); Float32(d) > 40.9f0 && (ifcdbh = 5)
    v = PN_FORMCL[ifor, sp, ifcdbh]; v == 0 ? 80 : v
end

# PN R6_EQN VOLEQ (forest 612 → FORNUM 12, validated vs FVSpn_clean). DF westside geosub '00'; else Behre.
function _pn_r6_eqn(fornum::Int, fia::Int)::String
    fia == 202 && return "F00FW2W202"           # DF → westside Flewelling (geosub 00 = base)
    fia == 263 && return "F03FW2W263"           # WH → westside Flewelling
    fia == 351 && return "NVBM240351"           # RA → Brackett (not reachable on pnt01)
    return "616BEHW" * lpad(string(fia), 3, '0')
end

# PN Behre per-tree volume — reuse the BM R6 machinery + PN form class.
function pn_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32)
    fclass = pn_formcl(sp, ifor, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    vol2 = 0f0; vol4 = 0f0
    v1 = if h <= 17.3f0
        0.00272708f0 * dbhib * dbhib * h
    else
        v = bm_r6vol3(d, dbtbh, fclass, h, 1)
        mtopp = 4.5f0 * bark
        xlogs, ld1 = bm_r6dibs(d, fclass, mtopp, h)
        lv1, lv4 = bm_r6vol1(d, fclass, xlogs, ld1)
        nlog = Int(floor(xlogs)); nacc = (xlogs - nlog) > 0f0 ? nlog + 1 : nlog
        for k in 1:nacc
            vol2 += bm_anint(lv1[k]); vol4 += bm_anint(lv4[k] * 10f0) / 10f0
        end
        v
    end
    return (max(v1, 0f0), max(vol4, 0f0), max(vol2, 0f0))
end

function compute_volumes_pn!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species
    ifor = Int(s.plot.forest_idx)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 39
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = wc_bratio(sd, sp, d)
        dbhmin = sp == 11 ? 6.0f0 : 7.0f0
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, bf = wc_fw2_westside_vol(eq, d, hv, bark)
        else
            tcf, mcf, bf = pn_behre_vol(sp, ifor, d, hv, bark)   # 616BEHW (NVBM RA deferred)
        end
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = d >= dbhmin ? max(bf, 0f0) : 0f0
    end
    return s
end
