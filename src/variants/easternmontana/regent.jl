# =============================================================================
# regent.jl (easternmontana) — EM small-tree growth (em/regent.f, NIVAR path).
#
# Small trees (dbh < XMAX[sp]) — for emt01, the 0.1" DF seedling. Mirrors KT's subcycle scaffold
# (REGYR=5 subcycles, per-subcycle density from the large trees) with EM's NIVAR forms:
#   HTGRL = CON + BH·ln(H1) + BCCF·RDJ + BBAL·BAL   (log-space; BH=0.3740, BCCF=−0.00391, BBAL=−0.22957)
#   H2    = H1 + EXP(HTGRL)·(KPER/REGYR)·XRHGRO     (em/regent.f:550 — increment is exp(HTGRL))
#   CON   = RHCON[sp] + HCOR[sp];  RHCON = REGCH + 1.0667 + RHHAB[MAPHAB[ITYPE]]  (em_regcons!)
#   REGCH = RHGL[IGL] + (RSAB0 + RSAB1·cos(ASP) + RSAB2·sin(ASP))·SLOPE
#   DG-dub (per subcycle, D<3 & H2>4.5): D1=DIAM[sp]+DADJ (or AX·(H1−4.5)^BX+DADJ if H1>4.5);
#     D2=AX·(H2−4.5)^BX+DADJ; DGJ=max(D2−D1,0)·XRDGRO; D2=D+DGJ. AX=0.0658, BX=1.3817.
#   DADJ = DELMAX·RELH²−2·DELMAX·RELH+0.65; DELMAX=min((AH/36)·(0.01232·RELDEN−1.75),0); RELH=(H1−4.5)/(AH−4.5).
# Final: HTGR1=WK3−H + ZZRAN·HSIGMA(0.59) (dgsd≥1); XWT blend d∈[XMIN,cap] with the large-tree HTG; size cap.
# NIVAR handles the EMVAR conifers + LL(5) — LL capped at D<3 (regent.f:858) so its D≥3 large trees stay on
# the pure large-tree path (validated: em_LLseed seedlings grow, em_LL large trees bit-exact).
# CRVAR(CO)/TTVAR(LM)/UTVAR(RM,AS,PB) small-tree branches deferred.
# =============================================================================

const EM_RG_XMAX = Float32[3,3,3,3,10,99,3,3,3,3,2,4,2,2,2,2,4,3,3]
const EM_RG_XMIN = Float32[1.5,1.5,1.5,1.5,2,90,1.5,1.5,1.5,1.5,0.5,2,0.5,0.5,0.5,0.5,2,1.5,1.5]
const EM_RG_DIAM = Float32[0.4,0.3,0.3,0.3,0.3,0.3,0.4,0.3,0.3,0.5,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
const EM_RG_RHGL = Float32[-0.2785, -0.0480, 0.0]                 # geo-location (IGL 1..3)
const EM_RG_RHHAB = Float32[-0.2146, -0.0941, -0.4916, -0.3582, 0.0]
const EM_RG_MAPHAB = Int[4,4,4,4,4,4,4,4,4,4,4,4,3,2,5,5,5,5,1,4,4,4,4,4,4,4,1,3,3,4]   # ITYPE→RHHAB idx
const _EM_RG_RSAB0 = -0.10987f0; const _EM_RG_RSAB1 = 0.22157f0; const _EM_RG_RSAB2 = -0.12432f0
const _EM_RG_BH = 0.3740f0; const _EM_RG_BCCF = -0.00391f0; const _EM_RG_BBAL = -0.22957f0
const _EM_RG_AX = 0.0658f0; const _EM_RG_BX = 1.3817f0
const _EM_RG_HSIGMA = 0.59f0; const _EM_RG_REGYR = 5.0f0

# Species handled by the current NIVAR-form regent: EMVAR conifers + LL(5) (NIVAR — same HTGRL form; CON
# RHCON+HCOR==RHCON·exp(HCOR) at HCOR=0; BH/BCCF/BBAL are the shared subalpine-fir values). LM/CO/RM/AS/PB deferred.
@inline _em_rg_nivar(sp::Int) = _em_orig_species(sp) || sp == 5
# Effective regent DIAMETER cap: NIVAR skips the regent for D≥3 (em/regent.f:858 GO TO 23). The conifers
# already have EM_RG_XMAX=3, but LL(5)'s XMAX=10 is only the height-blend range — its regent still caps at 3,
# so D≥3 LL stays on the pure large-tree path (else the XWT DG-override halves the large DG → em_LL breaks).
@inline _em_rg_cap(sp::Int) = sp == 5 ? 3.0f0 : EM_RG_XMAX[sp]

# em/regent.f RCON entry: RHCON[sp] = REGCH + 1.0667 + RHHAB[MAPHAB[ITYPE]]. Returns the RHCON vector.
function em_regcons!(s::StandState)
    p = s.plot
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    igl = Int(p.geo_location); (igl < 1 || igl > 3) && (igl = 1)
    irhhab = EM_RG_MAPHAB[itype]; (irhhab < 1 || irhhab > 5) && (irhhab = 5)
    regch = EM_RG_RHGL[igl] + (_EM_RG_RSAB0 + _EM_RG_RSAB1 * cos(p.aspect) + _EM_RG_RSAB2 * sin(p.aspect)) * p.slope
    base = regch + 1.0667f0 + EM_RG_RHHAB[irhhab]
    rhcon = fill(base, 19)                                        # site constant is species-independent here
    return rhcon
end

function small_tree_growth!(s::StandState, stash, ::EasternMontana; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    n = t.n; n == 0 && return s
    rhcon = em_regcons!(s)
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height; ah = avh
    dgsd = s.control.dg_sd; regyr = _EM_RG_REGYR
    # subcycle count/lengths (regent.f:186-203, reuse KT logic)
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    # per-subcycle density from the large trees (regent.f:236-256, KT form)
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            cc1 = em_tree_ccf(sp, d1); cc2 = em_tree_ccf(sp, d2)
            bi = (b2 - b1) / 10.0f0; ci = (cc2 - cc1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * 0.985f0^k
                rdnext[j] += k * ci / pr * pn; banext[j] += k * bi * pn
            end
        end
    end
    delmax = (avh / 36.0f0) * (0.01232f0 * relden - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    # per-tree accumulators over subcycles
    wk3 = Float32[t.height[i] for i in 1:n]
    dnow = Float32[t.dbh[i] for i in 1:n]
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]; kpj = Float32(kper[j])
        for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]
            (d >= _em_rg_cap(sp) || t.tpa[i] <= 0.0f0) && continue
            _em_rg_nivar(sp) || continue                          # EMVAR conifers + LL(5)=NIVAR
            con = rhcon[sp] + c.htg_cor_small[sp]                # + HCOR (0 until calibrated)
            h1 = wk3[i]; h1 <= 0f0 && (h1 = 0.1f0)
            pct = Float32(t.crown_ratio[i]); bal = baj * (100.0f0 - pct) * 0.01f0
            htgrl = con + _EM_RG_BH * log(h1) + _EM_RG_BCCF * rdj + _EM_RG_BBAL * bal
            h2 = h1 + exp(htgrl) * (kpj / regyr)                 # XRHGRO=1 (no MULTS)
            wk3[i] = h2
            # per-subcycle DG dub (D<3, H2>4.5)
            if !(j >= nper || d >= 3.0f0) && h2 > 4.5f0
                relh = (h1 - 4.5f0) / (ah - 4.5f0); relh > 1f0 && (relh = 1f0); relh < 0f0 && (relh = 0f0)
                dadj = delmax * relh * relh - 2.0f0 * delmax * relh + 0.65f0
                d1v = h1 > 4.5f0 ? _EM_RG_AX * (h1 - 4.5f0)^_EM_RG_BX + dadj : EM_RG_DIAM[sp] + dadj
                d2v = _EM_RG_AX * (h2 - 4.5f0)^_EM_RG_BX + dadj
                dgj = (d2v - d1v); dgj < 0f0 && (dgj = 0f0)
                dnow[i] += dgj
            end
        end
    end
    # final HTGR1 + ZZRAN + XWT blend + DG override (regent.f:473-560, KT form)
    _sp_order = sortperm(view(t.species, 1:n); alg = Base.Sort.MergeSort)
    @inbounds for oi in 1:n
        i = _sp_order[oi]
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= _em_rg_cap(sp) || t.tpa[i] <= 0.0f0) && continue
        _em_rg_nivar(sp) || continue
        h = t.height[i]; xmn = EM_RG_XMIN[sp]; xmx = _em_rg_cap(sp)
        htgr1 = wk3[i] - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 2.0f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = htgr1 + zzran * _EM_RG_HSIGMA; htgr < 0.15f0 && (htgr = 0.15f0)
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # DG override from the subcycle-accumulated diameter
        dgnew = dnow[i] - d; dgnew < 0f0 && (dgnew = 0f0)
        dgw = dgnew * (1.0f0 - xwt) + xwt * t.diam_growth[i]
        t.diam_growth[i] = dgw
    end
    return s
end
