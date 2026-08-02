# =============================================================================
# regent.jl (teton) — TT small-tree growth (tt/regent.f + smhtgf.f + smdgf.f). Chunk 6.
#
# small_tree_growth!(::Teton): subcycle (REGYR=5-yr) height + DBH for small trees, blended with the
# large-tree DG/HTG over DBH∈[XMIN,XMAX]. Structure mirrors EM's (NPER/KPER + BANEXT/RDNEXT density
# projection), but the per-subcycle increments come from TT's SMHTGF (height) + SMDGF (DBH):
#   SMHTGF DEFAULT: BETA1=exp(B0ACCF+B1ACCF·lnTPCCF); BETA2=exp(B0BCCF+B1BCCF·lnTPCCF);
#                   HTG1=BETA1+BETA2·CR; HTGRL=HTG1+ZRAND·HTG1·(B0ASTD+B1BSTD·CR).  (ZRAND ±2 per tree)
#   SMHTGF aspen(6): SITAGE=(H/26.9825)^(1/1.1752); HTGR=(26.9825·(SITAGE+5)^1.1752−26.9825·SITAGE^1.1752)/(2.54·12);
#                    HTGRL=(HTGR+ZRAND·0.1)·0.75·RSIMOD (RSIMOD from SITEAR(6)).
#   H2 = H1 + HTGRL·(kpj/REGYR).  SMDGF: SDIAM=SDHTCR+SDHPCF·H+SDCR·CR+SDHL4·RD (DF/BS/AS/LP/ES/AF) or a
#   HLESS4 form (WB/LM/OS); D2=max(SDIAM, DIAM).  Blend XWT=(d−XMIN)/(XMAX−XMIN).
# =============================================================================

# smhtgf coefficients (tt/blkdat.f)
const TT_B0ACCF = Float32[1.17527, 1.17527, -4.35709, 0.0, -0.55052, 0.0, -0.90086, -0.55052, -4.35709, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.17527, 0.0]
const TT_B1ACCF = Float32[-0.42124, -0.42124, 0.67307, 0.0, -0.02858, 0.0, 0.16996, -0.02858, 0.67307, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.42124, 0.0]
const TT_B0BCCF = Float32[-2.56002, -2.56002, -2.49682, 0.0, -2.26007, 0.0, -1.50963, -2.26007, -2.49682, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.56002, 0.0]
const TT_B1BCCF = Float32[-0.58642, -0.58642, -0.51938, 0.0, -0.67115, 0.0, -0.61825, -0.67115, -0.51938, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.58642, 0.0]
const TT_B0ASTD = Float32[1.0872, 1.0872, 1.13785, 0.0, 1.0973, 0.0, 1.00749, 1.0973, 1.13785, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0872, 0.0]
const TT_B1BSTD = Float32[-0.0023, -0.0023, -0.00185, 0.0, -0.0013, 0.0, -0.00435, -0.0013, -0.00185, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0023, 0.0]
# smdgf coefficients (tt/smdgf.f)
const TT_SDHTCR = Float32[0.000231, 0.000231, -0.28654, 0.0, 0.04125, -0.41227, -0.41227, 0.04125, -0.15906, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000231, 0.0]
const TT_SDHPCF = Float32[-5e-05, -5e-05, 0.13469, 0.0, 0.17486, 0.16944, 0.16944, 0.17486, 0.15323, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5e-05, 0.0]
const TT_SDCR   = Float32[0.001711, 0.001711, 0.002736, 0.0, -0.002371, 0.003191, 0.003191, -0.002371, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001711, 0.0]
const TT_SDHL4  = Float32[0.17023, 0.17023, 0.00036, 0.0, -0.0007, -0.0022, -0.0022, -0.0007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17023, 0.0]
# blend / floor (tt/regent.f)
const TT_RG_DIAM = Float32[0.4, 0.3, 0.3, 0.4, 0.3, 0.2, 0.4, 0.3, 0.3, 0.5, 0.3, 0.3, 0.2, 0.2, 0.3, 0.2, 0.2, 0.3]
const TT_RG_XMIN = Float32[1.5, 1.5, 1.5, 90.0, 1.5, 1.5, 1.5, 1.5, 1.5, 2.0, 90.0, 90.0, 90.0, 2.0, 0.5, 90.0, 1.5, 0.5]
const TT_RG_XMAX = Float32[3.0, 3.0, 3.0, 99.0, 3.0, 3.0, 3.0, 3.0, 3.0, 5.0, 99.0, 99.0, 99.0, 4.0, 2.0, 99.0, 3.0, 2.0]
const _TT_REGYR = 5.0f0
const _TT_BACON = 0.005454154f0

@inline _tt_smdg_alt(sp::Int) = sp == 3 || (5 <= sp <= 9)   # DF/BS/AS/LP/ES/AF use the alternate SMDGF form
@inline _tt_rg_default(sp::Int) = sp <= 3 || sp == 5 || sp == 6 || (7 <= sp <= 9) || sp == 17  # regent-handled (ttt01)

# tt/smdgf.f — small-tree DBH from height/CR/relative-density.
@inline function _tt_smdgf(sp::Int, h::Float32, cr::Float32, rd::Float32)::Float32
    if _tt_smdg_alt(sp)
        return TT_SDHTCR[sp] + TT_SDHPCF[sp] * h + TT_SDCR[sp] * cr + TT_SDHL4[sp] * rd
    else
        hl4 = h - 4.5f0
        return TT_SDHTCR[sp] * hl4 * cr + TT_SDHPCF[sp] * hl4 * rd + TT_SDCR[sp] * cr + TT_SDHL4[sp] * hl4 + 0.3f0
    end
end

# tt/smhtgf.f — small-tree height increment HTGRL (ZRAND passed in; drawn once per tree).
@inline function _tt_smhtgf(sp::Int, h::Float32, cr::Float32, tpccf::Float32, zrand::Float32, si6::Float32)::Float32
    if sp == 6                                   # quaking aspen (FINDAG closed-form inversion)
        sitage = (h / 26.9825f0)^(1f0 / 1.1752f0)
        hite1 = 26.9825f0 * sitage^1.1752f0
        hite2 = 26.9825f0 * (sitage + 5f0)^1.1752f0
        htgr = (hite2 - hite1) / (2.54f0 * 12f0)
        htgrl = (htgr + zrand * 0.1f0) * 0.75f0
        relsi = (si6 - 30f0) / 70f0; relsi > 1f0 && (relsi = 1f0); relsi < 0f0 && (relsi = 0f0)
        return htgrl * (0.5f0 * (1f0 + relsi))   # regent.f:521 RSIMOD for sp6
    else
        beta1 = exp(TT_B0ACCF[sp] + TT_B1ACCF[sp] * log(tpccf))
        beta2 = exp(TT_B0BCCF[sp] + TT_B1BCCF[sp] * log(tpccf))
        htg1 = beta1 + beta2 * cr
        stddev = htg1 * (TT_B0ASTD[sp] + TT_B1BSTD[sp] * cr)
        return htg1 + zrand * stddev
    end
end

function small_tree_growth!(s::StandState, stash, ::Teton; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    n = t.n; n == 0 && return s
    ba = p.basal_area; relden = p.relative_density
    dgsd = s.control.dg_sd; regyr = _TT_REGYR
    si6 = p.sp_site_index[6]
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    # density projection to each subcycle from the large trees (d≥3) — EM/KT form
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = _TT_BACON * d1 * d1; b2 = _TT_BACON * d2 * d2
            cc1 = tt_tree_ccf(sp, d1); cc2 = tt_tree_ccf(sp, d2)
            bi = (b2 - b1) / 10.0f0; ci = (cc2 - cc1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * 0.985f0^k
                rdnext[j] += k * ci / pr * pn; banext[j] += k * bi * pn
            end
        end
    end
    # per-tree ZRAND (BACHLO ±2), drawn once (small-tree ZZRAN)
    zrand = fill(0f0, n)
    @inbounds for i in 1:n
        (t.dbh[i] >= TT_RG_XMAX[Int(t.species[i])] || t.tpa[i] <= 0f0) && continue
        _tt_rg_default(Int(t.species[i])) || continue
        if dgsd >= 1.0f0
            z = 0f0
            while true; z = bachlo(s.rng, 0.0f0, 1.0f0); (-2f0 <= z <= 2f0) && break; end
            zrand[i] = z
        end
    end
    wk3 = Float32[t.height[i] for i in 1:n]         # subcycle height
    wk5 = Float32[t.dbh[i] for i in 1:n]            # subcycle DBH
    @inbounds for j in 1:nper
        rdj = rdnext[j]; kpj = Float32(kper[j])
        for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]
            (d >= TT_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
            _tt_rg_default(sp) || continue
            h1 = wk3[i]; cr = Float32(t.crown_pct[i])
            pt = Int(t.plot_id[i]); tpccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100f0
            tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
            htgrl = _tt_smhtgf(sp, h1, cr, tpccf, zrand[i], si6)
            h2 = h1 + htgrl * (kpj / regyr)
            wk3[i] = h2
            d2 = _tt_smdgf(sp, h2, cr, rdj)
            d2 < TT_RG_DIAM[sp] && (d2 = TT_RG_DIAM[sp])
            wk5[i] = d2
        end
    end
    # blend HTGR/DG over [XMIN,XMAX] with the large-tree prediction (regent.f:473-560)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= TT_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        _tt_rg_default(sp) || continue
        h = t.height[i]; xmn = TT_RG_XMIN[sp]; xmx = TT_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        htgr = wk3[i] - h; htgr < 0.0f0 && (htgr = 0.0f0)
        t.ht_growth[i] = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]
        bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
        dgnew = (wk5[i] - d) * bark; dgnew < 0f0 && (dgnew = 0f0)   # DG stored inside-bark (·bark)
        t.diam_growth[i] = dgnew * (1.0f0 - xwt) + xwt * t.diam_growth[i]
    end
    return s
end
