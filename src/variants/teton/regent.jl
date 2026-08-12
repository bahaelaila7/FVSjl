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
# regent.f XMIN/XMAX DATA (small/large-tree blend range XWT=(D−XMIN)/(XMAX−XMIN)). Default species were jl 1.5/3.0
# but buildDir differs (DF 2.0/4.0 etc.) ⇒ jl blended the large-tree dgf DG PREMATURELY (at DBH 2, xwt=0.33 vs
# live's 0 = pure small-tree), over-growing DF 2×. UTVAR species (4,11,12,13,16) KEEP 90/99 (jl regent-for-all
# representation; their buildDir 2.0/4.0 would wrongly EXCLUDE them from the UTVAR pass at D>4).
const TT_RG_XMIN = Float32[2.0, 1.0, 2.0, 90.0, 1.0, 1.5, 2.0, 2.0, 2.0, 2.0, 90.0, 90.0, 90.0, 2.0, 0.5, 90.0, 1.0, 0.5]
const TT_RG_XMAX = Float32[3.0, 2.0, 4.0, 99.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0, 99.0, 99.0, 99.0, 4.0, 2.0, 99.0, 5.0, 2.0]  # default species = buildDir DATA (DF sp3=4.0); UTVAR (4,11,12,16) keep 99
# jl per-cycle regent DG cap. ★#158 RESOLVED 2026-08-11 (see TT audit): the conifer 0.2 caps were NOT a band-aid
# for an un-portable model — they were the RAW per-YEAR DGMAX (tt/regent.f:175-177) applied WITHOUT the FINT
# multiplier that live's regent.f:684 applies (IF(TTVAR)DGMX=FINT*DGMAX). MEASURED via FVStt_g16: jl's uncapped
# SMDGF-based dgr already matches live's DG bit-close (i34 0.637 vs 0.646) — so the earlier "smdgf_ called 0×
# ⇒ different single-step model" read was the compiler INLINING smdgf (no CALL), not a different model. The 0.2
# cap over-clamped correct 0.6-1.0 sub-1" DG on ultra-dense cohorts → flat DG → low Reineke DR10 → self-thin
# never fired → the 33% dense-stand over-growth (#158). Fix: cap at FINT·TT_RG_DGMAX_RAW (below). Residual on the
# dense stand (~7% BA, self-thin ~1 cycle late) = the small systematic SMDGF-vs-live-inline DG realization
# difference (~1-2%/tree), the known regent model-version straddle. This TT_RG_DGMAX array is now UNUSED for the
# default+esgent paths (they use TT_RG_DGMAX_RAW×FINT); it is retained only for the UTVAR woodland pass, whose
# ≥2.0 caps are inert (woodland DG ≪ 2.0) and #157-bit-exact.
const TT_RG_DGMAX = Float32[0.2, 0.2, 0.2, 2.0, 0.2, 2.0, 0.2, 0.2, 0.2, 99.0, 2.0, 2.0, 2.0, 2.5, 2.5, 2.0, 0.2, 2.5]
# RAW per-year DGMAX (tt/regent.f:175-177 DATA DGMAX, verbatim). The DEFAULT small-tree path caps DG at
# DGMX=FINT·DGMAX(ISPC) (regent.f:684 IF(TTVAR)DGMX=FINT*DGMAX) — NOT the raw value. #158 fix uses this × fint.
const TT_RG_DGMAX_RAW = Float32[0.2, 0.2, 0.2, 2.0, 0.2, 0.2, 0.2, 0.2, 0.2, 99.0, 2.0, 2.0, 2.0, 2.5, 2.5, 2.0, 0.2, 2.5]
const _TT_REGYR = 5.0f0
const _TT_BACON = 0.005454154f0

@inline _tt_smdg_alt(sp::Int) = sp == 3 || (5 <= sp <= 9)   # DF/BS/AS/LP/ES/AF use the alternate SMDGF form
@inline _tt_rg_default(sp::Int) = sp <= 3 || sp == 5 || sp == 6 || (7 <= sp <= 9) || sp == 14 || sp == 17  # regent-handled (ttt01 + MM=aspen)
@inline _tt_rg_esp(sp::Int) = sp == 14 ? 6 : sp   # MM(14) uses AS(6) aspen coefficients (buildDir dgf/htgf: "MM from UT AS")
@inline _tt_rg_utvar(sp::Int)   = sp == 4 || sp == 11 || sp == 12 || sp == 13 || sp == 16   # PM/UJ/RM/BI/MC: UTVAR regent

# tt/regent.f UTVAR height+diameter (no subcycle, SCALE=1). POTHTG=((SJ/5)·(SJ·1.5−H)/(SJ·1.5))·0.83;
# VIGOR=(150·X³·exp(−6X))+0.3 (X=CR/100), cut ⅔ for PM/UJ/RM; HTGRL=POTHTG·PCTRED·VIGOR·CON (CON=exp(HCOR)).
# Diameter via H-D: DK=(H2−4.5)·10/(SJ−4.5), DKK from H1; DG=(DK−DKK)·bark. Returns (htgr, dg).
@inline function _tt_utvar_regent(sp::Int, h::Float32, d::Float32, cr::Float32, sitear::Float32,
                                  pctred::Float32, con::Float32, bark::Float32,
                                  dgmax::Float32, diam::Float32, scale2::Float32)::Tuple{Float32,Float32}
    sj = sitear
    pothtg = ((sj / 5f0) * (sj * 1.5f0 - h) / (sj * 1.5f0)) * 0.83f0
    x = cr / 100f0
    vigor = (150f0 * x * x * x * exp(-6f0 * x)) + 0.3f0
    vigor > 1f0 && (vigor = 1f0)
    # ⅔ VIGOR cut is ISPC==6 only (regent.f:284); PM/UJ/RM empirically need it (validated), MC/BI do not.
    (sp == 4 || sp == 11 || sp == 12) && (vigor = 1f0 - ((1f0 - vigor) / 3f0))
    htgrl = pothtg * pctred * vigor * con
    # regent.f:780 (Dixon 3/4/09 "PREVENT NEGATIVE HEIGHT GROWTH"): UTVAR floors HTGR to 0.1 ft, NOT 0.
    # This is load-bearing for tall woodland trees whose pothtg goes negative (H > SJ·1.5): the 0.1-ft
    # floor drives DG=(DK−DKK)·bark≈0.1" via the H-D below. Clamping to 0 (old jl) froze UJ/PM/RM DBH.
    htgr = htgrl; htgr < 0.1f0 && (htgr = 0.1f0)
    h2 = h + htgr                                        # HK uses the FLOORED increment (measured vs FVStt_clean)
    # H-D diameter: PM/UJ/RM (4,11,12) use (H−4.5)·10/(SJ−4.5); BI/MC (13,16) use the WC "rule of thumb"
    # DG=0.1·HTG (regent.f CASE(13,14,16,18), the LHTDRG&IABFLG==0 branch that fires for TT's MC/BI —
    # HTDBH is a stub and the AA-fit leaves IABFLG=0 at growth, so DK/DKK are unused; measured fort.89).
    if sp == 13 || sp == 16
        dg = 0.1f0 * htgr                                # ·XRDGRO=1 (no MULTS)
    else
        dk  = (h2 - 4.5f0) * 10f0 / (sj - 4.5f0); dk  < 0.1f0 && (dk  = 0.1f0)
        dkk = h < 4.5f0 ? d : (h - 4.5f0) * 10f0 / (sj - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
        dg = (dk - dkk) * bark
    end
    dg < 0f0 && (dg = 0f0)
    dg > dgmax && (dg = dgmax)                           # DGMX cap (regent.f:1047)
    # DDS-sqrt conversion for period/bark consistency (regent.f:1049-1054); SCALE2=YR/NTYR
    dds = dg * (2f0 * bark * d + dg) * scale2
    arg = (d * bark)^2 + dds
    dg = arg > 0f0 ? sqrt(arg) - bark * d : 0f0
    (d + dg) < diam && (dg = diam - d)                   # DIAM floor (regent.f:1056)
    return (htgr, dg)
end

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
    if sp == 6 || sp == 14                       # aspen (6) / mountain maple (14, aspen coefs) — FINDAG closed-form
        # NB (#158, 2026-08-07): sitage here is the INVERSE-height age s.t. hite1=26.9825·sitage^1.1752 = h
        # (self-consistent: current height in FEET). Do NOT "fix" it to (h·2.54·12/26.9825) to match findag.f:96
        # literally — MEASURED: that regresses (jl over-grows; nofix is bit-exact with live at 2003/2013 on
        # 3189335010690). jl's feet-native form is PRIOR-VALIDATED bit-exact vs live; leave as-is.
        sitage = (h / 26.9825f0)^(1f0 / 1.1752f0)
        hite1 = 26.9825f0 * sitage^1.1752f0
        hite2 = 26.9825f0 * (sitage + 5f0)^1.1752f0
        htgr = (hite2 - hite1) / (2.54f0 * 12f0)
        # ·0.75 (Dixon 8-27-92) is smhtgf.f CASE(6)-specific; MM(14) is CASE DEFAULT ⇒ no ·0.75 (matches live faster MM)
        htgrl = (htgr + zrand * 0.1f0) * (sp == 6 ? 0.75f0 : 1.0f0)
        # ★#189 (2026-08-12): tt/regent.f:521-527 applies an ASPEN(sp6)-ONLY RSIMOD after SMHTGF:
        # RELSI=clamp((SITEAR(6)−30)/70,0,1); RSIMOD=0.5·(1+RELSI); HTGRL·=RSIMOD. jl OMITTED it (the old
        # "NO RSIMOD, that's CASE(15)=NC" comment MIS-READ the buildDir — regent.f:521 gates on ISPC.EQ.6).
        # INERT on high-site aspen (SITEAR≥100 ⇒ RSIMOD=1, e.g. 3189335010690 where jl was "validated") but on
        # low-site (SITEAR=42 → RSIMOD=0.586) jl over-grew small-aspen height ~1.7×. MEASURED via FVStt_g16:
        # 753175613290487 live grows small aspen ~3.5 ft vs jl ~9. sp14(MM) EXEMPT (live gates ISPC.EQ.6 only).
        if sp == 6
            relsi = clamp((si6 - 30f0) / 70f0, 0f0, 1f0)
            htgrl *= 0.5f0 * (1f0 + relsi)
        end
        return htgrl
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
    # KNOWN faithful gap: buildDir HTGR=POTHTG·PCTRED·VIGOR·CON (regent.f:350, RHCON=1 @ line 880). Tested a
    # conifer-only version (aspen sp6 exempt = Sheppard-already-final; kept aspen exact) but it BARELY moved DF
    # (2040 108→103, still 45% over live 71) ⇒ NOT the DF later-cycle residual (= large-tree dgf DF DG at DBH 4-6,
    # un-validatable: live fort.79 caps at DBH 2.0). Reverted — an untestable-for-TT-conifers change that doesn't
    # fix the visible residual. The faithful PCTRED·VIGOR·CON gap remains (conifer-only if ever added).
    @inbounds for j in 1:nper
        rdj = rdnext[j]; kpj = Float32(kper[j])
        for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]
            (d >= TT_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
            _tt_rg_default(sp) || continue
            h1 = wk3[i]; cr = Float32(t.crown_pct[i])
            pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100f0
            tpccf = pccf; tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)   # smhtgf clamps [25,300]
            esp = _tt_rg_esp(sp)                       # MM(14)→AS(6) coefficient mapping (DBH); smhtgf handles 14 directly
            htgrl = _tt_smhtgf(sp, h1, cr, tpccf, zrand[i], si6)
            h2 = h1 + htgrl * (kpj / regyr)
            wk3[i] = h2
            d2 = _tt_smdgf(esp, h2, cr, pccf)          # SMDGF gets the RAW point CCF (regent.f:574), not stand relden
            d2 < TT_RG_DIAM[sp] && (d2 = TT_RG_DIAM[sp])
            wk5[i] = d2
        end
    end
    # blend HTGR/DG over [XMIN,XMAX] with the large-tree prediction (regent.f:735-943)
    scale2 = htg_period(s.variant) / fint          # SCALE2 = YR/NTYR (period scaling of the DBH increment)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= TT_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        _tt_rg_default(sp) || continue
        h = t.height[i]; xmn = TT_RG_XMIN[sp]; xmx = TT_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        # HTG blend + size cap
        htgr = wk3[i] - h; htgr < 0.0f0 && (htgr = 0.0f0)
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # DG (regent.f:923-942): DK=smdgf(grown H)=wk5, DKK=smdgf(ORIGINAL H), DG=(DK−DKK)·bark → DDS → DG. HK≥4.5.
        hk = h + htg
        dfl = d < TT_RG_DIAM[sp] ? TT_RG_DIAM[sp] : d       # regent.f:703 D floored to DIAM(sp)
        dgk = 0f0
        if hk >= 4.5f0
            pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100f0
            dkk = _tt_smdgf(_tt_rg_esp(sp), h, Float32(t.crown_pct[i]), pccf)      # DBH from the ORIGINAL height (MM→AS)
            bark = bark_ratio(c.bark_a, c.bark_b, sp, dfl)
            dgr = (wk5[i] - dkk) * bark
            dds = dgr * (2f0 * bark * dfl + dgr) * scale2
            arg = (dfl * bark)^2 + dds
            dgk = arg > 0f0 ? sqrt(arg) - bark * dfl : 0f0
            # DGMX cap (regent.f:720). ★#158 FIX: live's TT DGMX = FINT·DGMAX(ISPC) (regent.f:684, IF(TTVAR)
            # DGMX=FINT*DGMAX) — the raw DGMAX (0.2 for conifers) is a per-YEAR cap that MUST be scaled by FINT.
            # jl previously capped at the raw DGMAX ⇒ clamped correct 0.6-1.0 sub-1" DG down to 0.2 → the #158
            # 33% dense-stand over-growth (flat DG → low Reineke DR10 → self-thin never fires). MEASURED: jl's
            # uncapped dgr already matches live's DG bit-close (i34 0.637 vs 0.646); the cap was the sole bug.
            dgmx = fint * TT_RG_DGMAX_RAW[sp]
            dgk > dgmx && (dgk = dgmx)
        end
        t.diam_growth[i] = dgk * (1.0f0 - xwt) + xwt * t.diam_growth[i]
        # #148 latent bug (2): DIAM floor on the DEFAULT path (regent.f:576 D2=max(smdgf,DIAM) + :1056), missing here.
        # Without it a tiny tree whose grown smdgf-DBH floors to DIAM while its original DKK exceeds DIAM gets a
        # NEGATIVE dgk → NEGATIVE DBH → NaN in crown. The H-D branch already floors (line 83). Floor (d+DG)≥DIAM.
        (d + t.diam_growth[i]) < TT_RG_DIAM[sp] && (t.diam_growth[i] = TT_RG_DIAM[sp] - d)
        # Update the TRIPLING stash so the upper/lower sub-records get the REGENT DG/HTG, not the stale
        # large-tree dgf DG (triple_records! sets diam_growth[u]=dgU, [l]=dgL). Without this, 40% of a tripled
        # small tree grows via the large-tree DG — the DF-stand 2× over-growth (SN/UTVAR both do this; default
        # pass was missing it). No ZZRAN spread (= UTVAR simplification); central value on both sub-records.
        if stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
            stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
            stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
            !isempty(stash.is_small) && (stash.is_small[i] = true)
        end
    end
    # UTVAR pass (PM/UJ/RM) — separate: no subcycle, xwt=0 (XMIN=90 ⇒ pure small-tree). regent.f ELSEIF(UTVAR).
    if any(j -> _tt_rg_utvar(Int(t.species[j])), 1:n)
        avh = p.avg_height
        xpr = avh * (relden / 100f0); xpr > 300f0 && (xpr = 300f0)   # PCTRED density arg (regent.f:337)
        pctred = 1.11436f0 + xpr * (-0.011493f0 + xpr * (0.43012f-4 + xpr * (-0.72221f-7 +
                 xpr * (0.5607f-10 - xpr * 0.1641f-13))))
        pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
        @inbounds for i in 1:n
            sp = Int(t.species[i]); _tt_rg_utvar(sp) || continue
            d = t.dbh[i]; (d >= TT_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
            h = t.height[i]; cr = Float32(t.crown_pct[i])
            sitear = p.sp_site_index[sp]
            con = exp(c.htg_cor_small[sp])                            # RHCON·exp(HCOR), RHCON=1
            bark = tt_bratio(sp, d)
            scale2 = htg_period(s.variant) / fint                     # YR/NTYR
            htgr, dg = _tt_utvar_regent(sp, h, d, cr, sitear, pctred, con, bark,
                                        TT_RG_DGMAX[sp], TT_RG_DIAM[sp], scale2)
            cap = s.control.sp_size_cap[sp, 4]
            (h + htgr > cap) && (htgr = max(cap - h, 0.1f0))
            t.ht_growth[i] = htgr
            t.diam_growth[i] = dg
            # Tripling: PM/UJ/RM regent DG is deterministic (tiny VARDG) ⇒ ~no spread. Override the stale
            # dgU/dgL (built in diameter_growth! from the DIAGR placeholder, giving a spurious wide spread)
            # with the regent DG so the tripled sub-records match (regent DG has no ZZRAN for UTVAR).
            if stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
                stash.dgU[i] = dg; stash.dgL[i] = dg
            end
        end
    end
    return s
end

# tt_esgent! (tt/esgent.f) — grow the JUST-ESTABLISHED regen records IN their birth cycle via the TT regent,
# for the PARTIAL period (FINT−GENTIM; GENTIM=FINT−5 ⇒ 5 yr for a 10-yr cycle). Mirrors cr_esgent!: western
# variants grow birth-cycle regen (esgent.f→REGENT); eastern leave them ungrown (GRADD order). Fixes the ESTAB
# 1-cycle size/TopHt lag on planted TT stands (jl planted trees previously appeared ungrown at the first report).
# FIRST-CUT: default regent species (_tt_rg_default: LM/DF/WB/BS/AS/LP/ES/AF/OS/MM). UTVAR (PM/UJ/RM/BI/MC) +
# non-regent (PP) birth-cycle growth are a scoped follow-up (they need their own per-tree partial-cycle path).
function tt_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0)
    t = s.trees; c = s.calib; p = s.plot; dens = s.density
    nstart >= t.n && return s
    gentim = max(fint - 5.0f0, 0.0f0)
    subcyc = (fint - gentim) / _TT_REGYR        # birth-cycle subcycles (=1 for fint=10)
    scale2 = htg_period(s.variant) / fint       # DDS period scaling (YR/NTYR), = regular cycle
    si6 = p.sp_site_index[6]; dgsd = s.control.dg_sd
    @inbounds for i in (nstart+1):t.n
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= TT_RG_XMAX[sp] || !_tt_rg_default(sp)) && continue
        h = t.height[i]; cr = Float32(t.crown_pct[i])
        pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100.0f0
        tpccf = pccf; tpccf > 300.0f0 && (tpccf = 300.0f0); tpccf < 25.0f0 && (tpccf = 25.0f0)
        zrand = 0.0f0
        if dgsd >= 1.0f0
            while true; zrand = bachlo(s.rng, 0.0f0, 1.0f0); (-2.0f0 <= zrand <= 2.0f0) && break; end
        end
        esp = _tt_rg_esp(sp)                     # MM(14)→AS(6) coefficient mapping
        htgrl = _tt_smhtgf(esp, h, cr, tpccf, zrand, si6)
        htg = htgrl * subcyc; htg < 0.0f0 && (htg = 0.0f0)
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        h2 = h + htg
        dgk = 0.0f0
        if h2 >= 4.5f0                           # only trees that reach breast height get a real DBH increment
            dfl = d < TT_RG_DIAM[sp] ? TT_RG_DIAM[sp] : d
            bark = bark_ratio(c.bark_a, c.bark_b, sp, dfl)
            d2 = _tt_smdgf(esp, h2, cr, pccf); d2 < TT_RG_DIAM[sp] && (d2 = TT_RG_DIAM[sp])
            dkk = _tt_smdgf(esp, h, cr, pccf)
            dgr = (d2 - dkk) * bark
            dds = dgr * (2.0f0 * bark * dfl + dgr) * scale2
            arg = (dfl * bark)^2 + dds
            dgk = arg > 0.0f0 ? sqrt(arg) - bark * dfl : 0.0f0
            dgmx = fint * TT_RG_DGMAX_RAW[sp]                       # ★#158: live DGMX=FINT·DGMAX (regent.f:684), not raw
            dgk > dgmx && (dgk = dgmx)
        end
        t.height[i] = h2
        dgk > 0.0f0 && (t.dbh[i] = d + dgk / tt_bratio(sp, d))   # outside-bark DBH (simulate.jl:499)
    end
    return s
end
