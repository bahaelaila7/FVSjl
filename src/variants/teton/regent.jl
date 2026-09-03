# =============================================================================
# regent.jl (teton) — TT small-tree growth (tt/regent.f + smhtgf.f + smdgf.f). Chunk 6.
#
# small_tree_growth!(::Teton): subcycle (REGYR=5-yr) height + DBH for small trees, blended with the
# large-tree DG/HTG over DBH∈[XMIN,XMAX]. Structure mirrors EM's (NPER/KPER + BANEXT/RDNEXT density
# projection), but the per-subcycle increments come from TT's SMHTGF (height) + SMDGF (DBH):
#   SMHTGF DEFAULT: BETA1=exp(B0ACCF+B1ACCF·lnTPCCF); BETA2=exp(B0BCCF+B1BCCF·lnTPCCF);
#                   HTG1=BETA1+BETA2·CR; HTGRL=HTG1+ZRAND·HTG1·(B0ASTD+B1BSTD·CR).  (ZRAND ±2 per tree)
#   SMHTGF aspen(6): SITAGE=(H·2.54·12/26.9825)^(1/1.1752) [findag.f:96, metric]; HTGR=(26.9825·(SITAGE+5)^1.1752
#                    −26.9825·SITAGE^1.1752)/(2.54·12); HTGRL=(HTGR+ZRAND·0.1)·0.75·RSIMOD (RSIMOD from SITEAR(6)).
#   H2 = H1 + HTGRL·(kpj/REGYR).  ASPEN applies this on the FIRST subcycle ONLY (regent.f:415 label-16 gate).  SMDGF: SDIAM=SDHTCR+SDHPCF·H+SDCR·CR+SDHL4·RD (DF/BS/AS/LP/ES/AF) or a
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
# regent.f XMIN/XMAX DATA (small/large-tree blend range XWT=(D−XMIN)/(XMAX−XMIN)).
# VERBATIM tt/regent.f:161-171 DATA (corrected 2026-09-03). Prior jl used 2.0/4.0 for the TTVAR conifers — a
# mis-read buildDir. The authoritative DATA (and the relinked FVStt single-.o trace: XMN=1.5, XMX=3.0 for
# LM/LP/ES/AF) is 1.5/3.0. XMIN/XMAX govern the HEIGHT-increment blend window (XWT) only. Correcting the window
# to [1.5,3.0] alone over-grew ~half the M331D woodland sample because the OLD code also (wrongly) XWT-blended
# the DIAMETER increment — see the DG assignment below, which is fixed in the same commit to pure regent DG.
# UTVAR (4,11,12,13,16) 90/99; NC (15) 0.5/2.0; OS (18) 0.5/2.0.
const TT_RG_XMIN = Float32[1.5, 1.5, 1.5, 90.0, 1.5, 1.5, 1.5, 1.5, 1.5, 2.0, 90.0, 90.0, 90.0, 2.0, 0.5, 90.0, 1.5, 0.5]
const TT_RG_XMAX = Float32[3.0, 3.0, 3.0, 99.0, 3.0, 3.0, 3.0, 3.0, 3.0, 5.0, 99.0, 99.0, 99.0, 4.0, 2.0, 99.0, 3.0, 2.0]
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
@inline _tt_rg_utvar(sp::Int)   = sp == 4 || sp == 11 || sp == 12 || sp == 13 || sp == 15 || sp == 16 || sp == 18   # PM/UJ/RM/BI/NC/MC/OH: UTVAR regent (regent.f:386 CASE(4,11:16,18); MM(14) routed via default+aspen)
# regent.f DBH BREAK (blkdat BREAK) — D≥BREAK skips the small-tree DBH increment (large-tree DG kept), D<BREAK gets
# the regent H-D DBH. NC/OH (15,18) have BREAK=1 (they leave the small-tree DBH regime at 1"); all other UTVAR = 99.
const TT_RG_BREAK = Float32[3.,3.,3.,99.,3.,3.,3.,3.,3.,3.,99.,99.,99.,3.,1.,99.,3.,1.]
# tt/blkdat.f HT1(AX)/HT2(BX) — the NC/OH (15,18) Wykoff H-D coefficients used by the UTVAR DBH dub for D<4.5→>4.5
# crossers (regent.f:907 `DK=(BX/(ALOG(HK-4.5)-AX))-1`, IABFLG=1 ⇒ AX=HT1). Same for 15 and 18; other sp unused.
const TT_HT1 = Float32[4.1920,4.1920,4.5175,3.2000,4.5822,4.4625,4.4625,4.5822,4.3603,4.993,3.2000,3.2000,4.7000,4.4421,4.4421,5.1520,4.1920,4.4421]
const TT_HT2 = Float32[-5.1651,-5.1651,-6.5129,-5.0000,-6.4818,-5.2223,-5.2223,-6.4818,-5.2148,-12.430,-5.0000,-5.0000,-6.3260,-6.5405,-6.5405,-13.5760,-5.1651,-6.5405]

# tt/regent.f BI/MC (13,16) inventory Curtis-Arney HT-DBH (regent.f:874-901, the "SO originally from WC"
# equations). Since LHTDRG(13)=LHTDRG(16)=.FALSE. (tt/grinit.f), the `IF(.NOT.LHTDRG .OR. …)` gate at
# regent.f:872 is ALWAYS true ⇒ the inventory equation is used and the 0.1·HTG rule-of-thumb (regent.f:1009,
# gated LHTDRG.AND.IABFLG.EQ.0) NEVER fires. Returns the DBH predicted from a height HGT (≥4.5).
# P2/P3/P4 (regent.f:875-882): MC(16)=1709.7229/5.8887/−0.2286, BI(13)=76.5170/2.2107/−0.6365.
@inline function _tt_bimc_dk(hgt::Float32, hat3::Float32, p2::Float32, p3::Float32, p4::Float32)::Float32
    if hgt >= hat3
        # regent.f:885 DK=EXP(ALOG((ALOG(HK-4.5)-ALOG(P2))/(-1.*P3))*1./P4)  (·1./P4 ≡ /P4, L-to-R)
        return fexp(flog((flog(hgt - 4.5f0) - flog(p2)) / (-p3)) / p4)
    else
        # regent.f:888 DK=(((HK-4.51)*2.7)/(HAT3-4.51))+0.3 ; denom = 4.5+P2·exp(−P3·3^P4) − 4.51 = HAT3−4.51
        return (((hgt - 4.51f0) * 2.7f0) / (hat3 - 4.51f0)) + 0.3f0
    end
end

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
    # regent.f:817 `IF(.NOT.TTVAR .AND. HK.LT.4.5)`: a sub-breast-height UTVAR seedling gets DG=0 and DBH grows
    # only via the tiny +0.001·HK nudge — AND it SKIPS the DIAM floor (that floor lives inside the HK≥4.5 ELSE
    # branch, regent.f:1048). Applies to ALL UTVAR species (PM/UJ/RM/BI/NC/MC/OH). Without this, jl floored a
    # 0.1" woodland seedling's DG to DIAM(sp)−d (≈0.2) instead of 0 — masked on low-TPA junipers, but sp18/NC-OH
    # seedlings dominate woodland TPA (repro 1856054779290487: 3900 TPA of 0.1" OH grew to QMD 1.8 vs oracle 1.0).
    h2 < 4.5f0 && return (htgr, 0f0)
    # H-D diameter: PM/UJ/RM (4,11,12) use (H−4.5)·10/(SJ−4.5); BI/MC (13,16) use the inventory Curtis-Arney
    # HT-DBH equation (regent.f:874-901; NOT the 0.1·HTG rule-of-thumb — that requires LHTDRG.AND.IABFLG.EQ.0,
    # but LHTDRG(13/16)=.FALSE.); NC/OH (15,18) use the Wykoff DK=(BX/(ln(HK−4.5)−AX))−1 (regent.f:907).
    # ★TT M331D woodland over-growth ROOT FIX (2026-09-03): jl formerly used 0.1·HTG for BI/MC, which
    # over-grew Gambel-oak (FIA 322→BI/13) DBH ~1.66× (measured FVStt_dbg dense stand 51031230020004 cyc1:
    # jl dg=0.383 vs oracle DG=(DK−DKK)·bark=0.230, DK from inventory eqn=0.3324, DKK=D=0.1). Gambel oak
    # dominates the M331D woodland TPA, so this one-directional over-grow drives the whole dense-woodland cap.
    if sp == 13 || sp == 16
        p2, p3, p4 = sp == 16 ? (1709.7229f0, 5.8887f0, -0.2286f0) : (76.5170f0, 2.2107f0, -0.6365f0)
        hat3 = 4.5f0 + p2 * fexp(-p3 * fpow(3.0f0, p4))      # regent.f:883
        dk  = _tt_bimc_dk(h2, hat3, p2, p3, p4)
        dkk = h <= 4.5f0 ? d : _tt_bimc_dk(h, hat3, p2, p3, p4)   # regent.f:891 H≤4.5 ⇒ DKK=D
        # regent.f:1001 DK<0 or DKK<0 ⇒ 0.2·HTG rule-of-thumb, else (DK−DKK)·bark·XRDGRO (XRDGRO=1)
        dg = (dk < 0f0 || dkk < 0f0) ? htgr * 0.2f0 * bark : (dk - dkk) * bark
        dg < 0f0 && (dg = 0.1f0)                             # regent.f:1011 BI/MC floor is 0.1 (not 0)
    elseif sp == 15 || sp == 18
        ax = TT_HT1[sp]; bx = TT_HT2[sp]
        dk  = bx / (flog(h2 - 4.5f0) - ax) - 1f0; dk < 0.1f0 && (dk = 0.1f0)
        dkk = h <= 4.5f0 ? d : bx / (flog(h - 4.5f0) - ax) - 1f0
        # regent.f:1018 CASE(4,11,12,14,15,18): DK<0 or DKK<0 ⇒ the 0.2·HTG rule-of-thumb, else (DK−DKK)·bark
        dg = (dk < 0f0 || dkk < 0f0) ? htgr * 0.2f0 * bark : (dk - dkk) * bark
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
        # #205 (2026-08-13): sitage = FINDAG's aspen inverse-height age, findag.f:96 CASE(6,14):
        # SITAGE=(H·2.54·12/26.9825)^(1/1.1752) — H converted feet→cm (the Sheppard curve 26.9825·age^1.1752 is
        # METRIC/cm). jl formerly used (h/26.9825)^(1/1.1752) (feet, no conversion), which kept aspen on the flat
        # young part of the convex curve ⇒ a ~constant per-cycle increment instead of live's increasing one. The
        # old "don't add ·2.54·12, it regresses" note was a two-bug ARTIFACT: the SITAGE under-growth was masked by
        # the fint=10 subcycle DOUBLE-APPLY (now gated to J=1 below). With both fixes jl htgr is bit-close to live —
        # MEASURED FVStt_g16 SMHTGF: #205 seedling jl 2.76=live 4.77·RSIMOD; ttt01 seedling jl 5.06=live 5.09.
        sitage = (h * 2.54f0 * 12f0 / 26.9825f0)^(1f0 / 1.1752f0)
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

"""
    tt_regent_hcor_aspen_init!(s, isct, ind1, saved_dbh)

TT REGENT small-tree HEIGHT self-calibration for the SMHTGF-aspen species — aspen (sp6, TTVAR) and
Rocky-Mountain maple MM (sp14, jl models it on the aspen closed form). Ports tt/regent.f:1097-1362 (the
LSTART CORNEW pass, `CALL REGENT(.FALSE.,1)` from cratet.f:713): for each sub-5" aspen/MM record carrying a
measured height increment, the predicted increment EDH = NPER·SMHTGF(HT) (with the aspen RSIMOD folded in for
sp6 by `_tt_smhtgf`) accumulated over the NPER calibration subcycles — the aspen closed form is TPCCF/CR-
independent so every subcycle is identical, hence the ·NPER (verified vs FVStt_g16: stand 325585226489998
tree19 EDH=2·5.9956·RSIMOD=6.338). CORNEW = Σ(HTG·SCALE3·P)/Σ(EDH·P); HCOR_raw = ln(CORNEW), trapped to
[0.0821,12.1825] (the ±2.5σ ERRGRO trap). The RAW HCOR → `htg_cor_init`; `calibrate_diameter_growth!`'s shared
attenuation (`htg_cor_small = dg_cor_goal + cormlt_h·(htg_cor_init − dg_cor_goal)`, dgdriv.f:188-213) then
produces the per-cycle applied CON = exp(HCOR), which `small_tree_growth!` multiplies into the height growth.

Without it aspen/MM held CON=1 ⇒ ~3.6× small-tree HEIGHT over-growth on the M331D woodland cluster (stand
325585226489998: oracle CON cyc1=0.4476, jl 1.0; aspen i19 htg 1.69 vs oracle 0.756=1.69·0.4476). The
calibration's SMHTGF ZRAND draw is taken as 0 here — jl draws ZRAND lazily in the growth loop, one phase later
than FVS's cratet-time draw, so the streams are not aligned and replicating the draw here would perturb the
(separately-validated) TTVAR-conifer RNG without gaining the FVS values; the ±ZRAND·0.1 perturbation averages
across the calibration sample ⇒ CORNEW within ~0.4% of the oracle (measured CON 0.449 vs 0.4476).
"""
function tt_regent_hcor_aspen_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                                    saved_dbh::AbstractVector)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    si6 = p.sp_site_index[6]
    regyr = _TT_REGYR
    # SCALE3 = 10./FINTH (NOT REGYR/FINTH). regent.f:1106-1116 selects SCALE3 on `SELECT CASE(ISPC)`, but ISPC
    # there is STALE — it is MAXSP+1 (=19), left one past the end by the immediately-preceding `DO 45 ISPC=1,MAXSP`
    # HCOR-init loop — so 19 is out of every CASE range and every species falls to `CASE DEFAULT: SCALE3=10./FINTH`.
    # This single stand-level SCALE3 is then applied to ALL species (aspen's nominal CASE(...)=REGYR/FINTH branch
    # never runs). FINTH = the HTG measurement period (FVS_STANDINIT HTG_MEASURE / GROWTH kwd), default 10.
    # VERIFIED vs FVStt_g16 (stand 325585226489998): FINTH=10 ⇒ SCALE3=1.0, TERM=HTG·1 (RGTREE dump).
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : Float32(htg_period(s.variant))
    scale3 = 10f0 / finth
    ntyr = Int(round(htg_period(s.variant))); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    npf = Float32(nper)
    @inbounds for sp in (6, 14)
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        snx = 0f0; sny = 0f0; nh = 0
        for k in i1:i2
            i = ind1[k]
            saved_dbh[i] >= 5f0 && continue                   # DBH<5 excluded (regent.f:1228)
            htg = t.ht_growth[i]; htg < 0.001f0 && continue   # measured height increment > 0.001
            h = t.height[i]                                   # SMHTGF reads HT(I) — the current height
            (h - htg) < 0.01f0 && continue                    # backdated start-of-period H ≥ 0.01 (regent.f:1227-1228)
            cr = Float32(t.crown_pct[i])
            edh = npf * _tt_smhtgf(sp, h, cr, 100f0, 0f0, si6)   # NPER identical subcycles; tpccf/cr unused for aspen form
            pr = t.tpa[i]
            snx += edh * pr; sny += htg * scale3 * pr; nh += 1
        end
        nh < 5 && continue                                    # NCALHT (default 5)
        cornew = snx > 0f0 ? sny / snx : 1f0
        cornew <= 0f0 && (cornew = 1f-4)
        (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)   # ±2.5σ ln(C) trap (regent.f:1358)
        c.htg_cor_init[sp] = log(cornew)                      # raw HCOR; shared attenuation → htg_cor_small
    end
    return s
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
                # ★TT M331D woodland over-growth ROOT FIX (2026-09-03): tt/ccfcal.f returns CCFT=poly(D)·P (CCF ×
                # trees-per-acre). FVS's projection `RDNEXT+=K·CI/P·PN` (regent.f:274) nets ONE factor of P
                # (CI=(C2−C1)/10 already carries ·P; CI/P·PN=CI/P·P·0.985^K). jl's `tt_tree_ccf` returns the
                # per-tree CCF WITHOUT ·P, so dividing `ci/pr` dropped the density weighting entirely — rdnext
                # barely moved (28.7→28.9 vs oracle 28.7→32.7), so PPCCF≈1.007 vs oracle 1.137, so subcycle-2
                # SMHTGF saw too-low CCF ⇒ too-high height increment ⇒ small-tree DBH over-grew, ONE-DIRECTIONAL
                # and compounding across the woodland cluster. banext is already correct (BI=BACON·D² carries no
                # ·P, so `k·bi·pn` nets ·P). Match by dropping the `/pr`: `k·ci·pn` = k·(poly2−poly1)/10·P·0.985^k.
                rdnext[j] += k * ci * pn; banext[j] += k * bi * pn
            end
        end
    end
    # per-tree ZRAND (BACHLO ±2). tt/smhtgf.f draws ZRAND(I) ONCE per tree and PERSISTS it across
    # subcycles, cycles AND tripled sub-records (SMHTGF: `IF(ZRAND(I).NE.-999.) GO TO 20`; the reset
    # to -999 fires only at projection start / new tree / when a cycle's growth floored to 0.1 ft).
    # jl formerly re-drew a FRESH local zrand every cycle (and a distinct one per tripled record),
    # so the persisted-negative-draw stands over-grew and the shared-triple spread was lost — the
    # dominant TT seedling over-growth (CCF worst-col). Persist in t.tree_random (ZRAND), which the
    # tripling copy list already inherits to sub-records. Inventory default 0 ⇒ first-cycle draw.
    # ★#158 ZRAND draw ORDER: FVS draws ZRAND inside SMHTGF, called from REGENT's SPECIES-MAJOR subcycle loop
    # (`DO 17 J; DO 16 ISPC=1,MAXSP; DO 15 I3=I1,I2; I=IND1(I3); CALL SMHTGF`) — the first-draw of each tree
    # happens at J=1 in species order (sp6 AS before sp7 LP, etc.). jl drew in raw tree-index order ⇒ on a
    # multi-species seedling cohort the per-tree ZRAND deviates were mis-assigned (aspen index-3 drew 3rd, not
    # 1st) ⇒ conifer SMHTGF HTGRL=HTG1+ZRAND·STDDEV landed on the wrong deviate ⇒ dense over-growth. Draw
    # species-major (index order within species = FVS IND1). Same root as the crown-dub reorder above; both are
    # needed for the RNG stream to track FVS. TT-only (this method dispatches on ::Teton).
    zorder = sort(collect(1:n); by = ii -> (Int(t.species[ii]), ii))
    @inbounds for i in zorder
        (t.dbh[i] >= TT_RG_XMAX[Int(t.species[i])] || t.tpa[i] <= 0f0) && continue
        _tt_rg_default(Int(t.species[i])) || continue
        if dgsd >= 1.0f0
            zr = t.tree_random[i]
            if zr == 0f0 || zr == -999f0            # 0 = inventory default (first draw); -999 = reset
                z = 0f0
                while true; z = bachlo(s.rng, 0.0f0, 1.0f0); (-2f0 <= z <= 2f0) && break; end
                t.tree_random[i] = z
            end
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
        ky = 0; for m in 1:j; ky += kper[m]; end          # cumulative subcycle length (regent.f KY=KY+KPER(J))
        kymort = 0.985f0 ^ ky
        # ★TT M331D woodland over-growth FIX (2026-09-03): PPCCF — the subcycle proportional point-CCF
        # adjustment (regent.f:352-353 `PPCCF=1.0+(RDJ-RELDEN)/RELDEN`). FVS scales each tree's point CCF by
        # the PROJECTED stand-density increase for subcycle J before feeding it to SMHTGF as TPCCF. jl omitted
        # it — every subcycle used the RAW point CCF — so subcycle J≥2 saw a lower CCF than FVS (density always
        # grows ⇒ PPCCF>1), yielding a HIGHER SMHTGF height increment (BETA1/BETA2 fall with CCF) ⇒ over-grown
        # small-tree height ⇒ over-grown SMDGF DBH, ONE-DIRECTIONAL and compounding across the woodland cluster.
        # MEASURED vs FVStt (stand 335 cyc1 i6): oracle J=2 TPCCF=65.18=57.33·1.1369 HTGRL=1.079 vs jl raw
        # TPCCF=57.33 HTGRL=1.121. SMDGF (DBH, regent.f:574) keeps the RAW PCCF — only the height model uses PPCCF.
        ppccf = relden > 0f0 ? 1f0 + (rdj - relden) / relden : 0f0
        for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]; pr = t.tpa[i]
            (d >= TT_RG_XMAX[sp] || pr <= 0f0) && continue
            _tt_rg_default(sp) || continue
            h1 = wk3[i]; cr = Float32(t.crown_pct[i])
            pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100f0
            tpccf = pccf * ppccf; tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)   # PPCCF-adjusted; smhtgf clamps [25,300]
            esp = _tt_rg_esp(sp)                       # MM(14)→AS(6) coefficient mapping (DBH); smhtgf handles 14 directly
            # #205 (2026-08-13): tt/regent.f:415 label-16 gate — aspen(6)/CIVAR(10)/UTVAR(4,11:16,18) apply the
            # small-tree height+DBH increment ONLY on the FIRST subcycle (`(ISPC.EQ.6 .OR. UTVAR .OR. CIVAR) .AND.
            # J.GT.1 GO TO 16`); only TTVAR conifers (1:3,5:9,17) subcycle across all J. jl formerly subcycled aspen
            # every j ⇒ for a 10-yr cycle (nper=2) it DOUBLE-applied the Sheppard SMHTGF increment ⇒ aspen small-tree
            # height/DBH over-grew ~1.9× (11796095010690 +34% BA). At J=1 aspen SCALE=KPER(1)/REGYR=kpj/regyr matches.
            # (sp14 MM is UTVAR with SCALE=NTYR/YR — separate; jl handles it via aspen coefs, not fixed here.)
            (sp == 6 && j > 1) && continue
            htgrl = _tt_smhtgf(sp, h1, cr, tpccf, t.tree_random[i], si6)
            # smhtgf.f: if the estimated increment ≤ 0.1 ft, floor to 0.1 and reset ZRAND (redraw next
            # cycle). Reset only on the final subcycle so the persisted draw is not corrupted mid-cycle.
            if htgrl <= 0.1f0
                htgrl = 0.1f0
                (dgsd >= 1.0f0 && j == nper) && (t.tree_random[i] = -999f0)
            end
            # CON = RHCON·exp(HCOR) — the REGENT small-tree HEIGHT self-calibration (regent.f:420,552:
            # `H2=H1+HTGRL*SCALE*XRHGRO*CON`). For aspen(6)/MM(14) HCOR = the CORNEW calibration seeded in
            # tt_regent_hcor_aspen_init! and attenuated per-cycle in diameter_growth!; for the TTVAR conifers
            # htg_cor_small stays ~0 (no measured-HTG calibration ⇒ CON≈1). Without this factor aspen/MM held
            # CON=1 ⇒ ~3.6× small-tree HEIGHT over-growth on the M331D woodland cluster (stand 325585226489998:
            # oracle CON=0.4476, jl was 1.0 ⇒ aspen i19 htg 1.69 vs oracle 0.756). RHCON=1 (no REUSCORR).
            con = exp(c.htg_cor_small[sp])
            h2 = h1 + htgrl * (kpj / regyr) * con
            wk3[i] = h2
            d1s = wk5[i]                               # subcycle-START DBH (regent.f:461 D1=WK5(I)) for the density feedback
            d2 = _tt_smdgf(esp, h2, cr, pccf)          # SMDGF gets the RAW point CCF (regent.f:574), not stand relden
            d2 < TT_RG_DIAM[sp] && (d2 = TT_RG_DIAM[sp])
            wk5[i] = d2
            # ★TT M331D woodland over-growth ROOT FIX (2026-09-03): small-tree density FEEDBACK (regent.f:581-593).
            # During subcycle J<NPER, each grown small tree (D<3, grown H2>4.5) adds its CCF/BA INCREASE to the NEXT
            # subcycle's projected density RDNEXT(J+1)/BANEXT(J+1). For the dense-seedling M331D woodland this
            # DOMINATES the density rise (stand 335: large-tree DO6 gives Δ1.24, the small-tree feedback adds ~2.7
            # more → RDNEXT(2) 28.7→32.7, PPCCF 1.137). jl omitted it, so PPCCF stayed ≈1.04 ⇒ SMHTGF saw too-low
            # CCF ⇒ over-grown height/DBH. `pr` supplies the ·P that CCFCAL folds in (tt_tree_ccf returns per-tree).
            if j < nper && d < 3.0f0 && h2 > 4.5f0
                cc1f = tt_tree_ccf(sp, d1s); cc2f = tt_tree_ccf(sp, d2)
                rdnext[j+1] += ky * (cc2f - cc1f) / 10f0 * pr * kymort
                banext[j+1] += _TT_BACON * (d2 * d2 - d1s * d1s) * pr * kymort
            end
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
        # HTG blend + size cap. ★TT M331D woodland over-growth ROOT FIX (2026-09-03): regent.f:731 sets
        # HTG(I)=0.0 for EVERY TTVAR tree at the top of loop-2 (the `IF(TTVAR)…HTG(I)=0.0` block), BEFORE the
        # XWT height blend at regent.f:799 (`HTG(K)=HTGR*(1-XWT)+XWT*HTG(K)`). So the large-tree height increment
        # is DISCARDED for TTVAR and the blend collapses to HTG = HTGR·(1−XWT). jl was blending in the large-tree
        # t.ht_growth[i] (htgf.f fires for D≥1.5), which — for trees in the [XMIN,XMAX)=[1.5,3.0) window — added
        # up to xwt·(large-tree htg) of spurious HEIGHT growth. MEASURED FVStt_dbg (AF stand 388908802489998
        # cyc1 i10, D=1.9): oracle HTGlarge=0 ⇒ HTG=0.8222·0.7333=0.603, jl blended large-htg=1.53 ⇒ HTG=1.011
        # (+68%). The current-cycle DBH is unaffected (SMDGF uses the subcycle H2, not this blended HTG — DG was
        # already bit-exact 0.11805), but the over-grown height compounds: next cycle SMDGF(H) over-predicts DBH,
        # driving the whole M331D AF/aspen/conifer woodland over-growth. (sp14 MM is UTVAR in FVS with a separate
        # POTHTG path; jl approximates it here via aspen — zeroing the large-tree term matches the TTVAR aspen it
        # is modeled on.)
        htgr = wk3[i] - h; htgr < 0.0f0 && (htgr = 0.0f0)
        htg = htgr * (1.0f0 - xwt)
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # DG (regent.f:923-942): DK=smdgf(grown H)=wk5, DKK=smdgf(ORIGINAL H), DG=(DK−DKK)·bark → DDS → DG. HK≥4.5.
        hk = h + htg
        dfl = d < TT_RG_DIAM[sp] ? TT_RG_DIAM[sp] : d       # regent.f:703 D floored to DIAM(sp)
        dgk = 0f0
        if hk >= 4.5f0
            pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100f0
            # regent.f:925 TTVAR branch: DKK = SMDGF(HT(I)) — the ORIGINAL-height DBH — UNCONDITIONALLY
            # (there is NO `H<4.5 → DKK=D` guard here; that guard lives only in the CIVAR/UTVAR branches at
            # :823-838). The HK≥4.5 gate above already suppresses DG while the GROWN height is sub-breast-height,
            # so DKK never needs the `=D` fallback. jl's earlier #191 fix mis-imported the CIVAR rule to TTVAR:
            # for a seedling whose ORIGINAL H<4.5 but GROWN HK≥4.5 it used DKK=D (≪ smdgf(H)) ⇒ (wk5−DKK)
            # over-counted the DBH increment ⇒ the DBH jumped past the CCF BREAK=1" ⇒ the ~10× CCF blow-up
            # (worst-col CCF in the TT dig cluster). MEASURED vs FVStt_g16 (aspen 31353347010690 cyc2):
            # oracle DKK=smdgf(4.455)=0.582 (DG≈0.65, d→0.85) vs jl DKK=d=0.203 (DG≈1.0, d→1.2).
            dkk = _tt_smdgf(_tt_rg_esp(sp), h, Float32(t.crown_pct[i]), pccf)   # DBH from ORIGINAL height (MM→AS)
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
        # DG is NOT XWT-blended. FVS regent.f blends only the HEIGHT increment with XWT (regent.f:799
        # HTG(K)=HTGR*(1-XWT)+XWT*HTG(K)); the DIAMETER increment has no such blend — for D<BKPT the regent
        # small-tree DG (dgk) fully REPLACES the large-tree dgf DG (regent.f:935-939), and for D≥BKPT the tree
        # takes GO TO 23 (regent.f:809) keeping the pre-computed large-tree DG. jl previously did
        # dgk*(1-xwt)+xwt*large_tree_DG, which mixed in up to ~xwt of the (larger) dgf DG. On mature ref stands
        # the [XMIN,XMAX) trees are few so this rounded out, but on M331D woodland/seedling stands (many trees
        # in 1.5–3"), the XMIN/XMAX correction to [1.5,3.0] pushed xwt→~0.93 at the window top, over-blending
        # the large-tree DG and over-growing DBH/BA/QMD (measured net-regression vs baseline). BKPT=XMAX=3.0
        # for the TT default conifers, so within this loop (d<XMAX) D<BKPT always holds ⇒ pure regent DG.
        if d < TT_RG_BREAK[sp]
            t.diam_growth[i] = dgk
        end
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
            # regent.f:822 `IF(D.GE.BKPT) GO TO 23`: at/above the DBH breakpoint the small-tree DBH increment is
            # skipped and the large-tree DG (already in diam_growth) stands — height still grows via the UTVAR
            # POTHTG above. NC/OH (15,18) BREAK=1"; junipers/BI/MC BREAK=99 (always < ⇒ always regent DBH).
            d >= TT_RG_BREAK[sp] && continue
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
        # birth-cycle ZRAND: esgent.f → REGENT → SMHTGF draws ZRAND(I) for the new tree and stores it,
        # so the NEXT cycle's small_tree_growth persists the same deviate (t.tree_random inheritance).
        zrand = 0.0f0
        if dgsd >= 1.0f0
            zr0 = t.tree_random[i]
            if zr0 == 0f0 || zr0 == -999f0
                while true; zrand = bachlo(s.rng, 0.0f0, 1.0f0); (-2.0f0 <= zrand <= 2.0f0) && break; end
                t.tree_random[i] = zrand
            else
                zrand = zr0
            end
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
