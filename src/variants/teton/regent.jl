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
# difference (~1-2%/tree), the known regent model-version straddle. 2026-09-27: the FINT·DGMAX cap itself is
# ESTAB-only (regent.f:962 sits inside IF(LESTB)); the regular-cycle TTVAR DG has no DGMX cap (only DGBND), and the
# UTVAR cap is the raw DGMAX(ISPC) (regent.f:1039).
# RAW per-year DGMAX (tt/regent.f:175-177 DATA DGMAX, verbatim). TTVAR ESTAB growth caps DG at
# DGMX=FINT·DGMAX(ISPC) (regent.f:684, applied at :962); UTVAR caps at DGMAX(ISPC) itself (regent.f:1039).
const TT_RG_DGMAX_RAW = Float32[0.2, 0.2, 0.2, 2.0, 0.2, 0.2, 0.2, 0.2, 0.2, 99.0, 2.0, 2.0, 2.0, 2.5, 2.5, 2.0, 0.2, 2.5]
const _TT_REGYR = 5.0f0
const _TT_BACON = 0.005454154f0

@inline _tt_smdg_alt(sp::Int) = sp == 3 || (5 <= sp <= 9)   # DF/BS/AS/LP/ES/AF use the alternate SMDGF form
@inline _tt_rg_default(sp::Int) = sp <= 3 || sp == 5 || sp == 6 || (7 <= sp <= 9) || sp == 17  # TTVAR regent-handled (WB/LM/DF/BS/AS/LP/ES/AF/OS; regent.f:388 CASE 1:3,5:9,17)
@inline _tt_rg_esp(sp::Int) = sp == 14 ? 6 : sp   # (legacy) MM→AS aspen-coef map; unused now MM is UTVAR-routed
# ★TT M331D MM over-growth ROOT FIX (2026-09-03): MM(14) is UTVAR (regent.f:396 CASE 4,11:16,18), NOT the
# aspen/TTVAR path. It was mis-routed through the default+aspen path (SMHTGF CASE(6) +5-yr Sheppard SUBCYCLED,
# no ·0.75, no RSIMOD, aspen SMDGF DBH), which over-grew MM height ~1.5× (H2≈15 vs FVS 10.2) then compounded
# into DBH/BA. FVS MM = regent.f:464 FINDAG aspen-height (AG2=SITAGE+10, ·RSIMOD·0.75, applied ONCE) + Wykoff
# H-D DBH (regent.f:907, HT1/HT2(14)). Route it through the UTVAR pass with its own FINDAG-height branch.
@inline _tt_rg_utvar(sp::Int)   = sp == 4 || sp == 11 || sp == 12 || sp == 13 || sp == 14 || sp == 15 || sp == 16 || sp == 18   # PM/UJ/RM/BI/MM/NC/MC/OH: UTVAR regent (regent.f:396 CASE 4,11:16,18)
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

# tt/regent.f:820-1056 (non-LESTB, HK≥4.5) — the CIVAR (PP) and UTVAR (PM/UJ/RM/BI/MM/NC/MC/OH) DBH increment of
# record K from the central's D=DBH(I), H=HT(I) and the record's grown HK=H+HTG(K). `bark` is REGENT's carried BARK
# local: the UTVAR CASE(4,11,12,14,15,18) DG=(DK−DKK)·BARK·XRDGRO (regent.f:1024-1027) reads whatever the last
# assignment left there (the previous record's regent.f:1035 BRATIO, or pass 1's regent.f:461 one) — only CIVAR
# (regent.f:1037) and CASE(13,16) (regent.f:999) recompute it before use. Returns (DG(K), BARK after the record).
# DBH(K)=D for every copy: tt/dgdriv.f:270 `DBH(ITRIPU)=DBH(I)` tripled DBH before REGENT runs.
# sj = SITEAR(ISPC); lhtdrg/iabflg/aa = LHTDRG, IABFLG and the calibrated Wykoff intercept AA of the species.
function _tt_rg_cu_dg(sp::Int, civar::Bool, d::Float32, h::Float32, hk::Float32, htg::Float32,
                      bark::Float32, xrdgro::Float32, dgmx::Float32, scale2::Float32,
                      sj::Float32, lhtdrg::Bool, iabflg::Integer, aa::Float32)::Tuple{Float32,Float32}
    local dk::Float32, dkk::Float32
    if civar
        # regent.f:821-829 CI PP height-diameter
        dk = fexp(-1.10700f0 + 0.830144f0 * flog(hk))
        dkk = h <= 4.5f0 ? d : fexp(-1.10700f0 + 0.830144f0 * flog(h))
    elseif sp == 4 || sp == 11 || sp == 12
        # regent.f:833-838 PM/UJ/RM: DKK floors to 0.1 BEFORE the H<4.5 ⇒ DKK=D override
        dk = (hk - 4.5f0) * 10f0 / (sj - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
        dkk = (h - 4.5f0) * 10f0 / (sj - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
        h < 4.5f0 && (dkk = d)
    elseif sp == 13 || sp == 16
        # regent.f:853-901 BI/MC (SO/WC origin): linear dub, replaced by the inventory Curtis-Arney equation whenever
        # the Wykoff calibration is off or did not happen (grinit LHTDRG(13/16)=.FALSE. ⇒ always).
        dkk = 3.1020f0 + 0.0210f0 * h; dkk < 0f0 && (dkk = d)
        dk = 3.1020f0 + 0.0210f0 * hk; dk < dkk && (dk = dkk + 0.01f0)
        if !lhtdrg || iabflg == 1
            p2, p3, p4 = sp == 16 ? (1709.7229f0, 5.8887f0, -0.2286f0) : (76.5170f0, 2.2107f0, -0.6365f0)
            hat3 = 4.5f0 + p2 * fexp(-p3 * fpow(3.0f0, p4))
            dk = _tt_bimc_dk(hk, hat3, p2, p3, p4)
            dkk = h <= 4.5f0 ? d : _tt_bimc_dk(h, hat3, p2, p3, p4)
        end
    else
        # regent.f:903-913 MM/NC/OH Wykoff: BX=HT2, AX=HT1 (IABFLG=1) or the calibrated AA (IABFLG=0)
        bx = TT_HT2[sp]; ax = iabflg == 1 ? TT_HT1[sp] : aa
        dk = bx / (flog(hk - 4.5f0) - ax) - 1f0; dk < 0.1f0 && (dk = 0.1f0)
        dkk = h <= 4.5f0 ? d : bx / (flog(h - 4.5f0) - ax) - 1f0
    end
    local dg::Float32
    if civar                                        # regent.f:975-983
        h < 4.5f0 && (dkk = d)
        if dk < 0f0 || dkk < 0f0
            dg = htg * 0.2f0 * bark * xrdgro; dk = d + dg
        else
            dg = (dk - dkk) * bark * xrdgro
        end
        dg < 0f0 && (dg = 0f0)
    elseif sp == 13 || sp == 16                     # regent.f:997-1016
        h < 4.5f0 && (dkk = d)
        bark = tt_bratio(sp, d)
        if dk < 0f0 || dkk < 0f0
            dg = htg * 0.2f0 * bark * xrdgro; dk = d + dg
        else
            dg = (dk - dkk) * bark * xrdgro
        end
        (lhtdrg && iabflg == 0) && (dg = 0.1f0 * htg * xrdgro)
        dg < 0f0 && (dg = 0.1f0)
        dg > dgmx && (dg = dgmx)
    else                                            # regent.f:1018-1025 CASE(4,11,12,14,15,18) — carried BARK
        if dk < 0f0 || dkk < 0f0
            dg = htg * 0.2f0 * bark * xrdgro; dk = d + dg
        else
            dg = (dk - dkk) * bark * xrdgro
        end
    end
    dg < 0f0 && (dg = 0f0)
    bark = tt_bratio(sp, d)                         # regent.f:1035 BARK=BRATIO(ISPC,DBH(K),HT(K))
    if civar
        dg = (dk - dkk) * bark * xrdgro             # regent.f:1037 recomputed (the 0 floor above does not survive)
    else
        dg > dgmx && (dg = dgmx)                    # regent.f:1039 DGMX=DGMAX(ISPC) (not FINT-scaled for UTVAR)
    end
    dds = (dg * (2f0 * bark * d + dg)) * scale2     # regent.f:1041-1046 DDS period/bark conversion
    arg = (d * bark) * (d * bark) + dds
    dg = sqrt(max(arg, 0f0)) - bark * d
    return (dg, bark)
end

# tt/regent.f loop 1 UTVAR height increment HTGRL (regent.f:463-491, 531-545) for a record of height H, crown CR (%).
# MM(14): FINDAG aspen inverse-height age (findag.f CASE 6,14 — metric Sheppard curve), AG2=SITAGE+10,
# HTGRL=(H(AG2)−H(SITAGE))/(2.54·12)·RSIMOD·CON·0.75 with RSIMOD from the SITERANGE-clamped SI (regent.f:427-431).
# Others: POTHTG=((SJ/5)·(SJ·1.5−H)/(SJ·1.5))·0.83, VIGOR=(150·X³·exp(−6X))+0.3 (X=CR/100) cut by two-thirds for
# PM/UJ/RM, HTGRL=POTHTG·PCTRED·VIGOR·CON. No floor here — the 0.1-ft floor follows the loop-2 ZZRAN (regent.f:778).
@inline function _tt_utvar_htgrl(sp::Int, h::Float32, cr::Float32, sj::Float32, rsimod::Float32,
                                 pctred::Float32, con::Float32)::Float32
    if sp == 14
        sitage = (h * 2.54f0 * 12f0 / 26.9825f0)^(1f0 / 1.1752f0)
        hite1  = 26.9825f0 * sitage^1.1752f0
        hite2  = 26.9825f0 * (sitage + 10f0)^1.1752f0
        return (hite2 - hite1) / (2.54f0 * 12f0) * rsimod * con * 0.75f0
    end
    pothtg = ((sj / 5f0) * (sj * 1.5f0 - h) / (sj * 1.5f0)) * 0.83f0
    x = cr / 100f0
    vigor = (150f0 * x * x * x * exp(-6f0 * x)) + 0.3f0
    vigor > 1f0 && (vigor = 1f0)
    (sp == 4 || sp == 11 || sp == 12) && (vigor = 1f0 - ((1f0 - vigor) / 3f0))
    return pothtg * pctred * vigor * con
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

# tt/dgdriv.f DO 220 DG(I) — measured (capped at inside-bark DBH when IDG<2), 0 at HT<=4.5, else the dub from the
# second DGF(WK3) (COR final; dub_wk2/dub_wk3 stash). Read by the LSTART REGCAL DO 49.
@inline function tt_do220_dg(s::StandState, i::Int, dcur::Float32)::Float32
    t, c = s.trees, s.calib
    sp = Int(t.species[i])
    bark = tt_bratio(sp, dcur)
    if t.diam_growth[i] > 0f0 && t.height[i] > 4.5f0
        dg = t.diam_growth[i]
        (s.control.growth_idg < 2 && dg > dcur * bark) && (dg = dcur * bark)
        return dg
    elseif t.height[i] <= 4.5f0 || length(c.dub_wk2) < i
        return 0f0
    end
    sc = s.control.growth_fint / 10f0
    dd = c.dub_wk3[i] * bark
    dub = sqrt(dd * dd + fexp(c.dub_wk2[i] + t.old_random[i]) * sc) - dd
    dub > dd && (dub = dd)
    return dg_bound(nothing, nothing, sp, dcur, dub, s.control.sp_size_cap)
end

"""
    tt_regent_hcor_init!(s, isct, ind1, saved_dbh, temba, temccf, pccfv, avh)

tt/regent.f:1095-1380 — the LSTART small-tree HEIGHT calibration (REGENT(.FALSE.,1) from tt/cratet.f:713), for every
sub-model. jl had only a UTVAR-woodland block (PM/UJ/RM/BI/MC, with POTHTG pinned at H=0) and a deterministic
aspen block (ZRAND=0, no draws); the TTVAR conifers (1-3,5,7-9,17), CIVAR PP (10) and UTVAR 15/18 were never
calibrated, and the calibration SMHTGF ZRAND draws were missing from the main stream. Per species with >= NCALHT(5)
sub-5" records carrying a measured HTG, CORNEW = Σ(HTG·SCALE3·P)/Σ(EDH·P), HCOR = ln(CORNEW) trapped to [0.0821,12.1825]:
- TTVAR (1:3,5:9,17): tt/smhtgf.f over NPER subcycles around the persistent ZRAND (t.tree_random; drawn on the main
  stream when unset), floored at 0.1 (ZRAND reset); aspen (6) its Sheppard form on the CURRENT HT (FINDAG) ×0.75, then
  RSIMOD; EDH = (HK−H)·RHCON.
- CIVAR (10): 2.764559 − 0.009643·BA + 0.025303·RCR², with FVS's second in-loop H backdate (and drop at H<0.01), ×RHCON.
- UTVAR (4,11:16,18): POTHTG = (SJ/5)(SJ·1.5−H)/(SJ·1.5)·0.83 evaluated ONCE per species with FVS's STALE H — the value
  left by the last record visited before it (0 only if none) — ·PCTRED·VIGOR (cut for 4/11/12); MM (14) Sheppard; ×0.5.
SCALE3 = 10/FINTH: the SELECT CASE on ISPC reads the index DO 45 left at MAXSP+1 ⇒ CASE DEFAULT for every species.
NTYR = IFINTH; PPCCF on TEMCCF; TEMBA/TEMCCF/PCCF from tt/cratet.f:243's backdating DENSE; AVH = AVHT40 (:653);
t.dbh is the backdated WK3 here (DO 49 uses it with the measured DG).
"""
function tt_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector, temba::Float32, temccf::Float32,
                              pccfv::AbstractVector{Float32}, avh::Float32)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    ctl = s.control
    ctl.growth_ifinth == 0 && return s                  # regent.f:1104 IF(IFINTH.EQ.0) GOTO 100
    rhcon = ones(Float32, 18)                           # tt REGCON: 1.0, or RCOR2 under READCORR
    if ctl.regh_cor2_on
        @inbounds for sp in 1:min(18, length(ctl.regh_cor2)); ctl.regh_cor2[sp] > 0f0 && (rhcon[sp] = ctl.regh_cor2[sp]); end
    end
    temccf <= 0f0 && (temccf = p.relative_density); temba <= 0f0 && (temba = p.basal_area)
    isempty(pccfv) && (pccfv = s.density.point_ccf)
    relden = temccf
    finth = ctl.growth_finth > 0f0 ? ctl.growth_finth : 5f0
    scale3 = 10f0 / finth                               # CASE DEFAULT (stale ISPC = MAXSP+1)
    ntyr = Int(ctl.growth_ifinth); iyr = 5
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for k in 1:nper
        if nn == 1; kper[k] = itot; break; end
        kper[k] = itot ÷ nn; itot -= kper[k]; nn -= 1
    end
    banext = fill(temba, nper); rdnext = fill(temccf, nper)
    if nper > 1                                         # DO 49 (regent.f:1123-1141) on the DGDRIV DO-220 DG
        @inbounds for i in 1:t.n
            d1 = t.dbh[i]; sp = Int(t.species[i]); pr = t.tpa[i]
            d2 = d1 + tt_do220_dg(s, i, saved_dbh[i]) / tt_bratio(sp, d1)
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            c1 = tt_tree_ccf(sp, d1) * pr; c2 = tt_tree_ccf(sp, d2) * pr
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn; banext[j] += Float32(k) * bi * pn
            end
        end
    end
    x = avh * (relden / 100f0); x > 300f0 && (x = 300f0)
    pctred = 1.11436f0 + x*(-0.011493f0 + x*(0.43012f-4 + x*(-0.72221f-7 + x*(0.5607f-10 - x*0.1641f-13))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    ihtg = ctl.growth_ihtg
    hst = 0f0                                           # FVS's REGENT-local H (static, 0 before any record)
    si6 = p.sp_site_index[6]
    rs6 = 0.5f0 * (1f0 + clamp((si6 - 30f0) / 70f0, 0f0, 1f0))   # regent.f:1233-1237 aspen RSIMOD
    @inbounds for sp in 1:18
        ttvar = (1 <= sp <= 3) || (5 <= sp <= 9) || sp == 17
        civar = sp == 10
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        slo = TT_SITELO[sp]; shi = TT_SITEHI[sp]
        pothtg = 0f0; rsimod = 0f0
        if !(ttvar || civar)                            # UTVAR species constants (regent.f:1173-1182) on the STALE H
            si = p.sp_site_index[sp]
            si > shi && (si = shi); si <= slo && (si = slo + 0.5f0)
            rsimod = 0.5f0 * (1f0 + (si - slo) / (shi - slo))
            sj = p.sp_site_index[sp]
            pothtg = ((sj / 5f0) * (sj * 1.5f0 - hst) / (sj * 1.5f0)) * 0.83f0
        end
        snp = 0f0; snx = 0f0; sny = 0f0; nh = 0
        for k in i1:i2
            i = Int(ind1[k])
            hg = t.ht_growth[i]
            hst = t.height[i]; ihtg < 2 && (hst -= hg)
            (saved_dbh[i] >= 5f0 || hst < 0.01f0) && continue
            hg < 0.001f0 && continue
            h = hst; hk = h; cr = Float32(t.crown_pct[i]); edh = 0f0
            ipccf = Int(t.plot_id[i])
            pccf_i = (1 <= ipccf <= length(pccfv)) ? pccfv[ipccf] : 0f0
            dropped = false
            for j in 1:nper
                rdj = rdnext[j]
                ppccf = temccf <= 0f0 ? 0f0 : (rdj - temccf) / temccf
                tpccf = pccf_i * ppccf
                tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
                if civar
                    iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9)
                    rcr = Float32(iicr)
                    ihtg < 2 && (hst -= hg; h = hst)               # regent.f:1253 — backdated AGAIN, every J
                    (saved_dbh[i] >= 5f0 || h < 0.01f0) && (dropped = true; break)
                    edh = 2.764559f0 - 0.009643f0 * temba + 0.025303f0 * rcr * rcr
                elseif ttvar                                       # tt/smhtgf.f
                    zr = t.tree_random[i]
                    if zr == -999f0 || zr == 0f0
                        z = 0f0
                        while true; z = bachlo(s.rng, 0f0, 1f0); (-2f0 <= z <= 2f0) && break; end
                        t.tree_random[i] = z
                    end
                    if saved_dbh[i] <= 0f0
                        edh = 0f0                                  # D<=0 ⇒ 0, bypassing the 0.1 floor
                    else
                        if sp == 6                                 # CASE(6): FINDAG age from the CURRENT HT(I)
                            sitage = (t.height[i] * 2.54f0 * 12f0 / 26.9825f0)^(1f0 / 1.1752f0)
                            htgr = (26.9825f0 * (sitage + 5f0)^1.1752f0 - 26.9825f0 * sitage^1.1752f0) / (2.54f0 * 12f0)
                            edh = (htgr + t.tree_random[i] * 0.1f0) * 0.75f0
                        else
                            edh = _tt_smhtgf(sp, h, cr, tpccf, t.tree_random[i], si6)
                        end
                        edh <= 0.1f0 && (edh = 0.1f0; t.tree_random[i] = -999f0)
                        sp == 6 && (edh *= rs6)
                    end
                    hk += edh
                else                                               # UTVAR
                    xv = cr / 100f0
                    vigor = 150f0 * xv^3 * exp(-6f0 * xv) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
                    (sp == 4 || sp == 11 || sp == 12) && (vigor = 1f0 - (1f0 - vigor) / 3f0)
                    if sp == 14
                        ag1 = (h * 12f0 * 2.54f0 / 26.9825f0)^0.8509f0
                        h2 = (26.9825f0 * (ag1 + 10f0)^1.1752f0) / (2.54f0 * 12f0)
                        edh = (h2 - h) * rsimod * rhcon[sp] * 0.75f0
                    else
                        edh = pothtg * pctred * vigor * rhcon[sp]
                    end
                    edh *= 0.5f0
                end
            end
            dropped && continue
            civar && (edh *= rhcon[sp])
            ttvar && (edh = (hk - h) * rhcon[sp])
            pr = t.tpa[i]
            snp += pr; snx += edh * pr; sny += hg * scale3 * pr; nh += 1
        end
        nh < 5 && continue                                         # NCALHT
        snx /= snp; sny /= snp
        cornew = sny / snx
        cornew <= 0f0 && (cornew = 1f-4)
        c.htg_cor_init[sp] = (cornew < 0.0821f0 || cornew > 12.1825f0) ? 0f0 : log(cornew)
    end
    return s
end

function small_tree_growth!(s::StandState, stash, ::Teton; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    n = t.n; n == 0 && return s
    cw = clim_wk4(s, Float32(current_cycle_year(s)) + fint / 2f0)   # CLGMULT WK4 (tt/regent.f:770); nothing ⇒ 1
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
            bark = tt_bratio(sp, d1)                  # regent.f:260 BRATIO (the PP power model, not the linear a+b·d)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = _TT_BACON * d1 * d1; b2 = _TT_BACON * d2 * d2
            cc1 = tt_tree_ccf(sp, d1) * pr; cc2 = tt_tree_ccf(sp, d2) * pr   # CCFCAL folds in ·P
            bi = (b2 - b1) / 10.0f0; ci = (cc2 - cc1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                # ★TT M331D woodland over-growth ROOT FIX (2026-09-03): tt/ccfcal.f returns CCFT=poly(D)·P (CCF ×
                # trees-per-acre). FVS's projection `RDNEXT+=K·CI/P·PN` (regent.f:274) nets ONE factor of P
                # (CI=(C2−C1)/10 already carries ·P; CI/P·PN=CI/P·P·0.985^K). jl's `tt_tree_ccf` returns the
                # per-tree CCF WITHOUT ·P, so dividing `ci/pr` dropped the density weighting entirely — rdnext
                # barely moved (28.7→28.9 vs oracle 28.7→32.7), so PPCCF≈1.007 vs oracle 1.137, so subcycle-2
                # SMHTGF saw too-low CCF ⇒ too-high height increment ⇒ small-tree DBH over-grew, ONE-DIRECTIONAL
                # and compounding across the woodland cluster. banext is already correct (BI=BACON·D² carries no
                # ·P, so `k·bi·pn` nets ·P). CI now carries ·P as CCFCAL's does, so the FVS expression is used as is.
                rdnext[j] += Float32(k) * ci / pr * pn; banext[j] += k * bi * pn
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
    # ZRAND is drawn INSIDE the subcycle loop (below), exactly where tt/smhtgf.f:70-73 draws it: on any call that
    # finds ZRAND(I)=-999 — the first subcycle of a new tree, or ANY later subcycle after a floored (≤0.1 ft)
    # increment reset it — with NO DGSD gate (SMHTGF has none; only REGENT's ZZRAN/DUBSCR draws test DGSD).
    # IND1 is the SPESRT species-major order with the post-TRIPLE lineage order inside a species (a triple walks as
    # copy1, original, copy2 — measured on FIA 2783239010690 cycle 2: live SMHTGF order I=167,54,168), NOT raw index
    # order; with raw order the J=1 draws landed on the wrong member of every triple.
    zorder = species_major_order(s)                 # IND1 (ISCT row ranges index into it)
    isct = s.control.sp_count_tab
    nsp = length(TT_RG_XMAX)
    cur_year = current_cycle_year(s)
    yr = htg_period(s.variant)                      # /CONTRL/ YR
    wk3 = Float32[t.height[i] for i in 1:n]         # subcycle height
    wk5 = Float32[t.dbh[i] for i in 1:n]            # subcycle DBH
    # PCTRED (regent.f:333-341): the UTVAR POTHTG density modifier from AVH and stand CCF (non-ESTAB call).
    xpr = p.avg_height * (relden / 100f0); xpr > 300f0 && (xpr = 300f0)
    pctred = 1.11436f0 + xpr * (-0.011493f0 + xpr * (0.43012f-4 + xpr * (-0.72221f-7 +
             xpr * (0.5607f-10 - xpr * 0.1641f-13))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # REGENT locals that CARRY from the subcycle loop into the growth loop (statics under -fno-automatic): loop 2
    # never reassigns CON or XRHGRO, so the CIVAR/UTVAR per-record ZZRAN step (regent.f:753, 765-770) uses the values
    # of the LAST species loop 1 processed; BARK is reassigned only on some loop-2 paths (see _tt_rg_cu_dg).
    con = 1f0; xrhgro = 1f0; bark = 0f0
    # KNOWN faithful gap: buildDir HTGR=POTHTG·PCTRED·VIGOR·CON (regent.f:350, RHCON=1 @ line 880). Tested a
    # conifer-only version (aspen sp6 exempt = Sheppard-already-final; kept aspen exact) but it BARELY moved DF
    # (2040 108→103, still 45% over live 71) ⇒ NOT the DF later-cycle residual (= large-tree dgf DF DG at DBH 4-6,
    # un-validatable: live fort.79 caps at DBH 2.0). Reverted — an untestable-for-TT-conifers change that doesn't
    # fix the visible residual. The faithful PCTRED·VIGOR·CON gap remains (conifer-only if ever added).
    # Loop 1 (regent.f:348-641, DO 17 J / DO 16 ISPC / DO 15 I3=IND1): TTVAR (1:3,5:9,17) SMHTGF subcycles;
    # CIVAR (10, PP) and UTVAR (4,11:16,18) take one full-cycle step at J=1 (regent.f:407).
    @inbounds for j in 1:nper
        rdj = rdnext[j]; kpj = Float32(kper[j])
        ky = 0; for m in 1:j; ky += kper[m]; end          # cumulative subcycle length (regent.f KY=KY+KPER(J))
        kymort = fpowi(0.985f0, ky)            # REAL**INTEGER (libgcc __powisf2)
        # ★TT M331D woodland over-growth FIX (2026-09-03): PPCCF — the subcycle proportional point-CCF
        # adjustment (regent.f:352-353 `PPCCF=1.0+(RDJ-RELDEN)/RELDEN`). FVS scales each tree's point CCF by
        # the PROJECTED stand-density increase for subcycle J before feeding it to SMHTGF as TPCCF. jl omitted
        # it — every subcycle used the RAW point CCF — so subcycle J≥2 saw a lower CCF than FVS (density always
        # grows ⇒ PPCCF>1), yielding a HIGHER SMHTGF height increment (BETA1/BETA2 fall with CCF) ⇒ over-grown
        # small-tree height ⇒ over-grown SMDGF DBH, ONE-DIRECTIONAL and compounding across the woodland cluster.
        # MEASURED vs FVStt (stand 335 cyc1 i6): oracle J=2 TPCCF=65.18=57.33·1.1369 HTGRL=1.079 vs jl raw
        # TPCCF=57.33 HTGRL=1.121. SMDGF (DBH, regent.f:574) keeps the RAW PCCF — only the height model uses PPCCF.
        ppccf = relden > 0f0 ? 1f0 + (rdj - relden) / relden : 0f0
        for sp in 1:nsp
            isct[sp, 1] == 0 && continue
            ttvar = _tt_rg_default(sp); civar = sp == 10; utvar = _tt_rg_utvar(sp)
            # #205 (2026-08-13): tt/regent.f:407 — aspen(6)/CIVAR(10)/UTVAR apply the small-tree height increment
            # ONLY on the FIRST subcycle; only TTVAR conifers subcycle across all J. jl formerly subcycled aspen
            # every j ⇒ for a 10-yr cycle (nper=2) it DOUBLE-applied the Sheppard SMHTGF increment ⇒ aspen small-tree
            # height/DBH over-grew ~1.9× (11796095010690 +34% BA). At J=1 aspen SCALE=KPER(1)/REGYR=kpj/regyr matches.
            ((sp == 6 || utvar || civar) && j > 1) && continue
            # regent.f:409-423: species constants — set even when every record of the species is ≥XMAX.
            # CON = RHCON·exp(HCOR) (RHCON=1, no REUSCORR): the REGENT small-tree HEIGHT self-calibration. For
            # aspen(6)/MM(14) HCOR = the CORNEW calibration seeded in tt_regent_hcor_init! and attenuated per-cycle
            # in diameter_growth!; without it aspen/MM held CON=1 ⇒ ~3.6× small-tree HEIGHT over-growth on the M331D
            # woodland cluster (stand 325585226489998: oracle CON=0.4476).
            xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
            xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
            con = exp(c.htg_cor_small[sp])
            si = p.sp_site_index[sp]
            si > TT_SITEHI[sp] && (si = TT_SITEHI[sp]); si <= TT_SITELO[sp] && (si = TT_SITELO[sp] + 0.5f0)
            rsimod = 0.5f0 * (1f0 + (si - TT_SITELO[sp]) / (TT_SITEHI[sp] - TT_SITELO[sp]))
            sj = p.sp_site_index[sp]
            for i3 in Int(isct[sp, 1]):Int(isct[sp, 2])
                i = zorder[i3]
                d = t.dbh[i]; pr = t.tpa[i]
                (d >= TT_RG_XMAX[sp] || pr <= 0f0) && continue
                h = t.height[i]; h1 = wk3[i]; cr = Float32(t.crown_pct[i])
                bark = tt_bratio(sp, d)                        # regent.f:461 BARK=BRATIO(ISPC,D,H)
                if civar
                    # regent.f:505-517 CI PP: HTGRL=2.764559−0.009643·BA+0.025303·RCR², RCR=crown-ratio class 1-9;
                    # regent.f:551 H2=H1+HTGRL·SCALE·XRHGRO·CON with SCALE=NTYR/REGYR; DBH is left for loop 2 (:569-572).
                    iicr = div(Int(t.crown_pct[i]) - 1, 10) + 1; iicr > 9 && (iicr = 9)
                    rcr = Float32(iicr)
                    htgrl = 2.764559f0 - 0.009643f0 * ba + 0.025303f0 * rcr * rcr
                    wk3[i] = h1 + htgrl * (Float32(ntyr) / regyr) * xrhgro * con
                    continue
                elseif utvar
                    htgrl = _tt_utvar_htgrl(sp, h, cr, sj, rsimod, pctred, con)
                    h2 = h1 + htgrl * (Float32(ntyr) / yr)        # regent.f:555 H2=H1+HTGRL·SCALE, SCALE=NTYR/YR
                    wk3[i] = h2
                    # regent.f:600-611: NC/OH (and BI) at/above BKPT keep D; a sub-4.5' seedling gets the 0.001·H2 nudge.
                    ((sp == 13 || sp == 15 || sp == 18) && d >= TT_RG_BREAK[sp]) && continue
                    d2 = d
                    h2 <= 4.5f0 && (d2 = d + 0.001f0 * h2)
                    wk5[i] = d2
                    # regent.f:622-624: the UTVAR record's CCF/BA change feeds the next subcycle's density.
                    if j < nper
                        c1 = tt_tree_ccf(sp, d) * pr; c2 = tt_tree_ccf(sp, d2) * pr
                        rdnext[j+1] += Float32(ky) * (c2 - c1) / 10f0 * kymort
                        banext[j+1] += (_TT_BACON * d2 * d2 - _TT_BACON * d * d) * pr * kymort
                    end
                    continue
                end
                ttvar || continue
                pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100f0
                tpccf = pccf * ppccf; tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)   # PPCCF-adjusted; smhtgf clamps [25,300]
                esp = _tt_rg_esp(sp)
                # ★#158 ZRAND (tt/smhtgf.f:70-73): draw when ZRAND=-999 (0 = jl inventory default), persisted across
                # subcycles, cycles and tripled sub-records (the tripling copy list inherits t.tree_random). The draw
                # sits in the species-major subcycle loop, so a tree reset at J=1 redraws at J=2 in species order.
                zr = t.tree_random[i]
                if zr == 0f0 || zr == -999f0
                    z = 0f0
                    while true; z = bachlo(s.rng, 0.0f0, 1.0f0); (-2f0 <= z <= 2f0) && break; end
                    t.tree_random[i] = z
                end
                htgrl = _tt_smhtgf(sp, h1, cr, tpccf, t.tree_random[i], si6)
                # smhtgf.f:131-133: an increment ≤ 0.1 ft floors to 0.1 and resets ZRAND to -999 on EVERY call, so the
                # next SMHTGF call (the next subcycle, or next cycle) draws a fresh deviate.
                if htgrl <= 0.1f0
                    htgrl = 0.1f0
                    t.tree_random[i] = -999f0
                end
                h2 = h1 + htgrl * (kpj / regyr) * xrhgro * con   # regent.f:553 H2=H1+HTGRL·SCALE·XRHGRO·CON
                wk3[i] = h2
                d1s = wk5[i]                               # subcycle-START DBH (regent.f:461 D1=WK5(I)) for the density feedback
                d2 = _tt_smdgf(esp, h2, cr, pccf) * xrdgro # SMDGF gets the RAW point CCF (regent.f:574), not stand relden
                d2 < TT_RG_DIAM[sp] && (d2 = TT_RG_DIAM[sp])
                wk5[i] = d2
                # ★TT M331D woodland over-growth ROOT FIX (2026-09-03): small-tree density FEEDBACK (regent.f:581-593).
                # During subcycle J<NPER, each grown small tree (D<3, grown H2>4.5) adds its CCF/BA INCREASE to the NEXT
                # subcycle's projected density RDNEXT(J+1)/BANEXT(J+1). For the dense-seedling M331D woodland this
                # DOMINATES the density rise (stand 335: large-tree DO6 gives Δ1.24, the small-tree feedback adds ~2.7
                # more → RDNEXT(2) 28.7→32.7, PPCCF 1.137). `pr` supplies the ·P that CCFCAL folds in.
                if j < nper && d < 3.0f0 && h2 > 4.5f0
                    cc1f = tt_tree_ccf(sp, d1s); cc2f = tt_tree_ccf(sp, d2)
                    rdnext[j+1] += ky * (cc2f - cc1f) / 10f0 * pr * kymort
                    banext[j+1] += _TT_BACON * (d2 * d2 - d1s * d1s) * pr * kymort
                end
            end
        end
    end
    # Loop 2 (regent.f:650-1076, DO 30 ISPC / DO 25 I3=IND1): the cycle increment of every record, species-major.
    # TTVAR copies take the central's values (regent.f:1069-1081); CIVAR/UTVAR re-run label 918 per tripled copy
    # (L=1,2; K=ITRN+2I−2+L) with a fresh ZZRAN that COMPOUNDS into HTGR, the XWT blend with the copy's large-tree HTG,
    # SIZCAP, the copy's own DBH/DG and DGBND — the same per-copy shape as the UT/CI REGENT.
    scale2 = yr / fint                              # SCALE2 = YR/NTYR (period scaling of the DBH increment)
    ntyr10 = Float32(ntyr) / 10f0
    ltrip = stash !== nothing && !isempty(stash.htgU)
    @inbounds for sp in 1:nsp
        isct[sp, 1] == 0 && continue
        ttvar = _tt_rg_default(sp); civar = sp == 10; utvar = _tt_rg_utvar(sp)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        xmx = TT_RG_XMAX[sp]; xmn = TT_RG_XMIN[sp]; bkpt = TT_RG_BREAK[sp]; diam = TT_RG_DIAM[sp]
        dgmx = TT_RG_DGMAX_RAW[sp]                  # UTVAR DGMX=DGMAX(ISPC); TTVAR's FINT·DGMAX binds only under ESTAB
        cap = s.control.sp_size_cap[sp, 4]
        for i3 in Int(isct[sp, 1]):Int(isct[sp, 2])
            i = zorder[i3]
            d = t.dbh[i]
            (d >= xmx || t.tpa[i] <= 0f0) && continue
            h = t.height[i]
            if ttvar
                # HTG blend + size cap. ★TT M331D woodland over-growth ROOT FIX (2026-09-03): regent.f:731 sets
                # HTG(I)=0.0 for EVERY TTVAR tree before the XWT blend (regent.f:799), so the large-tree height increment
                # is DISCARDED and the blend collapses to HTG = HTGR·(1−XWT). MEASURED FVStt_dbg (AF stand
                # 388908802489998 cyc1 i10, D=1.9): oracle HTGlarge=0 ⇒ HTG=0.8222·0.7333=0.603.
                dfl = d < diam ? diam : d                   # regent.f:703 D floored to DIAM(sp)
                xwt = dfl <= xmn ? 0.0f0 : (dfl - xmn) / (xmx - xmn)
                htgr = wk3[i] - h
                htg = htgr * (1.0f0 - xwt)
                if h + htg > cap; htg = cap - h; htg < 0.1f0 && (htg = 0.1f0); end
                t.ht_growth[i] = htg
                hk = h + htg
                dgk = 0f0
                if hk >= 4.5f0
                    pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100f0
                    # regent.f:925 TTVAR: DKK = SMDGF(HT(I)) — the ORIGINAL-height DBH — UNCONDITIONALLY (the `H<4.5 →
                    # DKK=D` guard lives only in the CIVAR/UTVAR branches). MEASURED vs FVStt_g16 (aspen 31353347010690
                    # cyc2): oracle DKK=smdgf(4.455)=0.582 (DG≈0.65, d→0.85).
                    dkk = _tt_smdgf(_tt_rg_esp(sp), h, Float32(t.crown_pct[i]), pccf)
                    b1 = tt_bratio(sp, dfl)                 # regent.f:934 BARK=BRATIO(ISPC,D,HT(K)), D floored
                    dgk = (wk5[i] - dkk) * b1
                    dds = dgk * (2f0 * b1 * dfl + dgk) * scale2
                    dgk = sqrt(max((dfl * b1)^2 + dds, 0f0)) - b1 * dfl
                end
                # regent.f:1033-1047: the non-ESTAB tail runs for TTVAR too — floor at 0, BARK=BRATIO(ISPC,DBH(K),HT(K)),
                # and a SECOND DDS conversion (the identity when SCALE2=1). There is no DGMX cap here: TTVAR's
                # DGMX=FINT·DGMAX (regent.f:684) is applied only inside the ESTAB branch (regent.f:962).
                dgk < 0f0 && (dgk = 0f0)
                bark = tt_bratio(sp, d)
                dds = dgk * (2f0 * bark * dfl + dgk) * scale2
                dgk = sqrt(max((dfl * bark)^2 + dds, 0f0)) - bark * dfl
                (d + dgk) < diam && (dgk = diam - d)        # regent.f:1049 DIAM floor on DBH(K)+DG(K)
                dgk = dg_bound(nothing, nothing, sp, d, dgk, s.control.sp_size_cap)   # regent.f:1055 DGBND
                # DG is NOT XWT-blended: for D<BKPT the regent DG fully REPLACES the large-tree DG (BKPT=XMAX=3 for TTVAR).
                d < bkpt && (t.diam_growth[i] = dgk)
                # TTVAR copies (regent.f:1069-1081): DBH/DG/HT/HTG/ICR of the central — no per-copy spread.
                if ltrip && i <= length(stash.dgU)
                    stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
                    stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
                    stash.is_small[i] = true
                end
                continue
            end
            # CIVAR / UTVAR (regent.f:747-812 then 815-1056), record K = central then its two copies.
            htgr = wk3[i] - h                               # regent.f:700 HTGR=HK−H; compounds over the copies
            large_htg = t.ht_growth[i]
            xwt = d <= xmn ? 0f0 : (d - xmn) / (xmx - xmn)
            wk4i = cw === nothing ? 1f0 : cw[i]
            nrec = ltrip && i <= length(stash.htgU) ? 3 : 1
            for l in 0:(nrec - 1)
                zzran = 0f0
                if dgsd >= 1f0
                    while true
                        zzran = bachlo(s.rng, 0f0, 1f0)
                        (zzran <= 0.5f0 && zzran >= -2f0) && break
                    end
                end
                if civar
                    htgr = (htgr + zzran * 0.1f0 * (Float32(ntyr) / regyr)) * xrhgro * con   # regent.f:753
                else
                    # regent.f:763-783: CR-surrogate NC/OH take ZZRAN·0.2 and the CLGMULT WK4; floor 0.1 ft (Dixon 3/4/09)
                    htgr = (sp == 15 || sp == 18) ? (htgr + zzran * 0.2f0 * ntyr10) * xrhgro * wk4i :
                                                    (htgr + zzran * 0.1f0 * ntyr10) * xrhgro
                    htgr < 0.1f0 && (htgr = 0.1f0)
                end
                # regent.f:790-806: XWT blend with record K's large-tree HTG (htgf.f:744-760 gives each copy TEMHTG),
                # then SIZCAP.
                lh = l == 0 ? large_htg : (stash.htg_copy[i] ? (l == 1 ? stash.htgU[i] : stash.htgL[i]) : large_htg)
                htg = htgr * (1f0 - xwt) + xwt * lh
                if h + htg > cap; htg = cap - h; htg < 0.1f0 && (htg = 0.1f0); end
                # regent.f:814 IF(D.GE.BKPT) GO TO 23: the record keeps its large-tree DG (NC/OH BKPT=1, PP BKPT=3).
                below = d < bkpt
                dbhk = d; dgk = 0f0; direct = false
                if below
                    hk = h + htg
                    if hk < 4.5f0
                        # regent.f:817-819: DG(K)=0, DBH(K)=D+0.001·HK (no DIAM floor on this path)
                        dbhk = d + 0.001f0 * hk; direct = true
                    else
                        dgk, bark = _tt_rg_cu_dg(sp, civar, d, h, hk, htg, bark, xrdgro, dgmx, scale2,
                                                 p.sp_site_index[sp], s.control.ht_drag_sp[sp],
                                                 c.ht_dbh_iabflg[sp], c.ht_dbh_aa[sp])
                        (d + dgk) < diam && (dgk = diam - d)         # regent.f:1049 DIAM floor
                    end
                    dgk = dg_bound(nothing, nothing, sp, dbhk, dgk, s.control.sp_size_cap)   # regent.f:1055 DGBND
                end
                if l == 0
                    t.ht_growth[i] = htg
                    if below
                        direct && (t.dbh[i] = dbhk)
                        t.diam_growth[i] = dgk
                    end
                elseif l == 1
                    stash.htgU[i] = htg; stash.is_small[i] = true
                    below && (stash.dgU[i] = dgk; stash.dbhU[i] = dbhk)
                else
                    stash.htgL[i] = htg
                    below && (stash.dgL[i] = dgk; stash.dbhL[i] = dbhk)
                end
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
    si6 = p.sp_site_index[6]
    # REGENT(LESTB) walks the new records species-major (DO 16 ISPC; DO 15 I3) — the SMHTGF ZRAND draw order.
    @inbounds for i in filter(>(nstart), species_major_order(s))
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= TT_RG_XMAX[sp] || !_tt_rg_default(sp)) && continue
        h = t.height[i]; cr = Float32(t.crown_pct[i])
        pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 100.0f0
        tpccf = pccf; tpccf > 300.0f0 && (tpccf = 300.0f0); tpccf < 25.0f0 && (tpccf = 25.0f0)
        # birth-cycle ZRAND: esgent.f → REGENT → SMHTGF draws ZRAND(I) for the new tree and stores it,
        # so the NEXT cycle's small_tree_growth persists the same deviate (t.tree_random inheritance).
        # tt/smhtgf.f:70-73 has no DGSD gate; a ≤0.1-ft increment floors to 0.1 and resets ZRAND (:131-133).
        zrand = t.tree_random[i]
        if zrand == 0f0 || zrand == -999f0
            while true; zrand = bachlo(s.rng, 0.0f0, 1.0f0); (-2.0f0 <= zrand <= 2.0f0) && break; end
            t.tree_random[i] = zrand
        end
        esp = _tt_rg_esp(sp)                     # MM(14)→AS(6) coefficient mapping
        htgrl = _tt_smhtgf(esp, h, cr, tpccf, zrand, si6)
        htgrl <= 0.1f0 && (htgrl = 0.1f0; t.tree_random[i] = -999f0)
        # regent.f:552 H2=H1+HTGRL·SCALE·XRHGRO·CON, CON=RHCON·EXP(HCOR) (RHCON=1 for TT)
        htg = htgrl * subcyc * exp(c.htg_cor_small[sp]); htg < 0.0f0 && (htg = 0.0f0)
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        h2 = h + htg
        t.height[i] = h2
        if h2 >= 4.5f0
            # tt/regent.f:948-960 TTVAR REGENT(LESTB) sets DBH(K)=DK = SMDGF(HK) ABSOLUTELY. The former code
            # booked dbh = d + dgk/bratio, i.e. the LESTB=F growth-increment reconstruction (DG=(DK−DKK)*BARK,
            # DDS-rescaled ≈ d+DK−DKK) capped at DGMX=FINT·DGMAX. Since DKK=SMDGF(H) on the sub-breast-height
            # birth height ≠ D, that reconstruction sits below DK, and the DGMX cap clips it further, so the
            # synchronized PLANT cohort entered the next cycle under DK — a one-directional BA/QMD deficit.
            # TTVAR takes the plain DBH(K)=DK arm (no UTVAR/CIVAR DIAM-floor/+0.001*HK, regent.f:962-967); the
            # jl d2 floor to DIAM is retained (harmless, prevents a sub-DIAM birth record). Mirror of IE/UT/BM.
            d2 = _tt_smdgf(esp, h2, cr, pccf); d2 < TT_RG_DIAM[sp] && (d2 = TT_RG_DIAM[sp])
            t.dbh[i] = d2
        end
    end
    return s
end
