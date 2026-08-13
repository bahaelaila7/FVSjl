# =============================================================================
# diameter_growth.jl (southeastalaska) — AK large-tree DDS (ak/dgf.f). Chunk 3 (beachhead).
#
# Hooks:
#   ak_dgcons!(s)          — per-species per-stand constants (ak/dgf.f ENTRY DGCONS):
#                            ATTEN=OBSERV; DGCON=0 (+ ln(COR2) if READCORD/LDCOR2); AK bark.
#   dgf!(s, ::SoutheastAlaska) — per-tree WK2 = DDS (ak/dgf.f main body).
#
# The AK model is a SINGLE Wykoff equation for ALL species (no per-species DDS branches):
#   BASEDG = exp(b1 + b2·D² + b3·ln(D) + b4·PBAL + b5·PRD + b6·ln(CR)
#                 + b7·ELEV + b8·SLOPE + b9·SLOPE·cos(ASPECT) + b10·ln(SI))   [annual, DOB]
#   with ELEV=elevation·100, SLOPE=slope·100, PRD = point ZeideSDI/point maxSDI, CR = crown %.
# then a multiplicative PERMAFROST modifier (species 4:7,13,16:23):
#   PFMOD = exp(PFCON+[PFPRES if LPERM]+PFDSQ·D²+PFLD·ln(D)+PFDBAL·PBAL+PFRD·PRD+PFLNCR·ln(CR)
#               +PFEL·ELEV+PFSLOP·SLOPE+PFSASP·SLOPE·cos(ASPECT)) / BASEDG
#   LPERM ⇒ cap PFMOD≤1 (permafrost slows growth); !LPERM ⇒ floor PFMOD≥1. Else species PFMOD=1.
# then period + species-multiplier + inside-bark conversion:
#   DGPRED = YR·BASEDG·PFMOD·ak_dgmult(sp);  BRAT=ak_bratio(sp,D)
#   DDS = ln((D+DGPRED)²·BRAT² − D²·BRAT²) + COR(sp) + DGCON(sp);  floor −9.21;  WK2 = DDS.
#
# VALIDATED BIT-EXACT vs the live FVSak `DEBUG DGF` dump on akt01 (54 records, 0 mismatch,
# WK2 rel-err 0.0). akt01 is all non-permafrost species (SS/WH/YC/RC/LP/MH) with DGRD=0 ⇒
# PRD & PFMOD do not enter; the PERMAFROST + PRD point-Zeide path is a LATER run (the engine
# lacks the per-point Zeide-SDI density SDICAL/SDICLS machinery ak/dgf.f needs for PRD).
# =============================================================================

# ak/bratio.f — 3 bark-equation types (a=AK_BARK_A, b=AK_BARK_B, type=AK_BARK_TYPE).
function ak_bratio(sp::Int, d::Float32)
    d <= 0f0 && return 0.99f0
    typ = AK_BARK_TYPE[sp]; a = AK_BARK_A[sp]; b = AK_BARK_B[sp]
    br = if typ == 1
        dbt = a * d^b; (d - dbt) / d           # DBT=a·DOB^b, BRATIO=(D−DBT)/D
    elseif typ == 2
        (a + b * d) / d                        # BRATIO=(a+b·DOB)/DOB
    elseif typ == 3
        (a * d^b) / d                          # DIB=a·DOB^b, BRATIO=DIB/D
    else
        0.9f0
    end
    br > 0.99f0 && (br = 0.99f0)
    br < 0.80f0 && (br = 0.80f0)
    return br
end

# ak/dgf.f ENTRY DGCONS — ATTEN=OBSERV; DGCON=0 unless READCORD (LDCOR2) adds ln(COR2).
function ak_dgcons!(s::StandState)
    c = s.calib; ctl = s.control
    @inbounds for sp in 1:23
        c.atten[sp] = AK_OBSERV[sp]
        dgcon = 0f0
        (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) && (dgcon += log(ctl.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
        # AK bark is applied inside dgf! via ak_bratio; expose (a,b) for the engine's DBH-update
        # step as the type-appropriate linear surrogate is NOT valid for AK — bark_a/b unused here.
    end
    return s
end

# ak/dgf.f main body — per-tree WK2 = DDS.
function dgf!(s::StandState, ::SoutheastAlaska)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    wk2 = view(s.scratch.wk, 2, :)
    yr = s.control.year                       # YR (=10)
    temel  = p.elevation * 100f0              # TEMEL = ELEV·100
    temslp = p.slope * 100f0                  # TEMSLP = SLOPE·100
    temsasp = temslp * cos(p.aspect)          # TEMSASP = TEMSLP·cos(ASPECT)
    lperm = false                             # LPERM (ak PERMAFROST keyword) — not yet wired; default off (later run)
    # PRD = point ZeideSDI/point maxSDI (ak/dgf.f SDICAL→XMAXPT + SDICLS→ZRD, computed ONCE per DGF call).
    # dgf.f: SDICAL(IWHO=2) fills XMAXPT (IWHO-independent); the SDICLS loop over points (IWHO=1, JSPEC=0,
    # DLO=0/DHI=500) fills ZRD(pt). Inert for the DGRD=0 coastal species (AK_DGRD[sp]=0); load-bearing for
    # the interior/permafrost/hardwood species (4-7,13-23). Computed on the CURRENT t.dbh at call time
    # (backdated during LSTART calibration, current during growth), matching FVS's per-pass SDICAL/SDICLS.
    xmaxpt, zrd, _ = ak_point_zeide!(s)
    npt_prd = length(xmaxpt)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        d2 = d * d
        cr = Float32(t.crown_pct[i])          # CR = REAL(ICR(I)) — crown ratio in PERCENT
        pbal = dens.point_bal[i]              # PBAL = PTBALT(I)
        ip = Int(t.plot_id[i])                # PRD = ZRD(pt)/XMAXPT(pt); 0 if the point has no BA (XMAXPT≤0)
        prd = (1 <= ip <= npt_prd && xmaxpt[ip] > 0f0) ? zrd[ip] / xmaxpt[ip] : 0f0
        ssite = p.sp_site_index[sp]           # SSITE = SITEAR(ISPC)
        dgcomp1 = AK_DGEL[sp]*temel + AK_DGSLOP[sp]*temslp + AK_DGSASP[sp]*temsasp + AK_DGLNSI[sp]*log(ssite)
        dgcomp2 = AK_DGDISQ[sp]*d2 + AK_DGLD[sp]*log(d) + AK_DGDBAL[sp]*pbal + AK_DGRD[sp]*prd + AK_DGLNCR[sp]*log(cr)
        basedg = exp(AK_DGCONB1[sp] + dgcomp2 + dgcomp1)
        # permafrost modifier
        pfmod = 1.0f0
        if sp in AK_PERM_SP
            pfcomp1 = AK_PFEL*temel + AK_PFSLOP*temslp + AK_PFSASP*temsasp
            pfcomp2 = AK_PFDSQ*d2 + AK_PFLD[sp]*log(d) + AK_PFDBAL[sp]*pbal + AK_PFRD[sp]*prd + AK_PFLNCR*log(cr)
            if lperm
                pfmod = exp(AK_PFCON[sp] + AK_PFPRES + pfcomp2 + pfcomp1) / basedg
                pfmod > 1f0 && (pfmod = 1f0)
            else
                pfmod = exp(AK_PFCON[sp] + pfcomp2 + pfcomp1) / basedg
                pfmod < 1f0 && (pfmod = 1f0)
            end
        end
        dgpred = yr * basedg * pfmod * ak_dgmult(sp)
        brat = ak_bratio(sp, d)
        tempd1 = d * brat                     # current DIB
        tempd2 = (d + dgpred) * brat          # grown DIB
        dds = log(tempd2*tempd2 - tempd1*tempd1) + c.dg_cor[sp] + c.dg_const[sp]
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end

# PRD point-Zeide relative density is now wired directly in dgf! via ak_point_zeide! (crown.jl),
# which ports ak/sdical.f SDICAL(XMAXPT) + SDICLS(ZRD). PRD = ZRD(pt)/XMAXPT(pt) per tree's point.
