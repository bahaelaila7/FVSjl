# =============================================================================
# so/cwcalc.f — SouthCentralOregon (VARACD "SO", MAXSP=33) crown width (ft) for
# FMCBA's PERCOV. Bit-exact transcription of the CWEQN SELECT CASE dispatch keyed
# by SOMAP (so/cwcalc.f DATA SOMAP ~line 240), one 5-char code per species index:
#
#   1 WP 11905  2 SP 11705  3 DF 20205  4 WF 01505  5 MH 26403  6 IC 08105
#   7 LP 10805  8 ES 09305  9 SH 02105 10 PP 12205 11 WJ 06405 12 GF 01703
#  13 AF 01905 14 SF 01105 15 NF 02206 16 WB 10105 17 WL 07303 18 RC 24205
#  19 WH 26305 20 PY 23104 21 WA 31206 22 RA 35106 23 BM 31206 24 AS 74605
#  25 CW 74705 26 CH 35106 27 WO 81505 28 WI 31206 29 GC 63102 30 MC 47502
#  31 MB 47502 32 OS 12205 33 OH 31206
#
# Reference stand forest KODFOR=601 (DESCHUTES, R6). The R6 Crookston models apply
# a per-forest crown factor BF, keyed on FIASP = first 3 chars of the 5-char code
# (so/cwcalc.f CASE(601,799) ~line 501): '015'->1.044 '019'->0.936 '022'->1.301
# '081'->0.837 '073'->0.818 '117'->1.048 '122'->0.918 '202'->1.055 '263'->1.097,
# otherwise 1.0. BF is folded into the leading coefficient (e.g. 5.0312f0*1.044f0)
# and is applied ONLY inside the R6 Crookston model-2 blocks that literally use BF.
# The Donnelly (…06), R1 log (…03/…04-power), Bechtold-m2 (…02) and the MH piecewise
# forms carry NO BF in the Fortran, so their FIASP BF (e.g. NF '022'=1.301, WL
# '073'=0.818) is inert and NOT applied — matching so/cwcalc.f verbatim.
#
# Helpers: _cr_r6m2 (src/variants/centralrockies/crown.jl) = Crookston R6 model-2
# a·D^b·H^c·CL^dd·(BA+1)^e·EXP(EL)^f with EL floor/cap, OMIND=1 small-D scaling and
# CW cap; fexp/flog = gfortran-companion EXP/ALOG (src/core/fmath.jl).
# OMIND=1., MIND=5. (so/cwcalc.f DATA MIND/5./,OMIND/1./). Math faithful to gfortran.
# =============================================================================
function so_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32)::Float32
    (1 <= sp <= 33) || return 0f0
    cl  = cr * h * 0.01f0
    ba1 = barea + 1f0

    # Donnelly (R6) / single-power form: a·D^b, OMIND=1 small-D scaling, CW cap.
    _don = function (a::Float32, b::Float32, cap::Float32)
        dm = d >= 1f0 ? d : 1f0
        cw = a * fpow(dm, b)
        d < 1f0 && (cw *= d)
        cw > cap ? cap : cw
    end

    if sp == 1          # WP 11905 R6-m2 BF=1.0 (FIASP 119); no (BA+1) term
        return _cr_r6m2(5.3822f0, 0.57896f0,-0.19579f0,0.14875f0, 0f0,      -0.00685f0, d,h,cl,ba1,el, 10f0,75f0,35f0)
    elseif sp == 2      # SP 11705 R6-m2 BF=1.048 (FIASP 117)
        return _cr_r6m2(3.5930f0*1.048f0, 0.63503f0,-0.22766f0,0.17827f0, 0.04267f0,-0.00290f0, d,h,cl,ba1,el, 5f0,75f0,56f0)
    elseif sp == 3      # DF 20205 R6-m2 BF=1.055 (FIASP 202). Inlined: so/cwcalc.f uses
                        # literal 1.0 (not BF) in the D<OMIND branch — verbatim quirk.
        elc = el < 1f0 ? 1f0 : (el > 75f0 ? 75f0 : el)
        cw = if d >= 1f0
            6.0227f0*1.055f0*fpow(d,0.54361f0)*fpow(h,-0.20669f0)*fpow(cl,0.20395f0)*fpow(ba1,-0.00644f0)*fpow(fexp(elc),-0.00378f0)
        else
            (6.0227f0*1.0f0*fpow(1f0,0.54361f0)*fpow(h,-0.20669f0)*fpow(cl,0.20395f0)*fpow(ba1,-0.00644f0)*fpow(fexp(elc),-0.00378f0))*d
        end
        return cw > 80f0 ? 80f0 : cw
    elseif sp == 4      # WF 01505 R6-m2 BF=1.044 (FIASP 015)
        return _cr_r6m2(5.0312f0*1.044f0, 0.53680f0,-0.18957f0,0.16199f0, 0.04385f0,-0.00651f0, d,h,cl,ba1,el, 2f0,75f0,35f0)
    elseif sp == 5      # MH 26403 Crookston-R1 piecewise (no BF); height-driven small-tree logic
        cw = if h < 5f0
            (0.8f0*h*max(0.5f0,cr*0.01f0))*(1f0-(h-5f0)*0.1f0)*6.90396f0*fpow(d,0.55645f0)*fpow(h,-0.28509f0)*fpow(cl,0.20430f0)*(h-5f0)*0.1f0
        elseif h >= 15f0
            6.90396f0*fpow(d,0.55645f0)*fpow(h,-0.28509f0)*fpow(cl,0.20430f0)
        else
            0.8f0*h*max(0.5f0,cr*0.01f0)
        end
        return cw > 45f0 ? 45f0 : cw
    elseif sp == 6      # IC 08105 R6-m2 BF=0.837 (FIASP 081)
        return _cr_r6m2(5.0446f0*0.837f0, 0.47419f0,-0.13917f0,0.14230f0, 0.04838f0,-0.00616f0, d,h,cl,ba1,el, 5f0,62f0,78f0)
    elseif sp == 7      # LP 10805 R6-m2 BF=1.0 (FIASP 108)
        return _cr_r6m2(6.6941f0, 0.81980f0,-0.36992f0,0.17722f0,-0.01202f0,-0.00882f0, d,h,cl,ba1,el, 1f0,79f0,40f0)
    elseif sp == 8      # ES 09305 R6-m2 BF=1.0 (FIASP 093); no (BA+1) term
        return _cr_r6m2(6.7575f0, 0.55048f0,-0.25204f0,0.19002f0, 0f0,      -0.00313f0, d,h,cl,ba1,el, 1f0,85f0,40f0)
    elseif sp == 9      # SH 02105 R6-m2 BF=1.0 (FIASP 021); no EL term (f=0)
        return _cr_r6m2(2.3170f0, 0.47880f0,-0.06093f0,0.15482f0, 0.05182f0, 0f0,       d,h,cl,ba1,el, 1f0,999f0,65f0)
    elseif sp == 10     # PP 12205 R6-m2 BF=0.918 (FIASP 122)
        return _cr_r6m2(4.7762f0*0.918f0, 0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el, 13f0,75f0,50f0)
    elseif sp == 11     # WJ 06405 R6-m2 BF=1.0 (FIASP 064); no EL term (f=0)
        return _cr_r6m2(5.1486f0, 0.73636f0,-0.46927f0,0.39114f0,-0.05429f0, 0f0,       d,h,cl,ba1,el, 1f0,999f0,36f0)
    elseif sp == 12     # GF 01703 Crookston-R1 log form (no BF)
        dm = d >= 1f0 ? d : 1f0
        v = 1.0303f0 * fexp(1.14079f0 + 0.20904f0*flog(cl) + 0.38787f0*flog(dm))
        d < 1f0 && (v *= d)
        return v > 40f0 ? 40f0 : v
    elseif sp == 13     # AF 01905 R6-m2 BF=0.936 (FIASP 019)
        return _cr_r6m2(5.8827f0*0.936f0, 0.51479f0,-0.21501f0,0.17916f0, 0.03277f0,-0.00828f0, d,h,cl,ba1,el, 10f0,85f0,30f0)
    elseif sp == 14     # SF 01105 R6-m2 BF=1.0 (FIASP 011)
        return _cr_r6m2(4.4799f0, 0.45976f0,-0.10425f0,0.11866f0, 0.06762f0,-0.00715f0, d,h,cl,ba1,el, 4f0,72f0,33f0)
    elseif sp == 15     # NF 02206 Donnelly (no BF; FIASP 022 BF=1.301 inert)
        return _don(3.0614f0, 0.6276f0, 40f0)
    elseif sp == 16     # WB 10105 R6-m2 BF=1.0 (FIASP 101); no (BA+1) and no EL term
        return _cr_r6m2(2.2354f0, 0.66680f0,-0.11658f0,0.16927f0, 0f0,      0f0,        d,h,cl,ba1,el, 1f0,999f0,40f0)
    elseif sp == 17     # WL 07303 Crookston-R1 log form (no BF; FIASP 073 BF=0.818 inert)
        dm = d >= 1f0 ? d : 1f0
        v = 1.02478f0 * fexp(0.99889f0 + 0.19422f0*flog(cl) + 0.59423f0*flog(dm) -
                             0.09078f0*flog(h) - 0.02341f0*flog(barea))
        d < 1f0 && (v *= d)
        return v > 40f0 ? 40f0 : v
    elseif sp == 18     # RC 24205 R6-m2 BF=1.0 (FIASP 242)
        return _cr_r6m2(6.2382f0, 0.29517f0,-0.10673f0,0.23219f0, 0.05341f0,-0.00787f0, d,h,cl,ba1,el, 1f0,72f0,45f0)
    elseif sp == 19     # WH 26305 R6-m2 BF=1.097 (FIASP 263)
        return _cr_r6m2(6.0384f0*1.097f0, 0.51581f0,-0.21349f0,0.17468f0, 0.06143f0,-0.00571f0, d,h,cl,ba1,el, 1f0,72f0,54f0)
    elseif sp == 20     # PY 23104 single-power form 6.1297·D^0.45424 (no BF)
        return _don(6.1297f0, 0.45424f0, 30f0)
    elseif sp == 21     # WA->BM 31206 Donnelly (no BF)
        return _don(7.5183f0, 0.4461f0, 30f0)
    elseif sp == 22     # RA 35106 Donnelly (no BF)
        return _don(7.0806f0, 0.4771f0, 35f0)
    elseif sp == 23     # BM 31206 Donnelly (no BF)
        return _don(7.5183f0, 0.4461f0, 30f0)
    elseif sp == 24     # AS 74605 R6-m2 BF=1.0 (FIASP 746); no (BA+1) and no EL term
        return _cr_r6m2(4.7961f0, 0.64167f0,-0.18695f0,0.18581f0, 0f0,      0f0,        d,h,cl,ba1,el, 1f0,999f0,45f0)
    elseif sp == 25     # CW 74705 R6-m2 BF=1.0 (FIASP 747); no (BA+1) and no EL term
        return _cr_r6m2(4.4327f0, 0.41505f0,-0.23264f0,0.41477f0, 0f0,      0f0,        d,h,cl,ba1,el, 1f0,999f0,56f0)
    elseif sp == 26     # CH->RA 35106 Donnelly (no BF)
        return _don(7.0806f0, 0.4771f0, 35f0)
    elseif sp == 27     # WO 81505 R6-m2 BF=1.0 (FIASP 815); no H, no (BA+1), no EL term
        return _cr_r6m2(2.4857f0, 0.70862f0, 0f0,      0.10168f0, 0f0,      0f0,        d,h,cl,ba1,el, 1f0,999f0,39f0)
    elseif sp == 28     # WI->BM 31206 Donnelly (no BF)
        return _don(7.5183f0, 0.4461f0, 30f0)
    elseif sp == 29     # GC->TO 63102 Bechtold-2004 model-2 (no BF); HI in [-55,15], MIND=5
        hv = hi < -55f0 ? -55f0 : (hi > 15f0 ? 15f0 : hi)
        cw = if d >= 5f0
            3.1150f0 + 0.7966f0*d + 0.0745f0*cr + (-0.0053f0)*barea + 0.0523f0*hv
        else
            (3.1150f0 + 0.7966f0*5f0 + 0.0745f0*cr + (-0.0053f0)*barea + 0.0523f0*hv)*(d/5f0)
        end
        return cw > 41f0 ? 41f0 : cw
    elseif sp == 30     # MC 47502 Bechtold-2004 model-2 (no BF); HI in [-37,27], MIND=5
        hv = hi < -37f0 ? -37f0 : (hi > 27f0 ? 27f0 : hi)
        cw = if d >= 5f0
            4.0105f0 + 0.8611f0*d + (-0.0431f0)*hv
        else
            (4.0105f0 + 0.8611f0*5f0 + (-0.0431f0)*hv)*(d/5f0)
        end
        return cw > 29f0 ? 29f0 : cw
    elseif sp == 31     # MB->MC 47502 Bechtold-2004 model-2 (no BF); HI in [-37,27], MIND=5
        hv = hi < -37f0 ? -37f0 : (hi > 27f0 ? 27f0 : hi)
        cw = if d >= 5f0
            4.0105f0 + 0.8611f0*d + (-0.0431f0)*hv
        else
            (4.0105f0 + 0.8611f0*5f0 + (-0.0431f0)*hv)*(d/5f0)
        end
        return cw > 29f0 ? 29f0 : cw
    elseif sp == 32     # OS 12205 R6-m2 BF=0.918 (FIASP 122) — same code as PP
        return _cr_r6m2(4.7762f0*0.918f0, 0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el, 13f0,75f0,50f0)
    else                # sp == 33  OH 31206 Donnelly (no BF)
        return _don(7.5183f0, 0.4461f0, 30f0)
    end
end
