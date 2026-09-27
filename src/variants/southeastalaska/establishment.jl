# =============================================================================
# establishment.jl (southeastalaska) — the AK regeneration establishment model.
#
# Ports what FVSak links: estb/esnutr.f (the ESTAB-call scheduler: user TALLY/TALLYONE/TALLYTWO, the LAUTAL
# automatic tally after a removal, the ≤19-yr continuation, LINGRW ingrowth, the PLANT/NATURAL catch-all) and
# ak/estab.f (2020 refit, forest-type based) with its ak/ helpers: estock, estpp, esnspe, espadv/espsub/espxcs
# (identical logistic occurrence models), esdlay, esadvh, essubh, esxcsh, estime, esprep, and ak/esgent.f
# (regent.jl ak_esgent!). Natural regeneration is predicted per replicated regen plot (MINREP=50 plots spread
# over the stockable inventory points) from the forest-type group IFT (1-14) of the current (else inventory)
# forest type; PLANT/NATURAL keyword trees are added per plot at ESSUBH (HTCALC site-curve) heights.
#
# Every ESRANN draw of estab.f is reproduced in order (the :estab stream); the main RANN stream is used only by
# ESGENT's crown/ZZRAN draws.
# =============================================================================

const AK_ES_MAXTPP = Int[17,30,10,22,21,21,21,10,10,15,15,17,5,0]
const AK_ES_MAXSPP = Int[4,4,4,5,6,5,6,4,3,4,4,4,3,0]
const AK_ES_MAXING = Int[3,3,3,3,3,3,3,3,3,3,3,3,3,0]
const AK_ES_MYHABG = Int[1,1,1,1,2,2,2,2,3,4,5,5,5,5,5,5]
const AK_ES_MAXPLT = 500
const AK_ES_REGNBK = 1.0f0            # ak/blkdat.f PDEN REGNBK
const AK_ES_IDCMP1 = 10000000

# ak/blkdat.f OCURFT(23 species, 14 forest-type groups) — species occurrence by forest type (0/1).
const AK_ES_OCURFT = let m = zeros(Float32, 23, 14)
    rows = (
        (0,0,0,1,1,0,1,0,0,0,0,1,0, 0,0,1,0,1,1,0,0,0,0),   # 122 white spruce
        (0,0,0,1,1,0,1,0,0,0,0,0,0, 0,0,1,0,1,1,0,0,0,0),   # 125 black spruce
        (0,0,1,0,1,0,0,1,1,0,1,1,0, 0,0,0,0,0,0,0,0,0,0),   # 270 mountain hemlock
        (0,0,1,0,0,0,0,1,1,1,1,1,0, 0,1,0,0,0,0,0,0,0,0),   # 271 Alaska cedar
        (0,0,1,0,0,0,0,1,1,1,1,1,0, 0,1,0,0,0,0,0,0,0,0),   # 281 lodgepole pine
        (0,0,1,0,0,0,0,1,1,1,1,1,0, 0,1,0,0,0,0,0,0,0,0),   # 301 western hemlock
        (0,0,1,0,0,0,0,1,1,1,1,1,0, 0,1,0,0,0,0,0,0,0,0),   # 304 western redcedar
        (0,0,1,0,0,0,0,1,0,1,1,1,0, 0,1,0,0,0,0,0,0,0,0),   # 305 Sitka spruce
        (0,0,0,0,1,0,1,1,0,0,1,0,0, 0,0,1,0,0,0,1,0,0,0),   # 703 cottonwood
        (0,0,0,1,1,0,1,0,0,0,0,0,0, 0,0,1,0,1,1,0,0,0,0),   # 901 aspen
        (0,0,0,1,1,0,1,1,0,0,0,1,0, 0,0,1,0,1,1,0,0,0,0),   # 902 paper birch
        (0,0,0,1,1,0,1,0,0,0,0,1,0, 0,0,1,0,1,1,0,0,0,0),   # 904 balsam poplar
        (0,0,0,0,0,0,0,1,0,0,1,0,0, 0,1,0,0,0,0,0,0,0,0),   # 911 red alder
        (0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0))   # other
    for (j, r) in enumerate(rows), i in 1:23
        m[i, j] = Float32(r[i])
    end
    m
end

# ak/espadv.f (= espsub.f = espxcs.f) coefficients: PNFT(23,14) + per-species PRDA/BA-share/ELEV/SLOPE/LAT terms.
const AK_ES_PNFT = let m = zeros(Float32, 23, 14)
    cols = (
        (0.0,0.0,0.0,-4.335977,-0.768725, 0.0,-13.023483,0.0,0.0,0.0, 0.0,-14.636196,0.0,0.0,0.0,
         4.543671,0.0,-0.129163,-10.517939,0.0, 0.0,0.0,0.0),
        (0.0,0.0,0.0,-4.335977,-1.793084, 0.0,-10.509385,0.0,0.0,0.0, 0.0,-15.787063,0.0,0.0,0.0,
         3.846221,0.0,-4.067106,-11.048545,0.0, 0.0,0.0,0.0),
        (0.0,0.0,17.052464,0.0,-3.104481, 0.0,0.0,-1.727407,13.303874,77.557037, 23.185698,-10.383977,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0),
        (0.0,0.0,17.998179,0.0,0.0, 0.0,0.0,-1.900489,14.426165,79.436554, 23.405695,-9.457501,0.0,0.0,35.036460,
         0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0),
        (0.0,0.0,18.116574,0.0,0.0, 0.0,0.0,-3.220587,14.50099,79.452263, 22.762266,-10.122027,0.0,0.0,35.036460,
         0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0),
        (0.0,0.0,15.963877,0.0,0.0, 0.0,0.0,-1.224259,12.070894,78.000985, 24.451408,-11.165193,0.0,0.0,35.036460,
         0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0),
        (0.0,0.0,17.094312,0.0,0.0, 0.0,0.0,-1.434335,13.561585,79.462069, 23.941341,-10.180572,0.0,0.0,35.036460,
         0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0),
        (0.0,0.0,14.997227,0.0,0.0, 0.0,0.0,-0.563591,11.217778,78.007480, 22.942334,-11.404114,0.0,0.0,35.036460,
         -0.265927,0.0,0.0,0.0,0.0, 0.0,0.0,0.0),
        (0.0,0.0,0.0,0.0,0.231061, 0.0,-14.506330,-2.473404,0.0,0.0, 21.694368,-14.517376,0.0,0.0,0.0,
         3.879223,0.0,0.0,-10.664679,-2.676858, 0.0,0.0,0.0),
        (0.0,0.0,0.0,-4.335977,0.65022, 0.0,-12.961809,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,
         3.691336,0.0,0.107820,-11.026581,0.0, 0.0,0.0,0.0),
        (0.0,0.0,0.0,-4.335977,-0.090381, 0.0,-12.383669,-5.003714,0.0,0.0, 20.687934,-15.353202,0.0,0.0,0.0,
         5.542455,0.0,-0.743739,-10.160502,0.0, 0.0,0.0,0.0),
        (0.0,0.0,0.0,-4.335977,-0.768725, 0.0,-13.023483,0.0,0.0,0.0, 0.0,-14.636196,0.0,0.0,0.0,
         4.543671,0.0,-0.129163,-10.517939,0.0, 0.0,0.0,0.0),
        (0.0,0.0,0.0,0.0,0.0, 0.0,0.0,-0.042003,0.0,0.0, 21.940216,0.0,0.0,0.0,35.03646,
         0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0),
        ntuple(_ -> 0.0, 23))
    for (j, c) in enumerate(cols), i in 1:23
        m[i, j] = Float32(c[i])
    end
    m
end
const AK_ES_PNRDA  = Float32[0.0,0.0,-2.727570,0.0,0.774395, 0.0,1.660146,-0.838168,-3.839624,-1.443889,
                             1.739987,-0.620028,0.0,0.0,0.851525, -1.846654,0.0,-0.464799,-0.397800,-0.464799, 0.0,0.0,0.0]
const AK_ES_PBAPER = Float32[0.0,0.0,3.266453,0.0,2.917051, 0.0,19.545809,0.891633,2.472698,1.074455,
                             5.056352,3.062509,0.0,0.0,5.115524, 0.957204,0.0,5.337844,3.757893,5.337844, 0.0,0.0,0.0]
const AK_ES_PNELEV = Float32[0.0,0.0,0.000346,0.0,0.000580, 0.0,0.000381,-0.000293,-0.000699,-0.001340,
                             -0.000998,0.000542,0.0,0.0,-0.004470, -0.000526,0.0,-0.000588,-0.000492,-0.000588, 0.0,0.0,0.0]
const AK_ES_PNSLO  = Float32[0.0,0.0,-0.003355,0.0,-0.020237, 0.0,-0.044745,0.003523,-0.002279,0.006463,
                             0.007062,-0.007500,0.0,0.0,-0.027259, 0.021042,0.0,0.019300,0.029324,0.019300, 0.0,0.0,0.0]
const AK_ES_PNTLAT = Float32[0.0,0.0,-0.314786,0.0,-0.025084, 0.0,0.156958,0.023200,-0.253624,-1.427424,
                             -0.415727,0.162942,0.0,0.0,-0.707129, -0.077977,0.0,-0.035312,0.145927,-0.035312, 0.0,0.0,0.0]
# ak/estock.f
const AK_ES_STK_PNFT   = Float32[1.428325,3.401549,1.906730,3.481250,3.172262,2.816007,4.073819,1.255094,0.357560,
                                 2.756703,1.471673,1.428325,1.554008,1.428325]
const AK_ES_STK_PNELEV = Float32[0.022431,-0.045237,0.036110,0.006615,0.017794,0.072023,0.012058,0.066790,0.135631,
                                 0.012600,0.032712,0.022431,0.037521,0.022431]
const AK_ES_STK_PNXCOS = Float32[0.006748,0.000011,0.000738,-0.000432,-0.000715,-0.002044,0.000694,-0.003398,-0.013222,
                                 0.003788,-0.009518,0.006748,-0.001041,0.006748]
# ak/estpp.f
const AK_ES_TPP_BB = Float32[5.617512,13.273635,4.291578,9.770818,9.452008,8.075209,9.635466,3.632153,3.794353,
                             6.107958,5.251271,5.617512,2.695103,5.617512]
const AK_ES_TPP_CC = Float32[1.213626,1.405827,1.442975,1.496351,1.523877,1.274965,1.505655,1.414683,1.462014,
                             1.399479,1.268050,1.213626,1.709675,1.213626]
# ak/esnspe.f NSSPB1(6,14) / NSSPB2(6)
const AK_ES_NSSPB1 = Float32[
     1.514466 -1.130485 -4.345526 -6.061343 -9.282172  -9.263319;
     2.33692  -1.364289 -4.725742 -7.190956 -10.817028 -9.695317;
     0.787703 -0.822546 -2.824892 -5.172281 -8.815643  -9.039528;
    -0.227513 -0.743091 -1.941313 -3.391035 -6.505262  -9.421316;
    -0.041975 -0.643339 -1.961849 -4.009763 -8.205985  -5.996876;
     0.837373 -0.687024 -3.137052 -5.605608 -9.147572  -9.551078;
     0.295898 -0.757339 -2.411081 -3.820975 -6.728512  -7.271471;
     1.261125 -1.085806 -3.752193 -5.807558 -8.732909  -8.933396;
     1.535307 -1.288482 -3.233754 -5.765195 -8.679735  -8.831665;
     0.827381 -0.632341 -3.309655 -6.342733 -9.042995  -9.132168;
     0.955467 -0.845367 -3.046738 -6.182925 -9.180378  -9.272577;
     1.514466 -1.130485 -4.345526 -6.061343 -9.282172  -9.263319;
     0.856235 -0.937544 -2.764232 -5.426666 -8.64982   -8.767986;
     1.514466 -1.130485 -4.345526 -6.061343 -9.282172  -9.263319]   # [IFT, k]
const AK_ES_NSSPB2 = Float32[-0.095656, 0.027575, 0.082907, 0.111271, 0.196761, 0.09174]
# ak/esdlay.f
const AK_ES_BSUB = Float32[4.34376 6.55916 9.16226; 3.45725 6.34975 8.65545; 4.16284 7.52536 10.20937;
    3.52946 7.62339 12.79801; 3.52946 7.62339 12.79801; 3.52946 7.62339 12.79801; 3.52946 7.62339 12.79801;
    5.36466 7.44468 9.69507; 5.33757 6.78727 9.45827; 5.23792 7.38005 10.42350; 4.33094 6.30802 8.63060;
    4.17909 5.88262 8.49857; 4.33094 6.30802 8.63060;
    3.81610 5.74622 9.36345; 3.81610 5.74622 9.36345; 3.81610 5.74622 9.36345; 3.81610 5.74622 9.36345;
    3.81610 5.74622 9.36345; 3.81610 5.74622 9.36345; 3.81610 5.74622 9.36345; 3.81610 5.74622 9.36345;
    3.81610 5.74622 9.36345; 3.81610 5.74622 9.36345]                                   # [sp, IT]
const AK_ES_CSUB = Float32[2.33194 2.63560 2.21663; 2.06804 2.79933 2.28687; 2.03892 3.12279 2.68340;
    1.71621 2.72466 3.98359; 1.71621 2.72466 3.98359; 1.71621 2.72466 3.98359; 1.71621 2.72466 3.98359;
    2.89777 2.80504 3.27745; 4.16994 3.59937 2.39138; 3.11598 3.04038 2.87196; 1.97408 2.35053 2.00997;
    2.47058 2.25957 2.00065; 1.97408 2.35053 2.00997;
    3.01975 2.09376 1.80925; 3.01975 2.09376 1.80925; 3.01975 2.09376 1.80925; 3.01975 2.09376 1.80925;
    3.01975 2.09376 1.80925; 3.01975 2.09376 1.80925; 3.01975 2.09376 1.80925; 3.01975 2.09376 1.80925;
    3.01975 2.09376 1.80925; 3.01975 2.09376 1.80925]
const AK_ES_BADV = Float32[13.121021 22.186540; 21.962337 32.176727; 17.779381 24.241344;
    6.699826 13.431179; 6.699826 13.431179; 6.699826 13.431179; 6.699826 13.431179;
    13.990273 21.779362; 7.358880 32.955809; 9.768223 27.100242; 13.604594 19.605485;
    11.269182 18.245664; 13.604594 19.605485;
    8.986115 10.312660; 8.986115 10.312660; 8.986115 10.312660; 8.986115 10.312660; 8.986115 10.312660;
    8.986115 10.312660; 8.986115 10.312660; 8.986115 10.312660; 8.986115 10.312660; 8.986115 10.312660]  # [sp, IBAA]
const AK_ES_CADV = Float32[1.043254 1.215651; 1.101266 1.317022; 1.337217 1.540663;
    1.262533 1.279302; 1.262533 1.279302; 1.262533 1.279302; 1.262533 1.279302;
    1.051296 1.416222; 0.912295 1.230540; 1.152577 1.319692; 1.267057 1.287445;
    1.122139 1.082805; 1.267057 1.287445;
    1.068472 1.470405; 1.068472 1.470405; 1.068472 1.470405; 1.068472 1.470405; 1.068472 1.470405;
    1.068472 1.470405; 1.068472 1.470405; 1.068472 1.470405; 1.068472 1.470405; 1.068472 1.470405]

"estab.f IESFT→IFT: the forest-type group (1-13 named types, 14 = other)."
function ak_es_ift(iesft::Integer)::Int
    iesft == 122 && return 1
    iesft == 125 && return 2
    (iesft == 270 || iesft == 268) && return 3
    iesft == 271 && return 4
    iesft == 281 && return 5
    (iesft == 301 || iesft == 264) && return 6
    iesft == 304 && return 7
    iesft == 305 && return 8
    (iesft == 703 || iesft == 704 || iesft == 709) && return 9
    iesft == 901 && return 10
    iesft == 902 && return 11
    iesft == 904 && return 12
    iesft == 911 && return 13
    return 14
end
const _AK_ES_NAMED_FT = (122, 125, 264, 268, 270, 271, 281, 301, 304, 305, 703, 704, 709, 901, 902, 904, 911)

"SAVEd AK ESTAB state (commons + SAVEd locals, -fno-automatic)."
mutable struct AKEstabState
    ift0::Int
    esdraw::Float32
    kdtold::Int
    esb::Float32
    esa::Float32
    zmech::Float32; zburn::Float32
    pnone::Float32; pmech::Float32; pburn::Float32
    meth::Int; ialn::Vector{Int}
    load::Int
    esb1::Vector{Float32}; pnn::Vector{Float32}; plprob::Vector{Float32}; nstore::Vector{Int}
    prob1::Vector{Float32}; xstore::Vector{Float32}; ipprep::Vector{Int}
    delay::Float32            # SAVEd local DELAY (read stale by the excess-tree AGEXC)
    time::Float32             # ESCOMN TIME
    xtes::Float32             # ESHAP XTES — the removal ratio of the last LAUTAL evaluation
    xmaxpt::Vector{Float32}   # XMAXPT of the last SDICAL (morts.f:202) — the DENSE PRDA divisor
    xmax::Float32
end
AKEstabState() = AKEstabState(0, 0f0, -99, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0, zeros(Int, 3), 0,
                              zeros(Float32, AK_ES_MAXPLT), zeros(Float32, AK_ES_MAXPLT), zeros(Float32, AK_ES_MAXPLT),
                              zeros(Int, AK_ES_MAXPLT), zeros(Float32, AK_ES_MAXPLT), zeros(Float32, AK_ES_MAXPLT),
                              ones(Int, AK_ES_MAXPLT + 1), 0f0, 0f0, 0f0, Float32[], 1f0)

@inline _ak_es_state(s::StandState) = (s.estab.ak_state === nothing && (s.estab.ak_state = AKEstabState());
                                       s.estab.ak_state::AKEstabState)

"ESFLTR (fvs.f:201) for AK: the inventory per-point overstory BAAINV/TPAAINV (D≥REGNBK, ×PIX=PI−NONSTK)."
function ak_esfltr!(s::StandState)
    s.variant isa SoutheastAlaska || return s
    p, t = s.plot, s.trees
    npt = max(1, Int(p.points_inv))
    ba = zeros(Float32, max(npt, AK_ES_MAXPLT)); tp = zeros(Float32, max(npt, AK_ES_MAXPLT))
    pix = Float32(p.pi) - Float32(p.nonstockable)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d < AK_ES_REGNBK && continue
        n = Int(t.plot_id[i]); (1 <= n <= length(ba)) || continue
        ba[n] += 0.005454154f0 * d * d * t.tpa[i] * pix
        tp[n] += t.tpa[i] * pix
    end
    s.estab.inv_point_baaold = ba
    s.estab.esb_shift_pt = tp          # AK: TPAAINV (the IE ESB-shift vector is unused by AK)
    return s
end

# ak/estock.f — stocking logit.
@inline function ak_estock(elev::Float32, ift::Int, xba::Float32, xtpa::Float32, slo::Float32, xcos::Float32,
                           elevsq::Float32)::Float32
    xqmd = xtpa > 0f0 ? sqrt((xba / xtpa) / 0.005454f0) : 0f0
    adjslo = slo * 100f0; adjxc = xcos * 100f0
    return AK_ES_STK_PNFT[ift] + (-0.00212f0) * xba + (-0.194748f0) * xqmd + (-0.002364f0) * adjslo +
           AK_ES_STK_PNELEV[ift] * elev + (-0.000742f0) * elevsq + AK_ES_STK_PNXCOS[ift] * adjxc +
           0.000316f0 * xba * xqmd
end

# ak/espadv.f (= espsub.f = espxcs.f): per-species occurrence probability on plot NNID.
function ak_espadv!(out::AbstractVector{Float32}, ift::Int, prda::Float32, over::AbstractVector{Float32},
                    baaa::Float32, slo::Float32, elev::Float32, tlat::Float32, xesmlt)
    if ift == 14
        fill!(out, 0f0); return out
    end
    denom = baaa <= 0f0 ? 1f0 : baaa
    adjslo = slo * 100f0; adjelv = elev * 100f0
    @inbounds for j in 1:23
        if AK_ES_OCURFT[j, ift] == 0f0
            out[j] = 0f0
        else
            baper = over[j] / denom
            pn = AK_ES_PNFT[j, ift] + AK_ES_PNRDA[j] * prda + AK_ES_PBAPER[j] * baper + AK_ES_PNELEV[j] * adjelv +
                 AK_ES_PNSLO[j] * adjslo + AK_ES_PNTLAT[j] * tlat
            v = (fexp(pn) / (1f0 + fexp(pn))) * AK_ES_OCURFT[j, ift] * xesmlt[j]
            v < 0f0 && (v = 0f0); v > 1f0 && (v = 1f0)
            out[j] = v
        end
    end
    return out
end

# ak/esdlay.f
function ak_esdlay(isp::Int, ias::Int, draw::Float32, time::Float32, baa::Float32)::Float32
    it = 1
    (time > 7.5f0 && time < 12.5f0) && (it = 2)
    time > 12.5f0 && (it = 3)
    ibaa = baa > 25.5f0 ? 2 : 1
    local delay::Float32
    if ias == 1
        bb = AK_ES_BADV[isp, ibaa]; cc = AK_ES_CADV[isp, ibaa]
        delay = fpow(-(flog(1.0f0 - draw)), 1.0f0 / cc) * bb
        delay = (delay + 3.0f0) * (-1.0f0)
    else
        bb = AK_ES_BSUB[isp, it]; cc = AK_ES_CSUB[isp, it]
        delay = fpow(-(flog(1.0f0 - draw)), 1.0f0 / cc) * bb
        delay = delay - 4.0f0
    end
    if delay >= 10f0
        delay = 10f0
    elseif delay <= 0f0
        delay = 0f0
    end
    return delay
end

"ak/essubh.f — returns (HHT, DELAY', TRAGE') (DELAY/TRAGE are modified in place in Fortran)."
function ak_essubh(sp::Int, delay::Float32, gentim::Float32, trage::Float32, time::Float32, xsite::Float32)
    n = unsafe_trunc(Int, delay + 0.5f0)
    n < -3 && (n = -3)
    delay = Float32(n)
    itime = unsafe_trunc(Int, time + 0.5f0)
    n > itime && (delay = time)
    age = time - delay - gentim
    age = age + trage
    age < 1.0f0 && (age = 1.0f0)
    trage = time - delay
    hht = ak_htcalc_htmax(sp, xsite) - 0f0 <= 1f0 ? 0f0 : ak_htcalc_height(sp, xsite, age)
    hht < AK_ES_XMIN[sp] && (hht = AK_ES_XMIN[sp])
    return hht, delay, trage
end

"ak/esadvh.f — advance-regen height; returns (HHT, DELAY', TRAGE')."
function ak_esadvh(sp::Int, delay::Float32, gentim::Float32, xsite::Float32)
    n = unsafe_trunc(Int, delay + 0.5f0)
    n > 2 && (n = 1)
    delay = Float32(n)
    trage = 3.0f0 - delay
    age = 3.0f0 - delay - gentim
    age < 1.0f0 && (age = 1.0f0)
    hht = ak_htcalc_htmax(sp, xsite) <= 1f0 ? 0f0 : ak_htcalc_height(sp, xsite, age)
    hht < AK_ES_XMIN[sp] && (hht = AK_ES_XMIN[sp])
    return hht, delay, trage
end

"ak/esxcsh.f — excess-tree height = the site-curve height at age 1."
function ak_esxcsh(sp::Int, xsite::Float32)::Float32
    hht = ak_htcalc_htmax(sp, xsite) <= 1f0 ? 0f0 : ak_htcalc_height(sp, xsite, 1.0f0)
    hht < AK_ES_XMIN[sp] && (hht = AK_ES_XMIN[sp])
    return hht
end


_ak_act_date(s::StandState, a) = (0 < Int(a.year) < 1000) ? Int(cycle_year_at(s.control, Int(a.year) - 1)) : Int(a.year)
# OPFIND "this cycle": dates in [IY(ICYC), IY(ICYC+1)); cycle 1 also collects any earlier (pre-inventory) date.
@inline _ak_in_cycle(d::Integer, icyc::Integer, year::Integer, next_year::Integer) =
    (icyc == 1 ? d < next_year : year <= d < next_year)

"morts.f:202 SDICAL(0,SDIMAX) — stash XMAXPT/XMAX; the next DENSE (gradd.f:192) divides PRDA by them."
function ak_es_stash_xmaxpt!(s::StandState, xmaxpt::Vector{Float32}, xmax::Float32)
    st = _ak_es_state(s)
    st.xmaxpt = copy(xmaxpt); st.xmax = xmax
    return s
end

"""
    ak_esnutr!(s; fint) -> Int

estb/esnutr.f for AK (the ADDTREES/sprout steps run elsewhere): decide whether and how ESTAB is called this cycle
(user TALLYONE/TALLYTWO, TALLY, the LAUTAL automatic tally, the ≤19-yr continuation, LINGRW ingrowth, the
PLANT/NATURAL catch-all), run ak/estab.f (+ ESGENT), and reset NTALLY after a LONE call. Returns the number of
records added.
"""
function ak_esnutr!(s::StandState; fint::Float32 = 10.0f0)::Int
    s.variant isa SoutheastAlaska || return 0
    est = s.estab; ctl = s.control; st = _ak_es_state(s)
    per = round(Int, fint)
    icyc = Int(ctl.cycle) + 1
    year = Int(current_cycle_year(s))
    next_year = year + per
    inv_year = Int(ctl.cycle_year[1])
    incyc(a) = _ak_in_cycle(_ak_act_date(s, a), icyc, year, next_year)
    # STOCKADJ (440) due this cycle (esnutr.f:63-110; SPECMULT/HTADJ are applied at keyword time in jl).
    for a in ctl.schedule
        (a.icflag == Int32(440) && incyc(a)) && (est.stoadj = a.params[1])
    end
    est.idsdat == Int32(-9999) && (est.idsdat = Int32(inv_year - 20))
    kdt = next_year - 1
    lone = false; fire = false; cancel = false
    # esnutr.f:163-212 TALLYONE (428), else TALLYTWO (429), due this cycle (the last one found).
    for code in (Int32(428), Int32(429))
        fire && break
        due = [a for a in ctl.schedule if a.icflag == code && incyc(a)]
        isempty(due) && continue
        a = due[end]
        est.ntally = code - Int32(427)
        idsdat = round(Int, a.params[1]); (1 <= idsdat < 1000) && (idsdat = Int(cycle_year_at(ctl, idsdat - 1)))
        est.idsdat = Int32(idsdat)
        lone = true
        if kdt + 1 - idsdat > 20
            cancel = true
        else
            if est.ntally == 2 && !any(b -> b.icflag == Int32(428) && idsdat <= _ak_act_date(s, b) <= kdt, ctl.schedule)
                est.ntally = Int32(1)
            end
            fire = true
        end
        break
    end
    # esnutr.f:218-252 TALLY (427): the latest disturbance date.
    if !fire && !cancel
        due = [a for a in ctl.schedule if a.icflag == Int32(427) && incyc(a)]
        if !isempty(due)
            idsdat = -1
            for a in due
                i = round(Int, a.params[1]); i > idsdat && (idsdat = i)
            end
            (1 <= idsdat < 1000) && (idsdat = Int(cycle_year_at(ctl, idsdat - 1)))
            est.idsdat = Int32(idsdat)
            if kdt + 1 - idsdat > 20
                cancel = true
            else
                est.ntally = Int32(1); fire = true
            end
        end
    end
    if !fire && !cancel
        # esnutr.f:264-289 LAUTAL automatic tally after a removal.
        if est.lautal
            st.xtes = est.last_xtes
            lone = (st.xtes >= est.thres1 && st.xtes < est.thres2)
            if lone || st.xtes >= est.thres2
                est.idsdat = Int32(year); est.ntally = Int32(1); fire = true
            end
        end
    end
    if !fire && !cancel && kdt - Int(est.idsdat) <= 19 && est.ntally > 0     # esnutr.f:298 continuation
        est.ntally += Int32(1); fire = true
    end
    if !fire && !cancel && est.lingrw                                          # esnutr.f:313 ingrowth
        keep = per <= 20 && any(a -> (a.icflag == Int32(427) || a.icflag == Int32(428) || a.icflag == Int32(429)) &&
                                     next_year <= _ak_act_date(s, a) < next_year + per, ctl.schedule)
        if keep
            cancel = true            # GOTO 400 (no ESTAB, no deletions)
        elseif (s.trees.n == 0 && icyc == 1) || (next_year - Int(est.idsdat) >= 40)
            est.ntally = Int32(99); est.idsdat = Int32(next_year - 20); lone = true; fire = true
        end
    end
    if !fire && !cancel && any(a -> (a.icflag == Int32(430) || a.icflag == Int32(431)) && incyc(a), ctl.schedule)
        est.ntally = Int32(99); est.idsdat = Int32(next_year - 20); fire = true   # esnutr.f:345 catch-all
    end
    est.last_xtes = 0f0
    fire || return 0
    nadd = ak_estab!(s, kdt; fint = fint)
    lone && (est.ntally = Int32(0))
    return nadd
end

"""
    ak_estab!(s, kdt; fint) -> Int

ak/estab.f (called from ESNUTR with the tally-year KDT): the regeneration tally over NPTIDS·IDUP replicated regen
plots — natural stocking (ESTOCK), trees per plot (ESTPP), number of species (ESNSPE), species (ESPADV/ESPSUB),
advance/subsequent heights (ESDLAY/ESADVH/ESSUBH), excess trees (ESPXCS/ESXCSH), the PLANT/NATURAL keyword trees,
the best-tree choice — then the new records, ESGENT and ABIRTH+=GENTIM. Returns the number of records added.
"""
function ak_estab!(s::StandState, kdt::Integer; fint::Float32 = 10.0f0)::Int
    p, t, est, ctl = s.plot, s.trees, s.estab, s.control
    st = _ak_es_state(s)
    rng = s.rng
    nofspe = 23
    icyc = Int(ctl.cycle) + 1
    year = Int(current_cycle_year(s)); next_year = year + round(Int, fint)
    # STEP 1
    itrnin = t.n + 1
    ingro = 0
    zharv = Float32(est.idsdat)
    flokdt = Float32(kdt)
    st.time = flokdt - zharv + 1.0f0
    nptids = max(1, Int(p.points_inv) - Int(p.nonstockable))
    flonpt = Float32(nptids)
    elev = p.elevation; elevsq = elev * elev
    itmp = icyc * 10000
    mode = 0
    if kdt + 1 - Int(ctl.cycle_year[1]) > 20
        est.inadv = true; st.load = 0
    end
    pasmax = est.pasmax
    note = zeros(Int, 99); agadsb = zeros(Float32, 99); agepl = zeros(Float32, 99); agexc = zeros(Float32, 99)
    tottpa = zeros(Float32, nofspe); aveht = zeros(Float32, nofspe); sumpx = zeros(Float32, nofspe)
    sumpi = zeros(Float32, nofspe); bestpa = zeros(Float32, nofspe); pastpa = zeros(Float32, nofspe)
    iplant = zeros(Int, nofspe)
    tcrop1 = 0f0; tcrop2 = 0f0; ttottp = 0f0
    ifortp = Int(p.forest_type); iifortp = Int(p.inv_forest_type)
    iesft = ifortp in _AK_ES_NAMED_FT ? ifortp : iifortp
    ift = ak_es_ift(iesft)
    est.ntally == 1 && (st.ift0 = ift)
    est.ntally < 99 && (ift = st.ift0)
    minrep = Int(est.minrep)
    dup = 0f0
    for i in 1:minrep
        dup += 1f0
        n = nptids * i
        nptids * (i + 1) > AK_ES_MAXPLT && break
        n >= minrep && break
    end
    idup = unsafe_trunc(Int, dup + 0.5f0)
    dupnpt = flonpt * dup
    expand = ift in (1, 2, 10, 11, 12) ? 250f0 : 300f0
    # STEP 2
    if est.ntally > 1
        estab_prep_cancel_cycle!(s)                      # OPFIND(2,MYACTS(3)) + OPDEL1 — this cycle's BURN/MECH preps
    end
    match = (est.ntally == 1 && kdt + 1 - Int(ctl.cycle_year[1]) > 20) ? 1 : 0
    shorty = 0f0
    if est.ntally == 99
        pnone = 1.0f0; st.pmech = 0f0; st.pburn = 0f0; st.meth = 0; st.load = 0; st.ialn[2] = 1
        ingro = 1
        pasmax = 15.0f0
        itemp = unsafe_trunc(Int, ((st.xtes / est.thres1) * fint) + 0.5f0)
        itemp < 1 && (itemp = 1)
        shorty = Float32(itemp)
        zharv = Float32(kdt + 1) - shorty
        st.zmech = zharv; st.zburn = zharv
        st.pnone = pnone
        est.ntally = Int32(1)
        @goto l44
    end
    est.ntally != 1 && @goto l276
    # ESETPR (esetpr.f): the last MECHPREP/BURNPREP in [IDSDAT,KDT]; a percent-less one sets METH and stops.
    let
        st.meth = 0; st.zmech = 0f0; st.zburn = 0f0; st.pnone = 0f0; st.pmech = 0f0; st.pburn = 0f0; fill!(st.ialn, 0)
        for (code, slot) in ((Int32(491), 3), (Int32(493), 2))
            for a in ctl.schedule
                a.icflag == code || continue
                d = _ak_act_date(s, a)
                (Int(est.idsdat) <= d <= kdt) || continue
                if slot == 3; st.zburn = Float32(d); else; st.zmech = Float32(d); end
                st.load = 0
                pct = a.params[1]
                if pct <= 0f0
                    st.meth = slot; break
                end
                if slot == 3; st.pburn = pct / 100f0; else; st.pmech = pct / 100f0; end
                st.ialn[slot] = 1
            end
            slot == 3 && st.meth == 3 && break
        end
        esetpr_mark!(s, Int(est.idsdat), Int(kdt))
    end
    st.zmech < zharv && (st.zmech = zharv)
    st.zburn < zharv && (st.zburn = zharv)
    @label l44
    st.kdtold = unsafe_trunc(Int, zharv - 0.5f0)
    sumpre = zeros(Float32, 5)
    let ndp = nptids * idup
        @inbounds for i in 1:ndp
            st.esb1[i] = 0f0; st.pnn[i] = 0f0; st.plprob[i] = 0f0; st.nstore[i] = 0
        end
    end
    st.esb = 0f0
    @label l276
    if est.ntally == 1
        draw = esrann!(rng)
        st.esdraw = Float32(unsafe_trunc(Int32, draw * 100000.0f0 + 0.5f0))
    end
    (st.esdraw % 2f0 == 0f0) && (st.esdraw += 1f0)
    rng.es0 = Float64(st.esdraw)                          # ESRNSD(.TRUE.,ESDRAW)
    # LOAD REGENERATION VECTORS FROM CURRENT INVENTORY (D<REGNBK)
    sum0 = 0f0
    @inbounds for i in 1:t.n
        ftemp = t.dbh[i]
        ftemp >= AK_ES_REGNBK && continue
        nn = Int(t.species[i]); (1 <= nn <= nofspe) || continue
        ftemp2 = t.tpa[i]
        pastpa[nn] += ftemp2; tcrop2 += ftemp2; sum0 += ftemp2
        n = Int(t.plot_id[i])
        (1 <= n <= AK_ES_MAXPLT) && (st.plprob[n] += ftemp2 / dup)
    end
    if !(est.inadv || est.ntally != 1)
        tpacre = sum0 < 1f0 ? 1f0 : sum0
        st.esa = 1.0f0 / (1.0f0 + fexp(3.190754f0 - 0.59894f0 * flog(tpacre)))
        ftemp = st.esa
        ftemp < 0.10f0 && (ftemp = 0.10f0); ftemp > 0.90f0 && (ftemp = 0.90f0)
        st.esb = -1.0f0 * flog(1.0f0 / ftemp - 1.0f0)
    end
    wk6 = zeros(Float32, max(idup * nptids, 2 * nofspe, 2 * maxint(AK_ES_MAXTPP), 6) + 1)
    @inbounds for i in 1:(idup * nptids)
        wk6[i] = esrann!(rng)
    end
    if est.ntally == 1
        if st.load == 1
            for i in 1:nptids; st.ipprep[i] = st.ipprep[i]; end
            if idup >= 2
                for i in 1:idup, n in 1:nptids
                    j = i * nptids + n; j <= length(st.ipprep) && (st.ipprep[j] = st.ipprep[n])
                end
            end
        else
            if st.ialn[2] == 0 && st.ialn[3] == 0
                st.pnone = 1.0f0; st.pmech = 0f0; st.pburn = 0f0         # ak/esprep.f
            else
                sm = st.pmech + st.pburn
                if sm > 1.0f0
                    st.pmech = st.pmech / sm; st.pburn = st.pburn / sm
                end
                st.pnone = 1.0f0 - st.pmech - st.pburn
            end
            sm = st.pnone + st.pmech + st.pburn
            sumup3 = Float32[st.pnone / sm, st.pmech / sm, st.pburn / sm]
            n = 0
            for ii in 1:idup, nn in 1:nptids
                n += 1
                draw = wk6[n]
                draw = draw * (((dupnpt + 1.0f0) - Float32(n)) / dupnpt)
                s0 = 0f0; sel = 0
                for i in 1:2
                    s0 += sumup3[i]
                    draw > s0 && continue
                    sel = i; break
                end
                if sel == 0
                    sel = 3; sumup3[3] < 0f0 && (sel = 1)
                end
                st.ipprep[n] = sel
                sumup3[sel] = sumup3[sel] - (1.0f0 / dupnpt)
                sumup3[sel] < 0f0 && (sumup3[sel] = 0f0)
            end
        end
    end
    # PLANT/NATURAL this cycle (OPFIND(2,MYACTS)); zero-TPA/survival ones are deleted.
    plants = [a for a in ctl.schedule if (a.icflag == Int32(430) || a.icflag == Int32(431)) &&
                                         _ak_in_cycle(_ak_act_date(s, a), icyc, year, next_year)]
    filter!(a -> !(a.params[2] <= 0.001f0 || a.params[3] <= 0.001f0), plants)
    mode = isempty(plants) ? 0 : 1
    # shade sums (loops 245/2451) and the per-point overstory DENSE statistics (gradd.f:192 DENSE, IND1 order)
    baaa = zeros(Float32, AK_ES_MAXPLT); ptpaa = zeros(Float32, AK_ES_MAXPLT); prda = zeros(Float32, AK_ES_MAXPLT)
    pccf = zeros(Float32, AK_ES_MAXPLT); over = zeros(Float32, nofspe, AK_ES_MAXPLT)
    let pi_f = p.pi, gross = p.gross_space, xmp = st.xmaxpt
        @inbounds for i in _ind1_order(s)
            ip = Int(t.plot_id[i]); (1 <= ip <= AK_ES_MAXPLT) || continue
            d = t.dbh[i]; pr = t.tpa[i]; sp = Int(t.species[i])
            pccf[ip] += ak_tree_ccf(sp, d) * pr * pi_f / gross
            d < AK_ES_REGNBK && continue
            batree = 0.005454154f0 * (d * (d * pr))
            baaa[ip] += batree * pi_f / gross
            over[sp, ip] += batree * pi_f / gross
            ptpaa[ip] += pr * pi_f / gross
            xm = ip <= length(xmp) ? xmp[ip] : st.xmax
            xm > 0f0 && (prda[ip] += (pr * fpow(d / 10f0, 1.605f0) * pi_f / gross) / xm)
        end
    end
    sum1 = 0f0; sum2 = 0f0
    for nn in 1:nptids
        tb = baaa[nn]; tb < 1f0 && (tb = 1f0); sum1 += tb
    end
    for nn in 1:nptids
        tb = baaa[nn]; tb < 1f0 && (tb = 1f0); sum2 += sum1 / tb
    end
    # stand/point site
    xxslp = Float32(p.slope_raw) * 0.01f0; xxasp = Float32(p.aspect_deg) * 0.0174533f0
    # estab.f:809-810 RADIAN=PASP(NNID), SLO=PSLO(NNID) — PER-PLOT. DATABASE tree input calls INTREE with IRDPLV=2
    # (dbsin.f:361 ⇒ NPNVRS=5, IPINFO=2): esplt1.f:69-70 stores each plot's FIRST tree record SLOPE/ASPECT and
    # esplt2.f:218-227 scales them (×0.01, ×0.0174533; a NULL column reads 0 ⇒ 0) — the FIA reader's p.point_slope/
    # point_aspect. Plots appended for IPTINV−NONSTK > found (esplt2.f:184-193, PSLO=−1 ⇒ XXSLP) and a bare
    # (IPINFO=0) stand take the stand XXSLP/XXASP — the same resolution as the IE port (inlandempire/establishment.jl).
    # MEASURED FIA 10708179010497 (stand slope 16, first tree 33/190°): live ESTOCK ADJSLO=33 ADJXC=−32.4986, jl had
    # the stand 16 ⇒ PN 1.91895 vs 1.88599 ⇒ PROB1 1-2 ULP off ⇒ every 2017 ES record PROB off and a WH/RC species swap.
    pslo_es = p.point_slope; pasp_es = p.point_aspect
    if s.trees.n > 0 && !isempty(pslo_es)
        if any(k -> k > length(p.point_ids) || p.point_ids[k] == 0, 1:length(pslo_es))
            pslo_es = copy(pslo_es); pasp_es = copy(pasp_es)
            @inbounds for k in 1:length(pslo_es)
                (k > length(p.point_ids) || p.point_ids[k] == 0) && (pslo_es[k] = xxslp; pasp_es[k] = xxasp)
            end
        end
    end
    tlat = p.latitude
    xesmlt = Float32[get(est.spec_mult, Int32(j), 1f0) for j in 1:nofspe]
    htadj = Float32[get(est.ht_adj, Int32(j), 0f0) for j in 1:nofspe]
    stoadj = est.stoadj
    maxtpp = AK_ES_MAXTPP[ift]
    ncount = 0
    gentim = fint - 5.0f0; gentim < 0f0 && (gentim = 0f0)
    padv = zeros(Float32, nofspe); psub = zeros(Float32, nofspe); pxcs = zeros(Float32, nofspe)
    first2 = fill(0.1f0, nofspe)
    icode = zeros(Int, 99); height = zeros(Float32, 99); esprob = zeros(Float32, 99); iasep = zeros(Int, 99)
    htimlt = ones(Float32, 99)
    tall = zeros(Float32, nofspe); stomlt = ones(Float32, nofspe); sumup = zeros(Float32, nofspe)
    ibest = zeros(Int, nofspe); ichoi = zeros(Int, 2, nofspe); jnull = zeros(Int, nofspe); ilsp = zeros(Int, nofspe)
    pspe = zeros(Float32, 6)
    wk6b = zeros(Float32, max(2 * nofspe, 2 * maxtpp, 6) + 1)
    itpp = 0; newtpp = 0; itp = 0
    esave = 0f0; slo = 0f0; baa = 1f0; nnid = 1
    ntally_now = Int(est.ntally)
    for nn in 1:nptids
        nnprep = zeros(Int, 4)
        istart = nn * idup - idup + 1; iend = nn * idup
        for i in istart:iend
            n = st.ipprep[i]; (n < 1 || n > 4) && (n = 1)
            nnprep[n] += 1
        end
        for itypep in 1:4
            ipold = 0
            nnprep[itypep] < 1 && continue
            ntimes = nnprep[itypep]
            for irep in 1:ntimes
                ncount += 1
                iprep = itypep
                if iprep != ipold
                    ipold = iprep
                    nnid = nn                                    # IPTIDS(NN) = NN (IPINFO=0)
                    st.prob1[ncount] = 1.0f0
                    radian = nnid <= length(pasp_es) ? pasp_es[nnid] : xxasp
                    slo    = nnid <= length(pslo_es) ? pslo_es[nnid] : xxslp
                    xcosas = fcos(radian); xcos = xcosas * slo
                    baa = baaa[nnid]; baa < 1f0 && (baa = 1f0); baa > 400f0 && (baa = 400f0)
                    baaold = nnid > length(est.inv_point_baaold) ? 0f0 : est.inv_point_baaold[nnid]
                    pta = ptpaa[nnid]
                    tpaold = nnid > length(est.esb_shift_pt) ? 0f0 : est.esb_shift_pt[nnid]
                    baaold < 1f0 && (baaold = 1f0); baaold > 400f0 && (baaold = 400f0)
                    if stoadj >= 0.0001f0
                        if !(est.inadv || ntally_now != 1)
                            ftemp = Float32(ctl.cycle_year[1]); ftemp < zharv && (ftemp = zharv)
                            st.time = iprep == 2 ? ftemp - st.zmech : iprep == 3 ? ftemp - st.zburn : ftemp - zharv
                            st.time < 0f0 && (st.time = 0f0)
                            pn = ak_estock(elev, ift, baaold, tpaold, slo, xcos, elevsq)
                            st.esb1[ncount] = pn
                            ftemp = fexp(pn) / (1f0 + fexp(pn))
                            itp0 = unsafe_trunc(Int, (st.plprob[nnid] * dupnpt) / (ftemp * expand) + 0.5f0)
                            st.nstore[ncount] = itp0; st.pnn[ncount] = st.esa
                            for i in 1:(ntimes - 1)
                                st.nstore[ncount + i] = itp0; st.pnn[ncount + i] = st.esa
                            end
                        end
                        st.time = iprep == 2 ? flokdt + 1f0 - st.zmech : iprep == 3 ? flokdt + 1f0 - st.zburn :
                                  flokdt + 1f0 - zharv
                        pn = ak_estock(elev, ift, baa, pta, slo, xcos, elevsq)
                        stoadj < 0.001f0 && (stoadj = 0.001f0; est.stoadj = stoadj)
                        ftemp = 1.0f0 / (1.0f0 + fexp(-(pn + st.esb - st.esb1[ncount])))
                        ftemp = ftemp * stoadj
                        ftemp < 0.0001f0 && (ftemp = 0.0001f0); ftemp > 0.9990f0 && (ftemp = 0.9990f0)
                        ftemp < st.pnn[ncount] && (ftemp = st.pnn[ncount] + 0.0001f0)
                        st.prob1[ncount] = ftemp
                        itp1 = 0
                        if ingro == 1 || match == 1
                            itp1 = unsafe_trunc(Int, (st.plprob[nnid] * dupnpt) / (ftemp * expand) + 0.5f0)
                            st.nstore[ncount] = itp1
                        end
                        for i in 1:(ntimes - 1)
                            st.prob1[ncount + i] = ftemp
                            st.esb1[ncount + i] = st.esb1[ncount]
                            (ingro == 1 || match == 1) && (st.nstore[ncount + i] = itp1)
                        end
                    end
                end
                # label 137
                st.time = ingro == 1 ? shorty : 10.0f0
                itime = unsafe_trunc(Int, st.time + 0.5f0)
                fill!(padv, 0f0); fill!(psub, 0f0); fill!(pxcs, 0f0)
                ak_espadv!(padv, ift, prda[nnid], view(over, :, nnid), baaa[nnid], slo, elev, tlat, xesmlt)
                itime > 2 && ak_espadv!(psub, ift, prda[nnid], view(over, :, nnid), baaa[nnid], slo, elev, tlat, xesmlt)
                ak_espadv!(pxcs, ift, prda[nnid], view(over, :, nnid), baaa[nnid], slo, elev, tlat, xesmlt)
                st.time = iprep == 2 ? flokdt + 1f0 - st.zmech : iprep == 3 ? flokdt + 1f0 - st.zburn :
                          flokdt + 1f0 - zharv
                d1 = esrann!(rng)
                emsqr = d1 < 0.5f0 ? -1.0f0 : 1.0f0
                emsqr = emsqr * esrann!(rng)
                if stoadj < 0.0001f0
                    esrann!(rng)
                    for _ in 1:(6 + 6 + nofspe + 2 * nofspe + 2 * maxtpp); esrann!(rng); end
                    itpp = 0; newtpp = 0
                else
                    draw = esrann!(rng)
                    ntally_now == 1 && (st.xstore[ncount] = draw)
                    draw = st.xstore[ncount]
                    # ak/estpp.f
                    tpp = AK_ES_TPP_BB[ift] * fpow(-1f0 * flog(1f0 - draw), 1f0 / AK_ES_TPP_CC[ift])
                    ift == 14 && (tpp = 0f0)
                    itpp = unsafe_trunc(Int, tpp + 0.5f0)
                    itpp < 1 && (itpp = 1)
                    itpp > maxtpp && (itpp = maxtpp)
                    (ingro == 1 && itpp > AK_ES_MAXING[ift]) && (itpp = AK_ES_MAXING[ift])
                    newtpp = itpp - st.nstore[ncount]; newtpp < 0 && (newtpp = 0)
                    st.nstore[ncount] = itpp
                    # STEP 6: number of species
                    for i in 1:6
                        pspe[i] = 0f0
                        wk6b[i] = esrann!(rng)
                    end
                    numspe = 1
                    if itpp != 1
                        ritpp = Float32(itpp)
                        for k in 1:min(itpp, 6)
                            pspa = AK_ES_NSSPB1[ift, k] + AK_ES_NSSPB2[k] * ritpp
                            pspe[k] = fexp(pspa) / (1 + fexp(pspa))
                        end
                        s0 = pspe[1]
                        for i in 2:6; s0 += pspe[i]; end
                        sumup[1] = pspe[1] / s0
                        for i in 2:5; sumup[i] = sumup[i - 1] + pspe[i] / s0; end
                        numspe = 6
                        for i in 1:5
                            wk6b[i] > sumup[i] && continue
                            numspe = i; break
                        end
                    end
                    numspe > AK_ES_MAXSPP[ift] && (numspe = AK_ES_MAXSPP[ift])
                    # STEP 7: which species
                    s0 = 0f0; nspnz = 0
                    for i in 1:nofspe
                        ibest[i] = 0
                        ftemp = padv[i] + psub[i]
                        ftemp > 0.0001f0 && (nspnz += 1)
                        s0 += ftemp; sumup[i] = ftemp
                    end
                    for i in 1:nofspe
                        s0 < 0.0001f0 && continue
                        sumup[i] = sumup[i] / s0
                    end
                    numspe > nspnz && (numspe = nspnz)
                    for i in 1:6; wk6b[i] = esrann!(rng); end
                    for i in 1:numspe
                        s0 = 0f0; draw = wk6b[i]; picked = false
                        for j in 1:(nofspe - 1)
                            s0 += sumup[j]
                            draw > s0 && continue
                            sumup[j] = 0f0; ibest[j] = 1; picked = true; break
                        end
                        if !picked
                            ibest[nofspe] = 1; sumup[nofspe] = 0f0
                        end
                        s0 = 0f0
                        for k in 1:nofspe; s0 += sumup[k]; end
                        s0 < 0.0001f0 && continue
                        for k in 1:nofspe; sumup[k] = sumup[k] / s0; end
                    end
                    fill!(padv, 0f0); fill!(psub, 0f0)
                    ntally_now == 1 && ak_espadv!(padv, ift, prda[nnid], view(over, :, nnid), baaa[nnid], slo, elev, tlat, xesmlt)
                    ak_espadv!(psub, ift, prda[nnid], view(over, :, nnid), baaa[nnid], slo, elev, tlat, xesmlt)
                    fill!(ichoi, 0)
                    for i in 1:nofspe
                        draw = esrann!(rng)
                        ibest[i] != 1 && continue
                        s0 = padv[i] + psub[i]
                        s0 < 0.0001f0 && continue
                        ftemp = padv[i] / s0
                        j = draw > ftemp ? 2 : 1
                        ichoi[j, i] = 1
                    end
                    # STEP 8: heights
                    st.time = flokdt - Float32(st.kdtold)
                    for i in 1:(2 * nofspe); wk6b[i] = esrann!(rng); end
                    ndraw = 0
                    for i in 1:nofspe
                        ndraw += 1
                        tall[i] = 0.001f0; stomlt[i] = 1.0f0
                        xsite = p.sp_site_index[i]
                        if ichoi[1, i] == 1
                            draw = wk6b[ndraw]
                            st.delay = ak_esdlay(i, 1, draw, st.time, baa)
                            hht, st.delay, trage = ak_esadvh(i, st.delay, gentim, xsite)
                            tall[i] = hht; agadsb[i] = trage
                            ftemp = trage; ftemp > gentim && (ftemp = gentim)
                            stomlt[i] = ftemp / (gentim + 0.0001f0)
                        end
                        if ichoi[2, i] == 1
                            ndraw += 1
                            draw = wk6b[ndraw]
                            st.delay = ak_esdlay(i, 2, draw, st.time, baa)
                            hht, st.delay, trage = ak_essubh(i, st.delay, gentim, 0f0, st.time, xsite)
                            tall[i] = hht; agadsb[i] = trage
                            ftemp = trage; ftemp > gentim && (ftemp = gentim)
                            stomlt[i] = ftemp / (gentim + 0.0001f0)
                        end
                    end
                    for i in 1:nofspe
                        ibest[i] != 1 && continue
                        tall[i] = tall[i] + htadj[i]
                        tall[i] - AK_ES_XMIN[i] < 0.2f0 && (tall[i] = AK_ES_XMIN[i] + 0.2f0)
                        tall[i] > AK_ES_HHTMAX[i] && (tall[i] = AK_ES_HHTMAX[i])
                    end
                    wk6x = zeros(Float32, nofspe)
                    for i in 1:nofspe
                        wk6x[i] = pxcs[i] * Float32(ibest[i])
                        ibest[i] != 1 && continue
                        wk6x[i] < 0.0001f0 && (wk6x[i] = 0.0001f0)
                    end
                    s0 = 0f0
                    for i in 1:nofspe; s0 += wk6x[i]; end
                    if s0 <= 0f0
                        fill!(sumup, 0f0)
                    else
                        sumup[1] = wk6x[1] / s0
                        for i in 2:(nofspe - 1); sumup[i] = sumup[i - 1] + wk6x[i] / s0; end
                    end
                    itp = itpp
                    for i in 1:(itpp + nofspe)
                        i > 99 && break
                        note[i] = 0; icode[i] = 0; iasep[i] = 0; height[i] = 0f0; esprob[i] = 0f0; htimlt[i] = 1f0
                    end
                    n = 0
                    for i in 1:nofspe
                        jnull[i] = 0; ilsp[i] = 0
                        ibest[i] == 0 && continue
                        n += 1
                        icode[n] = i; iasep[n] = ichoi[1, i] == 1 ? 1 : 2
                        height[n] = tall[i]; htimlt[n] = stomlt[i]
                    end
                    for i in 1:(2 * maxtpp); wk6b[i] = esrann!(rng); end
                    ndraw = 0
                    for n2 in (numspe + 1):itpp
                        ndraw += 1
                        draw = wk6b[ndraw]
                        ii = nofspe; ilast = nofspe
                        for i in 1:(nofspe - 1)
                            if !(draw > sumup[i])
                                ii = i; ilast = i; break
                            end
                        end
                        icode[n2] = ii; iasep[n2] = 3
                        ndraw += 1
                        ftemp = flokdt - Float32(st.kdtold)
                        hht = ak_esxcsh(ii, p.sp_site_index[ii])
                        height[n2] = hht + htadj[ii]
                        agexc[n2] = ftemp - st.delay - gentim
                        height[n2] < AK_ES_XMIN[ii] && (height[n2] = AK_ES_XMIN[ii])
                        height[n2] > AK_ES_HHTMAX[ii] && (height[n2] = AK_ES_HHTMAX[ii])
                        htimlt[n2] = stomlt[ilast]
                    end
                    ftemp = st.prob1[ncount]
                    ftemp2 = itpp > 0 ? Float32(newtpp) / Float32(itpp) : 0f0
                    itemp = itpp - newtpp
                    for i in 1:itpp
                        esprob[i] = ftemp
                        i < itemp + 1 && (esprob[i] = ftemp - st.pnn[ncount])
                        (ingro == 1 || match == 1) && (esprob[i] = ftemp * ftemp2)
                        esprob[i] < 0.0001f0 && (esprob[i] = 0.0001f0)
                    end
                    st.pnn[ncount] = ftemp
                end
                # label 163
                esave = Float32(unsafe_trunc(Int32, esrann!(rng) * 100000.0f0 + 0.5f0))
                itp = itpp
                if mode == 1
                    itodo = 0
                    for a in plants
                        itodo += 1
                        itpp + itodo > 99 && (itodo -= 1; continue)
                        ipnspe = round(Int, a.params[1])
                        ptree = (a.params[2] * (a.params[3] / 100.0f0)) / dupnpt
                        trage = a.params[4]
                        trage < 0.5f0 && (trage = 1.0f0)
                        trage > 2.0f0 && (trage = 2.0f0)
                        agepl[itodo] = trage
                        ipyear = _ak_act_date(s, a)
                        st.delay = Float32(ipyear) - Float32(kdt + 1 - fint)
                        treeht = a.params[5]
                        ishade = unsafe_trunc(Int, a.params[6])
                        ftemp = baa; ftemp < 1.0f0 && (ftemp = 1.0f0)
                        if ishade == 1
                            ptree = ptree * nptids * ftemp / sum1
                        elseif ishade == 2
                            ptree = ptree * nptids * (sum1 / ftemp) / sum2
                        end
                        tmtime = st.time
                        st.time = fint
                        hht, st.delay, trage = ak_essubh(ipnspe, st.delay, gentim, trage, st.time, p.sp_site_index[ipnspe])
                        st.time = tmtime
                        if treeht >= 0.1f0
                            hht = treeht
                            xh = flog(hht)
                            xxh = 0f0
                            while true
                                xxh = fexp(bachlo(rng, xh, 0.5f0; stream = :estab))
                                (xxh < 0.5f0 * hht || xxh > 2.0f0 * hht) || break
                            end
                            hht = xxh + htadj[ipnspe]
                            hht < 0.05f0 && (hht = 0.05f0)
                        else
                            hht = hht + htadj[ipnspe]
                            hht < AK_ES_XMIN[ipnspe] && (hht = AK_ES_XMIN[ipnspe])
                        end
                        hht > AK_ES_HHTMAX[ipnspe] && (hht = AK_ES_HHTMAX[ipnspe])
                        k = itpp + itodo
                        icode[k] = ipnspe
                        iasep[k] = a.icflag == Int32(430) ? 3 : 4
                        height[k] = hht; esprob[k] = ptree
                        ftemp = trage; ftemp > gentim && (ftemp = gentim)
                        htimlt[k] = ftemp / (gentim + 0.0001f0)
                        agepl[k] = fint - st.delay - gentim + agepl[itodo]
                    end
                    itp = itpp + itodo
                end
                # label 321
                (esave % 2f0 == 0f0) && (esave += 1f0)
                rng.es0 = Float64(esave)
                if stoadj >= 0.0001f0
                    if itp < 5
                        for i in 1:itp; note[i] = 1; end
                    else
                        istart2 = st.prob1[ncount] < 0.00011f0 ? itpp + 1 : 1
                        nbest = 0; num2pk = itp - istart2 + 1
                        itemp2 = 0
                        done = false
                        for _j in 1:2
                            nbest >= num2pk && (done = true; break)
                            ftemp = 0.001f0
                            for i in istart2:itp
                                note[i] == 1 && continue
                                height[i] < ftemp && continue
                                ftemp = height[i]; itemp2 = i
                            end
                            jnull[icode[itemp2]] = 1
                            note[itemp2] = 1
                            nbest += 1
                        end
                        if !done
                            for j in 1:nofspe, i in istart2:itp
                                icode[i] != j && continue
                                note[i] == 1 && continue
                                height[i] < tall[j] && continue
                                tall[j] = height[i]; ilsp[j] = i
                            end
                            for j in 1:nofspe
                                jnull[j] == 1 && continue
                                itemp2 = ilsp[j]; itemp2 == 0 && continue
                                note[itemp2] = 1; nbest += 1
                            end
                            if !(nbest >= 4 || nbest >= num2pk)
                                while true
                                    ftemp = 0.001f0
                                    for j in istart2:itp
                                        note[j] == 1 && continue
                                        height[j] < ftemp && continue
                                        ftemp = height[j]; itemp2 = j
                                    end
                                    nbest += 1; note[itemp2] = 1
                                    nbest >= num2pk && break
                                    nbest < 4 || break
                                end
                            end
                        end
                    end
                end
                # label 229: records
                excess = zeros(Float32, nofspe); sumhts = zeros(Float32, nofspe); sumesp = zeros(Float32, nofspe)
                for n2 in 1:itpp
                    i = icode[n2]
                    hht = height[n2]
                    if esprob[n2] >= 0.00011f0 && note[n2] == 1
                        k = _ak_new_record!(s, i, nnid, hht, (esprob[n2] * expand) / dupnpt, 1, agadsb[n2],
                                            htimlt[n2], Int32(est.idsdat + 20), itmp)
                        k == 0 && break
                    else
                        excess[i] += 1f0; sumhts[i] += hht; sumesp[i] += esprob[n2]
                    end
                end
                for i in 1:nofspe
                    excess[i] < 0.5f0 && continue
                    if sumesp[i] < 0.00011f0 && i < nofspe
                        sumesp[i + 1] += sumesp[i]; sumesp[i] = 0f0; continue
                    end
                    ftemp2 = sumesp[i] / (excess[i] + 0.00001f0)
                    ibrkup = unsafe_trunc(Int, excess[i] / 5.0f0 + 1.0f0); brkup = Float32(ibrkup)
                    for _ in 1:ibrkup
                        hht = sumhts[i] / excess[i]
                        xcsmax = excess[i] / brkup
                        ftemp = pasmax / brkup
                        xcsmax > ftemp && (xcsmax = ftemp)
                        pr = (ftemp2 * expand * xcsmax) / dupnpt
                        k = _ak_new_record!(s, i, nnid, hht, pr, 2, agexc[i], stomlt[i], Int32(0), itmp)
                        k == 0 && break
                        pastpa[i] += pr; tcrop2 += pr
                    end
                end
                if mode == 1
                    jcnt = itp - itpp
                    for jj in 1:jcnt
                        i = icode[itpp + jj]
                        ftemp = esprob[itpp + jj]
                        if ftemp < 0.00011f0 && jj < jcnt
                            esprob[itpp + jj + 1] += ftemp; esprob[itpp + jj] = 0f0; continue
                        end
                        ibrkup = unsafe_trunc(Int, ftemp / 10.0f0 + 1.0f0); brkup = Float32(ibrkup)
                        tottpa[i] += ftemp; pastpa[i] += ftemp; bestpa[i] += ftemp
                        tcrop1 += ftemp; tcrop2 += ftemp; ttottp += ftemp
                        imc = (note[itpp + jj] != 1 && stoadj > 0f0) ? 2 : 1
                        for _ in 1:ibrkup
                            k = _ak_new_record!(s, i, nnid, height[itpp + jj], ftemp / brkup, imc, agepl[itpp + jj],
                                                htimlt[itpp + jj], Int32(0), itmp)
                            k == 0 && break
                        end
                    end
                end
            end
        end
    end
    nadd = t.n - itrnin + 1
    if nadd > 0
        ak_esgent!(s, itrnin - 1; fint = fint, pccf_pre = pccf)
        @inbounds for i in itrnin:t.n
            t.birth_age[i] = t.birth_age[i] + gentim           # estab.f:2330 ABIRTH=ABIRTH+GENTIM
        end
    end
    st.kdtold = Int(kdt)
    return nadd
end

@inline maxint(v) = maximum(v)

"Append one ESTAB record (estab.f:2002-2063 / 2120-2177 / 2218-2278). Returns its index (0 if the list is full)."
function _ak_new_record!(s::StandState, sp::Int, nnid::Int, hht::Float32, prob::Float32, imc::Int, abirth::Float32,
                         wk4::Float32, iestat::Int32, itmp::Int)::Int
    t = s.trees
    n = t.n + 1
    n + Int(t.ndead) > length(t.dbh) && return 0
    t.n = n
    t.species[n] = Int32(sp); t.plot_id[n] = Int32(nnid)
    t.mort_code[n] = Int32(imc)
    t.cuft_vol[n] = 0f0; t.merch_cuft_vol[n] = 0f0; t.saw_cuft_vol[n] = 0f0; t.bdft_vol[n] = 0f0
    t.cull[n] = 0f0; t.decay_code[n] = Int32(0); t.woodland_stems[n] = Int32(0)
    t.trunc[n] = Int32(0); t.defect[n] = Int32(0); t.special[n] = Int32(0); t.norm_ht[n] = Int32(0)
    t.merch_top_bf[n] = 0f0; t.merch_top_cf[n] = 0f0
    t.mort_pa[n] = 0f0                                   # estab.f:2050/2169/2266 WK2(ITRN)=0 (TreeList MortPA)
    t.tpa[n] = prob
    t.dbh[n] = 0.1f0
    t.height[n] = hht
    t.birth_age[n] = abirth
    t.zrand[n] = -999f0
    t.crown_pct[n] = Int32(0); t.crown_ratio[n] = 0f0
    t.diam_growth[n] = 0f0; t.ht_growth[n] = 0f0
    t.old_crown_pct[n] = 0f0; t.old_random[n] = 0f0
    t.dg_prev[n] = 0f0
    t.htimlt[n] = wk4
    t.iestat[n] = iestat
    t.tree_id[n] = Int32(AK_ES_IDCMP1 + itmp + n)
    t.sort_key[n] = Float64(n)
    t.history[n] = Int32(0); t.cut_code[n] = Int32(0)
    return n
end
