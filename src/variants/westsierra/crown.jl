# =============================================================================
# crown.jl (westsierra) — WS crown ratio (ws/crown.f + ws/ccfcal.f). Chunk 5.
#
# Rank-based Weibull crown ratio, SAME framework as SO/EC/WC/PN (crown_ratio_update! mirrors so_crown), with
# WS-specific coefficients and — the key MEASURE catch — a SPECIES-GROUPED SCALE (ws/crown.f:364-371):
#   CASE(9:10,12,14:17,19:20,25:27) SCALE = 1.5 − RELSDI
#   CASE(1:8,11,13,18,22:24,28:40,42:43) SCALE = 1 − 0.00333·(RELDEN−50)   ← MAIN group (SP/DF/WF/RF …)
#   CASE(41) SCALE = 1 − 0.00167·(RELDEN−100)
# (NOT SO's single 1−0.00167·(RELDEN−100).) B floor 3, C floor 2; rank X=ISORT/ITRN·SCALE clamp[.05,.95];
# CRNEW=(A+B·(−ln(1−X))^(1/C))·10; ±1%/yr change limit; CRMAX cap; NORMHT topkill. GB(21) = CL=−0.59373+
# 0.67703·HF form; GS/RW(4,23) = logistic on ln(HDR)/PRD/(D/QMDPLT). DBH<1"@LSTART → ws/dubscr.f (STUB, chunk
# 5b — bypassed on wst01 whose crowns are all present). RELDEN = stand CCF via ws_tree_ccf (ws/ccfcal.f).
# MEASURED vs FVSws_g16 crown per-tree ICR.
# =============================================================================

# ── WS per-species Weibull crown coefficients (ws/crown.f CRCONS DATA, index = ISPC).
const WS_CROWN_WEIBA = Float32[
  0,0,0,0,0, 2.0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0]
const WS_CROWN_WEIBB0 = Float32[
  0.32957,0.39996,0.17606,0.32957,0.155, -1.2458,0.16601,0.20199,-0.13121,-0.13121,
  0.32957,0.16267,0.17606,0.16267,0.16267, 0.16267,0.16267,0.20199,0.16267,0.16267,
  0,0.39996,0.32957,0.32957,0.16267, 0.16267,0.16267,-0.14217,-0.14217,-0.14217,
  -0.14217,-0.14217,-0.14217,-0.14217,-0.14217, -0.14217,-0.14217,-0.14217,-0.14217,-0.14217,
  -0.23830,-0.098,-0.14217]
const WS_CROWN_WEIBB1 = Float32[
  1.03916,1.03150,1.07984,1.03916,1.08747, 0.94476,1.08150,1.07198,1.15976,1.15976,
  1.03916,1.07340,1.07984,1.07340,1.07340, 1.07340,1.07340,1.07198,1.07340,1.07340,
  0,1.03150,1.03916,1.03916,1.07340, 1.07340,1.07340,1.15448,1.15448,1.15448,
  1.15448,1.15448,1.15448,1.15448,1.15448, 1.15448,1.15448,1.15448,1.15448,1.15448,
  1.18016,1.11809,1.15448]
const WS_CROWN_WEIBC0 = Float32[
  -0.83314,-0.98287,-0.89140,-0.83314,0.85877, -10.54490,0.91420,0.75409,2.59824,2.59824,
  -0.83314,3.28850,-0.89140,3.28850,3.28850, 3.28850,3.28850,0.75409,3.28850,3.28850,
  0,-0.98287,-0.83314,-0.83314,3.28850, 3.28850,3.28850,0.59185,0.59185,0.59185,
  0.59185,0.59185,0.59185,0.59185,0.59185, 0.59185,0.59185,0.59185,0.59185,0.59185,
  3.04,4.05181,0.59185]
const WS_CROWN_WEIBC1 = Float32[
  0.91493,0.88449,0.76518,0.91493,0.40125, 2.45822,0.45768,0.52191,0,0,
  0.91493,0,0.76518,0,0, 0,0,0.52191,0,0,
  0,0.88449,0.91493,0.91493,0, 0,0,0.37245,0.37245,0.37245,
  0.37245,0.37245,0.37245,0.37245,0.37245, 0.37245,0.37245,0.37245,0.37245,0.37245,
  0,0,0.37245]
const WS_CROWN_C0 = Float32[
  7.12189,5.91609,6.86237,7.12189,6.32336, 7.33055,6.14578,6.15172,4.89032,4.89032,
  7.12189,6.48494,6.86237,6.48494,6.48494, 6.48494,6.48494,6.15172,6.48494,6.48494,
  0,5.91609,7.12189,7.12189,6.48494, 6.48494,6.48494,4.00579,4.00579,4.00579,
  4.00579,4.00579,4.00579,4.00579,4.00579, 4.00579,4.00579,4.00579,4.00579,4.00579,
  4.62512,6.35669,4.00579]
const WS_CROWN_C1 = Float32[
  -0.02817,-0.00943,-0.03278,-0.02817,-0.02987, -0.01539,-0.02781,-0.01994,-0.01884,-0.01884,
  -0.02817,-0.02325,-0.03278,-0.02325,-0.02325, -0.02325,-0.02325,-0.01994,-0.02325,-0.02325,
  0,-0.00943,-0.02817,-0.02817,-0.02325, -0.02325,-0.02325,-0.00522,-0.00522,-0.00522,
  -0.00522,-0.00522,-0.00522,-0.00522,-0.00522, -0.00522,-0.00522,-0.00522,-0.00522,-0.00522,
  -0.01604,-0.00846,-0.00522]

# ws/ccfcal.f DATA — per-species tree-CCF (crown-width² form for the WS-native species).
const WS_CCF_RD1 = Float32[
  6.74,6.81,5.82,5.82,7.11, 3.10,6.71,5.13,0,0, 6.74,0,5.82,0,0, 0,0,5.13,0,0,
  0,6.81,5.82,6.74,0, 0,0,10,10,10, 10,10,10,10,10, 10,10,10,10,10, 0,7.07,10]
const WS_CCF_RD2 = Float32[
  0.623,0.732,0.591,0.591,0.470, 0.839,0.421,0.658,0,0, 0.623,0,0.591,0,0, 0,0,0.658,0,0,
  0,0.732,0.591,0.623,0, 0,0,1.200,1.200,1.200, 1.200,1.200,1.200,1.050,1.050,
  1.050,1.050,1.050,1.050,1.200, 0,0.551,1.200]
const WS_CCF_RDA = Float32[
  0.007244,0.017299,0.015248,0.011109,0.008915, 0.007875,0.011402,0.007813,0,0,
  0.007244,0,0.015248,0,0, 0,0,0.007813,0,0, 0,0.017299,0.011109,0.007244,0,
  0,0,0.009187,0.009187,0.009187, 0.009187,0.009187,0.009187,0.011109,0.011109,
  0.011109,0.011109,0.011109,0.011109,0.009187, 0,0.009884,0.009187]
const WS_CCF_RDB = Float32[
  1.8182,1.5571,1.7333,1.7250,1.7800, 1.7360,1.7560,1.7780,0,0, 1.8182,0,1.7333,0,0, 0,0,1.7780,0,0,
  0,1.5571,1.7250,1.8182,0, 0,0,1.7600,1.7600,1.7600, 1.7600,1.7600,1.7600,1.7250,1.7250,
  1.7250,1.7250,1.7250,1.7250,1.7600, 0,1.6667,1.7600]

# CA-variant surrogate species (need R5CRWD crown width for CCF/crown SCALE — chunk 5b stub).
const WS_CROWN_CA_SURR = Set{Int}([9,10,12,14,15,16,17,19,20,25,26,27])

# ws/ccfcal.f MODE=1 CCFT (per tree, before ×P). WS-native species use crown-width²; GB(21)/MC(41) specials.
@inline function ws_tree_ccf(sp::Integer, d::Real, h::Real)::Float32
    (sp < 1 || sp > 43) && return 0f0
    D = Float32(d)
    if sp == 41                                            # MC (SO surrogate)
        D < 1f0 && return D * (0.0204f0 + 0.0246f0 + 0.0074f0)
        return 0.0204f0 + 0.0246f0 * D + 0.0074f0 * D * D
    elseif sp == 21                                        # GB (UT surrogate)
        D >= 10f0 && return 0.01925f0 + 0.01676f0 * D + 0.00365f0 * D * D
        D > 0.1f0 && return 0.009187f0 * fpow(D, 1.7600f0)
        return 0.001f0
    elseif sp in WS_CROWN_CA_SURR                          # CA-surrogate (R5CRWD crown width — chunk 5b)
        return 0.001f0                                     # stub (absent from wst01)
    end
    # CASE DEFAULT (WS-native, incl SP/DF/WF/RF): crown-width² (RD1+D·RD2)² · 0.001803
    D >= 1f0 && return ((WS_CCF_RD1[sp] + D * WS_CCF_RD2[sp])^2) * 0.001803f0
    D > 0.1f0 && return WS_CCF_RDA[sp] * fpow(D, WS_CCF_RDB[sp])
    return 0.001f0
end

# ws/crown.f — rank-based Weibull crown ratio. sdiac = SDIAC (stand SDI); RELDEN = stand CCF.
function crown_ratio_update!(s::StandState, ::WestSierra; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    relden = p.relative_density; sdiac = crown_sdi
    dens = s.density
    # rank trees by projected DBH (ws/crown.f ISORT via RDPSRT on D+DG/BARK), descending → ISORT[i] ∈ 1..n
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = ws_bratio(sd, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 43) && continue
        (lstart && t.crown_pct[i] > 0) && continue         # crown present → bypass
        icr = Int(t.crown_pct[i])
        # RELSDI (per-species SDImax = sp_sdi_def, chunk 2b)
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = WS_CROWN_C0[sp] + WS_CROWN_C1[sp] * relsdi * 100f0
        A = WS_CROWN_WEIBA[sp]
        B = WS_CROWN_WEIBB0[sp] + WS_CROWN_WEIBB1[sp] * acrnew; B < 3f0 && (B = 3f0)
        C = WS_CROWN_WEIBC0[sp] + WS_CROWN_WEIBC1[sp] * acrnew; C < 2f0 && (C = 2f0)

        if d < 1f0 && lstart                               # ws/crown.f:327 → DUBSCR (chunk 5b stub) / GB CL
            icr != 0 && continue
            if sp == 21
                hf = h + t.ht_growth[i]
                cl = -0.59373f0 + 0.67703f0 * hf
                cl < 1f0 && (cl = 1f0); cl > hf && (cl = hf)
                icri = trunc(Int, (cl / hf) * 100f0 + 0.5f0)
            else
                # ws_dubscr STUB (chunk 5b) — bypassed on wst01 (crowns present). Faithful floor.
                icri = 10
            end
            icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri); continue
        end

        if sp == 21                                        # GB — CL form
            hf = h + t.ht_growth[i]
            cl = -0.59373f0 + 0.67703f0 * hf
            cl < 1f0 && (cl = 1f0); cl > hf && (cl = hf)
            crnew = (cl / hf) * 100f0
        elseif sp == 4 || sp == 23                         # GS/RW — logistic (PRD/QMDPLT/HDR)
            pt_i = Int(t.plot_id[i])
            baplt = (1 <= pt_i <= length(dens.point_ba)) ? dens.point_ba[pt_i] : 0f0
            tpaplt = (1 <= pt_i <= length(dens.point_tpa)) ? dens.point_tpa[pt_i] : 0f0
            qmdplt = tpaplt > 0f0 ? sqrt((baplt / tpaplt) / 0.005454f0) : 1f0
            qmdplt < 1f0 && (qmdplt = 1f0)
            prd = ws_point_prd(s, pt_i)
            hdr = (h * 12f0) / d
            x = -1.021064f0 + 0.309296f0 * log(hdr) + 0.869720f0 * prd - 0.116274f0 * (d / qmdplt)
            x = 1f0 / (1f0 + exp(x))
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = x * 100f0
        else                                               # CASE DEFAULT — rank-Weibull w/ grouped SCALE
            scale = if sp in WS_CROWN_CA_SURR
                1.5f0 - relsdi
            elseif sp == 41
                1f0 - 0.00167f0 * (relden - 100f0)
            else                                           # MAIN group: SP/DF/WF/RF …
                1f0 - 0.00333f0 * (relden - 50f0)
            end
            scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
            x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
        end

        if !(lstart || icr == 0)                           # crown CHANGE, ±1%/yr limit
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg
        end
        icri = trunc(Int, crnew + 0.5f0)
        if !(lstart || icr == 0)                           # CRMAX crown-length cap
            htg = t.ht_growth[i]; crln = h * Float32(icr) / 100f0
            crmax = (crln + htg) / (h + htg) * 100f0
            (icri < 10) && (icri = trunc(Int, crmax + 0.5f0))
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        if lstart && t.trunc[i] != 0                       # topkill NORMHT reduction
            hn = Float32(t.norm_ht[i]) / 100f0; hd = hn - Float32(t.trunc[i]) / 100f0
            cl = (Float32(icri) / 100f0) * hn - hd
            icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
        end
        icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# ---------------------------------------------------------------------------
# WS FFE crown-biomass species map (ws/fmcrow.f:108 DATA ISPMAP) — the Jenkins/FMCROWE
# crown-biomass group per species. ws/fmcrow.f:165 routes CASE(36,39,41)=AS/DG/MC → FMCROWE
# (eastern Jenkins TOTABV), all others → FMCROWW (western crown-width, shared cr_crownw). WS FFE chunk-0.
# ---------------------------------------------------------------------------
const WS_ISPMAP = Int[
  15, 3, 4, 19, 20, 15, 4, 13, 11, 14,
  15, 12, 4, 11, 11, 11, 11, 13, 11, 11,
   9, 3, 19, 24, 16, 16, 16, 17, 17, 21,
  21, 21, 17, 17, 17, 41, 17, 10, 56, 5,
  41, 11, 21]
@inline ws_uses_fmcrowe(sp::Integer) = (sp == 36 || sp == 39 || sp == 41)
