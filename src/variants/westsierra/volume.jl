# =============================================================================
# volume.jl (westsierra) — WS volume (ws VEQNNC). Chunk 8. R5 California, forest 500.
#
# WS VEQNNC (dumped from FVSws_g16 wst01 sitset VOLEQDEF): conifers = "500WO2W<fia>" (R5 Wensel/Krumland
# taper, r5tap.f) ; hardwoods/woodland = "500DVEW<fia>" (R5 California hardwood D²H, r5harv.f). These are the
# SAME two profile kernels NC/Klamath ported (nc_r5tap_dib taper + nc_r5harv_vol). WS reuses the shared
# R5TAP taper (nc_r5tap_dib) + the R6/R3 log-bucking helpers (nc_wo2w_merch, _fw2_board, _fw2_tcubic).
#
# ★★ MEASURE CATCH (fvsvol.f:168-205 + grinit.f:160-170, PROVEN vs FVSws_g16 per-tree MCF/BBFV dump):
#   the R5 merch TOP DIAMETERS are NOT a fixed 6" — they are TOPD/BFTOPD (DOB) × BRATIO (species bark, DIB):
#     • cubic merch VOL(4):  top DIB = TOPD(=4.5)  · ws_bratio(sp,d)   (fvsvol.f:195 TOPDIAM=TOPD·BARK)
#     • board  Scribner VOL(2): top DIB = BFTOPD(=6.0) · ws_bratio(sp,d)  (fvsvol.f:663 TDIBB=BFTOPD·BRATIO)
#   e.g. SP(bark 0.8863): cubic 4.5·0.8863=3.99, board 6.0·0.8863=5.318. Using a flat 6" top under-counts
#   merch ~50% on small trees. Merch bucking rules = R3/R5 (mrules REGN 5): OPT=22, EVOD=2, MAXLEN=16,
#   MINLEN=2, TRIM=0.5, MERCHL=8, STUMP=1. Gates (grinit.f): DBHMIN=7 (MCF=0 below), BFMIND=10 (board=0).
#   SPFLG=0 on wst01 ⇒ MCF=VOL(4) only (no VOL(7) topwood). MEASURED bit-exact vs FVSws_g16 gated MCF/BBFV.
# =============================================================================

# ws/sitset.f VEQNNC (dumped from FVSws_g16), index = ISPC.
const WS_VOL_EQ = String[
  "500WO2W117","500WO2W202","500WO2W015","500DVEW212","500WO2W081",   # SP DF WF GS IC
  "500WO2W116","500WO2W020","500WO2W122","500WO2W108","500WO2W108",   # JP RF PP LP WB
  "500WO2W117","500WO2W116","500WO2W015","500WO2W108","500WO2W108",   # WP PM SF KP FP
  "500WO2W108","500WO2W108","500WO2W108","500WO2W108","500WO2W117",   # CP LM MP GP WE
  "500WO2W108","500WO2W202","500WO2W211","500WO2W015","500DVEW060",   # GB BD RW MH WJ
  "500DVEW060","500DVEW060","500DVEW801","500DVEW805","500DVEW807",   # UJ CJ LO CY BL
  "500DVEW818","500DVEW821","500DVEW839","500DVEW631","500DVEW431",   # BO VO IO TO GC
  "500DVEW818","500DVEW981","500DVEW361","500DVEW807","500DVEW312",   # AS CL MA DG BM
  "500DVEW801","500WO2W108","500DVEW821"]                             # MC OS OH

const WS_VOL_TOPD  = 4.5f0     # grinit.f:160 TOPD  (cubic merch top DOB)
const WS_VOL_BFTOPD = 6.0f0    # grinit.f:166 BFTOPD (board top DOB)
const WS_VOL_DBHMIN = 7.0f0    # grinit.f:161 DBHMIN (cubic merch gate)
const WS_VOL_BFMIND = 10.0f0   # grinit.f:167 BFMIND (board gate)

function compute_volumes_ws!(s::StandState)
    t = s.trees; sd = s.coef.species
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 43
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = WS_VOL_EQ[sp]; mdl = eq[4:6]
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        bark = ws_bratio(sd, sp, d)                       # DIB/DOB (fvsvol BARK=BRATIO)
        if mdl == "WO2"
            s5 = _nc_r5tap_sp(eq[8:10])
            if s5 == 0 || hv < 5f0
                t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
                t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
            end
            dibat = ht -> nc_r5tap_dib(s5, d, hv, Float32(ht))
            t.cuft_vol[i] = Float32(round(_fw2_tcubic(dibat, hv) * 10.0)) / 10f0     # VOL(1)
            t.merch_cuft_vol[i] = d >= WS_VOL_DBHMIN ?
                nc_wo2w_merch(dibat, hv; mtopp = WS_VOL_TOPD * bark,
                              stump = 1f0, minlen = 2f0, merchl = 8f0) : 0f0          # VOL(4), top=TOPD·bark
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= WS_VOL_BFMIND ?
                _fw2_board(dibat, hv, WS_VOL_BFTOPD * bark, 1f0, 2f0, 8f0) : 0f0      # VOL(2), top=BFTOPD·bark
        else                                              # DVE — California hardwood D²H (r5harv.f)
            tcf, mcf, _ = nc_r5harv_vol(eq, d, hv, WS_VOL_TOPD * bark)
            t.cuft_vol[i] = tcf
            t.merch_cuft_vol[i] = d >= WS_VOL_DBHMIN ? mcf : 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0   # DVE board deferred (as NC)
        end
    end
    return s
end
