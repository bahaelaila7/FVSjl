# =============================================================================
# volume.jl (teton) — TT volume (chunk 8, COMPLETE). VOLEQ = 400MATW<code> (Region-4 Matney = r4vol) for
# conifers + 400DVEW/401DVEW (Chojnacky INT-339 D2H = R4D2H) for the woodland species PM/RM/MC/OH/UJ.
# Conifer cubic (total CF0 + merch CFGRS via r4vol_volumes) + board-foot Scribner (r4vol_board) validated
# bit-exact-or-cornered. Woodland R4D2H (r4d2h_vol1 below) validated per-tree BIT-EXACT vs FVStt: .sum
# MC/RM/OH bit-exact, PM/UJ ±1-2 cuft cornered (Float32 tpa/accumulation tail); NO Behre/CFTOPK trim.
# =============================================================================

# TT VOLEQ per species (ttt01.out; VOLCODE ≠ FIA — e.g. WB uses 108). geocode 400/401.
const TT_VOL_EQ = String[
    "400MATW108", "400MATW108", "400MATW202", "400DVEW133", "400MATW093", "400MATW746",
    "400MATW108", "400MATW093", "400MATW019", "400MATW122", "401DVEW065", "400DVEW066",
    "400MATW108", "400MATW108", "400MATW108", "400DVEW475", "400MATW108", "400DVEW998"]

# R4D2H (volume/NVEL/r4d2h.f): Chojnacky INT-339 D2H woodland cubic. VOL(1)=(a+b·D2H^⅓+c·MSTEM)³,
# D2H=DBH²·HTTOT (DRC=0 in the FVS call ⇒ DBH used), MSTEM=0 (FCLASS=0 ⇒ c-term drops). UJ/PM get a
# 0.1 floor for DBH<3. Returns VOL(1) (= the entire-cubic TCF; VOL(4)=VOL(1)). eq = 10-char VOLEQ.
@inline function r4d2h_vol1(eq::AbstractString, d::Float32, h::Float32)
    (d <= 0f0 || h <= 0f0) && return 0f0
    d2h = d * d * h
    c = fpow(d2h, 1f0 / 3f0)                                 # D2H**(1./3.) via gfortran powf (doctrine #8)
    # NB: Fortran writes `(...)**3.` — a REAL exponent ⇒ powf(x,3.0), NOT x*x*x; match via fpow(.,3f0).
    code = strip(eq)[8:10]
    if code == "064"          # Western Juniper (WJ) — r4d2h.f:63-65 (no DBH<3 floor)
        return fpow(-0.22048f0 + 0.125468f0 * c, 3f0)
    elseif code == "106"      # Pinyon Pine (PI) — r4d2h.f:101-103 (no DBH<3 floor)
        return fpow(-0.20296f0 + 0.150283f0 * c, 3f0)
    elseif code == "066"      # Rocky Mountain Juniper (RM)
        return fpow(0.02434f0 + 0.119106f0 * c, 3f0)
    elseif code == "065"      # Utah Juniper (UJ) — TT uses 401 ⇒ VOLEQ(2:3)="01" (W.CO/E.UT/WY)
        v = strip(eq)[2:3] == "01" ? fpow(-0.08728f0 + 0.135420f0 * c, 3f0) :
            strip(eq)[2:3] == "02" ? fpow(-0.03655f0 + 0.135689f0 * c, 3f0) :
            strip(eq)[2:3] == "03" ? fpow( 0.04829f0 + 0.114358f0 * c, 3f0) :
                                     fpow(-0.13386f0 + 0.133726f0 * c, 3f0)
        return d < 3f0 ? 0.1f0 : v
    elseif code == "133"      # Single-leaf Pinyon Pine (PM)
        v = fpow(-0.14240f0 + 0.148190f0 * c, 3f0)
        return d < 3f0 ? 0.1f0 : v
    elseif code == "475"      # Curl-leaf Mtn-Mahogany (MC)
        return fpow(-0.13363f0 + 0.128222f0 * c, 3f0)
    elseif code == "998"      # Other Hardwoods (OH)
        return fpow(-0.13822f0 + 0.121850f0 * c, 3f0)
    end
    return 0f0
end

function compute_volumes_tt!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; c = s.control
    merch = (stmp = c.sp_stump_ht, topd = c.sp_top_diam, scfstmp = c.sp_scf_stump,
             scftop = c.sp_scf_topd, bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        # FVS volumes only trees with DBH ≥ 1 (measured: the FVStt volume loop never calls r4d2h for a
        # sub-inch record — a DBH-0.5 PM never reaches r4d2h, while DBH-2.0/2.5 do and take the 0.1 floor).
        # So the DVEW 0.1 floor (r4d2h.f:92/98, DBH<3) applies only to volumed trees, i.e. 1 ≤ DBH < 3.
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; mdl = length(strip(eq)) >= 7 ? strip(eq)[4:6] : "   "
        if mdl == "MAT"
            # TT merch cubic top = TOPD = 6.0 UNCONDITIONALLY (NOT 6·bark). TT's fvsvol.f (line 192) is a newer
            # version: "MTOPS AND TOPDIAM BOTH EQUAL THE TOPD" — no bark multiply. This DIFFERS from CI/UT whose
            # fvsvol uses MTOPS=TOPD*BARK for NFS equations (so CI/UT keep 6·bark; both bit-exact). Using 6·bark
            # here made TT small-tree MAT merch ~2-4% high (WB D=11.8: 9.70 vs live 9.30). DBHMIN=8 (sp7=7).
            mtopp = 6.0f0
            dbhmin = sp == 7 ? 7.0f0 : 8.0f0
            tcf, mcf = r4vol_volumes(eq, d, h, mtopp, 0f0)   # (total CF0, merch CFGRS) — bit-exact
            mcf = d >= dbhmin ? max(mcf, 0f0) : 0f0
            bfmind = sp == 7 ? 7.0f0 : 8.0f0                 # board DBHMIN (BFMIND, tt/grinit.f)
            bf = d >= bfmind ? r4vol_board(eq, d, h, mtopp, 0f0) : 0f0   # BFGRS Scribner (M=1)
            # Broken/killed-top reduction (r4_topkill → CFTOPK/BFTOPK): a top-killed tree over-volumes without
            # it — this WAS the TT forest-405 .sum residual (one AF idx-11 broke at 40ft: jl 53.2 vs live 44.0
            # ⇒ +55 TCuFt; TCuFt now bit-exact). Latent in CI/UT too (same MAT/FW2 path); no-op for un-killed.
            tcf, mcf, bf = r4_topkill(t, i, sp, d, h, tt_bratio(sp, d), max(tcf, 0f0), mcf, bf, merch)
            t.cuft_vol[i] = max(tcf, 0f0)
            t.merch_cuft_vol[i] = max(mcf, 0f0)
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        else
            # DVEW (PM/UJ/RM/MC/OH) — TT woodland volume = R4D2H (Chojnacky INT-339 D2H regression,
            # volume/NVEL/r4d2h.f), routed via GROSSVOL→DVEST (region 4). VOL(1)=VOL(4)=entire cubic;
            # NO Behre/CFTOPK trim in the DVE path (grossvol.f:172 + dvest.f:90). fvsvol.f NATCRS then maps
            # TCF=VOL(1), MCF=(D≥DBHMIN)·VOL(4) [VOL(7)=0], SCF=0 (region-4, not R8/R9), BdFt=0 (DVE has no
            # board-foot: fvsvol.f:411 skips DVE, BFPFLG=0 for R4). DRCOB=0 passed ⇒ D2H uses DBH not DRC;
            # FCLASS=0 ⇒ MSTEM=0 (the c-coefficient term drops). Earlier CR region-2 Chojnacky was ~5× high.
            vol1 = r4d2h_vol1(eq, d, h)
            dbhmin = 8.0f0                                   # tt/grinit.f DBHMIN(I)=8 (woodland: no override)
            t.cuft_vol[i] = max(vol1, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(vol1, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end
