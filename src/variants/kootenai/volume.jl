# =============================================================================
# volume.jl (kootenai) — KT volume (chunk 8). Region-1 Flewelling FW2 (NVEL), the SAME kernel CR uses
# (cr_fw2_vol) — KT's VOLEQ "I00FW2W<FIA>" maps through _fw2_jsp's 'I'/INGY geocode branch (already ported).
# Merch standards from kt/grinit.f: cubic TOPD=4.5/DBHMIN=7 (sp7 lodgepole=6); board BFTOPD=4.5/BFMIND=7
# (sp7=6); stump=1. KT bark = bark_ratio(KT_BKRAT) (kt/bratio.f), NOT cr_bratio.
# =============================================================================
# kt/ktfctr.f — merch-to-total cubic ratio for top-killed trees (0 ⇒ fall back to the Behre trim).
const KT_RCF1 = Float32[0.620, 1.133, 0.709, 0.592, 0.0, 0.0, 0.688, 0.0, 0.0, 1.047, 0.0]
const KT_RCF2 = Float32[3.358, 3.561, 3.475, 3.595, 0.0, 0.0, 3.580, 0.0, 0.0, 3.450, 0.0]
const KT_RCF3 = Float32[3.137, 3.418, 3.229, 3.329, 0.0, 0.0, 3.405, 0.0, 0.0, 3.290, 0.0]
kt_ktfctr(sp::Int, dtop::Float32, d::Float32)::Float32 =
    (sp <= 4 || sp == 7 || sp == 10) ? 1f0 - (KT_RCF1[sp] * fpow(dtop, KT_RCF2[sp]) / fpow(d, KT_RCF3[sp])) : 0f0

function compute_volumes_kt!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    iregn = 1                                    # KT = Northern (Region 1)
    topd = 4.5f0; bftopd = 4.5f0; stump = 1.0f0
    # vols.f:86-90 zeroes HT2TD for every record; FVSVOL/NATCRS then fills the FW2 merch-top heights.
    fill!(t.merch_top_cf, 0f0); fill!(t.merch_top_bf, 0f0)
    htb = zeros(Float32, 2)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        # ie/vols.f:144-145 (KT compiles ie/vols.f) — top-killed: H = NORMHT/100 for the full volume + CFTOPK trim.
        (h >= 4.5f0 && t.trunc[i] > 0 && t.norm_ht[i] > 0) && (h = Float32(t.norm_ht[i]) / 100f0)
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        dbhmin = sp == 7 ? 6f0 : 7f0             # kt/grinit.f (lodgepole sp7 = 6)
        bfmind = sp == 7 ? 6f0 : 7f0
        # vols.f:132,150-151: BARK=BRATIO(ISPC,DBH_start,H) before `D=D+DG(I)/BARK` ⇒ projected cycles use the stashed
        # start-of-cycle bark (t.vol_bark) for the merch tops / DBTBH / CFTOPK; grown-DBH bark at cycle 0 / dead records.
        bark = (i <= t.n && t.vol_bark[i] > 0f0) ? t.vol_bark[i] : KT_BKRAT[sp]   # kt/bratio.f BRATIO = BKRAT(IS)
        v = startswith(veq[sp], "I") ?
            cr_fw2_vol(veq[sp], d, h; bark = bark, topd = topd, bftopd = bftopd, stump = stump, iregn = iregn, sf_hs = true,
                       ht2td = htb) :
            zeros(Float32, 15)
        if startswith(veq[sp], "I")                       # fvsvol.f:337-339 cubic / :484-487 board HT1PRD → HT2TD
            d >= dbhmin && (t.merch_top_cf[i] = htb[1])
            d >= bfmind && (t.merch_top_bf[i] = htb[2])
        end
        tcf = max(v[1], 0f0)
        mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
        bf  = d >= bfmind ? v[2] : 0f0
        if t.trunc[i] > 0 && tcf > 0f0 && h >= 4.5f0     # broken-top reduction (Behre taper)
            vmax = tcf
            mcf0 = mcf
            tcf, mcf = cr_cftopk(tcf, mcf, d, h, vmax, bark, Int(t.trunc[i]), stump, topd)
            # kt/cftopk.f:52-58 — KT's own merch ratio for a top-killed tree: KTFCTR's RCF (LP/WL/DF/GF/WH/…
            # sp 1-4,7,10) replaces the Behre merch trim, MCF=RCF·TCF(trimmed). ktt01 WL/DF broken tops: MCF 7.000
            # vs live 6.767 without it.
            if mcf0 > 0f0
                rcf = kt_ktfctr(sp, topd, d)
                rcf > 0f0 && (mcf = rcf * tcf)
            end
            bf = cr_bftopk(bf, d, h, vmax, bark, Int(t.trunc[i]), stump, bftopd)
        end
        t.cuft_vol[i] = tcf; t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
    end
    return s
end
