# =============================================================================
# volume.jl (inlandempire) — IE volume (chunk 8). Region-1 NVEL. Merch standards IDENTICAL to KT
# (ie/grinit.f: TOPD=4.5/DBHMIN=7, sp7=6; BFTOPD=4.5/BFMIND=7, sp7=6; stump=1). sp1-14,23 use Flewelling
# FW2 (I00FW2W<fia>, reuse cr_fw2_vol); sp15-22 (PM/RM/PY/AS/CO/MM/PB/OH) use DVE/Behre (reuse cr_dve_vol;
# Behre sp17 deferred). VOLEQ strings DUMPED from the live binary (VOLEQDEF VAR='IE', forest 118).
# NOTE: VOLEQDEF is forest-keyed; this is the forest-118 (iet01) default — a full voleqdef port covers all forests.
# =============================================================================

# IE per-species VOLEQ (dumped from live FVSie sitset VEQNNC, forest 118).
const IE_VOL_EQ = String[
    "I00FW2W119",  # sp1
    "I00FW2W073",  # sp2
    "I00FW2W202",  # sp3
    "I00FW2W017",  # sp4
    "I00FW2W260",  # sp5
    "I00FW2W242",  # sp6
    "I00FW2W108",  # sp7
    "I00FW2W093",  # sp8
    "I00FW2W019",  # sp9
    "I00FW2W122",  # sp10
    "I00FW2W260",  # sp11
    "I00FW2W012",  # sp12
    "I00FW2W073",  # sp13
    "I00FW2W019",  # sp14
    "102DVEW106",  # sp15
    "102DVEW106",  # sp16
    "616BEHW231",  # sp17
    "102DVEW746",  # sp18
    "102DVEW740",  # sp19
    "200DVEW746",  # sp20
    "101DVEW375",  # sp21
    "200DVEW746",  # sp22
    "I00FW2W260",  # sp23
]

# ie/formcl.f COLVFC — Colville NF (IFOR=5, forest 621, a Region-6 forest) Girard form-class table,
# [IFCDBH 1..5, ISPC 1..23]. Column-major Fortran DATA transcribed to [ifcdbh, sp]. Region-1 IE forests
# ignore form class (their FW2/DVE routines don't use it) and get FC=80; Colville feeds it to the R6 Behre
# routine. FRMCLS keyword override (ie/formcl.f) not modeled — FIA stands carry none.
const IE_COLVFC = permutedims(Float32[
    78 78 78 76 76 64 80 77 78 78 75 81 75 78 56 56 56 77 76 70 70 70 75
    80 78 76 78 78 65 82 79 76 80 78 82 75 76 56 56 60 77 78 70 70 70 78
    80 80 75 77 80 66 82 80 74 80 79 82 75 74 56 56 60 77 78 70 70 70 79
    82 80 74 76 80 66 80 80 74 82 79 80 74 74 56 56 60 77 78 70 70 70 79
    80 80 74 76 82 66 80 81 74 80 78 80 74 74 56 56 60 77 78 70 70 70 78])

# ie/formcl.f FORMCL: Region-6 Colville (IFOR=5) uses COLVFC; all other IE forests default FC=80.
@inline function ie_formcl(sp::Integer, ifor::Int, d::Real)::Int
    (ifor == 5 && sp >= 1 && sp <= 23) || return 80
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1); Float32(d) > 40.9f0 && (ifcdbh = 5)
    return Int(IE_COLVFC[sp, ifcdbh])
end

# IE sp17 (PY, VEQNNC 616BEHW231) region-6 Behre volume — reuses the validated BM/SO R6 machinery
# (bm_r6vol3/bm_r6dibs/bm_r6vol1) with the IE form class. IE TOPD=BFTOPD=4.5 ⇒ MTOPP=4.5·BARK (fvsvol.f).
function ie_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32)
    fclass = ie_formcl(sp, ifor, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    vol2 = 0f0; vol4 = 0f0
    v1 = if h <= 17.3f0                                    # R6VOL short-tree guard (TTH≤FC_HT): cylinder VOL(1)
        0.00272708f0 * dbhib * dbhib * h
    else
        v = bm_r6vol3(d, dbtbh, fclass, h, 1)             # ZONE 1 total cubic → VOL(1)
        mtopp = 4.5f0 * bark                              # TOPDIAM = TOPD·BARK
        xlogs, ld1 = bm_r6dibs(d, fclass, mtopp, h)       # log bucking → small-end diameters
        lv1, lv4 = bm_r6vol1(d, fclass, xlogs, ld1)       # per-log Scribner (VOL2) + merch cubic (VOL4)
        nlog = Int(floor(xlogs)); nacc = (xlogs - nlog) > 0f0 ? nlog + 1 : nlog
        for k in 1:nacc
            vol2 += bm_anint(lv1[k]); vol4 += bm_anint(lv4[k] * 10f0) / 10f0
        end
        v
    end
    return (max(v1, 0f0), max(vol4, 0f0), max(vol2, 0f0))
end

function compute_volumes!(s::StandState, ::InlandEmpire)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq
    ifor = Int(s.plot.forest_idx)
    topd = 4.5f0; bftopd = 4.5f0; stump = 1.0f0; iregn = 1
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        dbhmin = sp == 7 ? 6f0 : 7f0
        bfmind = sp == 7 ? 6f0 : 7f0
        bark = ie_bratio(sp, d)
        eq = veq[sp]
        if occursin("BEH", eq)                            # region-6 Behre (sp17 PY, 616BEHW231)
            # Broken/dead-top trees (trunc = break-ht·100, norm_ht = predicted full ht·100): FVS builds the
            # full-bole volume at NORMHT then trims to the break via the Behre form-class taper (cftopk/bftopk),
            # exactly as the FW2 variants do. (Validated vs FVSie_g16 treelist: D=11.7 broken tree 9.16/7.44/29.5
            # vs live 9.2/7.4/29.5.)
            broken = t.trunc[i] > 0 && t.norm_ht[i] > 0
            hbase = broken ? Float32(t.norm_ht[i]) / 100f0 : h
            tcf, mcf, bf = ie_behre_vol(sp, ifor, d, hbase, bark)
            if broken && tcf > 0f0 && hbase >= 4.5f0
                vmax = tcf
                tcf, mcf = cr_cftopk(tcf, mcf, d, hbase, vmax, bark, Int(t.trunc[i]), 1f0, 4.5f0)
                bf = cr_bftopk(bf, d, hbase, vmax, bark, Int(t.trunc[i]), 1f0, 4.5f0)
            end
            t.cuft_vol[i] = max(tcf, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= bfmind ? max(bf, 0f0) : 0f0
            continue
        end
        v = if startswith(eq, "I")                       # Flewelling FW2 (sp1-14,23)
            cr_fw2_vol(eq, d, h; bark = bark, topd = topd, bftopd = bftopd, stump = stump, iregn = iregn)
        else                                              # Gevorkiantz DVE (sp15,16,18,19,20,21,22)
            cr_dve_vol(eq, d, h)
        end
        tcf = max(v[1], 0f0)
        mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
        bf  = d >= bfmind ? v[2] : 0f0
        t.cuft_vol[i] = tcf; t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
    end
    return s
end
