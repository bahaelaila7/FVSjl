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

function compute_volumes!(s::StandState, ::InlandEmpire)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq
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
        v = if startswith(eq, "I")                       # Flewelling FW2 (sp1-14,23)
            cr_fw2_vol(eq, d, h; bark = bark, topd = topd, bftopd = bftopd, stump = stump, iregn = iregn)
        elseif occursin("DVE", eq)                        # Gevorkiantz DVE (sp15,16,18,19,20,21,22)
            cr_dve_vol(eq, d, h)
        else                                              # Behre (sp17 PY) — deferred
            zeros(Float32, 15)
        end
        tcf = max(v[1], 0f0)
        mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
        bf  = d >= bfmind ? v[2] : 0f0
        t.cuft_vol[i] = tcf; t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
    end
    return s
end
