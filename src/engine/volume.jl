# =============================================================================
# volume.jl — per-tree volume driver (VOLS / CFVOL, the SN R8 Clark path)
#
# Ported from: base/vols.jl (VOLS) + base/fvsvol.jl (CFVOL).
#
# For each live tree, select the merchantability product (sawtimber when
# DBH ≥ SCFMIND, else pulpwood), call the pure `_R8CLARK_VOL` taper model, and
# load the four per-tree volumes the .sum reports:
#   cuft_vol  (CFV)  = vol[1]                     — total cubic
#   merch_cuft_vol (MCFV) = vol[4]+vol[7] if D≥DBHMIN — merch cubic
#   saw_cuft_vol  (SCFV) = vol[4]      if D≥SCFMIND   — sawtimber cubic
#   bdft_vol  (BFV)  = vol[10]                    — board feet
# (snt01 has no defect; defect correction + the two-pass dead handling come later.)
# =============================================================================

"""
    dub_missing_heights!(state)

HTDBH (htdbh.f, mode 0): assign a height to every tree whose input height is ≤ 0,
from the species height/DBH curve `h = 4.5 + p2·exp(-p3·D^p4)` (linear below 3").
Run before volume/growth so missing-height trees still get cubic volume (CRATET).
"""
function dub_missing_heights!(s::StandState)
    t = s.trees; sd = s.coef.species
    p2 = sd[:htdbh_p2]; p3 = sd[:htdbh_p3]; p4 = sd[:htdbh_p4]; db = sd[:htdbh_db]
    @inbounds for i in 1:t.n
        t.height[i] > 0f0 && continue
        d = t.dbh[i]; d <= 0f0 && continue
        sp = t.species[i]
        if d >= 3f0
            t.height[i] = 4.5f0 + p2[sp] * exp(-p3[sp] * d ^ p4[sp])
        else
            hat3 = 4.5f0 + p2[sp] * exp(-p3[sp] * 3f0 ^ p4[sp])
            t.height[i] = (hat3 - 4.51f0) * (d - db[sp]) / (3f0 - db[sp]) + 4.51f0
        end
    end
    return s
end

"""
    compute_volumes!(state)

Fill `trees.{cuft_vol,merch_cuft_vol,saw_cuft_vol,bdft_vol}` for every live tree
from the R8 Clark taper model and the per-species merch specs. Needs
`setup_volume_equations!` to have set `species.vol_eq`.
"""
function compute_volumes!(s::StandState)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species
    scfmin = sd[:scf_min_dbh]; scftop = sd[:scf_top_dib]; topd = sd[:top_dib]
    stmp = sd[:stump]; scfstmp = sd[:scf_stump]; dbhmin = sd[:dbh_min]
    @inbounds for i in 1:t.n
        d = t.dbh[i]; h = t.height[i]; sp = t.species[i]
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
            continue
        end
        if d >= scfmin[sp]
            prod = "01"; stump = scfstmp[sp]; mtopp = scftop[sp]
        else
            prod = "02"; stump = stmp[sp]; mtopp = topd[sp]
        end
        mtops = topd[sp]
        v, _, _ = _R8CLARK_VOL(veq[sp], d, h, mtopp, mtops, stump, prod)
        t.cuft_vol[i]       = v[1]
        t.merch_cuft_vol[i] = d >= dbhmin[sp] ? v[4] + v[7] : 0f0
        t.saw_cuft_vol[i]   = d >= scfmin[sp] ? v[4] : 0f0
        t.bdft_vol[i]       = v[10]
    end
    return s
end
