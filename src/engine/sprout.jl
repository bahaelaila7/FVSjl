# =============================================================================
# sprout.jl — stump-sprout regeneration sub-routines (NSPREC / SPRTHT / ESSPRT)
#
# Ported from: bin/FVSsn_buildDir/essprt.f (SN-variant SELECT CASE blocks).
# These are the three pure helpers consumed by ESUCKR (esuckr.f) when a harvested
# stump regenerates as stump/root sprouts:
#
#   nsprec_sn  — number of sprouts produced per stump   (species + stump DBH)
#   sprtht_sn  — sprout height at a given sprout age     (species + site index)
#   essprt_sn  — per-record survival multiplier on the carried sprout TPA (PREM)
#
# Each is bit-faithful to the Fortran SN block. The ESSPRT per-species coefficient
# blob lives in data/southern/sprout_essprt.csv (loaded as the per-species columns
# essprt_kind / essprt_p1 / essprt_p2 / essprt_fsp). NSPREC and SPRTHT are tiny
# piecewise rules, so they stay inline rather than in a CSV.
#
# Status: pure functions, unit-tested in isolation. Wired into the ESUCKR
# generation loop in Chunk C (until then they have no .sum effect).
# =============================================================================

"""
    nsprec_sn(ispc, dstmp) -> Int

Number of sprouts produced by one cut stump (`NSPREC`, essprt.f:1102-1120, SN).
`dstmp` is the stump diameter (in). Most species yield a single sprout; species
5 yields one only below 7 in, and the oak/sweetgum group {33,61,80,82} ramps
1→3 over the 5–10 in stump-DBH range.
"""
function nsprec_sn(ispc::Integer, dstmp::Float32)::Int
    if ispc == 5
        return dstmp < 7f0 ? 1 : 0
    elseif ispc == 33 || ispc == 61 || ispc == 80 || ispc == 82
        if dstmp < 5f0
            return 1
        elseif dstmp <= 10f0                       # 5.0 ≤ DSTMP ≤ 10.0
            return Int(nint(-1f0 + 0.4f0 * dstmp)) # NINT (ties away from zero)
        else
            return 3
        end
    else
        return 1
    end
end

"SPRTHT SN species set that uses the `(0.1 + SI/50)·age` curve (essprt.f:1389)."
@inline _sprtht_sn_curve(ispc::Integer) =
    ispc == 5 || ispc == 15 || ispc == 16 ||
    (18 <= ispc <= 57) || (59 <= ispc <= 87)

"""
    sprtht_sn(ispc, si, iag) -> Float32

Sprout height (ft) at sprout age `iag` for site index `si` (`SPRTHT`,
essprt.f:1387-1393, SN). The sprouting hardwoods use `(0.1 + SI/50)·age`;
everything else falls back to the original NI regen rule `0.5 + 0.5·age`.
"""
@inline function sprtht_sn(ispc::Integer, si::Float32, iag::Real)::Float32
    a = Float32(iag)
    return _sprtht_sn_curve(ispc) ? (0.1f0 + si / 50f0) * a : 0.5f0 + 0.5f0 * a
end

"""
    sprout_dbh(coef, ispc, ht) -> Float32

Sprout DBH (in) from sprout height `ht` (ft) via the ESUCKR Wykoff height–diameter
inverse (esuckr.f:296-307): `DBH = HT2/(ln(HT−4.5) − AX) − 1`, floored at 0.1.

`AX = HT1` (the default coefficient) in the standard case: `IABFLG` defaults 1
(grinit.f:105) and CRATET only re-fits it (`IABFLG=0`, `AX=AA`) when the species'
height–diameter regression is enabled — `LHTDRG`, off by default (grinit.f:104),
flipped only by an explicit HT-DBH keyword. So for an ordinary SPROUT stand the
per-stand AA fit never runs and the Wykoff defaults `HT1`/`HT2` (sitset.f, in
`sprout_htdbh_wykoff.csv`) apply directly. Heights ≤ 4.5 ft get the 0.1-in floor.
"""
@inline function sprout_dbh(coef::SpeciesCoefficients, ispc::Integer, ht::Float32)::Float32
    ht > 4.5f0 || return 0.1f0
    ax = coef_col(coef, :wykoff_ht1)[ispc]
    bx = coef_col(coef, :wykoff_ht2)[ispc]
    d = bx / (log(ht - 4.5f0) - ax) - 1f0
    return d < 0.1f0 ? 0.1f0 : d
end

"""
    sprtht_ne(ispc, si, iag) -> Float32

NE sprout height (essprt.f ENTRY SPRTHT, `CASE('NE')`, lines 1302-1308): the sprouting
hardwoods {26-70, 72-97, 99-108} use the site-driven curve `(0.1 + SI/50)·age`; every
other (sproutable) species uses the flat `0.5 + 0.5·age`. (Same formula as SN's branch,
different species set.) `iag` = sprout age (ISHAG).
"""
@inline function sprtht_ne(ispc::Integer, si::Float32, iag::Real)::Float32
    a = Float32(iag)
    tall = (26 <= ispc <= 70) || (72 <= ispc <= 97) || (99 <= ispc <= 108)
    return tall ? (0.1f0 + si / 50f0) * a : 0.5f0 + 0.5f0 * a
end

"""
    ne_sprout_dbh(coef, ispc, ht) -> Float32

NE sprout DBH (esuckr.f:294-301) — the same Wykoff-form inverse as `sprout_dbh` but reading
NE's `:htdbh_ht1`/`:htdbh_ht2` columns: `DBH = HT2/(ln(HT−4.5) − HT1) − 1`, floored 0.1, with
HT≤4.5 ⇒ 0.1. This is the IABFLG=1 path (the default; the IABFLG=0/AA branch needs the LHTDRG
per-stand HT-DBH re-fit, which a no-measured-height SPROUT stand never triggers — cratet.f:311).
"""
@inline function ne_sprout_dbh(coef::SpeciesCoefficients, ispc::Integer, ht::Float32)::Float32
    ht > 4.5f0 || return 0.1f0
    ax = coef_col(coef, :htdbh_ht1)[ispc]
    bx = coef_col(coef, :htdbh_ht2)[ispc]
    d = bx / (log(ht - 4.5f0) - ax) - 1f0
    return d < 0.1f0 ? 0.1f0 : d
end

"""
    nsprec_ne(issp, dstmp) -> Int

NSPREC NE branch (essprt.f:1006, `CASE('NE')`): number of sprout records produced by one cut stump,
by species and stump DBH (DSTMP). NINT = round-half-away (all branch values are ≥0 ⇒ `floor(x+0.5)`).
"""
@inline function nsprec_ne(issp::Integer, dstmp::Float32)::Int
    nint(x) = floor(Int, x + 0.5f0)
    if issp == 49
        return 2
    elseif issp == 51
        return dstmp < 25f0 ? 1 : 0
    elseif issp == 53
        return dstmp < 12f0 ? 1 : 0
    elseif issp in (72, 73, 75, 78)
        return dstmp < 8f0 ? 1 : 0
    elseif issp in (26, 27, 28, 29, 43, 45, 59, 60, 61, 67, 68, 70, 86, 90, 102, 104)
        return dstmp < 5f0 ? 1 : (dstmp <= 10f0 ? nint(0.2f0 * dstmp) : 2)
    elseif issp in (40, 46, 50, 52, 82, 87, 92, 93, 94, 101)
        return dstmp < 5f0 ? 1 : (dstmp <= 10f0 ? nint(-1f0 + 0.4f0 * dstmp) : 3)
    else
        return 1
    end
end

"""
    essprt_ne(issp, prem, dstmp) -> Float32

ESSPRT NE branch (essprt.f:362, `CASE('NE')`): the per-stump sprout survival multiplier applied to the
removed TPA (PREM), by species and stump DBH. Logistic forms kept in the exact FVS expression
(`exp(z)/(1+exp(z))` or `1/(1+exp(−w))`) for bit-exactness.
"""
@inline function essprt_ne(issp::Integer, prem::Float32, dstmp::Float32)::Float32
    logi(a, b, x) = (e = exp(a + b * x); e / (1f0 + e))
    if issp in (26, 29, 41, 42, 43, 44, 45, 46, 54)
        return prem * (dstmp < 12f0 ? 0.80f0 : 0.50f0)
    elseif issp == 27 || issp == 28
        return dstmp < 34.1f0 ? prem * ((89.191f0 - 2.611f0 * dstmp) / 100f0) : 0f0
    elseif issp in (30, 78, 100)
        return prem * 0.3f0
    elseif issp in (31, 32, 34, 47, 48, 57, 62, 76, 88, 95, 96, 97, 102, 103, 105, 107, 108)
        return prem * 0.70f0
    elseif issp in (33, 84, 85, 86, 90, 91)
        return prem * 0.90f0
    elseif issp in (35, 38, 39)
        return prem * (dstmp < 24f0 ? 0.95f0 : 0.60f0)
    elseif issp in (36, 37)
        return prem * (dstmp < 24f0 ? 0.75f0 : 0.50f0)
    elseif issp == 40
        return prem * 0.93f0
    elseif issp == 50
        return prem * (dstmp < 25f0 ? 0.80f0 : 0.50f0)
    elseif issp == 51 || issp == 53
        return prem * 0.40f0
    elseif issp in (52, 56, 63, 80, 82, 87, 92, 93, 94, 101, 106)
        return prem * 0.80f0
    elseif issp == 55
        return prem * logi(1.6134f0, -0.0184f0, ((dstmp / 0.7788f0) - 0.4403f0) * 2.54f0)
    elseif issp == 58
        return prem * (1f0 / (1f0 + exp(-(-2.8058f0 + 22.6839f0 * (1f0 / ((dstmp / 0.7788f0) - 0.4403f0))))))
    elseif issp in (59, 60, 61, 67, 70)
        return prem * ((57.3f0 - 0.0032f0 * dstmp^3) / 100f0)
    elseif issp == 64 || issp == 66
        return prem * logi(6.4205f0, -0.1097f0, ((dstmp / 0.8188f0) - 0.23065f0) * 2.54f0)
    elseif issp == 68 || issp == 89
        return prem * (dstmp < 10f0 ? 0.80f0 : 0.50f0)
    elseif issp == 69
        return prem * logi(6.0065f0, -0.0777f0, (dstmp / 0.7801f0) * 2.54f0)
    elseif issp in (72, 73, 75)
        return prem * 0.40f0
    elseif issp in (74, 77, 83)
        return prem * 0.50f0
    elseif issp == 79
        return prem * (dstmp < 8f0 ? 0.80f0 : 0.50f0)
    elseif issp == 81
        return prem * (1f0 / (1f0 + exp(-(2.7386f0 + (-0.1076f0 * dstmp)))))
    elseif issp == 99
        return prem * (dstmp < 15f0 ? 0.60f0 : 0.30f0)
    elseif issp == 104
        return prem * (dstmp < 8f0 ? 0.70f0 : 0.90f0)
    else
        return prem
    end
end

"Special-establishment forests (R8/R9 NFs) that trigger the ESSPRT overrides."
@inline _es_special_forest(isefor::Integer) =
    isefor == 809 || isefor == 810 || isefor == 905 || isefor == 908

"""
    essprt_sn(coef, ispc, prem, dstmp, isefor) -> Float32

Apply the per-record sprout-survival multiplier to `prem` (carried sprout TPA)
for species `ispc` and stump diameter `dstmp` (`ESSPRT`, essprt.f:514-590, SN).

Most species use either a constant multiplier or a logistic in stump DBH,
`1/(1 + exp(-(a + b·DSTMP)))`; both forms (plus a per-species flag) are read
from `sprout_essprt.csv`. Five species (64/66/70/75/77) carry a distinct
special-forest variant (forests 809/810/905/908) handled explicitly here; for
all other forests their CSV row already holds the common-forest (ELSE) form.
"""
function essprt_sn(coef::SpeciesCoefficients, ispc::Integer, prem::Float32,
                   dstmp::Float32, isefor::Integer)::Float32
    if coef_col(coef, :essprt_fsp)[ispc] == 1f0 && _es_special_forest(isefor)
        d = dstmp
        m = if ispc == 64 || ispc == 66 || ispc == 75
                (57.3f0 - 0.0032f0 * d^3) / 100f0          # essprt.f:547/554/571
            elseif ispc == 70
                1f0 / (1f0 + exp(-(2.3656f0 - 0.2781f0 * (d / 0.7801f0))))  # :561
            else # ispc == 77
                1f0 / (1f0 + exp(-(-2.8058f0 + 22.6839f0 *
                                    (1f0 / ((d / 0.7788f0) - 0.4403f0)))))  # :578
            end
        return prem * Float32(m)
    end
    kind = coef_col(coef, :essprt_kind)[ispc]
    p1 = coef_col(coef, :essprt_p1)[ispc]
    if kind == 0f0
        return prem * p1                                    # constant multiplier
    end
    p2 = coef_col(coef, :essprt_p2)[ispc]
    return prem * (1f0 / (1f0 + exp(-(p1 + p2 * dstmp))))   # logistic in DSTMP
end

"""
    esuckr!(s; fint) -> Bool

Stump/root-sprout regeneration (ESUCKR, esuckr.f:156-349). For each cut record logged
by `_log_cut!` this cycle (a sprouting species, in removal order) create its sprouts:
`nsprec_sn` gives the number of sprout records, `essprt_sn` the survival-adjusted
carried TPA, `sprtht_sn` the height (× HMULT, plus a clamped `BACHLO(0,0.5)` `:estab`
deviation scaled by HT/5.5), and `sprout_dbh` the Wykoff DBH. New records carry
`IMC=2`, `ICR=70`, `ABIRTH=ISHAG`, and `PROB = PREM·SMULT`.

Runs in the ESNUTR phase (GRADD order, esnutr.f:113) *before* `establish!` — both draw
the `:estab` stream, and ESUCKR consumes it first. Gated on LSPRUT; a no-op when
sprouting is off or nothing was cut, so the default (LSPRUT off) path is untouched.
SMULT/HMULT use the stand-level SPROUT keyword values (the per-species/DBH-range OPGET
table, esuckr.f:96-150, is a later refinement — the common SPROUT form applies one pair
to all sprouting species).
"""
function esuckr!(s::StandState; fint::Float32 = 5f0)::Bool
    (s.control.lsprut && !isempty(s.control.cut_log)) || return false
    ovr = s.control.sprout_overrides                   # per-species SPROUT activity-450 table (esuckr.f:96-205)
    t = s.trees; coef = s.coef
    # IFORDI = KODFOR÷100 (forkod.f:183) — the 3-digit district index the special-forest gate (809/810/905/908)
    # keys on, NOT the 5-digit KODFOR. (forkod's SELECT CASE remaps a few non-canonical aliases, e.g. 7207→809;
    # those rare cases would need the full forkod port — see INDEX.) snt01=80106 ⇒ IFORDI=801 (not special).
    isefor = Int(s.plot.user_forest_code) ÷ 100
    icyc = Int(s.control.cycle)
    ne = s.variant isa Northeast     # NE ESUCKR (1 record/stump, SPRTHT/Wykoff DBH) vs SN ESSPRT model
    created = false
    @inbounds for rec in s.control.cut_log
        prem = rec.prem
        prem < 0.001f0 && continue                     # esuckr.f:170
        issp = Int(rec.species); dstmp = rec.dstmp
        ishag = Int(rec.ishag); iplot = Int(rec.plot)
        # SPROUT keyword multipliers, looked up by the PARENT species + stump DBH (esuckr.f:197-205, DO 450):
        # default 1/1, then each matching activity (species match AND DSTMP ∈ [dmin,dmax)) overwrites; last wins.
        smult = 1f0; hmult = 1f0
        @inbounds for (code, sm, hm, dmn, dmx) in ovr
            matched = code > 0f0 ? (Int(code) == issp) :
                      code == 0f0 ? true :
                      (let g = -Int(code); 1 <= g <= length(s.control.sp_groups) && (Int32(issp) in s.control.sp_groups[g]) end)
            (matched && dstmp >= dmn && dstmp < dmx) && (smult = sm; hmult = hm)
        end
        smult <= 0f0 && continue                       # esuckr.f:211 — SMULT≤0 ⇒ no sprouts for this stump
        # NE and SN share the ESUCKR structure: NUMSPR=NSPREC records, each PREM reduced by ESSPRT survival
        # (both VARACD-branched in essprt.f). NE uses its own CASE('NE') tables (nsprec_ne / essprt_ne); SN
        # uses nsprec_sn / essprt_sn. ⚠ NE aspen suckering (ESASID(NE)=49 → ASSPTN) is still TODO — for a cut
        # sp49 record NE would call ASSPTN to reset PREM before ESSPRT; absent it, sp49 uses the plain PREM.
        numspr = ne ? nsprec_ne(issp, dstmp) : nsprec_sn(issp, dstmp)
        prem = ne ? essprt_ne(issp, prem, dstmp) : essprt_sn(coef, issp, prem, dstmp, isefor)
        prem < 0.001f0 && continue                     # esuckr.f:170/244
        si = s.plot.sp_site_index[issp]                # SITEAR(ISSP)
        sp2 = s.species.class_codes[issp, 1][1:2]      # 2-char alpha code (for CWCALC)
        prob = prem * smult                            # PROB(ITRN)
        for _ in 1:numspr
            n = t.n + 1; n > length(t.dbh) && break    # no ESCPRS compression — list-overflow guard
            t.n = n
            # height: SPRTHT × HMULT + clamped BACHLO(0,0.5,ESRANN) deviation (× HT/5.5). NE & SN share the
            # SPRTHT formula but differ in the per-variant sprouting-species set (and the DBH coef columns).
            ht = (ne ? sprtht_ne(issp, si, ishag) : sprtht_sn(issp, si, ishag)) * hmult
            randev = 0f0
            while true
                randev = bachlo(s.rng, 0f0, 0.5f0; stream = :estab)
                -1f0 <= randev <= 1f0 && break
            end
            ht += randev * ht / 5.5f0
            dbh = ne ? ne_sprout_dbh(coef, issp, ht) : sprout_dbh(coef, issp, ht)
            # CWCALC's CR arg is the dummy CRDUM=1.0 (esuckr.f:317), NOT the record's ICR=70 (that is the
            # discarded 6th arg IICR, cwcalc.f). Passing 70 inflated sprout CrWidth by cr_coef·69 for Bechtold spp.
            cw = crown_width(coef, sp2, dbh, ht, 1f0, 0,
                             s.plot.latitude, s.plot.longitude, s.plot.elevation)
            # tree-record initialization (esuckr.f:258-343)
            t.mort_code[n]   = Int32(2)                # IMC = 2 (sprout regeneration)
            t.species[n]     = Int32(issp)
            t.plot_id[n]     = Int32(iplot)
            t.tpa[n]         = prob
            t.dbh[n]         = dbh
            t.height[n]      = ht
            t.crown_pct[n]   = Int32(70)               # ICR
            t.crown_ratio[n] = 70f0                    # regen convention (see establish!)
            t.crown_width[n] = cw
            t.norm_ht[n]     = Int32(0)
            t.trunc[n]       = Int32(0)
            t.birth_age[n]   = Float32(ishag)          # ABIRTH
            t.tree_id[n]     = Int32(10000000 + icyc * 10000 + n)  # IDTREE
            t.tree_random[n] = -999f0                  # ZRAND
            t.sort_key[n]    = Float64(n)
            # zero the carried-over fields (the slot may hold a previously-deleted record)
            t.diam_growth[n] = 0f0; t.ht_growth[n] = 0f0
            t.cuft_vol[n]    = 0f0; t.merch_cuft_vol[n] = 0f0
            t.saw_cuft_vol[n] = 0f0; t.bdft_vol[n] = 0f0
            t.defect[n]      = Int32(0); t.special[n] = Int32(0)
            t.cull[n]        = 0f0; t.decay_code[n] = Int32(0); t.woodland_stems[n] = Int32(0)
            t.old_random[n]  = 0f0; t.old_crown_pct[n] = 0f0
            created = true
        end
    end
    created && compute_density!(s)
    return created
end
