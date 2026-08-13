# =============================================================================
# organon_volume.jl — OC (Oregon Coast) BLM cubic volume (chunk C10a).
#
# Ported from the NVEL BLM volume library (volume/NVEL/blmvol.f + blmtap.f) + oc/formcl.f. OC's
# VOLEQDEF equations are BLM Behre (`B00/B01 BEHW`+FIA, oracle VOLEQ table); `volinit.f:366` routes
# `VOLEQ(1:1)=='B'` to BLMVOL. The TOTAL cubic volume is Smalian's over 4-ft sections using the BLM
# Behre-hyperbola stem taper (`BLMTAP` feet path); small trees use the closed-form `0.00272708·DIB²·H`.
#   • oc_formcl        — oc/formcl.f FORMCL (Medford IFOR=9 ⇒ BLM711 per-species form class, INT-truncated)
#   • oc_blmtapeq      — blmvol.f BLMTAPEQ (VOLEQ → stem PROFILE 1..10 + TAPEQU)
#   • oc_double_bark   — blmvol.f DOUBLE_BARK (DBHOB → DBHIB, per TAPEQU)
#   • oc_blmtap        — blmtap.f BLMTAP (feet path: DIB at height via Behre hyperbola)
#   • oc_blmtcub       — blmvol.f BLMTCUB (Smalian total cubic)
#   • oc_tree_cuft     — blmvol.f BLMVOL total-cubic driver (small-tree closed form vs BLMTCUB)
#
# MEASURED vs FVSoc_clean TREELIST (ocmin cyc0): per-tree total cuft (e.g. DF D12.7/H67 → 30.1,
# GF D6.2/H38 → 6.7, LP D11.5/H73 → 16.9). Board-foot (SCRIB) is C10b. See docs/OC_VARIANT_PORT_AUDIT.md.
# =============================================================================

# oracle VOLEQDEF table (VAR='OC', REGN 7, FORST 11) — per FVS species index 1..50.
const OC_VOLEQ = String[
    "B00BEHW081","B00BEHW081","B00BEHW242","B00BEHW017","B00BEHW021","B00BEHW021","B01BEHW202",
    "B00BEHW263","B00BEHW260","B00BEHW119","B00BEHW108","B00BEHW108","B00BEHW108","B00BEHW108",
    "B00BEHW116","B00BEHW117","B00BEHW119","B00BEHW122","B00BEHW108","B00BEHW108","B00BEHW242",
    "B00BEHW093","B00BEHW211","B00BEHW231","B00BEHW999","B00BEHW800","B00BEHW800","B00BEHW800",
    "B00BEHW800","B00BEHW800","B00BEHW800","B00BEHW800","B00BEHW800","B00BEHW312","B00BEHW800",
    "B00BEHW351","B00BEHW361","B00BEHW431","B00BEHW999","B00BEHW312","B00BEHW999","B00BEHW631",
    "B00BEHW800","B00BEHW999","B00BEHW747","B00BEHW999","B00BEHW231","B00BEHW631","B00BEHW999",
    "B00BEHW211"]

# oc/formcl.f BLM711 (Medford) per-species form class.
const OC_FORMCL_BLM711 = Float32[
    70.,78.,78.,91.,78.,78.,87.,91.,70.,73.,68.,68.,68.,68.,70.,84.,76.,85.,68.,68.,
    70.,74.,70.,76.,70.,80.,80.,80.,80.,80.,80.,80.,80.,84.,80.,88.,81.,83.,83.,84.,
    83.,84.,80.,72.,72.,75.,76.,84.,70.,75.]

# blmtap.f BLMTHT(4,10): Behre-hyperbola coefficients B0,B1,B2,B3 per stem PROFILE 1..10.
const OC_BLMTHT_B0 = Float32[0.6448,0.6096,0.31385,0.4779,0.5455,0.45648,0.6014,0.54568,0.4606,0.6200]
const OC_BLMTHT_B1 = Float32[-0.00196,-0.00196,0.0,0.0,-0.00196,0.00289,0.0,0.0,0.0,0.0]
const OC_BLMTHT_B2 = Float32[0.0,0.0,0.002985,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OC_BLMTHT_B3 = Float32[0.0,0.0,-0.00003386,0.0,0.0,0.0,0.0,0.00000546,0.0,0.0]

"oc/formcl.f FORMCL (IFOR=9 Medford): per-species BLM711 form class, INT-truncated."
@inline oc_formcl(sp::Int) = trunc(Int, OC_FORMCL_BLM711[sp])

"blmvol.f BLMTAPEQ — VOLEQ → stem PROFILE (1..10). Uses FIA=VOLEQ[8:10], geosub=VOLEQ[1:3]."
function oc_blmtapeq(veq::AbstractString)
    fia = veq[8:10]; geo = veq[1:3]
    if fia == "202"
        geo == "B01" && return 1
        geo == "B02" && return 2
        return 10
    end
    fia == "211" && return 10
    fia == "122" && return 3
    fia == "116" && return 10
    fia == "117" && return 4
    fia == "119" && return 5
    fia == "108" && return 10
    fia == "231" && return 10
    fia == "631" && return 10
    fia == "351" && return 10
    fia == "998" && return 10
    fia == "312" && return 10
    fia == "361" && return 10
    fia == "431" && return 10
    fia == "542" && return 10
    fia == "747" && return 10
    fia == "800" && return 10
    fia == "015" && return (geo == "B01" ? 6 : 10)
    fia == "021" && return 7
    fia == "017" && return 6
    fia == "011" && return 7
    fia == "022" && return 7
    fia == "093" && return 10
    fia == "098" && return 10
    (fia == "260" || fia == "263") && return 8
    fia == "081" && return 9
    fia == "042" && return 10
    fia == "041" && return 10
    fia == "242" && return 9
    fia == "073" && return 9
    return 10
end

"blmvol.f BLMTAPEQ — VOLEQ → TAPEQU (for DOUBLE_BARK)."
function oc_blmtapeq_tapequ(veq::AbstractString)
    fia = veq[8:10]; geo = veq[1:3]
    fia == "202" && return (geo == "B01" ? 1 : geo == "B02" ? 2 : geo == "B04" ? 4 : geo == "B05" ? 6 : 3)
    fia == "211" && return 5
    fia == "122" && return (geo == "B01" ? 10 : 11)
    fia == "116" && return 12
    fia == "117" && return 13
    fia == "119" && return 14
    fia == "108" && return 15
    fia == "231" && return 20
    fia == "631" && return 21
    fia == "351" && return 22
    fia == "998" && return 23
    fia == "312" && return 24
    fia == "361" && return 25
    fia == "431" && return 26
    fia == "542" && return 27
    fia == "747" && return 28
    fia == "800" && return 29
    fia == "015" && return (geo == "B01" ? 30 : 31)
    fia == "021" && return 32
    fia == "017" && return 33
    fia == "011" && return 34
    fia == "022" && return 35
    fia == "093" && return 41
    fia == "098" && return 42
    (fia == "260" || fia == "263") && return 48
    fia == "081" && return 51
    fia == "042" && return 52
    fia == "041" && return 53
    fia == "242" && return 54
    fia == "073" && return 55
    return 56
end

"blmvol.f DOUBLE_BARK — DBHOB → DBHIB (large-end butt-log DIB) per TAPEQU."
function oc_double_bark(tapequ::Int, d::Float32)
    if tapequ in (1,2,3,5,35)
        return 0.903563f0*fpow(d, 0.989388f0)
    elseif tapequ == 11 || tapequ == 12
        return 0.809427f0*fpow(d, 1.016866f0)
    elseif tapequ == 13 || tapequ == 14
        return 0.859045f0*d
    elseif tapequ == 15
        return d - (0.3147f0 + 0.0274f0*d)
    elseif tapequ == 20 || tapequ == 25
        return -0.03425f0 + 0.98155f0*d
    elseif tapequ == 21
        return -4.36852f0 + 0.95354f0*d + 0.18307f0*4.5f0
    elseif tapequ in (22,23,24,26,27)
        return 0.39534f0 + 0.90182f0*d
    elseif tapequ == 28 || tapequ == 29
        return -0.78034f0 + 0.95956f0*d
    elseif tapequ == 31 || tapequ == 33
        return 0.904973f0*d
    elseif tapequ == 32 || tapequ == 34
        return 0.86951f0*fpow(d, 1.00983f0)
    elseif tapequ == 41 || tapequ == 42
        return d - (0.2113f0 + 0.0445f0*d)
    elseif tapequ == 48 || tapequ == 56
        return d/1.071f0
    elseif tapequ == 52 || tapequ == 54
        return d/1.053f0
    elseif tapequ == 51 || tapequ == 53
        return 0.837291f0*d
    elseif tapequ == 55
        return d - (0.1231f0 + 0.1306f0*d)
    else
        return d
    end
end

"blmtap.f BLMTAP (feet path, TLH=0) — DIB at height `htup` for stem `profile`, DBH `dbhob`, total ht `tth`, D17."
@inline function oc_blmtap(profile::Int, dbhob::Float32, tth::Float32, htup::Float32, d17::Float32, xlen::Float32)
    hbutt = tth - (xlen + 1.5f0)
    hbutt <= 0f0 && return 0f0
    htdib = tth - htup
    a = OC_BLMTHT_B0[profile] + OC_BLMTHT_B1[profile]*dbhob + OC_BLMTHT_B2[profile]*tth +
        OC_BLMTHT_B3[profile]*dbhob*tth
    b = 1f0 - a
    r = htdib/hbutt
    return d17*(r/(a*r + b))
end

"blmvol.f BLMTCUB — Smalian total cubic over 4-ft sections (stump cylinder + tip)."
function oc_blmtcub(profile::Int, dbhob::Float32, tth::Float32, d17::Float32, htlog::Float32)
    tth <= 0f0 && return 0f0
    htloop = trunc(Int, (tth + 0.5f0 - 1.0f0)/4.0f0)
    hgt2 = 1.0f0
    d2 = oc_blmtap(profile, dbhob, tth, hgt2, d17, htlog)
    r = d2/2.0f0
    tcvol = (3.1416f0*r*r)/144.0f0                    # 1-ft stump cylinder
    @inbounds for _ in 1:htloop
        d2old = d2
        hgt2 += 4.0f0
        d2 = oc_blmtap(profile, dbhob, tth, hgt2, d17, htlog)
        tcvol += 0.00272708f0*(d2old*d2old + d2*d2)*4.0f0     # Smalian
    end
    (tth - hgt2) > 0f0 && (tcvol += 0.00272708f0*(d2*d2)*(tth - hgt2))   # tip
    return tcvol
end

"""
    oc_tree_cuft(sp, dbh, ht) -> total cubic feet

blmvol.f BLMVOL total-cubic driver for OC. Small trees (total ht ≤ 17.8 or SMD_17 < MTOPP) use the
closed form `0.00272708·DBHIB²·TTH`; larger trees use `BLMTCUB` (Behre-taper Smalian). TTH = HT+1.5.
"""
function oc_tree_cuft(sp::Int, dbh::Float32, ht::Float32)
    dbh <= 0f0 && return 0f0
    veq = OC_VOLEQ[sp]
    profile = oc_blmtapeq(veq)
    tapequ = oc_blmtapeq_tapequ(veq)
    dbhib = oc_double_bark(tapequ, dbh)
    dbhib <= 0.0001f0 && return 0f0
    mtopp = 4.5f0 * oc_bratio(sp, dbh)                # fvsvol MTOPP = TOPD(=4.5, sitset.f:244 IFOR 6-10) · BARK
    tth = ht + 1.5f0
    fclass = Float32(oc_formcl(sp))
    if tth > 0f0
        if tth <= 17.8f0
            return 0.00272708f0*(dbhib*dbhib)*tth
        end
        smd_17 = trunc(sqrt(dbhib*dbhib - (dbhib*dbhib)*17.3f0/tth) + 0.5f0)
        if smd_17 < mtopp
            return 0.00272708f0*(dbhib*dbhib)*tth
        end
    end
    d17 = round((dbh*fclass)/100.0f0, RoundNearestTiesAway)
    return oc_blmtcub(profile, dbh, tth, d17, 16.3f0)
end

"""
    compute_volumes_oc!(s)

OC total cubic volume (BLM Behre taper). Fills `t.cuft_vol` = `oc_tree_cuft` (total cuft) for every
live + dead record. Merch-cubic and board-foot (SCRIB) are chunk C10b; the `.sum` TCuFt aggregates
`cuft_vol·tpa/GROSPC` in the shared reporter.
"""
function compute_volumes_oc!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees
    @inbounds for i in 1:(t.n + t.ndead)
        sp = Int(t.species[i])
        v = (1 <= sp <= 50) ? oc_tree_cuft(sp, t.dbh[i], t.height[i]) : 0f0
        t.cuft_vol[i] = v
        t.merch_cuft_vol[i] = 0f0   # C10b
        t.bdft_vol[i] = 0f0         # C10b
        t.saw_cuft_vol[i] = 0f0
    end
    return s
end
