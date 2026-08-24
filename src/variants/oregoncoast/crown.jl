# =============================================================================
# crown.jl (oregoncoast) — OC crown-ratio CHANGE for the FVS-native (IORG=0) trees (oc/crown.f).
#
# The valid-ORGANON (IORG=1) trees take CRNEW = ANINT(CR2·100) from the ORGANON crown-recession
# (organon_crngro.jl), applied inside the growth hook (organon_hook.jl). The NON-ORGANON (IORG=0)
# trees — species ORGANON does not model (e.g. LP, and every FIA-surrogate species) — keep the
# FVS-native rank-based Weibull crown-change model (oc/crown.f CASE DEFAULT), and it must run EVERY
# cycle. Before this port `crown_ratio_update!(::OregonCoast)` was a no-op, so IORG=0 crown ratios
# were FROZEN at their cycle-0 inventory value; that starved the crown-vigor modifier CRMOD in the
# FVS-native height-growth kernel (oc/htgf.f: HTG = POTHTG·1.016605·(1−exp(−4.26558·CR))·RHMOD), so
# LP under-grew ~0.15 ft/cycle (oracle CR 0.35→0.53 over 10 cycles; jl stuck at 0.35). See
# scratchpad/lp_height/VALIDATION.md for the per-tree A/B vs FVSoc_clean.
#
# Coefficients (oc/crown.f CRCONS DATA) indexed by the 17 OC crown groups; species→group via IMAP.
# oc/crown.f uses LORGANON=.TRUE. for OC ⇒ the `.NOT.LORGANON`-guarded ICRI<10 bumps (crown.f:386,430)
# NEVER fire for OC, so they are omitted here. SCALE = clamp(1.5 − RELSDI, 0.30, 1.0) is the OC form
# (crown.f:312) — DISTINCT from the WestCascades/Olympic RELDEN form. ISORT is the whole-stand GROWN-
# DBH descending rank (OC grows DBH inline in the hook, so t.dbh is already the grown value here).
# =============================================================================

# oc/crown.f DATA IMAP — 50 species → 17 crown groups.
const OC_CROWN_IMAP = Int[
    6, 6, 6, 4, 9, 9, 3, 12, 12,
    13, 13, 17, 13, 13, 10, 2, 2, 10, 10,
    10, 1, 1, 1, 1, 3, 7, 7, 7, 7,
    7, 7, 7, 7, 14, 16, 15, 5, 16, 16,
    16, 16, 8, 16, 16, 16, 16, 16, 16, 16,
    1]

# oc/crown.f CRCONS 17-group coefficients (DATA WEIBA/WEIBB0/WEIBB1/WEIBC0/WEIBC1/C0/C1).
const OC_WEIBA  = Float32[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0]
const OC_WEIBB0 = Float32[0.52909,0.25115,0.52909,0.48464,0.08402,0.29964,0.06607,0.25667,0.16601,0.03685,0.25667,0.49085,0.16267,-0.81881,-1.11274,-0.23830,-0.13121]
const OC_WEIBB1 = Float32[1.00677,1.05987,1.00677,1.01272,1.10297,1.05398,1.10705,1.06474,1.08150,1.09499,1.06474,1.01414,1.07340,1.05418,1.12314,1.18016,1.15976]
const OC_WEIBC0 = Float32[-3.48211,0.33383,-3.48211,-2.78353,0.91078,-1.09270,2.04714,0.11729,0.91420,4.01340,0.11729,3.16456,3.28850,-2.36611,2.53316,3.04413,2.59824]
const OC_WEIBC1 = Float32[1.38780,0.63833,1.38780,1.27283,0.45819,0.80687,0.15070,0.61681,0.45768,0.04946,0.61681,0.00000,0.00000,1.20241,0.00000,0.00000,0.00000]
const OC_CRC0   = Float32[7.48846,6.92893,7.48846,7.44422,3.64292,5.12357,6.82187,5.95912,6.14578,6.04928,5.95912,5.48853,6.48494,4.42000,4.12048,4.62512,4.89032]
const OC_CRC1   = Float32[-0.02899,-0.04053,-0.02899,-0.04779,-0.00317,-0.01042,-0.02247,-0.01812,-0.02781,-0.01091,-0.01812,-0.00717,-0.02325,-0.01066,-0.00636,-0.01604,-0.01884]

"""
    crown_ratio_update!(s, ::OregonCoast; fint, crown_sdi) -> s

oc/crown.f CASE DEFAULT — the rank-based Weibull crown-ratio CHANGE for the FVS-native (IORG=0)
trees, run at GRADD-CROWN (simulate.jl) each cycle. IORG=1 trees are skipped (their crown ratio is
already set from the ORGANON CR2 in the growth hook). `crown_sdi` is the pre-growth Reineke SDI
(SDIBC) FVS's SDICAL feeds into RELSDI. LSTART inventory dubbing is NOT here — OC dubs missing
inventory crowns in `oc_organon_prepare!`/`oc/cratet.f` at setup, so this hook is cycling-only.
"""
function crown_ratio_update!(s::StandState, ::OregonCoast; fint::Float32 = 5.0f0,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sdiac = crown_sdi
    nsd = length(p.sp_sdi_def)
    # ISORT: whole-stand GROWN-DBH descending rank (oc/crown.f:183 — ISORT(IND(JJ))=ITRN−JJ+1, so the
    # largest tree → n, smallest → 1). OC applies DBH growth inline in the hook, so t.dbh is grown.
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n; key[i] = t.dbh[i]; idx[i] = Int32(i); end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    # IORG stashed by the growth hook (op_iorg is the shared ORGANON per-tree flag). A tree beyond the
    # stashed length (fresh regen added after the hook) is FVS-native ⇒ treated as IORG=0.
    org_ran = length(c.op_iorg) >= n
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        (org_ran && c.op_iorg[i] == 1) && continue        # ORGANON crown already set in the hook
        sp = Int(t.species[i]); (sp < 1 || sp > 50) && continue
        d = t.dbh[i]; h = t.height[i]; icr = Int(t.crown_pct[i])
        grp = OC_CROWN_IMAP[sp]
        relsdi = (sp <= nsd && p.sp_sdi_def[sp] > 0f0) ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = OC_CRC0[grp] + OC_CRC1[grp] * relsdi * 100f0
        A = OC_WEIBA[grp]
        B = OC_WEIBB0[grp] + OC_WEIBB1[grp] * acrnew; B < 3f0 && (B = 3f0)
        C = OC_WEIBC0[grp] + OC_WEIBC1[grp] * acrnew; C < 2f0 && (C = 2f0)
        scale = 1.5f0 - relsdi                            # oc/crown.f:312 (OC form, NOT RELDEN)
        scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
        x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale   # RANN path (d=0) inert
        x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
        crnew = A + B * (-flog(1f0 - x))^(1f0 / C)
        crnew *= 10f0
        if icr != 0                                       # ±1%/yr limit (crown.f:346-349), CRNMLT=1
            chg = crnew - Float32(icr)
            pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg                    # statement 41 (DLOW/DHI cover all)
        end
        icri = trunc(Int, crnew + 0.5f0)                  # 9052
        if icr != 0                                       # CRMAX cap (crown.f:368-380), cycling ICR>0
            htg = t.ht_growth[i]
            crln = h * Float32(icr) / 100f0
            crmax = (crln + htg) / (h + htg) * 100f0
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        icri > 95 && (icri = 95)                          # statement 59 (LORGANON ⇒ no <10 bump for OC)
        icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
