# =============================================================================
# site_index.jl (centralidaho) — CI habtyp/forkod/sitset (chunk 2). Traced ci/habtyp.f + ci/sitset.f.
#
# CI has TWO habitat indices (ci/habtyp.f): ICINDX (1..130, CI-type, bracket ICITYP) and ITYPE
# (1..30, NI, = NIHMAP(I-1)). The DG chunk reads ICHBCL(ICINDX,sp) (dgf.f:546) ⇒ we stash ICINDX
# in p.habitat_input (what dgf! will read). SITEAR default = site-RANGE interpolation (SITELO/HI in
# species_coefficients.csv) — NO MAPSIT (unlike IE). SDImax forest-dependent: IFOR<2 (forest 117)
# = Stage BAMAX-derived; IFOR≥2 = ZEIDE SDIDEF(sp)=R4SDI(sp). cit01 forest 412→IFOR 4 ⇒ R4SDI.
# =============================================================================

# ci/blkdat.f ICITYP(130): CI habitat-type code brackets.
const CI_ICITYP = Int32[50,60,70,80,100,120,130,140,160,161,162,170,190,195,200,210,220,221,222,250,
    260,262,264,265,280,290,310,313,315,320,323,324,325,330,331,332,334,340,341,343,344,360,370,371,
    372,375,380,385,390,392,393,395,396,397,398,400,410,440,490,493,500,505,510,511,515,520,525,526,
    527,580,585,590,591,592,593,600,605,620,621,625,635,636,637,638,640,645,650,651,652,654,655,660,
    661,662,663,670,671,672,690,691,692,694,700,705,720,721,723,730,731,732,734,740,745,750,780,790,
    791,793,810,830,831,833,850,870,900,905,920,940,955,999]
# ci/habtyp.f NIHMAP(130): CI-type (ICINDX) → NI 30 habitat type (ITYPE). Expanded from repeat-syntax.
const CI_NIHMAP = Int32[
    fill(1,11); fill(2,3); fill(3,3); fill(2,2); 3; fill(4,4); 5; 6; fill(7,3); fill(8,4);
    fill(9,24); 10; fill(11,2); fill(12,5); 13; 12; fill(13,2); fill(17,6); fill(18,2); fill(19,7);
    fill(20,7); fill(21,4); fill(22,3); fill(23,3); fill(24,2); fill(25,2); fill(26,2); fill(27,11);
    fill(28,4); fill(29,7); 30]
# ci/sitset.f BAMAXA(130): BA-max by ICINDX (used only when IFOR<2).
const CI_BAMAXA = Float32[fill(100,7); 130; fill(95,3); 220; 213; 140; 175; 154; 137; 122; 152; 100;
    fill(220,2); 233; 220; 175; 250; fill(190,2); 226; 140; 200; 210; 140; 160; 200; 160; 200; 160;
    180; 160; fill(200,4); 160; 225; fill(160,2); fill(170,2); 215; 210; 225; 210; 220; 300; 225; 350;
    280; 225; 280; 342; 200; 230; 282; 350; 280; 300; 280; 150; 253; fill(242,2); 250; 275; 225; 300;
    220; 270; 220; fill(225,2); fill(250,2); 150; 230; 140; 175; 240; 140; 210; 200; 250; 230; 200;
    170; 200; 170; 160; 230; 220; 160; fill(200,2); fill(180,2); 190; 120; 140; 156; 120; 210; 140;
    160; 214; 130; 154; 120; 300; fill(130,2); 200; 100; 75; fill(100,2); 110; fill(150,3)]
# ci/sitset.f R4SDI(MAXSP): per-species ZEIDE SDImax (IFOR≥2).
const CI_R4SDI = Float32[529,423,570,562,682,762,679,620,602,446,621,576,562,272,501,409,452,409,452]
# ci/forkod.f JFOR(6)/KFOR(6): Central Idaho national forests + geographic-location class.
const CI_JFOR = Int[117,402,406,412,413,414]
const CI_KFOR = Int[1,2,2,3,2,2]

"ci/habtyp.f numeric path: habitat code KODTYP → (ICINDX 1..130, ITYPE 1..30). Bracket search ICITYP."
function ci_habtyp(kodtyp::Integer)
    (kodtyp < 10 || kodtyp > 999) && return (1, 1)      # ERRGRO(14) → I=1 path
    ii = findfirst(i -> kodtyp < CI_ICITYP[i], 1:130)
    i = ii === nothing ? 131 : ii
    i == 1 && return (1, 1)
    return (i - 1, Int(CI_NIHMAP[i-1]))                  # (ICINDX, ITYPE)
end

"ci/forkod.f: KODFOR → (IFOR 1..6, IGL=KFOR). Not found ⇒ IFOR=1, no IGL."
function ci_forkod!(p)
    kodfor = Int(p.user_forest_code)
    if kodfor == 7721 || kodfor == 8107                   # ci/forkod.f reservation pseudo-codes → IFOR=3
        p.forest_idx = Int32(3); p.geo_location = Int32(CI_KFOR[3])
        return 3
    end
    idx = findfirst(==(kodfor), CI_JFOR)
    if idx === nothing
        # ci/forkod.f CASE DEFAULT: an unrecognized code keeps the grinit.f:196 default IFOR=4 (Payette) and
        # does NOT update IGL (USEIGL=.FALSE.) ⇒ geo_location stays at its init. IFOR=4≥2 ⇒ Zeide R4SDI SDImax.
        p.forest_idx = Int32(4)
    else
        p.forest_idx = Int32(idx)
        p.geo_location = Int32(CI_KFOR[idx])
    end
    return Int(p.forest_idx)
end

"ci/sitset.f: SITEAR site-range default + forest-dependent SDImax (IFOR<2 Stage / IFOR≥2 Zeide R4SDI)."
function ci_sitset!(s::StandState, icindx::Int, ifor::Int)
    p = s.plot; sd = s.coef.species
    slo_a = sd[:site_lo]; shi_a = sd[:site_hi]
    isisp = Int(p.site_species)
    tem = 50.0f0
    (isisp > 0 && p.sp_site_index[isisp] > 0f0) && (tem = p.sp_site_index[isisp])
    isisp == 0 && (isisp = 3)
    slossp = slo_a[isisp]; shissp = shi_a[isisp]         # SITERANGE(isisp)
    @inbounds for i in 1:19
        tem < slossp && (tem = slossp)
        slo = slo_a[i]; shi = shi_a[i]                   # SITERANGE(i)
        if p.sp_site_index[i] <= 0f0
            p.sp_site_index[i] = slo + (tem - slossp) / (shissp - slossp) * (shi - slo)
        end
    end
    # SDImax (SDIDEF). PMSDIU is a PERCENT here (CI uses PMSDIU/100).
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 85.0f0
    bamax = s.control.ba_max
    if bamax > 0f0
        @inbounds for i in 1:19
            p.sp_sdi_def[i] <= 0f0 && (p.sp_sdi_def[i] = bamax / (0.5454154f0 * (pmsdiu / 100f0)))
        end
    else
        bamax = CI_BAMAXA[icindx]
        @inbounds for i in 1:19
            if ifor < 2
                p.sp_sdi_def[i] <= 0f0 && (p.sp_sdi_def[i] = bamax / (0.5454154f0 * (pmsdiu / 100f0)))
            else
                p.sp_sdi_def[i] <= 0f0 && (p.sp_sdi_def[i] = CI_R4SDI[i])   # ZEIDE
            end
        end
    end
    return s
end

# ci/cratet.f:127-168 — adjust SITEAR to a 50-YEAR age base for WB/LM/PY (11,12,16),
# whose growth eqns were fit on a 50-yr-base site index (Alexander-Tackle-Dahms RM-29).
# CRATET does this once at init, after SITSET, using stand CCF (TEMCCF, floored 125; DBH-only
# open-grown CCF, valid at site_setup!). Inert unless a WB/LM/PY tree/site-species is present.
function ci_cratet_site_adjust!(s::StandState)
    p, t = s.plot, s.trees
    temccf = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        temccf += ci_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
    end
    temccf < 125f0 && (temccf = 125f0)
    @inbounds for sp in (11, 12, 16)
        si = p.sp_site_index[sp]
        p.sp_site_index[sp] = 9.89311f0 - 0.19177f0 * 50f0 + 0.00124f0 * 50f0^2 -
            0.00082f0 * (temccf - 125f0) * si + 0.01387f0 * 50f0 * si -
            0.0000455f0 * 50f0^2 * si
    end
    return s
end

function ci_site_index_setup!(s::StandState)
    p = s.plot
    ifor = ci_forkod!(p)
    kodtyp = Int(p.habitat_code)
    icindx, itype = kodtyp > 0 ? ci_habtyp(kodtyp) : (max(Int(p.habitat_input), 1), max(Int(p.habitat_input), 1))
    (icindx < 1 || icindx > 130) && (icindx = 1)
    p.habitat_input = Int32(icindx)                       # DG reads ICHBCL(ICINDX) ← stash ICINDX here
    ci_sitset!(s, icindx, ifor)
    ci_cratet_site_adjust!(s)                             # CRATET 50-yr-base site adjust (WB/LM/PY)
    return s
end

site_setup!(s::StandState, ::CentralIdaho) = ci_site_index_setup!(s)
