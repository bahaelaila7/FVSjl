# =============================================================================
# fia_database.jl — FVS DBS database INPUT path (DATABASE / DSNIN / StandSQL / TreeSQL)
#
# Faithful port of the FVS "FVS-ready" database reader, so FVSjl consumes the SAME
# keyfile live FVS reads (a DATABASE input block) and simulates an FIA stand
# byte-identically:
#   vdbsqlite/dbsstandin.f  — FVS_STANDINIT_{COND,PLOT} columns → stand/plot state
#   dbsqlite/dbstreesin.f   — FVS_TREEINIT_{COND,PLOT}  columns → tree records
#   base/notre.f            — TPA expansion of the raw TREE_COUNT (already `notre!`)
#
# The stand SQL row plays the role of the STDINFO/SITECODE/DESIGN/GROWTH cards; the
# tree SQL rows play the role of TREEDATA. TREE_COUNT is the raw per-plot count
# (PROB) — `notre!` expands it to trees/acre from BAF/FPA/PI exactly as FVS does.
# =============================================================================

using SQLite
using DBInterface

# Run `sql` (with %StandID% / %Stand_CN% substituted by `sid`) and return every row
# as a Dict{String,Any} keyed by UPPERCASE column name; SQLite NULL stays `missing`.
function _fia_rows(db::SQLite.DB, sql::AbstractString, sid::AbstractString)
    q = sql
    for tok in ("%StandID%", "%Stand_ID%", "%StandCN%", "%Stand_CN%", "%STANDID%")
        q = replace(q, tok => sid)
    end
    out = Dict{String,Any}[]
    for r in DBInterface.execute(db, q)
        d = Dict{String,Any}()
        for nm in propertynames(r)
            d[uppercase(String(nm))] = getproperty(r, nm)
        end
        push!(out, d)
    end
    return out
end

_fia_present(d, k) = haskey(d, k) && d[k] !== missing && d[k] !== nothing
# Some FVS-ready columns are TEXT-typed but hold numbers (e.g. TREEINIT.SEVERITY3 = "2.0"); live FVS parses
# them. `_fia_num` accepts a real OR a numeric string (tryparse → nothing if unparseable → caller uses dv).
_fia_num(x) = x isa AbstractString ? tryparse(Float64, strip(x)) : Float64(x)
_fia_f32(d, k, dv) = (_fia_present(d, k) && (v = _fia_num(d[k])) !== nothing) ? Float32(v) : dv
_fia_int(d, k, dv) = (_fia_present(d, k) && (v = _fia_num(d[k])) !== nothing) ? round(Int, v) : dv
_fia_str(d, k, dv) = _fia_present(d, k) ? strip(string(d[k])) : dv

# FIA numeric species codes arrive un-padded from the DB (e.g. "71"), but the FVS species
# tables key on the 3-digit FIA code ("071"). Zero-pad a purely-numeric 1–2 char code to 3
# digits (dbstreesin.f:353-360 auto-pads a 2-char numeric code with a leading '0'); leave
# alpha codes untouched. Without this, every FIA code < 100 mis-resolves to Other-Hardwood.
_fia_spcode(c::AbstractString) = (s = strip(c); (!isempty(s) && all(isdigit, s) && length(s) <= 2) ? lpad(s, 3, '0') : String(s))

"""
    apply_fia_stand!(s, d)

Map one FVS_STANDINIT_COND/PLOT row `d` (uppercase-keyed Dict) into the stand's
plot/control state — the STDINFO/SITECODE/DESIGN/GROWTH-card equivalent (dbsstandin.f).
"""
function apply_fia_stand!(s::StandState, d::Dict{String,Any})
    p = s.plot; c = s.control
    # INV_YEAR (IY(1)) → first-cycle (start) year
    iy = _fia_int(d, "INV_YEAR", 0); iy > 0 && (c.cycle_year[1] = Int32(iy))
    # LOCATION (KODFOR): direct if present, else composite REGION*100 + FOREST (dbsstandin.f:569).
    # ⚠ FAITHFULNESS NOTE (latent, inert on all FIA data): dbsstandin.f uses a WIDER composite for a couple of
    # variants — KT builds RFFDDCCC (region*10^6 + forest*10^3 + district*... ) and WS/SN build RFFDD
    # (region*10^4 + forest*100 + district). We deliberately DO NOT port those encodings: the REGION*100+FOREST
    # fallback only fires when LOCATION is absent, which never happens on the western FIA tables (LOCATION is
    # populated on every stand); KT has ZERO FIA stands (so its encoding can never be oracle-validated); and WS's
    # composite equals its LOCATION anyway. Documenting rather than adding unexercised, un-testable encoding code.
    loc = _fia_int(d, "LOCATION", 0)
    (loc == 0 && _fia_present(d, "REGION")) && (loc = _fia_int(d, "REGION", 0) * 100 + _fia_int(d, "FOREST", 0))
    loc != 0 && (p.user_forest_code = Int32(loc))
    # Fort Bragg (forkod.f CASE 701): a region-7/forest-1 FIA code (composite LOCATION=701) is remapped to
    # NC Uwharrie 81110 (region 8) — forkod runs for every stand, incl. DB input. Without it VOLEQDEF sees
    # region 7 ⇒ no R8 Clark equation ⇒ zero volume. (Shared with kw_stdinfo!.)
    s.variant isa Southern && sn_forkod_remap!(p)
    # AGE (IAGE)
    _fia_present(d, "AGE") && (p.stand_age = Int32(_fia_int(d, "AGE", 0)))
    # ASPECT degrees → radians (TRNASP ×0.0174533); SLOPE percent → fraction (PLOT.F77)
    _fia_present(d, "ASPECT") && (p.aspect = _fia_f32(d, "ASPECT", 0f0) * 0.0174533f0)
    _fia_present(d, "ASPECT") && (p.aspect_deg = trunc(Int32, _fia_f32(d, "ASPECT", 0f0)))   # IASPEC=IFIX(ASPECT), initre.f:436
    # SLOPE: grinit.f:226 defaults a MISSING/NULL slope to 5.0% (→0.05 fraction) BEFORE the DB
    # overrides it — all 4 variants (sn/ne/cs/grinit.f:221-226 SLOPE=5.0). jl previously left a
    # missing slope at the 0.0 constructor default, which zeroed the DGF slope/aspect DGCON term
    # (TANS·SLOPE + FCOS·SLOPE·cos(ASP) + FSIN·SLOPE·sin(ASP), dgf.f:1125-1127) → over-grew species
    # with large slope coefficients (e.g. sp39 loblolly-bay FCOS=-10.15: a 0.05 slope = -0.68 in
    # ln(DDS) ⇒ ~2× DBH growth). Apply the grinit default so a missing slope matches live FVS.
    p.slope = _fia_present(d, "SLOPE") ? _fia_f32(d, "SLOPE", 0f0) / 100f0 : 5f0 / 100f0
    p.slope_raw = _fia_present(d, "SLOPE") ? trunc(Int32, _fia_f32(d, "SLOPE", 0f0)) : Int32(5)   # ISLOP=IFIX(SLOPE), initre.f:437
    # ELEVATION in hundreds of feet; ELEVFT is feet → ×0.01 (dbsstandin.f:710).
    # ⚠ METRIC DBs (BC/ON) store ELEVATION in METRES: the metric dbsstandin.f (FVSbc_buildDir, header
    # "METRIC-VDBSQLITE") does RSTANDDATA(9) = ELEVATION * MtoFt / 100 (:351) — metres→hundreds-of-feet —
    # while the western VDBSQLITE reader (FVSem_buildDir:303) uses ELEVATION raw (already hundreds-of-ft).
    # jl elevation feeds the BC crown model (CRCON EL·elev + EL2·elev²) and AUTOES, so a raw-metre store is
    # ~30.48× too large. ELEVFT is feet in both readers (metric reader still ×0.01, no MtoFt).
    # ⚠ Gate on a NUMERIC ELEVFT, not mere presence: this FVS-ready DB stores an unset numeric column as an
    # EMPTY STRING (e.g. LD3001 ElevFt="" while Elevation=518 m), which _fia_present() reports as present.
    # A blank/non-numeric ELEVFT would then set elevation=0 and short-circuit the ELEVATION fallback, dropping
    # the real metres value. FVS's dbsstandin treats a null ELEVFT as absent and reads ELEVATION — mirror that.
    if _fia_present(d, "ELEVFT") && _fia_num(d["ELEVFT"]) !== nothing
        p.elevation = _fia_f32(d, "ELEVFT", 0f0) * 0.01f0
    elseif _fia_present(d, "ELEVATION") && _fia_num(d["ELEVATION"]) !== nothing
        if s.variant isa BritishColumbia || s.variant isa Ontario
            p.elevation = _fia_f32(d, "ELEVATION", 0f0) * 3.280839895f0 / 100f0   # m → hundreds of ft (BC/ON metric DB)
        else
            p.elevation = _fia_f32(d, "ELEVATION", p.elevation)
        end
    end
    # Western variants set a per-variant DEFAULT elevation (hundreds of ft) in grinit.f (EM 55, BM 45, IE 38,
    # KT 35, UT 83, TT 65, CI 50); the DB overrides ONLY when >0 (dbsstandin.f:647), so a NULL/≤0 ELEVATION/
    # ELEVFT keeps that default (a present-but-NULL ELEVFT column reads as 0 above). Without it the ESTOCK and
    # DG elevation terms use 0 ⇒ wrong AUTOES PROB1 (over-produces: measured EM stand 5332701010661 NULL elev →
    # live grinit ELEV=55 → ESTOCK PN −1.490, PROB1 0.184; jl elev 0 → PN −0.001, PROB1 0.4996 = 2.7× the trees)
    # AND wrong large-tree DG (EM_DGEL·elev + EM_DGEL2·elev²). CR grinit ELEV=0 (no default); SN uses forest_location.
    if p.elevation <= 0f0
        edf = s.variant isa EasternMontana ? 55f0 : s.variant isa BlueMountains ? 45f0 :
              s.variant isa InlandEmpire   ? 38f0 : s.variant isa Kootenai      ? 35f0 :
              s.variant isa Utah           ? 83f0 : s.variant isa Teton         ? 65f0 :
              s.variant isa CentralIdaho   ? 50f0 : 0f0
        edf > 0f0 && (p.elevation = edf)
    end
    # LATITUDE/LONGITUDE (TLAT/TLONG, dbsstandin.f:254-259) — feed the Hopkins bioclimatic
    # index in the eastern crown-width models.
    _fia_present(d, "LATITUDE")  && (p.latitude  = _fia_f32(d, "LATITUDE", p.latitude))
    _fia_present(d, "LONGITUDE") && (p.longitude = _fia_f32(d, "LONGITUDE", p.longitude))
    # ECOREGION (ecological unit / EUT, e.g. "223Db") → eco_unit. FVS reads it from STANDINIT and adds a
    # per-species ecological-unit DG term (dgf.f EUT categorical coefficients dg_phys_*), plus it drives the
    # montane site/height/estab branches (eco_unit[1]=='M'). Without it the whole EUT DG term is dropped for
    # FIA-DB stands — the FIA/FVS campaign's slice-1 divergence (yellow-poplar large-tree DG ~20% low on a
    # 223Db stand: jl omitted dg_phys_p222 ≈ +0.255 of the −0.344 ln(DDS) deficit). SN-gated: resolve_eco_unit
    # / SNECU is the SOUTHERN ecological-unit table; NE/CS/LS use their own eco-unit handling (left blank as
    # before, a documented follow-up if their FIA differentials show an analogous gap).
    if s.variant isa Southern && _fia_present(d, "ECOREGION")
        p.eco_unit = rpad(resolve_eco_unit(_fia_str(d, "ECOREGION", ""), 0), 10)
    elseif s.variant isa BritishColumbia && _fia_present(d, "ECOREGION")
        # BC's ECOREGION is the BEC string (e.g. 'CAR-IDFdk3/01'), which bc_habtyp/bc_becset parse to the
        # zone/subzone/series driving the V3 DGCON coefficients (dgf.f ZNKONST/SSKONST). Store it RAW — do
        # NOT route through the SN resolve_eco_unit (which mangles it to a garbage EUT code). Without this the
        # DB path leaves eco_unit blank ⇒ bc_becset DEFAULTS to ICHmw2/01 ⇒ WRONG DGCON for every non-ICH DB
        # stand (e.g. YSM029-271 is IDFdk3: ICH cedar-hemlock coeffs over-predict Pl DG on a dry IDF stand).
        p.eco_unit = rpad(_fia_str(d, "ECOREGION", ""), 10)
    end
    # CPVREF (dbsstandin.f:568-574): PV_REF_CODE > 0 ⇒ WRITE(CPVREF,'(I10)') (PVREF1/6 compare it ADJUSTL'd), ≤0 ⇒ blank.
    if _fia_present(d, "PV_REF_CODE")
        r = _fia_f32(d, "PV_REF_CODE", 0f0)
        p.pv_ref = r > 0f0 ? string(Int(round(r))) : ""
    end
    # dbsstandin.f:579-592 — HABTYP runs for a non-blank PV_CODE (SN with an ECOREGION skips it); only its ERRGRO
    # branches (FVS14/32/33/34) are mirrored here, the habitat itself is resolved below.
    if _fia_present(d, "PV_CODE") && !isempty(strip(_fia_str(d, "PV_CODE", ""))) &&
       !(s.variant isa Southern && _fia_present(d, "ECOREGION"))
        chab = strip(_fia_str(d, "PV_CODE", ""))
        habtyp_errors!(s, chab, strip(p.pv_ref), something(tryparse(Int, chab), 0))   # READ(CHAB,'(I10)',ERR=40)
        s.control.habtyp_done = true
    end
    # PV_CODE (potential-vegetation / habitat-type code, e.g. 531) → habitat_code (KODTYP), the input to
    # habtyp. KT/western: the DG habitat term (KKTYPE→MAPHAB→DGHAB) needs it; site_setup!(::Kootenai) maps
    # KODTYP→KKTYPE. Eastern variants key DG off forest type (habitat-input left a documented gap there), so
    # KT-gated. Live FVSkt reads PV_CODE from the DB automatically; this matches it on the jl side.
    if s.variant isa Kootenai && _fia_present(d, "PV_CODE")
        p.habitat_code = Int32(round(_fia_f32(d, "PV_CODE", 0f0)))
    end
    # BM (region-6): PV_CODE is the ALPHA plant-community code (e.g. "CWF312"), decoded to the KODTYP index
    # into BM_PCOML by habtyp/HBDECD. Without it, bm_sitset! gets no ECOCLS row ⇒ SDIDEF=0 ⇒ stand_sdimax=0 ⇒
    # bm/morts.f's "SDIMAX<5 ⇒ kill ALL trees" fires and the stand COLLAPSES to 0 TPA at cycle 1 (cycle-0-only
    # sweeps never caught this). Match live FVSbm, which reads PV_CODE and HBDECDs it to the PCOML index.
    if s.variant isa BlueMountains
        pv = _fia_present(d, "PV_CODE") ? String(strip(_fia_str(d, "PV_CODE", ""))) : ""
        # bm/habtyp.f: when PV_REF_CODE is present, PVREF6 crosswalks (PV_CODE, PV_REF_CODE) → the canonical
        # PCOML code BEFORE the HBDECD/PCOML match (e.g. "CJG111"/622 → "CPG111"). A raw FIA PV code not in
        # PCOML would otherwise fall through to the CWG113 default (SDIMAX 395 vs 166) ⇒ under-thinning (#140).
        # No PV_REF_CODE ⇒ the raw PV code is matched directly (PVREF6 is not called).
        # ★ pvref6.f BLANKS KARD2 on entry (line 2577) and sets it to HABPVR only on a FULL (pv AND ref) match
        # (EXIT at 2586); ANY unmatched (pv,ref) pair returns BLANK, NOT the raw pv. So when PVREF6 is invoked
        # (pvr>0) and no full match exists, the habitat is UNRESOLVED → sitset's CWG113 default (SDIMAX≈433 for a
        # PP/DF stand), NOT the raw code's ecoclass. jl previously KEPT the raw pv on no-match: e.g. CDG111/653 has
        # no PVREF6 row (CDG111 exists only for refs 604/622/639/640), so live blanks it → CWG113 433, but jl kept
        # CDG111 → SDIMAX 278/351 ⇒ jl OVER-thins (#140 over-thin case, measured: 504530915126144 −52% TPA by 2097
        # vs live; live ZZPVREF out=BLANK, ZZMORTS SDIMAX=433). Assign unconditionally — mapped is "" on no-match.
        if !isempty(pv) && _fia_present(d, "PV_REF_CODE")
            pvr = Int(round(_fia_f32(d, "PV_REF_CODE", 0f0)))
            if pvr > 0
                pv = bm_pvref6(pv, string(pvr))   # PVREF6 blanks on no-match (mapped=="")
            end
        end
        hc = 0
        if !isempty(pv)
            idx = findfirst(==(pv), BM_PCOML)
            idx !== nothing && (hc = Int(idx))
        end
        # bm/habtyp.f:67-69 — a MISSING or UNRESOLVED habitat defaults to CWG113 = KODTYP 79 (measured live
        # ICL5=79 on habitat-less FIA stands). jl previously left habitat_code at 1 here ⇒ the DF small-tree
        # SMCON picked SMHAB(2,2)=−0.337 instead of the neutral SMHAB(1,·)=0 ⇒ ~6% DF DG under-shoot ⇒ the #140
        # self-thin under-kill (9:2 skew). Setting 79 makes the DF small-tree DDS bit-exact vs live.
        p.habitat_code = Int32(hc == 0 ? 79 : hc)
    end
    # CI (ci/habtyp.f:44-75 via dbsstandin.f:579-592): with a PV_REF_CODE, PVREF4 crosswalks the (PV_CODE text,
    # PV_REF_CODE) pair to the CI habitat KODTYP; a pair PVREF4 does not know (FVS34/33/32) keeps the grinit default
    # ICINDX=21 (habitat 260). Without a reference code KODTYP = the PV_CODE's integer value. jl used "PV_CODE mod
    # 1000" instead — the 9999999 "no habitat" sentinel became habitat 999 and unknown pairs (45101/494) a real
    # code (MEASURED FVSci_g16: 3369538010690 9999999/491 and 12276084010690 45101/494 both "MAPPED TO 260").
    if s.variant isa CentralIdaho
        pvs = _fia_present(d, "PV_CODE") ? String(strip(_fia_str(d, "PV_CODE", ""))) : ""
        if !isempty(pvs)
            kod, _, lpvxxx = ci_habitat_kodtyp(pvs, strip(p.pv_ref), something(tryparse(Int, pvs), 0))
            p.habitat_code = Int32(lpvxxx ? 0 : kod)
        end
    end
    # EC (region-6, Wykoff-DDS): PV_CODE is the ALPHA plant-association code (e.g. "CDG131"). The STDINFO
    # KEYWORD path was fixed in 70535053 (ec_hbdecd), but the FIA-DB PV_CODE column had no reader branch ⇒
    # habitat_code stayed 0 ⇒ ec_sitset! grew every FIA stand on the poor default site (measured live vs jl,
    # CN 1143082551290487: live SDI MAX 530..629 per species vs jl uniform 331). ec/habtyp.f PVREF6-crosswalks
    # (PV_CODE,PV_REF_CODE)→HABPVR then HBDECDs it to the KODTYP index into EC_PCOML. PV_CODE priority.
    if s.variant isa EastCascades
        pv = _fia_present(d, "PV_CODE") ? String(strip(_fia_str(d, "PV_CODE", ""))) : ""
        pvref = ""
        if _fia_present(d, "PV_REF_CODE")
            r = _fia_f32(d, "PV_REF_CODE", 0f0); r > 0f0 && (pvref = string(Int(round(r))))
        end
        if !isempty(pv)
            hc = ec_habitat_kodtyp(pv, pvref)
            hc != 0 && (p.habitat_code = Int32(hc))
        end
    end
    # CA (region-5/6, Wykoff-DDS): PV_CODE is the ALPHA plant-association / ecoclass code (e.g. "CDH524").
    # For R6 forests (KODFOR≥600 — ~99.7% of CA's FIA population, LOCATION 610/611) ca/habtyp.f PVREF6-
    # crosswalks (PV_CODE,PV_REF_CODE)→HABPVR then HBDECDs it to the KODTYP index into CA_PCOML that
    # ca_sitset! consumes. Without it habitat_code stays 0 ⇒ ca_sitset! defaults to CWC221 (SDImx 815)
    # instead of the stand's ecoclass (CDH524/641 → CDS511 → 635), inflating the BA-weighted SDIMAX ~1.3×
    # so the density self-thin under-fires ⇒ dense stands keep too much TPA (measured live vs jl; cyc0/ref
    # stands never caught it — cat01's habitat resolves to the benign CWC221 default). PV_CODE priority.
    # (R5 forests <600 use the ca/pvref5.f + R5HABT path, which jl's R6-only ca_sitset! does not implement —
    #  only ~6 CA FIA stands, left at the CWC221 default as before.)
    if s.variant isa CentralCalifornia && Int(p.user_forest_code) >= 600
        pv = _fia_present(d, "PV_CODE") ? String(strip(_fia_str(d, "PV_CODE", ""))) : ""
        pvref = ""
        if _fia_present(d, "PV_REF_CODE")
            r = _fia_f32(d, "PV_REF_CODE", 0f0); r > 0f0 && (pvref = string(Int(round(r))))
        end
        if !isempty(pv)
            hc = ca_habitat_kodtyp(pv, pvref)
            hc != 0 && (p.habitat_code = Int32(hc))
        end
    end
    # NC (Klamath, region-6/NE-California Wykoff-DDS): PV_CODE is the ALPHA plant-association code (e.g.
    # "HTS121"). For NC's R6 forests (611 Siskiyou = IFOR 4, 712 BLM Coos Bay = IFOR 7) nc/habtyp.f
    # PVREF6-crosswalks (PV_CODE,PV_REF_CODE)→HABPVR (no ref ⇒ the raw code) then HBDECDs it to the KODTYP
    # index into NC_PCOML that nc_sitset! consumes to seed the ECOCLS per-species SDImax. Without it
    # habitat_code stays 0 ⇒ nc_sitset! rode a provisional uniform default (~720) instead of the ecoclass
    # (CWC221 default → DF 815 + C6 fan capped at 850 = live's "SDI MAX 850 850 815 850..."; a real PA → its
    # own RSDI). MEASURED live vs jl, CN 1127525637290487. PV_CODE priority (nc/habtyp.f).
    if s.variant isa Klamath && Int(p.user_forest_code) in (611, 712, 8103, 8105)   # 8103/8105 → 611 (nc/forkod.f)
        pv = _fia_present(d, "PV_CODE") ? String(strip(_fia_str(d, "PV_CODE", ""))) : ""
        pvref = ""
        if _fia_present(d, "PV_REF_CODE")
            r = _fia_f32(d, "PV_REF_CODE", 0f0); r > 0f0 && (pvref = string(Int(round(r))))
        end
        if !isempty(pv)
            hc = nc_habitat_kodtyp(pv, pvref)
            hc != 0 && (p.habitat_code = Int32(hc))
        end
    end
    # SO (region-6/NE-California, Wykoff-DDS): PV_CODE is the ALPHA plant-association code (e.g. "CWS313").
    # For SO's R6 forests (601/602/620 = IFOR 1-3, 799 = IFOR 10; reservation 7710/7711 → 602) so/habtyp.f
    # PVREF6-crosswalks (PV_CODE,PV_REF_CODE)→HABPVR (no ref ⇒ the raw code) then HBDECDs it to the KODTYP
    # index into SO_PCOML that so_sitset! consumes to seed the site species' ecoclass SDImax. Without it
    # habitat_code stays 0 ⇒ so_sitset! rode the CPS111 default RSDI (285) for EVERY stand instead of the
    # stand's own PA (CWS313 → 810), so on a non-default PA the SDIMAX was ~2.8× wrong (measured live vs jl,
    # CN 24397775010900) — most SO alpha stands crosswalk to CPS111 (=default) so were already right; the
    # minority on other PAs self-thinned at the wrong level. PV_CODE priority (SO habtyp.f:206).
    if s.variant isa SouthCentralOregon && Int(p.user_forest_code) in (601, 602, 620, 799, 7710, 7711)
        pv = _fia_present(d, "PV_CODE") ? String(strip(_fia_str(d, "PV_CODE", ""))) : ""
        pvref = ""
        if _fia_present(d, "PV_REF_CODE")
            r = _fia_f32(d, "PV_REF_CODE", 0f0); r > 0f0 && (pvref = string(Int(round(r))))
        end
        if !isempty(pv)
            hc = so_habitat_kodtyp(pv, pvref)
            hc != 0 && (p.habitat_code = Int32(hc))
        end
    end
    # PN/WC (region-6, ORGANON vwc/morts.f): PV_CODE is the ALPHA plant-association code (e.g. "CHS512"),
    # decoded to the KODTYP index into PN_PCOML/WC_PCOML by habtyp/PVREF6/HBDECD. Without it habitat_code
    # stays 0 ⇒ pn_sitset!/wc_sitset! fall back to the PA default (PN CHS133 SDIDEF 1606→FORMAX 950; WC
    # CFS551 815) instead of the stand's ecoclass SDIDEF (e.g. CHS512 → 485). SDICAL's weighted SDIMAX then
    # runs ~2× high, so vwc/morts.f's SDI density self-thin (the integer-PASS loop, morts.f:456-500) never
    # fires and the stand keeps ~3-4× too much TPA at cycle 1 (measured live vs jl on dense over-max-SDI FIA
    # stands; cyc0/ref-stand sweeps never caught it — the ref stands' habitat resolves to a benign default).
    # Match live FVSpn/FVSwc, which read PV_CODE and PVREF6/HBDECD it to the PCOML index.
    if s.variant isa PacificNorthwest || s.variant isa WestCascades
        pv = _fia_present(d, "PV_CODE") ? String(strip(_fia_str(d, "PV_CODE", ""))) : ""
        pvref = ""
        if _fia_present(d, "PV_REF_CODE")
            r = _fia_f32(d, "PV_REF_CODE", 0f0)
            r > 0f0 && (pvref = string(Int(round(r))))
        end
        if !isempty(pv) || !isempty(pvref)
            hc = s.variant isa PacificNorthwest ? pn_habitat_kodtyp(pv, pvref) : wc_habitat_kodtyp(pv, pvref)
            p.habitat_code = Int32(hc)
        end
    end
    # EM/UT/TT/IE (western, habitat-type-group DG): read PV_CODE (KODTYP) like CI/KT/BM. These variants'
    # DG constant + mortality ITYPE key off a habtyp(KODTYP) map (em_habtyp/ut_habtyp/tt_habtyp/ie_habtyp),
    # but the FIA reader never set habitat_code for them ⇒ it defaulted to 0/1 ⇒ WRONG habitat vs live.
    # Live e.g. FVSem prints "HABITAT TYPE MAPPED TO 470 ... PV_REF_CODE WAS IGNORED" (PV_CODE priority).
    # A 5-digit PV_CODE (e.g. 41780 = 2-digit state prefix + 3-digit habitat) → 780 = 41780 mod 1000; the
    # 3-digit code (470/310) is used directly. PV_REF_CODE is the fallback. Each variant's site_setup! then
    # maps habitat_code → habitat_input via its own habtyp. Confirmed live bug on EM stand 12344705010690
    # (live habitat 470, jl was defaulting to 1).
    if s.variant isa EasternMontana || s.variant isa Utah || s.variant isa Teton || s.variant isa InlandEmpire
        hc = 0
        isie = s.variant isa InlandEmpire
        pvref = _fia_present(d, "PV_REF_CODE") ? Int(round(_fia_f32(d, "PV_REF_CODE", 0f0))) : 0
        # ★#143 (IE): ie/pvref1.f crosswalks the (PV_CODE, PV_REF_CODE) PAIR → HABPVR and takes PRECEDENCE whenever
        # a reference code is present (habtyp.f:91). A full match yields the habitat; a present-but-UNRECOGNIZED pair
        # (e.g. PV_CODE "ABR8" + ref 639, 4/5 IE sweep stands) defaults to 260 — NOT the raw reference code. jl's old
        # bug fell back to the raw PV_REF_CODE (639) ⇒ ESTOCK ihab 11 (GF-dominant) vs live's ihab 3 (DF/PP) ⇒ AUTOES
        # over-establishment. MEASURED vs live esplt2.f/habtyp.f (ihab/MYGRUP/OCURHT/OCURNF/PADV all match at ihab 3).
        # ★ IE: when a PV_REF_CODE is present (CPVREF≠blank), ie/habtyp.f calls PVREF1 which BLANKS KARD2 and zeros
        # KODTYP UNLESS it finds a full (PV_CODE,PV_REF) match. On a non-match the subsequent HBDECD therefore sees a
        # BLANK code ⇒ default 260 — the oracle does NOT re-try the bare PCOML string match NOR the raw numeric
        # PV_CODE. So gate BOTH fallbacks below on `!(isie && pvref>0)`: for an IE stand with a ref, only ie_pvref1's
        # result counts, else default 260. (CWS422/626 → jl formerly fell back to ie_pa_habitat_code=520; 41691/401 →
        # jl formerly %1000-stripped to 691. MEASURED vs FVSie_g16 habtyp.f+pvref1.f: both default to 260.)
        ie_ref_present = isie && pvref > 0
        # EM: em/habtyp.f calls PVREF1 whenever a PV reference code is present (CPVREF≠blank); only a full
        # (PV_CODE, PV_REF_CODE) match in em/pvref1.f resolves the habitat, anything else defaults to 260 —
        # the raw PV code is NOT used (jl fed 130 straight to em_habtyp ⇒ IEMTYP ≠ live 29 ⇒ wrong EMMD fuel
        # models, 474157097489998 fire flame 5.03 vs live 1.95).
        em_ref_present = s.variant isa EasternMontana && pvref > 0 && _fia_present(d, "PV_CODE")
        if em_ref_present
            hc = em_pvref1_kodtyp(_fia_str(d, "PV_CODE", ""), pvref)
            # dbsstandin.f:590 ICL5 = KODTYP = IFIX(numeric PV_CODE) BEFORE HABTYP remaps KODTYP — the DGCONS
            # habitat (rcon.f IDTYPE) keys off this RAW code, not the PVREF1/default IEMTYP (live 474157097489998:
            # KODTYP 130 → 260 for FMCFMD, but DG on IHCODE(130)).
        end
        # dbsstandin.f:590-593 (EVERY stand with a PV_CODE, ref or not): ICL5 = KODTYP = IFIX(numeric PV_CODE) BEFORE
        # HABTYP translates KODTYP; em/habtyp.f never touches ICL5, and `IF(ICL5.LE.0) ICL5=KODTYP` only fires for a
        # non-numeric code (READ …'(I10)' ERR ⇒ 0). ICL5 is "the actual input habitat code" (grinit.f:200): rcon.f's
        # DGCONS IDTYPE and esplt2.f's AUTOES habitat bracket (IHTYPE) both key off it, not the translated KODTYP.
        # MEASURED FVSem_g16: 3006831010690 PV_CODE 9999999 (→ habitat 260) runs ESTAB with ISER 5 (bracket(9999999)
        # = group 16, the AF series), 2999215010690 PV_CODE 356 (→ 250) with bracket(356) = group 4.
        if s.variant isa EasternMontana && _fia_present(d, "PV_CODE")
            raw = tryparse(Int, strip(_fia_str(d, "PV_CODE", "")))      # READ(CHAB,'(I10)'): integer text only
            (raw !== nothing && 0 < raw < typemax(Int32)) && (s.control.icl5 = Int32(raw))
        end
        # ★ ESTAB ICL5 fix: for IE a 6-char plant-association PV_CODE is decoded by ie/habtyp.f's PCOML block
        # (dbsstandin.f DB path keeps KARD2 as the alpha code ⇒ HBDECD→NITYPE→ICL5=JTYPE(NITYPE)), INDEPENDENT
        # of PV_REF_CODE. That North-Idaho code (NOT pvref1's HABPVR, NOT the MTYPE growth code) is what
        # esplt2.f/estab.f key the AUTOES habitat bracket off. MEASURED FVSie_g16 645170890126144 (PV_CODE
        # CWS821, ref 627): ICL5=590=JTYPE(55) (growth KODTYP 510), ESTOCK IHAB=5. jl formerly stored the MTYPE/
        # pvref1 code (510) ⇒ estab bracket group 6 (SHAB −0.06) vs the oracle's group 5 (SHAB +0.54), a −0.60
        # stocking logit that under-produced AUTOES ingrowth ~100 TPA/cycle. Growth is invariant: ie_habtyp(590)
        # == ie_habtyp(510) == ITYPE 12 (MTYPE codes are the canonical JTYPE representative of each ITYPE group),
        # so the DG/site path re-derives the identical habitat_input. Take precedence for any recognized 6-char code.
        if isie
            ni = ie_pa_ni_code(_fia_str(d, "PV_CODE", ""))
            ni != 0 && (hc = ni)
        end
        if hc == 0 && ie_ref_present
            # ★ ESTAB ICL5 (numeric+ref path): ie/pvref1.f gives the HABPVR *intermediate* KODTYP (e.g.
            # PV_CODE 508 / ref 110 → 506). For a NUMERIC PV_CODE the oracle's dbsstandin.f then stores
            # ICL5 = the FINAL MTYPE-mapped KODTYP (dbsstandin.f:590-593: ICL5=KODTYP; CALL HABTYP zeros ICL5
            # on the numeric path since habtyp.f:125 sets ICL5 only inside the 6-char/PCOML alpha block; then
            # `IF(ICL5.LE.0) ICL5=KODTYP` ⇒ ICL5 = the projection code 510 the .out reports "MAPPED TO 510").
            # esplt2.f keys the AUTOES habitat bracket off ICL5, so the estab IHAB comes from 510 (IHAB 6), NOT
            # the pvref1 intermediate 506 (IHAB 8). jl stored 506 ⇒ estab bracketed to IHAB 8 ⇒ PP occupancy
            # (OCURHT) opened at 8 that is closed by the correct habitat's species-probs ⇒ a spurious fast-growing
            # ponderosa cohort on bare stands (measured FVSie_g16 estab-debug: over-stand ICL5→IHAB=6; jl→8). Apply
            # ie_estab_kodtyp (= habtyp.f:134-147 MTYPE remap). GROWTH IS PROVABLY INVARIANT: ie_estab_kodtyp(x) =
            # MTYPE[ie_habtyp(x)] is the canonical JTYPE representative of x's ITYPE ⇒ ie_habtyp(ie_estab_kodtyp(x))
            # == ie_habtyp(x) (ie_habtyp(506)==ie_habtyp(510)==ITYPE 12), so the DG/site path is unchanged. The 6-char
            # ALPHA path (ie_pa_ni_code above) already stores ICL5=JTYPE (habtyp.f:125) and must NOT be remapped.
            pv = ie_pvref1(_fia_str(d, "PV_CODE", ""), pvref)
            hc = pv == 0 ? 0 : ie_estab_kodtyp(pv)
        end
        # PV_CODE as a 6-char plant-association string given directly (ie/habtyp.f PCOML path, e.g. "CDS715"→260).
        # NOT run when a ref was present (PVREF1 already blanked KARD2 ⇒ no PCOML retry).
        if hc == 0 && isie && !ie_ref_present && _fia_present(d, "PV_CODE")
            hc = ie_pa_habitat_code(_fia_str(d, "PV_CODE", ""))
        end
        # numeric PV_CODE (5-digit 41780 → 780; 3-digit used directly). The FIA "not collected" sentinel
        # 9999999 (and any ≥6-digit value) is NOT a habitat code: ie/habtyp.f HBDECD fails to recognize it ⇒
        # the code is unrecognized ⇒ the default habitat (below), NOT 9999999%1000=999 (ITYPE 30). The %1000
        # strip is only for a genuine 5-digit state-prefixed code (41780 → 780).
        # NOT run for an IE stand with a ref present (PVREF1 zeroed KODTYP ⇒ the raw numeric value is discarded).
        if hc == 0 && !ie_ref_present && !em_ref_present && _fia_present(d, "PV_CODE")
            pvc = Int(round(_fia_f32(d, "PV_CODE", 0f0)))
            if pvc < 100000
                # EM/UT/TT strip a 2-digit state prefix (41780 → 780); IE does NOT — ie/habtyp.f keeps the raw
                # KODTYP, and a >999 value is out of range ⇒ default 260 (MEASURED 41691 → 260, NOT 691).
                !isie && pvc > 999 && (pvc = pvc % 1000)
                # IE: mirror dbsstandin.f/habtyp.f — the ESTAB ICL5 for a NUMERIC PV_CODE is the FINAL
                # MTYPE-mapped KODTYP (see the pvref1 branch above), not the raw code. Growth-invariant.
                (10 <= pvc <= 999) && (hc = isie ? ie_estab_kodtyp(pvc) : pvc)
            end
        end
        # Unresolved habitat: IE ⇒ live default habitat 260 (ie/habtyp.f: an out-of-range/unrecognized KODTYP
        # leaves ITYPE at the grinit default 4, then KODTYP=MTYPE(4)=260 — the ".out" "HABITAT TYPE WILL BE
        # MAPPED TO 260"). This default is UNCONDITIONAL for IE (matches FVS, which always ends KODTYP=MTYPE(ITYPE));
        # without it a no-PV/sentinel-PV IE stand kept habitat_code=999 (ITYPE 30) or 0, so AUTOES (ie_estab_kodtyp)
        # picked the wrong climax species (ES/AF vs DF/PP) ⇒ half-height regen establishment ⇒ the D2 growth deficit.
        # EM/UT/TT ⇒ the raw PV_REF_CODE (their readers use it directly — validated separately; unchanged).
        if hc == 0
            hc = isie ? 260 : (pvref > 0 ? ((10 <= pvref <= 999) ? pvref : 0) : 0)
        end
        hc != 0 && (p.habitat_code = Int32(hc))
    end
    # FORKOD phase-3 default (forkod.f:540-546, mirrored from kw_stdinfo!): fill any geo field the
    # DB left at 0 from the national-forest table. FVS runs forkod BEFORE the DB overrides, and the
    # DB overrides elevation ONLY when >0 (dbsstandin.f:647) — so a null/≤0 ELEVATION keeps the
    # forest default (e.g. LOCATION 80215 → forest 802 → 12.0 hundred-ft, live-confirmed). That
    # elevation drives the Hopkins index for hardwood open-grown crowns, so without it HI (and the
    # reported CCF) drift. Southern-gated: the SN forest_location table is keyed by KODFOR÷100 the
    # same way kw_stdinfo! keys it; NE/CS/LS use a different forkod keying (left as a follow-up).
    # ON (canada/on/forkod.f, called by dbsstandin.f:618 BEFORE the DB LATITUDE/LONGITUDE overrides): a code outside
    # JFOR keeps the grinit.f IFOR=9 (⇒ KODFOR=JFOR(9)=915); KODFOR 915/916 sets TLAT=46.78, TLONG=92.11, ELEV=16 where
    # still 0. The FVSDataHardwood LD3001 stand (Region 9, Forest 15) has LATITUDE but no LONGITUDE ⇒ live TLONG=92.11
    # (Hopkins index of the open-grown hardwood crown width ⇒ CCF); jl left 0.
    if s.variant isa Ontario && (Int(p.user_forest_code) in (915, 916) ||
                                 !(Int(p.user_forest_code) in (902, 903, 904, 906, 907, 909, 910, 913, 924)))
        p.latitude  == 0f0 && (p.latitude  = 46.78f0)
        p.longitude == 0f0 && (p.longitude = 92.11f0)
        p.elevation == 0f0 && (p.elevation = 16f0)
    end
    if s.variant isa Southern
        lat0, long0, elev0 = forest_location(s.coef, div(Int(p.user_forest_code), 100))
        p.latitude  == 0f0 && (p.latitude  = lat0)
        p.longitude == 0f0 && (p.longitude = long0)
        p.elevation == 0f0 && (p.elevation = elev0)
    end
    # Sampling design (DESIGN card): BAF / FPA / BRK / IPTINV / NONSTK / SAMWT / GROSPC
    # ON (canada/on dbsstandin.f:357-376, METRIC reader): BASAL_AREA_FACTOR <0 (inverse fixed-plot size, per ha) is
    # /HAtoACR, ≥0 (m²/ha) is *M2pHAtoFT2pACR; INV_PLOT_SIZE /HAtoACR; BRK_DBH (cm) *CMtoIN — so notre.f expands the
    # per-plot TREE_COUNT to trees/ACRE itself (PROB = count·(−BAF)/PI). (TREE_COUNT is NOT converted.)
    on_db = s.variant isa Ontario
    if _fia_present(d, "BASAL_AREA_FACTOR")
        b = _fia_f32(d, "BASAL_AREA_FACTOR", 0f0)
        on_db && (b = b < 0f0 ? b / 2.471f0 : b * 4.3560773f0)
        p.baf = b
    end
    _fia_present(d, "INV_PLOT_SIZE")     && (p.fixed_plot_inv = _fia_f32(d, "INV_PLOT_SIZE", 0f0) / (on_db ? 2.471f0 : 1f0))
    _fia_present(d, "BRK_DBH")           && (p.min_dbh_var_plot = on_db ? _fia_f32(d, "BRK_DBH", 0f0) * 0.3937f0 :
                                                                          _fia_f32(d, "BRK_DBH", p.min_dbh_var_plot))
    _fia_present(d, "NUM_PLOTS")         && (p.points_inv = Int32(_fia_int(d, "NUM_PLOTS", 1)))
    _fia_present(d, "NONSTK_PLOTS")      && (p.nonstockable = Int32(_fia_int(d, "NONSTK_PLOTS", 0)))
    _fia_present(d, "SAM_WT")            && (p.sample_weight = _fia_f32(d, "SAM_WT", p.sample_weight))
    # GROSPC: STK_PCNT (1..100 → ÷100) if given, else (IPTINV − NONSTK)/IPTINV (dbsstandin.f:740)
    if _fia_present(d, "STK_PCNT")
        g = _fia_f32(d, "STK_PCNT", 1f0)
        (g > 1f0 && g <= 100f0) && (g *= 0.01f0)
        (g > 0f0 && g <= 1f0) && (p.gross_space = g)
    else
        ip = max(1, Int(p.points_inv))
        p.gross_space = Float32(ip - Int(p.nonstockable)) / Float32(ip)
    end
    # Growth calibration transition/measurement (GROWTH card: IDG/FINT/IHTG/FINTH/FINTM).
    # DG_TRANS=1 ⇒ the DG field is a PAST diameter (not an increment) measured DG_MEASURE yrs ago.
    _fia_present(d, "DG_TRANS")     && (c.growth_idg   = Int32(_fia_int(d, "DG_TRANS", 0)))
    # dbsstandin.f:700-703 / :709-711: DG_MEASURE / HTG_MEASURE replace FINT / FINTH only when > 0 (else the
    # variant default stays); MORT_MEASURE always replaces FINTM but ≤ 0 ⇒ 5 (:715-718).
    if _fia_present(d, "DG_MEASURE")
        v = _fia_f32(d, "DG_MEASURE", 0f0); v > 0f0 && (c.growth_fint = v)
        c.dbs_ifint = trunc(Int32, c.growth_fint)          # dbsstandin.f:702 IFINT = IFIX(FINT)
    end
    # The FIA-DB DG_TRANS/DG_MEASURE pair IS a GROWTH card — mark growth_dg_set so the DG calibration
    # NORMALIZES the observed increment by YR/FINT (simulate.jl:47 gates dgscale on growth_dg_set). Without it,
    # a non-native FINT (e.g. the 9-yr FIA remeasurement) is NOT scaled ⇒ the DGSCOR self-calibration over-fits
    # (loblolly COR 0.98→0.34, matching FVS's fort.13 raw scale 1.411 once set). Only when a measured-DG col present.
    (_fia_present(d, "DG_TRANS") || _fia_present(d, "DG_MEASURE")) && (c.growth_dg_set = true)
    _fia_present(d, "HTG_TRANS")    && (c.growth_ihtg  = Int32(_fia_int(d, "HTG_TRANS", 0)))
    if _fia_present(d, "HTG_MEASURE")
        v = _fia_f32(d, "HTG_MEASURE", 0f0); v > 0f0 && (c.growth_finth = v)
        c.growth_ifinth = trunc(Int32, c.growth_finth)          # dbsstandin.f:711 IFINTH = IFIX(FINTH)
    end
    if _fia_present(d, "MORT_MEASURE")
        v = _fia_f32(d, "MORT_MEASURE", 5f0); c.growth_fintm = v <= 0f0 ? 5f0 : v
    end
    # SITE_SPECIES (ISISP) + SITE_INDEX (SITEAR): assign to the site species only if given, else to all
    # species (dbsstandin.f:841). A SITE_INDEX ≤ 7 is a DUNNING site-CLASS code, NOT a site index in feet
    # (dbsstandin.f:763 `IF (RSTANDDATA(35).LE.7.) ... DUNNING CODE`); FVS calls DUNN to convert it, but the
    # SN/NE/CS/LS DUNN routines are all DUMMIES ⇒ the code is not usable as an SI ⇒ FVS falls through to the
    # SITSET default. So we record the site species but leave sp_site_index=0 for SITSET to default; using the
    # code literally (e.g. 5) cripples growth (frozen TopHt / suppressed DBH — audit 43bl).
    if _fia_present(d, "SITE_INDEX")
        si = _fia_f32(d, "SITE_INDEX", 0f0)
        # ON (canada/on dbsstandin.f:419-422, METRIC reader): RSTANDDATA(35)=SITE_INDEX*MtoFt — the DB site index is
        # metres and is converted BEFORE the ≤7 Dunning test. (The inline .tre/SITECODE path is converted in initre.)
        s.variant isa Ontario && (si *= 3.28084f0)
        isp = 0
        if _fia_present(d, "SITE_SPECIES")
            code = _fia_spcode(_fia_str(d, "SITE_SPECIES", ""))
            if !isempty(code)
                # FVS matches the SITE species STRICTLY against the variant's OWN
                # alpha/FIA/PLANTS codes only (dbsstandin.f:741-748) — NOT the SPCTRN
                # regional crosswalk that resolve_species falls through to for trees.
                # An unrecognized site species ⇒ ISISP=0 ⇒ SITE_INDEX applied to ALL
                # species (dbsstandin.f:776-779), which the `isp<1` branch below does.
                # NUMERIC codes: `_fia_spcode` zero-pads the input to 3 chars ("93"→"093")
                # but the variant coef stores FIA codes UNPADDED ("93"), so a strict `==`
                # here MISSED a 2-digit FIA site species (e.g. ES=093 vs table "93") ⇒ ISISP
                # fell to 0 ⇒ SITE_INDEX wrongly filled ALL species, freezing SITSET's
                # per-species interpolation (aspen SITEAR 86 vs live's interpolated 60.67 ⇒
                # RSIMOD 0.9 vs 0.719 ⇒ TT small-tree aspen height over-grows in dense
                # woodland). Normalize leading zeros the same way resolve_species does for
                # trees (but WITHOUT its SPCTRN fallthrough — the strict own-code match stands).
                cc = uppercase(code); sp = s.species
                ccn = all(isdigit, cc) ? lstrip(cc, '0') : ""
                @inbounds for j in 1:nspecies(s.variant)
                    fj = strip(sp.fia[j])
                    if strip(sp.alpha[j]) == cc || fj == cc || strip(sp.plants[j]) == cc ||
                       (!isempty(ccn) && all(isdigit, fj) && lstrip(fj, '0') == ccn)
                        isp = j; break
                    end
                end
            end
        end
        if si > 7f0                                # a real site index in feet
            if isp >= 1
                p.sp_site_index[isp] = si; p.site_species = Int32(isp); p.site_index = si
            else
                fill!(p.sp_site_index, si); p.site_index = si
            end
        else                                       # Dunning class code (≤7): record site species, SITSET defaults SI
            isp >= 1 && (p.site_species = Int32(isp))
        end
    end
    # FFE initial dead surface fuels (FUINI override) from the FVS_STANDINIT FUEL_* columns
    # (dbsstandin.f:396-458 → FUELINIT/FUELSOFT → fmcba.f:318-362). Measured FIA down-woody-material tons/ac by
    # size class; missing ⇒ -1 sentinel (fmcba! keeps the FUINI-table default for that class). Class order =
    # jl stfuel_hard (state.jl / _fuelinit!): <.25/.25-1/1-3/3-6/6-12/12-20/20-35/35-50/>50/litter/duff.
    _fuelcol(nm::Vararg{String}) = begin
        v = -1f0
        for n in nm; _fia_present(d, n) && (v = _fia_f32(d, n, -1f0); break); end
        v
    end
    hard = Float32[
        _fuelcol("FUEL_0_25_H","FUEL_0_25"),   _fuelcol("FUEL_25_1_H","FUEL_25_1"),
        _fuelcol("FUEL_1_3_H","FUEL_1_3"),     _fuelcol("FUEL_3_6_H","FUEL_3_6"),
        _fuelcol("FUEL_6_12_H","FUEL_6_12"),   _fuelcol("FUEL_12_20_H","FUEL_12_20","FUEL_GT_12"),
        _fuelcol("FUEL_20_35_H","FUEL_20_35"), _fuelcol("FUEL_35_50_H","FUEL_35_50"),
        _fuelcol("FUEL_GT_50_H","FUEL_GT_50"), _fuelcol("FUEL_LITTER"), _fuelcol("FUEL_DUFF")]
    # lumped <1" (FUEL_0_1) splits into classes 1&2 when the split columns are absent (fmcba.f:329-340)
    c01 = _fuelcol("FUEL_0_1","FUEL_0_1_H")
    if c01 >= 0f0
        if hard[1] < 0f0 && hard[2] < 0f0
            hard[1] = c01 * 0.5f0; hard[2] = c01 * 0.5f0
        elseif hard[1] < 0f0
            hard[1] = max(c01 - hard[2], 0f0)
        elseif hard[2] < 0f0
            hard[2] = max(c01 - hard[1], 0f0)
        end
    end
    soft = Float32[
        _fuelcol("FUEL_0_25_S"), _fuelcol("FUEL_25_1_S"), _fuelcol("FUEL_1_3_S"),
        _fuelcol("FUEL_3_6_S"), _fuelcol("FUEL_6_12_S"), _fuelcol("FUEL_12_20_S"),
        _fuelcol("FUEL_20_35_S"), _fuelcol("FUEL_35_50_S"), _fuelcol("FUEL_GT_50_S"), -1f0, -1f0]
    any(x -> x >= 0f0, hard) && (p.ffe_fuel_hard = hard)
    any(x -> x >= 0f0, soft) && (p.ffe_fuel_soft = soft)
    # FVS_StandInit FUEL_MODEL (dbsstandin.f:460): a fixed surface fuel model that FVS injects as a FUELMODL
    # activity (fmusrfm.f LUSRFM), overriding FMCFMD auto-selection. 0/missing ⇒ auto-select (unchanged).
    if _fia_present(d, "FUEL_MODEL")
        fm = _fia_int(d, "FUEL_MODEL", 0)
        fm > 0 && (p.fuel_model_input = Int32(fm))
    end
    return s
end

"""
    apply_fia_trees!(s, rows) -> Int

Map FVS_TREEINIT_COND/PLOT rows into `TreeRecord`s (dbstreesin.f column reads) and
ingest them via the shared `ingest_tree_records!` (same fixups as the .tre loader).
TREE_COUNT is stored raw as `tpa` (PROB); `notre!` later expands it to trees/acre.
"""
function apply_fia_trees!(s::StandState, rows::Vector{Dict{String,Any}})
    recs = TreeRecord[]
    for d in rows
        dbh = _fia_present(d, "DIAMETER") ? _fia_f32(d, "DIAMETER", 0f0) : _fia_f32(d, "DBH", 0f0)
        dmg = (Int32(_fia_int(d, "DAMAGE1", 0)), Int32(_fia_int(d, "SEVERITY1", 0)),
               Int32(_fia_int(d, "DAMAGE2", 0)), Int32(_fia_int(d, "SEVERITY2", 0)),
               Int32(_fia_int(d, "DAMAGE3", 0)), Int32(_fia_int(d, "SEVERITY3", 0)))
        rec = TreeRecord(
            Int32(_fia_int(d, "PLOT_ID", 1)),           # plot (ITREI) → subplot/IPVEC
            Int32(_fia_int(d, "TREE_ID", 0)),           # id (IDTREE)
            _fia_f32(d, "TREE_COUNT", 1f0),             # tpa (raw PROB; notre! expands)
            Int32(_fia_int(d, "HISTORY", 1)),           # history (ITH; 1 = live default)
            _fia_spcode(_fia_str(d, "SPECIES", "OT")),  # species_code (FIA 3-digit / alpha / PLANTS)
            dbh,                                         # dbh
            _fia_f32(d, "DG", 0f0),                     # diam_growth (PAST dbh when IDG=1)
            _fia_f32(d, "HT", 0f0),                     # height
            _fia_f32(d, "HTTOPK", 0f0),                 # top_height (broken/dead)
            _fia_f32(d, "HTG", 0f0),                    # ht_growth
            Int32(_fia_int(d, "CRRATIO", 0)),           # crown_pct (ICR)
            dmg,
            Int32(_fia_int(d, "TREEVALUE", 0)),         # mort_code (IMC1)
            Int32(_fia_int(d, "PRESCRIPTION", 0)),      # cut_code (KUTKOD)
            (Int32(0), Int32(0), Int32(0), Int32(0), Int32(0)),  # pest_vars
            _fia_f32(d, "AGE", 0f0),                    # birth_age (ABIRTH)
        )
        push!(recs, rec)
    end
    # Per-plot topography (PSLO/PASP, esplt2.f): each FVS_TREEINIT row carries its plot's SLOPE(%)/ASPECT(deg).
    # Keyed by the RAW DB PLOT_ID here; ingest_tree_records! REMAPS plot numbers to a 1..NPTS internal index
    # (first-appearance order) and records the raw→internal map in p.point_ids. So collect raw per-plot topo now,
    # then reindex to the internal point order AFTER ingest. A missing/NULL SLOPE ⇒ 0 (matches live's PSLO default)
    # — the #143 slope-dependent under-establishment root (jl formerly used the uniform stand slope for all points).
    raw_slo = Dict{Int,Float32}(); raw_asp = Dict{Int,Float32}()
    for d in rows
        pid = _fia_int(d, "PLOT_ID", 1)
        haskey(raw_slo, pid) && continue
        raw_slo[pid] = _fia_present(d, "SLOPE")  ? _fia_f32(d, "SLOPE", 0f0) * 0.01f0     : 0f0
        raw_asp[pid] = _fia_present(d, "ASPECT") ? _fia_f32(d, "ASPECT", 0f0) * 0.0174533f0 : 0f0
    end
    # Metric-variant DATABASE input (BC + ON): the FVS_TreeInit DB is METRIC (cm DBH, m HT, trees/ha).
    # Convert cm→in / m→ft on ingest exactly as the inline/.tre path (treeinput.jl:82); the trees/ha→
    # trees/acre expansion is applied just below (metric-DB only). Gated on the metric variants ⇒ US-FIA
    # sweeps unchanged. (ON was missing here — an Ontario DB previously ingested cm as inches, ~2.5× off,
    # the same class as the BC bug 42eb555.)
    metric_db = s.variant isa BritishColumbia || s.variant isa Ontario
    res = ingest_tree_records!(s, recs; metric = metric_db)
    if metric_db && !(s.variant isa Ontario)     # ON: the design factors are converted instead (see BASAL_AREA_FACTOR)
        # DB TREE_COUNT (PROB) is per-HECTARE; FVS's metric expansion yields per-acre (via the metric
        # plot area). notre! here multiplies by the design factor only, so pre-scale the raw PROB
        # per-ha→per-acre (× ACRtoHA) so the expanded internal TPA is per-acre like every other variant.
        @inbounds for i in 1:s.trees.n
            s.trees.tpa[i] *= 0.40468564f0
        end
    end
    p = s.plot
    npt_trees = s.trees.n > 0 ? maximum(Int(pp) for pp in @view s.trees.plot_id[1:s.trees.n]) : 0
    # ★ Size point_slope/point_aspect to the FULL inventory-point count (IPTINV-NONSTK = nptids) — but ONLY when
    # the stand carries ≥1 tree record. esplt2.f keys the per-plot slope on IPINFO:
    #   • IPINFO=2 (tree records present, carrying plot site data): FVS builds all IPTINV points; a point with NO
    #     tree record keeps PSLO(NNID)=0 (BLKDAT default → *0.01 → 0), while the tree-bearing points get their slope.
    #     jl formerly sized this to max(tree plot_id), so those empty points fell OUT of the array and
    #     ie_autoes_tally/ie_autoes_run used the empty→STAND-slope fallback — feeding the stand slope into ESTPP/
    #     ESTOCK for points the oracle drives with slope 0. Since ESTPP's BB carries −2.76·SLO, a spurious nonzero
    #     slope on the empty points UNDER-counts trees-per-plot (ITPP) ⇒ AUTOES ingrowth under-produces (MEASURED
    #     FVSie_g16 4733242010690 @2013: oracle per-point PSLO=[0.68,0,0,0], ITPP 92 → 100 TPA; jl fed 0.68 to all
    #     4 points ⇒ ITPP 62 → 64 TPA). Extending to npt_inv with empty points → 0 (get(...,0f0)) matches PSLO=0.
    #   • IPINFO=0 (fully BARE stand, ZERO tree records): FVS takes the NM=1/GOTO 150 branch and sets EVERY point to
    #     the STAND slope XXSLP (MEASURED FVSie_g16 1856009969290487, 0 trees, 4 points: PSLO=[0.44,0.44,0.44,0.44]).
    #     So for a bare stand leave point_slope EMPTY (npt=0) — ie_autoes_tally's empty→stand-slope fallback then
    #     drives every point with the stand slope, matching IPINFO=0. (Forcing empty points to 0 here OVER-produced:
    #     1856009969290487 189→470, 1629553740290487 252→826.)
    #   The IPINFO discriminator is the RAW tree-record count (length(recs)), NOT s.trees.n: a stand can carry a
    #   tree record that ingest DROPS (e.g. a HISTORY=6 dead / non-projectable stem) yet FVS still processed that
    #   record in esplt1 (IPINFO=2) — so its plot keeps its slope and the empty points get 0. A truly-bare stand
    #   (zero tree rows) is the only IPINFO=0 case. (MEASURED: 4733242010690 has 1 DB tree row but s.trees.n=0.)
    npt_inv = max(1, Int(p.points_inv) - Int(p.nonstockable))          # IPTINV-NONSTK (esplt2.f nptids)
    npt = length(recs) == 0 ? npt_trees : max(npt_trees, npt_inv)
    # esplt1.f:68 IPINFO=2 for tree records carrying plot site data ⇒ esplt2.f:274 LOAD=1 (the establishment site
    # prep comes from the plot data, not ESPREP's default proportions) — MEASURED FVSem_g16 196378260020004 thinbba
    # disturbance tally: "PNONE= 0 PMECH= 0 PBURN= 0", no ESPREP draw, every plot IPREP 1.
    length(recs) > 0 && (s.estab.load = true)
    # Always populate for DATABASE input (this reader is DATABASE-only; TREEDATA uses a different path so iet01 etc.
    # keep the empty→stand-slope fallback). A NULL DB slope ⇒ 0, which is the CORRECT establishment slope live uses.
    if npt >= 1
        ps = zeros(Float32, npt); pa = zeros(Float32, npt)
        @inbounds for k in 1:npt
            rawp = k <= length(p.point_ids) ? Int(p.point_ids[k]) : k    # internal k → raw plot number (IPVEC)
            ps[k] = get(raw_slo, rawp, 0f0); pa[k] = get(raw_asp, rawp, 0f0)
        end
        p.point_slope = ps; p.point_aspect = pa
    end
    return res
end

"""
    load_fia_stand!(s, dbpath, standsql, treesql) -> StandState

Populate stand `s` from an FIA "FVS-ready" SQLite database: run `standsql` (one row →
stand/plot state) and `treesql` (tree records), substituting `%StandID%` with the
stand's id (from STDIDENT). Mirrors the FVS DATABASE/DSNIN input block.
"""
function load_fia_stand!(s::StandState, dbpath::AbstractString,
                         standsql::AbstractString, treesql::AbstractString)
    sid = String(strip(s.plot.stand_id))
    # Open READ-ONLY (URI mode=ro, immutable=1): FVSjl only ever SELECTs from an FIA
    # database — it must never create/modify/journal the source file.
    db = SQLite.DB(startswith(dbpath, "file:") ? dbpath : "file:$(dbpath)?mode=ro&immutable=1")
    try
        srows = _fia_rows(db, standsql, sid)
        isempty(srows) && error("FIA database: no FVS_STANDINIT row for stand '$sid'")
        apply_fia_stand!(s, srows[1])
        isempty(treesql) || apply_fia_trees!(s, _fia_rows(db, treesql, sid))
    finally
        SQLite.close(db)
    end
    return s
end
