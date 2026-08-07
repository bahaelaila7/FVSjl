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
    # LOCATION (KODFOR): direct if present, else composite REGION*100 + FOREST (dbsstandin.f:569)
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
    # SLOPE: grinit.f:226 defaults a MISSING/NULL slope to 5.0% (→0.05 fraction) BEFORE the DB
    # overrides it — all 4 variants (sn/ne/cs/grinit.f:221-226 SLOPE=5.0). jl previously left a
    # missing slope at the 0.0 constructor default, which zeroed the DGF slope/aspect DGCON term
    # (TANS·SLOPE + FCOS·SLOPE·cos(ASP) + FSIN·SLOPE·sin(ASP), dgf.f:1125-1127) → over-grew species
    # with large slope coefficients (e.g. sp39 loblolly-bay FCOS=-10.15: a 0.05 slope = -0.68 in
    # ln(DDS) ⇒ ~2× DBH growth). Apply the grinit default so a missing slope matches live FVS.
    p.slope = _fia_present(d, "SLOPE") ? _fia_f32(d, "SLOPE", 0f0) / 100f0 : 5f0 / 100f0
    # ELEVATION in hundreds of feet; ELEVFT is feet → ×0.01 (dbsstandin.f:710)
    if _fia_present(d, "ELEVFT")
        p.elevation = _fia_f32(d, "ELEVFT", 0f0) * 0.01f0
    elseif _fia_present(d, "ELEVATION")
        p.elevation = _fia_f32(d, "ELEVATION", p.elevation)
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
    # CI: the habitat KODTYP fed to ci_habtyp is the 3-digit NI code in PV_REF_CODE (e.g. 401); PV_CODE holds the
    # 5-digit FIA code (out of ci_habtyp's 10-999 range). Without it habitat_code=0 ⇒ habitat_input defaults to 1
    # ⇒ wrong DG (DGHAB via ICINDX) + wrong mortality ITYPE ⇒ multi-cycle divergence.
    # FVS uses PV_CODE (the FIA habitat code) with PRIORITY over PV_REF_CODE — the .out prints "PV_CODE WAS
    # USED, PV_REF_CODE WAS IGNORED". A 5-digit PV_CODE (e.g. 41732) is the 2-digit state prefix + 3-digit
    # habitat (→ 732 = 41732 mod 1000); ci_habtyp then maps 732→(ICINDX 110, ITYPE 27). PV_REF_CODE (401) is
    # only the fallback. Using PV_REF_CODE first (jl's old behavior) picked the WRONG habitat ⇒ wrong DGHAB
    # (dg_const) + ITYPE (htgf) + mortality on stands where the two codes differ.
    if s.variant isa CentralIdaho
        hc = 0
        if _fia_present(d, "PV_CODE")
            pvc = Int(round(_fia_f32(d, "PV_CODE", 0f0)))
            pvc > 999 && (pvc = pvc % 1000)               # strip the 2-digit state prefix (41732 → 732)
            (10 <= pvc <= 999) && (hc = pvc)
        end
        if hc == 0 && _fia_present(d, "PV_REF_CODE")       # fallback
            pvr = Int(round(_fia_f32(d, "PV_REF_CODE", 0f0)))
            (10 <= pvr <= 999) && (hc = pvr)
        end
        hc != 0 && (p.habitat_code = Int32(hc))
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
        if _fia_present(d, "PV_CODE")
            pvc = Int(round(_fia_f32(d, "PV_CODE", 0f0)))
            pvc > 999 && (pvc = pvc % 1000)               # strip the 2-digit state prefix (41780 → 780)
            (10 <= pvc <= 999) && (hc = pvc)
        end
        if hc == 0 && _fia_present(d, "PV_REF_CODE")       # fallback
            pvr = Int(round(_fia_f32(d, "PV_REF_CODE", 0f0)))
            (10 <= pvr <= 999) && (hc = pvr)
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
    if s.variant isa Southern
        lat0, long0, elev0 = forest_location(s.coef, div(Int(p.user_forest_code), 100))
        p.latitude  == 0f0 && (p.latitude  = lat0)
        p.longitude == 0f0 && (p.longitude = long0)
        p.elevation == 0f0 && (p.elevation = elev0)
    end
    # Sampling design (DESIGN card): BAF / FPA / BRK / IPTINV / NONSTK / SAMWT / GROSPC
    _fia_present(d, "BASAL_AREA_FACTOR") && (p.baf = _fia_f32(d, "BASAL_AREA_FACTOR", 0f0))
    _fia_present(d, "INV_PLOT_SIZE")     && (p.fixed_plot_inv = _fia_f32(d, "INV_PLOT_SIZE", 0f0))
    _fia_present(d, "BRK_DBH")           && (p.min_dbh_var_plot = _fia_f32(d, "BRK_DBH", p.min_dbh_var_plot))
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
    _fia_present(d, "DG_MEASURE")   && (c.growth_fint  = _fia_f32(d, "DG_MEASURE", 5f0))
    # The FIA-DB DG_TRANS/DG_MEASURE pair IS a GROWTH card — mark growth_dg_set so the DG calibration
    # NORMALIZES the observed increment by YR/FINT (simulate.jl:47 gates dgscale on growth_dg_set). Without it,
    # a non-native FINT (e.g. the 9-yr FIA remeasurement) is NOT scaled ⇒ the DGSCOR self-calibration over-fits
    # (loblolly COR 0.98→0.34, matching FVS's fort.13 raw scale 1.411 once set). Only when a measured-DG col present.
    (_fia_present(d, "DG_TRANS") || _fia_present(d, "DG_MEASURE")) && (c.growth_dg_set = true)
    _fia_present(d, "HTG_TRANS")    && (c.growth_ihtg  = Int32(_fia_int(d, "HTG_TRANS", 0)))
    _fia_present(d, "HTG_MEASURE")  && (c.growth_finth = _fia_f32(d, "HTG_MEASURE", 5f0))
    _fia_present(d, "MORT_MEASURE") && (c.growth_fintm = _fia_f32(d, "MORT_MEASURE", 5f0))
    # SITE_SPECIES (ISISP) + SITE_INDEX (SITEAR): assign to the site species only if given, else to all
    # species (dbsstandin.f:841). A SITE_INDEX ≤ 7 is a DUNNING site-CLASS code, NOT a site index in feet
    # (dbsstandin.f:763 `IF (RSTANDDATA(35).LE.7.) ... DUNNING CODE`); FVS calls DUNN to convert it, but the
    # SN/NE/CS/LS DUNN routines are all DUMMIES ⇒ the code is not usable as an SI ⇒ FVS falls through to the
    # SITSET default. So we record the site species but leave sp_site_index=0 for SITSET to default; using the
    # code literally (e.g. 5) cripples growth (frozen TopHt / suppressed DBH — audit 43bl).
    if _fia_present(d, "SITE_INDEX")
        si = _fia_f32(d, "SITE_INDEX", 0f0)
        isp = 0
        if _fia_present(d, "SITE_SPECIES")
            code = _fia_spcode(_fia_str(d, "SITE_SPECIES", ""))
            if !isempty(code)
                # FVS matches the SITE species STRICTLY against the variant's OWN
                # alpha/FIA/PLANTS codes only (dbsstandin.f:741-748) — NOT the SPCTRN
                # regional crosswalk that resolve_species falls through to for trees.
                # An unrecognized site species ⇒ ISISP=0 ⇒ SITE_INDEX applied to ALL
                # species (dbsstandin.f:776-779), which the `isp<1` branch below does.
                cc = uppercase(code); sp = s.species
                @inbounds for j in 1:nspecies(s.variant)
                    if strip(sp.alpha[j]) == cc || strip(sp.fia[j]) == cc || strip(sp.plants[j]) == cc
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
    return ingest_tree_records!(s, recs)
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
