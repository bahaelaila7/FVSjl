"""
    FVSjl

An idiomatic, maintainable, thread-safe Julia reimplementation of the USFS Forest
Vegetation Simulator. It is a drop-in replacement for the live Fortran FVS — Southern
(`FVSsn`) and Northeast (`FVSne`) are complete and validated, Central States (`FVScs`)
in progress — same `.key`/`.tre` inputs, same SQLite / `.sum` outputs, and (in the
default `faithful=true` mode) bit-exact results.

Design (see docs/ARCHITECTURE.md):
  * all simulation state lives on an explicit `StandState` — no globals;
  * one state per stand/thread → safe parallelism;
  * pure numeric kernels + stateful orchestration → testable;
  * variants are dispatched via `AbstractVariant` singletons.

Ported from the Fortran sources under /workspace/ForestVegetationSimulator and
validated against the live Fortran builds (SN additionally vs the faithful port at
/workspace/FVSjulia; NE/CS have no faithful port — live binary is the sole oracle).
"""
module FVSjl

using Printf
using SQLite
using DBInterface

# --- core (order matters: parameters → rng/units/trees → variant → state) ----
include("core/fmath.jl")                       # gfortran-identical exp/log/pow (doctrine #8)
using .FMath: fexp, flog, fpow, fexp_julia, flog_julia, fpow_julia
include("core/parameters.jl")
include("core/rng.jl")
include("core/units.jl")
include("core/trees.jl")
include("variants/variant.jl")
include("core/coefficients.jl")
include("core/state.jl")

# --- variants ---------------------------------------------------------------
include("variants/southern/southern.jl")
include("variants/southern/coefficients.jl")
include("variants/southern/species.jl")
include("variants/southern/forest_location.jl")
include("variants/southern/habitat.jl")
include("variants/southern/site_index.jl")
include("variants/southern/bark_and_bounds.jl")
include("variants/southern/serial_correlation.jl")
include("variants/southern/diameter_growth.jl")
include("variants/southern/height_growth.jl")
include("variants/southern/small_tree_growth.jl")
include("variants/southern/crown_ratio.jl")
include("variants/southern/mortality.jl")

# CS singleton is registered BEFORE the NE methods because NE's shared TWIGS crown method
# dispatches on `Union{Northeast,CentralStates}` (CS reuses it). The CS *methods* land after.
include("variants/centralstates/centralstates.jl")  # CS singleton + registration (MAXSP 96)
include("variants/lakestates/lakestates.jl")         # LS singleton + registration (MAXSP 68)
include("variants/centralrockies/centralrockies.jl") # CR singleton + registration (MAXSP 38) — first WESTERN variant
include("variants/centralrockies/species.jl")
include("variants/centralrockies/site_index.jl")
include("variants/centralrockies/diameter_growth.jl")
include("variants/centralrockies/height_growth.jl")
include("variants/centralrockies/crown.jl")
include("variants/centralrockies/small_tree_growth.jl")
include("variants/centralrockies/mortality.jl")
include("variants/centralrockies/dwarf_mistletoe_model.jl")  # CR dwarf mistletoe per-cycle spread (mistoe.f)

# --- kootenai (KT) — first of the western Rockies cluster; scaffold, equations chunk by chunk ---
include("variants/kootenai/kootenai.jl")             # KT singleton + registration (MAXSP 11, western Wykoff DDS)
include("variants/kootenai/species.jl")              # KT blkdat init (11 species, seed 55329, Stage SDI) + SPCTRN col 12
include("variants/kootenai/habitat_tables.jl")       # KT habtyp/sitset lookup tables (KOTHAB/JTYPE/KTYPE/MTYPE/BAMAXA)
include("variants/kootenai/site_index.jl")           # KT habtyp 2-level mapping + sitset SDImax (site_setup!)
include("variants/kootenai/dg_coefficients.jl")      # KT large-tree DDS coefficient arrays (kt/dgf.f DATA)
include("variants/kootenai/crown.jl")                # KT per-tree CCF (ccfcal MODE=1) → RELDEN
include("variants/kootenai/diameter_growth.jl")      # KT dgf! (Wykoff DDS) + kt_dgcons!
include("variants/kootenai/height_growth.jl")        # KT htgf! (exp-form height increment)
include("variants/kootenai/regent.jl")               # KT regent (chunk 6): kt_regcons! + small_tree_growth!
include("variants/kootenai/mortality.jl")            # KT mortality (chunk 7): Hamilton MORTS
include("variants/kootenai/volume.jl")               # KT volume (chunk 8): Region-1 Flewelling FW2

# --- inland empire (IE) — 3rd western variant, reuses the KT engine at a discount ---
include("variants/britishcolumbia/britishcolumbia.jl") # BC singleton + registration (MAXSP 15) — chunk 0 scaffold
include("variants/britishcolumbia/species.jl")          # BC species block-data init (bc/blkdat.f) — chunk 1
include("variants/britishcolumbia/becset.jl")           # BC BEC-zone parse (habtyp) + sitset BAMAX table — chunk 2
include("variants/britishcolumbia/site_index.jl")       # BC site/habitat (chunk 2): KODTYP→ITYPE + Stage SDImax
include("variants/britishcolumbia/dg_coefficients.jl")  # BC V3 DG coeffs (ZNKONST/SSKONST) + matcher/DDS/bark — chunk 3
include("variants/britishcolumbia/diameter_growth.jl")  # BC large-tree DDS hooks (bc_dgcons! + dgf!) — chunk 3 (not yet in simulate dispatch)
include("variants/britishcolumbia/v2_diameter_growth.jl") # BC V2-regime large-tree DDS formula (ESSF/MS/PP) — chunk 2b (task #133)
include("variants/britishcolumbia/height_growth_coefficients.jl") # BC V3 HTG coeffs (LTHG/LTSP) + bc_v3_htg — chunk 4
include("variants/britishcolumbia/height_growth.jl")             # BC height_growth!(::BritishColumbia) — chunk 4 (V3)
include("variants/britishcolumbia/crown_coefficients.jl")       # BC V3 crown-ratio coeffs (CRKONST ICH/IDF) + CRNMD — chunk 5b
include("variants/britishcolumbia/crown.jl")                    # BC per-tree CCF + stand RELDEN (ccfcal.f) + crown_ratio_update! — chunk 5
include("variants/britishcolumbia/regent_coefficients.jl")     # BC small-tree V3 coeffs (ST_COEF) + bc_v3_sthg/bc_st_dbh — chunk 6
include("variants/britishcolumbia/regent.jl")                   # BC small_tree_growth!(::BritishColumbia) — chunk 6
include("variants/britishcolumbia/mortality.jl")                # BC mortality!(::BritishColumbia): V3 FMRT tabular + LMRT logistic — chunk 7
include("variants/britishcolumbia/volume.jl")                   # BC total cubic volume (Kozak taper: CFVOL→MIN→LOG) — chunk 8
include("variants/britishcolumbia/newspred_shd.jl")             # BC NEWSPRED Shd1/ShdPtr encoded trajectory table (#196 C4)
include("variants/britishcolumbia/newspred.jl")                 # BC NEWSPRED/NISI spatial dwarf-mistletoe (#196) — C1 state+consts (not yet engine-wired)
include("variants/centralidaho/centralidaho.jl")     # CI singleton + registration (MAXSP 19) — chunk 0 scaffold
include("variants/centralidaho/species.jl")          # CI species block-data init (ci/blkdat.f + grinit.f) — chunk 1
include("variants/centralidaho/site_index.jl")       # CI habtyp/forkod/sitset (chunk 2): ICINDX/ITYPE, R4SDI SDImax
include("variants/centralidaho/dg_coefficients.jl")  # CI large-tree DDS coefficient arrays (ci/dgf.f) — chunk 3
include("variants/centralidaho/diameter_growth.jl")  # CI large-tree DDS hooks (ci_dgcons! + dgf!) — chunk 3
include("variants/centralidaho/crown.jl")            # CI per-tree CCF (ci/ccfcal.f MODE=1) → RELDEN — chunk 5
include("variants/centralidaho/height_growth.jl")    # CI large-tree height growth (ci/htgf.f) — chunk 4
include("variants/centralidaho/regent.jl")           # CI small-tree growth (ci/regent.f) — chunk 6
include("variants/centralidaho/establishment.jl")    # CI ESSUBH planted/subsequent base height (ci/essubh.f) — #154
include("variants/centralidaho/mortality.jl")        # CI mortality (ci/morts.f) — Hamilton, chunk 7
include("variants/centralidaho/volume.jl")           # CI volume (ci VEQNNC): MATW/FW2W/DVEW — chunk 8
include("variants/inlandempire/inlandempire.jl")     # IE singleton + registration (MAXSP 23)
include("variants/inlandempire/species.jl")          # IE species block-data init (ie/blkdat.f + grinit.f)
include("variants/inlandempire/site_index.jl")       # IE habtyp/forkod/sitset (chunk 2): ITYPE/IFOR/SITEAR
include("variants/inlandempire/dg_coefficients.jl")  # IE large-tree DDS coefficient arrays (ie/dgf.f)
include("variants/inlandempire/diameter_growth.jl")  # IE large-tree DDS hooks (ie_dgcons! + dgf!)
include("variants/inlandempire/height_growth.jl")    # IE large-tree height growth (ie/htgf.f)
include("variants/inlandempire/crown_coefficients.jl") # IE crown-ratio tables (ie/crown.f, dumped from live)
include("variants/inlandempire/crown.jl")            # IE per-tree CCF (ie/ccfcal.f) + crown ratio (ie/crown.f)
include("variants/inlandempire/mort_coefficients.jl")# IE mortality Hamilton coefficients (ie/morts.f, dumped)
include("variants/inlandempire/mortality.jl")        # IE mortality (ie/morts.f) — Hamilton, reuses KT form
include("variants/inlandempire/mistoe_coefficients.jl") # IE MISTOE (dwarf mistletoe) effect coefficients (mistoe/misintie.f)
include("variants/inlandempire/pvref1_data.jl")        # IE PVREF1 (PVCODE,PVREF)->HABPVR crosswalk table (ie/pvref1.f)
include("variants/inlandempire/establishment.jl")      # IE ESSUBH subsequent/planted-tree height model (ie/essubh.f)
include("variants/inlandempire/volume.jl")           # IE volume (ie/sitset.f VOLEQ): FW2 + DVE, reuse CR kernels
include("variants/inlandempire/regent_coefficients.jl") # IE REGENT small-tree coefficients (ie/regent.f, dumped)
include("variants/inlandempire/regent.jl")           # IE REGENT small-tree growth (ie/regent.f): ie_regcons! + hook

# --- easternmontana (EM) — next western Rockies cluster variant (KT-engine discount); chunk 0 scaffold ---
include("variants/easternmontana/easternmontana.jl")  # EM singleton + registration (MAXSP 19, western Wykoff DDS)
include("variants/easternmontana/species.jl")         # EM blkdat init (19 species, seed 55329, Stage SDI) + SPCTRN col 10
include("variants/easternmontana/site_index.jl")      # EM habtyp (JTYPE/NIHMAP→ITYPE) + sitset (site index/SDIDEF) — site_setup!
include("variants/easternmontana/crown.jl")           # EM per-tree CCF (em/ccfcal.f MODE=1 polynomial) — em_tree_ccf
include("variants/easternmontana/dg_coefficients.jl") # EM large-tree DG coefficient arrays (em/dgf.f DATA, 1D+2D)
include("variants/easternmontana/diameter_growth.jl") # EM large-tree DDS (em/dgf.f): em_dgcons! + dgf! (Wykoff main)
include("variants/easternmontana/height_growth.jl") # EM large-tree height growth (em/htgf.f+pothtg.f): POTHTG+RALPH+modifier
include("variants/easternmontana/mortality.jl") # EM mortality (em/morts.f): RI-background + SDI-trend RN
include("variants/easternmontana/regent.jl") # EM small-tree growth (em/regent.f NIVAR): em_regcons! + small_tree_growth!
include("variants/easternmontana/volume.jl") # EM volume: FW2 conifers (cr_fw2_vol) + DVEW woodland (R1KEMP cubic)
include("variants/easternmontana/establishment.jl") # EM ESSUBH base height (em/essubh.f) — em_essubh_hht + em_ihtser

# --- teton (TT) — next western Rockies cluster variant (EM-DG + CR-Zeide-density discount) ---
include("variants/teton/teton.jl")            # TT singleton + registration (MAXSP 18, western Wykoff DDS, oracle verified)
include("variants/teton/species.jl")          # TT blkdat init (18 species, seed 55329, ZEIDE SDI) + SPCTRN
include("variants/teton/habtyp_table.jl")      # TT R4HABT(363) habitat code → ITYPE (tt/habtyp.f)
include("variants/teton/site_index.jl")        # TT sitset (SITEAR interp + SDIDEF) + forkod (chunk 2)
include("variants/teton/crown.jl")             # TT per-tree CCF (tt/ccfcal.f MODE=1) — RELDEN+PCCF for DG (chunk 3)
include("variants/teton/dg_coefficients.jl")   # TT large-tree DG coefficient arrays (tt/dgf.f DATA, generated)
include("variants/teton/diameter_growth.jl")   # TT large-tree DDS (tt/dgf.f): tt_dgcons! + dgf! (MAIN + ASPEN)
include("variants/teton/htgf_coefficients.jl") # TT Schreuder-Hafley SBB height coefficients (tt/htgf.f, generated)
include("variants/teton/height_growth.jl")     # TT large-tree height growth (tt/htgf.f DEFAULT SBB) — chunk 4
include("variants/teton/regent.jl")            # TT small-tree growth (tt/regent.f+smhtgf.f+smdgf.f) — chunk 6
include("variants/teton/mortality.jl")         # TT mortality (tt/morts.f = EM Hamilton form, TT coeffs) — chunk 7
include("variants/teton/volume.jl")             # TT volume (chunk 8): R4VOL Matney cubic + DVE
include("variants/teton/establishment.jl")      # TT establishment clamps (tt/blkdat.f ESCOMN XMIN/HHTMAX)

# --- utah (UT) — Region-4 western Wykoff DDS (Zeide SDI), 24 species, heavy PJ/woodland; chunk 0 scaffold ---
include("variants/utah/utah.jl")              # UT singleton + registration (MAXSP 24, western Wykoff DDS, Zeide SDI)
include("variants/utah/species.jl")           # UT blkdat init (24 species, seed 55329, Zeide SDI) + SPCTRN col 4
include("variants/utah/site_index.jl")         # UT site_setup! (chunk 2): forkod + habtyp (reuse TT R4HABT) + sitset
include("variants/utah/crown.jl")              # UT per-tree CCF (ut/ccfcal.f MODE=1) → RELDEN + PCCF (chunk 5 partial)
include("variants/utah/dg_coefficients.jl")    # UT large-tree DG coefficient arrays (ut/dgf.f DATA, CSV-loaded)
include("variants/utah/diameter_growth.jl")    # UT large-tree DDS (chunk 3): ut_dgcons! + dgf! (5-branch)
include("variants/utah/htgf_coefficients.jl")  # UT Schreuder-Hafley SBB height coefficients (ut/htgf.f, CSV)
include("variants/utah/height_growth.jl")      # UT large-tree height (chunk 4): SBB (conifer/aspen)
include("variants/utah/regent.jl")             # UT small-tree growth (chunk 6): ut/regent.f POTHTG + ht_dbh DG
include("variants/utah/mortality.jl")          # UT mortality (chunk 7): ut/morts.f uniform Hamilton + SDI self-thin
include("variants/utah/volume.jl")             # UT volume (chunk 8): MATW r4vol + FW2 + DVEW (reuse)
include("variants/utah/establishment.jl")      # UT ESSUBH planted/subsequent base height (ut/essubh.f) — #154

# --- bluemountains (BM) — Region-6 western Wykoff DDS (Stage SDI), 18 species; chunk 0 scaffold ---
include("variants/bluemountains/bluemountains.jl")  # BM singleton + registration (MAXSP 18, western Wykoff DDS)
include("variants/bluemountains/species.jl")        # BM blkdat init (18 species, seed 55329, Stage SDI, DGSD=1.5)
include("variants/bluemountains/site_index.jl")     # BM chunk 2: forkod/habtyp(PCOML)/ecocls/sichg/htcalc/sitset
include("variants/bluemountains/dg_coefficients.jl")# BM chunk 3: DG coefficient arrays + BM_PSIGSQ + bm_bratio
include("variants/bluemountains/crown.jl")          # BM chunk 5 (partial): bm_tree_ccf (CCF/RELDEN, needed by DG)
include("variants/bluemountains/diameter_growth.jl")# BM chunk 3: bm_dgcons! + dgf!(::BlueMountains) MSS spline
include("variants/bluemountains/height_growth.jl")  # BM chunk 4: bm_findag + height_growth!(::BlueMountains)
include("variants/bluemountains/regent.jl")         # BM chunk 6: bm_smhtgf + small_tree_growth!(::BlueMountains)
include("variants/bluemountains/mortality.jl")      # BM chunk 7: mortality!(::BlueMountains) Hamilton + SDI self-thin
include("variants/bluemountains/volume.jl")         # BM chunk 8: compute_volumes_bm! (FW2W Flewelling conifers)

# --- variants: northeast (NE) — skeleton; equations + data ported chunk by chunk ---
include("variants/northeast/northeast.jl")
include("variants/northeast/species.jl")
include("variants/northeast/site_index.jl")
include("variants/northeast/diameter_growth.jl")    # NE large-tree DG (B1/B2/B3 + BAL competition)
include("variants/northeast/height_growth.jl")      # NE height growth (NC-128 curve + BAL modifier)
include("variants/northeast/mortality.jl")          # NE VARMRT efficiency (relative-height + VARADJ)
include("variants/northeast/crown_ratio.jl")        # NE crown ratio (TWIGS BA-based model)
include("variants/northeast/small_tree_growth.jl")  # NE small-tree growth (REGENT, dbh<5)

# --- centralstates (CS) methods (singleton already registered above, before NE) -------------
include("variants/centralstates/species.jl")        # CS blkdat init (codes/RNG/YR=10/Zeide SDI)
include("variants/centralstates/site_index.jl")     # CS SITSET (per-species ASITE/BSITE linear)
include("variants/centralstates/diameter_growth.jl")# CS DG: cs_dgf! ln(DDS) regression + cs_dgcons!
include("variants/centralstates/height_growth.jl")  # CS HTGF: NC-128 (MAPCS) + cs_balmod (reuses LTBHEC)
include("variants/centralstates/small_tree_growth.jl") # CS REGENT (d<5): NE shape, XMIN=3, cs_balmod
# CS HT-DBH rides the Southern Curtis-Arney+Wykoff path (cs/htdbh.f ≡ sn/htdbh.f); CS crown
# rides NE's TWIGS method (see northeast/crown_ratio.jl) — both are coefficient-driven, no CS code.

# --- lakestates (LS) methods (singleton already registered above, before NE) -----------------
include("variants/lakestates/species.jl")           # LS blkdat init (codes/RNG/YR=10/Zeide SDI) + spctrn col 5
include("variants/lakestates/site_index.jl")         # LS SITSET (68×68 SICOEF fan-out + SDICON) + FORKOD
include("variants/lakestates/diameter_growth.jl")    # LS DG: ls_dgf! ln(DDS) (== CS model, LS caps) + ls_dgcons!
include("variants/lakestates/height_growth.jl")      # LS HTGF: NC-128 (MAPLS/IVAR=1) + ls_balmod competition
include("variants/lakestates/small_tree_growth.jl")  # LS REGENT (d<5): NE/CS shape, XMIN=3, MAPLS + ls_balmod
include("variants/klamath/klamath.jl")               # NC (Klamath Mtns) singleton + registration (MAXSP 12, Zeide) — chunk 0 scaffold
include("variants/klamath/species.jl")               # NC species block-data init (nc/blkdat.f + grinit.f) — chunk 1
include("variants/klamath/site_index.jl")            # NC site index + SDImax (nc/sitcind.f + ecocls.f) — chunk 2
include("variants/klamath/diameter_growth.jl")       # NC large-tree DDS (nc/dgf.f): nc_bratio + nc_dgcons! + dgf! — chunk 3
include("variants/klamath/height_growth.jl")         # NC height growth (nc/htgf.f+findag.f+htcalc.f) — chunk 4
include("variants/klamath/regent.jl")                # NC small-tree growth (nc/regent.f+htgr5.f+htdbh.f) — chunk 6
include("variants/klamath/mortality.jl")             # NC mortality (nc/morts.f — reuses EM/UT Zeide form) — chunk 7
include("variants/klamath/crown.jl")                 # NC crown ratio (nc/crown.f — Weibull) — chunk 5
include("variants/klamath/volume.jl")                # NC volume (nc VEQNNC): WO2W R5TAP taper + DVEW r5harv — chunk 8
include("variants/oregoncoast/oregoncoast.jl")       # OC (Oregon Coast) singleton + registration (MAXSP 50, ORGANON SWO) — chunk 0 foundation
include("variants/oregoncoast/species.jl")           # OC species block-data init (oc/blkdat.f + grinit.f) — chunk 1; growth = UNPORTED ORGANON
include("variants/oregoncoast/organon_interface.jl") # OC FVS↔ORGANON boundary marshalling (orgspc/IORG/big-6/buffer) — chunk C1
include("variants/oregoncoast/organon_setup.jl")     # OC ORGANON PREPARE setup calibration (ACALIB/TMPCAL + HT/CR dubbing, SWO) — chunk C2
include("variants/oregoncoast/organon_diamgro.jl")   # OC ORGANON SWO diameter growth (DG_SWO/DIAMGRO_RUN + SSTATS/GET_BAL/SUBMAX/DGCALIB) — chunk C3
include("variants/oregoncoast/organon_htgro.jl")     # OC ORGANON SWO height growth (HG_SWO/HTGRO1 + HS_HG/CRNCLO/LCW/HLCW/CW/LIMIT) — chunk C4
include("variants/oregoncoast/organon_crngro.jl")    # OC ORGANON SWO crown growth (CROWGRO + HCB_SWO/MAXHCB_SWO/OLDGRO) — chunk C5
include("variants/oregoncoast/organon_mortality.jl") # OC ORGANON SWO mortality (MORTAL_RUN + PM_SWO/PM_FERT/QUAD1) — chunk C6
include("variants/oregoncoast/organon_execute.jl")   # OC ORGANON EXECUTE/GROW orchestration + FVS bark/DDS copy-back — chunk C7
include("variants/oregoncoast/organon_hook.jl")      # OC ORGANON live growth hook + StandState copy-back (diameter_growth!) — chunk C7
include("variants/oregoncoast/organon_cratet.jl")    # OC setup height dubbing (oc/cratet.f HTDBH Curtis-Arney) — chunk C8
include("variants/oregoncoast/organon_dgf.jl")       # OC FVS-native Wykoff DGF for IORG=0 trees (oc/dgf.f + DGCONS + forkod) — chunk C9
include("variants/oregoncoast/organon_htgf.jl")      # OC FVS-native HTGF for IORG=0 trees (oc/htgf.f + findag + htcalc) — chunk C9 step 4
include("variants/oregoncoast/organon_volume.jl")    # OC BLM cubic volume (blmvol/blmtap/formcl BLM711) — chunk C10a

# --- southeastalaska (AK) — Region-10 Wykoff DDS (Zeide SDI), 23 species, permafrost DG modifier; chunk 0 beachhead ---
include("variants/southeastalaska/southeastalaska.jl") # AK singleton + registration (MAXSP 23, Wykoff DDS, Zeide SDI)
include("variants/southeastalaska/species.jl")         # AK blkdat init (23 species, seed 55329, Zeide SDI, LHTDRG=false) + SPCTRN col 4
include("variants/southeastalaska/dg_coefficients.jl") # AK large-tree DG + permafrost + bark coefficient arrays (ak/dgf.f, ak/bratio.f)
include("variants/southeastalaska/diameter_growth.jl") # AK large-tree DDS (chunk 3): ak_bratio + ak_dgcons! + dgf! — VALIDATED bit-exact
include("variants/southeastalaska/site_index.jl")      # AK site index (ak/sitset.f SITEAR) + SDImax + forkod — chunk 2 (SITEAR validated)
include("variants/southeastalaska/height_growth.jl")   # AK large-tree height growth (ak/htgf.f) — chunk 4 (validated bit-exact cyc0)
include("variants/southeastalaska/crown.jl")           # AK crown ratio (ak/crown.f logistic + dubscr) + point-Zeide (sdical) — chunk 5
include("variants/southeastalaska/regent.jl")          # AK small-tree growth (ak/regent.f) — chunk 6 STUB (no-op)
include("variants/southeastalaska/mortality.jl")       # AK mortality (ak/morts.f logistic survival + SDI/BA iterative pass) — chunk 7
include("../data/southeastalaska/volume_coefficients.jl") # AK R10 F32 Flewelling profile coeffs (SHP_AK/FDBT_AK)
include("variants/southeastalaska/volume.jl")          # AK Region-10 volume (VOLEQDEF→NVEL): F32 Flewelling (chunk 8)

# --- westcascades (WC) — R6 westside Prognosis (Reineke SDI), 39 species; chunk 0 foundation + DGF validated ---
include("variants/westcascades/westcascades.jl")     # WC (West Cascades) singleton + registration (MAXSP 39, Reineke) — chunk 0 foundation
include("variants/westcascades/species.jl")          # WC species-coefficient table binding + blkdat init (bark/site/crown/vol) — chunk 1
include("variants/westcascades/diameter_growth.jl")  # WC large-tree DDS (wc/dgf.f): wc_dgcons! + dgf! (19-group Wykoff) — chunk 3 (DGF validated)
include("variants/westcascades/height_growth.jl")    # WC large-tree HTG (wc/htgf.f + findag.f + htcalc.f): wc_findag/wc_htcalc + height_growth! — chunk 4
include("variants/westcascades/site_index.jl")       # WC site index + Reineke SDImax (forkod/habtyp/ecocls/sichg/sitset) — chunk 2
include("variants/westcascades/crown.jl")            # WC crown ratio (wc/crown.f + wc/dubscr.f): Weibull CR + RW logistic + DUBSCR — chunk 5
include("variants/westcascades/regent.jl")           # WC small-tree growth (wc/regent.f + vwc/smhgdg.f + wc/htdbh.f + dgbnd.f): SMHGDG + small_tree_growth! — chunk 6
include("variants/westcascades/volume.jl")           # WC volume (R6 NVEL): westside Flewelling SHP_W3/W4/W5 + INGY + Behre + wc_formcl — chunk 8
include("variants/westcascades/mortality.jl")        # WC mortality — chunk 7 PLACEHOLDER no-op (vwc/morts.f ORGANON RIP not yet ported; unblocks cyc0 .sum)

# --- pacificnorthwest (PN) — R6 westside Prognosis near-clone of WC (coefficient swap on the shared engine) ---
include("variants/pacificnorthwest/pacificnorthwest.jl")  # PN singleton + registration (MAXSP 39, Reineke, slot-6 SS) — chunk 0
include("variants/pacificnorthwest/species.jl")           # PN species-coefficient table binding + blkdat init — chunk 1
include("variants/pacificnorthwest/diameter_growth.jl")   # PN large-tree DDS (pn/dgf.f): pn_dgcons! + dgf! (20-group Wykoff) — chunk 3
include("variants/pacificnorthwest/height_growth.jl")     # PN HTG (pn/htcalc.f Farr SS/RC + King DF; findag/htgf = WC) — chunk 4
include("variants/pacificnorthwest/site_index.jl")        # PN site (pn/forkod/habtyp/ecocls/sichg/sitset) — chunk 2
include("variants/pacificnorthwest/crown.jl")             # PN crown (pn/crown.f 17-group Weibull + dubscr) + CCF (pn/ccfcal.f 19-group) — chunk 5 + density
include("variants/pacificnorthwest/regent.jl")            # PN small-tree growth (WC SMHGDG no-King + PN htdbh) — chunk 6
include("variants/pacificnorthwest/mortality.jl")         # PN mortality (WC ORGANON arrays, raw SITEAR) — chunk 7
include("variants/pacificnorthwest/volume.jl")            # PN volume (WC westside F00 + PN Behre/formcl) — chunk 8

# --- eastcascades (EC) — R6 westside Wykoff-DDS variant (32-species, own uncompressed coefficient DATA) ---
include("variants/eastcascades/eastcascades.jl")          # EC singleton + registration (MAXSP 32, Reineke) — chunk 0
include("variants/eastcascades/species.jl")               # EC species-coefficient table binding + blkdat init — chunk 1
include("variants/eastcascades/diameter_growth.jl")       # EC large-tree DDS (ec/dgf.f): ec_dgcons! + dgf! (32-sp, 3-branch) — chunk 3
include("variants/eastcascades/height_growth.jl")         # EC potential-height curves (ec/htcalc.f) — chunk 4a
include("variants/eastcascades/site_index.jl")            # EC site (ec/forkod/habtyp/ecocls/sichg/sitset) — chunk 2
include("variants/eastcascades/crown.jl")                 # EC crown (ec/crown.f per-species Weibull + dubscr) + CCF (ec/ccfcal.f) — chunk 5
include("variants/eastcascades/regent.jl")                # EC small-tree growth (ec/regent.f + smhtgf.f) + htdbh — chunk 6
include("variants/eastcascades/mortality.jl")             # EC mortality (ec/morts.f Stage/Reineke via shared driver + ec/varmrt.f) — chunk 7
include("variants/eastcascades/volume.jl")                # EC volume (ec/voleqdef R6_EQN forest-8 INGY + Behre) — chunk 8

# --- centralcalifornia (CA) — Inland California westside Wykoff-DDS (50-species, 13-group-compressed, Zeide SDI) ---
include("variants/centralcalifornia/centralcalifornia.jl")  # CA singleton + registration (MAXSP 50, Zeide) — chunk 0
include("variants/centralcalifornia/species.jl")            # CA species-coefficient table binding + blkdat init — chunk 1
include("variants/centralcalifornia/diameter_growth.jl")    # CA large-tree DDS (ca/dgf.f): ca_dgcons! + dgf! (13-group, 3-branch) — chunk 3
include("variants/centralcalifornia/site_index.jl")         # CA site (ca/forkod/habtyp/ecocls/sitset, R6ADJ fan) — chunk 2
include("variants/centralcalifornia/htdbh.jl")              # CA height↔DBH Curtis-Arney missing-height dub (ca/htdbh.f) — chunk 4a
include("variants/centralcalifornia/crown.jl")              # CA crown ratio (ca/crown.f 17-group Weibull + ca/dubscr.f) — chunk 5
include("variants/centralcalifornia/volume.jl")             # CA volume (R6 NVEL: Behre + FW2 Flewelling, ca/formcl.f) — chunk 8
include("variants/centralcalifornia/height_growth.jl")      # CA large-tree HTG (ca/htgf.f + findag.f + htcalc.f) — chunk 4b
include("variants/centralcalifornia/small_tree_growth.jl")  # CA small-tree growth (ca/regent.f + smhtgf.f) — chunk 6
include("variants/centralcalifornia/mortality.jl")          # CA mortality (ca/morts.f + varmrt.f): shared driver + ri_scale 0.5 — chunk 7
include("variants/southcentraloregon/southcentraloregon.jl")# SO (SORNEC-33) singleton + registration (MAXSP 33, DGSD 2.0, seed 55329) — chunk 0
include("variants/southcentraloregon/species.jl")           # SO species table binding + so_bratio (so/bratio.f 3-path) — chunk 1
include("variants/southcentraloregon/site_index.jl")        # SO site (forkod/sichg/htcalc/sitset): SITEAR fan + LZEIDE reset — chunk 2
include("variants/southcentraloregon/diameter_growth.jl")   # SO large-tree DDS (so/dgf.f): so_dgcons! + dgf! (33-sp, 4-branch, RMAI/RELDEN) — chunk 3
include("variants/southcentraloregon/htdbh.jl")             # SO height↔DBH (so/htdbh.f): forest-fanned Curtis DESCHT/FREMNT/WINEMA — chunk 4a
include("variants/southcentraloregon/height_growth.jl")     # SO large-tree HTG (so/htgf.f + so/findag.f): potential-ht + Hoerl/CR + Ritchie-Hann + Johnson-SBB — chunk 4b
include("variants/southcentraloregon/crown.jl")             # SO crown ratio (so/crown.f + so/dubscr.f): rank-Weibull + logistic/linear DUBSCR — chunk 5
include("variants/southcentraloregon/mortality.jl")         # SO mortality VARMRT (so/scomrt.f): PEFF·VARADJ·0.01 — chunk 7
include("variants/southcentraloregon/small_tree_growth.jl") # SO small-tree HEIGHT (so/smhtgf.f) so_smhtgf — chunk 6 step 1 (driver TODO)
include("variants/southcentraloregon/volume.jl")           # SO volume (so/formcl.f + VOLEQDEF): Behre 616BEHW + INGY FW2 — chunk 8
include("variants/westsierra/westsierra.jl")                # WS (Western Sierra Nevada) singleton + registration (MAXSP 43, DGSD 2.0, Zeide, seed 55329) — chunk 0
include("variants/westsierra/species.jl")                   # WS species table binding + ws_bratio (ws/bratio.f IMAP dispatch) — chunk 1
include("variants/westsierra/site_index.jl")                # WS site (forkod/sichg/htcalc/sitset): SITEAR fan (Zeide) — chunk 2
include("variants/westsierra/diameter_growth.jl")           # WS large-tree DDS (ws/dgf.f): ws_dgcons! + dgf! (43-sp, 6-branch) — chunk 3
include("variants/westsierra/height_growth.jl")             # WS large-tree HTG (ws/htgf.f): ws_htcons! + height_growth! (linear CASE DEFAULT) — chunk 4a
include("variants/westsierra/htdbh.jl")                     # WS Curtis-Arney H→D (ws/htdbh.f MODE=1) for regent CA-surrogate DK/DKK — chunk 4b
include("variants/westsierra/crown.jl")                     # WS crown ratio (ws/crown.f + ccfcal.f): ws_tree_ccf + crown_ratio_update! (rank-Weibull, grouped SCALE) — chunk 5
include("variants/westsierra/small_tree_growth.jl")        # WS small-tree HTG+DBH (ws/regent.f + smhtgf.f): ws_smhtgf + small_tree_growth! — chunk 6
include("variants/westsierra/mortality.jl")                   # WS mortality VARMRT efficiency (ws/varmrt.f): _varmrt_efftr! (×0.1 default, GB ×0.01) — chunk 7
include("variants/westsierra/volume.jl")                      # WS volume (ws VEQNNC): 500WO2W R5TAP + 500DVEW r5harv (reuse NC kernels) — chunk 8

# --- olympic (OP) — ORGANON NWO (VERSION=2) follow-on to OC, MAXSP=39; foundation + DG_NWO core ---
include("variants/olympic/olympic.jl")                        # OP (Olympic) singleton + registration (MAXSP 39, ORGANON NWO) — foundation
include("variants/olympic/species.jl")                        # OP species map (orgspc OSPMAP ∘ SPGROUP_RUN SCODE2) → NWO group — validated
include("variants/olympic/organon_diamgro_nwo.jl")            # OP ORGANON NWO diameter-growth core (DG_NWO_RUN) — VALIDATED bit-exact vs FVSop_clean
include("variants/olympic/organon_nwo.jl")                    # OP ORGANON NWO engine: bark + MCW/LCW/HLCW/CW/HCB/MAXHCB + HG_NWO/B_HG/HD_NWO/LIMIT + PM_NWO + SUBMAX + SSTATS/CRNCLO + DG/HG/CR/MORT passes
include("variants/olympic/diameter_growth.jl")                # OP FVS-native large-tree DGF (op/dgf.f Wykoff ln(DDS)) for the non-ORGANON species — VALIDATED per-tree vs FVSop_clean
include("variants/olympic/height_growth.jl")                  # OP FVS-native large-tree HTGF (op/htgf.f+findag+htcalc) for the non-ORGANON species — VALIDATED per-tree vs FVSop_clean
include("variants/olympic/site_index.jl")                     # OP site-index fan (op/forkod+habtyp+ecocls+sichg+htcalc+sitset) — chunk 2; SITEAR validated vs FVSop_clean
include("variants/olympic/volume.jl")                         # OP BLM Behre-taper cubic+board volume (blmvol/blmtap, reuse OC) — chunk 2
include("variants/olympic/organon_hook_op.jl")                # OP COOPERATING driver: op_organon_prepare! (LSTART) + diameter_growth!(::Olympic) (DGDRIV: dgf!+ORGANON fold→WK2→DG) — multi-cycle
include("variants/olympic/crown.jl")                          # OP crown ratio (op/crown.f Weibull + op/dubscr.f) cooperating with ORGANON CR2 — multi-cycle
include("variants/olympic/mortality.jl")                      # OP mortality (op/morts.f: ORGANON MORTEXP-for-all when ORGANON ran) — multi-cycle
include("variants/olympic/small_tree_growth.jl")              # OP small-tree (op/regent.f) guarded no-op (S248112 all large trees) — multi-cycle

# --- io ---------------------------------------------------------------------
include("io/treedata.jl")
include("io/keyword.jl")
include("io/csv_trees.jl")
include("io/yaml_keywords.jl")
include("io/yaml_stand.jl")
include("io/input.jl")

# --- engine -----------------------------------------------------------------
include("engine/species_translation.jl")
include("engine/treeinput.jl")
include("engine/keyword_dispatch.jl")
include("engine/root_disease.jl")        # Western Root Disease (WRD) — Chunk −1 infra + reader + inert seam
include("engine/climate.jl")
include("io/fia_database.jl")            # DATABASE/DSNIN input: FIA "FVS-ready" SQLite → stand
include("io/fia_translate.jl")           # raw FIADB (PLOT/COND/TREE/…) → FVS-ready records
include("engine/init.jl")
include("engine/crown_width.jl")
include("engine/forest_type.jl")
include("engine/r8clark_vol.jl")
include("engine/r9clark_vol.jl")        # NE: NVEL Region-9 Clark profile volume
include("engine/cr_dve_vol.jl")         # CR: NVEL R3 D2H (DVE) volume
include("engine/cr_nvb_vol.jl")         # CR: NVEL NSVB (NVB) volume
include("engine/cr_fw2_vol.jl")         # CR: NVEL Flewelling (FW2) stem-profile volume
include("engine/r4vol.jl")              # TT: NVEL Region-4 Matney (MATW) cubic volume
include("engine/r9vol_gevorkiantz.jl")  # CS: NVEL Region-9 Gevorkiantz '900DVEE' (VOLUME METHC=5)
include("engine/volume_equations.jl")
include("engine/volume.jl")
include("engine/standstats.jl")
include("engine/event_monitor.jl")
include("engine/cuts.jl")
include("engine/quickersort.jl")        # RDPSRT/IQRSRT (Scowen 1965) — faithful sorts for COMPRESS
include("engine/compress.jl")            # COMPRESS (act 250) — PC-score tree-record clustering
include("engine/structure_stage.jl")    # SSTAGE — stand structural-stage class (1-6)
include("../data/centralrockies/establishment/estab_coefs.jl")   # CR establishment coefs (XMIN/HHTMAX/ESSUBH-HHT)
include("../data/centralrockies/dwarf_mistletoe.jl")             # CR dwarf mistletoe coefs + pure kernels (mistoe/misintcr.f)
include("engine/establishment.jl")
include("engine/sprout.jl")
include("engine/fire/biomass.jl")        # FFE F1 — Jenkins tree biomass
include("engine/fire/crown_biomass.jl")  # FFE F2 — crown biomass by size class (FMCROWE)
include("engine/fire/cr_crown_biomass.jl") # FFE F2 — CR western crown biomass (FMCROWW)
include("../data/centralrockies/fire/ffe_live_fuel.jl")  # CR FFE FULIVE/FULIVI live-fuel data
include("../data/inlandempire/fire/ffe_fuel.jl")         # IE FFE FULIVE/FUINIE fuel data (reuses _cr_algslp2)
include("../data/easternmontana/fire/ffe_fuel.jl")       # EM FFE FULIVE/FUINIE + MD1/MD2 fuel-model map
include("../data/centralidaho/fire/ffe_fuel.jl")         # CI FFE FULIVE/FUINIE + MAPPVG/MAPS9B fuel-model map
include("../data/bluemountains/fire/ffe_fuel.jl")        # BM FFE FULIVE/FUINIE + ISPMAP/cwcalc-remap (reuses _cr_algslp2, cr_cwcalc, cr_crownw)
include("../data/teton/fire/ffe_fuel.jl")                # TT FFE FULIVE/FUINIE (reuses cr_select fmcfmd)
include("../data/utah/fire/ffe_fuel.jl")                 # UT FFE FULIVE/FUINIE (reuses cr_select fmcfmd)
include("../data/klamath/fire/ffe_fuel.jl")              # NC FFE FULIVE/FULIVI/FUINIE/FUINII (top-2 cover-type; reuses _cr_algslp2)
include("../data/westsierra/fire/ffe_fuel.jl")           # WS FFE FULIVE/FULIVI/FUINIE/FUINII (43-species top-2 cover-type; reuses _cr_algslp2)
include("../data/centralcalifornia/fire/ffe_fuel.jl")    # CA FFE FULIVE/FULIVI/FUINIE/FUINII (50-species top-2 cover; reuses _cr_algslp2)
include("../data/westcascades/fire/ffe_fuel.jl")         # WC FFE FULIVE/FULIVI/FUINIE/FUINII (39-species SINGLE cover-type) + wc_cwcalc
include("../data/pacificnorthwest/fire/ffe_fuel.jl")     # PN FFE FULIVE/FULIVI/FUINIE/FUINII (39-species SINGLE cover-type) + pn_cwcalc (forest-612 BF)
include("../data/eastcascades/fire/ffe_fuel.jl")         # EC FFE FULIVE/FULIVI/FUINIE/FUINII (32-species SINGLE cover-type) + ec_cwcalc (forest-608 BF) + ec_moist
include("../data/southcentraloregon/fire/so_cwcalc.jl")  # SO crown width (so/cwcalc.f SOMAP, forest-601 DESCHUTES BF)
include("../data/southcentraloregon/fire/ffe_fuel.jl")   # SO FFE FCCS/Ottmar fuel loading (so/fmcba.f COVRINI/FUELINI, FMSSTAGE-keyed)
include("engine/fire/fuel_loading.jl")   # FFE F3 — initial surface fuel loading (FMCBA)
include("engine/fire/fmcba.jl")          # FFE F3 — per-cycle fuel & cover-type update (FMCBA)
include("engine/fire/fuel_decay.jl")     # FFE F3 — per-cycle surface-fuel decay (FMCWD)
include("engine/fire/fuel_additions.jl") # FFE F3 — annual fuel additions / litterfall (FMCADD)
include("engine/fire/fire_effects.jl")   # FFE F6 — fire-caused tree mortality (FMEFF/FMBRKT)
include("engine/fire/rothermel.jl")      # FFE F5 — Rothermel surface fire behavior (FMFINT)
include("engine/fire/fuel_moisture.jl")  # FFE F5b — fuel-moisture scenario + wind reduction (FMMOIS)
include("engine/fire/fuel_model.jl")     # FFE F5b — dynamic fuel-model construction (FMCFMD3)
include("engine/fire/nc_fuel_model.jl")  # FFE F5b — NC California-CWHR fuel-model selection (nc/fmcfmd.f + cwhr.f)
include("engine/fire/ws_fuel_model.jl")  # FFE F5c — WS California-CWHR fuel-model selection (ws/fmcfmd.f; reuses _ca_cwhr)
include("engine/fire/ca_fuel_model.jl")  # FFE F5d — CA California-CWHR fuel-model selection (ca/fmcfmd.f; reuses _ca_cwhr)
include("engine/fire/wc_fuel_model.jl")  # FFE F4b — WC FIRE-VPN cover-metagroup fuel-model selection (wc/fmcfmd.f)
include("engine/fire/ec_fuel_model.jl")  # FFE F4  — EC FMDYN dynamic cover-metagroup fuel-model selection (ec/fmcfmd.f)
include("engine/fire/so_fuel_model.jl")  # FFE F4  — SO FMDYN Oregon 8-plant-group fuel-model selection (so/fmcfmd.f)
include("engine/fire/fmburn.jl")         # FFE F5b — fire event driver (FMBURN/FMEFF) → kill TPA
include("engine/fire/carbon.jl")         # FFE F8 — standing live-tree carbon pools (FMCRBOUT)
include("engine/fire/snag.jl")           # FFE F7 — snag falldown + decay dynamics (FMSFALL)
include("engine/fire/consumption.jl")    # FFE F7/F8 — fire fuel consumption + carbon release (FMCONS)
include("engine/econ.jl")                # C8 — ECON economic-analysis core (eccalc.f)
include("io/summary.jl")
include("io/dbs_output.jl")
include("engine/simulate.jl")
include("engine/svs.jl")                 # SVS (Stand Visualization System) data path — chunk 0

# --- more engine, extensions, cli are added in later chunks -----------------
# include("engine/...")    # C2–C5
# include("extensions/...")# C6–C8
# include("cli.jl")        # C8

export StandState, Southern, Northeast, CentralStates, LakeStates, CentralRockies, Kootenai, InlandEmpire, EasternMontana, Teton, Utah, BlueMountains, CentralIdaho, BritishColumbia, Klamath, OregonCoast, SoutheastAlaska, WestCascades, PacificNorthwest, EastCascades, AbstractVariant, variant_code, variant_from_code
export load_species_coefficients!, init_blockdata!
export resolve_species, translate_species
export FVSRng, rann!, esrann!, bachlo, TreeList, ntrees
export parse_tree_format, parse_tree_record, read_tree_file, TreeRecord, DEFAULT_TREE_FORMAT
export KeywordReader, read_keyword!, KeywordRecord, KeywordStatus, KW_OK, KW_EOF, KW_STOP, KW_PARMS
export read_tree_records, read_trees_csv, write_trees_csv, convert_tre_to_csv, TREE_CSV_HEADER
export read_keyword_records, read_keywords_yaml, write_keywords_yaml, convert_key_to_yaml, read_keyfile_records
export write_keyfile, write_tree_file, convert_yaml_to_key, convert_csv_to_tre, translate_io
export write_sum_csv, SUM_CSV_HEADER, yaml_variant_code, yaml_output_format
export initialize, initialize!, each_stand, run_keyfile, process_keywords!, load_trees!, strip_key_ext
export dgcons!, dgf!, notre!, stand_tpa, stand_ba, stand_qmd, stand_sdi, stand_ccf, stand_top_height

end # module FVSjl
