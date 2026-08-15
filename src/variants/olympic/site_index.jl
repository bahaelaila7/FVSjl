# =============================================================================
# site_index.jl (olympic) — OP site-index fan + forest code + SDImax. Chunk 2.
# (op/forkod.f, op/habtyp.f, op/ecocls.f (buildDir), op/sichg.f, op/htcalc.f, op/sitset.f)
#
# Region-6 chain, MEASURED from the live FVSop_clean run on stand S248112 (forest 708 BLM Salem,
# habitat 40 → plant association CHS133):
#   • forkod: JFOR=[609,612,800,708,709,712]; 708 → IFOR=4 (BLM Salem), IGL=1. Reservation
#     pseudo-codes 8101..8129 cross-walk to the NF/BLM forest (op/forkod.f).
#   • habtyp: KODTYP indexes the 75-entry PCOML plant-association list; numeric habitat N (1..75)
#     → PCOML[N]; anything else → default PA CHS133 (index 40) (op/habtyp.f).
#   • ecocls: PA → site species / site index / default max-SDI (op/ecocls.f 75-entry table). The PA
#     order is IDENTICAL to habtyp's PCOML, so KODTYP indexes it directly.
#   • sichg: SIAGE(I) reference-age fan (op/sichg.f A/B/REFAGE/REFLOC); MH(20) is the metric species.
#   • sitset: ISISP/SITEAR(ISISP) from ecocls (Region-6 default site sp = DF 16, SITEAR 100 if unset);
#     DF(16)↔WH(19) Nigh-1995 conversion; the fan HTCALC(SINDX, ISISP, SIAGE(ISPC)) uses the SITE-
#     species (DF) curve for EVERY species (op_htcalc), then per-species reductions + WO(28) King
#     transform + MH(20) metric. SDIDEF ← min(RSDI, FORMAX=950); non-site species inherit SDIDEF(ISISP).
#
# VALIDATED vs the live FVSop_clean SITECODE dump (S248112, forest 708, habitat 40):
#   SITEAR fan reproduces WF=97.98189, ES=139.15092, LP=98.22826, SP=139.15092, PP=139.15092,
#   DF=98.00 (and the rounded whole-fan SF/AF=139, GF/RF/SS/RC/WH=98, YC/IC/JP/WP=139 …).
# =============================================================================

const OP_FORMAX = 950.0f0                    # op/sitset.f DATA FORMAX/950./
const OP_PMSDIU = 85.0f0                     # op/grinit.f PMSDIU

# op/forkod.f JFOR / KFOR (all 1). Reservation pseudo-codes → NF/BLM forest index.
const OP_JFOR = Int[609, 612, 800, 708, 709, 712]
const OP_FORKOD_RES = Dict(
    8101=>2, 8102=>2, 8103=>2, 8104=>6, 8105=>6,
    8110=>1, 8111=>1, 8113=>1, 8114=>1, 8115=>1, 8116=>1, 8119=>1, 8120=>1,
    8121=>1, 8122=>1, 8123=>1, 8125=>1, 8126=>1, 8127=>1, 8128=>1, 8129=>1)

function op_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 1; useigl = true
    if haskey(OP_FORKOD_RES, kodfor)
        ifor = OP_FORKOD_RES[kodfor]
    else
        idx = findfirst(==(kodfor), OP_JFOR)
        idx === nothing ? (useigl = false; ifor = 1) : (ifor = idx)
    end
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(1))            # KFOR all = 1
    p.user_forest_code = Int32(OP_JFOR[ifor])
    return ifor
end

# op/ecocls.f 75-entry plant-association table (PA order == op/habtyp.f PCOML). Every entry is a single
# site species (NUMBR=1, IFLAG=1); we store (site_index, max_sdi, fvsseq). Index = KODTYP.
const OP_ECOCLS = [
    (54f0,750f0,16),(62f0,955f0,16),(33f0,600f0,16),(50f0,367f0,4),(65f0,535f0,4),
    (91f0,955f0,4),(31f0,560f0,4),(150f0,1050f0,1),(84f0,950f0,16),(83f0,1093f0,1),
    (145f0,995f0,1),(154f0,845f0,1),(84f0,861f0,16),(149f0,1015f0,1),(83f0,1050f0,1),
    (127f0,1090f0,1),(108f0,835f0,1),(84f0,1090f0,16),(101f0,1101f0,1),(136f0,1055f0,1),
    (111f0,1080f0,1),(115f0,955f0,1),(118f0,920f0,1),(107f0,996f0,1),(96f0,470f0,1),
    (104f0,780f0,19),(110f0,960f0,19),(114f0,925f0,19),(94f0,950f0,19),(116f0,1010f0,16),
    (108f0,1040f0,16),(50f0,696f0,19),(105f0,1165f0,19),(52f0,760f0,22),(118f0,985f0,16),
    (114f0,820f0,19),(112f0,1210f0,19),(78f0,1050f0,19),(67f0,880f0,16),(98f0,1606f0,16),
    (78f0,810f0,16),(84f0,895f0,19),(108f0,975f0,16),(71f0,1095f0,16),(119f0,955f0,16),
    (122f0,825f0,16),(128f0,875f0,16),(98f0,1107f0,19),(114f0,840f0,16),(84f0,945f0,19),
    (80f0,865f0,19),(56f0,1145f0,19),(56f0,610f0,16),(80f0,1065f0,16),(66f0,810f0,16),
    (88f0,845f0,16),(110f0,675f0,19),(94f0,660f0,19),(119f0,600f0,16),(134f0,485f0,16),
    (114f0,375f0,19),(118f0,935f0,19),(98f0,1025f0,19),(70f0,610f0,16),(94f0,570f0,19),
    (92f0,915f0,19),(14f0,1021f0,20),(120f0,930f0,6),(115f0,930f0,6),(120f0,930f0,6),
    (125f0,1000f0,6),(117f0,615f0,6),(123f0,545f0,6),(111f0,535f0,6),(121f0,1000f0,6)]
const OP_ECOCLS_DEFAULT = 40                  # op/habtyp.f default ITYPE (CHS133)

"op/habtyp.f — numeric habitat KODTYP → ecocls index (default CHS133=40 when out of 1..75)."
@inline op_habtyp_index(kodtyp::Integer) = (1 <= kodtyp <= 75) ? Int(kodtyp) : OP_ECOCLS_DEFAULT

# op/sichg.f DATA (A, B, REFAGE, REFLOC[0=B,1=T]) — 39 species.
const OP_SICHG_A = Float32[
    21.35000,9.01840,9.01840,33.72545,14.81367, 6.01954,7.938059,11.56252,8.00000,33.72545,
    10.65724,8.000000,21.02281,21.02281,8.00000, 4.94166,11.56252,6.01954,6.15767,45.24969,
    11.56252,3.28241,11.56252,11.56252,11.56252, 11.56252,11.56252,4.94166,11.56252,8.11668,
    11.56252,11.56252,11.56252,11.56252,11.56252, 11.56252,11.56252,11.56252,11.56252]
const OP_SICHG_B = Float32[
    -0.10290,-0.05700,-0.05700,-0.27451,-0.11174, -0.03636,-0.02873,-0.05586,-0.04286,-0.27451,
    -0.10667,-0.04259,-0.15978,-0.15978,-0.04286, -0.02419,-0.05586,-0.03636,-0.03596,-1.28885,
    -0.05586,-0.02987,-0.05586,-0.05586,-0.05586, -0.05586,-0.05586,-0.02419,-0.05586,-0.05661,
    -0.05586,-0.05586,-0.05586,-0.05586,-0.05586, -0.05586,-0.05586,-0.05586,-0.05586]
const OP_SICHG_REFAGE = Int[
    100,50,50,100,50, 50,100,100,100,100, 50,100,100,100,100, 50,100,50,50,100,
    100,20,100,100,100, 100,100,50,100,50, 100,100,100,100,100, 100,100,100,100]
const OP_SICHG_REFLOC = Int[                    # 0='B', 1='T'  (10*B, T, 10*B, T, 17*B)
    0,0,0,0,0, 0,0,0,0,0, 1,0,0,0,0, 0,0,0,0,0, 1,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0]

"op/sichg.f — SIAGE(I) reference-age fan for site species `isisp`, site index `ssite`. MH(20) metric."
function op_sichg(isisp::Int, ssite::Float32)
    isiloc = OP_SICHG_REFLOC[isisp]
    siage = Vector{Float32}(undef, 39)
    @inbounds for i in 1:39
        diff = 0
        (isiloc == 1 && OP_SICHG_REFLOC[i] == 0) && (diff = -1)
        (isiloc == 0 && OP_SICHG_REFLOC[i] == 1) && (diff = 1)
        age2bh = 0f0
        if diff != 0
            age2bh = OP_SICHG_A[i] + OP_SICHG_B[i]*ssite
            (isisp != 20 && i == 20) && (age2bh = OP_SICHG_A[i] + OP_SICHG_B[i]*(ssite/3.281f0))
            (isisp == 20 && i != 20) && (age2bh = OP_SICHG_A[i] + OP_SICHG_B[i]*(ssite*3.281f0))
        end
        siage[i] = Float32(OP_SICHG_REFAGE[i]) + age2bh*Float32(diff)
    end
    return siage
end

# op/sitset.f per-species site reductions (misc hardwoods; ISPC≠ISISP): ISPC → factor.
const OP_SITE_REDUX = Dict(21=>0.75f0, 23=>0.65f0, 24=>0.70f0, 25=>0.70f0, 26=>0.75f0,
                           27=>0.85f0, 29=>0.23f0, 31=>0.70f0, 33=>0.25f0, 34=>0.60f0,
                           35=>0.25f0, 36=>0.50f0, 37=>0.50f0)

function op_sitset!(s::StandState)
    p = s.plot
    maxsp = 39
    si = p.sp_site_index; sdi = p.sp_sdi_def
    nsiset = count(>(0f0), @view si[1:maxsp])

    # ECOCLS for the stand's plant association (single entry).
    idx = op_habtyp_index(Int(p.habitat_code))
    rsi, rsdimx, iseq = OP_ECOCLS[idx]
    rsdi = min(rsdimx, OP_FORMAX)

    isisp = (1 <= Int(p.site_species) <= maxsp) ? Int(p.site_species) : 0
    if iseq > 0
        (isisp <= 0) && (isisp = iseq)                         # IFLAG=1 site species
        (si[iseq] <= 0f0 && nsiset == 0) && (si[iseq] = rsi)
        (sdi[iseq] <= 0f0) && (sdi[iseq] = rsdi)
    end
    isisp <= 0 && (isisp = 16)                                 # Region-6 global default site sp = DF
    si[isisp] <= 0f0 && (si[isisp] = 100f0)

    sindx0 = si[isisp]
    siage = op_sichg(isisp, sindx0)

    # DF(16)↔WH(19) ORGANON site conversion, Nigh (1995) — op/sitset.f:148-154.
    if si[16] > 0f0 || si[19] > 0f0
        if si[16] <= 0f0
            si[16] = 0.480f0 + 1.110f0*si[19]
        elseif si[19] <= 0f0
            si[19] = -0.432f0 + 0.899f0*si[16]
        end
    end

    sindx = si[isisp]
    @inbounds for ispc in 1:maxsp
        si[ispc] > 0f0 && continue
        v = op_htcalc(sindx, isisp, siage[ispc])               # site-SPECIES curve (op/htcalc.f)
        if ispc == 20 && isisp != 20
            v = v/3.281f0; v > 28f0 && (v = 28f0)              # MH metric, cap 28
        end
        haskey(OP_SITE_REDUX, ispc) && isisp != ispc && (v *= OP_SITE_REDUX[ispc])
        (ispc == 28 && isisp != 28) &&
            (v = 114.2f0*(1f0 - exp(-0.0266f0*v))^2.26f0)      # WO King-from-DF (op/sitset.f:196)
        si[ispc] = v
    end

    # SDIDEF: non-set species inherit the site-species default.
    @inbounds for i in 1:maxsp
        sdi[i] <= 0f0 && (sdi[i] = sdi[isisp])
    end
    s.control.ba_max <= 0f0 &&
        (s.control.ba_max = sdi[isisp]*(OP_PMSDIU/100f0)*0.54542f0)   # BAMAX (op/sitset.f:207)
    p.site_species = Int32(isisp)
    return s
end

function op_site_index_setup!(s::StandState)
    op_forkod!(s.plot)
    op_sitset!(s)
    return s
end

site_setup!(s::StandState, ::Olympic) = op_site_index_setup!(s)
