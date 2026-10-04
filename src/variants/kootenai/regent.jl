# =============================================================================
# regent.jl (kootenai) — KT small-tree growth (kt/regent.f). CHUNK 6.
# kt_regcons! = REGCON entry (RHCONS site constant, analog of kt_dgcons!).
# small_tree_growth!(::Kootenai) = the multi-subcycle height+DG model (PENDING — this commit
# lands the validated RHCON setup + all coefficients; the growth loop is the next step).
# NOTE: regent MAPLOC/MAPHAB DIFFER from dgf.f (measured) — own tables below.
# =============================================================================

const KT_RG_HTSLOP = Float32[-0.739841f0, -1.150486f0, -0.74956f0, -0.730825f0, 0.047743f0, -0.034439f0, -0.171346f0, -0.140742f0, 0.459091f0, -1.923583f0, 0f0]
const KT_RG_HTSLSQ = Float32[-0.347447f0, 0f0, 0.726317f0, 0.755344f0, 0f0, 0f0, -1.70513f0, 0f0, -1.0158f0, 0f0, 0f0]
const KT_RG_HTEL   = Float32[-0.017447f0, -0.001659f0, -0.009145f0, -0.009281f0, -0.008651f0, 0.001325f0, -0.007956f0, -0.006651f0, -0.005558f0, 0.005946f0, 0f0]
const KT_RG_HTCASP = Float32[-0.319768f0, 0.506689f0, 0.02547f0, -0.019985f0, -0.200859f0, -0.148969f0, 0.238473f0, -0.174545f0, 0.047339f0, -0.562835f0, 0f0]
const KT_RG_HTSASP = Float32[-0.22656f0, -0.019319f0, -0.061316f0, -0.002931f0, 0.035112f0, -0.086932f0, 0.171464f0, -0.083216f0, -0.016025f0, 0.56742f0, -0f0]
const KT_RG_RHSC   = Float32[0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0.8953f0]
const KT_RG_RHGL   = Float32[-0.2785f0, -0.048f0, 0f0]   # sp11 IGL 1..3
const KT_RG_RSAB   = Float32[-0.010987f0, 0.22157f0, -0.12432f0]   # sp11 RSAB0/1/2
const KT_RG_HTFOR  = Float32[
    2.269835f0 2.10643f0 2.64339f0 0f0 0f0;
    2.887315f0 2.138032f0 2.577398f0 3.316322f0 0f0;
    2.77609f0 2.731853f0 2.898513f0 0f0 0f0;
    1.508043f0 1.336607f0 1.653463f0 0f0 0f0;
    1.856635f0 1.882195f0 0f0 0f0 0f0;
    1.102706f0 1.227212f0 1.508519f0 0f0 0f0;
    3.900741f0 3.814237f0 4.274753f0 0f0 0f0;
    1.605062f0 1.39565f0 2.304788f0 2.051997f0 1.297702f0;
    1.401828f0 1.625757f0 1.99335f0 1.66862f0 0f0;
    1.370572f0 1.604815f0 1.083899f0 0f0 0f0;
    0f0 0f0 0f0 0f0 0f0]  # [sp,irhfor 1..5]
const KT_RG_RHHAB  = Float32[
    0f0 -0.288729f0 0.222331f0 0f0 0f0 0f0;
    0f0 0.095681f0 0f0 0f0 0f0 0f0;
    0f0 0.259585f0 0.143218f0 0.107012f0 -0.136078f0 0f0;
    0f0 0.111443f0 0.18247f0 0f0 0f0 0f0;
    0f0 0.128458f0 -0.115672f0 0.271298f0 0f0 0f0;
    0f0 -0.107416f0 0f0 0f0 0f0 0f0;
    0f0 -0.637727f0 -0.808173f0 -0.39611f0 0.666094f0 0f0;
    0f0 0.245211f0 0.726173f0 0.341746f0 0.553888f0 0.427478f0;
    0f0 0.252628f0 0.326399f0 0.4012f0 0.568246f0 0.785996f0;
    0f0 -0.121653f0 -0.374652f0 0f0 0f0 0f0;
    -0.2146f0 -0.0941f0 -0.3738f0 0f0 0f0 0f0]  # [sp,irhhab 1..6]
const KT_RG_MAPLOC = Int32[
    1 1 1 1 1 1 1 1 1 1 1;
    2 1 2 2 2 1 1 2 2 2 1;
    1 2 2 1 2 2 2 2 2 1 1;
    2 1 2 2 1 1 3 2 2 3 1;
    2 2 2 1 1 1 1 5 4 2 1;
    2 2 1 1 2 2 1 4 3 1 1;
    1 3 3 1 1 1 1 3 4 1 1;
    3 1 3 3 1 3 2 1 4 3 1;
    3 3 3 2 2 1 2 5 4 2 1;
    3 4 3 3 1 2 3 1 4 2 1]  # [kotfor 1..10, sp]
const KT_RG_MAPHAB = Int32[
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 2 1 1 1 1 2 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 2 1 1 1 2 1 1 2 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 2 3;
    1 2 2 1 1 1 2 1 1 2 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 2 1 1 1 1 1 1 2 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 2 3;
    2 2 2 2 1 1 3 1 1 3 3;
    1 2 3 1 1 1 2 1 1 3 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 2 2 1 1 2 1 1 3 3;
    2 2 2 2 1 1 2 2 2 2 3;
    2 2 3 2 1 1 4 2 2 1 3;
    2 2 5 2 1 1 3 2 2 3 3;
    2 2 4 2 1 1 3 2 2 3 3;
    1 2 3 2 1 1 2 2 4 2 3;
    2 2 4 2 1 1 2 2 2 2 3;
    2 2 4 2 1 1 2 2 2 2 3;
    1 2 3 1 1 1 2 2 2 2 3;
    2 2 4 2 1 1 2 2 2 2 3;
    1 2 3 2 1 1 2 2 2 1 3;
    2 2 3 2 1 1 3 2 2 1 3;
    2 2 4 2 1 1 2 2 2 2 3;
    2 2 4 2 1 1 2 2 2 2 3;
    2 2 3 2 1 1 2 2 2 2 3;
    1 2 3 1 1 1 2 1 1 2 3;
    2 2 4 2 1 1 3 2 2 2 3;
    2 2 5 1 1 1 3 2 2 2 3;
    1 1 1 1 1 1 1 1 1 1 3;
    2 2 4 2 1 1 4 2 2 1 3;
    1 2 3 2 1 1 2 2 2 3 3;
    2 2 4 2 1 1 2 2 2 2 3;
    2 2 4 2 1 1 2 2 2 1 3;
    2 2 2 2 1 1 4 1 2 3 3;
    1 2 3 1 1 1 3 1 2 2 3;
    1 2 4 1 1 1 1 1 1 1 3;
    1 2 3 1 1 1 3 1 1 1 3;
    2 2 4 1 1 1 2 1 1 2 3;
    1 2 2 1 1 1 2 1 2 2 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 1 1 1 1 2 1 2 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 4 1 1 1 5 5 2 1 3;
    1 2 3 1 1 1 4 2 2 1 3;
    1 2 3 1 1 1 3 2 2 1 3;
    1 2 3 1 1 1 3 2 2 1 3;
    1 2 1 1 1 1 2 2 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 3 1 1 1 2 2 2 1 3;
    1 1 1 1 1 1 2 1 1 1 3;
    1 2 3 1 1 1 2 2 2 1 3;
    1 2 3 1 1 1 2 5 2 1 3;
    1 1 1 1 1 1 4 1 1 1 3;
    1 1 3 3 1 1 1 1 1 3 3;
    1 2 3 1 4 1 4 4 4 3 3;
    1 2 2 2 2 1 4 3 3 2 3;
    1 1 3 1 1 1 4 1 1 3 3;
    1 2 2 3 1 1 1 1 1 1 3;
    1 2 3 3 1 1 4 1 1 3 3;
    1 2 3 3 4 1 2 3 4 3 3;
    1 2 2 3 1 1 4 1 1 3 3;
    1 2 3 3 1 1 1 1 1 1 3;
    1 2 3 3 4 1 4 3 4 2 1;
    1 2 4 3 4 1 4 4 4 2 1;
    1 2 4 3 1 1 4 3 6 2 1;
    1 2 3 3 1 1 2 4 4 2 1;
    1 2 3 3 4 1 4 4 1 2 1;
    1 2 4 1 1 1 1 1 1 1 1;
    1 2 2 3 2 1 4 4 4 3 2;
    1 2 4 1 2 1 4 4 4 3 2;
    1 2 2 3 2 1 4 4 4 3 2;
    1 2 4 3 2 1 2 4 4 1 2;
    1 2 4 3 2 1 4 4 5 1 2;
    1 2 2 1 4 1 4 3 5 1 2;
    1 2 3 1 2 1 2 3 5 1 2;
    1 2 3 1 1 1 1 3 1 1 2;
    1 2 4 1 4 1 4 3 5 3 4;
    1 2 2 2 2 1 2 3 4 1 4;
    3 2 2 2 2 1 4 4 3 2 4;
    3 2 3 2 2 1 2 2 4 2 4;
    1 2 2 2 2 1 4 3 3 2 4;
    1 2 2 2 2 1 4 2 4 2 4;
    1 2 2 2 2 1 2 4 3 1 4;
    1 2 2 2 2 1 2 3 4 1 4;
    1 2 1 1 1 1 1 2 3 1 4;
    1 2 2 2 2 1 2 4 3 2 4;
    1 2 2 2 2 1 2 4 3 1 4;
    1 1 1 1 2 1 1 1 1 1 4;
    1 2 3 1 1 1 1 1 1 1 4;
    1 2 4 3 1 1 4 3 2 1 4;
    1 2 3 3 1 1 4 3 2 1 4;
    1 2 4 1 1 1 2 4 2 1 4;
    1 2 3 1 3 1 5 3 6 1 4;
    3 2 3 3 1 1 5 4 5 1 1;
    1 2 4 3 1 1 3 4 4 1 1;
    1 2 3 3 3 1 2 4 5 1 1;
    3 2 4 3 1 1 2 6 3 1 1;
    3 2 3 3 3 1 3 6 3 1 1;
    3 2 3 3 1 1 4 5 5 1 1;
    1 2 3 1 3 1 5 4 5 1 1;
    1 1 1 1 1 1 1 1 1 1 1;
    1 2 3 1 3 1 1 3 6 1 1;
    3 2 1 1 3 1 5 4 3 1 1;
    1 2 3 3 1 1 3 4 5 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 1 1 3 1 2 3 6 1 3;
    1 2 3 1 1 1 3 1 2 1 3;
    1 1 1 1 1 1 2 1 1 1 3;
    1 2 1 1 1 1 3 3 1 1 3;
    1 2 1 1 1 1 2 4 3 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    3 2 1 3 1 1 4 6 2 1 3;
    1 2 1 3 1 1 3 4 3 1 3;
    1 2 3 3 3 1 3 6 3 1 3;
    3 2 1 3 3 1 3 4 2 1 3;
    3 2 3 3 3 1 2 6 3 1 2;
    3 2 1 1 1 1 1 4 3 1 2;
    1 1 1 1 1 1 1 4 5 1 2;
    1 2 3 1 1 1 1 4 3 1 2;
    1 2 1 1 1 2 4 3 3 1 2;
    1 2 3 3 1 1 4 6 3 1 3;
    3 2 3 3 1 2 3 5 3 1 3;
    1 2 3 3 1 2 2 5 3 1 3;
    3 2 4 3 1 2 2 6 3 1 3;
    1 1 4 1 1 1 3 5 3 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    3 2 3 3 1 2 4 4 6 1 3;
    1 1 1 1 1 1 5 3 1 1 3;
    3 2 3 3 1 1 2 4 2 1 3;
    1 2 3 3 1 2 3 4 2 1 3;
    3 2 4 3 1 2 3 6 3 1 3;
    3 2 1 1 1 2 1 4 2 1 3;
    1 2 1 1 1 1 1 4 2 1 3;
    1 2 1 1 1 1 2 1 1 1 3;
    3 2 3 3 1 2 2 5 3 1 3;
    3 2 1 1 1 2 4 5 3 1 3;
    1 2 1 1 1 1 2 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    2 2 1 3 1 2 2 4 2 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 1 1 1 1 1 6 1 1 3;
    1 1 1 1 1 1 2 1 1 1 3;
    1 1 1 1 1 1 3 1 6 1 3;
    1 2 1 1 1 1 1 5 5 1 3;
    1 2 3 1 1 1 1 6 6 1 3;
    1 2 1 1 1 1 1 6 3 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 1 1 1 1 1 1 1 1 1 3;
    1 2 1 1 1 1 1 1 3 1 3;
    1 1 1 1 1 1 5 1 1 1 3]  # [kktype 1..175, sp]
const KT_RG_DIAM = Float32[0.4f0, 0.3f0, 0.3f0, 0.3f0, 0.2f0, 0.2f0, 0.4f0, 0.3f0, 0.3f0, 0.5f0, 0.2f0]
const KT_RG_DCON = Float32[-0.145707f0, -0.111929f0, -0.12769f0, -0.125198f0, -0.178162f0, -0.233566f0, -0.010818f0, -0.186874f0, -0.15887f0, -0.185011f0, 0f0]
const KT_RG_HCON = Float32[0.125404f0, 0.098451f0, 0.115066f0, 0.12879f0, 0.129191f0, 0.143211f0, 0.090927f0, 0.131209f0, 0.128174f0, 0.141907f0, 0f0]
const KT_RG_RHBAL = Float32[-0.002569f0, 0.002745f0, 0.002367f0, 0.000683f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, -0.25349f0]
const KT_RG_RHLH = Float32[0.233302f0, 0.187195f0, 0.185496f0, 0.192256f0, 0.126835f0, 0.083622f0, 0.175631f0, 0.160429f0, 0.169834f0, 0.387617f0, 0.2354f0]
const KT_RG_RHCCF = Float32[0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, -0.00391f0]
const KT_RG_HTH2 = Float32[-0.007483f0, -0.00557f0, -0.005958f0, -0.006796f0, -0.004391f0, -0.00245f0, -0.004546f0, -0.005592f0, -0.005955f0, -0.015505f0, 0f0]
const KT_RG_HT2MOD = Float32[0.5f0, 0.4f0, 0.5f0, 0.5f0, 0.6f0, 0.6f0, 0.7f0, 0.5f0, 0.5f0, 0.5f0, 1f0]
const KT_RG_RHBA = Float32[-0.225869f0, -0.612994f0, -0.580046f0, -0.303384f0, -0.310048f0, -0.188838f0, -0.523915f0, -0.270415f0, -0.337342f0, -0.25931f0, 0f0]
const KT_RG_HTCR = Float32[-0.037519f0, 1.793489f0, 0.41613f0, 1.384217f0, 0.802352f0, 0.70806f0, 0.771466f0, -0.520312f0, -0.274359f0, 2.976688f0, 0f0]
const KT_RG_HTCR2 = Float32[1.895808f0, 1.713339f0, 1.070272f0, 0f0, 0.343107f0, 0.231768f0, 2.223783f0, 1.859532f0, 1.212637f0, -0.803694f0, 0f0]
const KT_RG_HTPCC1 = Float32[0.002726f0, 0.001935f0, 0.002364f0, 0.001374f0, 0.000325f0, 0.000516f0, 0.002926f0, 0.001409f0, 0.001964f0, 0.003596f0, 0f0]
const KT_RG_XMAX = Float32[10f0, 10f0, 10f0, 10f0, 10f0, 10f0, 5f0, 10f0, 10f0, 10f0, 10f0]
const KT_RG_XMIN = Float32[2f0, 2f0, 2f0, 2f0, 2f0, 2f0, 1f0, 2f0, 2f0, 2f0, 2f0]
const KT_RG_HSIGMA = 0.59f0
const KT_RG_REGYR = 5.0f0
"""
    kt_regcons!(s) -> Vector{Float32}

KT REGCON entry (kt/regent.f:1017): the per-species small-tree height-model site constant RHCONS.
sp≠11: RHCON = HTFOR(irhfor,sp) + (HTSLOP + HTSLSQ·SLOPE + HTCASP·cosASP + HTSASP·sinASP)·SLOPE
       + HTEL·ELEV + RHSC(sp) + RHHAB(irhhab,sp);  irhfor = MAPLOC(KOTFOR,sp), irhhab = MAPHAB(KKTYPE,sp).
sp=11 (mtn hemlock): RHGL(IGL) + (RSAB0+RSAB1·cosASP+RSAB2·sinASP)·SLOPE + RHSC(11) + RHHAB(irhhab,11).
(regent's MAPLOC/MAPHAB are its OWN tables, distinct from the DG ones — verified by measurement.)
"""
function kt_regcons!(s::StandState)
    p = s.plot
    kotfor = Int(p.forest_idx); kktype = Int(p.habitat_code)
    slope = p.slope; elev = p.elevation; asp = p.aspect
    ca = fcos(asp); sa = fsin(asp)
    igl = Int(p.geo_location)                       # IGL = KFOR(IFOR) from forkod (sp11 only); 0 → default 3
    (igl < 1 || igl > 3) && (igl = 3)
    rhcon = zeros(Float32, 11)
    @inbounds for sp in 1:11
        irhfor = (1 <= kotfor <= 10) ? Int(KT_RG_MAPLOC[kotfor, sp]) : 1
        irhhab = (1 <= kktype <= 175) ? Int(KT_RG_MAPHAB[kktype, sp]) : 1
        (irhfor < 1 || irhfor > 5) && (irhfor = 1)
        (irhhab < 1 || irhhab > 6) && (irhhab = 1)
        if sp == 11
            regch = KT_RG_RHGL[igl] + (KT_RG_RSAB[1] + KT_RG_RSAB[2]*ca + KT_RG_RSAB[3]*sa) * slope
            rhcon[sp] = regch + KT_RG_RHSC[sp] + KT_RG_RHHAB[sp, irhhab]
        else
            rhcon[sp] = KT_RG_HTFOR[sp, irhfor] +
                (KT_RG_HTSLOP[sp] + KT_RG_HTSLSQ[sp]*slope + KT_RG_HTCASP[sp]*ca + KT_RG_HTSASP[sp]*sa) * slope +
                KT_RG_HTEL[sp]*elev + KT_RG_RHSC[sp] + KT_RG_RHHAB[sp, irhhab]
        end
    end
    return rhcon
end

"""
    small_tree_growth!(s, stash, ::Kootenai; fint)

KT REGENT (kt/regent.f) small-tree height+diameter growth. Multi-subcycle (NPER subcycles ≤5yr each). For each
small tree (D<XMAX): advance height through subcycles (density RDNEXT/BANEXT grow from the large trees), add the
ZZRAN stochastic increment, then BLEND with the large-tree htgf HTG via XWT=(D−XMIN)/(XMAX−XMIN). Diameter for
D<3 is dubbed from height (HCON·H+DCON+DADJ). Writes t.ht_growth (central) + the tripling stash.
CON = exp(HCOR) small-tree height calibration reuses c.htg_cor_small (shared REGENT calib; 0 until a KT branch
computes it). NOTE: density subcycle feedback from OTHER small trees (regent.f:444) is omitted in this first
implementation (only the large-tree density contribution is included) — refine if the differential needs it.
"""
function small_tree_growth!(s::StandState, stash, ::Kootenai; fint::Float32 = 10.0f0,
                            lestb::Bool = false, itrnin::Int = 1,
                            atba::Float32 = -1f0, atccf::Float32 = -1f0, atavh::Float32 = -1f0,
                            ba_now::Float32 = -1f0, relden_now::Float32 = -1f0, pccf_now::Vector{Float32} = Float32[])
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    t.n == 0 && return s
    n = t.n
    rhcon = kt_regcons!(s)
    # ESTAB mode (REGENT(.TRUE.,ITRNIN) from esgent.f): BA/RELDEN/PCCF are the gradd.f:192 DENSE's (post-growth,
    # PRE-ESNUTR — snapshotted by the caller before the new cohort exists); the new records (≥ITRNIN) have PCT=0.
    ba = (lestb && ba_now >= 0f0) ? ba_now : p.basal_area
    relden = (lestb && relden_now >= 0f0) ? relden_now : p.relative_density
    avh = p.avg_height
    _pccf(pt) = (lestb && 1 <= pt <= length(pccf_now)) ? pccf_now[pt] :
                (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
    _pct(i) = (lestb && i >= itrnin) ? 0f0 : t.crown_ratio[i]
    managed = p.managed == Int32(1)
    dgsd = s.control.dg_sd
    regyr = KT_RG_REGYR                                    # 5.0
    # subcycle count + lengths (regent.f:186-203)
    ntyr = Int(round(fint)); iyr = Int(regyr)
    lskiph = false
    if lestb                                               # regent.f:187-188 ESTAB: the cycle after year 5
        ntyr -= 5; lskiph = ntyr <= 0
    end
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); (!lestb && nper < 1) && (nper = 1)
    lestb && nper < 0 && (nper = 0)
    kper = zeros(Int, max(nper, 1)); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    nper > 0 && (kper[nper] = itot == ntyr ? kper[nper] : itot)          # KPER(NPER)=ITOT (regent.f:203)
    # recompute ITOT-based last (mirror FVS exactly): after the loop ITOT holds the remainder
    # (the loop already assigned kper[nper]=itot when nn hit 1)
    # per-subcycle stand density from the LARGE trees (regent.f:236-256)
    banext = fill(ba, max(nper, 1)); rdnext = fill(relden, max(nper, 1))
    # regent.f:206-211 after-thin TEMBA/TEMCCF/TEMAHT (ATBA/ATCCF ≤0 ⇒ the current BA/RELDEN)
    temba = (atba > 0f0) ? atba : ba; temccf = (atccf > 0f0) ? atccf : relden; temaht = atavh
    if lestb
        # regent.f:269-281 — ESTAB: interpolate each subcycle's start density from the after-thin and the current
        # (post-growth) values; BAYR/CCFYR divide by ITOT (the remainder after the KPER split), 0 when LSKIPH.
        bayr = 0f0; ccfyr = 0f0
        if !lskiph
            bayr = (ba - temba) / Float32(itot); ccfyr = (relden - temccf) / Float32(itot)
        end
        nyr = 5
        @inbounds for j in 1:nper
            rdnext[j] = temccf + Float32(nyr) * ccfyr
            banext[j] = temba + Float32(nyr) * bayr
            nyr += kper[j]
        end
        # regent.f:287-300 DO 13 — dub the crown of each new record, STORAGE order, before any growth draw
        @inbounds for i in itrnin:n
            pccf = _pccf(Int(t.plot_id[i]))
            cr0 = 0.89722f0 - 0.0000461f0 * pccf
            ran = 0f0
            while true
                ran = bachlo(s.rng, 0f0, 1f0); (-1f0 <= ran <= 1f0) && break
            end
            cr0 = cr0 + 0.07985f0 * ran
            cr0 > 0.90f0 && (cr0 = 0.90f0); cr0 < 0.20f0 && (cr0 = 0.20f0)
            icr0 = unsafe_trunc(Int32, cr0 * 100f0 + 0.5f0)
            t.crown_pct[i] = icr0; t.crown_ratio[i] = Float32(icr0)
        end
    end
    if !lestb && nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            # CCFCAL returns CCFT·P (kt/ccfcal.f MODE=1), so RDNEXT(J)+K·CI/P·PN nets one P; kt_tree_ccf is per tree.
            c1 = kt_tree_ccf(sp, d1) * pr; c2 = kt_tree_ccf(sp, d2) * pr
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]
                pn = pr * fpowi(0.985f0, k)
                rdnext[j] += k * ci / pr * pn
                banext[j] += k * bi * pn
            end
        end
    end
    # DELMAX (regent.f:313), R=RELDEN, AH=AVH
    ah = lestb ? temaht : avh                               # regent.f:306-312 ESTAB: R=TEMCCF, AH=TEMAHT
    delmax = (ah / 36.0f0) * (0.01232f0 * (lestb ? temccf : relden) - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    # per-tree height accumulator (WK3), starts at HT
    wk3 = Float32[t.height[i] for i in 1:n]
    # ---- subcycle height loop (regent.f:321-460) ----
    ind1 = species_major_order(s)                     # DO 16 ISPC / DO 15 I3=I1,I2 / I=IND1(I3)
    ky = 0
    @inbounds for j in 1:nper
        ky += kper[j]                                  # kt/regent.f:319 KY=KY+KPER(J)
        baj = banext[j]; rdj = rdnext[j]
        alba = baj > 0.0f0 ? flog(baj) : 0.0f0
        kpj = Float32(kper[j])
        for i in ind1
            sp = Int(t.species[i]); d = t.dbh[i]
            d >= KT_RG_XMAX[sp] && continue
            (lestb && i < itrnin) && continue                 # regent.f:376
            t.tpa[i] <= 0.0f0 && continue
            # kt/regent.f:329-331 CON = RCOR2 (READCORR, LRCOR2) · EXP(HCOR)
            con = ((s.control.lrcor2 && s.control.sp_rcor2[sp] > 0f0) ? s.control.sp_rcor2[sp] : 1f0) *
                  fexp(c.htg_cor_small[sp])
            h1 = wk3[i]
            pct = _pct(i)
            bal = baj * (100.0f0 - pct) * 0.01f0
            balmh = baj * (100.0f0 - pct) * 0.0001f0
            cr = Float32(t.crown_pct[i]) / 100.0f0
            pt = Int(t.plot_id[i]); pccf1 = _pccf(pt)
            htpc1 = managed ? KT_RG_HTPCC1[sp] : 0.0f0
            xrhgro = active_multiplier(s.control, :regh, sp, current_cycle_year(s))
            local htgrl::Float32
            if sp != 11
                htgrl = rhcon[sp] + KT_RG_RHLH[sp]*h1 + (KT_RG_HTH2[sp]*KT_RG_HT2MOD[sp])*h1*h1 +
                        KT_RG_RHBAL[sp]*bal + KT_RG_RHBA[sp]*alba + htpc1*pccf1 +
                        (KT_RG_HTCR[sp] + KT_RG_HTCR2[sp]*cr)*cr
                htgrl < 0.0f0 && (htgrl = 0.0f0)
                wk3[i] = h1 + htgrl * (kpj / regyr) * xrhgro * con
            else
                htgrl = fexp(rhcon[sp] + KT_RG_RHLH[sp]*flog(h1) + KT_RG_RHCCF[sp]*rdj + KT_RG_RHBAL[sp]*balmh + c.htg_cor_small[sp])
                htgrl < 0.0f0 && (htgrl = 0.0f0)
                wk3[i] = h1 + htgrl * (kpj / regyr) * xrhgro
            end
            # kt/regent.f:411-447 — UPDATE DENSITY FOR NEXT SUBCYCLE: a small tree (D<3) whose subcycle height passed
            # 4.5 ft adds its dubbed-diameter BA/CCF increase (1.5%/yr mortality) to BANEXT/RDNEXT(J+1). jl skipped it,
            # so subcycle 2 saw a thinner stand (ktt01 BAJ 97.254 / RDJ 103.87 vs live 97.548 / 116.27) ⇒ HTGRL high.
            h2 = wk3[i]
            if j < nper && d < 3.0f0 && h2 > 4.5f0
                pr = t.tpa[i]
                relh = (h1 - 4.5f0) / (ah - 4.5f0); relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
                dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
                d1 = KT_RG_DIAM[sp] + dadj
                local d2::Float32
                if sp == 11
                    h1 > 4.5f0 && (d1 = 0.0729f0*fpow(h1 - 4.5f0, 1.1988f0) + dadj)
                    d2 = 0.0729f0*fpow(h2 - 4.5f0, 1.1988f0) + dadj
                else
                    h1 > 4.5f0 && (d1 = KT_RG_HCON[sp]*h1 + KT_RG_DCON[sp] + dadj)
                    d2 = KT_RG_HCON[sp]*h2 + KT_RG_DCON[sp] + dadj
                end
                xrdgro = active_multiplier(s.control, :regd, sp, current_cycle_year(s))
                dgj = (d2 - d1) * xrdgro; dgj < 0.0f0 && (dgj = 0.0f0)
                d2 = d + dgj
                c1 = kt_tree_ccf(sp, d) * pr; c2 = kt_tree_ccf(sp, d2) * pr     # CCFCAL(…,P,…) = CCFT·P
                b1 = 0.005454154f0 * d * d
                f = fpowi(0.985f0, ky)
                rdnext[j+1] += Float32(ky) * (c2 - c1) / 10.0f0 * f
                banext[j+1] += (0.005454154f0 * d2 * d2 - b1) * pr * f
            end
        end
    end
    # ---- final: HTGR1 + ZZRAN + XWT blend + DG dub (kt/regent.f DO 30, species-sorted). A FRESH ZZRAN per
    #      TRIPLED record (L-loop L=0 central→tree, L=1/2→the two copies). WITHOUT the per-record stash writes
    #      the tripled small-tree records inherit the EXPLOSIVE large-tree DG/HTG from diameter_growth! (the DDS
    #      gemdg balloons tiny DBH) ⇒ regen BA over-predicted 30-60%. The DG dub is the FAITHFUL non-ESTAB path
    #      (kt/regent.f:591-604): DG(K)=(DK−D1)·XRDGRO·BARK on the DDS scale, DBH grows via GRADD — the old code
    #      did a raw (DK−D1)·XRDGRO with NO bark/DDS/SIZCAP scaling. DBH-direct is the HK<4.5 tiny edge only. ----
    if lestb
        # regent.f:473-651 under LESTB, the new records only (I<ITRNIN skipped), species-major: LSKIPH ⇒ HTG=0 (no
        # draw); else HTGR = (WK3−H) + ZZRAN·HSIGMA (ZZRAN∈[−1.5,1]), ≥0.15, XWT=0, SIZCAP. D<3: HK<4.5 ⇒ DBH =
        # 0.1+DIAM·0.01+0.001·HK, DG=0; else the height dub DK (≥DIAM, +0.001·HK) is BOTH DBH and DG (:582-584).
        # DGBND; no tripling, no DUBSCR.
        dgsd_e = s.control.dg_sd
        @inbounds for i in species_major_order(s)
            i < itrnin && continue
            sp = Int(t.species[i]); d = t.dbh[i]
            d >= KT_RG_XMAX[sp] && continue
            h = t.height[i]
            local htg::Float32
            if lskiph
                htg = 0f0
            else
                htgr1 = wk3[i] - h; htgr1 < 0f0 && (htgr1 = 0f0)
                zzran = 0f0
                if dgsd_e >= 1f0
                    while true
                        zzran = bachlo(s.rng, 0f0, 1f0)
                        (zzran <= 1f0 && zzran >= -1.5f0) && break
                    end
                end
                htgr = htgr1 + zzran * KT_RG_HSIGMA; htgr < 0.15f0 && (htgr = 0.15f0)
                htg = htgr
                cap = s.control.sp_size_cap[sp, 4]
                if h + htg > cap
                    htg = cap - h; htg < 0.1f0 && (htg = 0.1f0)
                end
            end
            t.ht_growth[i] = htg
            if d < 3f0
                relh = (h - 4.5f0) / (ah - 4.5f0); relh > 1f0 && (relh = 1f0); relh < 0f0 && (relh = 0f0)
                dadj = delmax*relh*relh - 2f0*delmax*relh + 0.65f0
                hk = h + htg
                local dbhk::Float32, dgk::Float32
                if hk < 4.5f0
                    dbhk = 0.1f0 + KT_RG_DIAM[sp]*0.01f0 + hk*0.001f0; dgk = 0f0
                else
                    dk = sp == 11 ? 0.0729f0*fpow(hk - 4.5f0, 1.1988f0) + dadj : KT_RG_HCON[sp]*hk + KT_RG_DCON[sp] + dadj
                    dk < KT_RG_DIAM[sp] && (dk = KT_RG_DIAM[sp])
                    dk = dk + hk*0.001f0
                    dbhk = dk; dgk = dk
                    (dbhk + dgk) < KT_RG_DIAM[sp] && (dgk = KT_RG_DIAM[sp] - dbhk)
                end
                dgk = dg_bound(nothing, nothing, sp, dbhk, dgk, s.control.sp_size_cap)
                t.dbh[i] = dbhk; t.diam_growth[i] = dgk
            end
        end
        return s
    end
    scale = fint > 0.0f0 ? 10.0f0 / fint : 1.0f0           # SCALE=YR/FINT (kt/regent.f:220), YR=10
    _sp_order = species_major_order(s)   # IND1: SPESRT lineage order within a species (post-TRIPLE copy1, original, copy2)
    @inbounds for i in _sp_order
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= KT_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        h = t.height[i]
        xmn = KT_RG_XMIN[sp]; xmx = KT_RG_XMAX[sp]
        htgr1 = wk3[i] - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        cap = s.control.sp_size_cap[sp, 4]
        large_htg = t.ht_growth[i]                          # large-tree htgf value for the blend
        xrdgro = active_multiplier(s.control, :regd, sp, current_cycle_year(s))
        bark = bark_ratio(c.bark_a, c.bark_b, sp, d)        # BRATIO(DBH(K)=D, HT(K)=H)
        small_d = d < 3.0f0
        # deterministic dub base D1 (kt/regent.f:544-555)
        relh = (h - 4.5f0) / (ah - 4.5f0); relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
        dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
        d1 = KT_RG_DIAM[sp] + dadj
        if h > 4.5f0
            d1 = sp == 11 ? 0.0729f0*fpow(h - 4.5f0, 1.1988f0) + dadj : KT_RG_HCON[sp]*h + KT_RG_DCON[sp] + dadj
        end
        nrec = stash !== nothing ? 3 : 1
        central_dbh = d
        for l in 0:(nrec - 1)
            zzran = 0.0f0
            if dgsd >= 1.0f0
                while true
                    zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                    (zzran <= 1.0f0 && zzran >= -1.5f0) && break
                end
            end
            htgr = htgr1 + zzran * KT_RG_HSIGMA; htgr < 0.15f0 && (htgr = 0.15f0)
            # kt/regent.f:522 HTG(K)=HTGR*(1-XWT)+XWT*HTG(K): K is the copy's own slot, whose large-tree HTG
            # kt_triple_htg! (htgf.f:139-161) put in htgU/htgL — read before this loop overwrites it below.
            lh = l == 0 ? large_htg : (stash.htg_copy[i] ? (l == 1 ? stash.htgU[i] : stash.htgL[i]) : large_htg)
            htg = htgr * (1.0f0 - xwt) + xwt * lh
            (h + htg > cap) && (htg = max(cap - h, 0.1f0))
            dg_inc = 0.0f0; dbh_dir = -1.0f0
            if small_d
                hk = h + htg
                if hk < 4.5f0
                    dbh_dir = 0.1f0 + KT_RG_DIAM[sp]*0.01f0 + hk*0.001f0   # kt/regent.f:562 (DBH set, DG=0)
                else
                    dk = sp == 11 ? 0.0729f0*fpow(hk - 4.5f0, 1.1988f0) + dadj :
                                    KT_RG_HCON[sp]*hk + KT_RG_DCON[sp] + dadj
                    dk < KT_RG_DIAM[sp] && (dk = KT_RG_DIAM[sp])           # kt/regent.f:582
                    dk = dk + hk*0.001f0
                    dgk = (dk - d1) * xrdgro; dgk < 0.0f0 && (dgk = 0.0f0) # kt/regent.f:591 (DK−D1)·XRDGRO
                    dg0 = dgk * bark                                        # DG(K)=DGK·BARK (:599)
                    dds = dg0 * (2.0f0*bark*d + dg0) * scale                # :600
                    dg_inc = sqrt((d*bark)^2 + dds) - bark*d              # :601
                    (d + dg_inc) < KT_RG_DIAM[sp] && (dg_inc = KT_RG_DIAM[sp] - d)  # :603 DIAM floor
                    dg_inc = dg_bound(nothing, nothing, sp, d, dg_inc, s.control.sp_size_cap)  # DGBND :610
                end
            end
            if l == 0
                t.ht_growth[i] = htg
                if small_d
                    if dbh_dir >= 0.0f0
                        t.dbh[i] = dbh_dir; t.diam_growth[i] = 0.0f0; central_dbh = dbh_dir
                    else
                        t.diam_growth[i] = dg_inc                           # increment; DBH via GRADD
                    end
                end
            elseif l == 1
                stash.htgU[i] = htg; stash.is_small[i] = true
                small_d && (stash.dgU[i] = dbh_dir >= 0.0f0 ? (dbh_dir - central_dbh)*bark : dg_inc)
            else
                stash.htgL[i] = htg
                small_d && (stash.dgL[i] = dbh_dir >= 0.0f0 ? (dbh_dir - central_dbh)*bark : dg_inc)
            end
            # kt/regent.f:614-650 — a SMALL record (D<3; `IF(D.GE.3.0) GO TO 23` skips the rest) that reaches DBH≥3
            # this cycle (DNEW=D+DG on the cycle's DDS scale) gets a fresh DUBSCR crown (one FCR draw), capped at the
            # all-HTG-to-crown CRMAX. The tripled copies run it too: kt/dgdriv.f:253/261 already put the central DBH in
            # their slots, so they draw as well — but their ICR(K) is overwritten by TRIPLE's ICR(ITFN)=ICR(I), so for
            # them only the RNG draw survives. (Measured on ktt01 cycle 2: live dubs K=50,180,181,12,104,105,51,…)
            if small_d
                dK = dbh_dir >= 0.0f0 ? dbh_dir : central_dbh           # DBH(K): the HK<4.5 direct set, else DBH(I)
                dgK = dbh_dir >= 0.0f0 ? 0.0f0 : dg_inc
                barkK = bark_ratio(c.bark_a, c.bark_b, sp, dK)
                dds2 = dgK * (2.0f0 * barkK * dK + dgK) * (fint / 10.0f0)       # SCALE2=FINT/YR
                dg2 = sqrt((dK * barkK)^2 + dds2) - barkK * dK; dg2 < 0.0f0 && (dg2 = 0.0f0)
                if dK + dg2 >= 3.0f0
                    crf = kt_dubscr_cr(sp, dK + dg2, h + htg, s.plot.basal_area, dgsd, s.rng)
                    if l == 0
                        temcr = crf * 100.0f0 + 0.5f0
                        icrK = Int(t.crown_pct[i])
                        if icrK != 0
                            crln = h * Float32(icrK) / 100.0f0
                            crmax = (crln + htg) / (h + htg) * 100.0f0
                            temcr > crmax && (temcr = crmax)
                        end
                        t.crown_pct[i] = Int32(trunc(Int, temcr))
                    end
                end
            end
        end
    end
    return s
end

# kt/dgdriv.f DO 220 (:707-731) DG(I) — measured (capped at inside-bark DBH when IDG<2), 0 at HT<=4.5, else the dub
# from the second DGF(WK3) (:705, COR final; dub_wk2/dub_wk3 stash). Read by the LSTART REGCAL DO 49.
@inline function kt_do220_dg(s::StandState, i::Int, dcur::Float32)::Float32
    t, c = s.trees, s.calib
    sp = Int(t.species[i])
    bark = bark_ratio(c.bark_a, c.bark_b, sp, dcur)
    if t.diam_growth[i] > 0f0 && t.height[i] > 4.5f0
        dg = t.diam_growth[i]
        (s.control.growth_idg < 2 && dg > dcur * bark) && (dg = dcur * bark)
        return dg
    elseif t.height[i] <= 4.5f0 || length(c.dub_wk2) < i
        return 0f0
    end
    sc = s.control.growth_fint / 10f0
    dd = c.dub_wk3[i] * bark
    dub = sqrt(dd * dd + fexp(c.dub_wk2[i] + t.old_random[i]) * sc) - dd
    dub > dd && (dub = dd)
    return dg_bound(nothing, nothing, sp, dcur, dub, s.control.sp_size_cap)
end

"""
    kt_regent_hcor_init!(s, isct, ind1, saved_dbh)

kt/regent.f:674-848 — the LSTART small-tree HEIGHT calibration (REGENT(.FALSE.,1) from kt/cratet.f:648), which jl had
never ported (HCOR_init 0 for every KT species). Per species with >= NCALHT(5) sub-5" records carrying a measured
HTG, EDH = HK−H grown over the NPER subcycles by the KT small-tree model:
- species 1-10: EDH = BHAB + BLH·HK + HTH2·HT2MOD·HK² + RHBA·ln(BAJ) + RHBAL·BAL + HTPCC1·MANAGD·PCCF + (HTCR+HTCR2·CR)·CR,
  floored at 0, HK += EDH·CON (CON = RCOR2 under READCORR, else 1); BAL = BAJ·(100−PCT)·0.01;
- species 11 (the NI form): EDH = exp(BHAB + RHLH·ln HK + RHCCF·RDJ + RHBAL·BALMH), BALMH = BAJ·(100−PCT)·0.0001.
CORNEW = Σ(HTG·SCALE3·P)/Σ(EDH·P) (SCALE3 = REGYR/FINTH), HCOR = ln(CORNEW) trapped to [0.0821, 12.1825]. NTYR = IFINTH.
Stand values (TEMBA/TEMCCF, PCCF, PCT) come from the kt/cratet.f:184 backdating DENSE snapshot. ALBA is reset per
species and only updated when BAJ > 0 (as FVS). t.dbh is the backdated WK3 here.
"""
function kt_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    ctl = s.control
    ctl.growth_ifinth == 0 && return s                  # regent.f:681 IF(IFINTH.EQ.0) GOTO 100
    rhcon = kt_regcons!(s)
    if ctl.regh_cor2_on                                 # REGCON: sp11 + ln RCOR2 (1-10 carry RCOR2 as CON below)
        length(ctl.regh_cor2) >= 11 && ctl.regh_cor2[11] > 0f0 && (rhcon[11] += flog(ctl.regh_cor2[11]))
    end
    snap = c.cratet_relden > 0f0
    temba = snap ? c.cratet_ba : p.basal_area; temccf = snap ? c.cratet_relden : p.relative_density
    pccfv = (snap && !isempty(c.cratet_pccf)) ? c.cratet_pccf : s.density.point_ccf
    pctv = (snap && length(c.cratet_pct) >= t.n) ? c.cratet_pct : t.crown_ratio
    finth = ctl.growth_finth > 0f0 ? ctl.growth_finth : 5f0
    scale3 = KT_RG_REGYR / finth
    ntyr = Int(ctl.growth_ifinth); iyr = Int(KT_RG_REGYR)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for k in 1:nper
        if nn == 1; kper[k] = itot; break; end
        kper[k] = itot ÷ nn; itot -= kper[k]; nn -= 1
    end
    banext = fill(temba, nper); rdnext = fill(temccf, nper)
    if nper > 1                                         # DO 49 (regent.f:699-717)
        @inbounds for i in 1:t.n
            d1 = t.dbh[i]; sp = Int(t.species[i]); pr = t.tpa[i]
            d2 = d1 + kt_do220_dg(s, i, saved_dbh[i]) / bark_ratio(c.bark_a, c.bark_b, sp, d1)
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            c1 = kt_tree_ccf(sp, d1) * pr; c2 = kt_tree_ccf(sp, d2) * pr
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn; banext[j] += Float32(k) * bi * pn
            end
        end
    end
    ihtg = ctl.growth_ihtg
    dum1 = p.managed == Int32(1) ? 1f0 : 0f0
    @inbounds for sp in 1:11
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        alba = 0f0
        bhab = rhcon[sp]; blh = KT_RG_RHLH[sp]; bccf = KT_RG_RHCCF[sp]
        hths2 = KT_RG_HTH2[sp] * KT_RG_HT2MOD[sp]; bba = KT_RG_RHBA[sp]
        htcrs = KT_RG_HTCR[sp]; htcrs2 = KT_RG_HTCR2[sp]
        bbal = KT_RG_RHBAL[sp]; bbalmh = 0f0
        sp == 11 && (bbalmh = bbal; bbal = 0f0)
        con = (ctl.regh_cor2_on && sp <= length(ctl.regh_cor2) && ctl.regh_cor2[sp] > 0f0) ? ctl.regh_cor2[sp] : 1f0
        htpc1 = KT_RG_HTPCC1[sp] * dum1
        snp = 0f0; snx = 0f0; sny = 0f0; nh = 0
        for k in i1:i2
            i = Int(ind1[k])
            hg = t.ht_growth[i]
            h = t.height[i]; ihtg < 2 && (h -= hg)
            (saved_dbh[i] >= 5f0 || h < 0.01f0) && continue
            hg < 0.001f0 && continue
            ipccf = Int(t.plot_id[i])
            pccf1 = (1 <= ipccf <= length(pccfv)) ? pccfv[ipccf] : 0f0
            cr = Float32(t.crown_pct[i]) / 100f0
            hk = h
            for j in 1:nper
                baj = banext[j]; rdj = rdnext[j]
                baj > 0f0 && (alba = flog(baj))
                bal = baj * (100f0 - pctv[i]) * 0.01f0
                balmh = baj * (100f0 - pctv[i]) * 0.0001f0
                if sp != 11
                    edh = bhab + blh * hk + hths2 * hk * hk + bba * alba + bbal * bal + htpc1 * pccf1 +
                          (htcrs + htcrs2 * cr) * cr
                    edh < 0f0 && (edh = 0f0)
                    hk += edh * con
                else
                    edh = fexp(bhab + blh * flog(hk) + bccf * rdj + bbalmh * balmh)
                    edh < 0f0 && (edh = 0f0)
                    hk += edh
                end
            end
            edh = hk - h
            pr = t.tpa[i]
            snp += pr; snx += edh * pr; sny += hg * scale3 * pr; nh += 1
        end
        nh < 5 && continue                                         # NCALHT
        snx /= snp; sny /= snp
        cornew = sny / snx
        cornew <= 0f0 && (cornew = 1f-4)
        c.htg_cor_init[sp] = (cornew < 0.0821f0 || cornew > 12.1825f0) ? 0f0 : flog(cornew)
    end
    return s
end

"""
    kt_esgent!(s, nstart; fint, atavh, atba, atrelden, relden_pre, ba_pre, pccf_pre)

estb/esgent.f (KT) — SPESRT, REGENT(.TRUE.,ITRNIN) (small_tree_growth! in ESTAB mode) for the records ESTAB created
this cycle (nstart+1:n), then DO 100: HTEMP=HT+HTG; HTG=HTG·WK4; HT=HT+HTG; WK4<1 ⇒ HT<4.5: DBH=0.1+0.001·HT, DG=0,
else DBH and DG scaled by HT/HTEMP; HT>HHTMAX ⇒ HT=HHTMAX, DBH=2.95. estab.f:1504 then adds GENTIM to ABIRTH.
"""
function kt_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atavh::Float32 = -1.0f0, atba::Float32 = -1.0f0, atrelden::Float32 = -1.0f0,
                    relden_pre::Float32 = -1.0f0, ba_pre::Float32 = -1.0f0, pccf_pre::Vector{Float32} = Float32[])
    t = s.trees
    nstart >= t.n && return s
    species_sort!(s)                                   # esgent.f CALL SPESRT
    small_tree_growth!(s, nothing, s.variant; fint = fint, lestb = true, itrnin = nstart + 1,
                       atba = atba, atccf = atrelden, atavh = atavh, ba_now = ba_pre, relden_now = relden_pre,
                       pccf_now = pccf_pre)
    @inbounds for i in (nstart+1):t.n
        sp = Int(t.species[i])
        htemp = t.height[i] + t.ht_growth[i]
        wk4 = t.htimlt[i]
        t.ht_growth[i] = t.ht_growth[i] * wk4
        t.height[i] = t.height[i] + t.ht_growth[i]
        if wk4 < 1f0
            if t.height[i] < 4.5f0
                t.dbh[i] = 0.1f0 + 0.001f0 * t.height[i]; t.diam_growth[i] = 0f0
            else
                t.dbh[i] = t.dbh[i] * (t.height[i] / htemp)
                t.diam_growth[i] = t.diam_growth[i] * (t.height[i] / htemp)
            end
        end
        if t.height[i] > _KT_ES_HHTMAX[sp]
            t.height[i] = _KT_ES_HHTMAX[sp]; t.dbh[i] = 2.95f0
        end
    end
    esgent_add_gentim!(s, nstart, fint)                # estab.f:1504 ABIRTH += GENTIM (after ESGENT)
    return s
end
