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
    ca = cos(asp); sa = sin(asp)
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
function small_tree_growth!(s::StandState, stash, ::Kootenai; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    t.n == 0 && return s
    n = t.n
    rhcon = kt_regcons!(s)
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    managed = p.managed == Int32(1)
    dgsd = s.control.dg_sd
    regyr = KT_RG_REGYR                                    # 5.0
    # subcycle count + lengths (regent.f:186-203)
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    kper[nper] = itot == ntyr ? kper[nper] : itot          # KPER(NPER)=ITOT (regent.f:203)
    # recompute ITOT-based last (mirror FVS exactly): after the loop ITOT holds the remainder
    # (the loop already assigned kper[nper]=itot when nn hit 1)
    # per-subcycle stand density from the LARGE trees (regent.f:236-256)
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            c1 = kt_tree_ccf(sp, d1); c2 = kt_tree_ccf(sp, d2)   # CCFCAL per-tree (no tpa) — matches C1/C2 in regent
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]
                pn = pr * 0.985f0^k
                rdnext[j] += k * ci / pr * pn
                banext[j] += k * bi * pn
            end
        end
    end
    # DELMAX (regent.f:313), R=RELDEN, AH=AVH
    delmax = (avh / 36.0f0) * (0.01232f0 * relden - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    # per-tree height accumulator (WK3), starts at HT
    wk3 = Float32[t.height[i] for i in 1:n]
    ah = avh
    # ---- subcycle height loop (regent.f:321-460) ----
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]
        alba = baj > 0.0f0 ? log(baj) : 0.0f0
        kpj = Float32(kper[j])
        for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]
            d >= KT_RG_XMAX[sp] && continue
            t.tpa[i] <= 0.0f0 && continue
            con = exp(c.htg_cor_small[sp])                # RHCON·EXP(HCOR); HCOR=0 until KT calib branch
            h1 = wk3[i]
            pct = t.crown_ratio[i]
            bal = baj * (100.0f0 - pct) * 0.01f0
            balmh = baj * (100.0f0 - pct) * 0.0001f0
            cr = Float32(t.crown_pct[i]) / 100.0f0
            pt = Int(t.plot_id[i]); pccf1 = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
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
                htgrl = exp(rhcon[sp] + KT_RG_RHLH[sp]*log(h1) + KT_RG_RHCCF[sp]*rdj + KT_RG_RHBAL[sp]*balmh)
                htgrl < 0.0f0 && (htgrl = 0.0f0)
                wk3[i] = h1 + htgrl * (kpj / regyr) * xrhgro
            end
        end
    end
    # ---- final: HTGR1 + ZZRAN + XWT blend + DG dub (regent.f:473-560) ----
    _sp_order = sortperm(view(t.species, 1:n); alg = Base.Sort.MergeSort)
    @inbounds for oi in 1:n
        i = _sp_order[oi]
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= KT_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        h = t.height[i]
        xmn = KT_RG_XMIN[sp]; xmx = KT_RG_XMAX[sp]
        htgr1 = wk3[i] - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 1.0f0 && zzran >= -1.5f0) && break
            end
        end
        htgr = htgr1 + zzran * KT_RG_HSIGMA; htgr < 0.15f0 && (htgr = 0.15f0)
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # diameter dub for D<3 (regent.f:534-560)
        if d < 3.0f0
            relh = (h - 4.5f0) / (ah - 4.5f0); relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
            dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
            d1 = KT_RG_DIAM[sp] + dadj
            if sp == 11
                h > 4.5f0 && (d1 = 0.0729f0*(h - 4.5f0)^1.1988f0 + dadj)
            else
                h > 4.5f0 && (d1 = KT_RG_HCON[sp]*h + KT_RG_DCON[sp] + dadj)
            end
            hk = h + htg
            local d2::Float32
            if sp == 11
                d2 = hk > 4.5f0 ? 0.0729f0*(hk - 4.5f0)^1.1988f0 + dadj : KT_RG_DIAM[sp] + dadj
            else
                d2 = hk > 4.5f0 ? KT_RG_HCON[sp]*hk + KT_RG_DCON[sp] + dadj : KT_RG_DIAM[sp] + dadj
            end
            xrdgro = active_multiplier(s.control, :regd, sp, current_cycle_year(s))
            dg = (d2 - d1) * xrdgro; dg < 0.0f0 && (dg = 0.0f0)
            t.diam_growth[i] = dg
        end
    end
    return s
end
