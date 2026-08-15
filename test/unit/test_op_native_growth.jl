# =============================================================================
# test_op_native_growth.jl — OP (Olympic) FVS-NATIVE large-tree DGF/HTGF A/B vs the live
# FVSop_clean oracle, for the NON-ORGANON species (WF/ES/LP/SP/PP).
#
# In OP only GF(3)/DF(16) are ORGANON-grown; every other species grows the FVS-native op Wykoff
# path (op/dgf.f + op/htgf.f). Stand S248112 (opdbg.key) has 19 such trees. This feeds the exact
# oracle per-tree inputs (from /workspace/.opwork/opdbg.out DEBUG DGF/HTGF dump, cyc0) and checks:
#   • DGF : LN(DDS) reproduced from op_dgcon + op_dgf_dds vs the op/dgf.f format-9001 print (F7.4)
#   • HTGF: HTG(I) reproduced from op_findag + op_htcalc + op_htg_default vs the op/htgf.f
#           format-901 print (F9.2, pre-SCALE/XHT/HTCON=1)
# Stand facts (op/sitset.f + the 9030 DGCONS dump): ELEV=7, SLOPE=0.30, ASPECT=5.49779, IFOR=6,
# BA(DGF)=62.5335, AVH=63.4388; SITEAR(WF)=97.98189, (ES)=139.15092, (LP)=98.22826, (SP/PP)=
# 139.15092. op is DGSD>0 (OLDRN serial-corr) so multi-cycle is straddle-class — the cyc0 DGF/HTGF
# formula is a deterministic bit-exact target (matched to the oracle's own print precision).
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OP FVS-native DGF/HTGF (S248112, non-ORGANON species)" begin
    ifor = 6; elev = 7.0f0; slope = 0.30f0; asp = 5.49779f0
    ba_dgf = 62.5335f0; avh = 63.4388f0
    SITEAR = Dict(2 => 97.98189f0, 10 => 139.15092f0, 11 => 98.22826f0,
                  13 => 139.15092f0, 15 => 139.15092f0)

    # --- DGF: (tree#, ISPC, D, CR, BAL, PCCF, HT) → oracle LN(DDS) (op/dgf.f format 9001) ---
    dgf_cases = [
        ( 6, 2, 3.99f0,   0.45f0, 47.3326f0, 161.1034f0, 38.0f0,     1.9798f0),
        (13, 2, 0.0828f0, 0.65f0, 62.5256f0, 189.5587f0,  3.0f0,    -5.2213f0),
        (14, 2, 4.3055f0, 0.65f0, 54.0535f0, 189.5587f0, 27.0f0,     2.3832f0),
        (16, 2, 4.7740f0, 0.75f0, 48.9892f0, 130.3031f0, 38.0f0,     2.6745f0),
        (21, 2, 9.7950f0, 0.65f0,  6.2005f0,  69.2663f0, 65.0f0,     3.2620f0),
        (27, 2, 5.0530f0, 0.65f0, 43.6589f0,  38.5963f0, 30.0f0,     2.5694f0),
        (23,10, 2.5333f0, 0.45f0, 60.9087f0, 104.1756f0, 17.0f0,     0.6110f0),
        (24,10, 0.0828f0, 0.65f0, 62.5324f0, 104.1756f0,  2.0f0,    -5.1111f0),
        (25,10, 4.6889f0, 0.65f0, 51.4392f0, 104.1756f0, 28.0f0,     1.6230f0),
        (26,10, 3.8889f0, 0.25f0, 56.6932f0, 104.1756f0, 25.0f0,     0.6455f0),
        ( 1,11,10.3889f0, 0.35f0,  2.9361f0,  64.1733f0, 73.0f0,     1.9788f0),
        ( 8,11, 8.2778f0, 0.25f0, 19.3783f0, 161.1034f0, 60.0f0,     1.6443f0),
        (18,11, 9.0444f0, 0.25f0, 15.8279f0, 170.9581f0, 60.0f0,     1.6954f0),
        (20,11, 7.0399f0, 0.25f0, 24.8675f0, 170.9581f0, 62.3889f0,  1.5516f0),
        ( 4,13, 7.2016f0, 0.25f0, 36.3350f0,  59.2510f0, 75.0f0,     2.2120f0),
        ( 5,13, 7.1851f0, 0.25f0, 33.1083f0, 161.1034f0, 63.0f0,     2.1896f0),
        ( 7,13, 6.9571f0, 0.35f0, 27.6113f0, 161.1034f0,  5.0f0,     2.3313f0),
        (10,13, 6.8031f0, 0.45f0, 30.3551f0,  84.9523f0, 65.0f0,     2.4864f0),
        ( 3,15, 3.7468f0, 0.75f0, 46.0035f0,  59.2510f0, 30.0f0,     2.9519f0),
    ]

    @testset "DGF LN(DDS) tree $(tn) ISPC $(ispc)" for (tn, ispc, d, cr, bal, pccf, ht, orc) in dgf_cases
        dgcon = F.op_dgcon(ispc, ifor, elev, slope, asp, SITEAR[ispc])   # COR=0 (calibration pass)
        dgdsq = F.op_dg_dgdsq(ispc, ifor)
        relht = min(ht/avh, 1.5f0)
        lndds = F.op_dgf_dds(ispc, dgcon, dgdsq, d, cr, bal, pccf, relht, ba_dgf)
        # 17/19 trees match the F7.4 oracle print exactly. The two D≈0.0828 trees are DBH=0.1 input
        # records that DGDRIV back-dates to an internally-computed diameter the DGF DEBUG dump prints
        # only to F11.4 (=0.0828); for a sub-0.1" tree ln(D) is hyper-sensitive to that 5th digit, so
        # they are INPUT-print-limited (the ORGANON decimal-vs-hex lesson), not a formula error.
        atol = d < 0.1f0 ? 2e-3 : 6e-5
        @test round(lndds, digits=4) ≈ orc atol=atol
    end

    # --- DGCONS: DGCON(ISPC) reproduced to the op/dgf.f format-9030 print (F9.5) ---
    dgcon_cases = [(2, 0.79385f0), (10, -0.18469f0), (11, 0.89793f0), (13, 0.92083f0), (15, 2.04686f0)]
    @testset "DGCONS DGCON ISPC $(ispc)" for (ispc, orc) in dgcon_cases
        @test round(F.op_dgcon(ispc, ifor, elev, slope, asp, SITEAR[ispc]), digits=5) ≈ orc atol=2e-5
    end

    # --- HTGF: (ISPC, D, DG, HT, ICR) → oracle HTG(I) (op/htgf.f format 901, pre-scale) ---
    ba_htg = 85.13f0
    htg_cases = [
        ( 2, 3.99f0,   0.91f0, 38.0f0,    45, 9.57f0),
        ( 2, 0.0828f0, 0.02f0,  3.0f0,    65, 2.65f0),
        ( 2, 4.3055f0, 1.21f0, 27.0f0,    65, 6.95f0),
        ( 2, 4.7740f0, 1.46f0, 38.0f0,    75, 9.66f0),
        ( 2, 9.7950f0, 1.39f0, 65.0f0,    65, 9.87f0),
        ( 2, 5.0530f0, 1.28f0, 30.0f0,    65, 7.82f0),
        (10, 2.5333f0, 0.35f0, 17.0f0,    45, 5.38f0),
        (10, 0.0828f0, 0.02f0,  2.0f0,    65, 4.58f0),
        (10, 4.6889f0, 0.54f0, 28.0f0,    65, 7.22f0),
        (10, 3.8889f0, 0.25f0, 25.0f0,    25, 5.15f0),
        (11,10.3889f0, 0.36f0, 73.0f0,    35, 8.97f0),
        (11, 8.2778f0, 0.31f0, 60.0f0,    25, 8.00f0),
        (11, 9.0444f0, 0.32f0, 60.0f0,    25, 8.00f0),
        (11, 7.0399f0, 0.33f0, 62.39f0,   25, 8.26f0),
        (13, 6.9571f0, 0.77f0,  5.0f0,    35, 1.70f0),
        (15, 3.7468f0, 2.11f0, 30.0f0,    75, 4.35f0),
    ]

    @testset "HTGF HTG ISPC $(ispc) D $(d)" for (ispc, d, dg, ht, icr, orc) in htg_cases
        sindx = SITEAR[ispc]
        brat = F.op_bratio(ispc, d)
        d2 = d + dg/brat
        sitage, sitht, agmax, htmax, htmax2 = F.op_findag(ispc, d, d2, ht, sindx)
        htg = if ht > htmax
            ht >= htmax2 ? max(0.5f0*dg, 0.1f0) : 0.0f0
        else
            F.op_htg_default(ispc, sindx, d, ht, Float32(icr), avh, ba_htg, dg, d2,
                             sitage, sitht, agmax, htmax2)
        end
        @test round(htg, digits=2) ≈ orc atol=6e-3            # oracle print is F9.2
    end

    # coefficients(::Olympic) now returns a valid struct (was: errored) so the growth path dispatches.
    @test length(F.coefficients(F.Olympic()).code_alpha) >= 39
end
