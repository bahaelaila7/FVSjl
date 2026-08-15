# =============================================================================
# test_op_organon_nwo.jl — OP (Olympic) ORGANON NWO engine A/B vs the live FVSop_clean oracle.
#
# Stand S248112 (opdbg.key, /workspace/.opwork/opdbg.out DGDRIV/HTGF/CROWN/MORTS DEBUG dump).
# Feeds the exact /ORGANON/ FOR EXECUTE buffer (27 trees) + SI_1/SI_2 + ACALIB + MSDI and checks the
# engine outputs against the oracle per-tree dump:
#   • C2 PREPARE : ACALIB(1,1)=0.795231164 (height), ACALIB(2,1)=0.667298734 (crown)  — BIT-EXACT
#   • C3 DG_NWO  : DGRO for the 7 IORG=1 (DF) trees                                    — BIT-EXACT
#   • C4 HG_NWO  : HGRO for all 27 buffer trees                       — bit-exact-or-1-ULP (libm floor)
#   • C5 CROWGRO : CR2  for all 27 buffer trees (NWO CALIB(2) path)   — bit-exact-or-≤2-ULP
#   • C6 PM_NWO  : MORTEXP (DEADEXP) for all 27 buffer trees                            — BIT-EXACT
# ORGANON is DGSD=0 (deterministic) ⇒ these are genuine bit-exact targets; residuals are the
# documented irreducible gfortran/Julia libm exp/log/pow ULP.
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OP ORGANON NWO engine (S248112)" begin
    # (ISP, DBH1, HT1OR, CR1, EXPAN1, IORG) — FOR EXECUTE + ISP/IORG dump
    rows = [
     (11,11.5f0,73.0f0,0.35f0,5.54545212f0,0),(16,0.1f0,4.59999990f0,0.55f0,90.0000076f0,0),
     (15,6.5f0,30.0f0,0.75f0,17.3582478f0,0),(13,7.9f0,75.0f0,0.25f0,11.7510977f0,0),
     (13,8.0f0,63.0f0,0.25f0,11.4591560f0,0),(2,6.19999981f0,38.0f0,0.45f0,19.0787220f0,0),
     (13,8.39999962f0,5.0f0,0.35f0,10.3937941f0,0),(11,9.5f0,60.0f0,0.25f0,8.12616062f0,0),
     (16,4.0f0,20.0f0,0.25f0,30.0000019f0,1),(13,8.19999981f0,65.0f0,0.45f0,10.9069901f0,0),
     (16,1.20000005f0,11.0f0,0.55f0,30.0000019f0,1),(16,1.89999998f0,13.0f0,0.45f0,30.0000019f0,1),
     (2,0.1f0,4.59999990f0,0.65f0,90.0000076f0,0),(2,5.30000019f0,27.0f0,0.65f0,26.1084366f0,0),
     (16,10.0f0,65.0f0,0.35f0,7.33385992f0,1),(2,6.09999990f0,38.0f0,0.75f0,19.7093792f0,0),
     (16,12.6999998f0,67.0f0,0.35f0,4.54700232f0,1),(11,9.60000038f0,60.0f0,0.25f0,7.95774698f0,0),
     (16,10.3999996f0,55.0f0,0.45f0,6.78056669f0,1),(11,8.5f0,62.3888664f0,0.25f0,10.1506720f0,0),
     (2,10.8999996f0,65.0f0,0.65f0,6.17276382f0,0),(16,9.39999962f0,60.0f0,0.35f0,8.29997826f0,1),
     (10,3.20000005f0,17.0f0,0.45f0,30.0000019f0,0),(10,0.1f0,4.59999990f0,0.65f0,30.0000019f0,0),
     (10,5.80000019f0,28.0f0,0.65f0,21.8010101f0,0),(10,5.0f0,25.0f0,0.25f0,29.3354397f0,0),
     (2,6.59999990f0,30.0f0,0.65f0,16.8362274f0,0)]
    n = length(rows)
    buf = F.OrganonBuffer(n); buf.ntrees = n
    isp_fvs = Int[]
    for (i,r) in enumerate(rows)
        buf.species[i] = F.op_organon_fia(r[1])
        buf.dbh1[i]=r[2]; buf.ht1or[i]=r[3]; buf.cr1[i]=r[4]; buf.expan1[i]=r[5]; buf.iorg[i]=r[6]
        push!(isp_fvs, r[1])
    end
    SITE1=98.0f0; SITE2=87.6699982f0; si_1=SITE1-4.5f0; si_2=SITE2-4.5f0; MSDI=950.0f0
    calib1=ones(Float32,11); calib1[1]=0.795231164f0
    calib2=ones(Float32,11); calib2[1]=0.667298734f0
    g = F.op_execute_nwo(buf, isp_fvs; si_1=si_1, si_2=si_2, msdi=MSDI, calib1=calib1, calib2=calib2)

    ref_dgro = Dict(9=>0.453597367f0,11=>0.608625233f0,12=>0.569343150f0,15=>0.825030863f0,
                    17=>0.869497478f0,19=>1.01872468f0,22=>0.806305468f0)
    ref_hgro = Float32[7.16223001,3.27322006,5.61525536,7.10141420,6.43143177,2.90402436,3.37348509,
     6.13480234,5.06713915,7.12914515,4.09005976,4.33374834,1.25820363,3.80295682,6.84004736,
     4.87350750,6.93274069,6.13480234,5.88690805,6.38035202,8.18715763,6.35791016,4.77119637,
     3.27322006,5.53115034,5.46211624,3.92671394]
    ref_cr2 = Float32[0.333584845,0.439911902,0.670865059,0.232286513,0.236058652,0.418279767,
     0.401604533,0.237311244,0.238987386,0.424361765,0.465543568,0.394213021,0.569035053,
     0.598949075,0.335176408,0.703863263,0.339111507,0.244272947,0.430125713,0.233839989,
     0.643511236,0.344468594,0.397868514,0.500091732,0.577729225,0.237465680,0.618761182]
    ref_mort = Float32[2.11330671e-2,18.6668568,0.249865741,0.139908895,0.128600851,0.206817865,
     9.47695896e-2,5.62895387e-2,1.19826448,0.103270188,2.91477108,2.22952390,5.39207745,
     0.339106530,4.07665037e-2,0.218046814,1.36980172e-2,5.20351306e-2,3.24037969e-2,9.18116644e-2,
     2.60991566e-2,5.78977503e-2,1.41518009,6.07447767,0.424153298,0.870218158,0.168391213]

    for (i,r) in ref_dgro; @test g.dgro[i] == r; end                       # DGRO bit-exact
    @test maximum(abs(g.hgro[i]-ref_hgro[i]) for i in 1:n) <= 1f-6         # HGRO ≤1 ULP
    @test maximum(abs(g.cr2[i]-ref_cr2[i])   for i in 1:n) <= 2f-7         # CR2  ≤2 ULP
    for i in 1:n; @test g.deadexp[i] == ref_mort[i]; end                    # MORTEXP bit-exact

    # C2 PREPARE — ACALIB(1,1)/ACALIB(2,1) bit-exact (HT floored to 4.6 as cratet.f:234; tree20 blank)
    prows = [(r[1], r[2], (r[1] in (16,10) && r[3]<=4.6f0 ? 4.59999990f0 : r[3]), r[4], r[5]) for r in rows]
    ht = Float32[p[3] for p in prows]; ht[20] = 0.0f0                       # tree20 originally blank
    species = Int32[F.op_organon_fia(p[1]) for p in prows]
    dbh = Float32[p[2] for p in prows]; cr = Float32[p[4] for p in prows]; ex = Float32[p[5] for p in prows]
    res = F.op_prepare_nwo(species, dbh, ht, cr, ex, n, 1, 60, 54, 98.0f0, 87.6699982f0)
    @test res.acalib[1,1] == 0.795231164f0
    @test res.acalib[2,1] == 0.667298734f0
end
