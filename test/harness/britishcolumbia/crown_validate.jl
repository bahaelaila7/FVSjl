# BC chunk-5b crown-ratio validation: replay oracle-dumped crown inputs (BCCRSP per-species site coeffs +
# BCCRTREE per-tree) through bc_crnmd + the crown-change update, compare EXPPCR/EXPDCR/ICRI bit-exactly.
using FVSjl
const OUT = "/workspace/.bcwork/mrun/all_BC.out"
f32(x) = parse(Float32, x)

function main()
    lines = readlines(OUT)
    # per (cycle,sp) resolved coeffs + stand terms
    sp = Dict{Tuple{Int,Int},NamedTuple}()
    for l in lines
        startswith(l, "BCCRSP") || continue
        p = split(l)
        sp[(parse(Int,p[2]), parse(Int,p[3]))] = (crcon=f32(p[4]), crhtdbh=f32(p[5]), crht=f32(p[6]),
            crdbh2=f32(p[7]), crbal=f32(p[8]), crlnccf=f32(p[9]), xcrcon=f32(p[10]), dcrcon=f32(p[11]),
            relden=f32(p[12]), rdm1=f32(p[13]), ba=f32(p[14]), oba=f32(p[15]))
    end
    # assert my extracted coefficient constants match the oracle's resolved slopes (sp14, ICH)
    z = FVSjl._BC_CR_ICH; col = findfirst(==(14), z.seq)
    r14 = first(v for ((cy,s),v) in sp if s==14)
    okc = z.htdbh[col]==r14.crhtdbh && z.ht[col]==r14.crht && z.dbh2[col]==r14.crdbh2 &&
          z.bal[col]==r14.crbal && z.lnccf[col]==r14.crlnccf
    println("sp14 constant coeffs (htdbh/ht/dbh2/bal/lnccf) match oracle: ", okc)

    nE=0; exP=0; exD=0; exI=0; mxP=0f0; mxD=0f0
    for l in lines
        startswith(l, "BCCRTREE") || continue
        p = split(l)
        cyc=parse(Int,p[2]); s=parse(Int,p[4]); dbh=f32(p[5]); ht=f32(p[6]); dg=f32(p[7]); htg=f32(p[8])
        pct=f32(p[9]); oldpct=f32(p[10]); icr=parse(Int,p[11]); ppO=f32(p[12]); pdO=f32(p[13]); icriO=parse(Int,p[14])
        st = sp[(cyc,s)]
        bark = FVSjl.bc_bratio(s)
        # FORWARD (EXPPCR): use oracle XCRCON + my crnmd
        pp = pct < 0.01f0 ? 0.01f0 : pct
        balf = (1f0 - pp/100f0) * st.ba
        exppcr = FVSjl.bc_crnmd(st.xcrcon, st.crhtdbh, st.crht, st.crdbh2, st.crbal, dbh, ht, balf)
        # BACKWARD (EXPDCR): backdate D,H,P + oracle DCRCON/OBA
        db = dbh - dg/bark; db <= 0f0 && (db = dbh)
        hb = ht - htg;      hb <= 0f0 && (hb = ht)
        pb = oldpct <= 0f0 ? pct : oldpct; pb < 0.01f0 && (pb = 0.01f0)
        balb = (1f0 - pb/100f0) * st.oba
        expdcr = FVSjl.bc_crnmd(st.xcrcon, st.crhtdbh, st.crht, st.crdbh2, st.crbal, db, hb, balb)  # V3 uses XCRCON, not DCRCON
        chg = exppcr - expdcr
        if icr > 0
            pdifpy = chg / Float32(icr) / 10f0 * 100f0
            pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * 10f0 / 100f0)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * 10f0 / 100f0)
        end
        icri = trunc(Int, Float32(icr) + chg*100f0 + 0.50005f0)
        icri > 95 && (icri = 95); icri < 5 && (icri = 5)
        nE += 1
        exppcr == ppO && (exP += 1); mxP = max(mxP, abs(exppcr-ppO))
        expdcr == pdO && (exD += 1); mxD = max(mxD, abs(expdcr-pdO))
        icri == icriO && (exI += 1)
    end
    println("BC crown-ratio replay (all_BC, all cycles, V3 CRNMD):")
    println("  EXPPCR: $exP/$nE bit-exact ($(round(100*exP/nE,digits=1))%)  maxerr=$mxP")
    println("  EXPDCR: $exD/$nE bit-exact ($(round(100*exD/nE,digits=1))%)  maxerr=$mxD")
    println("  ICRI  : $exI/$nE bit-exact ($(round(100*exI/nE,digits=1))%)")
end
main()
