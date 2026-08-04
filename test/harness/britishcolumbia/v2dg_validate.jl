# BC V2 (LV2ATV) large-tree DG formula validation (dgf.f:1940 + 1992-1996). Replays the oracle BCV2DG dump
# (ESSFdk/01) through bc_v2_dds using the extracted per-species coeffs + the oracle's CONSPP/DGDSQ (the
# habitat-class CONSPP resolution is a separate step). Validates the V2 DDS formula + DGLD/DGCR/DGCRSQ/DGBAL/DGDBAL.
const OUT = "/tmp/essf.out"
f32(x) = parse(Float32, x)

# V2 per-species DG coefficients (dgf.f DATA, 15 species)
const DGLD   = Float32[0.56445,0.54140,0.56888,0.68810,0.68712,0.58705,0.89503,0.73045,0.86240,0.66101,0.89778,0.89778,0.89778,0.56888,0.89778]
const DGCR   = Float32[1.08338,1.03478,2.06850,1.93969,1.64133,1.29360,1.85558,1.54643,0.52044,1.31618,1.28403,1.28403,1.28403,2.06850,1.28403]
const DGCRSQ = Float32[0.0,0.07509,-0.62361,-0.78258,-0.27244,0.0,-0.36393,-0.26635,0.86236,0.0,0.0,0.0,0.0,-0.62361,0.0]
const DGBAL  = Float32[0.42112,0.43637,0.50202,0.45142,0.0,0.74596,-0.03665,0.25639,0.0,0.0,0.0,0.0,0.0,0.50202,0.0]
const DGDBAL = Float32[-2.08272,-2.03256,-2.11590,-1.76812,-0.80918,-2.28375,-0.43329,-1.18218,-0.51270,-1.25881,-0.66110,-0.66110,-0.66110,-2.11590,-0.66110]

"""V2 large-tree DDS (dgf.f:1940 + 1992-1995). `d` inches, `conspp`/`dgdsq` from the (separate) CONSPP resolution."""
function bc_v2_dds(sp, d, bal, cr, conspp, dgdsq)
    ald = log(d)
    dds = conspp + DGLD[sp]*ald + DGBAL[sp]*bal + cr*(DGCR[sp] + cr*DGCRSQ[sp]) +
          dgdsq*d*d + DGDBAL[sp]*bal/log(d + 1f0)
    return max(-9.21f0, log(max(0.001f0, dds)))
end

function main()
    rows = [split(l) for l in readlines(OUT) if startswith(l, "BCV2DG ")]
    n = 0; ex = 0; mx = 0f0
    for r in rows
        sp = parse(Int, r[2]); d = f32(r[3]); bal = f32(r[4]); cr = f32(r[5]); conspp = f32(r[6])
        dgdsq = f32(r[11]); wk2o = f32(r[13])
        my = bc_v2_dds(sp, d, bal, cr, conspp, dgdsq)
        n += 1; my == wk2o && (ex += 1); mx = max(mx, abs(my - wk2o))
    end
    println("BC V2 DG formula replay (ESSF, sp14): $ex/$n bit-exact ($(round(100ex/n,digits=1))%)  maxerr=$mx")
end
main()
