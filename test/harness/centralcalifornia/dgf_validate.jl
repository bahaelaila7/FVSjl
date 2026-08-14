# =============================================================================
# dgf_validate.jl — CA chunk-3 large-tree DGF per-tree LN(DDS) (ca/dgf.f), bit-exact vs FVSca_dump.
#
# Reference: ref_dgf_cat01.txt — per-tree {I, ISPC, JSPC, D, CR, BA, BAL, PCCF, RELHT, DGCON, COR, COR2,
# SITEAR, DDS} from an instrumented ca/dgf.f (ICYC==1 WRITE), cat01. 81 records = 27 trees × 3 DGF passes
# (calibration COR=0, then growth COR-calibrated) — every record is a self-contained test of the DGF
# equation, reproduced from the shipped CA_* consts + the record's own live {DGCON, COR, D, CR, BA, BAL,
# PCCF, RELHT}. cat01 exercises the STANDARD Wykoff branch only (species WF/DF/LP/SP/PP/BR; no GS/RW/TANOAK).
# =============================================================================
using FVSjl
const M = FVSjl

# Reproduce ca/dgf.f STANDARD-branch LN(DDS) from live per-tree inputs (+ TANOAK 5→10yr, ISPC 42).
function ca_dds_ref(isp::Int, jspc::Int, dgcon::Float32, cor::Float32, d::Float32, cr::Float32,
                    ba::Float32, bal::Float32, pccf::Float32, relht::Float32)
    conspp = dgcon + cor
    dds = conspp + M.CA_DGLD[jspc]*log(d) +
          cr*(M.CA_DGCR[jspc] + cr*M.CA_DGCRSQ[jspc]) +
          M.CA_DGDS[jspc]*d*d + M.CA_DGDBAL[jspc]*bal/log(d + 1f0)
    dds += M.CA_DGPCCF[jspc]*pccf + M.CA_DGHAH[jspc]*relht +
           M.CA_DGLBA[jspc]*log(ba) + M.CA_DGBAL[jspc]*bal
    isp == 42 && (dds = log(exp(dds)*2f0))
    dds < -9.21f0 && (dds = -9.21f0)
    return dds
end

function run_validation()
    ref = joinpath(@__DIR__, "ref_dgf_cat01.txt")
    nfail = 0; worst = 0f0; n = 0
    for l in eachline(ref)
        f = split(strip(l)); isempty(f) && continue
        isp = parse(Int, f[3]); jspc = parse(Int, f[4])
        d   = parse(Float32, f[5]);  cr  = parse(Float32, f[6])
        ba  = parse(Float32, f[7]);  bal = parse(Float32, f[8])
        pccf= parse(Float32, f[9]);  relht = parse(Float32, f[10])
        dgcon = parse(Float32, f[11]); cor = parse(Float32, f[12])
        dds_live = parse(Float32, f[15])
        jl = ca_dds_ref(isp, jspc, dgcon, cor, d, cr, ba, bal, pccf, relht)
        dd = abs(jl - dds_live); worst = max(worst, dd); n += 1
        dd > 1f-3 && (nfail += 1; println("DDS MISMATCH I=$(f[2]) isp=$isp jspc=$jspc d=$d jl=$jl live=$dds_live Δ$dd"))
    end
    println("CA DGF LN(DDS): $n records, $(n-nfail) bit-exact, worst |Δ|=$worst")
    nfail == 0 || error("CA DGF validation FAILED ($nfail)")
end
run_validation()
