# =============================================================================
# volume.jl (easternmontana) — EM volume. Conifers = Region-1 Flewelling FW2 (I00FW2W<fia>, same kernel as
# KT via cr_fw2_vol). Non-conifers = DVEW woodland (em VEQNNC): region-1 codes (102/101DVEW) route through
# NVEL dvest.f → R1KEMP (VOLEQ 2:3='02') / R1ALLEN ('01'); region-2 (200DVEW) → R2OLDV/cr_dve_vol.
# This chunk ports the R1KEMP cubic (r1kemp.f) for the non-conifer FIA codes present in the EM set.
# =============================================================================

# R1KEMP cubic coefficients CBVOLE(ISPEC, 1:8) (volume/NVEL/r1kemp.f DATA CBVOLE), for the EM DVEW species.
# ISPEC map (r1kemp.f:178-204): 746→1 (aspen), 740→2 (cottonwood).
const _EM_R1KEMP_CB = Dict{String,NTuple{8,Float32}}(
    "746" => (0.3482f0, -0.0384f0, 0.001427f0, -0.842503f0, 0.224f0, -0.343f0, 0.217f0, 1.071f0),
    "740" => (0.1064f0, -0.00778f0, 0.000176f0, -0.265342f0, 0.204f0, -0.749f0, 0.194f0, 4.285f0),
)
# R1KEMP board-foot BFVOL(ISPEC,1,1:4) (r1kemp.f:39-40), JTAB=1 table (the '02' species).
const _EM_R1KEMP_BF = Dict{String,NTuple{4,Float32}}(
    "746" => (1.197f0, -18.544f0, 1.216f0, -21.309f0),
    "740" => (1.046f0, -15.966f0, 1.140f0, -46.735f0),
)

# R1KEMP gross board-foot (r1kemp.f:300, JTAB=1): BFGRS = BFVOL·D2H100 + intercept (DBH split at 21).
@inline function _em_r1kemp_board(fia::AbstractString, d::Float32, h::Float32)::Float32
    bf = get(_EM_R1KEMP_BF, fia, nothing); bf === nothing && return 0f0
    d2h100 = d * d * h / 100f0
    bfgrs = d < 21f0 ? bf[1] * d2h100 + bf[2] : bf[3] * d2h100 + bf[4]
    return bfgrs < 0f0 ? 0f0 : bfgrs
end

# R1KEMP gross cubic (r1kemp.f:345-384, ISPEC≠14,15; KLASS=1 default). D2H100=DBH²·HT/100.
@inline function _em_r1kemp_cubic(fia::AbstractString, d::Float32, h::Float32)::Float32
    cb = get(_EM_R1KEMP_CB, fia, nothing); cb === nothing && return 0f0
    d2h100 = d * d * h / 100f0
    cbgrs = if d < 5f0
        0f0                         # coeffs 9-11 are 0 for these species ⇒ CBGRS=0 for DBH<5 (→ min 1.6)
    elseif d <= 9.5f0
        d2h100 * (cb[1] * d + cb[2] * d * d + cb[3] * d * d * d + cb[4])
    elseif d <= 20.5f0
        cb[5] * d2h100 + cb[6]
    else
        cb[7] * d2h100 + cb[8]
    end
    cbgrs < 1.6f0 && (cbgrs = 1.6f0)      # KLASS≤2 gross cubic minimum (r1kemp.f:378-379)
    return cbgrs
end

function compute_volumes_em!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    topd = 4.5f0; bftopd = 4.5f0; stump = 1.0f0
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]
        if startswith(eq, "I")                      # conifers — Region-1 FW2 (same as KT)
            dbhmin = sp == 7 ? 6f0 : 7f0
            bfmind = sp == 7 ? 6f0 : 7f0
            bark = bark_ratio(ba_a, ba_b, sp, d)
            v = cr_fw2_vol(eq, d, h; bark = bark, topd = topd, bftopd = bftopd, stump = stump, iregn = 1)
            tcf = max(v[1], 0f0)
            mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
            bf  = d >= bfmind ? v[2] : 0f0
            if t.trunc[i] > 0 && tcf > 0f0 && h >= 4.5f0
                vmax = tcf
                tcf, mcf = cr_cftopk(tcf, mcf, d, h, vmax, bark, Int(t.trunc[i]), stump, topd)
                bf = cr_bftopk(bf, d, h, vmax, bark, Int(t.trunc[i]), stump, bftopd)
            end
            t.cuft_vol[i] = tcf; t.merch_cuft_vol[i] = mcf
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        elseif length(eq) >= 6 && eq[4:6] == "DVE"   # non-conifers — DVEW woodland
            fia = strip(eq)[8:10]
            # dvest.f region-1: VOLEQ(2:3)='02'→R1KEMP; '01'→R1ALLEN (TODO); region-2 (VOLEQ(1:1)='2')→R2OLDV (TODO).
            r1kemp = eq[1] == '1' && eq[2:3] == "02"
            tcf = r1kemp ? _em_r1kemp_cubic(fia, d, h) : 0f0
            t.cuft_vol[i] = max(tcf, 0f0)
            t.merch_cuft_vol[i] = max(tcf, 0f0)      # R1KEMP: VOL(1)=VOL(4) (no separate merch trim)
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = (r1kemp && d >= 7f0) ? _em_r1kemp_board(fia, d, h) : 0f0   # VOL(2)=BFGRS, BFMIND=7
        else
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end
