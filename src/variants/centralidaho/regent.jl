# =============================================================================
# regent.jl (centralidaho) — CI small-tree growth (ci/regent.f). Chunk 6. Overrides DG/HTG for
# small trees (D < XMAX). NIVAR path (CI conifers): per-subcycle height via the HTGRL regression,
# accumulated, then DBH from HTDBH; DG = new−old DBH. Coeffs from data/centralidaho/regent_coeffs_1d.csv.
#   HTGRL = CNST + HTBA·RELHT·PTBAA + PBAL·PTBALI + STBA·BA + RLHT·RELHT + CRSQ·RCR² + CR·RCR
#         + BAL·TBAL + HDM1·RHDM1 + HDM2·RHDM2 + PTBA·PTBAA   (ci/regent.f:641)
#   CON = RHCON·exp(HCOR)  (RHCON≈1, HCOR = c.htg_cor_small); H2 = H1 + HTGRL·SCALE·XRHGRO·CON
#   DK = exp(DHCN + DHHT·ln(HK) + DHCR·ln(RCR))  (HTDBH); DG = DK − D_old
# SMHTGF (stochastic BETA/FINDAG) is the establishment/aspen path — deferred; cit01 small trees are conifers.
# =============================================================================

let
    path = joinpath(CI_DATADIR, "regent_coeffs_1d.csv")
    hdr = split(strip(readlines(path)[1]), ',')
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    ci(n) = findfirst(==(n), hdr)
    col(n) = Float32[parse(Float32, strip(rows[sp][ci(n)])) for sp in 1:19]
    for n in ("CNST","CR","CRSQ","DGMAX","DHCN","DHCR","DHHT","HDM1","HDM2","HTBA","RLHT","STBA",
              "XMAX","XMIN","BAL","PBAL","PTBA")
        @eval global const $(Symbol("CI_RG_" * n)) = $(col(n))
    end
end
const CI_RG_REGYR = 5.0f0
const CI_RG_DIAM = Float32[0.4,0.3,0.3,0.3,0.2,0.2,0.4,0.3,0.3,0.5,0.3,0.3,0.2,0.3,0.2,0.3,0.2,0.2,0.2]  # ci/regent.f DATA DIAM (min DBH)

function small_tree_growth!(s::StandState, stash, ::CentralIdaho; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species
    t.n == 0 && return s
    n = t.n
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    kodtyp = Int(p.habitat_code)
    rhdm1 = (500 <= kodtyp < 600) ? 1.0f0 : 0.0f0
    rhdm2 = (600 <= kodtyp < 700) ? 1.0f0 : 0.0f0
    regyr = CI_RG_REGYR
    ntyr = Int(round(fint))                               # NTYR (grow cycle; LSTART/ESTAB partial-cycle deferred)
    scale_h = Float32(ntyr) / regyr                       # ci/regent.f:518 CIVAR SCALE = NTYR/REGYR (single pass)
    scale2 = 1.0f0                                        # SCALE2 = YR/NTYR = 1 for a full cycle ⇒ DDS round-trip identity
    cur_year = current_cycle_year(s)
    # CIVAR (sp 1-10,18) small-tree height+DBH — single pass (regent.f:528 J>1 skips). Verified vs the PRISTINE
    # FVSci_clean via the DEBUG keyword: seedling I=24 H=1.01→HK=5.267, DK=0.5179, DKK=0.1 ⇒ DG=0.387 (matches).
    # ZZRAN random height (regent.f:934) DEFERRED = deterministic mean HTGR=HTGR1 (grow-phase RNG desync class).
    # TT(11,12,16)=SMHTGF and UT(13,14,15,17,19)=aspen/juniper get the CIVAR form as a fallback (rare/deferred).
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d0 = t.dbh[i]
        d0 >= CI_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        (sp < 1 || sp > 19) && continue
        h0 = t.height[i]
        pt = Int(t.plot_id[i])
        ptba = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : ba
        pct = t.crown_ratio[i]
        relht = avh > 0f0 ? h0 / avh : 1.5f0; relht > 1.5f0 && (relht = 1.5f0)
        ptbali = ptba * (1.0f0 - pct / 100.0f0)
        tbal = (1.0f0 - pct / 100.0f0) * ba
        iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9); iicr < 1 && (iicr = 1)
        rcr = Float32(iicr)
        htgrl = CI_RG_CNST[sp] + CI_RG_HTBA[sp]*relht*ptba + CI_RG_PBAL[sp]*ptbali +
                CI_RG_STBA[sp]*ba + CI_RG_RLHT[sp]*relht + CI_RG_CRSQ[sp]*rcr*rcr +
                CI_RG_CR[sp]*rcr + CI_RG_BAL[sp]*tbal + CI_RG_HDM1[sp]*rhdm1 +
                CI_RG_HDM2[sp]*rhdm2 + CI_RG_PTBA[sp]*ptba
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        con = exp(c.htg_cor_small[sp])                    # CON = RHCON·exp(HCOR); RHCON=1
        h2 = h0 + htgrl * scale_h * xrhgro * con
        htgr1 = h2 - h0; htgr1 < 0.0f0 && (htgr1 = 0.0f0) # HTGR1 = HK−H (regent.f:907-908)
        htgr = htgr1                                      # + ZZRAN*0.1*SCALE deferred
        htgr < 0.1f0 && (htgr = 0.1f0)                    # Dixon 3/4/09 floor (regent.f:961)
        htg = htgr                                        # XWT=0 for D≤XMN small trees ⇒ HTG=HTGR
        t.ht_growth[i] = htg
        hk = h0 + htg                                     # HK = H + HTG (regent.f:1002)
        if hk < 4.5f0
            t.diam_growth[i] = 0.0f0
            continue
        end
        dhcn = CI_RG_DHCN[sp]; dhht = CI_RG_DHHT[sp]; dhcr = CI_RG_DHCR[sp]
        dk = exp(dhcn + dhht * log(hk) + dhcr * log(rcr))
        dkk = h0 < 4.5f0 ? d0 : exp(dhcn + dhht * log(h0) + dhcr * log(rcr))   # regent.f:1009-1012
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        bark = ci_bratio(sd, sp, d0)                      # BRATIO(ISPC,DBH,HT)
        dg = (dk - dkk) * bark * xrdgro; dg < 0.0f0 && (dg = 0.0f0)            # regent.f:1245
        dds = dg * (2.0f0 * bark * d0 + dg) * scale2                          # regent.f:1251
        dg = sqrt((d0 * bark)^2 + dds) - bark * d0                            # regent.f:1253
        dg < 0.0f0 && (dg = 0.0f0)
        (d0 + dg) < CI_RG_DIAM[sp] && (dg = CI_RG_DIAM[sp] - d0)              # regent.f:1256 min-DBH floor
        t.diam_growth[i] = dg
    end
    return s
end
