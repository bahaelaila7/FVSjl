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

function small_tree_growth!(s::StandState, stash, ::CentralIdaho; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species
    t.n == 0 && return s
    n = t.n
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    kodtyp = Int(p.habitat_code)
    rhdm1 = (500 <= kodtyp < 600) ? 1.0f0 : 0.0f0
    rhdm2 = (600 <= kodtyp < 700) ? 1.0f0 : 0.0f0
    yr = s.control.year
    regyr = CI_RG_REGYR
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    cur_year = current_cycle_year(s)
    wk3 = Float32[t.height[i] for i in 1:n]                # accumulated height (H)
    @inbounds for j in 1:nper
        scale = Float32(kper[j]) / regyr
        for i in 1:n
            sp = Int(t.species[i]); d0 = t.dbh[i]
            d0 >= CI_RG_XMAX[sp] && continue
            t.tpa[i] <= 0.0f0 && continue
            (sp < 1 || sp > 19) && continue
            h1 = wk3[i]
            pt = Int(t.plot_id[i])
            ptba = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : ba
            pct = t.crown_ratio[i]
            relht = avh > 0f0 ? h1 / avh : 1.5f0; relht > 1.5f0 && (relht = 1.5f0)
            ptbali = ptba * (1.0f0 - pct / 100.0f0)
            tbal = (1.0f0 - pct / 100.0f0) * ba
            iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9); iicr < 1 && (iicr = 1)
            rcr = Float32(iicr)
            htgrl = CI_RG_CNST[sp] + CI_RG_HTBA[sp]*relht*ptba + CI_RG_PBAL[sp]*ptbali +
                    CI_RG_STBA[sp]*ba + CI_RG_RLHT[sp]*relht + CI_RG_CRSQ[sp]*rcr*rcr +
                    CI_RG_CR[sp]*rcr + CI_RG_BAL[sp]*tbal + CI_RG_HDM1[sp]*rhdm1 +
                    CI_RG_HDM2[sp]*rhdm2 + CI_RG_PTBA[sp]*ptba
            xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
            con = exp(c.htg_cor_small[sp])                 # RHCON≈1
            h2 = h1 + htgrl * scale * xrhgro * con
            h2 < h1 && (h2 = h1)
            wk3[i] = h2
        end
    end
    # HTDBH: final DBH from accumulated height; DG = new − old. Apply to diam_growth / ht_growth.
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d0 = t.dbh[i]
        d0 >= CI_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        (sp < 1 || sp > 19) && continue
        hk = wk3[i]; hk <= 4.5f0 && continue
        iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9); iicr < 1 && (iicr = 1)
        rcr = Float32(iicr)
        dk = exp(CI_RG_DHCN[sp] + CI_RG_DHHT[sp] * log(hk) + CI_RG_DHCR[sp] * log(rcr))
        dk < CI_RG_XMIN[sp] && (dk = CI_RG_XMIN[sp])
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        dg = (dk - d0) * xrdgro; dg < 0.0f0 && (dg = 0.0f0)
        t.diam_growth[i] = dg
        t.ht_growth[i] = wk3[i] - t.height[i]
    end
    return s
end
