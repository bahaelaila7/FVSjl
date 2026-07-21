# =============================================================================
# small_tree_growth.jl (centralrockies) — CR small-tree/regen growth (cr/regent.f)
#
# For trees with DBH < XMAX[sp]: a height increment from a site/density/vigor model
# (POTHTG·PCTRED·VIGOR·CON, or the Sheppard aspen/birch curve), weighted (XWT) with
# the large-tree htgf estimate; then a diameter increment from the height change via
# a species-specific inverse height–DBH curve, DDS-scaled to the cycle.
#
# Math: ALOG->flog, EXP->fexp, **realexp->fpow, SQRT->Base. Stochastic ZZRAN
# (DGSD≥1, ASYMMETRIC reject window [-2, 0.5]) injected by the caller.
# =============================================================================

"""
    _cr_regent_tree(...) -> (htg, dg)

Per-tree small-tree height + diameter increment (cr/regent.f growth path). `htg_large` = the htgf
large-tree height increment for the XWT blend (0 for pure seedlings). `pothtg`/`pctred`/`rsimod` are
stand/species terms from the caller. `zzran` is the injected stochastic deviate (0 ⇒ deterministic).
"""
function _cr_regent_tree(sp::Int, d::Float32, h::Float32, icr::Int, abirth::Float32,
                         rsimod::Float32, pothtg::Float32, pctred::Float32, con::Float32,
                         xrhgro::Float32, xrdgro::Float32, scale::Float32, scale2::Float32,
                         wk4::Float32, htg_large::Float32, sizcap4::Float32, sitear::Float32,
                         bark::Float32, ivflag::Bool, lskiph::Bool, zzran::Float32,
                         dgmax::Float32, break_sp::Float32, xmn::Float32, xmx::Float32,
                         diam_sp::Float32, ht1::Float32, ht2::Float32)
    dgmx = dgmax * scale
    # ---- HEIGHT increment ----
    local htg::Float32
    if lskiph
        htg = 0.0f0
    else
        if sp == 20 || sp == 28                            # aspen / paper birch (Sheppard curve)
            ag1 = abirth < 5.0f0 ? 5.0f0 : abirth
            hite1 = 26.9825f0 * fpow(ag1, 1.1752f0)
            ag2 = ag1 + 10.0f0
            hite2 = 26.9825f0 * fpow(ag2, 1.1752f0)
            htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 0.75f0
        else
            x = Float32(icr) / 100.0f0
            vigor = 150.0f0 * fpow(x, 3.0f0) * fexp(-6.0f0 * x) + 0.3f0
            vigor > 1.0f0 && (vigor = 1.0f0)
            ivflag && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)
            htgr = pothtg * pctred * vigor * con
        end
        htgr = (htgr + zzran * 0.2f0) * xrhgro * scale * wk4
        xwt = (d - xmn) / (xmx - xmn)
        d <= xmn && (xwt = 0.0f0)
        htg = htgr * (1.0f0 - xwt) + xwt * htg_large
        htg < 0.1f0 && (htg = 0.1f0)
        if h + htg > sizcap4
            htg = sizcap4 - h; htg < 0.1f0 && (htg = 0.1f0)
        end
    end
    # ---- DIAMETER increment (only for D < BREAK[sp]) ----
    dg = 0.0f0
    if d < break_sp
        hk = h + htg
        if hk <= 4.5f0
            dg = 0.0f0                                     # DBH set to D+.001*HK by caller if needed
        else
            local dk::Float32, dkk::Float32
            if sp == 13 || sp == 36                        # ponderosa / Chihuahua pine
                dk = (hk - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0; dk < 0.1f0 && (dk = 0.1f0)
                dkk = (h - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0; dkk < 0.1f0 && (dkk = 0.1f0)
                h < 4.5f0 && (dkk = d)
            elseif ivflag                                  # pinyon/juniper/oak/bristlecone
                dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                dkk = (h - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                h < 4.5f0 && (dkk = d)
            else                                           # all other species (IABFLG=1 ⇒ AX=HT1)
                bx = ht2; ax = ht1
                dk = (bx / (flog(hk - 4.5f0) - ax)) - 1.0f0; dk < 0.1f0 && (dk = 0.1f0)
                dkk = h <= 4.5f0 ? d : (bx / (flog(h - 4.5f0) - ax)) - 1.0f0
            end
            if dk < 0.0f0 || dkk < 0.0f0
                dg = htg * 0.2f0 * bark * xrdgro
                dk = d + dg
            else
                dg = (dk - dkk) * bark * xrdgro
            end
            dg < 0.0f0 && (dg = 0.0f0)
            dg > dgmx && (dg = dgmx)
            dds = dg * (2.0f0 * bark * d + dg) * scale2
            dg = sqrt(fpow(d * bark, 2.0f0) + dds) - bark * d
        end
        (d + dg) < diam_sp && (dg = diam_sp - d)           # DBH(K)≈D for cycling (no estab reassign)
    end
    return htg, dg
end
