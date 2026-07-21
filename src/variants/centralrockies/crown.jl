# =============================================================================
# crown.jl (centralrockies) — CR crown ratio = GENGYM (cr/crown.f + gemcr.f)
#
# gemcr.f predicts crown LENGTH CL (linear in HF/DF/BAT/BAU) and converts to a
# ratio CR=CL/HF, except MH/RC/WL (CASE 6/7/8) which predict CR directly (ICRFLG=1).
# crown.f is the DUB/change wrapper (crown_ratio! hook).
#
# Math: ALOG->flog; RELDEN**2 (int) -> RELDEN*RELDEN.
# =============================================================================

"""
    cr_gemcr(imodty, is, bau, bat, hf, df, h, relden, pcti) -> CR

GENGYM crown ratio (cr/gemcr.f). Returns crown ratio (0-1). `HF`=future height (H+HTG),
`DF`=future dbh, `BAT`=stand BA, `BAU`=BA-above-class, `PCTI`=BA percentile.
"""
function cr_gemcr(imodty::Int, is::Int, bau::Float32, bat::Float32, hf::Float32, df::Float32,
                  h::Float32, relden::Float32, pcti::Float32)::Float32
    icrflg = 0
    cl = 0.0f0; cr = 0.0f0
    if is == 1 || is == 2
        cl = imodty <= 2 ? 0.50706f0 + 0.73070f0 * hf : 0.36135f0 + 0.57085f0 * hf
    elseif is == 3
        cl = 6.47479f0 + 0.54482f0 * df + 0.50703f0 * hf - 0.03326f0 * bat
    elseif is == 4 || is == 5
        cl = 6.22959f0 + 0.67587f0 * hf - 0.03098f0 * bat
    elseif is == 6
        cr = 0.3450f0 - 0.00264f0 * bat + 0.00000512f0 * relden * relden -
             0.25138f0 * flog(h) + 0.05140f0 * flog(pcti)
        icrflg = 1
    elseif is == 7
        cr = -1.6053f0 + 0.17479f0 * flog(bat) - 0.00183f0 * relden -
             0.00560f0 * df + 0.11050f0 * flog(pcti)
        icrflg = 1
    elseif is == 8
        cr = 0.03441f0 - 0.00204f0 * bat + 0.30066f0 * flog(df) - 0.59302f0 * flog(h)
        icrflg = 1
    elseif is == 11
        cl = 5.00215f0 + 0.06334f0 * hf + 0.88236f0 * df - 0.03821f0 * bau
    elseif is == 9 || is == 10 || is == 12 || is == 14 || is == 16 ||
           (23 <= is <= 27) || (29 <= is <= 35) || is == 37
        cl = -0.59373f0 + 0.67703f0 * hf
    elseif is == 13 || is == 36
        if imodty == 1
            cl = 5.63367f0 + 0.56252f0 * hf - 0.06411f0 * bat
        elseif imodty == 2
            cl = 4.35671f0 + 0.84714f0 * df + 0.32549f0 * hf - 0.03802f0 * bat
        else
            cl = 3.49178f0 + 0.80767f0 * df + 0.17421f0 * hf - 0.03272f0 * bat
        end
    elseif is == 15
        cl = 3.03832f0 + 0.65587f0 * hf - 0.01792f0 * bat
    elseif is == 17
        cl = 3.61635f0 + 0.93639f0 * df + 0.61547f0 * hf - 0.02360f0 * bat
    elseif is == 18
        cl = imodty <= 2 ? 1.05857f0 + 0.68442f0 * hf : 3.22244f0 + 0.44315f0 * hf + 0.44755f0 * df
    elseif is == 19
        cl = 0.15768f0 + 0.74697f0 * hf
    elseif is == 20 || is == 21 || is == 22 || is == 28 || is == 38
        cl = 5.17281f0 + 0.32552f0 * hf - 0.01675f0 * bat
    else
        cl = 1.0f0
    end
    if icrflg == 0
        cl < 1.0f0 && (cl = 1.0f0)
        cl > hf && (cl = hf)
        cr = cl / hf
    end
    return cr
end
