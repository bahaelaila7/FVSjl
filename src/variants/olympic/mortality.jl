# =============================================================================
# mortality.jl (olympic) — OP mortality (op/morts.f).
#
# op/morts.f is NOT a cooperating split: whenever ORGANON ran this cycle (SMORMT = ΣMORTEXP > 0),
# the ORGANON per-record mortality MORTEXP is used for EVERY tree record — IORG=0 and IORG=1 alike
# (op/morts.f:13-19,331-337): `WKI = MORTEXP(I)·(FINT/5)`, capped at PROB. Only when NO valid ORGANON
# tree exists (SMORMT=0) does op fall back to its FVS-native RIP mortality (the per-species BM0..BM5
# logistic + Gould–Harrington small-tree model). For S248112 (has DF big-6) ORGANON always runs, so
# the MORTEXP-for-all path is the one exercised and validated. The native RIP path errors loudly until
# a no-big-6 OP stand needs it (doctrine #5).
#
# The op/morts.f SDIMAX<5 kill-all and the iterative SDI-cap over-density pass (morts.f:508-547) are
# inert for a stand ORGANON already thinned below SDImax; they are omitted here (PASS=1) and will be
# ported if a dense stand shows a residual. Deterministic (DGSD=0). FINT=5 on OP ⇒ FINT/5 = 1.
# =============================================================================

function mortality!(s::StandState, ::Olympic; fint::Float32 = 5.0f0, book_snags::Bool = true)
    t, c = s.trees, s.calib
    n = t.n; n == 0 && return s
    killed = zeros(Float32, n)
    if c.op_org_ran && length(c.op_mortexp) == n
        # ORGANON MORTEXP for ALL records (op/morts.f:331-336): WKI = MORTEXP·(FINT/5), capped at PROB.
        fscale = fint / 5f0
        @inbounds for i in 1:n
            p = t.tpa[i]; p <= 0f0 && continue
            wki = c.op_mortexp[i] * fscale
            wki > p && (wki = p)
            killed[i] = wki
        end
    else
        error("OP mortality: FVS-native RIP path (no valid ORGANON tree, SMORMT=0) is UNPORTED. " *
              "S248112 has a DF big-6 so ORGANON always runs; port op/morts.f BM0..BM5 + Gould–Harrington " *
              "small-tree logistic before running a no-big-6 OP stand.")
    end
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
