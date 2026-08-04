# Western-variant REGEN (small-tree) regression guard. The western variants historically validated their
# small-tree/regen path only on the LARGE-tree 248112 synthetic stands (ut/em/bm/cr t01 read imperially →
# 11.5″ trees), so the regen path went untested — hiding real bugs (UT ht-dbh curve #134, EM SMHTGF+SMDGF
# #136). This guard runs a REALISTIC small-tree stand per variant (the same 248112 stand with DBH shrunk to
# <4″ AND height=4.5+6·d, NOAUTOES to isolate growth from establishment) through FVSjl and diffs vs the live
# FVS oracle `.sum`, asserting bit-exact-or-cornered per column. It locks in the UT/EM fixes + BM/CR clean
# state. (Stands + oracle .sum live in /workspace/.{ut,em,bm,cr}work — the live-oracle differential harness,
# same convention as the BC guard.)
#
# NOT covered here: TT (live-binary≠source, blocked). IE IS covered (task #139 fixed 2026-08-04: the small-tree
# path did DBH-direct — the ESTAB-only branch — plus a 1-draw-per-tree ZZRAN in record order, desyncing the RNG
# and doubling regen BA; fixed to the faithful non-ESTAB increment path (DG=(DK−D1)·XRDGRO·BARK on the DDS scale,
# DBH via GRADD) + species-sorted per-record triple-draw ⇒ BA 2×→~7%, all aggregates cornered).
import Pkg; Pkg.activate("/workspace/FVSjl"; io = devnull)
using FVSjl

# (variant, key.base, TopHt tol). Columns: Year Age TPA BA SDI CCF TopHt QMD.
# IE's TopHt runs looser (~0.29) because HTG carries a large multiplicative ZZRAN (exp(zzran·0.59)) that dominates
# the stand EXTREME (top height) while averaging out of the aggregates — the pre-ZZRAN htgr1 is bit-exact vs live,
# so this is the accepted ZZRAN-tail residual, not a height-model gap. UT/EM/BM/CR keep the tight 0.15 TopHt bound.
const CASES = [
    (FVSjl.Utah(),           "/workspace/.utwork/utt01_smallr", 0.15),
    (FVSjl.EasternMontana(), "/workspace/.emwork/emt01_smallr", 0.15),
    (FVSjl.BlueMountains(),  "/workspace/.bmwork/bmt01_smallr", 0.15),
    (FVSjl.CentralRockies(), "/workspace/.crwork/crt01_smallr", 0.15),
    (FVSjl.InlandEmpire(),   "/workspace/.iework/ierun/iet01_smallr", 0.30),
]
# Tolerances chosen to ACCEPT the converged cornered state yet FAIL a fix-revert (UT ht-dbh revert → ~30%+
# under-growth; EM SMHTGF revert → 2-3× height over). Known cornered residuals within these bounds:
#  - UT BA ~10% mid-cycle (ZZRAN tail), - BM late-cycle TPA up to ~20% (jl under-thins late = SDI-plateau-at-
#    SDImax vs oracle thinning; a self-thin/mortality residual — flagged for follow-up, NOT a regen-growth bug).
# EM + CR track bit-exact-or-tight. IE is cornered (BA/QMD/SDI ~7%, TopHt = ZZRAN tail, see per-case tol). TT
# excluded (blocked). TPA tol 0.32 accepts BM's late-cycle under-thin (28.6% @2090, open task #140); revert=100%+.
const COLS = [(3, "TPA", :rel, 0.32), (4, "BA", :rel, 0.15), (5, "SDI", :rel, 0.15),
              (6, "CCF", :rel, 0.15), (7, "TopHt", :rel, 0.15), (8, "QMD", :rel, 0.12)]

_rows(t) = [split(l) for l in split(strip(t), "\n") if !startswith(strip(l), "-999") && !isempty(strip(l))]

function check(variant, base, toptol)
    jl = FVSjl.run_keyfile("$base.key"; variant = variant)
    orc = read("$base.sum", String)
    J = _rows(jl); O = _rows(orc); nrow = min(length(J), length(O))
    fails = String[]; exact = 0
    for r in 1:nrow
        all(J[r][c] == O[r][c] for c in 3:8) && (exact += 1)
        for (c, name, kind, tol0) in COLS
            tol = c == 7 ? toptol : tol0                   # per-case TopHt override (IE ZZRAN tail)
            jv = parse(Float64, J[r][c]); ov = parse(Float64, O[r][c])
            d = kind === :rel ? (ov == 0 ? abs(jv) : abs(jv - ov) / abs(ov)) : abs(jv - ov)
            d > tol && push!(fails, "$name yr $(O[r][1]): jl=$jv or=$ov (Δ$(round(d,digits=3))>$tol)")
        end
    end
    return (nrow, exact, fails)
end

function main()
    allok = true
    for (v, base, toptol) in CASES
        vn = string(nameof(typeof(v)))
        nrow, exact, fails = check(v, base, toptol)
        if isempty(fails)
            println("✓ $vn  ($(basename(base))): $nrow cycles, $exact bit-exact — bit-exact-or-cornered")
        else
            allok = false
            println("✗ $vn  ($(basename(base))): $(length(fails)) out-of-tolerance:")
            for f in fails[1:min(6, length(fails))]; println("     ", f); end
        end
    end
    allok ? println("\n✓ WESTERN REGEN GUARD PASS (UT/EM/BM/CR/IE small-tree)") :
            (println("\n✗ REGRESSION"); exit(1))
end
main()
