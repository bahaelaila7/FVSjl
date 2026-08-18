# =============================================================================
# ontario/small_tree_growth.jl — ON small-tree REGENT (canada/on/regent.f).
#
# CHUNK 8a. `small_tree_growth!(s, stash, ::Ontario; fint)` mirrors canada/on/regent.f.
#
# ON's REGENT (Penner 2010) applies the small-tree HEIGHT + DBH increment ONLY to records with
#   DBH < XMAX(ISPC).  For the whole 72-species table XMIN/XMAX are the CONSTANTS 3.15"/4.72"
#   (= 8 cm / 12 cm; regent.f:78 `DATA XMIN /MAXSP*3.15/, XMAX/MAXSP*4.72/`).  A record with
#   D >= XMAX branches to `GO TO 25` (regent.f:139) — the large-tree DG/HTG already computed by
#   `diameter_growth!`/`height_growth!` is kept unchanged, and no tripled small-tree record is
#   touched (the triple loop sits AFTER the D>=XMX skip).  So for an all-large-tree stand REGENT
#   is a genuine no-op.
#
# The small-tree branch itself (D < XMAX) is the NC128 height-age curve (htcalc.f) blended with
# the large-tree HTG over [XMIN,XMAX], plus the Wykoff HT->DBH dubbing (htont.f ONHTDBH / htdbh.f)
# and DGBND size-cap — ~760 lines of Fortran with NO ON validation stand exercising it (the only
# oracle stand, ont01, is 8 large trees whose internal DBH stays >= 5.9" the entire projection).
# Per doctrine we do NOT stage that unvalidated math: a record that is actually small errors
# loudly rather than silently producing wrong numbers.  When an ON small-tree/regen validation
# stand exists this becomes a full port (htcalc + ONHTDBH + htdbh + DGBND) — see the handoff.
#
# VALIDATED: on ont01 every live record has D >= XMX every cycle, so this reproduces REGENT's
# no-op exactly (the tree list — dbh/diam_growth/ht_growth and the tripling stash — is unchanged
# across the call, matching FVSon_g16 which likewise skips every record).
# =============================================================================

const ON_REG_XMAX = 4.72f0    # regent.f:78 XMAX = 12 cm (uniform across all 72 species)
const ON_REG_XMIN = 3.15f0    # regent.f:78 XMIN = 8 cm  (uniform)

"""
    small_tree_growth!(s, stash, ::Ontario; fint) — ON REGENT (canada/on/regent.f).

Applies the small-tree height/DBH increment to records with DBH < XMAX (4.72" = 12 cm). For a
stand of only large trees (D >= XMAX for every record — the ont01 case) this is a validated no-op:
the large-tree DG/HTG stand, unchanged. A record below XMAX errors (the ON small-tree model is a
separate, un-validated chunk; no oracle stand exercises it).
"""
function small_tree_growth!(s::StandState, stash, ::Ontario; fint::Float32 = 10.0f0)
    t = s.trees
    n = t.n
    n == 0 && return s
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        d = t.dbh[i]
        d >= ON_REG_XMAX && continue          # regent.f:139 D>=XMX -> GO TO 25 (keep large-tree DG/HTG)
        error("ON small-tree REGENT (canada/on/regent.f htcalc/ONHTDBH/htdbh path) is not yet " *
              "ported: record $i has DBH=$(d)\" < XMAX=$(ON_REG_XMAX)\" (12 cm). No ON validation " *
              "stand exercises the small-tree model; port it (with a small-tree oracle stand) " *
              "before running a stand with sub-12 cm trees or active regeneration.")
    end
    return s
end
