# =============================================================================
# species_translation.jl — unknown-species translation (SPCTRN)
#
# Ported from: base/spctrn.f  (the 562-row ASPT crosswalk table).
#
# When a tree's species code isn't one of the variant's own codes, FVS maps it to
# the nearest variant species via this crosswalk: columns are
#   (alpha, FIA, PLANTS, CS-target, LS-target, NE-target, SN-target).
# The 562-row table now lives in data/southern/species_translation.csv and is loaded
# into `coef.translation::Vector{NTuple{7,String}}`. The variant picks its own target
# column (spctrn_column) and "other" catch-all species (other_species), so the table
# is shared and the engine stays variant-agnostic (requirement #10).
# =============================================================================

"Column of the translation table holding this variant's target code (1=alpha…7=SN)."
function spctrn_column end

"This variant's catch-all 'other' species index."
function other_species end

spctrn_column(::Southern) = 7
other_species(::Southern) = Int32(90)        # OT (other)

"""
    translate_species(code, variant, sp, coef) -> (index::Int32, format::Int32)

Map an unrecognized species `code` to a variant species index via the SPCTRN
crosswalk in `coef.translation`. `format` is 1=alpha, 2=FIA, 3=PLANTS (the column
the code matched). Pure. Ported from spctrn.f.
"""
function translate_species(code::AbstractString, variant::AbstractVariant,
                           sp::SpeciesData, coef::SpeciesCoefficients)
    c   = rstrip(code)
    col = spctrn_column(variant)
    target = "XX"
    fmt = Int32(3)
    for row in coef.translation
        a = rstrip(row[1]); f = rstrip(row[2]); p = rstrip(row[3])
        if !isempty(a) && c == a
            fmt = Int32(1); target = row[col]; break
        elseif !isempty(f) && c == f
            fmt = Int32(2); target = row[col]; break
        elseif !isempty(p) && c == p
            fmt = Int32(3); target = row[col]; break
        end
    end
    rstrip(target) == "XX" && return (other_species(variant), fmt)
    # find target (first 2 chars) in the variant species class codes (spctrn.f:48)
    t2 = rstrip(target)
    for j in 1:MAXSP
        cc = sp.class_codes[j, 1]
        if length(cc) >= 2 && t2 == cc[1:2]
            return (Int32(j), fmt)
        end
    end
    return (other_species(variant), fmt)
end

"""
    resolve_species(code, variant, sp, coef) -> (index::Int32, format::Int32)

Resolve a raw tree-record species code to a species index: first try a direct
match against the variant's own alpha/FIA/PLANTS codes, then fall back to the
SPCTRN crosswalk. Ported from intree.f:240-263.
"""
function resolve_species(code::AbstractString, variant::AbstractVariant,
                         sp::SpeciesData, coef::SpeciesCoefficients)
    c = uppercase(strip(code))
    isempty(c) && (c = "OT")
    @inbounds for j in 1:MAXSP
        strip(sp.alpha[j])  == c && return (Int32(j), Int32(1))
        strip(sp.fia[j])    == c && return (Int32(j), Int32(2))
        strip(sp.plants[j]) == c && return (Int32(j), Int32(3))
    end
    return translate_species(c, variant, sp, coef)
end
