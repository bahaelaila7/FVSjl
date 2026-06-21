# habitat.jl — HABTYP/HBDECD: ecological-unit (PCOM) code decoder (Southern variant).
#
# Ported from: habtyp.f + hbdecd.f. The STDINFO habitat field is read BOTH as a
# number (ARRAY2) and as text (KARD2). When numeric it is a 1-based INDEX into the
# 320-entry Southern ecological-unit table (SNECU); when alphabetic it is matched
# (uppercased, trimmed) against that same table. The resolved table entry — not the
# raw field — is the PCOM string the diameter-growth physiography term keys off.
#
# Example: STDINFO field "232" is index 232 → "234AA" (p234 group), NOT the literal
# "232" (which would be p232). The default when nothing matches is entry 122,
# "231DD" (habtyp.f:6).

# 320 Southern ecological-unit codes (NPA=320), DATA-initialized (habtyp.f).
const SNECU = String[
    "221DB","221DD","221DE","221EB","221EG","221EJ","221EN","221HA","221HB","221HC",
    "221HD","221HE","221JA","221JB","221JC","222AB","222AG","222AH","222AL","222AM",
    "222AN","222CB","222CC","222CD","222CE","222CF","222CG","222CH","222DA","222DB",
    "222DC","222DD","222DE","222DG","222DI","222DJ","222EA","222EB","222EC","222ED",
    "222EE","222EF","222EG","222EH","222EI","222EJ","222EK","222EN","222EO","222FA",
    "222FB","222FC","222FD","222FF","223AB","223AG","223AH","223AM","223AN","223BB",
    "223BC","223BD","223DA","223DB","223DC","223DD","223DE","223DG","223DI","223DJ",
    "223EA","223EB","223EC","223ED","223EE","223EF","223EG","223EH","223FA","223FB",
    "223FC","223FD","223FF","231AA","231AB","231AC","231AD","231AE","231AF","231AG",
    "231AH","231AI","231AJ","231AK","231AL","231AM","231AN","231AO","231AP","231BA",
    "231BB","231BC","231BD","231BE","231BF","231BG","231BH","231BI","231BJ","231BK",
    "231BL","231CA","231CB","231CC","231CD","231CE","231CF","231CG","231DA","231DB",
    "231DC","231DD","231DE","231EA","231EB","231EC","231ED","231EE","231EF","231EG",
    "231EH","231EI","231EJ","231EK","231EL","231EM","231EN","231EO","231FA","231FB",
    "231GA","231GB","231GC","231HA","231HB","231HC","231HD","231HE","231HF","231HH",
    "231HI","231IA","231IB","231IC","231ID","231IE","231IF","231IG","232AD","232BA",
    "232BB","232BC","232BD","232BE","232BF","232BG","232BH","232BI","232BJ","232BK",
    "232BL","232BM","232BN","232BO","232BP","232BQ","232BR","232BS","232BT","232BU",
    "232BV","232BX","232BZ","232CA","232CB","232CC","232CD","232CE","232CF","232CG",
    "232CH","232CI","232CJ","232DA","232DB","232DC","232DD","232DE","232EA","232EB",
    "232EC","232ED","232EE","232EF","232FA","232FB","232FC","232FD","232FE","232FF",
    "232GA","232GB","232GC","232GD","232HA","232HB","232HC","232IA","232IB","232JA",
    "232JB","232JC","232JD","232JE","232JF","232JG","232KA","232KB","232LA","232LB",
    "232LC","234AA","234AB","234AC","234AD","234AE","234AF","234AG","234AH","234AI",
    "234AJ","234AK","234AL","234AM","234AN","234CA","234CB","234CC","234CD","234DA",
    "234DB","234DC","234DO","234EA","234EB","234EC","251EA","251EC","251ED","251FB",
    "251FC","255AA","255AB","255AC","255AD","255AE","255AF","255AG","255AH","255AI",
    "255AJ","255AK","255AM","255BA","255CA","255CC","255CD","255CE","255CF","255CG",
    "255CH","255DA","255DB","255DC","255DD","255EA","255EB","255EC","255ED","255EE",
    "411AA","411AB","411AC","411AD","411AE","411AF","411AG","M221AA","M221AB","M221AC",
    "M221BA","M221BD","M221BE","M221CA","M221CB","M221CC","M221CD","M221CE","M221DA","M221DB",
    "M221DC","M221DD","M222AA","M222AB","M223AA","M223AB","M231AA","M231AB","M231AC","M231AD",
]

"Default ecological-unit index when the habitat field matches nothing (habtyp.f:6)."
const ECU_DEFAULT_IDX = 122   # "231DD"

"""
    resolve_eco_unit(field, value) -> String

HABTYP/HBDECD: decode the STDINFO habitat field into its PCOM ecological-unit code.
`value` is the field's numeric parse (0 when non-numeric); `field` is the raw text.
A positive numeric value indexes SNECU directly; otherwise the uppercased, trimmed
text is matched against SNECU. Anything unresolved falls back to "231DD".
"""
function resolve_eco_unit(field::AbstractString, value::Real)
    npa = length(SNECU)
    ihb = floor(Int, value)
    if ihb >= 1 && ihb <= npa
        return SNECU[ihb]                       # numeric index into the table
    elseif ihb == 0
        s = strip(field)
        (isempty(s) || s[1] == '0') && return SNECU[ECU_DEFAULT_IDX]
        temp = first(uppercase(s), 8)
        idx = findfirst(==(temp), SNECU)
        return idx === nothing ? SNECU[ECU_DEFAULT_IDX] : SNECU[idx]
    end
    return SNECU[ECU_DEFAULT_IDX]               # out-of-range numeric
end
