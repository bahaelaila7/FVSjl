# R6 snag dynamics shared by the Region-6 FFE variants (BM/EC/SO/PN/OP/WC; tables in r6_snag_tables.jl):
#   FMR6SDCY (fmr6sdcy.f) — species group, snag size class SML, years-to-soft JYRSOFT and the decay adjustment JADJ;
#   FMR6FALL (fmr6fall.f) — the annual fall FRACTION BASE = SNFL(SPG, MOIS, OTSH, OTRT, JADJ, SML);
#   FMR6HTLS (fmr6htls.f) — the random snag height loss (a RANN draw per FMSNGHT call).
# Unlike the SN/CS FMSFALL (a CONSTANT modrate·ORIGDEN stems/yr with a last-5% ramp), the R6 form falls a fraction of
# the CURRENT density: DFALLN = BASE·FALLX·DENTTL (ec/fmsfall.f:44-46, == bm/wc/pn/op/ak; so/fmsfall.f Oregon branch).

"""
    r6_ffe_code(v, kodfor) -> Symbol

Which FMR6SDCY/FMR6FALL table set the variant reads (`:EC`, `:BM`, `:SO`, `:PN` (PN and OP share it), `:WC`), or
`:none` for a variant outside the R6 snag model. SO is R6 only on its Oregon forests; each SO routine has its own
test (fall: so/fmsfall.f KODFOR 5xx/701 = California; decay time: fmsngdk.f 601/602/620/799 = Oregon; height loss:
fmsnght.f 505/506/509/511/701/514 = California), so this returns `:SO` and the callers apply their own test.
"""
@inline function r6_ffe_code(v)::Symbol
    v isa EastCascades && return :EC
    v isa BlueMountains && return :BM
    v isa SouthCentralOregon && return :SO
    (v isa PacificNorthwest || v isa Olympic) && return :PN
    v isa WestCascades && return :WC
    v isa SoutheastAlaska && return :AK   # vbase fmr6sdcy/fmr6fall CASE('AK'); FVSak's fmsfall.f/fmsngdk.f take the R6 path
    return :none
end

@inline _r6_clamp(i::Integer, n::Integer) = clamp(Int(i), 1, n)

"""
    r6_sdcy(code, ksp, dbh, itype) -> (yrsoft, jadj, sml)

FMR6SDCY (fmr6sdcy.f:356-431): the snag size class SML from the species' DBH breakpoints (cm; the 2nd is 50 when the
1st is 20, else 75 — fmr6sdcy.f:163-196), the habitat temperature/moisture class, and from the YRSOFT/DCYADJ tables
(PN/OP/AK → PN*, WC → WC*, else ES*) the years-to-soft and the fall-adjustment class.
"""
function r6_sdcy(code::Symbol, ksp::Integer, dbh::Float32, itype::Integer)
    dbhcm = dbh * 2.54f0
    if code === :EC
        it = _r6_clamp(itype, 155); temp = _R6SD_ECHMC[it]; mois = _R6SD_ECWMD[it]
        spg = _R6SD_ECSPEC[ksp]; d1 = _R6SD_ECDBH1[ksp]
    elseif code === :BM
        it = _r6_clamp(itype, 92); temp = _R6SD_BMHMC[it]; mois = _R6SD_BMWMD[it]
        spg = _R6SD_BMSPEC[ksp]; d1 = _R6SD_BMDBH1[ksp]
    elseif code === :SO
        it = _r6_clamp(itype, 92); temp = _R6SD_SOHMC[it]; mois = _R6SD_SOWMD[it]
        spg = _R6SD_SOSPEC[ksp]; d1 = _R6SD_SODBH1[ksp]
    elseif code === :PN
        it = _r6_clamp(itype, 75); temp = _R6SD_PNWMC[it]; mois = _R6SD_PNWMD[it]
        spg = _R6SD_WSSPEC[ksp]; d1 = _R6SD_WSDBH1[ksp]
    elseif code === :AK
        temp = 3; mois = 1                                   # fmr6sdcy.f 'AK': assume cold, assume wet
        spg = _R6SD_AKSPEC[ksp]; d1 = _R6SD_AKDBH1[ksp]
    else  # :WC
        it = _r6_clamp(itype, 139); temp = _R6SD_WCWMC[it]; mois = _R6SD_WCWMD[it]
        spg = _R6SD_WSSPEC[ksp]; d1 = _R6SD_WSDBH1[ksp]
    end
    d1f = Float32(d1); d2f = d1 == 20 ? 50f0 : 75f0
    sml = dbhcm < d1f ? 1 : (dbhcm < d2f ? 2 : 3)
    if code === :PN || code === :AK                          # fmr6sdcy.f CASE('PN','AK','OP')
        return (Int(_R6SD_PNYRSOFT[spg, temp, mois, sml]), Int(_R6SD_PNDCYADJ[spg, temp, mois, sml]), sml)
    elseif code === :WC
        return (Int(_R6SD_WCYRSOFT[spg, temp, mois, sml]), Int(_R6SD_WCDCYADJ[spg, temp, mois, sml]), sml)
    end
    return (Int(_R6SD_ESYRSOFT[spg, temp, mois, sml]), Int(_R6SD_ESDCYADJ[spg, temp, mois, sml]), sml)
end

# fmr6fall.f:799-806 LOREG: these forests are in WASHINGTON (species-group column J=2); every other KODFOR is Oregon.
@inline _r6_washington(kodfor::Integer) =
    kodfor in (609, 603, 605, 608, 617, 699, 613, 701) || kodfor > 999 || kodfor == 800

"""
    r6_fall_base(code, ksp, dbh, itype, kodfor) -> Float32

FMR6FALL BASE (fmr6fall.f:817-850): SNFL(SPG, MOIS, OTSH, OTRT, JADJ, SML) with JADJ/SML from `r6_sdcy`, the species
group from the variant's (Oregon|Washington) column, and the habitat moisture / overstory-shade / root-rot classes.
"""
function r6_fall_base(code::Symbol, ksp::Integer, dbh::Float32, itype::Integer, kodfor::Integer)::Float32
    _, jadj, sml = r6_sdcy(code, ksp, dbh, itype)
    j = _r6_washington(kodfor) ? 2 : 1
    if code === :EC
        it = _r6_clamp(itype, 155)
        spg = _R6FL_ECSPEC[ksp, j]; mois = _R6FL_ECWMD[it]; otsh = _R6FL_ECOTSH[it]; otrt = _R6FL_ECOTRT[it]
    elseif code === :BM
        it = _r6_clamp(itype, 92)
        spg = _R6FL_BMSPEC[ksp]; mois = _R6FL_BMWMD[it]; otsh = _R6FL_BMOTSH[it]; otrt = _R6FL_BMOTRT[it]
    elseif code === :SO
        it = _r6_clamp(itype, 92)
        spg = _R6FL_SOSPEC[ksp]; mois = _R6FL_SOWMD[it]; otsh = _R6FL_SOOTSH[it]; otrt = _R6FL_SOOTRT[it]
    elseif code === :PN
        it = _r6_clamp(itype, 75)
        spg = _R6FL_WSSPEC[ksp, j]; mois = _R6FL_PNWMD[it]; otsh = _R6FL_PNOTSH[it]; otrt = _R6FL_PNOTRT[it]
    elseif code === :AK                                       # fmr6fall.f CASE('AK'): MOIS=OTSH=OTRT=1
        spg = _R6FL_AKSPEC[ksp]; mois = 1; otsh = 1; otrt = 1
    else  # :WC
        it = _r6_clamp(itype, 139)
        spg = _R6FL_WSSPEC[ksp, j]; mois = _R6FL_WCWMD[it]; otsh = _R6FL_WCOTSH[it]; otrt = _R6FL_WCOTRT[it]
    end
    return _R6FL_SNFL[spg, mois, otsh, otrt, jadj, sml]
end

"BM snag fall BASE (kept for callers/tests): FMR6FALL CASE('BM') — BM has one species-group column, so no KODFOR."
bm_r6_fall_base(ksp::Integer, dbh::Float32, itype::Integer)::Float32 = r6_fall_base(:BM, ksp, dbh, itype, 0)

"""
    r6_htls(code, ksp, y) -> Float32

FMR6HTLS (fmr6htls.f): the fraction of height a snag loses this year given the uniform draw `y` = RANN — SNHTLS(SPG)
when y ≤ PRHTLS(SPG), else 0. Species group: ECSPEC / BMSPEC / SOSPEC, the westside WSSPEC for every other variant.
"""
@inline function r6_htls(code::Symbol, ksp::Integer, y::Float32)::Float32
    spg = code === :EC ? _R6HL_ECSPEC[ksp] : code === :BM ? _R6HL_BMSPEC[ksp] :
          code === :SO ? _R6HL_SOSPEC[ksp] : _R6HL_WSSPEC[ksp]
    return y <= _R6HL_PRHTLS[spg] ? _R6HL_SNHTLS[spg] : 0f0
end

"FMR6SDCY/FMR6FALL habitat index ITYPE: BM decodes its own (bm_itype); SO defaults to CPS111 = 49 (so/habtyp.f); EC/WC/PN/OP carry it in habitat_input."
function _r6_itype(s::StandState)::Int
    v = s.variant
    v isa BlueMountains && return bm_itype(Int(s.plot.habitat_code))
    hi = Int(s.plot.habitat_input)
    v isa SouthCentralOregon && hi <= 0 && return 49
    return hi
end

# so/fmsfall.f:48-49 — the California forests (KODFOR 5xx, 701) keep the western 18" linear fall; Oregon is R6.
@inline _so_california_fall(kodfor::Integer) = (500 <= kodfor < 600) || kodfor == 701
# so/fmsngdk.f — only these Oregon forests use the R6 JYRSOFT decay time; every other SO forest the default formula.
@inline _so_oregon_dk(kodfor::Integer) = kodfor in (601, 602, 620, 799)
# so/fmsnght.f — these California forests use the HTX height-loss formula; every other SO forest FMR6HTLS.
@inline _so_california_ht(kodfor::Integer) = kodfor in (505, 506, 509, 511, 701, 514)

"""
Variants whose FMSFALL is the western 18"-linear form ({v}/fmsfall.f, md5 fd611f… = IE/KT/NC/WS/CA/OC; CI/TT/UT/CR/EM
variants of it): the IE/KT/CI/TT/UT/EM/CR layer plus NC/WS/CA/OC and SO's California forests (R6 handled first).
"""
_ffe_west_fallform(v) = _ffe_west_vol(v) || v isa Klamath || v isa WestSierra || v isa CentralCalifornia ||
                        v isa OregonCoast || v isa SouthCentralOregon

"""
    _snag_dktime(s, sp, d, dcx) -> Float32

FMSNGDK (fmsnag.f:279-285) years-since-death for a snag to turn soft. PN/WC/BM/EC/OP (and SO's Oregon forests
601/602/620/799): DKTIME = JYRSOFT·DECAYX from FMR6SDCY. Everything else: (1.24·DECAYX·D) + (13.82·DECAYX), in that
exact Float32 order (fmsngdk.f:80 CASE DEFAULT; XMOD = 1).
"""
function _snag_dktime(s::StandState, sp::Int, d::Float32, dcx::Float32)::Float32
    r6 = r6_ffe_code(s.variant)
    (r6 === :SO && !_so_oregon_dk(Int(s.plot.user_forest_code))) && (r6 = :none)
    if r6 !== :none
        yrsoft, _, _ = r6_sdcy(r6, sp, d, _r6_itype(s))
        return Float32(yrsoft) * dcx
    end
    return (1.24f0 * dcx * d) + (13.82f0 * dcx)
end
