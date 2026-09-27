# CRATET "ESTIMATE MISSING TOTAL TREE AGES" (every western <v>/cratet.f, just before CALL DGDRIV):
#
#       DO I=1,ITRN
#       IF(ABIRTH(I) .LE. 0.)THEN
#         ... D1=DBH(I); H=HT(I); D2=0.0
#         CALL FINDAG(I,ISPC,D1,D2,H,SITAGE,SITHT,AGMAX,HTMAX,HTMAX2,DEBUG)
#         IF(SITAGE .GT. 0.)ABIRTH(I)=SITAGE
#
# The variants below never read ABIRTH in their own growth, so this dub feeds only Climate-FVS: clgmult.f
# BIRTHYR=THISYR−ABIRTH (the Leites transfer distance behind TREEMULT) and clmorts.f:170 (the DE* dClimMort).
# Without it every inventory tree had BIRTHYR=THISYR ⇒ zero transfer distance ⇒ GrowthMult/dClimMort wrong.
# UT/TT/CR/BM/IE/EM dub in their own setup (their growth reads ABIRTH); KT's findag.f leaves SITAGE=0 and AK's
# cratet only clamps ABIRTH<0 to 0, so neither dubs. Each FINDAG is that variant's own findag.f (the jl ports
# the height-growth models already use), called with D1=DBH, D2=0, H=HT and the SITEAR of the species.
cratet_findag_dub!(s::StandState) = cratet_findag_dub!(s, s.variant)
cratet_findag_dub!(s::StandState, ::AbstractVariant) = s

@inline _sheppard_age(h::Float32) = fpow(h * 2.54f0 * 12.0f0 / 26.9825f0, 1.0f0 / 1.1752f0)   # aspen FINDAG

function _findag_dub_each!(f, s::StandState)
    t, p = s.trees, s.plot
    @inbounds for i in 1:t.n
        t.birth_age[i] > 0f0 && continue
        sp = Int(t.species[i])
        sitage = Float32(f(sp, t.dbh[i], t.height[i], p.sp_site_index[sp]))
        sitage > 0f0 && (t.birth_age[i] = sitage)
    end
    return s
end

# ci/findag.f: AS (13) Sheppard; MC (15) iterates the Alexander-form curve (AGMAX1=50, HTMAX1=20; INCRNG is set
# once before label 75, not per pass); every other species SITAGE=0 (no dub).
function _ci_findag(sp::Int, h::Float32, sindx::Float32)::Float32
    sp == 13 && return _sheppard_age(h)
    sp == 15 || return 0f0
    agmax1 = 50f0; htmax1 = 20f0; toler = 2f0
    h >= htmax1 && return agmax1 + (h - htmax1) / 0.10f0
    ag = 2f0; incrng = 0; hguess = 0f0
    while true
        oldhg = hguess
        a14 = fpow(ag, -1.4f0)
        hguess = (sindx - 4.5f0) / (0.6192f0 - 5.3394f0 / (sindx - 4.5f0) + 240.29f0 * a14 +
                                    (3368.9f0 / (sindx - 4.5f0)) * a14)
        hguess = hguess + 4.5f0
        if hguess >= 1f0
            (abs(hguess - h) <= toler || h < hguess) && return ag
            diff = hguess - oldhg
            (oldhg != 0f0 && diff >= 0.05f0) && (incrng = 1)
            (incrng == 1 && diff < 0.05f0) && return ag
        end
        ag += 2f0
        ag > agmax1 && return agmax1
    end
end

# canada/on/cratet.f FINDAG (canada/on/findag.f → HTCALC MODE0=0, HTMAX−1.1 retry) with SITEAR(ISPC): ON's volont.f
# Mowraski net-merch cull reads the dubbed ABIRTH (volume.jl on_mowraski), aged +FINT per cycle after UPDATE.
cratet_findag_dub!(s::StandState, ::Ontario) =
    _findag_dub_each!((sp, d, h, si) -> on_tree_age(sp, h, si), s)
cratet_findag_dub!(s::StandState, ::CentralIdaho) =
    _findag_dub_each!((sp, d, h, si) -> _ci_findag(sp, h, si), s)
cratet_findag_dub!(s::StandState, ::Klamath) =
    _findag_dub_each!((sp, d, h, si) -> nc_findag(h, sp, si)[1], s)
cratet_findag_dub!(s::StandState, ::PacificNorthwest) =
    _findag_dub_each!((sp, d, h, si) -> pn_findag(sp, d, 0f0, h, si)[1], s)
cratet_findag_dub!(s::StandState, ::WestCascades) =
    _findag_dub_each!((sp, d, h, si) -> wc_findag(sp, d, 0f0, h, si)[1], s)
cratet_findag_dub!(s::StandState, ::Olympic) =
    _findag_dub_each!((sp, d, h, si) -> op_findag(sp, d, 0f0, h, si)[1], s)
cratet_findag_dub!(s::StandState, ::EastCascades) =
    _findag_dub_each!((sp, d, h, si) -> ec_findag(sp, d, 0f0, h, si)[1], s)
function cratet_findag_dub!(s::StandState, ::SouthCentralOregon)
    ifor = Int(s.plot.forest_idx)
    _findag_dub_each!((sp, d, h, si) -> so_findag(sp, ifor, h, si)[1], s)
end
function cratet_findag_dub!(s::StandState, ::CentralCalifornia)
    ifor = Int(s.plot.forest_idx)
    _findag_dub_each!((sp, d, h, si) -> ca_findag(sp, h, si, ifor)[1], s)
end
cratet_findag_dub!(s::StandState, ::OregonCoast) =
    _findag_dub_each!((sp, d, h, si) -> oc_findag(sp, si, h)[1], s)
# ws/findag.f GB (21) reads BAU(ICLS)/BA; WS never calls BADIST, so the GGCOM BAU stays 0 ⇒ BAUTBA=0.
function cratet_findag_dub!(s::StandState, ::WestSierra)
    ifor = Int(s.plot.forest_idx)
    _findag_dub_each!((sp, d, h, si) -> ws_findag(sp, h, si, ifor, d, 0f0)[1], s)
end
# bc/cratet.f dubs only SELECT CASE(ISPC) 11,12,13,15, through the IE findag.f (Sheppard for every species).
cratet_findag_dub!(s::StandState, ::BritishColumbia) =
    _findag_dub_each!((sp, d, h, si) -> (sp == 11 || sp == 12 || sp == 13 || sp == 15) ? _sheppard_age(h) : 0f0, s)
