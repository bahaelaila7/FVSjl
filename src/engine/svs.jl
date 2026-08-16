# =============================================================================
# svs.jl — SVS (Stand Visualization System) data path — CHUNK 0
#
# Ported from:
#   base/svkey.f    (SVKEY)    — the SVS keyword flag + _index.svs open
#   vbase/svstart.f (SVSTART)  — the cycle-0 inventory-picture seam (fvs.f:333)
#   base/svgtpl.f   (SVGTPL)   — deterministic subplot layout (IPLGEM=0 branch)
#   base/svestb.f   (SVESTB)   — TPA→objects placement (integer path)
#   base/svgtpt.f   (SVGTPT)   — random point in a rectangle (IPLGEM<2 branch)
#   base/svobol.f   (SVOBOL)   + base/svcrol.f (SVCROL) — circle/circle overlap
#   base/svrann.f   (SVRANN)   — the 2nd Park–Miller stream (divisor 2^31) — in rng.jl
#   base/cwcalc.f   (CWCALC)   — KT western forest-grown crown width (Crookston R1/R6)
#   vbase/svout.f   (SVOUT)    — header (:243-258) + live-tree object record (:429-440)
#
# SCOPE (chunk 0): cycle-0 inventory picture only; IPLGEM=0 (one 208.71-ft square acre,
# subplot ids ignored); INTEGER TPA (FRACAD≡0 ⇒ the fractional lottery is skipped, so the
# ONLY svrann draws are SVGTPT's two-per-object); live green trees (IOBJTP=1). Western
# #TREEFORM = WEST.TRF. Later chunks add: fractional lottery + overlap retries (1),
# multi-subplot layout (2), later cycles (3), cut/mortality removal (4), snags/CWD (5).
#
# Validated bit-exact vs the live relinked FVSkt on stand S248112 (kt0.key/kt0.tre):
# _001.svs is byte-identical (see test/integration/test_svs_chunk0.jl).
# =============================================================================

const _SVS_SIDE_IMPERIAL = 208.7103f0   # sqrt(43560) — square-acre side, feet (svgtpl.f:46)
const _SVS_DIAI_FACTOR   = 0.04166667f0 # DBH(in) → stem radius(ft): 1/12 × 0.5 (svestb.f:308)

# KT crown-width equation code per KT species sequence number (cwcalc.f KTMAP, ISPC 1..11):
#   1 WP 2 WL 3 DF 4 GF 5 WH 6 RC 7 LP 8 ES 9 AF 10 PP 11 OT/MH
const _KT_CWEQN = ("11903","07303","20203","01703","26303","24203",
                   "10803","09303","01903","12203","26405")

"""
    kt_crown_width(ispc, D, H, CR, BA, EL) -> Float32

Forest-grown crown width (ft) for KT species sequence `ispc`, DBH `D`, height `H`, crown
ratio percent `CR`, stand basal area `BA`, elevation `EL` (100s ft). Mirrors CWIDTH→CWCALC
(cwcalc.f) for the Kootenai variant: KT is FS region 1 (KODFOR<601) ⇒ the R6 forest factor
BF=1 and the R6-forest adjustment block is skipped. `CL = CR·H·0.01` is crown length;
`BAREA = max(BA, 1)`. All arithmetic in Float32 to match the Fortran REAL single precision.
"""
function kt_crown_width(ispc::Integer, D::Float32, H::Float32, CR::Float32,
                        BA::Float32, EL::Float32)::Float32
    cl    = CR * H * 0.01f0
    barea = BA <= 1f0 ? 1f0 : BA
    eq    = _KT_CWEQN[ispc]
    cw    = 0f0
    if eq == "11903"           # WP — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.0405f0 * exp(1.2799f0 + 0.11941f0*log(cl) + 0.42745f0*log(b) - 0.07182f0*log(barea))
        D < 1f0 && (cw *= D / 1f0); cw > 35f0 && (cw = 35f0)
    elseif eq == "07303"       # WL — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.02478f0 * exp(0.99889f0 + 0.19422f0*log(cl) + 0.59423f0*log(b) -
                             0.09078f0*log(H) - 0.02341f0*log(barea))
        D < 1f0 && (cw *= D / 1f0); cw > 40f0 && (cw = 40f0)
    elseif eq == "20203"       # DF — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.01685f0 * exp(1.48372f0 + 0.27378f0*log(cl) + 0.49646f0*log(b) -
                             0.18669f0*log(H) - 0.01509f0*log(barea))
        D < 1f0 && (cw *= D / 1f0); cw > 80f0 && (cw = 80f0)
    elseif eq == "01703"       # GF — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.0303f0 * exp(1.14079f0 + 0.20904f0*log(cl) + 0.38787f0*log(b))
        D < 1f0 && (cw *= D / 1f0); cw > 40f0 && (cw = 40f0)
    elseif eq == "26303"       # WH — Crookston (R1)
        b = D >= 0.1f0 ? D : 0.1f0
        cw = 1.02460f0 * exp(1.3522f0 + 0.24844f0*log(cl) + 0.412117f0*log(b) -
                             0.104357f0*log(H) + 0.03538f0*log(barea))
        D < 0.1f0 && (cw *= D / 0.1f0); cw > 54f0 && (cw = 54f0)
    elseif eq == "24203"       # RC — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.03597f0 * exp(1.46111f0 + 0.26289f0*log(cl) + 0.18779f0*log(b))
        D < 1f0 && (cw *= D / 1f0); cw > 45f0 && (cw = 45f0)
    elseif eq == "10803"       # LP — Crookston (R1)
        b = D >= 0.7f0 ? D : 0.7f0
        cw = 1.03992f0 * exp(1.58777f0 + 0.30812f0*log(cl) + 0.64934f0*log(b) - 0.38964f0*log(H))
        D < 0.7f0 && (cw *= D / 0.7f0); cw > 40f0 && (cw = 40f0)
    elseif eq == "09303"       # ES — Crookston (R1)
        b = D >= 0.1f0 ? D : 0.1f0
        cw = 1.02687f0 * exp(1.28027f0 + 0.2249f0*log(cl) + 0.47075f0*log(b) - 0.15911f0*log(H))
        D < 0.1f0 && (cw *= D / 0.1f0); cw > 40f0 && (cw = 40f0)
    elseif eq == "01903"       # AF — Crookston (R1)
        b = D >= 0.1f0 ? D : 0.1f0
        cw = 1.02886f0 * exp(1.01255f0 + 0.30374f0*log(cl) + 0.37093f0*log(b) - 0.13731f0*log(H))
        D < 0.1f0 && (cw *= D / 0.1f0); cw > 30f0 && (cw = 30f0)
    elseif eq == "12203"       # PP — Crookston (R1)
        b = D >= 2f0 ? D : 2f0
        cw = 1.02687f0 * exp(1.49085f0 + 0.1862f0*log(cl) + 0.68272f0*log(b) - 0.28242f0*log(H))
        D < 2f0 && (cw *= D / 2f0); cw > 46f0 && (cw = 46f0)
    elseif eq == "26405"       # OT/MH — Crookston (R6) model 2 (BF=1 for KT region 1)
        el = EL < 10f0 ? 10f0 : (EL > 79f0 ? 79f0 : EL)
        b  = D >= 1f0 ? D : 1f0
        cw = 3.7854f0 * (b^0.54684f0) * (H^(-0.12954f0)) * (cl^0.16151f0) *
             ((barea + 1f0)^0.03047f0) * (exp(el)^(-0.00561f0))
        D < 1f0 && (cw *= D / 1f0); cw > 45f0 && (cw = 45f0)
    end
    return cw
end

# --- SVOBOL/SVCROL: two stem circles overlap iff (r1+r2)² > centre-distance² (svcrol.f:15-22) ---
@inline function _svs_circles_overlap(x1::Float32, y1::Float32, r1::Float32,
                                      x2::Float32, y2::Float32, r2::Float32)::Bool
    d2 = (x2 - x1)^2 + (y2 - y1)^2
    rr = r1 + r2
    return rr * rr > d2
end

# --- SVGTPT rectangle branch (IPLGEM<2): two svrann draws map to (x,y) in [x1,x2]×[y1,y2] ---
@inline function _svgtpt_rect(rng::FVSRng, x1::Float32, x2::Float32, y1::Float32, y2::Float32)
    xr = svrann!(rng)
    x  = x1 + (x2 - x1) * xr
    yr = svrann!(rng)
    y  = y1 + (y2 - y1) * yr
    return (x, y)
end

# --- SVS object (SVDATA.F77 /SVOBJ/): one visualization stem. IS2F is the tree-record pointer
#     (remapped by SVTRIP on tripling / SVRMOV on removal); (x,y)=XSLOC/YSLOC persist across cycles.
mutable struct SVSObj
    species::Int          # ISP snapshot (fixed; the remapped record is always the same species)
    is2f::Int             # IS2F: tree-record index this object represents (0 ⇒ removed)
    x::Float32            # XSLOC(ISVOBJ)
    y::Float32            # YSLOC(ISVOBJ)
end

"""
    svs_place_objects(s) -> Vector{SVSObj}

SVSTART inventory placement (cycle-0): SVGTPL(IPLGEM=0)+SVESTB(integer)+SVGTPT, advancing the SVS
random stream (seeded here from the main stream at the SVSTART seam). Returns the persistent SVS
object list (IOBJTP=1 live green trees) in placement order — the array whose (x,y) survives every
later cycle while IS2F is remapped by tripling/removal. Chunk assumptions (asserted): IPLGEM=0,
integer TPA. `s` must be at the inventory state (post `setup_growth!`/`compute_volumes!`, pre-growth).
"""
function svs_place_objects(s::StandState)::Vector{SVSObj}
    t = s.trees
    @assert Int(s.control.svs_iplgem) == 0 "svs: only IPLGEM=0 supported"

    # --- SVSTART seed: SVS stream starts at the MAIN stream's current s0 (svstart.f:50) ---
    svs_seed!(s.rng)

    side = _SVS_SIDE_IMPERIAL          # imperial (IMETRIC=0); one square acre per subplot

    objs = SVSObj[]
    # --- SVESTB (integer path): NTOADD(i) = ifix(PROB(i)); FRACAD≡0 ⇒ no lottery START draw ---
    for i in 1:t.n
        prob = t.tpa[i]
        ntoadd = trunc(Int, prob)                         # IFIX(PROB) with NOBTS=0
        # integer-TPA guard. A fractional remainder means the lottery path (chunk 1) is needed.
        @assert abs(prob - ntoadd) < 1f-4 "svs: non-integer TPA on record $i (PROB=$prob) needs the lottery (chunk 1)"
        diai = t.dbh[i] * _SVS_DIAI_FACTOR                # stem radius (ft)
        itre = t.plot_id[i]
        for _ in 1:ntoadd
            x = 0f0; y = 0f0
            ncks = 0
            while true
                (x, y) = _svgtpt_rect(s.rng, 0f0, side, 0f0, side)
                # SVOBOL overlap reject vs already-placed live objects on the SAME subplot (svestb.f:323-348).
                overlap = false
                for o in objs
                    t.plot_id[o.is2f] == itre || continue
                    rj = t.dbh[o.is2f] * _SVS_DIAI_FACTOR
                    if _svs_circles_overlap(x, y, diai, o.x, o.y, rj)
                        overlap = true
                        break
                    end
                end
                if overlap
                    ncks += 1
                    ncks > 40 && break                    # give up after 40 tries; place anyway (svestb.f:343)
                    continue
                end
                break
            end
            push!(objs, SVSObj(Int(t.species[i]), i, x, y))
        end
    end
    return objs
end

"""
    svs_svtrip!(objs, nlive)

SVTRIP (svtrip.f): update the object→record pointers for record tripling. FVS TRIPLE splits every
live base record `i` (1..`nlive`) into (central=i, upper=nlive+2i-1, lower=nlive+2i); this mirrors
svtrip.f, distributing the NR objects that pointed at base `i` as ⌊0.6·NR+0.5⌋ to central, then
⌊0.625·rem+0.5⌋ to upper, the rest to lower (identical to triple_records!'s .60/.25/.15 TPA split).
(x,y) are untouched. No-op for a base with ≤1 object (svtrip.f:44). Faithful to the physical
`nlive+2i-1 / nlive+2i` append order in triple_records! (diameter_growth.jl).
"""
function svs_svtrip!(objs::Vector{SVSObj}, nlive::Integer)
    for i in 1:nlive
        itr = (i, nlive + 2i - 1, nlive + 2i)             # (central, upper, lower) — TRIPLE layout
        # count NR = live objects currently pointing at base record i (svtrip.f:40-42)
        nr = 0
        for o in objs
            o.is2f == i && (nr += 1)
        end
        nr <= 1 && continue
        nrt1 = trunc(Int, nr * 0.6f0 + 0.5f0)             # IFIX(FLOAT(NR)*.6+.5)
        nrt2 = max(0, trunc(Int, (nr - nrt1) * 0.625f0 + 0.5f0))
        nrt3 = max(0, nr - nrt1 - nrt2)
        nrt2 <= 0 && continue                             # svtrip.f:52
        nrt = (nrt1, nrt2, nrt3)
        slot = 1; ntoc = nrt1
        for o in objs
            o.is2f == i || continue
            o.is2f = itr[slot]
            ntoc -= 1
            if ntoc == 0
                slot += 1
                slot > 3 && break
                ntoc = nrt[slot]
                ntoc == 0 && break
            end
        end
    end
    return objs
end

"""
    svs_picture(objs, s; year, msg) -> String

Format one SVS picture (SVOUT header svout.f:243-258 + live-tree object records svout.f:429-440)
from the persistent object list `objs` and the CURRENT tree state `s.trees`. Each object is emitted
with its (remapped) record's DBH/HT/ICR/crown-width and its persistent (x,y). Removed objects
(IS2F=0) are skipped (svout.f:391 `I.GT.0`).
"""
function svs_picture(objs::Vector{SVSObj}, s::StandState; year::Integer,
                     msg::AbstractString)::String
    t = s.trees; p = s.plot
    ba = p.basal_area; el = p.elevation
    io = IOBuffer()
    stand = strip(String(p.stand_id))
    print(io, "#TITLE Stand=", stand, " Year=", _svs_i4(year), " ", msg, "\n")
    print(io, "#TREEFORM WEST.TRF\n")     # western cluster
    print(io, "#FORMAT 2\n")
    print(io, "#PLOTSIZE 208.71 208.71\n")
    print(io, "#UNITS ENGLISH\n")
    print(io, ";                  trcl  stus             fang\n")
    print(io, ";species        tr#  |crcl|   dbh   ht lang |edia crd  cr    crd  cr    crd  cr    crd  cr ex mk  xloc    yloc  z\n")
    for o in objs
        o.is2f == 0 && continue           # svout.f:391 I=IS2F(ISVOBJ); IF (I.GT.0)
        rec  = o.is2f
        sp2  = rpad(rstrip(String(s.species.code2[o.species])), 2)   # SPCD: 2-char, left-justified
        icr  = abs(t.crown_pct[rec])
        xicr = Float32(icr) * 0.01f0
        cw   = kt_crown_width(o.species, t.dbh[rec], t.height[rec],
                              Float32(t.crown_pct[rec]), ba, el)      # CW=CRWDTH(I)
        crad = cw / 2f0
        _svs_write_tree!(io, sp2, rec, t.dbh[rec], t.height[rec], crad, xicr, o.x, o.y)
    end
    return String(take!(io))
end

"""
    svs_render_cycle0(s; msg="Inventory conditions") -> String

Cycle-0 inventory picture (SVSTART seam): place objects then format. Kept for the chunk-0 test.
"""
function svs_render_cycle0(s::StandState; msg::AbstractString = "Inventory conditions")::String
    objs = svs_place_objects(s)
    return svs_picture(objs, s; year = current_cycle_year(s), msg = msg)
end

# --- SVOUT format 30 for IOBJTP=1 (svout.f:429-431/439):
#     (A,T16,I5,I3,2I2,F6.1,F6.0,I2,I4,I2,4(F6.1,1X,F4.2),2I2,2F8.2,I2)
#     A=SPCD  I5=tree#  I3=class(0)  I2 I2=(0,IPS=1)  F6.1=DBH  F6.0=HT  I2=lean(0)  I4=dir(0)
#     I2=edia(0)  4×(F6.1 crownRad, 1X, F4.2 crownRatio)  I2 I2=(ex=1,mk=0)  F8.2 F8.2=xloc,yloc  I2=z(0)
function _svs_write_tree!(io::IO, sp2::AbstractString, rec::Integer, dbh::Float32, ht::Float32,
                          crad::Float32, xicr::Float32, x::Float32, y::Float32)
    print(io, rpad(sp2, 15))                       # A + T16 (SPCD in cols 1-2, next field at col 16)
    print(io, _svs_i(rec, 5))                      # I5 tree#
    print(io, _svs_i(0, 3))                        # I3 tree class
    print(io, _svs_i(0, 2), _svs_i(1, 2))          # crown class, plant status (IPS=1)
    print(io, _svs_f(dbh, 6, 1))                   # F6.1 dbh
    print(io, _svs_f0(ht, 6))                      # F6.0 ht (trailing '.')
    print(io, _svs_i(0, 2), _svs_i(0, 4), _svs_i(0, 2))   # lean, felling dir, small-end dia
    for _ in 1:4
        print(io, _svs_f(crad, 6, 1), " ", _svs_f(xicr, 4, 2))   # 4×(crown radius, crown ratio)
    end
    print(io, _svs_i(1, 2), _svs_i(0, 2))          # expansion factor (1), marking status (0)
    print(io, _svs_f(x, 8, 2), _svs_f(y, 8, 2))    # xloc, yloc
    print(io, _svs_i(0, 2))                        # z
    print(io, "\n")
end

# --- minimal Fortran edit-descriptor emitters (all fixed-width, blank-filled) ---
_svs_i(v::Integer, w::Integer) = lpad(string(v), w)
_svs_i4(v::Integer) = lpad(string(v), 4, '0')                       # I4.4
_svs_f(v::Real, w::Integer, d::Integer) =
    lpad(_fmt_fixed(Float64(v), d), w)
"Fortran F<w>.0 — round to integer, print WITH a trailing decimal point (e.g. 73.0 → \"   73.\")."
_svs_f0(v::Real, w::Integer) = lpad(string(round(Int, Float64(v))) * ".", w)

"""
    svs_write_cycle0_files(stem, s; imageno=1, msg="Inventory conditions") -> String

The SVSTART seam (fvs.f:333, cycle-0 inventory picture): write the two SVS files FVS emits with
the default JSVPIC=91 — `<stem>_NNN.svs` (the picture: SVOUT header + per-tree object records) and
`<stem>_index.svs` (`#TREELISTINDEX` + one `"Stand=… Year=…" "<pic>"` index line, svout.f:200/211).
Returns the picture-file path. Gated by the caller on `control.svs_on` (SVS keyword seen).
"""
function svs_write_cycle0_files(stem::AbstractString, s::StandState;
                                imageno::Integer = 1, msg::AbstractString = "Inventory conditions")
    body    = svs_render_cycle0(s; msg = msg)
    picfile = string(stem, "_", lpad(string(imageno), 3, '0'), ".svs")
    open(picfile, "w") do io; write(io, body); end
    stand = strip(String(s.plot.stand_id))
    year  = current_cycle_year(s)
    open(string(stem, "_index.svs"), "w") do io
        print(io, "#TREELISTINDEX\n")
        print(io, "\"Stand=", stand, " Year=", _svs_i4(year), " ", msg, "\" \"", basename(picfile), "\"\n")
    end
    return picfile
end

"""
    svs_project!(stem, s; fint=10f0) -> Vector{String}

Full multi-cycle SVS data path for the standard projection: mirror the engine's SVS seams —
SVSTART (fvs.f:333, inventory picture), GRINCR (grincr.f:277, `IF ICYC>1` "Beginning of cycle"
picture) each later cycle, and MAIN (fvs.f:453, "End of projection") — driving `grow_cycle!`
between them and remapping the persistent object list through SVTRIP on every tripling cycle.

Writes `<stem>_NNN.svs` for each picture and the accumulated `<stem>_index.svs` (`#TREELISTINDEX`
+ one line per picture). Returns the picture-file paths. `s` must be at the inventory state
(post `setup_growth!`/`compute_volumes!`, pre-growth). Cycle count comes from the keyword schedule
(`s.control`); tripling fires only for the first ICL4 cycles, exactly as the engine decides it.
"""
function svs_project!(stem::AbstractString, s::StandState; fint::Float32 = 10f0)
    t = s.trees
    stand = strip(String(s.plot.stand_id))
    ncyc  = Int(s.control.ncycle)

    # --- SVSTART (cyc0 inventory picture) ---
    objs   = svs_place_objects(s)
    index  = Tuple{Int,String,String}[]     # (year, msg, picfile-basename)
    pics   = String[]
    imageno = 1
    picfile = string(stem, "_", lpad(string(imageno), 3, '0'), ".svs")
    open(picfile, "w") do io; write(io, svs_picture(objs, s; year = current_cycle_year(s),
                                                     msg = "Inventory conditions")); end
    push!(pics, picfile)
    push!(index, (current_cycle_year(s), "Inventory conditions", basename(picfile)))

    # --- projection loop (mirrors fvs.f: for ICYC=1..NUMCYCLE) ---
    for icyc in 1:ncyc
        # GRINCR seam (grincr.f:277): IF ICYC>1 emit "Beginning of cycle" BEFORE growing this cycle.
        if icyc > 1
            imageno += 1
            picfile = string(stem, "_", lpad(string(imageno), 3, '0'), ".svs")
            open(picfile, "w") do io; write(io, svs_picture(objs, s; year = current_cycle_year(s),
                                                             msg = "Beginning of cycle")); end
            push!(pics, picfile)
            push!(index, (current_cycle_year(s), "Beginning of cycle", basename(picfile)))
        end
        # grow one cycle; TRIPLE (if it fired) split every live base record → remap the objects.
        nlive_pre = t.n
        grow_cycle!(s; fint = fint)
        compute_volumes!(s)
        t.n > nlive_pre && svs_svtrip!(objs, nlive_pre)   # tripling appended records ⇒ SVTRIP
    end

    # --- MAIN seam (fvs.f:453): "End of projection" picture ---
    imageno += 1
    picfile = string(stem, "_", lpad(string(imageno), 3, '0'), ".svs")
    open(picfile, "w") do io; write(io, svs_picture(objs, s; year = current_cycle_year(s),
                                                     msg = "End of projection")); end
    push!(pics, picfile)
    push!(index, (current_cycle_year(s), "End of projection", basename(picfile)))

    # --- accumulated _index.svs ---
    open(string(stem, "_index.svs"), "w") do io
        print(io, "#TREELISTINDEX\n")
        for (yr, msg, pf) in index
            print(io, "\"Stand=", stand, " Year=", _svs_i4(yr), " ", msg, "\" \"", pf, "\"\n")
        end
    end
    return pics
end

"Format a Float64 with exactly `d` fractional digits (round-half-away, like Fortran F edit)."
function _fmt_fixed(v::Float64, d::Integer)::String
    neg = v < 0
    a = abs(v)
    scale = 10.0^d
    n = floor(Int, a * scale + 0.5)                # round half up (matches Fortran F for our data)
    ip = n ÷ Int(scale)
    fp = n % Int(scale)
    s = d == 0 ? string(ip) : string(ip, ".", lpad(string(fp), d, '0'))
    return neg ? "-" * s : s
end
