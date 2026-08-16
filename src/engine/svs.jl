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

"""
    svs_render_cycle0(s; msg="Inventory conditions") -> String

Produce the byte content of the cycle-0 SVS picture file (`<stem>_001.svs`) for the current
(inventory) stand: run SVGTPL(IPLGEM=0)+SVESTB(integer)+SVGTPT placement — advancing the SVS
random stream (seeded here from the main stream at the SVSTART seam) — then write the SVOUT
header and the live-tree (IOBJTP=1) object records. Pure w.r.t. files; the caller writes it.

Requires `s` at the inventory state (post `setup_growth!`/`compute_volumes!`, pre-growth) so
DBH/HT/PCT(ICR)/PROB/ITRE/BA match the live SVSTART seam. Chunk-0 assumptions (asserted):
IPLGEM=0, integer TPA. Non-integer PROB or IPLGEM≠0 → later chunks.
"""
function svs_render_cycle0(s::StandState; msg::AbstractString = "Inventory conditions")::String
    t = s.trees
    p = s.plot
    @assert Int(s.control.svs_iplgem) == 0 "svs chunk 0: only IPLGEM=0 supported"

    # --- SVSTART seed: SVS stream starts at the MAIN stream's current s0 (svstart.f:50) ---
    svs_seed!(s.rng)

    # --- SVGTPL (IPLGEM=0): one square acre per subplot; subplot ids ignored for placement ---
    side = _SVS_SIDE_IMPERIAL          # imperial (IMETRIC=0)

    ba  = p.basal_area
    el  = p.elevation

    # object list: (species_index, record_index, x, y) in placement order
    objs = Tuple{Int,Int,Float32,Float32}[]

    # --- SVESTB (integer path): NTOADD(i) = ifix(PROB(i)); FRACAD≡0 ⇒ no lottery START draw ---
    for i in 1:t.n
        prob = t.tpa[i]
        ntoadd = trunc(Int, prob)                         # IFIX(PROB) with NOBTS=0
        # chunk-0 guard: integer TPA. A fractional remainder means the lottery path (chunk 1) is needed.
        @assert abs(prob - ntoadd) < 1f-4 "svs chunk 0: non-integer TPA on record $i (PROB=$prob) needs the lottery (chunk 1)"
        diai = t.dbh[i] * _SVS_DIAI_FACTOR                # stem radius (ft)
        itre = t.plot_id[i]
        for _ in 1:ntoadd
            x = 0f0; y = 0f0
            ncks = 0
            while true
                (x, y) = _svgtpt_rect(s.rng, 0f0, side, 0f0, side)
                # SVOBOL overlap reject vs already-placed live objects on the SAME subplot (svestb.f:323-348).
                # Chunk-0 stems are sub-foot radii on a 208-ft square ⇒ this never rejects (0 retries,
                # verified: every object consumes exactly 2 svrann draws). Ported faithfully for later reuse.
                overlap = false
                dobj = 0
                for (spj, recj, xj, yj) in objs
                    t.plot_id[recj] == itre || continue
                    rj = t.dbh[recj] * _SVS_DIAI_FACTOR
                    if _svs_circles_overlap(x, y, diai, xj, yj, rj)
                        overlap = true
                        break
                    end
                    dobj += 1
                end
                if overlap
                    ncks += 1
                    ncks > 40 && break                    # give up after 40 tries; place anyway (svestb.f:343)
                    continue
                end
                break
            end
            push!(objs, (Int(t.species[i]), i, x, y))
        end
    end

    # --- SVOUT: header (svout.f:243-258) + live-tree object records (svout.f:429-440) ---
    io = IOBuffer()
    stand = strip(String(p.stand_id))
    year  = current_cycle_year(s)
    print(io, "#TITLE Stand=", stand, " Year=", _svs_i4(year), " ", msg, "\n")
    print(io, "#TREEFORM WEST.TRF\n")     # western cluster
    print(io, "#FORMAT 2\n")
    print(io, "#PLOTSIZE 208.71 208.71\n")
    print(io, "#UNITS ENGLISH\n")
    print(io, ";                  trcl  stus             fang\n")
    print(io, ";species        tr#  |crcl|   dbh   ht lang |edia crd  cr    crd  cr    crd  cr    crd  cr ex mk  xloc    yloc  z\n")

    for (sp, rec, x, y) in objs
        sp2  = rpad(rstrip(String(s.species.code2[sp])), 2)   # SPCD: 2-char, left-justified
        icr  = abs(t.crown_pct[rec])
        xicr = Float32(icr) * 0.01f0
        cw   = kt_crown_width(sp, t.dbh[rec], t.height[rec], Float32(t.crown_pct[rec]), ba, el)
        crad = cw / 2f0
        _svs_write_tree!(io, sp2, rec, t.dbh[rec], t.height[rec], crad, xicr, x, y)
    end
    return String(take!(io))
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
