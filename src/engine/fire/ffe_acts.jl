# ffe_acts.jl — OPFIND/OPGET/OPDONE semantics for the FMIN activities (opfind.f, opmerg.f, opexpn.f, opcycl.f).
#
# FVS schedules every dated FFE keyword as an activity (OPNEW) and each FFE routine asks OPFIND for the activities of ITS code
# in the CURRENT cycle — NOT the ones whose date equals the year it is called in (the `IF (JYR .NE. IYR)` tests in fmburn.f are
# commented out). A date ≤ MAXCYC is a CYCLE number converted to that cycle's start year (opexpn.f:40-44); a calendar date is
# bucketed into the cycle that contains it (opcycl.f:58-64); a date before the first cycle belongs to cycle 1. Activities are
# returned in date order, ties in keyword order (OPSORT), and a done one (OPDONE) is not returned again.

"True when the activity date `date` (calendar year, or a cycle number < 1000) falls in the current cycle of `s`."
function act_due(s::StandState, date::Integer)::Bool
    c = s.control; cyc = Int(c.cycle)
    cs = cycle_year_at(c, cyc); ce = cycle_year_at(c, cyc + 1); ce <= cs && (ce = cs + 1)
    d = Int(date)
    d <= 0 && return true                                      # a date ≤ 0 is expanded to every cycle (opexpn.f:50-80)
    (0 < d < 1000) && return d == cyc + 1                      # cycle-number form (IY1.LE.MAXCYC)
    return (cs <= d < ce) || (d < cs && cyc == 0)
end

"The activity's date as a calendar year (a cycle-number date → that cycle's start year)."
act_year(s::StandState, date::Integer)::Int =
    (0 < Int(date) < 1000) ? cycle_year_at(s.control, Int(date) - 1) : Int(date)

"OPNEW for an FFE keyword: schedule activity `code` at `date` with parameter vector `params`."
function ffe_schedule!(s::StandState, code::Integer, date::Integer, params::AbstractVector{<:Real})
    u = s.fire.user
    push!(u.acts, FFEAct(Int32(code), Int32(date), Int32(length(u.acts) + 1), Int32(-1), Float32.(params)))
    return nothing
end

"OPFIND(1, code, NTODO) + the OPGET loop: this cycle's not-done activities of `code`, in date order (ties: keyword order)."
function ffe_find(s::StandState, code::Integer)::Vector{FFEAct}
    fs = s.fire; out = FFEAct[]
    (fs === nothing) && return out
    cyc = Int32(s.control.cycle)
    for a in fs.user.acts
        (a.code == code && a.dcyc != cyc && (a.dcyc < 0 || a.date <= 0) && act_due(s, a.date)) && push!(out, a)
    end
    sort!(out; by = a -> (act_year(s, a.date), a.seq))
    return out
end

"OPDONE: mark `a` done for this cycle."
ffe_done!(s::StandState, a::FFEAct) = (a.dcyc = Int32(s.control.cycle); nothing)

# ---------------------------------------------------------------------------------------------------------------------------
# FMUSRFM (fmusrfm.f): FIRECALC (2549), FUELMODL (2538) and FUELTRET (2525) activities of the current cycle → the persistent
# FMCOM state (IFLOGIC/…, LUSRFM + FMDUSR/FWTUSR, IFTYR/FTREAT/HARTYP/DPMOD). FMUSRFM runs at every FMCFMD call; every cycle's
# call comes through `ffe_cycle_acts!`, so a keyword whose cycle has no fire still takes effect (and the activity is never
# looked at again once its cycle is past).

# fmusrfm.f:60-94 — the model numbers FMUSRFM refuses (no such standard model): LUSRFM stays FALSE
_usrfm_bad(m::Int) = (31 <= m <= 100) || (110 <= m <= 120) || (125 <= m <= 140) || (150 <= m <= 160) ||
                     (166 <= m <= 180) || (190 <= m <= 200) || (205 <= m <= 256)
const _DPMULT = ((1.0f0, 1.3f0, 1.6f0), (0.83f0, 0.83f0, 0.83f0), (0.75f0, 0.75f0, 0.75f0))   # DPMULT(HARTYP, FTREAT+1) by FTREAT

function _fmusrfm_update!(s::StandState)
    fs = s.fire; fs === nothing && return
    u = fs.user
    for a in ffe_find(s, 2549)                                   # FIRECALC (fmusrfm.f:37-52): the last one of the cycle stands
        p = a.params
        u.iflogic = Int32(trunc(p[1])); u.ifmset = Int32(trunc(p[2]))
        u.usav = (p[3], p[4], p[5]); u.ubd = (p[6], p[7]); u.ulhv = p[8]
        ffe_done!(s, a)
    end
    for a in ffe_find(s, 2538)                                   # FUELMODL
        p = a.params; np = length(p)
        if np <= 0 || p[1] == 0f0
            u.lusrfm = false; ffe_done!(s, a)
        elseif any(p[i] > 256f0 for i in (1, 3, 5, 7) if i <= np) ||
               any(_usrfm_bad(Int(trunc(p[i]))) for i in (1, 3, 5, 7) if i <= np)
            u.lusrfm = false                                     # refused: not marked done, LUSRFM off
        else
            ffe_done!(s, a); u.lusrfm = true
            q = copy(p)
            sm = 0f0
            for i in 2:2:np; q[i] <= 0f0 && (q[i] = 1f0); sm += q[i]; end
            sm = 1f0 / sm
            for i in 2:2:np; q[i] = q[i] * sm; end
            fill!(u.fmdusr, Int32(0)); fill!(u.fwtusr, 0f0)
            j = 0
            for i in 1:2:np-1
                j += 1; j > 4 && break
                u.fmdusr[j] = Int32(trunc(q[i])); u.fwtusr[j] = q[i + 1]
            end
        end
    end
    for a in ffe_find(s, 2525)                                   # FUELTRET (fmusrfm.f:187-214)
        p = a.params
        lok = !(p[1] < 0f0 || p[1] > 2f0 || p[2] < 1f0 || p[2] > 3f0)
        if lok
            u.ftreat = Int32(trunc(p[1])); u.hartyp = Int32(trunc(p[2]))
            u.dpmod = p[3]
            u.dpmod < 0f0 && (u.dpmod = _DPMULT[Int(u.ftreat) + 1][Int(u.hartyp)])
            u.iftyr = Int32(current_cycle_year(s))
        end
        ffe_done!(s, a)                                          # OPDONE, or OPDEL1 for an invalid card
    end
    if u.iftyr > 0
        if Int(current_cycle_year(s)) - Int(u.iftyr) > 5
            u.iftyr = Int32(0); u.ftreat = Int32(0); u.hartyp = Int32(0); u.dpmod = 1f0
        end
    else
        u.ftreat = Int32(0); u.hartyp = Int32(0); u.dpmod = 1f0
    end
    return
end

# FMBURN's per-cycle activity handling that runs whether or not a fire is scheduled (fmburn.f:88-249): DROUGHT (2529) → IDRYB/IDRYE
# (the FIRST one of the cycle only, fmburn.f:123 GOTO 406) and DEFULMOD (2539) → the SURFVL/FMLOAD/FMDEP/MOISEX table entries.
function _fmburn_acts!(s::StandState)
    fs = s.fire; fs === nothing && return
    u = fs.user
    for a in ffe_find(s, 2529)
        d = act_year(s, a.date)
        u.idryb = Int32(d); u.idrye = Int32(trunc(Float32(d) + a.params[1] - 1f0))
        ffe_done!(s, a)
        break
    end
    for a in ffe_find(s, 2539)
        p = a.params                                             # PRMS(1:13): model, SAV(3), live-woody SAV, loads, depth, mext, herb SAV/load
        ifmd = Int(trunc(p[1]))
        ffe_done!(s, a)
        (1 <= ifmd <= size(s.coef.ffe_fuel_models, 1)) || continue
        load0, sav0, depth0, mext0 = fuel_model_resolved(s, ifmd)
        xsur = zeros(Float32, 2, 3); xfml = zeros(Float32, 2, 3)
        for i in 2:4; xsur[1, i - 1] = p[i] < 0f0 ? sav0[1, i - 1] : p[i]; end
        xsur[2, 1] = p[5] < 0f0 ? sav0[2, 1] : p[5]
        xsur[2, 2] = p[12] < 0f0 ? sav0[2, 2] : p[12]
        sum(xsur) <= 0.0001f0 && continue                        # LOK = FALSE ⇒ OPDEL1
        for i in 6:8; xfml[1, i - 5] = p[i] < 0f0 ? load0[1, i - 5] : p[i]; end
        xfml[2, 1] = p[9] < 0f0 ? load0[2, 1] : p[9]
        xfml[2, 2] = p[13] < 0f0 ? load0[2, 2] : p[13]
        indd = count(>(0f0), xfml[1, :]); inl = count(>(0f0), xfml[2, :])
        (indd <= 0 && inl <= 0) && continue
        xdep = p[10] < 0f0 ? depth0 : p[10]; xdep < 0f0 && continue
        xext = p[11] < 0f0 ? mext0 : p[11]; (xext < 0f0 || xext > 1f0) && continue
        load = copy(load0); sav = copy(sav0)
        for j in 1:3; sav[1, j] = Float32(trunc(Int, xsur[1, j])); load[1, j] = xfml[1, j]; end    # SURFVL is INTEGER (fmburn.f:225)
        sav[2, 1] = Float32(trunc(Int, xsur[2, 1])); sav[2, 2] = Float32(trunc(Int, xsur[2, 2]))
        load[2, 1] = xfml[2, 1]; load[2, 2] = xfml[2, 2]
        fs.defulmod[Int32(ifmd)] = (load, sav, xdep, xext)
    end
    return
end

"Per-cycle FFE keyword activities (FMBURN/FMUSRFM) — call once per cycle before the fire."
function ffe_cycle_acts!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return
    _fmburn_acts!(s)
    _fmusrfm_update!(s)
    return
end

# ---------------------------------------------------------------------------------------------------------------------------
# Post-activity fuel models for the eastern builds that carry the AFWT block (NE, LS; CS/SN have none).
# FMFMOV (fmtret.f:372-400, label 506 — after the FUELMOVE moves, every FMMAIN): OLARGE/OSMALL = the previous SMALL/LARGE, the new
# SMALL/LARGE over BOTH pile pools, SLCHNG = 100·(LARGE+SMALL−X)/X (X = OLARGE+OSMALL, 0 when X ≤ 1E-6). fmcfmd.f then opens the <=5-yr
# post-harvest window (FM11) when SLCHNG ≥ SLCRIT (10, fmvinit.f:1252) or LATFUEL is set. (west-tail implements the same state for the
# western builds as fs.slchng/_ffe_afwt!; keep ONE implementation when the two branches meet.)
function fmfmov_slchng!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active || !fs.fuels_init) && return s
    u = fs.user; yr = Int32(current_cycle_year(s))
    u.fmtret_year == yr && return s
    u.fmtret_year = yr
    sm, lg = _small_large_fuel(fs)
    x = u.o_large + u.o_small
    u.slchng = x > 1f-6 ? 100f0 * (lg + sm - x) / x : 0f0
    u.o_large = lg; u.o_small = sm
    return s
end

"fmcfmd.f AFWT/LATFUEL block (NE: ne/fmcfmd.f:136-142; LS shares it for the JP/MW/OA cover types): returns AFWT, sets EQWT(11)."
function _east_afwt!(s::StandState, eqwt::AbstractVector{Float32})::Float32
    u = s.fire.user
    afwt = max(0f0, 1f0 - Float32(Int(current_cycle_year(s)) - Int(s.fire.harvyr)) / 5f0)
    if u.slchng >= 10f0 || u.latfuel
        u.latfuel = true
        eqwt[11] = afwt
        afwt <= 0f0 && (u.latfuel = false)
    end
    u.latfuel || (afwt = 0f0)
    return afwt
end

# ---------------------------------------------------------------------------------------------------------------------------
# FUELFOTO photo series (fmphotoval.f, byte-identical in every FVS build): FOTOVAL(1:11) hard + FOTOVALS(1:9) soft tons/ac per (reference, code),
# -1 = not listed in the guide. Data generated from the oracle source (data/fire/photo_series.csv, test/harness/fmphoto/README.txt).
function _load_photo(path::AbstractString)
    d = Dict{Tuple{Int,Int},Tuple{NTuple{11,Float32},NTuple{9,Float32}}}()
    for line in eachline(path)
        (isempty(line) || line[1] == '#' || startswith(line, "ref,")) && continue
        f = split(strip(line), ',')
        length(f) >= 22 || continue
        v = Float32[parse(Float32, f[i]) for i in 3:22]
        d[(parse(Int, f[1]), parse(Int, f[2]))] = (ntuple(i -> v[i], 11), ntuple(i -> v[11 + i], 9))
    end
    d
end
const _PHOTO = _load_photo(joinpath(@__DIR__, "..", "..", "..", "data", "fire", "photo_series.csv"))
