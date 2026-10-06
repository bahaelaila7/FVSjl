# fmcfmd2.jl — FIRECALC (fmusrfm.f / fmcfmd2.f): IFLOGIC 1 = the "new fuel model logic" (FMCFMD2: Scott's rules pick the two closest of the 13
# Anderson / 40 Scott-Burgan models by a departure index), IFLOGIC 2 = modelled loads (FMCFMD3 fills the custom FM 89 from the current down
# wood / live fuel; FMUSRFM forces FMD=89). Both replace the variant's FMCFMD selection; FMBURN calls FMCFMD3 at its start (fmburn.f:246)
# and FMMAIN again after a burn (fmmain.f:190) before FMPOFL, so the burn and the PotFire report each read the stand as it is at that call.

# fmcfmd2.f:205-262 FMBURN pick lists (LARID = the ARID variants; the humid ones are AK/CS/LS/NE/OP/PN/SN/WC/ON)
_fm_larid(v) = !(v isa SoutheastAlaska || v isa CentralStates || v isa LakeStates || v isa Northeast || v isa Olympic ||
                 v isa PacificNorthwest || v isa Southern || v isa WestCascades || v isa Ontario)

function _fm_picklist(larid::Bool, ifcft::Int)::Vector{Int}
    if larid
        ifcft == 1 && return [1, 2, 3, 101, 102, 104, 107]
        ifcft == 2 && return [1, 2, 3, 5, 102, 104, 121, 122, 141, 142]
        ifcft == 3 && return [2, 4, 5, 7, 10, 141, 142, 145, 147, 161, 164, 165]
        return vcat([8, 9], 11:13, 181:189, 201:204)
    else
        ifcft == 1 && return [1, 2, 3, 101, 103, 105, 106, 108, 109]
        ifcft == 2 && return [1, 2, 3, 7, 103, 105, 106, 123, 124, 141, 143, 144]
        ifcft == 3 && return vcat([2, 4, 5, 7, 10, 143, 144, 146, 148, 149], 161:163)
        return vcat([8, 9], 11:13, 181:189, 201:204)
    end
end

# The (load[2,4], sav[2,4], depth, mext) of model `m` for the departure index: defined models resolve through the tables, an undefined number is empty.
function _fm_model_or_empty(s::StandState, m::Int)
    fs = s.fire
    (fs !== nothing && haskey(fs.defulmod, Int32(m))) && return fuel_model_resolved(s, m)
    (m <= size(s.coef.ffe_fuel_models, 1) || haskey(_SB_FUEL_MODELS, m)) && return fuel_model_resolved(s, m)
    return (zeros(Float32, 2, 4), fill(1500f0, 2, 4), 0f0, 0f0)
end

# CURRCWD(J), HERB, WOODY of fmcfmd2.f:132-166 / fmcfmd.f's FMCFMD3 (tons/ac; ×0.04591 later for FMCFMD3's lb/ft²)
function _fm_surface_inputs(s::StandState)
    fs = s.fire; t = s.trees
    cur = zeros(Float32, 11)
    @inbounds for j in 1:11
        for k in 1:2, l in 1:4; cur[j] = cur[j] + fs.cwd[j, k, l]; end
    end
    herb = fs.flive[1]
    woody = 0f0
    pr = isempty(fs.prob_prefire) ? nothing : fs.prob_prefire     # FMPROB: the density the FFE year saw (before this cycle's burn kill)
    @inbounds for i in 1:t.n
        t.height[i] <= fs.user.canmht || continue
        sp = Int(t.species[i])
        xv = _ffe_crownw(s, i, sp, t.dbh[i], t.height[i], Int(round(t.crown_pct[i])))
        p = (pr !== nothing && i <= length(pr)) ? pr[i] : t.tpa[i]
        woody = woody + (xv[1] + 0.5f0 * xv[2]) * p * _FM_P2T
    end
    woody = woody + fs.flive[2]
    return cur, herb, woody
end

"""
    fmcfmd2_select(s) -> Vector{Tuple{Int,Float32}}

FMCFMD2 (fmcfmd2.f:65-472): the fuel models the FIRECALC "new logic" (IFLOGIC=1) burns on. IFCFT (fire-carrying fuel type) from the live fraction /
herb fraction / herb:woody ratio picks the arid/humid list, IFMSET (0 = the 13 originals, 1 = the 40 new, 2 = both) trims it, FMODLIST (act 2550)
forces models on/off, then the two lowest departure indices
  D = .25((FBSAV-FMSAV)/405.2)² + .25((FBBD-FMBD)/.3992)² + .5((FBFFL-FMFFL)/3.051)²
(fuel bed vs model SAV, bulk density, fine-fuel load) are weighted 1/D (exact match / one candidate ⇒ weight 1). No candidate ⇒ FM 8 (+ warning).
"""
function fmcfmd2_select(s::StandState)::Vector{Tuple{Int,Float32}}
    fs = s.fire; u = fs.user
    for a in ffe_find(s, 2550)                                  # FMODLIST: IFUELMON(model) = 0 (force on) / 1 (force off)
        m = Int(trunc(a.params[1])); (1 <= m <= 256) && (u.ifuelmon[m] = Int32(trunc(a.params[2])))
        ffe_done!(s, a)
    end
    cur, herb, woody = _fm_surface_inputs(s)
    fbffl = cur[1] + cur[10] + herb + woody
    ifcft = 4
    if fbffl > 0f0
        livefrac = (herb + woody) / fbffl; herbfrac = herb / fbffl
        herbrat = woody > 0f0 ? herb / woody : 10f0
        if livefrac <= 0.20f0
            ifcft = 4
        elseif herbfrac >= 0.75f0
            ifcft = 1
        elseif herbrat > 2.0f0
            ifcft = 1
        elseif herbrat > 0.25f0
            ifcft = 2
        else
            ifcft = 3
        end
    end
    lfuelmon = falses(256)
    for m in _fm_picklist(_fm_larid(s.variant), ifcft); lfuelmon[m] = true; end
    if u.ifmset == 0
        for i in 101:256; lfuelmon[i] = false; end
    elseif u.ifmset == 1
        for i in 1:13; lfuelmon[i] = false; end
    end
    usav = u.usav; ubd = u.ubd
    # fuel bed bulk density and SAV (fmcfmd2.f:271-318)
    fdfl = cur[1] + cur[10]
    fbbd = 0f0
    if fbffl > 0f0
        wf = fdfl / fbffl
        fbbd = ubd[1] + (wf * (ubd[2] - ubd[1]))
    end
    sa = zeros(Float32, 5)
    sa[1] = usav[1] * fdfl; sa[2] = 109f0 * cur[2]; sa[3] = 30f0 * cur[3]; sa[4] = usav[3] * woody; sa[5] = usav[2] * herb
    sadead = sa[1] + sa[2] + sa[3]
    salive = max(0.0000001f0, sa[4] + sa[5])
    wt = zeros(Float32, 5)
    for i in 1:3
        wt[i] = sadead > 0f0 ? sa[i] / sadead : 0f0
        (salive > 0f0 && i < 3) ? (wt[i + 3] = sa[i + 3] / salive) : (i < 3 && (wt[i + 3] = 0f0))
    end
    fbsav = 0f0
    if !(salive <= 0f0 && sadead <= 0f0)
        wd = sadead / (salive + sadead); wl = 1f0 - wd
        savdead = usav[1] * wt[1] + 109f0 * wt[2] + 30f0 * wt[3]
        savlive = usav[3] * wt[4] + usav[2] * wt[5]
        fbsav = savlive * wl + savdead * wd
    end
    # per-model SAV / bulk density / fine-fuel load and the departure index (fmcfmd2.f:322-402)
    dindex = fill(9999f0, 256)
    for j in 1:256
        ((!lfuelmon[j]) && u.ifuelmon[j] != 0) && continue
        load, sav, depth, _ = _fm_model_or_empty(s, j)
        sav11 = sav[1, 1]; sav12 = sav[1, 2]; sav13 = sav[1, 3]; sav21 = sav[2, 1]; sav22 = sav[2, 2]
        msa = zeros(Float32, 5)
        msa[1] = sav11 * load[1, 1] / 0.04591f0; msa[2] = sav12 * load[1, 2] / 0.04591f0; msa[3] = sav13 * load[1, 3] / 0.04591f0
        msa[4] = sav21 * load[2, 1] / 0.04591f0; msa[5] = sav22 * load[2, 2] / 0.04591f0
        msd = msa[1] + msa[2] + msa[3]; msl = max(0.0000001f0, msa[4] + msa[5])
        w2 = zeros(Float32, 5)
        for i in 1:3
            w2[i] = msd > 0f0 ? msa[i] / msd : 0f0
            (msl > 0f0 && i < 3) ? (w2[i + 3] = msa[i + 3] / msl) : (i < 3 && (w2[i + 3] = 0f0))
        end
        fmsav = 0f0
        if !(msl <= 0f0 && msd <= 0f0)
            wd2 = msd / (msl + msd); wl2 = 1f0 - wd2
            sd = sav11 * w2[1] + sav12 * w2[2] + sav13 * w2[3]
            sl = sav21 * w2[4] + sav22 * w2[5]
            fmsav = sl * wl2 + sd * wd2
        end
        fmffl = (load[1, 1] + load[2, 1] + load[2, 2]) / 0.04591f0
        fmbd = (load[1, 1] + load[1, 2] + load[1, 3] + load[2, 1] + load[2, 2]) / depth
        dindex[j] = 0.25f0 * ((fbsav - fmsav) / 405.2f0)^2 + 0.25f0 * ((fbbd - fmbd) / 0.3992f0)^2 + 0.50f0 * ((fbffl - fmffl) / 3.051f0)^2
    end
    neg = -dindex
    ifm = zeros(Int32, 256)
    rdpsrt!(256, neg, ifm, true)
    lfm = Int[]; lowdi = Float32[]
    for i in 1:256
        m = Int(ifm[i])
        if lfuelmon[m] || u.ifuelmon[m] == 0
            if u.ifuelmon[m] != 1
                push!(lfm, m); push!(lowdi, dindex[m])
            end
        end
    end
    j = min(length(lfm), 2)
    out = Tuple{Int,Float32}[]
    if j == 0
        push!(out, (8, 1f0))                                    # no model can be picked: FM 8 (RCDSET(2) warning)
    elseif lowdi[1] <= 0f0 || j == 1
        push!(out, (lfm[1], 1f0))
    else
        bot = 0f0
        for i in 1:j; bot = bot + 1f0 / lowdi[i]; end
        sumwt = 0f0
        for i in 1:(j - 1)
            w = (1f0 / lowdi[i]) / bot
            push!(out, (lfm[i], w)); sumwt = sumwt + w
        end
        push!(out, (lfm[j], 1f0 - sumwt))
    end
    u.ldynfm || (out = [(out[1][1], 1f0)])
    return out
end

"""
    fmcfmd3_load_fm89!(s)

FMCFMD3 with IFLOGIC=2 (fmcfmd2.f:520-612): load the custom fuel model 89 from the stand's own down wood / live fuel — loads in lb/ft² (tons/ac ×
0.04591: 0-.25" + litter, .25-1", 1-3", live woody = foliage + half the 0-.25" branchwood of trees ≤ CANMHT + shrubs, herb), SAV = USAV / 109 / 30 (INTEGER),
bed depth = (FFL + CURRCWD(2) + CURRCWD(3))/BDavg and moisture of extinction (12 + 480·BDavg/32)/100 with BDavg = UBD(1)+WF·(UBD(2)−UBD(1)). The
table is only replaced when the parameters pass the FMBURN LOK tests.
"""
function fmcfmd3_load_fm89!(s::StandState)
    fs = s.fire; u = fs.user
    cur0, herb0, woody0 = _fm_surface_inputs(s)
    cur = cur0 .* 0.04591f0
    # the Fortran multiplies each CWD cell before summing (CWD·0.04591): re-sum in that order
    c2 = zeros(Float32, 11)
    @inbounds for j in 1:11, k in 1:2, l in 1:4; c2[j] = c2[j] + fs.cwd[j, k, l] * 0.04591f0; end
    cur = c2
    herb = fs.flive[1] * 0.04591f0
    woody = woody0 * 0.04591f0
    xfml = zeros(Float32, 2, 3)
    xfml[1, 1] = max(0f0, cur[1] + cur[10]); xfml[1, 2] = max(0f0, cur[2]); xfml[1, 3] = max(0f0, cur[3])
    xfml[2, 1] = max(0f0, woody); xfml[2, 2] = max(0f0, herb)
    fdfl = cur[1] + cur[10]
    ffl = cur[1] + cur[10] + herb + woody
    wf = fdfl / ffl
    bdavg = u.ubd[1] + (wf * (u.ubd[2] - u.ubd[1]))
    xdep = (ffl + cur[2] + cur[3]) / bdavg
    xext = (12f0 + 480f0 * bdavg / 32f0) / 100f0
    lok = !(xdep < 0f0) && !(xext < 0f0 || xext > 1f0)
    indd = count(>(0f0), xfml[1, :]); inl = count(>(0f0), xfml[2, :])
    (indd <= 0 && inl <= 0) && (lok = false)
    isnan(xdep) && (lok = false)
    if lok
        load = zeros(Float32, 2, 4); sav = zeros(Float32, 2, 4)
        xsur = Float32[u.usav[1] 109f0 30f0; u.usav[3] u.usav[2] 0f0]
        for i in 1:2, j in 1:3
            sav[i, j] = Float32(trunc(Int, xsur[i, j])); load[i, j] = xfml[i, j]
        end
        sav[1, 4] = sav[2, 2]                                  # dead-herb SAV (FMGFMV copies it again)
        fs.defulmod[Int32(89)] = (load, sav, xdep, xext)
    elseif !haskey(fs.defulmod, Int32(89))
        fs.defulmod[Int32(89)] = (zeros(Float32, 2, 4), fill(1500f0, 2, 4), 0f0, 0.15f0)
    end
    return s
end

# FMUSRFM's IFLOGIC=2 forcing (fmusrfm.f:204-215) and IFLOGIC=1's FMCFMD2, for select_fuel_models
function select_iflogic_models(s::StandState)
    fs = s.fire; u = fs.user
    if u.iflogic == 2
        fmcfmd3_load_fm89!(s)
        return [(89, 1f0)]
    end
    return fmcfmd2_select(s)
end

# fmfint.f:108-115 heat content: model 106 (GR6) 9000, else ULHV under IFLOGIC=2, else 8000 BTU/lb
_fm_lhv(s::StandState, fm::Integer)::Float32 =
    fm == 106 ? 9000f0 : (s.fire !== nothing && s.fire.user.iflogic == 2 ? s.fire.user.ulhv : 8000f0)
