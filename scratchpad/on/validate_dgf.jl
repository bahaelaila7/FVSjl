# Validate the WIRED Ontario dgf! end-to-end against production FVSon_g16 WK2.
# Feeds the oracle's EXACT internal per-tree state (DIAM/HT/PCT/IMC + stand BA/RMSQD/SITEAR,
# dumped as Float32-hex from the instrumented dgf.f, unit 773) into a StandState, runs
# species_sort! + dgf!(::Ontario), and compares each tree's WK2 (scratch.wk[2,:]) bit-exact
# to the production WK2 (unit 772, col 6) — the ground truth that USES the leftover BARK.
using FVSjl
using FVSjl: StandState, Southern, Ontario, species_sort!, dgf!, on_bratio, MAXSP

h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
f2h(x::Float32) = uppercase(string(reinterpret(UInt32, x); base = 16, pad = 8))

const INSTR = "/workspace/FVSjl/.claude/worktrees/agent-a66cf29fdfd18b39e/scratchpad/on/instr/rundir"

# --- parse fort.773 (I ISPC IMC | DIAM HT PCT SITEAR BA RMSQD) ---
rows773 = Dict{Int,NamedTuple}()
for ln in eachline(joinpath(INSTR, "fort.773"))
    f = split(strip(ln)); isempty(f) && continue
    I = parse(Int, f[1]); ISPC = parse(Int, f[2]); IMC = parse(Int, f[3])
    hx = f[4]                      # 6×Z8.8 packed into one token by the format? split handles spaces
    # the FORMAT is (3I5,1X,6Z8.8) -> the 6 hexes are concatenated (48 chars) in f[4]
    diam = h2f(hx[1:8]); ht = h2f(hx[9:16]); pct = h2f(hx[17:24])
    sitear = h2f(hx[25:32]); ba = h2f(hx[33:40]); rmsqd = h2f(hx[41:48])
    rows773[I] = (; ISPC, IMC, diam, ht, pct, sitear, ba, rmsqd)
end

# --- parse fort.772 (I ISPC ICYC | DIAM DIAGR BARK DDS COR WK2) ---
wk2_oracle = Dict{Int,Float32}(); diagr_oracle = Dict{Int,Float32}()
barklast_oracle = Dict{Int,Float32}(); dds_oracle = Dict{Int,Float32}()
for ln in eachline(joinpath(INSTR, "fort.772"))
    f = split(strip(ln)); isempty(f) && continue
    I = parse(Int, f[1]); hx = f[4]
    diagr_oracle[I] = h2f(hx[9:16]); barklast_oracle[I] = h2f(hx[17:24])  # LEFTOVER bark (0.94 all)
    dds_oracle[I] = h2f(hx[25:32]); wk2_oracle[I] = h2f(hx[41:48])
end

# --- per-tree BARK ground truth = original per-tree dump ont01_dgfdump.txt (col 8) ---
bark_oracle = Dict{Int,Float32}()
for ln in eachline("/workspace/FVSjl/scratchpad/on/ont01_dgfdump.txt")
    f = split(strip(ln)); isempty(f) && continue
    I = parse(Int, f[1]); bark_oracle[I] = h2f(f[13])   # DIAM DBHM SIM BAM QMDM BALM HTM BARK ...
end

Is = sort(collect(keys(rows773)))
n = length(Is)

# --- build a StandState container and populate the loaded stand exactly as the oracle had it ---
s = StandState(Southern())          # container only; dgf! dispatches on the Ontario() arg below
t = s.trees; p = s.plot; c = s.calib
t.n = n
ba_v = rows773[Is[1]].ba; rmsqd_v = rows773[Is[1]].rmsqd
p.basal_area = ba_v
p.qmd = rmsqd_v
for I in Is
    r = rows773[I]
    t.species[I] = Int32(r.ISPC)
    t.dbh[I] = r.diam
    t.height[I] = r.ht
    t.crown_ratio[I] = r.pct          # FVS PCT
    t.mort_code[I] = Int32(r.IMC)
    t.tpa[I] = 1f0
    p.sp_site_index[r.ISPC] = r.sitear
    c.dg_cor[r.ISPC] = 0f0            # ground-truth COR = 0 at cyc0
    # sort key must be set for species_sort! within-species ordering
    t.sort_key[I] = Int32(I)
end

species_sort!(s)
dgf!(s, Ontario())
wk2 = view(s.scratch.wk, 2, :)

# --- standalone check: on_bratio reproduces the PER-TREE bark (drives DIAGR) ---
println("== on_bratio per-tree (grown D) vs oracle DIAGR/bark ==")
global okbark = 0
for I in Is
    r = rows773[I]
    # recompute grown D via the Penner core to feed on_bratio (matches dgf! internal)
    dbhm_f, _, _ = FVSjl.on_penner_dds(FVSjl.ON_OSPMAP[r.ISPC] + ((r.ISPC in FVSjl.ON_LQUAL_SPP && r.IMC==1) ? 1 : 0),
        (r.ISPC in FVSjl.ON_LQUAL_SPP && r.IMC==1) ? 1 : 0, r.diam,
        r.sitear*FVSjl.ON_FTtoM, r.ba*FVSjl.ON_FT2pACRtoM2pHA, r.rmsqd*FVSjl.ON_INtoCM,
        (1f0-r.pct/100f0)*r.ba*FVSjl.ON_FT2pACRtoM2pHA, r.ht*FVSjl.ON_FTtoM, 1f0)
    bark_j = on_bratio(r.ISPC, dbhm_f*FVSjl.ON_CMtoIN, r.ht)
    m = f2h(bark_j) == f2h(bark_oracle[I])
    global okbark += m
    println("  I=$I ISPC=$(r.ISPC)  bark jl=$(f2h(bark_j)) or=$(f2h(bark_oracle[I]))  $(m ? "OK" : "MISMATCH ($(bark_j) vs $(bark_oracle[I]))")")
end
println("  bark per-tree: $okbark/$n bit-exact\n")

# --- MAIN go/no-go: WK2 (with production leftover-bark DDS) ---
println("== dgf! WK2 vs production FVSon_g16 (fort.772) ==")
global okwk = 0
for I in Is
    m = f2h(wk2[I]) == f2h(wk2_oracle[I])
    global okwk += m
    println("  I=$I ISPC=$(rows773[I].ISPC)  WK2 jl=$(f2h(wk2[I])) or=$(f2h(wk2_oracle[I]))  $(m ? "OK" : "MISMATCH ($(wk2[I]) vs $(wk2_oracle[I]))")")
end
println("\nRESULT: bark $okbark/$n  WK2 $okwk/$n bit-exact vs FVSon_g16")
println(okwk == n && okbark == n ? "PASS" : "FAIL")
