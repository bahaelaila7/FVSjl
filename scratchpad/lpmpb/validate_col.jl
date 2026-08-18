# Standalone dump-replay validation of the LPMPB COL (Cole 1983) deterministic
# mortality path against FVSie_lpmpb goldens (lp_dbg.stderr).
# Ports: COLIND, COLDBH(START), COLMOD (default IBOUSE=0), COLMRT (PRKILL/XT).

hex2f32(s) = reinterpret(Float32, parse(UInt32, s; base=16))
f32hex(x::Float32) = uppercase(string(reinterpret(UInt32, x); base=16, pad=8))

# --- MPBINT constants (mpbint.f) ---
const ZINMOR = Float32[0.0,0.0,0.0038,0.0128,0.0206,0.0353,0.0500,0.1429,0.1500,0.1500]
const PRNOIN = Float32[1.0,1.0,0.9935,0.982,0.965,0.909,0.743,0.309,0.285,0.285]

# COLIND (colind.f): DBH -> size class 1..10
function colind(dbh::Float32)::Int
    ival = trunc(Int, dbh)                 # IFIX
    idx  = div(ival, 2) + mod(ival, 2)
    idx > 10 && (idx = 10)
    idx < 1  && (idx = 1)
    return idx
end

# COLMOD default path (colmod.f, ICYC=1, LCURMR=LINVMR=false, IBOUSE=0).
# Returns GREEN(1..NUMYRS, 1..10).
function colmod(start::Vector{Float32}, numyrs::Int)
    DEAD  = zeros(Float32, numyrs, 10)
    GREEN = zeros(Float32, numyrs, 10)
    TDEAD  = zeros(Float32, numyrs)
    TGREEN = zeros(Float32, numyrs)
    for i in 1:10
        DEAD[1,i]  = start[i] * ZINMOR[i]
        GREEN[1,i] = start[i] - DEAD[1,i]
        TDEAD[1]  += DEAD[1,i]
        TGREEN[1] += GREEN[1,i]
    end
    for k in 2:numyrs
        if TDEAD[k-1] <= 0.0005f0 || TGREEN[k-1] <= 0.0005f0
            break
        end
        for j in 1:10
            pnew = PRNOIN[j]                                    # IBOUSE=0
            DEAD[k,j]  = GREEN[k-1,j] * (1.0f0 - (pnew ^ DEAD[k-1,j]))
            TDEAD[k]  += DEAD[k,j]
            GREEN[k,j] = GREEN[k-1,j] - DEAD[k,j]
            TGREEN[k] += GREEN[k,j]
        end
        TDEAD[k] += TDEAD[k-1]
    end
    return GREEN
end

# --- parse goldens ---
recs = Tuple{Int,Int,Float32,Float32,Float32}[]     # (J, IDX, DBH, PROB, XT)
gold_start  = zeros(Float32, 10)
gold_green  = zeros(Float32, 10)   # GREEN(NUMYRS, class)
gold_prkill = zeros(Float32, 10)
numyrs = 0
for ln in eachline("lp_dbg.stderr")
    if startswith(ln, "DBGCM_HDR")
        m = match(r"NUMYRS\s+(\d+)", ln); global numyrs = parse(Int, m.captures[1])
    elseif startswith(ln, "DBGCM_CLS")
        t = split(ln)
        c = parse(Int, t[2])
        gold_start[c]  = hex2f32(t[4])
        gold_green[c]  = hex2f32(t[6])
        gold_prkill[c] = hex2f32(t[8])
    elseif startswith(ln, "DBGCM_REC")
        t = split(ln)
        push!(recs, (parse(Int,t[3]), parse(Int,t[5]),
                     hex2f32(t[7]), hex2f32(t[9]), hex2f32(t[11])))
    end
end

# COLDBH: rebuild START from the per-record DBH/PROB (in oracle IND1/IPT order).
start = zeros(Float32, 10)
for (_, _, dbh, prob, _) in recs
    start[colind(dbh)] += prob
end

# COLMOD + PRKILL (colmrt.f): PRKILL(i) = (START-GREEN(NUMYRS))/START.
GREEN = colmod(start, numyrs)
prkill = zeros(Float32, 10)
for i in 1:10
    prkill[i] = start[i] <= 0.0f0 ? 0.0f0 : (start[i] - GREEN[numyrs,i]) / start[i]
end

# report
ok_start = ok_green = ok_prkill = ok_idx = ok_xt = 0
n_start = n_green = n_prkill = 0
println("== COLDBH START / COLMOD GREEN(NUMYRS) / COLMRT PRKILL ==")
for i in 1:10
    gs = f32hex(gold_start[i]) == f32hex(start[i]);        global ok_start += gs
    gg = f32hex(gold_green[i]) == f32hex(GREEN[numyrs,i]); global ok_green += gg
    gp = f32hex(gold_prkill[i]) == f32hex(prkill[i]);      global ok_prkill += gp
    global n_start += 1; global n_green += 1; global n_prkill += 1
    ulp_g = Int(reinterpret(Int32,GREEN[numyrs,i]) - reinterpret(Int32,gold_green[i]))
    ulp_p = Int(reinterpret(Int32,prkill[i]) - reinterpret(Int32,gold_prkill[i]))
    println("cls $i  START ", f32hex(start[i]), gs ? " ok" : " X("*f32hex(gold_start[i])*")",
            "  GREEN ", f32hex(GREEN[numyrs,i]), gg ? " ok" : " Δ$ulp_g",
            "  PRKILL ", f32hex(prkill[i]), gp ? " ok" : " Δ$ulp_p")
end
println("== COLMRT per-record XT = PRKILL(idx)*PROB ==")
for (j, idx, dbh, prob, xt) in recs
    myidx = colind(dbh)
    gi = myidx == idx; global ok_idx += gi
    myxt = prkill[myidx] * prob
    gx = f32hex(myxt) == f32hex(xt); global ok_xt += gx
    gi && gx && continue
    println("REC J$j idx $myidx(g$idx) XT ", f32hex(myxt), gx ? "" : " X g="*f32hex(xt))
end
println()
println("START  $ok_start/$n_start   GREEN(NUMYRS) $ok_green/$n_green   PRKILL $ok_prkill/$n_prkill")
println("COLIND $ok_idx/$(length(recs))   XT $ok_xt/$(length(recs))")
