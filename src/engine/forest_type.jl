# =============================================================================
# forest_type.jl — FIA forest-type classification (FORTYP + STKVAL)
#
# Ported from: base/fortyp.jl + base/stkval.jl (FVSjulia, validated). The decision
# tree is reproduced faithfully (it only READS the per-group stocking array `s`
# and a few stand scalars); the data tables (TAB2 stocking coeffs, TAB3 FIA→group
# /equation maps, the valid FIA type codes) come from data/<variant>/*.csv via
# `coef.stock_b0/b1`, `coef.fia_group/fia_stock_eq`, `coef.forest_type_codes`.
#
# `compute_forest_type!` sets `plot.forest_type` (IFORTP), which the diameter-growth
# model reads for its forest-type coefficient term. Stand size/stocking class
# (ISZCL/ISTCL) is not needed for the type code and is omitted.
# =============================================================================

# Normalize an FIA species-code string to its integer code (blank/0/out-of-range
# fall back to 998/999, matching stkval.f's _get_fia_and_coeffs).
@inline function _fia_code(fs::AbstractString)
    z = strip(fs)
    ifia = isempty(z) ? 998 : something(tryparse(Int, z), 998)
    ifia == 0 && (ifia = 999)
    (ifia < 1 || ifia > 1000) && (ifia = 999)
    return ifia
end

# STKVAL: per-ITG-group stocking array s[1:210] from the live tree list. Stocking
# per tree = b0·D^b1·PROB·cf (cf scales sub-5" trees; a "future stand" correction
# lifts D to 5" when large-tree stocking < 20). s[g] sums species mapped to group g.
function _stkval_stocking(st::StandState)
    coef = st.coef
    isct = st.control.sp_count_tab; ind1 = st.scratch.idx1
    dbh = st.trees.dbh; prob = st.trees.tpa; fiajsp = coef.code_fia
    b0t = coef.stock_b0; b1t = coef.stock_b1
    eq999 = get(coef.fia_stock_eq, 999, 26)
    @inline function _coeffs(ispc::Int)
        ifia = _fia_code(fiajsp[ispc])
        eq = get(coef.fia_stock_eq, ifia, eq999)
        (eq < 1 || eq > 36) && (eq = eq999)
        b0 = b0t[eq]; b1 = b1t[eq]
        (b0 == 0f0 || b1 == 0f0) && (eq = eq999; b0 = b0t[eq]; b1 = b1t[eq])
        return ifia, b0, b1
    end
    dmxss = 0f0
    @inbounds for ispc in 1:MAXSP
        isct[ispc,1] == 0 && continue
        for i3 in isct[ispc,1]:isct[ispc,2]
            d = dbh[Int(ind1[i3])]
            (d < 5f0 && d > dmxss) && (dmxss = d)
        end
    end
    ttst51 = 0f0
    @inbounds for ispc in 1:MAXSP
        isct[ispc,1] == 0 && continue
        _, b0, b1 = _coeffs(ispc)
        for i3 in isct[ispc,1]:isct[ispc,2]
            d = dbh[Int(ind1[i3])]
            (d > 0f0 && d >= 5f0) && (ttst51 += b0 * d^b1 * prob[Int(ind1[i3])])
        end
    end
    _cf(d) = d >= 5f0 ? 1f0 : (dmax = ttst51 >= 10f0 ? 5f0 : dmxss; dmax > 0f0 ? d / dmax : 0f0)
    ttst52 = 0f0
    @inbounds for ispc in 1:MAXSP
        isct[ispc,1] == 0 && continue
        _, b0, b1 = _coeffs(ispc)
        for i3 in isct[ispc,1]:isct[ispc,2]
            i = Int(ind1[i3]); d = dbh[i]; d <= 0f0 && continue
            d >= 5f0 && (ttst52 += b0 * d^b1 * prob[i] * _cf(d))
        end
    end
    s = zeros(Float32, 210)
    totstk = 0f0; szcl1 = 0f0; szcl2 = 0f0; szcl3 = 0f0   # for size/stocking class
    @inbounds for ispc in 1:MAXSP
        isct[ispc,1] == 0 && continue
        ifia, b0, b1 = _coeffs(ispc)
        grp = get(coef.fia_group, ifia, 0)
        acc = 0f0
        for i3 in isct[ispc,1]:isct[ispc,2]
            i = Int(ind1[i3]); d = dbh[i]; d <= 0f0 && continue
            d_adj = (ttst52 < 20f0 && d < 5f0) ? 5f0 : d
            stk = b0 * d_adj^b1 * prob[i] * _cf(d)
            acc += stk
            totstk += stk
            # size-class bin by ORIGINAL dbh (stkval.f:466-477); softwood (FIA<300)
            # large ≥9", hardwood (≥300) large ≥11".
            if d < 5f0;                              szcl1 += stk
            elseif (ifia < 300 ? d < 9f0 : d < 11f0); szcl2 += stk
            else;                                     szcl3 += stk
            end
        end
        (1 <= grp <= 210) && (s[grp] += acc)
    end
    # stocking class (stkval.f:498) and size class (stkval.f:487-496)
    istcl = totstk > 100f0 ? 1 : totstk >= 60f0 ? 2 : totstk >= 35f0 ? 3 : totstk >= 10f0 ? 4 : 5
    iszcl = totstk < 10f0 ? 5 : szcl1 > totstk * 0.5f0 ? 3 : szcl2 > szcl3 ? 2 : 1
    return s, Int32(iszcl), Int32(istcl)
end

"""
    compute_forest_type!(state) -> Int32

FORTYP: classify the stand into an FIA forest-type code (e.g. 520 = mixed upland
hardwoods) from the per-group stocking array, and store it in `plot.forest_type`.
"""
function compute_forest_type!(st::StandState)
    KODFOR = Int(st.plot.user_forest_code)
    VARACD = "SN"
    south  = 0f0
    species_sort!(st)
    s, iszcl, istcl = _stkval_stocking(st)
    st.plot.size_class = iszcl; st.plot.stocking_class = istcl
    sftwds = sum(s[1:58]) + s[60] + sum(s[62:79]) + s[161] + s[162] + s[170]
    trfirspr = sum(s[1:5]) + s[7] + s[9] + s[14] + s[15] + s[28] + s[34] + s[35]
    spsafir  = s[4] + s[14] + s[15]
    engsafir = s[4] + s[14]
    salpfir  = s[4]
    engspr   = s[14]
    bluspr   = s[15]
    whmlcks  = s[34] + s[35]
    whmlck   = s[34]
    mtnhmlk  = s[35]
    trufir   = sum(s[1:5]) + s[7]
    pslvrfir = s[1]; whtfir = s[2]; grndfir = s[3]
    redfir   = s[5]; noblfir = s[7]
    akylcdr  = s[9]; wwhpin = s[28]
    dflrwp   = s[8]+s[10]+s[11]+s[13]+s[23]+s[24]+s[26]+s[27]+s[30]+s[31]+s[36]
    dgfrlr   = s[11]+s[13]+s[31]
    dgfir    = s[31]; wlrch = s[13]; wrcedar = s[11]
    dgfrpin  = s[8]+s[10]+s[23]+s[24]+s[26]+s[27]+s[30]+s[31]+s[36]
    pndrosa  = s[26]+s[36]; lodgpole = s[23]
    porfcdr  = s[8]; sugrpin = s[27]; inscdr = s[10]
    jfbcdgfr = s[24]+s[30]+s[36]
    lrchpin  = s[13]+s[23]+s[26]+s[36]
    stksphm  = s[11]+s[18]+s[34]
    stkaspr  = s[18]
    othwpin  = s[6]+s[12]+s[19]+s[20]+s[21]+s[22]+s[25]+s[29]+s[40]
    knbpin   = s[21]; swhtpin = s[22]; bshppin = s[6]
    mntrypin = s[29]; fxtlpin = s[20]; lmbrpin = s[25]; whbrkpin = s[19]
    miscwsfw = s[12]+s[40]
    redwds   = s[31]+s[32]+s[33]; redwood = s[32]; gntseq = s[33]
    epine    = sum(s[41:54]) + s[66]
    rwjpin   = s[41]+s[42]+s[53]+s[66]
    whphem   = s[53]+s[66]; ewhpin = s[53]; ehmlck = s[66]
    redpin   = s[42]; jackpin = s[41]
    lnglfslh = s[46]+s[48]; lnglfp = s[48]; slashp = s[46]
    lobshrtp = s[44]+s[45]+s[47]+s[49]+s[50]+s[51]+s[52]+s[54]
    lobp     = s[52]; shrtlfp = s[45]; virgp = s[54]
    sandp    = s[44]; tblmtnp = s[49]; pondp = s[51]; pitchp = s[50]
    sprucep  = s[47]
    pnynjn   = s[38]+s[63]+s[64]+s[161]+s[162]
    ercedar  = s[64]; rmjunipr = s[63]; wjunipr = s[38]
    juniprw  = s[161]; pinjunw = s[161]+s[162]
    esprfir  = s[16]+s[17]+s[55]+s[58]+s[60]+s[65]
    usprfir  = s[16]+s[55]+s[58]; bfirspr = s[55]+s[58]
    blsmfir  = s[55]; redspr = s[58]; whspr = s[16]
    lsprfir  = s[17]+s[60]+s[65]; blkspr = s[17]; tamrck = s[65]; nwhcdr = s[60]
    exsftwd  = s[70]+s[71]+s[72]; scotchp = s[71]; japblkp = s[72]; othexswd = s[70]
    othsft   = s[170]

    # --- HARDWOODS
    hrdwds = sum(s[81:153]) + sum(s[156:160]) + sum(s[201:210]) +
             s[59] + s[61] + s[163] + s[180] + s[190]
    oakpin = sum(s[41:54]) + s[64]
    ercedar = s[64]; shrtlfp = s[45]; ewhpin = s[53]; lnglfp = s[48]
    virgp = s[54]; lobp = s[52]; slashp = s[46]; jackpin = s[41]; redpin = s[42]
    sandp = s[44]; sprucep = s[47]; tblmtnp = s[49]; pitchp = s[50]; pondp = s[51]
    oakhck = sum(s[81:86]) + s[88]+s[89]+s[92]+s[93]+s[101]+s[108]+s[110]+
             s[120]+s[122]+s[202]+s[206]+s[207]
    whoak = s[81]; buroak = s[83]; chstoak = s[84]; nroak = s[85]
    scrltoak = s[82]; yp = s[110]; bwalnut = s[108]; blcst = s[122]
    redmapl = s[95]
    pstbljko = s[86]+s[206]
    chblsc = s[82]+s[84]+s[120]
    worohk = s[81]+s[85]+s[92]+s[94]+s[120]+s[207]
    ypworo = s[81]+s[85]+s[110]
    scruboak = s[89]+s[203]+s[206]
    swgyp = s[109]+s[110]
    sasfprsm = s[93]
    mxdhwd = s[83]+s[88]+s[94]+s[101]+s[106]+s[108]+s[113]+s[122]+s[125] +
             sum(s[201:204])
    okgmcyp = s[59]+s[61]+s[87]+s[90]+s[111]+s[112]+s[114]+s[127]+s[128]+s[143]
    atlwcdr = s[59]
    sbstrm = s[95]+s[111]+s[113]+s[114]+s[127]
    schchbo = s[87]
    swgwilo = s[109]+s[125]+s[143]+s[201]+s[203]+s[204]
    cypwtup = s[61]+s[112]
    ovrcupwh = s[90]+s[128]
    elmashcw = s[91]+s[97]+s[100]+s[104]+s[115]+s[116]+s[118]+s[123]+
               s[129]+s[135]+s[137]+s[208]
    ctnwd = s[118]+s[137]; willow = s[123]
    rbrchsyc = s[108]+s[116]+s[123]+s[129]
    sycpcelm = s[91]+s[94]+s[109]+s[116]
    baelmmap = s[104]
    sgbelmga = s[94]+s[100]+s[105]+s[115]+s[208]+s[209]
    slvmaelm = s[94]+s[97]
    ctnwwilo = s[118]+s[123]+s[130]+s[131]+s[137]
    orgnash = s[135]
    mbb = s[66]+s[96]+s[98]+s[107]+s[110]+s[122]+s[124]
    blch = s[121]
    bcwayp = s[103]+s[110]+s[121]
    mplbaswd = s[96]+s[124]
    elmashlo = s[94]+s[105]+s[122]
    mplbchyb = s[66]+s[94]+s[95]+s[96]+s[98]+s[102]+s[105]+s[107]+s[108]
    aspbrch = s[99]+s[117]+s[119]; aspen = s[119]; blsmpop = s[117]; paprbrch = s[99]
    aldrmapl = s[130]+s[131]; redaldr = s[131]; bglfmpl = s[130]
    wstoak = s[134]+s[138]+s[139]+s[140]+s[142]+s[158]+s[163]+s[210]
    clfbo = s[139]; orgwho = s[140]; blueoak = s[134]; graypine = s[163]
    costlo = s[138]; canylo = s[142]; doakwdld = s[158]; eoakwdld = s[210]
    tanoklrl = s[133]+s[136]+s[141]; tanoak = s[136]; clflrel = s[141]; gntchnk = s[133]
    othwhwd = s[132]+s[156]+s[157]+s[159]+s[160]
    pacmdrn = s[132]; mesquitw = s[157]; mtnbrshw = s[156]
    intmaplw = s[159]; miscwhwd = s[160]
    trophwds = s[147]+s[149]; palm = s[147]; mangrv = s[149]
    exhwds = s[144]+s[145]+s[146]+s[148]
    palonia = s[144]; meluca = s[145]; euclpt = s[148]; othexhwd = s[146]

    # --- Special groups
    uplooak = s[125]+s[201]+s[203]+s[204]
    uplohwd = s[95]+s[103]+s[105]
    sroak = s[88]; aelm = s[94]; wcelm = s[209]; slvrmpl = s[97]
    whash = s[103]; ectnwd = s[118]
    blkgum = s[113]; beech = s[102]; holly = s[106]; sweetgum = s[109]; pinoak = s[205]
    othall = s[170]+s[180]+s[190]; othsft = s[170]; othhrd = s[180]; othspc = s[190]

    # Total stocking
    totstk = sum(s)

    # Upland/lowland category totals
    upmpbcbr = s[66]+s[85]+s[94]+s[95]+s[96]+s[98]+s[102]+s[103]+
               s[105]+s[107]+s[108]+s[110]+s[121]+s[122]+s[124]
    upoakhic = s[81]+s[82]+s[83]+s[84]+s[85]+s[86]+s[88]+s[89]+s[92]+
               s[93]+s[94]+s[95]+s[101]+s[102]+s[103]+s[105]+s[106]+
               s[108]+s[109]+s[110]+s[113]+s[120]+s[121]+s[122]+s[125]+
               sum(s[201:204])+s[206]+s[207]+s[209]
    lwelascw = s[91]+s[94]+s[95]+s[97]+s[100]+s[103]+s[104]+s[105]+
               s[108]+s[109]+s[115]+s[116]+s[118]+s[123]+s[129]+
               s[130]+s[131]+s[135]+s[137]+s[208]+s[209]
    lwokgmcy = s[59]+s[61]+s[87]+s[90]+s[94]+s[95]+s[102]+s[103]+
               s[105]+s[106]+s[109]+s[111]+s[112]+s[113]+s[114]+
               s[125]+s[127]+s[128]+s[143]+s[201]+s[203]+s[204]+s[205]

    # Proportion thresholds
    pcond = Float32(1)
    half   = Float32(0.50) * totstk
    quartr = Float32(0.25) * totstk
    five   = Float32(0.05) * totstk
    ten    = Float32(0.10) * totstk
    fiftn  = Float32(0.15) * totstk
    twenty = Float32(0.20) * totstk
    eighty = Float32(0.80) * totstk
    ninety  = Int(floor(Float32(0.90) * totstk))
    ninty5  = Int(floor(Float32(0.95) * totstk))

    # Forest service region from KODFOR
    irgn = Int(0)
    if KODFOR ÷ 10000000 >= 1
        irgn = KODFOR ÷ 10000000
    elseif KODFOR ÷ 10000 >= 1
        irgn = KODFOR ÷ 10000
    elseif KODFOR ÷ 100 >= 1
        irgn = KODFOR ÷ 100
    end

    # ---------------------------------------------------------------------------
    # Main decision tree
    # ---------------------------------------------------------------------------
    ift = Int32(0)

    if totstk / pcond < Float32(10)
        st.plot.forest_type = Int32(999)   # NONSTOCKED — skip the decision tree
        return Int32(999)
    end


    if sftwds >= hrdwds  # SOFTWOODS
        ift = Int32(996)
        p1 = max(trfirspr,stksphm,redwds,dflrwp,othwpin,
                 epine,esprfir,pnynjn,exsftwd,othsft) - Float32(0.001)
        if p1 > Float32(0)
            if p1 < trfirspr                   # TRUE FIRS AND SPRUCE
                p2 = max(spsafir,whmlcks,trufir,akylcdr,wwhpin) - Float32(0.001)
                if p2 > Float32(0)
                    if p2 < spsafir
                        if engsafir >= bluspr
                            if (salpfir >= five && salpfir < half) && (engspr >= five && engspr < half)
                                ift = Int32(266)   # ENGELMANN SPRUCE-SUBALPINE FIR
                            elseif engspr >= salpfir
                                ift = Int32(265)   # ENGELMANN SPRUCE
                            elseif salpfir > Float32(0)
                                ift = Int32(268)   # SUBALPINE FIR
                            end
                        elseif bluspr > Float32(0)
                            ift = Int32(269)       # BLUE SPRUCE
                        end
                    elseif p2 < whmlcks
                        if whmlck >= mtnhmlk
                            ift = Int32(301)       # WESTERN HEMLOCK
                        elseif mtnhmlk > Float32(0)
                            ift = Int32(270)       # MOUNTAIN HEMLOCK
                        end
                    elseif p2 < trufir
                        p3 = max(pslvrfir,whtfir,grndfir,salpfir,redfir,noblfir) - Float32(0.001)
                        if p3 > Float32(0)
                            if p3 < pslvrfir
                                ift = Int32(264)   # PACIFIC SILVER FIR
                            elseif p3 < whtfir
                                ift = Int32(261)   # WHITE FIR
                            elseif p3 < grndfir
                                ift = Int32(267)   # GRAND FIR
                            elseif p3 < salpfir
                                ift = Int32(268)   # SUBALPINE FIR
                            elseif p3 < redfir
                                ift = Int32(262)   # RED FIR
                            elseif noblfir > Float32(0)
                                ift = Int32(263)   # NOBLE FIR
                            end
                        end
                    elseif p2 < akylcdr
                        ift = Int32(271)           # ALASKA YELLOW CEDAR
                    elseif wwhpin > Float32(0)
                        ift = Int32(241)           # WESTERN WHITE PINE
                    end
                end
            elseif p1 < stksphm                # SITKA SPRUCE-HEMLOCK
                p2 = max(whmlck,stkaspr,wrcedar) - Float32(0.001)
                if p2 > Float32(0)
                    if p2 < whmlck
                        ift = Int32(301)           # WESTERN HEMLOCK
                    elseif p2 < stkaspr
                        ift = Int32(305)           # SITKA SPRUCE
                    elseif wrcedar > Float32(0)
                        ift = Int32(304)           # WESTERN REDCEDAR
                    end
                end
            elseif p1 < redwds                 # REDWOODS
                p2 = max(redwood,gntseq,dgfir) - Float32(0.001)
                if p2 > Float32(0)
                    if p2 < redwood
                        ift = Int32(341)           # REDWOOD
                    elseif p2 < gntseq
                        ift = Int32(342)           # GIANT SEQUOIA
                    elseif dgfir > Float32(0)
                        ift = Int32(201)           # DOUGLAS-FIR
                    end
                end
            elseif p1 < dflrwp                 # DF-LARCH-W.WHITE PINES
                p2 = max(dgfrlr,dgfrpin,lrchpin) - Float32(0.001)
                if p2 > Float32(0)
                    if p2 < dgfrlr
                        p3 = max(dgfir,wlrch,wrcedar) - Float32(0.001)
                        if p3 > Float32(0)
                            if p3 < dgfir
                                ift = Int32(201)   # DOUGLAS-FIR
                            elseif p3 < wlrch
                                ift = Int32(321)   # WESTERN LARCH
                            elseif wrcedar > Float32(0)
                                ift = Int32(304)   # WESTERN REDCEDAR
                            end
                        end
                    elseif p2 < dgfrpin
                        p3 = max(dgfir,pndrosa,lodgpole,porfcdr,sugrpin,inscdr,jfbcdgfr) - Float32(0.001)
                        if p3 > Float32(0)
                            if p3 < dgfir
                                ift = Int32(201)   # DOUGLAS-FIR
                            elseif p3 < pndrosa
                                ift = Int32(221)   # PONDEROSA PINE
                            elseif p3 < lodgpole
                                ift = Int32(281)   # LODGEPOLE PINE
                            elseif p3 < porfcdr
                                ift = Int32(202)   # PORT ORFORD CEDAR
                            elseif p3 < sugrpin
                                ift = Int32(224)   # SUGAR PINE
                            elseif p3 < inscdr
                                ift = Int32(222)   # INCENSE CEDAR
                            elseif jfbcdgfr > Float32(0)
                                ift = Int32(223)   # JEFFREY-COULTER PINE-BIGCONE DF
                            end
                        end
                    elseif p2 < lrchpin
                        p3 = max(wlrch,pndrosa,lodgpole) - Float32(0.001)
                        if p3 > Float32(0)
                            if p3 < wlrch
                                ift = Int32(321)   # WESTERN LARCH
                            elseif p3 < pndrosa
                                ift = Int32(221)   # PONDEROSA PINE
                            elseif lodgpole > Float32(0)
                                ift = Int32(281)   # LODGEPOLE PINE
                            end
                        end
                    end
                end
            elseif p1 < othwpin                # OTHER WESTERN PINES
                p2 = max(knbpin,swhtpin,bshppin,mntrypin,fxtlpin,lmbrpin,whbrkpin,miscwsfw) - Float32(0.001)
                if p2 > Float32(0)
                    if p2 < knbpin
                        ift = Int32(361)           # KNOBCONE PINE
                    elseif p2 < swhtpin
                        ift = Int32(362)           # SW WHITE PINE
                    elseif p2 < bshppin
                        ift = Int32(363)           # BISHOP PINE
                    elseif p2 < mntrypin
                        ift = Int32(364)           # MONTEREY PINE
                    elseif p2 < fxtlpin
                        ift = Int32(365)           # FOXTAIL-BRISTLECONE PINE
                    elseif p2 < lmbrpin
                        ift = Int32(366)           # LIMBER PINE
                    elseif p2 < whbrkpin
                        ift = Int32(367)           # WHITEBARK PINE
                    elseif miscwsfw > Float32(0)
                        ift = Int32(368)           # MISC. WESTERN SOFTWOODS
                    end
                end
            elseif p1 < epine                  # EASTERN PINES
                p2 = max(rwjpin,lnglfslh,lobshrtp) - Float32(0.001)
                if p2 > Float32(0)
                    if p2 < rwjpin
                        p3 = Float32(0)
                        if (whphem >= half) && (ewhpin >= five && ewhpin < half) &&
                           (ehmlck >= five && ehmlck < half)
                            ift = Int32(104)       # E. WHITE PINE-HEMLOCK
                        else
                            p3 = max(ewhpin,redpin,jackpin,ehmlck) - Float32(0.001)
                        end
                        if p3 > Float32(0)
                            if p3 < ewhpin
                                ift = Int32(103)   # E. WHITE PINE
                            elseif p3 < redpin
                                ift = Int32(102)   # RED PINE
                            elseif p3 < jackpin
                                ift = Int32(101)   # JACK PINE
                            elseif ehmlck > Float32(0)
                                ift = Int32(105)   # EASTERN HEMLOCK
                            end
                        end
                    elseif p2 < lnglfslh
                        if lnglfp >= slashp
                            ift = Int32(141)       # LONGLEAF PINE
                        elseif slashp > Float32(0)
                            ift = Int32(142)       # SLASH PINE
                        end
                    elseif lobshrtp > Float32(0)
                        p3 = max(lobp,shrtlfp,virgp,sandp,tblmtnp,pondp,pitchp,sprucep) - Float32(0.001)
                        if p3 > Float32(0)
                            if p3 < lobp
                                ift = Int32(161)   # LOBLOLLY PINE
                            elseif p3 < shrtlfp
                                ift = Int32(162)   # SHORTLEAF PINE
                            elseif p3 < virgp
                                ift = Int32(163)   # VIRGINIA PINE
                            elseif p3 < sandp
                                ift = Int32(164)   # SAND PINE
                            elseif p3 < tblmtnp
                                ift = Int32(165)   # TABLE MOUNTAIN PINE
                            elseif p3 < pondp
                                ift = Int32(166)   # POND PINE
                            elseif p3 < pitchp
                                ift = Int32(167)   # PITCH PINE
                            elseif sprucep > Float32(0)
                                ift = Int32(168)   # SPRUCE PINE
                            end
                        end
                    end
                end
            elseif p1 < esprfir                # EASTERN SPRUCE-FIR
                if usprfir >= lsprfir           # UPLAND E. SPRUCE-FIR
                    if (bfirspr >= half) && (blsmfir >= five && blsmfir < half) &&
                       (redspr >= five && redspr < half)
                        ift = Int32(124)           # RED SPRUCE-BALSAM FIR
                    else
                        p2 = max(blsmfir,whspr,redspr) - Float32(0.001)
                        if p2 > Float32(0)
                            if p2 < blsmfir
                                ift = Int32(121)   # BALSAM FIR
                            elseif p2 < whspr
                                ift = Int32(122)   # WHITE SPRUCE
                            elseif redspr > Float32(0)
                                ift = Int32(123)   # RED SPRUCE
                            end
                        end
                    end
                else                            # LOWLAND E. SPRUCE-FIR
                    p2 = max(blkspr,tamrck,nwhcdr) - Float32(0.001)
                    if p2 > Float32(0)
                        if p2 < blkspr
                            ift = Int32(125)       # BLACK SPRUCE
                        elseif p2 < tamrck
                            ift = Int32(126)       # TAMARACK
                        elseif nwhcdr > Float32(0)
                            ift = Int32(127)       # NORTHERN WHITE-CEDAR
                        end
                    end
                end
            elseif p1 < pnynjn                 # PINYON-JUNIPER
                p2 = max(ercedar,rmjunipr,wjunipr,juniprw,pinjunw) - Float32(0.001)
                if p2 > Float32(0)
                    if p2 < ercedar
                        ift = Int32(181)           # EASTERN REDCEDAR
                    elseif p2 < rmjunipr
                        ift = Int32(182)           # ROCKY MTN. JUNIPER
                    elseif p2 < wjunipr
                        ift = Int32(183)           # WESTERN JUNIPER
                    elseif p2 < juniprw
                        ift = Int32(184)           # JUNIPER WOODLAND
                    elseif pinjunw > Float32(0)
                        ift = Int32(185)           # PINYON-JUNIPER WOODLAND
                    end
                end
            elseif p1 < exsftwd                # EXOTIC SOFTWOODS
                p2 = max(scotchp,japblkp,othexswd) - Float32(0.001)
                if p2 > Float32(0)
                    if p2 < scotchp
                        ift = Int32(381)           # SCOTCH PINE
                    elseif p2 < japblkp
                        ift = Int32(383)           # JAPANESE BLACK PINE
                    elseif othexswd > Float32(0)
                        ift = Int32(383)           # OTHER EXOTIC SOFTWOODS
                    end
                end
            elseif p1 < othsft                 # FVS OTHER SOFTWOODS
                if othsft > Float32(0)
                    ift = Int32(996)
                end
            end
        end

    else  # HARDWOODS
        ift = Int32(997)

        if oakpin >= quartr                    # OAK-PINE
            p1 = max(ercedar,shrtlfp,ewhpin,lnglfp,virgp,lobp,slashp,
                     jackpin,redpin,sandp,tblmtnp,pitchp,pondp,sprucep) - Float32(0.001)
            if p1 > Float32(0)
                if p1 < ercedar
                    ift = Int32(402)               # E. REDCEDAR-HARDWOOD
                elseif p1 < shrtlfp
                    ift = Int32(404)               # SHORTLEAF PINE-OAK
                elseif p1 < ewhpin
                    ift = Int32(401)               # E. WHITE PINE-RED OAK-WHITE ASH
                elseif p1 < lnglfp
                    ift = Int32(403)               # LONGLEAF PINE-OAK
                elseif p1 < virgp
                    ift = Int32(405)               # VIRGINIA PINE-S. RED OAK
                elseif p1 < lobp
                    ift = Int32(406)               # LOBLOLLY PINE-HARDWOOD
                elseif p1 < slashp
                    ift = Int32(407)               # SLASH PINE-HARDWOOD
                elseif jackpin+redpin+sandp+tblmtnp+pitchp+pondp+sprucep > Float32(0)
                    ift = Int32(409)               # OTHER PINE-HARDWOOD
                end
            end
        else
            # Adjust group totals before upland/lowland decision
            if pstbljko > Float32(0.01); pstbljko = pstbljko + sroak; end
            if chstoak   > Float32(0.01); chstoak   = chstoak   + sroak; end
            if baelmmap  > Float32(0.01)
                baelmmap = baelmmap + aelm + redmapl + slvrmpl + whash + ectnwd
            end
            upland = Float32(0)
            p1 = max(lwelascw,lwokgmcy,upmpbcbr,upoakhic) - Float32(0.001)
            if p1 > Float32(0)
                if p1 < lwelascw
                    upland = Float32(0)
                elseif p1 < lwokgmcy
                    upland = Float32(0)
                elseif p1 < upmpbcbr
                    upland = Float32(1)
                elseif upoakhic > Float32(0)
                    upland = Float32(1)
                end
            end

            if upland == Float32(1)             # UPLAND site
                oakhck = oakhck + sweetgum + uplooak + aelm + wcelm + holly + blkgum
                south = VARACD == "SN" ? Float32(1) : Float32(0)
                if south == Float32(1)
                    oakhck = oakhck + uplohwd + beech + blch
                else
                    if redmapl >= half
                        oakhck = oakhck + uplohwd + beech + blch
                    elseif (uplohwd + beech + blch) < half
                        if oakhck > five
                            oakhck = oakhck + uplohwd + beech + blch
                        end
                    end
                    mbb = mbb + uplohwd + beech + blch + aelm
                    if mbb > five
                        if (nroak + yp + bwalnut) < half
                            mbb = mbb + nroak + yp + bwalnut
                        end
                    end
                end
            else                                # LOWLAND site
                upland = Float32(0)
                okgmcyp = okgmcyp + uplooak
                if okgmcyp > Float32(0.01)
                    okgmcyp = okgmcyp + uplohwd + sweetgum + aelm + beech + holly + blkgum + pinoak
                end
                elmashcw = elmashcw + aelm + wcelm
                if elmashcw > Float32(0.01)
                    elmashcw = elmashcw + uplohwd + sweetgum + bwalnut
                elseif okgmcyp < Float32(0.01)
                    if oakhck > mbb
                        oakhck = oakhck + uplohwd + blch + beech
                    else
                        mbb = mbb + uplohwd + blch + beech
                    end
                end
            end


            p1 = max(mbb,elmashcw,oakhck,okgmcyp,aspbrch,aldrmapl,
                     wstoak,tanoklrl,othwhwd,trophwds,exhwds,
                     othhrd,othspc) - Float32(0.001)
            if p1 > Float32(0)
                if p1 < mbb                     # MAPLE-BEECH-BIRCH
                    if blch >= half
                        ift = Int32(802)           # BLACK CHERRY
                    elseif redmapl >= half
                        ift = Int32(809)           # RED MAPLE-UPLAND
                    else
                        p2 = max(bcwayp,mplbaswd,elmashlo,mplbchyb) - Float32(0.001)
                        if p2 > Float32(0)
                            if p2 < bcwayp
                                ift = Int32(803)   # B.CHERRY-W.ASH-Y.POPLAR
                            elseif p2 < mplbaswd
                                ift = Int32(805)   # HARD MAPLE-BASSWOOD
                            elseif p2 < elmashlo
                                ift = Int32(807)   # ELM-ASH-LOCUST
                            elseif mplbchyb > Float32(0)
                                ift = Int32(801)   # MAPLE-BEECH-YELLOW BIRCH
                            end
                        end
                    end
                elseif p1 < elmashcw            # ELM-ASH-COTTONWOOD
                    if ctnwd >= half
                        ift = Int32(703)           # COTTONWOOD
                    elseif willow >= half
                        ift = Int32(704)           # WILLOW
                    elseif redmapl >= half
                        ift = Int32(708)           # RED MAPLE-LOWLAND
                    else
                        p2 = max(rbrchsyc,sycpcelm,baelmmap,sgbelmga,
                                 slvmaelm,ctnwwilo,orgnash) - Float32(0.001)
                        if p2 > Float32(0)
                            if p2 < rbrchsyc
                                ift = Int32(702)   # RIVER BIRCH-SYCAMORE
                            elseif p2 < sycpcelm
                                ift = Int32(705)   # SYCAMORE-PECAN-ELM
                            elseif p2 < baelmmap
                                ift = Int32(701)   # B.ASH-ELM-MAPLE
                            elseif p2 < sgbelmga
                                ift = Int32(706)   # SUGARBERRY-HACKBERRY-ELM-GREEN ASH
                            elseif p2 < slvmaelm
                                ift = Int32(707)   # SILVER MAPLE-ELM
                            elseif p2 < ctnwwilo
                                ift = Int32(709)   # COTTONWOOD-WILLOW
                            elseif orgnash > Float32(0)
                                ift = Int32(722)   # OREGON ASH
                            end
                        end
                    end
                elseif p1 < oakhck              # OAK-HICKORY
                    if whoak >= half
                        ift = Int32(504)           # WHITE OAK
                    elseif buroak >= half
                        ift = Int32(509)           # BUR OAK
                    elseif chstoak >= half
                        ift = Int32(502)           # CHESTNUT OAK
                    elseif nroak >= half
                        ift = Int32(505)           # N. RED OAK
                    elseif scrltoak >= half
                        ift = Int32(510)           # SCARLET OAK
                    elseif yp >= half
                        ift = Int32(511)           # YELLOW POPLAR
                    elseif bwalnut >= half
                        ift = Int32(512)           # BLACK WALNUT
                    elseif blcst >= half
                        ift = Int32(513)           # BLACK LOCUST
                    elseif redmapl >= half
                        ift = Int32(519)           # RED MAPLE-OAK
                    else
                        qoak2 = Float32(0.25) * oakhck
                        p2 = max(pstbljko,chblsc,worohk,ypworo,scruboak,swgyp,sasfprsm,mxdhwd)
                        if p2 >= qoak2
                            p2 = p2 - Float32(0.001)
                            if p2 < pstbljko
                                ift = Int32(501)   # POST-BLACKJACK OAK
                            elseif p2 < chblsc
                                ift = Int32(515)   # CHESTNUT-BLACK-SCARLET OAK
                            elseif p2 < worohk
                                ift = Int32(503)   # W.OAK-R.OAK-HICKORY
                            elseif p2 < ypworo
                                ift = Int32(506)   # Y.POPLAR-W.OAK-RED OAK
                            elseif p2 < scruboak
                                ift = Int32(514)   # SOUTHERN SCRUB OAK
                            elseif p2 < swgyp
                                ift = Int32(508)   # SWEETGUM-Y.POP
                            elseif p2 < sasfprsm
                                ift = Int32(507)   # SASSAFRAS-PERSIMMON
                            elseif mxdhwd > Float32(0)
                                ift = Int32(520)   # MIXED UPLAND HARDWOODS
                            end
                        else
                            if mxdhwd > Float32(0)
                                ift = Int32(520)   # MIXED UPLAND HARDWOODS
                            end
                        end
                    end
                elseif p1 < okgmcyp             # OAK-GUM-CYPRESS
                    if atlwcdr >= half
                        ift = Int32(606)           # ATLANTIC WHITE CEDAR
                    else
                        p2 = max(sbstrm,schchbo,swgwilo,cypwtup,ovrcupwh) - Float32(0.001)
                        if p2 > Float32(0)
                            if p2 < sbstrm
                                ift = Int32(608)   # SWEETBAY-SWAMP TUPELO-RED MAPLE
                            elseif p2 < schchbo
                                ift = Int32(601)   # SWAMP CHESTNUT-CHERRYBARK OAK
                            elseif p2 < swgwilo
                                ift = Int32(602)   # SWEETGUM-NUTTALL-WILLOW OAK
                            elseif p2 < cypwtup
                                ift = Int32(607)   # CYPRESS-WATER TUPELO
                            elseif ovrcupwh > Float32(0)
                                ift = Int32(605)   # OVERCUP OAK-WATER HICKORY
                            end
                        end
                    end
                elseif p1 < aspbrch             # ASPEN-BIRCH
                    p2 = max(aspen,blsmpop,paprbrch) - Float32(0.001)
                    if p2 > Float32(0)
                        if p2 < aspen
                            ift = Int32(901)       # ASPEN
                        elseif p2 < blsmpop
                            ift = Int32(904)       # BALSAM POPLAR
                        elseif paprbrch > Float32(0)
                            ift = Int32(902)       # PAPER BIRCH
                        end
                    end
                elseif p1 < aldrmapl            # RED ALDER-MAPLE
                    p2 = max(redaldr,bglfmpl) - Float32(0.001)
                    if p2 < redaldr
                        ift = Int32(911)           # RED ALDER
                    elseif bglfmpl > Float32(0)
                        ift = Int32(912)           # BIGLEAF MAPLE
                    end
                elseif p1 < wstoak              # WESTERN OAKS
                    p2 = max(clfbo,orgwho,blueoak,graypine,costlo,canylo,doakwdld,eoakwdld) - Float32(0.001)
                    if p2 > Float32(0)
                        if p2 < clfbo
                            ift = Int32(922)       # CALIFORNIA BLACK OAK
                        elseif p2 < orgwho
                            ift = Int32(923)       # OREGON WHITE OAK
                        elseif p2 < blueoak
                            ift = Int32(924)       # BLUE OAK
                        elseif p2 < graypine
                            ift = Int32(921)       # GRAY PINE
                        elseif p2 < costlo
                            ift = Int32(931)       # COAST LIVE OAK
                        elseif p2 < canylo
                            ift = Int32(932)       # CANYON-INTERIOR LIVE OAK
                        elseif p2 < doakwdld
                            ift = Int32(925)       # DECIDUOUS OAK WOODLAND
                        elseif eoakwdld > Float32(0)
                            ift = Int32(926)       # EVERGREEN OAK WOODLAND
                        end
                    end
                elseif p1 < tanoklrl            # TAN OAK-LAUREL
                    p2 = max(tanoak,clflrel,gntchnk) - Float32(0.001)
                    if p2 > Float32(0)
                        if p2 < tanoak
                            ift = Int32(941)       # TAN OAK
                        elseif p2 < clflrel
                            ift = Int32(942)       # CALIFORNIA LAUREL
                        elseif gntchnk > Float32(0)
                            ift = Int32(943)       # GIANT CHINKAPIN
                        end
                    end
                elseif p1 < othwhwd             # OTHER WESTERN HARDWOODS
                    p2 = max(pacmdrn,mesquitw,mtnbrshw,intmaplw,miscwhwd) - Float32(0.001)
                    if p2 > Float32(0)
                        if p2 < pacmdrn
                            ift = Int32(951)       # PACIFIC MADRONE
                        elseif p2 < mesquitw
                            ift = Int32(952)       # MESQUITE WOODLAND
                        elseif p2 < mtnbrshw
                            ift = Int32(953)       # MTN. BRUSH WOODLAND
                        elseif p2 < intmaplw
                            ift = Int32(954)       # INT. MTN. MAPLE WOODLAND
                        elseif miscwhwd > Float32(0)
                            ift = Int32(955)       # MISC. W. HARDWOODS
                        end
                    end
                elseif p1 < trophwds            # TROPICAL HARDWOODS
                    p2 = max(palm,mangrv) - Float32(0.001)
                    if p2 > Float32(0)
                        if p2 < palm
                            ift = Int32(981)       # SABAL PALM
                        elseif mangrv > Float32(0)
                            ift = Int32(982)       # MANGROVE
                        end
                    end
                elseif p1 < exhwds              # EXOTIC HARDWOODS
                    p2 = max(palonia,meluca,euclpt,othexhwd) - Float32(0.001)
                    if p2 > Float32(0)
                        if p2 < palonia
                            ift = Int32(991)       # ROYAL PAULOWNIA
                        elseif p2 < meluca
                            ift = Int32(992)       # MELALUCA
                        elseif p2 < euclpt
                            ift = Int32(993)       # EUCALYPTUS
                        elseif othexhwd > Float32(0)
                            ift = Int32(995)       # OTHER EXOTIC HARDWOODS
                        end
                    end
                elseif p1 < othhrd              # FVS OTHER HARDWOODS
                    if othhrd > Float32(0)
                        ift = Int32(997)
                    end
                elseif p1 < othspc              # FVS ALL OTHER SPECIES
                    if othspc > Float32(0)
                        ift = Int32(998)
                    end
                end
            end
        end
    end

    # --- California special case
    st.plot.forest_type = Int32(ift)
    return Int32(ift)
end
