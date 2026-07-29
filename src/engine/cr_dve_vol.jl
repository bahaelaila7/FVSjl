# =============================================================================
# cr_dve_vol.jl — CR volume, DVE method (NVEL volume/NVEL/r3d2hv.f, "R3 D2H").
#
# Volume from D2H = DBHOB²·HTTOT via per-(region+species) polynomials. Two families:
#   TIMBER (093,113, 300-122, 301-015, 301-202): ENTIRE (total cuft), GCUFT6 (to 6″),
#     GCUFT4 (to 4″), TWVOL (4-6″ topwood), INTBDFT/SCBDFT (board feet).
#   WOODLAND (060,106,800,999): DRC/D2HA-based, GCUFT4 only (no board), FCLASS branch.
# Returns VOL[1..15] (r3d2hv.f:505-546 fill): VOL1=ENTIRE, VOL4=GCUFT6(UNT1)|GCUFT4(UNT3),
# VOL7=TWVOL, VOL10=INTBDFT, VOL2=SCBDFT. `unt` = 1 for prod "01" (sawtimber), 3 for "02".
# CR uses only the 9 eq ids above (no HANN_PP — CR takes 300DVEW122, not 301/302).
#
# Math: **int→x*x/x^3, **frac→fpow, D^-n→1/D^n (openlibm ok — polynomial, no RNG).
# =============================================================================

@inline _k3h(k, h, dp) = (k * k * k * Float32(h)) / dp    # (k³·HTTOT)/DBHOB^p, p folded into dp arg

"CR DVE per-tree volume (r3d2hv.f). `d`=DBHOB, `h`=HTTOT, `drc`=root-collar dia (0⇒use d), `fclass`=1 single/else
multistem, `unt`=1|3 (prod 01|02). Returns a 15-vec; caller reads [1]/[4]/[7]/[10] like the Clark path."
function cr_dve_vol(voleq::AbstractString, d::Float32, h::Float32; unt::Int = 1,
                    drc::Float32 = 0.0f0, fclass::Int = 0, httfll::Float32 = 0.0f0, ht1prd::Float32 = 0.0f0)
    vol = zeros(Float32, 15)
    (d < 1.0f0 && drc < 1.0f0) && return vol
    reg = voleq[1:3]; spc = voleq[8:10]
    d2h = d * d * h
    # Region-2 DVE (VOLEQ 1:1=='2') → R2OLDV (volume/NVEL/r2oldv.f). CR's region-2 DVE set is all WOODLAND/
    # hardwood (065/066/069/106/814/823/998 + 475): Chojnacky INT-339 cubes TCUFT=(a+b·(D²H)^⅓+c·MSTEM)³, D²H
    # on DRC if present, MSTEM=1 single-stem (FCLASS==1) else 0. GCUFT=TCUFT (VOL1=VOL4=TCUFT, no board). The
    # `**3.` is gfortran-folded to X·X·X (integer cube — negatives cube to <0 → floored to 0, not NaN).
    if reg[1] == '2'
        h <= 0.0f0 && return vol
        dd2h = drc > 0.0f0 ? drc * drc * h : d2h
        cr3 = fpow(dd2h, 1.0f0 / 3.0f0)
        ms = fclass == 1 ? 1.0f0 : 0.0f0
        b = if spc == "065";      -0.08728f0 + 0.135420f0*cr3 - 0.019587f0*ms   # Utah juniper
            elseif spc == "066";   0.02434f0 + 0.119106f0*cr3                    # Rocky Mtn juniper
            elseif spc == "069";  -0.19321f0 + 0.136101f0*cr3 + 0.038187f0*ms   # oneseed juniper
            elseif spc == "106";  -0.20296f0 + 0.150283f0*cr3 + 0.054178f0*ms   # pinyon pine (R2)
            elseif spc == "475";  -0.13363f0 + 0.128222f0*cr3 + 0.080208f0*ms   # mountain mahogany
            elseif spc == "814";  -0.13600f0 + 0.145743f0*cr3                    # Gambel oak
            elseif spc == "823";   0.12853f0 + 0.105885f0*cr3                    # bur oak
            elseif spc == "998";  -0.13822f0 + 0.121850f0*cr3                    # other hardwoods
            else; return vol                                                     # region-2 conifer — not in CR forest table
            end
        tcuft = b * b * b
        spc == "475" && drc < 3.0f0 && d < 3.0f0 && (tcuft = 0.1f0)              # FIA <3" floor (r2oldv.f:447)
        tcuft < 0.0f0 && (tcuft = 0.0f0)
        vol[1] = tcuft; vol[4] = tcuft                                           # VOL1 total, VOL4 gross-merch
        return vol
    end
    entire = 0.0f0; gcuft6 = 0.0f0; gcuft4 = 0.0f0; twvol = 0.0f0
    intbdft = 0.0f0; scbdft = 0.0f0; um4 = 0.0f0; um6 = 0.0f0
    d1 = 1.0f0 / d; d15 = fpow(d, 1.5f0)                       # DBHOB^-1, DBHOB^1.5

    if spc == "093"                                            # Engelmann spruce, etc.
        entire = 0.225466f0 + 0.002170f0 * d2h
        um6 = -0.2664752f0 + 0.006129f0 * ((6f0^3 * h) / d15) + 0.007431f0 * d * d
        gcuft6 = max(entire - um6, 0f0)
        intbdft = gcuft6 * (5.987363f0 - 9.847918f0 * d1 - (-300.812808f0) * d1 * d1 - 2855.342454f0 * d1 * d1 * d1)
        scbdft = intbdft * (0.878454f0 - 15.998458f0 * d1 * d1)
        um4 = -0.2664752f0 + 0.006129f0 * ((4f0^3 * h) / d15) + 0.007431f0 * d * d
        gcuft4 = entire - um4; twvol = um6 - um4
    elseif spc == "113"                                        # SW white/bristlecone/limber/foxtail pine
        entire = 0.160889f0 + 0.002032f0 * d2h
        um6 = -0.213005f0 + 0.004912f0 * ((6f0^3 * h) / d15) + 0.006061f0 * d * d
        gcuft6 = max(entire - um6, 0f0)
        intbdft = gcuft6 * (6.691967f0 - 7.520114f0 * d1 - 216.348366f0 * d1 * d1)
        scbdft = intbdft * (1.006086f0 - 2.384660f0 * d1)
        um4 = -0.213005f0 + 0.004912f0 * ((4f0^3 * h) / d15) + 0.006061f0 * d * d
        gcuft4 = entire - um4; twvol = um6 - um4
    elseif spc == "122" && reg == "300"                        # ponderosa pine (R3 self-contained)
        if d2h <= 31629.91964f0
            scbdft = -1.786f0 + 0.00098814f0 * d2h
        else
            htf = httfll <= 0f0 ? 25.0f0 : httfll
            scbdft = -52.897f0 + 0.12826f0 * htf + 0.0017678f0 * d2h + 879120.0f0 / d2h
        end
        scbdft *= 10.0f0
        gcuft6 = d2h <= 33590.92207f0 ? (-1.7751f0 + 0.0018897f0 * d2h) : (-13.542f0 + 0.00224f0 * d2h)
        entire = 0.081072f0 + 0.001984f0 * d2h
        um4 = -0.125349f0 + 0.003604f0 * ((4f0^3 * h) / d15) + 0.005406f0 * d * d
        gcuft4 = entire - um4
        um6 = -0.125349f0 + 0.003604f0 * ((6f0^3 * h) / d15) + 0.005406f0 * d * d
        twvol = um6 - um4
        # NOTE: 122+300 leaves INTBDFT=0 (r3d2hv fills VOL(10)=INTBDFT=0; board = SCBDFT only)
    elseif spc == "015" && reg == "301"                        # white fir (R3-01)
        entire = 0.210904f0 + 0.001840f0 * d2h
        um6 = -0.182700f0 + 0.001248f0 * ((6f0^3 * h) / d) + 0.006245f0 * d * d       # NOTE: /D (^1.0), not /D^1.5
        gcuft6 = max(entire - um6, 0f0)
        intbdft = gcuft6 * (6.246875f0 - 7.019940f0 * d1 - 201.958728f0 * d1 * d1)
        scbdft = intbdft * (1.0f0 - 1.888144f0 * d1 - 8.851449f0 * d1 * d1)
        um4 = -0.182700f0 + 0.001248f0 * ((4f0^3 * h) / d) + 0.006245f0 * d * d
        gcuft4 = entire - um4; twvol = um6 - um4
    elseif spc == "202" && reg == "301"                        # Douglas-fir (R3-01)
        entire = 0.438374f0 + 0.001756f0 * d2h
        um6 = -0.083149f0 + 0.001219f0 * ((6f0^3 * h) / d) + 0.005417f0 * d * d       # /D (^1.0)
        gcuft6 = max(entire - um6, 0f0)
        intbdft = gcuft6 * (6.587353f0 - 0.892716f0 * d1 - 243.514909f0 * d1 * d1)
        scbdft = intbdft * (1.000897f0 - 4.100072f0 * fpow(d, -1.177748f0))
        um4 = -0.083149f0 + 0.001219f0 * ((4f0^3 * h) / d) + 0.005417f0 * d * d
        gcuft4 = entire - um4; twvol = um6 - um4
        scbdft *= 0.932f0
    elseif spc == "060" || spc == "106" || spc == "800" || spc == "999"   # WOODLAND (VOLEQU(2:3)="00" ⇒ INT-391)
        if d > 3.0f0 || drc > 3.0f0
            drc > 0f0 && (d2h = drc * drc * h)
            d2ha = d2h / 1000.0f0
            if spc == "060"
                bp = 6.0f0
                gcuft4 = fclass != 1 ?
                    (d2ha <= bp ? -0.129f0 + 2.0255f0 * d2ha + 0.1011f0 * d2ha^2 : 10.786f0 + 2.0255f0 * d2ha - 43.663f0 / d2ha) :
                    (d2ha <= bp ? -0.032f0 + 2.1076f0 * d2ha + 0.1454f0 * d2ha^2 : 15.675f0 + 2.1076f0 * d2ha - 62.827f0 / d2ha)
            elseif spc == "106"                                # no FCLASS branch, bp=3
                gcuft4 = d2ha <= 3.0f0 ? -0.060f0 + 2.5139f0 * d2ha + 0.1466f0 * d2ha^2 :
                                          3.898f0 + 2.5139f0 * d2ha - 7.917f0 / d2ha
            elseif spc == "800"
                bp = 4.0f0
                gcuft4 = fclass != 1 ?
                    (d2ha <= bp ? -0.028f0 + 1.9545f0 * d2ha + 0.1400f0 * d2ha^2 : 6.691f0 + 1.9545f0 * d2ha - 17.918f0 / d2ha) :
                    (d2ha <= bp ? -0.068f0 + 2.4048f0 * d2ha + 0.1383f0 * d2ha^2 : 6.571f0 + 2.4048f0 * d2ha - 17.704f0 / d2ha)
            else                                               # 999 (mesquite/default), bp=2
                gcuft4 = fclass != 1 ?
                    (d2ha <= 2.0f0 ? 0.020f0 + 1.8972f0 * d2ha + 0.5756f0 * d2ha^2 : 6.927f0 + 1.8972f0 * d2ha - 9.210f0 / d2ha) :
                    (d2ha <= 2.0f0 ? -0.043f0 + 2.3378f0 * d2ha + 0.8024f0 * d2ha^2 : 9.586f0 + 2.3378f0 * d2ha - 12.839f0 / d2ha)
            end
        else
            gcuft4 = spc == "800" ? 0.1f0 : 0.0f0              # 800 floors to 0.1 (FIA match), others 0
        end
        entire = gcuft4
        unt == 1 && (gcuft6 = gcuft4)
    else                                                       # Smalian fallback (r3d2hv.f:495)
        gcuft4 = (d * d * 0.005454f0 + 16.0f0 * 0.005454f0) / 2.0f0 * ht1prd
        entire = gcuft4
        unt == 1 && (gcuft6 = gcuft4)
    end

    d < 9.0f0 && (scbdft = 0f0; intbdft = 0f0)                 # r3d2hv.f:508 DBH<9 ⇒ no board
    vol[1] = entire
    if unt == 1
        vol[2] = scbdft; vol[10] = intbdft; vol[4] = gcuft6; vol[7] = twvol
    elseif unt == 3
        vol[4] = gcuft4; vol[6] = vol[4] / 79.0f0
    end
    vol[15] = um4
    @inbounds for k in (1, 2, 4, 6, 7, 10, 15)
        vol[k] < 0f0 && (vol[k] = 0f0)
    end
    return vol
end

# compute_volumes_cr! — CR per-tree volume (vols.f → NATCRS method dispatch). Dispatches each species' eq id
# (VEQNNC method chars 4:6, or "NVB" prefix) to the ported NVEL method, then applies the fvsvol.f return
# mapping (fvsvol.f:510-529): TCF=TVOL(1)≥0; MCF=TVOL(4)+TVOL(7) gated D≥DBHMIN; SCF=TVOL(4) ONLY for
# region 8/9/FIA-NVB (so 0 for CR region 2/3); BF gated D≥BFMIND. Merch specs are the CR sitset defaults
# (cr/sitset.f:520-555, by IMODTY): DBHMIN/TOPD/STUMP. DVE/NVB/FW2 all ported (cr_dve_vol/cr_nvb_vol/cr_fw2_vol);
# woodland FCLASS=0 (FIA multi-stem) fixed. Residuals: FW2 ~1.5% precision + size-dependent NVB/DVE taper on a few
# conifer stands (partly entangled with the dense self-thinning/DGSCOR tail); board-foot for NVB still partial.
function compute_volumes_cr!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; c = s.control; sd = s.coef.species
    scfmin = c.sp_scf_dbhmin
    imodty = Int(s.plot.model_type)
    # CR volume merch standards (cr/sitset.f): IMODTY 3 (Black-Hills PP) vs the rest.
    is3 = imodty == 3
    dbhmin = is3 ? 9f0 : 5f0
    topd   = is3 ? 6f0 : 4f0
    stump  = 1f0
    # Board-foot min DBH (cr/sitset.f): IMODTY 3 → 9; else 7 (IFOR<IGFOR=13) / 9. Board top DOB = 6.
    ifor   = Int(s.plot.forest_idx)
    bfmind = is3 ? 9f0 : ((ifor > 0 && ifor < 13) ? 7f0 : 9f0)
    iregn  = Int(s.plot.user_forest_code) ÷ 100    # stand region (MRULES keys merch bucking on REGN)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]
        mdl = length(eq) >= 6 ? eq[4:6] : "   "
        nvb = startswith(eq, "NVB")
        v = if mdl == "DVE"
            cr_dve_vol(eq, d, h; unt = d >= scfmin[sp] ? 1 : 3)
        elseif nvb
            bark = cr_bratio(sd, sp, d, imodty)
            cr_nvb_vol(eq, d, h; bark = bark, topd = topd, stump = stump, iregn = iregn)   # TCF+MCF+board
        elseif mdl == "FW2"
            cr_fw2_vol(eq, d, h; bark = cr_bratio(sd, sp, d, imodty), topd = topd, stump = stump, iregn = iregn)   # TCF+MCF+board
        else
            zeros(Float32, 15)
        end
        tcf = max(v[1], 0f0)
        mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
        scf = 0f0                              # CR is region 2/3: fvsvol.f sets SCF only for region 8/9
        # BdFt = BBFV = TVOL(2) Scribner for CR (METHB=6≠9), gated D≥BFMIND. DVE/NVB/FW2 all fill VOL(2).
        bf  = d >= bfmind ? v[2] : 0f0
        t.cuft_vol[i] = tcf; t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i] = scf; t.bdft_vol[i] = bf
    end
    return s
end
