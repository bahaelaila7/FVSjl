# =============================================================================
# io/errgro.jl — ERRGRO error/warning messages (errgro.f) and the FVS_Error DBS table (dbserror.f)
#
# errgro.f builds CMSG (CHARACTER*256) for each error number and CALLs DBSERROR(NPLT,CMSG), which writes one
# FVS_Error row (CaseID, StandID, Message) — the message stored as the full 256-char blank-padded CMSG. DBSERROR
# forces the case (DBSCASE(1)), so the row is written whether or not any other output table was requested.
# The messages are collected per stand (`s.control.error_msgs`) and written with the stand's other DBS tables.
# =============================================================================

# keywds.f BLOCK DATA TABLE (identical in every variant build): the base keywords FNDKEY accepts. Anything else read
# by the base keyword reader (initre.f:140-170) is ERRGRO(.TRUE.,1) "FVS01 INVALID KEYWORD".
const FVS_BASE_KEYWORDS = Set([
    "ADDFILE", "AGPLABEL", "ALSOTRY", "ANIN", "ATRTLIST", "BAIMULT", "BAMAX", "BFDEFECT", "BFFDLN", "BFVOLEQU",
    "BFVOLUME", "BMIN", "BRUST", "CALBSTAT", "CCADJ", "CFVOLEQU", "CHEAPO", "CLIMATE", "CLOSE", "COMPRESS",
    "COMPUTE", "COVER", "CRNMULT", "CUTEFF", "CUTLIST", "CWEQN", "CYCLEAT", "DATABASE", "DATASCRN", "DEBUG",
    "DEFECT", "DELOTAB", "DESIGN", "DFB", "DFTM", "DGSTDEV", "ECHO", "ECHOSUM", "ECON", "ENDFILE", "ENDIF",
    "ESTAB", "FERTILIZ", "FIAVBC", "FIXCW", "FIXDG", "FIXHTG", "FIXMORT", "FMIN", "FVSSTAND", "GROWTH",
    "HTGMULT", "HTGSTOP", "IF", "INVYEAR", "LOCATE", "MANAGED", "MCDEFECT", "MCFDLN", "MGMTID", "MINHARV",
    "MISTOE", "MODTYPE", "MORTMSB", "MORTMULT", "MPB", "NOAUTOES", "NOCALIB", "NODEBUG", "NOECHO", "NOHTDREG",
    "NOSCREEN", "NOSUM", "NOTREES", "NOTRIPLE", "NUMCYCLE", "NUMTRIP", "OPEN", "ORGANON", "POINTREF", "PRMFROST",
    "PROCESS", "PRUNE", "PTGROUP", "RANNSEED", "RDIN", "READCORD", "READCORH", "READCORR", "REGDMULT", "REGHMULT",
    "RESETAGE", "REUSCORD", "REUSCORH", "REUSCORR", "REWIND", "RRIN", "SCREEN", "SDICALC", "SDIMAX", "SERLCORR",
    "SETPTHIN", "SETSITE", "SITECODE", "SPCODES", "SPECPREF", "SPGROUP", "SPLABEL", "SPLEAVE", "STANDCN", "STATS",
    "STDIDENT", "STDINFO", "STRCLASS", "SVS", "TCONDMLT", "TFIXAREA", "THEN", "THINABA", "THINATA", "THINAUTO",
    "THINBBA", "THINBTA", "THINCC", "THINDBH", "THINHT", "THINMIST", "THINPRSC", "THINPT", "THINQFA", "THINRDEN",
    "THINRDSL", "THINSDI", "TIMEINT", "TOPKILL", "TREEDATA", "TREEFMT", "TREELIST", "TREESZCP", "VOLEQNUM",
    "VOLUME", "WSBW", "YARDLOSS", "XSALVAGE"])

# Extension IN keywords whose variant build links the ex*.f stub that raises ERRGRO(.TRUE.,11) (exclim.f CLIN,
# excov.f CVIN, exmist.f MISIN, exrd.f RDIN — per bin/FVS<v>_buildDir). The stub reads nothing, so the block's
# sub-keywords and its END reach the base reader (FVS01 each). (exbm/exbrus/exdfb/exdftm/exmpb are stubs in these
# builds too, but jl runs those extensions — not switched here.)
ext_stub_keywords(::Union{Southern,CentralStates,LakeStates,Northeast,Ontario}) = ("CLIMATE", "COVER", "MISTOE", "RDIN")
ext_stub_keywords(::SoutheastAlaska) = ("CLIMATE", "RDIN")
ext_stub_keywords(::BritishColumbia) = ("COVER",)
ext_stub_keywords(::Union{CentralCalifornia,OregonCoast,Olympic}) = ("RDIN",)
ext_stub_keywords(::AbstractVariant) = ()

# Fortran Iw edit descriptor (right-justified; w asterisks on overflow).
_fI(n::Integer, w::Int) = (t = string(n); length(t) > w ? "*"^w : lpad(t, w))

const _ERRGRO_TEXT = Dict{Int,String}(
    3 => "FVS03 WARNING:  FOREST CODE INDICATES THE GEOGRAPHIC LOCATION IS OUTSIDE THE RANGE OF THE MODEL.  " *
         "DEFAULT CODE IS USED.",
    11 => "FVS11 ERROR:  REQUESTED EXTENSION IS NOT PART OF THIS PROGRAM.",
    14 => "FVS14 WARNING:  HABITAT/PLANT ASSOCIATION/ECOREGION CODE WAS NOT RECOGNIZED; HABITAT/PLANT " *
          "ASSOCIATION/ECOREGION SET TO DEFAULT CODE.",
    32 => "FVS32 WARNING:  PV REFERENCE CODE WAS NOT RECOGNIZED; HABITAT/PLANT ASSOCIATION/ECOREGION SET TO " *
          "DEFAULT CODE.",
    33 => "FVS33 WARNING:  PV CODE WAS NOT RECOGNIZED; HABITAT/PLANT ASSOCIATION/ECOREGION SET TO DEFAULT CODE.",
    34 => "FVS34 WARNING:  PV CODE/PV REFERENCE CODE COMBINATION WAS NOT RECOGNIZED; HABITAT/PLANT " *
          "ASSOCIATION/ECOREGION SET TO DEFAULT CODE.",
    40 => "FVS40 WARNING:  TREE RECORD REPRESENTING GREATER THAN 1000 TPA ENCOUNTERED. MAY CAUSE MATHEMATICAL ERRORS.",
    41 => "FVS41 WARNING:  INITIAL STAND STOCKING IS MORE THAN 5% ABOVE LIMIT, SDI MAXIMUM RESET.")

"""
    errgro!(s, ierrn; irecnt, irec1, irecrd)

ERRGRO(.TRUE.,IERRN): record the CMSG text (errgro.f FORMAT 1111/1811/2111/…) for the stand's FVS_Error table.
FVS01 carries the keyword record count IRECNT (I4); FVS08 the projectable/total tree-record counts (I2/I4) and the
stand id (A26).
"""
function errgro!(s::StandState, ierrn::Integer; irecnt::Integer = 0, irec1::Integer = 0, irecrd::Integer = 0)
    msg = if ierrn == 1
        "FVS01 ERROR:  INVALID KEYWORD WAS SPECIFIED.  RECORDS READ=" * _fI(irecnt, 4)            # errgro.f:1111
    elseif ierrn == 8
        "FVS08 WARNING:  TOO FEW PROJECTABLE TREE RECORDS.  PROJECTABLE RECORDS:" * _fI(irec1, 2) *   # errgro.f:1811
            "; TREE RECORDS:" * _fI(irecrd, 4) * "; STAND ID: " * rpad(first(s.plot.stand_id, 26), 26)
    else
        _ERRGRO_TEXT[Int(ierrn)]
    end
    push!(s.control.error_msgs, msg)
    return s
end

# --- HABTYP error branches (ERRGRO 14/32/33/34) ---------------------------------------------------------------
# The (PVCODE, PVREF) crosswalk tables as the full row lists (PVREF1/PVREF6 set LPVCOD/LPVREF from ANY row, including
# rows whose HABPVR is blank — so BM reads data/bluemountains/pvref6_all.csv, all 3144 bm/pvref6.f DATA rows).
_pvtable_rows(raw::AbstractString) = [(String(f[1]), String(f[2])) for f in (split(l) for l in eachline(IOBuffer(raw)))
                                      if length(f) >= 2]
const _PVSETS = Dict{Symbol,Any}()
function _pv_sets(key::Symbol)
    get!(_PVSETS, key) do
        rows = key === :ie ? _pvtable_rows(_IE_PVREF1_RAW) :
               key === :em ? _pvtable_rows(_EM_PVREF1_RAW) :
               [(String(strip(f[1])), String(strip(f[2]))) for f in
                (split(l, ',') for l in readlines(joinpath(BM_DATADIR, "pvref6_all.csv"))[2:end]) if length(f) >= 2]
        (Set(first.(rows)), Set(last.(rows)))
    end
end

# PVREF1/PVREF6 flag scan for a NON-full-match (a full match returns LPVCOD=LPVREF=.TRUE. with a resolved code).
function _pvref_flags(key::Symbol, code::AbstractString, cpvref::AbstractString)
    codes, refs = _pv_sets(key)
    return (code in codes, cpvref in refs)
end

# habtyp.f KARD2 preprocessing (IE/EM): ADJUSTL; a ≤2-char all-digit code gets a leading '0' (habtyp.f:84-94).
function _habtyp_kard2(pv::AbstractString)
    k = String(strip(pv))
    (length(k) == 2 && all(isdigit, k)) && (k = "0" * k)
    return k
end
# PVREF1 drops everything from the last '.' (ie/em pvref1.f DO I=10,1,-1).
_pvref1_key(k::AbstractString) = (i = findlast('.', k); i === nothing ? k : k[1:i-1])

"""
    habtyp_errors!(s, kard2, cpvref, kodtyp)

The ERRGRO calls HABTYP makes for this habitat input (KARD2 = `kard2`, CPVREF = `cpvref` (blank = ""), and the
caller's KODTYP = IFIX(ARRAY2) — the numeric value of the code, 0 when it is not an integer), for the
variants whose habtyp.f branches are ported here (IE ie/habtyp.f:84-135, EM em/habtyp.f:84-121, BM bm/habtyp.f:
90-163). The habitat RESOLUTION itself is done by the FIA reader / site setup; this only mirrors the error branches.
"""
habtyp_errors!(::StandState, ::AbstractString, ::AbstractString, ::Integer) = nothing
function habtyp_errors!(s::StandState{InlandEmpire}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    k = _habtyp_kard2(pv); lpvxxx = false; kodtyp = Int(kodtyp)
    if !isempty(cpvref)
        key = _pvref1_key(k)
        if haskey(_IE_PVREF1, (key, cpvref))                  # full match ⇒ LPVCOD=LPVREF=.TRUE., KARD2=HABPVR
            h = _IE_PVREF1[(key, cpvref)]
            kodtyp = all(isdigit, h) ? parse(Int, h) : ie_pa_habitat_code(h)
        else
            lc, lr = _pvref_flags(:ie, key, cpvref)
            kodtyp = 0
            if lc && lr
                errgro!(s, 34); lpvxxx = true                  # KODTYP out of range and KARD2 blank
            elseif !lc && !lr
                errgro!(s, 33); errgro!(s, 32); lpvxxx = true
            elseif lc
                errgro!(s, 32); lpvxxx = true
            else
                errgro!(s, 33); lpvxxx = true
            end
        end
    elseif length(k) == 6
        kodtyp = ie_pa_habitat_code(k)                          # HBDECD on the plant-association code
    end
    (kodtyp < 10 || kodtyp > 999) && !lpvxxx && errgro!(s, 14)  # habtyp.f:133-135
    return nothing
end
function _em_jtype_ok(kodtyp::Int)
    (kodtyp == 0 || kodtyp < EM_JTYPE[1]) && return false
    return kodtyp < EM_JTYPE[end]                              # DO 20 J=1,118 … falls through ⇒ label 21
end
function habtyp_errors!(s::StandState{EasternMontana}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    k = _habtyp_kard2(pv); kodtyp = Int(kodtyp)
    if !isempty(cpvref)
        key = _pvref1_key(k)
        kodtyp = get(_EM_PVREF1, (key, cpvref), 0)
        if kodtyp <= 0
            lc, lr = _pvref_flags(:em, key, cpvref)
            if lc && lr
                errgro!(s, 34); return nothing                  # em/habtyp.f:94-97 (GOTO 21 with LPVXXX)
            elseif !lc && !lr
                errgro!(s, 33); errgro!(s, 32); return nothing
            elseif lc
                errgro!(s, 32); return nothing
            else
                errgro!(s, 33); return nothing
            end
        end
    end
    _em_jtype_ok(kodtyp) || errgro!(s, 14)                   # em/habtyp.f:112-121
    return nothing
end
function habtyp_errors!(s::StandState{BlueMountains}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    k = String(strip(pv))
    if !isempty(cpvref)
        h = bm_pvref6(k, cpvref)
        if isempty(h)
            lc, lr = _pvref_flags(:bm, k, cpvref)
            if lc && lr
                errgro!(s, 34)                                  # full match absent (or blank HABPVR) ⇒ KARD2 blank
            elseif !lc && !lr
                errgro!(s, 33); errgro!(s, 32)
            elseif lc
                errgro!(s, 32)
            else
                errgro!(s, 33)
            end
            return nothing
        end
        k = h
    end
    # HBDECD against PCOML, then the sequence-number fallback IHB=INT(ARRAY2) (bm/habtyp.f:136-158)
    findfirst(==(k), BM_PCOML) !== nothing && return nothing
    ihb = isempty(cpvref) ? Int(kodtyp) : 0                    # IHB = INT(ARRAY2); PVREF6 zeroes ARRAY2
    (1 <= ihb <= length(BM_PCOML)) && return nothing
    errgro!(s, 14)
    return nothing
end

function habtyp_errors!(s::StandState{CentralIdaho}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    _, errs, _ = ci_habitat_kodtyp(pv, cpvref, kodtyp)                       # ci/habtyp.f:54-75
    foreach(e -> errgro!(s, e), errs)
    return nothing
end

function habtyp_errors!(s::Union{StandState{Teton},StandState{Utah}}, pv::AbstractString, cpvref::AbstractString,
                        kodtyp::Integer)
    _, errs = r4_habitat_itype(pv, cpvref, kodtyp)                     # tt|ut/habtyp.f:152-228
    foreach(e -> errgro!(s, e), errs)
    return nothing
end

# WC/PN/EC habtyp.f (Region-6 plant associations, identical shape; defaults ITYPE 52 / 40 / 114): with a reference code
# PVREF6 crosswalks (PV_CODE, ref) → HABPVR (blank / no full match ⇒ FVS34 / FVS33+FVS32 / FVS32 / FVS33 and LPVXXX);
# then HBDECD matches the code against PCOML, else the IHB = IFIX(ARRAY2) sequence number in 1..NPA; unresolved and
# not LPVXXX ⇒ FVS14. The PVREF6 rows (incl. blank HABPVR) come from each build's pvref6.f (data/<v>/pvref6_all.csv).
const _R6PV = Dict{Symbol,Any}()
function _r6pv(key::Symbol, dir::AbstractString)
    get!(_R6PV, key) do
        m = Dict{Tuple{String,String},String}(); codes = Set{String}(); refs = Set{String}()
        for l in readlines(joinpath(dir, "pvref6_all.csv"))[2:end]
            f = split(l, ','); length(f) >= 3 || continue
            c = String(strip(f[1])); r = String(strip(f[2])); h = String(strip(f[3]))
            push!(codes, c); push!(refs, r); haskey(m, (c, r)) || (m[(c, r)] = h)
        end
        (m, codes, refs)
    end
end
function _r6_habtyp_errs(pv::AbstractString, cpvref::AbstractString, kodtyp::Integer, pcoml, key::Symbol, dir)
    errs = Int[]; k = String(strip(pv)); r = String(strip(cpvref)); a2 = Int(kodtyp)
    if !isempty(r)
        m, codes, refs = _r6pv(key, dir)
        h = get(m, (k, r), nothing)
        if h === nothing
            lc = k in codes; lr = r in refs; k = ""
        else
            lc = true; lr = true; k = h; a2 = 0                             # full match: KARD2=HABPVR, ARRAY2=0
        end
        if lc && lr && isempty(k); push!(errs, 34); return errs
        elseif !lc && !lr;         push!(errs, 33, 32); return errs
        elseif !lr && lc;          push!(errs, 32); return errs
        elseif !lc && lr;          push!(errs, 33); return errs
        end
    end
    npa = length(pcoml)
    ihb = a2
    resolved = if 0 <= ihb <= npa                                          # HBDECD
        ihb > 0 ? true : (!isempty(k) && k[1] != '0' && findfirst(==(uppercase(first(k, 8))), pcoml) !== nothing)
    else
        false
    end
    resolved || (1 <= a2 <= npa) || push!(errs, 14)                        # label 20 sequence number, else 14
    return errs
end
for (V, pc, key, dir) in ((:WestCascades, :WC_PCOML, :wc, :WC_DATADIR), (:PacificNorthwest, :PN_PCOML, :pn, :PN_DATADIR),
                          (:EastCascades, :EC_PCOML, :ec, :EC_DATADIR))
    @eval function habtyp_errors!(s::StandState{$V}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
        foreach(e -> errgro!(s, e), _r6_habtyp_errs(pv, cpvref, kodtyp, $pc, $(QuoteNode(key)), $dir))
        return nothing
    end
end

# Region-5 / Region-5+6 habtyp.f (WS; NC/SO/CA). PVREF5 zeroes ARRAY2 on every call, PVREF6 only on a full match;
# both blank KARD2 unless a full (PV_CODE, ref) row matches. Tables: data/<v>/pvref5_all.csv, pvref6_all.csv,
# r5habt.csv (generated from each build's pvref5.f / pvref6.f / habtyp.f R5HABT).
const _PVALL = Dict{String,Any}()
function _pvall(path::AbstractString)
    get!(_PVALL, path) do
        m = Dict{Tuple{String,String},String}(); codes = Set{String}(); refs = Set{String}()
        for l in readlines(path)[2:end]
            f = split(l, ','); length(f) >= 3 || continue
            c = String(strip(f[1])); r = String(strip(f[2])); h = String(strip(f[3]))
            push!(codes, c); push!(refs, r); haskey(m, (c, r)) || (m[(c, r)] = h)
        end
        (m, codes, refs)
    end
end
const _R5HABT = Dict{String,Vector{String}}()
_r5habt(dir::AbstractString) = get!(() -> [String(strip(l)) for l in readlines(joinpath(dir, "r5habt.csv"))[2:end]],
                                    _R5HABT, dir)
_crdecd(k::AbstractString, tab) = isempty(k) ? 0 : something(findfirst(==(uppercase(first(k, 8))), tab), 0)

# ws/habtyp.f (R5HABT only; the TT/UT shape with PVREF5): a ≥2-char code with a reference code → PVREF5 (FVS34/33/32),
# CRDECD vs R5HABT, then IHB = IFIX(ARRAY2) ≤ NR5; unresolved and not LPVXXX ⇒ FVS14.
function _ws_habtyp_errs(pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    errs = Int[]; lpvxxx = false; k0 = String(strip(pv)); k = k0; r = String(strip(cpvref)); a2 = Int(kodtyp); kod = 0
    r5 = _r5habt(WS_DATADIR)
    if length(k) >= 2
        if !isempty(r)
            m, codes, refs = _pvall(joinpath(WS_DATADIR, "pvref5_all.csv"))
            h = get(m, (k, r), nothing); a2 = 0
            if h === nothing
                lc = k0 in codes; lr = r in refs; k = ""
            else
                lc = true; lr = true; k = String(strip(h))
            end
            if lc && lr && isempty(k);  push!(errs, 34); lpvxxx = true
            elseif !lc && !lr;          push!(errs, 33, 32); lpvxxx = true
            elseif !lr && lc;           push!(errs, 32); lpvxxx = true
            elseif !lc && lr;           push!(errs, 33); lpvxxx = true
            end
        end
        kod = _crdecd(k, r5)
    end
    kod == 0 && a2 <= length(r5) && (kod = a2)
    kod <= 0 && !lpvxxx && push!(errs, 14)
    return errs
end

# nc/so/ca habtyp.f (R5HABT on the Region-5 forests, PCOML on the Region-6 ones — IR5/IR6 by the FORKOD-mapped
# KODFOR). The PVREF branches only fire while KODTYP ≤ 0; an unresolved code (the R5HABT/PCOML text loop and the
# IHB sequence number both fail) is FVS14 unless LPVXXX.
function _r56_habtyp_errs(pv::AbstractString, cpvref::AbstractString, kodtyp::Integer, kodfor::Int, ir5::Bool,
                          ir6::Bool, dir::AbstractString, pcoml)
    errs = Int[]; lpvxxx = false; k0 = String(strip(pv)); k = k0; r = String(strip(cpvref))
    a2 = Int(kodtyp); kod = Int(kodtyp)
    r5 = _r5habt(dir); nr5 = length(r5); npa = length(pcoml)
    if !isempty(r)
        lc = false; lr = false
        if ir5 || ir6
            m, codes, refs = _pvall(joinpath(dir, ir5 ? "pvref5_all.csv" : "pvref6_all.csv"))
            kod = 0; ir5 && (a2 = 0)
            h = get(m, (k0, r), nothing)
            if h === nothing
                lc = k0 in codes; lr = r in refs; k = ""
            else
                lc = true; lr = true; k = String(strip(h)); a2 = 0
            end
        end
        if kod <= 0
            if lc && lr && isempty(k);  push!(errs, 34); lpvxxx = true
            elseif !lc && !lr;          push!(errs, 33, 32); lpvxxx = true
            elseif !lr && lc;           push!(errs, 32); lpvxxx = true
            elseif !lc && lr;           push!(errs, 33); lpvxxx = true
            end
        end
        lpvxxx && return errs                                           # GO TO 300 with KODTYP 0, LPVXXX
    end
    if ir5
        kod = _crdecd(k, r5)
    elseif ir6                                                          # HBDECD vs PCOML
        ihb = a2
        kod = (0 <= ihb <= npa) ? (ihb > 0 ? ihb + nr5 :
              ((isempty(k) || k[1] == '0') ? 0 : (i = findfirst(==(uppercase(first(k, 8))), pcoml); i === nothing ? 0 : i + nr5))) : 0
    end
    if kodfor == 0 || kod == 0
        kk = first(k, 8); found = false
        if kodfor == 0 || ir5
            i = findfirst(==(kk), r5); i !== nothing && (kod = i; found = true)
        end
        if !found
            i = findfirst(==(kk), pcoml); i !== nothing && (kod = i + nr5; found = true)
        end
        if !found
            kod = (a2 <= nr5 && ir5) ? a2 : (nr5 < a2 <= nr5 + npa && ir6) ? a2 : 0
        end
    end
    kod == 0 && !lpvxxx && push!(errs, 14)
    return errs
end
_mapped_kodfor(s::StandState, forkod!) = (q = deepcopy(s.plot); forkod!(q); Int(q.user_forest_code))

function habtyp_errors!(s::StandState{WestSierra}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    foreach(e -> errgro!(s, e), _ws_habtyp_errs(pv, cpvref, kodtyp))
    return nothing
end
function habtyp_errors!(s::StandState{Klamath}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    kf = _mapped_kodfor(s, nc_forkod!)
    ir5 = kf > 0 && (kf < 600 || kf == 705 || kf == 800); ir6 = !ir5 && kf > 0 && (kf == 611 || kf == 712)
    foreach(e -> errgro!(s, e), _r56_habtyp_errs(pv, cpvref, kodtyp, kf, ir5, ir6, NC_DATADIR, NC_PCOML))
    return nothing
end
function habtyp_errors!(s::StandState{SouthCentralOregon}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    kf = _mapped_kodfor(s, so_forkod!)
    ir5 = kf in (505, 506, 509, 511, 701, 514); ir6 = kf in (601, 602, 620, 799)
    foreach(e -> errgro!(s, e), _r56_habtyp_errs(pv, cpvref, kodtyp, kf, ir5, ir6, SO_DATADIR, SO_PCOML))
    return nothing
end
function habtyp_errors!(s::StandState{CentralCalifornia}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    kf = _mapped_kodfor(s, ca_forkod!)
    ir5 = 0 < kf < 600; ir6 = kf >= 600
    foreach(e -> errgro!(s, e), _r56_habtyp_errs(pv, cpvref, kodtyp, kf, ir5, ir6, CA_DATADIR, CA_PCOML))
    return nothing
end

# --- FVS_Error DBS table (dbserror.f) ------------------------------------------------------------------------------
const _FVS_ERROR_CREATE = """
CREATE TABLE IF NOT EXISTS FVS_Error(
  CaseID text not null, StandID text not null, Message text)"""

"""
    write_dbs_error!(dbpath, caseid, standid, msgs) -> dbpath

One FVS_Error row per ERRGRO message (dbserror.f), Message = the 256-char blank-padded CMSG.
"""
function write_dbs_error!(dbpath, caseid::AbstractString, standid::AbstractString, msgs::AbstractVector{<:AbstractString})
    isempty(msgs) && return dbpath
    db = SQLite.DB(dbpath)
    try
        _ensure_table!(db, _FVS_ERROR_CREATE)
        stmt = DBInterface.prepare(db, "INSERT INTO FVS_Error VALUES (?,?,?)")
        for m in msgs
            DBInterface.execute(stmt, (caseid, standid, rpad(first(m, 256), 256)))
        end
    finally
        SQLite.close(db)
    end
    return dbpath
end
