# ls/habtyp.f LSNPC (63 plant associations) and ls/pvref9.f (PVCODE, PVREF, HABPVR) — generated from the LS buildDir.
const LS_LSNPC = String["FDN12", "FDN22", "FDN32", "FDN33", "FDN43", "FDC12", "FDC23", "FDC24", "FDC25", "FDC34", "MHN35", "MHN44", "MHN45", "MHN46", "MHN47", "MHC26", "MHC36", "MHC37", "MHC47", "FFN57", "FFN67", "WFN53", "WFN55", "WFN64", "WFS57", "WFW54", "FPN62", "FPN63", "FPN71", "FPN72", "FPN73", "FPN81", "FPN82", "FPS63", "FPW63", "APN80", "APN81", "APN90", "APN91", "CTN11", "CTN12", "CTN24", "CTN32", "CTN42", "CTU22", "RON12", "RON23", "LKI32", "LKI43", "LKI54", "LKU32", "LKU43", "RVX32", "RVX43", "RVX54", "OPN81", "OPN91", "OPN92", "OPN93", "WMN82", "MRN83", "MRN93", "MRU94"]
const LS_PVREF9 = Tuple{String,String,String}[("APN80", "904", "APN80"), ("APN81", "904", "APN81"), ("APN90", "904", "APN90"), ("APN91", "904", "APN91"), ("CTN11", "904", "CTN11"), ("CTN12", "904", "CTN12"), ("CTN24", "904", "CTN24"), ("CTN32", "904", "CTN32"), ("CTN42", "904", "CTN42"), ("CTU22", "904", "CTU22"), ("FDC12", "904", "FDC12"), ("FDC23", "904", "FDC23"), ("FDC24", "904", "FDC24"), ("FDC25", "904", "FDC25"), ("FDC34", "904", "FDC34"), ("FDN12", "904", "FDN12"), ("FDN22", "904", "FDN22"), ("FDN32", "904", "FDN32"), ("FDN33", "904", "FDN33"), ("FDN43", "904", "FDN43"), ("FFN57", "904", "FFN57"), ("FFN67", "904", "FFN67"), ("FPN62", "904", "FPN62"), ("FPN63", "904", "FPN63"), ("FPN71", "904", "FPN71"), ("FPN72", "904", "FPN72"), ("FPN73", "904", "FPN73"), ("FPN81", "904", "FPN81"), ("FPN82", "904", "FPN82"), ("FPS63", "904", "FPS63"), ("FPW63", "904", "FPW63"), ("LKI32", "904", "LKI32"), ("LKI43", "904", "LKI43"), ("LKI54", "904", "LKI54"), ("LKU32", "904", "LKU32"), ("LKU43", "904", "LKU43"), ("MHC26", "904", "MHC26"), ("MHC36", "904", "MHC36"), ("MHC37", "904", "MHC37"), ("MHC47", "904", "MHC47"), ("MHN35", "904", "MHN35"), ("MHN44", "904", "MHN44"), ("MHN45", "904", "MHN45"), ("MHN46", "904", "MHN46"), ("MHN47", "904", "MHN47"), ("MRN83", "904", "MRN83"), ("MRN93", "904", "MRN93"), ("MRU94", "904", "MRU94"), ("OPN81", "904", "OPN81"), ("OPN91", "904", "OPN91"), ("OPN92", "904", "OPN92"), ("OPN93", "904", "OPN93"), ("RON12", "904", "RON12"), ("RON23", "904", "RON23"), ("RVX32", "904", "RVX32"), ("RVX43", "904", "RVX43"), ("RVX54", "904", "RVX54"), ("WFN53", "904", "WFN53"), ("WFN55", "904", "WFN55"), ("WFN64", "904", "WFN64"), ("WFS57", "904", "WFS57"), ("WFW54", "904", "WFW54"), ("WMN82", "904", "WMN82"), ("10", "903", ""), ("11", "903", ""), ("110", "903", ""), ("120", "903", ""), ("21", "903", ""), ("211", "903", ""), ("212", "903", ""), ("22", "903", ""), ("221", "903", ""), ("222", "903", ""), ("23", "903", ""), ("230", "903", ""), ("24", "903", ""), ("31", "903", ""), ("311", "903", ""), ("312", "903", ""), ("32", "903", ""), ("321", "903", ""), ("322", "903", ""), ("33", "903", ""), ("330", "903", ""), ("34", "903", ""), ("411", "903", ""), ("412", "903", ""), ("62", "902", ""), ("63", "902", ""), ("64", "902", ""), ("72", "902", ""), ("73", "902", ""), ("74", "902", ""), ("80", "902", ""), ("81", "902", ""), ("82", "902", ""), ("AOC", "901", ""), ("AQVac", "901", ""), ("AQVib", "901", ""), ("ATD", "901", ""), ("AVO", "901", ""), ("DSH", "902", ""), ("FE", "901", ""), ("FI", "901", ""), ("FMC", "901", ""), ("HPM", "902", ""), ("HPM42", "902", ""), ("HPM43", "902", ""), ("HRM", "902", ""), ("HRM47", "902", ""), ("MSH", "902", ""), ("MSH37", "902", ""), ("OWP", "902", ""), ("OWP1", "902", ""), ("PCS", "901", ""), ("PO", "901", ""), ("PVC", "901", "FDC23"), ("PVD", "901", "FDC12"), ("QAE", "901", ""), ("TAM", "901", ""), ("TM", "901", ""), ("TMC", "901", ""), ("TMV", "901", ""), ("TTL", "901", ""), ("TTM", "901", ""), ("TTP", "901", ""), ("TTS", "901", "")]

"""
    habtyp_errors!(s::StandState{LakeStates}, pv, cpvref, kodtyp)

ls/habtyp.f's ERRGRO branches. With a PV reference code PVREF9 crosswalks (PV_CODE, PVREF) → HABPVR (FVS34 when both
codes are known but the pair is not / maps to blank, FVS33/FVS32 when either is unknown — those set LPVXXX and skip
FVS14); then HBDECD against LSNPC, then the sequence number IHB = IFIX(ARRAY2) ∈ 1..63; anything else is FVS14. The
PROCESS-time default call (initre.f:384-387, KODTYP never set — every FIA stand without a PV_CODE) has a blank code
⇒ FVS14, as live FVSls writes for every LS tiered stand.
"""
function habtyp_errors!(s::StandState{LakeStates}, pv::AbstractString, cpvref::AbstractString, kodtyp::Integer)
    k = String(strip(pv)); ihb = Int(kodtyp)
    if !isempty(cpvref)
        kard2 = ""; lc = false; lr = false
        for (code, ref, hab) in LS_PVREF9
            if code == k && ref == cpvref
                kard2 = hab; lc = true; lr = true; break
            end
            code == k && (lc = true)
            ref == cpvref && (lr = true)
        end
        if lc && lr && isempty(kard2)
            errgro!(s, 34); return nothing
        elseif !lc && !lr
            errgro!(s, 33); errgro!(s, 32); return nothing
        elseif lc && !lr
            errgro!(s, 32); return nothing
        elseif !lc && lr
            errgro!(s, 33); return nothing
        end
        k = kard2; ihb = 0                                  # PVREF9 zeroes ARRAY2
    end
    findfirst(==(k), LS_LSNPC) !== nothing && return nothing  # HBDECD
    (1 <= ihb <= length(LS_LSNPC)) && return nothing
    errgro!(s, 14)
    return nothing
end

