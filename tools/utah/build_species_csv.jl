# Build a full-schema UT species_coefficients.csv to UNBLOCK stand loading (chunk-0 scaffold).
# PLACEHOLDER, mirroring tools/easternmontana/build_species_csv.jl: non-DG columns (dbh_max,
# st_*, ht*, mort, vol, htdbh) are borrowed from CR's western rows (matched by FIA, fallback
# listed per species) — to be REPLACED with real UT values (ut/blkdat.f + bratio/htdbh/varvol)
# for bit-exact .sum validation in later chunks. The REAL UT bark (ut/bratio.f), regent site
# range (ut/siterange.f SLO/SHI), and DG residual SD (ut/blkdat.f SIGMAR) are injected here so
# the ported DG/height/regent FORMS use the correct per-species constants. DG per-tree WK2 is
# UNAFFECTED (reads dg_coefficients.jl + dg_*.csv, ported in the DG chunk).
const CRCSV = "data/centralrockies/species_coefficients.csv"
hdr = split(strip(readlines(CRCSV)[1]), ',')
crrows = Dict{Int,Vector{String}}()
for ln in readlines(CRCSV)[2:end]
    f = String.(split(strip(ln), ','))
    fia = tryparse(Int, strip(f[3]))
    fia !== nothing && (crrows[fia] = f)
end
# UT species in order: (alpha, fia, plants, CR-template-FIA). Fallbacks for FIA absent in CR:
# WJ 064 → RM juniper 66; GB 142 → 143; FC 748 → NC 749; MC 475 → mtn-shrub 810; BI 322/BE 313 → OH 998.
ut = [("WB","101","PIAL",101),("LM","113","PIFL2",113),("DF","202","PSME",202),
      ("WF","015","ABCO",15),("BS","096","PIPU",96),("AS","746","POTR5",746),
      ("LP","108","PICO",108),("ES","093","PIEN",93),("AF","019","ABLA",19),
      ("PP","122","PIPO",122),("PI","106","PIED",106),("WJ","064","JUOC",66),
      ("GO","814","QUGA",814),("PM","133","PIMO",133),("RM","066","JUSC2",66),
      ("UJ","065","JUOS",65),("GB","142","PILO",143),("NC","749","POAN3",749),
      ("FC","748","POFR2",749),("MC","475","CELE3",810),("BI","322","ACGR3",998),
      ("BE","313","ACNE2",998),("OS","299","2TN",299),("OH","998","2TB",998)]
# REAL UT per-species constants (ut/bratio.f BARK1/BARK2/IMAP; ut/blkdat.f SIGMAR; ut/siterange.f SLO/SHI):
const UT_BARK1 = Float32[0.9625,0.9625,0.867,0.890,0.9502,0.950,0.9625,0.9502,0.890,0.8967,0.0,0.0,0.93789,0.0,0.0,0.0,0.9625,0.892,0.892,0.9,0.94782,0.892,0.9625,0.93789]
const UT_BARK2 = Float32[-0.1141,-0.1141,0.0,0.0,-0.2528,0.0,-0.1141,-0.2528,0.0,-0.4448,0.0,0.0,-0.24096,0.0,0.0,0.0,-0.1141,-0.086,-0.086,0.0,0.0836,-0.086,-0.1141,-0.24096]
const UT_IMAP  = Int[3,3,2,2,3,2,3,3,2,3,1,1,3,1,1,1,3,3,3,2,3,3,3,3]
const UT_SIGMAR = Float32[0.46710,0.46710,0.34418,0.24060,0.35168,0.37500,0.28860,0.35168,0.28005,0.27338,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5357,0.5107,0.2,0.46710,0.2]
const UT_SLO = Float32[15,15,30,20,20,30,27,20,27,23,5,5,5,5,5,5,20,30,30,5,5,15,15,5]
const UT_SHI = Float32[32,32,70,50,58,70,53,65,58,50,20,15,20,20,15,15,60,120,90,15,30,40,32,20]
ci = Dict(String(h) => i for (i, h) in enumerate(hdr))
open("data/utah/species_coefficients.csv", "w") do io
    println(io, join(hdr, ","))
    for (idx, (al, fi, pl, tf)) in enumerate(ut)
        tmpl = get(crrows, tf, crrows[202])
        row = copy(tmpl)
        row[ci["species_index"]] = string(idx)
        row[ci["code_alpha"]]    = al
        row[ci["code_fia"]]      = fi
        row[ci["code_plants"]]   = pl
        row[ci["bark1"]]     = string(UT_BARK1[idx])   # REAL UT bark (rest of the row still CR placeholder)
        row[ci["bark2"]]     = string(UT_BARK2[idx])
        row[ci["bark_imap"]] = string(UT_IMAP[idx])
        row[ci["dg_resid_sd"]] = string(UT_SIGMAR[idx])
        row[ci["site_lo"]] = string(UT_SLO[idx]); row[ci["site_hi"]] = string(UT_SHI[idx])
        println(io, join(row, ","))
    end
end
println("wrote data/utah/species_coefficients.csv (24 rows, CR-placeholder non-DG cols)")
