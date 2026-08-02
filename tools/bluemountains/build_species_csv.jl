# Build a full-schema BM species_coefficients.csv to UNBLOCK stand loading (chunk-0 scaffold), mirroring
# tools/utah/build_species_csv.jl. Non-DG columns borrowed from CR templates (matched by FIA, fallback per
# species) — to be REPLACED with real BM values in later chunks. REAL injected: bark BARK1/BARK2 (bm/bratio.f
# POWER model DIB=BARK1·D^BARK2), SITELO/SITEHI (bm/siterange.f), SIGMAR dg_resid_sd (bm/blkdat.f).
const CRCSV = "data/centralrockies/species_coefficients.csv"
hdr = split(strip(readlines(CRCSV)[1]), ',')
crrows = Dict{Int,Vector{String}}()
for ln in readlines(CRCSV)[2:end]
    f = String.(split(strip(ln), ',')); fia = tryparse(Int, strip(f[3])); fia !== nothing && (crrows[fia] = f)
end
# BM species (alpha, fia, plants, CR-template-FIA). Fallbacks: WP119→108, WJ064→66, PY231/YC042→998, CW747→746.
bm = [("WP","119","PIMO3",108),("WL","073","LAOC",73),("DF","202","PSME",202),("GF","017","ABGR",17),
      ("MH","264","TSME",264),("WJ","064","JUOC",66),("LP","108","PICO",108),("ES","093","PIEN",93),
      ("AF","019","ABLA",19),("PP","122","PIPO",122),("WB","101","PIAL",101),("LM","113","PIFL2",113),
      ("PY","231","TABR2",998),("YC","042","CANO9",998),("AS","746","POTR5",746),("CW","747","POBAT",746),
      ("OS","299","2TN",299),("OH","998","2TB",998)]
# REAL BM per-species constants (bm/bratio.f BARK1/BARK2 POWER; bm/siterange.f SITELO/SITEHI; bm/blkdat.f SIGMAR):
const BM_BARK1 = Float32[0.859045,0.859045,0.903563,0.904973,0.903563,0.0,0.9,0.9,0.904973,0.809427,0.969,0.9625,0.933290,0.837291,0.950,0.075256,0.809427,0.9000]
const BM_BARK2 = Float32[1.0,1.0,0.989388,1.0,0.989388,0.0,1.0,1.0,1.0,1.016866,0.0,-0.1141,1.0,1.0,0.0,0.949670,1.016866,1.0]
const BM_SIGMAR = Float32[0.5670,0.3383,0.2580,0.2548,0.5571,0.34663,0.2820,0.3348,0.3249,0.2745,0.34663,0.46710,0.4842,0.3931,0.34663,0.5357,0.2745,0.5357]
const BM_SITELO = Float32[20,50,50,50,15,5,30,40,50,70,20,20,5,50,30,10,70,5]
const BM_SITEHI = Float32[80,110,110,110,30,40,70,120,150,140,65,50,75,110,66,191,140,125]
ci = Dict(String(h) => i for (i, h) in enumerate(hdr))
open("data/bluemountains/species_coefficients.csv", "w") do io
    println(io, join(hdr, ","))
    for (idx, (al, fi, pl, tf)) in enumerate(bm)
        tmpl = get(crrows, tf, crrows[202]); row = copy(tmpl)
        row[ci["species_index"]] = string(idx); row[ci["code_alpha"]] = al
        row[ci["code_fia"]] = fi; row[ci["code_plants"]] = pl
        row[ci["bark1"]] = string(BM_BARK1[idx]); row[ci["bark2"]] = string(BM_BARK2[idx])
        row[ci["bark_imap"]] = "4"   # BM power-model bark (DIB=BARK1·D^BARK2) — bm_bratio handles it; 4=sentinel
        row[ci["dg_resid_sd"]] = string(BM_SIGMAR[idx])
        row[ci["site_lo"]] = string(BM_SITELO[idx]); row[ci["site_hi"]] = string(BM_SITEHI[idx])
        println(io, join(row, ","))
    end
end
println("wrote data/bluemountains/species_coefficients.csv (18 rows, CR-placeholder non-DG cols)")
