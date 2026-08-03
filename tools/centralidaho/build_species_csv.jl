# Build data/centralidaho/species_coefficients.csv (chunk 1) to UNBLOCK cit01 loading.
# Non-DG columns borrowed from IE template (FIA-matched, fallback per species) — REAL CI values
# injected: bark BARK1/BARK2 (ci/bratio.f POWER model DIB=BARK1·D^BARK2), SITELO/SITEHI
# (ci/siterange.f), SIGMAR dg_resid_sd (ci/blkdat.f). To be refined chunk-by-chunk vs FVSci_clean.
const IECSV = "data/inlandempire/species_coefficients.csv"
hdr = split(strip(readlines(IECSV)[1]), ',')
ierows = Dict{Int,Vector{String}}()
for ln in readlines(IECSV)[2:end]
    f = String.(split(strip(ln), ',')); fia = tryparse(Int, strip(f[3])); fia !== nothing && (ierows[fia] = f)
end
# CI species (alpha, fia, plants, IE-template-FIA). Fallbacks: WJ064→66(RM juniper), MC475→998, CW747→746(AS).
ci_sp = [("WP","119","PIMO3",119),("WL","073","LAOC",73),("DF","202","PSME",202),("GF","017","ABGR",17),
         ("WH","263","TSHE",263),("RC","242","THPL",242),("LP","108","PICO",108),("ES","093","PIEN",93),
         ("AF","019","ABLA",19),("PP","122","PIPO",122),("WB","101","PIAL",101),("PY","231","TABR2",231),
         ("AS","746","POTR5",746),("WJ","064","JUOC",66),("MC","475","CELE3",998),("LM","113","PIFL2",113),
         ("CW","747","POBAT",746),("OS","299","2TN",299),("OH","998","2TB",998)]
# REAL CI per-species constants (ci/bratio.f BARK1/BARK2; ci/siterange.f SITELO/SITEHI; ci/blkdat.f SIGMAR):
const CI_BARK1 = Float32[.859045,.900000,.903563,.904973,.903563,.837291,.900000,.900000,.903563,.809427,
                         .969,.969,.950,0.0,0.9,.969,.892,.900000,.892]
const CI_BARK2 = Float32[0.0,0.0,.989388,0.0,.989388,0.0,0.0,0.0,.989388,1.016866,
                         0.0,0.0,0.0,0.0,0.0,0.0,-0.086,0.0,-0.086]
const CI_SITELO = Float32[20,50,30,50,6,29,20,40,40,40,25,25,30,5,5,25,30,30,30]
const CI_SITEHI = Float32[80,110,70,110,203,152,100,100,90,80,50,50,70,15,15,50,120,70,120]
const CI_SIGMAR = Float32[0.230,0.206,0.267,0.260,0.260,0.206,0.203,0.232,0.245,0.230,
                          0.4671,0.4671,0.3750,0.2,0.5357,0.4671,0.2000,0.203,0.2000]
col = Dict(String(h) => i for (i, h) in enumerate(hdr))
open("data/centralidaho/species_coefficients.csv", "w") do io
    println(io, join(hdr, ","))
    for (idx, (al, fi, pl, tf)) in enumerate(ci_sp)
        tmpl = get(ierows, tf, ierows[202]); row = copy(tmpl)
        row[col["species_index"]] = string(idx); row[col["code_alpha"]] = al
        row[col["code_fia"]] = fi; row[col["code_plants"]] = pl
        row[col["bark1"]] = string(CI_BARK1[idx]); row[col["bark2"]] = string(CI_BARK2[idx])
        row[col["bark_imap"]] = "4"   # POWER-model bark (DIB=BARK1·D^BARK2) — same as BM; ci_bratio handles it
        row[col["dg_resid_sd"]] = string(CI_SIGMAR[idx])
        row[col["site_lo"]] = string(CI_SITELO[idx]); row[col["site_hi"]] = string(CI_SITEHI[idx])
        println(io, join(row, ","))
    end
end
println("wrote data/centralidaho/species_coefficients.csv (19 rows, IE-template non-DG cols)")
