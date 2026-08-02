# Build a full-schema EM species_coefficients.csv to UNBLOCK stand loading (the original was never
# committed and is lost). PLACEHOLDER SCAFFOLD, mirroring .sweep_work/build_kt_species_csv.jl: non-DG
# columns (site/bark/height/mort/vol/htdbh) are borrowed from CR's western rows (matched by FIA code;
# fallback listed per species) — to be REPLACED with real EM values (em/blkdat.f + bratio/htdbh/varvol)
# for bit-exact .sum validation. DG per-tree WK2 is UNAFFECTED (reads dg_coefficients.jl + dg_*.csv);
# the ported non-conifer DG + height FORMS can be smoke-tested against live FVSem with this in place.
const CRCSV = "data/centralrockies/species_coefficients.csv"
hdr = split(strip(readlines(CRCSV)[1]), ',')
crrows = Dict{String,Vector{String}}()
for ln in readlines(CRCSV)[2:end]
    f = String.(split(strip(ln), ','))
    crrows[f[3]] = f
end
# EM species in order: (alpha, fia, plants, CR-template-FIA-to-borrow). Fallbacks for FIA absent in CR:
# LL 072 → WL 73 (larch); GA 544 → OH 998; CW 747 → AS 746; BA 741 → PW 745 (poplars).
em = [("WB","101","PIAL","101"),("WL","073","LAOC","73"),("DF","202","PSME","202"),
      ("LM","113","PIFL2","113"),("LL","072","LALY","73"),("RM","066","JUSC2","66"),
      ("LP","108","PICO","108"),("ES","093","PIEN","93"),("AF","019","ABLA","19"),
      ("PP","122","PIPO","122"),("GA","544","FRPE","998"),("AS","746","POTR5","746"),
      ("CW","747","POBAT","746"),("BA","741","POBA2","745"),("PW","745","PODEM","745"),
      ("NC","749","POAN3","749"),("PB","375","BEPA","375"),("OS","299","2TN","299"),
      ("OH","998","2TB","998")]
ci = Dict(String(h) => i for (i, h) in enumerate(hdr))
open("data/easternmontana/species_coefficients.csv", "w") do io
    println(io, join(hdr, ","))
    for (idx, (al, fi, pl, tf)) in enumerate(em)
        tmpl = get(crrows, tf, crrows["202"])
        row = copy(tmpl)
        row[ci["species_index"]] = string(idx)
        row[ci["code_alpha"]]    = al
        row[ci["code_fia"]]      = fi
        row[ci["code_plants"]]   = pl
        println(io, join(row, ","))
    end
end
println("wrote data/easternmontana/species_coefficients.csv (19 rows, CR-placeholder non-DG cols)")
